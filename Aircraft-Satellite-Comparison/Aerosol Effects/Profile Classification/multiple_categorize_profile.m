function [ OutStruct, Diagnostics ] = multiple_categorize_profile( varargin )
%multiple_categorize_profile Run categorize_aerosol_profile 3 times and take the most common result.
%   The general idea of categorizing the relative position of NO2 and
%   aerosol profiles by comparing the altitude below which x% of the column
%   can be found works well, but can sometimes ignore the structure of the
%   profile at lower altitudes.  By checking the categorization for
%   multiple critical fractions and taking the most common result from
%   them, this should overcome that problem.
%
%   This requires either 1, 2 or 4 inputs. In the 1 or 2 argument forms,
%   pass a campaign name valid for merge_field_names as the first input (as
%   a string). Pass the string 'green' as the second argument to force this
%   function to use the original green wavelength aerosol data, instead of
%   the blue data provided by Lee Thornhill (which better reflects the
%   aerosol at wavelengths important to NO2 retrievals). or pass an NO2
%   profile, its corresponding altitudes, an aerosol extinction profile,
%   and its corresponding altitudes as the four inputs. The first mode will
%   categorize all profiles in that campaign; the second will categorize
%   the given profile.
%
%   Josh Laughner <joshlaugh5@gmail.com> 7 Apr 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%
%%%%% INPUT %%%%%
%%%%%%%%%%%%%%%%%

narginchk(1,4);
if ischar(varargin{1})
    campaign_name = varargin{1};
    wavelength = 'blue';
    if numel(varargin) > 1
        if ischar(varargin{2}) && strcmpi(varargin{2},'green')
            wavelength = 'green';
        else
            E.badinput('The second input, if present and if the first is a campaign name, must be the string "green"');
        end
    end
elseif nargin == 4 && all(iscellcontents(varargin,'isnumeric'))
    no2_in = varargin{1};
    no2_alt_in = varargin{2};
    aer_in = varargin{3};
    aer_alt_in = varargin{4};
    if ~all(size(no2_in) == size(no2_alt_in))
        E.badinput('no2_in and no2_alt_in (first two arguments) must be the same size');
    elseif ~all(size(aer_in) == size(aer_alt_in))
        E.badinput('aer_in and aer_alt_in (third and fourth arguments) must be the same size');
    end
else
    E.badinput('multiple_categorize_profile requires either a valid campaign name as a string or 4 inputs: NO2 profile & altitudes, aerosol extinction and altitudes');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL CATEGORIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

crit_fracs = [0.5, 0.75, 0.9];
cat_structs = cell(1,numel(crit_fracs));
diag_structs = cell(1,numel(crit_fracs));
for a=1:numel(crit_fracs)
    if nargin < 4
        [cat_structs{a}, diag_structs{a}] = categorize_aerosol_profile('campaign_name',campaign_name,'crit_frac',crit_fracs(a),'wavelength',wavelength,'DEBUG_LEVEL',1);
        out_fxn = @add_cat;
    else
        [cat_structs{a}, diag_structs{a}] = categorize_aerosol_profile(no2_in, no2_alt_in, aer_in, aer_alt_in,'crit_frac',crit_fracs(a),'DEBUG_LEVEL',1);
        out_fxn = @add_cat_simple;
    end
end

if nargin < 4
    % This is the case where we are dealing with an entire campaign's worth
    % of data.
    % Since tabulate_profile_categories already nicely arranges everything but
    % profile we'll go ahead and use that.
    [cat_tab, crit_tab] = tabulate_profile_categories(cat_structs{:});
    
    cat_data = table2cell(cat_tab);
    % Column 1 = profile #, Column 2 = date string, Column 3+ = category
else
    % If, on the other hand, we're just looking at a specific profile, then
    % we've already got a cell array of categorizations.
    cat_data = [cell(1,2), cat_structs];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% META-CATEGORIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prep the output structure: it will have a similar form to the usual
% categorize_aerosol_prof output structures, although it will lose some of
% the extra fields beyond profile ID and date.

if nargin < 4
    OutStruct.Column_Fraction_Criterion = crit_tab{'Column_Fraction_Criterion',:};
    OutStruct.Magnitude_Criterion_Type = crit_tab{'Magnitude_Criterion_Type',:};
    OutStruct.Magnitude_Critical_Value = crit_tab{'Magnitude_Critical_Value',:};
    OutStruct.dz = crit_tab{'dz',:};
end

OutStruct.CoincidentLow = [];
OutStruct.CoincidentLowDates = {};
OutStruct.CoincidentHigh = [];
OutStruct.CoincidentHighDates = {};
OutStruct.AerosolAboveLow = [];
OutStruct.AerosolAboveLowDates = {};
OutStruct.AerosolAboveHigh = [];
OutStruct.AerosolAboveHighDates = {};
OutStruct.NO2AboveLow = [];
OutStruct.NO2AboveLowDates = {};
OutStruct.NO2AboveHigh = [];
OutStruct.NO2AboveHighDates = {};

% Rules: if two or three categorizations agree on the location
% (coinc/aer above/no2 above), that categorization is chosen.

diag_reasons = zeros(size(cat_data,1),2);

for a=1:size(cat_data,1)
    cats = cat_data(a,:);
    % Look for the categories - allow there to be a space in between
    % aerosol & above or no2 & above.
    xx_coinc = ~iscellcontents(regexpi(cats(3:end),'Coincident'),'isempty');
    xx_aer = ~iscellcontents(regexpi(cats(3:end),'Aerosol ?Above'),'isempty');
    xx_no2 = ~iscellcontents(regexpi(cats(3:end),'NO2 ?Above'),'isempty');
    
    if sum(xx_coinc) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_coinc);
    elseif sum(xx_aer) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_aer);
    elseif sum(xx_no2) > 1
        OutStruct = out_fxn(OutStruct, cats, xx_no2);
    else
        diag_reasons(a,1) = cat_data{a,1}(1); % This is the profile number. We explicitly take the first element to handle UTC ranges.
        diag_reasons(a,2) = meta_diag_reason(xx_coinc, xx_aer, xx_no2);
    end
    
    
end

% Make the Diagnostic output structure. It will contain each critical
% fraction's diagnostics plus the meta-categorization diagnostics.
Diagnostics = struct;
for a=1:numel(crit_fracs)
    fn = sprintf('CritFrac%d',crit_fracs(a)*100);
    Diagnostics.(fn) = diag_structs{a};
end
Diagnostics.IndividualTables.categories = cat_tab;
Diagnostics.IndividualTables.criteria = crit_tab;
Diagnostics.MultipleCat.reasons = {'101 - all runs of categorize_aerosol_profile failed to return a categorization';...
                                   '102 - all but one run of categorize_aerosol_profile failed to return a categorization';...
                                   '200 - no majority category';...
                                   '999 - other'};

dd = diag_reasons(:,2) ~= 0;                               
Diagnostics.MultipleCat.profiles = diag_reasons(dd,:);


end

function OutStruct = add_cat(OutStruct, cat_line, xx)
    % find the categorization name - we want to read it from the table to
    % capture whether it is "low" or "high" AOD
    E = JLLErrors;
    cats = cat_line(3:end);
    cat_name = unique(cats(xx));
    if numel(cat_name) > 1 % both names should be the same. if not, don't continue.
        E.badvartype(cat_name,'1x1 cell array');
    end
    
    catdate_name = strcat(cat_name,'Dates');
    catprof = cat_line{1};
    catdate = cat_line(2); % we want the actual cell, not its contents, to concatenate
    
    OutStruct.(cat_name{1}) = cat(1,OutStruct.(cat_name{1}),catprof);
    OutStruct.(catdate_name{1}) = cat(1,OutStruct.(catdate_name{1}),catdate);
end

function cat_name = add_cat_simple(~, cat_line, xx)
    E = JLLErrors;
    cats = cat_line(3:end);
    cat_name = unique(cats(xx));
    if numel(cat_name) > 1 % both names should be the same. if not, don't continue.
        E.badvartype(cat_name,'1x1 cell array');
    end
end

function reason = meta_diag_reason(xx_coinc, xx_aer, xx_no2)
    % Identify the reason a profile was not meta-categorized. Codes are:
    %   100 - All but one critical fraction (102), or all critical
    %   fractions (101) are uncategorized.
    %   200 - There are two or more categorized fractions, but they all
    %   disagree
    sum_all = sum(xx_coinc) + sum(xx_aer) + sum(xx_no2);
    if sum_all == 0
        reason = 101;
    elseif sum_all == 1
        reason = 102;
    elseif sum(xx_coinc) < 2 && sum(xx_aer) < 2 && sum(xx_no2) < 2
        reason = 200;
    else
        reason = 999;
    end
end

