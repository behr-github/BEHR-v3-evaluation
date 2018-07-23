function [  ] = iterate_vertical_profile_plot( iterable_field, iterable_ranges, start_date, end_date, binwidth, data_field_1, varargin )
%Creates successive plots that show vertical profiles within different criteria specified by "iterable_field" and "iterable range"
%   This function repeatedly calls "avg_vertical_profiles" and
%   "plot_vertical_profiles."  The parameter values "other_crit" and
%   "other_crit_ranges" allow you to specify restrictions on the vertical
%   profiles in the same way as "avg_vertical_profiles."  The
%   "iterable_field" and "iterable_ranges" similarly specify restrictions,
%   but allow you to iterate over several values for the field specified.
%
%   iterable_field (required): A string specifying the field in the
%   aircraft merge files that you wish to make several various plots for.
%
%   iterable_ranges (required): A n x 1 or n x 2 matrix specifying each
%   value or range to restrict 'iterable_field' to. Each row represents a
%   unique value or range, see examples below for clarification. Note that
%   the fields discoveraqSiteFlag1sec and ProfileSequenceNum can be called
%   with the shorthand strings 'siteflag' and 'profnum' respectively. There
%   are also special fields, 'starttime-lt' and 'starttime-utc' which
%   collect profiles with starting times within the range specified.
%
%   start_date (required): The beginning date of the period to average
%   over, passed as a string.
%
%   end_date (required): The end date of the period to average over, passed
%   as a string.
%
%   binwidth (required): The width of altitude bins to group the data in,
%   in kilometers.
%
%   data_field_1 (required): The data to plot for the vertical profile.
%   Must match a field in the Merge file.
%
%   data_field_2 (optional): This optional 7th argument is used to plot a
%   second vertical profile on the same figure as the first, with its own
%   x-axis.
%
%   other_crit (parameter): A string or cell array of strings specifying
%   addition criteria by which to restrict the vertical profiles.  See
%   notes for iterable_field for shorthand field names.
%
%   other_crit_ranges (parameter): A 1 x 2 matrix or cell array of said
%   matrices that correspond to the fields given in other_crit and describe
%   the ranges to which to restrict the other values.
%
%   Examples:
%       (a) Plot NO2 profiles for all of Baltimore/DC campaign (Jul 2011),
%       dividing profile data by local solar time in 2-hr increments:
%           iterate_vertical_profiles('LOCAL_SUN_TIME',[7 9; 9 11; 11 13;
%           13 15], '07/01/2011','07/31/2011',0.3,'NO2_LIF')
%
%       (b) Plot NO2 profiles against potential temperature for each site
%       in the Baltimore/DC campaign, restricted to morning profiles.
%           iterate_vertical_profiles('siteflag',[1;2;3;4;5;6],
%           '07/01/2011','07/31/2011','NO2_LIF','THETA','other_crit','LOCAL_SUN_TIME',
%           [0 12])
%
%       (c) Like (b), but also restricted CO mixing ratio to be < 100 ppbv.
%           iterate_vertical_profiles('('siteflag',[1;2;3;4;5;6],
%           '07/01/2011','07/31/2011','NO2_LIF','THETA','other_crit',{'LOCAL_SUN_TIME',
%           'Carbon_Monoxide_mixing_ratio'},{[0 12],[0 100]})
%
%   Josh Laughner <joshlaugh5@gmail.com> 28 May 2014

p = inputParser;
p.addRequired('iterable_field',@isstr);
p.addRequired('iterable_ranges',@ismatrix);
p.addRequired('start_date',@isstr);
p.addRequired('end_date',@isstr);
p.addRequired('binwidth',@isscalar);
p.addRequired('data_field_1',@isstr);
p.addOptional('data_field_2','',@isstr);
p.addParamValue('other_crit',{},@(x) (ischar(x) || iscell(x)));
p.addParamValue('other_crit_ranges',{},@(x) (isnumeric(x) || iscell(x)));

p.parse(iterable_field, iterable_ranges, start_date, end_date, binwidth, data_field_1, varargin{:});
pout = p.Results;

iterable_field = pout.iterable_field;
iterable_ranges = pout.iterable_ranges;
start_date = pout.start_date;
end_date = pout.end_date;
binwidth = pout.binwidth;
data_field_1 = pout.data_field_1;
data_field_2 = pout.data_field_2;
other_crit = pout.other_crit;
other_ranges = pout.other_crit_ranges;

% If the user input a row vector for the iterable_ranges, switch it to a
% column vector so that each row corresponds to a value.
if isrow(iterable_ranges); iterable_ranges = iterable_ranges'; end

% Ensure the other criteria names and their ranges are cell structures;
% this simplifies their implementation later.
if ~iscell(other_crit); other_crit = {other_crit}; end
if ~iscell(other_ranges); other_ranges = {other_ranges}; end

% Combine the iterated and other criteria
if ~isempty(other_crit)
    all_crit = [{iterable_field}, other_crit];
else
    all_crit = iterable_field;
end

% Parse the shorthand entries for profile number and site flag if we are
% considering aircraft data.  Allow the user to capitalize the first letter
all_crit = regexprep(all_crit,'[pP]rofnum','ProfileSequenceNum');
all_crit = regexprep(all_crit,'[sS]iteflag','discoveraqSiteFlag1sec');

s = size(iterable_ranges);

for a = 1:s(1)
    % Combine the current range for the iterated value with the additional
    % criteria ranges.
    this_range = iterable_ranges(a,:);
    if ~isempty(other_ranges)
        curr_ranges = [{this_range}, other_ranges];
    else
        curr_ranges = {this_range};
    end
    
    % If a second data field was passed, bin both and plot both.
    % Otherwise, plot a single profile.
    if ~isempty(data_field_2)
        [val, mid, err] = avg_vertical_profiles(data_field_1,'Aircraft',binwidth,start_date,end_date,all_crit,curr_ranges);
        [val2, mid2, err2] = avg_vertical_profiles(data_field_2,'Aircraft',binwidth,start_date,end_date,all_crit,curr_ranges);
        plot_vertical_profile_bins(val, mid, err, val2, mid2, err2);
    else
        [val, mid, err] = avg_vertical_profiles(data_field_1,'Aircraft',binwidth,start_date,end_date,all_crit,curr_ranges);
        plot_vertical_profile_bins(val, mid, err);
    end
end

end

