function [ prof_table, crit_table ] = tabulate_profile_categories( varargin )
%tabulate_profile_categories Returns profiles sorted by date and number
%   To compare the changes in categorization with different criteria, this
%   will take any number of structures returns from
%   categorize_aerosol_profile in and sort the profiles by date and number,
%   appending the categorization given in each structure.  It will return a
%   table with all of this information. The second table will contain a
%   summary of the criteria set for each experiment.
%
%   If you pass a structure of user-defined categorizations, make it the
%   last one and make the final argument a scalar value > 0.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,Inf);

% Figure out if there is a user-defined categorization structure - it'll
% need handled specially
if ~isstruct(varargin{end}) && isscalar(varargin{end}) && varargin{end} > 0
    user_struct = varargin{end-1};
    varargin = varargin(1:end-2);
    user_struct_bool = true;
else 
    user_struct_bool = false;
end

% Check that only structures were passed
for a=1:numel(varargin)
    if ~isstruct(varargin{a})
        error(E.badinput('sort_profile_categories only accepts structures as inputs'));
    end
end


% The "Info" structure will contain important stats about the first
% structure passed, to check against following structure to ensure that
% they are consistent.

fns = fieldnames(varargin{1});
Info.fieldnames = fns;
xx = ~iscellcontents(regexp(fns,'(Coincident|AerosolAbove|NO2Above)'),'isempty');

Info.data_fields = fns(xx);
Info.info_fields = fns(~xx);
Info.profnum_fields = {'CoincidentLow','CoincidentHigh','AerosolAboveLow','AerosolAboveHigh','NO2AboveLow','NO2AboveHigh'};
Info.date_fields = strcat(Info.profnum_fields,'Dates');

if user_struct_bool
    % A user defined structure from categorize_aerosol_profiles will ONLY
    % have profnum and date fields, so we just need to separate those out
    user_fns = fieldnames(user_struct);
    xx = ~iscellcontents(regexpi(user_fns,'Dates'),'isempty');
    UserInfo.profnum_fields = user_fns(~xx);
    UserInfo.date_fields = user_fns(xx);
end

% Check that the second and later structures have the same fields as the
% first one. We need to use two ismember() calls, once to ensure that all
% the fields in the new structure are in the first, and once to ensure that
% all the fields in the first are in the new.
for a=2:numel(varargin)
    next_fns = fieldnames(varargin{a});
    if ~all(ismember(Info.profnum_fields,next_fns))
        pnfs = strjoin(Info.profnum_fields{:});
        E.badinput('All input structures must contain the fields %s',pnfs);
    end
    if ~all(ismember(Info.date_fields,next_fns))
        dfs = strjoin(Info.date_fields{:});
        E.badinput('All input structures must contain the fields %s',dfs);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Get all the profile numbers and dates from each structure, identify
% unique profiles, and sort them.  If a user created structure was passed,
% it has been removed from the varargin and saved as user_struct, so we can
% loop over the varargin as automatic structures and just handle the
% user_struct if needed.
Out.all_profnums = [];
Out.all_dates = {};
for a=1:numel(varargin)
    curr_struct = varargin{a};
    for f=1:numel(Info.profnum_fields)
        if ~isempty(curr_struct.(Info.profnum_fields{f}))
            Out.all_profnums = cat(1,Out.all_profnums,curr_struct.(Info.profnum_fields{f}));
            Out.all_dates = cat(1,Out.all_dates,curr_struct.(Info.date_fields{f}));
        end
    end
end
if user_struct_bool
    for f=1:numel(UserInfo.profnum_fields)
        if ~isempty(user_struct.(UserInfo.profnum_fields{f}))
            Out.all_profnums = cat(1,Out.all_profnums,user_struct.(UserInfo.profnum_fields{f}));
            Out.all_dates = cat(1,Out.all_dates,user_struct.(UserInfo.date_fields{f}));
        end
    end
end

% Next we want to arrange the lines of the table in ascending order and
% remove duplicates - the idea is that multiple categorizations will have
% some profiles in common (profiles that were able to be categorized with
% both sets of criteria) and some that will not overlap (those that failed
% to be categorized for some reason in one, but not the other). Once we
% have this list of unique profiles, later we will match them up with the
% categorizations from each run.

% We're using sortrows here because it will behave correctly for either
% profile numbers or UTC ranges.
[Out.all_profnums, sort_perm] = sortrows(Out.all_profnums);
Out.all_dates = Out.all_dates(sort_perm);
all_datenums = datenum(Out.all_dates);

keep_prof = true(size(Out.all_dates));

% This is where we actually decide figure out which profiles are duplicates
% and remove them. We use the all() statement to compare the start and end
% of the UTC ranges, but when dealing with profile numbers, there's only a
% single profile identifier to compare.
for b=2:size(Out.all_profnums,1)
    if all(Out.all_profnums(b,:) == Out.all_profnums(b-1,:)) && all_datenums(b) == all_datenums(b-1)
        keep_prof(b) = false;
    end
end

Out.all_profnums = Out.all_profnums(keep_prof,:);
Out.all_dates = Out.all_dates(keep_prof);
all_datenums = all_datenums(keep_prof);


% Now loop through each structure and make the table. First column is the
% dates, second is prof nums, third and on will be the categorization of
% each profile. Each field in the structure will become a column in the
% table, so we will have a field for each categorization carried out.

if user_struct_bool
    Out.UserDefined = cell(size(Out.all_profnums));
    for f=1:numel(UserInfo.profnum_fields)
        xx = ismember(Out.all_profnums,user_struct.(UserInfo.profnum_fields{f})) & ismember(all_datenums,datenum(user_struct.(UserInfo.date_fields{f})));
        Out.UserDefined(xx) = UserInfo.profnum_fields(f);
    end
end

for a=1:numel(varargin)
    curr_struct = varargin{a};
    curr_table_col = sprintf('InputStct%02d',a);
    Out.(curr_table_col) = cell(size(Out.all_dates));
    for f=1:numel(Info.profnum_fields)
        xx = all(ismember(Out.all_profnums,curr_struct.(Info.profnum_fields{f})),2) & ismember(all_datenums,datenum(curr_struct.(Info.date_fields{f})));
        Out.(curr_table_col)(xx) = Info.profnum_fields(f);
    end
    
    Criteria.(curr_table_col) = cell(numel(Info.info_fields),1);
    for f=1:numel(Info.info_fields)
        curr_info = curr_struct.(Info.info_fields{f});
        if ~iscell(curr_info)
            Criteria.(curr_table_col){f} = curr_info;
        else
            % Handle the possibility of a structure output from the
            % multiple categorization function that has multiple criteria:
            % in that case, the criteria will be cell arrays.
            if all(iscellcontents(curr_info,'isnumeric'))
               Criteria.(curr_table_col){f} = mat2str(cell2mat(curr_info));
            else
                u_curr_info = unique(curr_info);
                if numel(u_curr_info) == 1
                    Criteria.(curr_table_col){f} = u_curr_info{1};
                else
                    Criteria.(curr_table_col){f} = 'various';
                end
            end
        end
    end
end

% Convert both structures to tables, and in the case of the Criteria table,
% assign the row variable names.

prof_table = struct2table(Out);
crit_table = struct2table(Criteria,'RowNames',Info.info_fields);


end

