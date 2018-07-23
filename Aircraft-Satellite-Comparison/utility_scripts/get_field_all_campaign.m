function [ data, dates ] = get_field_all_campaign( campaign_name, field )
%get_field_all_campaign Get all values of "field" for the given campaign
%   Iterates through all merge fields from the given campaign and
%   concatenates all values from the field passed as the second argument.
%   Requires a campaign name understood by merge_field_names as the first
%   argument and a valid field name from that campaign's Merge.Data
%   structure. Returns a vector of the data with fill values replaced with
%   NaNs; optionally also returns a cell array of dates in yyyy-mm-dd
%   format.

E = JLLErrors;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(campaign_name)
    E.badinput('campaign_name must be a string');
elseif ~ischar(field)
    E.badinput('field must be a string');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,merge_dir] = merge_field_names(campaign_name);
F = dir(fullfile(merge_dir,'*.mat'));
n = numel(F);

data = [];
dates = {};

for a=1:n
    load(fullfile(merge_dir, F(a).name),'Merge'); % loads the variable "Merge" into the workspace
    this_date = Merge.metadata.date;
    
    % Give the user intelligible errors: if there's a problem loading the
    % field, consider if it's the first file or not. If it is, assume that
    % that means the user entered the wrong field. If not, then somehow the
    % Merge files loaded have inconsistent fields.
    try
        this_days_data = remove_merge_fills(Merge,field);
    catch err
        if strcmp(err.identifier, 'MATLAB:nonExistentField')
            if a==1
                E.badinput('%s is not a valid field for the campaign %s',field,campaign_name);
            else
                E.callError('inconsistentFields','%s was found in previous Merge files, but is not present for the file on %s',field,this_date);
            end
        else
            rethrow(err)
        end
    end
    
    data = cat(2,data,this_days_data);
    if nargout > 1
        date_cell = repmat({this_date},size(data));
        dates = cat(2,date_cell);
    end
end

end

