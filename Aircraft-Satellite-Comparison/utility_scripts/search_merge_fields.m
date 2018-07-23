function [ fields_out ] = search_merge_fields( Merge, field_string, varargin )
%search_merge_fields(Merge, field_string, [case sensitive]) Returns all field names containing field_string
%
%   Finds all fields in Merge.Data that contain the string field_string.
%   By default it is case insensitive, pass a 1 as the optional third
%   argument to require it match case as well.  The field string can be a
%   regular expression if desired.

fields = fieldnames(Merge.Data);
if nargin > 2 && varargin{1}
    matches = regexp(fields,field_string);
else
    matches = regexpi(fields,field_string);
end

fields_out = cell(numel(matches),1);
N=0;
for a=1:numel(matches)
    if ~isempty(matches{a})
        N=N+1;
        fields_out{N} = fields{a};
    end
end
fields_out = fields_out(1:N);
end