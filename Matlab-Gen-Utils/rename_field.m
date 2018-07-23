function [ S ] = rename_field( S, oldfieldname, newfieldname )
%rename_field Renames a field in the structure S

narginchk(3,3);
if nargout < 1;
    warning('rename_field must be assigned back into structure');
end

% Find the position of the old field
fields = fieldnames(S);
xx = find(strcmp(fields,oldfieldname));

% Create a new field that copies the old field
S.(newfieldname) = S.(oldfieldname);

% Delete the old field 
S = rmfield(S,oldfieldname);

% Put the new field (which was appended to the end of the structure) in the
% same place as the old one
n = numel(fields);
perm = [1:(xx-1),n,xx:(n-1)];
S = orderfields(S,perm);

end

