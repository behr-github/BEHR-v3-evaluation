function [ field_value ] = find_substruct_field(S, field_name)
%FIND_SUBSTRUCT_FIELD Finds a field arbitrarily deep in a structure
%   FIELD_VALUE = FIND_SUBSTRUCT_FIELD( S, FIELD_NAME ) Searchs S and any
%   substructures within S for fields with the name FIELD_NAME. S and any
%   substructures must, at present, be scalar. Once found, returns the
%   value of that field in FIELD_VALUE.
%
%   This function only returns the first match. If multiple substructures
%   contain the desired field, only the first one is returned. This
%   function does not continue checking and so currently offers no
%   guarantee that there were not other instances of that field.  It
%   searches in a quasi-depth first manner; if the requested field is
%   present in S, that value is returned. If not, it recurses into the
%   field field of S that is a structure. It checks all the fields in that
%   substructure, if still not found it recurses into the first structure
%   within that substructure.

E = JLLErrors;

if ~isstruct(S)
    E.badinput('S must be a structure');
elseif ~isscalar(S)
    E.notimplemented('FIND_SUBSTRUCT_FIELD not yet designed for non-scalar structures');
end

if ~ischar(field_name)
    E.badinput('FIELD_NAME must be a character array')
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

if isfield(S, field_name)
    field_value = S.(field_name);
    return
else
    fns = fieldnames(S);
    for i_fn = 1:numel(fns)
        if isstruct(S.(fns{i_fn}))
            try
                field_value = find_substruct_field(S.(fns{i_fn}), field_name);
                return
            catch err
                if ~strcmp(err.identifier, 'MATLAB:nonExistentField')
                    rethrow(err)
                end
            end
        end
    end
    
    error('MATLAB:nonExistentField', 'Could not find the field "%s" in S or any of its substructures', field_name);
end

end

