function new_struct = make_struct_from_field_values(fieldnames, fieldvalues)
%MAKE_STRUCT_FROM_FIELD_VALUES Generate a structure from a list of fieldnames and the values they should take on.
%   S = MAKE_STRUCT_FROM_FIELD_VALUES( FIELDNAMES, FIELDVALUES ) Creates a
%   struct S with fields FIELDNAMES and values FIELDVALUES. FIELDNAMES must
%   be a cell array of chars or a string array; FIELDVALUES must be an
%   array of some sort with an equal number of elements to that of
%   FIELDNAMES. 

E = JLLErrors;

if isstring(fieldnames)
    fieldnames = cellstr(fieldnames);
elseif ~iscellstr(fieldnames)
    E.badinput('FIELDNAMES must be a string array or cell array of chars');
end

if numel(fieldvalues) ~= numel(fieldnames)
    E.badinput('Number of elements in FIELDVALUES and FIELDNAMES must be equal');
end

arg_cell = cell(1,2*numel(fieldnames));
for i_field = 1:numel(fieldnames)
    i_arg = (i_field-1)*2 + 1;
    arg_cell{i_arg} = fieldnames{i_field};
    if iscell(fieldvalues)
        arg_cell{i_arg + 1} = fieldvalues{i_field};
    else
        arg_cell{i_arg + 1} = fieldvalues(i_field);
    end
end

new_struct = struct(arg_cell{:});

end

