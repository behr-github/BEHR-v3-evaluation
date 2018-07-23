function C = struct2cell2(S)
%STRUCT2CELL2 Convert structure to cell, preserving field names
%   C = STRUCT2CELL2(S) If S is a scalar structure, C will be a cell array
%   alternating field name and field value. If S is nonscalar, then the
%   values will themselves be cell arrays the same shape as S. So if S is a
%   2x2 structure with field 'alpha', and 
%
%       S(1,1).alpha = 1;
%       S(1,2).alpha = 2;
%       S(2,1).alpha = 3;
%       S(2,2).alpha = 4;
%
%   then C will be
%
%   {'alpha', {1, 2; 3, 4}}

fns = fieldnames(S);
sz = size(S);
C = cell(1, 2*numel(fns));

tmp_cell = struct2cell(S);
for i_fn = 1:numel(fns)
    i_cell = 2*(i_fn - 1) + 1;
    C{i_cell} = fns{i_fn};
    C{i_cell + 1} = format_field_cell(tmp_cell(i_fn,:), sz);
end

end

function c = format_field_cell(c_in, struct_sz)
if isscalar(c_in)
    c = c_in{1};
else
    c = reshape(c_in, struct_sz);
end
end