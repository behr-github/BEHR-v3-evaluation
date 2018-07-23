function [ C ] = mat2cell_simple( M )
% C = MAT2CELL_SIMPLE( M ) Converts array M to a cell array C on a 1-by-1 basis
%   NOTE: Deprecated, should be replaced by built-in NUM2CELL.
%
%   Unlike the built-in MAT2CELL this function assumes that you wish each
%   element of M to inhabit its own cell in C. Since the syntax of MAT2CELL
%   makes that annoying, this convenience function saves you the trouble.

warning('mat2cell_simple deprecated, use builtin NUM2CELL instead')

sz = cell(1,ndims(M));
for i=1:numel(sz)
    sz{i} = ones(1, size(M,i));
end

C = mat2cell(M, sz{:});


end

