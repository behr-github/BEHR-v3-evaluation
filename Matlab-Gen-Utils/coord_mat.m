%Returns a matrix whose values represent the xyz coordinates (or xy
%coordinates for a 2D matrix). matrices of greater than length 9 will be
%weird.

function matrix_out = coord_mat(varargin)%x, y, z)

dims = [varargin{:}];
if numel(dims) == 1
    dims(2) = 1;
end
matrix_out = zeros(prod(dims),1);
place = 10.^(0:numel(dims)-1);
for a=1:prod(dims)
    C = cell(1,numel(dims));
    [C{:}] = ind2sub(dims,a);
    C = [C{:}];
    matrix_out(a) = sum(C .* place);
end
matrix_out = reshape(matrix_out,dims);