%This function will return a matrix of the specified dim'ns with each entry
%labeled with its sequential index.

function matrix_out = seq_mat(varargin)

matrix_out = zeros(varargin{:});
for a=1:numel(matrix_out)
    matrix_out(a) = a;
end

end