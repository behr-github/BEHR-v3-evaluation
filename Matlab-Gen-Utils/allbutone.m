function [ logi ] = allbutone( mat_in )
%allbutone - Returns TRUE if only one of the elements of the input matrix is false, FALSE otherwise.

% The logical output variable will be set to 1 if any of the matrices
% resulting from removing one element of the matrix in pass the "all" test.
logi = 0;

for a = 1:numel(mat_in)
    tmp_mat = mat_in;
    tmp_mat(a) = [];
    if all(tmp_mat);
        logi = 1;
        break
    end
end

end

