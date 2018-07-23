function [ log_mat ] = test_str( str_in, test_in )
%test_str(str_in, test_in): Tests if str_in exists in test_in.
%   Uses regexp to test if str_in exists in test_in. If test_in is a string
%   as well, this function will return a scalar value. If test_in is a cell
%   structure, this function will return a logical matrix.

log_1 = regexp(test_in, str_in);
log_mat = true(size(log_1));

for a=1:numel(log_1)
    if isempty(log_1{a}); log_mat(a) = 0;
end

end

