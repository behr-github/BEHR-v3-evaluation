function [ output_args ] = n_max( X, n )
%n_max Returns the n largest values and the positions

E = JLLErrors;
warning('n_max will be deprecated, replace using findbysize')
narginchk(1,2);

% This function is currently only working for vectors, so it will operate
% along the vector's length.
if nargin > 2
    warning('Specified dimension will not be used, that function has not been incorporated yet')
end

if ~isvector(X)
    error(E.callError('X_not_vector','This function is currently not able to accept non-vector inputs'));
end

[M,I] = sort(X);
Xs = M(end-n+1:end);
Is = I(end-n+1:end);



end

