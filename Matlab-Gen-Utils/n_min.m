function [ Xs, Is ] = n_min( X, n, dim )
%n_min Finds the smallest n values in X and their indices

E = JLLErrors;
warning('n_min will be deprecated, replace using findbysize')

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
Xs = M(1:n);
Is = I(1:n);

% This is the start of my idea to handle matrices, but I'm not sure how to
% implement it totally (JLL 18 Sept 2014)
% [M,I] = sort(X,dim);
% M = permute(M,[dim,1,3:ndims(M)]);
% I = permute(I,[dim,1,3:ndims(I)]);
% s = size(M);
% Xs = reshape(M(1:n,:),s(2:end));


end

