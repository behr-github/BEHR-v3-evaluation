function [blockinds, blockvals] = findBlock(M, dim)
% FINDBLOCK Finds blocks of contiguous values in M
%   [ INDS, VALS ] = FINDBLOCK( M ) If M is a vector, FINDBLOCK operates
%   along the non-singleton dimension and defines blocks as a series of
%   adjacent elements with the same value, i.e. in [1 2 2 2 4 3 3 2 2],
%   blocks would be defined as the chunks separated by | below:
%
%    [ 1 | 2 2 2 | 4 | 3 3 | 2 2 ]
%
%   If M is an array, this syntax will operate along the first dimension
%   and defines blocks as adjacent slices which are identical:
%
%    [ 1 2
%      1 2
%    -------
%      2 2
%    -------
%      3 2
%      3 2
%      3 2 ]
%
%   If M has n elements along the operative dimension, then INDS will be an
%   n-by-2 matrix where the first column defines the start indicies of each
%   block and the second column the end indices. In the array example
%   above, INDS would be [1 2; 3 3; 4 6]. 
%
%   VALS will be n-by-PROD(SIZE(M,2:end)) and contain the values of each
%   block, be it a single value if M is a vector or a slice (flattened) if
%   M is an array.
%
%   [ ___ ] = FINDBLOCK( M, DIM ) operates along dimension DIM instead of
%   the first dimension. Only affects arrays (not vectors), a warning to
%   this effect is issued if used with a vector. VALS will now by n-by-x,
%   where x is the product of the lengths of the other dimensions.

E = JLLErrors;

% Check the input
narginchk(1,2);
if ~ismatrix(M)
    warning('findBlock has not been tested on arrays with >2 dimensions')
end
if ~exist('dim','var')
    dim = 1;
elseif ~isnumeric(dim) || ~isscalar(dim) || dim < 0 || dim > ndims(M) || mod(dim,1) ~= 0
    E.badinput('DIM must be a valid dimension of M')
elseif isvector(M)
    warning('DIM has no effect if M is a vector')
end



% If M is a vector, then we just need to make sure it is a column.
% Otherwise, permute the matrix as necessary to put the operative dimension
% first
if isvector(M)
    M = M(:);
else
    permvec = 1:ndims(M);
    permvec(dim) = [];
    permvec = [dim, permvec];
    M = permute(M, permvec);
end

% Set up blocks with the maximum possible number of blocks.  We'll remove
% extra entries at the end. It needs to have two columns for start and 
sz = size(M);
n = sz(1);
blockinds = nan(n, 2);
blockvals = nan(n, prod(sz(2:end)));

% Loop through M. For each entry, go through the following values
% until a different one is found.

i = 1;
b = 1;
while true 
    start_ind = i;
    chk_val = M(i,:);
    while true
        i = i+1;
        if i > n || any(M(i,:) ~= chk_val);
            i = i-1;
            break
        end
    end
    last_ind = i;
    blockinds(b,:) = [start_ind, last_ind];
    blockvals(b,:) = chk_val;
    b = b+1;
    i = i+1;
    if i > n
        break
    end
end

blockinds(b:end,:) = [];
blockvals(b:end,:) = [];

end