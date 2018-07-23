function [ A ] = bitopmat( A, op, dim )
%bitopmat Applies a bitwise AND or OR to a matrix or vector
%   As sum() or mean() can operate along various dimensions of a matrix,
%   this function can apply a bitwise AND, OR, or XOR to a vector or along
%   one dimension of a matrix. Takes 2 required, 1 optional argument:
%
%       A: the matrix or vector to operate on. It must be of integer type.
%
%       op: a string ('and', 'or', 'xor') that specifies what operation to
%       carry out.
%
%       dim (optional): what dimension to operate along.  Vectors will
%       always be operated along their longest dimension, matrices will by
%       default be operated on along the first dimension, i.e. a 3 x 3
%       matrix will become a 1 x 3.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,3);

if ~isinteger(A)
    E.badinput('A must be a matrix or vector of integers (int or uint type)');
end

if ~ischar(op)
    E.badinput('op must be a string');
end
% Make op lowercase to simplify comparisons
op = lower(op);
if ~ismember(op,{'and','or','xor'})
    E.badinput('op must be one of the strings: ''and'', ''or'', ''xor''');
end
op = str2func(sprintf('bit%s',op));

% Figure out what dimension to operate along. If the user gave a dimension,
% make sure its valid; if not, operate along the first dimension of a
% matrix or the long dimension of a vector.

if nargin > 2
    if ~isnumeric(dim)
        E.badinput('dim must be a number')
    elseif dim < 1 || dim > ndims(A)
        E.badinput('dim must be a valid dimension of A')
    end
else
    if isvector(A)
        [~,dim] = max(size(A));
    else
        dim = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% The way we'll do this is to permute A so that the dimension to operate
% along is the first dimension, do the operation, then unpermute it.

perm_vec = 1:ndims(A);
perm_vec(1) = dim;
perm_vec(dim) = 1;

A = permute(A,perm_vec);
sz = size(A);
n = sz(1);
sz(1) = 1;

% B will be the matrix we store our results in
B = nan(sz);

% We need to iterate over the second, third, fourth, etc. dimensions of A.
% In Matlab 2014b, specifying only 2 coordinates for a 3 or higher
% dimensional matrix "linearizes" the 2nd and higher dimensions. For
% example, if A is 2x2x2, then A(1,2) = A(1,2,1) and A(1,3) = A(1,1,2). If
% this behavior changes, that will need to be addressed.

for d=1:prod(sz)
    bit = A(1,d);
    for a=2:n
        bit = op(bit,A(a,d));
    end
    B(d) = bit;
end

% Unpermute B to get A back in the right shape
A = permute(B,perm_vec);


end

