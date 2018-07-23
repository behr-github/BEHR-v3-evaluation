function [ A ] = clipmat( A, varargin )
%CLIPMAT Restrict an array to a range of values
%   A = CLIPMAT( A, X, Y ) clamps any value of A < min(X,Y) to min(X,Y) and
%   any value of A > max(X,Y) to max(X,Y). Therefore the order of X, Y does
%   not matter.
%
%   A = CLIPMAT( A, [X, Y] ) alternate syntax, but same behavior.

E = JLLErrors;

if nargin == 3
    limits = [varargin{1}, varargin{2}];
elseif nargin == 2
    limits = varargin{1};
else
    E.badinput('CLIPMAT requires 2 or 3 inputs')
end

A = min(A, max(limits));
A = max(A, min(limits));

end

