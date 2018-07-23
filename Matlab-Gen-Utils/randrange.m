function [ r ] = randrange( a, b, varargin )
%RANDRANGE Produces a range of random numbers between two values
%   R = RANDRANGE( A, B, N ) produces an NxN matrix of random floating
%   point values between A and B.
%
%   R = RANDRANGE( A, B, N1, N2, ... , Nn ) Produces an N1 x N2 x ... Nn
%   n-dimensional array of random values.

E = JLLErrors;

if ~isnumeric(a) || ~isscalar(a)
    E.badinput('A must be a scalar number');
elseif ~isnumeric(a) || ~isscalar(a)
    E.badinput('A must be a scalar number');
elseif numel(varargin) < 1
    E.badinput('At least three arguments must be supplied')
elseif any(~iscellcontents(varargin, 'isnumeric')) || any(~iscellcontents(varargin, 'isscalar')) || any(iscellcontents(varargin, @(x) x < 1)) || any(iscellcontents(varargin, @(x) mod(x,1) ~= 0))
    E.badinput('All the inputs from the third position on must be scalar, positive, whole numbers');
end

r = (b-a)*rand(varargin{:})+a;

end

