function [ M ] = range_squeeze( M, varargin )
%RANGE_SQUEEZE Scale M so that its values fit in the range given
%   M = range_squeeze(M, A, B) Scales array M so that its values are mapped
%   between A and B, i.e. min(M(:)) --> A and max(M(:)) --> B. M may be an
%   numeric array of any size; A and B must be scalar values with A < B.
%
%   M = RANGE_SQUEEZE(M, [A, B]) Alternative syntax to above.

E = JLLErrors;

if nargin == 2
    rng = varargin{1};
elseif nargin == 3
    rng = cell2mat(varargin(1:2));
else
    E.badinput('RANGE_SQUEEZE takes two or three input arguments')
end

if ~isnumeric(rng) || numel(rng) ~= 2 || rng(1) > rng(2)
    E.badinput('The range must either be given as two scalar numeric inputs with the first less than the second, or a two element numeric vector with the first element less than the second')
end

if ~isnumeric(M)
    E.badinput('M must be numeric')
end

M = M - min(M(:)) + rng(1);
M = M .* rng(2) ./ max(M(:));

end

