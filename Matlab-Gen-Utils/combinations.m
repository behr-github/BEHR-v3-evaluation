function [ M ] = combinations( varargin )
%COMBINATIONS Return all combinations of the values in the given vectors
%   M = COMBINATIONS( V1, V2, V3... ) Returns M which is a matrix containing
%   all combinations of the values in the input vectors V1, V2, V3, etc.
%   The values from each vector end up in the corresponding column of M.
%   Each row is one combination of values.

n = zeros(size(varargin));
for a = 1:numel(varargin)
    n(a) = numel(varargin{a});
end

M = nan(prod(n), numel(varargin));

inds = num2cell(zeros(size(n)));
for a = 1:prod(n)
    [inds{:}] = ind2sub(n, a);
    for b = 1:numel(varargin)
        M(a,b) = varargin{b}(inds{b});
    end
end

end

