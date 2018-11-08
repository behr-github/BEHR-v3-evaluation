function [mu, sigma] = weighted_mean(x, w, varargin)
%WEIGHTED_MEAN Compute a weight mean and standard deviation.
%   [MU, SIGMA] = WEIGHTED_MEAN( X, W ) Given values X and weights W,
%   computes the weighted mean (MU) and std. deviation (SIGMA) of X.
%   Operates along the first dimension, unless X and W are row vectors,
%   then operates along the vector dimension. X and W must be the same
%   size.
%
%   [ ___ ] = WEIGHTED_MEAN( ___, DIM ) Allows you to specify a different
%   dimension to operate along.

E = JLLErrors;
p = advInputParser;
p.addOptional('dim', 0);
p.parse(varargin{:});
pout = p.Results;

if ~isequal(size(x), size(w))
    E.badinput('X and W must be the same size');
elseif ~isnumeric(x) || ~isnumeric(w)
    E.badinput('X and W must be numeric');
end

dim = pout.dim;
if dim == 0
    if isrow(x)
        dim = 2;
    else
        dim = 1;
    end
end

if ~isnumeric(dim) || ~isscalar(dim) || dim < 1
    E.badinput('DIM must be a scalar number >= 1');
end

w_sum = nansum2(w, dim);
x_sum = nansum2(x .* w, dim);
mu = x_sum ./ w_sum;

% Following the formula for a standard deviation of a weighted mean from https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
weighted_ssr = nansum2(w .* (x - mu).^2, dim);
nz_wts = sum(~isnan(w) & w > 0, dim);
m_ratio = (nz_wts-1)./nz_wts;
sigma = sqrt(weighted_ssr ./ (m_ratio .* w_sum));

end

