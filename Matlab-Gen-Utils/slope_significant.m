function is_sig = slope_significant(P, x, y, varargin)
% SLOPE_SIGNIFICANT Test is a slope is significant
%
%   IS_SIG = SLOPE_SIGNIFICANT( P, X, Y ) Given the fit polynomial
%   P (a two element vector: [slope intercept]) and the X and Y
%   values fit by it, returns true if the slope is different from 0
%   at the 95% confidence level.
%
%   IS_SIG = SLOPE_SIGNIFICANT( P, X, Y, M ) Tests if the slope
%   is significantly different from M, instead of 0.
%
%   Parameters:
%
%       'confidence' - the confidence level to test at. Default is
%       0.95.

% TODO: add citation from the stats book
E = JLLErrors;
p = inputParser;
p.addOptional('test_slope', 0, @(x) isscalar(x) && isnumeric(x));
p.addParameter('confidence', 0.95);

p.parse(varargin{:});
pout = p.Results;

test_slope = pout.test_slope;
confidence = pout.confidence;

if ~isnumeric(P) || numel(P) ~= 2
    E.badinput('P must be a two element numeric vector')
elseif ~isnumeric(x) || ~isnumeric(y)
    E.badinput('x and y must both be numeric')
elseif ~isequal(size(x), size(y))
    E.badinput('x and y must be the same size');
end

if ~isnumeric(confidence) || ~isscalar(confidence) || confidence < 0 || confidence > 1
    E.badinput('"confidence" must be a scalar number between 0 and 1')
end

% tinv gives 1-sided t values. For a two sided t-test (i.e. one where
% the values may be above or below each other) the same percent chance
% of being different must be split onto either side. Concretely, for a
% 95% two-sided confidence interval, we need to compare against the
% 97.5% t-value returned by tinv.
confidence = 1 - ( (1 - confidence)/2 );
n_vals = sum(~isnan(x) & ~isnan(y));

sigma = sqrt(nansum2((y - P(1) .* x + P(2)).^2)./(n_vals - 2));
stderr_slope = sigma ./ sqrt(nansum2((x - nanmean(x)).^2));
t_slope = (P(1) - test_slope) ./ stderr_slope;
t_table = tinv(confidence, n_vals - 2);

is_sig = t_slope >= t_table;

end

