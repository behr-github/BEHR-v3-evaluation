function [ val_lims ] = calc_plot_limits( values, varargin )
%CALC_PLOT_LIMITS Compute ideal limits for plotting
%   VAL_LIMS = CALC_PLOT_LIMITS( VALUES ) Returns the 1x2 vector VAL_LIMS
%   that contains the minimum and maximum values in VALUES, rounded to the
%   next lower and next higher integer, respectively.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( VALUES, ROUND_TO ) Rounds the limits to
%   the next multiple of ROUND_TO instead of 1. ROUND_TO must be a scalar
%   number >= 0. If ROUND_TO == 0, the limits will not be rounded at all.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( ___, 'zero' ) will cause one side of the
%   limits to be set to 0, depending on whether the values are positive or
%   negative: if positive, the lower limit will be 0, if negative, the
%   upper limit. If the values are both positive and negative, a warning
%   will be given and the behavior will be as if 'zero' was not specified.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( ___, 'difference' ) will force the upper
%   and lower limits to have the same magnitude, but opposite sign. This is
%   useful when ploting differences as color where having 0 be the exact
%   middle of the color scale makes it easier to see (e.g. if using a
%   colormap where red should be positive and blue negative). This is
%   mutually exclusive with 'zero'.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( ___, 'diff' ) identical to the previous
%   syntax, just a shorter flag name.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( ___, 'pow10' ) scales ROUND_TO by:
%       10^(floor(log10(max(abs(values(:)))))-1)
%   i.e. ROUND_TO is taken to be one order of magnitude less than the
%   maximum value. This can help make a plot's limits more generic if you
%   are plotting different data in the same function and the data varies by
%   many orders of magitude; e.g. NO2 column densities are ~10^15, while
%   surface reflectance is ~1, if we want to plot these in a generic
%   function, rounding the limits to some multiple of 5, using 0.05 makes
%   sense for reflectance but 5e14 makes sense for NO2.
%
%   VAL_LIMS = CALC_PLOT_LIMITS( ___, 'max', [MINVAL MAXVAL] ) will clip
%   the limits to the range [MINVAL MAXVAL] if they exceed those values.

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

E = JLLErrors;

if ~isnumeric(values)
    E.badinput('VALUES must be numeric');
end

p = advInputParser;
p.addOptional('round_to', 1, @isnumeric);
p.addFlag('diff');
p.addFlag('difference');
p.addFlag('zero');
p.addFlag('pow10');
p.addParameter('max', [-Inf Inf]);
p.addParameter('sigma', 0);

p.parse(varargin{:});
pout = p.AdvResults;

round_to = pout.round_to;
make_equal_pos_neg = pout.diff || pout.difference;
force_zero = pout.zero;
scale_by_log = pout.pow10;
max_lims = pout.max;
sigma_lims = pout.sigma;

if sigma_lims > 0
    values(abs(values) > sigma_lims * nanstd(values(:))) = [];
end


if make_equal_pos_neg && force_zero
    E.badinput('''diff'' or ''difference'' and ''zero'' are mutually exclusive options');
end
    
if make_equal_pos_neg
    limit = max(abs(values(:)));
    val_lims = [-limit, limit];
elseif force_zero
    val_lims = [min(values(:)), max(values(:))];
    if all(val_lims >= 0)
        val_lims(1) = 0;
    elseif all(val_lims <= 0)
        val_lims(2) = 0;
    else
        warning('Values are both positive and negative, will not set one limit to 0');
    end
else
    val_lims = [min(values(:)), max(values(:))];
end

if scale_by_log
    max_pow_10 = floor(log10(max(abs(values(:)))));
    round_to = round_to * 10^(max_pow_10-1);
end

if round_to > 0
    val_lims(1) = floor(min(val_lims)/round_to)*round_to;
    val_lims(2) = ceil(max(val_lims)/round_to)*round_to;
end

% Apply any user specified max limits. If no limits specified, these are
% -Inf and Inf, so any finite limit will be permitted.

val_lims(1) = max(val_lims(1), max_lims(1));
val_lims(2) = min(val_lims(2), max_lims(2));

end

