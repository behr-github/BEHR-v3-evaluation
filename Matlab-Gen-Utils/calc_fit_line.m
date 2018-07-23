function [ line_x, line_y, legend_str, LineData ] = calc_fit_line( x, y, varargin )
%CALC_FIT_LINE Calculates a fit line for data in X and Y
%   [ LINE_X, LINE_Y, LEGEND_STR, LINEDATA ] = CALC_FIT_LINE( X, Y )
%   Calculates a line of best fit for the given data after removing NaNs.
%   Returns vectors of x & y points to plot the line (LINE_X and LINE_Y), a
%   string formatted with information on the fit line for a legend
%   (LEGEND_STR), and a structure containing the slope/intercept vector P
%   (compatible with polyval), the R^2 values, std. deviations of slope and
%   intercept, and p-value (LINEDATA).
%
%   This makes use of the various linear regression functions by Edward
%   Peitzer (http://www.mbari.org/staff/etp3/regress/index.htm). Which
%   regression is used can be set by the 'regression' parameter.
%       
%       'y-resid' - (default) minimizes differences between the points and
%       the fit line in the y-direction only.  Best when x is considered to
%       have no error.
%
%       'x-resid' - minimized differences between the points and the fit
%       line in the x-direction only.
%
%       'majoraxis' - minimized the normal distance between the fit line
%       and the points. Said to fit a line along the major axis of an
%       ellipse that represents the variance.  Should only be used if x and
%       y have the same units and scale, and if both x & y have
%       non-negligible error.
%
%       'orth-origin' - minimized the mean square orthogonal difference
%       between the points and the slope, forcing the slope through the
%       origin.
%
%       'RMA' - takes a geometric mean of the y-on-x and x-on-y regression.
%       Effectively minimizes the area of a triangle between each point and
%       the fit line. Most computationally expensive, but good when x & y
%       may not have the same scale or units and both have error terms.
%
%   http://www.unc.edu/courses/2007spring/biol/145/001/docs/lectures/Nov5.html
%   has a good explanation of the various sorts of regressions.
%
%   You can alter what x-coordinates for the line are returned with the
%   parameter 'xcoord' which accepts the following values:
%
%       'figxlim' (default) - line will have the x coordinates as the
%       limits of the current figure. If no figure open, defaults to 0 to
%       1.
%
%       'dataxlim' - line will have x coordinates as the min and max of the
%       data.
%
%       [xcoords] - giving a monotonically increasing numeric vector will
%       generate the line using those values as the x-coordinates.

E = JLLErrors;

default_x_coord = 0:1;

p = inputParser;
p.addRequired('x',@isnumeric);
p.addRequired('y',@isnumeric);
p.addParameter('pvalue',false);
p.addParameter('significance',false);
p.addParameter('regression','y-resid');
p.addParameter('xcoord','figxlim');

p.parse(x,y,varargin{:});
pout = p.Results;
x = pout.x;
y = pout.y;
regression = lower(pout.regression); % explicitly make the regression string lower case to ease comparison in the switch-case statement
pvalue = pout.pvalue;
xcoord = lower(pout.xcoord);
mark_significance = pout.significance;

allowed_regressions = {'y-resid','x-resid','majoraxis','RMA','orth-origin'};
if ~any(strcmpi(regression, allowed_regressions))
    E.badinput('The parameter ''REGRESSION'' must be one of %s', strjoin(allowed_regressions, ', '));
end
allowed_xcoords = {'figxlim','dataxlim'};
if isnumeric(xcoord)
    if numel(xcoord) < 2 || ~isvector(xcoord) || any(diff(xcoord) < 0)
        E.badinput('If giving the parameter ''XCOORD'' as numeric values, it must be a monotonically increasing vector with at least two elements');
    end
elseif ischar(xcoord) && ~any(strcmpi(xcoord, allowed_xcoords)) || ~ischar(xcoord)
    E.badinput('The parameter ''XCOORD'' must be a numerical vector or one of the strings %s', strjoin(allowed_xcoords, ', ')); 
end

if any(isnan(x)) || any(isnan(y))
    warning('fitting:nans_removed','NaNs detected, removing any points with a value of NaN for either coordinate');
    nans = isnan(x) | isnan(y);
    x = x(~nans); y = y(~nans);
end

% Make x and y into vectors
x = x(:);
y = y(:);

switch regression
    case 'y-resid'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfity(x,y);
        R = R^2; % All of the lsqfit** functions do not square R before returning
        
    case 'x-resid'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitx(x,y);
        
    case 'majoraxis'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitma(x,y);
        R = R^2; 
        
    case 'rma'
        [P(1), P(2), R, sigma_m, sigma_b] = lsqfitgm(x,y);
        R = R^2;
        
    case 'orth-origin'
        [P(1), sigma_m] = lsqfitnmorg(x,y);
        P(2) = 0;
        R = nan;
        sigma_b = nan;
end

p_val = p_val_slope(P(1), sigma_m, sum(~isnan(x) & ~isnan(y)));
is_significant = slope_significant(P, x, y);

LineData.P = P;
LineData.R2 = R;
LineData.StdDevM = sigma_m;
LineData.StdDevB = sigma_b;
LineData.p_value = p_val;
LineData.is_significant = is_significant;
LineData.num_pts = sum(~isnan(x) & ~isnan(y));
   
if isnumeric(xcoord)
    line_x = xcoord;
elseif strcmp(xcoord, 'figxlim') && ~isempty(get(0,'children'))
    line_x = get(gca,'xlim');
elseif strcmp(xcoord, 'dataxlim')
    line_x = [min(x(:)), max(x(:))];
else
    line_x = default_x_coord;
end
line_y = P(1) * line_x;
if ~isnan(P(2))
    line_y = line_y + P(2);
end


if mark_significance && is_significant
    sig_asterisk = '*';
else
    sig_asterisk = '';
end

legend_str = sprintf('Fit: %.4fx',P(1));
if ~isnan(P(2))
    legend_str = [ legend_str, sprintf(' + %.2g', P(2)) ];
end
if ~isnan(R)
    if pvalue
        legend_str = [ legend_str, sprintf(' \nR^2 = %.4f (p = %.2f)',R,p_val) ];
    else
        legend_str = [ legend_str, sprintf(' \nR^2 = %.4f%s',R,sig_asterisk) ];
    end
end

end

