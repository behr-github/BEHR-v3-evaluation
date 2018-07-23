function [ LineData ] = plot_fit_line( x, varargin )
%plot_fit_line Draws a fit line and 1:1 line on the current figure.
%   Automatically draws a fit line on the current plot (given the data on
%   the plot as input variables - i.e. after plot(x,y) call
%   plot_fit_line(x,y) to make the lines). Also creates a legend with the
%   linear fit info.  By default also draws a 1:1 line.  This can be
%   disabled by setting the parameter 'one2one' to 0.
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
%   This function has two other parameters:
%
%       one2one - a boolean indicating whether to plot a 1:1 line. Defaults
%       to true.
%
%       pvalue - a boolean indicating whether to plot the p value in the
%       legend box.
%
%       sigma_m - a number indicating how many standard deviations of the
%       slope to plot as an envelope.  Defaults to 0, i.e. no envelope will
%       be plotted.
%
%       sigma_b - a number indicating how many standard deviations of the
%       intercept to plot as an envelope.  Defaults to 0, i.e. no envelope
%       will be plotted.

E = JLLErrors;

p = inputParser;
p.addRequired('x');
p.addOptional('y',[]);
p.addParameter('one2one',true);
p.addParameter('pvalue',false);
p.addParameter('regression','y-resid');
p.addParameter('sigma_m',0,@isscalar);
p.addParameter('sigma_b',0,@isscalar);

p.parse(x,varargin{:});
pout = p.Results;
x = pout.x;
y = pout.y;
one2one = pout.one2one;
regression = pout.regression;
pvalue = pout.pvalue;
mult_sigma_m = pout.sigma_m;
mult_sigma_b = pout.sigma_b;

% Input validation
if isa(x,'matlab.graphics.axis.Axes') 
    [x,y] = get_child_data(x);
else
    if ~isnumeric(x) || ~isnumeric(y) && ~isempty(y)
        E.badinput('X and Y must be numeric if not giving a set of axes')
    end
end
        

% Make x & y into vectors
x = x(:);
y = y(:);

if any(isnan(x)) || any(isnan(y))
    warning('NaNs detected, removing any points with a value of NaN for either coordinate');
    nans = isnan(x) | isnan(y);
    x = x(~nans); y = y(~nans);
end

[xline, yline, lstr, LineData] = calc_fit_line(x, y, 'regression', regression,'pvalue',pvalue);

% Plot the line of best fit and format the string for the legend with the
% slope/intercept info and the R^2 value.
h(1) = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legendcell = {lstr};

% Plot the 1:1 line
if one2one
    h(end+1) = line(xline,xline,'color','r','linestyle',':','linewidth',2);
    legendcell{end+1} = '1:1';
end


% Plot an envelope showing the variance in the slope, intercept, or both.
% This will only plot a single envelope, showing the maximum variance,
% adding the slope and intercept in the same direction.  The std.
% deviations will be multiplied by the parameters sigma_m and sigma_b, so
% as to allow the user to plot a 2-sigma envelope.
if mult_sigma_m > 0 || mult_sigma_b > 0
    slope_plus_sigma = LineData.P(1) + LineData.sigma_m * mult_sigma_m;
    int_plus_sigma = LineData.P(2) + LineData.sigma_b * mult_sigma_b;
    
    slope_minus_sigma = LineData.P(1) - LineData.sigma_m * mult_sigma_m;
    int_minus_sigma = LineData.P(2) - LineData.sigma_b * mult_sigma_b;
    
    P_upper = [slope_plus_sigma, int_plus_sigma];
    P_lower = [slope_minus_sigma, int_minus_sigma];
    
    env_y_upper = polyval(P_upper,xline);
    env_y_lower = polyval(P_lower,xline);
    
    env_x = [xline, fliplr(xline)];
    env_y = [env_y_upper, fliplr(env_y_lower)];
    
    h(end+1) = patch(env_x,env_y,'b','facealpha',0.5,'linestyle','-','linewidth',1);
    
    % The format string to fill in the standard deviations and multipliers
    env_str = 'Var: %2$d%5$s_m = %2$d * %1$.2f\n     %4$d%5$s_b = %4$d * %3$.2G';
    legendcell{end+1} = sprintf(env_str,sigma_m,mult_sigma_m,sigma_b,mult_sigma_b,char(hex2dec('03C3')));
end

% legend req. that the vector of handles be a column vector, hence we
% transpose it.
legend(h',legendcell{:})
end

function [x,y] = get_child_data(ax)
E = JLLErrors;
if ~isscalar(ax)
    E.badinput('PLOT_FIT_LINE can only take a single set of axes, not an array')
end

ch = cat(1,findobj(ax,'type','line'),findobj(ax,'type','scatter'));
if isempty(ch)
    E.badinput('Input axes have no data')
elseif ~isscalar(ch)
    E.notimplemented('%s','Multiple child data detected, fitting multiple series not implemented')
end

x = get(ch,'XData');
y = get(ch,'YData');
end
