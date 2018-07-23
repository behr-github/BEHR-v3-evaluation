function [ fighandle ] = plot_error_envelope_x( y_in, x_lower, x_upper, varargin )
%plot_error_envelope_x(y_in, x_lower, x_upper, [opt] colorspec, [opt] figure number): Uses the fill function to plot an error envelope for error in x value.
%   This plots an error envelope on your plot for error in the x_values.
%   This is set up to allow the use of standard deviation/std. error or
%   quartiles, hence the lower and upper bounds are input separately.  The
%   default color of the envelope is light grey, but the color can be
%   specified as an optional argument.  Most likely this will need to be
%   called before you plot your actual data, as the fill will probably
%   cover the line if plotted afterwards.
%
%   Returns the figure number it used to plot.

p = inputParser;
p.addRequired('y',@isnumeric);
p.addRequired('x_lower',@isnumeric);
p.addRequired('x_upper',@isnumeric);
p.addParameter('colorspec',[0.8 0.8 0.8]);
p.addParameter('facealpha',1,@isscalar);
p.addParameter('ax',[]);

p.parse(y_in, x_lower, x_upper, varargin{:});
pout = p.Results;

y = pout.y;
xl = pout.x_lower;
xu = pout.x_upper;

colspec = pout.colorspec;
ax = pout.ax;
facealpha = pout.facealpha;

% If no figure handle is passed, create one. If one is, switch to that
% figure and turn hold on so that we don't delete existing data.
if isempty(ax)
    ax = gca;
end
hold on

X = [xl(:); flipud(xu(:))];
Y = [y(:); flipud(y(:))];

fill(ax, X,Y,colspec,'edgecolor','none','FaceAlpha',facealpha);
hold off

end

