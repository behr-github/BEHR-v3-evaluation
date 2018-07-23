function [ fighandle ] = plot_error_envelope_y( x_in, y_lower, y_upper, varargin )
%plot_error_envelope_y(x_in, y_lower, y_upper, [opt] colorspec, [opt] figure handle): Uses the fill function to plot an error envelope for error in x value.
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
p.addOptional('colorspec',[0.8 0.8 0.8], @(t) (isnumeric(t) || ischar(t)));
p.addOptional('fignum',-1, @(t) (ishandle(t) && strcmp(get(t,'type'),'figure')));

p.parse(x_in, y_lower, y_upper, varargin{:});
pout = p.Results;

x = pout.y;
yl = pout.x_lower;
yu = pout.x_upper;

%Convert all to column vectors
x = x(:);
yl = yl(:);
yu = yu(:);

colspec = pout.colorspec;
nfig = pout.fignum;

% If no figure handle is passed, create one. If one is, switch to that
% figure and turn hold on so that we don't delete existing data.
if nfig == -1
    fighandle = figure();
else
    fighandle = figure(nfig);
    hold on
end


Y = [yl; flipud(yu)];
X = [x; flipud(x)];

fill(X,Y,colspec,'edgecolor','none');
hold off

end

