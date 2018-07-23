function [ varargout ] = draw_circle( center_x, center_y, radius, varargin )
%DRAW_CIRCLE Draws a circle on a plot.
%   DRAW_CIRCLE( CENTER_X, CENTER_Y, RADIUS ) draws a circle centered on
%   CENTER_X and CENTER_Y with given RADIUS. This implicitly assumes that
%   the x and y axes of the plot are in the same units.
%
%   DRAW_CIRCLE( CENTER_X, CENTER_Y, RADIUS, lineopts ) You can also pass
%   any parameter key-value pairs that you can pass on to the built-in
%   function LINE.
%
%   [ GOBJECT ] = DRAW_CIRCLE( ___ ) returns the graphic object
%   representing the circle.


theta = 0:0.01:6.28;
x = cos(theta)*radius + center_x;
y = sin(theta)*radius + center_y;

l=line(x,y,varargin{:});

if nargout > 0
    varargout{1} = l;
end

end

