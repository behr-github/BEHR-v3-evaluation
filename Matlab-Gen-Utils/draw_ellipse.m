function [ varargout ] = draw_ellipse( center_x, center_y, x_radius, y_radius, varargin )
%DRAW_ELLIPSE Draws an ellipse on a plot.
%   DRAW_ELLIPSE( CENTER_X, CENTER_Y, X_RADIUS, Y_RADIUS ) draws an ellipse
%   centered on CENTER_X and CENTER_Y with given radii in the x and y
%   direction. This implicitly assumes that the x and y axes of the plot
%   are in the same units.
%
%   DRAW_ELLIPSE( CENTER_X, CENTER_Y, X_RADIUS, Y_RADIUS, lineopts ) You can
%   also pass any parameter key-value pairs that you can pass on to the
%   built-in function LINE.
%
%   [ GOBJECT ] = DRAW_ELLIPSE( ___ ) returns the graphic object
%   representing the circle.
%
% Answer at https://www.mathworks.com/matlabcentral/answers/86615-how-to-plot-an-ellipse
% by Azzi Abdelmalek

t=-pi:0.01:pi;
x=center_x + x_radius .* cos(t);
y=center_y + y_radius .* sin(t);
l = line(x,y,varargin{:});

if nargout > 0
    varargout{1} = l;
end

end

