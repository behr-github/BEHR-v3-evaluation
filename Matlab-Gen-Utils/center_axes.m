function [ ] = center_axes(ax_to_center, left_ax, right_ax)
%CENTER_AXES Center an axes object between two other objects
%   CENTER_AXES( AX_TO_CENTER, LEFT_AX, RIGHT_AX ) moves the horizontal
%   position of AX_TO_CENTER so that it's midpoint is approximately halfway
%   between the left edge of LEFT_AX and the right edge of RIGHT_AX. All
%   three must be handles to graphics objects with the Position property.
%
%   AX_TO_CENTER will not be resized, but the calculation does depend on
%   its horizontal size, so any resizing should be done before calling
%   CENTER_AXES().

left_edge = left_ax.Position(1);
right_edge = right_ax.Position(1) + right_ax.Position(3);
centered_x = (left_edge + right_edge)/2 - ax_to_center.Position(3)/2;

ax_to_center.Position(1) = centered_x;

end

