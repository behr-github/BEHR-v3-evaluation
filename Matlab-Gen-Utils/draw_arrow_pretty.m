function draw_arrow_pretty(x, y, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = advInputParser;
p.addOptional('headfrac', 0.025);
p.addParameter('headlength', []);
p.addParameter('headwidth', []);
p.addParameter('linewidth', 1);
p.addParameter('linestyle', '-');
p.addParameter('color', 'k');
p.addParameter('parent', gca);

p.parse(varargin{:});
pout = p.Results;


arrow_head_fraction = pout.headfrac;
arrow_head_length = pout.headlength;
arrow_head_width = pout.headwidth;
line_width = pout.linewidth;
line_style = pout.linestyle;
arrow_color = pout.color;
parent = pout.parent;

[line_x, line_y, head_x, head_y] = calculate_coordinates(x, y, arrow_head_fraction, arrow_head_length, arrow_head_width);
make_arrow(line_x, line_y, head_x, head_y);

    function h = make_arrow(line_x, line_y, head_x, head_y)
        h = hggroup(parent);
        lhandle = line(h, line_x, line_y, 'color', arrow_color, 'linewidth', line_width, 'linestyle', line_style);
        phandle = patch(h, head_x, head_y, arrow_color, 'edgecolor', arrow_color);
    end

end

function [line_x, line_y, head_x, head_y] = calculate_coordinates(x, y, head_frac, head_len, head_width)

% First, if the head length or width are not specified, calculate them
% based on the arrow head fraction. Assume an equilateral triangle.
xdel = x(2) - x(1);
ydel = y(2) - y(1);
vector_length = sqrt(xdel.^2 + ydel.^2);
vector_angle = atan2(ydel, xdel);

head_size = vector_length * head_frac;
if isempty(head_len)
    head_len = head_size;
end
if isempty(head_width)
    head_width = head_size/2;
end

% Now calculate where the line will end. Basically we need to move the end
% back by the length of the arrow head
line_x(1) = x(1);
%line_x(2) = x(1) + (xdel - head_len * cos(vector_angle));
line_y(1) = y(1);
%line_y(2) = y(1) + (ydel - head_len * sin(vector_angle));
[line_x(2), line_y(2)] = coord_along_vector(x(2), y(2), -head_len, vector_angle);

% Next we will calculate the arrow points. Fortunately, patch() doesn't
% seem terribly picky about the points being clockwise or counterclockwise.
% The tip of the arrow will be at [x(2) y(2)], but then the other points
% will be head_width out to the side of the end of the line part.
head_x(1) = x(2);
head_y(1) = y(2);
[head_x(2), head_y(2)] = coord_along_vector(line_x(2), line_y(2), head_width, vector_angle - pi/2);
[head_x(3), head_y(3)] = coord_along_vector(line_x(2), line_y(2), head_width, vector_angle + pi/2);

end



function [x,y] = coord_along_vector(xstart, ystart, r, theta)
x = xstart + r .* cos(theta);
y = ystart + r .* sin(theta);
end