function [ varargout ] = label_axis_with_letter( label, varargin )
%LABEL_AXIS_WITH_LETTER Label a given axis with a string
%   LABEL_AXIS_WITH_LETTER( LABEL ) Places the string LABEL at the top left
%   corner of the current axes. If LABEL is a scalar number, it is
%   converted to a lowercase letter matching its place in the alphabet and
%   that letter, enclosed in parentheses, is used as the label.
%
%   Parameter arguments:
%       axis - give a handle to an axis to label instead of the current
%       one.
%
%       xshift - additional displacement from the axes in the x direction
%       as a fraction of the width of the plot. Default is 0.12.
%
%       yshift - additional displacement in the y direction as xshift is in
%       the x direction. Default is 0.
%
%       fontweight - any valid method of specifying a font weight to the
%       TEXT function. Default is 'b' (i.e. bold).
%
%       fontcolor - any valid color specification; changes the color of the
%       text. Default is 'k' (black).
%
%       fontsize - number specifying the font size in points. Default is
%       16.

E = JLLErrors;
if isnumeric(label) && isscalar(label)
    if label > 26
        warning('Numeric value of LABEL exceeds number of letters in alphabet');
    end
    label = sprintf('(%s)',char(label+96));
end
if ~ischar(label)
    E.badinput('LABEL must be a string or a scalar number')
end

p = inputParser;
p.addParameter('axis', gobjects(0));
p.addParameter('xshift', 0.12);
p.addParameter('yshift', 0);
p.addParameter('fontweight', 'b');
p.addParameter('fontcolor', 'k');
p.addParameter('fontsize', 16);

p.parse(varargin{:});
pout = p.Results;

ax = pout.axis;
shift_percent_x = pout.xshift;
shift_percent_y = pout.yshift;
fontweight = pout.fontweight;
fontcolor = pout.fontcolor;
fontsize = pout.fontsize;

if ~isempty(ax)
    if ~ishandle(ax) || ~strcmp(get(ax,'type'), 'axes')
        E.badinput('AX (if given) must be a handle to axes')
    end
else
    ax = gca;
end

is3D = ~isempty(get(ax,'ZTickLabel'));

if ~is3D
    text_handle = label2D();
else
    text_handle = label3D();
end

if nargout > 0
    varargout{1} = text_handle;
end

    function text_handle = label2D()
        xl = get(ax, 'xlim');
        if strcmp(get(ax, 'xdir'), 'normal')
            x_left = min(xl);
        else
            x_left = max(xl);
        end
        
        yl = get(ax, 'ylim');
        if strcmp(get(ax, 'ydir'), 'normal')
            y_top = max(yl);
        else
            y_top = min(yl);
        end
        
        x_left = x_left - shift_percent_x * diff(xl);
        y_top = y_top + shift_percent_y * diff(yl);
        text_handle = text(double(x_left), double(y_top), label, 'color', fontcolor, 'fontweight', fontweight, 'fontsize', fontsize, 'parent', ax);
    end


    function text_handle = label3D()
        % Determining the upper-left corner of the graph from the
        % perspective of the camera is more complicated when it is a 3D
        % graph. By taking the cross product of the camera's up direction
        % and the ray from the camera to it's focal target, we get a vector
        % that points left from the camera's point of view (as long as we
        % are in a right-handed coordinate system).
        %
        % If an odd number of the axes are reversed, than this graph is
        % a left-handed coordinate system, so we would need to flip the
        % cross product vector to get it pointed in the right direction.
        left_handed = mod(strcmp(get(ax, 'xdir'), 'reverse') + strcmp(get(ax, 'ydir'), 'reverse') + strcmp(get(ax,'zdir'),'reverse'),2) == 1;
        
        camera_up = get(ax,'CameraUpVector');
        camera_out = get(ax,'CameraTarget') - get(ax,'CameraPosition');
        left_vec = cross(camera_up, camera_out);
        left_vec = (left_vec / norm(left_vec)) .* shift_percent_x;
        if left_handed
            left_vec = -left_vec;
        end
        
        % Once we have the left-pointing vector, we want to choose the x
        % and y limits that are in the direction pointed by that vector, so
        % if the x component is negative, we need the lower x-limit and
        % vice versa. For the z-direction, we use the same trick with the
        % camera's up vector (since that is usually [0 0 1] or [0 0 -1],
        % the left pointing vector usually does not have a nonzero z
        % component).
        
        xl = get(ax, 'xlim');
        if left_vec(1) > 0
            x_left = max(xl);
        else
            x_left = min(xl);
        end
        
        yl = get(ax, 'ylim');
        if left_vec(2) > 0
            y_left = max(yl);
        else
            y_left = min(yl);
        end
        
        zl = get(ax,'zlim');
        if camera_up(3) > 0
            z_top = max(zl);
        else
            z_top = min(zl);
        end
        
        text_pos = [x_left, y_left, z_top] + left_vec + camera_up .* shift_percent_y;
        text_handle = text(double(text_pos(1)), double(text_pos(2)), double(text_pos(3)), label, 'color', fontcolor, 'fontweight', fontweight, 'fontsize', fontsize, 'parent', ax);
    end
end




