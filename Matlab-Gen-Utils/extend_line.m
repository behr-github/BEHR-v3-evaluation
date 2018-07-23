function [  ] = extend_line( line_handle, x_lims )
%EXTEND_LINE Lengthen a plotted line
%   EXTEND_LINE( LINE_HANDLE, X_LIMS ) Given a handle to a line,
%   LINE_HANDLE, this will attempt to extend the line to new X_LIMS, which
%   should be a 2-element vector specifying the minimum and maximum x
%   values. This function uses linear extrapolation to extend the line to
%   reach those limits. It will not add new points to the interior of the
%   line if any value of X_LIMS is within the range of x values present in
%   the line.

E = JLLErrors;

if ~ishandle(line_handle) || ~isscalar(line_handle)
    E.badinput('LINE_HANDLE must be a scalar handle to a line')
end
if ~isnumeric(x_lims) || numel(x_lims) ~= 2
    E.badinput('X_LIMS must be a 2 element numeric vector');
end


line_x = get(line_handle,'XData');
line_y = get(line_handle,'YData');

% Avoid adding new points to the interior of the line
xx = x_lims < min(line_x(:)) | x_lims > max(line_x(:));

new_x = sort(cat(2, line_x, x_lims(xx)));
new_y = interp1(line_x, line_y, new_x, 'linear', 'extrap');

set(line_handle, 'XData', new_x, 'YData', new_y);

end

