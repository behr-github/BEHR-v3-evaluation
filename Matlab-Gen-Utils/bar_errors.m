function [  ] = bar_errors( bar_vals, error_vals, varargin )
%BAR_ERRORS Places error bars on the simplest bar graph possible.
%   BAR_ERRORS( BAR_VALS, ERROR_VALS ) will generate a bar graph using the
%   syntax BAR( BAR_VALS ). Symmetrical error bars will be placed on the
%   bars ends with spread specified by ERROR_VALS. Error bars will be
%   black.
%
%   BAR_ERRORS( BAR_VALS, ERROR_LOWER, ERROR_UPPER ) will put asymmetrical
%   error bars on the plot.


E = JLLErrors;
if ~isequal(size(bar_vals), size(error_vals))
    E.badinput('bar_vals and error_vals must be the same size')
end

B = bar(bar_vals);

%error_vals = reshape(error_vals,1,[]);
error_x = nan(size(error_vals));
error_y = nan(size(error_vals));

for a=1:numel(B)
    error_x(:,a) = B(a).XData + B(a).XOffset;
    error_y(:,a) = B(a).YData;
end

scatter_errorbars(error_x, error_y, error_vals, varargin{:}, 'color', 'k', 'linewidth', 2, 'tipscale', 2);

end

