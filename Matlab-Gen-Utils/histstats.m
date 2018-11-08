function histstats(varargin)
%HISTSTATS Plot a histogram with statistics included on it
%   HISTSTATS( ___ ) Plot a histogram of values, with the mean marked with
%   a vertical line and the mean, median, standard deviation, and upper and
%   lower quartiles marked on the plot as a text object. The text will
%   automatically try to find the more open side of the plot to be on. This
%   can accept any call format that HIST can, and will return N and X from
%   HIST if outputs exist. Note that VALS must be a vector for this
%   function, unlike HIST.

% Allow for the case where the user passes a handle as the first argument.
if isnumeric(varargin{1})
    vals = varargin{1};
elseif isnumeric(varargin{2})
    vals = varargin{2};
else
    E.badinput('The values must be given as the first argument, or the second if the first argument is an axis handle.')
end

if ~isvector(vals)
    E.badinput('The values must be given as a vector. Plotting a matrix of values is not supported.')
end

% Extract our specific options
[ our_args, hist_args ] = parse_args(varargin{:});

% Prepare the stats text.
stats_text = sprintf('Mean = %.4g\n1\\sigma = %.4g\nMedian = %.4g\nQuartiles = %.4g, %.4g', nanmean(vals), nanstd(vals), nanmedian(vals), quantile(vals, 0.25), quantile(vals, 0.75));

hist(hist_args{:});

avg_diff = nanmean(vals);

% Currently, I assume that all the extra arguments to this function are
% formatting arguments.
fns = fieldnames(our_args);
for i_fn = 1:numel(fns)
    if ~isempty(our_args.(fns{i_fn}))
        set(gca, fns{i_fn}, our_args.(fns{i_fn}));
    end
end

% Get these after changing the font size because that can
% change the limits (in order to make the ticks fit better
% visually)
x_lim_vals = get(gca,'xlim');
y_lim_vals = get(gca,'ylim');

l=line([avg_diff, avg_diff], y_lim_vals, 'color', 'r', 'linestyle', '--', 'linewidth', 2);
legend(l,{'Mean'});

% Put the text near the top left or right corner with a little
% space between it and the edge of the plot. Choose the corner
% based on which side is further from zero, on the assumption
% that the histogram will peak near 0 so going further from 0
% gives the text more room.
if abs(x_lim_vals(1)) > abs(x_lim_vals(1))
    text_x = x_lim_vals(1) + 0.02 * diff(x_lim_vals);
    text_y = y_lim_vals(2) - 0.1 * diff(y_lim_vals);
    h_align = 'left';
else
    text_x = x_lim_vals(2) - 0.02 * diff(x_lim_vals);
    % leave a little extra space for the legend on the right
    % side
    text_y = y_lim_vals(2) - 0.2 * diff(y_lim_vals);
    h_align = 'right';
end
text(text_x, text_y, stats_text, 'verticalalignment', 'top', 'horizontalalignment', h_align, 'fontsize', 10);

end

function [ our_args, hist_args] = parse_args(varargin)
% Have to manually parse arguments because inputParser chokes if only
% expecting parameters but gets positional args
param_struct = struct('fontsize',[]);
fns = fieldnames(param_struct);
our_args = struct();
our_args_idx = false(size(varargin));
for i_fn = 1:numel(fns)
    xx = find(strcmpi(fns{i_fn}, varargin), 1, 'last');
    if ~isempty(xx)
        our_args.(fns{i_fn}) = varargin{xx+1};
        our_args_idx(xx:xx+1) = true;
    end
end

hist_args = varargin(~our_args_idx);

end
