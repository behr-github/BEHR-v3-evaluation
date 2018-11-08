function varargout = scatter_grouped(values, varargin)
%SCATTER_GROUPED Plots clusters of values
%   SCATTER_GROUPED( VALUES ) Given a M-by-N matrix VALUES, plot M clusters
%   of N points, with each point slightly offset from its fellows. This
%   behaves essentially similar to BAR( VALUES ) except using markers
%   instead of bars.
%
%   SCATTER_GROUPED( VALUES, ERRORS )
%   SCATTER_GROUPED( VALUES, LOWER_ERRORS, UPPER_ERRORS ) will
%   automatically add error bars to the points, either symmetric (first
%   case) or asymmetric (second case).
%
%   Additional parameters:
%
%       'group_fmt' - a structure that defines the appearance of each
%       group. The field names must be parameters that LINE() understands
%       to format the series. The i-th group will take on the format given
%       by group_fmt(i). E.g., group_fmt = struct('marker', {'<', '>', '^',
%       'v'}) would cause the groups in each cluster to use the four
%       different triangle markers. Note that 'linestyle' is 'none' by
%       default, unless overridded. Also note that if a the structure given
%       has fewer elements than there are groups, it will be repeated as
%       necessary.
%
%       'error_bar_fmt' - a structure that defines the appearance of the
%       error bars. Its field names must be arguments that
%       SCATTER_ERRORBARS() understands for formatting. Note that the
%       parameter 'direction' of SCATTER_ERRORBARS is set to the same value
%       as 'data_axis' (see below) unless overridden.
%
%       'data_axis' - which axis to plot VALUES along. Default is 'y', i.e.
%       the clusters are spaced out along the x-axis and the values given
%       are plotted along the y-axis. Can also be 'x', which reverses this.

E = JLLErrors;

p = advInputParser;
p.addOptional('lower_errors',[]);
p.addOptional('upper_errors',[]);
p.addParameter('group_fmt', []);
p.addParameter('error_bar_fmt', []);
p.addParameter('data_axis','y');
p.addParameter('width', 0.25);

p.parse(varargin{:});
pout = p.Results;

lower_errors = pout.lower_errors;
upper_errors = pout.upper_errors;
group_fmt = pout.group_fmt;
error_bar_fmt = pout.error_bar_fmt;
data_axis = pout.data_axis;
cluster_width = pout.width;

if isempty(group_fmt)
    group_fmt = struct('marker', {'o','x','^','v','*'});
end
if isempty(error_bar_fmt)
    error_bar_fmt = struct('color', 'k');
end

if ~isscalar(cluster_width) || ~isnumeric(cluster_width)
    E.badinput('"cluster_width" must be a scalar number');
elseif cluster_width <= 0 || cluster_width >= 0.5
    E.badinput('"cluster_width" must be between 0 and 0.5 (exclusive)');
end

% Repeat the group format structure enough to cover all groups needed
n_clusters = size(values,1);
n_groups = size(values,2);

n_reps = ceil(n_groups/numel(group_fmt));
group_fmt = repmat(group_fmt(:), n_reps, 1);
n_reps_eb = ceil(n_groups/numel(error_bar_fmt));
error_bar_fmt = repmat(error_bar_fmt(:), n_reps_eb, 1);

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

% First calculate the coordinates for the groups, assuming that each
% cluster of points will be centered around 1, 2, 3, etc. and will spread
% out to +/- 0.25 (by default) to either side so that there is a gap between clusters.

coords = repmat((1:n_clusters)',1,n_groups);
offsets = linspace(-cluster_width, cluster_width, n_groups);
offsets = repmat(offsets,n_clusters,1);

if strcmpi(data_axis, 'y')
    X = coords + offsets;
    Y = values;
elseif strcmpi(data_axis, 'x')
    X = values;
    Y = coords + offsets;
else
    E.badinput('"data_axis" must be "x" or "y"');
end

% Each group will be it's own series
l = gobjects(n_groups,1);
for i_group = 1:n_groups
    fmt_args = struct2cell2(group_fmt(i_group));
    l(i_group) = line(X(:,i_group), Y(:,i_group), 'linestyle', 'none', fmt_args{:});
    
    fmt_args_eb = struct2cell2(error_bar_fmt(i_group));
    if ~isempty(upper_errors)
        scatter_errorbars(X(:,i_group), Y(:,i_group), lower_errors(:,i_group), upper_errors(:,i_group), 'direction', data_axis, fmt_args_eb{:});
    elseif ~isempty(lower_errors)
        scatter_errorbars(X(:,i_group), Y(:,i_group), lower_errors(:,i_group), 'direction', data_axis, fmt_args_eb{:});
    end
end

% Force the ticks to be present for every cluster so the user can rename
% them easily
set(gca,'xtick', 1:n_clusters);

if nargout > 0
    varargout = {l};
end

end

