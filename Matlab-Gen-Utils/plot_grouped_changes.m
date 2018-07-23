function [ varargout ] = plot_grouped_changes(changes, varargin)
%PLOT_GROUPED_CHANGES Plot changes between initial and final states in one or more groups
%   PLOT_GROUPED_CHANGES( CHANGES ) Given CHANGES, an N-by-2 array or cell
%   array of N-by-2 arrays, plot the changes connected by arrows. Each row
%   of CHANGES is two points that should be connected, the rows will be
%   plotted successively along the x-axis. Examples:
%
%       plot_grouped_changes([1 3; 2 4; 3 5])
%
%   would plot points (1,1) and (1,3) connected by an arrow, then points
%   (2, 2) and (2, 4), then (3, 3) and (3, 5). If CHANGES is a cell array
%   of arrays, then all the first rows of every array in the cell array
%   would be clustered together, then all the second rows, and so on.
%
%   There are a number of parameters that control the appearance of the
%   points. 
%
%       'common_fmt' - format that should be applied to all points.
%
%       'group_fmt' - formats that should be applied to different groups.
%
%       'series1_fmt' and 'series2_fmt' - formats applied to the two series
%       (before and after, i.e. CHANGES(:,1) and CHANGES(:,2) if CHANGES is
%       an array itself and not a cell array).
%
%       'arrow_fmt' - the format of the arrows.
%
%   All of these parameters are expected to be structures that if passed as
%   the final argument to LINE() would set the visual characteristics of
%   the line. For example, to make all points be black circles with no
%   connecting line, then pass
%
%       struct('color', 'k', 'marker', 'o', 'linestyle', 'none')
%
%   as the value for 'common_fmt'. Essentially, the field names are the
%   parameter names and the field values the parameter values. All of the
%   values are expected to be scalar structures except 'group_fmt', which
%   can be non-scalar. For 'group_fmt', group_fmt(i) is applied as the
%   format for CHANGES{i}, i.e. 'group_fmt' gives you a way to specify how
%   each group's formatting differs. If there are more arrays in CHANGES
%   than there are indicies in 'group_fmt', then 'group_fmt' is repeated;
%   so if 'group_fmt' is a 1x3 structure and CHANGES is a 1x4 cell array,
%   CHANGES{4} will be formatted using group_fmt(1).
%
%   These are applied in the order common_fmt, group_fmt, series1_fmt or
%   series2_fmt. That is, any properties set in group_fmt override those
%   set in common_fmt and any set in series1_fmt or series2_fmt override
%   both common_fmt and group_fmt. The exception is for the default values
%   (see below). arrow_fmt is independent of the other three.
%
%   The default values for each are:
%       common_fmt = struct('linestyle', 'none')
%       group_fmt = struct('marker', {'o','x','^','p','v'});
%       series1_fmt = struct('color','b');
%       series2_fmt = struct('color', 'r');
%       arrow_fmt = struct('color','k')
%   i.e. by default, there are no lines connecting the markers, group 1
%   will use circles, group 2 X's, etc., and the first series will be red
%   and the second blue. For group_fmt, remember that giving a cell array
%   as the value in struct() creates a non-scalar structure.
%
%   However, the defaults do not override user preferences, so if you pass
%   struct('linestyle','none','color','k') for common_fmt, then both series
%   formats will not specify color. Only if you pass a structure with the
%   field 'color' to either group_fmt, series1_fmt, or series2_fmt will the
%   color in common_fmt be overridden.
%
%
%   There are further parameters to control labelling:
%
%       'tick_labels' - sets the labels for the X ticks. The X ticks are
%       automatically set to be the center of each cluster of changes. Must
%       be a cell array of strings.
%
%       'group_labels' - how to label each group in the legend. Must be a
%       cell array of strings.
%
%       'series_labels' - how to label each series in the legend. Must be a
%       1-by-2 cell array of strings.
%
%   If one or both of group_labels and series_labels is/are given, then a
%   legend is automatically added. Dummy lines are used to create legend
%   entries that are common_fmt + group_fmt (one per group) and common_fmt
%   + series_fmt (one per series).

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

E = JLLErrors;

p = inputParser;
p.addParameter('common_fmt', struct('linestyle','none'));
p.addParameter('group_fmt',struct());
p.addParameter('series1_fmt',struct());
p.addParameter('series2_fmt',struct());
p.addParameter('arrow_fmt', struct('color', 'k'));
p.addParameter('tick_labels', {});
p.addParameter('group_labels', {});
p.addParameter('series_labels', {});
p.addParameter('intra_space', 0.5);
p.addParameter('inter_space', 1);

p.parse(varargin{:});
pout = p.Results;

common_fmt = pout.common_fmt;
group_fmt = pout.group_fmt;
series1_fmt = pout.series1_fmt;
series2_fmt = pout.series2_fmt;
arrow_fmt = struct2opts(pout.arrow_fmt);

tick_labels = pout.tick_labels;
group_labels = pout.group_labels;
series_labels = pout.series_labels;

intra_group_space = pout.intra_space;
inter_group_space = pout.inter_space;

if isnumeric(changes)
    changes = {changes};
elseif ~iscell(changes)
    E.badinput('CHANGES must be an n-by-2 numeric array or a cell array of such arrays')
end

changes_check = find(~cellfun(@(x) isnumeric(x) && size(x,2) == 2, changes));
if ~isempty(changes_check)
    E.badinput('The arrays at indices %s of CHANGES are either not numeric or not n-by-2', mat2str(changes_check));
end

% To avoid overwriting user input, we need to make sure that the user
% doesn't pass the color or marker options in common_fmt or the color
% option in group_fmt before we set the defaults.
if numel(fieldnames(group_fmt)) == 0 && ~isfield(common_fmt, 'marker')
        group_fmt = struct('marker', {'o','x','^','p','v'});
end
if ~isfield(common_fmt, 'color') && ~isfield(group_fmt, 'color')
    if numel(fieldnames(series1_fmt)) == 0
        series1_fmt = struct('color','b');
    end
    if numel(fieldnames(series2_fmt)) == 0
        series2_fmt = struct('color', 'r');
    end
    
end

[format_array, legend_groups, legend_series] = design_format_array(numel(changes), common_fmt, group_fmt, series1_fmt, series2_fmt);

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

% First figure out how to plot the groups.
[group_ticks, central_ticks] = calculate_x_coords(changes, intra_group_space, inter_group_space);

n_groups = numel(changes);

% Now plot the individual series;
fig = figure;

for i_group = 1:n_groups
    line(group_ticks{i_group}, changes{i_group}(:,1), format_array(i_group, 1));
    line(group_ticks{i_group}, changes{i_group}(:,2), format_array(i_group, 2));
    for i_tick = 1:numel(group_ticks{i_group})
        draw_arrow(repmat(group_ticks{i_group}(i_tick),1,2), changes{i_group}(i_tick,:), arrow_fmt{:});
    end
end

xmin = min(cellfun(@(x) min(x(:)), group_ticks)) - inter_group_space;
xmax = max(cellfun(@(x) max(x(:)), group_ticks)) + inter_group_space;
set(gca,'xtick',central_ticks,'xlim', [xmin, xmax]);
if ~isempty(tick_labels)
    set(gca,'xticklabel',tick_labels);
end


% Finally make the legend. Use dummy lines that won't show up anywhere to
% create legend entries that only have the relevant pieces to the groups
% vs. series.
if ~isempty(group_labels)
    group_handles = gobjects(n_groups,1);
    for i_group = 1:n_groups
        group_handles(i_group) = line(nan,nan,legend_groups(i_group));
    end
    
    % Ensure the legend labels are a row
    legend_labels = group_labels(:)';
    
    % If we also need to put series labels in, use this to add a blank line
    % in the legend
    if ~isempty(series_labels)
        group_handles(end+1) = line(nan,nan,'linestyle','none','marker','none');
        legend_labels{end+1} = '';
    end
else
    group_handles = gobjects(0);
    legend_labels = cell(0);
end

% If only one group is passed, then adding the extra "blank" group handle
% will make a row vector instead of a column vector.
if isrow(group_handles)
    group_handles = group_handles';
end

if ~isempty(series_labels)
    series_handles = gobjects(2,1);
    for i_series = 1:2
        series_handles(i_series) = line(nan,nan,legend_series(i_series));
    end
    
    group_handles = cat(1, group_handles, series_handles);
    legend_labels = cat(2, legend_labels, series_labels);
end

if ~isempty(group_handles)
    legend(group_handles, legend_labels);
end

if nargout > 0
    varargout{1} = fig;
end

end

function [group_ticks, central_ticks] = calculate_x_coords(changes, intra_group_space, inter_group_space)
n_groups = numel(changes);
max_ticks = max(cellfun(@(x) size(x,1), changes));

% Each group will by default have 0.5 X units between its members, and the
% edge members of different groups will be separated by 1 unit by default.
% For even numbers of groups, the central tick will go in the space between
% the middle two groups; for odd numbers it will be aligned with the middle
% group.

group_width = intra_group_space * (n_groups - 1);
xtick_spacing = group_width + inter_group_space;

central_ticks = 1:xtick_spacing:(1 + xtick_spacing*(max_ticks-1));
group_ticks = cell(n_groups, 1);
offsets = -(n_groups-1)/2*intra_group_space:intra_group_space:(n_groups-1)/2;
for i_group = 1:n_groups
    group_ticks{i_group} = central_ticks + offsets(i_group);
end

end

function [format_array, group_examples, series_examples] = design_format_array(n_groups, common_fmt, group_fmt, series1_fmt, series2_fmt)
% Create a struct array (n_groups-by-2) that gives the formatting options for each
% data point based on it's group and which series it is in (before or after
% the change). Series formatting overrides group formatting.

E = JLLErrors;

if ~isstruct(common_fmt) || ~isscalar(common_fmt)
    E.badinput('The value for "common_fmt" must be a struct');
end
if ~isstruct(group_fmt)
    E.badinput('The value for "group_fmt" must be a struct');
end
if ~isstruct(series1_fmt) || ~isscalar(series1_fmt)
    E.badinput('The value for "series1_fmt" must be a struct');
end
if ~isstruct(series2_fmt) || ~isscalar(series2_fmt)
    E.badinput('The value for "series2_fmt" must be a struct');
end

all_fields = unique(veccat(fieldnames(common_fmt), fieldnames(group_fmt), fieldnames(series1_fmt), fieldnames(series2_fmt)));
format_array = repmat(make_empty_struct_from_cell(all_fields), n_groups, 2);

do_err_empty_field = true;

format_array(:) = copy_fields(common_fmt, format_array(1), do_err_empty_field, 'The field "%s" in common_fmt is empty');

% Prepare example series that only have the relevant characteristics for
% each group and the two series for use in the legend.

group_fields = unique(veccat(fieldnames(common_fmt), fieldnames(group_fmt)));
series_fields = unique(veccat(fieldnames(common_fmt), fieldnames(series1_fmt), fieldnames(series2_fmt)));

group_examples = repmat(make_empty_struct_from_cell(group_fields),n_groups,1);
group_examples(:) = copy_fields(common_fmt, group_examples(1));

series_examples = repmat(make_empty_struct_from_cell(series_fields),2,1);
series_examples(:) = copy_fields(common_fmt, series_examples(1));
series_examples(1) = copy_fields(series1_fmt, series_examples(1));
series_examples(2) = copy_fields(series2_fmt, series_examples(2));

for i_group = 1:n_groups
    i_g_wrapped = mod(i_group - 1, numel(group_fmt)) + 1;
    format_array(i_group, :) = copy_fields(group_fmt(i_g_wrapped), format_array(i_group, 1), do_err_empty_field, sprintf('The field "%%s" in group_fmt(%d) is empty', i_g_wrapped));
    
    format_array(i_group, 1) = copy_fields(series1_fmt, format_array(i_group, 1), do_err_empty_field, 'The field "%s" in series1_fmt is empty');
    format_array(i_group, 2) = copy_fields(series2_fmt, format_array(i_group, 2), do_err_empty_field, 'The field "%s" in series2_fmt is empty');
    
    group_examples(i_group) = copy_fields(group_fmt(i_g_wrapped), group_examples(i_group));
end

group_examples = ensure_visible(group_examples);
series_examples = ensure_visible(series_examples);

end

function opts = ensure_visible(opts)
% Make sure that the given options structure will result in a visible
% series
for i=1:numel(opts)
    if isfield(opts,'linestyle') && strcmp(opts(i).linestyle, 'none')
        if ~isfield(opts,'marker') || strcmp(opts(i).marker, 'none')
            opts(i).linestyle = '-';
        end
    end
end
end

function dest = copy_fields(source, dest, error_if_empty_field, err_msg)
% If source and dest are structs with different fields, then source(i) =
% dest will fail. This should avoid that, so long as all fields in dest are
% already in source, which is how format_array is constructed in
% design_format_array.

E = JLLErrors;

if nargin < 3
    error_if_empty_field = false;
end
if nargin < 4
    err_msg = 'The field "%s" in SOURCE was empty';
end

fns = fieldnames(source);
for f=1:numel(fns)
    if isempty(source.(fns{f}))
        if error_if_empty_field
            E.badinput(err_msg, fns{f});
        end
    else
        dest.(fns{f}) = source.(fns{f});
    end
end
end

function opts = struct2opts(S)
% Unfortunately, QUIVER() cannot handle options as a structure so we need
% to convert the options structure to name-value pairs
fns = fieldnames(S);
opts = cell(1,2*numel(fns));
for f = 1:2:numel(opts)
    opts{f} = fns{f};
    opts{f+1} = S.(fns{f});
end
end