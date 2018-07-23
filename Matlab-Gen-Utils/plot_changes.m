function [ varargout ] = plot_changes(x, y, varargin)
%PLOT_CHANGES Plot a series of changes connected by lines
%   PLOT_CHANGES( X, Y ) plot each column of X and Y as a separate series
%   with different marker. Each row of X and Y are connected by lines to
%   represent how a pair of values changes over time (or any other
%   coordinate).
%
%   SERIES = PLOT_CHANGES( ___ ) With one output, this returns an array of
%   graphics objects corresponding to each series (i.e. each column of X
%   and Y).
%
%   Additional parameters:
%
%       'group_fmts' - a structure with field names that can be passed as
%       formatting options to LINE. The structure must have the same number
%       of elements as X and Y have columns in order to give each series
%       its own format. If 'linestyle' is not a field of the structure,
%       then it defaults to 'none' (this does not affect the connecting
%       lines).
%
%       'parent' - a axes object to plot into. If omitted, the current axes
%       are used. If no figure is open, a new set of axes is created.

E = JLLErrors;

p = inputParser;
p.addParameter('group_fmts',[]);
p.addParameter('connector_fmt', struct('color','k'));
p.addParameter('parent', gobjects(0));

p.parse(varargin{:});
pout = p.Results;

group_fmts = pout.group_fmts;
connector_fmt = pout.connector_fmt;
ax = pout.parent;

if ~isnumeric(x) || ~isnumeric(y)
    E.badinput('X and Y must be numeric');
elseif ~isequal(size(x), size(y))
    E.badinput('X and Y must be the same size');
end

n_groups = size(x,2);

if isempty(group_fmts)
    group_fmts = make_default_formats(n_groups);
else
    group_fmts = validate_formats(group_fmts, n_groups);
end

if ~isstruct(connector_fmt)
    E.badinput('"connector_fmt" must be a structure');
end

if isempty(ax)
    ax = gca;
elseif ~isvalid(ax) 
    E.badinput('The value given to parameter "ax" must be a valid handle');
elseif ~strcmp(ax.Type, 'axes')
    E.badinput('The value given to parameter "ax" must be a handle to axes');
end

% Create groups to hold the individual series and the arrows
series_group = hggroup(ax);
arrows_group = hggroup(ax);

% Plot the series 
for i_grp=1:n_groups
    line(x(:,i_grp), y(:,i_grp), group_fmts(i_grp), 'parent', series_group);
end

% Add the connecting lines
for i_change = 1:size(x,1)
    line(x(i_change, :), y(i_change, :), connector_fmt, 'parent', arrows_group);
end

% Return the handles if requested. Graphic objects are added in
% last-in-first-out order
if nargout > 0
    varargout{1} = flipud(series_group.Children);
end

end

function fmt = make_default_formats(n_groups)
markers = {'o','x','^','v','*'};
linestyle = 'none';

rep_size = ceil(numel(markers)/n_groups);
markers = repmat(markers, 1, rep_size);
markers = markers(1:n_groups);

fmt = struct('marker', markers, 'linestyle', linestyle);

end

function fmt = validate_formats(fmt, n_groups)
E = JLLErrors;
if ~isstruct(fmt)
    E.badinput('The parameter "group_fmts" must be a structure');
elseif numel(fmt) ~= n_groups
    E.badinput('The parameter "group_fmts" must have a number of elements equal to the number of columns in X and Y');
end

% If linestyle is not part of the format, set it to "none" so that lines
% are only given connecting rows of X and Y.
if ~isfield(fmt, 'linestyle')
    for i_group = 1:n_groups
        fmt(i_group).linestyle = 'none';
    end
end
end