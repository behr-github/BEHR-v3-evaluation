function [  ] = xylims( varargin )
%XYLIMS Set X and Y limits simultaneously
%   XYLIMS() Makes the X and Y limits of the current axes the same. The
%   smaller lower bound and greater upper bound are used.
%
%   XYLIMS( LIMITS ) Makes the X and Y limits of the current axes equal to
%   LIMITS
%
%   XYLIMS( AX, ___ ) Apply either of the previous syntaxes to the axes
%   specified by AX instead of the current axis. 
%
%   XYLIMS( FIG, ___ ) Finds the child axes of FIG and acts on those.

E = JLLErrors;

ax = gca;
lims = [];
axes_in = 0;
if numel(varargin) > 0
    if isa(varargin{1},'matlab.graphics.axis.Axes')
        ax = varargin{1};
        axes_in = 1;
        if numel(ax) > 1
            E.badinput('AX must be scalar');
        end
    elseif isa(varargin{1},'matlab.ui.Figure')
        ax = findobj(varargin{1},'type','axes');
        axes_in = 1;
        if numel(ax) > 1
            E.badinput('XLIMS cannot handle subplots (or figures with multiple axes in general) yet')
        end
    end
    if numel(varargin) > 1 || axes_in == 0
        lims = varargin{1+axes_in};
    end
end

if ~isnumeric(lims) || numel(lims) ~= 2 && ~isempty(lims)
    E.badinput('LIMS must be a two-element vector')
elseif ~isempty(lims) && lims(2) < lims(1)
    E.badinput('LIMS(1) must be less than LIMS(2)')
end

if isempty(lims)
    % figure out the right limits to use
    x = get(ax,'xlim');
    y = get(ax,'ylim');
    lims = [min([x,y]), max([x,y])];
end

set(ax,'xlim',lims);
set(ax,'ylim',lims);

end