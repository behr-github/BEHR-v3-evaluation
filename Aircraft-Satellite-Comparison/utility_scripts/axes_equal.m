function [  ] = axes_equal( varargin )
%AXES_EQUAL Set axis limits equal
%   AXES_EQUAL() will set the current axes to have the same limits. Those
%   limits will be calculated automatically as the nearest multiple of a
%   power of 10 (the lower limit will always be 0)
%
%   AXES_EQUAL(val) will set both axes in the current axes to the upper
%   limit given as val.  If val is a 2-element vector, that will be used as
%   the upper and lower limits.
%
%   AXES_EQUAL(ax,___) will use the axis handle specified. val may or may
%   not be passed.

%%%%% INPUT PARSING %%%%%
E = JLLErrors;

% If 2 arguments in, assign both inputs to the appropriate variables.
% If 1 argument, check whether it is a number; if it is, it's the limits if
% not, it's an axis handle.
% If 0 arguments, set both automatically.
if nargin == 2
    ax = varargin{1};
    val = varargin{2};
    if numel(val)==1
        val = [0, val];
    end
elseif nargin == 1
    if isnumeric(varargin{1})
        ax = gca;
        val = varargin{1};
        if numel(val)==1
            val = [0, val];
        end
    else
        ax = varargin{1};
        val = setLimits(ax);
    end
elseif nargin == 0
    ax = gca;
    val = setLimits(ax);
end

% Check that both args are what we expect: val should be a 2-element vector
% and ax should be an axis object
if ~isgraphics(ax)
    error(E.badinput('ax must be an Axes object. It must be the first argument.'));
elseif ~isnumeric(val) || numel(val) ~= 2
    error(E.badinput('val must be a number, and must be 1 or 2 elements only'));
end

% Now set both axes to the desired limits
ax.XLim = val;
ax.YLim = val;

end


function val = setLimits(ax)
    upper_lim = max([ax.XLim, ax.YLim]);
    pow = floor(log10(upper_lim));
    num = ceil(upper_lim / 10^pow);
    val = [0, num * 10^pow];
end

