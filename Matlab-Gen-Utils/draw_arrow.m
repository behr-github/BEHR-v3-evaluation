function [ ] = draw_arrow( x, y, varargin )
%DRAW_ARROW( X, Y ) Draws an arrow on the current axes
%   Uses quiver to draw the arrow starting at X(1) Y(1) and ending at X(2)
%   Y(2).
%
%   DRAW_ARROW( X, Y, name-value pairs ) allows you to pass name value
%   pairs through to QUIVER to alter the appearance of the arrows.
%
%   Credit to http://stackoverflow.com/questions/25729784/how-to-draw-an-arrow-in-matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if numel(x) ~= 2 || numel(y) ~= 2 || ~isnumeric(x) || ~isnumeric(y)
    E.badinput('X and Y must be 2-element numeric vectors');
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if plot is held currently; will reset this state afterwards
pltstate = get(gca,'nextplot');
finishup = onCleanup(@() myCleanup(pltstate));
set(gca,'nextplot','add'); % do not overwrite the current plot
quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0,varargin{:});

end

function myCleanup(pltstate)
% Ensure that the plot state is reset even if this script errors
set(gca,'nextplot',pltstate);
end