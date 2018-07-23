function [  ] = reset_fig_size( fighnd )
%reset_fig_size Resets the figure to the default 8" x 6".  If not passed a
%figure handle, use the current figure (gcf). The string 'all' will reset
%all currently open figures.


E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the input is a figure handle, if one is given. Otherwise we will
% use the current figure
if nargin < 1
    fig = gcf;
elseif ischar(fighnd) && strcmpi('all',fighnd)
    fig = get(0,'children');
else
    if ~strcmp(class(fighnd),'matlab.ui.Figure');
        error(E.badinput('The function reset_fig_size only accepts figure handles'));
    end
    
    fig = fighnd;
end

% We'll reset the units for the figure to its original setting at the end,
% but to resize it, we want them in inches.  So we'll save the original
% unit so that we can reset it after we've resized the picture.

origunits = get(fig,'Units');
set(fig,'Units','inches');



% Size is saved as the last two values in the position vector.  So we want
% to replace those with 8 and 6 respectively (8" high, 6" wide is the
% default matlab figure size).  We'll keep its position on the screen the
% same, unless, resizing will move it off the screen.

newHeight = 8;
newWidth = 6;

for a = 1:numel(fig)
    oldpos = get(fig(a),'Position');
    newpos = [oldpos(1), oldpos(2), newHeight, newWidth];
    
    % Get the screen size to check if the figure will be offscreen
    origscrunits = get(0,'Units');
    set(0,'Units','inches');
    scrsize = get(0,'ScreenSize');
    set(0,'Units',origscrunits);
    
    if newpos(3) + newHeight > scrsize(3)
        newpos(1) = scrsize(3) - newHeight;
    end
    
    if newpos(4) + newWidth > scrsize(4)
        newpos(2) = scrsize(4) - newWidth;
    end
    
    set(fig(a),'Position',newpos);
    if iscell(origunits)
        set(fig(a),'Units',origunits{a});
    else
        set(fig(a),'Units',origunits);
    end
end
end

