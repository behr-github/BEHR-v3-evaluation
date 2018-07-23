function [ col ] = map2colmap( val, minval, maxval, map )
%MAP2COLMAP Maps a values to a color
%   Maps a values to a color in a colormap similarly to how Matlab
%   internally does it (not guaranteed to be identical). Specifically, it
%   takes a range of values (default 0 to 1, but can be specified) and
%   makes the bottom of the range the first color in the colormap, the top
%   of the range the last value in the colormap, and linearly interpolates
%   in between. Values out of range are clipped to the colormap.
%
%   COL = MAP2COLMAP( VAL ) Maps VAL to the parula colormap, assuming that
%   VAL lies on the range 0 to 1.
%
%   COL = MAP2COLMAP( VAL, MAXVAL ) Maps VAL to the parula colormap,
%   assuming that VAL lies on the range 0 to MAXVAL (if MAXVAL > 0) or
%   MAXVAL to 0 (if MAXVAL < 0).
%
%   COL = MAP2COLMAP( VAL, MINVAL, MAXVAL ) Maps VAL to the parula colormap
%   for VAL on the range MINVAL to MAXVAL.
%
%   COL = MAP2COLMAP( ___, MAP ) Uses colormap MAP instead of parula with
%   any of the previous syntaxes.

narginchk(1,4);
E = JLLErrors;
[val, minval, maxval, map] = parse_input(nargin);

figtmp = figure; % colormap always opens a figure
cmap = colormap(map);
close(figtmp);

ncolors = size(cmap,1);
ind = round((val - minval)/(maxval - minval) * (ncolors-1) + 1);
ind = max(ind, 1);
ind = min(ind, ncolors);

col = cmap(ind,:);

    function [val_in, minval_in, maxval_in, map_in] = parse_input(nargs)
        def_map = 'parula'; 
        
        val_in = val;
        if nargs == 1
            minval_in = 0;
            maxval_in = 1;
            map_in = def_map;
        elseif nargs == 2
            minval_in = 0;
            if is_colormap(minval)
                maxval_in = 1;
                map_in = minval;
            elseif isnumeric(minval) && isscalar(minval);
                maxval_in = minval;
                map_in = def_map;
            else
                E.badinput('With 2 inputs, the second must be a scalar number (maximum value), an n-by-3 array (color map) or string (color map)');
            end
        elseif nargs == 3
            if ~isnumeric(minval) || ~isscalar(minval)
                E.badinput('With 3 inputs, the second must be a scalar number')
            end
            
            if is_colormap(maxval)
                map_in = maxval;
                minval_in = 0;
                maxval_in = minval;
            elseif isnumeric(maxval) && isscalar(maxval)
                map_in = def_map;
                minval_in = minval;
                maxval_in = maxval;
            else
                E.badinput('With 3 inputs, the third must be a scalar number (maximum value), an n-by-3 array (color map) or string (color map)');
            end
        elseif nargs == 4
            minval_in = minval;
            maxval_in = maxval;
            map_in = map;
        else
            E.notimplemented('nargin == %d', nargs);
        end
        
        if ~isnumeric(val_in) || ~isscalar(val_in)
            E.badinput('VAL must be a scalar number');
        elseif ~isnumeric(minval_in) || ~isscalar(minval_in)
            E.badinput('Failed to parse input: MINVAL is not a scalar number')
        elseif ~isnumeric(maxval_in) || ~isscalar(maxval_in)
            E.badinput('Failed to parse input: MAXVAL is not a scalar number')
        elseif ~is_colormap(map_in)
            E.badinput('Failed to parse input: MAP is not a colormap (string or n-by-3 array')
        end
        
        if minval_in > maxval_in
            tmp = minval_in;
            minval_in = maxval_in;
            maxval_in = tmp;
        end
    end

end

function b = is_colormap(map)
b = ischar(map) || (isnumeric(map) && size(map,2) == 3);
end