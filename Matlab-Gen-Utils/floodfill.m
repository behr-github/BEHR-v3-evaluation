function [in, edge] = floodfill(val, x, y, threshold, center_x, center_y)
% FLOODFILL Find a contigous section of a matrix that meets some threshold
%   [IN, EDGE] = FLOODFILL( VAL, X, Y, THRESHOLD, CENTER_X, CENTER_Y )
%   FLOODFILL with explore the values given in the VAL matrix and find
%   points that meet the given THRESHOLD. X and Y define the x and y
%   coordinates of VAL and so must have the same size. The algorithm will
%   start in the box closest to CENTER_X and CENTER_Y and work outwards,
%   stopping when no new neighbor values meet the given THRESHOLD.
%   THRESHOLD can either be a scalar number, in which case values meet the
%   threshold if they are greater than or equal to the value given as
%   THRESHOLD. THRESHOLD can also be a function handle, which is passed the
%   difference between the current value and the center one (defined by
%   center_x and center_y). It must return a scalar logical value, true if
%   the current pixel meets the threshold.
%
%   The return value IN is a logical matrix the same size as VAL that is
%   true for all pixels part of a contigous block including the center
%   pixel that meet the threshold. EDGE is a logical matrix of the same
%   size, but only true for values on the edge of the block (i.e. pixels
%   that have at least one neighbor that does not meet the threshold).
%   Pixels are considered to have 8 neighbors.

%%%%%%%%%%%%%%%%%%
% INPUT CHECKING %
%%%%%%%%%%%%%%%%%%
E = JLLErrors;

if ~isnumeric(val) || ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(center_x) || ~isnumeric(center_y)
    E.badinput('VAL, X, Y, CENTER_X, and CENTER_Y must all be numeric')
elseif ~ismatrix(val) || ~ismatrix(x) || ~ismatrix(y)
    E.badinput('VAL, X, and Y must be, at most 2D matrices. Higher dimensional arrays are not supported.')
elseif ~isequal(size(val), size(x)) || ~isequal(size(val), size(y))
    E.badinput('VAL, X, and Y must all be the same size')
elseif ~isscalar(center_x) || ~isscalar(center_y)
    E.badinput('CENTER_X and CENTER_Y must be scalar values')
end

if isnumeric(threshold)
    if ~isscalar(threshold)
        E.badinput('If given as a number, THRESHOLD must be scalar')
    end
elseif isa(threshold, 'function_handle')
    try
        threshold_check = threshold(0);
    catch err
        E.badinput('The function given as THRESHOLD failed when called with an input of 0. The error message is: "%s"', err.message)
    end
    if ~islogical(threshold_check) || ~isscalar(threshold_check)
        E.badinput('THRESHOLD must return a scalar logical value when given a scalar numeric input')
    end
else
    E.badinput('THRESHOLD must be a number or function handle')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find the pixel corresponding to the center x & y %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the difference between all of the x/y values and the center
% x/y, then find the smallest total difference
lat_dif = abs(y - center_y);
lon_dif = abs(x - center_x);
combined_dif = lat_dif + lon_dif;

%The linear index of the center pixel
[~, center_pix] = min(combined_dif(:));

%if the threshold is numeric, redefine it as a anonymous function
%As in the original code, the 'if' statement in line 57 is
%matrix_in(n_row(b), n_col(b)) >= threshold. To keep it as the same format,
%the threshold is defined as @(t) t>=th-matrix_in(center_pix);
if isnumeric(threshold)
    th = threshold;
    threshold = @(t) t>=th-val(center_pix);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialize our starting condition and loop until the plume is located %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plume_mat = zeros(size(val));
plume_mat(center_pix) = 1;

while true
    [row, col] = find(plume_mat == 1);
    if isempty(row)
        %If no pixels with a value of 1 are left, assume that the full
        %plume has been found and exit.
        break; 
    end 
    for a = 1:length(row)
        %Take each of the eight neighbors and test: (a) if they are already
        %part of the plume and (b) if they are higher than the threshold
        n_row = [row(a) - 1, row(a) - 1, row(a) - 1, row(a), row(a) row(a) + 1, row(a) + 1, row(a) + 1];
        n_col = [col(a) - 1, col(a), col(a) + 1, col(a) - 1, col(a) + 1, col(a) - 1, col(a), col(a) + 1];
        for b = 1:8
            if n_row(b) < 1 || n_row(b) > size(plume_mat,1) || n_col(b) < 1 || n_col(b) > size(plume_mat,2) 
                % edge case. Cannot evaluate b/c requested grid cell does
                % not exist
            elseif plume_mat(n_row(b), n_col(b)) > 0 
                %A 0 indicates this cell has not yet been visited.  If the
                %value is >0, skip this cell because it has already been
                %considered
            elseif threshold(val(n_row(b), n_col(b))-val(center_pix))
                %A value of 1 indicates that this pixel is part of the
                %plume but its neighbors have not yet been checked.
                plume_mat(n_row(b), n_col(b)) = 1; 
                %By setting this to 1, the loop will catch it next time and
                %check its neighbors.
            end                              
        end
        %A value of 2 indicates that the neighbors of this pixel have been
        %evaluated.
        plume_mat(row(a),col(a)) = 2; 
    end
end

in = plume_mat == 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Now find the edge pixels, i.e. any with a neighbor marked 0 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row2, col2] = find(plume_mat == 2);
edge_pixels_ind = zeros(1,numel(row2)); edge_pixels_ind(:)=NaN;
edge_ind = 1;
for a = 1:length(row2)
    n_row = [row2(a) - 1, row2(a) - 1, row2(a) - 1, row2(a), row2(a) row2(a) + 1, row2(a) + 1, row2(a) + 1];
    n_col = [col2(a) - 1, col2(a), col2(a) + 1, col2(a) - 1, col2(a) + 1, col2(a) - 1, col2(a), col2(a) + 1];
    % Edge cases - restrict to valid indices
    n_row = max(min(n_row,size(plume_mat,1)),1);
    n_col = max(min(n_col,size(plume_mat,2)),1);
    plume_mat_sub = plume_mat(n_row,n_col);
    
    for b=1:8
        if any(plume_mat_sub(:)==0)
            lin_ind = sub2ind(size(plume_mat),row2(a),col2(a));
            edge_pixels_ind(edge_ind) = lin_ind;
            edge_ind = edge_ind + 1;
        end
    end
end
edge_pixels_ind(isnan(edge_pixels_ind)) = [];
edge = false(size(plume_mat));
edge(edge_pixels_ind) = true;
end

