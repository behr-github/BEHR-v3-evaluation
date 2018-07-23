function [ X_corner, Y_corner ] = pcolor_pixel_corners( X, Y )
%pcolor_pixel_corners Converts from pixel centers to lower left corner points.
%   The pcolor function assumes that the points passed to it define the
%   corner of a pixel, rather than its center.  With most data
%   applications, this distinction may not matter, but when the actual
%   position of the pixel edges relative to other data (e.g. geography)
%   matters, this function will make a best guess as to the position of the
%   corners that will put the points passed in as X and Y in the center of
%   the pixels.
%
%   Josh Laughner <joshlaugh5@gmail.com> 6 Aug 2014

X_corner = zeros(size(X));
Y_corner = zeros(size(Y));

% For all except the first row and column, the corner for pixel (i,j) will
% be the average of points (i,j), (i-1,j), (i,j-1), and (i-1,j-1).

tmp(:,:,1) = X(2:end,2:end);
tmp(:,:,2) = X(1:end-1,2:end);
tmp(:,:,3) = X(2:end,1:end-1);
tmp(:,:,4) = X(1:end-1,1:end-1);

X_corner(2:end,2:end) = mean(tmp,3);
clear ('tmp');

tmp(:,:,1) = Y(2:end,2:end);
tmp(:,:,2) = Y(1:end-1,2:end);
tmp(:,:,3) = Y(2:end,1:end-1);
tmp(:,:,4) = Y(1:end-1,1:end-1);

Y_corner(2:end,2:end) = mean(tmp,3);

% For the first row and column, take the difference between the 2nd and 3rd
% row or column as the best guess for the difference between the 1st and
% 2nd.

X_corner(1,:) = X_corner(2,:) + (X_corner(2,:) - X_corner(3,:));
X_corner(:,1) = X_corner(:,2) + (X_corner(:,2) - X_corner(:,3));
Y_corner(1,:) = Y_corner(2,:) + (Y_corner(2,:) - Y_corner(3,:));
Y_corner(:,1) = Y_corner(:,2) + (Y_corner(:,2) - Y_corner(:,3));

end

