function [lon_corners, lat_corners] = geos_latlon_bounds_to_corners(lon_bounds, lat_bounds, varargin)
% Convert GEOS-Chem longitude/latitude bounds to more standard formats
%
%   In the GEOS-Chem ND51 output diagnostic (at least for v9-02), latitude
%   and longitude bounds for grid cells are output in 2-by-nlat and
%   2-by-nlon arrays, respectively, where latitude_bounds(:,1) are the
%   south and north limits for the first row of grid cells and
%   longitude_bounds(:,1) are the west and east bounds for the first column
%   of grid cells. ("row" and "column" refer to the geospatial orientation;
%   whether this corresponds to actual rows and columns in the data array
%   probably depends on whether you read them in as column-major or
%   row-major). This format for lat/lon corners is different from any
%   format that the rest of my code is used to, so this function converts
%   these bounds into a more standard format.
%
%   [ LON_CORNERS, LAT_CORNERS ] = GEOS_LATLON_BOUNDS_TO_CORNERS(
%   LON_BOUNDS, LAT_BOUNDS ) Takes the 2-by-n arrays LON_BOUNDS and
%   LAT_BOUNDS from the ND51 file and returns the corner arrays LON_CORNERS
%   and LAT_CORNERS. These will be (nlon+1)-by-(nlat+1) arrays, so that the
%   corners for grid cell (i,j) are LON_CORNERS(i:i+1, j:j+1) and
%   LAT_CORNERS(i:i+1, j:j+1).
%
%   Additional parameters:
%
%       'format' - string that controls the output format. Default is
%       'grid' which gives the (nlon+1)-by-(nlat+1) arrays. The other
%       option is 'pixel' which returns 4-by-nlon-by-nlat 3D arrays, so
%       that the corners for pixel (i,j) are LON_CORNERS(:, i, j) and
%       LAT_CORNERS(:, i, j). This mimics the structure of pixel corners in
%       the OMPIXCOR product.

E = JLLErrors;

p = inputParser;
p.addParameter('format','grid');
p.parse(varargin{:});
pout = p.Results;

output_format = pout.format;
allowed_formats = {'grid','pixel'};
if ~ismember(output_format, allowed_formats)
    E.badinput('The parameter ''format'' must be one of: %s', strjoin(allowed_formats, ', '));
end

if size(lon_bounds,1) ~= 2 || size(lat_bounds,1) ~= 2
    E.badinput('LON_BOUNDS and LAT_BOUNDS are expected to have length 2 in the first dimension');
end

if strcmpi(output_format, 'grid')
    [lon_corners, lat_corners] = make_grid_format(lon_bounds, lat_bounds);
elseif strcmpi(output_format, 'pixel')
    [lon_corners, lat_corners] = make_pixel_format(lon_bounds, lat_bounds);
else
    E.notimplemented('No method for format = "%s"', output_format);
end

end

function [loncorn, latcorn] = make_grid_format(lon_bounds, lat_bounds)
% What I call "grid" format is where the corners are given by 2D arrays
% that for N grid cells in the lon direction and M in the lat direction,
% are (N+1) by (M+1). So the corners for grid cell (i,j) will be given by
% loncorn(i:i+1, j:j+1) and latcorn(i:i+1, j:j+1).

E = JLLErrors;

% First verify that the adjacent grid cells do indeed have the same value
% for the common corner.
if ~check_corners_abut(lon_bounds) || ~check_corners_abut(lat_bounds)
    E.badinput('In "grid" format, the grid cells must abut, so bounds(2,i) must equal bounds(1,i+1)');
end

% Given that the bounds abut, we can collapse them to a single vector and
% then replicate it. For sanity, unlike MESHGRID we will put longitude on
% the first dimension and latitude on the second. We checked that
% bounds(2,i) == bounds(1,i+1), so we take the first row then append the
% last element of the second row onto it.
lon_bounds = [lon_bounds(1,:), lon_bounds(2,end)];
lat_bounds = [lat_bounds(1,:), lat_bounds(2,end)];

% need lon_bounds as column vector, lat_bounds as row vector
loncorn = repmat(lon_bounds(:), 1, numel(lat_bounds));
latcorn = repmat(lat_bounds, numel(lon_bounds), 1);

end

function [loncorn, latcorn] = make_pixel_format(lon_bounds, lat_bounds)
% In this case, it doesn't matter if the pixels abut. We just need to make
% sure that the corners go counterclockwise. I'll define the NW corner as
% the first.
loncorn = nan(4, size(lon_bounds,2), size(lat_bounds,2));
latcorn = nan(4, size(lon_bounds,2), size(lat_bounds,2));
for i_lon = 1:size(lon_bounds,2)
    for i_lat = 1:size(lat_bounds,2)
        % Remember, lon_bounds(1,:) = west and lon_bounds(2,:) = east while
        % lat_bounds(1,:) = south and lat_bounds(2,:) = north
        loncorn(:,i_lon,i_lat) = lon_bounds([1 1 2 2], i_lon);
        latcorn(:,i_lon,i_lat) = lat_bounds([2 1 1 2], i_lat);
    end
end
end

%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%

function chk = check_corners_abut(gc_bounds)
% The bounds have the form
%   [ left/bottom ... ]
%   [ right/top ... ]
% for each pixel in sequence, so we want gc_bounds(2,i) == gc_bounds(1,i+1)
% The end bounds don't wrap around in an easily checkable way.

x1 = gc_bounds(2,1:end-1);
x2 = gc_bounds(1,2:end);

% I chose a somewhat loose tolerance because these are being read in as
% single precision values.
chk = all(abs(x1 - x2) < 1e-5);

end