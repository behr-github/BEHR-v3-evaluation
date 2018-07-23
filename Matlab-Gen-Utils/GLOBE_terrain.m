function [ GLOBE_LON_GRID, GLOBE_LAT_GRID, GLOBE_ALT_METERS ] = GLOBE_terrain( lonbdy, latbdy, varargin )
%GLOBE_terrain Returns a grid of globe terrain heights in meters.
%   Pass the longitude and latitude boundaries of the area you want the
%   globe terrain height for, optionally, pass the factor to reduce the
%   resolution by.

if numel(varargin) < 1;
    ptskip = 1;
else 
    ptskip = varargin{1};
end

globe_dir = '/Volumes/share-sat/SAT/BEHR/GLOBE_Database/';
[GLOBE_ALT_METERS, refvec] = globedem(globe_dir,ptskip,latbdy,lonbdy);

cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(GLOBE_ALT_METERS,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
GLOBE_LAT_GRID = repmat(globe_lat_matrix,1,size(GLOBE_ALT_METERS,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(GLOBE_ALT_METERS,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
GLOBE_LON_GRID = repmat(globe_lon_matrix,size(GLOBE_ALT_METERS,1),1); 

end

