function [ mcd43_lon, mcd43_lat ] = mcd43_lon_lats( )
%mcd43_lon_lats Returns the center longitude and latitude for an MCD43C3 dataset


%MODIS albedo is given in 0.05 degree cells and a single file covers the
%full globe, so figure out the lat/lon of the middle of the grid cells as:
mcd43_lat=-90+0.05/2:0.05:90-0.05/2; mcd43_lat=mcd43_lat'; mcd43_lat=repmat(mcd43_lat,1,7200);
mcd43_lon=-180+0.05/2:0.05:180-0.05/2; mcd43_lon=repmat(mcd43_lon,3600,1);

end

