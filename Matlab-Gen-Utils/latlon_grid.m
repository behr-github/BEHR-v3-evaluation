function [ LON, LAT ] = latlon_grid( res, lonlim, latlim )
%LATLON_GRID Produces lat/lon grids with the specified resolution.
%   Only required input is the resolution in degrees. lonlim and latlim
%   will assume that it is global unless they are given. LON and LAT are
%   returned as matrices with longitude along the first dimension and
%   latitude along the second.

if ~exist('lonlim','var')
    lonlim = [-180 180];
end
if ~exist('latlim','var')
    latlim = [-90 90];
end

lonvec = min(lonlim)+res/2:res:max(lonlim)-res/2;
latvec = min(latlim)+res/2:res:max(latlim)-res/2;

[LON, LAT] = meshgrid(lonvec, latvec);
LON = LON';
LAT = LAT';


end

