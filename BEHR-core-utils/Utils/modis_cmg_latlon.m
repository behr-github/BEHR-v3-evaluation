function [ lon, lat, xx, yy ] = modis_cmg_latlon( resolution, varargin )
%MODIS_CMG_LATLON Generate the CMG grid for MODIS data
%   [ LON, LAT ] = MODIS_CMG_LATLON( RESOLUTION ) will yield longitude and
%   latitude vectors LON and LAT that describe the CMG grid with RESOLUTION
%   spacing in degrees.
%
%   [ ___ ] = MODIS_CMG_LATLON( RESOLUTION, LONLIM, LATLIM ) will restrict
%   the lon and lat vectors to the limits specified in LONLIM and LATLIM,
%   which must be two element numeric vectors.
%
%   [ LON, LAT, XX, YY ] = MODIS_CMG_LATLON( RESOLUTION, LONLIM, LATLIM )
%   provides logical vectors XX and YY that can be used to subset a MODIS
%   CMG array at the given resolution to the specified longitude and
%   latitude limits. XX is the vector for the longitudinal direction, YY
%   for the latitudinal direction. For a typical MODIS CMG array M, which
%   is usually lat x lon, M(yy,xx) will cut it down to the longitude and
%   latitude limits.
%
%   [ ___ ] = MODIS_CMG_LATLON( ___, 'grid' ) with any of the previous
%   syntaxes, this will output LON and LAT as a meshgrid, instead of
%   vectors.

E = JLLErrors;

xx = strcmpi(varargin, 'grid');
as_grid = any(xx);
varargin = varargin(~xx);

if numel(varargin) < 1
    lonlim = [-180 180];
else
    lonlim = varargin{1};
end
if numel(varargin) < 2
    latlim = [-90 90];
else
    latlim = varargin{2};
end

if ~isnumeric(resolution) || ~isscalar(resolution)
    E.badinput('RESOLUTION must be a scalar number')
elseif ~isnumeric(lonlim) || numel(lonlim) ~= 2
    E.badinput('LONLIM must be a two element numeric vector')
elseif ~isnumeric(latlim) || numel(latlim) ~= 2
    E.badinput('LATLIM must be a two element numeric vector')
end

lon = (-180 + resolution/2):resolution:(180 - resolution/2);
lat = (90 - resolution/2):-resolution:(-90 + resolution/2);

xx = lon >= min(lonlim) & lon <= max(lonlim);
yy = lat >= min(latlim) & lat <= max(latlim);

lon = lon(xx);
lat = lat(yy);

if as_grid
    [lon, lat] = meshgrid(lon, lat);
end

end

