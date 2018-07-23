function [ is_ocean, lw_lon, lw_lat ] = get_modis_ocean_mask( lonlim, latlim )
%GET_MODIS_OCEAN_MASK Load the 7-category land use mask
%   [ IS_OCEAN, LW_LON, LW_LAT ] = GET_MODIS_OCEAN_MASK( LONLIM, LATLIM )
%   Reads the MODIS land-water mask file at the path specified by
%   behr_paths.modis_land_mask and returns IS_OCEAN, a boolean array that
%   is true over ocean at the points defined by LW_LON and LW_LAT. LONLIM
%   and LATLIM should be 2-element vectors giving the min/max longitude and
%   latitude to read in.
%
%   The land-water mask was available from
%   ftp://rsftp.eeos.umb.edu/data02/Gapfilled/ as the file
%   Land_Water_Mask_7Classes_UMD as of 21 Sept 2017.

% The land mask seems to be given at 30 arc second resolution
[lw_lon, lw_lat, lw_xx, lw_yy] = modis_cmg_latlon(1/120, lonlim, latlim);
lw_lat = lw_lat';

hdfi = hdfinfo(behr_paths.modis_land_mask);
lw_mask = hdfreadmodis(hdfi.Filename, hdfdsetname(hdfi, 1, 1, 'LW_MASK_UMD'), 'log_index', {lw_yy, lw_xx});

% According to the attribute "LW_Label" on the LW_MASK_UMD dataset, of the
% 7 labels, 0 = shallow ocean, 6 = moderate or continental ocean, and 7 =
% deep ocean. 1-5 are land, coastlines, lakes, or ephemeral wate.
is_ocean = lw_mask >= 6 | lw_mask == 0;
end

