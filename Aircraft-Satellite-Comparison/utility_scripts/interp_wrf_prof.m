function [ wrf_profile, wrf_pressures ] = interp_wrf_prof( PROFILE, lat, lon)
%interp_wrf_prof Given a WRF profile structure and a lat/lon, will
%interpolate the profile to those coordinate.

wrf_lonvec = PROFILE.Longitude(1,:);
wrf_latvec = PROFILE.Latitude(:,1);
wrf_pressures = fliplr(PROFILE.Pressure);
wrf_no2 = flipdim(PROFILE.NO2_profile,1);
[Y,Z,X] = meshgrid(wrf_latvec,wrf_pressures,wrf_lonvec);

% Interpolate the WRF profiles to the location of the aircraft profile. The
% order of the arguments may look funny, but it seems to work to match the
% ordering of the NO2 profile matrix.
wrf_profile = interp3(Y,Z,X,wrf_no2,lat,wrf_pressures,lon);
wrf_profile = flipud(wrf_profile);

% Flip the pressure around again (and make them a column) to match the
% flipped profile for output
wrf_pressures = flipud(wrf_pressures');

end

