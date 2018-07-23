function [ sza ] = solar_zenith_from_time(date_in, lon, lat, varargin )
%SOLAR_ZENITH_FROM_TIME Calculate a solar zenith angle for given time and location
%   SZA = SOLAR_ZENITH_FROM_TIME( DATE_IN, LON, LAT ) computes the solar
%   zenith angle for the time and position on the Earth's surface specified
%   by DATE_IN, LON, and LAT. DATE_IN may be any Matlab recognized date
%   identifier (such as a date number or string). LON and LAT are to be
%   given in degrees, for LON, west is negative.
%
%   SZA = SOLAR_ZENITH_FROM_TIME( ___, 'noon' ) will replace the time part
%   of DATE_IN with the estimated time for solar noon, i.e. when the sun is
%   highest in the sky. The date component of DATE_IN is still used.

date_in = as_utc(date_in, lon);
gamma = fractional_year(date_in);

if ~ismember('noon', varargin)
    eqtime = equation_of_time(gamma);
    toff = time_offset(eqtime, lon);
    tst = true_solar_time(date_in, toff);
else
    % Solar noon should be the exact middle of the day
    tst = (24*60)/2;
end

ha = solar_hour_angle(tst);
decl = solar_declination(gamma);

sza = solar_zenith_angle(lat, decl, ha);
end


function date_in = as_utc(date_in, longitude)
date_in = datenum(date_in) - longitude ./ 15;
end


%https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
function gamma = fractional_year(date_in)
% Convert a date into radians around the sun
if leapyear(year(date_in))
    n_days = 366;
else
    n_days = 365;
end

gamma = 2*pi/n_days * (modis_date_to_day(date_in) - 1 + (hour(date_in) - 12)/24);
end

function eqtime = equation_of_time(gamma)
% Calculate the difference between the mean and apparent solar time in
% minutes
eqtime = 229.18 .* (0.000075 + 0.001868 .* cos(gamma) - 0.032077 .* sin(gamma) - 0.014615 .* cos(2 .* gamma) - 0.040849 .* sin(2 .* gamma) );
end

function toff = time_offset(eqtime, longitude)
% Since I convert the time to UTC, no timezone offset should be needed.
toff = eqtime + 4 .* longitude;
end

function tst = true_solar_time(date_in, toff)
tst = hour(date_in)*60 + minute(date_in) + second(date_in)/60 + toff;
end

function ha = solar_hour_angle(tst)
% Returns the solar angle in degrees give the true solar time
ha = (tst ./ 4) - 180;
end

function decl = solar_declination(gamma)
decl = 0.006918 - 0.399912 .* cos(gamma) + 0.070257 .* sin(gamma) - 0.006758 .* cos(2 .* gamma)...
    + 0.000907 .* sin(2 .* gamma)- 0.002697 .* cos(3 .* gamma) + 0.00148 .* sin(3 .* gamma);
end

function sza = solar_zenith_angle(lat, decl, ha)

% Convert solar hour angle and latitude from degrees to radians
lat = lat .* pi ./ 180;
ha = ha .* pi ./ 180;

cos_sza = sin(lat) .* sin(decl) + cos(lat) .* cos(decl) .* cos(ha);
sza = acosd(cos_sza);
end