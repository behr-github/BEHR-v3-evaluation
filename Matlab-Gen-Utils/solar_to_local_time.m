function [ std_time ] = solar_to_local_time( solar_time, dates, lon )
%SOLAR_TO_LOCAL_TIME Converts solar time to local (standard) time
%   Solar time is defined such that noon is when the sun is highest in the
%   sky.  This does not always correspond to local standard time.  This
%   function will convert solar time to standard time.  Requires the solar
%   time (in decimal hours) and the longitude. dates and longitude can be
%   vectors or scalar, but if they are vector they must be the same size as
%   solar_time.
%
%   Equations from http://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-time
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 Nov 2015

E = JLLErrors;

%%%%% INPUT %%%%%
if iscell(dates)
    E.badinput('Input multiple dates as a vector of datenums or a character array. Cell arrays tend to behave poorly for conversion to datenums.')
elseif ischar(dates)
    dates = datenum(dates);
end
if isscalar(dates)
    lon = repmat(lon, size(solar_time));
elseif ndims(dates) ~= ndims(solar_time) || any(size(dates) ~= size(solar_time))
    E.badinput('Size of dates and solar_time must be the same, if lon is not a scalar')
end

if isscalar(lon)
    lon = repmat(lon, size(solar_time));
elseif ndims(lon) ~= ndims(solar_time) || any(size(lon) ~= size(solar_time))
    E.badinput('Size of lon and solar_time must be the same, if lon is not a scalar')
end


%%%%% CALCULATION %%%%%

% The local standard time meridian is the meridian closest to the point in
% question that is a multiple of 15.
lstm = round(lon/15)*15;

std_time = nan(size(solar_time));

for a=1:numel(solar_time)
    % The equation of time corrects for eccentricity in the Earth's orbit.
    d = modis_date_to_day(dates(a));
    B = 360/365 * (d - 81);
    eot = 9.87 * sind(2*B) - 7.53 * cos(B) - 1.5 * sin(B);
    
    % Time correction accounts for longitudinal variation
    tc = 4*(lon(a) - lstm(a)) + eot;
    
    % Convert from solar time
    std_time(a) = solar_time(a) - tc/60;
end

end

