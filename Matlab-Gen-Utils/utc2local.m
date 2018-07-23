function [ local ] = utc2local( utcsec, timezone)
%utc2local: Converts UTC time in seconds after midnight to local time.
%   Pass this function a time in seconds after midnight UTC and a time
%   zone, it will return a string with the local time in military time.
%
%   Timezones:
%       EST = Eastern Std.      EDT = Eastern Daylight
%       CST = Central Std.      CDT = Central Daylight
%       MST = Mountain Std.     MDT = Mountain Daylight
%       PST = Pacific Std.      PDT = Pacific Daylight

zones = {'est','edt','cst','cdt','mst','mdt','pst','pdt'};
offset = [-5, -4, -6, -5, -7, -6, -8, -7];

% If the user passes a string, it's assumed to be one of the US timezone
% abbreviations, otherwise it is assumed to be the numerical offset.
if ischar(timezone)
    t = strcmpi(timezone,zones);
    hours_added = offset(t);
else
    hours_added = timezone;
end
h = floor(utcsec/3600) + hours_added;
m = floor(mod(utcsec,3600)/60);

local = sprintf('%2d:%02d',h,m);

end

