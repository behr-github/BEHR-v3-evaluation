function [ utcsec ] = local2utc( localstr, timezone)
%local2utc: Converts a local time (as string, using military time) into
%seconds after midnight utc.
%
%   Timezones:
%       EST = Eastern Std.      EDT = Eastern Daylight
%       CST = Central Std.      CDT = Central Daylight
%       MST = Mountain Std.     MDT = Mountain Daylight
%       PST = Pacific Std.      PDT = Pacific Daylight

zones = {'est','edt','cst','cdt','mst','mdt','pst','pdt'};
offset = [-5, -4, -6, -5, -7, -6, -8, -7];

% if the user passes a US timezone abbreviation, figure out the offset
if ischar(timezone)
    t = strcmpi(timezone,zones);
    hours_added = offset(t);
% otherwise, assume the user passed a numerical definition of timezone
% (e.g. +/- hours from UTC).
else
    hours_added = timezone;
end
colon_pos = strfind(localstr,':');
h = str2double(localstr(1:colon_pos-1)) - hours_added;
m = str2double(localstr(colon_pos+1:end));

utcsec = h*3600+m*60;

end

