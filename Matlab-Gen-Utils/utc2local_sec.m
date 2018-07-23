function [ local ] = utc2local_sec( utcsec, timezone )
%utc2local: Converts UTC time in seconds after midnight to local time.
%   Pass this function a vector of times in seconds after midnight UTC and
%   a time zone or numerical UTC offset, it will return a vector of times
%   in seconds after local midnight.  The time zone must either be a single
%   string, a cell of strings of the same size as the vector of utc times,
%   or a vector of numeric offsets of the same size as the vector of times.
%
%   Timezones:
%       EST = Eastern Std.      EDT = Eastern Daylight
%       CST = Central Std.      CDT = Central Daylight
%       MST = Mountain Std.     MDT = Mountain Daylight
%       PST = Pacific Std.      PDT = Pacific Daylight

E = JLLErrors;
zones = {'est','edt','cst','cdt','mst','mdt','pst','pdt'};
offset = [-5, -4, -6, -5, -7, -6, -8, -7];

% Check the inputs
if iscell(timezone) || (~isscalar(timezone) && ismatrix(timezone) && isnumeric(timezone))
    if numel(timezone) ~= numel(utcsec)
        error(E.badinput('Input vector/cell array timezone must have the same number of elements as the vector utcsec or be a scalar'));
    end
end

% If the user passes a string, it's assumed to be one of the US timezone
% abbreviations, otherwise it is assumed to be the numerical offset.
if ischar(timezone)
    t = strcmpi(timezone,zones);
    hours_added = offset(t);
elseif iscell(timezone)
    hours_added = zeros(size(timezone));
    for a=1:numel(timezone)
        t = strcmpi(timezone{a},zones);
        hours_added(a) = offset(t);
    end
else
    hours_added = timezone;
end

local = utcsec + (hours_added * 3600);

end

