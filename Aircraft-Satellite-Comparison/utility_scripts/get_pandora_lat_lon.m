function [lon, lat] = get_pandora_lat_lon(Merge)
% The Pandora's latitude and longitude should be given in the
% Merge.metadata LOCATION field. We'll allow any capitalization of LOCATION
% and then search for Lat #, Long # in the string. 
E = JLLErrors;

fns = fieldnames(Merge.metadata);
xx = strcmpi('location', fns);
if sum(xx) ~= 1
    E.callError('location_metadata', 'Could not identify the location metadata')
end

loc_string = Merge.metadata.(fns{xx});

% An example location string from the Aldino site in Maryland:
% ground- Aldino (Maryland), USA, Lat 39.563333°, Long -76.203889°, alt= 128m
% So we want to search for the first occurence of numbers after Lat and
% Long respectively. To be safe, I will search for Lat and Lon and allow
% the to be any characters between them and the beginning of the number,
% which the first character will either be a - or a digit.
lat_string = regexpi(loc_string, '(?<=lat).+?\-?\d+\.?\d*', 'match', 'once');
lon_string = regexpi(loc_string, '(?<=lon).+?\-?\d+\.?\d*', 'match', 'once');
% This match will include the characters before the beginning of the number
% after "lat" or "lon" because you cannot have a variable length
% look-behind assertion
lat = str2double(regexp(lat_string, '\-?\d.*', 'match', 'once'));
lon = str2double(regexp(lon_string, '\-?\d.*', 'match', 'once'));

end