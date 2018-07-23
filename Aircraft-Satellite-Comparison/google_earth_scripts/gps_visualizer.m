function [  ] = gps_visualizer( lon, lat, alt, varargin )
%gps_visualizer Produces a text file to upload for Google Earth visualization
%   www.gpsvisualizer.com accepts text input files to produce flight tracks
%   in a .kmz or .kml file.  It needs a .csv file with latitude, longitude,
%   and altitude.  Pass these into the function, which will open a dialogue
%   to select where to save the file.  You can also pass name/value pairs
%   of additional datasets to use as custom fields on the visualizer (not
%   fully implemented).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = JLLErrors;
n = numel(varargin)/2;
if n ~= round(n)
    E.callError('bad_input','Inputs after altitude must be pairs of arguments: field name as string, field data as vector');
end

if numel(lon) ~= numel(lat) || numel(lon) ~= numel(alt)
    E.numelMismatch('lon','lat','alt');
end

% Check that the extra data has the same numel of elements as lon and save
% it to a matrix as column vector


[fname, dirname] = uiputfile('.txt','Choose where to save the output text file');

file = fullfile(dirname,fname);
fid = fopen(file,'w');


titles = cell(1,n);
varspec = [repmat('%s,',1,n-1), '%s'];
for a=0:n-1
    titles{a+1} = varargin{a*2+1};
end

%titlespec = sprintf('longitude,latitude,altitude,%s\n',varspec);
%fprintf(fid,titlespec,titles{:});
fprintf(fid,'longitude,latitude,altitude\n');

for a=1:numel(lon)
    fprintf(fid,'%f,%f,%f\n',lon(a),lat(a),alt(a)*1000); %Generally, the altitude is given in km, but the visualizer expects meters
end

fclose(fid);

end

