function [ data, varargout ] = remove_merge_fills( Merge, field, varargin )
% REMOVE_MERGE_FILLS - Return data from a Merge structure with fill value exchanged for NaNs
%   DATA = REMOVE_MERGE_FILLS( MERGE, FIELD ) Returns the values of
%   MERGE.Data.(FIELD) with fill values replaced with NaNs, where MERGE is
%   a Merge structure created by reading an ICART file with
%   read_icartt_file and FIELD is a string that gives a field of
%   MERGE.Data.
%
%   [ DATA, UTC, ALT, LON, LAT ] = REMOVE_MERGE_FILLS( MERGE, FIELD ) Still
%   returns DATA with fill values removed, but also does the same thing for
%   the UTC, ALTP (pressure altitude), LONGITUDE, and LATITUDE fields. (The
%   UTC field is assumed to have no fill values). This syntax was created
%   because it assumes that you'd often like to know when and where a
%   measurement occurred. If the MERGE structure uses different names for
%   these fields, you can override the names with the 'time', 'alt', 'lon',
%   and 'lat' parameters respectively.
%
%   [ ... ] = REMOVE_MERGE_FILLS( MERGE, FIELD, 'unit', UNIT ) Used with
%   either prior syntax, this will try to convert the data field to the
%   unit given by the string UNIT by using convert_units.m. This requires
%   that the unit given in the merge structure for that field be understood
%   by convert_units and of the same type as the unit you wish to convert
%   to.

p = advInputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('field',@isstr);
p.addParameter('time','UTC',@isstr);
p.addParameter('alt','ALTP',@isstr);
p.addParameter('lat','LATITUDE',@isstr);
p.addParameter('lon','LONGITUDE',@isstr);
p.addParameter('unit','');
p.addFlag('utc2datenum');
p.addParameter('DEBUG_LEVEL',1,@(x) (isscalar(x) && isnumeric(x)));

p.parse(Merge,field,varargin{:});
pout = p.AdvResults;
Merge = pout.Merge;
field = pout.field;
utc_field = pout.time;
alt_field = pout.alt;
lon_field = pout.lon;
lat_field = pout.lat;
unit_out = pout.unit;
do_utc_to_datenum = pout.utc2datenum;

E = JLLErrors;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

data = Merge.Data.(field).Values;
fill_val = Merge.Data.(field).Fill;
try
    ulod = Merge.metadata.upper_lod_flag;
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField')
        ulod = fill_val;
        if DEBUG_LEVEL > 0; fprintf('remove_merge_fills: ulod field does not exist\n'); end
    end
end

try
    llod = Merge.metadata.lower_lod_flag;
catch err
    if strcmp(err.identifier,'MATLAB:nonExistentField')
        llod = fill_val;
        if DEBUG_LEVEL > 0; fprintf('remove_merge_fills: llod field does not exist\n'); end
    end
end

if isnumeric(fill_val) && isscalar(fill_val)
    fills = data == fill_val | data == ulod | data == llod;
    data(fills) = NaN;
elseif strcmpi(fill_val,'n/a')
    % do nothing - this means that there should be absolutely no data
    % missing because no fill value was assigned in the Merge file. Usually
    % only happens for UTC.
else
    E.badvar('fill_val','Is not a numeric scalar value or "N/A"');
end


% Let the user specify longitude as the main input, in which case we need
% to correct the sign
if ~isempty(regexpi(field, 'longitude', 'ONCE'))
    data = lon_fix(data);
end

% Let the user specify UTC as the main input, in which case check if we
% should convert it to a date number
if regcmpi(field, 'utc') && do_utc_to_datenum
    data = convert_utc_to_datenum(data, Merge.metadata.date);
end

% If a unit conversion is requested, try to do so
if ~isempty(unit_out)
    data = convert_units(data, Merge.Data.(field).Unit, unit_out);
end

if nargout > 5; varargout{5} = fills; end

if nargout > 1
    utc = Merge.Data.(utc_field).Values; 
    if do_utc_to_datenum
        utc = convert_utc_to_datenum(utc, Merge.metadata.date);
    end
    %utcfills = eval(sprintf('Merge.Data.%s.Fill',utcfield));
    %utc(utc==utcfills) = NaN; % There shouldn't ever be fill values in the UTC field
    varargout{1} = utc;
end
if nargout > 2
    alt = Merge.Data.(alt_field).Values; 
    altfills = Merge.Data.(alt_field).Fill;
    alt(alt==altfills) = NaN;
    varargout{2} = alt;
end
if nargout > 4
    lat = Merge.Data.(lat_field).Values;
    latfills = Merge.Data.(lat_field).Fill;
    lat(lat==latfills) = NaN;
    varargout{4} = lat;
end
if nargout > 3
    lon = Merge.Data.(lon_field).Values;
    lonfills = Merge.Data.(lon_field).Fill;
    lon(lon==lonfills) = NaN;
    lon = lon_fix(lon);
    varargout{3} = lon;
end



end

function lon = lon_fix(lon)
lon(lon>180) = lon(lon>180) - 360;
end

function utc = convert_utc_to_datenum(utc, merge_date)
utc = utc / 86400 + datenum(merge_date);
end