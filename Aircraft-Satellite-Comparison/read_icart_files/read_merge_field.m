function [ varargout ] = read_merge_field( Merge, fieldname, varargin )
%READ_MERGE_FIELD Read a fill from a Merge structure, replacing fills with NaNs
%
%   [ VALUE ] = READ_MERGE_FIELD( MERGE, FIELDNAME ) will read the Value
%   field from MERGE.Data.(FIELDNAME), replace any fill, upper LOD, or
%   lower LOD values with NaNs, and return that as VALUE.
%
%   [ VALUE, UTC, ALT, LON, LAT ] = READ_MERGE_FIELD( MERGE, FIELDNAME )
%   will simultaneously return the UTC time, altitude, longitude, and
%   latitude, assuming that their field names are "UTC", "ALTP",
%   "LONGITUDE", and "LATITUDE".
%
%   [ ___ ] = READ_MERGE_FIELD( ___, 'unit', UNIT ) will use CONVERT_UNITS
%   to convert the data given in FIELDNAME from the units listed in the
%   Merge structure to UNIT.


error('deprecated:superior_function', 'This function is no longer supported; use REMOVE_MERGE_FILLS INSTEAD');


end

