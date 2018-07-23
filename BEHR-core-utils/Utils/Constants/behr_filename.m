function [ fname ] = behr_filename( date_in, prof_mode, region, ext, behr_version )
%BEHR_FILENAME Create the file name for a BEHR file of a given date
%   FNAME = BEHR_FILENAME( DATE_IN ) Given the date, DATE_IN, as a date
%   number or a string which datestr can parse without a specified format,
%   will construct the proper BEHR .mat filename for the current BEHR
%   version with wildcards ('*') for the region and profile type.
%
%   FNAME = BEHR_FILENAME( DATE_IN, PROF_MODE ) will insert the string
%   PROF_MODE in the file name in the proper place (in upper case) for
%   files with that profile type (usually 'daily' or 'monthly'). The region
%   remains a wildcard.
%
%   FNAME = BEHR_FILENAME( DATE_IN, PROF_MODE, REGION ) will insert the
%   string REGION in the proper place (in upper case) for files retrieved
%   over that region. 
%
%   FNAME = BEHR_FILENAME( DATE_IN, PROF_MODE, REGION, EXT ) uses the
%   extension EXT instead of .mat. The leading . may be included or not. 
%
%   FNAME = BEHR_FILENAME( DATE_IN, PROF_MODE, REGION, EXT, BEHR_VERSION )
%   will put the string BEHR_VERSION in place of the current version string
%   in the filename.
%
%   FNAME = BEHR_FILENAME( DATE_IN, PROF_MODE, REGION, EXT, true ) is an
%   old syntax to put a wildcard in for the version string. It is included
%   for backwards compatibility.
%
%   The idea of putting wildcards in the file name is that it can be passed
%   to DIR() in order to get a list of files that match for the given
%   DATE_IN, which can help you write functions that work for any region,
%   profile mode, or version.
%
%   An alternate syntax exists to generate older filenames that do not
%   include the region and profile mode:
%
%   FNAME = BEHR_FILENAME( DATE_IN, 'prefix', PREFIX, EXT, BEHR_VERSION )
%   will use PREFIX as the string to insert before the version string, i.e.
%   FNAME will behr PREFIX_VERSION_DATESTR.EXT

E = JLLErrors;

if strcmpi(prof_mode, 'prefix')
    prefix = region;
else
    if ~exist('prof_mode', 'var')
        prof_mode = '*';
    end
    
    if ~exist('region', 'var')
        region = '*';
    end
    
    prefix = sprintf('OMI_BEHR-%s_%s', upper(prof_mode), upper(region));
end
    
if ~exist('ext','var')
    ext = 'mat';
else
    % Remove a leading "."; it is included already
    ext = regexprep(ext, '^\.', '');
end

if ~exist('behr_version', 'var')
    behr_version = false;
elseif ~ischar(behr_version) && ((~islogical(behr_version) && ~isnumeric(behr_version)) || ~isscalar(behr_version))
    E.badinput('BEHR_VERSION must be a string or a scalar number or boolean (if given)')
end

if ischar(behr_version)
    ver_str = behr_version;
elseif behr_version
    ver_str = '*';
else
    ver_str = BEHR_version();
end
    


if ischar(date_in) && strcmp(date_in, '*')
    date_string = date_in;
else
    date_string = datestr(date_in, 'yyyymmdd');
end
fname = sprintf('%s_%s_%s.%s', prefix, ver_str, date_string, ext);

end

