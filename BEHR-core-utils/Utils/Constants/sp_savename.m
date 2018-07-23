function [ fname ] = sp_savename( date_in, region, ext, behr_version )
%SP_SAVENAME Create the file name for an SP file of a given date
%   FNAME = SP_SAVENAME( DATE_IN ) Given the date, DATE_IN, as a date
%   number or a string which datestr can parse without a specified format,
%   will construct the proper SP .mat filename for the file which contains
%   the SP, pixel corner, modis, and globe data. This will include a
%   wildcard for the region.
%
%   FNAME = SP_SAVENAME( DATE_IN, REGION ) as before, but includes REGION
%   in the proper place in the name (in upper case).
%
%   FNAME = SP_SAVENAME( DATE_IN, REGION, EXT ) uses the extension EXT
%   instead of .mat. EXT may include the leading . or not. You can pass a
%   '*' for REGION to allow any region. An empty string for region will
%   omit the region part of the file name, as in the version 2 naming
%   convention.
%
%   FNAME = SP_SAVENAME( DATE_IN, REGION, EXT, BEHR_VERSION ) will put the
%   string BEHR_VERSION in the filename as the version string. You can pass
%   a '*' for BEHR_VERSION to allow any version.
%
%   FNAME = SP_SAVENAME( DATE_IN, REGION, EXT, true ) will put a wildcard
%   in place of the version string. Included for backwards compatibility.
%
%   The idea of putting wildcards in the file name is that it can be passed
%   to DIR() in order to get a list of files that match for the given
%   DATE_IN, which can help you write functions that work for any region,
%   profile mode, or version.

E = JLLErrors;

validate_date(date_in);

if ~exist('region', 'var')
    region = '*';
end

if ~exist('ext','var')
    ext = 'mat';
else
    if ~ischar(ext)
        E.badinput('EXT (if given) must be a string');
    end
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

if isempty(region)
    fname = sprintf('OMI_SP_%s_%s.%s', ver_str, datestr(date_in, 'yyyymmdd'), ext);
else
    fname = sprintf('OMI_SP_%s_%s_%s.%s', upper(region), ver_str, datestr(date_in, 'yyyymmdd'), ext);
end

end

