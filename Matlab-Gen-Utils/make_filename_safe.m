function [ filename ] = make_filename_safe ( filename, subchar, varargin )
%MAKE_FILENAME_SAFE Returns a sanitized version of the given filename
%   FILENAME = MAKE_FILENAME_SAFE( FILENAME ) returns a sanitized version
%   of FILENAME. The following characters are removed:
%
%       ! @ # $ ^ & * ( ) [ ] { } < > ? ' " : ; + = | ~ `
%
%     Percent signs are treated specially in that they get replaced with
%     the word 'percent' (capitalized if not preceeded by whitespace). This
%     is because this function is geared towards santizing figure titles
%     for file names, and '% difference' for example seems like a likely
%     plot title.
%
%   FILENAME = MAKE_FILENAME_SAFE( FILENAME, SUBCHAR ) replaces unsafe
%   characters with SUBCHAR instead of removing them.
%
%   The special behaviour of percent signs can be disabled by setting the
%   value of the parameter 'percent' to false.

p = inputParser;
p.addOptional('subchar', '', @ischar);
p.addParameter('percent', true);

if ~exist('subchar','var')
    subchar = '';
end
p.parse(subchar, varargin{:});
pout = p.Results;

subchar = pout.subchar;
percent_special = pout.percent;

if ~ischar(filename)
    E.badinput('FILENAME must be a string')
end

if ~ischar(subchar)
    E.badinput('SUBCHAR must be a string')
end

if ~isscalar(percent_special)
    E.badinput('The parameter PERCENT must be a scalar value')
end


re = '[!@#$^&*()\[\]{}<>?''":;+=|~`]';
filename = regexprep(filename, re, subchar);
if percent_special
    filename = regexprep(filename,'(?<=\w\s+)\%','percent'); % preceded by a letter before a whitespace - lowercase percent
    filename = regexprep(filename, '%', 'Percent'); % otherwise capital
else
    filename = strrep(filename, '%', subchar);
end

end

