function [ files_cell, bad_format ] = files_input( files, varargin )
%FILES_INPUT Convert and validate input of file names
%   Often it's convenient to be able to accept a list of files in several
%   formats, converting that list internally to a consistent type in order
%   to make the rest of the code easier to maintain. This function accepts
%   an input of a list of files in any one of three formats:
%       1) A string (for a single file)
%       2) A cell array of strings
%       3) A structure output from DIR() (specifically, one with the field
%       "name")
%   and converts any of those three to a cell array of strings.
%
%   FILES_CELL = FILES_INPUT( FILES ) Given a list of files in any of the
%   three inputs mentioned above, converts them to a cell array
%   (FILES_CELL) or throws an error if the list is not in one of the
%   appropriate formats. 
%
%   [ FILES_CELL, BAD_FORMAT ] = FILES_INPUT( FILES, 'softfail', true )
%   will not throw an error if FILES is not in one of the expected formats,
%   instead the return value BAD_FORMAT will be true, an the user is
%   expected to handle whatever recovery is necessary.
%
%   Other parameters:
%
%       'argname' - changes the name of the input argument in the error
%       message. Default is "FILES", i.e. the message will is by default
%       'FILES must be a string, cell array of strings, or a structure with
%       the field "name"'. This value must be a string.
%
%       'checkexist' - default false, if true, then FILES_INPUT will check
%       whether each file listed exists as a file or a directory.
%       Currently, this just throws an error if one doesn't exist,
%       regardless of the value of 'softfail'.

E = JLLErrors;

p = inputParser;
p.addParameter('argname', 'FILES');
p.addParameter('softfail', false);
p.addParameter('checkexist', false);

p.parse(varargin{:});
pout = p.Results;

arg_name = pout.argname;
soft_fail = pout.softfail;
check_exist = pout.checkexist;

if ~ischar(arg_name)
    E.badinput('The value for "argname" must be a string');
elseif ~isscalar(soft_fail) || (~islogical(soft_fail) && ~isnumeric(soft_fail))
    E.badinput('The value for "softfail" must be a scalar number or logical');
elseif ~isscalar(check_exist) || (~islogical(check_exist) && ~isnumeric(check_exist))
    E.badinput('The value for "checkexist" must be a scalar number or logical');
end

bad_format = false;

if ischar(files)
    files_cell = {files};
elseif isstruct(files) && isfield(files, 'name')
    files_cell = {files.name};
elseif iscellstr(files)
    files_cell = files;
else
    if ~soft_fail
        E.badinput('%s must be a string, cell array of strings, or a structure with the field "name"', arg_name);
    end
    bad_format = true;
end

if check_exist
    for a=1:numel(files)
        if ~exist(files{a}, 'file') && ~exist(files{a}, 'dir')
            E.badinput('%s is neither a file nor a directory', files{a});
        end
    end
end

end

