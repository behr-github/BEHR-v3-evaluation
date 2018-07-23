function [ S, files ] = load_by_glob( filepattern, multifile )
%LOAD_BY_GLOB Load single file that matches wildcard pattern
%   S = LOAD_BY_GLOB( FILEPATTERN ) Takes a pattern FILEPATTERN which can
%   use the * wildcard to match one file. If multiple files are matched, an
%   error will occur. S will be a structure that is the result of executing
%   S = load(file).
%
%   S = LOAD_BY_GLOB( FILEPATTERN, TRUE ) will allow multiple files to be
%   loaded as long as they have the same variables, otherwise an error will
%   occur. S will now be a multi-element structure, with each element
%   corresponding to one file.
%
%   [S, FILENAMES] = LOAD_BY_GLOB( ___ ) returns the resolved filenames as
%   a cell array with either previous syntax.

E = JLLErrors;

if ~ischar(filepattern)
    E.badinput('FILEPATTERN must be a string')
end

if ~exist('multifile', 'var')
    multifile = false;
elseif (~isnumeric(multifile) && ~islogical(multifile)) || ~isscalar(multifile)
    E.badinput('MULTIFILE (if given) must be a scalar number or logical');
end

[file_dir, ~, file_ext] = fileparts(filepattern);

if isempty(file_dir)
    file_dir = '.';
end

if ~exist(file_dir, 'dir')
    E.dir_dne(file_dir);
elseif ~strcmp(file_ext, '.mat')
    warning('FILEPATTERN does not specify the .mat extension, load may fail or behave strangely if pattern matches non .mat file')
end

F = dir(filepattern);
if numel(F) < 1
    E.filenotfound(filepattern);
elseif ~multifile && numel(F) > 1
    E.toomanyfiles(filepattern);
else
    files = cell(1, numel(F));
    for a=1:numel(F)
        files{a} = fullfile(file_dir, F(a).name);
        try
            S(a) = load(fullfile(file_dir, F(a).name));
        catch err
            if strcmp(err.identifier, 'MATLAB:heterogeneousStrucAssignment')
                E.callError('File %s has different variables than previous files; cannot concatenate', F(a).name);
            else
                rethrow(err)
            end
        end
    end
end

end

