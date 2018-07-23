function [ F ] = rfind( pattern, varargin )
%RFIND Finds all files matching a pattern recursively under the given directory
%   F = RFIND( PATTERN ) Returns a structure F similar to the output of the
%   Matlab built-in DIR() for all files that match PATTERN in the directory
%   given in PATTERN and all its subdirectories. Specifically, given:
%
%       [top_dir, p1, p2] = fileparts(PATTERN)
%
%   any file matching [p1, p2] in top_dir or any subdirectory therein is
%   returned. The path given for the name of each element of F will be the
%   path to the file relative to top_dir. As with DIR() the path given as
%   part of PATTERN can be absolute or relative to the current directory.
%
%   F = RFIND( PATTERN, 'fullpath' ) prepends top_dir (see above) to each
%   file name returned. 



E = JLLErrors;
xx = strcmpi('fullpath', varargin);
if any(xx)
    return_full = true;
    varargin(xx) = [];
else
    return_full = false;
end

[top_dir, fname, fext] = fileparts(pattern);
fpattern = [fname, fext];
if isempty(top_dir)
    top_dir = '.';
end

% First find files in this directory that match the given pattern
F = dir(pattern);

% Now find directories and recurse into them
F_all = dir(top_dir);
for a=3:numel(F_all) % start from 3 to skip '.' and '..'
    if F_all(a).isdir
        F_tmp = rfind(fullfile(top_dir, F_all(a).name, fpattern));
        % The paths returned will be relative to fullfile(top_dir,
        % F_all(a).name), add F_all(a).name so that they are relative to
        % top_dir
        for b=1:numel(F_tmp)
            F_tmp(b).name = fullfile(F_all(a).name, F_tmp(b).name);
        end
        
        F = veccat(F, F_tmp);
    end
end

% If full paths requested, add them now
if return_full
    for a=1:numel(F)
        F(a).name = fullfile(top_dir, F(a).name);
    end
end
end

