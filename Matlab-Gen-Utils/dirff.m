function [ files ] = dirff( path_in )
%DIRFF Return directory structure with full paths to each file
%   FILES = DIRFF( PATH_IN ) Returns a structure of files matching the
%   pattern PATH_IN, where FILES is the value of DIR( PATH_IN ) with each
%   entries name prepended with the directory part of PATH_IN.

files = dir(path_in);
dirname = fileparts(path_in);

for a=1:numel(files)
    files(a).name = fullfile(dirname, files(a).name);
end

end

