function [  ] = python_make_pickle( obj, filename )
%PYTHON_MAKE_PICKLE Save a Python object as a pickle file
%   PYTHON_MAKE_PICKLE( OBJ, FILENAME ) Saves the Python object, OBJ as a
%   pickle file, given by the path FILENAME. 

if isempty(regexp(class(obj),'py', 'once'))
    error('pyinterface:badinput','OBJ must be a Python type')
end
file_path = fileparts(filename);
if ~isempty(file_path) && ~exist(file_path, 'dir')
    error('pyinterface:dir_dne', 'The directory "%s" does not exist', file_path)
end

my_path = fileparts(mfilename('fullpath'));
py_util_dir = fullfile(my_path,'Utils');
if count(py.sys.path, py_util_dir) == 0
    insert(py.sys.path, int32(0), py_util_dir);
end

py.pyside.make_pickle(obj, filename);

end

