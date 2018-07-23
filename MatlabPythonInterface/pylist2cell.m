function [ cell_out ] = pylist2cell( list )
%PYLIST2CELL Converts a Python list into a cell array.
%   CELL_OUT = PYLIST2CELL( LIST ) Given a py.list object LIST, converts it
%   to a Matlab cell array CELL_OUT by calling PYTHON2MATLAB on each
%   element of the list.

if ~isa(list, 'py.list') && ~isa(list, 'py.tuple')
    error('pyinterface:bad_input', 'LIST must be of type "py.list" or "py.tuple"')
end

cell_out = cell(list);
for a=1:numel(cell_out)
    cell_out{a} = python2matlab(cell_out{a});
end

end

