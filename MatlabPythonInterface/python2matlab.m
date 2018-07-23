function [ matval ] = python2matlab( val )
%PYTHON2MATLAB Convert a Python object to corresponding Matlab type
%   MATVAL = PYTHON2MATLAB( VAL ) Convert a python object VAL into a Matlab
%   type MATVAL. Calls the appropriate conversion function based on the
%   type:
%       numeric -> numeric, direct conversion, result will be typed and may
%           not be a double
%       py.str -> string, via char()
%       py.list -> cell array, via pylist2cell
%       py.numpy.ndarray -> matlab array
%       py.dict -> structure (pydict2struct)

if isnumeric(val) || islogical(val)
    matval = val;
elseif isa(val, 'py.int')
    matval = int64(val);
elseif isa(val, 'py.str')
    matval = char(val);
elseif isa(val,'py.list') || isa(val, 'py.tuple')
    matval = pylist2cell(val);
elseif isa(val, 'py.numpy.ndarray')
    matval = numpyarray2matarray(val);
elseif isa(val, 'py.dict')
    matval = pydict2struct(val);
else
    error('pyinterface:not_implemented','Cannot convert value of type %s', class(val));
end

end

