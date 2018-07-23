function [ pyval ] = matlab2python( val, varargin )
%MATLAB2PYTHON Convert Matlab type into Python type.
%   PYVAL = MATLAB2PYTHON( VAL ) Converts the Matlab type VAL into a Python
%   type PYVAL. Calls the appropriate conversion function based on the type
%   of the input. Several of the other functions in this repo use this
%   recursively to convert nested cell arrays or structs.
%
%   Current type conversions are:
%       Scalar numbers -> scalar numbers (direct conversion)
%       Numeric arrays -> numpy ndarrays (matarray2numpyarray)
%       Strings -> strings (direct conversion)
%       Cell arrays -> lists (cell2pylist)
%       Structures -> dictionaries (struct2pydict)
%
%   PYVAL = MATLAB2PYTHON( VAL, 'array1' ) forces any scalar values to be
%   converted to 1D numpy arrays, rather than normal Python literals. (In
%   the future, 'array0' might be an option to convert to Numpy 0D arrays.
%   But it's not implemented yet, that's just why it's "array1" instead of
%   "array".)
%
%   PYVAL = MATLAB2PYTHON( VAL, VEC_AS_MAT ) controls how Matlab vectors
%   are converted to Numpy arrays. Options are 'never', 'row', 'column',
%   'always' (default is 'never'). See LIST_RECURSION for information on
%   each choice.

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

xx = ismember(varargin, {'array1','scalar'});
if sum(xx) == 0
    force_array = 'scalar';
else
    xxf = find(xx,1,'last'); % find the last argument, so that if both 'scalar' and 'array1' are passed, the later one is used
    force_array = varargin{xxf};
end

xx = ismember(varargin, {'never','row','column','always'});
if sum(xx) == 0
    vec_as_mat = 'never';
else
    xxf = find(xx,1,'last');
    vec_as_mat = varargin{xxf};
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

if isnumeric(val) || islogical(val)
    if ~isscalar(val) || strcmpi(force_array, 'array1')
        pyval = matarray2numpyarray(val, [], [], vec_as_mat);
    else
        % Scalar numeric values can be directly converted.
        pyval = val;
    end
elseif ischar(val)
    % strings can be directly converted
    pyval = val;
elseif iscell(val)
    pyval = cell2pylist(val, '', force_array, vec_as_mat);
elseif isstruct(val)
    pyval = struct2pydict(val, force_array, vec_as_mat);
else
    error('pyinterface:not_implemented','Unable to convert field of type "%s" into appropriate Python type', class(val));
end
end

