function l = list_recursion(A, varargin)
% LIST_RECURSION Turns a multidimension array into nested Python lists
%   L = LIST_RECURSION(A) Takes a Matlab array, A, and converts it into
%   nested Python lists. Each slice along the first dimension will be the
%   inner-most lists, and successive dimensions will be the next level out.
%   A may be a numeric or cell array. If a cell array, MATLAB2PYTHON will
%   be called on each element to convert it to the appropriate python type.
%
%   Examples:
%       list_recursion([1 2; 3 4])
%       = [[1.0, 3.0], [2.0, 4.0]]
%
%       A = cat(3, [1 2 3; 4 5 6], [10 20 30; 40 50 60])
%       list_recursion(A)
%       = [[[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]], 
%          [[10.0, 40.0], [20.0, 50.0], [30.0, 60.0]]]
%
%   L = LIST_RECURSION( A, BASE_DIM, BASE_SIZE ) should not usually be used
%   directly, rather this is usually used during the internal recursion.
%   BASE_DIM is a scalar number indicating which dimension LIST_RECURSION
%   is operating on, it reduces by 1 in each recursive call. BASE_SIZE is
%   the size that the final Numpy array should have, it must be equivalent
%   to the size of A. When these are not passed, they are initialized to
%   BASE_DIM = ndims(A) and BASE_SIZE = size(A).
%
%   Additional parameters:
%       'vec_as_mat' - controls how vectors passed into this function are
%       treated. Unlike Matlab vectors, which are just 2D arrays with one
%       dimension having length 1, Numpy arrays can actually be strictly
%       one dimensional. This parameter can have values 'never', 'column',
%       'row', or 'always'. 'never' means that any vector passed in will be
%       converted to a 1D numpy array (default). 'row' and 'column' mean
%       that row or column vectors, respectively will be converted a 2D
%       array with length 1 in one dimension, the other will become a 1D
%       numpy array. 'always' means that any vector passed in becomes a 2D
%       numpy array with length 1 in one dimension.
%
%       'force_array' - passed through to MATLAB2PYTHON. Pass the string
%       'array1' as the value to force scalar values to be created as size
%       1, 1D numpy arrays.

p = inputParser;
p.addOptional('base_dim', []);
p.addOptional('base_size', []);
p.addParameter('vec_as_mat', 'never', @ischar);
p.addParameter('force_array', 'scalar', @ischar);

p.parse(varargin{:});
pout = p.Results;

base_size = pout.base_size;
base_dim = pout.base_dim;
vector_as_matrix = pout.vec_as_mat;
force_array = pout.force_array;

allowed_vam = {'never', 'column', 'row', 'always'};
if ~ischar(vector_as_matrix) || ~ismember(vector_as_matrix, allowed_vam)
    error('pyinterface:bad_input', 'The parameter "vec_as_mat" must be one of the strings: %s', strjoin(allowed_vam, ', '));
end

base_inputs = ~isempty(base_size) + ~isempty(base_dim);
if base_inputs == 1
    error('pyinterface:bad_input', 'LIST_RECURSION requires both or neither of BASE_DIM and BASE_SIZE to be specified; you cannot specify only one');
elseif base_inputs < 1
    base_size = size(A);
    base_dim = numel(base_size);
else
    if ~isscalar(base_dim) || ~isnumeric(base_dim) || base_dim < 1 || base_dim > numel(base_size)
        error('pyinterface:bad_input', 'BASE_DIM must be a scalar number between 1 and numel(BASE_SIZE)');
    elseif ~isnumeric(base_size) || ~isvector(base_size)
        error('pyinterface:bad_input', 'BASE_SIZE must be a numeric vector');
    end
end

if isnumeric(A)
    A = num2cell(A);
elseif ~iscell(A)
    error('pyinterface:bad_input', 'LIST_RECURSION is only defined for A as a numeric array or cell array');
end

if isvector(A) && (base_dim == 1 || strcmpi(vector_as_matrix, 'never') || (strcmpi(vector_as_matrix, 'column') && iscolumn(A)) || (strcmpi(vector_as_matrix, 'row') && isrow(A)))
    % We should enter this branch under the following conditions:
    %   A) The input array is a vector and either:
    %       i) It is the last dimension of the original input array
    %       ii) We are at the top level, so the original input was a
    %       vector, and the user has not requested that vectors be retained
    %       as matrices.
    % In either case, we need to actually put the elements of A into a
    % list, so we will iterate over the elements of A, converting them into
    % Python types, and then converting the resulting cell array to a
    % Python list.
    for i=1:numel(A)
        A{i} = matlab2python(A{i}, force_array);
    end
    
    if iscolumn(A)
        % py.list requires A to be a row vector
        A = A';
    end
    l = py.list(A);
else
    % We should enter this branch if we want to treat A as an array with
    % more than one dimension (i.e. not a vector). This distinction is
    % important because in Matlab all values have at least 2 dimensions;
    % vectors just have one of those as a length 1 dimension. In Numpy, an
    % array can have only 1 (or even 0!) dimensions.
    sz = base_size(1:base_dim);
    A2 = reshape(A, prod(sz(1:end-1)), sz(end));
    l = py.list;
    for i=1:size(A2,2)
        % So what we need to do is iterate over A, taking out each subarray
        % along the last dimension, i.e. for a 3D array, A(:,:,i). Using
        % this reshaping approach allows us to work with an array with an
        % arbitrary number of dimensions. (We work in this order because
        % Numpy arrays use C-style ordering, while Matlab arrays use
        % Fortran style ordering; working "backwards" from our perspective
        % makes the dimensions match up).
        if base_dim == 2
            target_sz = [sz(1), 1];
        else
            target_sz = sz(1:end-1);
        end
        Aslice = reshape(A2(:,i), target_sz);
        % We recursively pass this subarrays into list_recursion so that it
        % continues to extract subarrays until it gets to the last
        % dimension, at which point it actually starts creating lists from
        % the individual vectors. We tell it to always treat vectors as
        % matrices in order to preserve singleton dimensions. I.e., if we
        % didn't do that, an array of size (4,1,3) would get squashed into
        % a (4,3) numpy array. 
        l.append(list_recursion(Aslice, base_dim-1, base_size, 'vec_as_mat', 'always', 'force_array', force_array));
    end
end

end