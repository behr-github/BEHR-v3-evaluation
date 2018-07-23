function [ A ] = numpyarray2matarray( N, dimorder )
%NUMPYARRAY2MATARRAY Convert a Matlab array into a Python numpy array
%   A = NUMPYARRAY2MATARRAY( N ) Converts the Numpy array, N, into a Numpy
%   array A. Dimension order is treated so that A(:,1,1) == N[:,0,0].
%
%   A = NUMPYARRAY2MATARRAY( N, 'match' ) behaves the same as the first
%   method.
%
%   A = NUMPYARRAY2MATARRAY( N, 'native' ) Retains the native Python
%   dimension order in A, such that A(:,1,1) == N[0,0,:].

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(N, 'py.numpy.ndarray')
    error('pyinterface:badinput','N should be a Numpy ndarray')
end

if ~exist('dimorder','var')
    dimorder = 'match';
else
    allowed_orders = {'match', 'native'};
    if ~any(strcmpi(dimorder, allowed_orders))
        error('pyinterface:badinput','DIMORDER (if given) must be one of %s', strjoin(allowed_orders, ', '));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Credit to https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double
% David Laurenson for most of this, though the final permutation is handled
% differently.

dsize = cellfun(@int64, cell(N.shape));
A = double(py.array.array('d', py.numpy.nditer(N, pyargs('order', 'C'))));
if numel(dsize) > 1
    % Python shapes (unlike Matlab's sizes) will actually only have 1
    % element if the array is 1D. So if A is a 10x1 vector, while Matlab's
    % size(A) would return [10, 1], Python's A.shape would only return 10.
    % At this point, A will always be a vector, so we only need to reshape
    % if there is actual >1D shape to it. (Also RESHAPE cannot accept a
    % size vector with one element.)
    A = reshape(A, fliplr(dsize));
end

if strcmpi(dimorder, 'match')
    A = permute(A, ndims(A):-1:1);
end

if N.dtype == py.numpy.dtype('bool')
    A = logical(A);
end

end

