function [ A ] = nansum2( A, varargin )
%nansum2 Version of nansum that returns NaN if all elements along the
%summation are NaN
%   The built in nansum function will return 0 if all elements summed are
%   NaNs.  This function will instead return a NaN if all elements summed
%   are NaNs; this is useful in e.g. a running average of a quantity where
%   NaNs are used as fill values. Takes a matrix or vector and optionally a
%   dimension to sum along.  Vectors are always summed along their
%   non-singleton dimension, matrices default to summation along the first
%   dimension.

% Figure out what dimension to sum along
if nargin > 1
    dim = varargin{1};
elseif isvector(A)
    [~,dim] = max(size(A));
else
    dim = 1;
end

% Find instances where the terms summed are all NaNs
xx = all(isnan(A),dim);

% Do the sum; we'll adjust for all nans in a bit
A = nansum(A,dim);

% Set the terms where all the summed terms were NaNs to NaN themselves
A(xx) = nan;


end