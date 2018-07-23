function [ varargout ] = find_submatrix2( submatx, submaty, parmatx, parmaty, tolerance )
%FIND_SUBMATRIX2 Identifies the indicies where a smaller matrix exists in a larger one
%   This function searches parmat for the location of submat. This version
%   assumes that the matrices describe an x and y coordinate, thus there
%   needs to be two matrices input for the submat and two for the parmat
%   (parent matrix). An example of the application of this function is if
%   you have a subset of satellite pixels (which have a latitude and
%   longitude) that you want to find where they came from in the full array
%   of pixels. This assumes that the submatrices are contiguous in the
%   parent matrices, so it finds the location of the first and last points
%   of the submatrices in the parent matrices and assumes that those
%   describe a rectangular block that is the submatrix location.
%
%   The fifth input (tolerance) is optional. It describes the allowable
%   percent difference between the two values to be considered the same. It
%   defaults to 0.1%.
%
%   The output will be a 2x2 matrix where the first row describes the
%   beginning and ending indices in the first dimension, and the second row
%   the same for the second dimension. It will be empty if no match is
%   found.
%
%   Josh Laughner <joshlaugh5@gmail.com> 29 Jan 2016

E=JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if any([ndims(submatx), ndims(submaty), ndims(parmatx), ndims(parmaty)] ~= 2)
    E.badinput('All input matrices are expected to be 2D')
end
if any(size(submatx) ~= size(submaty))
    E.badinput('submatx and submaty must be the same size')
end
if any(size(parmatx) ~= size(parmaty))
    E.badinput('parmatx and parmaty must be the same size')
end 
if any(size(submatx) > size(parmatx))
    E.badinput('The submats must be smaller than the parmats')
end

if any(isnan([submatx(1); submatx(end); submaty(1); submaty(end)]))
    E.badinput('The first and last points in both SUBMATX and SUBMATY cannot be NaNs')
end

if ~exist('tolerance','var')
    tolerance = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Find where the sub matrix starts.
xx1 = abs((parmatx - submatx(1))/submatx(1)) * 100 < tolerance;
yy1 = abs((parmaty - submaty(1))/submaty(1)) * 100 < tolerance;

zz1 = xx1 & yy1;

% Find where it ends
xx2 = abs((parmatx - submatx(end))/submatx(end)) * 100 < tolerance;
yy2 = abs((parmaty - submaty(end))/submaty(end)) * 100 < tolerance;

zz2 = xx2 & yy2;

if sum(zz1(:)) == 0 || sum(zz2(:)) == 0
    nargoutchk(0,2);
    if nargout <= 1
        varargout{1} = [];
    elseif nargout == 2
        varargout = {[], []};
    end
    return
elseif sum(zz1(:)) > 1 || sum(zz2(:)) > 1
    E.callError('bad_solution','The starting and ending coordinates are not unique. Try a tighter tolerance than %f',tolerance);
end

[i1,j1] = find(zz1);
[i2,j2] = find(zz2);

parmatx_sub = parmatx(i1:i2, j1:j2);
parmaty_sub = parmaty(i1:i2, j1:j2);

if any(size(parmatx_sub) ~= size(submatx)) || any(size(parmaty_sub) ~= size(submaty))
    E.callError('bad_solution','The identified subset in parmat is not the same size as the submat')
elseif any(parmatx_sub(:) - submatx(:) > tolerance) || any(parmaty_sub(:) - submaty(:) > tolerance)
    E.callError('bad_solution','The identified subset in parmats is the same size, but it''s elements do not match to within %f%%',tolerance);
end

nargoutchk(0,2);
if nargout <= 1
    varargout{1} = [i1, i2; j1, j2];
elseif nargout == 2
    varargout = {i1:i2; j1:j2};
end
end

