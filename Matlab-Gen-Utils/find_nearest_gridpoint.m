function [ inds ] = find_nearest_gridpoint( x, y, Xgrid, Ygrid )
%FIND_NEAREST_GRIDPOINT( x, y, Xgrid, Ygrid)
%   Finds the indices in Xgrid, Ygrid closest to x, y. Returns those
%   indices in a single 2x1 vector.

E = JLLErrors;

%%%%% INPUT CHECKING %%%%%

if ndims(Xgrid) ~= ndims(Ygrid) || any(size(Xgrid) ~= size(Ygrid))
    E.badinput('Xgrid and Ygrid must be the same size')
end
if ~isscalar(x) || ~isscalar(y)
    E.badinput('x and y are expected to be scalar')
end
if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(Xgrid) || ~isnumeric(Ygrid)
    E.badinput('All inputs are expected to be numeric')
end

%%%%% MAIN FUNCTION %%%%%
r2 = (Xgrid - x).^2 + (Ygrid - y).^2;
inds = cell(ndims(Xgrid),1);
[inds{:}] = find(r2 == min(r2(:)));
inds = cell2mat(inds);

end

