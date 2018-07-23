function [q_in, q_on] = inpolygon_mg(xq, yq, xv, yv)
%INPOLYGON_MG Faster inpolygon() if XQ and YQ are on a mesh grid
%   [ Q_IN, Q_ON ] = INPOLYGON_MG( XQ, YQ, XV, YV ) Returns
%   logical arrays Q_IN and Q_ON the same size as XQ and YQ
%   that a true for points defined by XQ and YQ that are in
%   and on the border, respectively, of the polygon defined 
%   by XV and YV. Internally, this uses INPOLYGON, but assumes
%   that XQ and YQ define a meshgrid, i.e. 
%
%       [XQ, YQ] = meshgrid(X, Y)
%
%   where X and Y are vectors. Specifically, this function first
%   cuts down XQ and YQ by comparing the first row of XQ and 
%   first column of YQ to the range of XV and YV, respectively.
%   Only the submatrices of XQ and YQ that have points inside
%   that range are passed to INPOLYGON. In simple tests, when
%   the polygon is an order of magnitude smaller than the query
%   field, this function is an order of magnitude faster than
%   INPOLYGON. However, passing the arrays seems to incur some 
%   overhead, so this is usually slower than cutting down the 
%   arrays in the calling function, then passing the smaller
%   arrays to INPOLYGON.

E = JLLErrors;

if ~isequal(size(xq), size(yq))
    E.badinput('XQ and YQ must be the same size');
end

if any(diff(xq(:,1),[],1) ~= 0) || any(diff(yq(1,:),[],2) ~= 0)
    E.badinput('XQ and YQ do not appear to be a mesh grid (XQ is not constant in the first dim/YQ is not constant in the second');
end
    

if any(isnan(xv)) || any(isnan(yv))
    warning('INPOLYGON_MG has not been tested for polygons defined with inner and outer loops');
end

xx = xq(1,:) >= min(xv) & xq(1,:) <= max(xv);
yy = yq(:,1) >= min(yv) & yq(:,1) <= max(yv);

q_in = false(size(xq));
q_on = false(size(xq));

[q_in(yy,xx), q_on(yy,xx)] = inpolygon(xq(yy,xx), yq(yy,xx), xv, yv);

end
