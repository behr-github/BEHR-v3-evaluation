function [ vout ] = griddata_nointerp( x, y, v, xq, yq )
%griddata_nointerp Grids data to xq and yq without interpolation
%   The built-in function griddata will grid data but interpolates missing
%   values to give a complete output value grid.  This function will grid
%   data without interpolation, so that missing values are preserved.  For
%   this to work well, the query grids must be such that change only occurs
%   along one dimension; e.g. this is okay:
%       xq = [1 2 3     and     yq = [1 1 1
%             1 2 3                   2 2 2
%             1 2 3]                  3 3 3]
%
%   but this is not:
%       xq = [1 2 3             yq = [1 2 3
%             2 3 4                   2 3 4
%             3 4 5]                  3 4 5] 
%
%   This function is only meant for 2D matrices
%
%   This function will assign values to the nearest vertex in (xq,yq).
%
%   Also, the inputs xq and yq must be matrices; this function will not
%   automatically generate a mesh grid based on input vectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that xq and yq are matrices of the same size and that x, y, and v
% are the same size
if isvector(xq) || isvector(yq);
    error('griddata_nointerp:query_input', 'xq and yq must be matrices, not vectors');
elseif ~all(size(xq)==size(yq));
    error('griddata_nointerp:query_input', 'xq and yq must have the same dimensions');
elseif any(size(x) ~= size(y)) || any(size(x) ~= size(v))
    error('griddata_nointerp:data_input', 'x, y, and v must have the same dimensions');
end

%Check that one dimension in xq and yq doesn't change - this is the best I
%can do towards checking that the matrices are meshgrids
dx1 = diff(xq,1,1); dx2 = diff(xq,1,2);
dx1 = sum(dx1,2); dx2 = sum(dx2,1);
dy1 = diff(yq,1,1); dy2 = diff(yq,1,2);
dy1 = sum(dy1,2); dy2 = sum(dy2,1);

if ~(all(dx1==0) || all(dx2==0)) || ~(all(dy1==0) || all(dy2==0))
    error('griddata_nointerp:query_input','One dimension of xq and yq must be constant');
end

% Check that xq and yq have the change occurring in different dimensions
test = [all(dx1==0), all(dx2==0); all(dy1==0), all(dy2==0)];
if any(all(test,1))
    error('griddata_nointerp:query_input','xq and yq should have changing values in different dimensions');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTRUCTING CALUCATION MATRICES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vals = nan(size(xq));
count = zeros(size(xq));

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

% Remove any values that lie outside of the query matrices
xmin = min(xq(:)); xmax = max(xq(:));
ymin = min(yq(:)); ymax = max(yq(:));

xx = x >= xmin & x < xmax;
yy = y >= ymin & y < ymax;

x = x(xx & yy);
y = y(xx & yy);
v = v(xx & yy);

% Iterate through each element in (x,y,v) and assign the value to the
% nearest vertex in (xq,yq)

for a=1:numel(v)
    delta = sqrt((xq - x(a)).^2 + (yq - y(a)).^2);
    [~,ind] = min(delta(:)); %find the smallest distance 
    
    oc = count(ind); %"old count" used for running average
    count(ind) = count(ind)+1;
    vals(ind) = (nansum([vals(ind)*oc, v(a)]))/(oc+1);
end

vout = vals;

end

