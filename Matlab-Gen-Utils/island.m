function [ isLand ] = island( lon, lat, varargin )
%island Test is points defined by a pair of lon/lat vectors are over land
%   Returns a logical vector of the same length as the input lon/lat
%   vectors that will be 1 is the point is over land and 0 if the point is
%   over the ocean.  Modified from code posted at 
%
%   http://www.mathworks.com/matlabcentral/answers/1065-determining-whether-a-point-on-earth-given-latitude-and-longitude-is-on-land-or-ocean
%
%   by Brett Shoelson.
%
%   By default, this will do a 0.1 degree grid around the continental US.
%   To change this, use the following parameters:
%
%       density = # of points per degree latitude and longitude. Must be
%       scalar, defaults to 10.
%
%       latlim, lonlim = sets the borders of the matrix used to test
%       whether something is land or ocean.  Defaults to [20, 55] and
%       [-130, -60].
%
%   Making either the density or area smaller will speed up the
%   calculation.
%
%   Josh Laughner <joshlaugh5@gmail.com> 24 Sept 2014


p = inputParser;
p.addParameter('density',10,@isscalar);
p.addParameter('lonlim',[-130,-60],@(x) numel(x)==2);
p.addParameter('latlim',[20,55], @(x) numel(x)==2);
p.parse(varargin{:});
pout=p.Results;
density = pout.density;
lonlim = pout.lonlim;
latlim = pout.latlim;

E = JLLErrors;

% Input validation
narginchk(2,8);

if ~all(size(lon)==size(lat))
    E.dimMismatch('lon','lat')
end

if density < 1;
    E.callError('bad_input','''density'' must be a scalar value of 1 or greater.');
elseif lonlim(1) > lonlim(2) || latlim(1) > latlim(2)
    E.callError('bad_input','latlim and lonlim must have the smaller value first');
end

% Calculate
coast = load('coast.mat');
[Z, R] = vec2mtx(coast.lat, coast.long, density, latlim, lonlim, 'filled');
val = ltln2val(Z, R, lat, lon);
isLand = val == 0;

end

