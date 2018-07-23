function [ varargout ] = scatter_errorbars( x,y,l,varargin )
%xerrorbars - Plots error bars acceptable for a scatter plot in x or y
%directions.
%   Valid inputs are (x,y,l,u) or (x,y,e)   
%
%   This function takes the x & y coordinates of the data points, along
%   with either upper and lower bounds for the error (in distance from the
%   point) or a single value for symmetric error.  It will plot just error
%   bars, with no connecting lines, and will not overwrite an existing
%   plot.
%
%   Parameters:
%       'direction' = 'x' or 'y' to indicate which axis the error bars
%       should be plotted on. Defaults to 'y'.
%
%       'color' = Any valid colorspec, will recolor the error bar lines.
%
%       'linewidth' = The width of the error bar lines in pixels.
%
%       'tipscale' = a scale factor for the cross bars on the tips of the
%       errorbars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     INPUT PARSING     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addRequired('x',@ismatrix);
p.addRequired('y',@ismatrix);
p.addRequired('l',@ismatrix);
p.addOptional('u',[],@isnumeric);
p.addParameter('direction','y',@(x) any(strcmpi(x,{'x','y'})));
p.addParameter('color','b');
p.addParameter('linewidth',0.5);
p.addParameter('parent',gca);
p.addParameter('tipscale',1,@isscalar);

p.parse(x,y,l,varargin{:});
pout = p.Results;
x = pout.x;
y = pout.y;
l = pout.l;
u = pout.u;
direction = pout.direction;
color = pout.color;
width = pout.linewidth;
parent = pout.parent;
tipscale = pout.tipscale;

if isempty(u); u=l; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     PLOTTING     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(direction,'x')
    err_axis = x; ind_axis = y;
elseif strcmpi(direction,'y')
    err_axis = y; ind_axis = x;
end

% These values will be necessary for plotting the errorbars
tee = (max(ind_axis(:))-min(ind_axis(:)))/100 * tipscale;  % make tee .02 x-distance for error bars
tip_l = ind_axis - tee;
tip_r = ind_axis + tee;
emax = err_axis + u;
emin = err_axis - l;
n = numel(err_axis);

% Following the MATLAB errorbar.m function, make an n by 9 matrix of
% NaN-separated line segments that will draw out the error bars; [1:2] are
% the main bar, [4:5] the lower end cross bar, and [7:8] the upper end
% cross bar.

errpts = zeros(9,n);
errpts(1,:) = emin(:);
errpts(2,:) = emax(:);
errpts(3,:) = NaN;
errpts(4,:) = emin(:);
errpts(5,:) = emin(:);
errpts(6,:) = NaN;
errpts(7,:) = emax(:);
errpts(8,:) = emax(:);
errpts(9,:) = NaN;

indpts = zeros(9,n);
indpts(1,:) = ind_axis(:);
indpts(2,:) = ind_axis(:);
indpts(3,:) = NaN;
indpts(4,:) = tip_l(:);
indpts(5,:) = tip_r(:);
indpts(6,:) = NaN;
indpts(7,:) = tip_l(:);
indpts(8,:) = tip_r(:);
indpts(9,:) = NaN;

if strcmpi(direction,'x');
    l=line(errpts(:), indpts(:), 'color', color, 'linewidth', width, 'parent', parent);
elseif strcmpi(direction, 'y');
    l=line(indpts(:), errpts(:), 'color', color, 'linewidth', width, 'parent', parent);
end

if nargout > 0
    varargout{1} = l;
end
end

