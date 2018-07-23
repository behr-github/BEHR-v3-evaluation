function [  ] = draw_tangent( x, y, slope, varargin )
%DRAW_TANGENT Draws a tangent line on the current graph
%   DRAW_TANGENT( X, Y, SLOPE ) will draw a line on the current plot at X,
%   Y with SLOPE as its slope. If X, Y, and SLOPE all have multiple values,
%   then each one will be drawn.
%
%   DRAW_TANGENT( X, Y, SLOPE, 'LENGTH', L ) the parameter 'LENGTH' can eb
%   used to alter how long the slope drawn will be; this is 1 by default
%   (0.5 units to either side of X).
%
%   DRAW_TANGENT( X, Y, SLOPE, NAME-VALUE PAIRS ) You may also pass name
%   value pairs that can be given to LINE to alter the appearance of the
%   tangent lines.

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

if numel(x) ~= numel(y) || numel(x) ~= numel(slope)
    E.badinput('X, Y, and SLOPE must all have the same number of elements')
end

if mod(numel(varargin),2) ~= 0
    E.badinput('One or more parameters are missing values')
end

xx = find(strcmpi('length', varargin),1,'last');
if isempty(xx)
    len = 1;
else
    len = varargin{xx+1};
    varargin(xx:xx+1) = [];
    if ~isscalar(len) || len <= 0
        E.badinput('The value of the parameter "LENGTH" must be a scalar and > 0')
    end
end

%%%%% MAIN FUNCTION %%%%%

for i=1:numel(x)
    x1 = x-len/2;
    x2 = x+len/2;
    y1 = y-slope*len/2;
    y2 = y+slope*len/2;
    
    line([x1, x2], [y1, y2], varargin{:});
end

end

