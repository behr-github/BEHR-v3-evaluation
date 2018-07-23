function [  ] = highlight_plot( ranges, varargin  )
%HIGHLIGHT_PLOT Highlight ranges of x or y values on a plot.
%   HIGHLIGHT_PLOT(ranges) will draw a patch object for each row in the
%   matrix ranges. The first column defines the bottom of the range, the
%   second the top. The patches are drawn with an alpha of 0.5 by default.
%
%   HIGHLIGHT_PLOT(ranges, name-value property pairs) allows you to define
%   the properties of the patches. Any patch property will work, plus two
%   additionally defined properties:
%
%       'color' - defines the overall color of the patches, i.e. it is
%       passed as C in patch(X,Y,C). Defaults to red.
%
%       'axis' - 'x' or 'y', indicates if the values in ranges are x or y
%       values. Defaults to x.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 June 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(ranges)
    E.badinput('ranges must be a numeric matrix')
elseif size(ranges, 2) ~= 2
    E.badinput('ranges must be an n-by-2 matrix')
elseif any(ranges(:,1) >= ranges(:,2))
    E.badinput('The left column of ranges is expected to be less than the right')
end

% varargin should have an even number of elements because it should
% consist of name-value pairs only.

if mod(numel(varargin),2) ~= 0
    E.badinput('One or more input properties is missing a value')
end

% Look for the two properties we need to handle separately
properties_in = varargin;
xx = strcmpi('color',properties_in);
if any(xx)
    f = find(xx);
    patch_color = properties_in{f+1};
    properties_in(f:f+1) = [];
    % Not validating this, we'll let patch() do it
else
    patch_color = 'r';
end

xx = strcmpi('axis',properties_in);
if any(xx)
    f = find(xx);
    range_axis = lower(properties_in{f+1});
    if ~ismember(range_axis,{'x','y'})
        E.badinput('The property ''range_axis'' must be ''x'' or ''y''');
    end
    properties_in(f:f+1) = [];
else 
    range_axis = 'x';
end

% We'll need to know if the user is setting a different parent set of axes
% (to get the right limits to draw the patches) and if the user wants a
% different FaceAlpha, since it defaults to 0.5.
xx = strcmpi('parent',properties_in);
if any(xx)
    f = find(xx);
    parent = properties_in{f+1};
    % do pass this along to the patch() function
else
    parent = gca;
end

xx = strcmpi('FaceAlpha',properties_in);
if any(xx)
    f = find(xx);
    facealpha = properties_in{f+1};
    properties_in(f:f+1) = [];
    % this we'll be passing explicitly, so remove it
else 
    facealpha = 0.5;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DRAW PATCHES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
switch range_axis
    case 'x'
        lims = get(parent,'ylim');
        py = [lims(1), lims(1), lims(2), lims(2)];
    case 'y'
        lims = get(parent, 'xlim');
        px = [lims(1), lims(1), lims(2), lims(2)];
end
    
for a=1:size(ranges,1)
	r = ranges(a,:);
    switch range_axis
        case 'x'
            px = [r(1), r(2), r(2), r(1)];
        case 'y'
            py = [r(1), r(2), r(2), r(1)];
    end
    
    patch(px, py, patch_color, 'facealpha', facealpha, properties_in{:});
    
end

end

