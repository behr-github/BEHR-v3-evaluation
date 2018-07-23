function [ ] = reset_axis_pos( varargin )
%reset_axis_pos([opt] fighandle): Fixes a second axis so that it is scaled correctly 
%   This function looks for all non-legend axes in a figure and
%   finds the one that has its x-axis on the bottom and its y-axis on the
%   left.  This is assumed to be the first axis for a figure, the one that
%   scales correctly when a figure window is resized.  All other axes are
%   then resized to match this one.  Optionally, a figure handle can be
%   passed, if none is, the current figure is assumed.

p = inputParser;
p.addOptional('fighandle',gcf);
p.parse(varargin{:});
pout = p.Results;
fighandle = pout.fighandle;

haxes = findall(fighandle,'Type','Axes');
hax = haxes(~ismember(get(haxes,'Tag'),{'Legend'}));

xaxloc = get(hax,'XAxisLocation');
yaxloc = get(hax,'YAxisLocation');

xaxbool = strcmpi('bottom',xaxloc);
yaxbool = strcmpi('left',yaxloc);

hax_first = hax(xaxbool & yaxbool);
pos = get(hax_first,'Position');

for a=1:numel(hax)
    if hax(a) ~= hax_first
        set(hax(a),'Position',pos);
    end
end

end

