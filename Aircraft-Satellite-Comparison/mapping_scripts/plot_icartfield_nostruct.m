function [ cb ] = plot_icartfield_nostruct( lon, lat, data, varargin )
%plot_flightpath(Merge_str): Plots the flight path from lon and lat data passed, also coloring the line by the data passed.
%   This function will automatically plot the flight path based on lon and
%   lat data given.  This is intended to more easily plot subsets of ICART
%   files, by allowing you to extract the data you want from an ICART merge
%   and only pass the points you wish to see plotted.  For full day
%   plotting, the regular plot_icartfield function accepts Merge data
%   structure produced by read_merge_data.
%
%   The first two arguments are the longitude and latitude values.
%
%   The third argument is the data you want to plot.
%
%   The optional fourth argument will reduce the number of points plotted,
%   e.g. (lon, lat, data, 4) would only plot every 4th point.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 May 2014

p = inputParser;
p.addRequired('lon');
p.addRequired('lat');
p.addRequired('data');
p.addOptional('skip',1,@isscalar);
p.addParamValue('lonbdy',[],@ismatrix);
p.addParamValue('latbdy',[],@ismatrix);

p.parse(lon,lat,data, varargin{:})
pout = p.Results;
skip = pout.skip;
latbdy = [min(pout.latbdy), max(pout.latbdy)];
lonbdy = [min(pout.lonbdy), max(pout.lonbdy)];

% Read in the data, lat, and lon data:
p_alt = pout.data;
lat = pout.lat;
lon = pout.lon;
if any(lon>180)
    lon = lon - 360; %Since m_map uses the convention that W is < 0, we must make the conversion if the user inputs longitude values in the original ICART style
end

% Get the current colormap; we'll need it to color the line
cm = colormap;
cml = size(cm,1); 

% Reduce the size of the plotted vectors
tmp = true(size(lat));
tmp(1:skip:end) = 0;

lat(tmp) = [];
lon(tmp) = [];
p_alt(tmp) = [];

% Calculate the lat/lon boundaries with some padding
if isempty(latbdy); latbdy = [floor(min(lat))-0.5, ceil(max(lat))+0.5]; end
if isempty(lonbdy); lonbdy = [floor(min(lon))-0.5, ceil(max(lon))+0.5]; end

% Calculate an integer version of pressure altitude that has its lowest
% value at 1 and its highest as the first dimension size of the colormap.
% This is necessary so that these can serve as indices to the colormap.
p_alt64 = ceil((p_alt - min(p_alt)) / max(p_alt - min(p_alt)) * (cml-1)) + 1;

% Plot things
m_proj('Mercator','long',lonbdy,'lat',latbdy);
m_coast('color','k');
m_states('k');
m_grid('linestyle','none')


for a = 1:numel(lat)
    m_line(lon(a), lat(a), 'color',cm(p_alt64(a),:)); 
end

% Now, unfortunately since we've done a custom implementation of coloring,
% we'll need to create our own tick marks. We'll want each tick mark to
% represent the highest pressure that would have that color.

cb = colorbar;
ticks = get(cb,'Ticks'); nTicks = numel(ticks);

% We're just going to relabel our existing ticks, so first we find our
% colormap indices that will correspond to those ticks
tick_indices = floor(linspace(1,cml,nTicks));

% Then we want to find the highest pressure that corresponds to each of the
% indices
pressures = cell(size(tick_indices));
for b = 1:numel(tick_indices)
    tmp = p_alt(p_alt64 == tick_indices(b));
    pressures{b} = num2str(round(max(tmp)));
end
%pressures = fliplr(pressures); % This will put the highest pressure as the lowest tick on the colorbar
set(cb,'TickLabels',pressures);

end

