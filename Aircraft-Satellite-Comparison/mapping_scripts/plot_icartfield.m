function [ cb ] = plot_icartfield( Merge_str, data_field, varargin )
%plot_flightpath(Merge_str): Plots the data (as color) and flight path from an ICART merge structue
%   This function will automatically plot the flight path of an ICART merge
%   file, coloring the line by the value of whatever field is passed.  
%
%   The first argument must be a Merge data structure created by the
%   read_merge_data script.
%
%   The second argument is a data field in the Merge data structure, passed
%   as a string.
%
%   The optional second argument will reduce the number of points plotted,
%   e.g. (Merge, 4) would only plot every 4th point.
%
%   The parameter value 'cbrange' allows you to force the colorbar to a
%   specific range.  Since the color mapping is done custom and cannot be
%   changed with caxis, this is the only way to override the default
%   mapping.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 May 2014

p = inputParser;
p.addRequired('Merge_str');
p.addRequired('data_field',@isstr);
p.addOptional('skip',1,@isscalar);
p.addParamValue('cbrange',[],@ismatrix);

p.parse(Merge_str, data_field, varargin{:})
pout = p.Results;
Merge= pout.Merge_str;
field = pout.data_field;
skip = pout.skip;
cbrange = pout.cbrange;

if ~isfield(Merge.Data,field)
    error('plot_icart:not_field','%s is not a data field in the Merge data structure.',field);
end

% Read in the data, lat, and lon data:
p_alt = eval(sprintf('Merge.Data.%s.Values',field));
lat = Merge.Data.LATITUDE.Values;
lon = Merge.Data.LONGITUDE.Values;
lon = lon - 360; %Since m_map uses the convention that W is < 0, we must make the conversion.

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
latbdy = [floor(min(lat))-0.5, ceil(max(lat))+0.5];
lonbdy = [floor(min(lon))-0.5, ceil(max(lon))+0.5];

% Calculate an integer version of pressure altitude that has its lowest
% value at 1 and its highest as the first dimension size of the colormap.
% This is necessary so that these can serve as indices to the colormap.
if isempty(cbrange)
    min_palt = min(p_alt); max_palt = max(p_alt);
else
    min_palt = min(cbrange); max_palt = max(cbrange);
end

p_alt64 = ceil((p_alt - min_palt) / max_palt * (cml-1)) + 1;
p_alt64 = max(p_alt64,1); p_alt64 = min(p_alt64,cml); % Clamp the values of p_alt64 such that they will always be a valid index of the colormap
% Plot things
m_proj('Mercator','long',lonbdy,'lat',latbdy);
m_coast('color','k');
m_states('k');
m_grid('linestyle','none')
title(sprintf('Flight path on %s ', Merge.metadata.date),'fontsize',18)

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
cblabel = sprintf('%s (%s)',field,eval(sprintf('Merge.Data.%s.Unit',field)));
ylabel(cb,cblabel,'fontsize',16)
end

