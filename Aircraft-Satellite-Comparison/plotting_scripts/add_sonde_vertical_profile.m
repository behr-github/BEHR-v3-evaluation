function [hax1, hax2] = add_sonde_vertical_profile(Merge_in, data_field, varargin)
% add_vertical_profile(Merge, data_field): Adds a vertical profile of "data_field" in "Merge" to the current figure.
%	This function will add a second y-axis to an existing vertical profile
%	plot.  Requires a merge data structure created using "read_merge_file.m" and a
%	data field from that structure passed as a string. As an optional argument, pass
%	the figure handle of the figure to add the second y axis to.  Parameter
%	arguments: 
%
%       'binwidth' = width of the altitude bins in km.  If not specified,
%       all points will be plotted with no binning.
%       
%       'binmode' = sets if the binned value will be a 'median' with 25th and
%       75th quartile or 'mean' with std. error. Median is default.

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('data_field',@isstr);
p.addOptional('fighandle',-1,@isscalar);
p.addParamValue('binwidth',-1,@isscalar);
p.addParamValue('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(Merge_in,data_field,varargin{:});
pout = p.Results;

Merge = pout.Merge;
field = pout.data_field;
fignum = pout.fighandle;
binwidth = pout.binwidth;
binmode = pout.binmode;

% If the user passes a valid figure handle, select that figure.
if fignum > 0; 
	figure(fignum);
end

% Find the entries that have the desired site flag and profile number

pressures = Merge.Data.Pressure_hPa.Values;
data_vals = eval(sprintf('Merge.Data.%s.Values',field));

% Replace any fill values, upper LOD, or lower LOD values with NaNs
fill_val = eval(sprintf('Merge.Data.%s.Fill',field));
pres_fill_val = Merge.Data.Pressure_hPa.Fill;
ULOD_val = Merge.metadata.upper_lod_flag;
LLOD_val = Merge.metadata.lower_lod_flag;

data_vals(data_vals == fill_val) = NaN;
data_vals(data_vals == ULOD_val) = NaN;
data_vals(data_vals == LLOD_val) = NaN;

xx = find(pressures == pres_fill_val);
data_vals(xx) = []; pressures(xx) = [];

% Plot the values converting pressure to altitude, using P0 = 1013.25 hPa
% and scale height H = 7.4 km

altitude = -log(pressures ./ 1013.25) .* 7.4;

% Get the existing y-limits
old_ylim = get(gca,'YLim');

% Create the second set of axes, with the xaxis at the top of the graph
hold on;
hax1 = gca;
ax_pos = get(gca,'Position');
set(gca,'box','off') %Turn of the box so that tick marks don't appear on both sides.
hax2 = axes('Position', ax_pos,'XAxisLocation','top','YLim',old_ylim,'YTick',[],'Color','none');

% If binwidth is negative (i.e. not specified by the user), plot every data
% point.  If a bin width is specified, average together all values within
% each bin slice, then plot that average at the midpoint altitude.
if binwidth < 0;
    scatter(hax2, data_vals, altitude);
else
    top_alt = max(altitude(:));
    bin_edges = 0:binwidth:top_alt;
    
    % The colon operator will only include the end value if it falls along
    % the interval given; i.e. 0:2:3 will give [0 2] but 0:2:4 gives [0 2 4]
    % So we need to check what the last value of our bin_edges matrix is;
    % if it is not our top altitude, we want to append the previous
    % entry plus the binwidth.
    if bin_edges(end) ~= top_alt;
        bin_edges(end+1) = bin_edges(end) + binwidth;
    end
    
    % Calculate the midpoints
    bin_midpoints = bin_edges(1:end-1) + diff(bin_edges);
    
    % Do the binning. Points whose altitude falls exactly on a bin edge are
    % assigned to the upper bin.
    bin_avgs = zeros(size(bin_midpoints));
    if strcmpi(binmode,'mean') % For now, I'm keeping the bin_error variable in until I have a way to plot x-errorbars
        bin_error = zeros(size(bin_midpoints));
    else
        bin_error = zeros(length(bin_midpoints),2);
    end
    
    for a=1:numel(bin_avgs)
        bin_data_vals = data_vals(altitude > bin_edges(a) & altitude <= bin_edges(a+1));
        if strcmpi(binmode,'mean')
            bin_avgs(a) = nanmean(bin_data_vals(:));
            bin_error(a) = nanstd(bin_data_vals(:)) / sqrt(numel(bin_data_vals));
        else
            bin_avgs(a) = nanmedian(bin_data_vals(:));
            bin_error(a,:) = quantile(bin_data_vals(:),[0.25,0.75]);
        end     
    end
    
    % Remove any NaNs in the average
    xx = find(isnan(bin_avgs));
    bin_avgs(xx) = []; bin_midpoints(xx) = []; bin_error(xx,:) = [];
    
    % Make a line plot of the vertical profile with the standard error
    % plotted as an envelope
    
    line(bin_avgs, bin_midpoints,'Parent',hax2,'Color','r');
end
end