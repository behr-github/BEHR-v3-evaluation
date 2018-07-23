function [ bin_values, bin_midpoints, bin_error ] = bin_vertical_profile( altitude, data_values, binwidth, varargin )
%[bin_values, bin_midpoints, bin_error] = bin_vertical_profile(altitude, data_values, binwidth, [opt] binmode): Returns values of "data_values" by "altitude" binned to "binwidth" 
%   This function carries out binning of vertical profile data.  "Altitude"
%   should be a vector of altitudes in kilometers that correspond to the
%   values in "data_values."  "Binwidth" is the width of the bins in
%   kilometers.  By default, this function finds the bin value as the
%   median, with error presented as 25th and 75th quantiles.  To use mean
%   and std. deviation instead, enter 'mean' as the optional fourth
%   argument.
%
%   This function will return 3 variables: bin_values, the data values
%   median or mean for the bin; bin_midpoints, the midpoint altitude of the
%   bin; and bin_error, a 2 x n matrix of quartiles if running in median
%   mode, or a 1 x n row vector of std. dev. in mean mode.
%
%   Josh Laughner <joshlaugh5@gmail.com> 27 May 2014

p = inputParser;
p.addRequired('altitude',@isnumeric);
p.addRequired('data_values',@isnumeric);
p.addRequired('binwidth',@isscalar);
p.addOptional('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(altitude, data_values, binwidth, varargin{:});
pout = p.Results;
altitude = pout.altitude;
data_vals = pout.data_values;
binwidth = pout.binwidth;
binmode = pout.binmode;

% Define the edges of the bins

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
bin_midpoints = bin_midpoints;

% Do the binning. Points whose altitude falls exactly on a bin edge are
% assigned to the upper bin.
bin_values = zeros(size(bin_midpoints));
if strcmpi(binmode,'mean')
    bin_error = zeros(size(bin_midpoints));
else
    bin_error = zeros(length(bin_midpoints),2);
end

for a=1:numel(bin_values)
    bin_data_vals = data_vals(altitude > bin_edges(a) & altitude <= bin_edges(a+1));
    if strcmpi(binmode,'mean')
        bin_values(a) = nanmean(bin_data_vals(:));
        bin_error(a) = nanstd(bin_data_vals(:)) / sqrt(numel(bin_data_vals));
    else
        bin_values(a) = nanmedian(bin_data_vals(:));
        bin_error(a,:) = quantile(bin_data_vals(:),[0.25,0.75]);
    end
end

if size(bin_error,1)>size(bin_error,2); bin_error = bin_error'; end % Make bin_error into a row rather than column matrix
end

