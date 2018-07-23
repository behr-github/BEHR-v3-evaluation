function [ bin_values, bin_midpoints, bin_errors, bin_edges ] = bin_omisp_pressure( pressure, data_values, varargin )
%bin_omisp_pressure(): Returns 'data' binned to bins corresponding to OMI pressure bins.
%   This function aims to replicate the binning described in Bucsela et.
%   al. (J. Geophys. Res. 2008, 113, D16S31), which bins data by pressure
%   to match the 25 lowest pressures that BAMFs are calculated at for the
%   OMI algorithm (see OMNO2 Readme File, Document Version 6.7, 10 Jan
%   2013, p. 6).
%
%   It requires only the data to be binned and the corresponding PRESSURE
%   (not altitude) values.  It will return the center pressure of each bin,
%   the median (or mean) value for data, and a measure of error in that
%   bin. 
%
%   By default, the data value for each bin is a median value and the
%   errors will be given as 25th and 75th quartiles.  To output mean and
%   std. deviation instead, pass 'mean' as the optional third argument.

p = inputParser;
p.addRequired('pressure',@isnumeric);
p.addRequired('data_values',@isnumeric);
p.addOptional('binmode','median',@(x) any(strcmpi(x,{'median','mean','binonly'})));


p.parse(pressure, data_values, varargin{:});
pout = p.Results;
pressure = pout.pressure;
data_vals = pout.data_values;
binmode = pout.binmode;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Construct the bins   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Before we used the pressure levels defined in the OMNO2 readme. When we
% updated BEHR to use a variable tropopause, we needed to extend the levels
% above our previous 200 hPa end point, so at that time, it made sense to
% switch to using then same pressure levels as BEHR. The difference should
% be minor, with most levels still present.
bin_midpoints = behr_pres_levels();

% Bucsela 2008 states that "bin widths [are] equal to the level spacing;"
% so I will define bin edges as the center +/- 0.5*D, where D is the larger
% of the spacings above or below.

delta = diff(bin_midpoints);
bins = zeros(numel(bin_midpoints,3));
for a=1:numel(bin_midpoints)
    % The first and last bins need special handling, since diff produces a
    % vector one shorter than its input
    if a==1;
        D = delta(1);
    elseif a==numel(bin_midpoints)
        D = delta(end);
    else
        D = max(delta((a-1):a));
    end
    % bins will have three columns: the bottom of the bin, the bin
    % midpoint, and the top of the top.
    bins(a,1) = bin_midpoints(a) - D/2;
    bins(a,2) = bin_midpoints(a);
    bins(a,3) = bin_midpoints(a) + D/2;
end

% Do the binning. 
if strcmpi(binmode, 'binonly')
    bin_values = cell(1,size(bins,1));
else
    bin_values = zeros(1,size(bins,1)); % output as row
end

if strcmpi(binmode,'mean') % For means, we output the std. dev.
    bin_errors = zeros(1,size(bins,1));
elseif strcmpi(binmode,'median') % For medians, we use the 25th and 75th quartiles
    bin_errors = zeros(2,size(bins,1));
else % For just binning, there is no need for bin_errors
    bin_errors = [];
end

for a=1:numel(bin_values)
    bin_data_vals = data_vals(pressure > bins(a,3) & pressure <= bins(a,1));
    if strcmpi(binmode,'mean')
        bin_values(a) = nanmean(bin_data_vals(:));
        bin_errors(a) = nanstd(bin_data_vals(:)) / sqrt(numel(bin_data_vals));
    elseif strcmpi(binmode, 'median')
        bin_values(a) = nanmedian(bin_data_vals(:));
        bin_errors(:,a) = quantile(bin_data_vals(:),[0.25,0.75]);
    else
        bin_values{a} = bin_data_vals(:);
    end
end 

% Output the bin edges
bin_edges = bins(:,[1,3]);

end

