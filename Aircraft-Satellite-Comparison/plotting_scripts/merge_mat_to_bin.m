function [ bin_values, bin_midpoints, bin_error ] = merge_mat_to_bin( Merge_in, data_field, binwidth, varargin )
%[bin_values, bin_midpoints, bin_error] = merge_mat_to_bin(Merge, data_field, binwidth): An interactive function that passes a merge data structure to the bin_vertical_profile function.
%   This function takes a merge data structure from "read_merge_file.m," a
%   data field name as a string, a binwidth, and optionally the binning
%   mode ('median' = median + quartiles (default); 'mean' = mean + std.
%   dev.)  This is intended primarily as an interactive function to allow
%   quick binning of whole merge files.  For better control over what is
%   binned, use "bin_vertical_profile.m" in a custom script to pass
%   specific data.  Parameter values:
%
%       'siteflag' = A single site flag value or a two-element matrix
%       specifying a range of site flag values.  If used, only data points
%       matching that site flag will be used in the profiles. If both this
%       and the 'profnum' flag are specified, only data points that match
%       both flags will be used.
%
%       'profnum' = A single Profile Sequence Number or two-element matrix
%       specifying a range of these numbers.  Works identically to
%       'siteflag' except using the profile number instead of the site
%       flag.
%       
%       'binmode' = sets if the binned value will be a 'median' with 25th and
%       75th quartile or 'mean' with std. error. Median is default.
%
%   Josh Laughner <joshlaugh5@gmail.com> 27 May 2014

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('data_field',@isstr);
p.addRequired('binwidth',@isscalar);
p.addParamValue('siteflag',[-1e3,1e3],@isnumeric);
p.addParamValue('profnum',[-1e10,1e10],@isnumeric);
p.addParamValue('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(Merge_in,data_field,binwidth,varargin{:});
pout = p.Results;

Merge = pout.Merge;
field = pout.data_field;
siteflag = [min(pout.siteflag), max(pout.siteflag)];
profnum = [min(pout.profnum), max(pout.profnum)];
binwidth = pout.binwidth;
binmode = pout.binmode;

% Find the entries that have the desired site flag and profile number
site_logical = (Merge.Data.discoveraqSiteFlag1sec.Values >= siteflag(1) & Merge.Data.discoveraqSiteFlag1sec.Values <= siteflag(2));
profnum_logical = (Merge.Data.ProfileSequenceNum.Values >= profnum(1) & Merge.Data.ProfileSequenceNum.Values <= profnum(2));

pressures = Merge.Data.PRESSURE.Values(site_logical & profnum_logical);
data_vals = eval(sprintf('Merge.Data.%s.Values(site_logical & profnum_logical)',field));

% Replace any fill values, upper LOD, or lower LOD values with NaNs
fill_val = eval(sprintf('Merge.Data.%s.Fill',field));
ULOD_val = Merge.metadata.upper_lod_flag;
LLOD_val = Merge.metadata.lower_lod_flag;

data_vals(data_vals == fill_val) = NaN;
data_vals(data_vals == ULOD_val) = NaN;
data_vals(data_vals == LLOD_val) = NaN;

% Plot the values converting pressure to altitude, using P0 = 1013.25 hPa
% and scale height H = 7.4 km

altitude = -log(pressures ./ 1013.25) .* 7.4;

% Carry out the binning
[bin_values, bin_midpoints, bin_error] = bin_vertical_profile(altitude,data_vals,binwidth,binmode);
end

