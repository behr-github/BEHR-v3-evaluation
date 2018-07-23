function [ bin_avgs, bin_midpoints ] = plot_vertical_profile( Merge_in, data_field, varargin )
%plot_vertical_profile(Merge, data_field): Plots a vertical profile of the data in "data_field" against pressure altitude.
%   This function accepts a Merge data structure prepared by
%   "read_merge_files" and makes a vertical profile of the data specified
%   in data_field.  This function is built to work for aircraft merge files.
%   Parameter arguments:
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
%       'binwidth' = width of the altitude bins in km.  If not specified,
%       all points will be plotted with no binning.
%
%       'binmode' = sets if the binned value will be a 'median' with 25th and
%       75th quartile or 'mean' with std. error. Median is default.
%
%   Josh Laughner <joshlaugh5@gmail.com> 23 May 2014


p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('data_field',@isstr);
p.addParamValue('siteflag',[-1e3,1e3],@isnumeric);
p.addParamValue('profnum',[-1e10,1e10],@isnumeric);
p.addParamValue('binwidth',-1,@isscalar);
p.addParamValue('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));
p.addParamValue('opmode','std',@(x) any(strcmpi(x,{'std','std.','standard','rolling'})));

p.parse(Merge_in,data_field,varargin{:});
pout = p.Results;

Merge = pout.Merge;
field = pout.data_field;
siteflag = [min(pout.siteflag), max(pout.siteflag)];
profnum = [min(pout.profnum), max(pout.profnum)];
binwidth = pout.binwidth;
binmode = pout.binmode;
opmode = pout.opmode;

% Find the entries that have the desired site flag and profile number, if
% we are using discover-aq data.  Otherwise, just import the data.
clear profnum_logical site_logical %Clear these just in case they are hanging around from a previous function call
if isfield(Merge.Data,'ProfileSequenceNum');
    profnum_logical = (Merge.Data.ProfileSequenceNum.Values >= profnum(1) & Merge.Data.ProfileSequenceNum.Values <= profnum(2));
elseif isfield(Merge.Data,'ProfileNumber');
    profnum_logical = (Merge.Data.ProfileNumber.Values >= profnum(1) & Merge.Data.ProfileNumber.Values <= profnum(2));
end

if isfield(Merge.Data,'discoveraqSiteFlag1sec')
    site_logical = (Merge.Data.discoveraqSiteFlag1sec.Values >= siteflag(1) & Merge.Data.discoveraqSiteFlag1sec.Values <= siteflag(2));
elseif isfield(Merge.Data, 'SiteSeqNumber')
    site_logical = (Merge.Data.SiteSeqNumber.Values >= siteflag(1) & Merge.Data.SiteSeqNumber.Values <= siteflag(2));
end

if exist('profnum_logical','var') && exist('site_logical','var')
    pressures = Merge.Data.PRESSURE.Values(site_logical & profnum_logical);
    data_vals = eval(sprintf('Merge.Data.%s.Values(site_logical & profnum_logical)',field));
    local_times = Merge.Data.LOCAL_SUN_TIME.Values(site_logical & profnum_logical);
else
    pressures = Merge.Data.PRESSURE.Values;
    data_vals = eval(sprintf('Merge.Data.%s.Values',field));
    local_times = Merge.Data.LOCAL_SUN_TIME.Values;
end

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

% If binwidth is negative (i.e. not specified by the user), plot every data
% point.  If a bin width is specified, average together all values within
% each bin slice, then plot that average at the midpoint altitude.
if binwidth < 0;
    scatter(data_vals,altitude);
else
    if strcmpi(opmode,'rolling')
        [bin_avgs, bin_midpoints, bin_error] = bin_rolling_vertical_profile(altitude,data_vals,binwidth,0.1,binmode);
        xx = find(isnan(bin_avgs)); %Removes nans to ensure that the error envelope plots. 
        bin_avgs(xx) = []; bin_midpoints(xx) = []; bin_error(:,xx) = [];
    else
        [bin_avgs, bin_midpoints, bin_error] = bin_vertical_profile(altitude,data_vals,binwidth,binmode);
    end
    
    % Make a line plot of the vertical profile with the standard error
    % plotted as an envelope
    
    if strcmpi(binmode,'mean')
        plot_error_envelope_x(bin_midpoints, bin_avgs - bin_error, bin_avgs + bin_error);
    else
        plot_error_envelope_x(bin_midpoints, bin_error(1,:), bin_error(2,:));
    end
    hold on
    plot(bin_avgs, bin_midpoints);
    
    printfield = regexprep(field,'_',' ');
    starttime = min(local_times);
    starttime_str = sprintf('%d:%02d',floor(starttime),floor(mod(starttime,1)*60));
    endtime = max(local_times);
    endtime_str = sprintf('%d:%02d',floor(endtime),floor(mod(endtime,1)*60));
    titlestr = sprintf('%s %s vertical profile for %s\n LT %s-%s\n',binmode,printfield,Merge.metadata.date,starttime_str,endtime_str);
    if ~any(profnum < -100); titlestr = [titlestr, sprintf('Profile num: %s  ',mat2str(profnum))];
    else titlestr = [titlestr, 'Profile number: all  '];
    end
    
    if ~any(siteflag < -100); titlestr = [titlestr, sprintf('Site flags: %s  ', mat2str(siteflag))];
    else titlestr = [titlestr, 'Site flags: all  '];
    end
    title(titlestr,'fontsize',18)
    
    xlabelstr = sprintf('%s (%s)',printfield, eval(sprintf('Merge.Data.%s.Unit  ',field)));
    if binwidth > 0; xlabelstr = [xlabelstr, sprintf('  %2.2f km bins ',binwidth)]; end
    xlabel(xlabelstr,'fontsize',14);
    ylabel('Altitude (km)','fontsize',14);
    ylim([(min(altitude)), ceil(max(altitude))]);
    hold off
end

end

