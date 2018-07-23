function [ bin_avgs, bin_midpoints ] = plot_sonde_vertical_profile( Merge_in, data_field, varargin )
%plot_vertical_profile(Merge, data_field): Plots a vertical profile of the data in "data_field" against pressure altitude.
%   This function accepts a Merge data structure prepared by
%   "read_merge_files" and makes a vertical profile of the data specified
%   in data_field.  This function is built to work for aircraft merge files.
%   Parameter arguments:
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
p.addParamValue('binwidth',-1,@isscalar);
p.addParamValue('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(Merge_in,data_field,varargin{:});
pout = p.Results;

Merge = pout.Merge;
field = pout.data_field;
binwidth = pout.binwidth;
binmode = pout.binmode;

pressures = Merge.Data.Pressure_hPa.Values;
data_vals = eval(sprintf('Merge.Data.%s.Values',field));
utc_times = Merge.Data.Start_UTC.Values;

% Replace any fill values, upper LOD, or lower LOD values with NaNs
data_fill_val = eval(sprintf('Merge.Data.%s.Fill',field));
pres_fill_val = Merge.Data.Pressure_hPa.Fill;
ULOD_val = Merge.metadata.upper_lod_flag;
LLOD_val = Merge.metadata.lower_lod_flag;

data_vals(data_vals == data_fill_val) = NaN;
data_vals(data_vals == ULOD_val) = NaN;
data_vals(data_vals == LLOD_val) = NaN;

xx = find(pressures == pres_fill_val);
data_vals(xx) = []; pressures(xx) = []; utc_times(xx) = [];

% Plot the values converting pressure to altitude, using P0 = 1013.25 hPa
% and scale height H = 7.4 km

altitude = -log(pressures ./ 1013.25) .* 7.4;

% If binwidth is negative (i.e. not specified by the user), plot every data
% point.  If a bin width is specified, average together all values within
% each bin slice, then plot that average at the midpoint altitude.
if binwidth < 0;
    scatter(data_vals,altitude);
else
    [bin_avgs, bin_midpoints, bin_error] = bin_vertical_profile(altitude,data_vals,binwidth,binmode);
    
    % Make a line plot of the vertical profile with the standard error
    % plotted as an envelope
    
    if strcmpi(binmode,'mean')
        plot_error_envelope_x(bin_midpoints, bin_avgs - bin_error, bin_avgs + bin_error);
    else
        plot_error_envelope_x(bin_midpoints, bin_error(:,1)', bin_error(:,2)');
    end
    hold on
    plot(bin_avgs, bin_midpoints);
    
    printfield = regexprep(field,'_',' ');
    starttime = min(utc_times);
    starttime_str = sprintf('%d:%02d',floor(starttime/3600),floor(mod(starttime,3600)/60));
    endtime = max(utc_times);
    endtime_str = sprintf('%d:%02d',floor(endtime/3600),floor(mod(endtime,3600)/60));
    titlestr = sprintf('%s %s vertical profile for %s\n UTC %s-%s',binmode,printfield,Merge.metadata.date,starttime_str,endtime_str);
    
    title(titlestr,'fontsize',18)
    
    xlabelstr = sprintf('%s (%s)',printfield, eval(sprintf('Merge.Data.%s.Unit  ',field)));
    if binwidth > 0; xlabelstr = [xlabelstr, sprintf('  %2.2f km bins ',binwidth)]; end
    xlabel(xlabelstr,'fontsize',14);
    ylabel('Altitude (km)','fontsize',14);
    ylim([(min(altitude)), ceil(max(altitude))]);
    hold off
end

end

