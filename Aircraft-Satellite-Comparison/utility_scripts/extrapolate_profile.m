function [ data_bins, pres_bins, flag ] = extrapolate_profile( data_in, pres_in, varargin )
%extrapolate_profile Extrapolate an aircraft profile to fill in all pressure bins
%   This function will bin aircraft data by omi_sp_bins and also
%   automatically extrapolate the profile to fill in any empty bins unless
%   specifically told not to. It assumes a default surface pressure of 1020
%   hPa, this can be overridded.
%
%   Parameters:
%       shape - can be 'exp' or 'linear' (linear is default).  If set to
%           exp, interpolation will be carried out on ln(data) and ln(pres)
%           rather than their normal values.  This should be used when data
%           has an exponential variation with pressure/altitude.
%
%       surfacePressure - the surface pressure in hPa. Defaults to 1020,
%           i.e. all OMI SP model pressure bins are reported.
%
%       bottom - can be 'median' (default), 'fit', 'ground', or 'none'.
%
%       top - Determines how the top of the profile is extrapolated. Can be
%           'median' (default), 'fit', 'wrf', 'wrf-scaled', or 'none'.  Median
%           simply takes the top npoints data measurements and sets the top
%           bin equal to their median. Fit takes a best fit of the npoints
%           top and bottom non-nan measurments and using that fit to
%           calculate the concentration at each empty bin center. This uses
%           a regression that only minimized residuals in data, not
%           pressure.  'wrf' and 'wrf-scaled' both substitude the 12-km
%           monthly WRF profile nearest to the profile, wrf will just
%           directly append the shape, while wrf-scaled will multiply the
%           profile shape by a constant such that the WRF shape factor
%           corresponding to the last bin of the in-situ profile equals the
%           shape factor of the in-situ profile for that bin.
%
%       npoints - how many points to use for either the median or fit
%           extrapolation approaches.
%
%       city - a string with the city name.  Needed to load the ground site
%           super merge file for the 'ground' bottom extrapolation option.
%
%       date - the date string or date number that the profile is from;
%           needed for the ground extrapolation option for the bottom
%           profile.
%
%       utc - the starting or average utc time for the profile. Needed for
%           ground extrapolation.
%
%       sitenum - the site number for the profile.  Needed for the ground
%           extrapolation method for the bottom of the profile.
%
%       month - the numeric month of the data, needed for either of the WRF
%           extrapolation options. Can be passed as a scalar or a string,
%           but must be e.g. 7 or '7', not 'Jul' or 'July'.
%
%       conversion - A conversion of the in-situ data from its native units
%           to a straight mixing ratio.  Defaults to 1e-12, e.g. ppt to
%           part-per-part.  Needed only for the WRF options.
%
%       lat - the latitude of the profile. Can be a vector of latitude
%           values (with fills removed) or already averaged. Required for
%           either WRF extrapolation method.
%
%       lon - the longitude of the profile, same formats as lat allowed.
%           Also required for either WRF extrapolation method.
%
%   Dependencies: bin_omisp_pressure, findbysize
%
%   Josh Laughner <joshlaugh5@gmail.com> Aug 2014

E = JLLErrors;

p = inputParser;
p.addParameter('shape','exp',@(x) any(strcmpi(x,{'exp','linear'})));
p.addParameter('surfacePressure',1020,@isscalar);
p.addParameter('bottom','median',@(x) any(strcmp(x,{'median','fit','ground','none'})));
p.addParameter('top','median',@(x) any(strcmp(x,{'median','fit','wrf','wrf-scaled','none'})));
p.addParameter('npoints',10,@isscalar);
p.addParameter('city','',@isstr);
p.addParameter('date','',@(x) ischar(x) || isscalar(x));
p.addParameter('sitenum',[],@isscalar);
p.addParameter('utc',[],@isscalar);
p.addParameter('ground_convert',[],@(x) (isempty(x) || isscalar(x)));
p.addParameter('month','',@(x) (ischar(x) || isscalar(x)));
p.addParameter('conversion',1e-12,@isscalar);
p.addParameter('lat',[],@isnumeric);
p.addParameter('lon',[],@isnumeric);
p.parse(varargin{:});
pout = p.Results;

npoints = pout.npoints;
prof_shape = pout.shape;
surfpres = pout.surfacePressure;
bottom = pout.bottom;
top = pout.top;
city = pout.city;
date = pout.date;
sitenum = pout.sitenum;
prof_utc = pout.utc;
ground_convert = pout.ground_convert;
month = pout.month; 
conversion = pout.conversion;
lat = pout.lat;
lon = pout.lon;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% The path where the super merge files for the ground stations are stored.
% Will need to be changed if they move.
groundpath = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Ground/SuperMerges/';

% Check that any additional parameters needed for the bottom or top
% extrapolation methods have been passed.
if any(strcmp(top,{'wrf','wrf-scaled'})) 
    if isempty(month)
        error('extrp_profile:month_input','WRF profile extrapolation selected. A month must be passed');
    end
    if isempty(lat) || isempty(lon)
        error('extrap_profile:lat_lon_input','WRF profile extrapolation selected; latitude and longitude of the profile must be passed.');
    end
end
if strcmp(bottom,'ground')
    if isempty(city)
        error('extrap_profile:city_input','Ground site extrapolation selected, city name required.');
    end
    if isempty(date)
        error('extrap_profile:city_input','Ground site extrapolation selected, date required.');
    end
    if isempty(sitenum)
        error('extrap_profile:city_input','Ground site extrapolation selected, site number required.');
    end
    if isempty(prof_utc)
        error('extrap_profile:utc_input','Ground site extrapolation selected, profile starting or average UTC required.');
    end
end

% Convert month to a number, and check that it is an integer between 1 and 12
if any(strcmp(top,{'wrf','wrf-scaled'}))
    if ischar(month); month = str2double(month); end
    if mod(month,1)~=0 || month < 1 || month > 12;
        error('extrap_profile:month_input','''Month'' must be a whole number between 1 and 12');
    end
    % Use the month to generate the correct file name
    wrf_profile_name = sprintf('m%02d_NO2_profile.mat',month);
    
    % Store the WRF-Profile file path.  The directory will need to be
    % recoded if it moves.
    wrf_profile_file = sprintf('/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles/%s',wrf_profile_name);
end

% Flag that will be set to 1 if any problems arise.
flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[data_bins, pres_bins] = bin_omisp_pressure(pres_in, data_in);
above_surface = pres_bins <= surfpres;
data_bins = data_bins(above_surface); pres_bins = pres_bins(above_surface);

% Convert everything to logs if the profile is exponential
if strcmpi(prof_shape, 'exp');
    data_bins = log(data_bins); pres_bins = log(pres_bins);
    data_in = log(data_in); pres_in = log(pres_in);
end

nans = isnan(data_in) | isnan(pres_in);
data_in(nans) = []; pres_in(nans) = [];

if numel(data_in) == 0;
    flag = 2;
    return
end

% Check that the top bin is a nan, otherwise we don't need to do any
% extrapolation
if ~isnan(data_bins(end));
    fprintf('!!! Top bin already has a value. No extrapolation carried out. !!!\n')
else
    switch top
        case 'none'
            %Do nothing
        case 'median'
            % This case doesn't really work when there's no existing data
            % points (i.e. all are below the surface) so we need to error
            % out if thats the case
            if all(isnan(data_bins))
                E.nodata('data_bins');
            end
            
            % Find the top npoints data point
            top_points = findbysize(pres_in,npoints,'smallest');
            
            % Set the last bin to equal the median of these points.  The
            % intervening bins will be filled in when any internal nans are
            % filled by fill_nans below
            data_bins(end) = nanmedian(data_in(top_points));
            
        case 'fit'
            % Find the top npoints data point
            top_points = findbysize(pres_in,npoints,'smallest');
            
            % Remove any points > 2 std. dev. from the mean
            pres_top = pres_in(top_points); data_top = data_in(top_points);
            top_std = nanstd(data_top); top_mean = nanmean(data_top);
            two_sd = data_top > top_mean - 2*top_std & data_top < top_mean + 2*top_std;
            fit_top = polyfit(pres_top(two_sd), data_top(two_sd), 1);
            
            % Extrapolate the top data points
            lastbin = find(~isnan(data_bins),1,'last');
            if isempty(lastbin)
                lastbin = 0;
            end
            data_bins((lastbin+1):end) = fit_top(1) .* pres_bins((lastbin+1):end) + fit_top(2);
            
        case {'wrf','wrf-scaled'}
            %Both wrf cases are very similar, except that the "scaled" version
            %needs to match to the top valid bin of the in-situ profile
            
            load(wrf_profile_file); % Loads the PROFILE variable
            
            % Find the closest WRF profile (geographically) to this one
            lat = nanmean(lat); lon = nanmean(lon);
            delta = abs(PROFILE.Latitude - lat) + abs(PROFILE.Longitude - lon);
            [~, ind] = min(delta(:)); [x,y] = ind2sub(size(PROFILE.Latitude), ind);
            
            % Read the WRF profile and scale it to the units of the in-situ profile
            wrf_prof = (1e-6 / conversion) .* PROFILE.NO2_profile(:,x,y);
            
            % If we are doing an exponential fit, then convert the WRF
            % pressures and profile mixing ratios to natural logs
            if strcmpi(prof_shape,'exp')
                wrf_prof = log(wrf_prof);
                wrf_pres = log(PROFILE.Pressure);
            else
                wrf_pres = PROFILE.Pressure;
            end
            
            % Intepolate the WRF profile to the pressure of the in-situ bins
            wrf_prof_interp = interp1(wrf_pres, wrf_prof, pres_bins);
            
            % Find the last non-nan bin of the in-situ profile.  If we want
            % the scaled wrf profile, scale the corresponding bin to match
            % and multiply the rest of the profile by that factor.  In
            % either case, replace all following in-situ bins with their
            % (scaled or unscaled) WRF counterparts.
            last_bin = find(~isnan(data_bins),1,'last');
            if isempty(last_bin)
                last_bin = 0;
            end
            
            if strcmp(top,'wrf-scaled') && last_bin > 0;
                scale_factor = data_bins(last_bin) / wrf_prof_interp(last_bin);
                wrf_prof_interp = wrf_prof_interp .* scale_factor;
            end
            
            data_bins((last_bin+1):end) = wrf_prof_interp((last_bin+1):end);
    end
end

if ~isnan(data_bins(1))
    fprintf('!!! Bottom bin has a value. No extrapolation carried out. !!!\n');
else
    switch bottom
        case 'none'
            % Do nothing
        case 'median'
            % Find the top npoints data point
            bottom_points = findbysize(pres_in,npoints,'largest');
            
            % Set the first bin to equal the median of these points.  The
            % intervening bins will be filled in when any internal nans are
            % filled by fill_nans below
%             if strcmpi(prof_shape,'exp')
%                 data_bins(1) = nanmedian(exp(data_in(bottom_points)));
%             else
                data_bins(1) = nanmedian(data_in(bottom_points));
%             end
        case 'fit'
            % Find the bottom npoint data points
            bottom_points = findbysize(pres_in, npoints, 'largest');
            
            % Remove any points > 2 std. dev. from the mean
            pres_bottom = pres_in(bottom_points); data_bottom = data_in(bottom_points);
            bottom_std = nanstd(data_bottom); bottom_mean = nanmean(data_bottom);
            two_sd = data_bottom > bottom_mean - 2*bottom_std & data_bottom < bottom_mean + 2*top_std;
            
            % Find the fit
            fit_bottom = polyfit(pres_bottom(two_sd), data_bottom(two_sd), 1);
            
            % Find any bins on the bottom that have nans and extrapolate their value
            % using the bottom fit
            firstbin = find(~isnan(data_bins),1,'first');
            firstnans = 1:(firstbin-1);
            data_bins(firstnans) = fit_bottom(1) .* pres_bins(firstnans) + fit_bottom(2);
        case 'ground'
            supermerge_files = dir(fullfile(groundpath,sprintf('%s*.mat',city)));
            if numel(supermerge_files) == 1; 
                load(fullfile(groundpath,supermerge_files(1).name)); % Loads the SuperMerge variable
            else
                error('extrap_profile:load_supermerge','Too many or no supermerge files found for the given city');
            end
            
            % Format the site number; Baltimore site numbers are a single
            % digit while CA and Houston numbers are six digits, where the
            % 2nd and 3rd actually represent the site number.
            if sitenum < 10; % then it is a Baltimore type site number
                sitenum_str = sprintf('%02d',sitenum);
            else % otherwise assume that it's a CA or Houston number
                sitenum_str = num2str(sitenum);
                sitenum_str = sitenum_str(2:3);
            end
            
            % Create the field names to access the super merge structure;
            % it should be arranged SuperMerge --> Date --> Site, where the
            % date fields have the format Date_yyyymmdd and the site fields
            % are Site_nn
            datefield = sprintf('Date_%s',datestr(date,'yyyymmdd'));
            sitefield = sprintf('Site_%s',sitenum_str);
            
            % Check if the require fields exist. If the date field doesn't
            % exist, it is likely that something incorrect was passed, so
            % throw an error.  If the site field doesn't exist, it may be
            % one of the sites that 
            if ~isfield(SuperMerge, datefield)
                error('extrap_profile:datefield','The field ''%s'' cannot be found in the SuperMerge structure.',datefield);
            elseif ~isfield(SuperMerge.(datefield),sitefield)
                warning('The field ''%s'' cannot be found within SuperMerge.%s. Site %s may not have data for NO2. Defaulting to median extrapolation.',sitefield,datefield,sitenum_str);
                flag = 1;
                bottom_points = findbysize(pres_in,npoints,'largest');
%                 if strcmpi(prof_shape,'exp')
%                     bottom_bin = nanmedian(exp(data_in(bottom_points)));
%                 else
                    bottom_bin = nanmedian(data_in(bottom_points));
%                 end
            else
                unit = SuperMerge.(datefield).(sitefield).NO2.Unit;
                ground_no2 = SuperMerge.(datefield).(sitefield).NO2.Values;
                ground_start_utc = SuperMerge.(datefield).(sitefield).UTC.Values;
                xx = find(ground_start_utc <= prof_utc,1,'last');
                
                % Figure out the unit if it is ppm, ppb or ppt.
                if isempty(ground_convert)
                    ind = regexpi(unit,'pp[m,b,t]');
                    u = lower(unit(ind+2)); % convert the character that indicates the part per million/billion/trillion to lower case for the switch statement
                    switch u
                        case 'm'
                            ground_convert = 1e-6;
                        case 'b'
                            ground_convert = 1e-9;
                        case 't'
                            ground_convert = 1e-12;
                        otherwise
                            error('extrap_profile:ground_unit','Unit for ground site (%s) not recognized.  Pass conversion factor as parameter ''ground_convert''',unit);
                    end
                end
                
                % Find the period of measurement that contains the profile
                % time, then convert it to the same units as our aircraft
                % profile and make it the bottom bin.
                bottom_bin = ground_no2(xx) * ground_convert / conversion;
                if bottom_bin < 0;
                    warning('The ground measurement for NO2 is < 0. Assuming fill value and defaulting to median extrapolation.');
                    flag = 1;
                    bottom_points = findbysize(pres_in,npoints,'largest');
                    if strcmpi(prof_shape,'exp')
                        bottom_bin = nanmedian(exp(data_in(bottom_points)));
                    else
                        bottom_bin = nanmedian(data_in(bottom_points));
                    end
                end
            end
            if strcmpi(prof_shape,'exp'); bottom_bin = log(bottom_bin); end
            data_bins(1) = bottom_bin;
    end
end

% Fill in internal nans
[pres_bins, data_bins] = fill_nans(pres_bins, data_bins,'noclip');

% If assuming an exponential profile, exponentiate the bins to undo the log
if strcmpi(prof_shape, 'exp');
    data_bins = exp(data_bins); pres_bins = exp(pres_bins);
end


end

