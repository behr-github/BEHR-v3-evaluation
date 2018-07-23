function [omi_lon, omi_lat, OMI_NO2, BEHR_NO2, DC8_NO2, varianceNO2, corner_lon, corner_lat, choices, concNO2debug, BLHdebug] = boundary_layer_verification( Merge, Data, tz, varargin )
%[ pix_lon, pix_lat, OMI_NO2, BEHR_NO2, aircraft_NO2]
%
%boundary_layer_verification: Compare satellite and aircraft measurements
%using the boundary layer method. Returns pixel latitude, longitude, BEHR
%tropospheric NO2, aircraft trop. NO2, and OMI trop. NO2.
%   This function compares satellite and aircraft NO2 columns using the
%   boundary layer method described in Russell et. al. Atmos. Chem. Phys.
%   11, 8543-8554, 2011. As inputs, it requires a Merge data structure (the
%   result of reading an ICART file using read_icart_file.m), a BEHR Data
%   data structure (one of the outputs of BEHR_main.m), and a time zone
%   abbreviation:
%
%       EST = Eastern Std.
%       CST = Central Std.
%       MST = Mountain Std.
%       PST = Pacific Std.
%
%   Since Aura satellite overpass is ~1:45 p in standard time, do not use
%   daylight savings time values.  Further, you must pass only one element
%   of Data; the one corresponding to the swath of interest.
%
% Parameters
%   timerange: By default, this method will restrict flight data to between
%   12:00p and 3:00p local, as per Russell et. al.  This can be overridden
%   using this parameter, which requires a cell array with the start and
%   end times as strings in military time format.
%
%   cloud_product: Which cloud product (omi or modis) to use in rejecting
%   pixels.  Defaults to omi.
%
%   cloud_frac_max: The maximum allowed geometric cloud fraction.  Defaults
%   to 0.2; recommended value for use with MODIS cloud product is 0.
%
%   rowanomaly: The method of rejecting pixels based on row anomaly.
%   Defaults to 'AlwaysByRow'.  See omi_pixel_reject or omi_rowanomaly for
%   more information on the possible choices ('AlwaysByRow', 'RowsByTime',
%   'XTrackFlags', and 'XTrackFlagsLight').
%
%   DEBUG_LEVEL: The level of output messages to write; 0 = none, 1 =
%   normal; 2 = verbose
%
%   Josh Laughner <joshlaugh5@gmail.com> June 2014

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst',})));
p.addOptional('choices',struct('no2',[],'h2o',[],'theta',[]),@isstruct)
p.addParamValue('timerange',{'12:00','15:00'},@iscell)
p.addParamValue('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis'})));
p.addParamValue('cloud_frac_max',0.2, @isscalar);
p.addParamValue('rowanomaly','AlwaysByRow',@(x) strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}));
p.addParamValue('DEBUG_LEVEL',1,@isscalar);

p.parse(Merge,Data,tz,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
choices_in = pout.choices;
timerange = pout.timerange;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

% Load the aircraft data
if DEBUG_LEVEL > 0; fprintf(' Getting data and restricting to 12:00-15:00 %s\n',tz); end
no2_raw = Merge.Data.NO2_UCB.Values;
alt_raw = Merge.Data.ALTP.Values;
radar_alt_raw = Merge.Data.Radar_Altitude.Values;
utc_raw = Merge.Data.UTC.Values;
lat_raw = Merge.Data.LATITUDE.Values;
lon_raw = Merge.Data.LONGITUDE.Values - 360;
pressure_raw = Merge.Data.PRESSURE.Values;
temperature_raw = Merge.Data.TEMPERATURE.Values;
theta_raw = Merge.Data.THETA.Values;
h2o_raw = Merge.Data.H2Ov.Values;

% Get the date
merge_date = Merge.metadata.date;

% Load the file containing the altitude ranges in which to look for
% boundary layer crossings. Loads the variable "Ranges," a structure.

load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/ARCTAS-CA Altitude Ranges Exclusive 2.mat');
% Store the date corresponding to each index for easy searching later
range_dates = {datestr(Ranges(1).Date,29), datestr(Ranges(2).Date,29), datestr(Ranges(3).Date,29), datestr(Ranges(4).Date,29)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%       LOAD IN DATA       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data for NO2, water, and potential temperature, replacing fill
% values with nans.
[no2, utc, alt, lon, lat] = remove_merge_fills(Merge,'NO2_UCB');
[h2o, utc_h2o, alt_h2o] = remove_merge_fills(Merge,'H2Ov');
[theta, utc_theta, alt_theta] = remove_merge_fills(Merge,'THETA');


% Now for temperature - actually remove these values
[temperature, utc_T, alt_T, ~, ~, xx_T] = remove_merge_fills(Merge,'TEMPERATURE');
alt_T = alt_T(~xx_T);
utc_T = utc_T(~xx_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  FIND BDY LAYER HEIGHTS  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DEBUG_LEVEL > 0; fprintf(' Calculating & interpolating boundary layer heights.\n'); end

% Get the range for today;
d = find(datenum(range_dates) == datenum(merge_date));

[heights_no2, times, choices.no2] = findall_no2_bdy_layer_heights(utc, no2, alt, Ranges(d).Ranges,choices_in.no2);
[heights_h2o, ~, choices.h2o] = findall_no2_bdy_layer_heights(utc_h2o, h2o, alt_h2o, Ranges(d).Ranges,choices_in.h2o);
[heights_theta, ~, choices.theta] = findall_no2_bdy_layer_heights(utc_theta, theta, alt_theta, Ranges(d).Ranges,choices_in.theta,'blmode','theta');

tmp = [heights_no2, heights_h2o, heights_theta];
heights = nanmean(tmp,2);
nans = isnan(heights);
heights = heights(~nans); times = times(~nans);
if numel(heights) < 2
    error('bdy_layer_verify:findBoundaryLayer','Could not find 2 boundary layer heights - needed for interpolation')
end

interp_height = interp1(times, heights, utc_raw); % Linearly interpolate the boundary layer height to every value of UTC.

% For values outside of the range of "times," assume that the boundary
% layer height equals the closest value.
first_height = find(~isnan(interp_height),1,'first');
last_height = find(~isnan(interp_height),1,'last');

if first_height > 1;
    interp_height(1:first_height-1) = interp_height(first_height);
end
if last_height < numel(interp_height)
    interp_height(last_height+1:end) = interp_height(last_height);
end

% Remove any values outside the time range (12:00-3:00p local standard by
% default)
utcstart = local2utc(timerange{1}, tz); utcend = local2utc(timerange{2}, tz);
time_logical = utc_raw >= utcstart & utc_raw <= utcend;
no2 = no2(time_logical);
alt = alt(time_logical);
utc = utc(time_logical);
lat = lat(time_logical);
lon = lon(time_logical);

% Also clip the interpolated height to the time range 
interp_height = interp_height(time_logical);
%figure; plot(utc,alt); line(utc,interp_height,'color','r'); title(merge_date);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PREP TEMP AND PRESSURE DATA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract a temperature profile; we will need this to convert mixing
% ratios to number densities
[T_bins, T_bin_alt] = bin_vertical_profile(alt_T, temperature, 0.3);
% The temperature profiles often have an inversion in the BL, so split the
% temperature data into slope <0 (remnant PBL/free troposphere) and >=0
% (mixed layer)
ML_ind = find(diff(T_bins) < 0,1,'first')+1; % The +1 accounts for the fact that length(diff(x)) is length(x) - 1
[T_poly_mixed, T_R2_mixed] = polyfit_R2(T_bin_alt(1:ML_ind),T_bins(1:ML_ind),1);
[T_poly_free, T_R2_free] = polyfit_R2(T_bin_alt(ML_ind:end),T_bins(ML_ind:end),1);
% We'll calculate the mixed layer height as the intersection of the two
% fits
A = [-T_poly_free(1), 1; -T_poly_mixed(1), 1]; b = [T_poly_free(2); T_poly_mixed(2)];
X = A\b; ML_height = X(1)*1000;

% We need the P0 and scale height (H) used in the conversion of pressure to
% altitude.  Fit the data to a line such that ln(P) = (-1/H)z + ln(P0),
% thus H = -1/slope and P0 = e^(intercept)
P = polyfit(alt_raw,log(pressure_raw),1);
P0 = exp(P(2)); H = -1/P(1);


% Also, define Avogadro's number and gas constant;
Av = 6.022e23; % molec. / mol
R = 8.314e4; % (hPa * cm^3) / (mol * K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  COLOCATE AIRCRAFT SEGMENTS WITH PIXELS   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find pixels that fail the criteria used for BEHR mapping.
% omi_pixel_reject.m works by setting Areaweight to 0 for any pixels that
% fail the criteria. Hence, we need to create an areaweight field that
% assumes all pixels are valid.

Data.Areaweight = ones(size(Data.Longitude));
if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
Data2 = omi_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
xx = Data2.Areaweight > 0;

% Read in corner lat/lon and BEHR NO2 columns. Latcorn/loncorn are 4 x n
% fields, where the first dimension is each corner for a pixel.
omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
corner_lat = Data2.Latcorn(:,xx); corner_lon = Data2.Loncorn(:,xx);
BEHR_NO2 = Data2.BEHRColumnAmountNO2Trop(xx); GLOBETerpres = Data2.GLOBETerpres(xx);
OMI_NO2 = Data2.ColumnAmountNO2Trop(xx); TropopausePres = Data2.TropopausePressure(xx);
vza = Data2.ViewingZenithAngle(xx); 

% Remove from consideration all pixels whose corners fall entirely outside
% the range of lat/lons that the aircraft flew
if DEBUG_LEVEL > 0; fprintf('   Removing pixels outside aircraft track\n'); end
lon_logical = corner_lon > min(lon) & corner_lon < max(lon);
lat_logical = corner_lat > min(lat) & corner_lat < max(lat);
latlon_logical = any(lon_logical) & any(lat_logical); % any() operates on matrices along the first dimension; here that is each corner for pixel n

omi_lat = omi_lat(latlon_logical);
omi_lon = omi_lon(latlon_logical);
corner_lat = corner_lat(:,latlon_logical);
corner_lon = corner_lon(:,latlon_logical);
BEHR_NO2 = BEHR_NO2(latlon_logical);
OMI_NO2 = OMI_NO2(latlon_logical);
GLOBETerpres = GLOBETerpres(latlon_logical);
TropopausePres = TropopausePres(latlon_logical);
vza = vza(latlon_logical);

% Make the matrix to save the integrated aircraft columns to
DC8_NO2 = -9e9*ones(size(BEHR_NO2));
concNO2debug = -9e9*ones(size(BEHR_NO2));
BLHdebug = -9e9*ones(size(BEHR_NO2));
varianceNO2 = -9e9*ones(size(BEHR_NO2));

% Now reject any no2 measurments more than 2 sigma from the mean
% mean_no2 = nanmean(no2); std_no2 = nanstd(no2);
% no2lb = mean_no2 - 2*std_no2; no2ub = mean_no2 + 2*std_no2;
% in_std = no2 >= no2lb & no2 <= no2ub;
% no2(~in_std) = NaN;
% For each pixel left, find all aircraft measurments that fall within that
% pixel.
vza_rej = 0;
for a=1:numel(BEHR_NO2)
    if DEBUG_LEVEL > 1; fprintf('   Pixel %d of %d\n',a,numel(BEHR_NO2)); end
    IN = inpolygon(lon,lat,corner_lon(:,a),corner_lat(:,a));
    lon_pix = lon(IN); lat_pix = lat(IN);
    no2_pix = no2(IN); interp_height_pix = interp_height(IN);
    alt_pix = alt(IN);
    
    % Remove measurements above the boundary layer and that have an no2
    % measurement more than 2 std. dev. away from the mean
    in_BL = alt_pix <= interp_height_pix;
    
    mean_no2 = nanmean(no2_pix); std_no2 = nanstd(no2_pix);
    no2lb = mean_no2 - 2*std_no2; no2ub = mean_no2 + 2*std_no2;
    in_std = no2_pix >= no2lb & no2_pix <= no2ub;
    no2_pix(~in_std) = NaN;
    
    no2_pix = no2_pix(in_BL); interp_height_pix = interp_height_pix(in_BL);
    
    if sum(~isnan(no2_pix)) < 20 % If there are not 20 measurments within this pixel, skip it
        continue
    elseif vza(a) > 60; % Do not compare the pixel if the viewing zenith angle is >60 deg. These pixels get large.
        vza_rej = vza_rej+1;
        continue
    else % Otherwise, average the NO2 measurement and integrate from ground height to tropopause of the pixel
        mean_no2_pix = nanmean(no2_pix); concNO2debug(a) = mean_no2_pix;
        BLH = nanmedian(interp_height_pix)*1000; BLHdebug(a) = BLH;
        
        % Prepare variables from the BEHR Data Structure
        % Standard values for scale height and P0 used here because these
        % values are not derived from the aircraft data
        surface_alt = round(-log(GLOBETerpres(a)/1013.25)*7400);
        TP_pres = TropopausePres(a);
        TP_alt = round(-log(TP_pres/1013.25)*7400);
        
        % Integrate the column in dz meter increments
        no2_column = 0;
        dz = 1;
        
        if surface_alt < BLH
            % Surface to boundary layer
            for h=surface_alt:dz:BLH
                P = P0 * exp(-h/(H*1000));
                if h <= ML_height; % Account for the temperature inversion in the mixed layer
                    %T = polyval(T_poly_mixed,h/1000);
                    T = T_poly_mixed(1)*(h/1000) + T_poly_mixed(2);
                elseif h > ML_height;
                    %T = polyval(T_poly_free,h/1000);
                    T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                end
                conc_NO2 = (Av * P * mean_no2_pix*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
            % Boundary layer to tropopause
            for h=(BLH+1):dz:TP_alt
                P = P0 * exp(-h/(H*1000));
                T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                conc_NO2 = (Av * P * 40*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv.  Assume 40 pptv in the free troposphere
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
        else
            % Surface to tropopause
            for h=surface_alt:dz:TP_alt
                P = P0 * exp(-h/(H*1000));
                T = T_poly_free(1)*(h/1000) + T_poly_free(2);
                conc_NO2 = (Av * P * 40*1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv.  Assume 40 pptv in the free troposphere
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
        end
        DC8_NO2(a) = no2_column;
        varianceNO2(a) = (nanstd(no2_pix))/sqrt(numel(no2_pix));
    end
end
fprintf('Num. pixels rejected for VZA = %d\n',vza_rej);
fills = DC8_NO2 == -9e9;
omi_lon = omi_lon(~fills);
omi_lat = omi_lat(~fills);
corner_lat = corner_lat(:,~fills);
corner_lon = corner_lon(:,~fills);
DC8_NO2 = DC8_NO2(~fills);
OMI_NO2 = OMI_NO2(~fills);
BEHR_NO2 = BEHR_NO2(~fills);
concNO2debug = concNO2debug(~fills);
BLHdebug = BLHdebug(~fills);
varianceNO2 = varianceNO2(~fills);
end