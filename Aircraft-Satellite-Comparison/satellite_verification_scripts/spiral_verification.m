function [ omi_lon_out, omi_lat_out, omi_no2_out, behr_no2_out, air_no2_out, db ] = spiral_verification( Merge, Data, timezone, varargin )
%[lon, lat, omi, behr, air, cov_frac] = spiral_verification(Merge,Data,timezone) Compare OMI pixel NO2 values to aircraft spirals.
%
%   This function is based off of the method described in Hains et. al. (J.
%   Geophys. Res. 2010, 115, D05301 doi 10.1029/2009JD012399) and Bucesela
%   et. al. (J. Geophys. Res. 2008, 113, D16S31 doi:10.1029/2007/D008838)
%   to compare OMI NO2 columns to aircraft measurements using full spiral
%   profiles as the means to integrate NO2 concentrations.  The profile
%   data is binned by pressure bins as in the Bucsela article.
%   Extrapolation above and below the spiral is usually handled as in Hains
%   et. al.; the median lowest/highest 10 NO2 measurements are extrapolated
%   as a constant.  As in Hains, the surface altitude is determined as the
%   difference between the median of the lowest 10 GPS altitudes and radar
%   altitudes. (Pressure altitude was not used because of potential
%   variability with local conditions).  Unlike Hains, the tropopause
%   pressure was not fixed, but taken from the OMI pixel being validated
%   against.
%
%   This requires 3 inputs: Merge - a data structure (produced by
%   read_merge_data.m) that contains all the data for one day of flight
%   campaigns. Data is one top-level index of one of the Data data
%   structures output by BEHR_main.m. Timezone must be one of the four
%   standard US timezones (est, cst, mst, or pst) and is needed to compare
%   UTC time to OMI overpass time.
%
%   This outputs pixel longitude and latitude, OMI and BEHR satellite NO2
%   columns for those pixels, the corresponding aircraft inferred columns,
%   and the fraction of spiral measurments that fall within the pixel
%   boundary.
%
%   The sixth output is a 16-bit integer that is a quality flag for the
%   column, similar to the vcdQualityFlag and XTrackFlag in the OMNO2
%   product.  Use bitget() to check the status of individual bits; the
%   meaning of each bit is specified here:
%       1st: Summary, set to 1 if any flags are set.
%       2nd: Reserved as a second summary bit against future need.
%       3rd: Indicates the column top was derived from a daily composite
%           rather than extrapolation of the top median NO2 values.
%       4th: Indicates that < 10% of the data points in the profile had NO2
%           data
%       5th: Indicates that the composite profile had < 10% of the data
%           points valid
%       6th: Indicates that radar values from higher in the column were
%           used to calculate the surface pressure.
%       7th: Indicates that the GLOBE terrain database was used to find the
%           surface pressure because no radar data was available
%       8th: No data points fall within the time frame of +/- 3 hr from OMI
%           overpass
%       9th: No BEHR data for this swath (probably used an OMI_SP only
%           file)
%       10-15: Unused
%       16th: Set if the column was skipped due to < 1% valid
%           NO2, pressure, or temperature data
%
%   Parameters:
%
%   profiles: Allows the user to pass either (1) the name of
%   the field containing profile ID numbers in the Merge structure, or(2)
%   an (n x 2) matrix containing the start and stop times (in seconds after
%   midnight UTC) of the periods during the flight when the aircraft is
%   sprialing.  The first is useful if the profile numbers field is not
%   recognized by this function; the second is useful for campaigns such as
%   ARCTAS-CA that do not identify spirals.
%
%   starttime: Profiles must have a start time later that this. Pass as a
%   string using military time, e.g. 16:00 instead of 4:00 pm.  This is
%   always in local standard time.  Set to 10:45 by default.  Allows user
%   to restrict aircraft data to times near satellite overpass.
%
%   endtime: Profiles must have a start time before this. See starttime for
%   more details.  Set to 16:45 by default.
%
%   no2field: Defaults to 'NO2_LIF', if this is not the NO2 field, use this
%   parameter to override that.
%
%   altfield: Defaults to ALTP, which is the correct field for pressure
%   altitude for DISCOVER campaigns.
%
%   radarfield: Defaults to the correct radar altitude field name for a
%   DISCOVER campaign based on the date of the merge file.  Use this field
%   to override that selection.
%
%   presfield: Defaults to PRESSURE.
%
%   tempfield: Defaults to TEMPERATURE.
%
%   cloud_product: Which cloud product (omi, modis, or rad) to use in rejecting
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
%   clean: Setting this to 0 will not remove any pixels with a fill value -
%   useful only for debugging why a pixel is rejected.
%
%   Josh Laughner <joshlaugh5@gmail.com> 4 Jul 2014

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst','auto'})));
p.addParamValue('starttime','10:45',@isstr);
p.addParamValue('endtime','16:45',@isstr);
p.addParamValue('profiles',[], @(x) size(x,2)==2 || ischar(x));
p.addParamValue('no2field','',@isstr);
p.addParamValue('altfield','ALTP',@isstr);
p.addParamValue('radarfield','',@isstr);
p.addParamValue('presfield','PRESSURE',@isstr);
p.addParamValue('tempfield','TEMPERATURE',@isstr)
p.addParamValue('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis','rad'})));
p.addParamValue('cloud_frac_max',0.2, @isscalar);
p.addParamValue('rowanomaly','AlwaysByRow',@(x) any(strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'})));
p.addParamValue('DEBUG_LEVEL',1,@isscalar);
p.addParamValue('clean',1,@isscalar);

p.parse(Merge,Data,timezone,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
starttime = pout.starttime;
endtime = pout.endtime;
profiles = pout.profiles;
no2field = pout.no2field;
altfield = pout.altfield;
radarfield = pout.radarfield;
presfield = pout.presfield;
Tfield = pout.tempfield;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
DEBUG_LEVEL = pout.DEBUG_LEVEL;
clean_bool = pout.clean;

merge_datenum = datenum(Merge.metadata.date);

if isempty(profiles)
    spiral_mode = 'profnum';
    if merge_datenum >= datenum('07/01/2011') && merge_datenum <= datenum('07/31/2011');
        profnum = Merge.Data.ProfileSequenceNum.Values; % DISCOVER-AQ in Baltimore
    elseif merge_datenum >= datenum('01/16/2013') && merge_datenum <= datenum('02/06/2013');
        profnum = Merge.Data.ProfileNumber.Values; % DISCOVER-AQ in CA
    elseif merge_datenum >= datenum('09/04/2013') && merge_datenum <= datenum('09/29/2013');
        profnum = Merge.Data.ProfileNumber.Values; % DISCOVER-AQ in Texas
    else
        error('sprial_ver:profnum','Profiles not identified. Pass a profile number field name or UTC time ranges as the parameter ''profiles''');
    end
elseif ischar(profiles)
    spiral_mode = 'profnum';
    profnum = eval(sprintf('Merge.Data.%s.Values',profiles));
elseif ismatrix(profiles);
    spiral_mode = 'utcranges';
    Ranges = profiles;
end

if isempty(no2field)
    if merge_datenum >= datenum('07/01/2011') && merge_datenum <= datenum('07/31/2011');
        no2field = 'NO2_LIF';
    elseif merge_datenum >= datenum('01/16/2013') && merge_datenum <= datenum('02/06/2013');
        no2field = 'NO2_MixingRatio_LIF';
    elseif merge_datenum >= datenum('09/04/2013') && merge_datenum <= datenum('09/29/2013');
        no2field = 'NO2_MixingRatio_LIF';
    end
end

if isempty(radarfield)
    if merge_datenum >= datenum('07/01/2011') && merge_datenum <= datenum('07/31/2011');
        radarfield = 'A_RadarAlt';
    elseif merge_datenum >= datenum('01/16/2013') && merge_datenum <= datenum('02/06/2013');
        radarfield = 'Radar_Altitude';
    elseif merge_datenum >= datenum('09/04/2013') && merge_datenum <= datenum('09/29/2013');
        radarfield = 'Radar_Altitude';
    else
        error('sprial_ver:profnum','Radar Alt not identified. Pass a radar altitude field name as the parameter ''radarfield''');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     LOAD DATA     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First handle the aircraft data
[no2, utc, alt, lon, lat] = remove_merge_fills(Merge,no2field,'alt',altfield);
no2(no2<0) = NaN; % Must remove any negative values from consideration because they will return imaginary components during the log-log interpolation
llfill = Merge.Data.LATITUDE.Fill; xxll = lon == llfill | lon == llfill-360 | lat == llfill; % Make sure there are no fill values in longitude or latitude;
lon(xxll) = NaN; lat(xxll) = NaN;
radar_alt = remove_merge_fills(Merge,radarfield,'alt',altfield);
pres = remove_merge_fills(Merge,presfield,'alt',altfield);
temperature = remove_merge_fills(Merge,Tfield,'alt',altfield);
altfill = eval(sprintf('Merge.Data.%s.Fill',altfield));
alt(alt==altfill) = NaN; % Switching to GPS altitude gave fill values for altitude.  These must be removed.


% If the timezone was set to "auto," calculate the difference from UTC
% based on the mean longitude
if strcmpi(tz,'auto')
    tz = round(nanmean(lon)/15);
end

% Check what percentage of values were fill values, if >99% are fills for
% data, temperature, or pressure, return NaNs and set the largest bit on
% the quality flag to 1 (as well as the summary flag). If <99% but >90%,
% set the 5th bit to 1 as a warning.
percent_nans = sum(isnan(no2))/numel(no2);
percent_nans_P = sum(isnan(pres))/numel(pres);
percent_nans_T = sum(isnan(temperature))/numel(temperature);
if percent_nans > 0.99 || percent_nans_P > 0.99 || percent_nans_T > 0.99;
    if DEBUG_LEVEL > 1 && percent_nans > 0.99;
        fprintf('%s had < 1%% of NO2 values valid\n',datestr(merge_datenum,2));
    end
    if DEBUG_LEVEL > 1 && percent_nans_P > 0.99;
        fprintf('%s had < 1%% of pressure values valid\n',datestr(merge_datenum,2));
    end
    if DEBUG_LEVEL > 1 && percent_nans_T > 0.99;
        fprintf('%s had < 1%% of temperature values valid\n',datestr(merge_datenum,2));
    end
    % We must skip this merge file if there is no NO2, pressure, or
    % temperature data
    omi_lon_out = NaN;
    omi_lat_out = NaN;
    omi_no2_out = NaN;
    behr_no2_out = NaN;
    air_no2_out = NaN;
    
    db.all_profiles = {NaN};
    db.quality_flags = {uint16(2^15+1)};
    db.coverage_fraction = {NaN};
    db.dist_vectors = {NaN};
    db.loncorn = {NaN(4,1)};
    db.latcorn = {NaN(4,1)};
    db.strat_NO2 = {NaN};
    db.modis_cloud = {NaN};
    if strcmp(spiral_mode,'profnum'); db.profnums = {NaN};
    else db.profnums = {NaN(1,2)};
    end
    db.reject = uint8(2);
else
    if percent_nans > 0.9;
        warning('Merge file for %s has %.0f%% NO2 values as NaNs',datestr(merge_datenum,2),percent_nans*100);
        q_base = uint16(bitset(0,5,1));
    else
        q_base = uint16(0);
    end
    
    
    % Make a composite profile for all data within 3 hours of OMI overpass.
    % Find the median of the 10 highest altitude NO2 values (along with their
    % associated temperature and pressure values) and append these as the top
    % of tropopause value.
    tt = utc >= local2utc('10:45',tz) & utc <= local2utc('16:45',tz);
    if sum(tt) == 0 % If no points fall within the time frame, return NaNs and exit
        omi_lon_out = NaN;
        omi_lat_out = NaN;
        omi_no2_out = NaN;
        behr_no2_out = NaN;
        air_no2_out = NaN;
        
        db.all_profiles = {NaN};
        q_base = bitset(q_base,1,1); db.quality_flags = {bitset(q_base,8,1)};
        db.coverage_fraction = {NaN};
        db.dist_vectors = {NaN};
        db.loncorn = {NaN(4,1)};
        db.latcorn = {NaN(4,1)};
        db.strat_NO2 = {NaN};
        db.modis_cloud = {NaN};
        if strcmp(spiral_mode,'profnum'); db.profnums = {NaN};
        else db.profnums = {NaN(1,2)};
        end
        db.reject = uint8(2);
        return
    end
    [no2_composite, pres_composite] = bin_omisp_pressure(pres(tt),no2(tt));
    temp_composite = bin_omisp_pressure(pres(tt),temperature(tt));
    
    % The top bin of bin_omisp_pressure (200 hPa) is right around the normal
    % boundary of the troposphere, 12 km.  If no NO2 data is available, (i.e.
    % that bin has a value of NaN), then we'll extrapolate the median of the
    % top ten NO2 measurements to that bin.  In the loop itself, we'll adjust
    % for the changing tropopause pressure.
    if isnan(no2_composite(end))
        M1 = sortrows([pres', no2', temperature']);
        top = find(~isnan(M1(:,2)),10,'first'); % Since lower pressure = higher altitude, we want the first 10 NO2 measurements when sorted by pressure.
        no2_comp_top_med = median(M1(top,2)); temp_comp_top_med = nanmedian(M1(top,3)); pres_comp_top_med = nanmedian(M1(top,1));
        if no2_comp_top_med < 3; 
            no2_comp_top_med = 1.5; 
            if DEBUG_LEVEL > 0; fprintf('Composite profile top < LoD\n'); end
        end
        no2_composite(end) = no2_comp_top_med;
        
        % Temperature should have a linear relationship to altitude
        % (i.e., log(pressure)) up to the tropopause; therefore extrapolate
        % based on the fit.
        nans = isnan(temp_composite);
        temp_composite(nans) = interp1(log(pres_composite(~nans)), temp_composite(~nans), log(pres_composite(nans)),'linear','extrap');
        
        % Do any interpolation in log-log space (per. Bucsela et. al. J.
        % Geophys. Res. 2008, 113, D16S31) since pressure has an exponential
        % relation to altitude.  Log-log space would assume that NO2 also has
        % an exponential relationship to altitude.
        [~, no2_tmp] = fill_nans(log(pres_composite),log(no2_composite),'noclip');
        no2_composite = exp(no2_tmp);
        
        
    end
    
    % Identify all spirals according to the 'profiles' input; reject any
    % without a start time between 10:45 and 4:45, that is, about 3 hours on
    % either side of the OMI overpass
    if strcmp(spiral_mode,'profnum')
        % Get all unique profile numbers and their start times
        unique_profnums = unique(profnum(profnum~=0));
        start_times = zeros(numel(unique_profnums),1);
        for a=1:numel(unique_profnums)
            xx = profnum == unique_profnums(a);
            start_times(a) = min(utc(xx));
        end
        
        % Remove from consideration any profiles with a start time before 10:45
        % am or after 4:45 pm local standard time
        yy = start_times >= local2utc(starttime,tz) & start_times <= local2utc(endtime,tz);
        unique_profnums = unique_profnums(yy); start_times = start_times(yy);
        
        % Save each profile's NO2, altitude, radar altitude, latitude, and
        % longitude as an entry in a cell array
        s = size(unique_profnums);
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s); profnum_array = cell(s);
        pres_array = cell(s); temp_array = cell(s);
        for a=1:numel(unique_profnums)
            xx = profnum == unique_profnums(a);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            profnum_array{a} = unique_profnums(a);
        end
    elseif strcmp(spiral_mode,'utcranges')
        % Find all the utc start times that are between 10:45 and 4:45 local
        % standard time
        yy = Ranges(:,1) >= local2utc(starttime,tz) & Ranges(:,1) <= local2utc(endtime,tz);
        ranges_in_time = Ranges(yy,:);
        s = [1,sum(yy)];
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s);
        pres_array = cell(s); temp_array = cell(s); profnum_array = cell(s);
        for a=1:s(2)
            xx = utc >= ranges_in_time(a,1) & utc <= ranges_in_time(a,2);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            profnum_array{a} = ranges_in_time(a,:);
        end
    end
    
    % Then get pixel information from BEHR: column NO2, pixel corners, etc. We
    % will not consider any pixels that were affected by the row anomaly, have
    % too great a cloud fraction, etc.
    
    Data.Areaweight = ones(size(Data.Longitude));
    if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
    try % Handles the case of both BEHR files and SP-only files
        Data2 = omi_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
    catch err;
        if strcmp(err.identifier,'MATLAB:nonExistentField');
            Data2 = omi_sp_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
        else
            rethrow(err)
        end
    end
    xx = Data2.Areaweight > 0;
    
    omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
    corner_lat = Data2.Latcorn(:,xx); corner_lon = Data2.Loncorn(:,xx);
    omi_no2 = Data2.ColumnAmountNO2Trop(xx);
    TropopausePres = Data2.TropopausePressure(xx);
    vza = Data2.ViewingZenithAngle(xx);
    TerrainPres = Data2.TerrainPressure(xx); %load terrain pressure (in hPa)
    % Try to load the BEHR column, if it is a field in the Data structure
    try
        behr_no2 = Data2.BEHRColumnAmountNO2Trop(xx);
    catch err
        % If not, fill the imported variable with fill values
        if strcmp(err.identifier,'MATLAB:nonExistentField')
            if DEBUG_LEVEL > 0; fprintf('    No BEHR data for this swath\n'); end
            q_base = bitset(q_base,9,1);
            behr_no2 = -127*ones(size(xx));
        else
            rethrow(err)
        end
    end
    % Try to load the MODIS cloud fraction - not needed within this script,
    % but included in the output structure "db" to examine if cloud
    % fraction has an impact on the retrieval
    try
        modis_cloud = Data2.MODISCloud(xx);
    catch err
        % If "MODISCloud" is not a field, or if it only contains a single
        % value of "0" or is empty (thus the field was created but never
        % had anything entered), fill with import variable with fill
        % values.
        if strcmp(err.identifier,'MATLAB:nonExistentField') || (strcmp(err.identifier,'MATLAB:badsubscript') && isempty(Data2.MODISCloud)) || (strcmp(err.identifier,'MATLAB:badsubscript') && numel(Data2.MODISCloud) == 1 && Data2.MODISCloud == 0)
            if DEBUG_LEVEL > 0; fprintf('    No MODIS data for this swath\n'); end
            q_base = bitset(q_base,9,1);
            modis_cloud = -127*ones(size(xx));
        else
            rethrow(err)
        end
    end
    % Extra fields carried through for curiosity; this is used to calculate
    % stratospheric NO2
    total_omi_no2 = Data2.ColumnAmountNO2(xx);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     MATCH PIXELS AND SPIRALS     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize matrices to hold the values for each pixel measured.  
    n = numel(no2_array);
    npix = numel(omi_lat);
    
    omi_lon_out = -9e9*ones(npix,1);
    omi_lat_out = -9e9*ones(npix,1);
    omi_no2_out = -9e9*ones(npix,1);
    behr_no2_out = -9e9*ones(npix,1);
    air_no2_out = -9e9*ones(npix,1);
    
    db.all_profiles = cell(npix,1);
    db.quality_flags = cell(npix,1);
    db.coverage_fraction = cell(npix,1);
    db.dist_vectors = cell(npix,1);
    db.loncorn = cell(npix,1);
    db.latcorn = cell(npix,1);
    db.strat_NO2 = cell(npix,1);
    db.modis_cloud = cell(npix,1);
    db.profnums = cell(npix,1);
    db.reject = cell(n,1);
    
    
    % Also, define Avogadro's number and gas constant;
    Av = 6.022e23; % molec. / mol
    R = 8.314e4; % (hPa * cm^3) / (mol * K)
    
    % Now iterate through each profile, and find all pixels which that profile
    % overlaps. As Hains did, we will require that the pixel contains
    % measurements between 0-3 km, specifically (in addition to what was
    % described in her paper), we will require that there be 20 measurements
    % below 3 km for that pixel, and that the pixel's VZA be less than 60
    % degrees.
    
    for pix=1:numel(omi_lat)
        omi_lon_p = omi_lon(pix); omi_lat_p = omi_lat(pix);
        loncorn_p = corner_lon(:,pix); latcorn_p = corner_lat(:,pix);
        omi_no2_p = omi_no2(pix); behr_no2_p = behr_no2(pix);
        tropopause_p = TropopausePres(pix); vza_p = vza(pix);
        modis_cloud_p = modis_cloud(pix); total_omi_no2_p = total_omi_no2(pix);
        
        pix_reject = uint8(0);
        % Skip pixels with corners near +/- 180 if the corners have
        % differing signs; this will confuse the algorithm that matches
        % spirals to pixels
        if any(abs(loncorn_p)>179) && abs(mean(sign(loncorn_p))) ~= 1
            pix_reject = bitset(pix_reject,4);
            db.reject{pix} = pix_reject;
            if ~clean_bool
                omi_lon_out(pix) = omi_lon_p;
                omi_lat_out(pix) = omi_lat_p;
            end
            continue
        end
        
        % Check the vza
        if vza_p > 60; 
            pix_reject = bitset(pix_reject,5);
            db.reject{pix} = pix_reject;
            if ~clean_bool
                omi_lon_out(pix) = omi_lon_p;
                omi_lat_out(pix) = omi_lat_p;
            end
            continue; 
        end
        
        % Initialize matrices to hold values that will be returned per
        % profile
        pix_aircraft_no2 = []; pix_coverage = []; pix_flags = []; pix_profnums = []; pix_distances = [];
        pix_q_flag = q_base;
        % Check if any profiles have any part within this pixel, if so,
        % integrate that profile and average the column value to the pixel
        for p=1:n
            % Now check whether the profile is inside the pixel using
            % logical operations, which are much faster than inpolygon().
            lontest = lon_array{p}; lontest(isnan(lontest)) = [];
            lattest = lat_array{p}; lattest(isnan(lattest)) = [];
            if all(lontest < min(loncorn_p)) || all(lontest > max(loncorn_p)) || all(lattest < min(latcorn_p)) || all(lattest > max(latcorn_p))
                continue
            % If the profile does not go below 500 m, skip it anyway (per
            % Hains, require for good BL sampling)
            elseif min(radar_array{p}) > 0.5
                pix_reject = bitset(pix_reject,1);
                continue
            % If the entire profile was fill values (i.e. instrument trouble) obviously we have to skip this profile.
            elseif all(isnan(no2_array{p}));
                pix_reject = bitset(pix_reject,2);
                continue
            end
            
            % Finally actually check if the profile falls in the pixel
            % using inpolygon(). Recall that we require there to be 20
            % valid measurements in the lowest 3 km.
            no2_3km = no2_array{p}(alt_array{p}<3);
            lon_3km = lon_array{p}(alt_array{p}<3);
            lat_3km = lat_array{p}(alt_array{p}<3);
            IN_3km = inpolygon(lon_3km, lat_3km, loncorn_p, latcorn_p);
            if sum(~isnan(no2_3km(IN_3km)))<20
                pix_reject = bitset(pix_reject,3);
                continue % Pixel must have 20 valid measurments between 0-3 km altitude (good sampling of boundary layer)
            end
            
            % Initialize the qualtity flag for this profile
            q_flag = pix_q_flag;
            
            % Bin the NO2 data by pressure.
            [no2bins, presbins] = bin_omisp_pressure(pres_array{p}, no2_array{p});
            % Get the standard error of the bins
            [~,~,no2stderr] = bin_omisp_pressure(pres_array{p}, no2_array{p}, 'binmode','mean');
            % Bin the temperature
            [tempbins, temp_presbins] = bin_omisp_pressure(pres_array{p}, temp_array{p});
            
            
            % Get the top 10 NO2 measurements; if their median value is < 100
            % pptv, assume that the profile has sampled the free troposphere
            % and can be extrapolated to the tropopause safely.  If not, then
            % we will need to use the composite profile to fill in the bins
            % above the profile top. At the same time, get the bottom 10
            % NO2 measurements; we'll need them to extrapolate to the
            % surface.
            M = sortrows([pres_array{p}', no2_array{p}', radar_array{p}', temp_array{p}', alt_array{p}']);
            xx = find(~isnan(M(:,2)),10,'last'); zz = find(~isnan(M(:,2)),10,'first');
            bottom_med_no2 = median(M(xx,2)); top_med_no2 = median(M(zz,2));
            
            %Hains substitutes 1.5 ppt for any median no2 mixing ratios
            %less than the LoD (3 ppt)
            if top_med_no2 < 3; 
                top_med_no2 = 1.5; 
                if DEBUG_LEVEL>0; fprintf('Profile top median NO2 < LoD\n'); end
            end
            
            bottom_med_temp = nanmedian(M(xx,4)); top_med_temp = nanmedian(M(zz,4));
            bottom_med_radar_alt = nanmedian(M(xx,3)); bottom_med_pres = nanmedian(M(xx,1));
            bottom_med_GPS_alt = nanmedian(M(xx,5));
            % There is a chance that the radar system wasn't working at the
            % same time as the NO2 measurments, so if there were no
            % corresponding radar measurements in the lowest part of the
            % column, take whatever lowest 10 are available.
            if isnan(bottom_med_radar_alt);
                yy = find(~isnan(M(:,3)),10,'last');
                bottom_med_radar_alt = nanmedian(M(yy,3)); bottom_med_GPS_alt = nanmedian(M(yy,5));
                q_flag = bitset(q_flag,6,1);
            end
            % If the radar altitude is now a valid number, use it to get
            % the surface pressure.  Otherwise, use the GLOBE database and
            % find the nearest surface altitude, which must be converted to
            % kilometers
            if ~isnan(bottom_med_radar_alt)
                surface_alt = bottom_med_GPS_alt-bottom_med_radar_alt; surface_pres = 1013*exp(-surface_alt/7.4);
            else
                if DEBUG_LEVEL > 0; fprintf('  Retrieving GLOBE surface altitude.  May take a second...\n'); end
                surface_alt = nearest_GLOBE_alt(nanmean(lon_array{p}), nanmean(lat_array{p}))/1000;
                surface_pres = 1013*exp(-surface_alt/7.4);
                q_flag = bitset(q_flag,7,1);
            end
            
            if sum(~isnan(no2_array{p}))/numel(no2_array{p}) < 0.1
                q_flag = bitset(q_flag,4,1); % Set the 4th bit as a flag if less than 10% of the data points are non-fill values
            end
            
            % If the top median no2 value is < 100 pptv, then it's safe to
            % assume that the profile has crossed the boundary layer and is
            % sampling the free troposphere.  If not, then extrapolating that
            % value to the tropopause will grossly overestimate the total
            % column.  In the latter case, append the composite profile on top
            % of the current one.
            if top_med_no2 < 100;
                no2bins(end) = top_med_no2;
                tempbins(end) = top_med_temp;
            else
                % Get the altitude of the highest bin in the profile and find
                % all bins in the composite profile above that
                xx = pres_composite < min(presbins(~isnan(no2bins)));
                no2bins(xx) = no2_composite(xx);
                tempbins(xx) = temp_composite(xx);
                
                % Set the 3rd bit of the quality flag to 1 to indicate that a
                % composite column was appended
                q_flag = bitset(q_flag,3,1);
            end
            
            % Fill in any NaNs with interpolated values
            [tmp_pres, tmp_no2] = fill_nans(log(presbins),log(no2bins),'noclip');
            no2bins = exp(tmp_no2); presbins = exp(tmp_pres);
            [~, tempbins] = fill_nans(log(temp_presbins),tempbins,'noclip');
            % Find the surface pressure, then compare it to the bottom bin with
            % NO2 measurements.  If it is above the bottom bin center, reset
            % that bin center to the surface pressure.  If below, then we'll
            % extrapolate the median lowest 10 NO2 measurements to surface
            % pressure.
            bb = find(~isnan(no2bins),1,'first');
            % Restrict the three bins to start from those that have NO2
            % values
            no2bins = no2bins(bb:end);
            presbins = presbins(bb:end);
            tempbins = tempbins(bb:end);
            % If the surface pressure is less (i.e. above) the second
            % remaining bin, check with the user to proceed.
            if surface_pres < presbins(2);
                queststring = sprintf('Surface P (%.4f) less than second bin (%.4f). \nLow alt radar nan flag is %d. \n Continue?',surface_pres, presbins(2), bitget(q_flag,6));
                choice = questdlg(queststring,'Surface pressure','Yes','No','Abort run','No');
                switch choice
                    case 'Yes'
                        cc = presbins < surface_pres;
                        no2bins = no2bins(cc);
                        presbins = presbins(cc);
                        tempbins = tempbins(cc);
                    case 'No'
                        continue
                    case 'Abort run'
                        error('spiral_ver:surface_pres','Surface pressure < second bin.');
                end
            elseif surface_pres < presbins(1); % if surface is above the bottom bin center...
                presbins(1) = surface_pres;
                nans = isnan(tempbins);
                tempbins(nans) = interp1(log(pres_composite(~nans)), tempbins(~nans), log(pres_composite(nans)),'linear','extrap'); % Just in case the temperature data doesn't cover all the remaining bins, interpolate it
            else
                no2nans = isnan(no2bins);
                no2bins = no2bins(~no2nans); tempbins = tempbins(~no2nans); presbins = presbins(~no2nans);
                no2bins = [bottom_med_no2, no2bins];
                tempbins = [bottom_med_temp, tempbins];
                presbins = [surface_pres, presbins];
            end
            
            % Insert the OMI tropopause pressure as the final pressure
            % bin center, then convert all pressure bins to altitude,
            % interpolate, and integrate.
            presbins(end) = tropopause_p;
            altbins = -log(presbins ./ 1013) * 7.4;
            
            % Interpolate the NO2, temperature, and pressure data
            dz = 1; % integration segments in meters
            alt_profile = altbins(1):(dz/1000):altbins(end);
            no2_profile = interp1(altbins,no2bins,alt_profile,'linear');
            temp_profile = interp1(altbins,tempbins,alt_profile,'linear');
            pres_profile = exp(interp1(altbins,log(presbins),alt_profile,'linear')); % Linearly interpolate ln(P) since that is what depends linearly on altitude
            
            % Carry out the numerical integration
            no2_column = 0;
            for z=1:numel(alt_profile)
                P_z = pres_profile(z); T = temp_profile(z); no2_z = no2_profile(z);
                conc_NO2 = (Av * P_z * no2_z * 1e-12)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
                no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
            end
            pix_aircraft_no2(end+1) = no2_column;
            
            % If any bits in the quality flag are set, set the summary
            % bit; then append the quality flag to all those for this pixel
            if any(q_flag); q_flag = bitset(q_flag,1,1); end
            pix_flags(end+1) = q_flag;
            
            % Calculate what percentage of the profile actually falls in
            % this pixel; append to all values for this pixel
            IN_all = inpolygon(lon_array{p}, lat_array{p}, loncorn_p, latcorn_p);
            pix_coverage(end+1) = sum(IN_all)/numel(no2_array{p});
            
            % Append the profnum number so that this can be used for
            % counting the number of profiles found
            pix_profnums = [pix_profnums; profnum_array{p}];
            
            % Calculate the distance from the mean of the boundary layer
            % part of the profile to the pixel center
            vec = [omi_lon_p - nanmean(lon_3km), omi_lat_p - nanmean(lat_3km)];
            pix_distances(end+1) = norm(vec);
        
        end % End the loop over all profiles
        
        % Save the results for this pixel
        
        if ~isempty(pix_aircraft_no2); 
            omi_lon_out(pix) = omi_lon_p;
            omi_lat_out(pix) = omi_lat_p;
            omi_no2_out(pix) = omi_no2_p;
            behr_no2_out(pix) = behr_no2_p;
            air_no2_out(pix) = nanmean(pix_aircraft_no2);
        end
        
        db.all_profiles{pix} = pix_aircraft_no2;
        db.quality_flags{pix} = pix_flags;
        db.coverage_fraction{pix} = pix_coverage;
        db.dist_vectors{pix} = pix_distances;
        db.latcorn{pix} = latcorn_p;
        db.loncorn{pix} = loncorn_p;
        db.strat_NO2{pix} = total_omi_no2_p - omi_no2_p;
        db.modis_cloud{pix} = modis_cloud_p;
        if ~isempty(pix_profnums); db.profnums{pix} = pix_profnums;
        else db.profnums{pix} = [];
        end
        db.reject{pix} = pix_reject;
    
    end % End the loop over all pixels
    
    % Clean up the output variables
    if clean_bool
        fills = omi_lon_out == -9e9;
        omi_lon_out = omi_lon_out(~fills);
        omi_lat_out = omi_lat_out(~fills);
        omi_no2_out = omi_no2_out(~fills);
        behr_no2_out = behr_no2_out(~fills);
        air_no2_out = air_no2_out(~fills);
        
        db.all_profiles = db.all_profiles(~fills);
        db.quality_flags = db.quality_flags(~fills);
        db.coverage_fraction = db.coverage_fraction(~fills);
        db.dist_vectors = db.dist_vectors(~fills);
        db.latcorn = db.latcorn(~fills);
        db.loncorn = db.loncorn(~fills);
        db.strat_NO2 = db.strat_NO2(~fills);
        db.modis_cloud = db.modis_cloud(~fills);
        db.profnums = db.profnums(~fills);
        db.reject = db.reject(~fills);
    end
end
end