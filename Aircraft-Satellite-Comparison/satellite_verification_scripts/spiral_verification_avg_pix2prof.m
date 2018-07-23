function [ prof_lon_out, prof_lat_out, omi_no2_out, behr_no2_out, air_no2_out, db ] = spiral_verification_avg_pix2prof( Merge, Data, timezone, varargin )
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
%   Unlike spiral_verification which averages all profiles to the pixel
%   they cross into, this function averages all measurements for the pixels
%   that cover a particular profile together.
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
%   The sixth output is a structure containing a wide range of ancillary
%   data that might be valuable for debugging, such as individual pixel
%   measurements, cloud fractions, etc. All fields contain cell arrays
%   where each cell corresponds to a profile. Most fields are self
%   explanatory, but two require some additional explanation.
%       The quality_flags field is a 16-bit integer that is a quality flag
%   for the column, similar to the vcdQualityFlag and XTrackFlag in the OMNO2
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
%       8th: No data points fall within the specified time frame
%       9th: No BEHR data for this swath (probably used an OMI_SP only
%           file)
%       10th: A warning that the flight path is considered to be in
%           multiple time zones.
%       11th: Data is present in the specified time range, but profiles or
%           ranges are not
%       12th: Indicates that ground site NO2 was used for this profile.
%       13-15: Unused
%       16th: Set if the column was skipped due to < 1% valid
%           NO2, pressure, or temperature data
%
%       The reject field contains more specific information on why a
% profile or pixel was rejected.  The difference between this and the
% quality_flags field is that the quality_flags field tries to represent
% problems with the underlying data (and generally involves the whole merge
% dataset), while the reject field describes reasons that a profile or the
% pixels impinging on it were not considered for comparison.  Generally, if
% a cell contains a matrix, each entry represents one pixel impinging on
% the profile. It is an 8-bit integer; the bits have the following
% meanings:
%       1st: Profile did not reach minRadarAlt of the ground (from Hains
%           et. al., recommended 0.5 km) 
%       2nd: Entire profile was NaNs 
%       3rd: Less than numBLpoints data points in the lowest 3 km (from Hains et.
%          al., recommended number is 20).
%       4th: No valid pixels (clouds, row anomaly) overlap profile
%       5th: VZA for this pixel > 60 deg.
%       6th: Profile height was too little (max(alt) - min(alt) <
%          min_height). The required height is set with the parameter
%          'min_height'
%
%   Parameters:
%
%   campaign_name: a string identifying the campaign, used to retrieve the
%   appropriate field names. See merge_field_names.m in the Utils/Constants
%   folder.  The string need only contain an identifiable campaign name,
%   note that the Discover campaigns need the string 'discover' plus the
%   state abbreviation.  An empty string will bypass this, note that
%   no2field, altfield, radarfield, and aerfield will need to be set then
%   or an error will be thrown.
%
%   profiles: Allows the user to pass either (1) the name of the field
%   containing profile ID numbers in the Merge structure, (2) an (n x 2)
%   matrix containing the start and stop times (in seconds after midnight
%   UTC) of the periods during the flight when the aircraft is sprialing,
%   or The first is useful if the profile numbers field is not recognized
%   by this function; the second is useful for campaigns such as ARCTAS-CA
%   that do not identify spirals.
%
%   profnums: A list of specific profile numbers to examine, usually used
%   to examine only certain profiles that have been pre-selected, for e.g.
%   their aerosol characteristics.  Defaults to an empty matrix, which
%   means all profiles are included.
%
%   behrfield: What field to use for BEHR NO2 data. Defaults to
%   BEHRColumnAmountNO2Trop, but can be reset to, for example, use
%   reprocessed columns using in-situ AMFs.
%
%   starttime: Profiles must have a start time later that this. Pass as a
%   string using military time, e.g. 16:00 instead of 4:00 pm.  This is
%   always in local standard time.  Set to 10:45 by default.  Allows user
%   to restrict aircraft data to times near satellite overpass.
%
%   endtime: Profiles must have a start time before this. See starttime for
%   more details.  Set to 16:45 by default.
%
%   no2field: If unspecified, passed an empty string, or set to 'lif' will
%   default to the LIF (our data) field, field name determined by calling
%   merge_field_names with the campaign name specified. If set to 'cl',
%   will use the chemiluminescence data - NCAR if available, NOAA
%   otherwise. Otherwise will use the field specified by the string given.
%
%   conv_fact: conversion factor that, when multiplied by NO2 data, will
%   convert that data into unscaled mixing ratio (parts-per-part). This is
%   not used unless CONVERT_UNITS() cannot convert the unit given in the
%   Merge structure, in which case, if this is not given, an error is
%   thrown.
%
%   altfield: If unspecified or set to 'gps' uses the GPS altitude field
%   for the campaign, if set to 'pressure', uses the pressure altitude
%   field. Otherwise will use the field specified by the string given.
%
%   radarfield: If unspecified uses the default radar altitude field for
%   the given campaign, otherwise uses the field specified by this input.
%
%   presfield: Defaults to PRESSURE.
%
%   tempfield: Defaults to TEMPERATURE.
%
%   aerfield: If set to an empty string (default) will try to figure out
%   the aerosol extinction field based on the campaign name. If set to 0,
%   will not return aerosol data. If set to a string, will use that string
%   as the field name.
%
%   ssafield: If set to an empty string (default) will try to figure out
%   the aerosol extinction field based on the campaign name. If set to a
%   string, will use that string as the field name. Note that this will not
%   be returned if the aerfield value is 0.
%
%   cloud_product: Which cloud product (omi or modis or rad) to use in
%   rejecting pixels.  Defaults to omi.
%
%   cloud_frac_max: The maximum allowed geometric cloud fraction.  Defaults
%   to 0.2; recommended value for use with MODIS cloud product is 0.
%
%   rowanomaly: The method of rejecting pixels based on row anomaly.
%   Defaults to 'AlwaysByRow'.  See omi_pixel_reject or omi_rowanomaly for
%   more information on the possible choices ('AlwaysByRow', 'RowsByTime',
%   'XTrackFlags', and 'XTrackFlagsLight').
%
%   min_height: The minimum difference between the lowest and highest
%   points in the profile in kilometers. Defaults to 0, i.e. any profile
%   height.
%
%   numBLpoints: The minimum number of data points (valid, not NaNs) in the
%   lowest 3 km of the atmosphere.  Hains et. al. recommends 20, which is
%   set as the default.
%
%   minRadarAlt: The altitude (in km) above the ground which the plane must
%   go below for the profile to be used. Hains et. al. recommend 0.5 km
%   which is the default.
%
%   useground: boolean whether or not to include ground site data where
%   available. Defaults to true.  If set to 0, the 12th bit of the quality
%   flag will still represent whether ground data was available or not, it
%   just won't be used.
%
%   usecomposite: boolean whether or not to force using the composite NO2
%   profile, rather than extrapolating from the median of the top 10
%   measurements.
%
%   surf_pres_choice: controls what happens if, in a profile, the surface
%   pressure derived from the aircraft's (altitude - radar altitude) is
%   above the second data bin with valid data. Default is 'ask', meaning
%   ask the user what to do in each case. Other options are 'yes' (proceed
%   with the profile), 'no' (do not use the profile, but continue the run),
%   or 'abort run' (stop the run).
%
%   DEBUG_LEVEL: The level of output messages to write; 0 = none, 1 =
%   normal; 2 = verbose, 3 = (reserved), 4 = plot NO2 profiles colored by
%   source - intrinsic, extrapolation, composite
%
%   clean: Setting this to 0 will not remove any pixels with a fill value -
%   useful only for debugging why a pixel is rejected.
%
%   loncorn, latcorn: Sets the field name to use for pixel corners.
%   Defaults are FoV75CornerLongitude and FoV75CornerLatitude (prior to 7
%   Jul 2017, i.e. version 3 of BEHR, these were Loncorn and Latcorn,
%   respectively)
%
% Note to anyone editing this file in the future: there is a subfunction
% "setReturnVar" that should be used to create the initial instance of the
% six output variables. It can set them to the null values (nans) or
% initial fills and empty cell arrays.  If you need to add any fields to
% db, you should do it in this function.  The reason it is done this way is
% if fields are added to db in different orders depending on the
% circumstances it is returned under, Matlab won't be able to concatenate
% two instances of db with different field orders, so the calling function,
% Run_Spiral_Verification, will fail.  By centralizing where db is
% generated for the first time, this shouldn't be a problem.
%
%   Josh Laughner <joshlaugh5@gmail.com> 25 Jul 2014

p = inputParser;
p.addRequired('Merge',@isstruct);
p.addRequired('Data',@isstruct);
p.addRequired('timezone', @(x) any(strcmpi(x,{'est','cst','mst','pst','auto'})));
p.addParameter('behrfield','BEHRColumnAmountNO2Trop',@isstr);
p.addParameter('starttime','10:45',@isstr);
p.addParameter('endtime','16:45',@isstr);
p.addParameter('campaign_name','',@isstr);
p.addParameter('profiles',[], @(x) size(x,2)==2 || ischar(x));
p.addParameter('profnums',[], @(x) (isvector(x) || isempty(x) || size(x,2)==2));
p.addParameter('no2field','',@isstr);
p.addParameter('conv_fact',NaN,@isscalar);
p.addParameter('altfield','',@isstr);
p.addParameter('radarfield','',@isstr);
p.addParameter('presfield','',@isstr);
p.addParameter('tempfield','',@isstr)
p.addParameter('aerfield','', @(x) (ischar(x) || x==1 || x==0) );
p.addParameter('ssafield','', @(x) (ischar(x) || x==0));
p.addParameter('cloud_product','omi',@(x) any(strcmpi(x,{'omi','modis','rad'})));
p.addParameter('cloud_frac_max',0.2, @isscalar);
p.addParameter('rowanomaly','AlwaysByRow',@(x) any(strcmp(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'})));
p.addParameter('min_height',0,@isscalar);
p.addParameter('numBLpoints',20,@isscalar);
p.addParameter('minRadarAlt',0.5,@isscalar);
p.addParameter('useground',1,@isscalar);
p.addParameter('usecomposite',0,@isscalar);
p.addParameter('loncorn','FoV75CornerLongitude',@ischar);
p.addParameter('latcorn','FoV75CornerLatitude',@ischar);
p.addParameter('DEBUG_LEVEL',1,@isscalar);
p.addParameter('clean',1,@isscalar);
p.addParameter('surf_pres_choice','ask',@ischar);

p.parse(Merge,Data,timezone,varargin{:});
pout = p.Results;

% Check that only one element of Data was passed
if numel(Data)>1; error('bdy_layer_verify:DataInput','Only pass one top-level element of the Data structure'); end

Merge = pout.Merge;
Data = pout.Data;
tz = pout.timezone;
behrfield = pout.behrfield;
starttime = pout.starttime;
endtime = pout.endtime;
profiles = pout.profiles;
user_profnums = pout.profnums;
campaign_name = pout.campaign_name;
no2field = pout.no2field;
conv_fact = pout.conv_fact;
altfield = pout.altfield;
radarfield = pout.radarfield;
presfield = pout.presfield;
aerfield = pout.aerfield;
ssafield = pout.ssafield;
Tfield = pout.tempfield;
cloud_prod = pout.cloud_product;
cloud_frac_max = pout.cloud_frac_max;
rowanomaly = pout.rowanomaly;
min_height = pout.min_height;
numBLpoints = pout.numBLpoints;
minRadarAlt = pout.minRadarAlt;
useground = pout.useground;
force_composite = pout.usecomposite;
loncorn_field = pout.loncorn;
latcorn_field = pout.latcorn;
DEBUG_LEVEL = pout.DEBUG_LEVEL;
clean_bool = pout.clean;
surface_pres_override = pout.surf_pres_choice;

merge_datenum = datenum(Merge.metadata.date);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Many of the errors in this script have not been updated yet to use the
% new error class.  
E = JLLErrors;

% This will check for the MATLAB_DISPLAY environmental variable.
% In my .bashrc file, I've set the alias startmatlab to open a command
% line Matlab with the -nodisplay argument. This alias also sets the 
% MATLAB_DISPLAY environmental variable to be 0 so that I can check
% for it in a script and disable any graphical components.  Here, it's
% the dialogue asking what to do if the second pressure bin is less
% than the surface pressure.
MATLAB_DISPLAY = str2double(getenv('MATLAB_DISPLAY'));
if isnan(MATLAB_DISPLAY)
    MATLAB_DISPLAY = 1;
end

% As long as the campaign name isn't blank, try to get the appropriate
% field names for the given campaign. Otherwise check that the merge
% fields we need are set and if not, error out.

if ~isempty(campaign_name);
    % get the information needed to read the merge files but do not ask
    % the user to choose a range file.
    [FieldNames,~,~,~,ground_site_dir] = merge_field_names(campaign_name, false);
elseif isempty(no2field) || isempty(altfield) || isempty(aerfield) || isempty(ssafield) || isempty(radarfield) || isempty(presfield) || isempty(Tfield)
    error(E.badinput('If no campaign is to be specified (using the parameter ''campaign_name'') then parameters no2field, altfield, presfield, tempfield, radarfield, aerfield, and ssafield must not be empty strings.'));
end

% Deal with the profiles or UTC range input. Set both the spiral_mode
% string (which determines whether to use profile numbers or UTC ranges
% later on in this function) and reads in the profile numbers or UTC ranges.
if isempty(profiles)
    spiral_mode = 'profnum';
    if isempty(FieldNames.profile_numbers)
        error(E.callError('profile_mode','The profile numbers field is not defined for this campaign; be sure to set the calling function to use UTC ranges'));
    end
    profnum = Merge.Data.(FieldNames.profile_numbers).Values;
elseif ischar(profiles)
    spiral_mode = 'profnum';
    profnum = Merge.Data.(profiles).Values;
elseif ismatrix(profiles);
    spiral_mode = 'utcranges';
    Ranges = profiles;
end

% Now for each field name not given by the user, set it from the
% merge field names retrieved.
if isempty(no2field) || strcmpi(no2field,'lif')
    no2field = FieldNames.no2_lif;
elseif strcmpi(no2field,'cl')
    no2field = FieldNames.no2_ncar;
end

if isempty(altfield) || strcmpi(altfield,'gps')
    altfield = FieldNames.gps_alt;
elseif strcmpi(altfield,'pressure')
    altfield = FieldNames.pressure_alt;
end

if isempty(radarfield)
    radarfield = FieldNames.radar_alt;
end

if isempty(aerfield) % This will not override a 0 passed as this field.
    aerfield = FieldNames.aerosol_extinction;
end

if isempty(ssafield)
    ssafield = FieldNames.aerosol_dry_ssa;
end

if isempty(presfield)
    presfield = FieldNames.pressure;
end

if isempty(Tfield)
    Tfield = FieldNames.temperature;
end

% Finally check that user_profnums fullfills the requirements for the
% profile mode, be it profile numbers or UTC ranges.
if ~isempty(user_profnums)
    if strcmp(spiral_mode,'profnum')
        if ~isvector(user_profnums)
            error(E.badinput('User defined profile numbers must be a vector.'));
        end
    elseif strcmp(spiral_mode,'utcranges')
        if size(user_profnums,2) ~= 2
            error(E.badinput('User defined UTC ranges must be in an n x 2 matrix'));
        end
    end
end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     LOAD DATA     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First handle the aircraft data. Convert NO2 into unscaled mixing ratio,
% which will simplify all our later calculations.
try
    [no2, utc, alt, lon, lat] = remove_merge_fills(Merge,no2field, 'alt', altfield, 'unit', 'ppp');
catch err
    if strcmp(err.identifier, 'convert_units:unit_not_found')
        if isnan(conv_fact)
            E.callError('unknown_units', 'Units for NO2 field (%s) unknown; give a non-NaN value for the "conv_fact" parameter to specify how to convert NO2 to unscaled mixing ratio', Merge.Data.(no2field).Unit)
        else
        [no2, utc, alt, lon, lat] = remove_merge_fills(Merge,no2field, 'alt', altfield);
        no2 = no2 * conv_fact;
        warning('read_no2:unknown_unit', 'Units for NO2 field (%s) unknown, converting NO2 to mixing ratio with given conversion factor (%g)', Merge.Data.(no2field).Unit, conv_fact);
        end
    else
        rethrow(err)
    end
end

no2(no2<0) = NaN; % Must remove any negative values from consideration because they will return imaginary components during the log-log interpolation
lon(isnan(lat))=NaN; lat(isnan(lon))=NaN; % Make any points that are NaNs in one fields also so in the other
radar_alt = remove_merge_fills(Merge,radarfield,'alt',altfield,'unit','kilometers');
pres = remove_merge_fills(Merge,presfield,'unit','hPa');
temperature = remove_merge_fills(Merge,Tfield);

altfill = Merge.Data.(altfield).Fill;
alt(alt==altfill) = NaN; % Switching to GPS altitude gave fill values for altitude.  These must be removed.
alt = convert_units(alt, Merge.Data.(altfield).Unit, 'kilometers');
if aerfield ~= 0
    aer_data = remove_merge_fills(Merge,aerfield);
    ssa_data = remove_merge_fills(Merge,ssafield);
end


% The variable that holds the quality flags
q_base = uint16(0);

% If the timezone was set to "auto," calculate the difference from UTC
% based on the mean longitude. This will produce a vector of the same
% length as the utc and lon variables
if strcmpi(tz,'auto')
    tz = round(lon/15);
    
    % If all elements in the timezone vector aren't the same, set the 10th
    % bit in the quality flag
    if ~all(tz==tz(1))
        q_base = bitset(q_base,10,1);
    end
end

% Now, use the timezone (entered or calculated) to produce a new vector of
% times that reflect the local time at each point
local_times = utc2local_sec(utc,tz);

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
    setReturnVar(uint16(2^15+1),spiral_mode,aerfield, 'null', 1, DEBUG_LEVEL);
    return;

else
    if percent_nans > 0.9;
        warning('Merge file for %s has %.0f%% NO2 values as NaNs',datestr(merge_datenum,2),percent_nans*100);
        q_base = uint16(bitset(q_base,5,1));
    end
    
    
    % Make a composite profile for all data within 3 hours of OMI overpass
    % - this will look for points based on their local time.
    % Find the median of the 10 highest altitude NO2 values (along with their
    % associated temperature and pressure values) and append these as the top
    % of tropopause value.
    tt = local_times >= local2utc('10:45',0) & local_times <= local2utc('16:45',0);
    if sum(tt) == 0 % If no points fall within the time frame, return NaNs and exit
        q_base = bitset(q_base,1,1); % Set the summary and specific quality flags
        setReturnVar(bitset(q_base,8,1), spiral_mode, aerfield, 'null', 1, DEBUG_LEVEL);
        return
    end
    [no2_composite, pres_composite] = bin_omisp_pressure(pres(tt),no2(tt));
    % Calculate the uncertainty in the composite bins as the std. error
    [~,~, no2_composite_stderr] = bin_omisp_pressure(pres(tt),no2(tt),'mean');
    temp_composite = bin_omisp_pressure(pres(tt),temperature(tt));
    
    % The top bin of bin_omisp_pressure (200 hPa) is right around the normal
    % boundary of the troposphere, 12 km.  If no NO2 data is available, (i.e.
    % that bin has a value of NaN), then we'll extrapolate the median of the
    % top ten NO2 measurements to that bin.  In the loop itself, we'll adjust
    % for the changing tropopause pressure.
    wasextrap = 0;
    if isnan(no2_composite(end))
        wasextrap = find(~isnan(temp_composite),1,'last');
        
        M1 = sortrows([pres', no2', temperature']);
        top = find(~isnan(M1(:,2)),10,'first'); % Since lower pressure = higher altitude, we want the first 10 NO2 measurements when sorted by pressure.
        no2_comp_top_med = median(M1(top,2)); temp_comp_top_med = nanmedian(M1(top,3)); pres_comp_top_med = nanmedian(M1(top,1));
        
        % Assume the error is the standard error of the top ten points. NO2
        % values should be given in mixing ratio (parts per part).
        no2_comp_top_med_stderr = std(M1(top,2))/sqrt(10);

        if no2_comp_top_med < 3e-12; 
            no2_comp_top_med = 1.5 * 1e-12; 
            if DEBUG_LEVEL > 0; fprintf('Composite profile top < LoD\n'); end
        end
        no2_composite(end) = no2_comp_top_med;
        no2_composite_stderr(end) = no2_comp_top_med_stderr;
        
        % For the free troposphere we'll use the full campaign composite.
        nans = isnan(temp_composite);
        lastbin = find(~nans,1,'last');
        nans(1:lastbin)=false;
        tmp = campaign_temperature_profile(campaign_name);
        temp_composite(nans) = tmp(nans);
        
        
        if isnan(temp_composite(1))
            nans2 = isnan(temp_composite);
            temp_composite(nans2) = interp1(log(pres_composite(~nans2)),temp_composite(~nans2),log(pres_composite(nans2)),'linear','extrap');
        end
        
        
        % Do any interpolation in log-log space (per. Bucsela et. al. J.
        % Geophys. Res. 2008, 113, D16S31) since pressure has an exponential
        % relation to altitude.  Log-log space would assume that NO2 also has
        % an exponential relationship to altitude.
        [~, no2_tmp, no2_tmperr] = fill_nans(log(pres_composite),log(no2_composite),log(no2_composite_stderr),'noclip');
        no2_composite = exp(no2_tmp);
        no2_composite_stderr = exp(no2_tmperr);    
        
        
    end
    
    % Identify all spirals according to the 'profiles' input; reject any
    % without a start time between 10:45 and 4:45, that is, about 3 hours on
    % either side of the OMI overpass
    if strcmp(spiral_mode,'profnum')
        % Get all unique profile numbers and their start times
        unique_profnums = unique(profnum(profnum~=0));
        start_times = zeros(numel(unique_profnums),1);
        end_times = zeros(numel(unique_profnums),1);
        prof_tzs = cell(numel(unique_profnums),1); % This is only set if using profnums b/c it's only purpose (now) is for ground sites - which only exist when there are profnums.
        for a=1:numel(unique_profnums)
            % If we're using a vector of timezones, get the most common
            % timezone from the profile - we'll treat the whole profile as
            % belonging to that timezone.  If timezone has been set
            % manually, then just use that time zone to convert the start
            % time of the profile. (obviously) Save the local start time.
            xx = profnum == unique_profnums(a);
            if ismatrix(tz) && isnumeric(tz)
                % Case where we're using a vector of timezones
                mct = mode(tz(xx));
                start_times(a) = utc2local_sec(min(utc(xx)),mct);
                end_times(a) = utc2local_sec(max(utc(xx)),mct);
                prof_tzs{a} = mct;
            elseif ischar(tz)
                start_times(a) = utc2local_sec(min(utc(xx)),tz);
                end_times(a) = utc2local_sec(max(utc(xx)),tz);
                prof_tzs{a} = tz;
            else
                error(E.callError('tz_not_recognized','Cannot understand the format of timezones'));
            end
        end
        
        % Remove from consideration any profiles with a start time before 10:45
        % am or after 4:45 pm local standard time
        yy = start_times >= local2utc(starttime,0) & start_times <= local2utc(endtime,0);
        
        % If no profiles are left, return null values and exit
        if sum(yy) == 0
            setReturnVar(bitset(q_base,11,1), spiral_mode, aerfield, 'null', 1, DEBUG_LEVEL);
            return
        end
        
        unique_profnums = unique_profnums(yy); 
        start_times = start_times(yy); 
        end_times = end_times(yy);
        prof_tzs = prof_tzs(yy);
        
        % If the user passed a list of profile numbers, remove any profile
        % numbers that do not match the list provided.
        if ~isempty(user_profnums)
            tmp = true(size(unique_profnums));
            for a = 1:numel(unique_profnums)
                tmp(a) = any(unique_profnums(a)==user_profnums);
            end
            unique_profnums(~tmp) = [];
            prof_tzs(~tmp) = [];
            start_times(~tmp) = [];
            end_times(~tmp) = [];
        end
        
        % Save each profile's NO2, altitude, radar altitude, latitude, and
        % longitude as an entry in a cell array
        s = size(unique_profnums);
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s); profnum_array = cell(s);
        pres_array = cell(s); temp_array = cell(s); aer_array = cell(s);
        ssa_array = cell(s);
        for a=1:numel(unique_profnums)
            xx = profnum == unique_profnums(a);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            if aerfield ~= 0; aer_array{a} = aer_data(xx); end
            if ssafield ~= 0; ssa_array{a} = ssa_data(xx); end
            profnum_array{a} = unique_profnums(a);
        end
    elseif strcmp(spiral_mode,'utcranges')
        % Find all the utc start times that are between within the
        % specified range of local times.  Go through each range, find the
        % data points that correspond to it, get the most common timezone,
        % use that to set whether to include that range or not. Also, check
        % the "user_profnums" variable which will have specific UTC ranges
        % to allow
        yy = false(size(Ranges,1),1);
        for a=1:size(Ranges,1)
            if ~isempty(user_profnums) && ~(ismember(Ranges(a,1),user_profnums(:,1)) && ismember(Ranges(a,2),user_profnums(:,2)))
                % If the UTC range is not one defined in the user_profnums
                % (which is really user_ranges in this case, but it was
                % called profnums first), skip the rest of this loop, which
                % will not set yy(a) to true thus skipping this range.
                continue
            end
            tz_ind = utc >= Ranges(a,1) & utc <= Ranges(a,2);
            mct = mode(tz(tz_ind));
            range_start_local = utc2local_sec(Ranges(a,1),mct);
            yy(a) = range_start_local >= local2utc(starttime,0) && range_start_local <= local2utc(endtime,0);
        end
        
        % If no profiles are left, return null values and exit
        if sum(yy) == 0
            setReturnVar(bitset(q_base,11,1), spiral_mode, aerfield, 'null',1,DEBUG_LEVEL);
            return
        end
        
        ranges_in_time = Ranges(yy,:);
        start_times = utc2local_sec(Ranges(yy,1)',0);
        end_times = utc2local_sec(Ranges(yy,2)',0);
        s = [1,sum(yy)];
        no2_array = cell(s); alt_array = cell(s); radar_array = cell(s);
        lat_array = cell(s); lon_array = cell(s); aer_array = cell(s);
        pres_array = cell(s); temp_array = cell(s); profnum_array = cell(s);
        ssa_array = cell(s);
        for a=1:s(2)
            xx = utc >= ranges_in_time(a,1) & utc <= ranges_in_time(a,2);
            no2_array{a} = no2(xx);
            alt_array{a} = alt(xx);
            radar_array{a} = radar_alt(xx);
            lat_array{a} = lat(xx);
            lon_array{a} = lon(xx);
            pres_array{a} = pres(xx);
            temp_array{a} = temperature(xx);
            if aerfield ~= 0; aer_array{a} = aer_data(xx); end
            if ssafield ~= 0; ssa_array{a} = ssa_data(xx); end
            profnum_array{a} = ranges_in_time(a,:);
        end
        
        % prof_tzs is used later to add ground site data for DISCOVER
        % campaigns. If we are using a UTC range file, we are more than
        % likely in a campaign that doesn't have ground site data. So,
        % check that we are not using ground data (error if we are) and
        % make a dummy prof_tz cell array.
        if ~useground
            prof_tzs = cell(size(profnum_array));
        else
            E.notimplemented('%s','Trying to use ground site data and a UTC range is not supported, since only DISCOVER campaigns have ground site data. Set the ''useground'' parameter to 0')
        end
    end
    
    % Then get pixel information from BEHR: column NO2, pixel corners, etc. We
    % will not consider any pixels that were affected by the row anomaly, have
    % too great a cloud fraction, etc.
    
    Data.Areaweight = ones(size(Data.Longitude));
    if isfield(Data, 'BEHRQualityFlags')
        % The BEHRQualityFlags field was added in v3.0A of BEHR, so as long
        % as it's present, we can use the newest version of
        % omi_pixel_reject.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
        reject_detail = struct('cloud_type', cloud_prod, 'cloud_frac', cloud_frac_max,...
            'row_anom_mode', rowanomaly, 'check_behr_amf', true);
        Data2 = omi_pixel_reject(Data,'detailed',reject_detail);
    elseif isfield(Data, 'BEHRColumnAmountNO2Trop')
        % If there's a BEHR NO2 column field but not the BEHR quality flags
        % field, we're probably working with version 2 data.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject_v2.m\n'); end
        Data2 = omi_pixel_reject_v2(Data, cloud_prod, cloud_frac_max, rowanomaly);
    else
        % If the BEHR NO2 column field isn't in Data, it's definitely not a
        % BEHR file, so use the SP reject field.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_sp_pixel_reject.m\n'); end
        Data2 = omi_sp_pixel_reject(Data,cloud_prod,cloud_frac_max,rowanomaly);
    end
    xx = Data2.Areaweight > 0;
    
    omi_lat = Data2.Latitude(xx); omi_lon = Data2.Longitude(xx);
    corner_lat = Data2.(latcorn_field)(:,xx); corner_lon = Data2.(loncorn_field)(:,xx);
    omi_no2 = Data2.ColumnAmountNO2Trop(xx);
    omi_cloudfrac = Data2.CloudFraction(xx);
    omi_cloudradfrac = Data2.CloudRadianceFraction(xx);
    TropopausePres = Data2.TropopausePressure(xx);
    vza = Data2.ViewingZenithAngle(xx);
    TerrainPres = Data2.TerrainPressure(xx); %load terrain pressure (in hPa)
    sza = Data2.SolarZenithAngle(xx);
    % Try to load the BEHR column if present. If not, it is probably a
    % standard product file.
    try
        behr_no2 = Data2.(behrfield)(xx);
        
    catch err
        % If not, fill the imported variable with fill values
        if strcmp(err.identifier,'MATLAB:nonExistentField')
            if DEBUG_LEVEL > 0; fprintf('    No BEHR data for this swath\n'); end
            q_base = bitset(q_base,9,1);
            behr_no2 = -127*ones(size(xx));
            alb = -127*ones(size(xx));
        else
            rethrow(err)
        end
    end
    % Try to load the MODIS cloud fraction and albedo - not needed within
    % this script, but included in the output structure "db" to examine if
    % cloud fraction has an impact on the retrieval or albedo. Prefer the
    % BRDF albedo, if present.
    try
        modis_cloud = Data2.MODISCloud(xx);
        if isfield(Data2,'MODISAlbedoBRDF')
            alb = Data2.MODISAlbedoBRDF(xx);
        else
            alb = Data2.MODISAlbedo(xx);
        end
    catch err
        % If "MODISCloud" is not a field, or if it only contains a single
        % value of "0" or is empty (thus the field was created but never
        % had anything entered), fill with import variable with fill
        % values.
        if strcmp(err.identifier,'MATLAB:nonExistentField') || (strcmp(err.identifier,'MATLAB:badsubscript') && isempty(Data2.MODISCloud)) || (strcmp(err.identifier,'MATLAB:badsubscript') && numel(Data2.MODISCloud) == 1 && Data2.MODISCloud == 0)
            if DEBUG_LEVEL > 0; fprintf('    No MODIS data for this swath\n'); end
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
    
    % Sets the initial fill values for prof_lon_out, prof_lat_out,
    % omi_no2_out, behr_no2_out, air_no2_out, and the structure db.  This
    % helps ensure that the fields for db are always in the right order.
    % Any fields to be added to db should be added in the setReturnVar
    % subfunction at the end of this file.
    setReturnVar(bitset(q_base,11,1), spiral_mode, aerfield, 'init', n, DEBUG_LEVEL);
    
    % All fields in db are set to an empty cell array with size (n,1), 
    db.start_times = num2cell(start_times);
    db.end_times = num2cell(end_times);
    
    
    % Also, define Avogadro's number and gas constant;
    Av = 6.022e23; % molec. / mol
    R = 8.314e4; % (hPa * cm^3) / (mol * K)
    
    % Now iterate through each profile, and find all pixels which that profile
    % overlaps. As Hains did, we will require that the pixel contains
    % measurements between 0-3 km, specifically (in addition to what was
    % described in her paper), we will require that there be 20 measurements
    % below 3 km for that pixel, and that the pixel's VZA be less than 60
    % degrees.
    
    for p=1:n
        % Initialize the qualtity flag for this profile
        q_flag = q_base;
        
        prof_reject = uint8(0);
        
        lontest = lon_array{p}; lontest(isnan(lontest)) = [];
        lattest = lat_array{p}; lattest(isnan(lattest)) = [];
        % If the profile does not go below 500 m, skip it anyway (per
        % Hains, require for good BL sampling)
        if min(radar_array{p}) > minRadarAlt
            prof_reject = bitset(prof_reject,1); 
            db.reject{p} = prof_reject;
            continue
            % If the entire profile was fill values (i.e. instrument trouble) obviously we have to skip this profile.
        elseif all(isnan(no2_array{p}));
            prof_reject = bitset(prof_reject,2); 
            db.reject{p} = prof_reject;
            continue
        elseif (max(alt_array{p}) - min(alt_array{p})) < min_height
            prof_reject = bitset(prof_reject,6);
            db.reject{p} = prof_reject;
            continue
        end
        % Find all the pixels that impinge on this profile
        lat_logical = max(lattest) > min(corner_lat,[],1) & min(lattest) < max(corner_lat,[],1);
        lon_logical = max(lontest) > min(corner_lon,[],1) & min(lontest) < max(corner_lon,[],1);
        latlon_logical = lat_logical & lon_logical;
        omi_lon_p = omi_lon(latlon_logical); omi_lat_p = omi_lat(latlon_logical);
        loncorn_p = corner_lon(:,latlon_logical); latcorn_p = corner_lat(:, latlon_logical);
        omi_no2_p = omi_no2(latlon_logical); behr_no2_p = behr_no2(latlon_logical);
        tropopause_p = TropopausePres(latlon_logical); vza_p = vza(latlon_logical);
        terrain_pres_p = TerrainPres(latlon_logical);
        omi_cloudfrac_p = omi_cloudfrac(latlon_logical);
        omi_cloudradfrac_p = omi_cloudradfrac(latlon_logical);
        modis_cloud_p = modis_cloud(latlon_logical); total_omi_no2_p = total_omi_no2(latlon_logical);
        sza_p = sza(latlon_logical); alb_p = alb(latlon_logical);
        
        % Check each pixel for rejection criteria.
        pix_xx = true(size(omi_no2_p)); pix_reject = prof_reject*uint8(ones(size(pix_xx)));
        if isempty(pix_xx); pix_reject = bitset(prof_reject,4,1); end
        
        % Skip pixels with corners if the corners have differing signs and
        % are not near 0; this will confuse the algorithm that matches
        % spirals to pixels.  Check along the first dimension of the corner
        % variables. First check that pix_xx isn't an empty matrix, which will happen
        % if latlon_logical is all 0's - because of how the loncorn is indexed, the
        % same latlon_logical will make omi_no2_p an empty matrix, and loncorn_p an
        % empty 4 x 0 matrix, which will confuse the & operation
        if ~isempty(pix_xx)
            pix_xx = pix_xx & ~((any(abs(loncorn_p)>20) & abs(mean(sign(loncorn_p),1))~=1)');
        end
        
        % Check the vza
        vza_xx = vza_p <= 60;
        pix_xx = pix_xx & vza_xx;
        pix_reject(~vza_xx) = bitset(pix_reject(~vza_xx),5*uint8(ones(size(find(~vza_xx)))),1);
        
        % Finally actually check if the profile falls in the pixel
        % using inpolygon(). Recall that we require there to be 20
        % valid measurements in the lowest 3 km.
        no2_3km = no2_array{p}(alt_array{p}<3);
        lon_3km = lon_array{p}(alt_array{p}<3);
        lat_3km = lat_array{p}(alt_array{p}<3);
        pix_coverage = zeros(size(omi_no2_p));
        for pix=1:size(loncorn_p,2)
            IN_3km = inpolygon(lon_3km, lat_3km, loncorn_p(:,pix), latcorn_p(:,pix));
            if sum(~isnan(no2_3km(IN_3km)))<numBLpoints
                pix_xx(pix) = false;
                pix_reject(pix) = bitset(pix_reject(pix),3,1);
                continue % Pixel must have 20 valid measurments between 0-3 km altitude (good sampling of boundary layer)
            end
            % Calculate what percentage of the profile actually falls in
            % this pixel; append to all values for this pixel
            IN_all = inpolygon(lon_array{p}, lat_array{p}, loncorn_p(:,pix), latcorn_p(:,pix));
            pix_coverage(pix) = sum(IN_all)/numel(no2_array{p});
        end
        dum=1;
        % If no valid pixels are left, skip this profile.
        if sum(pix_xx)==0;
            continue
        end
        
        %Otherwise, cut down the vectors representing pixels that do match
        %this profile to the pixels that should be considered
        omi_lon_p = omi_lon_p(pix_xx); omi_lat_p = omi_lat_p(pix_xx);
        loncorn_p = loncorn_p(:,pix_xx); latcorn_p = latcorn_p(:, pix_xx);
        omi_no2_p = omi_no2_p(pix_xx); behr_no2_p = behr_no2_p(pix_xx);
        tropopause_p = tropopause_p(pix_xx); vza_p = vza_p(pix_xx);
        omi_cloudfrac_p = omi_cloudfrac_p(pix_xx); omi_cloudradfrac_p = omi_cloudradfrac_p(pix_xx);
        modis_cloud_p = modis_cloud_p(pix_xx); total_omi_no2_p = total_omi_no2_p(pix_xx);
        pix_coverage = pix_coverage(pix_xx); terrain_pres_p = terrain_pres_p(pix_xx);
        sza_p = sza_p(pix_xx); alb_p = alb_p(pix_xx);
        
        % Calculate the distance vector between the mean lat/lon of the
        % lowest 3 km of the profile and each pixel left.  This can be used
        % later for weighting.
        dist_vectors = zeros(size(omi_lat_p));
        prof_latlon = [nanmean(lon_3km), nanmean(lat_3km)];
        for a=1:numel(dist_vectors)
            omi_latlon = [omi_lon_p(a), omi_lat_p(a)]; 
            dist_vectors(a) = norm((prof_latlon - omi_latlon));
        end
        
        % Bin the NO2 data by pressure.
        [no2bins, presbins] = bin_omisp_pressure(pres_array{p}, no2_array{p});
        % Get the standard error of the bins
        [~,~,no2stderr] = bin_omisp_pressure(pres_array{p}, no2_array{p}, 'binmode','mean');
        % Bin the temperature
        [tempbins, ~] = bin_omisp_pressure(pres_array{p}, temp_array{p});
        
        
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
        % Calculate the standard error for the top and bottom median points
        bottom_med_no2_stderr = std(M(xx,2))/sqrt(10);
        top_med_no2_stderr = std(M(zz,2))/sqrt(10);
        
        if DEBUG_LEVEL > 3;
            f2 = figure;
            plot(M(xx,2),M(xx,3),'color','b','linestyle','none','marker','o');
            line(M(zz,2),M(zz,3),'color','r','linestyle','none','marker','s');
            xlabel('NO2 mixing ratio'); ylabel('Radar altitude');
            title(sprintf('%s: Profile %d',Merge.metadata.date,profnum_array{p}));
            fprintf('\nPaused\n');
            pause;
            close(f2);
        end
        
        %Hains substitutes 1.5 ppt for any median no2 mixing ratios
        %less than the LoD (3 ppt). NO2 values should be given in mixing
        %ratio (parts-per-part)
        if top_med_no2 < 3e-12;
            top_med_no2 = 1.5e-12;
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
        
        if DEBUG_LEVEL > 3; f1 = figure; plot(no2bins, presbins, 'color','b','linewidth',10); end
        % If the top median no2 value is < 100 pptv, then it's safe to
        % assume that the profile has crossed the boundary layer and is
        % sampling the free troposphere.  If not, then extrapolating that
        % value to the tropopause will grossly overestimate the total
        % column.  In the latter case, append the composite profile on top
        % of the current one.
        if top_med_no2 < 100 * 1e-12 && ~force_composite;
            no2bins(end) = top_med_no2;
            no2stderr(end) = top_med_no2_stderr;
            topcol = [0 0.7 0];
            
        else
            % Get the altitude of the highest bin in the profile and find
            % all bins in the composite profile above that
            xx = pres_composite < min(presbins(~isnan(no2bins)));
            no2bins(xx) = no2_composite(xx);
            no2stderr(xx) = no2_composite_stderr(xx);
            topcol = [0.7 0 0.7];
            
            % Set the 3rd bit of the quality flag to 1 to indicate that a
            % composite column was appended
            q_flag = bitset(q_flag,3,1);
        end
        
        % Fill in any NaNs with interpolated values
        [tmp_pres, tmp_no2] = fill_nans(log(presbins),log(no2bins),'noclip');
        no2bins = exp(tmp_no2); presbins = exp(tmp_pres);
        
        % Always take temperature from the composite profile for free
        % troposphere bins
        nans = isnan(tempbins);
        lastbin = find(~isnan(tempbins),1,'last');
        nans(1:lastbin) = false;
        tempbins(nans) = temp_composite(nans);
        
        if DEBUG_LEVEL > 3; figure(f1); line(no2bins,presbins,'color',topcol,'linewidth',4); end
        
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
        no2stderr = no2stderr(bb:end);
        
        % If the surface pressure is less (i.e. above) the second
        % remaining bin, check with the user to proceed.
        if surface_pres < presbins(2);
            if strcmpi(surface_pres_override, 'ask')
                queststring = sprintf('Surface P (%.4f) less than second bin (%.4f). \nLow alt radar nan flag (radar and NO2 measurements not coincident) is %d. \n Continue?',surface_pres, presbins(2), bitget(q_flag,6));
                if MATLAB_DISPLAY
                    choice = questdlg(queststring,'Surface pressure','Yes','No','Abort run','No');
                else
                    queststring = sprintf('%s: [Yes] / No / Abort run',queststring);
                    choice = input(queststring,'s');
                end
                choice = lower(choice);
            else
                choice = lower(surface_pres_override);
            end
            switch choice
                case 'no'
                    if DEBUG_LEVEL > 0; fprintf('\tSkipping this column because surface pressure is smaller than the second bin with valid data.\n'); end 
                    continue
                case 'abort run'
                    error('spiral_ver:surface_pres','Surface pressure < second bin.');
                otherwise % Default to assuming yes.
                    cc = presbins < surface_pres;
                    no2bins = no2bins(cc);
                    presbins = presbins(cc);
                    tempbins = tempbins(cc);
                    no2stderr = no2stderr(cc);
            end
        elseif surface_pres < presbins(1); 
            % if surface is above the bottom bin center but below the
            % second bin, just reset the first bin pressure to be at the
            % surface.
            presbins(1) = surface_pres;
            nans = isnan(tempbins);
            tempbins(nans) = interp1(log(pres_composite(~nans)), tempbins(~nans), log(pres_composite(nans)),'linear','extrap'); % Just in case the temperature data doesn't cover all the remaining bins, interpolate it
        else
            % Otherwise, add a "surface bin" to interpolate to. Initially,
            % use the median of the bottom 10 NO2 measurements. If there is
            % ground site data available, we'll add that in a bit.
            no2nans = isnan(no2bins);
            no2bins = no2bins(~no2nans); tempbins = tempbins(~no2nans); presbins = presbins(~no2nans); no2stderr = no2stderr(~no2nans);
            no2bins = [bottom_med_no2, no2bins];
            no2stderr = [bottom_med_no2_stderr, no2stderr];
            
            tempbins = [bottom_med_temp, tempbins];
            presbins = [surface_pres, presbins];
        end
        
        % Will insert ground NO2 as the bottom bin, otherwise returns
        % no2bins and no2stderr unchanged.
        [no2bins, no2stderr, q_flag] = add_ground_no2(no2bins, no2stderr, q_flag, profnum_array{p}, prof_tzs{p}, spiral_mode, start_times(p), end_times(p), ground_site_dir, useground, Merge, FieldNames, DEBUG_LEVEL);
        
        
        if DEBUG_LEVEL > 3
            figure(f1); line(no2bins,presbins,'color','red');
            %scatter_errorbars(no2bins,presbins,no2stderr,'direction','x','color','k','linewidth',2);
            title(sprintf('%s: profile %d',Merge.metadata.date,profnum_array{p}));
            set(gca,'YDir','reverse');
            fprintf('\nPaused\n');
            pause;
            close(f1);
        end
            
        % Insert the average OMI tropopause pressure of the pixels
        % considered as the final pressure bin center, then convert all
        % pressure bins to altitude, interpolate, and integrate.
        presbins(end) = nanmean(tropopause_p);
        altbins = -log(presbins ./ 1013) * 7.4;
        
        % Interpolate the NO2, temperature, and pressure data
        dz = 1; % integration segments in meters
        alt_profile = altbins(1):(dz/1000):altbins(end);
        no2_profile = interp1(altbins,no2bins,alt_profile,'linear');
        temp_profile = interp1(altbins,tempbins,alt_profile,'linear');
        pres_profile = exp(interp1(altbins,log(presbins),alt_profile,'linear')); % Linearly interpolate ln(P) since that is what depends linearly on altitude
        
        % Fill in the standard error nans, but we don't need to interpolate
        % to every 100 cm altitude point.
        [~,~,no2stderr] = fill_nans(altbins,no2bins,no2stderr,'noclip');
        
        % Carry out the numerical integration, 
        no2_column = 0;
        for z=1:numel(alt_profile)
            P_z = pres_profile(z); T = temp_profile(z); no2_z = no2_profile(z);
            conc_NO2 = (Av * P_z * no2_z)/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
            no2_column = no2_column + conc_NO2 * dz * 100; % Integrating in 100dz cm increments
        end
        
        % The error propagation will be handled by considering this to be
        % effectively a trapezoid rule implementation.  The math behind
        % this is in my BEHR notebook from 4 Feb 2015 - Josh Laughner
        
        % Number density of air at each bin center
        Nair = (presbins .* Av)./(R .* tempbins);
        % Calculates the error for each trapezoid using vectors, then sum
        % up that vector and square root it to find the total column error.
        % Since the concentrations are defined in molec./cm^3, we need to
        % convert altitude bins from km to cm (hence the factor of 1e5).
        % Also, NO2 values are usually reported (by us at least) in pptv,
        % hence the conv_fact defaults to 1e-12.
        column_error = sqrt( sum( ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* no2stderr(2:end).^2 .* Nair(2:end).^2 +...
            ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* no2stderr(1:end-1).^2 .* Nair(1:end-1).^2 ) );
        % Double check that this number is a scalar. This will catch if a
        % matrix (rather than a vector) slips into this calculation.
        if ~isscalar(column_error)
            error(E.badvartype('column_error', 'scalar'));
        end
        
        
        % Save the output for this profile
        prof_lon_out(p) = nanmean(lon_array{p});
        prof_lat_out(p) = nanmean(lat_array{p});
        omi_no2_out(p) = nanmean(omi_no2_p);
        behr_no2_out(p) = nanmean(behr_no2_p);
        air_no2_out(p) = no2_column;
        
        % Bin the aerosol data for this profile and find the max extinction
        % value and total integrated aerosol extinction.
        if aerfield ~= 0;
            % Check the quality of the aerosol data in this profile. If the
            % entire profile is NaNs, skip it
            aer_quality = uint8(0);
            if all(isnan(aer_array{p}))
                % Set the first bit in aer_quality to 1 if all data is NaNs
                aer_quality = bitset(aer_quality,1);
                aer_max = -9e9;
                aer_int = -9e9;
                aer_prof_height = -9e9;
            else
                if sum(isnan(aer_array{p}))/numel(aer_array{p}) > 0.9
                    % Set the second bit if >90% of the profile has NaNs
                    % for aerosol data.
                    aer_quality = bitset(aer_quality,2);
                end
                [aerbins] = bin_omisp_pressure(pres_array{p}, aer_array{p});
                aer_max = max(aerbins);
                if isempty(aer_max);
                    aer_max = -9e9;
                end
                
                binwidth = 0.25; % km
                [aerintbins, aerbinmid] = bin_vertical_profile(alt_array{p}, aer_array{p}, binwidth);
                [aerbinmid, aerintbins] = fill_nans(aerbinmid,aerintbins);
                if numel(aerintbins) == 1
                    % Set the third bit if only one bin had actual data
                    aer_quality = bitset(aer_quality,3);
                end
                % Integrate, converting the aerosol extinction from Mm^-1
                % to m^-1 and the bin altitudes from km to m.
                aer_int = trapz(aerbinmid*1e3, aerintbins*1e-6);
                % Calculate the height of the aerosol profile as the
                % difference between the first and last non-NaN bin. This
                % will give us a criterion to gauge the comparability of
                % the AOD.
                aer_bottom = find(~isnan(aerintbins),1,'first');
                aer_top = find(~isnan(aerintbins),1,'last');
                aer_prof_height = aerbinmid(aer_top) - aerbinmid(aer_bottom);
            end
            
            % Check the quality of SSA data. If it is all nans, set the
            % fifth bit. If > 90% is nans, set the sixth bit. In the first
            % case, return fill values.
            if all(isnan(ssa_array{p}))
                aer_quality = bitset(aer_quality,5);
                aer_ssa = -9e9;
            elseif sum(isnan(ssa_array{p}))/numel(ssa_array{p}) > 0.9
                aer_quality = bitset(aer_quality,6);
                aer_ssa = nanmedian(ssa_array{p});
            else
                aer_ssa = nanmedian(ssa_array{p});
            end
        end
        
        % If any bits in the quality flag are set, set the summary
        % bit; then append the quality flag to all those for this pixel
        if any(q_flag); q_flag = bitset(q_flag,1,1); end
        
        db.all_omi{p} = omi_no2_p;
        db.all_behr{p} = behr_no2_p;
        db.sza_avg{p} = nanmean(sza_p);
        db.alb_avg{p} = nanmean(alb_p);
        db.quality_flags{p} = q_flag;
        db.coverage_fraction{p} = pix_coverage;
        db.dist_vectors{p} = dist_vectors;
        db.latcorn{p} = latcorn_p;
        db.loncorn{p} = loncorn_p;
        db.strat_NO2{p} = total_omi_no2_p - omi_no2_p;
        db.omi_cloudfrac{p} = omi_cloudfrac_p;
        db.omi_cloudradfrac{p} = omi_cloudradfrac_p;
        db.modis_cloud{p} = modis_cloud_p;
        db.profnums{p} = profnum_array{p};
        db.reject{p} = pix_reject;
        db.lon_3km{p} = lon_3km;
        db.lat_3km{p} = lat_3km;
        %db.start_times has alreadfy been handled
        if aerfield ~= 0; 
            db.aer_max_out{p} = aer_max; 
            db.aer_int_out{p} = aer_int;
            db.aer_quality{p} = aer_quality;
            db.aer_median_ssa{p} = aer_ssa; 
            db.aer_prof_height{p} = aer_prof_height;
        end
        db.column_error{p} = column_error;
        
    end % End the loop over all profiles

    % Clean up the output variables
    if clean_bool
        [prof_lon_out, prof_lat_out, omi_no2_out, behr_no2_out, air_no2_out, db] = cleanUpOutput(prof_lon_out, prof_lat_out, omi_no2_out, behr_no2_out, air_no2_out, db);
    end
end
end

function [lon_out, lat_out, omi_out, behr_out, air_out, db_out] = cleanUpOutput(lon_out, lat_out, omi_out, behr_out, air_out, db_out)
    fills = lon_out == -9e9 | lat_out == -9e9 | omi_out == -9e9 | behr_out == -9e9 | air_out == -9e9;
    lon_out = lon_out(~fills);
    lat_out = lat_out(~fills);
    omi_out = omi_out(~fills);
    behr_out = behr_out(~fills);
    air_out = air_out(~fills);
    
    db_fields = fieldnames(db_out);
    for f=1:numel(db_fields)
        db_out.(db_fields{f}) = db_out.(db_fields{f})(~fills);
    end
end

function setReturnVar(q_flag_val,spiral_mode,aerfield,null_or_init,n,DEBUG_LEVEL)
% Sets fill return values in the event that this function needs to return
% early. Needs the value for q_flag (since that represents why the function
% is returning early) and the value for aerfield, to determine whether to
% include aerosol data. null_or_init should be the string 'null' or 'init',
% which tells it whether to set the values to NaNs or the initial fill
% values. It needs to know how large to make the arrays if writing initial
% fill values, so pass n which should be the number of profiles to examine.
% If not using 'init', n can just be 1. Also uses DEBUG_LEVEL to print
% messages.

% Double check that the summary bit is set on the quality flag if needed
if q_flag_val > 0
    q_flag_val = bitset(q_flag_val,1,1);
end

if strcmpi(null_or_init,'null')
    if DEBUG_LEVEL>0; fprintf('*** Creating null return values ***\n'); end
    cellVal = {NaN};
    cellRangeVal = {NaN(1,2)};
    cornVal = {NaN(4,1)};
    regVal = NaN;
    rejectVal = {uint8(2)};
    qualVal = {q_flag_val};
elseif strcmpi(null_or_init,'init')
    if DEBUG_LEVEL>1; fprintf('*** Creating initial values ***\n'); end
    cellVal = cell(n,1);
    cellRangeVal = cell(n,1);
    cornVal = cell(n,1);
    regVal = -9e9*ones(n,1);
    rejectVal = cell(n,1);
    qualVal = cell(n,1);
end



prof_lon_out = regVal;
prof_lat_out = regVal;
omi_no2_out = regVal;
behr_no2_out = regVal;
air_no2_out = regVal;

db.all_omi = cellVal;
db.all_behr = cellVal;
db.sza_avg = cellVal;
db.alb_avg = cellVal;
db.quality_flags = qualVal;
db.coverage_fraction = cellVal;
db.dist_vectors = cellVal;
db.loncorn = cornVal;
db.latcorn = cornVal;
db.strat_NO2 = cellVal;
db.omi_cloudfrac = cellVal;
db.omi_cloudradfrac = cellVal;
db.modis_cloud = cellVal;
if strcmp(spiral_mode,'profnum'); db.profnums = cellVal;
else db.profnums = cellRangeVal;
end
db.reject = rejectVal;
db.lon_3km = cellVal;
db.lat_3km = cellVal;
db.start_times = cellVal;
db.end_times = cellVal;
if aerfield ~= 0;
    db.aer_max_out = cellVal;
    db.aer_int_out = cellVal;
    db.aer_quality = cellVal;
    db.aer_median_ssa = cellVal;
    db.aer_prof_height = cellVal;
end
db.column_error = cellVal;

assignin('caller','prof_lon_out',prof_lon_out);
assignin('caller','prof_lat_out',prof_lat_out);
assignin('caller','omi_no2_out',omi_no2_out);
assignin('caller','behr_no2_out',behr_no2_out);
assignin('caller','air_no2_out',air_no2_out);
assignin('caller','db',db);
end

function [no2bins, no2stderr, q_flag] = add_ground_no2(no2bins, no2stderr, q_flag, profnum, prof_tz, spiral_mode, start_time, end_time, ground_site_dir, useground, Merge, FieldNames, DEBUG_LEVEL)
    E = JLLErrors;
    
    if ~isempty(ground_site_dir) && strcmpi(spiral_mode,'profnum')
        sitenum = (profnum - mod(profnum,1000))/1000;
        % Handles the CA, TX campaigns where profile numbers start
        % in the hundred thousands.
        if sitenum > 100; sitenum = mod(sitenum,100); end
        % Check for any files for this site for this day
        ground_file = sprintf('*Ground%d*%s.mat',sitenum,regexprep(Merge.metadata.date,'-','_'));
        F = dir(fullfile(ground_site_dir,ground_file));
        if numel(F) > 1
            E.toomanyfile('ground site data');
        elseif numel(F) == 1 && useground
            % Load the merge file; must load it as another
            % structure because we've already got a Merge variable
            GM = load(fullfile(ground_site_dir,F(1).name));
            ground_utc = remove_merge_fills(GM.Merge,FieldNames.ground_utc{1,sitenum});
            ground_utcstop = remove_merge_fills(GM.Merge,FieldNames.ground_utc{3,sitenum});
            % If the ground NO2 is an unknown unit, this will fail.
            % Currently I don't have anything set up to allow you to pass
            % in a conversion factor for the ground NO2.
            ground_no2_vec = remove_merge_fills(GM.Merge,FieldNames.ground_no2{sitenum}, 'unit', 'ppp');


            % Convert the UTC values to local time (in seconds
            % after midnight). Use the timezone defined for this
            % profile.
            ground_utc = utc2local_sec(ground_utc,prof_tz);
            ground_utcstop = utc2local_sec(ground_utcstop,prof_tz);

            % Some sites have ~minute averaging, some have ~hourly.
            % To cover all cases, we first look for measurements
            % that have a start or end time inside the profile time
            % range, then if there's none of those, look for the
            % one that contains the start and end times.
            xx_start = ground_utc >= start_time & ground_utc < end_time;
            xx_stop = ground_utcstop > start_time & ground_utcstop <= end_time;
            xx_ground = xx_start | xx_stop;

            if sum(xx_ground) == 0
                xx_start = start_time > ground_utc & start_time < ground_utcstop;
                xx_stop = end_time > ground_utc & end_time < ground_utcstop;
                xx_ground = xx_start & xx_stop;
            end

            if sum(xx_ground) == 0
                if DEBUG_LEVEL > 0;
                    fprintf('\tNo ground measurements from site #%d overlap profile #%d\n',sitenum,profnum);
                end
            else
                ground_no2 = nanmedian(ground_no2_vec(xx_ground));
                ground_no2_stderr = nanstd(ground_no2_vec(xx_ground))/sqrt(sum(~isnan(ground_no2_vec(xx_ground))));
                if ~isnan(ground_no2) && ground_no2 > 0;
                    % Only overwrite the median values if there is a
                    % valid
                    if DEBUG_LEVEL > 1; fprintf('\tInserting ground site NO2 for site #%d\n',sitenum); end
                    no2bins(1) = ground_no2;
                    no2stderr(1) = ground_no2_stderr;
                    % Set the 12th quality bit flag here to indicate
                    % that ground NO2 was used in this profile
                    q_flag = bitset(q_flag,12,1);
                end
            end
        elseif numel(F) == 1 && ~useground
            % If we're not supposed to use the ground NO2 data but
            % we could have, still set the 12th quality bit flag.
            q_flag = bitset(q_flag,12,1);
        end
        % If there are no files, continue on with the median NO2
    end
end
