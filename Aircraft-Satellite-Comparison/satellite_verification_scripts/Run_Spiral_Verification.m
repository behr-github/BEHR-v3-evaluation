function varargout = Run_Spiral_Verification(varargin)
% RUN_SPIRAL_VERIFICATION - Driver function for profile based satellite verification
%   One method of verifying satellite retrieved NO2 columns is to calculate
%   an NO2 column from aircraft measurements of NO2 concentration by
%   integrating those concentrations over altitude (c.f. Bucsela et al., J.
%   Geophys. Res., 113, D16S31, 2008 and Hains et al., J. Geophys. Res.,
%   115, D05301, 2010). The function spiral_verification_avg_pix2prof.m
%   handles the mechanics of doing so for a single day, this function
%   handles looping over multiple days by loading the proper data etc.
%
%   This is intended to operate in one of two modes. First, all the options
%   that spiral_verification_avg_pix2prof.m needs are coded into this
%   function, with the intent that, if you were running interactively, you
%   could change those options here and not have to type them all out every
%   time. Alternately, the options can be specified as parameter inputs,
%   intended for cases where you call this from another function to prepare
%   the matched columns. If only some of the options are specified, a
%   warning is issued, since you are relying on hard-coded options that are
%   intended to be changed.
%
%   Output can be handled one of two ways:
%       RUN_SPIRAL_VERIFICATION( ... ) - with no output arguments, the
%       variables db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall,
%       airno2_iall, and dates_iall are directly placed in the base
%       workspace (description of variables below).
%
%       [LON, LAT, OMINO2, BEHRNO2, AIRNO2, DB, DATES] =
%       RUN_SPIRAL_VERIFICATION(...) - instead returns the values via the
%       typical function output mechanism.
%
%   Output values are:
%       LON - the mean longitude of the profile
%       LAT - the mean latitude of the profile
%       OMINO2 - the NO2 tVCD from the NASA standard product.
%       BEHRNO2 - the NO2 tVCD from the BEHR product.
%       AIRNO2 - the NO2 tVCD calculated from integrating the aircraft NO2
%           profile
%       DB - a structure containing extra information from the matching and
%           integrating process, which can be used to diagnose problems or
%           investigate other correlations within the data
%       DATES - a cell array that indicates the date each profile occurred
%           on
%
%   The parameter values available are:
%
%   **** Data selection variables ****
%   
%       'campaign' - a string indicating which campaign to compare against.
%       Must be understood by merge_field_names.m
%
%       'merge_dir' - the directory containing the Merge files to read
%       aircraft data from. Make an empty string to use the directory
%       specified in merge_field_names.m for the given campaign name.
%
%       'behr_dir' - the directory that the BEHR data may be found in. 
%
%       'behr_prefix' - the part of the BEHR file names before the date.
%       This does NOT understand wildcards at present.
%
%
%   **** Profile selection variables ****
%
%       'profnums' - this restrict which profiles to compare against (which
%       was used when I was looking at the effect of different aerosol
%       layers on the retrieval). This must either be a vector of profile
%       numbers (an empty one matches no profiles) or the string 'all' to
%       use all profiles.
%
%       'profile_input' - this can either be the name of the profile number
%       field in the Merge data, an empty string (which will use the
%       profile number field specified by merge_field_names.m for the given
%       campaign) or the string 'ranges' which indicates a range file
%       should be used.
%
%       'ask_ranges' - controls how merge_field_names behaves if multiple
%       range files are present. If given as a boolean value, it will force
%       merge_field_names to ask which file (true) or force it to just
%       return an empty string (false). If given a non-boolean input, it
%       will be given as the fourth input to merge_field_names with the
%       third input set to 'specified', so this allows control over which
%       Ranges file is chosen.
%
%       'startdate', 'enddate' - the start and end dates as strings, give
%       empty strings to use the values given in merge_field_names.m for
%       the campaign specified.
%
%       'starttime', 'endtime' - the earliest and latest time that profiles
%       should start between. Give as strings in 24-hr HH:MM format, in
%       local standard time.
%
%       'timezone' - the timezone to use to adjust the UTC time given in
%       the Merge file to local time. Set as 'auto' to guess the UTC offset
%       based on longitude, or see utc2local_sec for time zone
%       abbreviations.
%
%       'minheight' - the minimum height (in kilometers) between bottom and
%       top that a profile must have to be used. Set to 0 to ignore this
%       requirement.
%
%       'numBLpoints' - the minimum number of points below 3 km a profile
%       must have in order to have good boundary layer sampling. 20 is
%       recommended by Hains et al. 
%
%       'minagl' - a height above ground level (in kilometers) that a
%       profile must reach below to be used. I usually use 0.5 km.
%
%       'useground' - set to 1 to use ground site data to fill in the
%       bottom of the profile, 0 will use the median of the 10 lowest
%       measurements. Note that if ground data not available for a
%       campaign, setting this to 1 will cause an error.
%
%       'surf_pres_choice' - what to do if the surface pressure derived
%       from the aircraft (altitude - radar altitude) or GLOBE database is
%       above the second valid bin of the NO2 data. Can be 'ask', 'yes',
%       'no', or 'abort run'.
%
%
%   **** Merge field name control ****
%
%       'no2field' - the Merge.Data field containing the NO2 concentration.
%       An empty string or the string 'lif' will use the TD-LIF NO2
%       measurements, 'cl' will use chemiluminescence (specified in
%       merge_field_names). Any other string will be assumed to manually
%       specify the field name.
%
%       'no2conv' - the conversion factor for the NO2 measurements. Must be
%       a scalar number that when multiplied by the NO2 concentrations in
%       the Merge files returns a value in unscaled mixing ratio (i.e.
%       parts-per-part).
%
%       'aerfield' - the field to use for aerosol extinction. If given as
%       an empty string, will use the field specified by merge_field_names.
%
%       'ssafield' - the field to use for single-scattering albedo. If
%       given as an empty string, will use the field specified by
%       merge_field_names.
%
%       'altfield' - the field to use for altitude above sea level. Can be
%       'gps' or empty string to use the GPS altitude, 'pressure' to use
%       the pressure altitude, or a manually specified field.
%
%       'radarfield' - the field to use for altitude above ground. If given
%       as an empty string, will use the field specified by
%       merge_field_names.
%
%       'presfield' - the field to use for atmospheric pressure outside the
%       aircraft. If given as an empty string, will use the field specified
%       by merge_field_names.
%
%       'tempfield' - the field to use for ambient temperature outside the
%       aircraft. If given as an empty string, will use the field specified
%       by merge_field_names.
%
%
%   **** Satellite field control ****
%       
%       'cloudtype' - which cloud product to use to reject cloudy pixels.
%       Options are those accepted by omi_pixel_reject: 'omi' and 'modis'
%
%       'cloudfrac' - The maximum cloud fraction to allow. Usually 0.2 is
%       recommended.
%
%       'rowanomaly' - which method of rejecting pixels affected by the row
%       anomaly to use. Usually 'XTrackFlags' is recommended. See
%       omi_pixel_reject for possible values.
%
%       'behrfield' - which field in the BEHR Data structure to use for the
%       tVCD. Usually 'BEHRColumnAmountNO2Trop', or
%       'BEHRColumnAmountNO2TropVis', but in old files with columns
%       reprocessed with the in situ profiles, there is a
%       'BEHR_R_ColumnAmountNO2Trop' that used in situ profiles in the AMF
%       calculation.
%
%       'debug' - which debugging/verbosity level this and spiral_verif.
%       should use. Higher numbers print more information and generate
%       plots at debug > 2.
%
%       'clean' - boolean, controls whether profiles that have fill values
%       as the output for any of the main fields (lon, lat, behrno2,
%       omino2, or airno2) are removed (true) or kept (false). Usually
%       true, i.e. they are removed, is recommended.


E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSER %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

noinput = 'noinput';

p = inputParser;
p.addParameter('campaign', noinput);
p.addParameter('profnums', noinput);
p.addParameter('profile_input', noinput);
p.addParameter('ask_range', noinput);
p.addParameter('merge_dir', noinput);
p.addParameter('behr_dir', noinput);
p.addParameter('behr_prefix', noinput);
p.addParameter('startdate', noinput);
p.addParameter('enddate', noinput);
p.addParameter('starttime', noinput);
p.addParameter('endtime', noinput);
p.addParameter('timezone', noinput);
p.addParameter('no2field', noinput);
p.addParameter('no2conv', noinput);
p.addParameter('aerfield', noinput);
p.addParameter('ssafield', noinput);
p.addParameter('altfield', noinput);
p.addParameter('radarfield', noinput);
p.addParameter('presfield', noinput);
p.addParameter('tempfield', noinput);
p.addParameter('minheight', noinput);
p.addParameter('numBLpts', noinput);
p.addParameter('minagl', noinput);
p.addParameter('useground', noinput);
p.addParameter('cloudtype', noinput);
p.addParameter('cloudfrac', noinput);
p.addParameter('rowanomaly', noinput);
p.addParameter('behrfield', noinput);
p.addParameter('surf_pres_choice', noinput);
p.addParameter('DEBUG_LEVEL', noinput);
p.addParameter('clean', noinput);

p.parse(varargin{:});
pout = p.Results;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% The campaign name. Current valid strings are 'discover-md','discover-ca',
% 'discover-tx', 'seac4rs', 'dc3', 'arctas-b', 'arctas-carb'. This is used
% to automatically find the campaign dates, the campaign directory, and the
% data field names. If you don't want to retrieve this automatically, set
% this to an empty string.
campaign_name = set_input_value('soas', 'campaign'); % Which campaign this is for. Used to automatically find field names

% These variables are used to subset the aircraft data into profiles used
% to generate column data to compare against satellite data. 
%   'profiles' should either be set to: 
%       1) an empty string or the name of the profile number field 
%       2) the string 'ranges'. 
%   Option 1 only works with the DISCOVER campaigns so far, as those are the
% only campaigns that have specifically generated spirals labeled with
% profile numbers.  An empty string will automatically detect the correct
% field name for a campaign.  Using this option, the "profnums" variable
% can be set with an empty matrix to use all profiles, or to a matrix
% containing specific profile numbers to be examined. Such a matrix can be
% generated using categoriz_aerosol_profile in the Aerosol Effects/Profile
% Classification directory.
%   Option 2 requires a range file (a file containing the "Ranges"
% structure generated by "record_profile_ranges" in the utility_scripts
% folder) and for the path to that file to be specified.  This file
% contains a list of UTC ranges that have been specified as profiles of
% interest.  If range_file is left blank, the program will look at the
% range files returned from merge_field_names.  If there is one, that one
% will be used, otherwise the user is presented with his options.
profile_input = set_input_value('ranges', 'profile_input');
profnums = set_input_value('all', 'profnums'); % will only be used if the input profnums is 'noinput'

% Grab the dates and directory for the campaign unless the campaign name is
% empty.
ask_range = set_input_value(strcmp(profile_input, 'ranges'), 'ask_range');
if ~isempty(campaign_name)
    % If the user inputs a boolean value for ask_range, that controls
    % whether merge_field_names should ask for which file to use if
    % multiple are used or not. Otherwise, it is assumed to be a way to
    % identify which one of multiple range files to use.
    if islogical(ask_range)
        [~, dates, directory, range_file] = merge_field_names(campaign_name, ask_range);
    else
        [~, dates, directory, range_file] = merge_field_names(campaign_name, 'specified', ask_range);
    end
end

% The dates to run - leave as empty strings to automatically do the entire
% campaign, which is read from the merge_field_names function output above.
date_start = set_input_value('', 'startdate');
date_end = set_input_value('', 'enddate');

if isempty(date_start) || isempty(date_end)
    fprintf('Setting start and end dates based on the campaign name.\n');
    date_start = dates{1};
    date_end = dates{2};
end

% The directories where to find the data needed for the run.  The merge
% directory (where the campaign data is located) can be automatically
% determined from the campaign name (leave the string empty if this is
% desired), but for now the BEHR directory will be set manually.  The BEHR
% prefix is the part of the BEHR file name before the date; wildcards are
% allowed so long as there is enough of the prefix there to uniquely
% identify 1 file per date in the given directory.
merge_dir = set_input_value('', 'merge_dir');
behr_dir = set_input_value('/Volumes/share-sat/SAT/BEHR/AlbedoTestBRDF/BRDF', 'behr_dir');
behr_prefix = set_input_value('OMI_BEHR_v2-1C_', 'behr_prefix');


if isempty(merge_dir)
    fprintf('Setting merge file directory based on the campaign name.\n');
    merge_dir = directory;
end

% Start and end times (in military format) for which profiles to consider.
% General recommendation is +/-1.5 hr from overpass.
starttime = set_input_value('12:00', 'starttime');
endtime = set_input_value('15:00', 'endtime');

% Time zone (3 letter abbreviation). Set to 'auto' to determine based on
% the longitude of the data

tz = set_input_value('auto', 'timezone');


% Which fields from the merge files to use. Set them to empty strings ('')
% to automatically guess the correct field for the given campaign.  

no2field = set_input_value('lif', 'no2field'); % Which NO2 data field to use. Leave as empty string or 'lif' for our LIF data, set to 'cl' for chemiluminescence data, or any other string to override.
conv_fact = set_input_value(1e-12, 'no2conv'); % Conversion factor for NO2 data from part-per-whatever to part-per-part. Usually 1e-12, i.e. NO2 data is in pptv.
aerfield = set_input_value('', 'aerfield'); % Which aerosol extinction field to use.
ssafield = set_input_value('', 'ssafield'); % Which aerosol SSA field to use.
altfield = set_input_value('', 'altfield'); % Which altitude field to use. Can set to 'pressure' or 'gps' ('' defaults to gps), or override.
radarfield = set_input_value('', 'radarfield'); % The field for radar altitude.
presfield = set_input_value('', 'presfield'); % The field for atmospheric pressure
tempfield = set_input_value('', 'tempfield'); % The field for ambient temperature

% Variables to allow or disallow the use of a profile
min_height = set_input_value(0, 'minheight'); % The minimum difference between the top and bottom of a profile. Set to 0 to ignore (max - min > 0 always).
numBLpoints = set_input_value(20, 'numBLpts'); % The number of data points required in the bottom 3 km to ensure good BL sampling. Hains et. al. recommends 20.
minRadarAlt = set_input_value(0.5, 'minagl'); % Height above the surface (in km) a profile must be below to ensure good BL sampling. Hains et. al. recommends 0.5 km (500 m).

% Set to 1 to include ground site data, or 0 to use only aircraft data.
useground = set_input_value(0, 'useground');

% What to do if the surface pressure is above the second lowest bin with
% valid data
surf_pres_choice = set_input_value('ask', 'surf_pres_choice');


% Satellite variables
cloud_product = set_input_value('omi', 'cloudtype'); % Can be 'omi', 'modis', or 'rad'
cloud_frac_max = set_input_value(0.2, 'cloudfrac'); % Maximum cloud fraction to allow in a pixel. Recommended 0.2 for OMI, 0 for MODIS, and 0.5 for radiance.
row_anomaly = set_input_value('XTrackFlags', 'rowanomaly'); % How to reject for the row anomaly - can be 'AlwaysByRow', 'RowsByTime', 'XTrackFlags', and 'XTrackFlagsLight'
behrfield = set_input_value('BEHRColumnAmountNO2Trop', 'behrfield'); % The field in the Data structure with NO2 column data. 
                                          % Choices include 'ColumnAmountNO2Trop' (OMI SP column), 'BEHRColumnAmountNO2Trop' (BEHR column) 
                                          % and 'BEHR_R_ColumnAmountNO2Trop' only available in files where the column was reprocessed with
                                          % an AMF derived from in-situ measurements.

% Debugging variables
DEBUG_LEVEL = set_input_value(2, 'DEBUG_LEVEL'); % This will also be passed to the spiral verification function.
clean = set_input_value(1, 'clean'); % Set to 0 to keep all pixel comparisons, even those with fill values. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA PREPARATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If we are using ranges, load the range file and extract the date
% information - we'll need that to match up the correct set of ranges with
% the correct date.    
fprintf('Using %s as the range file\n',range_file);


if strcmpi(profile_input, 'ranges')
    R = load(range_file); % loads the Ranges variable, which is a data structure
    Ranges = R.Ranges;
    range_dates = cellstr(datestr({Ranges.Date},29));
    range_bool = true;
else
    range_bool = false;
end

% If no input for profile numbers is given, then use whatever is hard coded
% as profnums in this function.

if isempty(profnums)
    % An empty array for profnums input here means that no profiles
    % should be included. But to spiral_verif., that means include all
    % profiles. So instead we pass a number that will never be a
    % profile number, so no profiles are matched.
    if strcmpi(profile_input, 'ranges')
        profnums = [-127 -127];
    else
        profnums = -127;
    end
elseif strcmpi(profnums, 'all')
    % If all profiles are requested, then tell spiral_verif.
    profnums = [];
end


% This boolean will help identify whether the run succeeded, that is, if
% any files where actually loaded.
run_bool = false;


%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

dates = datenum(date_start):datenum(date_end);

S=0; %clear('dbs');
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('%s%s%s%s.mat',behr_prefix,year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    behr_files = dir(fullfile(behr_dir,behr_filename));
    if numel(behr_files)==1
        load(fullfile(behr_dir,behr_files(1).name),'Data')
    elseif isempty(behr_files)
        if DEBUG_LEVEL > 1; fprintf('No BEHR file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of BEHR files for %s is not 1 or 0',datestr(dates(d)));
    end

    % Find the UTC range data for this date, if we're using ranges instead
    % of profile numbers
    if range_bool
        xx = find(strcmp(curr_date,range_dates));
        if isempty(xx);
            if DEBUG_LEVEL > 0; fprintf('No UTC ranges found for %s, skipping\n',curr_date); end
            continue
        end
        profile_input = Ranges(xx).Ranges;
        
        % As we start considering very specific events (e.g. fires) not all
        % days will have ranges defined. 
        if isempty(profile_input); continue; end
    end
    
    % If no BEHR or Merge files are ever loaded, we won't get to this line,
    % so this will be false. This boolean will be used to call an error an
    % little later if still false, and that error will avoid the more
    % confusing error that dbs doesn't exist.
    if ~run_bool; run_bool = true; end
    
    for swath=1:numel(Data)
        S=S+1;
        [lon_i{S}, lat_i{S}, omino2_i{S}, behrno2_i{S}, airno2_i{S}, dbs(S)] = spiral_verification_avg_pix2prof(Merge,Data(swath),tz,...
            'behrfield',behrfield,...
            'starttime',starttime,... 
            'endtime',endtime,... 
            'profiles',profile_input,...
            'profnums',profnums,... 
            'campaign_name',campaign_name,...
            'no2field',no2field,... 
            'conv_fact',conv_fact,...
            'aerfield',aerfield,... 
            'ssafield',ssafield,...
            'altfield',altfield,... 
            'radarfield',radarfield,... 
            'presfield',presfield,... 
            'tempfield',tempfield,... 
            'cloud_product',cloud_product,... 
            'cloud_frac_max',cloud_frac_max,... 
            'rowanomaly',row_anomaly,... 
            'min_height',min_height,...
            'numBLpoints',numBLpoints,...
            'minRadarAlt',minRadarAlt,...
            'useground',useground,...
            'surf_pres_choice', surf_pres_choice,...
            'DEBUG_LEVEL',DEBUG_LEVEL,... 
            'clean',clean); 
        date_cell{S} = repmat({curr_date},numel(lon_i{S}),1);
    end
end

if ~run_bool;
    error(E.callError('run_failure','The loop never executed completely, meaning in no case was both a BEHR and Merge file loaded. Check the dates, file paths, and BEHR prefix.'));
end

% concatenate the output
[db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall, airno2_iall, dates_iall] = match_arrays2db(dbs, lon_i, lat_i, omino2_i, behrno2_i, airno2_i, date_cell);

% Save the inputs to the spiral_verification function; this way by saving
% the db_iall structure we can recreate this run later if necessary.
db_iall.run.date_start = date_start;
db_iall.run.date_end = date_end;
db_iall.run.merge_dir = merge_dir;
db_iall.run.behr_dir = behr_dir;
db_iall.run.behr_prefix = behr_prefix;
db_iall.run.behrfield = behrfield;
db_iall.run.range_file = range_file;
db_iall.run.starttime = starttime;
db_iall.run.endtime = endtime;
db_iall.run.timezone = tz;
if range_bool;
    db_iall.run.profile_input = 'ranges';
else
    db_iall.run.profile_input = profile_input;
end
db_iall.run.profnums = profnums;
db_iall.run.campaign_name = campaign_name;
db_iall.run.no2field = no2field;
db_iall.run.aerfield = aerfield;
db_iall.run.altfield = altfield;
db_iall.run.radarfield = radarfield;
db_iall.run.presfield = presfield;
db_iall.run.tempfield = tempfield;
db_iall.run.cloud_product = cloud_product;
db_iall.run.cloud_frac_max = cloud_frac_max;
db_iall.run.row_anomaly = row_anomaly;
db_iall.run.min_height = min_height;
db_iall.run.numBLpoints = numBLpoints;
db_iall.run.minRadarAlt = minRadarAlt;
db_iall.run.useground = useground;


if nargout == 0;
    putvar(db_iall, lon_iall, lat_iall, omino2_iall, behrno2_iall, airno2_iall, dates_iall);
else
    varargout{1} = lon_iall;
    varargout{2} = lat_iall;
    varargout{3} = omino2_iall;
    varargout{4} = behrno2_iall;
    varargout{5} = airno2_iall;
    varargout{6} = db_iall;
    varargout{7} = dates_iall;
end

    function val = set_input_value(user_value, parameter_name)
        if ~strcmpi(pout.(parameter_name), noinput)
            val = pout.(parameter_name);
        else
            val = user_value;
        end
    end

end
