function [ Names, dates, directory, range_file, ground_site_dir ] = merge_field_names( campaign_name, ask_ranges, range_file_in )
%merge_field_names Returns field names of key fields in merge files
%   Different field campaigns name the same data differently.  This
%   function will return a structure with all the appropriate field names
%   for use in any other functions.  This allows this function to be a
%   central clearing house for field names that any other function can make
%   use of.
%
%   [ NAMES, DATES, DIRECTORY ] = MERGE_FIELD_NAMES( CAMPAIGN_NAME )
%   returns the field names in NAMES, range of valid DATES as a two element
%   cell array of strings, and the Merge data DIRECTORY as a string for the
%   campaign.
%
%   Valid campaign designations are:
%       arctas
%       seac3rs/seacers
%       discover-md/ca/tx/co
%   This function will try to match the input string to one of these; it is
%   not especially picky, so for example discovermd, Discover-MD, and
%   Discover-AQ MD will all successfully indicate to use the field names
%   from the Maryland Discover-AQ campaign.  If no campaign name can be
%   matched, an error is thrown.
%
%   The second argument is optional. See below for its use.
%
%   Fields returned are:
%       pressure_alt - pressure derived altitude
%       gps_alt - GPS derived altitude
%       radar_alt - radar altitude, i.e. altitude above the ground
%       theta - potential temperature measurements
%       no2_lif - Our (Cohen group) NO2 measurements
%       no2_ncar - Andy Weinheimer's NO2 measurmenets
%       aerosol_extinction - Aerosol extinction measurments at green
%           wavelengths
%       aerosol_scattering - Aerosol scattering only measurements (no
%           absorption)
%       abs_angstr_exp - Angstrom exponent for absorption
%       scat_angst_exp - Angstrom exponent for scattering
%       aerosol_ssa - Aerosol single scattering albedo measurements
%       profile_numbers - The field with the number assigned to each
%           profile. Only a field in DISCOVER campaigns.
%   Any field that doesn't exist for a given campaign returns an empty
%   string.
%
%
%   [ ___, RANGE_FILE ] = MERGE_FIELD_NAMES( ___ ) will return the path to
%   a .mat file containing UTC ranges manually identified as profiles. If
%   there are more than one, it will ask you to choose one. If none are
%   available, and empty string is returned.
%
%   [ ___, RANGE_FILE, GROUND_SITE_DIR ] = MERGE_FIELD_NAMES( ___ ) will
%   also return the directory with ground site data for the specified
%   campaign. If not available, it is returned as an empty string.
%
%   [ ___, ~, GROUND_SITE_DIR ] = MERGE_FIELD_NAMES( ___, false ) returns
%   the ground site directory and under no circumstances asks which range
%   file to use.
%
%   [ ___ ] = MERGE_FIELD_NAMES( ___, 'specified' ) will return the first
%   range file, if one is available. 
%
%   [ ___ ] = MERGE_FIELD_NAMES( ___, 'specified', N ) will return the Nth
%   range file, if one is available. If N > number of range files
%   available, and error is thrown.
%
%   [ ___ ] = MERGE_FIELD_NAMES( ___, 'specified', REGEX ) will return the
%   file matching the regular expression REGEX. If none or > 1 file match,
%   an error is thrown.

E = JLLErrors;

% Check that the first input is a string. Default the second argument to
% true, and ensure that it is a logical variable or can be used as one.
if ~ischar(campaign_name)
    E.badinput('campaign_name must be a string')
end
if nargin < 2
    ask_ranges = true;
elseif (ischar(ask_ranges) && strcmpi(ask_ranges, 'specified')) || (islogical(ask_ranges) && ~ask_ranges)
    ask_ranges = false;
    if nargin < 3
        range_file_in = 1; 
    elseif ~(isscalar(range_file_in) && isnumeric(range_file_in) && range_file_in >= 1) && ~ischar(range_file_in)
        E.badinput('When specifying a range file, it must either be a scalar number >= 1 (to specify the index in merge_field_names) or the name of a file')
    end
elseif ~islogical(ask_ranges) && (~isnumeric(ask_ranges) || ~isscalar(ask_ranges))
    E.badinput('ask_ranges (if given) must be a logical or scalar numeric value')
end


% Setup the fields the output structure is expected to have - this will be
% validated against before the structure is returned.  That way, if I make
% a mistake adding a new campaign we can avoid instances of other functions
% expecting a field to be no2_lif and getting one that is NO2_LIF.  ADD ANY
% ADDITIONAL FIELDS TO RETURN HERE.

return_fields = {'utc','longitude','latitude','pressure_alt', 'gps_alt', 'radar_alt', 'temperature', 'pressure','theta', 'h2o', 'no2_lif', 'no2_ncar','acn','hcn','co'...
    'aerosol_extinction', 'aerosol_extinction_green', 'aerosol_scattering', 'aerosol_scattering_green', 'abs_angstr_exp','scat_angstr_exp','aerosol_ssa',...
    'aerosol_dry_ssa', 'profile_numbers','ground_no2','ground_utc'}';

% Initialize the return variables
for a=1:numel(return_fields)
    Names.(return_fields{a}) = '';
end

dates = cell(1,2);
directory = '';
ground_site_dir = '';
range_files = {''};

% All campaign data should be stored in a central directory, this is that
% directory
main_dir = fullfile('/Volumes','share2','USERS','LaughnerJ','CampaignMergeMats');

% Parse the campaign name and assign the fields

% DISCOVER-MD
if ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'md'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'A_RadarAlt';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_LIF';
    Names.no2_ncar = 'NO2_NCAR';
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXTamb532';
    Names.aerosol_scattering_green = 'SCamb532';
    Names.abs_angstr_exp = 'Angstrom_Exponent_of_Absorption_at_450and550nm';
    Names.scat_angstr_exp = 'Angstrom_Exponent_of_Scattering_at_450and550nm';
    Names.aerosol_ssa = 'SingleScatteringAlbedo_at_550nmambient';
    Names.aerosol_dry_ssa = 'SingleScatteringAlbedo_at_550nm';
    Names.profile_numbers = 'ProfileSequenceNum';
    Names.ground_no2 = {'','f_42c_NO2','NO2_conc_ppb','','NO2','NO2'};
    Names.ground_utc = {'','UTC','UTC','','UTC','UTC';...
        '','UTC_mid','Mid_UTC','','Mid_UTC','UTC_mid';...
        '','UTC_stop','Stop_UTC','','Stop_UTC','UTC_stop'};
    
    dates = {'2011-07-01','2011-07-31'};
    directory = fullfile(main_dir,'DISCOVER-AQ_MD','P3','1sec');
    ground_site_dir = fullfile(main_dir,'DISCOVER-AQ_MD','Ground','VariousTimePeriods');
    
    % DISCOVER-CA
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'ca'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'Radar_Altitude';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_MixingRatio_LIF';
    Names.no2_ncar = 'NO2_MixingRatio';
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXTamb532_TSI_PSAP';
    Names.aerosol_scattering_green = 'SCATamb532_TSI';
    Names.abs_angstr_exp = 'AngstromExponenetABS_470to532';
    Names.scat_angstr_exp = 'AngstromExponenetSCAT_450to550';
    Names.aerosol_ssa = 'SingleScatAlbedo550amb_TSIneph_PSAP';
    Names.aerosol_dry_ssa = 'SingleScatAlbedo550dry_TSIneph_PSAP';
    Names.profile_numbers = 'ProfileNumber';
    Names.ground_no2 = {'photo_NO2_ppbv','NO2_ppbv','NO2','','','NO2'};
    Names.ground_utc = {'start_secUTC','start_secUTC','UTC_start','','','UTC_start';...
        'mid_secUTC','mid_secUTC','UTC_mid','','','UTC_mid';...
        'stop_secUTC','stop_secUTC','UTC_stop','','','UTC_stop'};
    
    dates = {'2013-01-16','2013-02-06'};
    directory = fullfile(main_dir, 'DISCOVER-AQ_CA','P3','1sec');
    ground_site_dir = fullfile(main_dir,'DISCOVER-AQ_CA','Ground','VariousTimePeriods');
    
    % DISCOVER-TX
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'tx'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'Radar_Altitude';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.h2o = 'H2O_MixingRatio';
    Names.no2_lif = 'NO2_MixingRatio_LIF';
    Names.no2_ncar = 'NO2_MixingRatio';
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXT532nmamb_total_LARGE';
    Names.aerosol_scattering_green = 'SCAT550nmamb_total_LARGE';
    Names.abs_angstr_exp = 'AE_ABS_450to700nm_LARGE';
    Names.scat_angstr_exp = 'AE_SCAT_450to700nm_LARGE';
    Names.aerosol_ssa = 'SSA550nmamb_LARGE';
    Names.aerosol_dry_ssa = 'SSA550nmdry_LARGE';
    Names.profile_numbers = 'ProfileNumber';
    Names.ground_no2 = {'NO2ppbv', 'NO2_ppbv', 'CAPS_NO2_ppbv', 'NO2_099', 'photo_NO2_ppbv', 'NO2_099', '', 'NO2_099'};
    Names.ground_utc = {'UTC_start', 'StartTime_UTsec', 'start_sec-UTC', 'StartTime', 'start_sec-UTC', 'StartTime', '', 'StartTime';...
        'UTC_mid', 'MidTime_UTsec', 'mid_sec-UTC', '', 'mid_sec-UTC', '', '', '';...
        'UTC_stop', 'StopTime_UTsec', 'stop_sec-UTC', 'StopTime', 'stop_sec-UTC', 'StopTime', '', 'StopTime'};
    
    dates = {'2013-09-01','2013-09-30'};
    directory = fullfile(main_dir, 'DISCOVER-AQ_TX','P3','1sec');
    ground_site_dir = fullfile(main_dir,'DISCOVER-AQ_TX','Ground','VariousTimePeriods');
    
    % DISCOVER-CO
elseif ~isempty(regexpi(campaign_name,'discover')) && ~isempty(regexpi(campaign_name,'co'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'Radar_Altitude';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_LIF';
    Names.no2_ncar = 'NO2_MixingRatio';
    Names.acn = 'Acetonitrile_MixingRatio';
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXT532nmamb_total_LARGE';
    Names.aerosol_scattering_green = 'SCAT550nmamb_total_LARGE';
    Names.abs_angstr_exp = 'AE_ABSdry_450to700nm_LARGE';
    Names.scat_angstr_exp = 'AE_SCATamb_450to700nm_LARGE';
    Names.aerosol_ssa = 'SSA550nmamb_LARGE';
    Names.aerosol_dry_ssa = 'SSA550nmdry_LARGE';
    Names.profile_numbers = 'ProfileNumber';
    Names.ground_no2 = {'Photo_NO2_ppbv', 'NO2_ppbv', 'T500_NO2_ppbv', '', 'NO2', 'Photo_NO2_ppbv'};
    Names.ground_utc = {'start_sec-UTC', 'Start_UTC', 'start_sec-UTC', '', 'UTC_start', 'start_sec-UTC';...
        'mid_sec-UTC', 'Mid_UTC', 'mid_sec-UTC', '', 'UTC_mid', 'mid_sec-UTC';...
        'stop_sec-UTC', 'Stop_UTC', 'stop_sec-UTC', '', 'UTC_stop', 'stop_sec-UTC'};
    
    dates = {'2014-07-17','2014-08-10'};
    directory = fullfile(main_dir, 'DISCOVER-AQ_CO','P3','1sec/');
    ground_site_dir = fullfile(main_dir,'DISCOVER-AQ_CO','Ground','VariousTimePeriods');
    
    % SEAC4RS
elseif ~isempty(regexpi(campaign_name,'seac4rs')) || ~isempty(regexpi(campaign_name,'seacers'));
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'RadarAlt';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_TDLIF';
    Names.no2_ncar = 'NO2_ESRL'; % This is Ryerson's NO2, not sure if that's different from Weinheimer's
    Names.co = 'CO_DACOM';
    Names.acn = 'Acetonitrile'; % guessing this is Wisthaler's PTRMS
    Names.hcn = 'HCN_CIT'; % guessing this is Wennberg's CIT-CIMS
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXT532nmamb_total_LARGE';
    Names.aerosol_scattering_green = 'SCAT550nmamb_total_LARGE';
    Names.aerosol_ssa = 'SSA550nmamb_TSIandPSAP_LARGE';
    Names.aerosol_dry_ssa = 'SSA550nmdry_TSIandPSAP_LARGE';
    
    dates = {'2013-08-06','2013-09-23'};
    directory = fullfile(main_dir, 'SEAC4RS','DC8','1sec');
    
    range_files = {fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_Std.mat'),...
        fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_Porpoising.mat'),...
        fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_Curtaining.mat'),...
        fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_StdFires.mat'),...
        fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_CurtainingFires.mat'),...
        fullfile(main_dir, 'SEAC4RS','SEAC4RS_Profile_Ranges_Std-plus-Curtaining.mat')};
    
    % DC3 (not to be confused with the DC8 aircraft)
elseif ~isempty(regexpi(campaign_name,'dc3'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_ALT';
    Names.radar_alt = 'RadarAlt';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_TDLIF';
    Names.no2_ncar = 'NO2_ESRL'; % This is Ryerson's NO2, not sure if that's different from Weinheimer's
    Names.co = 'CO_DACOM';
    Names.acn = 'Acetonitrile_PTRMS'; % guessing this is Wisthaler's PTRMS
    Names.hcn = 'HCN_CIT';  % guessing this is Wennberg's CIT-CIMS
    Names.aerosol_extinction = 'Lee_ext450nm_amb';
    Names.aerosol_scattering = 'Lee_sc450nm_amb';
    Names.aerosol_extinction_green = 'EXTamb532nm_TSI_PSAP';
    Names.aerosol_scattering_green = 'SCATamb532nm_TSI';
    Names.aerosol_ssa = 'SingleScatAlbedo_amb550nm_TSIneph_PSAP';
    Names.aerosol_dry_ssa = 'SingleScatAlbedo_dry700nm_TSIneph_PSAP';
    
    dates = {'2012-05-18','2012-06-22'};
    directory = fullfile(main_dir, 'DC3','DC8','1sec');
    
    range_files = {fullfile(main_dir, 'DC3', 'DC3_Profile_Ranges_Std.mat'),...
        fullfile(main_dir, 'DC3', 'DC3_Profile_Ranges_Curtaining.mat')};
    
    % ARCTAS (-B and -CARB)
elseif ~isempty(regexpi(campaign_name,'arctas'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTP';
    Names.gps_alt = 'GPS_Altitude';
    Names.radar_alt = 'Radar_Altitude';
    Names.temperature = 'TEMPERATURE';
    Names.pressure = 'PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2_UCB';
    Names.no2_ncar = 'NO2_NCAR';
    Names.co = 'Carbon_Monoxide_mixing_ratio'; % Probably Glenn Diskin's DACOM CO
    Names.acn = 'Acetonitrile_PTRMS'; % I think the ACN on SEAC4RS was the PTRMS, so being consistent
    Names.hcn = 'HCN_CIT'; % Wennberg's CIT-CIMS
    Names.aerosol_extinction = 'TotalExtinctionGreen'; % I calculated this by adding scattering and abs
    Names.aerosol_scattering = 'Total_Scatter550_nm';
    % For Arctas, it is not specified whether these are taken at ambient or
    % dry conditions.  Since other later campaigns only have dry data at
    % multiple wavelengths, I'm inclined to guess that this is actually dry
    % data.
    Names.aerosol_ssa = 'Single_Scatter_Albedo_Green';
    Names.aerosol_dry_ssa = 'Single_Scatter_Albedo_Green';
    
    
    if ~isempty(regexpi(campaign_name,'carb'))
        dates = {'2008-06-18','2008-06-24'};
        directory = fullfile(main_dir,'ARCTAS-CARB','DC8','1sec');
        range_files = {fullfile(main_dir, 'ARCTAS-CARB','ARCTAS-CA Altitude Ranges Exclusive 3.mat')};
    elseif ~isempty(regexpi(campaign_name,'b'))
        dates = {'2008-06-29','2008-07-13'};
        directory = fullfile(main_dir,'ARCTAS-B','DC8','1sec');
    elseif nargout > 1
        E.badinput('Campaign is one of the ARCTAS segments, but could not ID which one (carb or b)')
    end
    
    % INTEX-B
elseif ~isempty(regexpi(campaign_name,'intex')) && ~isempty(regexpi(campaign_name,'b'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.pressure_alt = 'ALTITUDE_PRESSURE';
    Names.gps_alt = 'ALTITUDE_GPS';
    Names.radar_alt = 'ALTITUDE_RADAR';
    Names.temperature = 'TEMP_STAT_C';
    Names.pressure = 'STAT_PRESSURE';
    Names.theta = 'THETA';
    Names.no2_lif = 'NO2';
    Names.no2_ncar = '';
    Names.aerosol_extinction = 0;
    Names.aerosol_scattering = 0;
    Names.aerosol_ssa = 0;
    Names.aerosol_dry_ssa = 0;
    
    dates = {'2006-03-04','2006-05-15'};
    directory = fullfile(main_dir, 'INTEX-B','DC8','1sec');
    
    range_files = {fullfile(main_dir, 'INTEX-B','INTEXB_Profile_UTC_Ranges.mat'),...
        fullfile(main_dir, 'INTEX-B','INTEXB_Profile_UTC_Ranges_Inclusive.mat')};
    
    % TEXAQS
elseif ~isempty(regexpi(campaign_name,'texaqs','once'))
    directory = fullfile(main_dir, 'TexAQS2000','Electra','WAS');
    dates = {'2000-08-16','2000-09-13'};
    
    % SOAS
elseif ~isempty(regexpi(campaign_name, 'soas', 'once'))
    Names.utc = 'UTC';
    Names.longitude = 'LONGITUDE';
    Names.latitude = 'LATITUDE';
    Names.gps_alt = 'GpsAlt';
    Names.pressure_alt = 'PAlt';
    Names.radar_alt = 'RAlt';
    Names.pressure = 'StaticPrs';
    Names.temperature = 'AmbTemp';
    Names.theta = 'PotTemp';
    Names.no2_ncar = 'NO2_ppbv';
    
    directory = fullfile(main_dir, 'SOAS', 'P3', '1sec');
    dates = {'2013-05-31', '2013-07-10'};
    range_files = {fullfile(main_dir, 'SOAS', 'SOAS_Profile_UTC_Ranges.mat')};
else
    error(E.badinput('Could not parse the given campaign name - see help for this function for suggestions of proper campaign names.'));
end

% If there is only one range file, return it. If there are multiple
% options, ask the user which one to use (and don't continue until the
% input is valid).  If there is not a range file specified, return an empty
% string.  Only do this if the range file is actually output to a variable,
% otherwise, don't waste the user's time.

if nargout >= 4
    n = numel(range_files);
    if n == 1
        range_file = range_files{1};
    elseif n > 1
        if ask_ranges
            opts_str = cell(1,n);
            for a=1:n
                [~,opts_str{a}] = fileparts(range_files{a});
            end
            user_choice = ask_multichoice('Choose a range file', opts_str, 'list', true, 'index', true);
            range_file = range_files{user_choice};
        elseif ischar(range_file_in)
            % Try to match the filename given to one of the listed file
            % names
            xx = ~iscellcontents(regexp(range_files, range_file_in), 'isempty');
            if sum(xx) == 1
                range_file = range_files{xx};
            else
                E.callError('no_range_found', 'Found %d range files matching "%s" for %s:\n\t%s', sum(xx), range_file_in, campaign_name, strjoin(range_files(xx), '\n\t'));
            end
        elseif isnumeric(range_file_in)
            if range_file_in <= n
                range_file = range_files{range_file_in};
            else
                E.callError('bad_range_index', 'The given range file index (%d) exceeds the number of files available for %s (%d)', range_file_in, campaign_name, n);
            end
        else
            E.notimplemented('Unexpected type or value of range_file_in');
        end
    else
        range_file = '';
    end
else
    range_file = '';
end

% Check that all the fields of the output structure are what we expect
fields = fieldnames(Names);
if numel(fields) ~= numel(return_fields) || ~all(strcmp(fields,return_fields))
    error(E.callError('internal:fields_mismatch','Fields of output structure are not what is expected. Make sure any new fields are spelled correctly and that they have been added to ''return_fields'''));
end


end

