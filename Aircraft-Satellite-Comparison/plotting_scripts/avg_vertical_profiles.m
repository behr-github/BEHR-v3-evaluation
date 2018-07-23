function [ bin_vals, bin_midpts, bin_errors ] = avg_vertical_profiles( data_field, data_type, binwidth, start_date, end_date, varargin )
%[bin_vals, bin_midpts, bin_errors] = avg_vertical_profiles(data_field, data_type, binwidth, start_date, end_date): Averages together vertical profile data over several days, restricted by arbitrary fields.
%   This function will bin all data in Merge.Data.data_field for days
%   between start_date and end_date to altitude bins of width "binwidth"
%   (in kilometers).
%
%   This reads the Merge data structures from .mat files prepared by
%   read_merge_files.m.  Data type specifies where the data came from, e.g.
%   Aircraft or Sondes.  This must match the subdirectory name in
%   /Volumes/share/GROUPS/DISCOVER-AQ/Matlab Files.
%
%   Optionally (but recommended), you can pass an additional two arguments
%   specifying criteria by which to limit what is averaged together. The
%   first of these optional arguments is a string or cell array of strings
%   specifying the field name or names (in Merge.Data) and the range of
%   values of those fields that should be allowed into the average as a
%   matrix with form [min max].
%
%   As an example, if one wished to average NO2_LIF data from aircraft
%   between July 1st and 11th 2011, to 500 m (0.5 km) bins, for all vertical
%   profiles with profile numbers between 1000 and 1999, the correct
%   function call would be:
%   avg_vertical_profiles('NO2_LIF','Aircraft',0.5,'07/01/2011','07/11/2011','ProfileSequenceNum',[1000 1999])
%
%   To do the same range of data, but instead specify only vertical
%   profiles over site flag #1 and before 11:00 local sun time, do:
%   avg_vertical_profiles('NO2_LIF','Aircraft',0.5','07/01/2011', '07/11/2011',{'discoveraqSiteFlag1sec','LOCAL_SUN_TIME'},{[1 1], [0 11])
%
%   To ease typing, two shorthand descriptions are allowed in the criteria
%   field argument.  For aircraft data, the fields 'ProfileSequenceNum' and
%   'discoveraqSiteFlag1sec' can be called as 'profnum' or 'siteflag'.
%   This is only for the criteria field, not the first argument data field.
%   Additionally, some special arguments have been added: 'starttime-lt' or
%   'start_time-lt' will use profiles (IDed by profile number) that have a
%   start time in local sun time within the range specified.
%   'starttime-utc' or 'start_time-utc' do the same, but for the UTC field
%   (in seconds).
%
%   By default, this function will look for folders matching the data type
%   string at the directory '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files'.
%   This can be overridden using the parameter value 'mat_dir'.
%
%   Also by default, this function will bin using medians and 25th/75th
%   quantiles.  To use means and std. dev. instead, pass the parameter
%   argument 'binmode','mean'

p = inputParser;
p.addRequired('data_field',@isstr);
p.addRequired('data_type',@(x) any(strcmp(x,{'Aircraft','Sondes'})));
p.addRequired('binwidth',@isscalar);
p.addRequired('start_date',@isstr);
p.addRequired('end_date',@isstr);
p.addOptional('crit_field','',@(x) (ischar(x) || iscell(x)));
p.addOptional('crit_range',[],@(x) (ismatrix(x) || iscell(x)));
p.addParamValue('mat_dir','/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files', @(x) exist(x,'dir'));
p.addParamValue('pres_field','',@isstr);
p.addParamValue('binmode','median',@(x) any(strcmpi(x,{'median','mean'})));

p.parse(data_field,data_type,binwidth,start_date,end_date,varargin{:});
pout = p.Results;

field = pout.data_field;
dtype = pout.data_type;
binwidth = pout.binwidth;
start_date = pout.start_date;
end_date = pout.end_date;
crit_field = pout.crit_field;
crit_range = pout.crit_range;
mat_dir = pout.mat_dir;
pfield = pout.pres_field;
binmode = pout.binmode;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADDITIONAL INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if datenum(end_date) < datenum(start_date);
    error('avg_vert_prof:date_range_backwards','End date is before start date.');
end

% If either of the criteria inputs are not a cell array, make them a 1-long
% cell array. This will make things easier later.
if ~iscell(crit_field); crit_field = {crit_field}; end
if ~iscell(crit_range); crit_range = {crit_range}; end

% The length of both of the criteria cell arrays must be the same.
if numel(crit_field) ~= numel(crit_range)
    error('avg_vert_prof:criteria_input_size','The criteria field and criteria range cell arrays must have the same number of elements')
end

% Parse the shorthand entries for profile number and site flag if we are
% considering aircraft data.  Allow the user to capitalize the first letter
% Different campaigns use different field names, so check which campaign 
if strcmp(dtype,'Aircraft') && datenum(start_date) >= datenum('07/01/2011') && datenum(end_date) <= datenum('07/31/2011')
    crit_field = regexprep(crit_field,'[pP]rofnum','ProfileSequenceNum');
    crit_field = regexprep(crit_field,'[sS]iteflag','discoveraqSiteFlag1sec');
    if isempty(pfield); pfield = 'PRESSURE'; end
elseif strcmp(dtype,'Aircraft') && datenum(start_date) >= datenum('01/01/2013') && datenum(end_date) <= datenum('02/28/2013')
    crit_field = regexprep(crit_field,'[pP]rofnum','ProfileNumber');
    crit_field = regexprep(crit_field,'[sS]iteflag','SiteSeqNumber');
    if isempty(pfield); pfield = 'PRESSURE'; end
elseif strcmp(dtype,'Aircraft') && datenum(start_date) >= datenum('09/01/2013') && datenum(end_date) <= datenum('09/30/2013')
    crit_field = regexprep(crit_field,'[pP]rofnum','ProfileNumber');
    crit_field = regexprep(crit_field,'[sS]iteflag','SiteSeqNumber');
    if isempty(pfield); pfield = 'PRESSURE'; end
else
    warning('Start date does not conform to any known Discover campaigns, profnum and siteflag criteria will not work.');
end

% Prepare matrices for pressure, field values, and all the criteria. These
% will be appended to in the main loop.
pressures = [];
data_vals = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        MAIN PROGRAM       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the merge files
file_pattern = '*.mat';
mat_files = dir(fullfile(mat_dir,dtype, file_pattern));

for a=1:numel(mat_files);
    % Handle date ranges
    fname = mat_files(a).name;
    [s, e] = regexp(fname,'\d\d\d\d.*\d\d.*\d\d'); %Find the date in the file name, assuming it is of the form yyyymmdd with any or no single character separator
    fdate = datenum(fname(s:e));
    if ~isempty(start_date) && ~isempty(end_date) %If the user specified start and end dates, check that the file is within the date range
        % If the date of the file falls outside the range give, skip this file.
        if fdate < datenum(start_date) || fdate > datenum(end_date)
            continue
        end
    end
    
    % Load the file
    load(fullfile(mat_dir,dtype, mat_files(a).name), 'Merge');
    
    % Find the field that gives pressure that corresponds to altitude. Only
    % do this the first time through, assume that this field is always the
    % same in all following merge files.
    if a == 1 && isempty(pfield)
        pfield_names = {}; c=1;
        f = fieldnames(Merge.Data);
        for b=1:numel(f); %Find any fields that contain the phrase 'pressure'
            if ~isempty(regexpi(f{b},'pressure'))
                pfield_names{c} = f{b};
                c = c+1;
            end
        end
        if numel(pfield_names) == 0;
            error('avg_vert_prof:pressure_field','Pressure field could not be found. Set using the ''pres_field'' parameter');
        elseif numel(pfield_names) == 1;
            pfield = pfield_names{1};
        else
            x = zeros(1,numel(pfield_names));
            for b=1:numel(pfield_names)
                x(b) = length(pfield_names{b});
                xf = find(x(:) == min(x(:)));
                if xf == 1; %We would expect the pressure field we want to be early in the structure and have the shortest name. Shortest name overrides, but in that case write a warning to the user to double check that the right field was chosen.
                    pfield = pfield_names{1};
                else
                    pfield = pfield_names{xf};
                    warning('Pressure (altitude) field set to a %s. Confirm that the correct field was used. If not, override automatic selection using the ''pres_field'' parameter.',pfield);
                end
            end
        end
    end
    
    % Load in the data 
    x = true(1,eval(sprintf('numel(Merge.Data.%s.Values)',field))); % This will be our logical index matrix. We will successively restrict it for each criteria parameter.
    for b=1:numel(crit_field)
        if any(strcmpi(crit_field{b},{'starttime-lt','start_time-lt'}))
            unique_profnums = unique(Merge.Data.ProfileSequenceNum.Values);
            unique_profnums = unique_profnums(unique_profnums > 0);
            tmp = false(1,eval(sprintf('numel(Merge.Data.%s.Values)',field)));
            for c = 1:numel(unique_profnums)
                times = Merge.Data.LOCAL_SUN_TIME.Values(Merge.Data.ProfileSequenceNum.Values == unique_profnums(c));
                if min(times(:)) > min(crit_range{b}) && min(times(:)) < max(crit_range{b})
                    tmp = tmp | Merge.Data.ProfileSequenceNum.Values == unique_profnums(c);
                end
            end
        elseif any(strcmpi(crit_field{b},{'starttime-utc','start_time-utc'}));
            unique_profnums = unique(Merge.Data.ProfileSequenceNum.Values);
            unique_profnums = unique_profnums(unique_profnums > 0);
            tmp = false(1,eval(sprintf('numel(Merge.Data.%s.Values)',field)));
            for c = 1:numel(unique_profnums)
                times = Merge.Data.UTC.Values(Merge.Data.ProfileSequenceNum.Values == unique_profnums(c));
                if min(times(:)) > min(crit_range{b}) && min(times(:)) < max(crit_range{b})
                    tmp = tmp | Merge.Data.ProfileSequenceNum.Values == unique_profnums(c);
                end
            end
        else
            tmp = eval(sprintf('Merge.Data.%s.Values >= %f & Merge.Data.%s.Values <= %f',crit_field{b},min(crit_range{b}),crit_field{b},max(crit_range{b})));
            x = x & tmp;
        end
    end
    tmp_data = eval(sprintf('Merge.Data.%s.Values(x)',field));
    tmp_pres = eval(sprintf('Merge.Data.%s.Values(x)',pfield));
    
    % Get the fill values for the data a pressure fields, plus the ULOD and
    % LLOD markers
    data_fill = eval(sprintf('Merge.Data.%s.Fill',field));
    pres_fill = eval(sprintf('Merge.Data.%s.Fill',pfield));
    ulod = Merge.metadata.upper_lod_flag;
    llod = Merge.metadata.lower_lod_flag;
    
    % Remove any fill or LOD flags
    x2 = tmp_data == data_fill | tmp_pres == pres_fill | tmp_data == ulod | tmp_data == llod | tmp_pres == ulod | tmp_pres == llod;
    tmp_data(x2) = [];
    tmp_pres(x2) = [];
    
    data_vals = [data_vals, tmp_data];
    pressures = [pressures, tmp_pres];
end

% Plot the values converting pressure to altitude, using P0 = 1013.25 hPa
% and scale height H = 7.4 km

altitude = -log(pressures ./ 1013.25) .* 7.4;

% Get the fill

% Bin all that data
[bin_vals, bin_midpts, bin_errors] = bin_vertical_profile(altitude,data_vals,binwidth,binmode);
end
