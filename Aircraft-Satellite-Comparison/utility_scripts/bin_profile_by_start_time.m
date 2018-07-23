function [ BinStruct ] = bin_profile_by_start_time( campaign_name, start_date, end_date, varargin )
%bin_profile_by_start_time Bin DISCOVER-AQ profiles by start time
%   This function will create a structure containing profile data from the
%   DISCOVER-AQ campaigns binned to hour-wide time bins. The output
%   structure will have the form /Site#/StartTime. In addition to the
%   various sites' profiles, this structure will have an "AllData" field
%   that uses all the NO2 data from the entire file, not just the profiles.
%
%   For the individual sites, this will bin the profiles by their start
%   time, so a profile that begins at 10:59 AM will be included in the
%   10-11 bin.  For the AllData field, the data will simply be binned by
%   the time of acquisition. 
%
%   This requires three inputs: the campaign name that will be recognized
%   by merge_field_names, the first day to include in the average, and the
%   last day to include in the average (these should be given as standard
%   MATLAB date strings).
%
%   There are two optional parameters: start and end hour. Use these to
%   limit what hours will be binned, e.g. 8 and 15 will only produce the
%   bins 8-9, 9-10, 10-11, 12-1p, 1-2p, 2-3p, and 3-4p.
%
%   There are also two additional parameters that can be passed:
%   'topextrap' and 'bottomextrap'. These can take the values of the
%   parameters 'top' and 'bottom' respectively from extrapolate_profile. By
%   default, both will be set to 'none'. These parameters determine how the
%   site-by-site profiles will be extrapolated to the surface and the
%   tropopause. See the documentation from extrapolate_profile for more
%   information on the options.
%
%   Finally, note that this will only work for DISCOVER campaigns - other
%   campaigns will be rejected because they do not have profiles
%   identified, nor are the profiles present well constrained
%   geographically.

E = JLLErrors;
DEBUG_LEVEL = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(campaign_name)
    E.badinput('campaign_name must be a string')
elseif isempty(regexpi(campaign_name,'discover'))
    E.badinput('Only DISCOVER campaigns are valid for this function')
end

[Names, campaign_dates, merge_dir, merge_ground_dir] = merge_field_names(campaign_name);

% Check that the dates given are both in the right order and that they are
% within the dates for the given campaign
dstart = datenum(start_date);
dend = datenum(end_date);
dcamp = datenum(campaign_dates);
if dstart > dend
    E.badinput('start_date must be before end_date')
elseif dstart < dcamp(1) || dend > dcamp(2)
    E.badinput('Start date %s is not a valid date for the campaign %s',start_date, campaign_name)
elseif dend < dcamp(1) || dend > dcamp(2)
    E.badinput('End date %s is not a valid date for the campaign %s',end_date, campaign_name)
end

% Parse the arguments using an input parser instance
p = inputParser;
p.addOptional('start_hour',0,@(x) isnumeric(x) && isscalar(x) && x >=0 && x <= 23);
p.addOptional('end_hour',23,@(x) isnumeric(x) && isscalar(x) && x >=0 && x <= 23);
p.addParameter('topextrap', 'none', @(x) ismember(x,{'median','fit','wrf','wrf-scaled','none'}));
p.addParameter('bottomextrap', 'none', @(x) ismember(x, {'median','fit','ground','none'}));
p.parse(varargin{:});
pout = p.Results;
start_hour = pout.start_hour;
end_hour = pout.end_hour;
topextrap = pout.topextrap;
bottomextrap = pout.bottomextrap;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare a cell array to take the binned data, also go ahead and save the
% start and end times (in local standard time) of the bins;
bin_data = cell(1,24);
bin_start_times = start_hour:end_hour;
bin_end_times = bin_start_times + 1;

% Get the list of aircraft files. Find out which ones fall in the date
% range set.
Files = dir(fullfile(merge_dir,'*.mat'));
[date_ind(1), date_ind(2)] = regexp(Files(1).name, '\d\d\d\d[_-]?\d\d[_-]?\d\d'); % finds the date whether it is yyyymmdd, yyyy_mm_dd, or yyyy-mm-dd
dates = zeros(size(Files));
for a=1:numel(Files)
    file_date_str = regexprep(Files(a).name(date_ind(1):date_ind(2)),'[-_]','');
    dates(a) = datenum(file_date_str,'yyyymmdd');
end
xx = dates >= dstart & dates <= dend;
dates = dates(xx);
Files = Files(xx);

if isempty(Files)
    BinStruct = [];
    return
end

% Load the variable "Merge" from each file, then add its data to the bin.

BinStruct = struct;
BinStruct.AllSites.no2_time_bins = bin_data;
BinStruct.AllSites.pres_time_bins = bin_data;
BinStruct.AllSites.bin_start_hours = bin_start_times;
BinStruct.AllSites.bin_end_hours = bin_end_times;
BinStruct.AllSites.Lon = NaN;
BinStruct.AllSites.Lat = NaN;

for a=1:numel(Files)
    load(fullfile(merge_dir,Files(a).name), 'Merge');
    % Identify all the sites in the file. Create a field in the structure
    % if necessary
    pns = remove_merge_fills(Merge, Names.profile_numbers);
    pns = unique(pns(pns > 0));
    sitenums = unique(calc_sitenums(pns));
    
    for s=1:numel(sitenums)
        fldname = sprintf('Site%02d',sitenums(s));
        if ~isfield(BinStruct,fldname)
            BinStruct.(fldname).no2_time_bins = bin_data;
            BinStruct.(fldname).pres_time_bins = bin_data;
            BinStruct.(fldname).bin_start_hours = bin_start_times;
            BinStruct.(fldname).bin_end_hours = bin_end_times;
            BinStruct.(fldname).Lon = NaN;
            BinStruct.(fldname).Lat = NaN;
        end
        
        BinStruct.(fldname) = add2bin_site(BinStruct.(fldname), sitenums(s), Merge, Names);
    end
    BinStruct.AllSites = add2bin_all(BinStruct.AllSites, Merge, Names);
end

% Just to clean up the structure
BinStruct = orderfields(BinStruct);

% Need the month for WRF extrapolation. We'll take the most common month
% out of all the files used.
mcmonth = mode(month(dates));

% Finally, go through every field and bin each hour into OMI standard
% pressure bins, then add these back into the structure as a matrix where
% each column is an hour. Remove hours with no data first.
fns = fieldnames(BinStruct);
for a=1:numel(fns)
    % Note that the AllSites field will do WRF extrapolation based on the
    % average lat/lon of the entire campaign, so exactly how representative
    % the WRF profile will be is uncertain.  Given that it's the free
    % troposphere, it should be a reasonable approximation over the area a
    % DISCOVER campaign covered.
    BinStruct.(fns{a}) = bin_by_pressure(BinStruct.(fns{a}), mcmonth, topextrap, bottomextrap, DEBUG_LEVEL);
end

end

function sitenums = calc_sitenums(pns)
    % Some profile numbers are snnn, where s is the site number and nnn is
    % the sequential number of the profile over that site.  Others are
    % 10snnn, so this handles both.
    sitenums = pns - mod(pns,1000);
    sitenums = mod(sitenums,10000) / 1000;
end

function SiteStruct = add2bin_site(SiteStruct, site, Merge, Names)
% Load the lon (to calculate timezones), profile numbers, LIF NO2
% and GPS altitude
[no2,utc,~,lon,lat] = remove_merge_fills(Merge,Names.no2_lif);
no2(no2<0) = nan;
tz = round(lon/15);

pres = remove_merge_fills(Merge,Names.pressure);
pns = remove_merge_fills(Merge,Names.profile_numbers);
sitenums = calc_sitenums(pns);

% Find each profile that belongs to this site, get its start time,
% and use that to add it to the bin.
ss = sitenums == site;
site_pns = unique(pns(ss));
for a=1:numel(site_pns)
    xx = pns == site_pns(a);
    prof_start_time_utc = min(utc(xx));
    prof_tz = mode(tz(xx));
    prof_start_time = utc2local_sec(prof_start_time_utc, prof_tz)/3600; % profile start time, in hours local standard time
    tt = prof_start_time >= SiteStruct.bin_start_hours & prof_start_time < SiteStruct.bin_end_hours;
    if sum(tt) > 0
        SiteStruct.no2_time_bins{tt} = cat(2,SiteStruct.no2_time_bins{tt},no2(xx));
        SiteStruct.pres_time_bins{tt} = cat(2,SiteStruct.pres_time_bins{tt},pres(xx));
    end
end

SiteStruct.Lon = nanmean([nanmean(lon), SiteStruct.Lon]);
SiteStruct.Lat = nanmean([nanmean(lat), SiteStruct.Lat]);
end

function SiteStruct = add2bin_all(SiteStruct, Merge, Names)
[no2,utc,~,lon,lat] = remove_merge_fills(Merge,Names.no2_lif);
no2(no2<0)=nan;
tz = round(lon/15);
pres = remove_merge_fills(Merge,Names.pressure);

lst = utc2local_sec(utc,tz)/3600;

for a=1:numel(SiteStruct.bin_start_hours)
    xx = lst >= SiteStruct.bin_start_hours(a) & lst < SiteStruct.bin_end_hours(a);
    SiteStruct.no2_time_bins{a} = cat(2,SiteStruct.no2_time_bins{a},no2(xx));
    SiteStruct.pres_time_bins{a} = cat(2,SiteStruct.pres_time_bins{a},pres(xx));
end

SiteStruct.Lon = nanmean([SiteStruct.Lon, nanmean(lon)]);
SiteStruct.Lat = nanmean([SiteStruct.Lat, nanmean(lat)]);
end

function BinStruct = bin_by_pressure(BinStruct, month_in, topextrap, bottomextrap, DEBUG_LEVEL)
% Remove empty hours
xx = ~iscellcontents(BinStruct.no2_time_bins,'isempty');
BinStruct.no2_time_bins = BinStruct.no2_time_bins(xx);
BinStruct.pres_time_bins = BinStruct.pres_time_bins(xx);
BinStruct.bin_start_hours = BinStruct.bin_start_hours(xx);
BinStruct.bin_end_hours = BinStruct.bin_end_hours(xx);

% extrapolate_profile bins by pressure levels, so there will always be 25
% bins, so we can preinitialize the matrix to hold those bins.
BinStruct.final_no2_bins = nan(25, sum(xx));
BinStruct.final_pres_bins = nan(25, sum(xx));

% Bin each remaining hour - first time through, initialize the matrix
n = numel(BinStruct.no2_time_bins);
for a=1:n
    if sum(isnan(BinStruct.no2_time_bins{a}))/numel(BinStruct.no2_time_bins{a}) < 0.8 % require that a reasonable percentage of the profile be valid data
        [bin_no2, bin_pres] = extrapolate_profile(BinStruct.no2_time_bins{a}, BinStruct.pres_time_bins{a}, 'bottom', bottomextrap, 'top', topextrap, 'lon', BinStruct.Lon, 'lat', BinStruct.Lat, 'month', month_in);
        if a==1
            BinStruct.final_no2_bins = nan(length(bin_pres), n);
            BinStruct.final_pres_bins = nan(length(bin_pres), n);
        end
        
        if ~iscolumn(bin_no2); bin_no2 = bin_no2'; end
        if ~iscolumn(bin_pres); bin_pres = bin_pres'; end
        
        BinStruct.final_no2_bins(:,a) = bin_no2;
        BinStruct.final_pres_bins(:,a) = bin_pres;
    else
        if DEBUG_LEVEL > 1; fprintf('\t\tNot binning, > 80%% of profile data is NaNs\n'); end
        % Skip if too many nans
    end
end
end

