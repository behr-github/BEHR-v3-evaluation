function [ LT_table ] = find_aircraft_profile_times( start_date, end_date, varargin )
% find_aircraft_profile_time(start_date, end_date): Returns a list of
% profile numbers and their start and end time.  By default, the time is in
% local sun time decimal notiation, but this can be changed to 'lt-hhmm',
%'utc-sec' or 'utc-hhmm' as an optional third argument.
%
%   This function will iterate through all .mat files in the folder
%   /Volumes/share/GROUP/DISCOVER-AQ/Matlab Files within a range of dates
%   given and return a table of the profiles start and end times.
%
%   To override the default path, use the parameter value 'mat_dir'

p = inputParser;
p.addRequired('startdate',@isstr);
p.addRequired('enddate',@isstr);
p.addOptional('timemode', 'lt', @(x) any(strcmpi(x,{'lt','lt-hhmm','utc-sec','utc-hhmm'})));
p.addParamValue('mat_dir', '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft', @isstr);

p.parse(start_date, end_date, varargin{:})
pout = p.Results;

start_date = pout.startdate;
end_date = pout.enddate;
timemode = pout.timemode;
mat_dir = pout.mat_dir;

if ~exist(mat_dir,'dir');
    error('find_prof_times:mat_dir_DNE','Directory to matlab files does not exist.')
end

if any(strcmpi(timemode,{'lt','lt-hhmm'})); timefield = 'LOCAL_SUN_TIME';
elseif any(strcmpi(timemode,{'utc-sec','utc-hhmm'})); timefield = 'UTC';
end

file_pattern = '*.mat';
mat_files = dir(fullfile(mat_dir, file_pattern));

LT_table = {'Prof. Num.','Date','Start Time','End Time'};
tmp_table = [];
for a=1:numel(mat_files);
    % Handle date ranges
    fname = mat_files(a).name;
    [s, e] = regexp(fname,'\d\d\d\d.*\d\d.*\d\d'); %Find the date in the file name, assuming it is of the form yyyymmdd with any or no single character separator
    fdate = datenum(fname(s:e));
    if ~isempty(start_date) && ~isempty(end_date) %If the user specified start and end dates, check that the file is within the date range
        % If the date of the file falls outside the range give, skip this
        % file.
        if fdate < datenum(start_date) || fdate > datenum(end_date)
            continue
        end
    end
    
    % Load the file
    load(fullfile(mat_dir, mat_files(a).name), 'Merge');
    
    profnums = Merge.Data.ProfileSequenceNum.Values;
    unique_profnums = unique(profnums(profnums>0));
    for b = 1:numel(unique_profnums)
        times = eval(sprintf('Merge.Data.%s.Values',timefield));
        xx = profnums == unique_profnums(b);
        starttime = min(times(xx)); endtime = max(times(xx));
        s= size(tmp_table);
        tmp_table(s(1)+1,:) = [unique_profnums(b), fdate, starttime, endtime];
    end
end

tmp_table = sortrows(tmp_table);

for b = 1:size(tmp_table,1)
    profnum = tmp_table(b,1);
    date_str = datestr(tmp_table(b,2),2);
    starttime = tmp_table(b,3);
    endtime = tmp_table(b,4);
    switch timemode
        case 'lt'
            starttime_str = num2str(starttime);
            endtime_str = num2str(endtime);
        case 'lt-hhmm'
            starttime_str = sprintf('%2d:%02d',floor(starttime),floor(mod(starttime,1)*60));
            endtime_str = sprintf('%2d:%02d',floor(endtime),floor(mod(endtime,1)*60));
        case 'utc-sec'
            starttime_str = num2str(starttime);
            endtime_str = num2str(endtime);
        case 'utc-hhmm'
            starttime_str = sprintf('%2d:%02d',floor(starttime/3600),floor(mod(starttime,3600)/60));
            endtime_str = sprintf('%2d:%02d',floor(endtime/3600),floor(mod(endtime,3600)/60));
    end
    s = size(LT_table);
    LT_table(s(1)+1,:) = {profnum, date_str, starttime_str, endtime_str};
    
end

