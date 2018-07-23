%plot_all_INTEX_B_flights

date_start = '03/04/2006';
date_end = '05/15/2006';

latbdy = [10, 30];
lonbdy = [-120, -80];

colorbydate = false;

no2field = 'NO2';
altfield = 'ALTITUDE_GPS';
radarfield = 'ALTITUDE_RADAR';
tempfield = 'TEMP_STAT_C';
presfield = 'STAT_PRESSURE';

tz = 'auto'; %set to 'auto' to calculate the time zone based on the mean longitude of the flight

merge_dir = '/Volumes/share/GROUP/INTEX-B/Matlab files/';
range_file = '/Volumes/share/GROUP/INTEX-B/INTEXB_Profile_UTC_Ranges.mat';

DEBUG_LEVEL = 2;

load(range_file); range_dates = {Ranges.Date};
dates = datenum(date_start):datenum(date_end);
figure; worldmap(latbdy,lonbdy);
load coast;
plotm(lat,long,'color','k')
for d=1:numel(dates)
    % Load the merge files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    % Reduce the size of the plotted vectors
    lat = Merge.Data.LATITUDE.Values; lon = Merge.Data.LONGITUDE.Values-360; 
    tmp = true(size(lat));
    tmp(1:16:end) = 0;
    
    lat(tmp) = [];
    lon(tmp) = [];
    
    if colorbydate; scatterm(lat,lon,8,dates(d)*ones(size(lat)));
    else scatterm(lat,lon,8,'k');
    end
end

if colorbydate
    cb = colorbar;
    xt = get(cb,'Ticks');
    yt = datestr(xt,2);
    set(cb,'TickLabels',yt);
end