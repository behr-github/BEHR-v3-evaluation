%plot_flightpaths_by_day
%
%   Draws all flightpaths for the date range specified, colored by day.

date_start = '07/01/2011';
date_end = '07/31/2011';

states = {'md','pa','va','wv'}; %States to draw on the map

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

DEBUG_LEVEL = 2;

dates = datenum(date_start):datenum(date_end);

lon_array = cell(1,30);
lat_array = cell(1,30);
merge_dates = cell(1,30);
merge_datenums = zeros(1,30);
S=0;
for d=1:numel(dates)
    % Load the merge and BEHR files
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
    
    S=S+1;
    lon_array{S} = Merge.Data.LONGITUDE.Values - 360;
    lat_array{S} = Merge.Data.LATITUDE.Values;
    merge_dates{S} = Merge.metadata.date;
    merge_datenums(S) = datenum(Merge.metadata.date);
end

state_outlines(states{:});
for a=1:S
    datecolor = merge_datenums(a)*ones(size(lon_array{a}));
    scatter(lon_array{a},lat_array{a},8,datecolor);
end
title(sprintf('LIF vs. NCAR for %s to %s',date_start,date_end),'fontsize',16,'fontweight','bold');
cb = colorbar;
set(cb,'Ticks',merge_datenums(1:S));
set(cb,'TickLabels',merge_dates(1:S));