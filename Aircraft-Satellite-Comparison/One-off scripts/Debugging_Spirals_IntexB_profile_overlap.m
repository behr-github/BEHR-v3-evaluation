%Debugging_Spirals
%
%   Sequentially creates maps with BEHR NO2 maps, overlaid with the flight
%   path for that day and the pixel boundaries, plus a second figure where
%   the pixel centers are marked with the difference in columns

date_start = '03/16/2006';
date_end = '05/31/2006';
no2field = 'NO2';   %For Baltimore, this is NO2_LIF or NO2_NCAR
%For CA/TX, this is NO2_MixingRatio_LIF or NO2_MixingRatio
altfield = 'ALTITUDE_GPS';
radarfield = 'ALTITUDE_RADAR';
tempfield = 'TEMP_STAT_C';
presfield = 'STAT_PRESSURE';

tz = 'auto';

states={'tx'};

merge_dir = '/Volumes/share/GROUP/INTEX-B/Matlab files/';
behr_dir = '/Volumes/share-sat/SAT/OMI/Gridded_SP_Files/';
behr_prefix = 'OMI_SP';

range_file = '/Volumes/share/GROUP/INTEX-B/INTEXB_Profile_UTC_Ranges_Inclusive.mat';

DEBUG_LEVEL = 2;

load(range_file); range_dates = cellstr(datestr({Ranges.Date},29));
dates = datenum(date_start):datenum(date_end);

load('blue_red_cmap.mat'); load('coast'); coastlat = lat; coastlon = long;

for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('%s*%s%s%s.mat',behr_prefix,year,month,day);
    
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
    
    % Find the UTC range data for this date
    xx = find(strcmp(curr_date,range_dates));
    if isempty(xx);
        error('run_spiral:ranges','No UTC ranges found for %s',curr_date);
    end
    
    S=0;
    lon = cell(1,4); lat = cell(1,4); omino2 = cell(1,4); behrno2 = cell(1,4); airno2 = cell(1,4); clear('db');
    for swath=1:numel(Data)
        S=S+1;
        [lon{S}, lat{S}, omino2{S}, behrno2{S}, airno2{S}, db(S)] = spiral_verification_avg_pix2prof(Merge,Data(swath),tz,'DEBUG_LEVEL',DEBUG_LEVEL,'no2field',no2field,'profiles',Ranges(xx).Ranges,'radarfield',radarfield,'altfield',altfield,'presfield',presfield,'tempfield',tempfield);
    end
    
    f1 = figure;
    plot(coastlon, coastlat, 'color','k');
    
    for a=1:numel(airno2)
        if isempty(airno2{a}); continue; end
        for b=1:numel(airno2{a})
            figure(f1);
            if isnan(airno2{a}(b)); continue; end
            line(lon{a}(b),lat{a}(b),'color','r','marker','s');
            delta = omino2{a}(b) - airno2{a}(b); 
            text(lon{a}(b)+0.5 + b, lat{a}(b), sprintf('%+.1e',delta),'BackgroundColor',[0.7 0.7 0.7]);
            loncorn = cat(1,db(a).loncorn{b},db(a).loncorn{b}(1,:),nan(size(db(a).loncorn{b}(1,:))));
            latcorn = cat(1,db(a).latcorn{b},db(a).latcorn{b}(1,:),nan(size(db(a).latcorn{b}(1,:))));
            line(loncorn,latcorn,'marker','^');
            
            myrange = db(a).profnums{b};
            myutc = Merge.Data.UTC.Values > myrange(1) & Merge.Data.UTC.Values < myrange(2);
            mylat = Merge.Data.LATITUDE.Values(myutc); mylon = Merge.Data.LONGITUDE.Values(myutc);
            fills = mylat == Merge.Data.LONGITUDE.Fill | mylon == Merge.Data.LONGITUDE.Fill;
            mylat = mylat(~fills); mylon = mylon(~fills); mylon(mylon>180) = mylon(mylon>180)-360;
            
            line(mylon,mylat,'color',[0 0.6 0], 'linewidth',1);
            
            line(db(a).lon_3km{b}, db(a).lat_3km{b}, 'color',[0.5 0 0.5], 'linewidth',2);
            
            fprintf('\nPaused.\n');
            pause;
        end
    end
    close(f1);
    
end