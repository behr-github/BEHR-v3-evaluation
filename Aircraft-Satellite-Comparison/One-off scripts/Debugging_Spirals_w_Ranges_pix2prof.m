%Debugging_Spirals
%
%   Sequentially creates maps with BEHR NO2 maps, overlaid with the flight
%   path for that day and the pixel boundaries, plus a second figure where
%   the pixel centers are marked with the difference in columns

date_start = '05/04/2006';
date_end = '05/31/2006';
no2field = 'NO2';   %For Baltimore, this is NO2_LIF or NO2_NCAR
%For CA/TX, this is NO2_MixingRatio_LIF or NO2_MixingRatio
altfield = 'ALTITUDE_GPS';
radarfield = 'ALTITUDE_RADAR';
tempfield = 'TEMP_STAT_C';
presfield = 'STAT_PRESSURE';

tz = 'auto';

states={'ak'};

% Which plots to show
avgNO2 = true;
deltaNO2 = true; 
    deltatype = 'omi'; %set to "omi" or "behr"
    plotflightseg = true; % setting to true will plot the segments of the plane's flight that are used to get the pixel measurements
dailySatNO2 = false;
stratNO2 = false;

merge_dir = '/Volumes/share/GROUP/INTEX-B/Matlab files/';
behr_dir = '/Volumes/share-sat/SAT/OMI/Gridded_SP_Files/';
behr_prefix = 'OMI_SP_griddedNW_';

behr_map_file = '/Users/Josh/Documents/MATLAB/Figures/Sat Verification/INTEX B/OMI SP NO2 over NW May 06.mat';
range_file = '/Volumes/share/GROUP/INTEX-B/INTEXB_Profile_UTC_Ranges_Inclusive.mat';

DEBUG_LEVEL = 2;

load(range_file); range_dates = cellstr(datestr({Ranges.Date},29));
dates = datenum(date_start):datenum(date_end);

load('blue_red_cmap.mat'); load('coast'); coastlat = lat; coastlon = long;
load(behr_map_file);
latbdy = [min(LAT_GRID(:)), max(LAT_GRID(:))];
lonbdy = [min(LON_GRID(:)), max(LON_GRID(:))];
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
    
    lonall = cat(1,lon{:}); latall = cat(1,lat{:});
    if isempty(lonall)
        continue
    else
        
        % The satellite map
        if avgNO2
            fmap = figure; m_proj('Mercator','long',lonbdy,'lat',latbdy);
            m_pcolor(LON_GRID, LAT_GRID, NO2_GRID); shading flat
            m_states('w'); colorbar; caxis([0, 1e16]);
            m_coast('color','w');
            m_grid('linestyle','none');
            
            for a=1:numel(db)
                if ~isempty(airno2{a})
                    s = numel(db(a).loncorn);
                    for b=1:s
                        if ~isnan(airno2{a})
                            loncorn = db(a).loncorn{b}; latcorn = db(a).latcorn{b};
                            for c=1:size(loncorn,2)
                                m_line([loncorn(:,c); loncorn(1,c)], [latcorn(:,c); latcorn(1,c)],'color','w','linewidth',2)
                            end
                        end
                    end
                end
            end
            
            m_line(Merge.Data.LONGITUDE.Values-360, Merge.Data.LATITUDE.Values, 'color', 'r', 'linewidth',1);
            title(sprintf('%s',Merge.metadata.date),'fontsize',20)
            oldposmap = get(fmap,'Position'); set(fmap,'Position',[1,1, oldposmap(3:4)*1.5]);
        end
        
        % The difference map
        if deltaNO2
            fcent = state_outlines(states{:});
            line(coastlon, coastlat);
            if strcmpi(deltatype,'behr'); delta = cat(1,behrno2{:}) - cat(1,airno2{:});
            else delta = cat(1,omino2{:}) - cat(1,airno2{:});
            end
            for a=1:numel(db)
                if ~isempty(airno2{a}) 
                    s = numel(db(a).loncorn);
                    for b=1:s
                        if ~isnan(airno2{a})
                            loncorn = db(a).loncorn{b}; latcorn = db(a).latcorn{b};
                            for c=1:size(loncorn,2)
                                line([loncorn(:,c); loncorn(1,c)], [latcorn(:,c); latcorn(1,c)],'color','b','linewidth',2)
                            end
                        end
                    end
                end
            end
            hold on;
            % Plot the part of the plane's flight that contributed to
            % measurements
            rs = cat(1,db.profnums); rs = cat(1,rs{:});
            for r=1:size(rs,1)
                myutc = Merge.Data.UTC.Values > rs(r,1) & Merge.Data.UTC.Values < rs(r,2);
                mylat = Merge.Data.LATITUDE.Values(myutc); mylon = Merge.Data.LONGITUDE.Values(myutc);
                latlon_fills = mylat == Merge.Data.LATITUDE.Fill | mylon == Merge.Data.LONGITUDE.Fill;
                mylat = mylat(~latlon_fills); mylon = mylon(~latlon_fills); mylon(mylon>180) = mylon(mylon>180) - 360;
                line(mylon, mylat, 'color',[0 0.6 0],'linewidth',3);
            end
            % Make the plot of pixel/profile differences
            scatter(lonall, latall, 20, delta);
            dx = abs(lonbdy(2)-lonbdy(1))/8;
            for a=1:numel(delta)
                % Find if any previous texts are at the same place
                if a>1
                    samept = sum(abs(lonall(1:a-1)-lonall(a)) < 0.1 & abs(latall(1:a-1)-latall(a)) < 0.1);
                else
                    samept = 0;
                end
                text(lonall(a)+dx/2+samept*dx, latall(a), sprintf('%+.1e',delta(a)),'BackgroundColor',[0.6 0.6 0.6]);
            end
            %xlim(lonbdy); ylim(latbdy);
            colorbar; colormap(blue_red_cmap); caxis([-5e15 5e15]);
            title(sprintf('%s',Merge.metadata.date),'fontsize',20)
            oldposcent = get(fcent,'Position');
            if avgNO2
                set(fcent,'Position',[oldposmap(3)*2,1,oldposcent(3:4)*1.5]);
            end
            hold off
        end
        
        %A satellite map for that day
        if dailySatNO2
            cb = no2_column_map_2014(curr_date, curr_date, lonbdy, latbdy, 'rowanomaly', 'AlwaysByRow', 'projection', 'Mercator',...
                'behrdir',behr_dir,'fileprefix','OMI_BEHR_omiCloudAMF_','cbrange',[0 1e16]);
            fsatday = gcf;
        end
        
        %Plot the stratospheric columns
        if stratNO2
            fstrat = state_outlines(states{:});
            SNO2 = [];
            for a=1:numel(db)
                s = size(db(a).loncorn);
                if s(2) > 0;
                    SNO2 = cat(1,SNO2,db(a).db.strat_NO2);
                    for b=1:s(2)
                        line([db(a).loncorn(:,b); db(a).loncorn(1,b)], [db(a).latcorn(:,b); db(a).latcorn(1,b)],'color','b','linewidth',2)
                    end
                end
            end
            hold on;
            scatter(cat(1,lon{:}), cat(1,lat{:}), 20, SNO2); colorbar;
            xlim(lonbdy); ylim(latbdy);
            title(sprintf('%s Stratospheric NO2',Merge.metadata.date),'fontsize',20)
        end
        
        pause;
        if deltaNO2; close(fcent); end
        if avgNO2; close(fmap); end
        if dailySatNO2; close(fsatday); end
        if stratNO2; close(fstrat); end
    end
    
end