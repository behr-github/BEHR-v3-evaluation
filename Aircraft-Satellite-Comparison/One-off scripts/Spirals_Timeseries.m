%Debugging_Spirals
%
%   Sequentially creates maps with BEHR NO2 maps, overlaid with the flight
%   path for that day and the pixel boundaries, plus a second figure where
%   the pixel centers are marked with the difference in columns

date_start = '09/01/2013';
date_end = '09/30/2013';
no2field = 'NO2_MixingRatio_LIF';   %For Baltimore, this is NO2_LIF or NO2_NCAR
%For CA/TX, this is NO2_MixingRatio_LIF or NO2_MixingRatio

plottype = 'stratno2'; %can be deltabyday, stratno2

tz = 'auto';

merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';
behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';

DEBUG_LEVEL = 0;

dates = datenum(date_start):datenum(date_end);

load('blue_red_cmap.mat');

D=0;
clear('file_dates'); clear('file_data');
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('OMI_BEHR_*%s%s%s.mat',year,month,day);
    
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
    
    S=0;
    lon = cell(1,4); lat = cell(1,4); omino2 = cell(1,4); behrno2 = cell(1,4); airno2 = cell(1,4); clear('db');
    for swath=1:numel(Data)
        S=S+1;
        [lon{S}, lat{S}, omino2{S}, behrno2{S}, airno2{S}, db(S)] = spiral_verification(Merge,Data(swath),tz,'DEBUG_LEVEL',DEBUG_LEVEL,'no2field',no2field);
    end
    
    D=D+1;
    %%%%% Mutable code - change this to plot whatever is desired %%%%%
    
    % Plot the difference in columns (behr - aircraft) as a scatter plot by
    % day
    switch plottype
        case 'deltabyday'
            mydata{1} = cat(1,behrno2{:}) - cat(1,airno2{:});
            mydata{2} = D*ones(size(mydata{1}));
            
            file_dates{D} = datestr(curr_date,2);
            file_data{D} = mydata;
            
        case 'stratno2'
            mydata{1} = [];
            for a=1:numel(db);
                mydata{1} = [mydata{1}; db(a).strat_NO2];
            end
            mydata{2} = D*ones(size(mydata{1}));
            
            file_dates{D} = datestr(curr_date,2);
            file_data{D} = mydata;
    end
    %%%%% End Mutable Code                                      %%%%%
    
    
    
end

%%%%% Plotting code - change as needed %%%%%

switch plottype
    case 'deltabyday'
        figure; hold on
        for a=1:numel(file_data)
            scatter(file_data{a}{2}, file_data{a}{1});
        end
        set(gca,'xticklabel',file_dates)
        xlims = get(gca,'xlim'); xlim([xlims(1)-1, xlims(2)+1]);
        fdates(2:numel(file_dates)+1) = file_dates(:); fdates{1} = ' '; fdates{end+1} = '';
        set(gca,'xticklabel',fdates);
        ylims = get(gca,'ylim'); set(gca,'ytick',ylims(1):1e15:ylims(2));
        grid on
    case 'stratno2'
        figure; hold on
        for a=1:numel(file_data)
            scatter(file_data{a}{2}, file_data{a}{1});
        end
        xlims = get(gca,'xlim'); xlim([xlims(1)-1, xlims(2)+1]);
        fdates(2:numel(file_dates)+1) = file_dates(:); fdates{1} = ''; fdates{end+1} = '';
        set(gca,'xticklabel',fdates);
        ylims = get(gca,'ylim'); set(gca,'ytick',ylims(1):1e15:ylims(2));
        grid on
end

