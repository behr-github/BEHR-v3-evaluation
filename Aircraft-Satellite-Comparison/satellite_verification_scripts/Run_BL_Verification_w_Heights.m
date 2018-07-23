%Run_BL_Verification_w_Heights
%
%   Runs the boundary_layer_verification_presel_heights function for each
%   day specified.

date_start = '06/18/2008';
date_end = '06/24/2008';

no2field = 'NO2_UCB';
altfield = 'GPS_Altitude';
presfield = 'PRESSURE';

starttime = '12:00';
endtime = '15:00';

tz = 'auto';

merge_prefix = 'ARCTAS_CA_';
merge_dir = '/Volumes/share/GROUP/ARCTAS/Matlab Files/';

behr_prefix = 'OMI_BEHR_omi*';
behr_dir = '/Volumes/share-sat/SAT/BEHR/ARCTAS_BEHR/';
heights_file = '/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/ARCTAS-CA BL Heights Exclusive 3 Redone Edited 24 Sept 2014.mat';

DEBUG_LEVEL = 2;

load(heights_file); %Loads the variable "Heights"
dates = datenum(date_start):datenum(date_end);
if all(datenum({Heights(:).Date})<min(dates)) || all(datenum({Heights(:).Date})>max(dates));
    warning('All dates in the specified heights file fall outside the date range');
end

S=0; clear('db');
for d=1:numel(dates)
    % Load the merge and BEHR files
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    merge_filename = sprintf('%s%s_%s_%s.mat',merge_prefix,year,month,day);
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

    
    for swath=1:numel(Data)
        S=S+1;
        [lon_o{S}, lat_o{S}, omino2_o{S}, behrno2_o{S}, airno2_o{S}, db(S)] = ...
            boundary_layer_verification_presel_heights(Merge,Data(swath),Heights,tz,'DEBUG_LEVEL',DEBUG_LEVEL,...
            'no2field',no2field,'altfield',altfield,'presfield',presfield,'starttime',starttime,'endtime',endtime);
        date_cell{S} = repmat({curr_date},numel(lon_o{S}),1);
    end
end

% concatenate the output
[db_oall, lon_oall, lat_oall, omino2_oall, behrno2_oall, airno2_oall, dates_oall] = match_arrays2db(db, lon_o, lat_o, omino2_o, behrno2_o, airno2_o, date_cell);