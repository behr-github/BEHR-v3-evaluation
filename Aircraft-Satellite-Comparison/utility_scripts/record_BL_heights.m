% record_BL_heights
%
%   Script that automatically calls the select_BL_heights GUI for all
%   ranges in the specified range file (e.g. DISCOVER-Baltimore Altitude
%   Ranges.mat). It will load the correct merge file for each day as long
%   as the proper directory is given.  The resulting BL heights will be
%   saved in a data structure, Heights where each top-level index is a
%   given day and each day contains fields for BL height, median UTC times,
%   and a bit array that describes the accuracy of the BL height.
%
%   The first 9 bits of the quality flags are the only ones used currently.
%   Bits 1-3 describe the NO2 BL height, bits 4-6 the H2O BL height, and
%   bits 7-9 the potential temperature BL height.  The overall BL height is
%   an average of each individual height (NO2, H2O, Theta) unless a height
%   was rejected.  Bits 1, 4, and 7 will be set to 1 if their respective
%   species was not used.  Bits 2, 5, and 8 will be set to 1 if the BL
%   height was tweaked manually; the user can force the BL height to fall
%   within a certain range but not directly set it.  Bits 3, 6, and 9 will
%   be set if the BL height for that species is itself an average of
%   multiple heights, in particular if the remnant PBL creates a second
%   layer of lower concentration.
%
%   The heights structure will also contain the specific BL heights for
%   each species is a substructure.

merge_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignMergeMats/ARCTAS-CARB/DC8/1sec';
merge_prefix = 'ARCTAS_CA_';
ranges_file = '/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/ARCTAS-CA Altitude Ranges Exclusive 3.mat';

AltField = 'GPS_Altitude';
NO2field = 'NO2_UCB';
H2Ofield = 'H2Ov';
ThetaField = 'THETA';

DEBUG_LEVEL = 2;

% By default the GUI folder may not be on my Matlab path.
GUIdir = '/Users/Josh/Documents/MATLAB/NO2 Profiles/GUIs/';
addpath(GUIdir);

load(ranges_file);
%Make sure the prefix ends in an underscore
if ~strcmpi(merge_prefix(end),'_'); merge_prefix = strcat(merge_prefix,'_'); end

field = repmat({[]},1,numel(Ranges));

% Prepare the Heights structure to accept BL height data
Heights = struct('Date',field,'BL_Heights', field, 'UTC_Times', field, 'Quality_Flags', field);

for a=1:numel(Ranges)
    curr_date = datestr(Ranges(a).Date, 29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    if ~strcmp(curr_date,'2008-06-24'); continue; end
    
    if DEBUG_LEVEL > 0; fprintf('Finding heights for %s\n',curr_date); end
    
    merge_name = sprintf('%s%s_%s_%s.mat',merge_prefix,year,month,day);
    merge_file = fullfile(merge_dir,merge_name);
    load(merge_file,'Merge');
    
    bl_heights = nan(1,size(Ranges(a).Ranges,1));
    utc_times = nan(1,size(Ranges(a).Ranges,1));
    qual_flags = nan(1,size(Ranges(a).Ranges,1));
    bl_fields = struct([]);
    
    for b=1:size(Ranges(a).Ranges,1)
        if DEBUG_LEVEL > 1; fprintf('\tGetting range %d of %d\n',b,size(Ranges(a).Ranges,1)); end
        range = Ranges(a).Ranges(b,:);
        structout = select_BL_heights(Merge,range,AltField,NO2field,H2Ofield,ThetaField);
        bl_heights(b) = structout.Overall.height;
        utc_times(b) = structout.Overall.medianUTC;
        qual_flags(b) = structout.Overall.qualityFlag;
        bl_fields(b).no2 = structout.field1;
        bl_fields(b).h2o = structout.field2;
        bl_fields(b).theta = structout.field3;
    end
    
    Heights(a).Date = curr_date;
    Heights(a).BL_Heights = bl_heights;
    Heights(a).UTC_Times = utc_times;
    Heights(a).Quality_Flags = qual_flags;
    Heights(a).Fields = bl_fields;
end
