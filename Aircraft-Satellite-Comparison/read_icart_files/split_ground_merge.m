function [  ] = split_ground_merge(  )
%split_ground_merge Splits ground station data into individual days
%   The ground files for DISCOVER-AQ California combine all days into a
%   single file; since some of my code expects files for each day, I need
%   to separate these into individual days' files.

E = JLLErrors;
DEBUG_LEVEL = 1;

[loadfile, loaddir] = uigetfile('*.mat','Choose to file to split');
[savedir] = uigetdir('/Volumes', 'Select the directory to save the split files to');

if (isscalar(loadfile) && loadfile == 0) || (isscalar(savedir) && savedir == 0)
    E.userCancel;
end

M=load(fullfile(loaddir,loadfile),'Merge');

mergefields = fieldnames(M.Merge.Data);
[dayfield, ok] = listdlg('liststring',mergefields,'PromptString','Choose the field that indicates the day or date');

if ok == 0
    E.userCancel;
end

userans = questdlg('Is this a day-of-year or date field?','Day field type','Day-of-year','Date','Cancel','Cancel');

switch userans
    case 'Cancel'
        E.userCancel;
    case 'Day-of-year'
        useryear = inputdlg('You have selected day-of-year. What year is the data from (will be used in file names)?','Enter year');
        if isempty(useryear)
            E.userCancel;
        end
        
        doy = floor(M.Merge.Data.(mergefields{dayfield}).Values);
        udoy = unique(doy);
        udates = cell(size(udoy));
        for a=1:numel(udoy)
            udates{a} = datestr(modis_day_to_date(udoy(a),useryear{1}),29);
        end
    case 'Date'
        user_date_format = inputdlg('You have selected date. Describe the date format, using standard Matlab date format strings.','Date format');
        dates = M.Merge.Data.(mergefields{dayfield}).Values;
        if ~isa(dates,'double')
            % It should be that any date will be read in as a number, e.g.
            % 20130901. If not, this function will not be able to handle
            % that.
            E.badvartype(dates,'double');
        end
        % Make dates a column vector, then convert to a cell array - this
        % prevents num2str from putting spaces in there
        if ~iscolumn(dates); dates = dates'; end
        date_chararray = num2str(dates);
        tmp_datenums = datenum(date_chararray, user_date_format{1});
        dates = cellstr(datestr(tmp_datenums,29));
        
        % Make a doy vector for use in the main function
        doy = nan(size(dates));
        for a=1:numel(dates)
            if DEBUG_LEVEL > 1; 
                fprintf('\t Converting to doy: %d of %d\n',a,numel(dates));
            end
            doy(a) = modis_date_to_day(dates{a});
        end
        udoy = unique(doy);
        udates = cell(size(udoy));
        for a=1:numel(udoy)
            udates{a} = datestr(modis_day_to_date(udoy(a),year(dates{1})),29);
        end
        
end


% The new files will use a save name prefix that is the same as the old
% file, with any date removed (look for yyyy-mm-dd or yyyy_mm_dd or
% yyyymmdd)

save_prefix = regexprep(loadfile,'\d\d\d\d[-_]?\d\d[-_]?\d\d','');
save_prefix = regexprep(save_prefix,'\.\w+','');
if ~strcmp(save_prefix(end),'_')
    save_prefix = strcat(save_prefix,'_');
end

% Find the data that belongs to each day

for a=1:numel(udoy)
    xx = doy == udoy(a);
    Merge.metadata = M.Merge.metadata;
    Merge.Data = make_empty_struct_from_cell(mergefields);
    for b=1:numel(mergefields)
        Merge.Data.(mergefields{b}) = M.Merge.Data.(mergefields{b});
        if isfield(Merge.Data.(mergefields{b}),'Values')
            % Handle the fact that sometime the UTC field doesn't have
            % values because it's actually called something else.
            Merge.Data.(mergefields{b}).Values(~xx) = [];
        end
    end
    
    yr = udates{a}(1:4);
    mn = udates{a}(6:7);
    dy = udates{a}(9:10);
    
    savename = sprintf('%s%s%s%s.mat',save_prefix,yr,mn,dy);
    save(fullfile(savedir,savename),'Merge');
end


end

