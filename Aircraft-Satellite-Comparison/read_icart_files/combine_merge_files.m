function combine_merge_files
%
%   Allows manual combination of .mat files resulting from
%   read_merge_data.m. Combines two file sets at a time, to combine more
%   files, simply run this script again, combining the previous merge
%   output with the next file.  The first file chosen should probably be
%   the NAV data file.

lonfields = {'FMS_LON','Longitude'};
latfields = {'FMS_LAT','Latitude'};

save_unmatched = false;
DEBUG_LEVEL = 2;
mat_dir_1 = uigetdir('/Volumes','Select the directory with the first set of merge files in it');
fprintf('Will read .mat file set 1 from %s\n',mat_dir_1);
mat_dir_2 = uigetdir('/Volumes','Select the directory with the second set of merge files in it');
fprintf('Will read .mat file set 2 from %s\n',mat_dir_2);
save_path = uigetdir('/Volumes','Select the directory to save the resulting .mat files to');
fprintf('Will save mat files to %s\n',save_path);
save_bool = 1;
if all(mat_dir_1 == 0) || all(mat_dir_2 == 0) || all(save_path == 0); error('read_merge_data:user_cancel','User canceled run.'); end

savetitle = sprintf('Enter the save file name. This will preceed the data date.');
user_savename = inputdlg(savetitle); user_savename = user_savename{1};
if isempty(user_savename) 
    error('read_merge_data:user_cancel','User cancelled run');
elseif ~strcmp(user_savename(end),'_');
    user_savename = [user_savename,'_'];
end

% mat_dir_1 = '/Volumes/share/GROUP/INTEX-B/Unmerged files/Nav files/Matlab/';
% mat_dir_2 = '/Volumes/share/GROUP/INTEX-B/Unmerged files/Cohen-NO2 files/Matlab2/';
% save_path = '/Volumes/share/GROUP/INTEX-B/Matlab files/';
% user_savename = 'INTEXB_';

% Flag to know whether to fix the order of the fields
lonflag = 0;

merge_files_1 = dir(fullfile(mat_dir_1,'*.mat'));
merge_files_2 = dir(fullfile(mat_dir_2,'*.mat'));

% Make lists of the dates (as datenums) of each set of files, we will need
% this to match up the files by date.

datenums_1 = zeros(size(merge_files_1));
for a=1:numel(merge_files_1);
    ind = regexp(merge_files_1(a).name,'\d\d\d\d_\d\d_\d\d');
    ind2 = regexp(merge_files_1(a).name,'\d.mat');
    datenums_1(a) = datenum(merge_files_1(a).name(ind:ind2));
end


datenums_2 = zeros(size(merge_files_2));
for a=1:numel(merge_files_2);
    ind = regexp(merge_files_2(a).name,'\d\d\d\d_\d\d_\d\d');
    ind2 = regexp(merge_files_2(a).name,'\d.mat');
    datenums_2(a) = datenum(merge_files_2(a).name(ind:ind2));
end

% Go through each file fromt the first directory and match it with the
% corresponding file from the second.  If there is no corresponding file,
% either resave the first .mat file in the new save directory or just skip
% - this behavior is set by the 'save_unmatched' boolean.

first_time = true;
for a=1:numel(datenums_1);
    b = find(datenums_2 == datenums_1(a));
    curr_date = datestr(datenums_1(a),29);
    if DEBUG_LEVEL > 0; fprintf('Now combining files for %s.\n',curr_date); end
    year = curr_date(1:4); month = curr_date(6:7); day = curr_date(9:10);
    if isempty(b)
        if save_unmatched;
            Merge1 = load(fullfile(mat_dir_1,merge_files_1(a).name),'Merge');
            save(fullfile(save_path,sprintf('%s_%s_%s_%s.mat',user_savename,year,month,day)),'Merge1');
        else
            continue
        end
    end
    
    load(fullfile(mat_dir_1,merge_files_1(a).name),'Merge'); Merge1 = Merge; clear Merge
    load(fullfile(mat_dir_2,merge_files_2(b).name),'Merge'); Merge2 = Merge; clear Merge
    
    % Check that the UTC field actually has values in it; if not, ask the
    % user to pick which field to consider the UTC field.  Only do this
    % once.
    if first_time
        if ~isfield(Merge1.Data.UTC,'Values') || isempty(Merge1.Data.UTC.Values)
            fields1 = fieldnames(Merge1.Data);
            fields_wildcards = [repmat(' %s,',1,numel(fields1)-1),' %s.'];
            dlgstr = ['UTC field empty, enter the time field name.  Choices are:',fields_wildcards];
            fieldname1 = inputdlg(sprintf(dlgstr,fields1{:})); fieldname1 = fieldname1{1};
        else
            fieldname1 = 'UTC';
        end
        
        if ~isfield(Merge2.Data.UTC,'Values') || isempty(Merge2.Data.UTC.Values)
            fields2 = fieldnames(Merge2.Data);
            fields_wildcards = [repmat(' %s,',1,numel(fields2)-1),' %s.'];
            dlgstr = ['UTC field empty, enter the time field name.  Choices are:',fields_wildcards];
            fieldname2 = inputdlg(sprintf(dlgstr,fields2{:})); fieldname2 = fieldname2{1};
        else
            fieldname2 = 'UTC';
        end
        first_time = false;
    end
    UTC1 = eval(sprintf('Merge1.Data.%s.Values',fieldname1));
    UTC2 = eval(sprintf('Merge2.Data.%s.Values',fieldname2));
    UTCstart = min([min(UTC1(:)), min(UTC2(:))]);
    UTCend = max([max(UTC1(:)), max(UTC2(:))]);
    
    UTCall = UTCstart:UTCend;
    
    newMerge.metadata = Merge1.metadata;
    newMerge.Data.UTC.Unit = 'seconds after midnight UTC';
    newMerge.Data.UTC.Fill = NaN;
    newMerge.Data.UTC.Values = UTCall;
    
    fields1 = fieldnames(Merge1.Data);
    for n=1:numel(fields1)
        % Don't overwrite the UTC field in newMerge, and don't copy the old
        % time field
        if any(strcmp(fields1{n},{'UTC',fieldname1}));
            continue
        else
            eval(sprintf('newMerge.Data.%s = Merge1.Data.%s;',fields1{n},fields1{n}))
            eval(sprintf('newMerge.Data.%s.Values = zeros(size(UTCall));',fields1{n}));
            % Match up UTC times
            if DEBUG_LEVEL > 1; fprintf('   Matching %s from Merge1 to total UTC times.\n',fields1{n}); end
            for i=1:numel(UTCall)
                xx = find(UTC1 == UTCall(i));
                if ~isempty(xx)
                    eval(sprintf('newMerge.Data.%s.Values(i) = Merge1.Data.%s.Values(xx);',fields1{n},fields1{n}));
                else
                    eval(sprintf('newMerge.Data.%s.Values(i) = Merge1.Data.%s.Fill;',fields1{n},fields1{n}));
                end
            end
            unit = eval(sprintf('newMerge.Data.%s.Unit',fields1{n}));
            fills = eval(sprintf('newMerge.Data.%s.Values == newMerge.Data.%s.Fill',fields1{n},fields1{n}));
            % Convert from Celsius to Kelvin
            if strcmpi(unit,'degC') || strcmpi(unit,'deg C');
                if DEBUG_LEVEL > 1; fprintf('      Converting %s to units of Kelvin\n',fields1{n}); end
                eval(sprintf('newMerge.Data.%s.Values(~fills) = newMerge.Data.%s.Values(~fills) + 273;',fields1{n},fields1{n}))
                eval(sprintf('newMerge.Data.%s.Unit = ''K'';',fields1{n}));
            end
            % Convert from meters to kilometers
            if strcmpi(unit,'m')
                if DEBUG_LEVEL > 1; fprintf('      Converting %s to units of kilometers\n',fields1{n}); end
                eval(sprintf('newMerge.Data.%s.Values(~fills) = newMerge.Data.%s.Values(~fills)/1000;',fields1{n},fields1{n}))
                eval(sprintf('newMerge.Data.%s.Unit = ''km'';',fields1{n}));
            end
            % If the field is a latitude or longitude field rename it to
            % LATITUDE or LONGITUDE for consitency.  The any() statement is
            % in place to allow easy expansion of the fields to rename.
            if any(strcmp(lonfields,fields1{n}))
                lonflag=1;
                newMerge.Data.LONGITUDE = newMerge.Data.(fields1{n});
                newMerge.Data = rmfield(newMerge.Data,fields1{n});
            end
            if any(strcmp(latfields,fields1{n}))
                newMerge.Data.LATITUDE = newMerge.Data.(fields1{n});
                newMerge.Data = rmfield(newMerge.Data,fields1{n});
            end
        end
    end
    
    fields2 = fieldnames(Merge2.Data);
    for n=1:numel(fields2)
        myname = fields2{n};
        if any(strcmp(myname,{'UTC',fieldname2}))
            continue
        else
            % Check that the field name doesn't exist, if it does, append a
            % numerical index (and make sure that doesn't duplicate)
            while any(strcmp(myname,fields1));
                ind = str2double(myname(end));
                if isnan(ind)
                    myname = [myname,'_1'];
                else
                    myname(end)=num2str(ind+1);
                end
            end
            
            eval(sprintf('newMerge.Data.%s = Merge2.Data.%s;',myname,fields2{n}))
            eval(sprintf('newMerge.Data.%s.Values = zeros(size(UTCall));',myname));
            % Match up UTC times
            if DEBUG_LEVEL > 1; fprintf('   Matching %s from Merge2 to total UTC times as %s.\n',fields2{n},myname); end
            for i=1:numel(UTCall)
                xx = find(UTC2 == UTCall(i));
                if ~isempty(xx)
                    eval(sprintf('newMerge.Data.%s.Values(i) = Merge2.Data.%s.Values(xx);',myname,fields2{n}));
                else
                    eval(sprintf('newMerge.Data.%s.Values(i) = Merge2.Data.%s.Fill;',myname,fields2{n}));
                end
            end
            unit = eval(sprintf('newMerge.Data.%s.Unit',myname));
            fills = eval(sprintf('newMerge.Data.%s.Values == newMerge.Data.%s.Fill',myname,myname));
            % Convert from Celsius to Kelvin
            if strcmpi(unit,'degC') || strcmpi(unit,'deg C');
                if DEBUG_LEVEL > 1; fprintf('      Converting %s to units of Kelvin\n',myname); end
                eval(sprintf('newMerge.Data.%s.Values(~fills) = newMerge.Data.%s.Values(~fills) + 273;',myname,myname))
                eval(sprintf('newMerge.Data.%s.Unit = ''K'';',myname));
            end
            if strcmpi(unit,'m')
                if DEBUG_LEVEL > 1; fprintf('      Converting %s to units of kilometers\n',myname); end
                eval(sprintf('newMerge.Data.%s.Values(~fills) = newMerge.Data.%s.Values(~fills)/1000;',myname,myname))
                eval(sprintf('newMerge.Data.%s.Unit = ''km'';',myname));
            end
            % If the field is a latitude or longitude field rename it to
            % LATITUDE or LONGITUDE for consitency.  The any() statement is
            % in place to allow easy expansion of the fields to rename.
            if any(strcmp(lonfields,fields2{n}))
                lonflag=1;
                newMerge.Data.LONGITUDE = newMerge.Data.(fields2{n});
                newMerge.Data = rmfield(newMerge.Data,fields2{n});
            end
            if any(strcmp(latfields,fields2{n}))
                newMerge.Data.LATITUDE = newMerge.Data.(fields2{n});
                newMerge.Data = rmfield(newMerge.Data,fields2{n});
            end
        end
    end
    % If the lat/lon fields were renamed, reorder the fields so that they
    % are 2nd and 3rd (assuming UTC is the first field
    if lonflag
        fieldlist = fieldnames(newMerge.Data);
        xx = strcmp(fieldlist,'LATITUDE') | strcmp(fieldlist, 'LONGITUDE');
        fieldlist(xx) = [];
        fieldlist = [fieldlist(1); {'LATITUDE';'LONGITUDE'}; fieldlist(2:end)];
        newMerge.Data = orderfields(newMerge.Data,fieldlist);
    end
    
    % Convert longitude into the [0 360] representation
    if DEBUG_LEVEL > 1; fprintf('      Converting longitude to [0 360] representation\n'); end
    newMerge.Data.LONGITUDE.Values = mod(newMerge.Data.LONGITUDE.Values,360);
    
    % Save the new merge file as the variable 'Merge'
    Merge = newMerge;
    save(fullfile(save_path,sprintf('%s_%s_%s_%s.mat',user_savename,year,month,day)),'Merge');
    clear('newMerge'); clear('Merge');
end

end