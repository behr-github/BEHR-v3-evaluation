function ground_site_super_merge()

% ground_site_super_merge
%
%   Merges data from ground sites in the DISCOVER campaign to a single
%   super merge that is organized as Merge --> Date --> Site --> Data
%   fields.  All information needed from the user will be presented as
%   dialogue boxes, however, note that to work properly all files to be
%   merged together must be contained in a single directory with no sub
%   folders and the files should have Ground#_ somewhere in the name, where
%   # is the site number.
%
%   Josh Laughner <joshlaugh5@gmail.com> 11 Aug 2014

% Have the user select the matlab files that have the ground data to be
% merged.
filepath = uigetdir('/','Select the target folder');


% Find any hidden files (which on Macs start with a .) and remove them
files = dir(filepath);
filenames = {files(:).name};
xx = strncmp('.',filenames,1);
files = files(~xx);

fields = {};
SuperMerge = struct();

% Loop through each file
for a=1:numel(files)
    name = files(a).name;
    % Find the site number
    x1 = regexpi(name,'Ground\d+');
    x2 = regexp(name(x1:end),'_','ONCE');
    sitenum = str2double(name((x1+6):(x1+x2-2)));
    
    %Load the merge file
    load(fullfile(filepath, name),'Merge');
    
    % If the user has already defined which NO2 field to copy, we will read
    % that from the "fields" array, otherwise, have the user indicate which
    % NO2 field to use.  Automatically assign any field with the string
    % 'UTC' or 'Time' in it to be copied
    if numel(fields) >= sitenum && ~isempty(fields{sitenum})
    else
        mergefields = fieldnames(Merge.Data);
        selection = listdlg('ListString',mergefields, 'PromptString', 'Select the NO2 field to use');
        if isempty(selection)
            error('super_merge:user','User cancelled');
        else
            thisfields = {};
            for b=1:numel(mergefields)
                if ~isempty(strfind(mergefields{b},'UTC')) || ~isempty(strfind('Time',mergefields{b}))
                    thisfields{end+1} = mergefields{b};
                end
            end
            thisfields{end+1} = mergefields{selection};
            fields{sitenum} = thisfields;
        end
    end
    
    % Now copy all the fields listed, renaming the last one just to 'NO2'
    thisfields = fields{sitenum};
    datefield = sprintf('Date_%s',datestr(Merge.metadata.date,'yyyymmdd'));
    sitefield = sprintf('Site_%02d',sitenum);
    for b = 1:(numel(thisfields)-1)
        SuperMerge.(datefield).(sitefield).(thisfields{b}) = Merge.Data.(thisfields{b});
    end
    SuperMerge.(datefield).(sitefield).NO2 = Merge.Data.(thisfields{end});
    
    % Now include a field identifying the original file and the original
    % NO2 field name
    SuperMerge.(datefield).(sitefield).OriginalFile = fullfile(filepath,name);
    SuperMerge.(datefield).(sitefield).NO2.OriginalField = thisfields{end};
end

ind = regexp(name,'_','ONCE');
cityname = name(1:(ind-1));
savename = sprintf('%s_GroundSites');
uisave('SuperMerge',savename);

end