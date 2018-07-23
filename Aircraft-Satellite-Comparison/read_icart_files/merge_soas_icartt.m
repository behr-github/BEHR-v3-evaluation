function [  ] = merge_soas_icartt(  )
%MERGE_SOAS_ICARTT Merge disparate SOAS icartt Merge structures
%   SOAS aircraft data didn't have a single Merge file for all data,
%   so this function takes the individual instrument data and merges
%   it into a unified structure. This needed some customization compared
%   to the existing combine_merges function, so it was simpler to
%   write a wholly separate function.

E = JLLErrors;

soas_root = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw/SOAS';
save_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignMergeMats/SOAS/P3/1sec';
save_prefix = 'SOAS_P3B_1sec_';
% Assume that each directory in soas_root is a different set of data and
% find all the .ict files in it. Start with AircraftPos_All because without
% that data we can't do much.

Mdirs = dir(soas_root);
Mdirs = remove_hidden(Mdirs); % remove the '.' and '..' entries
Mdirs = {Mdirs.name};

pos_dir = 'AircraftPos_All';
xx = strcmp(Mdirs, pos_dir);
if ~any(xx)
    warning('No aircraft position data detected')
else
    Mdirs(xx) = [];
    Mdirs = veccat({pos_dir}, Mdirs);
end

ict_files = cell(size(Mdirs));
ict_dnums = cell(size(Mdirs));
for a=1:numel(Mdirs)
    ict_files{a} = rfind(fullfile(soas_root, Mdirs{a}, '*.ict'), 'fullpath');
    ict_dnums{a} = filenames2dates(ict_files{a});
    if a > 1 && ~isequal(ict_dnums{a}, ict_dnums{1})
        E.notimplemented('ict files with different sets of dates');
    end
end

ict_dnums = ict_dnums{1}; % we required that all the different data sets have the same dates available and in the same order

for a=1:numel(ict_dnums)
    Merges = cell(size(ict_files));
    for b=1:numel(ict_files)
        Merges{b} = read_icartt_file(ict_files{b}(a).name, 'existing'); % choose existing variable names - should rename the time variable to UTC
    end
    Merge = combine_merges(Merges{:});
    Merge.Data = rename_field(Merge.Data, 'GpsLon', 'LONGITUDE');
    Merge.Data = rename_field(Merge.Data, 'GpsLat', 'LATITUDE');
    Merge.Data.LONGITUDE.Values = mod(Merge.Data.LONGITUDE.Values, 360); % put the longitude values on the "degrees east" standard
    save_name = sprintf('%s%s.mat', save_prefix, datestr(ict_dnums(a), 'yyyy_mm_dd'));
    save(fullfile(save_dir, save_name), 'Merge');
end


end

function F = remove_hidden(F)
xx = true(size(F));
for a=1:numel(F)
    if strcmp(F(a).name(1), '.')
        xx(a) = false;
    end
end
F = F(xx);
end

function dnums = filenames2dates(fnames)
dnums = nan(size(fnames));
for a=1:numel(fnames)
    dstr = regexp(fnames(a).name, '\d\d\d\d\d\d\d\d', 'match', 'once');
    dnums(a) = datenum(dstr, 'yyyymmdd');
end
end
