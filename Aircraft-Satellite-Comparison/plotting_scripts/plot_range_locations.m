function plot_range_locations(campaign, range_file)
%PLOT_RANGE_LOCATIONS Plot the UTC range locations for non-DISCOVER campaigns
%   PLOT_RANGE_LOCATIONS( CAMPAIGN, RANGE_FILE ) Given the name of a
%   campaign (CAMPAIGN) recognized by merge_field_names() and the path to a
%   UTC range file, will plot the location of all ranges.

R = load(range_file);
Ranges = R.Ranges;

[Names, ~, campaign_dir] = merge_field_names(campaign);
merge_files = dirff(fullfile(campaign_dir, '*.mat'));
range_lon = [];
range_lat = [];
for i_range = 1:numel(Ranges)
    Merge = load_merge_for_date(Ranges(i_range).Date, merge_files);
    lon = remove_merge_fills(Merge, Names.longitude);
    lat = remove_merge_fills(Merge, Names.latitude);
    utc = remove_merge_fills(Merge, Names.utc);
    
    for j_range = 1:size(Ranges(i_range).Ranges,1)
        this_range = Ranges(i_range).Ranges(j_range, :);
        xx = utc >= this_range(1) & utc < this_range(2);
        range_lon = [range_lon; nanmean(lon(xx))];
        range_lat = [range_lat; nanmean(lat(xx))];
    end
    
end

figure; 
plot(range_lon, range_lat, 'ro');
state_outlines('k');

end

function Merge = load_merge_for_date(date_in, merge_files)
E = JLLErrors;
filenames = {merge_files.name};
date_regex = datestr(date_in, 'yyyy_mm_dd');
xx = regcmp(filenames, date_regex);
if sum(xx) ~= 1 
    E.filenotfound('Merge file for %s', date_in);
end

fprintf('Loading %s\n', filenames{xx});
M = load(filenames{xx}, 'Merge');
Merge = M.Merge;

end