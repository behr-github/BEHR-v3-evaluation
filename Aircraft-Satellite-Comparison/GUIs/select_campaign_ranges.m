% select_campaign_ranges - Wrapper script that will iterate through all the
% days specified, load any day with aircraft data into the select changing
% altitude GUI and save the returned UTC ranges in a structure for that
% campaign.

campaign_name = 'soas';
[~,dates,mat_dir] = merge_field_names(campaign_name);

start_date = dates{1};
end_date = dates{2};

R = 0;
datenums = datenum(start_date):datenum(end_date);
for d=1:numel(datenums)
    curr_date = datestr(datenums(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    fname = sprintf('*%s_%s_%s.mat',year,month,day);
    fpath = fullfile(mat_dir,fname);
    
    files = dir(fpath);
    if numel(files) < 1
        fprintf('%s not found\n',fpath);
        continue
    elseif numel(files) > 1
        fprintf('Many files found for %s\n.',fname);
        continue
    else
        filename = fullfile(mat_dir, files(1).name);
        r = select_changing_altitude(filename);
        R = R+1;
        Ranges(R).Date = datestr(curr_date,'mm/dd/yyyy');
        Ranges(R).Ranges = r;
    end
end