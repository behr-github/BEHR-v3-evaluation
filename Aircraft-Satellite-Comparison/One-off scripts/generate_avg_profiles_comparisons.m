% generate_avg_profile_comparisons
%
%   Script that runs compare_profile averaging all profiles within a
%   certain time range

start_date = '09/01/2013';
end_date = '09/29/2013';

no2field = 'NO2_MixingRatio';
profnumfield = 'ProfileNumber';

start_time = '12:00';
end_time = '15:00';
timezone = 'est';

DEBUG_LEVEL = 2;

profile_dir = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles/';
merge_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

dates = datenum(start_date):datenum(end_date);
oldprofilefile = '';

for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    % Load the profile file (if necessary) and the merge file for this day
    
    profile_file = sprintf('m%s_NO2_profile.mat',month);
    if ~strcmpi(profile_file, oldprofilefile) %If the profile file from the last loop is still in the same month, we don't need to reload it.
        load(fullfile(profile_dir,profile_file));
        oldprofilefile = profile_file;
    end
    
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('gen_prof:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    % Calculate what profile numbers have start times within the time
    % window chosen
    utcmin = local2utc(start_time,timezone);
    utcmax = local2utc(end_time,timezone);
    
    utc = remove_merge_fills(Merge,'UTC');
    profnums = remove_merge_fills(Merge,profnumfield);
    u_profnums = unique(profnums(profnums ~= 0 & ~isnan(profnums)));
    profnums_in_time = false(size(u_profnums));
    for a=1:numel(u_profnums)
        xx = find(profnums == u_profnums(a),1,'first');
        if utc(xx) >= utcmin && utc(xx) <= utcmax;
            profnums_in_time(a) = true;
        end
    end
    final_profnums = u_profnums(profnums_in_time);
    
    compare_profile(Merge,no2field,PROFILE,final_profnums,profnumfield);
    
end