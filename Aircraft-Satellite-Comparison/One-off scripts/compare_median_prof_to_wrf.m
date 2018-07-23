function [  ] = compare_median_prof_to_wrf( campaign_name, wrf_month )
%compare_median_prof_to_wrf Make plots of median campaign profiles vs. WRF
%   

if nargin < 1
    campaign_name = 'discover-tx';
end
if nargin < 2
    wrf_month = 9;
end

% Conversion from WRF mixing ratios to aircraft mixing ratios - WRF is in
% ppm. This will be multiplied by the WRF profile.
wrf_to_aircraft = 1e6;

% Load the WRF profiles used in BEHR for this month
wrf_dir = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles/';
wrf_file = sprintf('m%02d_NO2_profile.mat',wrf_month);
load(fullfile(wrf_dir,wrf_file));

[Names,~,merge_dir] = merge_field_names(campaign_name);

% Loop through all the merge files. For each day, take the median profile
% for each site for the day between 12 and 15 LST. We'll plot each day's
% profiles separately.

F = dir(fullfile(merge_dir,'*.mat'));
no2_profs_by_day = struct;
prof_lons = [];
prof_lats = [];

% The structure template that each site should follow
%site_struct = make_empty_struct_from_cell({'no2','alts','no2bins','altbins','lon','lat'});

for a=1:numel(F)
    load(fullfile(merge_dir,F(a).name),'Merge');
    [no2,utc,~,lon,lat] = remove_merge_fills(Merge,Names.no2_lif);
    alts = remove_merge_fills(Merge,Names.gps_alt);
    pns = remove_merge_fills(Merge,Names.profile_numbers);
    
    tz = round(nanmean(lon/15));
    xx = utc > local2utc('12:00',tz) & utc < local2utc('15:00',tz);
    
    pns(~xx) = [];
    no2(~xx) = [];
    alts(~xx) = [];
    lon(~xx) = [];
    lat(~xx) = [];
    
    upns = unique(pns(pns~=0));
    % Handle both the MD and CA type profile numbers, which are 4 or 6
    % digits respectively. In both cases, the site number is in the
    % thousands place.
    sites = (mod(pns,100000) - mod(pns,1000))/1000;
    usites = unique((mod(upns,100000) - mod(upns,1000))/1000);
    
    for b=1:numel(usites)
        % no2_profs_by_day will be a 1 x numel(F) structure, with fields
        % for each site number that contain binned NO2 and altitude, plus
        % the lat and lon
        site_field = sprintf('Site%02d',usites(b));
        ss = sites == usites(b);
        [no2bins, no2altbins] = bin_rolling_vertical_profile(alts(ss),no2(ss),0.5,0.1);
        % We'll use the raw NO2 values at the end to calculate a median
        % profile for the whole campaign and the individual profiles to
        % show the daily variation.
        no2_profs_by_day(a).(site_field).no2 = no2(ss);
        no2_profs_by_day(a).(site_field).alts = alts(ss);
        no2_profs_by_day(a).(site_field).no2bins = no2bins;
        no2_profs_by_day(a).(site_field).altbins = no2altbins;
        no2_profs_by_day(a).(site_field).lon = nanmean(lon(ss));
        no2_profs_by_day(a).(site_field).lat = nanmean(lat(ss));
    end
end

% Now, for each site, find the closest BEHR WRF profile and save it to its
% own structure
fns = fieldnames(no2_profs_by_day);
wrf_profs_by_site = struct;
% All the WRF profiles are on the same vertical coordinates. Convert from
% pressure to altitude in km
wrf_alt = squeeze(-log(PROFILE.Pressure/1013)*7.4);
for a=1:numel(fns)
    all_prof_lats = [];
    all_prof_lons = [];
    for b=1:numel(no2_profs_by_day)
        if ~isempty(no2_profs_by_day(b).(fns{a}))
            all_prof_lats = [all_prof_lats, no2_profs_by_day(b).(fns{a}).lat];
            all_prof_lons = [all_prof_lons, no2_profs_by_day(b).(fns{a}).lon];
        end
    end
    wrf_lat = nanmean(all_prof_lats);
    wrf_lon = nanmean(all_prof_lons);
    
    % Find the profile with the smallest distance to the given lat and lon
    dlat = PROFILE.Latitude - wrf_lat;
    dlon = PROFILE.Longitude - wrf_lon;
    dist = sqrt(dlat.^2 + dlon.^2);
    [~,min_dist_ind] = min(dist(:));
    [i,j] = ind2sub(size(dist),min_dist_ind);
    wrf_profs_by_site.(fns{a}).no2 = squeeze(PROFILE.NO2_profile(:,i,j))*wrf_to_aircraft;
    wrf_profs_by_site.(fns{a}).alt = wrf_alt;
    wrf_profs_by_site.(fns{a}).lat = PROFILE.Latitude(i,j);
    wrf_profs_by_site.(fns{a}).lon = PROFILE.Longitude(i,j);
end

% Finally make the plot for each site. Each day's profile will be draw in
% light grey with the median profile for the whole campaign in black. The
% WRF profile gets superimposed in red.
for a=1:numel(fns)
    % First collect the no2 and alts for this site for the whole campaign
    % and bin it (median)
    ThisSite = [no2_profs_by_day(:).(fns{a})];
    this_site_no2 = [ThisSite(:).no2];
    this_site_alts = [ThisSite(:).alts];
    [med_no2_prof, med_no2_alts] = bin_rolling_vertical_profile(this_site_alts,this_site_no2,0.5,0.1);
    
    figure;
    for b=1:numel(ThisSite)
        line(ThisSite(b).no2bins, ThisSite(b).altbins, 'color', [0.7 0.7 0.7]);
    end
    l(1) = line(med_no2_prof, med_no2_alts, 'color','k','linewidth',2);
    l(2) = line(wrf_profs_by_site.(fns{a}).no2, wrf_profs_by_site.(fns{a}).alt, 'color', 'r', 'linestyle', '--', 'linewidth',2);
    legend(l',{'Median campaign','WRF'});
    set(gca,'fontsize',16);
    xlabel('[NO2]');
    ylabel('Altitude (km)');
    titlestr = sprintf('%s: %s',upper(campaign_name),fns{a});
    title(titlestr);
end
end