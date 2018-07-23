% TOD & Site number
% Plot how vertical profiles changes for each of the DISCOVER sites by day,
% separating by time of day as well

matdir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

% Reading in the data and sorting by site
location = 'Baltimore'; numberofsites = 6;
files = dir(fullfile(matdir,sprintf('%s*.mat',location)));
clear profiles
for a=1:numel(files);
    load(fullfile(matdir,files(a).name),'Merge');
    [no2, utc_no2, alt_no2, lon_no2, lat_no2, xx_no2] = remove_merge_fills(Merge,'NO2_LIF');
    lt_no2 = Merge.Data.LOCAL_SUN_TIME.Values; lt_no2(xx_no2) = NaN;
    prof_no2 = Merge.Data.ProfileSequenceNum.Values; prof_no2(xx_no2) = NaN;
    u_prof_no2 = unique(prof_no2(prof_no2~=0 & ~isnan(prof_no2)));
    
    % Initialize the data structure
    profiles(a).date = Merge.metadata.date;
    profiles(a).starthour(24).vals = [];
    for b=1:numel(u_prof_no2)
        % Sort profile numbers according to what hour (solar time) the
        % spiral started in
        pp = prof_no2 == u_prof_no2(b);
        no2_p = no2(pp); alt_p = alt_no2(pp); lt_p = lt_no2(pp); utc_p = utc_no2(pp);
        starthour = floor(min(lt_p(:)));
        [val,alt,err,weight] = bin_rolling_vertical_profile(alt_p,no2_p,0.5,0.1);
        
        % Save data in structure
        cellindex = numel(profiles(a).starthour(starthour).vals)+1;
        profiles(a).starthour(starthour).vals{cellindex} = val;
        profiles(a).starthour(starthour).alts{cellindex} = alt;
        profiles(a).starthour(starthour).errs{cellindex} = err;
        profiles(a).starthour(starthour).weights{cellindex} = weight;
        profiles(a).starthour(starthour).localtimes{cellindex} = [min(lt_p(:)), max(lt_p(:))];
        profiles(a).starthour(starthour).utctimes{cellindex} = [min(utc_p(:)), max(utc_p(:))];
        profiles(a).starthour(starthour).siteflag{cellindex} = u_prof_no2(b);
    end
end

% Refine and plot the data; average any profiles from the same site/same
% hour, create plots that depict the variation by day in vertical profiles at each
% site at each hour


% First figure out which top level indices of "profile" actually have data:
% go through each hour of each day and test if there is data.  Once one day
% is found to have data for that hour, it will be plotted.
ptest = false(1,24);
for d=1:numel(profiles)
    for p=1:numel(profiles(d).starthour);
        ptest(p) = ~isempty(profiles(d).starthour(p).vals) || ptest(p); 
    end
end
xx = find(ptest);
n_xx = numel(xx);

colors = {'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],};
markers = {'s','s','s','s','s','^','^','^','^','^','o','o','o','o','o'};

for d=1:numel(profiles)
    for i=1:n_xx
        
    end
end