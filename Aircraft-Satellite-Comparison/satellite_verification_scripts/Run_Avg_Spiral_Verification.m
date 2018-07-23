function varargout = Run_Avg_Spiral_Verification(  )
%RUN_AVG_SPIRAL_VERIFICATION Compare campaign average profiles to satellite averages
%   Detailed explanation goes here


% First we'll use Run_Spiral_Verification to gather all the profile
% information for each of the four discover campaigns.

campaigns = {'discover-md','discover-ca','discover-tx','discover-co'};
blank_array = [];
prof_lon = blank_array;
prof_lat = blank_array;
prof_vcd = blank_array;
prof_all_lon = blank_array;
prof_all_lat = blank_array;
prof_nums = blank_array;
prof_dates = blank_array;
prof_campaign_ind = blank_array;

for a=1:numel(campaigns)
    [prof_lon_a, prof_lat_a, ~, ~, prof_vcd_a, prof_db_a, prof_dates_a] = Run_Spiral_Verification('all', campaigns{a});
    prof_lon = cat(1, prof_lon, prof_lon_a);
    prof_lat = cat(1, prof_lat, prof_lat_a);
    prof_vcd = cat(1, prof_vcd, prof_vcd_a);
    prof_dates = cat(1, prof_dates, prof_dates_a);
    prof_all_lon = cat(1, prof_all_lon, prof_db_a.lon_3km(:));
    prof_all_lat = cat(1, prof_all_lat, prof_db_a.lat_3km(:));
    prof_nums = cat(1, prof_nums, prof_db_a.profnums{:});
    prof_campaign_ind = cat(1, prof_campaign_ind, repmat(a, size(prof_lon_a)));
end

% Now we need to figure out which profiles are at the same site. The first
% step is to construct a campaign-site ID number.

prof_site_id = sitenum_from_profnum(prof_nums) + 100*prof_campaign_ind;
u_prof_site_id = unique(prof_site_id);
u_prof_site_id(isnan(u_prof_site_id)) = [];

% Now we can find the unique sites across all four campaigns. As we go
% through each, we will average the aircraft VCDs together. Then we need to
% put together the information necessary information to figure out which
% grid cells in the averaging methods match up with the aircraft data.
%
% The existing grid dumping a.k.a. constant value method (CVM) used by BEHR
% outputs the lat/lon as the center of the grid points, which are 0.05x0.05
% degrees. The PSM method is outputting mass-conserved values on a
% parabolic spline, so there isn't the same "grid box that values within
% are averaged," but we should be able to treat it similarly. I'm also
% gridding it to 0.05x0.05 degrees.
%
% Similar to Haines, and to the existing spiral verification method, we
% will require that there be 20 points below 3 km in a given grid cell on a
% given day to include it in the average.

blank_cell = cell(size(u_prof_site_id));
all_air_vcds = blank_cell;
all_cvm_vcds = blank_cell;
all_cvm_wts = blank_cell;
all_psm_vcds = blank_cell;
all_psm_wts = blank_cell;
all_air_lons = blank_cell;
all_air_lats = blank_cell;

behr_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/PSM-Comparison/BEHR_Files';

% Load one CVM file to get the lat/lon grid;
O = load(fullfile(behr_dir, behr_filename('2011-07-01')),'OMI');
[lon_edge, lat_edge] = compute_grid_edges(O.OMI(1).Longitude, O.OMI(1).Latitude);

for a=1:numel(u_prof_site_id)
    fprintf('Site %d\n', u_prof_site_id(a));
    site_inds = find(u_prof_site_id(a) == prof_site_id);
    site_dates = prof_dates(site_inds);
    
    all_air_vcds{a} = prof_vcd(site_inds);
    
    all_air_lons{a} = veccat(prof_all_lon{site_inds}); 
    all_air_lats{a} = veccat(prof_all_lat{site_inds});
    
    [all_cvm_vcds{a}, all_cvm_wts{a}, all_psm_vcds{a}, all_psm_wts{a}] = bin_sat_gridcells(prof_all_lon(site_inds), prof_all_lat(site_inds), site_dates, lon_edge, lat_edge);
end

% All thats left is the averaging
blank_array2 = nan(size(u_prof_site_id));
avg_air_vcds = blank_array2;
avg_cvm_vcds = blank_array2;
avg_psm_vcds = blank_array2;
avg_lons = blank_array2;
avg_lats = blank_array2;

fprintf('Averaging...\n');
for a=1:numel(u_prof_site_id)
    
    avg_lons(a) = nanmean(all_air_lons{a});
    avg_lats(a) = nanmean(all_air_lats{a});
    avg_air_vcds(a) = nanmean(all_air_vcds{a});
    avg_cvm_vcds(a) = nansum2(all_cvm_vcds{a} .* all_cvm_wts{a}) ./ nansum2(all_cvm_wts{a});
    avg_psm_vcds(a) = nansum2(all_psm_vcds{a} .* all_psm_wts{a}) ./ nansum2(all_psm_wts{a});
end

if nargout > 0
    varargout{1} = avg_air_vcds;
    varargout{2} = avg_cvm_vcds;
    varargout{3} = avg_psm_vcds;
    varargout{4} = avg_lons;
    varargout{5} = avg_lats;
    varargout{6} = u_prof_site_id;
else
    putvar(avg_air_vcds, avg_cvm_vcds, avg_psm_vcds, avg_lons, avg_lats, u_prof_site_id);
end
end

function [cvm_vcds, cvm_wts, psm_vcds, psm_wts] = bin_sat_gridcells(air_lons, air_lats, site_dates, lon_edge, lat_edge)
% We want any grid cell that has >= 20 data points inside it below 3 km
% (which are the only lat/lon points given) for that day. Assuming that CVM
% and PSM grid cells are on the same grid, the first step is to find the
% indicies of the grid cells for each date. Then we can grab the proper
% satellite VCDs, filter them, and add them to the list.

psm_behr_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/PSM-Comparison/BEHR_Files';

cvm_vcds = [];
cvm_wts = [];
psm_vcds = [];
psm_wts = [];

u_dates = unique(site_dates);
for a=1:numel(u_dates)
    fprintf('    Working on date %s\n', u_dates{a})
    day_inds = strcmp(site_dates, u_dates{a});
    day_lons = veccat(air_lons{day_inds});
    day_lats = veccat(air_lats{day_inds});
    
    include_bool = false(numel(lat_edge)-1, numel(lon_edge)-1);
    % Most of the country will not have the campaign occurring over it;
    % this will cut the loop down to the more relevant grid cells.
    y1 = find(lat_edge < min(day_lats(:)), 1, 'last') - 1;
    y2 = find(lat_edge > max(day_lats(:)), 1, 'first') + 1;
    x1 = find(lon_edge < min(day_lons(:)), 1, 'last') - 1;
    x2 = find(lon_edge > max(day_lons(:)), 1, 'first') + 1;
    
    for x=x1:x2
        for y=y1:y2
            lontest = lon_edge(x:x+1);
            lattest = lat_edge(y:y+1);
            pp = day_lons >= lontest(1) & day_lons < lontest(2) & day_lats >= lattest(1) & day_lats < lattest(2);
            include_bool(y,x) = sum(pp) >= 20;
        end
    end
    
    BEHR = load(fullfile(psm_behr_dir, behr_filename(u_dates{a})));
    
    [cvm_vcds_tmp, cvm_wts_tmp] = get_cvm_vcds(BEHR.OMI, include_bool);
    [psm_vcds_tmp, psm_wts_tmp] = get_psm_vcds(BEHR.Data, include_bool);
    cvm_vcds = cat(1, cvm_vcds, cvm_vcds_tmp(:));
    cvm_wts = cat(1, cvm_wts, cvm_wts_tmp(:));
    psm_vcds = cat(1, psm_vcds, psm_vcds_tmp(:));
    psm_wts = cat(1, psm_wts, psm_wts_tmp(:));
    
end

end

function [vcds, wts] = get_cvm_vcds(OMI, gridcell_bool)
vcds = [];
wts = [];

for a=1:numel(OMI)
    OMI(a) = omi_pixel_reject(OMI(a), 'omi', 0.2, 'XTrackFlags');
    xx = OMI(a).Areaweight > 0 & ~isnan(OMI(a).BEHRColumnAmountNO2Trop) & gridcell_bool;
    vcds = cat(1, vcds, OMI(a).BEHRColumnAmountNO2Trop(xx));
    wts = cat(1, wts, OMI(a).Areaweight(xx));
end
end

function [vcds, wts] = get_psm_vcds(Data, gridcell_bool)
    OMI = psm_wrapper(Data);
    vcds = OMI.BEHRColumnAmountNO2Trop(gridcell_bool);
    wts = OMI.Weights(gridcell_bool);
end

function [ lon_edge, lat_edge ] = compute_grid_edges(lon, lat)
if ~isvector(lon)
    lon = lon(1,:);
end
if ~isvector(lat)
    lat = lat(:,1)';
end

del_lon = mean(diff(lon));
del_lat = mean(diff(lat));

lon_edge = [lon - del_lon/2, lon(end) + del_lon/2];
lat_edge = [lat - del_lat/2, lat(end) + del_lat/2];
end