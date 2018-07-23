% show_bl_profiles.m
%
% Finds BLH heights and plots the part of the flight path that they were
% determined from, along with the vertical profiles of NO2, H2O, and theta.
%
%   Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2014

date = '622'; %month and day (mdd) with no separators
choice_ind = 6;
range_ind = 3;

load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/Scratch/ARCTAS Verification -  AMF using OMI cldfrac.mat');
load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/Scratch/BLH_choices modified 3.mat');
load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/ARCTAS-CA Altitude Ranges Exclusive 2.mat');
Merge = eval(sprintf('Merge%s',date));

[no2, utc, alt, lon, lat] = remove_merge_fills(Merge,'NO2_UCB');
[h2o, utc_h2o, alt_h2o] = remove_merge_fills(Merge,'H2Ov');
[theta, utc_theta, alt_theta] = remove_merge_fills(Merge,'THETA');

[heights_no2, times, ~] = findall_no2_bdy_layer_heights(utc, no2, alt, Ranges(range_ind).Ranges,blh_choices(choice_ind).no2);
[heights_h2o, ~, ~] = findall_no2_bdy_layer_heights(utc_h2o, h2o, alt_h2o, Ranges(range_ind).Ranges,blh_choices(choice_ind).h2o);
[heights_theta, ~, ~] = findall_no2_bdy_layer_heights(utc_theta, theta, alt_theta, Ranges(range_ind).Ranges,blh_choices(choice_ind).theta,'blmode','theta');

tmp = [heights_no2, heights_h2o, heights_theta];
heights = nanmean(tmp,2);
nans = isnan(heights);
heights = heights(~nans); times = times(~nans);
if numel(heights) < 2
    error('bdy_layer_verify:findBoundaryLayer','Could not find 2 boundary layer heights - needed for interpolation')
end

interp_height = interp1(times, heights, utc); % Linearly interpolate the boundary layer height to every value of UTC.

% For values outside of the range of "times," assume that the boundary
% layer height equals the closest value.
first_height = find(~isnan(interp_height),1,'first');
last_height = find(~isnan(interp_height),1,'last');

if first_height > 1;
    interp_height(1:first_height-1) = interp_height(first_height);
end
if last_height < numel(interp_height)
    interp_height(last_height+1:end) = interp_height(last_height);
end

utcstart = local2utc('12:00', 'pst'); utcend = local2utc('15:00', 'pst');
time_logical = utc >= utcstart & utc <= utcend;
utc_intime = utc(time_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot the vertical profiles with their location %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We'll be positioning plots, so we need the screen size
tmp = get(0,'screensize'); monx = tmp(3); mony = tmp(4);

% Plot the boundary layer heights on a map
state_outlines('ca'); map_fig = gcf; map_axes = gca;
set(map_fig,'Position',[monx - mony*0.44, mony*0.05, mony*0.4, mony*0.4]);
line(lon,lat); xlim([floor(min(lon)), ceil(max(lon))]); ylim([floor(min(lat)), ceil(max(lat))]);

% Plot the plane altitude against time
alt_fig = figure; alt_axes = gca; 
set(alt_fig,'Position',[monx - mony*0.044, mony*0.95, monx*0.33, mony*0.4]);
line(utc,alt); alt_ylim = get(alt_axes,'ylim'); y1=alt_ylim(1); y2 = alt_ylim(2);
line(utc,interp_height,'color','r');
hold on; fill([utcstart, utcstart, utcend, utcend],[y1,y2,y2,y1],'r','FaceAlpha',0.3); hold off

% Now, for each range, mark the location on the map where the profile comes
% from, highlight the range of times in the altitude plot, and plot the
% vertical profiles of NO2, H2O and potential temperature 
r = Ranges(range_ind).Ranges;
for a=1:size(r,1)
    xx = utc >= r(a,1) & utc <= r(a,2);
    if a==1; pre_xx = utc < r(a,1);
    else pre_xx = utc > r(a-1,2) & utc < r(a,1);
    end
    
    if a==size(r,1); post_xx = utc > r(a,2);
    else post_xx = utc > r(a,2) & utc < r(a+1,1);
    end
    
    % First the map   
    map_line_in = line(lon(xx), lat(xx), 'color',[0 0.7 0],'linewidth',2,'Parent',map_axes);
    map_line_pre = line(lon(pre_xx), lat(pre_xx), 'color','m','linewidth',2,'Parent',map_axes);
    map_line_post = line(lon(post_xx), lat(post_xx), 'color','r','linewidth',2,'Parent',map_axes);
    
    % Draw a second map showing the BLH for the segment before, during, and
    % after the range.
    state_outlines('ca'); map2_fig = gcf; 
    set(map2_fig,'Position',[monx - mony*0.88, mony*0.05, mony*0.4, mony*0.4]);
    xlim([floor(min(lon)), ceil(max(lon))]); ylim([floor(min(lat)), ceil(max(lat))]);
    scatter(lon(pre_xx | post_xx), lat(pre_xx | post_xx), 8, interp_height(pre_xx | post_xx));
    scatter(lon(xx), lat(xx), 16, interp_height(xx));
    cb = colorbar; ylabel(cb,'BL Height (km)','fontsize',14);
    
    % Now the altitude, we'll draw lines matching those on the first map
    alt_line_in = line(utc(xx), alt(xx), 'color',[0 0.7 0],'linewidth',4,'Parent',alt_axes);
    alt_line_pre = line(utc(pre_xx), alt(pre_xx), 'color','m','linewidth',4,'Parent',alt_axes);
    alt_line_post = line(utc(post_xx), alt(post_xx), 'color','r','linewidth',4,'Parent',alt_axes); 
    
    % Finally the profiles.
    % NO2
    [no2seg, altseg] = bin_rolling_vertical_profile(alt(xx),no2(xx),0.5,0.1);
    no2height = find_bdy_layer_height(no2seg, altseg,'max');
    no2_fig = figure; set(no2_fig,'Position',[0.01*monx,0.5*mony,0.45*mony,0.45*mony]);
    plot(no2seg,altseg,'marker','o','linewidth',2,'color',[0 0.5 0]); 
    line([min(no2seg),max(no2seg)],[no2height, no2height],'linestyle',':','color','k','linewidth',2);
    title(sprintf('NO2 #%d - %s',a,blh_choices(choice_ind).no2{a}),'fontsize',20);
    
    %H2O
    [h2oseg, altseg] = bin_rolling_vertical_profile(alt(xx),h2o(xx),0.5,0.1);
    h2oheight = find_bdy_layer_height(h2oseg, altseg,'max');
    h2o_fig = figure; set(h2o_fig,'Position',[0.3*monx,0.5*mony,0.45*mony,0.45*mony]);
    plot(h2oseg,altseg,'marker','o','linewidth',2); 
    line([min(h2oseg),max(h2oseg)],[h2oheight, h2oheight],'linestyle',':','color','k','linewidth',2);
    title(sprintf('H2O #%d - %s',a,blh_choices(choice_ind).h2o{a}),'fontsize',20);
    
    %Theta
    [thseg, altseg] = bin_rolling_vertical_profile(alt(xx),theta(xx),0.5,0.1);
    thheight = find_bdy_layer_height(thseg, altseg,'theta');
    th_fig = figure; set(th_fig,'Position',[0.15*monx,0.01*mony,0.45*mony,0.45*mony]);
    plot(thseg,altseg,'marker','o','linewidth',2,'color','r'); 
    line([min(thseg),max(thseg)],[thheight, thheight],'linestyle',':','color','k','linewidth',2);
    title(sprintf('THETA #%d - %s',a,blh_choices(choice_ind).theta{a}),'fontsize',20);
    
    fprintf('Pause - press a key\n'); pause
    delete(alt_line_in); delete(alt_line_post); delete(alt_line_pre);
    delete(map_line_in); delete(map_line_post); delete(map_line_pre);
    close(no2_fig); close(h2o_fig); close(th_fig); close(map2_fig);
end
close(map_fig); close(alt_fig);
