function [  ] = show_variability( airno2_iall, lat_iall, lon_iall, dates_iall, db_iall, campaign_name, ptskip, use_cl_instead )
%show_variability One-off function to make plots to show variability in aircraft data
%   Comparing OMI and BEHR columns to aircraft-derived column from
%   DISCOVER-MD shows an alarming lack of correlation. I expect that this
%   is because there's a good bit of small scale variation in NO2
%   concentration that gets averaged together in satellite pixels but that
%   the aircraft can resolve.
%
%   This function is meant to produce plots that show that. It needs output
%   variables from Run_Spiral_Verification and the campaign name to be
%   passed to merge_field_names.  Campaign name should be passed as, e.g.,
%   "discover-md" or "discover_md" and not "discovermd" since this will be
%   used to load the merge files.
%
%   Josh Laughner <joshlaugh5@gmail.com> 24 Mar 2015

if nargin < 7
    ptskip = 1;
end

if nargin < 8
    use_cl_instead = false;
end

if nargin < 9;
    tosave = false;
end

if tosave
    savedir = uigetdir('~/MATLAB/Figures','Choose where to save the figures to');
    if savedir == 0
        E.userCancel;
    end
end

scrsz = get(0,'screensize') + [0 80 0 -80];

[Names,~,file_dir] = merge_field_names(campaign_name);

dates_to_load = unique(dates_iall);
datenum_vec = datenum(dates_iall);
merge_name = [upper(regexprep(campaign_name,'\W','_')),'_P3_1sec_%s_%s_%s.mat'];

% The limits the colorbar should have
collim = [0 1e16];

for a=1:numel(dates_to_load)
    yr = dates_to_load{a}(1:4);
    mn = dates_to_load{a}(6:7);
    dy = dates_to_load{a}(9:10);
    fname = sprintf(merge_name,yr,mn,dy);
    M = load(fullfile(file_dir,fname),'Merge');
    Merge = M.Merge;
    
    if ~use_cl_instead
        no2 = remove_merge_fills(Merge,Names.no2_lif);
    else
        no2 = remove_merge_fills(Merge,Names.no2_ncar);
    end
    [profnums,utc,alt,lon,lat] = remove_merge_fills(Merge,Names.profile_numbers,'alt',Names.gps_alt);
    tz = round(mean(lon)/15);
    
    % Get the data for this day
    dd = datenum_vec == datenum(dates_to_load{a});
    airno2_dd = airno2_iall(dd);
    airlat_dd = lat_iall(dd);
    airlon_dd = lon_iall(dd);
    allbehr_dd = db_iall.all_behr(dd);
    latcorn_dd = db_iall.latcorn(dd);
    loncorn_dd = db_iall.loncorn(dd);
    profnum_dd = cell2mat(db_iall.profnums(dd));
    
    site_count_dd = mod(profnum_dd,1000);
    u_sc_dd = unique(site_count_dd);
    
    nd = numel(u_sc_dd);
    
    for d=1:nd
        % Each set of profiles will have a 2-panel subplot. On the left each
        % colocated satellite pixel will be plotted and the column density
        % of the aircraft column will be plotted. On the right, the
        % satellite pixels will be plotted as outlines only, and the mixing
        % ratio at each point along the path will be plotted. The size of
        % the points should correspond to their altitude.
        %
        % The "set" that we want to plot is a collection of profiles
        % occurring at unique sites.  Profile numbers for MD have the form
        % snnn, where s is the site number and nnn is a sequential ID of
        % the profiles at that site. CA and TX are 10snnn, so in either
        % case this line will identify the site.
        xx = site_count_dd == u_sc_dd(d);
        ff = find(xx);
        lsts = {};
        for f = 1:numel(ff)
            pp = profnums == profnum_dd(ff(f));
            lsts{f} = [utc2local(mean(utc(pp)),tz),' '];
        end
        lsts = cat(2,lsts{:});
        
        % Make the plots
        figure;
        subplot(1,2,1);
        
        shandle = PlotAircraftColDen(airno2_dd(xx),airlon_dd(xx),airlat_dd(xx),collim);
        
        subplot(1,2,1);
        set(gca,'fontsize',16);
        title(sprintf('%s: Profiles %s\n(LST: %s)\nColumn density',dates_to_load{a},mat2str(profnum_dd(xx)),lsts));
        
        subplot(1,2,2);
        PlotAircraftData(no2,lon,lat,alt,profnums,profnum_dd(xx),'[NO_2] (pptv)',ptskip);
        set(gca,'fontsize',16);
        title(sprintf('%s: Profiles %s\n(LST: %s)\nNO2/pptv',dates_to_load{a},mat2str(profnum_dd(xx)),lsts));
        
        
        this_behr = cat(1,allbehr_dd{xx});
        this_loncorn = cat(2,loncorn_dd{xx});
        this_latcorn = cat(2,latcorn_dd{xx});
        DrawSatPix(this_behr,this_loncorn,this_latcorn,collim);
        % Make sure that the scatter of aircraft column density is
        % drawn on top of the patches of satellite pixels (ensure
        % visibility)
        uistack(shandle,'top');
        
        set(gcf,'position',scrsz);
        
        if tosave
            savename = sprintf('%s variability - %s Profs %s',upper(regexprep(campaign_name,'[-_]',' ')),regexprep(dates_to_load{a},'\W','-'),regexprep(mat2str(profnum_dd(xx)),'[\[\]]',''));
            savepath = fullfile(savedir,savename);
            savefig(savepath);
        else
            fprintf('Paused...\n');
            pause
        end
        close(gcf);
    end
end


end

function s = PlotAircraftColDen(no2,lon,lat,clim)
    s = scatter(lon,lat,384,no2,'filled','marker','s','MarkerEdgeColor','k','linewidth',3);
    colorbar;
    caxis(clim);
end

function PlotAircraftData(data_in,lon,lat,alt,raw_profnums,matched_profnums,cblabel,ptskip)
    % First we want to find the NO2 data that corresponds to each of the
    % profiles given, plus the straight flight before and after.
    B = findBlock(raw_profnums);
    xx = false(size(data_in));
    for a=1:numel(matched_profnums)
        i = find(B(:,3) == matched_profnums(a));
        xx(B(i-1,1):B(i+1,2)) = true;
    end
    
    % Only plot every ptskip points
    ss = mod(1:numel(xx),ptskip) == 0;
    xx = xx & ss;
    
    % Then we plot, representing altitude by size and NO2 concentration by
    % color
    scatter3(lon(xx),lat(xx),alt(xx),12,data_in(xx));
    cb=colorbar;
    set(cb.Label,'String',cblabel,'fontsize',12);
end

function DrawSatPix(no2,loncorn,latcorn,clim)
    % Get the current colormap. We're going to manually assign colors to
    % the pixels, so we'll create a vector that linearly maps values to
    % colors and interpolate 
    cmap = colormap;
    cvec = (linspace(clim(1),clim(2),length(cmap)))';
    
    % Clamp the value to within the colormap range
    no2 = max(min(no2,clim(2)),clim(1));
    
    for a=1:numel(no2)
        col = zeros(1,3);
        col(1) = interp1(cvec,cmap(:,1),no2(a));
        col(2) = interp1(cvec,cmap(:,2),no2(a));
        col(3) = interp1(cvec,cmap(:,3),no2(a));
        
        subplot(1,2,1);
        phandle = patch(loncorn(:,a),latcorn(:,a),'w','edgecolor','k','facealpha',0.5,'linewidth',2);
        set(phandle,'facecolor',col)
        
        % Plotting in 3D space
        subplot(1,2,2);
        patch(loncorn(:,a),latcorn(:,a),zeros(4,1),'w','edgecolor',col,'linewidth',2);
        
    end
end

