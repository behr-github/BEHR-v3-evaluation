% TOD & Site number
% Plot how vertical profiles changes for each of the DISCOVER sites by day,
% separating by time of day as well

matdir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

% Reading in the data and sorting by site
location = 'Baltimore'; datafield = 'NO2_LIF'; datafield2 = 'HNO3_TDLIF'; numberofsites = 8; subplotsize = [2,4];

% Some cell arrays to help make a structure later
sitecell = cell(1,numberofsites); for a=1:numberofsites; sitecell{a} = a; end
valcell = cell(1,4); for a=1:4; valcell{a} = []; end
altcell = cell(1,4); for a=1:4; altcell{a} = []; end
errcell = cell(1,4); for a=1:4; errcell{a} = []; end
rangecell = cell(1,4); for a=1:4; rangecell{a} = []; end

files = dir(fullfile(matdir,sprintf('%s*.mat',location)));
clear profiles
for a=1:numel(files);
    load(fullfile(matdir,files(a).name),'Merge');
    [no2, utc_no2, alt_no2, lon_no2, lat_no2, xx_no2] = remove_merge_fills(Merge,datafield);
    lt_no2 = Merge.Data.LOCAL_SUN_TIME.Values; lt_no2(xx_no2) = NaN;
    prof_no2 = Merge.Data.ProfileSequenceNum.Values; prof_no2(xx_no2) = NaN;
    u_prof_no2 = unique(prof_no2(prof_no2~=0 & ~isnan(prof_no2)));
    
    [data2, utc_data2, alt_data2, lon_data2, lat_data2, xx_data2] = remove_merge_fills(Merge,datafield2);
    lt_data2 = Merge.Data.LOCAL_SUN_TIME.Values; lt_data2(xx_data2) = NaN;
    prof_data2 = Merge.Data.ProfileSequenceNum.Values; prof_data2(xx_data2) = NaN;
    u_prof_data2 = unique(prof_data2(prof_data2~=0 & ~isnan(prof_data2)));
    
    tmpno2 = cell(numberofsites,24); tmpdata2 = cell(numberofsites,24);
    tmpalt = cell(numberofsites,24); tmpalt2 = cell(numberofsites,24);
    
    for b=1:numel(u_prof_no2)
        profnum_b = mod(u_prof_no2(b),1e5); % Some of the DISCOVER profile numbers are in the 100,000 range, others are in the 1000s range.  This gets to the important part of the number.
        profnum_b = floor(profnum_b/1000); % We want the number in the 1000s place, since that corresponds to the site

        % Sort profile numbers according to what hour (solar time) the
        % spiral started in
        pp = prof_no2 == u_prof_no2(b);
        no2_p = no2(pp); alt_p = alt_no2(pp); lt_p = lt_no2(pp); utc_p = utc_no2(pp);
        starthour = floor(min(lt_p(:)));
        tmpno2{profnum_b,starthour} = [tmpno2{profnum_b,starthour}; no2_p];
        tmpalt{profnum_b,starthour} = [tmpalt{profnum_b,starthour}; alt_p];
        
        pp2 = prof_data2 == u_prof_data2(b);
        data2_p = data2(pp2); alt2_p = alt_data2(pp2); lt2_p = lt_data2(pp2); utc2_p = utc_data2(pp2);
        starthour = floor(min(lt2_p(:)));
        tmpdata2{profnum_b,starthour} = [tmpdata2{profnum_b,starthour}; data2_p];
        tmpalt2{profnum_b,starthour} = [tmpalt2{profnum_b,starthour}; alt2_p];
    end
    
    % Initialize the data structure
    profiles_no2(a).date = Merge.metadata.date;
    profiles_no2(a).site = struct('sitenum',sitecell,'starthour',struct('range',rangecell,'vals',valcell,'alts',altcell,'errs',errcell));
    
    profiles_data2(a).date = Merge.metadata.date;
    profiles_data2(a).site = struct('sitenum',sitecell,'starthour',struct('range',rangecell,'vals',valcell,'alts',altcell,'errs',errcell));
    
    % Find the vertical profiles
    ltranges = [8:9; 10:11; 12:13; 14:15];
    for h=1:4;
        r = ltranges(h,:);
        
        for s=1:numberofsites
            % Data field 1
            no2_r = [tmpno2{s,r}]; alt_r = [tmpalt{s,r}];
            profiles_no2(a).site(s).starthour(h).range = r;
            if isempty(no2_r) || isempty(alt_r)
                profiles_no2(a).site(s).starthour(h).vals = [];
                profiles_no2(a).site(s).starthour(h).alts = [];
                profiles_no2(a).site(s).starthour(h).errs = [];
            else
                [val, alt, err] = bin_rolling_vertical_profile(alt_r,no2_r,0.5,0.1);
                profiles_no2(a).site(s).starthour(h).vals = val;
                profiles_no2(a).site(s).starthour(h).alts = alt;
                profiles_no2(a).site(s).starthour(h).errs = err;
            end
            
            % Data field 2
            data2_r = [tmpdata2{s,r}]; alt2_r = [tmpalt2{s,r}];
            profiles_data2(a).site(s).starthour(h).range = r;
            if isempty(data2_r) || isempty(alt2_r)
                profiles_data2(a).site(s).starthour(h).vals = [];
                profiles_data2(a).site(s).starthour(h).alts = [];
                profiles_data2(a).site(s).starthour(h).errs = [];
            else
                [val, alt, err] = bin_rolling_vertical_profile(alt2_r,data2_r,0.5,0.1);
                profiles_data2(a).site(s).starthour(h).vals = val;
                profiles_data2(a).site(s).starthour(h).alts = alt;
                profiles_data2(a).site(s).starthour(h).errs = err;
            end
        end
    end
end

% Plot the data; average any profiles from the same site/same
% hour, create plots that depict the variation by day in vertical profiles at each
% site at each hour


colors = {'b','r',[0 0.5 0],[0.5 0 0.5]};
%markers = {'s','^','o','v'};

for d=1:numel(profiles_no2)
    l = gobjects(4,2);
    figure('Position',[100 100,subplotsize(2)*400,subplotsize(1)*400]);
    for h=1:4;
        dates{d} = profiles_no2(d).date;
        for s=1:numberofsites
            subplot(subplotsize(1),subplotsize(2),s)
            title(sprintf('Site #%d',s));
            this_no2 = profiles_no2(d).site(s).starthour(h).vals;
            if ~isempty(this_no2)
                % To compare the profile shape of two pieces of data that can
                % have very different ranges, we'll normalize the data by
                % dividing by its highest value.  In the case of temperature
                % and potential temperature, we also subtract the minimum
                % value since both are so offset from 0.
                
                if any(strcmp(datafield,{'THETA','TEMPERATURE'}));
                    this_no2 = (this_no2 - min(this_no2(:)))/(max(this_no2(:))-min(this_no2(:)));
                else
                    this_no2 = this_no2/max(this_no2(:));
                end
                this_alt = profiles_no2(d).site(s).starthour(h).alts;
                l(h,1) = line(this_no2, this_alt, 'marker', 's','color',colors{h},'linewidth',1,'markersize',9,'markerfacecolor',colors{h});
            end
            
            this_data2 = profiles_data2(d).site(s).starthour(h).vals;
            if ~isempty(this_data2)
                if any(strcmp(datafield2,{'THETA','TEMPERATURE'}));
                    this_data2 = (this_data2 - min(this_data2(:)))/(max(this_data2(:))-min(this_data2(:)));
                else
                    this_data2 = this_data2/max(this_data2(:));
                end
                this_alt2 = profiles_data2(d).site(s).starthour(h).alts;
                l(h,2) = line(this_data2, this_alt2, 'marker', '^','color',colors{h},'linewidth',1,'markersize',9);
            end
            
            if h==4;
                xlabel('Normalized profile');
            end
        end
    end
    l2 = [l(:,1);l(:,2)];
    ll = find(l2>0); 
    datafieldname = regexprep(datafield,'[^a-zA-Z0-9]',' ');
    datafieldname2 = regexprep(datafield2,'[^a-zA-Z0-9]',' ');
    rangestr = {sprintf('%s LT 800-1000',datafieldname),sprintf('%s LT 1000-1200',datafieldname),sprintf('%s LT 1200-1400',datafieldname),sprintf('%s LT 1400-1600',datafieldname)...
        sprintf('%s LT 800-1000',datafieldname2),sprintf('%s LT 1000-1200',datafieldname2),sprintf('%s LT 1200-1400',datafieldname2),sprintf('%s LT 1400-1600',datafieldname2)};
    leg = legend(l2(ll),rangestr(ll));
    set(leg,'Units','normalized');
    legpos = get(leg,'Position'); newpos = [1 - legpos(3)*1.1, 1 - legpos(4)*1.1, legpos(3), legpos(4)];
    set(leg,'Position',newpos);
    
    suptitle(sprintf('TOD variation for %s %s and %s, on %s',location, datafieldname, datafieldname2, profiles_no2(d).date));
end