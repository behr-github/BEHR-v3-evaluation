% TOD & Site number
% Plot how vertical profiles changes for each of the DISCOVER sites by day,
% separating by time of day as well

matdir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files/Aircraft/';

% Reading in the data and sorting by site
location = 'Baltimore'; datafield = 'NO2_LIF'; numberofsites = 8; subplotsize = [2,4];

% Some cell arrays to help make a structure later
sitecell = cell(1,numberofsites); for a=1:numberofsites; sitecell{a} = a; end
valcell = cell(1,24); for a=1:24; valcell{a} = []; end
altcell = cell(1,24); for a=1:24; altcell{a} = []; end
errcell = cell(1,24); for a=1:24; errcell{a} = []; end

files = dir(fullfile(matdir,sprintf('%s*.mat',location)));
clear profiles
for a=1:numel(files);
    load(fullfile(matdir,files(a).name),'Merge');
    [no2, utc_no2, alt_no2, lon_no2, lat_no2, xx_no2] = remove_merge_fills(Merge,datafield);
    lt_no2 = Merge.Data.LOCAL_SUN_TIME.Values; lt_no2(xx_no2) = NaN;
    prof_no2 = Merge.Data.ProfileSequenceNum.Values; prof_no2(xx_no2) = NaN;
    u_prof_no2 = unique(prof_no2(prof_no2~=0 & ~isnan(prof_no2)));
    
    tmpno2 = cell(numberofsites,24);
    tmpalt = cell(numberofsites,24);
    
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
    end
    
    % Initialize the data structure
    profiles(a).date = Merge.metadata.date;
    profiles(a).site = struct('sitenum',sitecell,'starthour',struct('vals',valcell,'alts',altcell,'errs',errcell));
    
    % Find the vertical profiles
    for h=1:24;
        for s=1:numberofsites
            if isempty(tmpalt{s,h}) || isempty(tmpno2{s,h})
                profiles(a).site(s).starthour(h).vals = [];
                profiles(a).site(s).starthour(h).alts = [];
                profiles(a).site(s).starthour(h).errs = [];
            else
                [val, alt, err] = bin_rolling_vertical_profile(tmpalt{s,h},tmpno2{s,h},0.5,0.1);
                profiles(a).site(s).starthour(h).vals = val;
                profiles(a).site(s).starthour(h).alts = alt;
                profiles(a).site(s).starthour(h).errs = err;
            end
        end
    end
end

% Plot the data; average any profiles from the same site/same
% hour, create plots that depict the variation by day in vertical profiles at each
% site at each hour


% First figure out which hours actually have data:
% go through each hour of each day and test if there is data.  Once one day
% is found to have data for that hour, it will be plotted.
ptest = false(1,24);
for d=1:numel(profiles)
    for s=1:numberofsites
        for p=1:numel(profiles(d).site(s).starthour);
            ptest(p) = ~isempty(profiles(d).site(s).starthour(p).vals) || ptest(p);
        end
    end
end
xx = find(ptest);
n_xx = numel(xx);

colors = {'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],'b','r',[0 0.5 0],[0.5 0 0.5],[0.5 0.5 0.5],};
markers = {'s','s','s','s','s','^','^','^','^','^','o','o','o','o','o'};

for i=1:n_xx
    figure('Position',[100 100,subplotsize(2)*400,subplotsize(1)*400]);
    for d=1:numel(profiles)
        dates{d} = profiles(d).date;
        for s=1:numberofsites
            subplot(subplotsize(1),subplotsize(2),s)
            title(sprintf('Site #%d',s));
            this_no2 = profiles(d).site(s).starthour(xx(i)).vals;
            if ~isempty(this_no2)
                this_alt = profiles(d).site(s).starthour(xx(i)).alts;
                l(d,i) = line(this_no2, this_alt, 'marker', markers{d},'color',colors{d},'linewidth',2,'markersize',9);
            end
        end
    end
    ll = find(l(:,i)>0);
    leg = legend(l(ll,i),dates(ll));
    set(leg,'Units','normalized');
    legpos = get(leg,'Position'); newpos = [1 - legpos(3)*1.1, 1 - legpos(4)*1.1, legpos(3), legpos(4)];
    set(leg,'Position',newpos);
    datafieldname = regexprep(datafield,'[^a-zA-Z0-9]',' ');
    suptitle(sprintf('Daily variation: %s %s Profiles starting solar time %d',location, datafieldname, xx(i)));
end