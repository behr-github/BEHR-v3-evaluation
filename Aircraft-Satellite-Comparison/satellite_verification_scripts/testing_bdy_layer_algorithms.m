%Testing various methods of finding all boundary layer crossings in
%ARCTAS-CA data, or other aircraft datasets in which the boundary layer is
%found by looking for sharp gradients in chemical concentrations

%% Loads data for an ARCTAS day
clear all
% Valid dates are 6/18/2008, 6/20/2008, 6/22/2008, 6/24/2008
file_date = '6/20/2008';
% Path to ARCTAS .mat files
file_path = '/Volumes/share/GROUP/ARCTAS/Matlab Files/';

file_date_str = datestr(file_date,29);
year = file_date_str(1:4); month = file_date_str(6:7); day = file_date_str(9:10);
filename = sprintf('ARCTAS_CA_%s_%s_%s.mat',year,month,day);
load(fullfile(file_path,filename),'Merge'); %Only load the Merge data structure

% Read altitude and utc data from the merge structure; assume only replace
% with NANs for fill values in altitude
utc_raw = Merge.Data.UTC.Values;
alt_raw = Merge.Data.ALTP.Values;

% Read NO2 data from the Merge structure, removing fill and LOD values
fill = Merge.Data.NO2_UCB.Fill; ulod = Merge.metadata.upper_lod_flag; llod = Merge.metadata.lower_lod_flag;
yy_no2 = Merge.Data.NO2_UCB.Values ~= fill & Merge.Data.NO2_UCB.Values ~= ulod & Merge.Data.NO2_UCB.Values ~= llod;
no2 = Merge.Data.NO2_UCB.Values; no2(~yy_no2) = NaN;
utc_no2 = Merge.Data.UTC.Values; utc_no2(~yy_no2) = NaN;
alt_no2 = Merge.Data.ALTP.Values; alt_no2(~yy_no2) = NaN;

% Read water data, doing the same
water_fill = Merge.Data.H2Ov.Fill;
yy_h2o = Merge.Data.H2Ov.Values ~= fill & Merge.Data.H2Ov.Values ~= ulod & Merge.Data.H2Ov.Values ~= llod;
h2o = Merge.Data.H2Ov.Values; h2o(~yy_h2o) = NaN;
utc_h2o = Merge.Data.UTC.Values; utc_h2o(~yy_h2o) = NaN;
alt_h2o = Merge.Data.ALTP.Values; alt_h2o(~yy_h2o) = NaN;

fprintf('Loading complete.\n');
%% Method 1: Find all points with a certain NO2 gradient that also change altitude
fprintf('Executing method 1...   ');

% Calcuate dNO2/dt and dz/dt.  Look for places where both exceed given
% criteria.

dno2 = diff(no2) ./ diff(utc_no2);
dalt = diff(alt_no2) ./ diff(utc_no2);

xx_no2 = dno2 > 1e3;
xx_alt = dalt > 0.01;
xx = xx_no2 & xx_alt;

% Plot NO2 and altitude vs. time and mark all the sites IDed as boundary
% layer crossings
fignum = max(findall(0,'Type','Figure'))+1; if isempty(fignum); fignum=1; end
figure(fignum);
[hax1, h1a, h1b] = plotyy(utc_no2,no2,utc_no2,alt_no2); line(utc_no2(xx),alt_no2(xx),'Parent',hax1(2),'linestyle','none','marker','^','markersize',16,'color','r');
fprintf('Method 1 complete\n');

%% Method 2: Locate parts of the flight with consistent altitude changes and find the largest dNO2 in each range.

% Find consistent changes in altitude, spanning more than 500 m.
firstns=find(~isnan(no2),3,'first');
n = firstns(3); last_n = fliplr(firstns(1:2));
ranges = [];
while n <= numel(utc_no2)
    start_n = n;
    while true && n <= numel(utc_no2)
        if isnan(no2(n)); n = n+1; continue; %If the value is a NaN, we'll need to skip this index
        elseif sign(alt_no2(n) - alt_no2(last_n(1))) == sign(alt_no2(last_n(1)) - alt_no2(last_n(2)));
        elseif sign(alt_no2(n) - alt_no2(last_n(1))) == 0;
        elseif abs(alt_no2(n) - alt_no2(last_n(1))) < 0.1*(alt_no2(last_n(1)) - alt_no2(last_n(2)))
        else
            if abs(alt_no2(start_n) - alt_no2(last_n(1))) > 0.5;
                ranges = [ranges; start_n, last_n(1)];%#ok<AGROW>
            end
            break;
        end
        
        last_n = [n, last_n(1)];
        n = n+1;
    end
    last_n = [n, last_n(1)];
    n = n+1;
end
fignum = max(findall(0,'Type','Figure'))+1; if isempty(fignum); fignum=1; end
figure(fignum);
[hax2, h2a, h2b] = plotyy(utc_no2,no2,utc_no2,alt_no2);
for a=1:size(ranges,1)
    line(utc_no2(ranges(a,:)),alt_no2(ranges(a,:)),'Parent',hax2(2),'linestyle','-','linewidth',4,'color','r','marker','o','markersize',12);
end
fprintf('Method 2 complete\n');

%% Method 3: Still trying to find contiguous changes in altitude
dalt = diff(alt_raw);
dalt_cont = false(size(dalt));
for a=2:numel(dalt_cont)
    if sign(dalt(a)) == sign(dalt(a-1));
        dalt_cont(a) = 1;
    end
end

dalt_zeros = find(dalt_cont==0);
for a=1:numel(dalt_zeros);
    if dalt_zeros(a) > 3 && dalt_zeros(a) < (numel(dalt_zeros)-3)
        if all(dalt_cont(dalt_zeros(a)-3:dalt_zeros(a)-1)==1) && all(dalt_cont(dalt_zeros(a)+1:dalt_zeros(a)+3)==1)
            dalt_cont(dalt_zeros(a)) = 1;
        end
    end
end
fprintf('Method 3 complete\n');

%% Method 4: I hate teaching computers to find semi-contiguous regions: Another altitude attempt
% Find all points where altitude is changing by more than 5 m/s
dalt = diff(alt_raw)./diff(utc_raw);
gtcrit = abs(dalt) > 0.005;

% Remove any singletons up to oct-letons, as these likely do not represent
% actual changes in airplane altitude for the purpose of measuring the BL
% height.

runs = contiguous(gtcrit,1);
ranges = runs{2};
gtcrit2 = gtcrit;
for a=1:size(ranges,1)
    if ranges(a,2) - ranges(a,1) < 8;
        gtcrit2(ranges(a,1):ranges(a,2)) = 0;
    end
end

% Now fill in any gaps of 300 s or less
runs = contiguous(gtcrit2,1);
ranges = runs{2};
gtcrit3 = gtcrit2;
for a=2:size(ranges,1)
    if utc_raw(ranges(a,1)) - utc_raw(ranges(a-1,2)) <= 300
        gtcrit3(ranges(a-1,2):ranges(a,1)) = 1;
    end
end

% Alternately, for the gaps, check if >75% of the differences there have
% the same sign as the bounding sections
runs = contiguous(gtcrit2,1);
ranges = runs{2};
gtcrit4 = gtcrit2;
for a=2:size(ranges,1)
    % First check that the surrounding ranges have the same median sign. If
    % not, skip this range.
    if sign(median(dalt(ranges(a-1,:)))) == sign(median(dalt(ranges(a,:))))
        % Next check that >75% of the measurements in the gap have the same
        % sign as the surrounding segments
        boundary_sign = sign(median(dalt(ranges(a-1,:))));
        tmp = sign(dalt(ranges(a-1,2):ranges(a,1)));
        gap_signs_logical = tmp == boundary_sign;
        sign_percent = sum(gap_signs_logical)/numel(gap_signs_logical);
        sign_zeros_percent = (sum(gap_signs_logical)+sum(tmp==0))/numel(gap_signs_logical);
        if sign_percent >= 0.75 || sign_zeros_percent >= 0.9
            gtcrit4(ranges(a-1,2):ranges(a,1)) = 1;
        end
    end
end

runs = contiguous(gtcrit4,1);
ranges = runs{2};
gtcrit5 = gtcrit4;
for a=1:size(ranges,1)
    if ranges(a,2) - ranges(a,1) < 25;
        gtcrit5(ranges(a,1):ranges(a,2)) = 0;
    end
end

fprintf('Method 4 complete\n');

%% Method 5: Running average of altitude changes
% Window size for the filter
ws = 200;

dalt = diff(alt_raw);
dalt_avg = filter(ones(1,ws)/ws,1,dalt);
gtcrit = abs(dalt) > 0.005;

% Check if >75% of the differences in the gaps have
% the same sign as the bounding sections
runs = contiguous(gtcrit,1);
ranges = runs{2};
gtcrit2 = gtcrit;
for a=2:size(ranges,1)
    % First, check that the range is not too long (>900 s)
    if utc_raw(ranges(a,1)) - utc_raw(ranges(a-1,2)) < 901;
        % Second, check that the surrounding ranges have the same median sign. If
        % not, skip this range.
        if sign(median(dalt(ranges(a-1,:)))) == sign(median(dalt(ranges(a,:))))
            % Next check that >75% of the measurements in the gap have the same
            % sign as the surrounding segments
            boundary_sign = sign(median(dalt(ranges(a-1,:))));
            tmp = sign(dalt(ranges(a-1,2):ranges(a,1)));
            gap_signs_logical = tmp == boundary_sign;
            sign_percent = sum(gap_signs_logical)/numel(gap_signs_logical);
            sign_zeros_percent = (sum(gap_signs_logical)+sum(tmp==0))/numel(gap_signs_logical);
            if sign_percent >= 0.75 || sign_zeros_percent >= 0.9
                gtcrit2(ranges(a-1,2):ranges(a,1)) = 1;
            end
        end
    end
end

fprintf('Method 5 complete\n');

%% Method 6: Use a GUI to manually find altitude ranges, then look in each of those ranges to find the greatest dNO2

% Call the GUI
addpath('/Users/Josh/Documents/MATLAB/NO2 Profiles/GUIs')
[ranges,~,~] = select_changing_altitude(fullfile(file_path,filename));
s=size(ranges);
blh = -99*ones(s(1),1);
utc = -99*ones(s(1),1);
maxes = -99*ones(s(1),1);
for a=1:s(1)
    startind = find(abs(utc_raw - ranges(a,1)) == min(abs(utc_raw - ranges(a,1))));
    endind = find(abs(utc_raw - ranges(a,2)) == min(abs(utc_raw - ranges(a,2))));
    dNO2 = abs(diff(no2(startind:endind)));
    maxes(a) = max(dNO2);
    if max(dNO2) > 800; 
        blh_ind = find(dNO2 == max(dNO2)); %Require that the change be greater than twice the average UT concentration
        blh(a) = alt_no2(startind+blh_ind);
        utc(a) = utc_no2(startind+blh_ind);
    end
end

hax = plotyy(utc_no2,no2,utc_raw,alt_raw);
hold on
line(utc,blh,'Parent',hax(2),'linestyle','none','marker','^','markersize',16,'color','r');
hold off

%% Method 7: Use the GUI plus the find_bdy_layer_height function
addpath('/Users/Josh/Documents/MATLAB/NO2 Profiles/GUIs')
[ranges,file,h] = select_changing_altitude(fullfile(file_path,filename));
s=size(ranges);
blh = -99*ones(s(1),1);
blh_water = -99*ones(s(1),1);
utc = -99*ones(s(1),1);
utc_water = -99*ones(s(1),1);
for a=1:s(1)
    startind = find(abs(utc_raw - ranges(a,1)) == min(abs(utc_raw - ranges(a,1))));
    endind = find(abs(utc_raw - ranges(a,2)) == min(abs(utc_raw - ranges(a,2))));
    no2seg = no2(startind:endind);
    altseg = alt_no2(startind:endind);
    utc(a) = nanmedian(utc_no2(startind:endind));
    
    h2oseg = h2o(startind:endind);
    alth2oseg = alt_h2o(startind:endind);
    utc_water(a) = nanmedian(utc_h2o(startind:endind));
    % if using the exponential boundary layer model, we need to bin the
    % data
    [no2seg, altseg] = bin_vertical_profile(altseg,no2seg,0.5);
    [h2oseg, alth2oseg] = bin_vertical_profile(alth2oseg,h2oseg,0.5);
    blh(a) = find_bdy_layer_height(no2seg, altseg,'exp2','debug',1);
    title(sprintf('Max NO2 - UTC: %.0f - %.0f',ranges(a,1),ranges(a,2)),'fontsize',16)
    %blh_water(a) = find_bdy_layer_height(h2oseg, alth2oseg,'exp2','debug',1);
    %title(sprintf('Water - UTC: %.0f - %.0f',ranges(a,1),ranges(a,2)),'fontsize',16)
end
nfig = nextfig; figure(nfig);
hax = plotyy(utc_no2,no2,utc_raw,alt_raw);
hold on
l1 = line(utc,blh,'Parent',hax(2),'linestyle','none','marker','^','markersize',16,'color','r');
l2 = line(utc,blh_water,'Parent',hax(2),'linestyle','none','marker','o','markersize',16,'color','c');
hold off
