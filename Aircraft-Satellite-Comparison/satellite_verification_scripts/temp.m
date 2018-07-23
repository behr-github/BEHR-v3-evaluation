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
    title(sprintf('NO2 - UTC: %.0f - %.0f',ranges(a,1),ranges(a,2)),'fontsize',16)
    blh_water(a) = find_bdy_layer_height(h2oseg, alth2oseg,'exp2','debug',1);
    title(sprintf('Water - UTC: %.0f - %.0f',ranges(a,1),ranges(a,2)),'fontsize',16)
end
nfig = nextfig; figure(nfig);
hax = plotyy(utc_no2,no2,utc_raw,alt_raw);
hold on
l1 = line(utc,blh,'Parent',hax(2),'linestyle','none','marker','^','markersize',16,'color','r');
l2 = line(utc,blh_water,'Parent',hax(2),'linestyle','none','marker','o','markersize',16,'color','c');
hold off