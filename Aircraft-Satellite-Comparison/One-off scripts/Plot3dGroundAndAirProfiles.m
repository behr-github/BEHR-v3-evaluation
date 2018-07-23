target_site = 5;
ground7_1 = SuperMerge.Date_20110701.Site_05;

[no2,utc,pres,lon,lat] = remove_merge_fills(Merge,'NO2_LIF','alt','PRESSURE');
sitenum = Merge.Data.discoveraqSiteFlag1sec.Values;

zz = ground7_1.Stop_UTC.Values > min(utc) & ground7_1.UTC.Values < max(utc);
gutc = ground7_1.UTC.Values(zz); gutcstop = ground7_1.Stop_UTC.Values(zz); gno2 = ground7_1.NO2.Values(zz);
gutcmid = ground7_1.Mid_UTC.Values(zz);
gno2(gno2<0) = NaN; gno2 = gno2*1e3;

no2vec = []; altvec = []; timevec = []; colvec = [];
for a=1:numel(gno2)
    no2vec(end+1) = gno2(a); altvec(end+1) = 1030; timevec(end+1) = gutcmid(a); colvec(end+1) = 0;
    xx = utc > gutc(a) & utc < gutcstop(a) & sitenum == target_site;
    profno2 = no2(xx); profpres = pres(xx);
    [no2bins, presbins] = bin_omisp_pressure(profpres, profno2);
    no2vec = [no2vec,no2bins]; altvec = [altvec, presbins];
    timevec = [timevec, gutcmid(a)*ones(size(no2bins))];
    colvec = [colvec, ones(size(no2bins))];
end

figure;
scatter3(no2vec,timevec,altvec,40,colvec)
set(gca,'zdir','reverse')
xlabel('NO2 (ppt)','fontsize',14)
ylabel('Time (UTC)','fontsize',14)
zlabel('Pressure (hPa)','fontsize',14)
title(['Difference in {\color{red}aircraft} and {\color{blue}ground} NO2 measurement at MD Site #',num2str(target_site)],'fontsize',14)
caxis([-0.2 1.2])
