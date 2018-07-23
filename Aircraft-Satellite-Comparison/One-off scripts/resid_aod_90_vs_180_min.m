function [  ] = resid_aod_90_vs_180_min( Comp90min, Comp180min, campaign_name )
%resid_aod_90_vs_180_min One off function to make plots to compare +/- 90
%vs +/- 180 min for residual vs AOD plots

f(1) = figure;
l180=line(Comp180min.coinc.aod,Comp180min.coinc.reldiff,'linestyle','none','marker','o','markersize',8,'color','b');
l90=line(Comp90min.coinc.aod,Comp90min.coinc.reldiff,'linestyle','none','marker','o','markersize',8,'color','b','markerfacecolor','b');
set(gca,'fontsize',16);
xlabel('AOD');
ylabel('% Difference VCDs');
legend([l90;l180], {'\pm 1.5 hr of overpass','\pm 3 hr of overpass'});
title(sprintf('%s coincident layers',upper(campaign_name)));

f(2) = figure;
l180=line(Comp180min.aerabove.aod,Comp180min.aerabove.reldiff,'linestyle','none','marker','o','markersize',8,'color','b');
l90=line(Comp90min.aerabove.aod,Comp90min.aerabove.reldiff,'linestyle','none','marker','o','markersize',8,'color','b','markerfacecolor','b');
set(gca,'fontsize',16);
xlabel('AOD');
ylabel('% Difference VCDs');
legend([l90;l180], {'\pm 1.5 hr of overpass','\pm 3 hr of overpass'});
title(sprintf('%s aerosol layer above',upper(campaign_name)));

f(3) = figure;
l180=line(Comp180min.no2above.aod,Comp180min.no2above.reldiff,'linestyle','none','marker','o','markersize',8,'color','b');
l90=line(Comp90min.no2above.aod,Comp90min.no2above.reldiff,'linestyle','none','marker','o','markersize',8,'color','b','markerfacecolor','b');
set(gca,'fontsize',16);
xlabel('AOD');
ylabel('% Difference VCDs');
legend([l90;l180], {'\pm 1.5 hr of overpass','\pm 3 hr of overpass'});
title(sprintf('%s NO_2 layer above',upper(campaign_name)));

for a=1:3;
    X(a,1:2) = f(a).Children(2).XLim;
    Y(a,1:2) = f(a).Children(2).YLim;
end

xmax = max(X(:,2));
ymax = max(abs(Y(:)));

for a=1:3
    f(a).Children(2).XLim = [0 xmax];
    f(a).Children(2).YLim = [-ymax, ymax];
end

end

