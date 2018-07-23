% Makes a categorical scatter plot of the slopes for the Leitao data

addCampaigns = true;

% load the LeitaoTables variable
load('~/Documents/MATLAB/NO2 Profiles/Workspaces/LeitaoTables.mat');

figure;
ydat = table2array(LeitaoTables.Coincident(:,5));
line(ones(size(ydat)),ydat,'color','k','linestyle','none','marker','x','markersize',8,'linewidth',2);
ydat = table2array(LeitaoTables.Coincident(:,8));
line(2*ones(size(ydat)),ydat,'color','k','linestyle','none','marker','o','markersize',8,'linewidth',2);

ydat = table2array(LeitaoTables.AerosolAbove(:,5));
line(3*ones(size(ydat)),ydat,'color','k','linestyle','none','marker','x','markersize',8,'linewidth',2);
ydat = table2array(LeitaoTables.AerosolAbove(:,8));
line(4*ones(size(ydat)),ydat,'color','k','linestyle','none','marker','o','markersize',8,'linewidth',2);

ydat = table2array(LeitaoTables.NO2Above(:,5));
line(5*ones(size(ydat)),ydat,'color','k','linestyle','none','marker','x','markersize',8,'linewidth',2);
ydat = table2array(LeitaoTables.NO2Above(:,8));
line(6*ones(size(ydat)),ydat,'color','k','linestyle','none','marker','o','markersize',8,'linewidth',2);

set(gca,'fontsize',16);
set(gca,'xtick',[1.5,3.5,5.5],'xticklabel',{'Coincident','Aerosol above','NO_2 above'});
ylabel('Slope A/A_0 vs. AOD');

xlim([0 7]);
%set(gca,'ylimmode','manual');
lc=line(-1,-1,'color','k','linestyle','none','marker','x','markersize',8,'linewidth',2);
lf=line(-1,-1,'color','k','linestyle','none','marker','o','markersize',8,'linewidth',2);
legend([lc;lf],{'Coarse aerosol','Fine aerosol'})

campaign_data = 'cl';
switch campaign_data
    % These are the fits obtained from doing a y-residual fit on the plots
    % made with LIF and Chemilum. data respectively
    case 'lif'
        mdvals = [-23.56, 31.48, nan];
        mderrs = [151, 177, nan];
        cavals = [-150, nan, nan];
        caerrs = [231, nan, nan];
        txvals = [-27.35, -42.079, -677.97];
        txerrs = [217, 176, Inf];
        covals = [-1805, 84.58, nan];
        coerrs = [549, 424, nan];
    case 'cl'
        mdvals = [-29.30, -78.26, -433.72];
        mderrs = [111, 79, 237];
        cavals = [-229.5, nan, nan];
        caerrs = [184, nan, nan];
        txvals = [755.57, -58.16, nan];
        txerrs = [804, 229, nan];
        covals = [-774.53, 56.96, nan];
        coerrs = [472, 457, nan];
    case 'eye'
        mdvals = [-78, -143, NaN];
        mderrs = nan(1,3);
        cavals = [-456, NaN, NaN];
        caerrs = nan(1,3);
        txvals = [-120, -39, NaN];
        txerrs = nan(1,3);
    case 'lif-avgdenom'
        mdvals = [];
        mderrs = [];
        cavals = [];
        caerrs = [];
        txvals = [];
        txerrs = [];
        covals = [];
        coerrs = [];
end


if addCampaigns
    xvals = [1.5, 3.5, 5.5];
    lmd = line(xvals-0.3,mdvals,'color','b','marker','s','markersize',8,'linewidth',2,'linestyle','none');
    scatter_errorbars(xvals-0.3,mdvals,mderrs*2,'color','b','tipscale',2,'linewidth',2);
    lca = line(xvals-0.15,cavals,'color','r','marker','s','markersize',8,'linewidth',2,'linestyle','none');
    scatter_errorbars(xvals-0.15,cavals,caerrs*2,'color','r','tipscale',2,'linewidth',2);
    ltx = line(xvals+0.15,txvals,'color',[0 0.5 0],'marker','s','markersize',8,'linewidth',2,'linestyle','none');
    scatter_errorbars(xvals+0.15,txvals,txerrs*2,'color',[0 0.5 0],'tipscale',2,'linewidth',2);
    lco = line(xvals+0.3,covals,'color','m','marker','s','markersize',8,'linewidth',2,'linestyle','none');
    scatter_errorbars(xvals+0.3,covals,coerrs*2,'color','m','tipscale',2,'linewidth',2);
    
    legend([lc;lf;lmd;lca;ltx;lco],{'Coarse aerosol','Fine aerosol','DISCOVER-MD','DISCOVER-CA','DISCOVER-TX','DISCOVER-CO'});
end