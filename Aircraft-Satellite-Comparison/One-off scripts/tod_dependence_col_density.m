% Sorts aircraft columns into bins based on their start times
% Run Run_Spiral_Verification first to get the variables needed.

% Campaign name, for the titles
campaign_name = 'Discover-CA';
uncert = true;
% Convert the output start times to numerical values. These are already in
% local times.
sec_after_midnight = cell2mat(db_iall.start_times);
timevec = sec_after_midnight / 3600;

% Define the bin edges as hour wide bins
bin_edges = [(0:23)',(1:24)'];
air_cols = cell(1,size(bin_edges,1));

for a=1:numel(timevec)
    if ~isnan(timevec(a))
        xx = timevec(a) >= bin_edges(:,1) & timevec(a) < bin_edges(:,2);
        air_cols{xx} = [air_cols{xx}; airno2_iall(a)];
    end
end

% Take the mean and median of each bin and plot

air_means = nan(size(air_cols));
air_std = nan(size(air_cols));
air_meds = nan(size(air_cols));
air_quart = nan(length(air_cols),2);
for a=1:numel(air_cols)
    air_means(a) = nanmean(air_cols{a});
    air_std(a) = nanstd(air_cols{a});
    air_meds(a) = nanmedian(air_cols{a});
    air_quart(a,:) = quantile(air_cols{a},[0.25,0.75]);
end

figure; subplot(2,1,1);
lmean = line(bin_edges(:,1)',air_means,'color','b','marker','x','linestyle','--');
if uncert; scatter_errorbars(bin_edges(:,1)',air_means,air_std,'color','b'); end
lmed = line(bin_edges(:,1)',air_meds,'color','r','marker','o','linestyle',':');
if uncert; scatter_errorbars(bin_edges(:,1)',air_means,air_quart(:,1)',air_quart(:,2)','color','r'); end
legend([lmean;lmed],{'Mean col density','Median col density'},'location','northwest');
set(gca,'fontsize',16);
xlabel('Bin start - hour of day local std. time');
ylabel('Column density');
ylim([-5e15,2e16]);
title(sprintf('%s col density trend:\n1 hr bins',campaign_name));

% Now do the same but with two-hour rolling windows
bin_edges = [(0:22)',(2:24)'];
air_cols = cell(1,size(bin_edges,1));

for a=1:numel(timevec)
    xx = find(timevec(a) >= bin_edges(:,1) & timevec(a) < bin_edges(:,2));
    for b=1:numel(xx)
        air_cols{xx(b)} = [air_cols{xx(b)}; airno2_iall(a)];
    end
end

% Take the mean and median of each bin and plot

air_means = nan(size(air_cols));
air_std = nan(size(air_cols));
air_meds = nan(size(air_cols));
air_quart = nan(length(air_cols),2);
for a=1:numel(air_cols)
    air_means(a) = nanmean(air_cols{a});
    air_std(a) = nanstd(air_cols{a});
    air_meds(a) = nanmedian(air_cols{a});
    air_quart(a,:) = quantile(air_cols{a},[0.25,0.75]);
end

subplot(2,1,2);
lmean = line(bin_edges(:,1)'+1,air_means,'color','b','marker','x','linestyle','--');
if uncert; scatter_errorbars(bin_edges(:,1)'+1,air_means,air_std,'color','b'); end
lmed = line(bin_edges(:,1)'+1,air_meds,'color','r','marker','o','linestyle',':');
if uncert; scatter_errorbars(bin_edges(:,1)'+1,air_means,air_quart(:,1)',air_quart(:,2)','color','r'); end
legend([lmean;lmed],{'Mean col density','Median col density'},'location','northwest');
set(gca,'fontsize',16);
xlabel('Bin center - hour of day local std. time');
ylabel('Column density');
ylim([-5e15,2e16]);
title(sprintf('%s col density trend:\n2 hr rolling bins',campaign_name));