% Script for accountability of data produced for DISCOVER -MD, -CA, -TX,
% and -CO plus SEAC4RS aerosol comparisons

%% SEAC4RS data
% categorize_aerosol_profiles run --> produced SEAC4RS_AerCat
% Run_all_aer_categories run --> produced Comparison

SEAC4RS.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
SEAC4RS.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

SEAC4RS.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
SEAC4RS.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
SEAC4RS.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
SEAC4RS.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
SEAC4RS.coinc.reldiff = (SEAC4RS.coinc.behrno2 - SEAC4RS.coinc.airno2)./SEAC4RS.coinc.airno2 * 100;

SEAC4RS.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
SEAC4RS.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
SEAC4RS.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
SEAC4RS.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
SEAC4RS.aerabove.reldiff = (SEAC4RS.aerabove.behrno2 - SEAC4RS.aerabove.airno2)./SEAC4RS.aerabove.airno2 * 100;

SEAC4RS.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
SEAC4RS.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
SEAC4RS.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
SEAC4RS.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
SEAC4RS.no2above.reldiff = (SEAC4RS.no2above.behrno2 - SEAC4RS.no2above.airno2)./SEAC4RS.no2above.airno2 * 100;

ComparisonSEAC4RS = Comparison;

%% SEAC4RS regressions

airval_c = SEAC4RS.coinc.reldiff;
airval_a = SEAC4RS.aerabove.reldiff;
airval_n = SEAC4RS.no2above.reldiff;
s_text = 'unscaled';


figure;

% Look for fill values and remove
xx = SEAC4RS.coinc.aod < -1;

xl = [min(min(SEAC4RS.coinc.aod(~xx)),0), max(SEAC4RS.coinc.aod(~xx))];

scatter(SEAC4RS.coinc.aod(~xx),airval_c(~xx));
set(gca,'xlim',xl)
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(SEAC4RS.coinc.aod(~xx),airval_c(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: Coincident layers\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = SEAC4RS.aerabove.aod < -1;

xl = [min(min(SEAC4RS.aerabove.aod(~xx)),0), max(SEAC4RS.aerabove.aod(~xx))];

scatter(SEAC4RS.aerabove.aod(~xx),airval_a(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(SEAC4RS.aerabove.aod(~xx),airval_a(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: Aerosol layer above\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = SEAC4RS.no2above.aod < -1;

xl = [min(min(SEAC4RS.no2above.aod(~xx)),0), max(SEAC4RS.no2above.aod(~xx))];

scatter(SEAC4RS.no2above.aod(~xx),airval_n(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(SEAC4RS.no2above.aod(~xx),airval_n(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: NO_2 layer above\n%s aircraft columns',s_text));

%% DISCOVER-MD
% categorize_aerosol_profiles run --> produced DISCOVERMD_AerCat
% Run_all_aer_categories run --> produced Comparison

db_tmp = cat(1,Comparison(:).db_iall);
stimes_all = cell2mat(cat(2,db_tmp(:).start_times))/3600;

DISCOVER_MD.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
if ~isempty(DISCOVER_MD.all.airno2)
    DISCOVER_MD.all.new_airno2 = scale_profiles_by_time(stimes_all, DISCOVER_MD.all.airno2, 'discovermd');
end
DISCOVER_MD.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

stimes_coinc = cell2mat(cat(2,db_tmp(1:2).start_times))/3600;

DISCOVER_MD.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
DISCOVER_MD.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
DISCOVER_MD.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
DISCOVER_MD.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
DISCOVER_MD.coinc.reldiff = (DISCOVER_MD.coinc.behrno2 - DISCOVER_MD.coinc.airno2)./DISCOVER_MD.coinc.airno2 * 100;
if ~isempty(DISCOVER_MD.coinc.airno2)
    DISCOVER_MD.coinc.new_airno2 = scale_profiles_by_time(stimes_coinc, DISCOVER_MD.coinc.airno2, 'discovermd');
    DISCOVER_MD.coinc.new_reldiff = (DISCOVER_MD.coinc.behrno2 - DISCOVER_MD.coinc.new_airno2)./DISCOVER_MD.coinc.new_airno2 * 100;
end

stimes_aer = cell2mat(cat(2,db_tmp(3:4).start_times))/3600;

DISCOVER_MD.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
DISCOVER_MD.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
DISCOVER_MD.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
DISCOVER_MD.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
DISCOVER_MD.aerabove.reldiff = (DISCOVER_MD.aerabove.behrno2 - DISCOVER_MD.aerabove.airno2)./DISCOVER_MD.aerabove.airno2 * 100;
if ~isempty(DISCOVER_MD.aerabove.airno2)
    DISCOVER_MD.aerabove.new_airno2 = scale_profiles_by_time(stimes_aer, DISCOVER_MD.aerabove.airno2, 'discovermd');
    DISCOVER_MD.aerabove.new_reldiff = (DISCOVER_MD.aerabove.behrno2 - DISCOVER_MD.aerabove.new_airno2)./DISCOVER_MD.aerabove.new_airno2 * 100;
end
stimes_no2 = cell2mat(cat(2,db_tmp(5:6).start_times))/3600;

DISCOVER_MD.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
DISCOVER_MD.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
DISCOVER_MD.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
DISCOVER_MD.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
DISCOVER_MD.no2above.reldiff = (DISCOVER_MD.no2above.behrno2 - DISCOVER_MD.no2above.airno2)./DISCOVER_MD.no2above.airno2 * 100;
if ~isempty(DISCOVER_MD.no2above.airno2)
    DISCOVER_MD.no2above.new_airno2 = scale_profiles_by_time(stimes_no2, DISCOVER_MD.no2above.airno2, 'discovermd');
    DISCOVER_MD.no2above.new_reldiff = (DISCOVER_MD.no2above.behrno2 - DISCOVER_MD.no2above.new_airno2)./DISCOVER_MD.no2above.new_airno2 * 100;
end

ComparisonDISCOVER_MD = Comparison;

%% DISCOVER-MD regressions
scaled = false;

if scaled 
    airval_c = DISCOVER_MD.coinc.new_reldiff;
    airval_a = DISCOVER_MD.aerabove.new_reldiff;
    airval_n = DISCOVER_MD.no2above.new_reldiff;
    s_text = 'scaled';
else
    airval_c = DISCOVER_MD.coinc.reldiff;
    airval_a = DISCOVER_MD.aerabove.reldiff;
    airval_n = DISCOVER_MD.no2above.reldiff;
    s_text = 'unscaled';
end

figure;

scatter(DISCOVER_MD.coinc.aod,airval_c);
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_MD.coinc.aod,airval_c,'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-MD: Coincident layers\n%s aircraft columns',s_text));

figure;

scatter(DISCOVER_MD.aerabove.aod,airval_a);
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_MD.aerabove.aod,airval_a,'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-MD: Aerosol layer above\n%s aircraft columns',s_text));

figure;

% There is an outlier here, find and remove it
sd = nanstd(airval_n);
xx = airval_n > (nanmean(airval_n) + 3*sd);

scatter(DISCOVER_MD.no2above.aod(~xx),airval_n(~xx));
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_MD.no2above.aod(~xx),airval_n(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-MD: NO_2 layer above\n%s aircraft columns',s_text));

%% DISCOVER-CA
% categorize_aerosol_profiles run --> produced DISCOVERCD_AerCat
% Run_all_aer_categories run --> produced Comparison

db_tmp = cat(1,Comparison(:).db_iall);
stimes_all = cell2mat(cat(2,db_tmp(:).start_times))/3600;

DISCOVER_CA.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
DISCOVER_CA.all.new_airno2 = scale_profiles_by_time(stimes_all, DISCOVER_CA.all.airno2, 'discoverca');
DISCOVER_CA.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

stimes_coinc = cell2mat(cat(2,db_tmp(1:2).start_times))/3600;

DISCOVER_CA.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
DISCOVER_CA.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
DISCOVER_CA.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
DISCOVER_CA.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
DISCOVER_CA.coinc.reldiff = (DISCOVER_CA.coinc.behrno2 - DISCOVER_CA.coinc.airno2)./DISCOVER_CA.coinc.airno2 * 100;
if ~isempty(DISCOVER_CA.coinc.airno2)
    DISCOVER_CA.coinc.new_airno2 = scale_profiles_by_time(stimes_coinc, DISCOVER_CA.coinc.airno2, 'discoverca');
    DISCOVER_CA.coinc.new_reldiff = (DISCOVER_CA.coinc.behrno2 - DISCOVER_CA.coinc.new_airno2)./DISCOVER_CA.coinc.new_airno2 * 100;
end

stimes_aer = cell2mat(cat(2,db_tmp(3:4).start_times))/3600;

DISCOVER_CA.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
DISCOVER_CA.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
DISCOVER_CA.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
DISCOVER_CA.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
DISCOVER_CA.aerabove.reldiff = (DISCOVER_CA.aerabove.behrno2 - DISCOVER_CA.aerabove.airno2)./DISCOVER_CA.aerabove.airno2 * 100;
if ~isempty(DISCOVER_CA.aerabove.airno2)
    DISCOVER_CA.aerabove.new_airno2 = scale_profiles_by_time(stimes_aer, DISCOVER_CA.aerabove.airno2, 'discoverca');
    DISCOVER_CA.aerabove.new_reldiff = (DISCOVER_CA.aerabove.behrno2 - DISCOVER_CA.aerabove.new_airno2)./DISCOVER_CA.aerabove.new_airno2 * 100;
end

stimes_no2 = cell2mat(cat(2,db_tmp(5:6).start_times))/3600;

DISCOVER_CA.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
DISCOVER_CA.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
DISCOVER_CA.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
DISCOVER_CA.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
DISCOVER_CA.no2above.reldiff = (DISCOVER_CA.no2above.behrno2 - DISCOVER_CA.no2above.airno2)./DISCOVER_CA.no2above.airno2 * 100;
if ~isempty(DISCOVER_CA.no2above.airno2)
    DISCOVER_CA.no2above.new_airno2 = scale_profiles_by_time(stimes_no2, DISCOVER_CA.no2above.airno2, 'discoverca');
    DISCOVER_CA.no2above.new_reldiff = (DISCOVER_CA.no2above.behrno2 - DISCOVER_CA.no2above.new_airno2)./DISCOVER_CA.no2above.new_airno2 * 100;
end

ComparisonDISCOVER_CA = Comparison;

%% DISCOVER-CA regressions
scaled = false;

if scaled 
    airval_c = DISCOVER_CA.coinc.new_reldiff;
    airval_a = DISCOVER_CA.aerabove.new_reldiff;
    airval_n = DISCOVER_CA.no2above.new_reldiff;
    s_text = 'scaled';
else
    airval_c = DISCOVER_CA.coinc.reldiff;
    airval_a = DISCOVER_CA.aerabove.reldiff;
    airval_n = DISCOVER_CA.no2above.reldiff;
    s_text = 'unscaled';
end

figure;

% Look for fill values and remove
xx = DISCOVER_CA.coinc.aod < -1;

scatter(DISCOVER_CA.coinc.aod(~xx),airval_c(~xx));
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CA.coinc.aod(~xx),airval_c(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CA: Coincident layers\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_CA.aerabove.aod < -1;

scatter(DISCOVER_CA.aerabove.aod(~xx),airval_a(~xx));
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CA.aerabove.aod(~xx),airval_a(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CA: Aerosol layer above\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_CA.no2above.aod < -1;

scatter(DISCOVER_CA.no2above.aod(~xx),airval_n(~xx));
set(gca,'xlimmode','manual');
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CA.no2above.aod(~xx),airval_n(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CA: NO_2 layer above\n%s aircraft columns',s_text));

%% DISCOVER-TX
% categorize_aerosol_profiles run --> produced DISCOVERCD_AerCat
% Run_all_aer_categories run --> produced Comparison

DISCOVER_TX.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
DISCOVER_TX.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

DISCOVER_TX.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
DISCOVER_TX.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
DISCOVER_TX.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
DISCOVER_TX.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
DISCOVER_TX.coinc.reldiff = (DISCOVER_TX.coinc.behrno2 - DISCOVER_TX.coinc.airno2)./DISCOVER_TX.coinc.airno2 * 100;

DISCOVER_TX.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
DISCOVER_TX.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
DISCOVER_TX.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
DISCOVER_TX.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
DISCOVER_TX.aerabove.reldiff = (DISCOVER_TX.aerabove.behrno2 - DISCOVER_TX.aerabove.airno2)./DISCOVER_TX.aerabove.airno2 * 100;

DISCOVER_TX.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
DISCOVER_TX.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
DISCOVER_TX.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
DISCOVER_TX.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
DISCOVER_TX.no2above.reldiff = (DISCOVER_TX.no2above.behrno2 - DISCOVER_TX.no2above.airno2)./DISCOVER_TX.no2above.airno2 * 100;

ComparisonDISCOVER_TX = Comparison;

%% DISCOVER-TX regressions

airval_c = DISCOVER_TX.coinc.reldiff;
airval_a = DISCOVER_TX.aerabove.reldiff;
airval_n = DISCOVER_TX.no2above.reldiff;
s_text = 'unscaled';


figure;

% Look for fill values and remove
xx = DISCOVER_TX.coinc.aod < -1;

xl = [min(min(DISCOVER_TX.coinc.aod(~xx)),0), max(DISCOVER_TX.coinc.aod(~xx))];

scatter(DISCOVER_TX.coinc.aod(~xx),airval_c(~xx));
set(gca,'xlim',xl)
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_TX.coinc.aod(~xx),airval_c(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: Coincident layers\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_TX.aerabove.aod < -1;

xl = [min(min(DISCOVER_TX.aerabove.aod(~xx)),0), max(DISCOVER_TX.aerabove.aod(~xx))];

scatter(DISCOVER_TX.aerabove.aod(~xx),airval_a(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_TX.aerabove.aod(~xx),airval_a(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: Aerosol layer above\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_TX.no2above.aod < -1;

xl = [min(min(DISCOVER_TX.no2above.aod(~xx)),0), max(DISCOVER_TX.no2above.aod(~xx))];

scatter(DISCOVER_TX.no2above.aod(~xx),airval_n(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_TX.no2above.aod(~xx),airval_n(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-TX: NO_2 layer above\n%s aircraft columns',s_text));

%% DISCOVER-CO
% categorize_aerosol_profiles run --> produced DISCOVERCO_AerCat
% Run_all_aer_categories run --> produced Comparison

DISCOVER_CO.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
DISCOVER_CO.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

DISCOVER_CO.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
DISCOVER_CO.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
DISCOVER_CO.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
DISCOVER_CO.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
DISCOVER_CO.coinc.reldiff = (DISCOVER_CO.coinc.behrno2 - DISCOVER_CO.coinc.airno2)./DISCOVER_CO.coinc.airno2 * 100;

DISCOVER_CO.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
DISCOVER_CO.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
DISCOVER_CO.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
DISCOVER_CO.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
DISCOVER_CO.aerabove.reldiff = (DISCOVER_CO.aerabove.behrno2 - DISCOVER_CO.aerabove.airno2)./DISCOVER_CO.aerabove.airno2 * 100;

DISCOVER_CO.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
DISCOVER_CO.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
DISCOVER_CO.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
DISCOVER_CO.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
DISCOVER_CO.no2above.reldiff = (DISCOVER_CO.no2above.behrno2 - DISCOVER_CO.no2above.airno2)./DISCOVER_CO.no2above.airno2 * 100;

ComparisonDISCOVER_CO = Comparison;

%% DISCOVER-CO regressions

airval_c = DISCOVER_CO.coinc.reldiff;
airval_a = DISCOVER_CO.aerabove.reldiff;
airval_n = DISCOVER_CO.no2above.reldiff;
s_text = 'unscaled';


figure;

% Look for fill values and remove
xx = DISCOVER_CO.coinc.aod < -1;

xl = [min(min(DISCOVER_CO.coinc.aod(~xx)),0), max(DISCOVER_CO.coinc.aod(~xx))];

scatter(DISCOVER_CO.coinc.aod(~xx),airval_c(~xx));
set(gca,'xlim',xl)
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CO.coinc.aod(~xx),airval_c(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CO: Coincident layers\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_CO.aerabove.aod < -1;

xl = [min(min(DISCOVER_CO.aerabove.aod(~xx)),0), max(DISCOVER_CO.aerabove.aod(~xx))];

scatter(DISCOVER_CO.aerabove.aod(~xx),airval_a(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CO.aerabove.aod(~xx),airval_a(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CO: Aerosol layer above\n%s aircraft columns',s_text));

figure;

% Look for fill values and remove
xx = DISCOVER_CO.no2above.aod < -1;

xl = [min(min(DISCOVER_CO.no2above.aod(~xx)),0), max(DISCOVER_CO.no2above.aod(~xx))];

scatter(DISCOVER_CO.no2above.aod(~xx),airval_n(~xx));
set(gca,'xlim',xl);
oylim = max(abs(get(gca,'ylim')));
set(gca,'ylim',[-oylim,oylim]);
set(gca,'fontsize',16)
[xline,yline,str] = calc_fit_line(DISCOVER_CO.no2above.aod(~xx),airval_n(~xx),'regression','rma');
l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
legend(l,{str});
xlabel('AOD');
ylabel('Percent difference in column');
title(sprintf('DISCOVER-CO: NO_2 layer above\n%s aircraft columns',s_text));
