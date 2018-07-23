function [  ] = ca_ground_site_intercomp(  )
%ca_ground_site_intercomp Plot NO2 values for different stations at same site
%   Detailed explanation goes here

% First do Bakersfield

p='/Volumes/share2/USERS/LaughnerJ/CampaignMergeMats/DISCOVER-AQ_CA/Ground-AllDays/VariousTimePeriods/';
Orig = load(fullfile(p,'California_Ground1_Bakersfield.mat'),'Merge');
Trevino = load(fullfile(p,'California_Groundx1_Bakersfield_Trevino_2013_01_06.mat'),'Merge');

o_doy = remove_merge_fills(Orig.Merge,'start_DOYUTC');
t_doy = floor(remove_merge_fills(Trevino.Merge,'Fractional_Julian_Day'));

o_udoy = unique(o_doy);
t_udoy = unique(t_doy);

o_utcmid = remove_merge_fills(Orig.Merge,'mid_secUTC');
t_utcmid = remove_merge_fills(Trevino.Merge,'UTC_mid');

o_no2 = remove_merge_fills(Orig.Merge,'photo_NO2_ppbv');
t_no2 = remove_merge_fills(Trevino.Merge,'NO2');

matched_no2_bakersfield = matchUTC(o_udoy, t_udoy, o_doy, t_doy, o_utcmid, t_utcmid, o_no2, t_no2);

figure; scatter(matched_no2_bakersfield.orig, matched_no2_bakersfield.trev,36,'k');
axes_equal;
set(gca,'fontsize',16)
xlabel('Russell Long - EPA: Teledyne NO2');
ylabel('Nathan Trevino - SJV UAPCD');
title('Bakersfield NO2');

p='/Volumes/share2/USERS/LaughnerJ/CampaignMergeMats/DISCOVER-AQ_CA/Ground-AllDays/VariousTimePeriods/';
Orig = load(fullfile(p,'California_Ground1_Bakersfield.mat'),'Merge');
Trevino1 = load(fullfile(p,'California_Groundx6_Fresno_Drummond_Trevino_2013_01_06.mat'),'Merge');
Trevino2 = load(fullfile(p,'California_Groundx6_Fresno_Skypark_Trevino_2013_01_06.mat'),'Merge');

o_doy = remove_merge_fills(Orig.Merge,'start_DOYUTC');
t1_doy = floor(remove_merge_fills(Trevino1.Merge,'Fractional_Julian_Day'));
t2_doy = floor(remove_merge_fills(Trevino2.Merge,'Fractional_Julian_Day'));

o_udoy = unique(o_doy);
t1_udoy = unique(t1_doy);
t2_udoy = unique(t2_doy);

o_utcmid = remove_merge_fills(Orig.Merge,'mid_secUTC');
t1_utcmid = remove_merge_fills(Trevino1.Merge,'UTC_mid');
t2_utcmid = remove_merge_fills(Trevino2.Merge,'UTC_mid');

o_no2 = remove_merge_fills(Orig.Merge,'photo_NO2_ppbv');
t1_no2 = remove_merge_fills(Trevino1.Merge,'NO2');
t2_no2 = remove_merge_fills(Trevino2.Merge,'NO2');

matched_no2_fresno_drummond = matchUTC(o_udoy, t1_udoy, o_doy, t1_doy, o_utcmid, t1_utcmid, o_no2, t1_no2);
matched_no2_fresno_skypark = matchUTC(o_udoy, t2_udoy, o_doy, t2_doy, o_utcmid, t2_utcmid, o_no2, t2_no2);

figure; s(1) = scatter(matched_no2_fresno_drummond.orig, matched_no2_fresno_drummond.trev,36,'b'); 
hold on
s(2) = scatter(matched_no2_fresno_skypark.orig, matched_no2_fresno_skypark.trev,36,'r');
axes_equal;
set(gca,'fontsize',16)
legend(s',{'Drummond SJV UAPCD','Skypark SJV UAPCD'})
xlabel('CARB NO2');
ylabel('Nathan Trevino - SJV UAPCD');
title('Fresno NO2');

end

function matched_no2 = matchUTC(o_udoy, t_udoy, o_doy, t_doy, o_utcmid, t_utcmid, o_no2, t_no2)

matched_no2.orig = [];
matched_no2.trev = [];

for a=1:numel(o_udoy)
    if containedin(o_udoy(a),t_udoy)
        xx_o = o_doy == o_udoy(a);
        xx_t = t_doy == t_udoy(a);
        
        this_o_utc = o_utcmid(xx_o);
        this_t_utc = t_utcmid(xx_t);
        this_o_no2 = o_no2(xx_o);
        this_t_no2 = t_no2(xx_t);
        
        for b=1:numel(this_o_utc)
            zz = this_t_utc == this_o_utc(b);
            if sum(zz) == 1
                matched_no2.orig = cat(1,matched_no2.orig,this_o_no2(b));
                matched_no2.trev = cat(1,matched_no2.trev,this_t_no2(zz));
            elseif sum(zz) > 1
                error('intercomp:too_many_points','More than one matching UTC');
            end
        end
    end
end

end

