function [  ] = frac_propene_texaqs(  )
%FRAC_PROPENE_TEXAQS Calculate the fraction of >=3 C alkenes that are propene
%   GEOS-Chem outputs all alkenes with 3 or more carbons as PRPE. WRF-Chem
%   (at least the R2SMH and RACM2 mechanisms) expects alkenes to be divided
%   into internal and terminal alkene. When converting from MOZARTS BIGENE
%   lumped species (which is alkenes with 4 or more carbons), it's just
%   split 50/50. So once we know how much of PRPE to assume is propene,
%   the simplest approach will be to make the WRF terminal alkenes that
%   fraction plus half of what's left.

% Fields beginning with wF are the most quantitative (run using the FID
% PLOT column). Fields beginning with wM somewhat quantitative (run with
% the HP-MSD column). Fields beginning with wV are not quantitative and
% should not be used. If a measurement is given for both wF and wM, wF is
% preferred.

% Define fields
propene = 'wF_propene';
internal_c4alkenes = {'wM_c_2_pentene',...
                      'wM_cyclopentene',...
                      'wM_t_2_pentene',...
                      'wF_t_2_butene',...
                      'wM_2_methyl_2_butene',...
                      'wM_c_2_butene',...
                      };
terminal_c4alkenes = {'wM_1_pentene',...
                      'wF_1_butene',...
                      'wF_1_i_butene',... 
                      'wF_i_butene',... % unclear how this an 1_i_butene differ
                      'wM_3_methyl_1_butene',...
                      };

propene_out = campaign_wide_ops('texaqs',propene,'cat','debug',0);
internal_out = campaign_wide_ops('texaqs',internal_c4alkenes,'cat','debug',0);
terminal_out = campaign_wide_ops('texaqs',terminal_c4alkenes,'cat','debug',0);

propene_sum = sum_fields(propene_out);
internal_sum = sum_fields(internal_out);
terminal_sum = sum_fields(terminal_out);
                  
propene_avg = campaign_wide_ops('texaqs',propene,'dayavg','debug',0);
internal_avg = campaign_wide_ops('texaqs',internal_c4alkenes,'dayavg','debug',0);
terminal_avg = campaign_wide_ops('texaqs',terminal_c4alkenes,'dayavg','debug',0);

propene_avg_sum = sum_fields(propene_avg); % should only be only field, but gets it converted to a vector rather than a structure
internal_avg_sum = sum_fields(internal_avg);
terminal_avg_sum = sum_fields(terminal_avg);

dnums = datenum(propene_avg.dates);
M = [propene_avg_sum, internal_avg_sum, terminal_avg_sum];

figure; 
bar(dnums, M,'stacked'); 
legend('Propene','Internal C4+C5','Terminal C4+C5');
datetick('x');

sums = [nansum2(propene_sum), nansum2(internal_sum), nansum2(terminal_sum)];

figure;
pie(sums,{'Propene','Internal C4+C5','Terminal C4+C5'});

fprintf('Average percent propene = %.2f\n', sums(1) / nansum2(sums) * 100);
fprintf('Average percent internal = %.2f\n', sums(2) / nansum2(sums) * 100);
fprintf('Average percent terminal = %.2f\n', sums(3) / nansum2(sums) * 100);



end

function s = sum_fields(Out)
fns = fieldnames(Out.data);
s = zeros(size(Out.data.(fns{1})));
for f=1:numel(fns)
    data = Out.data.(fns{1});
    data(isnan(data)) = 0;
    s = s + data;
end
end