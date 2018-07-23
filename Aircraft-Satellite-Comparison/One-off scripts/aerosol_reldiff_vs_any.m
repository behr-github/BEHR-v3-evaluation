function [  ] = aerosol_reldiff_vs_any( Comparison, quantity, campaign_name, avg_reldiff, colorby )
%AEROSOL_RELDIFF_VS_ANY Troubleshooting plots for aerosol project
%   I'm trying to figure out why I'm getting such large percent errors in
%   my percent-difference-vs-AOD plots, so this will let me plot the
%   percent different as a function of nearly any quantity in the
%   Comparison structure output in Run_all_aer_categories.
%
%   Josh Laughner <joshlaugh5@gmail.com> 18 Aug 2015

E = JLLErrors;
if ~exist('campaign_name','var')
    campaign_name = '';
end
if ~exist('avg_reldiff','var')
    avg_reldiff = true;
end
if ~exist('colorby','var')
    colorby = '';
end

% First organize everything into coincident, aerosol above, NO2 above.
% Extract values from the db_iall sub structure, averaging any cells with
% more than one value.

CompSorted = struct('Coincident',struct,'AerAbove',struct,'NO2Above',struct);

fns = fieldnames(Comparison);
for a=1:numel(fns)
    if strcmpi(fns{a}, 'dates_iall')
        CompSorted.Coincident.dates_iall = cat(1,datenum(Comparison(1).dates_iall),datenum(Comparison(2).dates_iall));
        CompSorted.AerAbove.dates_iall = cat(1,datenum(Comparison(3).dates_iall),datenum(Comparison(4).dates_iall));
        CompSorted.NO2Above.dates_iall = cat(1,datenum(Comparison(5).dates_iall),datenum(Comparison(6).dates_iall));
    elseif strcmpi(fns{a}, 'db_iall') || strcmpi(fns{a}, 'category')
        continue
        % we'll deal with db_iall next
    else
        CompSorted.Coincident.(fns{a}) = cat(1,Comparison(1:2).(fns{a}));
        CompSorted.AerAbove.(fns{a}) = cat(1,Comparison(3:4).(fns{a}));
        CompSorted.NO2Above.(fns{a}) = cat(1,Comparison(5:6).(fns{a}));
    end
end

if avg_reldiff
    CompSorted.Coincident.PerDiff = (CompSorted.Coincident.behrno2_iall - CompSorted.Coincident.airno2_iall) ./ mean([CompSorted.Coincident.airno2_iall, CompSorted.Coincident.behrno2_iall],2) * 100;
    CompSorted.AerAbove.PerDiff = (CompSorted.AerAbove.behrno2_iall - CompSorted.AerAbove.airno2_iall) ./ mean([CompSorted.AerAbove.airno2_iall, CompSorted.AerAbove.behrno2_iall],2) * 100;
    CompSorted.NO2Above.PerDiff = (CompSorted.NO2Above.behrno2_iall - CompSorted.NO2Above.airno2_iall) ./ mean([CompSorted.NO2Above.airno2_iall, CompSorted.NO2Above.behrno2_iall],2) * 100;
else
    CompSorted.Coincident.PerDiff = (CompSorted.Coincident.behrno2_iall - CompSorted.Coincident.airno2_iall) ./ CompSorted.Coincident.airno2_iall * 100;
    CompSorted.AerAbove.PerDiff = (CompSorted.AerAbove.behrno2_iall - CompSorted.AerAbove.airno2_iall) ./ CompSorted.AerAbove.airno2_iall * 100;
    CompSorted.NO2Above.PerDiff = (CompSorted.NO2Above.behrno2_iall - CompSorted.NO2Above.airno2_iall) ./ CompSorted.NO2Above.airno2_iall * 100;
end
    
fns = fieldnames(Comparison(1).db_iall);
for a=1:numel(fns)
    if ~iscell(Comparison(1).db_iall.(fns{a}))
        continue
        % all entries in db_iall should be a cell array except for the run
        % structure which just contains information about the parameters
        % used in the run.
    elseif ~isempty(regexpi(fns{a},'l..corn', 'once'))
        continue
        % skip the lat/lon corner fields
    end
    CompSorted.Coincident.(fns{a}) = nan(size(CompSorted.Coincident.dates_iall));
    db = cat(1, Comparison(1).db_iall.(fns{a})', Comparison(2).db_iall.(fns{a})');
    if numel(db) ~= numel(CompSorted.Coincident.dates_iall)
        E.badvar('db',sprintf('Different numbers of profiles in the debugging structure for field %s',fns{a}));
    end
    for b=1:numel(CompSorted.Coincident.(fns{a}))
        CompSorted.Coincident.(fns{a})(b) = nanmean(db{b});
    end
    
    CompSorted.AerAbove.(fns{a}) = nan(size(CompSorted.AerAbove.dates_iall));
    db = cat(1, Comparison(3).db_iall.(fns{a})', Comparison(4).db_iall.(fns{a})');
    if numel(db) ~= numel(CompSorted.AerAbove.dates_iall)
        E.badvar('db',sprintf('Different numbers of profiles in the debugging structure for field %s',fns{a}));
    end
    for b=1:numel(CompSorted.AerAbove.(fns{a}))
        CompSorted.AerAbove.(fns{a})(b) = nanmean(db{b});
    end
    
    CompSorted.NO2Above.(fns{a}) = nan(size(CompSorted.NO2Above.dates_iall));
    db = cat(1, Comparison(5).db_iall.(fns{a})', Comparison(6).db_iall.(fns{a})');
    if numel(db) ~= numel(CompSorted.NO2Above.dates_iall)
        E.badvar('db',sprintf('Different numbers of profiles in the debugging structure for field %s',fns{a}));
    end
    for b=1:numel(CompSorted.NO2Above.(fns{a}))
        CompSorted.NO2Above.(fns{a})(b) = nanmean(db{b});
    end
end

allfns = fieldnames(CompSorted.Coincident);
if ~ismember(quantity, allfns)
    E.badinput('quantity must be one of the following fields: %s', strjoin(allfns, ', '));
end

if ~isempty(colorby) && ~ismember(colorby,allfns)
    E.badinput('colorby (if given) must be one of the following fields: %s', strjoin(allfns, ', '));
end

figure; 
if ~isempty(colorby)
    scatter(CompSorted.Coincident.(quantity), CompSorted.Coincident.PerDiff, 36, CompSorted.Coincident.(colorby));
    cb=colorbar;
    cb.Label.String = regexprep(colorby, '_', ' ');
else
    scatter(CompSorted.Coincident.(quantity), CompSorted.Coincident.PerDiff);
end
xlabel(regexprep(quantity,'_',' '));
if avg_reldiff
    ylabel('Percent difference (sat - air)/(mean([air, sat])')
else
    ylabel('Percent difference (sat - air)/air')
end
set(gca,'fontsize',16)
title(sprintf('%s (coincident) - %% diff. vs. %s',upper(campaign_name), regexprep(quantity,'_',' ')));

figure;
if ~isempty(colorby)
    scatter(CompSorted.AerAbove.(quantity), CompSorted.AerAbove.PerDiff, 36, CompSorted.AerAbove.(colorby));
    cb=colorbar;
    cb.Label.String = regexprep(colorby, '_', ' ');
else
    scatter(CompSorted.AerAbove.(quantity), CompSorted.AerAbove.PerDiff);
end
xlabel(regexprep(quantity,'_',' '));
if avg_reldiff
    ylabel('Percent difference (sat - air)/(mean([air, sat])')
else
    ylabel('Percent difference (sat - air)/air')
end
set(gca,'fontsize',16)
title(sprintf('%s (aer. above) - %% diff. vs. %s',upper(campaign_name), regexprep(quantity,'_',' ')));

figure;
if ~isempty(colorby)
    scatter(CompSorted.NO2Above.(quantity), CompSorted.NO2Above.PerDiff, 36, CompSorted.NO2Above.(colorby));
    cb=colorbar;
    cb.Label.String = regexprep(colorby, '_', ' ');
else
    scatter(CompSorted.NO2Above.(quantity), CompSorted.NO2Above.PerDiff);
end
xlabel(regexprep(quantity,'_',' '));
if avg_reldiff
    ylabel('Percent difference (sat - air)/(mean([air, sat])')
else
    ylabel('Percent difference (sat - air)/air')
end
set(gca,'fontsize',16)
title(sprintf('%s (NO_2 Above) - %% diff. vs. %s',upper(campaign_name), regexprep(quantity,'_',' ')));
end

