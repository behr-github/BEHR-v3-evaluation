function [ ComparisonSorted ] = make_resid_vs_aod_plots( Comparison, reldiff_avg, campaign_name, colorby )
%make_resid_vs_aod_plots Function to ease the creation of residual vs AOD plots
%   I've been making lots of plots of the difference in satellite and
%   aircraft VCDs vs. AOD.  This function will take a Comparison structure
%   made by Run_all_aer_categories and parse it into a form where all the
%   relevant data is easier to get a hold of.
%
%   The second argument is optional and defaults to false. When true, it
%   uses the average of the BEHR and aircraft NO2 columns as the
%   denominator in the relative difference.  If false, it uses the aircraft
%   column.
%
%   If you include a third argument that is a string describing what
%   campaign this is for (does not need to be recognizable by
%   merge_field_names, it's only to go in plot titles) this will also go
%   ahead and make those plots.
%
%   A fourth argument (also optional) allows you to color the scatter plots
%   by various quantities. Currently implemented are: 'airno2', 'behrno2',
%   and 'ssa'. Make this an empty string if you need to enter a fifth
%   argument but want nothing to color points by.


E = JLLErrors;
% Check input
narginchk(1,4);
if ~isstruct(Comparison)
    E.badinput('The Comparison input must be a structure');
end
if nargin > 1 && ~islogical(reldiff_avg) && (~isnumeric(reldiff_avg) || ~isscalar(reldiff_avg))
    E.badinput('The (optional) input reldiff_avg must be a logical value or a scalar numeric value')
elseif nargin < 2
    reldiff_avg = false;
end
if nargin > 2 && ~ischar(campaign_name)
    E.badinput('campaign_name must be a string, if you choose to input it')
end
    

ComparisonSorted.all.airno2 = cat(1,Comparison(1:6).airno2_iall);
ComparisonSorted.all.behrno2 = cat(1,Comparison(1:6).behrno2_iall);

ComparisonSorted.coinc.airno2 = cat(1,Comparison(1:2).airno2_iall);
ComparisonSorted.coinc.behrno2 = cat(1,Comparison(1:2).behrno2_iall);
ComparisonSorted.coinc.aod = cell2mat(cat(2,Comparison(1).db_iall.aer_int_out,Comparison(2).db_iall.aer_int_out))';
ComparisonSorted.coinc.ssa = cell2mat(cat(2,Comparison(1).db_iall.aer_median_ssa,Comparison(2).db_iall.aer_median_ssa))';
if reldiff_avg
    ComparisonSorted.coinc.reldiff = (ComparisonSorted.coinc.behrno2 - ComparisonSorted.coinc.airno2)./mean([ComparisonSorted.coinc.airno2, ComparisonSorted.coinc.behrno2],2) * 100;
else
    ComparisonSorted.coinc.reldiff = (ComparisonSorted.coinc.behrno2 - ComparisonSorted.coinc.airno2)./ComparisonSorted.coinc.airno2 * 100;
end

ComparisonSorted.aerabove.airno2 = cat(1,Comparison(3:4).airno2_iall);
ComparisonSorted.aerabove.behrno2 = cat(1,Comparison(3:4).behrno2_iall);
ComparisonSorted.aerabove.aod = cell2mat(cat(2,Comparison(3).db_iall.aer_int_out,Comparison(4).db_iall.aer_int_out))';
ComparisonSorted.aerabove.ssa = cell2mat(cat(2,Comparison(3).db_iall.aer_median_ssa,Comparison(4).db_iall.aer_median_ssa))';
if reldiff_avg
    ComparisonSorted.aerabove.reldiff = (ComparisonSorted.aerabove.behrno2 - ComparisonSorted.aerabove.airno2)./mean([ComparisonSorted.aerabove.airno2, ComparisonSorted.aerabove.behrno2],2) * 100;
else
    ComparisonSorted.aerabove.reldiff = (ComparisonSorted.aerabove.behrno2 - ComparisonSorted.aerabove.airno2)./ComparisonSorted.aerabove.airno2 * 100;
end

ComparisonSorted.no2above.airno2 = cat(1,Comparison(5:6).airno2_iall);
ComparisonSorted.no2above.behrno2 = cat(1,Comparison(5:6).behrno2_iall);
ComparisonSorted.no2above.aod = cell2mat(cat(2,Comparison(5).db_iall.aer_int_out,Comparison(6).db_iall.aer_int_out))';
ComparisonSorted.no2above.ssa = cell2mat(cat(2,Comparison(5).db_iall.aer_median_ssa,Comparison(6).db_iall.aer_median_ssa))';
if reldiff_avg
    ComparisonSorted.no2above.reldiff = (ComparisonSorted.no2above.behrno2 - ComparisonSorted.no2above.airno2)./mean([ComparisonSorted.no2above.airno2, ComparisonSorted.no2above.behrno2],2) * 100;
else
    ComparisonSorted.no2above.reldiff = (ComparisonSorted.no2above.behrno2 - ComparisonSorted.no2above.airno2)./ComparisonSorted.no2above.airno2 * 100;
end
    

if nargin > 2
    if nargin < 4
        colorby = '';
    end

    make_plots(ComparisonSorted, campaign_name, colorby);
end

end

function make_plots(ComparisonSorted, campaign_name, colorby)
E = JLLErrors;
goodfields = unique(cat(1,fieldnames(ComparisonSorted.coinc), fieldnames(ComparisonSorted.aerabove), fieldnames(ComparisonSorted.no2above)));
E.addCustomError('colorby','bad_color_param',sprintf('The parameter (%%s) given for colorby is not available to use; available choices are %s', strjoin(goodfields, ', ')));

airval_c = ComparisonSorted.coinc.reldiff;
airval_a = ComparisonSorted.aerabove.reldiff;
airval_n = ComparisonSorted.no2above.reldiff;


aodval_c = ComparisonSorted.coinc.aod;
aodval_a = ComparisonSorted.aerabove.aod;
aodval_n = ComparisonSorted.no2above.aod;
x_text = 'AOD';


% Look for fill values and remove
xx = aodval_c < -1;

if sum(~xx) > 0
    figure;
    xl = [min(min(aodval_c(~xx)),0), max(aodval_c(~xx))];
    if isnan(xl(2)); xl(2) = 1; end
    
    if isempty(colorby)
        scatter(aodval_c(~xx),airval_c(~xx));
    else
        try
            scatter(aodval_c(~xx),airval_c(~xx), 36, ComparisonSorted.coinc.(colorby)(~xx));
            colormap('jet')
            cb=colorbar;
            cb.Label.String = colorby;
        catch err
            if strcmp(err.identifier, 'MATLAB:nonExistentField')
                E.callCustomError('colorby',colorby)
            else
                rethrow(err)
            end
        end
    end
    set(gca,'xlim',xl)
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(aodval_c(~xx),airval_c(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel(x_text);
    ylabel('Percent difference in column');
    title(sprintf('%s: Coincident layers',upper(campaign_name)));
end


% Look for fill values and remove
xx = aodval_a < -1;

if sum(~xx) > 0
    figure;
    xl = [min(min(aodval_a(~xx)),0), max(aodval_a(~xx))];
    if isnan(xl(2)); xl(2) = 1; end

    if isempty(colorby)
        scatter(aodval_a(~xx),airval_a(~xx));
    else
        try
            scatter(aodval_a(~xx),airval_a(~xx), 36, ComparisonSorted.aerabove.(colorby)(~xx));
            colormap('jet')
            cb=colorbar;
            cb.Label.String = colorby;
        catch err
            if strcmp(err.identifier, 'MATLAB:nonExistentField')
                E.callCustomError('colorby',colorby)
            else
                rethrow(err)
            end
        end
    end
    set(gca,'xlim',xl);
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(aodval_a(~xx),airval_a(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel(x_text);
    ylabel('Percent difference in column');
    title(sprintf('%s: Aerosol layer above',upper(campaign_name)));
end

% Look for fill values and remove
xx = aodval_n < -1;

if sum(~xx)
    figure;
    xl = [min(min(aodval_n(~xx)),0), max(aodval_n(~xx))];
    if isnan(xl(2)); xl(2) = 1; end

    if isempty(colorby)
        scatter(aodval_n(~xx),airval_n(~xx));
    else
        try
            scatter(aodval_n(~xx),airval_n(~xx), 36, ComparisonSorted.no2above.(colorby)(~xx));
            colormap('jet')
            cb=colorbar;
            cb.Label.String = colorby;
        catch err
            if strcmp(err.identifier, 'MATLAB:nonExistentField')
                E.callCustomError('colorby',colorby)
            else
                rethrow(err)
            end
        end
    end
    set(gca,'xlim',xl);
    oylim = max(abs(get(gca,'ylim')));
    set(gca,'ylim',[-oylim,oylim]);
    set(gca,'fontsize',16)
    [xline,yline,str] = calc_fit_line(aodval_n(~xx),airval_n(~xx),'regression','rma');
    l = line(xline,yline,'color','k','linestyle','--','linewidth',2);
    legend(l,{str});
    xlabel(x_text);
    ylabel('Percent difference in column');
    title(sprintf('%s: NO_2 layer above',upper(campaign_name)));
end
end
