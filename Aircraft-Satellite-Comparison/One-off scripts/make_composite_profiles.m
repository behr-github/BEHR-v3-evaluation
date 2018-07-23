% makes campaign composite temperature profiles

%campaigns = {'discover-md','discover-ca','discover-tx','seac4rs','intex-b'};
campaigns = {'dc3', 'arctas-carb'};

temps = [];
pres = [];

for a=1:numel(campaigns)
    [Names,~,fdir] = merge_field_names(campaigns{a});
    F = dir(fullfile(fdir,'*.mat'));
    for f=1:numel(F)
        load(fullfile(fdir,F(f).name),'Merge');
        temps = cat(2,temps,remove_merge_fills(Merge,Names.temperature));
        pres = cat(2,pres,remove_merge_fills(Merge,Names.pressure));
    end
    [tempbins,presbins] = bin_omisp_pressure(pres,temps);
    nans = isnan(tempbins);
    tempbins(nans) = interp1(log(presbins(~nans)),tempbins(~nans),log(presbins(nans)),'linear','extrap');
    cname = upper(regexprep(campaigns{a},'\W',''));
    eval(sprintf('%s_TEMP_COMP = tempbins;',cname));
end