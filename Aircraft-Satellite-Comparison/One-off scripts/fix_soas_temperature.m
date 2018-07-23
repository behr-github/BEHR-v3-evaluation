function [ output_args ] = fix_soas_temperature(  )
%FIX_SOAS_TEMPERATURE Convert the SOAS aircraft ambient temperature to K

save_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignMergeMats/SOAS/P3/1sec';
F = dir(fullfile(save_dir, '*.mat'));

for a=1:numel(F)
    fprintf('Correcting %s\n', F(a).name);
    load(fullfile(save_dir, F(a).name)); % Puts "merge" variable in workspace
    Merge.Data.AmbTemp.Values(Merge.Data.AmbTemp.Values > -1000) = Merge.Data.AmbTemp.Values(Merge.Data.AmbTemp.Values > -1000) + 273.15;
    Merge.Data.AmbTemp.Unit = 'K';
    save(fullfile(save_dir, F(a).name), 'Merge');
end

end

