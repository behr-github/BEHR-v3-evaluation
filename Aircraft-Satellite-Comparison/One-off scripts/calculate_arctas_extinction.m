function [  ] = calculate_arctas_extinction(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

campaign_name = 'arctas-b';

[~,~,merge_dir] = merge_field_names(campaign_name);
save_dir = '/Users/Josh/Documents/MATLAB/tmp';

F = dir(fullfile(merge_dir,'*.mat'));

for a=1:numel(F)
    load(fullfile(merge_dir, F(a).name));
    scattering_blue = remove_merge_fills(Merge, 'Total_Scatter450_nm');
    scattering_green = remove_merge_fills(Merge, 'Total_Scatter550_nm');
    scattering_red = remove_merge_fills(Merge, 'Total_Scatter700_nm');
    
    abs_blue = remove_merge_fills(Merge,'Total_Absorption470_nm');
    abs_green = remove_merge_fills(Merge, 'Total_Absorption532_nm');
    abs_red = remove_merge_fills(Merge, 'Total_Absorption660_nm');
    
    ext_blue = scattering_blue + abs_blue;
    ext_blue(isnan(ext_blue)) = Merge.Data.Total_Absorption470_nm.Fill;
    
    ext_green = scattering_green + abs_green;
    ext_green(isnan(ext_green)) = Merge.Data.Total_Absorption532_nm.Fill;
    
    ext_red = scattering_red + abs_red;
    ext_red(isnan(ext_red)) = Merge.Data.Total_Absorption660_nm.Fill;
    
    Merge.Data.TotalExtinctionBlue.Values = ext_blue;
    Merge.Data.TotalExtinctionBlue.Unit = Merge.Data.Total_Absorption470_nm.Unit;
    Merge.Data.TotalExtinctionBlue.Fill = Merge.Data.Total_Absorption470_nm.Fill;
    
    Merge.Data.TotalExtinctionGreen.Values = ext_green;
    Merge.Data.TotalExtinctionGreen.Unit = Merge.Data.Total_Absorption532_nm.Unit;
    Merge.Data.TotalExtinctionGreen.Fill = Merge.Data.Total_Absorption532_nm.Fill;
    
    Merge.Data.TotalExtinctionRed.Values = ext_red;
    Merge.Data.TotalExtinctionRed.Unit = Merge.Data.Total_Absorption660_nm.Unit;
    Merge.Data.TotalExtinctionRed.Fill = Merge.Data.Total_Absorption660_nm.Fill;
    
    DataTable = cat(2, DataTable, ext_blue', ext_green', ext_red');
    header = cat(2, header, {'TotalExtinctionBlue', 'TotalExtinctionGreen', 'TotalExtinctionRed'});
    
    % Save to a temporary directory so we can check them first
    save(fullfile(save_dir,F(a).name),'DataTable','header','Merge');
    clear DataTable header Merge
end

end

