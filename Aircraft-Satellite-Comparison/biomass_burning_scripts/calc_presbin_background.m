function [species_backgrounds, background_std_dev, bin_pressures, bin_edges] = calc_presbin_background(Merge, species_name, campaign_name)
% For the given species, calculate the background value as the mean of the
% bottom "bottom_percent" of the data. The species field must be present in
% Merge.Data. This will bin by pressure and return vectors
%
%   Josh Laughner <joshlaugh5@gmail.com> 12 June 2015

E=JLLErrors;

Names = merge_field_names(campaign_name);

if ~isfield(Names, species_name)
    E.badinput('%s is not a defined species for %s', species_name, campaign_name)
end

bottom_percent = 0.1;

species = remove_merge_fills(Merge, Names.(species_name));
pres = remove_merge_fills(Merge, Names.pressure);
nans = isnan(pres) | isnan(species);
species(nans) = [];
pres(nans) = [];
% We want to bin by the omi SP bins, but bin_omisp_pressure will bin and
% average/median, so instead we give it some fake data and just get the bin
% edges. Then we use those bin edges to bin (but not average) the data.
[~, bin_pressures, ~, bin_edges] = bin_omisp_pressure(1,1);
% Binning expects the matrix of bin edges to have the smaller value on the
% left, but this is backward with pressure, so we have to adjust.
bins = bin_data(pres, species, fliplr(bin_edges));

species_backgrounds = nan(size(bins));
background_std_dev = nan(size(bins));
for a=1:numel(bins)
    if isempty(bins{a})
        continue
    end
    species = sort(bins{a});
    fifth_percentile_ind = round(numel(species) * bottom_percent);
    
    species_subset = species(1:fifth_percentile_ind);
    species_backgrounds(a) = mean(species_subset);
    
    background_std_dev(a) = std(species_subset);
end

end

