function [ aitken_mean_mode_number, aitken_med_mode_number ] = sum_aitken_mode( matpath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
E = JLLErrors;
F = dir(fullfile(matpath, '*.mat'));
M = load(fullfile(matpath, F(1).name), 'Merge');

if ~isempty(strfind(matpath, 'DISCOVER-AQ_MD'))
    aitken_fns = aitken_fields_md(M.Merge);
elseif ~isempty(strfind(matpath, 'DISCOVER-AQ_CA'))
    aitken_fns = aitken_fields_ca(M.Merge);
elseif ~isempty(strfind(matpath, 'DISCOVER-AQ_TX'))
    aitken_fns = aitken_fields_tx(M.Merge);
else
    E.notimplemented('I don''t know what fields to look for given the path %s', matpath);
end

particle_counts = [];
n_points = 0;
for a=1:numel(F)
    M = load(fullfile(matpath, F(a).name), 'Merge');
    % All fields BETTER have the same number of points
    this_particle_count = zeros(size(M.Merge.Data.(aitken_fns{1}).Values));
    for b=1:numel(aitken_fns)
        % Add up the bins in the aitken mode
        val = remove_merge_fills(M.Merge, aitken_fns{b});
        this_particle_count = this_particle_count + val;
    end
    particle_counts = veccat(particle_counts, this_particle_count);
end

% By doing the average ourselves (rather than adding the mean of each field
% for each day above, i.e. aitken_mode_num = aitken_mode_num + nanmean(...)
% we avoid giving points on days with shorter flights more individidual
% weight.
%
% The numbers are reported in dN/dlogDp, and at least in CA and TX the bins
% are specified to be dlogDp = 0.05 wide, so to do the integration
% properly, it must be scaled by 0.05 (I think - this is assuming we can
% treat numerically integrating this distributions using a simple center
% point rule basically.)
aitken_mean_mode_number = nanmean(particle_counts)*0.05;
aitken_med_mode_number = nanmedian(particle_counts)*0.05;
end

function fields = aitken_fields_md(Merge)
% Maryland fields just start with a number in the ICARTT file, so they get
% "f_" prepended in the structure field name. They are also all given to 1
% decimal place, and we're looking for fields with centers between 10 and
% 100 nm. The decimal gets removed in the field name, so the numbers in the
% field name should be divided by 10 to get the actual midpoint value.

fns = search_merge_fields(Merge, 'f_');
fields = cell(1,100);
b=1;
for a=1:numel(fns)
    bin_mid = str2double(regexp(fns{a}, '\d+(?=nm)', 'once', 'match'))/10;
    if bin_mid >= 10 && bin_mid <= 100
        fields{b} = fns{a};
        b=b+1;
    end
end
fields(b:end) = [];
end

function fields = aitken_fields_ca(Merge)
% California fields begin with dNdlogDp. They do not include any decimal
% components in the bin midpoints.

fns = search_merge_fields(Merge, 'dNdlogDp_');
fields = cell(1,100);
b=1;
for a=1:numel(fns)
    bin_mid = str2double(regexp(fns{a}, '\d+(?=nm)', 'once', 'match'));
    if bin_mid >= 10 && bin_mid <= 100
        fields{b} = fns{a};
        b = b+1;
    end
end
fields(b:end) = [];
end

function fields = aitken_fields_tx(Merge)
% California fields begin with SMPS_Bin. They are labeled by bin #, the
% first 20 bins span 10-100 nm and the bin indices are zero-based

fns = search_merge_fields(Merge, 'SMPS_Bin');
fields = cell(1,20);
b=1;
for a=1:numel(fns)
    bin_num = str2double(regexp(fns{a}, '\d\d$', 'once', 'match'));
    if bin_num < 20
        fields{b} = fns{a};
        b = b+1;
    end
end
fields(b:end) = [];
end