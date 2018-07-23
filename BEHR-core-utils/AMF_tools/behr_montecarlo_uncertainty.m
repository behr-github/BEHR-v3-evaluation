function [ amf_var, amf_relvar, param_vals ] = behr_montecarlo_uncertainty( AMFS, uncert_param, uncert_range )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

params = {'SZAs', 'VZAs', 'RAAs', 'ALBs', 'SurfPs'};
n_params = numel(params);
if ~ismember(uncert_param, params)
    E.badinput('UNCERT_PARAM must be one of: %s', strjoin(uncert_param, ', '));
end

n_points = 10000;

% Get the range of values for each parameter and set up the random points
% within those ranges to sample. At the same time, extract the relevant
% vectors for the gridded interpolation, flipping them around to be
% strictly increasing if necessary.
param_ranges = nan(n_params, 2);
param_vals = nan(n_params, n_points);
param_vecs = cell(n_params, 1);
amfs = AMFS.AMFs;
for p=1:n_params
    if ndims(AMFS.(params{p})) ~= n_params
        E.badinput('AMFS.%s does not have the same number of dimensions as there are parameters', params{p});
    end
    
    this_range = [min(AMFS.(params{p})(:)), max(AMFS.(params{p})(:))];
    param_ranges(p,:) = this_range;
    param_vals(p, :) = randrange(this_range(1), this_range(2), 1, n_points);
    
    pvec = perm_vec(p, n_params);
    tmp = permute(AMFS.(params{p}), pvec);
    tmp_vec = tmp(:,1);
    if all(diff(tmp_vec) < 0)
        tmp_vec = flipud(tmp_vec);
        amfs = flip(amfs, p);
    elseif ~all(diff(tmp_vec) > 0)
        E.badinput('Dimension %d of AMFS.%s is not strictly increasing or decreasing', p, params{p});
    end
    
    param_vecs{p} = tmp_vec;
end


AInterp = griddedInterpolant(param_vecs, amfs);

i_uncert = strcmp(uncert_param, params);
amf_var = nan(1, n_points);

for a=1:n_points
    % For each sample point, vary the uncertain parameter by +/- the
    % uncertainty range and calculate the AMF value at the original point
    % and the +/- values
    test_vals = repmat(param_vals(:,p),1,3);
    test_vals(i_uncert, 2) = test_vals(i_uncert, 2) - uncert_range;
    test_vals(i_uncert, 3) = test_vals(i_uncert, 3) + uncert_range;
    test_vals(i_uncert, :) = clipmat(test_vals(i_uncert,:), param_ranges(i_uncert,:));
    test_vals = num2cell(test_vals);
    
    test_amfs = nan(1,3);
    for b=1:3
        test_amfs(b) = AInterp(test_vals(:,b));
    end
    
    amf_var(a) = (abs(test_amfs(2) - test_amfs(1)) + abs(test_amfs(3) - test_amfs(1)))/2;
    amf_relvar(a) = amf_var(a) / test_amfs(1);
end

end

function pvec = perm_vec(dim, n)
% Create a vector that, when used with PERMUTE() will permute an N
% dimensional array to have DIM along the first dimension.
pvec = 1:n;
pvec(pvec == dim) = [];
pvec = [dim, pvec];
end