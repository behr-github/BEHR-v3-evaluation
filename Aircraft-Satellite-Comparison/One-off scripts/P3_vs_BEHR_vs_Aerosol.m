function P3_vs_BEHR_vs_Aerosol

% This will run through all the permutations of DISCOVER-AQ sites,
% aerosol/NO2 distribution, etc. and generate the correlation coefficient
% and p-value matrices for:
%   - aircraft vs. behr vs. max aerosol extinction
%   - aircraft vs. behr vs. total integrated aerosol extinction
%   - (aircraft-behr) vs. max aerosol extinction
%   - (aircraft-behr) vs. total integrated aerosol extinction
%
% Yes I know the prolific use of eval statements is horrible programming,
% but this isn't intended to be used very frequently.

% Load the workspace that has all the variables prepared. This has a lot of
% variables.
load('/Users/Josh/Documents/MATLAB/NO2 Profiles/Workspaces/MD, CA, and TX data separated by aerosol-NO2 relative height.mat');

% Each variable has four parts to the name:
%   <data_type>_<site[rep]>_<no2/aer_dist>_<h/l>
%   data_type = air, behr, or db (debugging structure w/aerosol data in it)
%   site[rep] = md, ca, tx - the state where the flights occured. 'rep'
%       after indicates the columns were reprocessed using in-situ profiles
%   no2/aer_dist = coinc, aer, or no2 - where the no2 is relative to the
%      aerosol
%   h/l = high or low max extinction.

data_types = {'air','behr','aermax','aertot'}; nt = numel(data_types);
sites = {'md','mdrep','ca','carep','tx','txrep'}; ns = numel(sites);
distribution = {'coinc','aer','no2'}; nd = numel(distribution);

% First extract the aerosol max and total data from the db structure;
% transpose it to be a column
for s=1:ns
    for d=1:nd
        eval(sprintf('aermax_%1$s_%2$s_h = (cell2mat(db_%1$s_%2$s_h.aer_max_out))'';',sites{s},distribution{d}));
        eval(sprintf('aermax_%1$s_%2$s_l = (cell2mat(db_%1$s_%2$s_l.aer_max_out))'';',sites{s},distribution{d}));
        eval(sprintf('aertot_%1$s_%2$s_h = (cell2mat(db_%1$s_%2$s_h.aer_int_out))'';',sites{s},distribution{d}));
        eval(sprintf('aertot_%1$s_%2$s_l = (cell2mat(db_%1$s_%2$s_l.aer_int_out))'';',sites{s},distribution{d}));
    end
end


% Second concatenate the high and low max aerosol data
for t=1:nt;
    for s=1:ns;
        for d=1:nd;
            eval(sprintf('%1$s_%2$s_%3$s = cat(1,%1$s_%2$s_%3$s_l,%1$s_%2$s_%3$s_h);',data_types{t}, sites{s}, distribution{d}));
        end
    end
end

% Third make the correlation and p-value matrices for each data type, site,
% distribution comparing aircraft, BEHR, aerosol max, and total aerosol
% directly. Simultaneously run the correlation between the
% aircraft-satellite difference and the aerosol metrics.

corrStruct = struct;

for s=1:ns;
    for d=1:nd;
        % Correlation doesn't work if there are nans
        %xx_all = false(nt,numel(eval(sprintf('air_%1$s_%2$s',sites{s}, distribution{d}))));
        xx = any(isnan(eval(sprintf('[air_%1$s_%2$s, behr_%1$s_%2$s, aermax_%1$s_%2$s, aertot_%1$s_%2$s]',sites{s}, distribution{d}))),2);
        eval(sprintf('air_%1$s_%2$s(xx) = []',sites{s}, distribution{d}));
        eval(sprintf('behr_%1$s_%2$s(xx) = []',sites{s}, distribution{d}));
        eval(sprintf('aermax_%1$s_%2$s(xx) = []',sites{s}, distribution{d}));
        eval(sprintf('aertot_%1$s_%2$s(xx) = []',sites{s}, distribution{d}));
        [R,P] = corrcoef(eval(sprintf('[air_%1$s_%2$s, behr_%1$s_%2$s, aermax_%1$s_%2$s, aertot_%1$s_%2$s]',sites{s}, distribution{d})));
        corrStruct.all.(sites{s}).(distribution{d}).R2 = R;
        corrStruct.all.(sites{s}).(distribution{d}).Pvalues = P;
        
        delta = eval(sprintf('air_%1$s_%2$s - behr_%1$s_%2$s',sites{s}, distribution{d}));
        [R_2, P_2] = corrcoef(eval(sprintf('[delta, aermax_%1$s_%2$s, aertot_%1$s_%2$s]',sites{s}, distribution{d})));
        corrStruct.delta.(sites{s}).(distribution{d}).R2 = R_2;
        corrStruct.delta.(sites{s}).(distribution{d}).Pvalues = P_2;
    end
end

putvar(corrStruct);

end


