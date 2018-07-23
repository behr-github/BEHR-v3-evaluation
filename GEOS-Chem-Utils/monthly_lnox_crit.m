function [ lnox_gt60 ] = monthly_lnox_crit( Prod )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

E=JLLErrors;

tVec = Prod(1).tVec;

if any(year(tVec) ~= year(tVec(1)))
    E.notimplemented('different years in one tVec')
end

yr=year(tVec(1));

% Find the dates that correspond to each month
m{1} = tVec >= datenum(sprintf('%04d-01-01',yr)) & tVec < datenum(sprintf('%04d-02-01',yr));
m{2} = tVec >= datenum(sprintf('%04d-02-01',yr)) & tVec < datenum(sprintf('%04d-03-01',yr));
m{3} = tVec >= datenum(sprintf('%04d-03-01',yr)) & tVec < datenum(sprintf('%04d-04-01',yr));
m{4} = tVec >= datenum(sprintf('%04d-04-01',yr)) & tVec < datenum(sprintf('%04d-05-01',yr));
m{5} = tVec >= datenum(sprintf('%04d-05-01',yr)) & tVec < datenum(sprintf('%04d-06-01',yr));
m{6} = tVec >= datenum(sprintf('%04d-06-01',yr)) & tVec < datenum(sprintf('%04d-07-01',yr));
m{7} = tVec >= datenum(sprintf('%04d-07-01',yr)) & tVec < datenum(sprintf('%04d-08-01',yr));
m{8} = tVec >= datenum(sprintf('%04d-08-01',yr)) & tVec < datenum(sprintf('%04d-09-01',yr));
m{9} = tVec >= datenum(sprintf('%04d-09-01',yr)) & tVec < datenum(sprintf('%04d-10-01',yr));
m{10} = tVec >= datenum(sprintf('%04d-10-01',yr)) & tVec < datenum(sprintf('%04d-11-01',yr));
m{11} = tVec >= datenum(sprintf('%04d-11-01',yr)) & tVec < datenum(sprintf('%04d-12-01',yr));
m{12} = tVec >= datenum(sprintf('%04d-12-01',yr)) & tVec < datenum(sprintf('%04d-01-01',yr+1));

% For each month, calculate the total of anthropogenic, biomass burning,
% and lightning NO emissions and find locations where lightning is >60% of
% the emissions.

% Assumming anthropogenic is 2, biomass burning 3, and lightning 6 in the
% production structure:

anthro = Prod(2).dataBlock;
biob = Prod(3).dataBlock;
lnox = Prod(6).dataBlock;

sz = size(anthro);

lnox_gt60 = false(sz(1),sz(2),12);

for a=1:12
    month_anthro = nansum2(nansum2(anthro(:,:,:,m{a}),3),4);
    month_biob = nansum2(biob(:,:,m{a}),3);
    month_lnox = nansum2(nansum2(lnox(:,:,:,m{a}),3),4);
    lnox_gt60(:,:,a) = (month_lnox ./ (month_anthro + month_biob + month_lnox)) >= 0.6;
end

end

