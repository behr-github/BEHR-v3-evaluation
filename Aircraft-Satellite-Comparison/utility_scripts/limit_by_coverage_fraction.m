function [ behrno2 ] = limit_by_coverage_fraction( db_iall, frac_crit )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

narginchk(2,2);

behrno2 = nan(numel(db_iall.coverage_fraction),1);
for a=1:numel(db_iall.coverage_fraction)
    xx = db_iall.coverage_fraction{a} > frac_crit;
    behrno2(a) = nanmean(db_iall.all_behr{a}(xx));
end

end

