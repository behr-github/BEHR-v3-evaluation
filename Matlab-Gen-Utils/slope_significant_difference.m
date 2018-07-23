function [is_significant, t] = slope_significant_difference(P1, s1, n1, P2, s2, n2)
% SLOPE_SIGNIFICANT_DIFFERENCE Test if two slopes are different from each other
%
%   [ IS_SIG, T ] = SLOPE_SIGNIFICANT_DIFFERENCE( P1, S1, N1, P2, S2, N2 )
%   P is the two element [slope intercept] vector for each fit, S is the 
%   standard deviation of each slope, and N is the number of points fit
%   for each slope. IS_SIG will be true if the two slopes are significantly
%   different. T will be a vector [t_calc, t_table] of the calculated T
%   values and the T value for (N1 + N2 - 4) degress of freedom.

ssr1 = (s1.^2)*(n1-1);
ssr2 = (s2.^2)*(n2-1);
[t(1),t(2),is_significant] = two_sample_t_test(P1(1), ssr1, n1, P2(1), ssr2, n2, n1 + n2 - 4);


end

