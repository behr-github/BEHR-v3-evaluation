function [t_calc, t_table, is_different] = two_sample_t_test(val1, s1, n1, val2, s2, n2, n_dofs, varargin)
%TWO_SAMPLE_T_TEST Test whether two values are different based on the two value t-test.
%   [ T_CALC, T_TABLE, IS_DIFFERENT ] = TWO_SAMPLE_T_TEST( VAL1, SSR1, N1,
%   VAL2, SSR2, N2, N_DOFS ) Does a two-sample t-test on two values (VAL1
%   and VAL2) for their squared standard deviations time their degrees of
%   freedom (S1 and S2) and the number of data points involved in each
%   measurement (N1 and N2). N_DOFS is the number of degrees of freedom,
%   which is usually N1 + N2 minus the number of parameters fitted to the
%   data. T_CALC is the T value calculated for the two sample populations,
%   T_TABLE is the critical T value for, by default, a 95% confidence
%   two-sided t-test. IS_DIFFERENT is a logical value, indicating whether
%   the two samples are different at the 95% confidence level (i.e. t_calc
%   > t_table).
%
%   Additional parameters:
%       'confidence' - the confidence level, expressed as a decimal. 0.95 by default, i.e. 95%
%       confidence level.
%
%       'sided' - whether to use a one sided or two sided t-test. Default is 'two', must be either
%       'one' or 'two'.

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

E = JLLErrors;

p = inputParser;
p.addParameter('confidence', 0.95);
p.addParameter('sided', 'two');

p.parse(varargin{:});
pout = p.Results;

confidence_level = pout.confidence;
t_sidedness = pout.sided;

if ~isnumeric(val1) || ~isscalar(val1)
    E.badinput('VAL1 must be a numeric scalar')
end
if ~isnumeric(s1) || ~isscalar(s1)
    E.badinput('SSR1 must be a numeric scalar')
end

if ~isnumeric(n1) || ~isscalar(n1) || n1 < 1 || isnan(n1)
    E.badinput('N1 must be a positive, numeric scalar')
end


if ~isnumeric(val2) || ~isscalar(val2)
    E.badinput('VAL2 must be a numeric scalar')
end
if ~isnumeric(s2) || ~isscalar(s2)
    E.badinput('SSR2 must be a numeric scalar')
end
if ~isnumeric(n2) || ~isscalar(n2) || n2 < 1 || isnan(n2)
    E.badinput('N2 must be a positive, numeric scalar')
end


if ~isnumeric(n_dofs) || ~isscalar(n_dofs) || n_dofs < 0 || isnan(n_dofs)
   E.badinput('N_DOFS must be a positive, numeric scalar')
end


if ~isnumeric(confidence_level) || ~isscalar(confidence_level) || confidence_level < 0 || confidence_level > 1
   E.badinput('"confidence_level" must be a scalar number between 0 and 1');
end

allowed_sidedness = {'one','two'};
if ~ischar(t_sidedness) || ~ismember(t_sidedness, allowed_sidedness)
    E.badinput('"sided" must be one of: ''%s''', strjoin(allowed_sidedness, ''', '''));
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

s_pooled = sqrt( (s1 + s2) ./ n_dofs );
n_pooled = sqrt( n1 .* n2 / (n1 + n2) );
t_calc = abs( val1 - val2 ) ./ s_pooled .* n_pooled;

if strcmpi(t_sidedness, 'two')
    % tinv gives 1-sided t values. For a two sided t-test (i.e. one where
    % the values may be above or below each other) the same percent chance
    % of being different must be split onto either side. Concretely, for a
    % 95% two-sided confidence interval, we need to compare against the
    % 97.5% t-value returned by tinv.
    confidence_level = 1 - ( (1 - confidence_level)/2 );
end

t_table = tinv(confidence_level, n_dofs);

is_different = t_calc > t_table;

end

