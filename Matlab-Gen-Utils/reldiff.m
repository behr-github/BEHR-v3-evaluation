function [ varargout ] = reldiff( A, B, varargin )
%RDIFF = RELDIFF( A, B ) Computes the relative difference of A - B.
%   Convenience function to make it require less typing to compute a
%   relative difference between A and B as (A - B)./B. Multiply the result
%   by 100 to convert to percent difference.
%
%   RDIFF = RELDIFF( A, B, TRUE ) will reshape the output into a column
%   vector.
%
%   RDIFF = RELDIFF( A, B, 'avg' ) will use the average of A and B in the
%   denominator rather than B. This is good for percent differences where
%   neither A nor B is "right" or "first" (i.e. you're not doing a percent
%   change or percent error).
%
%   [RDIFF, SIGMA_DIFF] = RELDIFF( A, B, SIGMA_A, SIGMA_B, ___ ) will
%   propagate the uncertainty in A and B (SIGMA_A and SIGMA_B) through the
%   calculation. This works with either previous syntax, i.e. it can be
%   forced into a vector or used with the average in the denominator. This
%   assumes that the uncertainties can be treated using the formula s_f^2 =
%   (df/da s_a)^2 + (df/db s_a)^2.

E = JLLErrors;

forcevec = false;
denom_mode = 'first';
sigma_A = nan;
sigma_B = nan;

a=1;
while a <= numel(varargin)
    if isnumeric(varargin{a}) && isnumeric(varargin{a+1}) && isequal(size(varargin{a}), size(A)) && isequal(size(varargin{a+1}),size(B))
        sigma_A = varargin{a};
        sigma_B = varargin{a+1};
        a = a+1; % works with the later a=a+1 to advance by 2
    elseif (isnumeric(varargin{a}) || islogical(varargin{a})) && isscalar(varargin{a})
        forcevec = varargin{a};
    elseif ischar(varargin{a})
        denom_mode = varargin{a};
    else
        E.badinput('Optional arguments to RELDIFF must be either strings or scalar numbers/logicals.');
    end
    a = a+1;
end

if strcmpi(denom_mode, 'first')
    denom = B;
    uncert_fxn = @first_uncert;
elseif strcmpi(denom_mode, 'avg')
    denom = 0.5 * (A+B);
    uncert_fxn = @avg_uncert;
else
    E.badinput('The only denominator modes allowed for RELDIFF are ''first'' and ''avg''');
end

rdiff = (A - B)./denom;
if nargout > 1
    sigma_diff = uncert_fxn(A,B,sigma_A,sigma_B);
else
    sigma_diff = nan(size(rdiff));
end

if forcevec
    rdiff = rdiff(:);
    sigma_diff = sigma_diff(:);
end

varargout{1} = rdiff;
varargout{2} = sigma_diff;

end

function s = first_uncert(a,b,sa,sb)
% The uncertainty propagation for (a-b)/b is s^2 = [(1-b)/b * sa]^2 +
% {[b(a-1) - (a-b)]/b^2 * sb}^2.
s2 = ((1-b)./b .* sa).^2 + ((b.*(a-1) - (a-b))./b.^2 .* sb).^2;
s = sqrt(s2);
end

function s = avg_uncert(a,b,sa,sb)
% The uncertainty propagation for (a-b)/[0.5*(a+b)] is
% s^2 = [df/da sa]^2 + [df/db sb]^2, where
% df/da = 2*[(a+b)(1-b) - (a-b)(1+b)]/(a+b)^2
% df/db = 2*[(a+b)(a-1) - (a-b)(a+1)]/(a+b)^2
warning('This uncertainty propagation has not been verified');
df_da = 2 .* ((a+b).*(1-b) - (a-b).*(1+b)) ./ (a+b).^2;
df_db = 2 .* ((a+b).*(a-1) - (a-b).*(a+1)) ./ (a+b).^2;
s2 = (df_da .* sa).^2 + (df_db .* sb).^2;
s = sqrt(s2);
end

