function chk = difftol(A, B, varargin)
%DIFFTOL Check if two numeric arrays are the same to within tolerance
%   CHK = DIFFTOL(A, B) Returns true for each element of |A - B| that is
%   less than eps() (epsilon, i.e. machine precision). CHK will be the
%   same size as A and B.
%
%   Parameters:
%   
%       'tol' - modify the default tolerance. Can be either a positive
%       scalar number, or an empty array; the empty array will skip the
%       tolerance test.
%
%       'reltol' - relative tolerance. Default is an empty array, i.e. this
%       test is not used. Can be any positive scalar number. Relative
%       difference is calculated with ABS(RELDIFF(A(:), B(:), 'avg')).
%
%       'nans' - determines how nans are handled. Default is 'same',
%       meaning that NaNs must be in the same locations in A and B or the
%       test returns false. Currently this is the only option.

E = JLLErrors;

p = inputParser;
p.addParameter('tol', eps);
p.addParameter('reltol', []);
p.addParameter('nans', 'same');

p.parse(varargin{:});
pout = p.Results;

tolerance = pout.tol;
relative_tolerance = pout.reltol;
nan_treatment = pout.nans;

%%%%%%%%%%%%%%%%%%
% INPUT CHECKING %
%%%%%%%%%%%%%%%%%%

if ~isnumeric(A) || ~isnumeric(B)
    E.badinput('A and B must both be numeric')
elseif ~isequal(size(A), size(B))
    E.badinput('A and B must be the same size');
end

if ~isempty(tolerance) && (~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance < 0)
    E.badinput('"tolerance" must be a positive scalar number, an empty array, or the string "default"')
end

if ~isempty(relative_tolerance) && (~isnumeric(relative_tolerance) || ~isscalar(relative_tolerance) || relative_tolerance < 0)
    E.badinput('"relative_tolerance" must be a positive scalar number or an empty array')
end

%%%%%%%%
% MAIN %
%%%%%%%%

chk = true(size(A));

% First, verify that the nans in the two arrays match
switch lower(nan_treatment)
    case 'same'
        nans_pass = ~xor(isnan(A), isnan(B));
    otherwise
        E.badinput('The parameter "nans" must be one of: ''same''');
end

chk = chk & nans_pass;

% Only if neither value is a NaNs do we check that the numerical difference
% is within the specified relative and absolute tolerance, because NaNs
% will always fail this test.
neither_nans = ~isnan(A) & ~isnan(B);
if ~isempty(tolerance)
    delta = abs(A - B);
    chk(neither_nans) = chk(neither_nans) & delta(neither_nans) <= tolerance;
end
if ~isempty(relative_tolerance)
    delta = abs(reldiff(A, B, 'avg'));
    chk(neither_nans) = chk(neither_nans) & delta(neither_nans) <= relative_tolerance;
end

end

