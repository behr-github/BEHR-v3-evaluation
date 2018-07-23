function [ ka_kb ] = atkinson_branching( ncarbon, ndens_air, temperature )
%ATKINSON_BRANCHING Estimates a branching ratio for the number of carbons
%   GEOS-Chem's SMVGEARII solver seems to treat alkyl nitrate branching not
%   by a direct relation of the two rate constants, but by a ratio defined
%   by Eq. II in Atkinson, R., 1990, Atmos. Environ., 244, (1), 1-41.  This
%   uses the number of carbons (along with number density of air and
%   temperature) to determine the ratio of ka/kb, where ka is the formation
%   of ANs and kb is the production of RO and NO2.
%
%   ATKINSON_BRANCHING( NCARBON, NDENS_AIR, TEMPERATURE ) requires the
%   number of carbons, the number density of air in molec. cm^-3, and the
%   temperature in K.

% Input checking, and rename to shorten following lines
E = JLLErrors;

scalar_test = [~isscalar(ncarbon), ~isscalar(ndens_air), ~isscalar(temperature)];
if sum(scalar_test) > 1
    E.badinput('Only one of the inputs should be multiple valued, the others should be scalars')
end

if any([~isnumeric(ncarbon), ~isnumeric(ndens_air), ~isnumeric(temperature)])
    E.badinput('All inputs should be numeric')
end

T = temperature;
Cair = ndens_air;


% Calculation %
alpha = 1.94e-22;
beta = 0.97;
Y300_0 = alpha .* exp(beta .* ncarbon);
Y300 = 0.826;
m_0 = 0;
m_inf = 8.1;

X = Y300_0 .* Cair .* (T ./ 300) .^ -m_0;
Y = Y300 .* (T ./ 300) .^ -m_inf;

z = (1 + (log10(X ./ Y)) .^ 2) .^ -1;
F = 0.411;

temp = ( X ./ (1 + X ./ Y) ) .* F .^ z;

ka_kb = temp ./ (1 + temp);


end

