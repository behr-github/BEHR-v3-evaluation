function [ sf ] = shape_factor( pres, prof, lon, lat, month_in )
%SHAPE_FACTOR Calculates a shape factor using the Palmer 2001 definition
%   Palmer et al. 2001 (JGR Atmos. 106, p. 14539) specifies the shape
%   factor as the normalized number density (eq. 14). To compute shape
%   factors from mixing ratio profiles, that means we need to convert
%   mixing ratio to number density first, then divide by the VCD.
%
%   SF = SHAPE_FACTOR( PRES, PROF, LON, LAT, MONTH_IN ) calculates the
%   shape factor of PROF (given in mixing ratio, parts-per-part) defined at
%   pressures PRES assuming the temperature profile for MONTH_IN at LON and
%   LAT. PRES and PROF must be vectors of the same length; LON, LAT, and
%   MONTH_IN must all be scalar numbers.

E = JLLErrors;
if ~isvector(pres) || ~isnumeric(pres)
    E.badinput('PRES must be a numeric vector')
elseif ~isvector(prof) || ~isnumeric(prof)
    E.badinput('PROF must be a numeric vector')
elseif ~isscalar(lon) || ~isnumeric(lon)
    E.badinput('LON must be a numeric vector')
elseif ~isscalar(lat) || ~isnumeric(lat)
    E.badinput('LAT must be a numeric vector')
elseif ~isscalar(month_in) || ~isnumeric(month_in)
    E.badinput('MONTH_IN must be a numeric vector')
end

ndens = ndens_from_pres(pres, lon, lat, month_in);

% integPr2 takes mixing ratio profiles and integrates over pressure to
% return a VCD in molec. cm^-2.
vcd = integPr2(prof,pres,max(pres));
prof = ndens .* prof;
sf = prof ./ vcd;

end

