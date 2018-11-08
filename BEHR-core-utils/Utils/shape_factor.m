function [ sf, pres ] = shape_factor( pres, prof, varargin )
%SHAPE_FACTOR Calculates a shape factor using the Palmer 2001 definition
%   Palmer et al. 2001 (JGR Atmos. 106, p. 14539) specifies the shape
%   factor as the normalized number density (eq. 14). However, there are
%   other ways that we can calculate the shape factor.
%
%   SF = SHAPE_FACTOR( PRES, PROF ) calculates shape factor of PROF (given
%   in mixing ratio, parts-per-part) defined at pressures PRES by
%   converting PROF to number density assuming a temperature profiles for
%   the middle of the US in June.
%
%   SF = SHAPE_FACTOR( PRES, PROF, LON, LAT, MONTH_IN ) calculates the
%   shape factor of PROF (given in mixing ratio, parts-per-part) defined at
%   pressures PRES assuming the temperature profile for MONTH_IN at LON and
%   LAT. PRES and PROF must be vectors of the same length; LON, LAT, and
%   MONTH_IN must all be scalar numbers.
%
%   Additionally, you may specify a mode to calculate the shape factor
%   using the parameter 'mode':
%
%       'ndens' - the default, converts the mixing ratio to a number
%       density using the given pressure and a temperature profile looked
%       up from the NASA coarse temperature climatology that was used in
%       v2.1C and earlier of BEHR, then divides this number density by the
%       VCD to obtain the shape factor.
%
%       'mixing ratio' - simply divides the given mixing ratio by the VCD.
%       The units don't cancel, and technically this means that some
%       effects due to the different in number density of air from the
%       surface to the tropopause are neglected, but this is in many ways
%       the closest to how the profile is treated in BEHR.
%
%       'partial col' - this integrates partial columns from model level i
%       to i+1, for i = 1:length(prof)-1, then divides by the VCD. This
%       results in a shape factor one shorter than the input. A vector of
%       pressures is returned as the second output, calculated as the
%       midpoints between the original levels.

E = JLLErrors;
p = advInputParser;
p.addOptional('lon',-95);
p.addOptional('lat',37.5);
p.addOptional('month',6);
p.addParameter('mode','ndens');

p.parse(varargin{:});
pout = p.Results;
lon = pout.lon;
lat = pout.lat;
month_in = pout.month;
calc_mode = pout.mode;

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

notnans = ~isnan(pres);
% integPr2 takes mixing ratio profiles and integrates over pressure to
% return a VCD in molec. cm^-2.
vcd = integPr2(prof(notnans),pres(notnans),max(pres),min(pres));
switch lower(calc_mode)
    case 'mixing ratio'
        sf = prof ./ vcd;
    case 'ndens'
        ndens = ndens_from_pres(pres, lon, lat, month_in);
        prof = ndens .* prof;
        sf = prof ./ vcd;
    case 'partial col'
        pc_prof = nan(size(prof));
        for i = 1:length(prof)-1
            j = i + 1;
            if ~isnan(pres(i)) && ~isnan(pres(j))
                pc_prof(i) = integPr2(prof(i:j), pres(i:j), pres(i), pres(j));
            end
        end
        sf = pc_prof(1:end-1) ./ vcd;
        pres = (pres(1:end-1) + pres(2:end))/2;
end
end

