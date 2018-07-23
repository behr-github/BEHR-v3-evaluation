function [ ndens ] = ndens_from_pres( pres, varargin )
%NDENS_FROM_PRES Estimate the number density of air at a given pressure, time, and location
%   NDENS = NDENS_FROM_PRES( PRES, MONTH_IN ) Given an array of pressures
%   PRES and a location as LON and LAT and a month, MONTH_IN (numerical,
%   i.e. 1 = Jan, 2 = Feb, etc), estimate what the number density of air at
%   each pressure is. LON, LAT, and MONTH_IN must be scalar values which
%   will be used for all pressures. PRES must be in hPa.
%
%   Internally, this uses the temperature profiles from NASA SP v2 (used in
%   BEHR) to figure out what the temperature should be for each pressure.
%   Alternatively:
%
%   NDENS = NDENS_FROM_PRES( PRES, TEMP ) calculates number density from
%   given pressure and temperature. Pressure should be given in hPa,
%   temperature in Kelvin. This is mainly a convenience function so I can
%   stop looking up the gas constant.

E = JLLErrors;

if nargin == 4
    lon = varargin{1};
    lat = varargin{2};
    month_in = varargin{3};
    if ~isnumeric(pres)
        E.badinput('PRES must be numeric');
    elseif ~isnumeric(lon) || ~isscalar(lon)
        E.badinput('LON must be a scalar number');
    elseif ~isnumeric(lat) || ~isscalar(lat)
        E.badinput('LAT must be a scalar number');
    elseif ~isnumeric(month_in) || ~isscalar(month_in)
        E.badinput('MONTH_IN must be a scalar number');
    end
elseif nargin == 2
    temperature = varargin{1};
    if ~isequal(size(pres), size(temperature))
        E.badinput('With 2 arguments, TEMP must be the same size as PRES')
    end
    
    if any(temperature(:) < 200)
        warning('Values of TEMP < 200 suggests it is not given in Kelvin.')
    end
end

if any(pres(:) > 1100)
    warning('Pressure values > 1100 suggest that pressure is not in hPa, which this function assumes')
elseif all(pres(:) < 100)
    warning('Pressure values < 100 suggest that pressure is not in hPa, unless you are working with stratospheric pressures')
end

%%%%%%%%
% MAIN %
%%%%%%%%

if nargin == 4
    temperature_file = fullfile(behr_paths.amf_tools_dir, 'nmcTmpYr.txt');
    temperature = rNmcTmp2(temperature_file, pres(:), lon, lat, month_in);
    temperature = reshape(temperature, size(pres));
end
% From the ideal gas law, PV = nRT => n/V = P/RT. R = 8.3144598e3 cm^3 kPa
% K^-1 mol^-1 (wikipedia) so convert:
%
%   8.3144598e3 cm^3 kPa   10 hPa          mol
%   -------------------- * ------ * ------------------
%         K mol              kPa    6.022e23 molecules

R = 8.3144598e3 * 10 / 6.022e23;

ndens = pres ./ (R .* temperature);
end

