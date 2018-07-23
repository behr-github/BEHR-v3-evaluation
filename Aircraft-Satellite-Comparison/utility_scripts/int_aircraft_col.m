function [ gas_column, column_error ] = int_aircraft_col( altbins, gasbins, gasstderr, tempbins, presbins, varargin )
%int_aircraft_no2 Integrate trace gase measurements from aircraft
%   Detailed explanation goes here

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VARIABLE CHECK & SETUP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(4,8);

% 4 inputs: alt, gas, stderr, temperature. Pressures will be converted from
% altitude.
% 5 inputs: alt, gas, stderr, temperature, pressure. All used directly
% 6 inputs: alt, gas, stderr, lon, lat, month. Temperature will be read
% from the OMI temperature profiles for lon, lat, month. Pressures will be
% converted from altitude.
% 7 inputs: alt, gas, stderr, lon, lat, month, pressures. Temperature will be read
% from the OMI temperature profiles for lon, lat, month.

if nargin == 6 || nargin == 7;
    % Shift the variables around to make the names more sensible.  Clear
    % the original ones afterwards.
    lon = tempbins;
    tempbins = []; 
    lat = presbins;
    presbins = [];
    month = varargin{1};
    
    if ~isscalar(lon) || ~isscalar(lat) || ~isscalar(month)
        error(E.badinput('If not inputting a temperature profile, the fourth, fifth, and sixth arguments must be scalar values: lon, lat, and numeric month.'));
    end
end

if nargin == 7;
    presbins = varargin{2};
end

if nargin == 4 || nargin == 6
    % If no pressure bin centers passed, calculate the pressures from the
    % altitude values (in km), assuming a sea surface pressure of 1013 hPa
    % and a scale height of 7.4 km
    presbins = 1013 .* exp(altbins ./ 7/4);
end

% If we're getting temperature based on location and time, do so now.
if isempty(tempbins)
    [~,file_T] = amf_filepaths;
    tempbins = rNmcTmp2(file_T,presbins,lon,lat,month);
end

% If the user did not pass a std. error vector, make it a fill value
if isempty(gasstderr)
    gasstderr = -127*ones(size(gasbins));
end

% Check that the gas bins are in some sensible mixing ratio, i.e. less than
% a couple hundred ppm.  
if any(gasbins > 1e-3) || any(gasstderr > 1e-3)
    error(E.badinput('Trace gas profiles expected in mixing ratios of part/part. Please make the appropriate conversions'));
end

% Now check: the five expected vectors should be the same length
are_vec = [isvector(altbins), isvector(gasbins), isvector(gasstderr), isvector(tempbins), isvector(presbins)];
if any(~are_vec)
    error(E.badinput('The first five inputs must be vectors'));
end

vec_lens = [numel(altbins), numel(gasbins), numel(gasstderr), numel(tempbins), numel(presbins)];
if any(vec_lens~=vec_lens(1))
    error(E.badinput('The input vectors are expected to have the same number of elements'));
end


% Define Avogadro's #, the gas constant, and the integration segment
Av = 6.022e23; % molec. / mol
R = 8.314e4; % (hPa * cm^3) / (mol * K)
dz = 1; % integration segment in meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION - COLUMN DENSITY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate the NO2, temperature, and pressure data
alt_profile = altbins(1):(dz/1000):altbins(end);
gas_profile = interp1(altbins,gasbins,alt_profile,'linear');
temp_profile = interp1(altbins,tempbins,alt_profile,'linear');
pres_profile = exp(interp1(altbins,log(presbins),alt_profile,'linear')); % Linearly interpolate ln(P) since that is what depends linearly on altitude

% Fill in the standard error nans, but we don't need to interpolate
% to every 100 cm altitude point.
[~,~,gasstderr] = fill_nans(altbins,gasbins,gasstderr,'noclip');

% Carry out the numerical integration,
gas_column = 0;
for z=1:numel(alt_profile)
    P_z = pres_profile(z); T = temp_profile(z); gas_z = gas_profile(z);
    conc_gas = (Av * P_z * gas_z )/(R * T); % molec./cm^3, mixing ratio of NO2 is in pptv
    gas_column = gas_column + conc_gas * dz * 100; % Integrating in 100dz cm increments
end

% Double check that this number is a scalar. This will catch if a
% matrix (rather than a vector) slips into this calculation.
if ~isscalar(gas_column)
    error(E.badvartype(gas_column,'scalar'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION - ERROR PROPAGATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The error propagation will be handled by considering this to be
% effectively a trapezoid rule implementation.  The math behind
% this is in my BEHR notebook from 4 Feb 2015 - Josh Laughner

% Number density of air at each bin center
Nair = (presbins .* Av)./(R .* tempbins);
% Calculates the error for each trapezoid using vectors, then sum
% up that vector and square root it to find the total column error.
% Since the concentrations are defined in molec./cm^3, we need to
% convert altitude bins from km to cm (hence the factor of 1e5).
% Also, NO2 values are usually reported (by us at least) in pptv,
% hence the conv_fact defaults to 1e-12.
column_error = sqrt( sum( ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* (gasstderr(2:end)).^2 .* Nair(2:end).^2 +...
    ((altbins(2:end) - altbins(1:end-1))*1e5/2).^2 .* (gasstderr(1:end-1)).^2 .* Nair(1:end-1).^2 ) );
% Double check that this number is a scalar. This will catch if a
% matrix (rather than a vector) slips into this calculation.
if ~isscalar(column_error)
    error(E.badvartype(column_error,'scalar'));
end

end

