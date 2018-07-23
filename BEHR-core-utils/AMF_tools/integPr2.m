% INTEGPR2 Integrate mixing ratios over pressure
%
%   VCD = INTEGPR2(MIXINGRATIO, PRESSURE, PRESSURESURFACE) Integrates the
%   profile given by MIXINGRATIO (which must be unscaled, i.e. mol/mol or
%   volume/volume) defined at the vertical levels given by PRESSURE (which
%   must be monotonically decreasing and given in hPa). PRESSURESURFACE is
%   the lower level of integration, must be a scalar also given in hPa. The
%   upper limit is taken as the minimum value in PRESSURE.
%
%   VCD = INTEGPR2( ___, PRESSURETROPOPAUSE) will integrate up to
%   PRESSURETROPOPAUSE instead of the minimum value in PRESSURE.
%
%   [ VCD, P_OUT, F_OUT ] = INTEGPR2( ___, 'interp_pres', INTERP_PRES )
%   Allows you to specify one or more pressures as the vector INTERP_PRES
%   that the mixing ratio vector will be interpolated to. P_OUT will be
%   PRESSURES with the INTERP_PRES inserted and F_OUT will be the vector of
%   mixing ratios with extra entries corresponding to the INTERP_PRES
%   pressures. Note that the size of P_OUT and F_OUT are not guaranteed to
%   be numel(pressures) + numel(interp_pres); if any of the INTERP_PRES
%   pressures are already present in PRESSURES, then no extra level will be
%   added for that INTERP_PRES value.
%
%   [ ___ ] = INTEGPR2( ___, 'fatal_if_nans', true ) will cause INTEGPR2 to
%   throw an error if any of the mixing ratio values between the lower and
%   upper integration limits are NaNs. This allows you to ensure that your
%   profile is fully defined over all relevant pressures. If this parameter
%   is not given, a warning will still be issued if NaNs are detected.

% Legacy comment block from Ashley Russell:
%..........................................................................
% Integrates vector of mixing ratios (cm-3) above pressureSurface as function of pressure (hPa)
% to get vertical column densities (cm-2). Computes piecewise in layers between two pressures.
% Assumes exponential variation between mixing ratios that are positive or zero (zero is treated as 1.e-30).
% If one or both mixing ratios is negative, assumes constant (average=(f1+f2)/2) mixing ratio in the layer.
% The uncertainties are computed assuming constant mixing ratio in each layer (using exponential
% variation can lead to strange (large) results when f1*p1 is approx f2*p2. This sometimes happens
% when using averaging kernels).
%
% Arguments (vector indices are i = 0,1,2,...n-1):
%
%  mixingRatio(i)   = input volume mixing ratios (no units)
%  pressure(i)      = vector of input pressures from largest to smallest (hPa)
%  pressureSurface  = surface pressure (where integration starts) (hPa)
%  interpPres       = pressures to interpolate the output profile to (hPa) ( added by JLL for BEHR scattering weights )
%
% Restrictions:
%  Pressures must be greater than zero and monotonically decreasing
%
% Equation for exponential variation:
%  For a segment of the vertical column between ambient pressures p1 and p2,
%  at which the trace gas volume mixing ratios are f1 and f2, the vertical
%  column density VCD (assuming number density varies exponentially with p) is
%
%    VCD   =  (f1 * p1  -  f2 * p2) /  [(b + 1) * m * g]
%
%    where:   b = ALOG( f2/f1 ) / ALOG( p2/p1 )
%
%   Also, for interpolation of f between p1 and p2:
%
%   f  =  f1 * (p/p1)^b
%
% Note:  Can get large uncertainties when  f1*p1 = f2*p2
%
%..........................................................................

function [vcd, p_out, f_out] = integPr2(mixingRatio, pressure, pressureSurface, varargin)

E = JLLErrors;
p = inputParser;
p.addOptional('pressureTropopause', [], @(x) isscalar(x) && isnumeric(x) && (isnan(x) || x > 0));
p.addParameter('interp_pres', []);
p.addParameter('fatal_if_nans', false);

p.parse(varargin{:});
pout = p.Results;

pressureTropopause = pout.pressureTropopause;
interpPres = pout.interp_pres;
fatal_if_nans = pout.fatal_if_nans;

if any(pressure<0)
    E.badinput('PRESSURE must be all >= 0')
elseif any(diff(pressure)>0)
    E.badinput('PRESSURE must be monotonically decreasing')
end

if ~isscalar(pressureSurface)
    E.badinput('PRESSURESURFACE must be a scalar')
end

if isempty(pressureTropopause)
   pressureTropopause = min(pressure);
end

if isempty(interpPres)
    if nargout > 1
        E.callError('nargout','Without any interpPres values, p_out and f_out will not be set');
    end
else
    % Make sure the interpolation pressure is within the pressure vector
    % given.
    interpPres = clipmat(interpPres, min(pressure), max(pressure));
end


%   mean molecular mass (kg)  *  g (m/s2)  *  (Pa/hPa)   *   (m2/cm2)
mg  = (28.97/6.02E23)*1E-3    *   9.8      *    1E-2     *     1E4;

fmin = 1E-30;

vcd  = 0;
f    = mixingRatio;
p    = pressure;
p0   = max(pressureSurface, min(pressure));
n    = numel(p);

dvcd     = 0;
deltaVcd = zeros(numel(f),1); % Changed to make these vectors on 9/26/2014 JLL
df       = zeros(numel(f),1);



numIP = numel(interpPres);
if numIP > 0
    if ~iscolumn(p)
        p_out = p'; 
    else
        p_out = p;
    end
    f_out = f;
    for a=1:numIP
        if all(p~=interpPres(a))
            f_i = interpolate_surface_pressure(p,f,interpPres(a));
            bottom = p_out > interpPres(a);
            top = p_out < interpPres(a);
            p_out = [p_out(bottom); interpPres(a); p_out(top)];
            f_out = [f_out(bottom); f_i; f_out(top)];
        end
    end
end

% If the surface pressure is above the tropopause pressure, i.e. the lower
% integration limit is above the upper integration limit, return 0 b/c
% unlike abstract integration where reversing the integration limits
% reverses the sign, here, physically, if the surface pressure is above the
% top limit, then there is no column between the surface and the top limit.
% With this added here, we might want to remove the "pCld > pTropo" tests
% in omiAmfAK2.
if pressureSurface <= pressureTropopause
    vcd = 0;
    return
end

f0 = interpolate_surface_pressure(p,f,p0);

p(i0) = p0;
f(i0) = f0;
i_initial = i0;

if ~isnan(pressureTropopause)
    p_end   = max(pressureTropopause, min(pressure));
    f_end = interpolate_surface_pressure(p,f,p_end);
    p(i0+1) = p_end;
    f(i0+1) = f_end;
    i_end = i0;
else 
    i_end = n-1;
end


% Integrate................................................................
if isnan(pressureSurface)
    vcd = nan;
    return
end

for i = 1:n-1
    deltaVcd(i) = (f(i) + f(i + 1)) .* (p(i) - p(i + 1)) ./ (2 * mg);  %assume const mixing ratio in each layer
    b = (log(max(f(i + 1),fmin)) - log(max(f(i),fmin))) ./ (log(p(i + 1)) - log(p(i)));
    
    if f(i) >= 0 && f(i + 1)>=0 && abs(b + 1) >= 0.01;
        deltaVcd(i) = (f(i)*p(i) - f(i + 1)*p(i + 1)) ./ (b + 1) ./ mg;    %assume exponential variation in each layer
    end
    
end

if any(isnan(deltaVcd(i_initial:i_end)))
    if fatal_if_nans
        E.callError('nans_in_column', 'NaNs detected in partial columns.')
    else
        warning('integPr2:nans_in_column', 'NaNs detected in partial columns. They will not be added into the total column density.')
    end
end
vcd = nansum2(deltaVcd(i_initial:i_end));


    function [f0] = interpolate_surface_pressure(p_in,f_in,p0_in)
        pp=find(p_in>=p0_in); if isempty(pp); pp=0; end
        i0 = min(max((max(pp)),1),(n - 1));
        f0 = interp1(log([p_in(i0),p_in(i0 + 1)]), [f_in(i0),f_in(i0 + 1)], log(p0_in),'linear','extrap');   %assume linear variation in surface layer
        
        b = (log(max(f_in(i0 + 1),fmin)) - log(max(f_in(i0),fmin))) ./ (log(p_in(i0 + 1)) - log(p_in(i0)));
        if f_in(i0) >= 0 && f_in(i0 + 1) >= 0 && abs(b + 1) >= 0.01;
            f0 = f_in(i0) * ((p0_in./p_in(i0))^b);                                    %assume exponential variation in surface layer
        end
        
    end
end
