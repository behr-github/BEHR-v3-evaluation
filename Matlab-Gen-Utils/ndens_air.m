function [ N ] = ndens_air( T, P )
%NDENS_AIR Calculate the number density of air for given T and P
%   N = NDENS_AIR( T, P ) returns the number density of air in molec./cm^3
%   for temperature T (in Kelvin) and P (in hPa).

R = 8.314e3; % in cm^3 kPa K^-1 mol^-1
R = R*10;
Av = 6.022e23; % avogadro's number
N = P ./ (R .* T) * Av;

end

