%%rNmcTmp2
%%arr 07/23/2008

%..........................................................................
% Returns vector of temperatures on input pressure grid,
% for nearest longitude and latitude in geographic grid
%
% New, simplified version (without write capability) that 
% also interpolates profile onto input pressure grid (2008-05-30)
%
% Inputs:
%  fileTmp  = Complete directory/filename of temperature file
%  pressure = vector of pressures (hPa) monotonically decreasing
%  lon and lat are scalars (deg)
%
% Outputs (optional: Return these to avoid having to re-read fileTmp on each call):
%  presSAVE = pressure grid (vector) from fileTmp (hPA)
%  lonSAVE = vector of longitudes from fileTmp (deg)
%  latSAVE = vector of latitudes from fileTmp (deg)
%  monSAVE = vector of month numbers
%  tmpSAVE = large array of temperatures from fileTmp
%
%..........................................................................

function [temperature, tmpSAVE] = rNmcTmp2(fileTmp,  pressure, lon, lat, mon)

% Longitude in the temperature profile file is defined as degrees east of
% 0, so -90 in OMI data == 270 here, while +90 in OMI data == +90 here.
lon=mod(lon,360);

if any(diff(pressure) > 0);
    E.badinput('pressure must be a monotonically decreasing vector of pressures');
else
    
    if exist('tmpSAVE','var')==0;
    
        fid = fopen(fileTmp,'r');
        
        header = fgets(fid);

        nPres = fscanf(fid,'%g',[1 1]);
        
        LON = fscanf(fid,'%g',[1 3]);
        nLon = LON(1);  lonStart = LON(2);  lonEnd = LON(3);
        
        LAT = fscanf(fid,'%g',[1 3]);
        nLat = LAT(1);  latStart = LAT(2);  latEnd = LAT(3);
        
        MON = fscanf(fid,'%g',[1 3]);
        nMon = MON(1);  monStart = MON(2);  monEnd = MON(3);
         
        presSAVE = zeros(nPres,1);                  %pressure grid
        
        presSAVE = fscanf(fid,'%g',[1 18]);
 
        tmpSAVE = zeros(nPres, nLon, nLat, nMon);   %temperature array

        tmpSAVE_vec = fscanf(fid,'%g',inf);
        tmpSAVE = reshape(tmpSAVE_vec, [nPres nLon nLat nMon]);

        lonSAVE  = lonStart + ((0:nLon-1) + 0.5) * (lonEnd - lonStart) ./ nLon;  %%GOOD
        latSAVE  = latStart + ((0:nLat-1) + 0.5) * (latEnd - latStart) ./ nLat; %%GOOD
        monSAVE  = monStart + ((0:nMon-1) + 0) * (monEnd - monStart) ./ nMon; 
    end

    
    %Interpolate onto input pressure grid
    temperature_interp_horiz = nan([numel(presSAVE), size(lon)]);
    for pres_i=1:length(presSAVE)
        temperature_slice = reshape(tmpSAVE(pres_i,:,:,:), [nLon, nLat, nMon]);
        temperature_interp_horiz(pres_i,:,:)=(interpn(lonSAVE, latSAVE, monSAVE, temperature_slice, lon, lat, mon,'linear')); 
    end
    temperature=exp(interp1(log(presSAVE),log(temperature_interp_horiz), log(pressure),'linear','extrap'));
end

status = fclose(fid);
