function [ alt ] = nearest_GLOBE_alt( lon, lat )
%nearest_GLOBE_alt Returns the GLOBE altitude nearest the input lat/lon

%[globe_lon_matrix, globe_lat_matrix, terpres] = GLOBE_terrain([-125,-65],[25,50]);

% Retrieve globe data in a 1 degree box around the chosen lat/lon
[globe_lon_matrix, globe_lat_matrix, terpres] = GLOBE_terrain([lon-0.5,lon+0.5],[lat-0.5,lat+0.5]);

% The GLOBE database writes NaNs over the ocean, replace these with 0s
% (after all, the ocean should be at sea level)
terpres(isnan(terpres))=0;

% Find the point closest to the input lat/lon and return its altitude
dlon = abs(globe_lon_matrix - lon); dlat = abs(globe_lat_matrix - lat);
dtotal = dlon+dlat;
xx = (dtotal(:) == min(dtotal(:)));

alt = terpres(xx);

end

