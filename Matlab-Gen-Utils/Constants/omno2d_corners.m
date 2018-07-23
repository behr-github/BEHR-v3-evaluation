function [ omno2_lon, omno2_lat ] = omno2d_corners( isgrid )
%OMNO2D_CORNERS Outputs grids of corner points for OMNO2d gridded NO2
%   Set the input "isgrid" to false to only return vectors, since this is
%   on a grid. Otherwise, it will assume you need the full grid.

E=JLLErrors;

% INPUT CHECKING
if ~exist('isgrid','var')
    isgrid = true;
else
    try 
        isgrid=logical(isgrid);
    catch
        E.badinput('isgrid (if given) must be convertable to a logical value (true/false)')
    end
    if ~isscalar(isgrid)
        E.badinput('isgrid must be a scalar logical')
    end
end

omno2_lon = -180:0.25:180;
omno2_lat = -90:0.25:90;

if isgrid
    [omno2_lat, omno2_lon] = meshgrid(omno2_lat, omno2_lon);
end

end

