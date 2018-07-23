function [ vq ] = angle_interp1d( x, v, xq )
%ANGLE_INTERP1D Interpolates angles for 1D inputs
%   VQ = ANGLE_INTERP1D( X, V, XQ ) For a vector of angles V,
%   related to points X, interpolate to points XQ. Angles are
%   to be given in degrees, and will be returned in degrees
%   on the range (-180, 180].

v_rad = unwrap(v*pi/180);
vq = interp1(x, v_rad, xq);
vq = vq * 180/pi;
% If unwrap results in angles > 360 or < 0, this will put everything on [0,
% 360)
vq = mod(vq,360);
% Then transform back to quadrants III and IV being represented by negative
% angles
vq(vq>180) = vq(vq>180)-360;

end

