function [ del_theta ] = angle_diffd( theta_1, theta_2 )
%ANGLE_DIFF Computes the smallest difference between two arrays of angles
%   [ DEL_THETA ] = ANGLE_DIFF( THETA_1, THETA_2 ) DEL_THETA is an array
%   the same size as THETA_1 and THETA_2 where DEL_THETA(i) is the value
%   with the smallest magnitude out of:
%       THETA_1(i) - THETA_2(i)
%       THETA_1(i) - THETA_2(i) + 360
%       THETA_1(i) - THETA_2(i) - 360
%
%   THETA_2 may be a scalar, in which case that value is used for all
%   indices i.
%
%   The sign of the value is retained; in this way, the issue of dealing
%   with differences in angle near 0 or 360 is dealt with. Angles are to be
%   given in degrees between 0 and 360.

E = JLLErrors;

if isscalar(theta_2)
    theta_2 = repmat(theta_2, size(theta_1));
end

if ~isequal(size(theta_1), size(theta_2))
    E.badinput('THETA_1 and THETA_2 must be the same size, unless THETA_2 is a scalar');
elseif ~isnumeric(theta_1) || ~isnumeric(theta_2)
    E.badinput('THETA_1 and THETA_2 must be numeric');
end
check360 = any(theta_1(:) < 0 | theta_1(:) >= 360) || any(theta_2(:) < 0 | theta_2(:) >= 360);
check180 = any(theta_1(:) < -180 | theta_1(:) > 180) || any(theta_2(:) < -180 | theta_2(:) > 180);
if check180 && check360        
    E.badinput('All values of THETA_1 and THETA_2 must be >= 0 and < 360 or all must be >= -180 and <= 180');
end

d1 = theta_1 - theta_2;
d2 = theta_1 - theta_2 + 360;
d3 = theta_1 - theta_2 - 360;

del_theta = nan(size(d1));

for a=1:numel(del_theta)
    d_a = [d1(a), d2(a), d3(a)];
    [~,m] = min(abs(d_a));
    del_theta(a) = d_a(m);
end

end

