function [ mean_theta ] = angle_meand( theta, dim )
%ANGLE_MEAND Computes the mean of the angles in THETA
%   MEAN_THETA = ANGLE_MEAND( THETA ) computes the mean angle of all the
%   angles in the array THETA. If THETA is a vector, a single value is
%   returned. Otherwise this operates along the first non-singleton
%   dimension.
%
%   MEAN_THETA = ANGLE_MEAND( THETA, DIM ) allows you to specify which
%   dimension to operate along, DIM.
%
%   The average is computed as atan2d(meansin, meancos), where meansin =
%   nanmean(sind(theta)) and meancos = nanmean(cosd(theta)). See
%   https://rosettacode.org/wiki/Averages/Mean_angle

E = JLLErrors;


if ~isnumeric(theta) 
    E.badinput('THETA must be numeric');
end

if ~exist('dim','var')
    if isvector(theta)
        dim = find(size(theta) > 1);
    else 
        dim = 1;
    end
end

meansin = nanmean(sind(theta),dim);
meancos = nanmean(cosd(theta),dim);

mean_theta = atan2d(meansin, meancos);

end

