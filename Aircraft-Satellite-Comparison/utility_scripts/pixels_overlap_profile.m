function xx = pixels_overlap_profile(pixel_loncorn, pixel_latcorn, prof_lon, prof_lat, varargin)
% PIXELS_OVERLAP_PROFILE Return a logical array for which pixels overlap an aircraft profile
%
%   When verifying satellite columns by integrating an aircraft profile, we
%   need to know which pixels overlap the profile.
%
%   XX = PIXELS_OVERLAP_PROFILE(PIXEL_LONCORN, PIXEL_LATCORN, PROF_LON,
%   PROF_LAT) returns a logical array XX that has size [size(PIXEL_LONCORN,
%   2), size(PIXEL_LONCORN,3)] that is true for pixels that at least one
%   point falls in. PIXEL_LONCORN and PIXEL_LATCORN must be arrays that
%   define the pixel longitude and latitude corners, with the first
%   dimension having length 4 so that PIXEL_LONCORN(:,i) are the corner
%   longitudes for pixel i. PROF_LON and PROF_LON are numeric arrays that
%   give the longitude and latitude of each measurement.
%
%   XX = PIXELS_OVERLAP_PROFILE( ___, PROF_VALID ) PROF_VALID is an
%   optional array that must be the same size as PROF_LON and PROF_LAT that
%   is true for all points in the profile that are a valid measurement. If
%   a satellite pixel only overlaps pixels for which PROF_VALID is false,
%   it will not be included.
%
%   Additional parameters:
%
%       'min_frac_in_pixel' - sets the minimum fraction of valid points in
%       the profile that must fall within a pixel for the pixel to be
%       included. Default is 0, i.e. this criterion is disabled. If set to
%       a e.g. 0.25, will require that 25% of valid points in the profile
%       fall within a pixel for that pixels to be matched with the profile.
%
%       'min_points_in_pixel' - sets the minimum number of valid points in
%       the profile that must fall within a pixel for the pixel to be
%       included. Default is 1, i.e. any pixel with a valid profile
%       measurment in it is include. Both this and the minimum fraction
%       criterion must be satisfied for a pixel to be included, so if this
%       is set to 10 and min_frac is set to 0.2, then a pixel must contain
%       at least 10 valid profile measurements and those measurements must
%       constitude at least 20% of the valid measurements in the profile.
E = JLLErrors;

p = inputParser;
empty_logical = false(1,1);
empty_logical = empty_logical([]);
p.addOptional('prof_valid',empty_logical,@islogical);
p.addParameter('min_frac_in_pixel',0);
p.addParameter('min_points_in_pixel',1);
p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

prof_valid = pout.prof_valid;
min_frac = pout.min_frac_in_pixel;
min_points = pout.min_points_in_pixel;

if isempty(prof_valid)
    prof_valid = true(size(prof_lon));
end

if ~isnumeric(pixel_loncorn)
    E.badinput('PIXEL_LONCORN must be numeric')
end

if ~isnumeric(pixel_latcorn)
    E.badinput('PIXEL_LATCORN must be numeric')
end

if ~isnumeric(prof_lon)
    E.badinput('PROF_LON must be numeric')
end

if ~isnumeric(prof_lat)
    E.badinput('PROF_LAT must be numeric')
end

if ~islogical(prof_valid)
    E.badinput('PROF_VALID must be logical')
end

if ~isequal(size(prof_lon), size(prof_lat)) || ~isequal(size(prof_lon), size(prof_valid))
    E.badinput('PROF_LON, PROF_LAT, and (if given) PROF_VALID must all be the same size')
end

if ~isequal(size(pixel_loncorn), size(pixel_latcorn))
    E.badinput('PIXEL_LONCORN and PIXEL_LATCORN must be the same size')
elseif size(pixel_loncorn,1) ~= 4
    E.badinput('PIXEL_LONCORN and PIXEL_LATCORN must have a length of 4 in the first dimension')
end

if ~isnumeric(min_frac) || ~isscalar(min_frac) || min_frac < 0
    E.badinput('''min_frac'' must be a scalar, positive (or zero) number')
end

if ~isnumeric(min_points) || ~isscalar(min_points) || min_points < 1
    E.badinput('''min_points must be a scalar number >= 1');
end

xx = false(size(pixel_loncorn,2), size(pixel_loncorn,3));
% As far as I know, there's no way around the fact that we need to check
% each pixel individually in a loop, which would be slow if we went over
% every single pixel with inpolygon (the problem is we need to know which
% polygons contain the profile, rather than just what part of the profile
% lies within each polygon).
for i_pix = 1:numel(xx)
    % Will need to time this. Is it slow to call inpolygon for each pixel?
    % Is it faster to remove pixels definitely outside the profile?
    prof_in = inpolygon(prof_lon, prof_lat, pixel_loncorn(:,i_pix), pixel_latcorn(:,i_pix));
    
    % We have two criteria for whether a pixel overlaps the profile. First
    % is just the minimum number of points from the profile in the pixel.
    % Second is the fraction of the profile that lies within the pixel.
    % This way we can require that the pixel is reasonably well sampled by
    % the profile. We check both against only points in the profile with
    % valid NO2 measurements.
    frac_in = sum(prof_valid(prof_in))/sum(prof_valid);
    xx(i_pix) = sum(prof_valid(prof_in)) >= min_points && frac_in >= min_frac;
end
end