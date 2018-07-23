function [ vals ] = scale_to_range( vals, varargin )
%SCALE_TO_RANGE Linearly scales numeric data to a given range of values
%   VALS = SCALE_TO_RANGE( VALS, [MINVAL, MAXVAL] )
%   VALS = SCALE_TO_RANGE( VALS, MINVAL, MAXVAL ) In either form, VALS is
%   remapped to the range [MINVAL, MAXVAL] so that the old minimum value of
%   VALS is mapped to MINVAL and likewise for the old maximum value.

if nargin == 2
    if ~isnumeric(varargin{1}) || numel(varargin{1}) ~= 2
        E.badinput('In the two input form, RANGE must be a 2-element numeric vector');
    end
    minval = varargin{1}(1);
    maxval = varargin{1}(2);
elseif nargin == 3
    minval = varargin{1};
    maxval = varargin{2};
    if ~isnumeric(minval) || ~isscalar(minval)
        E.badinput('In the three input form, MINVAL must be a scalar number')
    end
    if ~isnumeric(maxval) || ~isscalar(maxval)
        E.badinput('In the three input form, MAXVAL must be a scalar number')
    end
end

vals = vals - min(vals(:));
vals = vals .* (maxval-minval)/max(vals(:));
vals = vals + minval;

end

