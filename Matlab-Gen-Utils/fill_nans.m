function [ varargout ] = fill_nans( x,y,varargin )
%fill_nans(x,y): Fills nans in 'y' unless they are at the end of the
%series, then they are removed
%   Linearly interpolates NaNs in 'y' that have non-NaN values around them.
%    NaNs at the beginning or end of 'y' are removed along with their
%    corresponding 'x' value.
%
%   This function can also interpolate error measurements assuming the
%   linear interpolation formula:
%
%       y = y0 + (y1 - y0) * (x - x0)/(x1 - x0)
%
%   Following the method of error propagation where:
%   
%       e_y = sqrt( sum_i (dy / dvi)^2 (e_vi)^2 )
%
%   e represents the error term and vi each variable with an error
%   component.  This function assumes that there is no error in the x
%   variables, so that the error propagation is:
%
%       e_y = sqrt( [1 - (x - x0)/(x1 - x0)]^2 * (e_y0)^2 + [(x-x0)/(x1-x0)]^2 * (e_y1)^2] )
%
%   Pass as a third argument a vector of equal length to x & y that
%   contains the known error terms.  This is optional.
%
%   If you wish to override the default behavior so that this function only
%   interpolates interior nans and leaves leading or trailing ones present,
%   pass 'noclip' as the optional final parameter
%
%   Josh Laughner <joshlaugh5@gmail.com> Updated 5 Feb 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION AND PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,4);
% First handle if there are four inputs
if nargin == 4
    % The first optional input should be the error vector
    if isnumeric(varargin{1})
        error_vec = varargin{1};
        error_bool = true;
    elseif ischar(varargin{1})
        error(E.badinput('''noclip'' must always be the final parameter passed'));
    else
        error(E.badinput('Option parameter not recognized - should be a numeric vector'));
    end
    % The second should be the string 'noclip.' Check that it is a string
    % and it is noclip, otherwise throw an error.
    if ischar(varargin{2})
        if strcmpi(varargin{2},'noclip')
            clipping_bool = false;
        else
            error(E.badinput(sprintf('The value %s for the final parameter is not understood. Should be ''noclip'' or nothing',varargin{2})));
        end
    else
        error(E.badinput('The final parameter if four arguments are passed should be a string.'));
    end
elseif nargin == 3
    % If there are 3 input, figure out if the optional one is a vector of
    % errors or 'noclip'. If it is neither, throw an error.
    if ischar(varargin{1})
        if strcmpi(varargin{1},'noclip')
            clipping_bool = false;
        else
            error(E.badinput('The only recognized string is ''noclip'''));
        end
    elseif isnumeric(varargin{1})
        error_vec = varargin{1};
        error_bool = true;
    else
        error(E.badinput('Option parameter not recognized - should be a numeric vector'));
    end
end

% Check if the booleans for error interpolation and clipping were defined.
% If not, set them to their default values (false for error - no error
% vector passed - and true for clipping).
if ~exist('clipping_bool','var')
    clipping_bool = true;
end
if ~exist('error_bool','var')
    error_bool = false;
end

% If an error vector was passed but only 2 outputs assigned, then warn the
% user
if error_bool && nargout < 3
    warning('An error vector has been passed but only 2 outputs assigned. The interpolated error vector cannot be returned.')
end

% Finally check that x & y are not all nans, not empty, and the same
% length.  Do the same for the error vector if it exists

if isempty(x)
    error(E.badinput('''x'' must not be empty'));
elseif all(isnan(x))
    error(E.badinput('''x'' must not be all nans'));
elseif isempty(y)
    error(E.badinput('''y'' must not be empty'));
elseif all(isnan(y))
    error(E.badinput('''y'' must not be all nans.'));
elseif ndims(x) ~= ndims(y) || ~all(size(x) == size(y))
    error(E.badinput('''x'' and ''y'' must have the same dimensions'));
end

if error_bool
    if isempty(error_vec)
        error(E.badinput('''error_vec'' must not be empty'));
    elseif all(isnan(error_vec))
        error(E.badinput('''x'' must not be all nans'));
    elseif ndims(error_vec) ~= ndims(x) || ~all(size(error_vec) == size(x))
        error(E.badinput('''error_vec'' must have the same dimensions as ''x'' and ''y'''));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Trim leading or trailing NaNs

if clipping_bool
    if isnan(y(1))
        first_val = find(~isnan(y),1,'first');
        y = y(first_val:end);
        x = x(first_val:end);
        if error_bool
            error_vec = error_vec(first_val:end);
        end
    end
    if isnan(y(end))
        last_val = find(~isnan(y),1,'last');
        y = y(1:last_val);
        x = x(1:last_val);
        if error_bool
            error_vec = error_vec(last_val:end);
        end
    end
end

% Interpolate to fill in the nans
nans = isnan(y);
y(nans) = interp1(x(~nans),y(~nans),x(nans));

% If handling an error vector as well, check that the locations of nans in
% it are the same as in y. If not, warn the user and find where the nans
% are in the error vector
if error_bool
    if ~all(isnan(error_vec)==nans)
        warning('NaNs in error_vec are not the same as those in y.');
        nans = isnan(error_vec);
    end
    first_val = find(~nans,1,'first');
    last_val = find(~nans,1,'last');
    
    % Now go through each "internal" NaN and calculate the error for the
    % interpolated point following the formula defined in the introduction
    for a=first_val:last_val
        if nans(a)
            % Find the previous and following original non-nan value: not an
            % interpolated one.
            prev = a;
            while true
                prev = prev-1;
                if ~nans(prev)
                    break
                elseif prev < 1
                    error(E.callError('index_not_found','Could not find a previous, original, non-NaN value'));
                end
            end
            
            next = a;
            while true
                next = next+1;
                if ~nans(next)
                    break
                elseif next > numel(nans)
                    error(E.callError('index_not_found','Could not find a following, original, non-NaN value'));
                end
            end
            
            x0 = x(prev);
            x1 = x(next);
            e_y0 = error_vec(prev);
            e_y1 = error_vec(next);
            
            error_vec(a) = sqrt( (1 - (x(a) - x0)/(x1 - x0))^2 * e_y0^2 + ((x(a) - x0)/(x1-x0))^2 * e_y1^2 );
        end
    end
end
% Set the output arguments
varargout{1} = x;
varargout{2} = y;
if error_bool
    varargout{3} = error_vec;
end

end