function [ scaled_column_densities ] = scale_profiles_by_time( prof_times, column_densities, campaign_name, varargin )
%scale_profiles_by_time Scale aircraft column densities to OMI overpass
%   Sally Pusede showed that, at least during the Balitmore campaign, there
%   is a definite change in NO2 profile shape with time-of-day.  Normally
%   with simple satellite-aircraft validation, this would not be a
%   significant problem because I would restrict the aircraft columns to be
%   used to +/- 1.5 hours on either side of OMI's overpass time, which
%   minimizes this TOD variation.  However, when considering other factors
%   (i.e. aerosol effects) it is useful to be able to use a broader
%   temporal range to increase the statistics available.
%
%   This function requires three inputs:
%       1) a vector of aircraft profile start times, in hours after
%       midnight local standard time.
%       2) a vector of column densities derived from aircraft measurements.
%       This must be the same size as the profile start times vector.
%       3) the campaign name (discover-md, discover-ca, discover-tx, or
%       discover-co). This will be used to select the scale factors to use
%
%   An optional fourth argument can be the string 'makenans' or 'keep'.
%   This determines how the function behaves when the time to calculate a
%   scale factor for is outside the range of times it has data for.
%   'makenans' is the default if no string is passed; this will return
%   these columns as NaNs. 'keep' returns these columns without scaling.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(3,4);

if ~isvector(prof_times) || ~isvector(column_densities)
    error(E.badinput('"prof_times" and "column_densities" are expected to be vectors'));
elseif ~ischar(campaign_name)
    error(E.badinput('"campaign_name" expected to be a string'));
end

if ~isnumeric(prof_times)
    error(E.badinput('"prof_times" expected to be a numeric vector'));
elseif any(prof_times > 24)
    error(E.badinput('"prof_times" should be given as hours after midnight local standard time; it looks like you passed seconds after mignight instead'));
end

if numel(prof_times) ~= numel(column_densities)
    error(E.badinput('"prof_times" and "column_densities" must have the same number of elements'));
end

if numel(varargin) > 0
    if ~any(strcmpi(varargin{1},{'keep','makenans'}))
        error(E.badinput('The optional 4th argument must be either the string ''keep'' or ''makenans'''));
    end
end

% Make both start times and column densities column vectors to save
% headaches
if ~iscolumn(prof_times); prof_times = prof_times'; end
if ~iscolumn(column_densities); column_densities = column_densities'; end

% Decide whether to leave scale factors outside the range of times
% available as NaNs (which will make the resulting scaled columns NaNs) or
% set them to 1, which will leave the columns unchanged.
if numel(varargin) == 0 || strcmpi(varargin{1},'makenans')
    set_to_one_bool = false;
elseif strcmpi(varargin{1},'keep')
    set_to_one_bool = true;
else
    error(E.unknownError('Unexpected case for optional arguments - something broke'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Simplify the campaign name to allow simpler comparisons
campaign_name = lower(regexprep(campaign_name,'[\W_]',''));

% This is the section to store the scale factors.  These are calculated by
% binning profile columns over the whole campaign by their start times
% (each bin being an hour wide), taking the median column density for the
% bin, and dividing all of those values by the value for the bin centered
% on 13:30 local standard time - i.e. the OMI overpass time.
switch campaign_name
    case 'discovermd'
        bintimes = [7.5000 9.5000 10.5000 11.5000 12.5000 13.5000 14.5000 15.5000 16.5000];
        binscale = [1.1648 0.9422 1.0607 1.0771 0.8195 1.0000 1.1866 1.1956 1.4429];
    case 'discoverca'
        bintimes = [8.5000 9.5000 10.5000 11.5000 12.5000 13.5000 14.5000];
        binscale = [1.1200 0.6503 0.5692 0.7499 3.0034 1.0000 1.0854];
    otherwise
        error(E.callError('unknown_campaign',sprintf('Could not identifty campaign "%s" (non-letter characters removed)',campaign_name)));
end

% Intepolate the scale factor to the start time of each column passed,
% these will be used individually to scale to columns to the OMI overpass
% time - this will attempt to correct for the time-of-day variation in
% column density.
scale_factors = interp1(bintimes,binscale,prof_times);

if set_to_one_bool
    scale_factors(isnan(scale_factors)) = 1;
end


% Dividing the measured columns by the scale factor scales them to OMI
% overpass (approximately)

scaled_column_densities = column_densities ./ scale_factors;
        

end

