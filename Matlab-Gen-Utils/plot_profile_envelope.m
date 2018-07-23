function [prof, prof_lb, prof_ub, l] = plot_profile_envelope(prof_x_ens, prof_y, envelope_mode, varargin)
%PLOT_PROFILE_ENVELOPE Plot the uncertainty in a profile as an patch envelope
%   [ PROF, PROF_LB, PROF_UB, L ] = PLOT_PROFILE_ENVELOPE( PROF_X_ENS,
%   PROF_Y, ENVELOPE_MODE ) Given an ensemble of profiles (PROF_X_ENS)
%   defined at y-coordinates PROF_Y, calculate the uncertainty given by the
%   variability of the profiles and plot it as a shaded patch object.
%   ENVELOPE_MODE can be:
%
%       'std' - envelope calculated as a multiple of the standard deviation
%       in profile values at each model level (default is 1 std.
%       deviation). The mean profile is plotted.
%
%       'quantile' - envelope calculated as an lower and upper quantile of
%       the variability (default is 0.25 and 0.75, respectively). The
%       median profile is plotted.
%
%       'range' - envelope calculated as the min and max values of the
%       ensemble at each level. The mean profile is plotted.
%
%   [ ___ ] = PLOT_PROFILE_ENVELOPE( ___, ENVELOPE_WIDTH ) ENVELOPE_WIDTH
%   affects how width the envelope should be. For ENVELOPE_MODE == 'std',
%   the width is the multiple of sigma that the envelope should be (default
%   1). For ENVELOPE_MODE == 'quantile', this must be a two element vector
%   giving the lower and upper quantile (in that order), default is [0.25
%   0.75]. This is not used for ENVELOPE_MODE == 'range'.
%
%   Additional parameters:
%
%       'axes' - which axis to plot in. Default is current.
%
%       'linecolor', 'linewidth', 'linestyle' - control the appearance of
%       the mean/median profile. Passed directly to the line() function.
%
%       'envalpha', 'envcolor' - control the appearance of the envelope.
%       'envalpha' controls the transparency (0 is fully transparent, 1 is
%       fully opaque, default is 1). 'envcolor' accepts any standard color
%       specification. Default is a light gray ([0.8 0.8 0.8]).

E = JLLErrors;

p = advInputParser;
p.addOptional('envelope_width', nan);

% Plot params
p.addParameter('axes', gca);

% Line params
p.addParameter('linecolor', 'k');
p.addParameter('linewidth', 2);
p.addParameter('linestyle', '-');

% Envelope params
p.addParameter('envalpha', 1);
p.addParameter('envcolor', [0.8 0.8 0.8]);
p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

allowed_env_modes = {'std', 'quantile', 'range'};

if size(prof_x_ens,1) ~= size(prof_y,1)
    E.badinput('PROF_X_ENS and PROF_Y must have the same length in the first dimension')
elseif size(prof_y,2) ~= 1
    E.badinput('PROF_Y must be a column vector')
elseif ~ismember(envelope_mode, allowed_env_modes)
    E.badinput('ENVELOPE_MODE must be one of: %s', strjoin(allowed_env_modes, ', '));
end

ax = pout.axes;

env_color = pout.envcolor;
env_alpha = pout.envalpha;

line_color = pout.linecolor;
line_width = pout.linewidth;
line_style = pout.linestyle;

env_width = pout.envelope_width;

if isnan(env_width)
    switch envelope_mode
        case 'std'
            env_width = 1;
        case 'quantile'
            env_width = [0.25 0.75];
        case 'range'
            % not used
            env_width = [];
        otherwise
            E.notimplemented('No default width defined for envelope type "%s"', envelope_mode);
    end
end

switch envelope_mode
    case 'std'
        if ~isscalar(env_width) || ~isnumeric(env_width) || env_width < 0
            E.badinput('ENVELOPE_WIDTH must be a scalar, positive number for ENVELOPE_MODE = "std"');
        end
        env_fxn = @calc_std_dev;
    case 'quantile'
        if numel(env_width) ~= 2 || ~isnumeric(env_width) || any(env_width < 0 | env_width > 1)
            E.badinput('ENVELOPE_WIDTH must be a two element vector with both elements between 0 and 1 for ENVELOPE_MODE = "quantile"');
        end
        env_fxn = @calc_quantiles;
    case 'range'
        env_fxn = @calc_range;
    otherwise
        E.notimplemented('No input checking implemented for ENVELOPE_MODE = "%s"', envelope_mode);
end


%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

% Calculate the upper and lower error bounds
[prof, prof_lb, prof_ub] = env_fxn(prof_x_ens, env_width);
plot_error_envelope_x(prof_y, prof_lb, prof_ub, 'colorspec', env_color, 'facealpha', env_alpha, 'ax', ax);
l = line(prof, prof_y, 'color', line_color, 'linestyle', line_style, 'linewidth', line_width);

end

function [y, lower_y, upper_y] = calc_std_dev(prof_ens, n_sigma)
y = nanmean(prof_ens, 2);
sigma_y = nanstd(prof_ens, 0, 2) * n_sigma;
lower_y = y - sigma_y;
upper_y = y + sigma_y;
end

function [y, lower_y, upper_y] = calc_quantiles(prof_ens, quantiles)
y = nanmedian(prof_ens, 2);
q = quantile(prof_ens, quantiles, 2);
lower_y = q(:,1);
upper_y = q(:,2);
end

function [y, lower_y, upper_y] = calc_range(prof_ens, ~)
y = nanmean(prof_ens, 2);
lower_y = nanmin(prof_ens, [], 2);
upper_y = nanmax(prof_ens, [], 2);
end
