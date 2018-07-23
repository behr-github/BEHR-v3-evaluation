function [ varargout ] = plot_vertical_profile_bins( data_vals, altitude, errors, varargin )
% plot_vertical_profile_bins(data_vals, altitude, error, [opt] data_vals2,
% altitude2, errors2): Plots vertical profiles using data already binned
% with one of the existing binning functions.
%   This function will plot vertical profiles for data passed to to.
%   "data_vals" must contain the values of the data that correspond to the
%   altitudes given in "altitude". "errors" will be a 2 x n or 1 x n matrix
%   with either 25th/75th quantiles or std. deviations.  These are all
%   returned by functions such as bin_vertical_profile and
%   avg_vertical_profiles.
%
%   Optionally, you can pass a second set of values, altitudes, and errors
%   to plot two vertical profiles.  In this case, the function
%   will return axis handles.
%
%   To make a legend, pass a cell array with identities of your first and
%   second data as strings using the parameter value 'legend'.
%   Other parameter values are 'xlabel','ylabel', and 'title', which are
%   useful if you wish to automate this function.

p = inputParser;
p.addRequired('vals',@ismatrix);
p.addRequired('altitude',@ismatrix);
p.addRequired('errors',@ismatrix);
p.addOptional('vals2',[],@ismatrix);
p.addOptional('alt2',[],@ismatrix);
p.addOptional('err2',[],@ismatrix);
p.addParamValue('legend',{},@iscell);
p.addParamValue('xlabel','',@isstr);
p.addParamValue('ylabel','',@isstr);
p.addParamValue('title','',@isstr);

p.parse(data_vals, altitude, errors, varargin{:});
pout = p.Results;

data_vals = pout.vals; if ~iscolumn(data_vals); data_vals = data_vals'; end
altitude = pout.altitude; if ~iscolumn(altitude); altitude = altitude'; end
errors = pout.errors; if size(errors,2)>2; errors = errors'; end

data_vals2 = pout.vals2; if ~iscolumn(data_vals2); data_vals2 = data_vals2'; end
altitude2 = pout.alt2; if ~iscolumn(altitude2); altitude2 = altitude2'; end
errors2 = pout.err2; if size(errors2,2)>2; errors2 = errors2'; end

if any([~isempty(data_vals2), ~isempty(altitude2), ~isempty(errors2)]) && ~all([~isempty(data_vals2), ~isempty(altitude2), ~isempty(errors2)])
    error('plot_vert_bins:second_field','Second profile variables must either be all passed or none passed')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    PLOT PROFILE     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plot the error envelopes first so that they don't cover the lines
% themselves
if size(errors,2) == 2;
    fh = plot_error_envelope_x(altitude',errors(:,1)',errors(:,2)','colorspec',[0.7 0.7 1]);
elseif size(1,errors) == 1;
    fh = plot_error_envelope_x(altitude,data_vals - errors, data_vals + errors, 'colorspec', [0.7 0.7 1]);
end

hold on

    hax1 = gca;
if ~isempty(data_vals2)
    ax_pos = get(hax1,'Position');
    old_ylim = get(hax1, 'YLim');
    set(hax1,'box','off') %Turn of the box so that tick marks don't appear on both sides.
    hax2 = axes('Position', ax_pos,'XAxisLocation','top','YLim',old_ylim,'YTick',[],'Color','none');
    
    if size(errors2,2) == 2;
        plot_error_envelope_x(altitude2,errors2(:,1),errors2(:,2),'colorspec',[1 0.7 0.7], 'FaceAlpha', 0.5, 'fignum', fh);
    elseif size(1,errors2) == 1;
        plot_error_envelope_x(altitude2,data_vals2 - errors2, data_vals2 + errors2, 'colorspec', [1 0.7 0.7], 'FaceAlpha',0.5, 'fignum', fh);
    end
    
end

hold on % error_envelope turns hold off

lh1 = line(data_vals, altitude, 'linestyle','-','color','b','linewidth',2,'Parent',hax1);

if ~isempty(data_vals2)
    lh2 = line(data_vals2, altitude2, 'linestyle','--','color','r','linewidth',2,'Parent',hax2);
    varargout = {[hax1, hax2], lh1, lh2};
end

if ~isempty(pout.legend)
    legend([lh1; lh2],pout.legend);
end

if ~isempty(pout.xlabel)
    xlabel(pout.xlabel,'fontsize',16,'Parent',hax1)
end

if ~isempty(pout.ylabel);
    ylabel(pout.ylabel,'fontsize',16,'Parent',hax1)
end

if ~isempty(pout.title)
    title(pout.title,'fontsize',18)
end

end

