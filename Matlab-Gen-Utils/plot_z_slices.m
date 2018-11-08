function plot_z_slices(x, y, z, varargin)
%PLOT_Z_SLICES Plot 3D data in a series of 2D slices
%
%   PLOT_Z_SLICES( X, Y, Z ) Given data Z defined on the coordinates X and Y,
%   plot 2D scatter plots of X and Y dividing the data into 9 bins based on the
%   Z values.
%
%   PLOT_Z_SLICES( X, Y, Z, C ) Similar to the previous syntax, except
%   additionally color the points by the data in C.
%
%   Additional parameters:
%
%       'plot' - which plot type to make. Options are 'scatter2' (default),
%       'surf', 'surfc', 'contourf', and 'pcolor', each of which make the plot
%       with the respective function. All but the scatter plots rely on the X
%       and Y binning, see next parameter.
%
%       'xbins' and 'ybins' - how many bins in the X and Y coordinates to make,
%       default is 10 for both. Used in the surface, contour, and pcolor plots
%       to average data along the X and Y coordinates to get it into a grid.
%
%       'zbins' - how many bins to use for the Z data; i.e. how many plots to
%       make.
%
%       'zname' - the name to use for the Z variable, which will be used in the
%       title for each figure to indicate the range of Z values. Default is 'z'.
%
%       'onefig' - controls whether to make subplots in one figure (true,
%       default) or multiple separate figures (false).
%
%       'ax_opts' - a cell array of options to set on each axis. Will be set by
%       calling SET(ax, ax_opts{:}). There are two special options: 'xlabelstr'
%       and 'ylabelstr' which allow you to set the X and Y axis labels directly
%       by passing the desired string, since 'xlabel' and 'ylabel' require a
%       graphics object.
%
%       'plot_opts' - additional options to pass to the plotting function, again
%       as a cell array: for example, if plotting a scatter plot, then it is
%       called as SCATTER(ax, X, Y, plot_opts{:}).
%
%       'cblabel' - a string to use as the labels for the colorbars on each
%       plot.
%
%       'clim' - color limits to set for each plot, as a two-element numeric
%       vector.

E = JLLErrors;
p = advInputParser;
p.addOptional('c', []);
p.addParameter('plot', 'scatter2');
p.addParameter('xbins', 10);
p.addParameter('ybins', 10);
p.addParameter('zbins', 9);
p.addParameter('zname','z');
p.addParameter('onefig', true);
p.addParameter('ax_opts', {});
p.addParameter('plot_opts', {});
p.addParameter('cblabel', '');
p.addParameter('clim',[]);

p.parse(varargin{:});
pout = p.Results;

c = pout.c;
plot_type = pout.plot;

nd_bins = {pout.xbins, pout.ybins};
zbins = setup_zbins(pout.zbins, z);
z_varname = pout.zname;

on_one_figure = pout.onefig;
ax_opts = pout.ax_opts;
plot_opts = pout.plot_opts;
cb_label = pout.cblabel;
c_limits = pout.clim;

x_limits = [min(x(:)), max(x(:))];
y_limits = [min(y(:)), max(y(:))];
if isempty(c_limits)
    c_limits = [min(c(:)), max(c(:))];
end

num_bins = size(zbins,1);

if on_one_figure
    % need to calculate the subplot x and y dims; try to keep as square as
    % possible
    subp_y = floor(sqrt(num_bins));
    subp_x = ceil(num_bins / subp_y);
    figs = figure;
else
    figs = gobjects(num_bins,1);
end


for i_z = 1:num_bins
    if on_one_figure
        curr_ax = subplot(subp_y, subp_x, i_z);
    else
        figs(i_z) = figure;
        curr_ax = gca;
    end
        
    
    z_idx = z >= zbins(i_z,1) & z < zbins(i_z,2);
    make_plot(curr_ax, z_idx);
end



    function make_plot(ax, z_idx)
        E = JLLErrors;
        switch lower(plot_type)
            case 'scatter2'
                make_scatter2(ax, z_idx);
            case 'surf'
                make_surf(ax, z_idx, @surf);
            case 'surfc'
                make_surf(ax, z_idx, @surfc);
            case 'contourf'
                make_surf(ax, z_idx, @contourf);
            case 'pcolor'
                make_surf(ax, z_idx, @pcolor);
                shading flat
            otherwise
                E.badinput('Plot type "%s" not recognized', plot_type)
        end
        xlim(x_limits);
        ylim(y_limits);
        title(sprintf('%s \\in = [%.4g, %.4g]', z_varname, nanmin(z(z_idx)), nanmax(z(z_idx))));
    end

    function make_scatter2(ax, z_idx)
        if ~isempty(c)
            scatter(ax, x(z_idx), y(z_idx), [], c(z_idx), plot_opts{:});
            add_colorbar(ax);
        else
            scatter(ax, x(z_idx), y(z_idx), plot_opts{:});
        end
        set_ax_opts(ax);
    end

    function make_scatter3(ax, z_idx)
        %TODO
    end

    function make_surf(ax, z_idx, surf_fxn)
        [bin_vals, bin_centers] = nd_averaging(c(z_idx), {x(z_idx), y(z_idx)}, nd_bins);
        [Y,X] = meshgrid(bin_centers{2}, bin_centers{1});
        surf_fxn(ax, X, Y, bin_vals, plot_opts{:});
        add_colorbar(ax);
        set_ax_opts(ax);
    end

    function add_colorbar(ax)
        cb = colorbar('peer', ax);
        cb.Label.String = cb_label;
        caxis(ax, c_limits);
    end

    function set_ax_opts(ax)
        %TODO: handle xlabel, ylabel, etc.
        if ~isempty(ax_opts)
            % Can't do set(ax,'xlabel',string) b/c it expects a graphics
            % object. So we treat the labels specially. Take advantage of
            % the fact that ax_opts must alternate option name/value and
            % turn it into a structure to make it easier to get specific
            % options.
            opts = struct(ax_opts{:});
            if isfield(opts, 'xlabelstr')
                xlabel(opts.xlabelstr)
            end
            if isfield(opts,'ylabelstr')
                ylabel(opts.ylabelstr);
            end
            
            tmp_opts = update_params('remove', ax_opts, 'xlabelstr', 'ylabelstr');
            set(ax, tmp_opts{:});
        end
    end
end

function zbins = setup_zbins(zbins_in, z)
if ~isnumeric(zbins_in)
    E.badinput('"zbins" must be numeric')
end

if isscalar(zbins_in)
    if zbins_in < 1
        E.badinput('If given as a scalar, "zbins" must be >= 1')
    end
    zbins_tmp = linspace(min(z), max(z)+eps(max(z)), zbins_in+1)';
    zbins = [zbins_tmp(1:end-1), zbins_tmp(2:end)];
elseif isvector(zbins_in)
    if any(diff(zbins_in) < 0)
        E.badinput('If given as a vector, "zbins" must be monotonically increasing')
    end
    zbins = [reshape(zbins_in(1:end-1),[],1), reshape(zbins_in(2:end),[],1)];
elseif ismatrix(zbins_in) && size(zbins_in,2) == 2
    zbins = zbins_in;
else
    E.badinput('"zbins" must be a scalar, vector, or n-by-2 matrix')
end

end
