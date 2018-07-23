function varargout = combine_plots( varargin )
%COMBINE_PLOTS Combine individual plots into a single figure
%   COMBINE_PLOTS() Combined all open figures into a single figure
%   keeping the number of axes along each dimension as close to
%   equal as possible.
%
%   COMBINE_PLOTS(FIGS) Combine only the figures specified by the
%   array of handles FIGS.
%
%   COMBINE_PLOTS( ___, 'dims', [y x] ) Combined with either previous
%   syntax, specify the dimensions as vertical x horizontal (as in 
%   subplot). If one of y or x are 0, it will be automatically calculated.
%   At most one may be zero.
%
%   COMBINE_PLOTS( ___, 'scale', 1 ) will automatically scale the combined
%   figure to try to keep all the sub plots at their original size. You may
%   change the factor of 1 to scale them to that fraction of their original
%   size. If this parameter is omitted, the combined plot is not resized at
%   all from the default figure size. This can also be a 1x2 vector if you
%   want to scale the vertical and horizontal directions differently. (The
%   vector would be [horizontal vertical].)

E = JLLErrors;
p = inputParser;

p.addOptional('figs',[]);
p.addParameter('dims',[]);
p.addParameter('scale',-1)
p.parse(varargin{:});
pout = p.Results;
figs = pout.figs;
subplot_dims = pout.dims;
scale_figure = pout.scale;

if isempty(figs)
    figs = get(0,'children');
elseif any(~isgraphics(figs,'figure'))
    E.badinput('FIGS must contain only figure handles');
end

% Calculate how many axes we'll need in the subplot
nfigs = 0;
figinds = [];
n_obj_types = 3; % how many object types we'll try to copy. Right now 3: axes, legend, and colorbar
inds = nan(1,n_obj_types+1);
for a=1:numel(figs)
    nfigs = nfigs + sum(isgraphics(figs(a).Children, 'axes'));
    inds(1) = a;
    % Go through the list of children for this figure. Assume that any
    % associated objects (e.g. colorbar) are listed before the associated
    % axis, so find each type of object we want to copy, then when we find
    % an axis instance, reset.
    for b=1:numel(figs(a).Children)
        if isgraphics(figs(a).Children(b), 'colorbar')
            inds(3) = b;
        elseif isgraphics(figs(a).Children(b), 'legend')
            inds(4) = b;
        elseif isgraphics(figs(a).Children(b), 'axes')
            inds(2) = b;
            figinds = cat(1, figinds, inds);
            inds = nan(1,n_obj_types+1);
            inds(1) = a;
        end
    end
end

% Make the subplots roughly even in each dimension, preferring width over
% height
if isempty(subplot_dims)
    x = ceil(sqrt(nfigs));
    y = ceil(nfigs/x);
else
    if ~isnumeric(subplot_dims) || numel(subplot_dims) > 2 || any(subplot_dims < 0)  || any(mod(subplot_dims,1) > 0)
        E.badinput('''dims'' must be a two-element positive vector of whole numbers')
    elseif prod(subplot_dims) < nfigs && ~any(subplot_dims == 0)
        E.badinput('''dims'' contains insufficient plots (%d axes to combine found)',nfigs)
    end
    
    xx_zero = subplot_dims == 0;
    if sum(xx_zero) == 1
        % Calculate how many plots in the flexible dimension we'll need
        subplot_dims(xx_zero) = ceil(nfigs/subplot_dims(~xx_zero));
    elseif sum(xx_zero) > 1
        E.badinput('''dims'' may have at most one 0 in it');
    end
    
    y = subplot_dims(1);
    x = subplot_dims(2);
end

spfig = figure;
for a=1:nfigs
    sp = subplot(y,x,a);
    sp.Units = 'Normalized';
    pos = sp.Position;
    delete(sp);
    
    f = figinds(a,1); 
    cc = ~isnan(figinds(a,:)); cc(1) = false; % never copy the whole figure (which is the handle index in the first column)
    c = figinds(a,cc);
    copy_ax = copyobj(figs(f).Children(c), spfig);
    new_ax = isgraphics(copy_ax, 'axes');
    copy_ax(new_ax).Position = pos;
    
    % Find the axis child and make sure it's x and y limits are unchanged.
    old_xlim = figs(f).Children(figinds(a,2)).XLim;
    old_ylim = figs(f).Children(figinds(a,2)).YLim;
    copy_ax(new_ax).XLim = old_xlim;
    copy_ax(new_ax).YLim = old_ylim;
end

if numel(scale_figure) > 1 || scale_figure > 0
    % Of course the order in the position vector is different. Position(3)
    % is the width of the figure and Position(4) is its height.
    spfig.Position(3:4) =  scale_figure .* [x y] .* spfig.Position(3:4);
end

if nargout > 0
    varargout{1} = spfig;
end


end

