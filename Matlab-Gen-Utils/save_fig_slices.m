function [  ] = save_fig_slices( M, varargin )
%SAVE_FIG_SLICES Saves each 2D slice of the matrix M as an image
%   SAVE_FIG_SLICES( M ) will take a matrix M which (currently) must have 3
%   non-singleton dimensions. It will iterate through each 2D slice of
%   squeeze(M) and save it as an image file, by default in the current
%   folder.
%
%   There are a number of parameters which control the behaviour of this
%   function:
%       'x' and 'y' - provide x and y coordinate data that corresponds to
%       M. The arrays given for x and y can either be 2D and the same size
%       as the first two dimensions of M, in which case they will be used
%       for each slice of M, or the same size as M (specifically,
%       squeeze(x), squeeze(y), and squeeze(M) must all be the same size).
%
%       'linex' and 'liney' - x and y data for lines to be drawn on each
%       plot. They must be vectors of the same length. 
%
%       'linespec' - a cell array of name-value pairs to pass to LINE to
%       draw LINEX and LINEY. Defaults to {'color','k'}.
%
%       'title' - the title string to use for each plot. Defaults to an
%       empty string, in which case the title will just be the index of the
%       third dimension slice being plotted.
%
%       'savedir' - the directory to save the files to. Can be relative or
%       absolute. Defaults to the current directory.
%
%       'format' - the format to save the file as. Must be a valid format
%       for SAVEAS (currently), or 'fig'. Defaults to 'png'.
%
%       'ind_format' - a format string suitable for sprintf, must contain
%       exactly 1 "d", "i", or "u" format specifier, e.g. the default is 'k
%       = %d'. This is put at the end of the title to specify the index of
%       the current 2D slice, e.g. using the default will give k = 1, k =
%       2...
%
%       'caxis' - a two element vector specifying the limits of the color
%       range for all the plots. If unspecified, defaults to [min(M(:)),
%       max(M(:))].
%
%       'colormap' - the colormap to use for the plots. Must be a string or
%       64x3 array. (Note that the built-in functions jet, parula, etc. and
%       the functions blue_red_cmap, four_color_cmap, etc. included in this
%       repo return 64x3 arrays, however the built-in functions usually
%       open a figure if one isn't open, so it's better to pass them as
%       strings).
%
%       'shading' - the shading to use for the pcolor plots. Must be
%       'faceted', 'flat', or 'interp'. Defaults to 'faceted'.
%
%       'cblabel' - a string to label the colorbar with. Default is
%       nothing.
%
%       'axisprop' - a cell array of property name value pairs to set for
%       each axes. Default is nothing.
%
%       'nosave' - a boolean (default false). If true, all the figures
%       created will not be saved and will be left open in the Matlab
%       window.

p = inputParser;
p.addParameter('x',[]);
p.addParameter('y',[]);
p.addParameter('linex',[]);
p.addParameter('liney',[]);
p.addParameter('linespec',{'color','k'});
p.addParameter('title', '');
p.addParameter('savedir', pwd);
p.addParameter('format','png');
p.addParameter('ind_format', 'k = %d');
p.addParameter('caxis',[min(M(:)), max(M(:))]);
p.addParameter('colormap','parula');
p.addParameter('shading', 'faceted');
p.addParameter('cblabel','');
p.addParameter('axisprop',{});
p.addParameter('nosave',false);


p.parse(varargin{:});
pout = p.Results;

X = pout.x;
Y = pout.y;
titlestr = pout.title;
savedir = pout.savedir;
saveformat = pout.format;
ind_format_spec = pout.ind_format;
linex = pout.linex;
liney = pout.liney;
linespec = pout.linespec;
clim = pout.caxis;
cmap = pout.colormap;
pcol_shading = pout.shading;
cblabel = pout.cblabel;
axisprop = pout.axisprop;
nosave = pout.nosave;

E = JLLErrors;

M = squeeze(M);
if ndims(M) ~= 3
    E.badinput('M must have exactly 3 non-singleton dimensions');
end

% Parse x and y coordinates
if xor(isempty(X), isempty(Y))
    E.badinput('If a value for X or Y is given, a value for both must be given.')
end
if ~isempty(X) 
    X = squeeze(X);
    Y = squeeze(Y);
    if ndims(X) ~= ndims(Y)
        E.badinput('X and Y should have the same number of non-singleton dimensions');
    end
    if ismatrix(X)
        X = repmat(X, 1, 1, size(M,3));
        Y = repmat(Y, 1, 1, size(M,3));
    elseif ~isequal(size(X), size(M)) || ~isequal(size(Y), size(M))
        E.badinput('The size of the arrays for X and Y must be the same as the size of M, if X and Y are not 2D');
    end
end

% Parse line inputs
if xor(isempty(linex), isempty(liney))
    E.badinput('Both or neither LINEX and LINEY must be given')
end
if ~isempty(linex)
    if ~isvector(linex) || ~isvector(liney)
        E.badinput('LINEX and LINEY must be vectors')
    end
    linex = linex(:);
    liney = liney(:);
    if length(linex) ~= length(liney)
        E.badinput('LINEX and LINEY must be the same length')
    end
end

if ~iscell(linespec) || mod(numel(linespec),2) ~= 0
    E.badinput('LINESPEC must be a cell with an even number of elements')
end

% Parse title inputs
if ~ischar(titlestr)
    E.badinput('TITLE must be a string');
end

if ~ischar(ind_format_spec)
    E.badinput('The value for IND_FORMAT must be a string');
end
dformats = regexp(ind_format_spec, '%[-+\o{40}0]*\d*[dui]');
other_formats = regexp(ind_format_spec, '%[-+\o{40}0]*\d*[oxXfeEgGcs]','once');
if numel(dformats) ~= 1
    E.badinput('The value for IND_FORMAT must contain exactly one "d", "i", or "u" format specifier')
elseif ~isempty(other_formats)
    E.badinput('The value for IND_FORMAT cannot contain o, x, X, f, e, E, g, G, c or s format specifiers');
end

if ~ischar(cblabel)
    E.badinput('CBLABEL must be a string')
end

% Parse save inputs
if ~ischar(savedir)
    E.badinput('SAVEDIR must be a string')
elseif ~exist(savedir,'dir') 
    E.badinput('SAVEDIR (%s) is not a directory')
end

if ~ischar(saveformat)
    E.badinput('SAVEFORMAT must be a string')
end

% Parse pcolor inputs
if ~isnumeric(clim) || numel(clim) ~= 2 || clim(1) > clim(2)
    E.badinput('CLIM must be a two element numeric vector with the second element >= the first')
end

if ischar(cmap)
    % all good
elseif isnumeric(cmap) && isequal(size(cmap), [64, 3])
    % all good
else
    E.badinput('CMAP must be a string or a 64-by-3 numeric array')
end

if ~ischar(pcol_shading) || ~ismember(pcol_shading, {'faceted','flat','interp'})
    E.badinput('SHADING must be one of the strings ''faceted'', ''flat'', or ''interp''');
end

% Other inputs
if ~iscell(axisprop) || mod(numel(axisprop),2) ~= 0
    E.badinput('AXISPROP must be a cell array with an even number of elements')
end

if ~isscalar(nosave)
    E.badinput('NOSAVE must be a scalar')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

for a=1:size(M,3)
    fig=figure;
    if isempty(X)
        pcolor(M(:,:,a));
    else
        pcolor(X(:,:,a), Y(:,:,a), M(:,:,a));
    end
    
    if ~isempty(linex)
        line(linex, liney, linespec{:})
    end
    
    caxis(clim);
    colormap(cmap);
    shading(pcol_shading);
    cb = colorbar;
    cb.Label.String = cblabel;
    
    a_str = sprintf(ind_format_spec, a);
    if isempty(titlestr)
        t_str = a_str;
    else
        t_str = sprintf('%s %s', titlestr, a_str);
    end
    title(t_str);
    
    if ~isempty(axisprop)
        set(gca, axisprop{:});
    end
    
    if ~nosave
        savename = make_filename_safe(t_str);
        saveas(fig, fullfile(savedir, savename), saveformat);
        close(fig);
    end
end

end

