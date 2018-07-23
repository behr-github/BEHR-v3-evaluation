function [ varargout ] = state_outlines( varargin )
%STATE_OUTLINES Draws state outlines using MATLAB shapefiles
%   STATE_OUTLINES() will draw all US states on the current figure in
%   black.
%
%   STATE_OUTLINES( FIGNUM ) will draw all US states on Figure number
%   FIGNUM in black.
%
%   STATE_OUTLINES( COLSPEC ) will draw all US stats on the current figure
%   using the color defined by COLSPEC. COLSPEC can be one of the Matlab
%   single character color specifications (see
%   http://www.mathworks.com/help/matlab/ref/colorspec.html) of a 1x3
%   vector specifying an RGB triplet.
%
%   STATE_OUTLINES( FIGNUM, COLSPEC ) combined the previous two syntaxes.
%
%   STATE_OUTLINES( ___, STATES ) can be used to plot a specific subset of
%   states, defined by their two letter postal abbreviation or the full
%   name. This can be combined with any of the previous syntaxes to specify
%   a figure number, color spec, both, or neither.
%
%   STATE_OUTLINES( ___, 'not', STATES ) can be used to print all states
%   EXCEPT those specified. Again, it can be combined with any of the first
%   four syntaxes to specify figure number, color spec, both, or neither.
%
%   FIG = STATE_OUTLINES( ___ ) will return the handle to the figure that
%   the outlines were drawn in.
%
%   [ SLON, SLAT ] = STATE_OUTLINES( ___ ) will not plot anything, but
%   instead return the lon and lat coordinates to draw the state outlines
%   yourself.
%
%   Examples:
%       state_outlines() will plot all US states on the current figure.
%
%       state_outlines(5) will plot all US states on Figure 5.
%
%       state_outlines('b') will plot all US states on the current figure
%       in blue.
%
%       state_outlines(2, 'r') will plot all US states on Figure 2 in red.
%
%       state_outlines('pa','oh','wv','va') will plot only Pennsylvania,
%       Ohio, West Virginia, and Virginia on the current figure in black.
%
%       state_outlines(3, 'b', 'not', 'ak', 'hi') will plot all states
%       EXCEPT Alaska and Hawaii in blue on Figure 3.
%
%   Josh Laughner <joshlaugh5@gmail.com> Jul 2014

argin = varargin;
fignum = -1;
if numel(argin) > 0
    if isnumeric(argin{1}) && numel(argin{1}) == 1;
        fignum = argin{1};
        argin(1) = [];
    elseif ~ischar(argin{1}) && ishandle(argin{1}) && strcmp(get(argin{1},'Type'),'figure')
        fignum = get(argin{1},'Number');
        argin(1) = [];
    elseif ~isempty(get(0,'children'))
        fignum = get(gcf,'Number');
    end
    
    if ischar(argin{1}) && length(argin{1}) == 1;
        colspec = argin{1};
        argin(1) = [];
    elseif isnumeric(argin{1}) && length(argin(1))==3;
        colspec = argin{1};
        argin(1) = [];
    else
        colspec = 'k';
    end
else
    if ~isempty(get(0,'children'))
        fignum = get(gcf,'Number');
    end
    colspec = 'k';
end

not_states = false;
if nargin == 0 || isempty(argin);
    states = 'all';
elseif strcmpi(argin{1},'not')
    states = argin(2:end);
    not_states = true;
else
    states = argin;
end

state_abbrev = {'al','ak','az','ar','ca','co','ct','de','fl','ga','hi','id','il','in','ia','ks','ky','la','me','md','ma','mi','mn','ms','mo','mt','ne','nv','nh','nj','nm','ny','nc','nd','oh','ok','or','pa','ri','sc','sd','tn','tx','ut','vt','va','wa','wv','wi','wy'};
state_names = {'alabama','alaska','arizona','arkansas','california','colorado','connecticut',...
    'delaware','florida','georgia','hawaii','idaho','illinois','indiana','iowa','kansas',...
    'kentucky','louisiana','maine','maryland','massachusetts','michigan','minnesota','mississippi',...
    'missouri','montana','nebraska','nevada','new_hampshire','new_jersey','new_mexico','new_york','north_carolina','north_dakota',...
    'ohio','oklahoma','oregon','pennsylvania','rhode_island','south_carolina','south_dakota','tennessee','texas','utah','vermont','virginia','washington','west_virginia',...
    'wisconsin','wyoming'};

% Read the states shapefile and plot
usa = shaperead('usastatehi.shp');




draw_lines = true;
if nargout <= 1
    if fignum < 1
        fig_handle = figure;
    else
        fig_handle = figure(fignum);
    end
    varargout{1} = fig_handle; 
elseif nargout == 2
    draw_lines = false;
    lon_out = [];
    lat_out = [];
end

% Figure out is the axes are a map or not. This may break in matlab
% versions before 2014b, when the graphics object was introduced.
im_a_map = false;
if draw_lines
    ch = get(figure(fignum),'children');
    for a=1:numel(ch)
        if strcmpi(get(ch(a),'Type'),'axes')
            im_a_map = ismap(ch(a));
            break
        end
    end
    
    % This will be necessary to keep the current x and y limits, for some
    % reason, setting xlimmode and ylimmode to 'manual' doesn't work. Save
    % the current limits, and reset them at the end.
    old_xlim = get(gca, 'xlim');
    old_ylim = get(gca, 'ylim');
end



if strcmpi(states,'all');
    for a=1:numel(usa)
        if im_a_map && draw_lines
            linem(usa(a).Y, usa(a).X, 'color', colspec);
        elseif draw_lines
            line(usa(a).X, usa(a).Y,'color',colspec);
        else
            lon_out = cat(1, lon_out, usa(a).X(:), nan);
            lat_out = cat(1, lat_out, usa(a).Y(:), nan);
        end
        hold on
    end
elseif not_states
    states_not_to_draw = cell(1,numel(states));
    for c=1:numel(states) 
        if length(states{c}) == 2
            xx = find(strcmpi(states{c},state_abbrev));
        else
            xx = find(strcmpi(states{c},state_names));
        end
        states_not_to_draw{c} = state_names{xx};
    end
    for b = 1:numel(state_names)
        if any(strcmpi(usa(b).Name,states_not_to_draw))
            continue
        else
            if im_a_map && draw_lines
                linem(usa(b).Y, usa(b).X, 'color', colspec);
            elseif draw_lines
                line(usa(b).X, usa(b).Y, 'color',colspec)
            else
                lon_out = cat(1, lon_out, usa(b).X(:), nan);
                lat_out = cat(1, lat_out, usa(b).Y(:), nan);
            end
            
        end
    end
else
    for b = 1:numel(states);
        if length(states{b}) == 2
            xx = find(strcmpi(states{b},state_abbrev));
        else
            xx = find(strcmpi(states{b},state_names));
        end
        if im_a_map && draw_lines
            linem(usa(xx).Y, usa(xx).X, 'color', colspec);
        elseif draw_lines
            line(usa(xx).X, usa(xx).Y,'color',colspec);
        else
            lon_out = cat(1, lon_out, usa(xx).X(:), nan);
            lat_out = cat(1, lat_out, usa(xx).Y(:), nan);
        end
    end
end

if ~draw_lines
    varargout{1} = lon_out;
    varargout{2} = lat_out;
else
    set(gca,'xlim',old_xlim,'ylim',old_ylim);
end
end

