classdef GlobeGrid < handle
    %GLOBEGRID Representation of lat/lon grids
    %   This class allows you to generate a representation of a lat/lon
    %   grid over a specified domain on the Earth's surface.
    %
    %   Constructor:
    %       G = GlobeGrid(gridres) - Constructs a global, equirectangular
    %       grid with the resolution GRIDRES in both the longitudinal and
    %       latitudinal direction.
    %
    %       G = GlobeGrid(lonres, latres) - Constructs a global,
    %       equirectangular grid with different resolutions in the
    %       longitudinal and latitudinal direction.
    %
    %       The constructor accepts two parameters:
    %           'projection' - specify a different map projection;
    %           currently, only 'equirectangular' is implemented.
    %
    %           'domain' - allows you to specify the domain. See the
    %           SetDomain method below for the allowed syntax.
    %
    %   Properties (publicly settable):
    %       LonRes, LatRes - resolution in the longitudinal and latitudinal
    %       directions, respectively.
    %
    %       Projection - the map projection to use. Currently, the only
    %       allowed value is 'equirectangular'.
    %
    %   Properties (publicly accessible, not settable):
    %       DomainLon, DomainLat - the longitude and latitude bounds of the
    %       domain as two-element vectors. Modify using the SetDomain
    %       method.
    %
    %       GridLon, GridLat - the center points of the grid cells for the
    %       grid. 
    %
    %       GridLoncorn, GridLatcorn - the corner points of the grid cells.
    %       Assumes that the cells are tiled, i.e. there is no space or
    %       overlap between the grid cells. This means that the corners for
    %       the (1,1) grid cell would the the four value (1:2, 1:2) in the
    %       GridLoncorn and GridLatcorn arrays, for the (3,4) grid cell,
    %       the corners would be (3:4, 4:5), and so on.
    %
    %   Methods (public)
    %      SetDomain(domainspec) - set the domain boundaries. DOMAINSPEC
    %      can either be a four-element vector, in which the first two
    %      elements are the numerical longitude bounds, and the second two
    %      are the latitude bounds, or a string specifying a region. When
    %      specifying the domain bounds numerically, note that the order
    %      does not matter, but longitude must be given on the range [-180,
    %      180] and latitude on the range [-90, 90]. If using a string, the
    %      following strings are recognized:
    %           'world' - longitude [-180, 180]; latitude [-90, 90]
    %           'us' - longitude [-125, -65]; latitude [25 50]
    %
    %       GridcellContains(lon, lat) - returns indicies of the gridcell
    %       which contains the given longitude and latitude. If a point is
    %       exactly on the border of two cells, it is placed in the cell to
    %       the north and east. If two outputs requested, the indices in
    %       the first and second dimension are output separately. If one or
    %       zero outputs requestion, they are returned as a vector.
    
    properties
        LonRes
        LatRes
        Projection
        Rotation
    end
    properties(SetAccess = protected)
        DomainLon = [-180 180];
        DomainLat = [-90 90];
        GridCenterLon = 0;
        GridCenterLat = 0;
    end
    properties(Dependent = true, SetAccess = private)
        GridLon
        GridLat
        GridLoncorn
        GridLatcorn
    end
    
    methods
        function obj = GlobeGrid(lon_res, varargin)
            p = inputParser;
            p.addOptional('lat_res',[]);
            p.addParameter('projection','equirectangular');
            p.addParameter('domain', 'world');
            p.addParameter('rotation',0);
            p.parse(varargin{:});
            
            lat_res = p.Results.lat_res;
            proj = p.Results.projection;
            domain = p.Results.domain;
            rotation = p.Results.rotation;
            
            if nargin < 1
                error('globegrid:bad_input','GlobeGrid requires at least one input (grid resolution) to the constructor')
            end
            obj.LonRes = lon_res;
            if ~isempty(lat_res)
                obj.LatRes = lat_res;
            else
                obj.LatRes = lon_res;
            end
            obj.Projection = proj;
            if rotation ~= 0 && ~strcmp(proj,'equirect-rotated')
                warning('Rotation only affects equirect-rotated projections')
            end
            obj.Rotation = rotation;
            obj.SetDomain(domain);
        end
        
        function glon = get.GridLon(obj)
            glon = obj.get_grid_centers();
        end
        
        function glat = get.GridLat(obj)
            [~, glat] = obj.get_grid_centers();
        end
        
        function gloncorn = get.GridLoncorn(obj)
            gloncorn = obj.get_grid_corners();
        end
        
        function glatcorn = get.GridLatcorn(obj)
            [~, glatcorn] = obj.get_grid_corners();
        end
        
        function set.LonRes(obj, val)
            obj.LonRes = val;
            obj.check_domain;
        end
        
        function set.LatRes(obj, val)
            obj.LatRes = val;
            obj.check_domain;
        end
        
        function set.Projection(obj, val)
            allowed_proj = {'equirectangular','equirect-rotated'};
            if ~ischar(val) || ~ismember(val, allowed_proj)
                % Use the same identifier as when setting an invalid
                % property for a graphics handle, should help consistency
                error('MATLAB:datatypes:InvalidEnumValue','When setting the Projection value of GlobeGrid:\n''%s'' is not a valid value. It must be one of %s', val, strjoin(allowed_proj, ','))
            end
            obj.Projection = val;
        end
        
        function SetDomain(obj, varargin)
            % Set the domain either by name, by edges, or by center point
            % and width. 
            %
            % To set by name, give as input one of the strings
            % 'world' or 'us' ('us' covers the continental US 125 to 65 W,
            % 25 to 50 N). 
            %
            % To set by edges, give a 4-element vector with
            % values [lonmin, lonmax, latmin, latmax]. 
            %
            % To set by center and width requires two inputs: the center in
            % degrees as [lon, lat] and the width in degrees as [size west,
            % size east, size south, size north].
            if numel(varargin) == 1 && ischar(varargin{1})
                obj.SetDomainByName(varargin{1});
            elseif numel(varargin) == 1 && isnumeric(varargin{1}) && numel(varargin{1}) == 4
                obj.SetDomainByEdges(varargin{1});
            elseif numel(varargin) == 2 && all(iscellcontents(varargin,'isnumeric')) && numel(varargin{1}) == 2 && numel(varargin{2}) == 4
                obj.SetDomainByCenterWidth(varargin{:});
            else
                error('globegrid:bad_input', 'To set a domain, pass either a string, a four element numeric vector ([lonmin lonmax latmin latmax]) or a center point and width')
            end
        end
        
        function SetDomainByName(obj, domain)
            switch lower(domain)
                case 'world'
                    obj.DomainLon = [-180 180];
                    obj.DomainLat = [-90 90];
                case 'us'
                    obj.DomainLon = [-125 -65];
                    obj.DomainLat = [25 50];
                case 'hk'
                    obj.DomainLon = [108 118];
                    obj.DomainLat = [19 26];
                otherwise
                    error('globegrid:bad_input','Domain name %s not recognized', domain)
            end
            obj.GridCenterLon = mean(obj.DomainLon);
            obj.GridCenterLat = mean(obj.DomainLat);
            obj.check_domain()
        end
        
        function SetDomainByEdges(obj, domain)
            if ~isvector(domain) || numel(domain) ~= 4
                error('globegrid:bad_input', 'Numeric domain limits must be given as a four element vector')
            elseif any(domain(1:2) < -180 | domain(1:2) > 180)
                error('globegrid:bad_input', 'The first two elements of a numeric domain limit are longitude and must lie between -180 and 180')
            elseif any(domain(3:4) < -90 | domain(3:4) > 90)
                error('globegrid:bad_input', 'The third and fourth elements of a numeric domain limit are latitude and must lie between -90 and 90')
            end
            
            obj.DomainLon = [min(domain(1:2)), max(domain(1:2))];
            obj.DomainLat = [min(domain(3:4)), max(domain(3:4))];
            obj.GridCenterLon = mean(obj.DomainLon);
            obj.GridCenterLat = mean(obj.DomainLat);
                
            obj.check_domain()
        end
        
        function SetDomainByCenterWidth(obj, center, width)
            % Width(1) is the size west of center, width(2) is the size
            % east of center, width(3) size south of center, width(4) size
            % north of center.
            
            if center(1) < -180 || center(1) > 180
                error('globegrid:bad_input', 'The first element of a domain center is longitude and must lie between -180 and 180')
            elseif center(2) < -90 || center(2) > 90
                error('globegrid:bad_input', 'The first element of a domain center is latitude and must lie between -90 and 90')
            elseif any(width < 0)
                error('globegrid:bad_input', 'All elements of a domain width must be positive')
            end
            
            obj.GridCenterLon = center(1);
            obj.GridCenterLat = center(2);
            
            res = [obj.LonRes, obj.LonRes, obj.LatRes, obj.LatRes];
            % Ensure that width is a multiple of resolution and add the
            % offset properly so that the center point is the center of a
            % grid cell
            width = ceil((width + res./2) ./ res) .* res - res./2;
            
            obj.DomainLon = [center(1) - width(1), center(1) + width(2)];
            obj.DomainLat = [center(2) - width(3), center(2) + width(4)];
            
            obj.check_domain()
        end
        
        function [varargout] = GridcellContains(obj, lon, lat, makebool)
            if ~exist('makebool', 'var')
                makebool = false;
            end
            switch obj.Projection
                case 'equirectangular'
                    [x, y] = obj.equirect_contains(lon, lat);
                case 'equirect-rotated'
                    [x, y] = obj.equirect_rot_contains(lon, lat);
                otherwise
                    error('globegrid:projection','Projection %s not implemented in GridcellContains method',obj.Projection)
            end
            
            if ~makebool
                x = find(x);
                y = find(y);
            end
            
            if nargout < 2
                varargout{1} = [x, y];
            elseif nargout == 2
                varargout{1} = x;
                varargout{2} = y;
            else
                error('MATLAB:TooManyOutputs', 'Too many output arguments')
            end
        end
        
        function grid_info = OmiGridInfo(obj)
            % Returns a representation of this grid as a Python omi.Grid
            % object.
            if ~strcmpi(obj.Projection, 'equirectangular')
                error('globegrid:omi_grid_format', 'Cannot convert a non-equirectangular grid to a Python omi.Grid instance');
            elseif obj.LonRes ~= obj.LatRes
                error('globegrid:omi_grid_format', 'Cannot convert a grid with different lon/lat resolutions to a Python omi.Grid instance');
            end
            
            % The omi.Grid constructor takes the lower left and upper right
            % corners. However, the upper right corner is exclusive (i.e.
            % it's not included) so we add an extra bit to that corner to
            % ensure the right- and top- most rows are included
            
            info_struct.llcrnrlon = min(obj.GridLon(:,1));
            info_struct.urcrnrlon = max(obj.GridLon(:,1)) + obj.LonRes / 2;
            info_struct.llcrnrlat = min(obj.GridLat(1,:));
            info_struct.urcrnrlat = max(obj.GridLat(1,:)) + obj.LatRes / 2;
            info_struct.resolution = obj.LonRes;
            
            grid_info = struct2pydict(info_struct);
        end
    end
    
    methods(Access=private)
        function check_domain(obj)
            % Check that the resolution divides the domain evenly, if we're just using the domain boundaries
            % (not a center coordinate)
            warncell = {};
            if mod(obj.DomainLon(2) - obj.DomainLon(1), obj.LonRes) ~= 0
                warncell{end+1} = 'longitudinal';
            end
            if mod(obj.DomainLat(2) - obj.DomainLat(1), obj.LatRes) ~= 0
                warncell{end+1} = 'latitudinal';
            end
            if ~isempty(warncell)
                warning('The domain is not a multiple of the resolution in the %s direction', strjoin(warncell, ' and '))
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Grid generative functions %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [grid_lon, grid_lat] = get_grid_centers(obj)
            switch obj.Projection
                case 'equirectangular'
                    [grid_lon, grid_lat] = obj.make_equirect_grid_centers();
                case 'equirect-rotated'
                    [grid_lon, grid_lat] = obj.make_equirect_rot_grid_centers();
                otherwise
                    error('globegrid:projection','Projection %s not implemented in get_grid_centers method',obj.Projection)
            end
        end
        
        function [grid_lon, grid_lat] = get_grid_corners(obj)
            switch obj.Projection
                case 'equirectangular'
                    [grid_lon, grid_lat] = obj.make_equirect_grid_corners();
                case 'equirect-rotated'
                    [grid_lon, grid_lat] = obj.make_equirect_rot_grid_corners();
                otherwise
                    error('globegrid:projection','Projection %s not implemented in get_grid_corners method',obj.Projection)
            end
        end

        function [grid_lon, grid_lat] = make_equirect_grid_centers(obj)
            lonvec = (obj.DomainLon(1)+obj.LonRes/2):obj.LonRes:(obj.DomainLon(2)-obj.LonRes/2);
            latvec = (obj.DomainLat(1)+obj.LatRes/2):obj.LatRes:(obj.DomainLat(2)-obj.LatRes/2);
            [grid_lat, grid_lon] = meshgrid(latvec, lonvec);
        end
        
        function [grid_lon, grid_lat] = make_equirect_grid_corners(obj)
            lonvec = obj.DomainLon(1):obj.LonRes:obj.DomainLon(2);
            latvec = obj.DomainLat(1):obj.LatRes:obj.DomainLat(2);
            [grid_lat, grid_lon] = meshgrid(latvec, lonvec);
        end
        
        function [grid_lon, grid_lat] = make_equirect_rot_grid_centers(obj)
            [grid_lon, grid_lat] = obj.make_equirect_grid_centers();
            [grid_lon, grid_lat] = obj.rotate_coordinates(grid_lon, grid_lat, obj.GridCenterLon, obj.GridCenterLat, obj.Rotation);
        end
        
        function [grid_lon, grid_lat] = make_equirect_rot_grid_corners(obj)
            [grid_lon, grid_lat] = obj.make_equirect_grid_corners();
            [grid_lon, grid_lat] = obj.rotate_coordinates(grid_lon, grid_lat, obj.GridCenterLon, obj.GridCenterLat, obj.Rotation);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Grid query functions %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function [xx, yy] = equirect_contains(obj, lon, lat)
            % Because each row is the same longitude and each column the
            % same latitude, we only have to do a logical comparison on the
            % first column and row of the grids
            xx = lon >= obj.GridLoncorn(1:end-1,1) & lon < obj.GridLoncorn(2:end,1);
            yy = lat >= obj.GridLatcorn(1,1:end-1) & lat < obj.GridLatcorn(1,2:end);
            if sum(xx) == 0 || sum(yy) == 0
                % Can have cases where one index is found, but not the
                % other. If one index can't find a grid cell in the domain,
                % don't let either one do so.
                xx(:) = false;
                yy(:) = false;
            end
        end
        
        function equirect_rot_contains(obj, lon, lat) %#ok<INUSD>
            error('globegrid:not_implemented','GridcellContains not implemented for the equirect_rotated projection')
        end
        
    end
    
        %%%%%%%%%%%%%%%%%%%%%
        % Utility functions %
        %%%%%%%%%%%%%%%%%%%%%
    methods(Static)
        function [lon, lat] = rotate_coordinates(lon, lat, center_lon, center_lat, theta)
            if ~isnumeric(lon) || ~isnumeric(lat)
                error('globegrid:internal:rotation','LON and LAT must both be numeric')
            elseif ~isequal(size(lon), size(lat))
                error('globegrid:internal:rotation','LON and LAT must be the same size')
            elseif ~isscalar(center_lon) || ~isscalar(center_lat) || ~isnumeric(center_lon) || ~isnumeric(center_lat)
                error('globegrid:internal:rotation', 'CENTER_LON and CENTER_LAT must be scalar numbers')
            elseif ~isscalar(theta)
                error('globegrid:internal:rotation', 'THETA must be a scalar number')
            end
            R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
            lon_tmp = lon - center_lon;
            lat_tmp = lat - center_lat;
            for a=1:numel(lon_tmp)
                v = R * [lon_tmp(a); lat_tmp(a)];
                lon(a) = v(1) + center_lon;
                lat(a) = v(2) + center_lat;
            end
        end
    end
    
end

