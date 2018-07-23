function [ varargout ] = misc_gc_plotting_fxns( plttype, varargin )
%GC_PLOTTING_FXNS Misc. plots to make for GC outputs
%   Rather than trying to do this all in the command window, or having a
%   zillion different .m files all over the place, I'm just going to
%   collect them all here.
%
%   Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2015

DEBUG_LEVEL = 2;
E = JLLErrors;
HOMEDIR = getenv('HOME');

plttype = lower(plttype);
switch plttype
    case 'subfrac_lnox'
        [varargout{1}, varargout{2}] = subfrac_lnox(varargin{1});
    case 'trop_mask'
        [varargout{1}] = trop_mask(varargin{1:2});
    case 'region_xxyy'
        [varargout{1}, varargout{2}] = region_xxyy(varargin{1});
        [~,~,varargout{3}] = define_regions(varargin{1});
    case 'timeser_column'
        plot_dcolumn_percentiles(varargin{:}, 'relative', false);
    case 'timeser_relcol'
        plot_dcolumn_percentiles(varargin{:}, 'relative', true);
    case 'timeser_plev'
        plot_level_dconc(varargin{:}, 'relative', false);
    case 'timeser_relplev'
        plot_level_dconc(varargin{:}, 'relative', true);
    case 'timeser_dlnox'
        plot_lnox_col_enhnc_diff(varargin{:}, 'relative', false);
    case 'timeser_reldlnox'
        plot_lnox_col_enhnc_diff(varargin{:}, 'relative', true);
    case 'meridional'
        plot_meridional_avg(varargin{:});
    case 'meridional-sat'
        add_sat_mer_avg_to_fig(varargin{:});
    case 'meridional_table'
        varargout{1} = avg_meridional_avg_table(varargin{:});
    case 'grid_mls'
        varargout{1} = grid_mls_data_to_gc(varargin{:});
    case 'grid_mls_daily'
        [varargout{1}, varargout{2}] = grid_daily_mls_data_to_gc(varargin{:});
    case 'grid_month_omno2'
        [varargout{1}, varargout{2}] = grid_omno2d_monthly(varargin{:});
    case 'grid_day_omno2'
        varargout{1} = grid_omno2d_daily(varargin{:});
    case 'dc3-diff'
        compare_cases_during_DC3();
    case 'wt-check'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = check_ak_vcd_domino_weights();
    case 'comp2sat'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}, varargout{8}] = compare_regions_to_sat(varargin{:});
    case 'plotprofs'
        plot_avg_gc_profs();
    case 'strat_no2'
        domino_sp_scatter();
    case 'compare-aks'
        compare_ak_vec();
    case 'compare-raw-aks'
        compare_raw_aks()
    case 'emissions'
        compute_emissions();
    otherwise
        fprintf('Did not recognize plot type\n');
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OTHER FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [prodfrac, lnoxgt60] = subfrac_lnox(prodfrac)
        % Calculates the fraction of NO emission due to lightning only
        % considering anthropogenic and biomass burning as other sources
        % Takes a production structure output from geos_frac_lnox
        subtot = prodfrac.NO_molec_cm2_sec.Anthropogenic_NO + prodfrac.NO_molec_cm2_sec.Biomass_NO + prodfrac.NO_molec_cm2_sec.Lightning_NO;
        subfrac = prodfrac.NO_molec_cm2_sec.Lightning_NO ./ subtot;
        prodfrac.NO_frac.Subfrac_LNOx = subfrac;
        lnoxgt60 = subfrac > 0.6;
    end

    function dataind = get_data_inds(casestruct)
        % Finds which indices correspond to pressure surface, box height,
        % number density, and tropopause level. Return as a struct with
        % fields bxhght, psurf, ndens, and tplev.
        dataind = struct('bxhght',nan,'psurf',nan,'ndens',nan,'tplev',nan);
        for n=1:numel(casestruct)
            if ~isempty(regexpi(casestruct(n).fullName,'bxheight'));
                dataind.bxhght = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'psurf|pressure'))
                dataind.psurf = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'number density'))
                dataind.ndens = n;
            elseif ~isempty(regexpi(casestruct(n).fullName,'tropopause'))
                dataind.tplev = n;
            end
        end
    end

    function [relative, ndens_bool, lonlim, latlim, timelim_bool] = parse_varargs(vargs)
        % Parse varargs. If one is the string 'relative', the next one is
        % the relative boolean. Also look for the string 'ndens', if there,
        % set it to use number density (and remove it). There should be 0
        % to 2 left after that, if 2, then those are the lat and lon
        % limits. If 1, it should be a region definition
        xx = strcmpi('relative',vargs);
        if sum(xx) == 1
            ff = find(xx);
            relative = vargs{ff+1};
            vargs(ff:(ff+1)) = [];
        elseif sum(xx) == 0
            relative = false;
        else
            E.badinput('The relative keyword appears more than once');
        end
        
        xx = strcmpi('ndens',vargs);
        if sum(xx) > 0
            ndens_bool = true;
            vargs(xx) = [];
        elseif sum(xx) == 0
            ndens_bool = false;
        end
        
        xx = strcmpi('timelim',vargs);
        if sum(xx) > 0
            timelim_bool = true;
        else
            timelim_bool = false;
        end
        
        if numel(vargs) == 0
            lonlim = [-200 200];
            latlim = [-100 100];
        elseif numel(vargs) == 2
            lonlim = vargs{1};
            latlim = vargs{2};
        elseif numel(vargs) == 1 && ischar(vargs{1}) && ismember(vargs{1}, {'na','sa','naf','saf','seas','neur'});
            [lonlim, latlim] = define_regions(vargs{1});
        else
            E.badinput('Parsing of optional arguments failed')
        end
    end

    function [lonlim, latlim, timeind] = define_regions(region)
        % If limiting by time, define the Northern Hemisphere to be
        % June-Aug and the Southern Hemisphere to be Jan, Feb, Dec.
        nh_times = modis_date_to_day('2012-06-01'):modis_date_to_day('2012-08-31');
        nh_mask = false(1,366);
        nh_mask(nh_times) = true;
        sh_times = [modis_date_to_day('2012-01-01'):modis_date_to_day('2012-02-29'), modis_date_to_day('2012-12-01'):modis_date_to_day('2012-12-31')];
        sh_mask = false(1,366);
        sh_mask(sh_times) = true;
        switch region
            case 'na'
                lonlim = [-120, -65];
                latlim = [20, 60];
                %latlim = [20, 50];
                timeind = nh_mask;
            case 'sa'
                lonlim = [-77, -39];
                latlim = [-35, 10];
                timeind = sh_mask;
            case 'naf'
                lonlim = [-15, 48];
                latlim = [3, 25];
                timeind = nh_mask;
            case 'saf'
                lonlim = [10, 48];
                latlim = [-30, 3];
                timeind = sh_mask;
            case 'seas'
                lonlim = [95, 146];
                latlim = [-9, 26];
                timeind = nh_mask;
            case 'neur'
                lonlim = [60, 130];
                latlim = [30, 68];
                timeind = nh_mask;
            case 'nh'
                lonlim = [-180,180];
                latlim = [0, 60];
                timeind = nh_mask;
            case 'sh'
                lonlim = [-180, 180];
                latlim = [-60, 60];
                timeind = sh_mask;
            case 'natl'
                lonlim = [-50 -25];
                latlim = [20 40];
                timeind = nh_mask;
            case 'spac'
                lonlim = [-160 -120];
                latlim = [-45 -10];
                timeind = sh_mask;
            otherwise
                lonlim = [-200, 200];
                latlim = [-100, 100];
                timeind = true(1,366);
        end
    end

    function [xx,yy] = region_xxyy(region)
        [lonlim, latlim] = define_regions(region);
        [glon, glat] = geos_chem_centers('2x25');
        xx = glon >= lonlim(1) & glon <= lonlim(2);
        yy = glat >= latlim(1) & glat <= latlim(2);
    end

    function [mask] = trop_mask(sz, TP)
        % Returns a mask (logical matrix) that is 1 for any grid box in the
        % troposphere (below the tropopause level) and 0 otherwise.
        % Requires the size of the mask desired (should match the size of
        % the dataBlock it will be applied to) and TP is a structure output
        % from read_geos_output with the tropopause level as the dataBlock.
        % The first, second, and final dimensions of the TP dataBlock must
        % be the same size as the corresponding dimensions in the mask
        % desired.
        
        if ~isvector(sz) || numel(sz) < 3 || numel(sz) > 4
            E.badinput('sz must be a size vector with 3 or 4 elements')
        elseif ~isstruct(TP) || ~isfield(TP,'dataBlock') || ~isfield(TP,'fullName') || ~isscalar(TP) 
            E.badinput('TP must be a scalar structure output from read_geos_output')
        elseif isempty(regexpi(TP.fullName,'tropopause'))
            E.badinput('TP does not seem to contain the tropopause level output')
        end
        
        sz_tp = size(TP.dataBlock);
        if numel(sz_tp) > 3 || numel(sz_tp) < 2
            E.badinput('The tropopause dataBlock is expected to have 2 or 3 dimensions');
        elseif any(sz_tp(1:2) ~= sz(1:2)) 
            E.badinput('The first two dimensions of the requested mask size must be the same as the first two dimensions of the TP dataBlock')
        elseif sz_tp(3) ~= sz(end)
            E.badinput('The final dimension of the requested mask is expected to be the same size as the third dimension of the TP dataBlock')
        end
        
        tplev = floor(TP.dataBlock);
        
        mask = true(sz);
        
        for t=1:sz_tp(3)
            for a=1:sz_tp(1)
                for b=1:sz_tp(2)
                    mask(a,b,tplev(a,b,t):end,t) = false;
                end
            end
        end
    end

    

    function [sat_lon, sat_avg] = sat_meridional_avg(specie, region, gridded_data, lnox_crit, noregrid_bool, timelim_bool)
        % Calculates a meridional average for satellite measurements of
        % NO2 (column), O3, or HNO3. Needs which species (as a string),
        % one of the regions (also a string). Hard coded to do 2012 right
        % now.
        [xx,yy] = region_xxyy(region);
        % Add one more point because we actually want to subset the corners
        % not the centers.
        xx_centers = xx; % we'll need this at the end.
        yy_centers = yy;
        xx(find(xx,1,'last')+1) = true;
        yy(find(yy,1,'last')+1) = true;
        [gloncorn, glatcorn] = geos_chem_corners;
        % geos_chem_corners outputs in lat x lon, so transpose just because
        % I'm more used (now) to having lon first.
        gloncorn = gloncorn(yy,xx)';
        glatcorn = glatcorn(yy,xx)';
        
        [lonlim, latlim, timemask] = define_regions(region);
        
        switch lower(specie)
            case 'no2'
                if ~exist('gridded_data','var') || isempty(gridded_data)
                    addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/OMI Utils');
                    gridded_data = omno2d_timeavg('2012-01-01','2012-12-31');
                end
                [omno2d_lon, omno2d_lat] = omno2d_centers;
                
                xx_omi = omno2d_lon(:,1) >= lonlim(1) & omno2d_lon(:,1) <= lonlim(2);
                yy_omi = omno2d_lat(1,:) >= latlim(1) & omno2d_lat(1,:) <= latlim(2);
                
                if exist('lnox_crit','var') && ~isempty(lnox_crit)
                    gridded_data(~lnox_crit) = nan;
                end
                
                if timelim_bool
                    if size(gridded_data,3) == 366
                        gridded_data(:,:,~timemask) = nan;
                    else
                        E.badinput('If you wish to filter by time, the input data must represent 366 days');
                    end
                end
                
                if size(gridded_data,3)>1
                    %gridded_data = nanmean(gridded_data,3);
                    gridded_data = nanmedian(gridded_data,3);
                end
                
                if size(gridded_data,1) * size(gridded_data,2) > 144*91
                    if noregrid_bool
                        gridded_data = gridded_data(xx_omi,yy_omi);
                    else
                        gridded_data = grid_omno2d_to_gc(gridded_data, omno2d_lon, omno2d_lat, gloncorn, glatcorn);
                    end
                else
                    gridded_data = gridded_data(xx_centers, yy_centers);
                end

                if noregrid_bool
                    sat_lon = omno2d_lon(xx_omi,1);
                else
                    sat_lon = geos_chem_centers('2x25');
                    sat_lon = sat_lon(xx_centers);
                end
            otherwise
                if ~exist('gridded_data','var') || isempty(gridded_data)
                    gridded_data = grid_mls_data_to_gc(gloncorn, glatcorn,specie);
                    if timelim_bool
                        warning('Time limits not applied; you must pass an existing matrix for that to work')
                    end
                else
                    if exist('lnox_crit','var') && ~isempty(lnox_crit)
                        gridded_data(~lnox_crit) = nan;
                    end
                    
                    if timelim_bool
                        if size(gridded_data,3) == 366
                            gridded_data(:,:,~timemask) = nan;
                        else
                            E.badinput('If you wish to filter by time, the input data must represent 366 days');
                        end
                    end
                    
                    %gridded_data = nanmean(gridded_data,3);
                    gridded_data = nanmedian(gridded_data,3);
                    gridded_data = gridded_data(xx_centers,yy_centers);
                end
                
                
                
                sat_lon = geos_chem_centers('2x25');
                sat_lon = sat_lon(xx_centers);
        end
        
        % Make the meridional average
        %sat_avg = nanmean(gridded_data,2)';
        sat_avg = nanmedian(gridded_data,2)';
        
    end

    function gridded_daily_omno2d = grid_omno2d_daily(gc_loncorn, gc_latcorn)
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/OMI Utils');
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        omno2d_path = '/Volumes/share-sat/SAT/OMI/OMNO2d';
        if ~exist(omno2d_path,'dir')
            E.dir_dne('omno2d_path')
        end
        [omno2d_lon, omno2d_lat] = omno2d_centers;
        
        fname = sprintf('OMI-Aura_L3-OMNO2d_%04dm%02d%02d_v003*.he5', year(dates(d)), month(dates(d)), day(dates(d)));
        F = dir(fullfile(omno2d_path, sprintf('%04d',year(dates(d))),fname));
        if isempty(F) && ~ignore_missing
            E.filenotfound(fname);
        elseif numel(F) > 1
            E.toomanyfiles(fname);
        end
        
        gridded_daily_omno2d = nan([size(gc_loncorn)-1, 366]);
        
        for d=1:366
            sdate = modis_day_to_date(d,2012);
            
            hinfo = h5info(fullfile(omno2d_path, sprintf('%04d',year(sdate)), F(1).name));
            todays_no2 = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2TropCloudScreened'));
            todays_weight = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'Weight'));
            
            % Ignore negative values - they are either fill values or unphysical
            todays_weight(todays_no2<0) = 0;
            
            gridded_daily_omno2d(:,:,d) = grid_omno2d_to_gc(todays_no2, omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn, todays_weight);
        end
    end

    function [gridded_monthly_omno2d, monthly_omno2d] = grid_omno2d_monthly(gc_loncorn, gc_latcorn)
        % Does monthly time averages of OMNO2d data then grids it to GC
        % resolution.
        [omno2d_lon, omno2d_lat] = omno2d_centers;
        monthly_omno2d = nan(size(omno2d_lon));
        gridded_monthly_omno2d = nan(size(gc_loncorn)-1);
        for m=1:12
            sdate = sprintf('2012-%02d-01',m);
            edate = sprintf('2012-%02d-%02d',m,eomday(2012,m));
            monthly_omno2d(:,:,m) = omno2d_timeavg(sdate,edate);
            gridded_monthly_omno2d(:,:,m) = grid_omno2d_to_gc(monthly_omno2d(:,:,m), omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn);
        end
        
    end

    function gridded_omno2d = grid_omno2d_to_gc(omno2d_no2, omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn, omno2d_weights)
        % Averages OMNO2d data to GEOS-Chem resolution. Needs OMNO2d data,
        % longitude, and latitude (all should be the same size) and the
        % GEOS-Chem corner points to average to.
        
        % Check input
        if ndims(omno2d_no2) ~= ndims(omno2d_lon) || ndims(omno2d_no2) ~= ndims(omno2d_lat) || any(size(omno2d_no2) ~= size(omno2d_lon)) || any(size(omno2d_no2) ~= size(omno2d_lat))
            E.sizeMismatch('omno2d_lon','omno2d_lat','omno2d_no2')
        end
        if ~exist(omno2d_weights)
            omno2d_weights = ones(size(omno2d_no2));
        end
        
        gridded_omno2d = nan(size(gc_loncorn)-1); % one smaller to switch from # of corners to # of grid cells
        
        % Loop over GC grid cells and find which OMNO2d data belongs in it
        for a=1:size(gridded_omno2d,1)
            if DEBUG_LEVEL > 1; fprintf('%d ',a); end
            for b=1:size(gridded_omno2d,2)
                xx = omno2d_lon >= gc_loncorn(a,b) & omno2d_lon < gc_loncorn(a+1,b+1) & omno2d_lat >= gc_latcorn(a,b) & omno2d_lat < gc_latcorn(a+1,b+1);
                gridded_omno2d(a,b) = nansum(omno2d_no2(xx) .* omno2d_weights(xx)) / nansum(omno2d_weights(xx));
            end
        end
        if DEBUG_LEVEL > 1; fprintf('\n'); end
    end

    function [gridded_mls, count] = grid_daily_mls_data_to_gc(gc_loncorn, gc_latcorn, specie, omi_overpass)
        % Calculates a running average of MLS data for the given GC
        % longitude and latitude corners.  Intended to be called from
        % sat_meridional_avg() or externally to pre-grid data and save
        % yourself time. (In that case, pass a full gc_loncorn/latcorn
        % matrix - don't cut it down into regions.)
        %
        % One optional argument, set to true to limit to approximately OMI
        % overpass time (will restrict to 11:00 - 15:00 solar time).
        
        if ~exist('omi_overpass','var')
            omi_overpass = false;
        end
        
        switch lower(specie)
            case 'o3'
                filepath = '/Volumes/share-sat/SAT/MLS/O3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [350 200];
            case 'hno3'
                filepath = '/Volumes/share-sat/SAT/MLS/HNO3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                %pressure_range = [250 200];
                %pressure_range = [250 140];
                pressure_range = [200 140];
        end
        

        gridded_mls = nan([size(gc_loncorn)-1, 366]);
        count = zeros([size(gc_loncorn)-1, 366]);

        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        
        fnames = {files.name};
        
        % Do the gridding and averaging
        d_ind = 0;
        for d=datenum('2012-01-01'):datenum('2012-12-31')
            if DEBUG_LEVEL > 1
                fprintf('Binning %s\n',datestr(d));
            end
            d_ind = d_ind + 1;
            fileglob = sprintf('MLS-Aura_L2GP-%s_v04-2[0-9]-c01_2012d%03d.he5',upper(specie),modis_date_to_day(d));
            filename = glob(fnames, fileglob);
            if isempty(filename)
                continue
            end
            Data = read_mls_data(fullfile(filepath,filename{1}));
            lon = Data.Longitude;
            lat = Data.Latitude;
            solar_time = Data.LocalSolarTime;
            
            % Find which levels of the MLS retrieval are in the desired
            % pressure range
            pp = Data.Pressure >= min(pressure_range) & Data.Pressure <= max(pressure_range);
            
            % Restrict (more or less) to OMI overpass time based on the
            % local solar time of the measurement. Mainly intended to
            % remove nighttime rather than completely limit to OMI
            % overpass.
            if omi_overpass
                tt = solar_time >= 11 & solar_time <= 15;
                Data.(upper(specie))(:,~tt) = nan;
            end
            
            for a = 1:size(gridded_mls,1)
                for b = 1:size(gridded_mls,2)
                    % Find the elements from the MLS data that fall in each
                    % GC grid cell and add them to the running average
                    xx = lon >= gc_loncorn(a,b) & lon < gc_loncorn(a+1,b+1) & lat >= gc_latcorn(a,b) & lat < gc_latcorn(a+1,b+1);
                    if sum(xx) < 1
                        continue
                    end
                    c = sum(xx);
                    data = nansum2(Data.(upper(specie))(pp,xx));
                    
                    gridded_mls(a, b, d_ind) = nansum2([gridded_mls(a,b) * count(a,b), data]) / (count(a,b) + c);
                    count(a, b, d_ind) = count(a,b) + c;
                end
            end
        end
    end

    function gridded_mls = grid_mls_data_to_gc(gc_loncorn, gc_latcorn, specie)
        % Calculates a running average of MLS data for the given GC
        % longitude and latitude corners.  Intended to be called from
        % sat_meridional_avg() or externally to pre-grid data and save
        % yourself time. (In that case, pass a full gc_loncorn/latcorn
        % matrix - don't cut it down into regions.)
        
        switch lower(specie)
            case 'o3'
                filepath = '/Volumes/share-sat/SAT/MLS/O3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [350 200];
            case 'hno3'
                filepath = '/Volumes/share-sat/SAT/MLS/HNO3v4-2/2012';
                files = dir(fullfile(filepath,'*.he5'));
                pressure_range = [250 200];
        end
        
        if temporal_avg
            gridded_mls = nan(size(gc_loncorn)-1);
            count = zeros(size(gc_loncorn)-1);
        else
            gridded_mls = nan([size(gc_loncorn)-1, numel(files)]);
            count = zeros([size(gc_loncorn)-1, numel(files)]);
        end
        addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/MLS');
        
        % Do the gridding and averaging
        for f=1:numel(files)
            if DEBUG_LEVEL > 1
                fprintf('Binning file %d of %d\n',f,numel(files));
            end
            Data = read_mls_data(fullfile(filepath,files(f).name));
            lon = Data.Longitude;
            lat = Data.Latitude;
            
            % Find which levels of the MLS retrieval are in the desired
            % pressure range
            pp = Data.Pressure >= min(pressure_range) & Data.Pressure <= max(pressure_range);
            
            for a = 1:size(gridded_mls,1)
                for b = 1:size(gridded_mls,2)
                    % Find the elements from the MLS data that fall in each
                    % GC grid cell and add them to the running average
                    xx = lon >= gc_loncorn(a,b) & lon < gc_loncorn(a+1,b+1) & lat >= gc_latcorn(a,b) & lat < gc_latcorn(a+1,b+1);
                    if sum(xx) < 1
                        continue
                    end
                    c = sum(xx);
                    d = nansum(Data.(upper(specie))(pp,xx));
                    
                    if temporal_avg
                        gridded_mls(a,b) = nansum2([gridded_mls(a,b) * count(a,b), d]) / (count(a,b) + c);
                        count(a,b) = count(a,b) + c;
                    else
                        gridded_mls(a,b,f) = nansum2([gridded_mls(a,b) * count(a,b), d]) / (count(a,b) + c);
                        count(a,b,f) = count(a,b) + c;
                    end
                end
            end
        end
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED PLOTTING FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_dcolumn_percentiles(case1sat, case1crit, case2sat, case2crit, struct_ind, varargin)
        % Plots a timeseries of the difference in total columns between two
        % cases. Requires structures with the information about the columns
        % as the first and third arguments, and logical matrices indicating
        % where lightning NOx is >60% of the emission as the second and
        % fourth. The next to last argument is the index for which field in
        % the structures to plot. The last argument is a boolean
        % determining whether to plot the absolute or relative difference
        % in columns. This is set by the internal function call and does
        % not need to be passed as a user argument to the main function.
        % The difference is defined as first - second or (first -
        % second)/second.
        
        % Parse varargs. If one is the string 'relative', the next one is
        % the relative boolean. There should be 0 or 2 left after that, if
        % 2, then those are the lat and lon limits.
        [relative, ~, lonlim, latlim] = parse_varargs(varargin);
        
        if ~ismember('Columns',fieldnames(case1sat)) %isfield wasn't behaving
            case1sat = integrate_geoschem_profile(case1sat,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2sat))
            case2sat = integrate_geoschem_profile(case2sat,1e-9);
        end
        
        [glon, glat] = geos_chem_centers(size(case1sat(1).dataBlock));
        xx = glon >= min(lonlim) & glon <= max(lonlim);
        yy = glat >= min(latlim) & glat <= max(latlim);
        
        case1col = case1sat(struct_ind).Columns;
        case2col = case2sat(struct_ind).Columns;
        
        case1col(~case1crit) = nan;
        case2col(~case2crit) = nan;
        
        case1col = case1col(xx,yy,:);
        case2col = case2col(xx,yy,:);
        if ~relative
            dcol = case1col - case2col;
        else
            dcol = (case1col - case2col) ./ case2col * 100;
        end
        
        plot_percentile_timeseries(dcol, case1sat(struct_ind).tVec);
        if ~relative
            ylabel('Absolute column differences (molec./cm^2)');
        else
            ylabel('Percent difference in column')
        end
    end

    function plot_level_dconc(case1sat, case1crit, case2sat, case2crit, struct_ind, plevels, varargin)
        % Most inputs the same as plot_dcolumn_percentiles, plevels is a
        % vector of level indicies to make plots for, 1 per level. 
        
        [relative, ndens_bool, lonlim, latlim] = parse_varargs(varargin);
               
        [glon, glat] = geos_chem_centers(size(case1sat(1).dataBlock));
        xx = glon >= min(lonlim) & glon <= max(lonlim);
        yy = glat >= min(latlim) & glat <= max(latlim);
        
        
        titlestr = input('Give a title for these plots with a %s for the pressure level: ', 's');
        caseinds = get_data_inds(case1sat);
        for a=1:numel(plevels)
            case1conc = squeeze(case1sat(struct_ind).dataBlock(xx,yy,plevels(a),2:end))*1e-9;
            case1ndens = squeeze(case1sat(caseinds.ndens).dataBlock(xx,yy,plevels(a),2:end))*1e-6;
            case2conc = squeeze(case2sat(struct_ind).dataBlock(xx,yy,plevels(a),2:end))*1e-9;
            case2ndens = squeeze(case2sat(caseinds.ndens).dataBlock(xx,yy,plevels(a),2:end))*1e-6;
            if ndens_bool
                case1conc = case1conc;% .* case1ndens;
                case2conc = case2conc;% .* case2ndens;
            end
            case1conc(~case1crit(xx,yy,2:end)) = nan;
            case2conc(~case2crit(xx,yy,2:end)) = nan;

            
            if ~relative
                dconc = case1conc - case2conc;
            else
                dconc = (case1conc - case2conc)./case2conc * 100;
            end
            
            % We don't want to include any results from above the
            % tropopause, so set those to NaN.
            tpbool = floor(case1sat(caseinds.tplev).dataBlock(xx,yy,2:end)) <= plevels(a);
            dconc(tpbool) = nan;
            
            plot_percentile_timeseries(dconc, case1sat(struct_ind).tVec(2:end));
            
            if ~relative && ndens_bool
                ylabel('Absolute concentration differences (molec./cm^3)');
            elseif ~relative && ~ndens_bool
                ylabel('Absolute VMR differences (ppp)');
            elseif relative && ndens_bool
                ylabel('Percent difference in num. density')
            elseif relative && ~ndens_bool
                ylabel('Percent difference in VMR');
            end
            
            titlestr = tex_in_printf(titlestr, 2);
            titlestr2 = sprintf(titlestr, '(trop only, Plev = %d, P range=[%.1f, %.1f])');
            sz = size(case1sat(caseinds.psurf).dataBlock(xx,yy,:,:));
            % List the range of pressures found for a given vertical
            % coordinate. However, we need to skip the first day because
            % any ND51 diagnostics that output not at midnight give weird
            % results on the first day
            minpres = min(reshape(case1sat(caseinds.psurf).dataBlock(xx,yy,plevels(a),2:end),sz(1)*sz(2)*(sz(4)-1),1));
            maxpres = max(reshape(case1sat(caseinds.psurf).dataBlock(xx,yy,plevels(a),2:end),sz(1)*sz(2)*(sz(4)-1),1));
            title(sprintf(titlestr2, plevels(a), minpres, maxpres));
        end
    end

    function plot_lnox_col_enhnc_diff(case1lnox, case1nolnox, case1crit, case2lnox, case2nolnox, case2crit, struct_ind, varargin)
        % Plots quantile differences in lightning enhancement over time.
        % case[12][lnox|nolnox] should be the same structures as
        % case[12]sat in the previous two functions, just one with and one
        % without lightning. case[12]crit likewise is the same logical
        % matrix of profiles with significant lightning emissions (e.g. 60%
        % of NO emissions from lightning, BB, anthro due to lightning).
        % struct_ind and relative are the same as before.
        
        [relative, ndens_bool, lonlim, latlim] = parse_varargs(varargin);
        
        if ~ismember('Columns',fieldnames(case1lnox)) %isfield wasn't behaving
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case1lnox\n'); end
            case1lnox = integrate_geoschem_profile(case1lnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case1nolnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case1nolnox\n'); end
            case1nolnox = integrate_geoschem_profile(case1nolnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2lnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case2lnox\n'); end
            case2lnox = integrate_geoschem_profile(case2lnox,1e-9);
        end
        if ~ismember('Columns',fieldnames(case2nolnox))
            if DEBUG_LEVEL > 0; fprintf('Calculating columns for case2nolnox\n'); end
            case2nolnox = integrate_geoschem_profile(case2nolnox,1e-9);
        end
        
        case1lnox_col = case1lnox(struct_ind).Columns;
        case1nolnox_col = case1nolnox(struct_ind).Columns;
        case2lnox_col = case2lnox(struct_ind).Columns;
        case2nolnox_col = case2nolnox(struct_ind).Columns;
        
        case1lnox_col(~case1crit) = nan;
        case1nolnox_col(~case1crit) = nan;
        case2lnox_col(~case2crit) = nan;
        case2nolnox_col(~case2crit) = nan;
        
        dcol = (case1lnox_col - case1nolnox_col) - (case2lnox_col - case2nolnox_col);
        if relative
            dcol = dcol ./ (case2lnox_col - case2nolnox_col) * 100;
        end
        
        plot_percentile_timeseries(dcol, case1lnox(struct_ind).tVec);
        if ~relative
            ylabel('Absolute enhancement differences (molec./cm^2)');
        else
            ylabel('Percent enhancement difference in column')
        end
    end

    function plot_percentile_timeseries(mat_in, tVec)
        % Assume the third dimension is time in the column matrices
        sz = size(mat_in);
        casequant = zeros(sz(3),5);
        
        for a=1:sz(3)
            sub = mat_in(:,:,a);
            casequant(a,:) = quantile(sub(:),[0.05, 0.25, 0.5, 0.75, 0.95]);
        end
        
        figure;
        plot(tVec, casequant);
        datetick('x');
        legend('0.05','0.25','0.50','0.75','0.95');
        set(gca,'ygrid','on');
    end

    function [avg_data, glon, region, p_max, column_bool, fig] = meridional_avg(dataStruct, criteria, varargin)

        
        % Parse the optional arguments
        fig = [];
        column_bool = false;
        timelim_bool = false;
        p_max = [];
        p_mat = [];
        region = '';
        tp_mask = [];
        for a=1:numel(varargin)
            if isnumeric(varargin{a}) && (isscalar(varargin{a}) || (isvector(varargin{a}) && numel(varargin{a}) == 2))
                p_max = varargin{a};
            elseif isstruct(varargin{a}) && ~isempty(regexpi(varargin{a}.fullName, 'PSURF', 'once'))
                p_mat = varargin{a}.dataBlock;
            elseif (isstruct(varargin{a}) && strcmpi(varargin{a}.fullName, 'Tropopause level')) || islogical(varargin{a})
                tp_mask = varargin{a};
            elseif ischar(varargin{a}) && strcmpi('columns',varargin{a})
                column_bool = true;
            elseif ~ischar(varargin{a}) && ishandle(varargin{a}) && strcmpi(get(varargin{a},'type'),'figure')
                fig = varargin{a};
            elseif ischar(varargin{a}) && ismember(varargin{a}, {'na','sa','naf','saf','seas','neur'})
                if ~isempty(region)
                    E.badinput('Multiple regions specified ("%s" and "%s"). Only specify one.', region, varargin{a});
                end
                region = varargin{a};
            elseif ischar(varargin{a}) && strcmpi(varargin{a}, 'timelim')
                timelim_bool = true;
            end
        end
        
        [lonlim, latlim, timemask] = define_regions(region);
        
        % Make sure the input structure is correct (has everything needed)
        req_fields = {'dataUnit','fullName'};
        if column_bool
            req_fields{end+1} = 'Columns';
        else
            req_fields{end+1} = 'dataBlock';
        end
        
        if ~isstruct(dataStruct) || ~isscalar(dataStruct)
            E.badinput('dataStruct must be a scalar structure')
        end
        flds = ismember(req_fields, fieldnames(dataStruct));
        if any(~flds)
            E.badinput('dataStruct must have the fields %s', strjoin(req_fields(~flds),', '));
        end
        
        % Read in Columns or dataBlock. If a level isn't specified and
        % we're reading in dataBlock, the user forgot.
        if column_bool
            data = dataStruct.Columns;
        else
            if isempty(p_max)
                E.badinput('Max pressure must be specified if plotting a concentration')
            elseif isempty(p_mat)
                E.badinput('A matrix of pressures must be give if plotting concentrations')
            elseif isempty(tp_mask) && isscalar(p_max)
                E.badinput('A troposphere mask must be given if plotting a concentration and only a max pressure given')
            elseif isstruct(tp_mask) && ~strcmpi(tp_mask.fullName,'Tropopause level')
                E.badinput('The tropopause mask is a structure but does not seem to refer to the tropopause level')
            end
            
            % Set stratospheric boxes to be NaN
            if isscalar(p_max)
                if isstruct(tp_mask)
                    tp_mask = trop_mask(size(dataStruct.dataBlock), tp_mask);
                end
                dataStruct.dataBlock(~tp_mask) = nan;
            end

            % Average the concentrations over the requested levels. Set
            % boxes below the max pressure to be NaN.
            if ndims(dataStruct.dataBlock) ~= ndims(p_mat) || any(size(dataStruct.dataBlock) ~= size(p_mat))
                E.badinput('Pressure and concentration matrices must be the same size')
            end
            
            if isscalar(p_max)
                pp = p_mat <= p_max;
            else
                pp = p_mat <= max(p_max) & p_mat >= min(p_max);
            end
            dataStruct.dataBlock(~pp) = nan;
            data = squeeze(nanmean(dataStruct.dataBlock,3));
        end
        
        % criteria must be a logical matrix with the same lat/lon/time
        % dimensions as dataBlock or Columns in the input structure.
        if any(size(data) ~= size(criteria))
            E.badinput('Size of input data and criteria do not match')
        end
        
        % Set times outside the time range to NaN so they are not included
        % in the average. Do this only if the 'timelim' argument has been
        % passed and if "data" is 366 elements long in the third dimension,
        % as the time limits assume so (i.e. 366 days of 2012).
        if timelim_bool
            if size(data,3) == 366
                data(:,:,~timemask) = nan;
            else
                E.badinput('To filter by time, the input matrix must represent 366 days')
            end
        end
        
        % Needed for x ticks and to limit to the region
        [glon, glat] = geos_chem_centers(size(data));
        
        % Average the data meridionally (i.e. along latitude). Ignore boxes
        % where the criterion is false. Martin seems to time average first,
        % then average latitudinally, so we'll do the same. Remove boxes
        % outside the area of interest first.
        
        data(~criteria) = nan;
        
        xx = glon >= lonlim(1) & glon <= lonlim(2);
        yy = glat >= latlim(1) & glat <= latlim(2);
        data = data(xx,yy,:);
        
        %avg_data = nanmean(nanmean(data,3),2);
        avg_data = nanmedian(nanmedian(data,3),2);
        glon = glon(xx);
    end

    function plot_meridional_avg(dataStruct, criteria, varargin)
        % Will produce a figure like Fig. 5 in Randall Martin's 2007
        % lightning NOx paper, where the NO2 columns or O3/HNO3
        % concentration are averaged meridionally for places that fit the
        % criteria. Req. 2 arguments: dataStruct must be one structure
        % output from read_geos_output or cleaned up from the ts_satellite
        % python reading code, and the criteria matrix (a logical matrix
        % true wherever the data should be included, i.e lightning >60% of
        % emissions). A troposphere mask - this can either be a logical
        % matrix the same size as the dataBlock or the tropopause level
        % output structure from read_geos_output - is required if doing
        % levels instead of columns.
        %
        % If plotting a level, you'll need to include a matrix of pressures
        % and the maximum pressure to plot.
        %
        % Optional arguments: a handle to a figure will plot on that figure
        % instead of creating a new one. The string 'columns' will cause it
        % to average columns rather than concentrations. Other strings will
        % specify what lon/lat limits to apply: 'na' (north america), 'sa'
        % (south america), 'naf' (north africa), 'saf' (south africa), and
        % 'seas' (southeast asia). The string 'timelim' will limit the
        % region to its three summer months (JJA for N. Hemisphere or DJF
        % for SH).
        
        [avg_data, glon, region, p_max, column_bool, fig] = meridional_avg(dataStruct, criteria, varargin{:});
        
        % Create or select the figure as necessary
        if isempty(fig)
            figure;
        else
            figure(fig);
            hold on
        end
        
        plot(glon', avg_data);
        
        if isempty(fig)
            species_name = regexprep(dataStruct.fullName,'[_^]',' ');
            xlabel('Longitude')
            if column_bool
                ylabel(sprintf('%s columns (molec. cm^{-2})', species_name));
            else
                ylabel(sprintf('%s (%s)', species_name, regexprep(dataStruct.dataUnit, '[_^]', ' ')));
            end
            if ~isempty(region)
                region_title = sprintf(' over %s',upper(region));
            else
                region_title = '';
            end
            if ~isempty(p_max)
                level_title = sprintf(' for P \\\\leq %d', p_max);
            else
                level_title = '';
            end
            species_name = strsplit(species_name);
            title(sprintf('%s Meridonal avg.%s%s', species_name{2}, region_title, level_title));
            set(gca,'fontsize',14);
        end
    end

    function add_sat_mer_avg_to_fig(fig, specie, region, gridded_data, lnox_crit, varargin)
        % Will call sat_meridional_avg to calculate the meridional average
        % of satellite data for the given region then add that to the
        % figure given by the handle fig.  The figure should have a legend
        % already.  If you have pregridded data, pass as the final argument
        % to save time gridding.  This WILL NOT check that you are adding
        % a like species to the figure, so it will happily say put NO2
        % columns on an HNO3 concentration plot or put data from N. Am. on
        % SE Asia. Just warning you...
        %
        % Some optional arguments: the string 'noregrid' will keep the data
        % in its native gridding, and the string 'timelim' will limit the
        % data for each region to that region's three summer months (JJA or
        % DJF).
        
        % Check & parse input
        if ~ishandle(fig) || ~strcmpi(fig.Type, 'figure')
            E.badinput('fig must be a handle to a figure')
        end
        if ~exist('gridded_data','var')
            gridded_data = []; % placeholder, sat_meridional_avg will recognize this as indicating gridding must be done.
        end
        if ~exist('lnox_crit','var')
            lnox_crit = [];
        end
        noregrid_bool = false;
        timelim_bool = false;
        if numel(varargin) > 0
            for a = 1:numel(varargin)
                if ischar(varargin{a}) && strcmpi(varargin{a}, 'noregrid')
                    noregrid_bool = true;
                elseif ischar(varargin{a}) && strcmpi(varargin{a}, 'timelim')
                    timelim_bool = true;
                end
            end
        end
        
        ax = findobj(fig,'Type','axes');
        leg = findobj(fig,'Type','legend');
        
        [sat_lon, sat_avg] = sat_meridional_avg(specie, region, gridded_data, lnox_crit, noregrid_bool, timelim_bool);
        
        line(sat_lon, sat_avg, 'color', 'k', 'linewidth', 1.5, 'linestyle', '--', 'parent', ax);
        %lstrings = leg.String;
        %lstrings{end+1} = 'Satellite';
        %legend(lstrings{:});
    end

    function [ avg_table ] = avg_meridional_avg_table(dataStruct_new, criteria_new, dataStruct_old, criteria_old, varargin)
       % Takes the meridional averages for all 6 regions and then averages
       % over the longitudes, giving a table of absolute and relative
       % average differences in meridional averages by region. See
       % plot_meridional_avg for a description of the inputs, with two
       % differences: 
       %    One, you must input two data structures (the "new" and "old"
       %    cases) as this is computing a difference.
       %
       %    Two, do not input a region as this function will reject it
       %    before calling meridional_avg, as this function will supply the
       %    regions itself.
       
       % This will be used in the table, so it should be a column vector
       regions = {'na';'sa';'naf';'saf';'seas';'neur'};
       
       for a=1:numel(varargin)
           if ischar(varargin{a})
               xx = ismember(varargin{a}, regions);
               
               if any(xx)
                   E.badinput('A region ("%s") was detected in input; remove this and try again. This function will iterate over all regions internally', varargin{a});
               end
           end
       end
       
       AbsoluteDifference = nan(size(regions));
       RelativeDifference = nan(size(regions));
       
       for a=1:numel(regions)
           avg_data_new = meridional_avg(dataStruct_new, criteria_new, varargin{:}, regions{a});
           avg_data_old = meridional_avg(dataStruct_old, criteria_old, varargin{:}, regions{a});
           
           AbsoluteDifference(a) = nanmean(avg_data_new - avg_data_old);
           RelativeDifference(a) = nanmean((avg_data_new - avg_data_old) ./ avg_data_old * 100);
       end
       
       avg_table = table(AbsoluteDifference, RelativeDifference, 'RowNames', regions);
    end

    function compare_cases_during_DC3
        newcase = ask_multichoice('Which case to use as the new case?',{'MPN','PNA','N2O5','HNO3','Updated'},'list',true);
        plotquant = ask_multichoice('Which quantity to plot the difference for?',{'NO2-tVCDs','UT-NOx','UT-HNO3'},'list',true);
        
        base_path = fullfile(HOMEDIR,'Documents','MATLAB','MPN Project','Workspaces','Pickering Parameterization','DailyDC3');
        base_case_path = fullfile(base_path ,'JPLnoMPN-p0');
        switch lower(newcase)
            case 'mpn'
                new_case_path = fullfile(base_path, 'JPLwMPN-p0');
            case 'pna'
                new_case_path = fullfile(base_path, 'JPLnoMPN-PNA-p0');
            case 'n2o5'
                new_case_path = fullfile(base_path, 'JPLnoMPN-N2O5-p0');
            case 'hno3'
                new_case_path = fullfile(base_path, 'HendwoMPN-p0');
            case 'updated'
                new_case_path = fullfile(base_path, 'HendwMPN-PNA-N2O5-p0');
        end
        
        % Load the NOx or HNO3 from both the old and new cases, along with
        % the pressure levels. Cut down to the region and time frame
        % desired, then restrict to 200-350 hPa.
        [xx,yy] = region_xxyy('na');
        dnumlims = [datenum('2012-05-01'), datenum('2012-06-29')];
        if ~isempty(regexpi(plotquant, 'tVCD'))
            % need to load from the OMI overpass instead
            E.notimplemented('tVCD')
        elseif ~isempty(regexpi(plotquant,'NOx'));
            [baseData, baseTVec] = loadNOx(base_case_path);
            [newData, newTVec] = loadNOx(new_case_path);
        elseif ~isempty(regexpi(plotquant,'HNO3'))
            [baseData, baseTVec] = loadGCVar(base_case_path,'.*HNO3.mat');
            [newData, newTVec] = loadGCVar(new_case_path, '.*HNO3.mat');
        end
        
        basePres = loadGCVar(base_case_path, '.*PSURF.mat');
        newPres = loadGCVar(new_case_path, '.*PSURF.mat');
    
        ttnew = newTVec >= dnumlims(1) & newTVec <= dnumlims(2);
        newData = newData(xx,yy,:,ttnew);
        newPres = newPres(xx,yy,:,ttnew);
        newData(newPres < 200 | newPres > 350) = nan;
        ttbase = baseTVec >= dnumlims(1) & baseTVec <= dnumlims(2);
        baseData = baseData(xx,yy,:,ttbase);
        basePres = basePres(xx,yy,:,ttbase);
        baseData(basePres < 200 | basePres > 350) = nan;
        
        del = newData - baseData;
        del = nanmean(nanmean(del,4),3);
        
        [glon, glat] = geos_chem_centers('2x25');
        [GLAT, GLON] = meshgrid(glat(yy),glon(xx));
        
        coastlines = load('coast');
        
        figure; pcolor(GLON,GLAT,del);
        colormap('jet')
        colorbar
        line(coastlines.long, coastlines.lat, 'color', 'k');
        title(sprintf('%s - base: %s',newcase,plotquant))
    end
    function [NOx, tvec] = loadNOx(file_path)
        [NO2, tvec] = loadGCVar(file_path,'.*NO2[-_]IJAVG.mat');
        
        NO = loadGCVar(file_path,'.*NO[-_]IJAVG.mat');
        
        NOx = NO + NO2;
    end
    function [var, tvec] = loadGCVar(directory, pattern)
        F = dir(fullfile(directory,'*.mat'));
        names = glob({F.name},pattern);
        if numel(names) ~= 1
            error('loadGCVars:multiple_files','Multiple files match the pattern %s', pattern);
        end
        filename = fullfile(directory, names{1});
        tmp = load(filename); % there should only be one file matching that pattern
        fn = fieldnames(tmp);
        fn = fn{1};
        var = tmp.(fn).dataBlock;
        tvec = tmp.(fn).tVec;
    end

    function [sa_wt_del, sa_wt_rdel, saf_wt_del, saf_wt_rdel] = check_ak_vcd_domino_weights()
        % Some of the DOMINO AK weights are ~10% lower that the satellite
        % VCD weights. This will give a rough estimate of how much of an
        % effect that would cause.
        
        % Load the AK weight test files that we have
        ak_path = '/Volumes/share2/USERS/LaughnerJ/MPN_Project/AKs/AKs-wttest/DOMINO';
        F_ak = dir(fullfile(ak_path, '*.mat'));
        start_dnum = datenum(regexp(F_ak(1).name,'\d\d\d\d\d\d\d\d','match','once'),'yyyymmdd');
        end_dnum = datenum(regexp(F_ak(end).name,'\d\d\d\d\d\d\d\d','match','once'),'yyyymmdd');
        total_weights = nan(144, 91, numel(F_ak));
        for f=1:numel(F_ak)
            AK = load(fullfile(ak_path, F_ak(f).name));
            for a=1:size(total_weights,1)
                for b=1:size(total_weights,2)
                    wt =  nansum2(AK.binned_weights{a,b});
                    if ~isempty(wt)
                        total_weights(a,b,f) = wt;
                    end
                end
            end
        end
        
        sat_weights = nan(size(total_weights));
        % Now load the corresponding satellite files
        filepat = 'OMI_DOMINO_%04d%02d%02d-newweight.mat';
        sat_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/DOMINO/2.5x2.0-avg-newweight';
        datevec=start_dnum:end_dnum;
        for a=1:numel(datevec)
            filename = sprintf(filepat, year(datevec(a)), month(datevec(a)), day(datevec(a)));
            sat = load(fullfile(sat_path, filename));
            sat_weights(:,:,a) = sat.GC_avg.Weight;
        end
        
        % Since we only have January outputs, just compare the two southern
        % hemisphere regions
        [xx,yy] = region_xxyy('sa');
        sa_wt_del = total_weights(xx,yy,:) - sat_weights(xx,yy,:);
        sa_wt_rdel = reldiff(total_weights(xx,yy,:), sat_weights(xx,yy,:));
        [xx,yy] = region_xxyy('saf');
        saf_wt_del = total_weights(xx,yy,:) - sat_weights(xx,yy,:);
        saf_wt_rdel = reldiff(total_weights(xx,yy,:), sat_weights(xx,yy,:));
    end

    function [omno2_no2, omno2_std, domino_no2, domino_std, omno2_gc, domino_gc, omno2_db, domino_db] = compare_regions_to_sat(varargin)
        % Will generate the (hopefully) final figure for Ben and my paper.
        % The current idea I have is to make two panels, one for OMNO2 one
        % for DOMINO that compare the average satellite columns for the
        % four "good" satellite regions (S.Am., S.Af., N.Af., and SE Asia)
        % to the AK-corrected base, HendwMPN-PNA-N2O5 and HendwMPN-PNA-N2O5
        % +33% cases. There will need to be two plots I think because the
        % AKs are different for each product. So each region will have 8
        % points, four per panel.
        
        plot_ratio_bool = strcmpi(ask_multichoice('Plot as ratios of model to satellite?',{'y','n'}),'y');
        plot_maps = ask_yn('Plot maps too?');
        if plot_maps
            map_product = ask_multichoice('Which product for the maps?', {'sp', 'domino'}, 'list', true);
        end
        
        if ismember('sdom',varargin)
            sdom = true;
            fprintf('Using std. dev. of mean\n');
        else
            sdom = false;
            fprintf('Using just std. dev.\n');
        end

        % Load the production based masks that we'll need to cut down the
        % GC and satellite matices to the lightning heavy observations. We
        % will use the lnox > 60% from the base runs for all runs so that
        % we are comparing the same grid cells in every case.
        M = load(fullfile('/Users/Josh/Documents/MATLAB/','MPN Project','Workspaces','Pickering Parameterization','Daily24hrs','All2012Daily','HendwMPN_PNA_N2O5-Pickp33-ProdFrac.mat'), 'hendwmpn_pna_n2o5_pickp33_lnoxgt60','anthro_lt15e10');
        lnox_mask = M.hendwmpn_pna_n2o5_pickp33_lnoxgt60;
        anthro_mask = M.anthro_lt15e10;
        
        % Then load the satellite data and average it down for each
        % regions.
        
        omno2_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/OMNO2/2.5x2.0-avg-newweight-allvcd';
        domino_path = '/Volumes/share2/USERS/LaughnerJ/DOMINO-OMNO2_comparision/DOMINO/2.5x2.0-avg-newweight-allvcd';
        
        fprintf('Loading OMNO2 data... Please be patient...   ');
        F_omno2 = dir(fullfile(omno2_path,'OMI*.mat'));
        omno2_data = nan(144,91,numel(F_omno2));
        omno2_weight = nan(144,91,numel(F_omno2));
        for a=1:numel(F_omno2)
            O = load(fullfile(omno2_path, F_omno2(a).name));
            if a == 1
                omno2_lon = O.GC_avg.Longitude;
                omno2_lat = O.GC_avg.Latitude;
            end
            omno2_data(:,:,a) = O.GC_avg.ColumnAmountNO2Trop;
            omno2_weight(:,:,a) = O.GC_avg.Weight;
            clear('O')
        end
        fprintf('Done.\n');
        omno2_mask = ~isnan(omno2_data);
        omno2_data(~lnox_mask | ~anthro_mask) = nan;
        
        fprintf('Loading DOMINO data... Please be patient...   ');
        F_domino = dir(fullfile(domino_path, 'OMI*.mat'));
        domino_data = nan(144,91,numel(F_domino));
        domino_weight = nan(144,91,numel(F_domino));
        for a=1:numel(F_domino)
            D = load(fullfile(domino_path, F_domino(a).name));
            if a == 1
                domino_lon = D.GC_avg.Longitude;
                domino_lat = D.GC_avg.Latitude;
            end
            domino_data(:,:,a) = D.GC_avg.TroposphericVerticalColumn;
            domino_weight(:,:,a) = D.GC_avg.Weight;
            clear('D')
        end
        domino_mask = ~isnan(domino_data);
        domino_data(~lnox_mask | ~anthro_mask) = nan;
        fprintf('Done.\n');
        
        fprintf('Dividing sat data into regions...   ');
        regions = {'sa','saf','naf','seas'}; region_xlabels = {'S. Am.','S. Af.','N. Af.', 'SE Asia'};
        %regions = {'sa','saf','naf','seas','natl','spac'}; region_xlabels = {'S. Am.','S. Af.','N. Af.', 'SE Asia', 'N. Atl.', 'S. Pac.'};
        %regions = {'na','sa','saf','naf','neur','seas'}; region_xlabels = {'N. Am.','S. Am.','S. Af.','N. Af.', 'N. Eur.', 'SE Asia'};
        omno2_no2 = make_empty_struct_from_cell(regions);
        omno2_std = make_empty_struct_from_cell(regions);
        domino_no2 = make_empty_struct_from_cell(regions);
        domino_std = make_empty_struct_from_cell(regions);
        
        omno2_db = make_empty_struct_from_cell(regions);
        domino_db = make_empty_struct_from_cell(regions);
        for r=1:numel(regions)
            [xx,yy] = region_xxyy(regions{r});
            [~,~,timeind] = define_regions(regions{r});
            omno2_reg = reshape(omno2_data(xx,yy,timeind),1,[]);
            omno2_wt = reshape(omno2_weight(xx,yy,timeind),1,[]);
            %omno2_no2.(regions{r}) = nanmean(omno2_reg);
            omno2_no2.(regions{r}) = nansum2(omno2_reg .* omno2_wt) / nansum2(omno2_wt);
            notnans = ~isnan(omno2_reg) & ~isnan(omno2_wt);
            %omno2_std.(regions{r}) = nanstd(omno2_reg);
            omno2_std.(regions{r}) = sqrt(var(omno2_reg(notnans),omno2_wt(notnans)));
            
            omno2_db.(regions{r}).columns = omno2_data(xx,yy,timeind);
            omno2_db.(regions{r}).weight = omno2_weight(xx,yy,timeind);
            omno2_db.(regions{r}).avg_vcd = nansum2(omno2_data(xx,yy,timeind) .* omno2_weight(xx,yy,timeind),3) ./ nansum2(omno2_weight(xx,yy,timeind),3);
            omno2_db.(regions{r}).lon = omno2_lon(xx,yy);
            omno2_db.(regions{r}).lat = omno2_lat(xx,yy);
            
            domino_reg = reshape(domino_data(xx,yy,timeind),1,[]);
            domino_wt = reshape(domino_weight(xx,yy,timeind),1,[]);
            %domino_no2.(regions{r}) = nanmean(domino_reg);
            domino_no2.(regions{r}) = nansum2(domino_reg .* domino_wt) / nansum2(domino_wt);
            notnans = ~isnan(domino_reg) & ~isnan(domino_wt);
            %domino_std.(regions{r}) = nanstd(domino_reg);
            domino_std.(regions{r}) = sqrt(var(domino_reg(notnans), domino_wt(notnans)));
            
            domino_db.(regions{r}).columns = domino_data(xx,yy,timeind);
            domino_db.(regions{r}).weight = domino_weight(xx,yy,timeind);
            domino_db.(regions{r}).avg_vcd = nansum2(domino_data(xx,yy,timeind) .* domino_weight(xx,yy,timeind),3) ./ nansum2(domino_weight(xx,yy,timeind),3);
            domino_db.(regions{r}).lon = domino_lon(xx,yy);
            domino_db.(regions{r}).lat = domino_lat(xx,yy);
            
        end
        fprintf('Done.\n');
        
        % Now load each of the 3 cases from GEOS-Chem for each of the
        % products' AKs and get the same regions.
        gc_path = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/DailyOMI/AllAKsDaily-newweight-allvcd-domtrop-ratioscale';
        chems = {'JPLnoMPN','HendwMPN-PNA-N2O5','HendwMPN-PNA-N2O5'};
        omno2_ak_files = {'JPLnoMPN-Pickp0-OMI_NO2_OMNO2_aks_newweight_allvcd_domtrop.mat','HendwMPN-PNA-N2O5-Pickp0-OMI_NO2_OMNO2_aks_newweight_allvcd_domtrop.mat','HendwMPN-PNA-N2O5-Pickp33-OMI_NO2_OMNO2_aks_newweight_allvcd_domtrop.mat'};
        dom_ak_files = {'JPLnoMPN-Pickp0-OMI_NO2_DOMINO_aks_newweight_allvcd_domtrop.mat','HendwMPN-PNA-N2O5-Pickp0-OMI_NO2_DOMINO_aks_newweight_allvcd_domtrop.mat','HendwMPN-PNA-N2O5-Pickp33-OMI_NO2_DOMINO_aks_newweight_allvcd_domtrop.mat'};
        fields = {'Base','Final','Final33'};
        gc_data = make_empty_struct_from_cell(fields);
        omno2_gc = make_empty_struct_from_cell(regions,gc_data);
        domino_gc = make_empty_struct_from_cell(regions,gc_data);
        for a=1:numel(omno2_ak_files)
            fprintf('Loading %s...   ', omno2_ak_files{a});
            O = load(fullfile(gc_path, chems{a}, omno2_ak_files{a}));
            fns = fieldnames(O);
            if numel(fns) > 1; warning('Multiple variables in %s, reading %s only',omno2_ak_files{a},fns{1}); end
            gc = O.(fns{1}).Columns;
            gc(~lnox_mask | ~anthro_mask | ~omno2_mask) = nan;
            for r=1:numel(regions)
                [xx,yy] = region_xxyy(regions{r});
                [~,~,timeind] = define_regions(regions{r});
                gc_reg = reshape(gc(xx,yy,timeind),1,[]);
                gc_wt = reshape(omno2_weight(xx,yy,timeind),1,[]);
                omno2_gc.(regions{r}).(fields{a}).mean = nansum2(gc_reg .* gc_wt)./nansum2(gc_wt);
                notnans = ~isnan(gc_reg) & ~isnan(gc_wt);
                omno2_gc.(regions{r}).(fields{a}).std = sqrt(var(gc_reg(notnans), gc_wt(notnans)));
                omno2_gc.(regions{r}).(fields{a}).vcds = nansum2(gc(xx,yy,timeind) .* omno2_weight(xx,yy,timeind),3) ./ nansum2(omno2_weight(xx,yy,timeind),3);
                %omno2_gc.(regions{r}).(fields{a}).mean = nanmean(reshape(gc(xx,yy,timeind),1,[]));
                %omno2_gc.(regions{r}).(fields{a}).std = nanstd(reshape(gc(xx,yy,timeind),1,[]));
                
                omno2_db.(regions{r}).(fields{a}) = gc(xx,yy,timeind);
            end
            fprintf('Done.\n');
            clear('O','gc')
            
            fprintf('Loading %s...   ', dom_ak_files{a});
            D = load(fullfile(gc_path, chems{a}, dom_ak_files{a}));
            fns = fieldnames(D);
            if numel(fns) > 1; warning('Multiple variables in %s, reading %s only',dom_ak_files{a},fns{1}); end
            gc = D.(fns{1}).Columns;
            gc(~lnox_mask | ~anthro_mask | ~domino_mask) = nan;
            for r=1:numel(regions)
                [xx,yy] = region_xxyy(regions{r});
                [~,~,timeind] = define_regions(regions{r});
                gc_reg = reshape(gc(xx,yy,timeind),1,[]);
                gc_wt = reshape(domino_weight(xx,yy,timeind),1,[]);
                domino_gc.(regions{r}).(fields{a}).mean = nansum2(gc_reg .* gc_wt)./nansum2(gc_wt);
                notnans = ~isnan(gc_reg) & ~isnan(gc_wt);
                domino_gc.(regions{r}).(fields{a}).std = sqrt(var(gc_reg(notnans), gc_wt(notnans)));
                domino_gc.(regions{r}).(fields{a}).vcds = nansum2(gc(xx,yy,timeind) .* domino_weight(xx,yy,timeind),3) ./ nansum2(domino_weight(xx,yy,timeind),3);
                %domino_gc.(regions{r}).(fields{a}).mean = nanmean(reshape(gc(xx,yy,timeind),1,[]));
                %domino_gc.(regions{r}).(fields{a}).std = nanstd(reshape(gc(xx,yy,timeind),1,[]));
                
                domino_db.(regions{r}).(fields{a}) = gc(xx,yy,timeind);
            end
            fprintf('Done.\n')
            clear('D','gc')
        end
        
        if sdom
            sqrtn_o = sqrt(sum(lnox_mask(:) & anthro_mask(:) & omno2_mask(:)));
            sqrtn_d = sqrt(sum(lnox_mask(:) & anthro_mask(:) & domino_mask(:)));
        else
            sqrtn_o = 1;
            sqrtn_d = 1;
        end


        % Maps plotting to actually look at the columns
        if plot_maps
            if strcmpi(map_product,'sp')
                S_sat = omno2_db;
                S_gc = omno2_gc;
            else
                S_sat = domino_db;
                S_gc = domino_gc;
            end
    
            for r=1:numel(regions)
                lon = S_sat.(regions{r}).lon;
                lat = S_sat.(regions{r}).lat;
                figure; 
                subplot(3,3,2);
                pcolor(lon, lat, S_sat.(regions{r}).avg_vcd); colorbar; title('Sat columns');
                for f=1:3
                    subplot(3,3,3+f);
                    pcolor(lon, lat, S_gc.(regions{r}).(fields{f}).vcds); colorbar; title(sprintf('%s columns', fields{f}));
                    subplot(3,3,6+f);
                    pcolor(lon, lat, S_gc.(regions{r}).(fields{f}).vcds - S_sat.(regions{r}).avg_vcd); colorbar; title(sprintf('%s - sat', fields{f}));
                end
            end
        end

        % Finally make the figures. Two panels: OMNO2 on top, DOMINO on
        % bottom.
        
        figure; 
        %subplot(2,1,1)
        
        cols = {'k','b','r'};
        marks = {'d','^','s'};
        
        for r=1:numel(regions)
            x = (r-1)*2+1;
            if ~plot_ratio_bool
                l(1) = line(x-0.5,omno2_no2.(regions{r}),'color','k','marker','o','linewidth',2);
                scatter_errorbars(x-0.5, omno2_no2.(regions{r}), omno2_std.(regions{r})/sqrtn_o, 'color','k','linewidth',2);
                for a=1:numel(fields)
                    dx = a*0.25 - 0.5;
                    l(a+1) = line(x+dx, omno2_gc.(regions{r}).(fields{a}).mean, 'color', cols{a}, 'marker', marks{a}, 'linewidth', 2);
                    scatter_errorbars(x+dx, omno2_gc.(regions{r}).(fields{a}).mean, omno2_gc.(regions{r}).(fields{a}).std/sqrtn_o, 'color', cols{a}, 'linewidth',2);
                end
            else
                % Plot upper and lower limits
%                 ul = 1 + omno2_no2.(regions{r})/(omno2_std.(regions{r})/sqrtn_o);
%                 ll = 1 - omno2_no2.(regions{r})/(omno2_std.(regions{r})/sqrtn_o);
%                 line(x + [-0.5, 0.5], [ul ul], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
%                 line(x + [-0.5, 0.5], [ll ll], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
                % The the ratio of each model case to the satellite
                for a=1:numel(fields)
                    dx = a*0.5 - 1;
                    mod_sat_ratio = omno2_gc.(regions{r}).(fields{a}).mean / omno2_no2.(regions{r});
                    l(a) = line(x+dx, mod_sat_ratio, 'color',cols{a},'marker',marks{1}, 'linestyle','none', 'linewidth', 2,'markersize',16);
                    
                    % Uncertainty propagation for f = m/y => 
                    % s_f^2 = s_m^2/y^2 + s_y^2 * m^2/y^4. Here m is the
                    % model values and y the satellite values.
                    m = omno2_gc.(regions{r}).(fields{a}).mean;
                    s_m = omno2_gc.(regions{r}).(fields{a}).std / sqrtn_o;
                    y = omno2_no2.(regions{r});
                    s_y = omno2_std.(regions{r}) / sqrtn_o;
                    s_f = sqrt((s_m.^2)/(y.^2) + (m.^2)./(y.^4).*s_y.^2);
                    
                    %s_f = (y + s_m)./y; % Ben's simple way (ignores uncertainty in the satellite)
                    %scatter_errorbars(x+dx, mod_sat_ratio, s_f, 'color', cols{a}, 'linewidth', 2);
                    
                    mod_sat_ratio = domino_gc.(regions{r}).(fields{a}).mean / domino_no2.(regions{r});
                    l(a+numel(fields)) = line(x+dx, mod_sat_ratio, 'color',cols{a},'marker',marks{2}, 'linestyle','none', 'linewidth', 2,'markersize',16);
                    
                    % Uncertainty propagation for f = x/y => 
                    % s_f^2 = s_m^2/y^2 + s_y^2 * m^2/y^4. Here m is the
                    % model values and y the satellite values.
                    m = domino_gc.(regions{r}).(fields{a}).mean;
                    s_m = domino_gc.(regions{r}).(fields{a}).std / sqrtn_d;
                    y = domino_no2.(regions{r});
                    s_y = domino_std.(regions{r}) / sqrtn_d;
                    s_f = sqrt((s_m.^2)/(y.^2) + (m.^2)./(y.^4).*s_y.^2);
                    
                    %s_f = (y + s_m)./y; % Ben's simple way (ignores uncertainty in the satellite)
                    %scatter_errorbars(x+dx+0.25, mod_sat_ratio, s_f, 'color', cols{a}, 'linewidth', 2);
                end
            end
        end
        %set(gca,'xtick',[1 3 5 7]);
        set(gca,'xtick',[1 3 5 7 9 11]);
        set(gca,'xticklabels',region_xlabels);
        set(gca,'fontsize',16,'ygrid','on');
        if ~plot_ratio_bool
            legend(l',{'OMNO2','Base case','Final case','Final case +33%'});
        else
            legend(l',{'Base case vs. SP','Final case vs. SP','Final case +33% vs. SP', 'Base case vs. DOM','Final case vs. DOM','Final case +33% vs. DOM'});
        end
        
        return;
        
        subplot(2,1,2);
        for r=1:numel(regions)
            x = (r-1)*2+1;
            if ~plot_ratio_bool
                l(1) = line(x-0.5,domino_no2.(regions{r}),'color','k','marker','o','linewidth',2);
                scatter_errorbars(x-0.5, domino_no2.(regions{r}), domino_std.(regions{r})/sqrtn_d, 'color','k','linewidth',2);
                for a=1:numel(fields)
                    dx = a*0.25 - 0.5;
                    l(a+1) = line(x+dx, domino_gc.(regions{r}).(fields{a}).mean, 'color', cols{a}, 'marker', marks{a}, 'linewidth', 2);
                    scatter_errorbars(x+dx, domino_gc.(regions{r}).(fields{a}).mean, domino_gc.(regions{r}).(fields{a}).std/sqrtn_d, 'color', cols{a}, 'linewidth',2);
                end
            else
                % Plot upper and lower limits
%                 ul = 1 + domino_no2.(regions{r})/(domino_std.(regions{r})/sqrtn_d);
%                 ll = 1 - domino_no2.(regions{r})/(domino_std.(regions{r})/sqrtn_d);
%                 line(x + [-0.5, 0.5], [ul ul], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
%                 line(x + [-0.5, 0.5], [ll ll], 'color', 'k', 'linewidth', 2, 'linestyle', '--');
                % The the ratio of each model case to the satellite
                for a=1:numel(fields)
                    dx = a*0.5 - 1;
                    mod_sat_ratio = domino_gc.(regions{r}).(fields{a}).mean / domino_no2.(regions{r});
                    l(a) = line(x+dx, mod_sat_ratio, 'color',cols{a},'marker',marks{a}, 'linestyle','none', 'linewidth', 2);
                    
                    % Uncertainty propagation for f = x/y => 
                    % s_f^2 = s_m^2/y^2 + s_y^2 * m^2/y^4. Here m is the
                    % model values and y the satellite values.
                    m = domino_gc.(regions{r}).(fields{a}).mean;
                    s_m = domino_gc.(regions{r}).(fields{a}).std / sqrtn_d;
                    y = domino_no2.(regions{r});
                    s_y = domino_std.(regions{r}) / sqrtn_d;
                    s_f = (s_m.^2)/(y.^2) + (m.^2)./(y.^4).*s_y;
                    scatter_errorbars(x+dx, mod_sat_ratio, s_f, 'color', cols{a}, 'linewidth', 2);
                end
            end
        end
        set(gca,'xtick',[1 3 5 7]);
        %set(gca,'xtick',[1 3 5 7 9 11]);
        set(gca,'xticklabels',{'S. Am.','S. Af.','N. Af.', 'SE Asia'});
        %set(gca,'xticklabels',{'N. Am.','S. Am.','S. Af.','N. Af.', 'N. Eur.', 'SE Asia'});
        set(gca,'fontsize',16);
        if ~plot_ratio_bool
            legend(l',{'DOMINO','Base case','Final case','Final case +33%'});
        else
            legend(l',{'Base case','Final case','Final case +33%'});
        end
    end

    function plot_avg_gc_profs()
        regions = ask_multiselect('Which regions to plot?', {'na','sa','naf','saf','neur','seas','natl','spac'},'all',true);
        cases = ask_multiselect('Which cases to plot?', {'base','update','update_p33'}, 'all', true);
        filter_emissions = ask_yn('Filter for lightning >0.6?');

        prof_dir = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/DailyOMI/All2012Daily';
        prof_files = struct('base', fullfile(prof_dir, 'JPLnoMPN-Pickp0-OMI_NO2.mat'),...
                            'update', fullfile(prof_dir, 'HendwMPN-PNA-N2O5-Pickp0-OMI-TIME-SER_NO2.mat'),...
                            'update_p33', fullfile(prof_dir, 'HendwMPN-PNA-N2O5-Pickp33-OMI-TIME-SER_NO2.mat'));

        emis_dir = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/Daily24hrs/All2012Daily';
        emis_files = struct('base', fullfile(emis_dir, 'JPLnoMPN-Pickp0-ProdFrac.mat'),...
                            'update', fullfile(emis_dir, 'HendwMPN_PNA_N2O5-Pickp0-ProdFrac.mat'),...
                            'update_p33', fullfile(emis_dir, 'HendwMPN_PNA_N2O5-Pickp33-ProdFrac.mat'));
        emis_vars = struct('base', 'jplnompn_pickp0_lnoxgt60',...
                           'update', 'hendwmpn_pna_n2o5_pickp0_lnoxgt60',...
                           'update_p33', 'hendwmpn_pna_n2o5_pickp33_lnoxgt60');

        % Load the pressure for vertical coordinate
        load(fullfile(prof_dir, 'JPLnoMPN-Pickp0-OMI_aux.mat'));
        JPLnoMPN_Pickp0_OMI_aux(2:end) = [];
        pres = JPLnoMPN_Pickp0_OMI_aux.dataBlock;

        for c=1:numel(cases)
            P = load(prof_files.(cases{c}));
            fn = fieldnames(P);
            no2 = P.(fn{1}).dataBlock;
            E = load(emis_files.(cases{c}));
            lnox_mask = repmat(E.(emis_vars.(cases{c})), 1, 1, size(no2,3));
            if filter_emissions
                no2(~lnox_mask) = nan;
            end

            for r=1:numel(regions)
                [xx,yy,tt] = misc_gc_plotting_fxns('region_xxyy', regions{r});
                no2slice = permute(no2(xx,yy,:,tt), [3 1 2 4]);
                presslice = permute(pres(xx,yy,:,tt), [3 1 2 4]);
                figure; plot(nanmean(no2slice(:,:),2), nanmean(presslice(:,:),2));
                title(sprintf('%s - %s', cases{c}, regions{r}));
            end
        end
    end

    function domino_sp_scatter()
        quantity = ask_multichoice('What would you like to compare?', {'Strat NO2','Trop AMFs','Trop VCDs','SCDs'}, 'list', true);
        switch lower(quantity)
            case 'strat no2'
                dom_field = 'AssimilatedStratosphericVerticalColumn';
                sp_field = 'ColumnAmountNO2Strat';
            case 'trop amfs'
                dom_field = 'AirMassFactorTropospheric';
                sp_field = 'AmfTrop';
            case 'trop vcds'
                dom_field = 'TroposphericVerticalColumn';
                sp_field = 'ColumnAmountNO2Trop';
            case 'scds'
                dom_field = 'SlantColumnAmountNO2';
                sp_field = 'SlantColumnAmountNO2';
        end
        regions = {'sa','saf','naf','seas'};
        dom_val = cell(size(regions));
        dom_val(:) = {nan(1,366)};
        dom_val_sd = cell(size(regions));
        dom_val_sd(:) = {nan(1,366)};
        sp_val = cell(size(regions));
        sp_val(:) = {nan(1,366)};
        sp_val_sd = cell(size(regions));
        sp_val_sd(:) = {nan(1,366)};
        
        dom_dir = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0/2012';
        dom_pattern = 'OMI-Aura_L2-OMDOMINO_%04dm%02d%02d*';
        sp_dir = '/Volumes/share-sat/SAT/OMI/OMNO2/2012';
        sp_pattern = 'OMI-Aura_L2-OMNO2_%04dm%02d%02d*';
        
        dvec = datenum('2012-01-01'):datenum('2012-12-31');
        [~,~,nh_tt] = define_regions('naf');
        [~,~,sh_tt] = define_regions('saf');
        either_tt = nh_tt | sh_tt;
        for d=1:numel(dvec)
            % Skip if this day isn't used in either hemisphere
            if ~either_tt(d)
                continue
            end
            fprintf('Working on %s\n', datestr(dvec(d)));
            % Find files for this day
            dom_fpat = fullfile(dom_dir,sprintf('%02d',month(dvec(d))), sprintf(dom_pattern, year(dvec(d)), month(dvec(d)), day(dvec(d))));
            dom_files = dir(dom_fpat);
            sp_fpat = fullfile(sp_dir, sprintf('%02d',month(dvec(d))), sprintf(sp_pattern, year(dvec(d)), month(dvec(d)), day(dvec(d))));
            sp_files = dir(sp_fpat);
            
            dom_day_no2 = cell(size(regions));
            sp_day_no2 = cell(size(regions));
            
            for a=1:numel(dom_files)
                % Load arrays
                dom = h5info(fullfile(dom_dir, sprintf('%02d',month(dvec(d))), dom_files(a).name));
                dom_stratno2_a = h5readomi(dom.Filename, h5dsetname(dom,1,2,1,1,dom_field));
                dom_cld_a = h5readomi(dom.Filename, h5dsetname(dom,1,2,1,1,'CloudFraction'));
                dom_alb_a = h5readomi(dom.Filename, h5dsetname(dom,1,2,1,1,'SurfaceAlbedo'));
                dom_flags_a = h5readomi(dom.Filename, h5dsetname(dom,1,2,1,1,'TroposphericColumnFlag'));
                %dom_stratno2_a(dom_cld_a > 0.3 | dom_alb_a > 0.3 | dom_flags_a < 0)=nan;
                dom_stratno2_a(dom_flags_a < 0)=nan;
                dom_lon = h5read(dom.Filename, h5dsetname(dom,1,2,1,2,'Longitude'));
                dom_lat = h5read(dom.Filename, h5dsetname(dom,1,2,1,2,'Latitude'));
                
                % Make a list of pixels that fall in each region
                
                for b=1:numel(regions)
                    [lonlim, latlim, timeind] = define_regions(regions{b});
                    % Skip this region if outside it's time period
                    if ~timeind(d)
                        continue
                    end
                    
                    xx = dom_lon >= lonlim(1) & dom_lon <= lonlim(2) & dom_lat >= latlim(1) & dom_lat <= latlim(2);
                    %fprintf('DOMINO b=%d, sum(xx) = %d\n',b,sum(xx(:)));
                    if sum(xx(:)) > 0
                        dom_day_no2{b} = cat(1, dom_day_no2{b}, dom_stratno2_a(xx));
                    end
                end
            end

            for a=1:numel(sp_files)
                sp = h5info(fullfile(sp_dir, sprintf('%02d',month(dvec(d))), sp_files(a).name));
                sp_stratno2_a = h5readomi(sp.Filename, h5dsetname(sp,1,2,1,1,sp_field));
                sp_cld_a = h5readomi(sp.Filename, h5dsetname(sp,1,2,1,1,'CloudFraction'));
                sp_alb_a = h5readomi(sp.Filename, h5dsetname(sp,1,2,1,1,'TerrainReflectivity'));
                sp_xtrack_a = h5readomi(sp.Filename, h5dsetname(sp,1,2,1,1,'XTrackQualityFlags'));
                sp_vcdflags_a = h5readomi(sp.Filename, h5dsetname(sp,1,2,1,1,'VcdQualityFlags'));
                %sp_stratno2_a(sp_cld_a > 0.3 | sp_alb_a > 0.3 | mod(sp_xtrack_a, 2) ~= 0 | mod(sp_vcdflags_a,2) ~= 0) = nan;
                sp_stratno2_a(mod(sp_xtrack_a, 2) ~= 0) = nan;
                sp_lon = h5read(sp.Filename, h5dsetname(sp,1,2,1,2,'Longitude'));
                sp_lat = h5read(sp.Filename, h5dsetname(sp,1,2,1,2,'Latitude'));
                
                for b=1:numel(regions)
                    [lonlim, latlim, timeind] = define_regions(regions{b});
                    % Skip this region if outside it's time period
                    if ~timeind(d)
                        continue
                    end
                    
                    xx = sp_lon >= lonlim(1) & sp_lon <= lonlim(2) & sp_lat >= latlim(1) & sp_lat <= latlim(2);
                    %fprintf('SP b=%d, sum(xx) = %d\n',b,sum(xx(:)));
                    if sum(xx(:)) > 0
                        sp_day_no2{b} = cat(1, sp_day_no2{b}, sp_stratno2_a(xx));
                    end
                end
            end
            
            for b=1:numel(regions)
                if ~isempty(dom_day_no2)
                    dom_val{b}(d) = nanmean(dom_day_no2{b});
                    dom_val_sd{b}(d) = nanstd(dom_day_no2{b});
                end
                if ~isempty(sp_day_no2)
                    sp_val{b}(d) = nanmean(sp_day_no2{b});
                    sp_val_sd{b}(d) = nanstd(sp_day_no2{b});
                end
            end
        end
        
        for b=1:numel(regions)
            figure;
            title(regions{b});
            line(dom_val{b}, sp_val{b}, 'color', 'k', 'marker', 'o', 'linestyle', 'none');
            scatter_errorbars(dom_val{b}, sp_val{b}, dom_val_sd{b}, 'direction', 'x');
            scatter_errorbars(dom_val{b}, sp_val{b}, sp_val_sd{b}, 'direction', 'y');
            xlabel(sprintf('DOMINO %s', dom_field))
            ylabel(sprintf('SP %s', sp_field))
        end
    end

    function compare_ak_vec()
        %fprintf('Loading file %d of %d\n',a,numel(F));
        region = ask_multichoice('Which region to look at?', {'sa','saf','naf','seas'}, 'list', true);
        ak_bool = ask_multichoice('Plot as scattering weights or AKs?', {'scattering weights', 'AKs'}, 'list', true);
        ak_bool = strcmpi(ak_bool, 'AKs');
        reldif_bool = ask_multichoice('Plot as relative difference?', {'y','n'});
        reldif_bool = strcmpi(reldif_bool,'y');
        
        dvec = datenum('2012-01-01'):datenum('2012-12-31');
        [lonlim,latlim,tt] = define_regions(region);
        dvec = dvec(tt);
        savenum=1;
        % Pick a random date and load the files from that date
        while true
            r = randi(numel(dvec),1);
            [dpath, dfiles] = domino_files(dvec(r));
            [dom_aks, dom_plevs, dom_lon, dom_lat] = load_domino_aks(dpath, dfiles, lonlim, latlim, ak_bool);
            [spath, sfiles] = sp_files(dvec(r));
            [sp_aks, sp_plevs, sp_lon, sp_lat] = load_sp_aks(spath, sfiles, lonlim, latlim, ak_bool);
            dom_nans = squeeze(all(isnan(dom_aks),1));
            sp_nans = squeeze(all(isnan(sp_aks),1));
            if all(dom_nans(:)) || all(sp_nans(:)) || numel(sp_lon) ~= numel(dom_lon)
                continue
            end
            dom_notin = dom_lon < lonlim(1) | dom_lon > lonlim(2) | dom_lat < latlim(1) | dom_lat > latlim(2);
            sp_notin = sp_lon < lonlim(1) | sp_lon > lonlim(2) | sp_lat < latlim(1) | sp_lat > latlim(2);
            % We only want to look at pixels that are valid in both
            % products
            dom_aks(:, dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_plevs(:,dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_lon(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_lat(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_aks(:, dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_lon(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_lat(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            
            % Now choose a random pixel and plot it
            while true
                p = randi(numel(sp_lon),1);
                fprintf('DOMINO pixel: %f, %f     SP pixel: %f, %f\n', dom_lon(p), dom_lat(p), sp_lon(p), sp_lat(p));
                fig=figure;
                if reldif_bool
                    dom_at_sp_plevs = interp1(dom_plevs(:,p), dom_aks(:,p), sp_plevs);
                    sp_at_dom_plevs = interp1(sp_plevs, sp_aks(:,p), dom_plevs(:,p));
                    plot(reldiff(dom_aks(:,p), sp_at_dom_plevs)*100, dom_plevs(:,p));
                    hold on
                    plot(reldiff(dom_at_sp_plevs, sp_aks(:,p))*100, sp_plevs);
                    legend('At DOMINO pres','At SP pres');
                    xlabel('% difference (DOMINO - SP)');
                else
                    plot(dom_aks(:,p), dom_plevs(:,p));
                    hold on
                    plot(sp_aks(:,p), sp_plevs);
                    legend('DOMINO', 'SP')
                end
                drawnow
                
                usersel = ask_multichoice('Enter n for a new profile, d for a new day, s to save the figure, or q to quit',{'n','d','s','q'});
                if strcmpi(usersel,'s')
                    if ak_bool
                        ak_str = 'AKs';
                    else
                        ak_str = 'SWs';
                    end
                    savename = sprintf('%s-%s_%04d.png',region,ak_str,savenum);
                    saveas(fig, savename);
                    savenum = savenum+1;
                    
                    fprintf('Figure saved as %s\n',savename);
                    usersel = ask_multichoice('Enter n for a new profile, d for a new day, or q to quit',{'n','d','q'});
                end
                close(fig);
                switch usersel
                    case 'd'
                        break
                    case 'q'
                        return
                end
                
            end
        end
        
    end

    function compare_raw_aks()
        region = ask_multichoice('Which region to look at?', {'sa','saf','naf','seas'}, 'list', true);
        ak_or_sw = ask_multichoice('Plot as scattering weights or AKs?', {'Scattering Weight', 'AKs'}, 'list', true);
        ak_bool = strcmpi(ak_or_sw, 'AKs');
        
        dvec = datenum('2012-01-01'):datenum('2012-12-31');
        [lonlim,latlim,tt] = define_regions(region);
        dvec = dvec(tt);
        
        figure; 
        title(sprintf('DOMINO ensemble - %s', upper(region)));
        
        day_incr = 15;
        for a=1:day_incr:numel(dvec)
            fprintf('Now on day %d of %d\n', a, floor(numel(dvec)/day_incr));
            [dpath, dfiles] = domino_files(dvec(a));
            [dom_aks, dom_plevs, dom_lon, dom_lat] = load_domino_aks(dpath, dfiles, lonlim, latlim, ak_bool);
            [spath, sfiles] = sp_files(dvec(a));
            [sp_aks, sp_plevs, sp_lon, sp_lat] = load_sp_aks(spath, sfiles, lonlim, latlim, ak_bool);
            dom_nans = squeeze(all(isnan(dom_aks),1));
            sp_nans = squeeze(all(isnan(sp_aks),1));
            if all(dom_nans(:)) || all(sp_nans(:)) || numel(sp_lon) ~= numel(dom_lon)
                continue
            end
            dom_notin = dom_lon < lonlim(1) | dom_lon > lonlim(2) | dom_lat < latlim(1) | dom_lat > latlim(2);
            sp_notin = sp_lon < lonlim(1) | sp_lon > lonlim(2) | sp_lat < latlim(1) | sp_lat > latlim(2);
            
            % We only want to look at pixels that are valid in both
            % products
            dom_aks(:, dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_plevs(:,dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_lon(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            dom_lat(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_aks(:, dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_lon(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            sp_lat(dom_nans | sp_nans | dom_notin | sp_notin) = [];
            
            % For both, we'll average them over time, but for DOMINO we
            % also want to plot the ensemble to prove that the
            % interpolation is behaving correctly (and hope that this
            % doesn't crash matlab b/c of too many lines)
            if a == 1
                sp_mean = nan(numel(sp_plevs),1);
                sp_count = 0;
                dom_mean = nan(numel(sp_plevs),1);
                dom_count = 0;
            end
            
            for b=1:size(dom_aks,2)
                line(dom_aks(:,b), dom_plevs(:,b), 'color', [0.5 0.5 0.5]);
                % We'll use the SP pressure levels to interpolate to b/c
                % they are nice and consistent
                interp_dom = interp1(dom_plevs(:,b), dom_aks(:,b), sp_plevs(:));
                dom_mean = nansum2(cat(2, dom_mean, interp_dom),2);
                dom_count = dom_count + 1;
            end
            
            for b=1:size(sp_aks,2)
                sp_mean = nansum2(cat(2, sp_mean, sp_aks(:,b)),2);
                sp_count = sp_count + 1;
            end
            
        end
        
        dom_mean = dom_mean / dom_count;
        sp_mean = sp_mean / sp_count;
        line(dom_mean, sp_plevs, 'color', 'b', 'linewidth', 2);
        set(gca,'fontsize',16,'ydir','reverse')
        xlabel(ak_or_sw); ylabel('Pressure (hPa)');
        
        figure;
        l = gobjects(2,1);
        l(1) = line(sp_mean, sp_plevs, 'color', 'r', 'linewidth', 2);
        l(2) = line(dom_mean, sp_plevs, 'color', 'b', 'linewidth', 2);
        legend(l, {sprintf('NASA SP (%d pixels)', sp_count), sprintf('DOMINO (%d pixels)', dom_count)}, 'location', 'NorthWest');
        xlabel(ak_or_sw); ylabel('Pressure (hPa)');
        set(gca,'fontsize',16,'ydir','reverse');
        ylim([0 1013]);
        title(sprintf('Avg. %s: %s', ak_or_sw, upper(region)));
    end

    function compute_emissions()
        emis_dir = '/Users/Josh/Documents/MATLAB/MPN Project/Workspaces/Pickering Parameterization/Daily24hrs/All2012Daily';
        area_file = fullfile(emis_dir, 'Area.mat');
        emis_files.base = fullfile(emis_dir, 'JPLnoMPN-Pickp0-Prod.mat');
        emis_files.updated = fullfile(emis_dir, 'HendwMPN-Pickp0-Prod.mat');
        emis_files.updated33 = fullfile(emis_dir, 'HendwMPN-Pickp33-Prod.mat');
        
        fns = fieldnames(emis_files);
        Afile = load(area_file);
        area_var = fieldnames(Afile);
        area = Afile.(area_var{1});
        
        for a=1:numel(fns);
            fprintf('Loading %s case...\n', fns{a});
            D = load(emis_files.(fns{a}));
            emis_var = fieldnames(D);
            emis = D.(emis_var{1});
            
            [total_prod, db] = geos_integrate_prod(area, emis(6));
            
            if strcmpi(fns{a}, 'base')
                base_prod = db.prod_grid;
            elseif strcmpi(fns{a}, 'updated33');
                new_prod = db.prod_grid;
            end
        
            % Convert from moles to Tg (* molec. mass N, /1e12)
            fprintf('   %s case total lightning NO production = %f Tg N yr^{-1}\n', fns{a}, total_prod * 14e-12);
        end
        
        [~,glat] = geos_chem_centers('2x25');
        n_hemisphere = glat >= 0;
        hybrid_prod = nan(size(base_prod));
        hybrid_prod(:,n_hemisphere) = new_prod(:,n_hemisphere);
        hybrid_prod(:,~n_hemisphere) = base_prod(:,~n_hemisphere);
        
        fprintf('\n    Production assuming +33%% only in northern hemisphere: %f Tg N yr^{-1}\n', nansum(hybrid_prod(:))*14e-12);
    end

%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%

    function [aks, plevs, lon, lat] = load_sp_aks(filepath, files, lonlim, latlim, ak_bool)
        aks = [];
        lon = [];
        lat = [];
        for a = 1:numel(files)
            sp = h5info(fullfile(filepath, files(a).name));
            
            % SP pressure levels are always the same
            if a == 1
                plevs = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'ScatteringWtPressure')));
            end
            
            this_lon = h5read(sp.Filename, h5dsetname(sp,1,2,1,2,'Longitude'));
            this_lat = h5read(sp.Filename, h5dsetname(sp,1,2,1,2,'Latitude'));
            xx = this_lon >= lonlim(1) & this_lon <= lonlim(2) & this_lat >= latlim(1) & this_lat <= latlim(2);
            if sum(xx) < 1
                continue
            end
            this_lon(~xx) = [];
            this_lat(~xx) = [];
            
            omi.sw = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'ScatteringWeight')));
            omi.sw(omi.sw < 1e-29) = nan;
            omi.amf = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'AmfTrop')));
            omi.amf(omi.amf < 1e-29) = nan;
            
            % These are fields needed to reject pixels
            omi.CloudFraction = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'CloudFraction')))*1e-3;
            omi.CloudFraction(omi.CloudFraction<0) = nan;
            omi.ColumnAmountNO2Trop = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'ColumnAmountNO2Trop')));
            omi.ColumnAmountNO2Trop(omi.ColumnAmountNO2Trop<-1e29) = nan;
            omi.vcdQualityFlags = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1, 'VcdQualityFlags')));
            omi.vcdQualityFlags(omi.vcdQualityFlags==65535)=nan;
            omi.XTrackQualityFlags = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'XTrackQualityFlags')));
            omi.XTrackQualityFlags(omi.XTrackQualityFlags==255) = nan;
            omi.TerrainReflectivity = double(h5read(sp.Filename, h5dsetname(sp,1,2,1,1,'TerrainReflectivity')));
            omi.TerrainReflectivity(omi.TerrainReflectivity==-32767) = nan;
            omi.TerrainReflectivity = omi.TerrainReflectivity * 1e-3;
            omi.Areaweight = ones(size(omi.amf));
            
            omi = omi_sp_pixel_reject(omi,'geo',0.3,'XTrackFlags');
            
            rejects = omi.Areaweight == 0 | isnan(omi.CloudFraction) | isnan(omi.ColumnAmountNO2Trop)...
                | isnan(omi.TerrainReflectivity) | omi.TerrainReflectivity > 0.3;
            omi.sw(:,rejects) = nan;
            omi.amf(rejects) = nan;
            
            omi.sw(:,~xx) = [];
            omi.amf(~xx) = [];
            
            if ak_bool
                this_aks = nan(size(omi.sw));
                for b=1:numel(omi.amf)
                    this_aks(:,b) = omi.sw(:,b) ./ omi.amf(b);
                end
            else
                this_aks = omi.sw;
            end
            
            
            aks = cat(2, aks, this_aks);
            lon = cat(2, lon, this_lon);
            lat = cat(2, lat, this_lat);
        end
    end

    function [aks, plevs, lon, lat] = load_domino_aks(filepath, files, lonlim, latlim, ak_bool)
        aks = [];
        plevs = [];
        lon = [];
        lat = [];
        for a = 1:numel(files)
            % Load in data and remove fill values
            hi = h5info(fullfile(filepath, files(a).name));
            
            this_lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
            this_lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
            xx = this_lon >= lonlim(1) & this_lon <= lonlim(2) & this_lat >= latlim(1) & this_lat <= latlim(2);
            if sum(xx) < 1
                continue
            end
            this_lon(~xx) = [];
            this_lat(~xx) = [];
            
            omi.aks = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AveragingKernel')))*0.001;
            omi.aks(omi.aks<-30) = nan;
            omi.amf = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AirMassFactor')));
            omi.amf(omi.amf<-1e29) = nan;
            omi.amftrop = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AirMassFactorTropospheric')));
            omi.amftrop(omi.amftrop<-1e29) = nan;
            omi.TM4PresA = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4PressurelevelA')))*0.01; % in Pa, want hPa.
            omi.TM4PresA(omi.TM4PresA<-1e30) = nan;
            omi.TM4PresB = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4PressurelevelB')));
            omi.TM4PresB(omi.TM4PresB<-1e30) = nan;
            omi.TM4SurfP = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4SurfacePressure')));
            omi.TM4SurfP(omi.TM4SurfP<-1e30) = nan;
            % Rearrange so that the vertical or corner coordinate is first
            omi.aks = permute(omi.aks, [3 1 2]);
            
            % Calculate the pressure levels for each pixel. 
            omi.PresLevs = nan(size(omi.aks));
            for x=1:size(omi.TM4SurfP,1)
                for y=1:size(omi.TM4SurfP,2)
                    omi.PresLevs(:,x,y) = omi.TM4PresA + omi.TM4SurfP(x,y) .* omi.TM4PresB;
                end
            end
            
            % These are fields needed to reject pixels
            cldScaleFactor = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,1,'CloudFraction'),'ScaleFactor'));
            omi.CloudFraction = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'CloudFraction')))*cldScaleFactor;
            omi.CloudFraction(omi.CloudFraction<0) = nan;
            vcdScaleFactor = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,1,'TroposphericVerticalColumn'),'ScaleFactor'));
            omi.ColumnAmountNO2Trop = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TroposphericVerticalColumn')))*vcdScaleFactor;
            omi.TropColumnFlag = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TroposphericColumnFlag')));
            albScaleFactor = double(h5readatt(hi.Filename, h5dsetname(hi,1,2,1,1,'SurfaceAlbedo'),'ScaleFactor'));
            omi.Albedo = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'SurfaceAlbedo')))*albScaleFactor;
            
            bad_vals = omi.CloudFraction > 0.3 | omi.ColumnAmountNO2Trop < 0 | omi.TropColumnFlag < 0 | omi.Albedo > 0.3;
            omi.aks(:,bad_vals) = nan;
            omi.PresLevs(:,bad_vals) = nan;
            omi.amf(bad_vals) = nan;
            if ak_bool
                for b=1:numel(omi.amf)
                    omi.aks(:,b) = omi.aks(:,b) * omi.amf(b) / omi.amftrop(b);
                end
                aks = cat(2, aks, omi.aks(:,xx));
            else
                omi.aks(:,~xx) = [];
                omi.amf(~xx) = [];
                this_sw = nan(size(omi.aks));
                for b = 1:numel(omi.amf)
                    this_sw(:,b) = omi.aks(:,b) .* omi.amf(b);
                end
                aks = cat(2, aks, this_sw);
            end
            plevs = cat(2, plevs, omi.PresLevs(:,xx));

            lon = cat(2, lon, this_lon);
            lat = cat(2, lat, this_lat);
        end
    end

    function [filepath, files] = domino_files(date_in)
        dom_dir = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0/2012';
        dom_pattern = 'OMI-Aura_L2-OMDOMINO_%04dm%02d%02d*';
        filepath = fullfile(dom_dir, sprintf('%02d', month(date_in)));
        filename = fullfile(filepath, sprintf(dom_pattern, year(date_in), month(date_in), day(date_in)));
        files = dir(filename);
    end

    function [filepath, files] = sp_files(date_in)
        sp_dir = '/Volumes/share-sat/SAT/OMI/OMNO2/version_3_2_1/2012';
        sp_pattern = 'OMI-Aura_L2-OMNO2_%04dm%02d%02d*';
        filepath = fullfile(sp_dir, sprintf('%02d', month(date_in)));
        filename = fullfile(filepath, sprintf(sp_pattern, year(date_in), month(date_in), day(date_in)));
        files = dir(filename);
    end
end

