function [ no2_bins, temp_bins, wrf_file, SurfPres, SurfPres_WRF, TropoPres, tropopause_interp_flag, pres_mode, temp_mode ] = rProfile_WRF( date_in, profile_mode, region, loncorns, latcorns, omi_time, globe_elevation, pressures, varargin )

%RPROFILE_WRF Reads WRF NO2 profiles and averages them to pixels.
%   This function is the successor to rProfile_US and serves essentially
%   the same purpose - read in WRF-Chem NO2 profiles to use as the a priori
%   for the BEHR OMI NO2 retrieval algorithm. The key difference is that
%   this no longer assumes every NO2 profile from WRF is on the same
%   pressure grid. Therefore, it recognizes that the profiles will need to
%   be interpolated to a uniform pressure grid. It will therefore
%   interpolate all NO2 profiles to be averaged together for a pixel to the
%   same pressures before averaging, and will ensure that one bin below the
%   surface is calculated; that way omiAmfAK2 and integPr2 will be able to
%   calculate the best approximation of the surface concentration for
%   integration.
%
%   Further this will read the profiles directly from the netCDF files
%   containing the WRF output, processed by the slurmrun_wrf_output.sh
%   utility. It will look for WRF_BEHR files. These should have variables
%   no2, XLAT, XLONG, and pres. (pres is calculated as the sum of the
%   original WRF variables P and PB, perturbation and base pressure. See
%   calculated_quantities.nco and read_wrf_output.sh in the WRF_Utils
%   folder.)
%
%   This can also use different profiles: monthly, daily, or hourly. These
%   should be saved under folders so named in the main BEHR WRF output
%   folder.  The averaging will have been done in the processing by
%   slurmrun_wrf_output.sh, so the only difference to the execution of this
%   function will be that if hourly is chosen, it will need to select the
%   hour closest to 1400 local standard time for each pixel - this tries to
%   get as close to overpass as possible, although until the fraction of
%   the orbit passed is imported as well, an exact calculation of overpass
%   time will not be possible.
%
%   The inputs to this function are:
%       date_in: The date being processed, used to find the right file.
%
%       profile_mode: Should be 'daily' or 'monthly'. Chooses which WRF
%       profiles are used, the daily or monthly outputs.
%
%       loncorns & latcorns: arrays containing the longitude and latitude
%       corners of the pixels. Can be 2- or 3-D, so long as the first
%       dimension has size 4 (i.e. loncorn(:,a) would give all 4 corners
%       for pixel a).
%
%       omi_time: the starting time of the OMI swath in TAI93 (the time
%       format given in the OMNO2 files). Used to match up daily profiles
%       to the OMI swath.
%
%       surfPres: a 1- or 2-D array containing the GLOBE surface pressures
%       for the pixels. This will be used to ensure enough bins are
%       calculated that omiAmfAK2 and integPr2 can interpolate to the
%       surface pressure.
%
%       pressure: the vector of pressures that the NO2 profiles should be
%       interpolated to.
%
%       wrf_output_path: optional, if provided, overrides which directory
%       this function will look for WRF output in. If passed an empty
%       string, the proper WRF directory will be found, just as if this
%       input was omitted.
%
%
%   Additional parameter inputs:
%
%       err_missing_att: controls whether an error is thrown if attributes
%       in WRF files cannot be found (true, which is the default) or
%       default values are assumed (false). Use with caution, as when set
%       to "false" there will be no error checking of the units in WRF
%       files.
%
%   Josh Laughner <joshlaugh5@gmail.com> 22 Jul 2015

E = JLLErrors;

parser = inputParser;
parser.addOptional('wrf_output_path', '', @ischar);
parser.addParameter('err_missing_att', true);
parser.addParameter('clip_at_int_limits', true);
parser.addParameter('DEBUG_LEVEL', 1);

parser.parse(varargin{:});
pout = parser.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
error_if_missing_attr = pout.err_missing_att;
wrf_output_path = pout.wrf_output_path;
clip_at_int_limits = pout.clip_at_int_limits;

if ~isnumeric(DEBUG_LEVEL) || ~isscalar(DEBUG_LEVEL)
    E.badinput('DEBUG_LEVEL must be a scalar number')
end
if ~islogical(error_if_missing_attr) || ~isscalar(error_if_missing_attr)
    E.badinput('The parameter "err_missing_att" must be a scalar logical');
end

% Defining custom errors
% Error for failing to find netCDF variable that is more descriptive about
% what likely went wrong. The error will take two additional arguments when
% called: the variable name and the file name.
E.addCustomError('ncvar_not_found','The variable %s is not defined in the file %s. Likely this file was not processed with (slurm)run_wrf_output.sh, or the processing failed before writing the calculated quantites.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CHECK DEPENDENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = GitChecker;
% Require that WRF_Utils has the version of find_wrf_tropopause that has
% been updated to identify jumps in the tropopause pressure and the general
% utils have been updated to include the floodfill algorithm needed by 
% find_wrf_tropopause.
G.addReqCommits(behr_paths.wrf_utils, '9272d08');
G.addReqCommits(behr_paths.utils, '51a0869');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(loncorns,1) ~= 4 || size(latcorns,1) ~= 4
    E.badinput('loncorns and latcorns must have the corners along the first dimension (i.e. size(loncorns,1) == 4')
elseif ndims(loncorns) ~= ndims(latcorns) || ~all(size(loncorns) == size(latcorns))
    E.badinput('loncorns and latcorns must have the same dimensions')
end

sz_corners = size(loncorns);
sz_surfPres = size(globe_elevation);

% Check that the corner arrays and the surfPres array represent the same
% number of pixels
if ndims(loncorns)-1 ~= ndims(globe_elevation) || ~all(sz_corners(2:end) == sz_surfPres(1:end))
    E.badinput('The size of the surfPres array must be the same as the corner arrays without their first dimension (size(surfPres,1) == size(loncorns,2, etc)')
end

% pressures should just be a vector, monotonically decreasing
if ~isvector(pressures) || any(diff(pressures)>0)
    E.badinput('pressures must be a monotonically decreasing vector')
end

% Make sure that the input for the date can be understood by MATLAB as a
% date
try
    date_num_in = datenum(date_in);
catch err
    if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
        E.badinput('%s could not be understood by MATLAB as a date')
    else
        rethrow(err);
    end
end

% Get the WRF output path - this function will itself throw an error if the
% profile mode is wrong or the path does not exist.
if isempty(wrf_output_path)
    wrf_output_path = find_wrf_path(region, profile_mode, date_in);
else
    if ~ischar(wrf_output_path)
        E.badinput('WRF_OUTPUT_PATH must be a string');
    elseif ~exist(wrf_output_path, 'dir')
        E.dir_dne(wrf_output_path);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD netCDF and READ VARIABLES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[wrf_no2, wrf_temp, wrf_pres, wrf_lon, wrf_lat, wrf_file, wrf_tropopres, wrf_surf_elev, pres_mode, temp_mode,tropopause_interp_indx] = load_wrf_vars();
num_profs = numel(wrf_lon);
prof_length = size(wrf_no2,3);

num_pix = numel(globe_elevation);
no2_bins = nan(length(pressures), size(globe_elevation,1), size(globe_elevation,2));
temp_bins = nan(length(pressures), size(globe_elevation,1), size(globe_elevation,2));
SurfPres = nan(size(globe_elevation));
SurfPres_WRF = nan(size(globe_elevation));
TropoPres = nan(size(globe_elevation));
tropopause_interp_flag = false(size(globe_elevation));

if any(size(wrf_lon) < 2) || any(size(wrf_lat) < 2)
    error('rProfile_WRF:wrf_dim','wrf_lon and wrf_lat should be 2D');
end
wrf_lon_bnds = [wrf_lon(1,1), wrf_lon(1,end), wrf_lon(end,end), wrf_lon(end,1)];
wrf_lat_bnds = [wrf_lat(1,1), wrf_lat(1,end), wrf_lat(end,end), wrf_lat(end,1)];


% If the WRF profiles are spaced at intervals larger than the smallest dimension
% of OMI pixels, interpolate instead of averaging b/c we will likely have at least
% some pixels with no profiles within them.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BIN PROFILES TO PIXELS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape the NO2 profiles and pressures such that the profiles are along
% the first dimension

% Reorder dimensions. This will make the perm_vec be [3, 1, 2, (4:end)]
perm_vec = 1:ndims(wrf_no2);
perm_vec(3) = [];
perm_vec = [3, perm_vec];

wrf_no2 = permute(wrf_no2, perm_vec);
wrf_temp = permute(wrf_temp, perm_vec);
wrf_pres = permute(wrf_pres, perm_vec);

wrf_no2 = reshape(wrf_no2, prof_length, num_profs);
wrf_temp = reshape(wrf_temp, prof_length, num_profs);
wrf_pres = reshape(wrf_pres, prof_length, num_profs);
wrf_lon = reshape(wrf_lon, 1, num_profs);
wrf_lat = reshape(wrf_lat, 1, num_profs);
wrf_tropopres = reshape(wrf_tropopres, 1, num_profs);

lons = squeeze(nanmean(loncorns,1));
lats = squeeze(nanmean(latcorns,1));
for p=1:num_pix
    if ~inpolygon(lons(p), lats(p), wrf_lon_bnds, wrf_lat_bnds)
        continue
    end

    [no2_bins(:,p), temp_bins(:,p), SurfPres(p), SurfPres_WRF(p), TropoPres(p), tropopause_interp_flag(p)] = avg_apriori();    
end


    function [no2_vec, temp_vec, pSurf, pSurf_WRF, pTropo, trop_interp_flag] = avg_apriori()
        xall = loncorns(:,p);
        xall(5) = xall(1);
        
        yall = latcorns(:,p);
        yall(5) = yall(1);
        
        % Try to speed this up by removing profiles outside a rectangle around
        % the pixel first, then deal with the fact that the pixel is angled
        % relative to lat/lon.
        
        xx = wrf_lon < max(xall) & wrf_lon > min(xall) & wrf_lat < max(yall) & wrf_lat > min(yall);
        tmp_no2 = wrf_no2(:,xx);
        tmp_temp = wrf_temp(:,xx);
        tmp_pres = wrf_pres(:,xx);
        tmp_lon = wrf_lon(xx);
        tmp_lat = wrf_lat(xx);
        tmp_elev = wrf_surf_elev(xx);
        tmp_pTropo = wrf_tropopres(xx);
        
        if any(xx(tropopause_interp_indx))
            trop_interp_flag = true;
        else
            trop_interp_flag = false;
        end
        
        yy = inpolygon(tmp_lon, tmp_lat, xall, yall);
        
        if sum(yy) < 1
            %E.callError('no_prof','WRF Profile not found for pixel near %.1, %.1f',mean(xall),mean(yall));
            no2_vec = nan(length(pressures),1);
            temp_vec = nan(length(pressures),1);
            pTropo = nan;
            pSurf = nan;
            pSurf_WRF = nan;
            return
        end
        
        tmp_no2(:,~yy) = [];
        tmp_temp(:,~yy) = [];
        tmp_pres(:,~yy) = [];
        tmp_elev(~yy) = [];
        tmp_pTropo(~yy) = [];
        
        % Scale the WRF surface pressure using the hypsometric equation and
        % the GLOBE terrain height.
        h_wrf = nanmean(tmp_elev);
        h_globe = globe_elevation(p);
        t_surf = nanmean(tmp_temp(1,:));
        R = 287; % J/kg/K, gas constant for dry air
        gamma = 6.5/1000; % K/m, lapse rate
        g = 9.8; % m/s^2, gravitational acceleration
        pSurf_WRF = nanmean(tmp_pres(1,:));
        pSurf = pSurf_WRF .* (t_surf ./ (t_surf + gamma .* (h_wrf - h_globe))) .^ (-g ./ (R .* gamma));
        
        
        % Interpolate all the NO2 and temperature profiles to the input
        % pressures, then average them. Extrapolate so that later we can be
        % sure to have one bin below the surface pressure for omiAmfAK2 and
        % integPr2. Interpolate NO2 in log-log space to account for the
        % exponential dependence of pressure on altitude and the often
        % exponential decrease of concentration with altitude. For
        % temperature, since a lapse rate implicitly assumes a linear
        % dependence of T on z, and since z is proportional to log(p):
        %
        %   z   = H * ln(p0/p)
        %       = H * ln(p0) - H * ln(p)
        %       = C - H * ln(p)
        %   z \propto ln(p)
        %
        % therefore, T should be proportional to ln(p)

        interp_no2 = nan(length(pressures), size(tmp_no2,2));
        interp_temp = nan(length(pressures), size(tmp_temp,2));
        
        if ~iscolumn(pressures); pressures = pressures'; end
        
        for a=1:size(tmp_no2,2)
            interp_no2(:,a) = interp1(log(tmp_pres(:,a)), log(tmp_no2(:,a)), log(pressures), 'linear', 'extrap');
            interp_temp(:,a) = interp1(log(tmp_pres(:,a)), tmp_temp(:,a), log(pressures), 'linear', 'extrap');
        end
        
        interp_no2 = exp(interp_no2);
        % do not need exp(interp_temp) since did not take the log of
        % tmp_temp
        
        pTropo = nanmean(tmp_pTropo);
        
        if clip_at_int_limits
            last_below_surf = find(pressures > pSurf,1,'last')-1;
            interp_no2(1:last_below_surf,:) = nan;
            interp_temp(1:last_below_surf,:) = nan;
            last_up_tropo = find(pressures < pTropo,1,'first')+1;
            interp_no2(last_up_tropo:end,:) = nan;
            interp_temp(last_up_tropo:end,:) = nan;
        end
        
        no2_vec = nanmean(interp_no2,2);
        temp_vec = nanmean(interp_temp,2);
        
    end

   function [wrf_no2, wrf_temp, wrf_pres, wrf_lon, wrf_lat, wrf_file, wrf_tropopres, wrf_elevation, pressure_mode, temperature_mode,tropopause_interp_indx] = load_wrf_vars()
       % Find the file for this day and the nearest hour May be "wrfout" or
        % "wrfout_subset"
        year_in = year(date_num_in);
        month_in = month(date_num_in);
        day_in = day(date_num_in);
        
        omi_utc_mean = omi_time_conv(nanmean(omi_time(:)));
        utc_hr = round(hour(omi_utc_mean) + minute(omi_utc_mean)/60);
        if strcmpi(profile_mode, 'daily')
            file_name = sprintf('wrfout_*_%04d-%02d-%02d_%02d-00-00', year_in, month_in, day_in, utc_hr);
            % Allow for the possibility that the filenames are "unsanitized" and have colons in them still
            file_name2 = sprintf('wrfout_*_%04d-%02d-%02d_%02d:00:00', year_in, month_in, day_in, utc_hr);
        elseif strcmpi(profile_mode, 'monthly')
            file_name = sprintf('WRF_BEHR_monthly_%02d.nc', month_in);
        end
        
        F = dir(fullfile(wrf_output_path,file_name));
        if numel(F) < 1 && strcmpi(profile_mode, 'daily')
            % For daily profiles, try the unsanitized name (with colons)
            % if we haven't found a file
            F = dir(fullfile(wrf_output_path, file_name2));
        end
        
        % Ensure we found exactly 1 file
        if numel(F) < 1
            E.filenotfound(file_name);
        elseif numel(F) > 1
            E.toomanyfiles(file_name);
        else
            wrf_info = ncinfo(fullfile(wrf_output_path,F(1).name));
        end
        
        % Load NO2 and check what unit it is - we'll use that to convert to
        % parts-per-part later. Make the location of units as general as possible
        try
            wrf_no2 = ncread(wrf_info.Filename, 'no2');
        catch err
            if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
                E.callCustomError('ncvar_not_found','no2',F(1).name);
            else
                rethrow(err);
            end
        end
        
        wrf_no2_units = ncreadatt_default(wrf_info.Filename, 'no2', 'units', 'ppm', 'fatal_if_missing', error_if_missing_attr);
        
        % Convert to be an unscaled mixing ratio (parts-per-part). Allow
        % the units to be different capitalization (i.e. since CMAQ seems
        % to output units of ppmV instead of ppm or ppmv).
        wrf_no2 = convert_units(wrf_no2, wrf_no2_units, 'ppp', 'case', false);
        
        wrf_vars = {wrf_info.Variables.Name};
        pres_precomputed = ismember('pres', wrf_vars);
        temp_precomputed = ismember('TT', wrf_vars);
        
        % Load the remaining variables
        try
            % TT is a custom variable that can be created using
            % calculated_quantities.nco or calculated_met_quantities.nco in
            % the WRF-nco-utils repo
            % (https://github.com/CohenBerkeleyLab/WRF-nco-tools). It
            % converts from WRF's perturbation to potential temperature to
            % an absolute temperature
            if temp_precomputed
                wrf_temp = ncread(wrf_info.Filename, 'TT');
                temp_units = ncreadatt_default(wrf_info.Filename, 'TT', 'units', 'K', 'fatal_if_missing', error_if_missing_attr);
                if ~strcmp(temp_units, 'K')
                    E.notimplemented('WRF temperature not in Kelvin');
                end
                temperature_mode = 'precomputed';
            else
                wrf_temp = convert_wrf_temperature(wrf_info.Filename, 'err_if_missing_units', error_if_missing_attr);
                temperature_mode = 'online';
            end
            
            if pres_precomputed % P and PB are already combined into the 'pres' variable in the monthly files
                varname = 'pres';
                p_tmp = ncread(wrf_info.Filename, varname);
                pb_tmp = 0; % Allows us to skip a second logical test later
                p_units = ncreadatt_default(wrf_info.Filename, 'pres', 'units','Pa', 'fatal_if_missing', error_if_missing_attr);
                pb_units = p_units;
                pressure_mode = 'precomputed';
            else
                varname = 'P';
                p_tmp = ncread(wrf_info.Filename, varname);
                p_units = ncreadatt_default(wrf_info.Filename, 'P', 'units', 'Pa', 'fatal_if_missing', error_if_missing_attr);
                varname = 'PB';
                pb_tmp = ncread(wrf_info.Filename, varname);
                pb_units = ncreadatt_default(wrf_info.Filename, 'PB', 'units', 'Pa', 'fatal_if_missing', error_if_missing_attr);
                pressure_mode = 'online';
            end
            
            wrf_elevation = read_wrf_preproc(wrf_info.Filename, 'elevation');
            wrf_elevation = wrf_elevation(:,:,1); % only need the surface elevation to adjust surface pressure
            
            varname = 'XLONG';
            wrf_lon = ncread(wrf_info.Filename, varname);
            varname = 'XLAT';
            wrf_lat = ncread(wrf_info.Filename, varname);
            [~,wrf_tropopres] = find_wrf_tropopause( wrf_info, 'error_if_missing_units', error_if_missing_attr );
        catch err
            if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
                E.callCustomError('ncvar_not_found',varname,F(1).name);
            else
                rethrow(err);
            end
        end

        % extrapolation wrf tropopause pressure when it's equal to 0
        tropopause_interp_indx = (wrf_tropopres == 0);

        if all(tropopause_interp_indx(:))
            E.notimplemented('No tropopause values found - no method to recover has been implemented')
        elseif any(tropopause_interp_indx(:)) 
            indx_nan = isnan(wrf_tropopres);
            wrf_tropopres(tropopause_interp_indx) = nan;
            
            %%% interplate the edge of wrf domain
            tp_nonbc = wrf_tropopres(2:end-1,2:end-1);
            wrf_tropopres(2:end-1,2:end-1) = nan;
            indx = ~isnan(wrf_tropopres);
            wrf_lon_good = wrf_lon(indx);
            wrf_lat_good = wrf_lat(indx);
            wrf_tropo_good = wrf_tropopres(indx);
            TropoScatObj = scatteredInterpolant(double(wrf_lon_good(:)),double(wrf_lat_good(:)),double(wrf_tropo_good(:)),'nearest');
            indx = isnan(wrf_tropopres);
            wrf_lon_bad = wrf_lon(indx);
            wrf_lat_bad = wrf_lat(indx);
            wrf_tropopres(indx) = TropoScatObj(double(wrf_lon_bad(:)),double(wrf_lat_bad(:)));
            wrf_tropopres(2:end-1,2:end-1) = tp_nonbc;
                
            %%% interpolate the rest of domain    

            indx = ~isnan(wrf_tropopres);
            wrf_lon_good = wrf_lon(indx);
            wrf_lat_good = wrf_lat(indx);
            wrf_tropo_good = wrf_tropopres(indx);
            TropoScatObj = scatteredInterpolant(double(wrf_lon_good(:)),double(wrf_lat_good(:)),double(wrf_tropo_good(:)));
            indx = isnan(wrf_tropopres);
            wrf_lon_bad = wrf_lon(indx);
            wrf_lat_bad = wrf_lat(indx);
            wrf_tropopres(indx) = TropoScatObj(double(wrf_lon_bad(:)),double(wrf_lat_bad(:)));
            wrf_tropopres(indx_nan) = nan;
        end

        if ~strcmp(pb_units, p_units)
            E.callError('unit_mismatch', 'Units for P and PB in %s do not match', wrf_info.Filename);
        end
        wrf_pres = convert_units(p_tmp + pb_tmp, p_units, 'hPa'); % convert from units in the wrfout files to hPa
        wrf_file = wrf_info.Filename;
    end
end


