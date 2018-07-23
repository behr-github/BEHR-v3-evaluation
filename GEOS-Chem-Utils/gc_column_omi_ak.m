function [ gc_no2 ] = gc_column_omi_ak( gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp, retrieval, no2_nox_ratio_bool )
%GC_COLUMN_OMI_AK Applies OMI averaging kernels to calculation of GEOS-Chem columns.
%   Averaging kernels are used to describe how a change in the "true" state
%   of a system is expressed as a change in the observations of a system.
%   In terms of satellite observations, a modeled NO2 column is usually
%   taken as "truth" (whether it is actually true is another matter
%   entirely) and the satellite column as observations. Thus, a vector of
%   averaging kernels is used to weight a modeled column for an
%   apples-to-apples comparison with a satellite observation.
%
%   For satellite retrievals, averaging kernels are defined as a vector of
%   box air mass factors or scattering weights scaled by the inverse of the
%   total AMF.  From Burrows, Platt, and Borrel (ed.), this is specifically
%   defined as:
%
%       AK_i = BAMF_i / AMF
%
%   and the box air mass factor (BAMF) is at one point defined as:
%
%       BAMF_i = SCD_i / VCD_i
%
%   therefore is the air mass factor for level i only. The observation VCD
%   can then be calculated from the true VCD as:
%
%       VCD_obs = sum_i(AK_i * VCD_i) = 1/AMF * sum_i(BAMF_i * VCD_i)
%
%   This is essentially calculating the observed slant column (the result
%   of BAMF_i * VCD_i summed over all levels) then dividing by the AMF to
%   simulate the retrieval process.
%
%   For OMI, the box air mass factors are given instead as scattering
%   weights, which reflect a similar quantity: basically the sensitivity of
%   the satellite to NO2 at that pressure level. From Eskes and Boersma
%   (2003), we know that the BAMF is independent of trace gas distribution
%   for optically thin absorbers, and NO2 is an optically thin absorber in
%   the atmosphere. Therefore, the scattering weights (which depend on
%   optical properties of the atmosphere & earth plus the geometric
%   arrangement of the sun and satellite) are sufficient for this purpose.
%
%   Reading:
%       Burrows, John P.; Platt, Ulrich; Borrell, Peter (eds.). The Remote
%       Sensing of Tropospheric Composition from Space. Springer, 2011. see
%       pp. 91-97.
%
%       Eskes, H.J. and Boersma, K.F. Averaging Kernels for DOAS
%       total-column satellite retrievals. Atmos. Chem. Phys. (2003), 3,
%       1285-1291.
%
%
%   There are three variables to be set, either as globals in a runscript or
%   manually in the code. omi_he5_dir must point to the base directory
%   for the satellite files.  This directory should contain the files
%   sorted into year and month subfolders (the months within each year folder).
%   run_mode determines if the function only tries to bin the AKs, apply them,
%   or both. It is one of the strings 'bin_aks','apply_aks','both'. Finally
%   ak_save_dir must point to the directory where the binned aks can be saved
%   or loaded. It should contain two subfolders, OMNO2 and DOMINO.
%
%   This function will apply OMNO2 averaging kernels to GEOS-Chem output.
%   It requires as input 5 output structures from GEOS-Chem created using
%   read_geos_output or by using read_gc_nd51.py and
%   convert_py_nd51_structure.m: the NO2 mixing ratios, the box heights, the matrix of
%   pressures, the number density of air, and the tropopause level. It will
%   assume that you want to apply the averaging kernels to the full time
%   extent of the NO2 structure.
%
%   Alternatively if running in "bin_aks" mode, the first two inputs should
%   be the start and end dates of the period to bin for.  This allows you to
%   break the year up into smaller chunks if necessary.
%
%   For each day, this function will load all 14 OMI orbits then average
%   the vectors of averaging kernels to the GEOS-Chem grid boxes.  The
%   averaging kernels 
%
%   Dependencies:
%       Classes/JLLErrors.m

E=JLLErrors;

% Get run directories
global onCluster
if isempty(onCluster)
    onCluster = false;
end

shareroot = getenv('SYNOMNT');

fprintf('onCluster = %d\n',onCluster);

global omi_he5_dir
if onCluster && isempty(omi_he5_dir)
    E.runscript_error('omi_he5_dir')
elseif ~onCluster && ~isempty(shareroot)
    omi_he5_dir = fullfile(shareroot,'share-sat','SAT','OMI','OMNO2');
elseif ~onCluster
    omi_he5_dir = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0';
end

global run_mode;
if onCluster && isempty(run_mode)
    E.runscript_error('run_mode');
elseif ~onCluster;
    run_mode = 'bin_aks';
    fprintf('Setting run_mode to %s\n',run_mode);
end

global ak_save_dir;
if onCluster && isempty(ak_save_dir)
    E.runscript_error('ak_save_dir');
elseif ~onCluster && ~isempty(shareroot)
    ak_save_dir = fullfile(shareroot,'share2','USERS','LaughnerJ','MPN_Project','AKs');
    fprintf('Setting ak_save_dir to %s\n',ak_save_dir);
elseif ~onCluster
    ak_save_dir = '/Volumes/share2/USERS/LaughnerJ';
    fprintf('Setting ak_save_dir to %s\n',ak_save_dir);
end

global overwrite
if onCluster && isempty(overwrite)
    E.runscript_error('overwrite')
elseif ~onCluster
    overwrite = true;
end

allowed_modes = {'bin_aks','apply_aks','both'};
if ~ismember(run_mode, allowed_modes)
    E.badinput('run_mode must be one of %s', strjoin(allowed_modes, ', '));
end

fprintf('OMI dir = %s\n',omi_he5_dir);
fprintf('Run mode = %s\n',run_mode);
fprintf('AK save dir = %s\n',ak_save_dir);

[gc_loncorn, gc_latcorn] = geos_chem_corners;

if isstruct(gc_no2)
    tVec = gc_no2.tVec;
    gc_no2_data = gc_no2.dataBlock;
    gc_bxhght_data = gc_bxhght.dataBlock;
    gc_pressure_data = gc_pressure.dataBlock;
    gc_ndens_air_data = gc_ndens_air.dataBlock;
    gc_tp_data = gc_tp.dataBlock;
    columns = nan(size(gc_tp_data));
    total_weights = nan(size(gc_tp_data));
elseif strcmpi(run_mode,'bin_aks')
    tVec = datenum(gc_no2):datenum(gc_bxhght);
    % parfor needs these to be slicable in the proper dimensions.
    gc_no2_data = nan(1,1,1,numel(tVec));
    gc_bxhght_data = nan(1,1,1,numel(tVec));
    gc_pressure_data = nan(1,1,1,numel(tVec));
    gc_ndens_air_data = nan(1,1,1,numel(tVec));
    gc_tp_data = nan(1,1,numel(tVec));
    columns = nan(144, 91, numel(tVec));
    total_weights = nan(1, 1, numel(tVec));
else
    E.badinput('If not running in bin only mode, then the GEOS-Chem output structures MUST be passed as arguments');
end

if ~exist('no2_nox_ratio_bool', 'var') && strcmpi(run_mode, 'apply_aks')
    E.badinput('In apply mode, NO2_NOX_RATIO_BOOL must be given');
elseif ~islogical(no2_nox_ratio_bool) || ~isscalar(no2_nox_ratio_bool)
    E.badinput('SCALE_NO2_NOX_RATIO must be a scalar logical');
elseif ~no2_nox_ratio_bool
    no2nox.dc3_ratio = [1 1];
    no2nox.gc_ratio = [1 1];
    no2nox.dc3_pres = [1013, 200];
    no2nox.gc_pres = [1013, 200];
else
    mydir = fileparts(mfilename('fullpath'));
    ratio_file = fullfile(mydir, 'Data', 'no2nox_ratios.txt');
    [ no2nox.dc3_pres, no2nox.dc3_ratio, no2nox.gc_pres, no2nox.gc_ratio ] = read_no2_nox_ratio_file( ratio_file );
end
    

sz_wt = size(total_weights);

parfor d=1:numel(tVec)
%for d=1:numel(tVec)
    fprintf('Loading OMI files for %s\n',datestr(tVec(d)));
    earth_ellip = referenceEllipsoid('wgs84','kilometer');
    if strcmpi(retrieval,'omno2')
        save_name = sprintf('OMNO2_AK_%04d%02d%02d.mat',year(tVec(d)),month(tVec(d)),day(tVec(d)));
        if strcmpi(run_mode,'bin_aks') && exist(fullfile(ak_save_dir,'OMNO2',save_name),'file') && ~overwrite
            fprintf('File %s exists, skipping\n',save_name);
            continue
        end
        if ismember(run_mode,{'bin_aks','both'})
            [omi_aks, omi_lon, omi_lat, omi_loncorn, omi_latcorn] = load_omi_files(year(tVec(d)), month(tVec(d)), day(tVec(d)), omi_he5_dir);
            omi_pixweight = calc_pix_areaweight(omi_loncorn, omi_latcorn);
            fprintf('Binnind OMI AKs for %s\n',datestr(tVec(d)));
            [binned_aks, binned_weights] = bin_omi_aks(gc_loncorn, gc_latcorn, omi_aks, omi_lon, omi_lat, omi_loncorn, omi_latcorn, omi_pixweight, earth_ellip);
            
            % binned_aks is output in ak_vec x gc_nlat x gc_nlon, rearrange to
            % gc_nlon x gc_nlat x ak_vec x time.
            binned_aks = permute(binned_aks, [3 2 1]);
            binned_weights = binned_weights';

            saveOmno2(fullfile(ak_save_dir,'OMNO2',save_name),binned_aks,binned_weights);
        elseif ismember(run_mode,{'apply_aks','both'})
            if strcmpi(run_mode,'apply_aks')
                AK = load(fullfile(ak_save_dir,'OMNO2',save_name));
                binned_aks = AK.binned_aks;
                binned_weights = AK.binned_weights;
            end
        
            columns(:,:,d) = integrate_omi_profile(gc_no2_data(:,:,:,d), gc_bxhght_data(:,:,:,d), gc_pressure_data(:,:,:,d), gc_ndens_air_data(:,:,:,d), gc_tp_data(:,:,d), binned_aks, no2nox);
            tmp_wts = nan(sz_wt(1:2));
            for i=1:sz_wt(1)
                for j=1:sz_wt(2)
                    if ~isempty(binned_weights{i,j})
                        tmp_wts(i,j) = nansum2(binned_weights{i,j});
                    end
                end
            end
            total_weights(:,:,d) = tmp_wts;
        end
    elseif strcmpi(retrieval,'domino')
        save_name = sprintf('DOMINO_AK_%04d%02d%02d.mat',year(tVec(d)),month(tVec(d)),day(tVec(d)));
        if strcmpi(run_mode,'bin_aks') && exist(fullfile(ak_save_dir,'DOMINO',save_name),'file') && ~overwrite
            fprintf('File %s exists, skipping\n',save_name);
            continue
        end
        if ismember(run_mode, {'bin_aks','both'})
            [omi_aks, omi_pres, omi_pres_edge, omi_lon, omi_lat, omi_loncorn, omi_latcorn] = load_domino_files(year(tVec(d)), month(tVec(d)), day(tVec(d)), omi_he5_dir);
            omi_pixweight = calc_pix_areaweight(omi_loncorn, omi_latcorn);
            fprintf('Binning OMI AKs for %s\n',datestr(tVec(d)));
            [binned_aks, binned_weights, binned_pres, binned_pres_edge, DB] = bin_omi_aks(gc_loncorn, gc_latcorn, omi_aks, omi_lon, omi_lat, omi_loncorn, omi_latcorn, omi_pixweight, earth_ellip, omi_pres, omi_pres_edge);
            
            % binned_aks and binned_pres are output in gc_nlat x gc_nlon with
            % the ak vectors in each cell, rearrange so that lon is the first
            % dimension
            binned_aks = binned_aks';
            binned_pres = binned_pres';
            binned_pres_edge = binned_pres_edge';
            binned_weights = binned_weights';
            
            saveDomino(fullfile(ak_save_dir,'DOMINO', save_name), binned_aks, binned_pres, binned_pres_edge, binned_weights);
        
        elseif ismember(run_mode, {'apply_aks','both'})
            if strcmpi(run_mode, 'apply_aks')
                AK = load(fullfile(ak_save_dir,'DOMINO',save_name));
                binned_aks = AK.binned_aks;
                binned_pres = AK.binned_pres;
                binned_pres_edge = AK.binned_pres_edge;
                binned_weights = AK.binned_weights;
            end
            columns(:,:,d) = integrate_domino_profile(gc_no2_data(:,:,:,d), gc_bxhght_data(:,:,:,d), gc_pressure_data(:,:,:,d), gc_ndens_air_data(:,:,:,d), gc_tp_data(:,:,d), binned_aks, no2nox, binned_pres, binned_pres_edge, binned_weights);
            tmp_wts = nan(sz_wt(1:2));
            for i=1:sz_wt(1)
                for j=1:sz_wt(2)
                    if ~isempty(binned_weights{i,j})
                        tmp_wts(i,j) = nansum2(binned_weights{i,j});
                    end
                end
            end
            total_weights(:,:,d) = tmp_wts;
        end
    else
        E.notimplemented(retrieval)
    end
    
    
    %tmp = integrate_geoschem_profile(cat(1,gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp), 1e-9, [], 0, binned_aks);
    %gc_no2.Columns = tmp(1).Columns;
    
    
end

if ~strcmpi(run_mode,'bin_aks')
    gc_no2.Columns = columns;
    gc_no2.TotalWeights = total_weights;
end

end

function saveOmno2(save_path, binned_aks, binned_weights)
        save(save_path,'binned_aks','binned_weights');
end

function saveDomino(save_path, binned_aks, binned_pres, binned_pres_edge, binned_weights)
    save(save_path, 'binned_aks', 'binned_pres', 'binned_pres_edge', 'binned_weights');
end

function weight = calc_pix_areaweight(loncorn, latcorn)
Lon1 = squeeze(loncorn(1,:,:));
Lon2 = squeeze(loncorn(2,:,:));
Lon3 = squeeze(loncorn(3,:,:));
Lon4 = squeeze(loncorn(4,:,:));

Lat1 = squeeze(latcorn(1,:,:));
Lat2 = squeeze(latcorn(2,:,:));
Lat3 = squeeze(latcorn(3,:,:));
Lat4 = squeeze(latcorn(4,:,:));

weight = nan(size(Lon1));

Amin = 312; % 13 km x 24 km
Amax = 7500; % approximate max seen from BEHR
for x=1:size(Lon1,1)
    for y=1:size(Lon1,2)
        pixelarea = (m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)]))*(m_lldist([Lon1(x,y)-180, Lon4(x,y)-180],[Lat1(x,y), Lat4(x,y)]));
        weight(x,y) = 1 - (pixelarea - Amin)/Amax;
    end 
end
end

function [omi_aks, omi_lon, omi_lat, omi_loncorn, omi_latcorn] = load_omi_files(yr, mn, dy, omi_he5_dir)
E = JLLErrors;

full_path = fullfile(omi_he5_dir,sprintf('%04d',yr),sprintf('%02d',mn));
file_pattern = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5',yr,mn,dy);
F = dir(fullfile(full_path, file_pattern));

fprintf('%s\n',full_path);
fprintf('%s\n',file_pattern);

if isempty(F)
    E.filenotfound(sprintf('OMI files for %04d-%02d-%02d',yr,mn,dy));
end


% We'll concatenate along the along-track dimension since that varys with
% each swath and the other two are constant.
omi_aks = [];
omi_lon = [];
omi_lat = [];
omi_loncorn = [];
omi_latcorn = [];

for a=1:numel(F)
    %fprintf('Loading file %d of %d\n',a,numel(F));
    hi = h5info(fullfile(full_path, F(a).name));
    omi.sw = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'ScatteringWeight')));
    omi.sw(omi.sw < 1e-29) = nan;
    omi.amf = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AmfTrop')));
    omi.amf(omi.amf < 1e-29) = nan;
    
    % These are fields needed to reject pixels
    omi.CloudFraction = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'CloudFraction')))*1e-3;
    omi.CloudFraction(omi.CloudFraction<0) = nan;
    omi.ColumnAmountNO2Trop = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'ColumnAmountNO2Trop')));
    omi.ColumnAmountNO2Trop(omi.ColumnAmountNO2Trop<-1e29) = nan;
    omi.vcdQualityFlags = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1, 'VcdQualityFlags')));
    omi.vcdQualityFlags(omi.vcdQualityFlags==65535)=nan;
    omi.XTrackQualityFlags = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'XTrackQualityFlags')));
    omi.XTrackQualityFlags(omi.XTrackQualityFlags==255) = nan;
    omi.TerrainReflectivity = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TerrainReflectivity')));
    omi.TerrainReflectivity(omi.TerrainReflectivity==-32767) = nan;
    omi.TerrainReflectivity = omi.TerrainReflectivity * 1e-3;
    omi.Areaweight = ones(size(omi.amf));
    
    omi = omi_sp_pixel_reject(omi,'geo',0.3,'XTrackFlags');
    
    rejects = omi.Areaweight == 0 | isnan(omi.CloudFraction) | isnan(omi.ColumnAmountNO2Trop)...
        | isnan(omi.TerrainReflectivity) | omi.TerrainReflectivity > 0.3;
    omi.sw(:,rejects) = nan;
    omi.amf(rejects) = nan;
    
    omi_sw = omi.sw;
    omi_amf = omi.amf;
    
    % Remove the large 'omi' structure since we have gotten what we need
    % out of it.
    clear('omi');
    
    this_aks = nan(size(omi_sw));
    for x=1:size(omi_amf,1)
        for y=1:size(omi_amf,2)
            this_aks(:,x,y) = omi_sw(:,x,y) ./ omi_amf(x,y);
        end
    end
    
    omi_aks = cat(3, omi_aks, this_aks);
    
    this_lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
    this_lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
    this_loncorn = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'FoV75CornerLongitude'));
    this_latcorn = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'FoV75CornerLatitude'));
    omi_lon = cat(2, omi_lon, this_lon);
    omi_lat = cat(2, omi_lat, this_lat);
    omi_loncorn = cat(3, omi_loncorn, this_loncorn);
    omi_latcorn = cat(3, omi_latcorn, this_latcorn);
    for p=1:numel(omi_lat)
        [omi_loncorn(:,p), omi_latcorn(:,p)] = uncross_corners(omi_loncorn(:,p), omi_latcorn(:,p));
    end
end

end

function [omi_aks, omi_pres, omi_pres_edge, omi_lon, omi_lat, omi_loncorn, omi_latcorn] = load_domino_files(yr, mn, dy, omi_he5_dir)
E = JLLErrors;

full_path = fullfile(omi_he5_dir,sprintf('%04d',yr),sprintf('%02d',mn));
file_pattern = sprintf('OMI-Aura_L2-OMDOMINO_%04dm%02d%02d*.he5',yr,mn,dy);
F = dir(fullfile(full_path, file_pattern));

fprintf('%s\n',full_path);
fprintf('%s\n',file_pattern);

if isempty(F)
    E.filenotfound(sprintf('OMI files for %04d-%02d-%02d',yr,mn,dy));
end


% We'll concatenate along the along-track dimension since that varys with
% each swath and the other two are constant.
omi_aks = [];
omi_pres = [];
omi_pres_edge = [];
omi_lon = [];
omi_lat = [];
omi_loncorn = [];
omi_latcorn = [];

for a=1:numel(F)
    % Load in data and remove fill values
    hi = h5info(fullfile(full_path, F(a).name));
    omi.aks = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AveragingKernel')))*0.001; 
    omi.aks(omi.aks<-30) = nan;
    omi.TM4PresA = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4PressurelevelA')))*0.01; % in Pa, want hPa.
    omi.TM4PresA(omi.TM4PresA<-1e30) = nan;
    omi.TM4PresB = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4PressurelevelB')));
    omi.TM4PresB(omi.TM4PresB<-1e30) = nan;
    omi.TM4SurfP = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'TM4SurfacePressure')));
    omi.TM4SurfP(omi.TM4SurfP<-1e30) = nan;
    % Rearrange so that the vertical or corner coordinate is first
    omi.aks = permute(omi.aks, [3 1 2]);

    % These are the aks vector using the total AMF. We need to convert to the tropospheric AMF
    amf = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AirMassFactor')));
    amf(amf < -1e30) = nan;
    amftrop = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,1,'AirMassFactorTropospheric')));
    amftrop(amftrop < -1e30) = nan;

    for x=1:size(amf,1)
        for y=1:size(amf,2)
            omi.aks(:,x,y) = omi.aks(:,x,y) * amf(x,y) / amftrop(x,y);
        end
    end
    
    % Calculate the pressure levels for each pixel. We'll also need the
    % bottom pressure edges to compare with box height in the integration
    % function. We'll assume that the edges are half way in between the
    % levels, except for the first one, which will be at the surface.
    omi.PresLevs = nan(size(omi.aks));
    omi.PresEdges = nan(size(omi.aks) + [1 0 0]); % need extra vertical coordinate for edge
    for x=1:size(omi.TM4SurfP,1)
        for y=1:size(omi.TM4SurfP,2)
            omi.PresLevs(:,x,y) = omi.TM4PresA + omi.TM4SurfP(x,y) .* omi.TM4PresB;
            omi.PresEdges(1,x,y) = omi.TM4SurfP(x,y);
            omi.PresEdges(2:end-1,x,y) = (omi.PresLevs(1:end-1,x,y) + omi.PresLevs(2:end,x,y)) / 2;
            omi.PresEdges(end,x,y) = 0.5 * omi.PresLevs(end,x,y); % set the final edge to 1/2 the box center - it's stratosphere, it won't matter anyway.
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
    
    bad_vals = omi.CloudFraction > 0.3 | omi.TropColumnFlag < 0 | omi.Albedo > 0.3 | isnan(omi.ColumnAmountNO2Trop) | isnan(omi.Albedo) | isnan(omi.CloudFraction) ;
    omi.aks(:,bad_vals) = nan;
    omi.PresLevs(:,bad_vals) = nan;
    
    omi_aks = cat(3,omi_aks,omi.aks);
    omi_pres = cat(3,omi_pres,omi.PresLevs);
    omi_pres_edge = cat(3, omi_pres_edge, omi.PresEdges);
    
    this_lon = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Longitude'));
    this_lat = h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'Latitude'));
    this_loncorn = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'LongitudeCornerpoints')));
    this_latcorn = double(h5read(hi.Filename, h5dsetname(hi,1,2,1,2,'LatitudeCornerpoints')));
    this_loncorn = permute(this_loncorn, [3 1 2]); 
    this_latcorn = permute(this_latcorn, [3 1 2]);
    omi_lon = cat(2, omi_lon, this_lon);
    omi_lat = cat(2, omi_lat, this_lat);
    omi_loncorn = cat(3, omi_loncorn, this_loncorn);
    omi_latcorn = cat(3, omi_latcorn, this_latcorn);
    for p=1:numel(omi_lat)
        [omi_loncorn(:,p), omi_latcorn(:,p)] = uncross_corners(omi_loncorn(:,p), omi_latcorn(:,p));
    end
end

end

function varargout = bin_omi_aks(gc_loncorn, gc_latcorn, omi_aks, omi_lon, omi_lat, omi_loncorn, omi_latcorn, omi_aw, earth_ellip, omi_pres, omi_pres_edge)
% This subfunction will handle the binning of OMNO2 averaging kernels to
% the GEOS-Chem grid. Inputs: matrices of GEOS-Chem corner points and
% matrices of OMI pixel AKs (with the different AK levels along the first
% dimension) and pixel lon/lats.  gc_loncorn and gc_latcorn are expected to
% be matrices that are (m+1)-by-(n+1) if there are m-by-n GC grid cells.

E = JLLErrors;
t = getCurrentTask();
if isempty(t)
    t.ID = 0;
end

fprintf('Checking input to bin_omi_aks\n');

if exist('omi_pres','var')
    retrieval = 'domino';
else
    retrieval = 'omno2';
end
    
if any(size(gc_loncorn) ~= size(gc_latcorn))
    E.sizeMismatch('gc_loncorn','gc_latcorn');
end

sz_aks = size(omi_aks);
sz_omilon = size(omi_lon);
sz_omilat = size(omi_lat);

if ndims(omi_aks) ~= 3
    E.badinput('omi_aks expected to be 3 dimensional')
elseif ~ismatrix(omi_lon)
    E.badinput('omi_lon expected to be 2 dimensional')
elseif ~ismatrix(omi_lat)
    E.badinput('omi_lat expected to be 2 dimensional')
end

% There should usually be fewer points in the latitude direction than the
% longitude direction of GEOS-Chem, use this to make sure the matrices have
% the right orientation.
if size(gc_loncorn,1) > size(gc_loncorn,2)
    E.badinput('The GC corner matrices are expected to be lat x lon')
end
if any(diff(gc_loncorn(1,:)) < 0) || any(diff(gc_latcorn(:,1)) < 0)
    E.badinput('The GC corner matrices are expected to have the SW corner at (1,1) and NE corner at (end,end)')
end


if any(sz_aks(2:3)~=sz_omilon) || any(sz_aks(2:3)~=sz_omilat)
    E.badinput('The second and third dimensions of omi_aks should be the same as the dimensions of omi_lon and omi_lat')
end

% Loop over all GEOS-Chem cells and bin the averaging kernels to them.
fprintf('Looping over all GEOS-Chem cells\n');
gc_nlat = size(gc_loncorn,1)-1;
gc_nlon = size(gc_loncorn,2)-1;
if strcmpi(retrieval,'omno2')
    aks = nan(size(omi_aks,1), gc_nlat, gc_nlon);
    weights = cell(gc_nlat, gc_nlon);
elseif strcmpi(retrieval,'domino')
    aks = cell(gc_nlat, gc_nlon);
    pres = cell(gc_nlat, gc_nlon);
    pres_edge = cell(gc_nlat, gc_nlon);
    weights = cell(gc_nlat, gc_nlon);
else
    E.notimplemented(retrieval)
end
DB=struct;
for a=1:gc_nlat
    fprintf('Binning %.1f%% complete\n',a/gc_nlat*100);
    for b=1:gc_nlon
        x1 = gc_loncorn(a,b);
        x2 = gc_loncorn(a,b+1);
        y1 = gc_latcorn(a,b);
        y2 = gc_latcorn(a+1,b);
        
        xx = omi_lon >= x1 & omi_lon < x2 & omi_lat >= y1 & omi_lat < y2;
        
        if strcmpi(retrieval,'omno2')
            % All OMI aks are given at the same pressures, so we can just
            % average across all kernels that belong to this grid cell.
            if sum(xx(:)) > 0
                this_ak = omi_aks(:,xx);
                this_aw = omi_aw(xx);
                this_loncorn = omi_loncorn(:,xx);
                this_latcorn = omi_latcorn(:,xx);
                Q = calc_overlap_weight([x1 x2], [y1 y2], this_loncorn, this_latcorn, earth_ellip);
                W = repmat(Q .* this_aw', size(this_ak,1), 1);
                W(:,all(isnan(this_ak),1)) = nan;
                this_ak_mean = nansum2(this_ak .* W, 2) ./ nansum2(W(1,:),2);
                aks(:,a,b) = this_ak_mean;
                weights{a,b} = W(1,:);
            else
                aks(:,a,b) = nan(size(omi_aks,1),1);
                weights{a,b} = nan;
            end
        elseif strcmpi(retrieval,'domino')
            % DOMINO AKs are NOT at the same pressure level, therefore when
            % we integrate, what we're going to have to do is apply each
            % averaging kernel in turn (interpolating the GC profile to
            % each set of pressure levels) then average the column at the
            % end.
            if sum(xx(:)) < 1; continue; end
            aks{a,b} = omi_aks(:,xx);
            pres{a,b} = omi_pres(:,xx);
            pres_edge{a,b} = omi_pres_edge(:,xx);
            
            this_loncorn = omi_loncorn(:,xx);
            this_latcorn = omi_latcorn(:,xx);
            Q = calc_overlap_weight([x1 x2], [y1 y2], this_loncorn, this_latcorn, earth_ellip);
            W = Q .* omi_aw(xx)';
            W(all(isnan(aks{a,b}),1)) = nan;
            %fprintf('W%d: (%d, %d): Size(Q) = %s, size(W) = %s, size(omi_aw(xx))'' = %s, size(aks{a,b}) = %s\n', t.ID, a, b, mat2str(size(Q)), mat2str(size(W)), mat2str(size(omi_aw(xx)')), mat2str(size(aks{a,b})));
            weights{a,b} = W;
            if b == 102 && a == 41
                DB.Q = Q;
                DB.AW = omi_aw(xx)';
                DB.W = W;
                DB.omilon = omi_lon(xx)';
                DB.omilat = omi_lat(xx)';
                DB.gclon = [x1, x2];
                DB.gclat = [y1, y2];
            end
        end
    end
end

varargout{1} = aks;
varargout{2} = weights;
if strcmpi(retrieval,'domino')
    varargout{3} = pres;
    varargout{4} = pres_edge;
    varargout{5} = DB;
end

end

function Q = calc_overlap_weight(gc_lon, gc_lat, pixloncorn, pixlatcorn, earth_ellip)
% get the grid cell corners
t = getCurrentTask;
%t.ID = 0;
Q = zeros(1,size(pixloncorn,2));
gc_xall = [gc_lon(1), gc_lon(1), gc_lon(2), gc_lon(2), gc_lon(1)];
gc_yall = [gc_lat(1), gc_lat(2), gc_lat(2), gc_lat(1), gc_lat(1)];
% ensure both are clockwise
[gc_xall, gc_yall] = poly2cw(gc_xall, gc_yall);
gc_area = areaint(gc_yall, gc_xall, earth_ellip);
for a=1:size(pixloncorn,2)
    if any(isnan(pixloncorn(:,a))) || any(isnan(pixlatcorn(:,a))) || any(pixloncorn(:,a) < -180) || any(pixloncorn(:,a) > 180) || any(pixlatcorn(:,a) < -90) || any(pixlatcorn(:,a) > 90) || (any(sign(pixloncorn(:,a))~=sign(pixloncorn(1,a))) && (any(abs(pixloncorn(:,a))>90) || any(abs(pixlatcorn(:,a))>80)))
        % the last test handles issues where a pixel straddles the international
        % date line. there's better ways to handle it (wrap the pixel corner around
        % to be the same sign as the GC corners) but I don't feel like doing that atm.
        continue;
    end
    pixloncorn_a = pixloncorn(:,a);
    pixlatcorn_a = pixlatcorn(:,a);
    %[pixloncorn_a, pixlatcorn_a] = uncross_corners(pixloncorn(:,a), pixlatcorn(:,a));
    [pixloncorn_a, pixlatcorn_a] = poly2cw(pixloncorn_a, pixlatcorn_a);
    % create a polygon that represents the area of overlap and calculate its area
    % in km.
    [xt, yt] = polybool('intersection',gc_xall,gc_yall,pixloncorn_a,pixlatcorn_a);
    if isempty(xt) || isempty(yt)
        continue
    elseif any(isnan(xt))
        error('load_and_grid_domino:calc_pix_grid_overlap','W%d: The pixel corners are wrong (%s, %s) at %d',t.ID,mat2str(pixloncorn_a),mat2str(pixlatcorn_a), a);
    end
    overlap_area = areaint(yt,xt,earth_ellip);
    Q(a) = overlap_area/gc_area;
end

end

function [x,y]= uncross_corners(x,y)
    m1 = (y(3) - y(2))/(x(3) - x(2));
    b1 = y(2) - m1*x(2);
    m2 = (y(4) - y(1))/(x(4) - x(1));
    b2 = y(1) - m2*x(1);
    flip_bool = false;
    if ~isinf(m1) && ~isinf(m2)
        % As long as neither slope is infinite, solve for the x-coordinate
        % of the intercept and see if it falls inside the polygon - if so,
        % the corners need flipped.
        inpt = (b2-b1)/(m1-m2);
        if inpt > min(x(2:3)) && inpt < max(x(2:3))
            flip_bool = true;
        end 
    elseif isinf(m1) && ~isinf(m2)
        % If one is infinite, fine the y-coord where the other one is at
        % it's x-coordinate and do the same test
        inpt = m2*x(2)+b2;
        if inpt > min(y(2:3)) && inpt < max(y(2:3))
            flip_bool = true;
        end 
    elseif isinf(m2) && ~isinf(m1)
        inpt = m1*x(1) + b1; 
        if inpt > min(y([1,4])) && inpt < max(y([1,4]))
            flip_bool = true;
        end 
        % If both are infinite, they are parallel and the corners do not
        % need flipped.
    end 
    if flip_bool
        tmp = x(4);
        x(4) = x(3);
        x(3) = tmp;
        tmp = y(4);
        y(4) = y(3);
        y(3) = tmp;
    end 
end

function no2_columns = integrate_omi_profile(gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp, binned_aks, no2nox_ratios)
% This will integrate the geos-chem columns weighted by the OMI AKs after
% interpolating to the OMI pressures.  The integration will be carried out
% after that in the same manner as integrate_geoschem_profile, that is,
% the integral will be calculated using the rectangular rule.
%
% In order to find the box height in meters, the height of each GC pressure
% level will be calculated as the cumulative sum of the box height values,
% which can then be intepolated to the OMI pressure levels, and the box
% height for level i taken as the difference between level i and i+1
% altitude in meters. Number density will also be interpolated to the new
% pressure levels, and the GC tropopause will be used to restrict the
% column.

% First, handle all the necessary intepolation to regrid the GC variables
% to the OMI pressure levels.
omi_pres = [1020 1010 1000 990 975 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200 120 60 35 20 12 8 5 3 1.5 0.8];
omi_pres_edge = [1025 1015 1005 995 982.5 967.5 952.5 935 912.5 887.5 862.5 837.5 812.5 785 755 720 680 635 585 525 475 425 375 315 240 150 72.5 42.5 24 14 9.5 6 3.75 1.85 1.15 0.45];

sz = size(gc_no2);
if numel(sz) < 4;
    sz = cat(2, sz, ones(1,4-numel(sz)));
end
gc_z = cumsum(gc_bxhght,3);
omi_bxhght = nan(sz(1), sz(2), numel(omi_pres), sz(4));
omi_ndens_air = nan(sz(1), sz(2), numel(omi_pres), sz(4));
omi_no2 = nan(sz(1), sz(2), numel(omi_pres), sz(4));
for a=1:sz(1)
    for b=1:sz(2)
        for t=1:sz(4)
            % For each column, interpolate box altitude to OMI pressures,
            % height is proportional to ln(p), so use linear interpolation
            % w.r.t. ln(p):
            %   z = -ln(p/p0)*H = H*log(p0) - H*log(p)
            % When GC outputs actual pressure edges the last edge is 0.01
            % which is not included in the satellite pressure values, we'll
            % include it to have edges to work with, then also add 0 to the
            % boxheight values so that we can take the difference.
            gc_pvec = cat(1,squeeze(gc_pressure(a,b,:,t)),0.01);
            gc_zvec = cat(1, 0, squeeze(gc_z(a,b,:,t)));
            omi_zvec = interp1(log(gc_pvec), gc_zvec, log(omi_pres_edge),'linear','extrap');
            omi_bxhght(a,b,:,t) = diff(omi_zvec);
            
            % Also interpolate number density to OMI pressure levels,
            % number density is directly proportional to p:
            %   pV = nRT => n/V = p/RT
            % so we'll directly interpolate.
            gc_nair_vec = squeeze(gc_ndens_air(a,b,:,t));
            omi_ndens_air(a,b,:,t) = interp1(gc_pvec(1:end-1), gc_nair_vec, omi_pres, 'linear','extrap');
            
            % Finally interpolate the GC NO2 mixing ratio values (scaling
            % to straight mixing ratios, i.e. parts-per-part). I'll do this
            % in log-log space, which I've seen before. At the same time,
            % scale the GC profiles to account for the discrepancy in
            % NO2:NOx ratio between GEOS-Chem and DC3. We will NOT
            % extrapolate this time, because I do not want there to be
            % values below the surface (so outside of the GC pressure
            % values, it will have the value of NaN)
            gc_no2vec = squeeze(gc_no2(a,b,:,t));
            gc_no2vec = scale_gc_prof(gc_no2vec, gc_pvec(1:end-1), no2nox_ratios);
            omi_no2(a,b,:,t) = exp(interp1(log(gc_pvec(1:end-1)), log(gc_no2vec), log(omi_pres)));
            strat = omi_pres < gc_pressure(a,b,floor(gc_tp(a,b,t)),t);
            omi_no2(a,b,strat,t) = nan;
        end
    end
end

% Finally the actual calculation of columns is fairly straightforward: sum
% over each box's partial column weighted by AK (after appropriate unit
% conversions)

omi_no2_ndens = (omi_no2 * 1e-9) .* (omi_ndens_air * 1e-6) .* (omi_bxhght * 100) .* binned_aks;
no2_columns = squeeze(nansum2(omi_no2_ndens,3));

end

function no2_columns = integrate_domino_profile(gc_no2, gc_bxhght, gc_pressure, gc_ndens_air, gc_tp, binned_aks, no2nox_ratios, binned_pres, binned_pres_edge, binned_weights)
% This will integrate the geos-chem columns weighted by the OMI AKs after
% interpolating to the OMI pressures.  The integration will be carried out
% after that in the same manner as integrate_geoschem_profile, that is,
% the integral will be calculated using the rectangular rule.
%
% In order to find the box height in meters, the height of each GC pressure
% level will be calculated as the cumulative sum of the box height values,
% which can then be intepolated to the OMI pressure levels, and the box
% height for level i taken as the difference between level i and i+1
% altitude in meters. Number density will also be interpolated to the new
% pressure levels, and the GC tropopause will be used to restrict the
% column.

tstr = getCurrentTask();
if isempty(tstr)
    tstr.ID = 0;
end

sz = size(gc_no2);
if numel(sz) < 4;
    sz = cat(2, sz, ones(1,4-numel(sz)));
end
gc_z = cumsum(gc_bxhght,3);
no2_columns = nan(sz(1), sz(2), sz(4));

for a=1:sz(1)
    for b=1:sz(2)
        for t=1:sz(4) % I think t is extraneous, but I'll leave it here just in case
            % For DOMINO, since each pixel has a different set of pressure
            % levels, we will interpolate GC output to the pressure levels
            % of each pixel, then calculate the VCD for that interpolation,
            % then average all the VCDs resulting from AKs in the grid cell
            % together to get the final, AK'ed VCD. This should be
            % mathematically equivalent (but likely less efficient
            % computationally) as the method I used for OMNO2 AKs.
            
            omi_aks = binned_aks{a,b};
            if size(omi_aks,2) == 0 || all(isnan(omi_aks(:)))
                continue
            end
            omi_pres = binned_pres{a,b};
            omi_pres_edge = binned_pres_edge{a,b};
            omi_weights = binned_weights{a,b};
            
            nans = all(isnan(omi_aks),1);
            omi_aks(:,nans) = [];
            omi_pres(:,nans) = [];
            omi_pres_edge(:,nans) = [];
            omi_weights(nans) = [];

            gc_pvec = cat(1,squeeze(gc_pressure(a,b,:,t)),0.01);
            gc_zvec = cat(1, 0, squeeze(gc_z(a,b,:,t)));
            gc_nair_vec = squeeze(gc_ndens_air(a,b,:,t));
            gc_no2vec = squeeze(gc_no2(a,b,:,t));
            gc_no2vec = scale_gc_prof(gc_no2vec, gc_pvec(1:end-1), no2nox_ratios);
            
            omi_no2_columns = nan(1,size(omi_aks,2));
            for p=1:size(omi_aks,2)
            
                % For each column, interpolate box altitude to OMI pressures,
                % height is proportional to ln(p), so use linear interpolation
                % w.r.t. ln(p):
                %   z = -ln(p/p0)*H = H*log(p0) - H*log(p)
                % When GC outputs actual pressure edges the last edge is 0.01
                % which is not included in the satellite pressure values, we'll
                % include it to have edges to work with, then also add 0 to the
                % boxheight values so that we can take the difference.
                
                omi_zvec_p = interp1(log(gc_pvec), gc_zvec, log(omi_pres_edge(:,p)),'linear','extrap');
                omi_bxhght_p = diff(omi_zvec_p);
                
                % Also interpolate number density to OMI pressure levels,
                % number density is directly proportional to p:
                %   pV = nRT => n/V = p/RT
                % so we'll directly interpolate.
                omi_ndens_air_p = interp1(gc_pvec(1:end-1), gc_nair_vec, omi_pres(:,p), 'linear','extrap');
                
                % Finally interpolate the GC NO2 mixing ratio values (scaling
                % to straight mixing ratios, i.e. parts-per-part). I'll do this
                % in log-log space, which I've seen before. We will NOT
                % extrapolate this time, because I do not want there to be
                % values below the surface (so outside of the GC pressure
                % values, it will have the value of NaN)
                omi_no2_p = exp(interp1(log(gc_pvec(1:end-1)), log(gc_no2vec), log(omi_pres(:,p))));
                strat = omi_pres(:,p) < gc_pressure(a,b,floor(gc_tp(a,b,t)),t);
                omi_no2_p(strat) = nan;
                
                % Finally the actual calculation of columns is fairly straightforward: sum
                % over each box's partial column weighted by AK (after appropriate unit
                % conversions)
                omi_no2_columns(p) = nansum2((omi_no2_p * 1e-9) .* (omi_ndens_air_p * 1e-6) .* (omi_bxhght_p * 100) .* omi_aks(:,p));

                if xor(all(isnan(omi_aks(:,p))), omi_weights(p) == 0 || isnan(omi_weights(p)))
                    %warning('For a = %d, b = %d, t = %d, p = %d, the AK is all NaNs (bool=%d) but the weight is positive (isnan=%d, ==0=%d) or vice versa',a,b,t,p, all(isnan(omi_aks(:,p))), isnan(omi_weights(p)), omi_weights(p)==0)
                end
            end
            %fprintf('W%d: (%d,%d): Size omi_no2_columns = %s, Size omi_weights = %s, size omi_aks = %s\n',tstr.ID,a,b,mat2str(size(omi_no2_columns)),mat2str(size(omi_weights)), mat2str(size(omi_aks)));
            no2_columns(a,b,t) = nansum2(omi_no2_columns .* omi_weights)./nansum2(omi_weights);
        end
    end
end

end

function gc_prof = scale_gc_prof(gc_prof, prof_pres, no2nox_ratios)
% The first step is to interpolate the ratios to the model pressure, in
% log-log space again. Values outside the ratios' pressure range are set to
% 0, so that when exponentiated, they become 1, i.e. don't scale the
% profile where we don't have the data to do so.
E = JLLErrors;
dc3_ratio = exp(interp1(log(no2nox_ratios.dc3_pres), log(no2nox_ratios.dc3_ratio), log(prof_pres), 'linear', nan));
gc_ratio = exp(interp1(log(no2nox_ratios.gc_pres), log(no2nox_ratios.gc_ratio), log(prof_pres), 'linear', nan));
% We can't have constant extrapolation with linear interpolation, so
% instead we do this manually. However, for levels where both the
% observations and model are not available, set to 1.
dc3_nans = isnan(dc3_ratio);
dc3_v1 = find(~dc3_nans,1,'first');
dc3_vlast = find(~dc3_nans,1,'last');
gc_nans = isnan(gc_ratio);
gc_v1 = find(~gc_nans,1,'first');
gc_vlast = find(~gc_nans,1,'last');
for a=1:numel(dc3_nans)
    if dc3_nans(a) && gc_nans(a)
        dc3_ratio(a) = 1;
        gc_ratio(a) = 1;
    else
        if dc3_nans(a) && ~gc_nans(a)
            if a <= dc3_v1
                dc3_ratio(a) = dc3_ratio(dc3_v1);
            elseif a >= dc3_vlast
                dc3_ratio(a) = dc3_ratio(dc3_vlast);
            else
                E.notimplemented('NaN in dc3_ratio between the first and last non-NaN values')
            end
        elseif gc_nans(a) && ~dc3_nans(a)
            if a <= gc_v1
                gc_ratio(a) = gc_ratio(gc_v1);
            elseif a >= gc_vlast
                gc_ratio(a) = gc_ratio(gc_vlast);
            else
                E.notimplemented('NaN in gc_ratio between the first and last non-NaN values')
            end
        else
            % nothing to do if neither is a NaN
        end
    end
end
gc_prof = gc_prof .* dc3_ratio ./ gc_ratio;
end
