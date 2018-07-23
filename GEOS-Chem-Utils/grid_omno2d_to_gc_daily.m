function gridded_daily_omno2d = grid_omno2d_to_gc_daily(gc_loncorn, gc_latcorn)
E = JLLErrors;

omno2d_path = '/global/home/users/laughner/myscratch/SAT/OMI/OMNO2d';
if ~exist(omno2d_path,'dir')
    E.dir_dne('omno2d_path')
end
[omno2d_lon, omno2d_lat] = omno2d_centers;



gridded_daily_omno2d = nan([size(gc_loncorn)-1, 366]);

parfor d=1:366
    E=JLLErrors;
    sdate = modis_day_to_date(d,2012);
    
    fname = sprintf('OMI-Aura_L3-OMNO2d_%04dm%02d%02d_v003*.he5', year(sdate), month(sdate), day(sdate));
    F = dir(fullfile(omno2d_path, sprintf('%04d',year(sdate)),fname));
    if isempty(F)
        error('gridded_omno2d:file_not_found','The file %s cannot be found',fname);
    elseif numel(F) > 1
        error('gridded_omno2d:too_many_files','The filespec %s returns multiple options',fname);
    end
    
    hinfo = h5info(fullfile(omno2d_path, sprintf('%04d',year(sdate)), F(1).name));
    todays_no2 = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2TropCloudScreened'));
    todays_weight = h5read(hinfo.Filename, h5dsetname(hinfo,1,2,1,1,'Weight'));
    
    % Ignore negative values - they are either fill values or unphysical
    todays_weight(todays_no2<0) = nan;
    
    gridded_daily_omno2d(:,:,d) = grid_omno2d_to_gc(todays_no2, omno2d_lon, omno2d_lat, gc_loncorn, gc_latcorn, todays_weight);
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
if ~exist('omno2d_weights','var')
    omno2d_weights = ones(size(omno2d_no2));
end

gridded_omno2d = nan(size(gc_loncorn)-1); % one smaller to switch from # of corners to # of grid cells

% Loop over GC grid cells and find which OMNO2d data belongs in it
for a=1:size(gridded_omno2d,1)
    for b=1:size(gridded_omno2d,2)
        xx = omno2d_lon >= gc_loncorn(a,b) & omno2d_lon < gc_loncorn(a+1,b+1) & omno2d_lat >= gc_latcorn(a,b) & omno2d_lat < gc_latcorn(a+1,b+1);
        gridded_omno2d(a,b) = nansum2(omno2d_no2(xx) .* omno2d_weights(xx)) / nansum2(omno2d_weights(xx));
    end
end
end
