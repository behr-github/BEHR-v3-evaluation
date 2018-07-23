function [  ] = gcstruct2ncdf( ncfile, varargin )
%GCSTRUCT2NCFILE Save a GEOS-Chem structure as a netCDF file.
%   Detailed explanation goes here
[glon,glat] = geos_chem_centers('2x25');
xx = find(strcmpi(varargin,'lon'));
if ~isempty(xx)
    glon = varargin{xx+1};
    varargin(xx:xx+1) = [];
end
xx = find(strcmpi(varargin,'lat'));
if ~isempty(xx)
    glat = varargin{xx+1};
    varargin(xx:xx+1) = [];
end
xx = find(strcmpi(varargin,'do_tedge'));
if ~isempty(xx)
    do_tedge = true;
    varargin(xx) = [];
else
    do_tedge = false;
end



xx = strcmpi(varargin,'overwrite');
if any(xx)
    overwrite = true;
    varargin(xx) = [];
else
    overwrite = false;
end

cmode = netcdf.getConstant('CLOBBER');
if exist(ncfile, 'file')
    if overwrite || ask_yn(sprintf('File %s exists. Overwrite?', ncfile))
        % nothing to do
    else
        error('io:file_exists', 'File %s exists, aborting', ncfile)
    end
end          

% Fields that should be the same in all input structures
common_fields = {'modelName', 'tVec', 'modelRes'};
first_struct = true;
ncid = -1;

cmode = bitor(cmode, netcdf.getConstant('NETCDF4'));
try
    for a=1:numel(varargin)
        for b=1:numel(varargin{a})
            GC = varargin{a}(b);
            if first_struct
                % Write the common fields, as global attributes or variables as
                % needed
                FirstGC = GC;
                ncid = netcdf.create(ncfile, cmode);
                
                DimIDs.lon = netcdf.defDim(ncid, 'lon', numel(glon));
                DimIDs.lat = netcdf.defDim(ncid, 'lat', numel(glat));
                DimIDs.model_level_47 = netcdf.defDim(ncid, 'model_level_47', 47);
                DimIDs.time = netcdf.defDim(ncid, 'time', numel(GC.tVec));
                if do_tedge
                    DimIDs.time_edge = netcdf.defDim(ncid, 'time_edge', numel(GC.tEdge));
                end
                
                write_var_with_atts(ncid,'lon',glon,make_dims(glon), 'degrees', 'Longitude of GEOS-Chem grid cell center (west is negative)','','');
                write_var_with_atts(ncid,'lat',glat,make_dims(glat), 'degrees', 'Latitude of GEOS-Chem grid cell center (south is negative)','','');
                write_var_with_atts(ncid,'time',GC.tVec - datenum('1985-01-01'),make_dims(GC.tVec), 'days since midnight 1 Jan 1985', 'Model date as a number of days since midnight 1 Jan 1985','','');
                if do_tedge
                    write_var_with_atts(ncid,'time_edge',GC.tEdge - datenum('1985-01-01'),make_dims(GC.tEdge),'days since midnight 1 Jan 1985','Edges of model averaging periods (assumes adjacent periods)','','');
                end
                netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'modelName',GC.modelName);
                netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'),'modelRes',GC.modelRes);
                
                first_struct = false;
            else
                for c=1:numel(common_fields)
                    if ~isequaln(GC.(common_fields{c}), FirstGC.(common_fields{c}))
                        error('common_field:not_equal', 'Common field %s not equal in structure %d index %d and first structure', common_fields{c}, a, b);
                    end
                end
            end
            gcvar = regexprep(GC.fullName, '\W', ''); % Remove everything but a-z, 0-9, and _
            if ~isempty(GC.fullCat)
                gccat = regexprep(GC.fullCat, '\W', '');
                gcvar = [gccat,'-',gcvar]; % and if a category is given, prepend it to the variable name (this deals with the production variables all being called "NO tracer")
            end
            dims_cell = make_dims(GC.dataBlock);
            write_var_with_atts(ncid, gcvar, GC.dataBlock, dims_cell, GC.dataUnit, GC.description, GC.fullCat, GC.fullName);
        end
    end
catch err
    if ncid >= 0
        netcdf.close(ncid);
    end
    rethrow(err);
end

netcdf.close(ncid);

    function dims = make_dims(dataBlock)
        dims = zeros(1,ndims(dataBlock));
        sz = size(dataBlock);
        sz_matches = {numel(glon), DimIDs.lon;
            numel(glat), DimIDs.lat;
            numel(FirstGC.tVec), DimIDs.time;
            47, DimIDs.model_level_47};
        if do_tedge
            sz_matches = cat(1, sz_matches, {numel(FirstGC.tEdge), DimIDs.time_edge});
        end
        dim_len_1 = sz == 1;
        for i=1:numel(sz)
            if dim_len_1(i)
                continue
            end
            xx = sz(i) == cat(1,sz_matches{:,1});
            if sum(xx) == 1
                dims(i) = sz_matches{xx,2};
            elseif i==3
                % The production arrays have variable numbers of levels. So
                % we assume that (1) singleton dimensions have been
                % removed, so that skipping dimensions of length 1 is still
                % appropriate, and (2) that if more than one level is
                % present, it will be in the third dimension. We'll define
                % a new dimension if necessary.
                fprintf('For variable of size %s, assuming that the third dimension is model level\n', mat2str(sz));
                fn = sprintf('model_level_%02d',sz(i));
                if ~isfield(DimIDs, fn)
                    DimIDs.(fn) = netcdf.defDim(ncid, fn, sz(i));
                end
                dims(i) = DimIDs.(fn);
            else
                error('ncdim:cannot_id','Cannot identify dimension length of %d (variable dimensions are %s)', sz(i), mat2str(sz));
            end
        end
        
        dims = dims(~dim_len_1);
    end

end

function write_var_with_atts(ncid, gcvar, val, dims, units, description, gc_cat, gc_name)
if ~isnumeric(val)
    error('ncval:not_numeric', 'VAL must be a numeric type');
end
val = single(val);

varid = netcdf.defVar(ncid, gcvar, 'NC_FLOAT', dims);
netcdf.defVarDeflate(ncid, varid, true, true, 1);
netcdf.putVar(ncid, varid, val);
netcdf.putAtt(ncid, varid, 'Units', units);
netcdf.putAtt(ncid, varid, 'Description', description);
netcdf.putAtt(ncid, varid, 'Diagnostic_Category', gc_cat);
netcdf.putAtt(ncid, varid, 'Tracer_Name', gc_name);

end
