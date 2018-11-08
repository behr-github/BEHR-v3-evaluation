function [  ] = merge_hong_kong_output( root_dir, start_date, end_date, varargin )
%MERGE_HONG_KONG_OUTPUT Merges CMAQ and WRF data into standard wrfout files
%   MERGE_HONG_KONG_OUTPUT( ROOT_DIR, START_DATE, END_DATE ) Looks for
%   folders dated 'dd mmm yyyy' in ROOT_DIR and merges the hourly pressure,
%   temperature, and NO2 files in those folders along with LON.nc and
%   LAT.nc in ROOT_DIR into wrfout_cmaq files that follow the normal format
%   for wrfout or wrfout_subset files read by rProfile_WRF. Only data
%   between START_DATE and END_DATE are operated on.
%
%   Parameters:
%
%       'overwrite' - boolean (default false), controls if existing
%       wrfout_cmaq files should be overwritten.
%
%       'save_dir' - string that overrides the default save directory,
%       which is fullfile(behr_paths.wrf_profiles{1}, 'hk'). Files are
%       always organized into year and month subfolders within this
%       save_dir, whether it is given or the default is used.
%
%       'DEBUG_LEVEL' - scalar number that controls how much information is
%       printed to the terminal. Default is 2, 0 == no output.

E = JLLErrors;

p = advInputParser;
p.addFlag('check_symlinks');
p.addParameter('overwrite', false);
p.addParameter('save_dir', '');
p.addParameter('hours', 1:24);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.AdvResults;

check_symlinks_bool = pout.check_symlinks;
overwrite = pout.overwrite;
hours = pout.hours;
save_dir = pout.save_dir;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~islogical(overwrite) || ~isscalar(overwrite)
    E.badinput('The parameter "overwrite" must be a scalar logical value')
end

if isempty(save_dir)
    save_dir = fullfile(behr_paths.wrf_profiles{1}, 'hk');
    if ~exist(save_dir, 'dir')
        E.dir_dne('The default root save directory (%s) does not exist', save_dir)
    end
elseif ~ischar(save_dir)
    E.badinput('The parameter "save_dir" must be a string')
elseif ~exist(save_dir, 'dir')
    E.badinput('The save directory given (%s) does not exist', save_dir);
end

if ~isnumeric(hours) || ~isrow(hours) || any(hours < 1 | hours > 24 | mod(hours,1) ~= 0)
    E.badinput('The parameter "hours" must be a numeric row vector with integer values between 1 and 24');
end

if ~isnumeric(DEBUG_LEVEL) || ~isscalar(DEBUG_LEVEL)
    E.badinput('The parameter "DEBUG_LEVEL" must be a scalar number')
end

% LON.nc and LAT.nc should be present in the root directory. Rename the
% variables in the info structures so that NCWRITESCHEMA names them
% properly. There should only be one variable per file. Also, these have
% the LAY (layer) dimension, which they really don't need, because they
% aren't defined on a 4D grid.
lon_info = ncinfo(fullfile(root_dir, 'LON.nc'));
lon_info.Variables.Name = 'XLONG';
lon_info.Dimensions(3) = [];
lon_info.Variables.Dimensions(3) = [];
lon_array = ncread(lon_info.Filename, 'LON');

lat_info = ncinfo(fullfile(root_dir, 'LAT.nc'));
lat_info.Variables.Name = 'XLAT';
lat_info.Dimensions(3) = [];
lat_info.Variables.Dimensions(3) = [];
lat_array = ncread(lat_info.Filename, 'LAT');

% Map the variable file names to the desired final named in the merged
% files. The field name of the structure should be both the file name
% (before the _hh) and the variable name in the file. LON and LAT will be
% automatically mapped to XLONG and XLAT, since their files are present
% separately.
var_mapping = struct('NO2', 'no2',...
    'T', 'T',...
    'P', 'P',...
    'PB', 'PB',...
    'PH', 'PH',...
    'PHB', 'PHB');

% List extra variables that are created in this function here.
extra_vars = {'z_center'};

% These give the size of the CMAQ and WRF grids in the initial week of data
% that Hugo sent. They will be used to check if the input data is from that
% week and so the WRF and CMAQ grids need reconciled
hk_wrf_sz = [222, 162, 38];
hk_cmaq_sz = [98, 74, 28];

% If the grids are the size for the Hong Kong focused run, cut them
% down according to the information I got from Hugo Mak on 1 Nov 2017:
% For your first question, "I'd assume that the CMAQ levels are just the bottom 26 WRF layers, but could you verify that?", the answer is as follows:
%
%    "The vertical layers in WRF runs are as follows, those in red
%    (starred) are the 26 vertical layers (the 3rd variable, z) used in
%    CMAQ runs:
%            1.0000*, 0.9979*, 0.9956*, 0.9931*,
%            0.9904*, 0.9875*, 0.9844*, 0.9807*,
%            0.9763*, 0.9711*, 0.9649*, 0.9575*,
%            0.9488*, 0.9385*, 0.9263*, 0.9120*,
%            0.8951*, 0.8753*, 0.8521*, 0.8251*,
%            0.7937*, 0.7597, 0.7229*, 0.6833,
%            0.6410*, 0.5960, 0.5484, 0.4985*,
%            0.4467, 0.3934, 0.3393, 0.2850*,
%            0.2316, 0.1801, 0.1324, 0.0903*,
%            0.0542, 0.0241
%    (Totally 26 labelled layers)
%
%   Basically, the 3rd dimension of all WRF files is 38, and only the
%   following ones are useful: 1-20, 21, 23, 25, 28, 32, 36, so it
%   becomes 26 vertical layers, i.e. the new dimension is 26, in line
%   with CMAQ.
%
%   Regarding your 2nd question: "how the T, P, and PB arrays match up
%   with the NO2, lat, and lon arrays in space?", the answer is as
%   follows:
%
%   Before feeding into CMAQ, WRF outputs are converted into
%   CMAQ-readable format by MCIP.  MCIP cuts the domains of WRF outputs
%   into corresponding size of CMAQ according to our configuration, so
%   basically CMAQ only reads some parts of WRF outputs (in domain
%   size).  For WRF output files, it is 222 x 162 for (x, y) variables,
%   and only the followings are useful:
%
%   For x variable, from 20th to 117th entry, so it becomes dimension
%   of 98 For y variable, from 35th to 108th entry, so it becomes
%   dimension of 74"
%
% The other problem is that the CMAQ profiles don't go high enough for the
% WRF temperature profiles to reach the tropopause (though they are close).
% So we actually need the last two layers from the WRF variables as well. 
% We can either fill in those layers in the NO2 profile by assuming that
% they are equal to the top CMAQ box or by using a monthly mean GEOS-Chem
% profile to get the shape.
hk_x_cut = 20:117;
hk_y_cut = 35:108;
hk_z_cut = [1:20, 21, 23, 25, 28, 32, 36:38];

if numel(hk_x_cut) ~= hk_cmaq_sz(1)
    E.callError('inconsistent_cmaq_cut', 'The x_cut and first CMAQ length are inconsistent')
elseif numel(hk_y_cut) ~= hk_cmaq_sz(2)
    E.callError('inconsistent_cmaq_cut', 'The y_cut and second CMAQ length are inconsistent')
elseif numel(hk_z_cut) ~= hk_cmaq_sz(3)
    E.callError('inconsistent_cmaq_cut', 'The z_cut and third CMAQ length are inconsistent')
end

datevec = datenum(start_date):datenum(end_date);
for d=1:numel(datevec)
    subdir = datestr(datevec(d), 'dd mmm yyyy');
    this_save_dir = fullfile(save_dir, datestr(datevec(d), 'yyyy'), datestr(datevec(d), 'mm'));
    if ~exist(this_save_dir, 'dir')
        mkdir(this_save_dir)
    end
    for h=hours
        outfile = sprintf('wrfout_cmaq_%s_%02d-00-00', datestr(datevec(d), 'yyyy-mm-dd'), h-1);
        full_outfile = fullfile(this_save_dir, outfile);
        if exist(full_outfile, 'file')
            if ~overwrite
                if DEBUG_LEVEL > 0
                    fprintf('%s already exists\n', full_outfile)
                end
                continue
            else
                delete(full_outfile);
            end
        end
        
        % Copy the variable names and attributes
        copy_schema(fullfile(root_dir, subdir), h, full_outfile, check_symlinks_bool, datevec(d));
        
        % Copy the actual data, cutting down the WRF variables if
        % necessary.
        copy_data(fullfile(root_dir, subdir), h, full_outfile);
        
    end
end

    function copy_schema(input_dir, file_hour, output_file, check_links_bool, curr_date)
        fns = fieldnames(var_mapping);
        input_info = make_empty_struct_from_cell(fns);
        
        for f=1:numel(fns)
            input_file = make_input_name(input_dir, fns{f}, file_hour);
            
            if check_links_bool
                check_symlink(input_file, curr_date);
            end
            
            input_info.(fns{f}) = ncinfo(input_file);
            input_info.(fns{f}).Variables.Name = var_mapping.(fns{f});
        end
        
        input_info.LON = lon_info;
        input_info.LAT = lat_info;
        
        input_info = match_cmaq_wrf_schema_grids(input_info);
        
        % Copy a non-staggered variable for z_center
        input_info.z_center = make_new_schema(input_info.P, 'z_center', {'description', 'altitude derived from geopotential at model box center', 'units', 'm'});
        
        ncwriteschema(output_file, input_info.LON);
        ncwriteschema(output_file, input_info.LAT);
        for f=1:numel(fns)
            ncwriteschema(output_file, input_info.(fns{f}));
        end
        
        for i_var=1:numel(extra_vars)
            ncwriteschema(output_file, input_info.(extra_vars{i_var}));
        end
    end

    function copy_data(input_dir, file_hour, output_file)
        % Load the data
        fns = fieldnames(var_mapping);
        data_struct = make_empty_struct_from_cell(fns);
        for f=1:numel(fns)
            input_file = make_input_name(input_dir, fns{f}, file_hour);
            data_struct.(fns{f}) = ncread(input_file, fns{f});
        end
        
        data_struct.LON = lon_array;
        data_struct.LAT = lat_array;
        
        % Add the z_center variable - since we are subsetting PH and PHB
        % next, if we want z at the center of grid cells, we need to
        % compute it while we have all the model layers
        z = calculate_wrf_altitude({data_struct.PH, data_struct.PHB});
        data_struct.z_center = (z(:,:,1:end-1,:) + z(:,:,2:end,:))/2;
        
        % Match up the CMAQ and WRF outputs
        data_struct = match_cmaq_wrf_data_grids(data_struct);
        
        % Write each array to the output file
        ncwrite(output_file, 'XLONG', data_struct.LON);
        ncwrite(output_file, 'XLAT', data_struct.LAT);
        for f=1:numel(fns)
            ncwrite(output_file, var_mapping.(fns{f}), data_struct.(fns{f}));
        end
        
        for i_var=1:numel(extra_vars)
            this_var = extra_vars{i_var};
            ncwrite(output_file, this_var, data_struct.(this_var));
        end
    end

    function data = match_cmaq_wrf_data_grids(data)
        
        data_fns = fieldnames(data);
        bt_stagger_size = hk_wrf_sz;
        bt_stagger_size(3) = bt_stagger_size(3) + 1;
        for i_fn = 1:numel(data_fns)
            this_fn = data_fns{i_fn};
            n_dim = ndims(data.(this_fn));
            data_size = size(data.(this_fn));
            if isequal(data_size, hk_wrf_sz(1:n_dim)) || isequal(data_size, bt_stagger_size(1:n_dim))
                % Allow for the possibility of bottom_top being a staggered
                % variable. We're cutting that down the same way Hugo did,
                % which is to just take the bottom edge of each box. It
                % would make more sense to average the top and bottom edges
                % and get the box center eventually, but that would
                % necessitate some modification to read_wrf_preproc to
                % handle a "z_center" variable or something.
                if DEBUG_LEVEL > 1
                    fprintf('  Cutting down "%s" from WRF to CMAQ size\n', this_fn);
                end
                data.(this_fn) = data.(this_fn)(hk_x_cut, hk_y_cut, hk_z_cut);
            elseif isequal(data_size, hk_cmaq_sz(1:n_dim)) || (strcmpi(this_fn, 'no2') && isequal(data_size(1:2), hk_cmaq_sz(1:2)))
                % If the variable is NO2 then we need to ignore the third
                % dimensions b/c we're adding the two extra necessary
                % layers
                if DEBUG_LEVEL > 1
                    fprintf('  "%s" is already the proper size for CMAQ\n', this_fn);
                end
            else
                E.callError('unknown_var_size', 'The "%s" variable is neither the size of a WRF variable nor a CMAQ variable', this_fn)
            end
        end
        
        % Regardless of the size otherwise, we need to extend the CMAQ NO2
        % profiles two more layers. Hugo said that WRF and CMAQ share the
        % same model top, so that means as far as CMAQ is concerned, the
        % top two WRF layers should have the same [NO2] as the top CMAQ
        % layer.
        n_extra_layers = numel(hk_z_cut) - size(data.NO2,3);
        extra_no2_layers = repmat(data.NO2(:,:,end,:), 1, 1, n_extra_layers, 1);
        data.NO2 = cat(3, data.NO2, extra_no2_layers);
        
     
    end

    function schema = match_cmaq_wrf_schema_grids(schema)
        % First, we need to make sure only one dimension is unlimited.
        % We'll use the WRF Time dimension
        unlimited_ind = [schema.P.Dimensions.Name];
        if sum(unlimited_ind) == 1
            wrf_unlimited_dim = schema.P.Dimensions(unlimited_ind);
        else
            wrf_unlimited_dim = [];
            warning('No unlimited dimension found!');
        end
        
        % The way that "schema" is set up, is that each variable is
        % represented by a schema for a single variable file. So, e.g.
        % schema.P is as if ncinfo was called on a file that only had one
        % variable, P. (That's why below we can reference
        % schema.(fns{f}).Variables.Dimensions without a multiple-reference
        % error). This whole block is dealing with different "Time"
        % dimensions between WRF and CMAQ. Basically, since version 3
        % netCDF files can only have one unlimited dimension, any unlimited
        % dimension gets overwritten by the WRF "Time" dimension so that
        % only one unlimited dimension is defined across all the schema.
        %
        % In the WRF files cut down to match CMAQ, the "Time" dimension is
        % missing. In that case, we'll just remove the unlimited dimension
        % entirely, since BEHR assumes that each file only contains one
        % time anyway.
        fns = fieldnames(schema);
        for f=1:numel(fns)
            for i=1:numel(schema.(fns{f}).Dimensions)
                if schema.(fns{f}).Dimensions(i).Unlimited
                    % Found an unlimited dimension - overwrite it with the
                    % WRF "Time" dimension, if that exists. 
                    if ~isempty(wrf_unlimited_dim)
                        schema.(fns{f}).Dimensions(i) = wrf_unlimited_dim;
                        if schema.(fns{f}).Variables.Dimensions(i).Unlimited
                            schema.(fns{f}).Variables.Dimensions(i) = wrf_unlimited_dim;
                        else
                            E.callError('time_dim', 'Unlimited dimension is different in the root and variable dimensions for %s', schema.(fns{f}).Filename);
                        end
                    else
                        % If WRF files do not have an unlimited dimension,
                        % we need to remove the unlimited dimension from
                        % the CMAQ file and variables.
                        
                        % Matlab seems to make a distinction between a
                        % literal empty array and one saved to a variable,
                        % so we have to explicitly use the empty brackets
                        % to remove the dimension.
                        schema.(fns{f}).Dimensions(i) = [];
                        if schema.(fns{f}).Variables.Dimensions(i).Unlimited
                            schema.(fns{f}).Variables.Dimensions(i) = [];
                        else
                            E.callError('time_dim', 'Unlimited dimension is different in the root and variable dimensions for %s', schema.(fns{f}).Filename);
                        end
                    end
                end
            end
            
            % Also make everything 64 bit
            schema.(fns{f}).Format = '64bit';
        end

        if all_sizes_equal(schema)
            return
        elseif isequal(size_from_schema(schema.NO2) + [0 0 2], hk_cmaq_sz)
            % As long as the NO2 is the expected CMAQ size plus 2 for extra
            % levels added to get the temperature profile up into the
            % tropopause, go ahead and set all the other variables to be
            % the CMAQ size
            schema_fns = fieldnames(schema);
            for i_fn = 1:numel(schema_fns)
                this_fn = schema_fns{i_fn};
                if strcmpi(this_fn, 'no2')
                    flag_args = {'force'};
                else
                    flag_args = {};
                end
                schema.(this_fn) = set_schema_dims(schema.(this_fn), hk_wrf_sz, hk_cmaq_sz, flag_args{:});
            end
        else
            E.callError('undef_grid', 'No subsetting defined for CMAQ/WRF grid sizes of %s vs. %s', mat2str(size_from_schema(schema.NO2)), mat2str(size_from_schema(schema.P)));
        end
    end

end

function input_file = make_input_name(input_dir, var_name, file_hour)
input_file = fullfile(input_dir, sprintf('%s_%d.nc', var_name, file_hour));
end

function chk = all_sizes_equal(data_in)
% Were we given a structure of the data fields, or of the schema?
if isfield(data_in.NO2,'Variables')
    sz_fxn = @(x) size_from_schema(x);
else
    sz_fxn = @(x) size(x);
end
% Assume the sizes are unequal until we prove otherwise
chk = false;
sz = sz_fxn(data_in.NO2); % Assume that all the 3D variables are centered in the grid cells, and so have the same size as NO2
fns = fieldnames(data_in);
% Check each field. If the size doesn't match, return false. If we get
% through to the end, all fields' sizes must have matched, so return true.
for f=1:numel(fns)
    if ismember(fns{f}, {'LON','LAT'})
        if ~isequal(sz_fxn(data_in.(fns{f})), sz(1:2))
            return
        end
    else
        if ~isequal(sz_fxn(data_in.(fns{f})), sz)
            return
        end
    end
end

chk = true;

end

function sz = size_from_schema(schema)
E = JLLErrors;
sz = [schema.Dimensions.Length];
if ~isequal(sz, [schema.Variables.Dimensions.Length])
    E.callError('dim_mismatch', 'Lenghts of root dimensions is different from the lengths of %s dimensions in %s', schema.Filename, schema.Variables.Name);
end

% Remove trailing singleton dimensions
xx = find(sz > 1, 1, 'last');
sz = sz(1:xx);
end

function schema = set_schema_dims(schema, starting_size, desired_size, varargin)
E = JLLErrors;
p = advInputParser;
p.addFlag('force');
p.parse(varargin{:});
pout = p.Results;

force_to_cmaq_dims = pout.force;

if numel(schema.Dimensions) ~= numel(schema.Variables.Dimensions)
    E.callError('dim_mismatch', 'Number of root dimensions is different from the number of %s dimensions in %s', schema.Filename, schema.Variables.Name);
end

n_dims = numel(schema.Dimensions);

for i=1:n_dims
    if i > numel(desired_size)
        % If this dimension does not have a desired size, set it to 1,
        % since in Matlab an extra dimension at the end of the array is
        % considered to have length 1.
        sz_i = 1;
    else
        sz_i = desired_size(i);
    end
    
    
    
    if i > numel(starting_size)
        chk_sz_i = 1;
    else
        chk_sz_i = starting_size(i);
    end
    
    if force_to_cmaq_dims || schema.Dimensions(i).Length == sz_i || schema.Dimensions(i).Length == chk_sz_i
        % Whether it is the expected starting (WRF) size or target (CMAQ)
        % size, we can set it - this avoids an extra elseif. But if it is
        % neither of those lengths, then something unexpected is going on,
        % and we don't want to make any assumptions. We can override this
        % to force this to set the schema dims to the CMAQ dims, which we
        % usually need to do for NO2 b/c we're adding two levels that
        % aren't in the CMAQ output.
        schema.Dimensions(i).Length = sz_i;
    elseif (schema.Dimensions(i).Length == sz_i + 1 || schema.Dimensions(i).Length == chk_sz_i + 1) && regcmp(schema.Dimensions(i).Name, '.*stag')
        % When Hugo sent me the first draft of PH and PHB variables, they
        % were subset in such a way that the lower edge of each box was
        % given, but not necessarily the upper edge. This mimics that.
        schema.Dimensions(i).Name = strrep(schema.Dimensions(i).Name, 'stag', 'lower_edge');
        schema.Dimensions(i).Length = sz_i;
    elseif schema.Dimensions(i).Length > 1 && i == n_dims
        % LAT and LON will have the third dimension be time, which should
        % be length 1, in general, a dimension doesn't match but is only 1
        % long it's probably okay if it is the last dimension because this
        % variable has fewer dims than the others
        E.callError('schema_global_dim_error', 'The dimension "%s" from file %s is length %d; neither the expected WRF (%d) nor CMAQ (%d) length',...
            schema.Dimensions(i).Name, schema.Filename, schema.Dimensions(i).Length, chk_sz_i, sz_i)
    end
    
    
    if force_to_cmaq_dims || schema.Variables.Dimensions(i).Length == sz_i || schema.Variables.Dimensions(i).Length == chk_sz_i
        schema.Variables.Dimensions(i).Length = sz_i;
    elseif (schema.Variables.Dimensions(i).Length == sz_i + 1 || schema.Variables.Dimensions(i).Length == chk_sz_i + 1) && regcmp(schema.Variables.Dimensions(i).Name, '.*stag')
        schema.Variables.Dimensions(i).Name = strrep(schema.Dimensions(i).Name, 'stag', 'lower_edge');
        schema.Variables.Dimensions(i).Length = sz_i;
    elseif schema.Variables.Dimensions(i).Length > 1
        E.callError('schema_global_dim_error', 'The dimension "%s" in variable "%s" from file %s is length %d; neither the expected WRF (%d) nor CMAQ (%d) length',...
            schema.Variables.Dimensions(i).Name, schema.Variables.Name, schema.Filename, schema.Variables.Dimensions(i).Length, chk_sz_i, sz_i)
    end
end
end

function base_schema = make_new_schema(base_schema, var_name, attribute_changes)
% Since Matlab does copy-on-write, we can just edit the given schema and
% return it
E = JLLErrors;
base_schema.Variables.Name = var_name;
if isempty(base_schema.Variables.Attributes)
    att_names = {};
else
    att_names = {base_schema.Variables.Attributes.Name};
end

for i_att = 1:2:numel(attribute_changes)
    xx_att = strcmp(att_names, attribute_changes{i_att});
    if sum(xx_att) < 1
        fprintf('No attribute named "%s" exists in the schema, creating\n', attribute_changes{i_att});
        base_schema.Variables.Attributes(end+1).Name = attribute_changes{i_att};
        base_schema.Variables.Attributes(end).Value = attribute_changes{i_att+1};
    elseif sum(xx_att) > 1
        E.callError('multiple_atts_matched', 'Multiple attributes named "%s" exist in the schema', attribute_changes{i_att});
    else
        base_schema.Variables.Attributes(xx_att).Value = attribute_changes{i_att+1};
    end
end
end

function check_symlink(link_file, path_date)
% Check that a symlink is pointing to a file for the current date. This
% assumes that the file pointed to by the link resides in a folder named as
% the date; i.e. the next to last path component is a folder with that
% date. 

E = JLLErrors;
% For some reason, readlink doesn't work if you use ~ for the home
% directory
homedir = getenv('HOME');
link_file = strrep(link_file, '~', homedir);
[stat, link_path] = system(sprintf('readlink "%s"', link_file));
if stat ~= 0
    return
end

tmp_path = strsplit(link_path, '/');
date_folder = tmp_path{end-1};

if datenum(date_folder, 'dd mmm yyyy') ~= datenum(path_date)
    E.callError('check_symlink:wrong_date', 'File (%s) linked to wrong date (%s instead of %s)', link_file, date_folder, datestr(path_date, 'dd mmm yyyy'));
else
    fprintf('Link %s correct\n', link_file);
end

end
