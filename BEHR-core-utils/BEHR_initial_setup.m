function [  ] = BEHR_initial_setup( )
%BEHR_INITIAL_SETUP Perform the actions necessary to set up BEHR to run
%   There are several steps to prepare the BEHR algorithm to run:
%       1) Clone the necessary repositories from GitHub
%       2) Identify the paths where OMI, MODIS, GLOBE, and WRF-Chem data
%       reside.
%       3) Add the necessary code paths to the Matlab search path
%       4) Build the "omi" Python package included in the BEHR-PSM-Gridding
%       repository.
%   If you're reading this, then you should have already done #1. This
%   function will automatically perform steps 2-4. You may be prompted for
%   each step, so you must run this function interactively (i.e. it cannot
%   be run as part of a batch job on a cluster).
%
%   NOTE TO DEVELOPERS: This function should only rely on built-in Matlab
%   functions and not any functions included as part of the BEHR code or
%   related repositories, since it is intended to be run before any of the
%   BEHR code has been added to the Matlab path.

[behr_utils_repo_path, myname] = fileparts(mfilename('fullpath'));
paths_filename = 'behr_paths.m';
template_filename = 'behr_paths_template.m';
constants_path = fullfile(behr_utils_repo_path,'Utils','Constants');
paths_file = fullfile(constants_path,paths_filename);
template_file = fullfile(constants_path,template_filename);
        
make_paths_file();
add_code_paths();
build_omi_python();

    function make_paths_file
        % Check that behr_paths.m does not already exist anywhere. If there's one
        % in the right place, then ask if we want to overwrite it. If there's one
        % somewhere else, the user will need to move or delete it so that we don't
        % get two behr_paths.m files.
        
        if exist(paths_file,'file')
            user_ans = input(sprintf('A %s file already exists in BEHR/Utils/Constants. Skip the paths setup and use existing?: ', paths_filename), 's');
            if strcmpi(user_ans,'Yes') || strcmpi(user_ans, 'y')
                % Previously I had this function delete behr_path.m.
                % However, since behr_paths is called later, that led to
                % issue where Matlab wouldn't recompile the class or
                % something, and try to use the old class. Having the user
                % manually delete it avoids this problem.
                return;                
            else
                error('behr_setup:file_exists', 'You must manually delete the existing behr_paths.m file at %s before this function will recreate that file', constants_path);
            end
        elseif ~isempty(which(paths_filename))
            error('behr_setup:file_exists','%1$s exists on your Matlab search path, but at %2$s, not at %3$s. Please delete or move that version to %3$s.',paths_filename, which(paths_filename), constants_path);
        end
        
        % This defines all the paths that should exist in that file. The
        % field name is the name the property will be given, the 'comment'
        % subfield is the explanatory comment, and the 'default' subfield
        % is what the default will be. If either field is omitted, that is
        % fine. A generic comment will be inserted and a blank path will be
        % inserted. Including the subfield "no_quote" at all will cause it
        % not to automatically quote the default path, which is useful if
        % you need to allow for multiple paths in a cell array. Including
        % the the field "isfile" (whatever value is set to it) will cause
        % behr_paths.ValidatePaths() to check if a file exists at that
        % path, rather than a directory. Including the field "is_code_dir"
        % will make behr_paths aware that the directory given contains
        % code, and should be added to the Matlab path on request. If its
        % value is the string 'norecurse', only that directory (and not its
        % subfolders) will be added to the path. Likewise, if "is_pypath"
        % is a field, it will be aware that this path contains Python code
        % and should be added to the Python path by SetPythonPath().
        
        sat_file_server = '128.32.208.13';
        wrf_file_server = 'cohenwrfnas.dyn.berkeley.edu';
        
        if ismac
            sat_folder = fullfile('/Volumes','share-sat');
            wrf1_folder = fullfile('/Volumes', 'share-wrf1','BEHR-WRF');
            wrf2_folder = fullfile('/Volumes', 'share-wrf2','BEHR-WRF');
        elseif isunix
            sat_folder = fullfile('/mnt','share-sat');
            wrf1_folder = fullfile('/mnt', 'share-wrf1','BEHR-WRF');
            wrf2_folder = fullfile('/mnt', 'share-wrf2','BEHR-WRF');
        elseif ispc
            drive_letter_request = 'Enter the drive letter (just the letter, not the ":\" that "%s" is mounted as';
            sat_folder = strcat(input(sprintf(drive_letter_request, 'share-sat'), 's'), ':');
            wrf1_folder = strcat(input(sprintf(drive_letter_request, 'share-wrf1','BEHR-WRF'), 's'), ':');
            wrf2_folder = strcat(input(sprintf(drive_letter_request, 'share-wrf2','BEHR-WRF'), 's'), ':');
        end
            
        
        % Local repos/folders
        paths.behr_core.comment = 'The directory of the BEHR-core repository. May be cloned from https://github.com/CohenBerkeleyLab/BEHR-core';
        paths.behr_core.default = fullfile(behr_utils_repo_path, '..', 'BEHR-core');
        paths.behr_core.is_code_dir = true;
        paths.behr_utils.comment = 'The directory of the BEHR-core-utils repository. May be cloned from https://github.com/CohenBerkeleyLab/BEHR-core-utils';
        paths.behr_utils.default = behr_utils_repo_path;
        paths.behr_utils.is_code_dir = true;
        paths.utils.comment = 'The directory of the general Matlab-Gen-Utils repository (not the BEHR-core-utils repo). May be cloned from https://github.com/CohenBerkeleyLab/Matlab-Gen-Utils';
        paths.utils.default = fullfile(behr_utils_repo_path, '..', 'Matlab-Gen-Utils');
        paths.utils.is_code_dir = true;
        paths.amf_tools_dir.comment = 'The AMF_tools directory in the BEHR-core-utils repository on your computer. It should contain the files damf.txt and nmcTmpYr.txt';
        paths.amf_tools_dir.default = fullfile(behr_utils_repo_path, 'AMF_tools');
        % do not need to explicitly add AMF_tools to path, since it is in
        % the behr_utils repo
        paths.psm_dir.comment = 'The PSM Gridding repository. It should contain the files PSM_Main.py and psm_wrapper.m. May be cloned from https://github.com/CohenBerkeleyLab/BEHR-PSM-Gridding';
        paths.psm_dir.default = fullfile(behr_utils_repo_path, '..', 'BEHR-PSM-Gridding');
        paths.psm_dir.is_code_dir = 'norecurse';
        paths.psm_dir.is_pypath = true;
        paths.python_interface.comment = 'The MatlabPythonInterface repository. May be cloned from https://github.com/CohenBerkeleyLab/MatlabPythonInterface';
        paths.python_interface.default = fullfile(behr_utils_repo_path, '..', 'MatlabPythonInterface');
        paths.python_interface.is_code_dir = true;
        paths.wrf_utils.comment = 'The WRF_Utils repository which should contain the function convert_wrf_temperature.m May be cloned from https://github.com/CohenBerkeleyLab/WRF_Utils';
        paths.wrf_utils.default = fullfile(behr_utils_repo_path, '..', 'WRF_Utils');
        paths.wrf_utils.is_code_dir = true;
        
        % Matlab file folders
        paths.sp_mat_dir.comment = sprintf('The default path where OMI_SP_*_yyyymmdd.mat files will be saved and read from. It should have subdirectories for each region to be produced (e.g. "us" - must be lower case). For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.sp_mat_dir.default = fullfile(sat_folder, 'SAT', 'BEHR', 'SP_Files');
        paths.behr_mat_dir.comment = sprintf('The default root path where OMI_BEHR_*_yyyymmdd.mat files will be saved and read from. It should have subdirectories for each region to be produced and within each region directories "daily" and "monthly". For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.behr_mat_dir.default = fullfile(sat_folder, 'SAT', 'BEHR', 'BEHR_Files');
        
        % OMI and ancillary data folders
        paths.omno2_dir.comment = sprintf('This should contain folders organized by year and month with OMI-Aura_L2-OMNO2 files in them. For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.omno2_dir.default = fullfile(sat_folder, 'SAT', 'OMI', 'OMNO2', 'version_3_3_0');
        paths.ompixcor_dir.comment = sprintf('This should contain folders organized by year and month with OMI-Aura_L2-OMNPIXCOR files in them. For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.ompixcor_dir.default = fullfile(sat_folder, 'SAT', 'OMI', 'OMPIXCOR', 'version_003');
        paths.myd06_dir.comment = sprintf('This should contain folders for each year with MYD06_L2 files in them. For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.myd06_dir.default = fullfile(sat_folder, 'SAT', 'MODIS', 'MYD06_L2');
        paths.mcd43d_dir.comment = sprintf('This should contain folders for each year with MCD43D* files in them. For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.mcd43d_dir.default = fullfile(sat_folder, 'SAT', 'MODIS', 'MCD43D');
        paths.modis_land_mask.comment = sprintf('This is the "Land_Water_Mask_7Classes_UMD file, available from ftp://rsftp.eeos.umb.edu/data02/Gapfilled/ (as of 21 Sept 2017). For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.modis_land_mask.default = fullfile(sat_folder, 'SAT', 'MODIS', 'Land_Water_Mask_7Classes_UMD.hdf');
        paths.modis_land_mask.isfile = true;
        paths.globe_dir.comment = sprintf('This is the folder with the GLOBE database, available from https://www.ngdc.noaa.gov/mgg/topo/gltiles.html (as of 21 Sept 2017). It should contain files a10g through p10g and their .hdr files. For UC Berkeley users, is it on the file server at %s which should be mounted on your computer.',sat_file_server);
        paths.globe_dir.default = fullfile(sat_folder, 'SAT', 'BEHR','GLOBE_Database');

        paths.website_staging_dir.comment = sprintf('The directory where data should be staged before being put on the folders visible to the website. Also on the file server at %s which should be mounted on your computer.', sat_file_server);
        paths.website_staging_dir.default = fullfile(sat_folder,'SAT','BEHR','WEBSITE','staging');
        
        % WRF data is spread across multiple volumes locally. It's just too
        % big.
        paths.wrf_monthly_profiles.comment = sprintf('The path that contains the WRF_BEHR*.nc monthly profile files. This should be on the file server at %s.', wrf_file_server);
        paths.wrf_monthly_profiles.default = fullfile(wrf1_folder,'MonthlyProfiles');
        
        paths.wrf_profiles.comment = sprintf('Add all paths that contain WRF profiles. These should be folders that are organized by year and month, with wrfout_d01 files in them. These will be found on the file server at %s, all volumes must be mounted on your computer.', wrf_file_server);
        paths.wrf_profiles.default = sprintf('{''%s'', ''%s''}', fullfile(wrf1_folder, 'Outputs'), fullfile(wrf2_folder, 'Outputs'));
        paths.wrf_profiles.no_quote = true;
        
        %%%%%%%%%%%%%%%%%%
        % Write the file %
        %%%%%%%%%%%%%%%%%%
        
        fid_template = fopen(template_file, 'r');
        fid_new = fopen(paths_file, 'w');
        
        tline = fgetl(fid_template);
        [~,paths_classname] = fileparts(paths_filename);
        [~,template_classname] = fileparts(template_filename);
        while ischar(tline)
            tline = strrep(tline, template_classname, paths_classname);
            if strcmp(strtrim(tline), '%#PATHS')
                write_paths(paths, fid_new);
            elseif strcmp(strtrim(tline), '%#ISFILE')
                write_is_file_struct(paths, fid_new);
            elseif strcmp(strtrim(tline), '%#ISCODEDIR')
                write_iscodedir_structs(paths, fid_new);
            else
                fprintf(fid_new, '%s\n', tline);
            end
            tline = fgetl(fid_template);
        end
        
        fclose(fid_template);
        fclose(fid_new);
        
        fprintf('\n!!! Defaults %s file created at %s. Review it and edit the paths as needed. !!!\n\n', paths_filename, paths_file);
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add BEHR paths to Matlab search path %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function add_code_paths
        fprintf('\n');
        fprintf('Trying to add BEHR code paths to Matlab path...\n');
        addpath(fullfile(behr_utils_repo_path, 'Utils', 'Constants'));
        try
            behr_paths.AddCodePaths();
            input('Press ENTER to continue','s');
        catch err
            if strcmpi(err.identifier, 'path_setup:bad_paths')
                fprintf('%s\n', err.message);
                fprintf('Correct the bad paths in %s and then run behr_paths.AddCodePaths()\n', paths_file);
                input('I will still try to finish the initial setup (Press ENTER to continue)', 's');
            else
                rethrow(err)
            end
        end
    end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build the omi Python package %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function build_omi_python
        fprintf('\n');
        % Get the version of python that Matlab will use
        [~, py_exe] = pyversion;
        build_omi_cmd = sprintf('%s setup.py build', py_exe);
        install_omi_cmd = sprintf('%s setup.py install --user', py_exe);
        psm_dir = fullfile(behr_paths.psm_dir, 'omi');
        
        if exist(psm_dir, 'dir')
            fprintf('BEHR relies on the omi Python package. This needs to be built and installed on your Python path.\n');
            fprintf('I can try to do this for you; I would run:\n\n\t%s\n\nfollowed by\n\n\t%s\n\nin %s\n\n', build_omi_cmd, install_omi_cmd, psm_dir);
            fprintf('I am using the Python executable that Matlab will use. You can change this with the pyversion() function.\n\n');
            fprintf('If you change Python versions for Matlab, you should rebuild the omi package. To do a completely clean build,\ndelete the "build" folder in %s\n\n', psm_dir);
            fprintf('Should I try to build the omi package? If not, you can do it manually later.\n');
            if strcmpi(input('  Enter y to build, anything else to skip: ', 's'), 'y')
                fprintf('Trying to build the omi package...\n');
                oldwd = cd(psm_dir);
                try
                    [build_status, build_result] = system(build_omi_cmd);
                catch err
                    cd(oldwd);
                    rethrow(err);
                end
                
                if build_status == 0
                    fprintf('%s\n\n', build_result)
                else
                    fprintf('Build failed. Output from build command:\n%s\n', build_result);   
                end
                
                try
                    [install_status, install_result] = system(install_omi_cmd);
                catch err
                    cd(oldwd);
                    rethrow(err);
                end
                cd(oldwd);
                
                if install_status == 0
                    fprintf('%s\n\n', install_result);
                    fprintf('Build appears to be successful.\n');
                else
                    fprintf('Build failed. Output from install command:\n%s\n', install_result);
                end
            else
                fprintf('\nTo build the omi package yourself, execute "%s" followed by "%s"\nin %s\n\n', build_omi_cmd, install_omi_cmd, psm_dir);
                fprintf('Note that your PATH and PYTHONPATH environmental variables may differ in the Terminal vs. the GUI Matlab.\n');
                fprintf('If you experience difficultly with the PSM package after building in the Terminal, try building from Matlab instead.\n');
            end
        else
            fprintf('The PSM directory (%s) in behr_paths is invalid. Fix that and rerun BEHR_initial_setup, or build the omi package manually.\n', behr_paths.psm_dir);
        end
    end

    
end

function str_out = wrap_comment(comment, indent_level)
% Assume each line is 76 characters long, allow 4 spaces for tabs and 2
% characters for the "% " at the beginning
nchars = 76 - 4*indent_level - 2;
tabs = repmat('\t', 1, indent_level);
i = 1;
str_out = '';
while i < length(comment)
    % Find the last space before the end of the line
    if length(comment) - i > nchars
    i2 = min(length(comment), i+nchars);
    j = strfind(comment(i:i2), ' ');
    j = j(end);
    else
        j = length(comment) - i + 1;
    end
    str_out = [str_out, tabs, '%% ', strtrim(comment(i:i+j-1)), '\n'];
    %str_out{end+1} = sprintf('%% %s\n',strtrim(comment(i:i+j-1)));
    i = i+j;
end
end

function write_paths(paths, fid)
fns = fieldnames(paths);
for a=1:numel(fns)
    if isfield(paths.(fns{a}), 'comment');
        this_comment = wrap_comment(paths.(fns{a}).comment, 2);
    else
        this_comment = wrap_comment(fns{a}, 2);
    end
    fprintf(fid, this_comment);
    
    if isfield(paths.(fns{a}), 'default')
        this_default = paths.(fns{a}).default;
    else
        this_default = '';
    end
    
    if isfield(paths.(fns{a}), 'no_quote')
        fprintf(fid, '\t\t%s = %s;\n\n', fns{a}, this_default);
    else
        fprintf(fid, '\t\t%s = ''%s'';\n\n', fns{a}, this_default);
    end
end
end

function write_is_file_struct(paths, fid)
is_file_fxn = @(path_struct) isfield(path_struct, 'isfile');
write_bool_structure(paths, fid, 'is_field_file', is_file_fxn);
end

function write_iscodedir_structs(paths, fid)
iscodedir_fxn = @(path_struct) isfield(path_struct, 'is_code_dir');
do_genpath_fxn = @(path_struct) isfield(path_struct, 'is_code_dir') && ~strcmpi(path_struct.is_code_dir, 'norecurse');
is_pypath_fxn = @(path_struct) isfield(path_struct, 'is_pypath');

write_bool_structure(paths, fid, 'is_code_dir', iscodedir_fxn);
fprintf(fid, '\n');
write_bool_structure(paths, fid, 'do_genpath', do_genpath_fxn);
fprintf(fid, '\n');
write_bool_structure(paths, fid, 'is_pypath', is_pypath_fxn);
end

function write_bool_structure(paths, fid, struct_name, bool_fxn)
% Write a structure containing each of the paths as fields with a boolean
% value. "paths" must be the paths structure, fid an open file identifier,
% struct_name the name you want the struct to have a bool_fxn a handle to a
% functions that, when given the value of a field in paths returns the
% boolean value you want stored in the struct for that field.
bool_struct = struct;
fns = fieldnames(paths);
for a=1:numel(fns)
    if bool_fxn(paths.(fns{a}))
        bool_struct.(fns{a}) = 'true';
    else
        bool_struct.(fns{a}) = 'false';
    end
end

fprintf(fid, '\t\t%s = struct(', struct_name);
for a=1:numel(fns)
    if a < numel(fns)
        fprintf(fid,'''%s'', %s,...\n\t\t\t', fns{a}, bool_struct.(fns{a}));
    else
        fprintf(fid,'''%s'', %s);\n', fns{a}, bool_struct.(fns{a}));
    end
end
end
