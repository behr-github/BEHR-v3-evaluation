classdef misc_behr_v3_validation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties(Constant = true)
        behr_v2_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_v2-1C';
        wrf_v2_dir = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
        gcas_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw/DISCOVER-AQ_TX/B200/GCAS-SAO';
        
        plot_colors = struct('aircraft', struct('raw', [0.5 0.5 0.5], 'avg', 'k'),... % black and grey for the aircraft data
            'v2', struct('raw', 'g', 'avg', [0 0.5 0]),... % greens for version 2
            'monthly', struct('raw', [1 0.75 0], 'avg', 'r'),... % red and orange for the monthly profiles
            'daily', struct('raw', 'c', 'avg', 'b')); % blue and cyan for the daily profiles
        
        plot_markers = struct('aircraft', struct('raw', '.', 'avg', 'o', 'filled', false),...
            'v2', struct('raw', '.', 'avg', '+', 'filled', false),...
            'monthly', struct('raw', '.', 'avg', 's', 'filled', true),...
            'daily', struct('raw', '.', 'avg', 'd', 'filled', false));
        
        plot_wrf_prof_legend_names = struct('aircraft', 'Aircraft', 'v2', 'V2', 'monthly', 'Monthly', 'daily', 'Daily');
        
        profile_extend_methods = {'wrf','geos','extrap'};
        
        discover_campaigns = {'discover_md','discover_ca','discover_tx','discover_co'};
    end
    
    
    
    methods(Static = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Property-like Methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        function value = my_dir()
            value = fileparts(mfilename('fullpath'));
        end
        
        function value = validation_root_dir()
            value = fullfile(misc_behr_v3_validation.my_dir, 'Workspaces', 'Validation');
        end
        
        function value = gc_data_path()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'GEOS-Chem-Monthly-Data');
        end
        
        function value = profile_comp_file(extend_method)
            E = JLLErrors;
            narginchk(1,1);
            extend_method = opt_ask_multichoice('Which profile extension method to use?', misc_behr_v3_validation.profile_extend_methods, extend_method, 'extend_method', 'list', true);
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', sprintf('profile-structs-%s.mat', extend_method));
        end
        
        function value = pandora_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'pandora-structs.mat');
        end
        
        function value = gcas_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'gcas-structs.mat');
        end
        
        function value = gcas_vec_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison', 'gcas-vec-structs.mat');
        end
        
        function value = scd_comp_dir()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'SCD-Wind-Comparisons');
            if ~exist(value, 'dir')
                mkdir(value)
            end
        end
        
        function value = scd_comp_file(is_partial)
            if ~exist('is_partial', 'var') || ~is_partial
                partial_str = '';
            else
                partial_str = '_partial';
            end
            value = fullfile(misc_behr_v3_validation.scd_comp_dir, sprintf('scd_daily_%s%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), partial_str));
        end
        
        function value = wrf_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'Profile-Comparison', 'wrf-match-structs.mat');
        end
        
        function value = pres_comp_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'Error-Analysis', 'pres-comp.mat');
        end
        
        function value = acarreta_cldpres_error_raw_file()
            value = fullfile(misc_behr_v3_validation.validation_root_dir, 'Error-Analysis', 'AcarrataCldPresUncert.txt');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Interactive methods for getting parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function extend_method = get_profile_extend_method(extend_method)
            allowed_methods = misc_behr_v3_validation.profile_extend_methods;
            if nargin < 1 || isempty(extend_method)
                extend_method = ask_multichoice('Which profile extension method to use?', allowed_methods, 'list', true);
            elseif ~ismember(extend_method, allowed_methods)
                E.badinput('EXTEND_METHOD must be one of: %s', strjoin(allowed_methods, ', '));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generative validation methods %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function insitu_comparison = generate_vcd_comparison_structs(varargin)
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('data_source', '');
            p.addParameter('extend_method', '');
            p.addParameter('v3_dir', '');
            p.addParameter('v2_dir', '');
            p.addParameter('do_save', nan);
            p.addParameter('overwrite', nan);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            data_source = pout.data_source;
            profile_extend_method = pout.extend_method;
            behr_v3_dir = pout.v3_dir;
            behr_v2_dir = pout.v2_dir;
            do_save = pout.do_save;
            do_overwrite = pout.overwrite;
            
            allowed_data_sources = {'aircraft', 'pandora'};
            if isempty(data_source)
                data_source = ask_multichoice('Which data source to use?', allowed_data_sources, 'list', true);
            elseif ~ismember(data_source, allowed_data_sources)
                E.badinput('DATA_SOURCE must be one of: %s', strjoin(allowed_data_sources, ', '));
            end
            
            if strcmpi(data_source, 'aircraft')
                profile_extend_method = misc_behr_v3_validation.get_profile_extend_method(profile_extend_method);
            end
            
            if isnan(do_save)
                do_save = ask_yn('Save the comparison files generated?');
            elseif ~isscalar(do_save) || ~islogical(do_save)
                E.badinput('DO_SAVE must be a scalar logical');
            end
            
            if do_save
                if strcmpi(data_source, 'aircraft')
                    comparison_save_name = misc_behr_v3_validation.profile_comp_file(profile_extend_method);
                elseif strcmpi(data_source, 'pandora')
                    comparison_save_name = misc_behr_v3_validation.pandora_comp_file;
                end
                
                if exist(comparison_save_name, 'file')
                    if isnan(do_overwrite)
                        do_overwrite = ask_yn(sprintf('%s exists. Overwrite?', comparison_save_name));
                    end
                    
                    if ~do_overwrite
                        insitu_comparison = [];
                        fprintf('%s exists and ''overwrite'' is false. Aborting.\n', comparison_save_name);
                        return
                    end
                end
            end
            
            versions = {'v2','v3'};
            regions = {'us'};
            prof_modes = {'monthly', 'daily'};
            insitu_comparison = struct();
            for v=1:numel(versions)
                if strcmpi(versions{v}, 'v3')
                    n_regions = numel(regions);
                    n_profs = numel(prof_modes);
                    % The proper version strings won't work as field names,
                    % so we map the "versions" version name to the proper
                    % one here.
                    behr_vers = 'v3-0B';
                else
                    n_regions = 1;
                    n_profs = 1;
                    behr_vers = 'v2-1C';
                end
                
                version_comp = struct();
                
                for r=1:n_regions
                    for p=1:n_profs
                        if strcmpi(versions{v}, 'v3')
                            if isempty(behr_v3_dir)
                                behr_dir = behr_paths.BEHRMatSubdir(regions{r}, prof_modes{p});
                            else
                                behr_dir = fullfile(behr_v3_dir, lower(prof_modes{p}));
                            end
                            behr_prefix = regexp(behr_filename(today, prof_modes{p}, regions{r}), 'OMI_BEHR.*(?=_v\d-\d[A-Z]_\d\d\d\d\d\d\d\d)', 'match', 'once');
                            wrf_prof_mode = ''; % an empty string tell verify_sat_vs_aircraft to determine the profile mode from the BEHR Data structure
                        else
                            if isempty(behr_v2_dir)
                                behr_dir = misc_behr_v3_validation.behr_v2_dir;
                            else
                                behr_dir = behr_v2_dir;
                            end
                            behr_prefix = 'OMI_BEHR';
                            wrf_prof_mode = 'monthly'; % verify_sat_vs_aircraft can read the BEHRProfileMode field in Data to determine which profiles to use, but the v2 files don't have that.
                        end
                        
                        if strcmpi(data_source, 'aircraft')
                            % Set up the necessary inputs for comparing
                            % against aircraft data
                            
                            % First we do the four DISCOVER campaigns, which are nice
                            % because they are geared towards satellite validation with
                            % spirals clearly marked in the data.
                            spiral_campaigns = {'discover_md', 'discover_ca', 'discover_tx', 'discover_co'};
                            other_campaigns = {'seac4rs', 'dc3', 'arctas_carb', 'intex_b', 'soas'};
                            % For each campaign that doesn't identify the
                            % profiles in the data, we need to specify which
                            % of the precreated range files that identify
                            % profiles by the UTC ranges to use.
                            range_files = {'SEAC4RS_Profile_Ranges_Std.mat',...
                                'DC3_Profile_Ranges_Std.mat',...
                                'ARCTAS-CA Altitude Ranges Exclusive 3.mat',...
                                'INTEXB_Profile_UTC_Ranges.mat',...
                                'SOAS_Profile_UTC_Ranges.mat'};
                            time_windows = [1.5, 3]; % how far away from satellite overpass (in hours) the profile is allowed to be
                            
                            behr_comp = struct();
                            all_campaigns = [spiral_campaigns, other_campaigns];
                            campaign_range_files = [repmat({''}, size(spiral_campaigns)), range_files];
                            % Specify whether to use the LIF or
                            % NCAR/chemiluminesnce/non-LIF NO2 measurement.
                            no2_fields = {'no2_lif','no2_lif','no2_lif','no2_lif','no2_lif','no2_lif','no2_lif','no2_lif','no2_ncar'};
                        elseif strcmpi(data_source, 'pandora')
                            % Right now I only have Pandora data for the
                            % four DISCOVER campaigns.
                            all_campaigns = {'discover_md','discover_ca','discover_tx','discover_co'};
                            time_windows = 1;
                        else
                            E.notimplemented('No input settings for data source = %s', data_source);
                        end
                        
                        for a=1:numel(all_campaigns)
                            campaign = struct();
                            for b=1:numel(time_windows)
                                [~, time_fieldname] = misc_behr_v3_validation.format_profile_time_range(time_windows(b));
                                if strcmpi(data_source, 'aircraft')
                                    all_profs_final = misc_behr_v3_validation.make_profile_vcd_comparison_for_one_campaign(all_campaigns{a}, prof_modes{p}, behr_dir, behr_prefix, behr_vers, campaign_range_files{a}, no2_fields{a}, time_windows(b), profile_extend_method, wrf_prof_mode);
                                elseif strcmpi(data_source, 'pandora')
                                    all_profs_final = misc_behr_v3_validation.make_pandora_vcd_comparison_for_one_campaign(all_campaigns{a}, prof_modes{p}, behr_dir, behr_prefix, behr_vers, time_windows(b));
                                else
                                    E.notimplemented('No validation method defined for data source = %s', data_source);
                                end
                                campaign.(time_fieldname) = all_profs_final;
                            end
                            behr_comp.(all_campaigns{a}) = campaign;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                    end
                end
                insitu_comparison.(versions{v}) = version_comp;
            end
            
            if do_save
                save(comparison_save_name, '-struct', 'insitu_comparison');
            end
        end
        
        function all_profs_final = make_pandora_vcd_comparison_for_one_campaign(campaign, prof_mode, behr_dir, behr_prefix, behr_vers, time_window)
            comparison_params = {'behr_dir', behr_dir,...
                'behr_prefix', behr_prefix,...
                'behr_version', behr_vers,...
                'behr_prof_mode', prof_mode,...
                'time_range', time_window};
            all_profs_final = run_pandora_verification(campaign, comparison_params{:});
        end
        
        function all_profs_final = make_profile_vcd_comparison_for_one_campaign(campaign, prof_mode, behr_dir, behr_prefix, behr_vers, campaign_range_file, no2_field, time_window, extension_mode, wrf_prof_mode)
            
            comparison_params = {'behr_dir', behr_dir,...
                'behr_prefix', behr_prefix,...
                'behr_version', behr_vers,...
                'utc_range_file', campaign_range_file,...
                'no2_field', no2_field,...
                'DEBUG_LEVEL',1,...
                'time_range', time_window,...
                'wrf_prof_mode', wrf_prof_mode,...
                'gc_data_dir', misc_behr_v3_validation.gc_data_path(),...
                'gc_data_year', 2012,...
                'gc_file_pattern', 'ts_12_14_satellite.%savg.nc',...
                'gc_file_date_fmt', 'yyyymm',...
                'prof_extension', extension_mode,...
                'match_bl_only', 3,...
                'min_height', 1,... % added b/c some SEAC4RS profiles are very short, causing issues with the extrapolation
                };
            try
                [all_profs_final, all_profs_detail] = run_insitu_verification(campaign, prof_mode, comparison_params{:});
                all_profs_final.details = all_profs_detail;
            catch err
                if strcmp(err.identifier, 'call_verify:file_not_found')
                    % The daily product will not have
                    % data for some days which causes
                    % an error in
                    % Run_Spiral_Verification.
                    all_profs_final = [];
                    fprintf(err.message);
                else
                    rethrow(err);
                end
            end
        end
        
        function [times, time_fieldname] = format_profile_time_range(win)
            base_time = datenum('2000-01-01 13:30:00');
            start_time = base_time - win/24;
            end_time = base_time + win/24;
            
            times = {datestr(start_time, 'HH:MM'), datestr(end_time, 'HH:MM')};
            time_fieldname = sprintf('t%s_%s', datestr(start_time, 'HHMM'), datestr(end_time, 'HHMM'));
        end
        
        function [gcas_comparison, gcas_comparison_vec] = generate_gcas_comparison_struct(do_save)
            E = JLLErrors;
            
            if ~exist('do_save', 'var')
                do_save = ask_yn('Save the comparison files generated?');
            elseif ~isscalar(do_save) || ~islogical(do_save)
                E.badinput('DO_SAVE must be a scalar logical');
            end
            
            versions = {'v2','v3'};
            regions = {'us'};
            prof_modes = {'monthly', 'daily'};
            sat_fields = {'BEHRColumnAmountNO2Trop', 'ColumnAmountNO2Trop'}; % can add VisOnly later if desired
            all_campaigns = {'discover_tx'};
            time_windows = [1.5, 3];
            gcas_comparison = struct();
            gcas_comparison_vec = struct();
            for v=1:numel(versions)
                if strcmpi(versions{v}, 'v3')
                    n_regions = numel(regions);
                    n_profs = numel(prof_modes);
                else
                    n_regions = 1;
                    n_profs = 1;
                end
                
                version_comp = struct();
                version_comp_vec = struct();
                
                for r=1:n_regions
                    for p=1:n_profs
                        if strcmpi(versions{v}, 'v3')
                            behr_dir = behr_paths.BEHRMatSubdir(regions{r}, prof_modes{p});
                        else
                            behr_dir = misc_behr_v3_validation.behr_v2_dir;
                        end
                        
                        behr_comp = struct();
                        behr_comp_vec = struct();
                        
                        for a=1:numel(all_campaigns)
                            campaign = struct();
                            campaign_vec = struct();
                            for b=1:numel(time_windows)
                                [~, time_fn] = misc_behr_v3_validation.format_profile_time_range(time_windows(b));
                                % Construct the options list
                                opts_list = {'campaign', all_campaigns{a},...
                                    'behr_dir', behr_dir,...
                                    'vectorize', false,...
                                    'cloud_prod', 'omi',...
                                    'cloud_frac_max', 0.2,...
                                    'row_anomaly', 'XTrackFlags',...
                                    'sat_fields', sat_fields,...
                                    'time_window', time_windows(b)};
                                try
                                    campaign.(time_fn) = run_gcas_verification(opts_list{:});
                                    campaign_vec.(time_fn) = vectorize_gcas_matches(campaign.(time_fn));
                                catch err
                                    if strcmp(err.identifier, 'load_behr_file_for_gcas:file_not_found')
                                        campaign.(time_fn) = [];
                                        campaign_vec.(time_fn) = [];
                                    else
                                        rethrow(err)
                                    end
                                end
                            end
                            
                            behr_comp.(all_campaigns{a}) = campaign;
                            behr_comp_vec.(all_campaigns{a}) = campaign_vec;
                        end
                        
                        fn = sprintf('%s_%s', regions{r}, prof_modes{p});
                        version_comp.(fn) = behr_comp;
                        version_comp_vec.(fn) = behr_comp_vec;
                    end
                end
                gcas_comparison.(versions{v}) = version_comp;
                gcas_comparison_vec.(versions{v}) = version_comp_vec;
            end
            
            if do_save
                save(misc_behr_v3_validation.gcas_comp_file, '-struct', 'gcas_comparison');
                save(misc_behr_v3_validation.gcas_vec_comp_file, '-struct', 'gcas_comparison_vec');
            end
        end
        
        function generate_profile_comparison_struct(varargin)
            % This will generate a structure containing WRF data matched to
            % individual DISCOVER profiles, as well as information about
            % the OMI overpass time vs. the profile time
            %
            % The hierarchy of the structure will be:
            %   profile type (daily, monthly, or version 2)
            %       campaign
            %           pXXXXXX (XXXXXX is the profile number)
            %               profdate (average profile UTC datetime)
            %               nearest OMI overpass (in space) time
            %               match
            
            p = inputParser;
            p.addParameter('do_save', nan);
            % By default, we want to restrict the campaign average profiles
            % to just the 1.5 hours before or after the standard OMI
            % overpass of 1:30 pm.
            p.addParameter('avg_prof_lst_range', [12, 15]);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            do_save = pout.do_save;
            avg_prof_lst_range = pout.avg_prof_lst_range;
            
            if isnan(do_save)
                do_save = ask_yn('Save the resulting match structures?');
            end
            
            profile_types = {'daily', 'monthly', 'v2'};
            campaigns = {'seac4rs', 'soas', 'dc3', 'discover_md', 'discover_ca', 'discover_tx', 'discover_co'};
            %campaigns = {'discover_ca'};
            
            % Loop through each campaign, load each merge file, find each
            % profile, extract the necessary fields to match up the
            % aircraft and WRF data, and pass it to match_wrf2aircraft to
            % do the actual matching
            profile_comp_struct = make_empty_struct_from_cell(profile_types);
            
            for a=1:numel(campaigns)
                % List the merge files
                [merge_names, ~, merge_dir] = merge_field_names(campaigns{a});
                merges = dirff(fullfile(merge_dir, '*.mat'));
                all_merge_wrf_dirs = make_empty_struct_from_cell(profile_types, {{}});
                for b=1:numel(merges)
                    M = load(merges(b).name);
                    
                    % Make raw structures for each profile that is listed
                    % in the file. This is only used when matching
                    % individual profiles, so if we're working on a
                    % campaign that does not include individual profiles
                    % (e.g. DC3) then can (and must) skip this part, since
                    % the Raws structure relies on the profile numbers
                    % field. We can also skip calculating the nearest OMI
                    % time, since that might be time consuming and it's
                    % only needed for matching/filtering individual
                    % profiles.
                    if ~isempty(merge_names.profile_numbers)
                        Raws = misc_behr_v3_validation.make_raw_struct(M.Merge, merge_names, campaigns{a});
                        profile_fns = fieldnames(Raws);
                        
                        
                        % Calculate the OMI overpass time from the v3 BEHR
                        % files, getting the closest overpass in space. We can
                        % use this to decide if a profile is important for the
                        % retrieval
                        OmiTimes = misc_behr_v3_validation.calc_nearest_omi_times(M.Merge.metadata.date, Raws);
                    end
                    
                    for p=1:numel(profile_types)
                        % Get the right WRF directory for the given profile
                        % type
                        if strcmpi(profile_types{p}, 'v2')
                            wrf_dir = misc_behr_v3_validation.wrf_v2_dir;
                        else
                            try
                                wrf_dir = find_wrf_path('us', profile_types{p}, M.Merge.metadata.date);
                            catch err
                                if strcmp(err.identifier, 'find_wrf_path:dir_does_not_exist')
                                    continue
                                else
                                    rethrow(err);
                                end
                            end
                        end
                        
                        if ~ismember(wrf_dir, all_merge_wrf_dirs.(profile_types{p}))
                            all_merge_wrf_dirs.(profile_types{p}){end+1} = wrf_dir;
                        end
                        
                        % Only DISCOVER campaigns have profile numbers.
                        % Other campaigns (e.g. DC3) didn't do specific
                        % satellite verification spirals, so we can't
                        % compare individual profiles (at least not as
                        % easily). So now that we've gathered the WRF
                        % directories that we need, skip the individual
                        % profile matching.
                        if isempty(merge_names.profile_numbers)
                            continue
                        end
                        
                        for f=1:numel(profile_fns)
                            Match = match_wrf2aircraft(Raws.(profile_fns{f}), wrf_dir, profile_types{p});
                            xlon = Match.wrf.xlon;
                            xlat = Match.wrf.xlat;
                            % These will be the same in every Match
                            % structure, including them greatly increases
                            % the size of the final file unnecessarily
                            Match.wrf = rmfield(Match.wrf, 'xlon');
                            Match.wrf = rmfield(Match.wrf, 'xlat');
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).match = Match;
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).prof_date = mean(Raws.(profile_fns{f}).dvec);
                            profile_comp_struct.(profile_types{p}).(campaigns{a}).(profile_fns{f}).omi_time = OmiTimes.(profile_fns{f});
                            profile_comp_struct.wrf_xlon = xlon;
                            profile_comp_struct.wrf_xlat = xlat;
                        end
                    end
                end
                
                % Match the WRF output to the entire P3 flight
                for p=1:numel(profile_types)
                    if ~isempty(all_merge_wrf_dirs.(profile_types{p}))
                        Match = match_wrf2campaigns.(campaigns{a})(profile_types{p}, all_merge_wrf_dirs.(profile_types{p}), 'lst_range', avg_prof_lst_range);
                        profile_comp_struct.(profile_types{p}).(campaigns{a}).All.match = Match;
                    end
                end
            end
            
            if do_save
                save(misc_behr_v3_validation.wrf_comp_file, '-struct', 'profile_comp_struct')
            end
        end
        
        function generate_behr_wrf_surfpres_comparison
            % This function will make monthly averages of the surface
            % pressure derived from GLOBE data and compare it to the
            % surface pressure in the monthly WRF files.
            
            test_year = 2012;
            for m = 1:12
                start_date = datenum(test_year, m, 1);
                end_date = datenum(test_year, m, eomday(test_year, m));
                % We don't need to reject any pixels because the terrain
                % pressure should be valid for all pixels
                [behr_pres, longrid, latgrid] = behr_time_average(start_date, end_date, 'avgfield', 'GLOBETerpres', 'rejectmode', 'none');
                
                wrf_info = ncinfo(find_wrf_path('us','monthly',start_date,'fullpath'));
                wrf_lon = double(ncread(wrf_info.Filename, 'XLONG'));
                wrf_lat = double(ncread(wrf_info.Filename, 'XLAT'));
                wrf_pres = double(ncread(wrf_info.Filename, 'pres'));
                wrf_pres = wrf_pres(:,:,1);
                
                behr_pres_interp = interp2(longrid, latgrid, behr_pres, wrf_lon, wrf_lat);
                
                if m == 1
                    all_wrf_pres = nan([size(wrf_pres), 12]);
                    all_behr_pres = nan([size(wrf_pres), 12]);
                    lon = wrf_lon;
                    lat = wrf_lat;
                end
                
                all_wrf_pres(:,:,m) = wrf_pres;
                all_behr_pres(:,:,m) = behr_pres_interp;
            end
            
            save(misc_behr_v3_validation.pres_comp_file, 'all_wrf_pres', 'all_behr_pres', 'lon', 'lat');
        end
        
        %%%%%%%%%%%%%%%%%%%%%
        % Utility functions %
        %%%%%%%%%%%%%%%%%%%%%
        
        function OmiTimes = calc_nearest_omi_times(this_date, Raws)
            Data = load_behr_file(this_date, 'monthly', 'us');
            % For each orbit, get the center longitude/latitude line. Also
            % get the average time
            center_lon = cell(size(Data));
            center_lat = cell(size(Data));
            avg_time = cell(size(Data));
            for a=1:numel(Data)
                n = round(size(Data(a).Longitude,2)/2);
                center_lon{a} = Data(a).Longitude(:,n);
                center_lat{a} = Data(a).Latitude(:,n);
                avg_time{a} = omi_time_conv(mean(Data(a).Time));
            end
            
            OmiTimes = make_empty_struct_from_cell(fieldnames(Raws));
            % For each profile, calculate the average lon/lat and figure
            % out which overpass comes the closest
            raw_fns = fieldnames(Raws);
            for f=1:numel(raw_fns)
                shortest_dist = Inf;
                omi_time = NaN;
                prof_avg_lon = nanmean(Raws.(raw_fns{a}).lon);
                prof_avg_lat = nanmean(Raws.(raw_fns{a}).lat);
                for a=1:numel(center_lon)
                    min_dist = min(sqrt( (prof_avg_lon - center_lon{a}).^2 + (prof_avg_lat - center_lat{a}).^2 ));
                    if min_dist < shortest_dist
                        shortest_dist = min_dist;
                        omi_time = avg_time{a};
                    end
                end
                OmiTimes.(raw_fns{f}) = omi_time;
            end
        end
        
        function Raws = make_raw_struct(Merge, merge_names, campaign)
            % The Raw structure for input into match_wrf2aircraft requires
            % the fields lon, lat, pres, dvec, and campaign. pres must be
            % pressure in hPa. dvec must be a vector that gives the date
            % and UTC time as a Matlab datenumber.
            
            
            
            % First get all the profile numbers
            profnums = remove_merge_fills(Merge, merge_names.profile_numbers);
            u_profnums = unique(profnums(profnums > 0));
            prof_fns = sprintfmulti('p%d', u_profnums);
            
            substruct = struct('lon', [], 'lat', [], 'pres', [], 'dvec', [], 'campaign', campaign);
            Raws = make_empty_struct_from_cell(prof_fns, substruct);
            
            % Read the other necessary variables. Read in UTC with NO2 to
            % ensure it has the same fill values.
            [AllRaw.no2, utc, ~, AllRaw.lon, AllRaw.lat] = remove_merge_fills(Merge, merge_names.no2_lif, 'unit', 'ppm',...
                'lon', merge_names.longitude, 'lat', merge_names.latitude);
            AllRaw.pres = remove_merge_fills(Merge, merge_names.pressure, 'unit', 'hPa');
            
            % UTC is given in seconds after midnight.
            if any(utc < 0)
                warning('Fill values in UTC time vector')
            end
            AllRaw.dvec = datenum(Merge.metadata.date) + utc ./ (24*60*60);
            
            data_fns = fieldnames(AllRaw);
            
            % Now assign the proper subsets to the individual profile
            % fields
            
            for a=1:numel(u_profnums)
                pp = profnums == u_profnums(a);
                for b=1:numel(data_fns)
                    Raws.(prof_fns{a}).(data_fns{b}) = AllRaw.(data_fns{b})(pp);
                end
            end
        end
        
        function [no2, xlon, xlat] = load_wrf_no2(prof_type, wrf_date, west_east_inds, north_south_inds, bottom_top_ind)
            if strcmpi(prof_type, 'v2')
                wrf_file_name = sprintf('m%02d_NO2_profile.mat', wrf_date);
                W = load(fullfile(misc_behr_v3_validation.wrf_v2_dir, wrf_file_name));
                no2 = permute(W.PROFILE.NO2_profile, [3 2 1]); % no2 in these files is arrange bottom_top, south_north, west_east; it needs to be the other way here
                xlon = W.PROFILE.Longitude'; % likewise lat and lon need transposed to have west_east in the first dimension
                xlat = W.PROFILE.Latitude';
            elseif strcmpi(prof_type, 'monthly')
                wrf_file_name = sprintf('WRF_BEHR_monthly_%02d.nc', wrf_date);
                [no2, xlon, xlat] = read_wrf_vars(find_wrf_path('us',prof_type,wrf_date), wrf_file_name, {'no2', 'XLONG', 'XLAT'});
            elseif strcmpi(prof_type, 'daily')
                wrf_file_name = sprintf('wrfout_d01_%s', datestr(wrf_date, 'yyyy-mm-dd_HH-MM-SS'));
                [no2, xlon, xlat] = read_wrf_vars(find_wrf_path('us',prof_type,wrf_date), wrf_file_name, {'no2', 'XLONG', 'XLAT'});
            else
                E.notimplemented('Loading WRF files for mode %s', prof_type);
            end
            
            no2 = no2(west_east_inds, north_south_inds, bottom_top_ind);
            xlon = xlon(west_east_inds, north_south_inds);
            xlat = xlat(west_east_inds, north_south_inds);
        end
        
        function locs = read_locs_file()
            locs_file = fullfile(misc_behr_v3_validation.my_dir, 'Workspaces', 'trend_locations.nc');
            ni = ncinfo(locs_file);
            ncvarnames = {ni.Variables.Name};
            struct_tmp_cell = cell(1, 2*numel(ncvarnames));
            for a=1:numel(ncvarnames)
                struct_tmp_cell{(a-1)*2 + 1} = ncvarnames{a};
                val = ncread(locs_file, ncvarnames{a});
                if ischar(val)
                    val = cellstr(val);
                else
                    val = num2cell(val);
                end
                
                struct_tmp_cell{a*2} = val;
            end
            
            % This makes it into a structure where each index is a location
            locs = struct(struct_tmp_cell{:});
        end
        
        function [xx, yy] = find_loc_indices(loc, lon, lat, radius)
            % LOC must be a scalar element of the locations structure, LON
            % and LAT must be 2D arrays of longitude and latitude
            % coordinates for an NO2 average or similar 2D field. RADIUS
            % must be a scalar number of grid cells in each direction to
            % get. If omitted, defaults to 0.
            E = JLLErrors;
            
            if ~exist('radius', 'var')
                radius = 0;
            end
            
            sz = size(lon);
            
            r = sqrt((lon - loc.Longitude).^2 + (lat - loc.Latitude).^2);
            [~, i_min] = min(r(:));
            [xx, yy] = ind2sub(size(lon), i_min);
            
            xx = (xx - radius):(xx + radius);
            xx = xx(xx > 0 & xx <= sz(1));
            yy = (yy - radius):(yy + radius);
            yy = yy(yy > 0 & yy <= sz(2));
        end
        
        function [behr_daily, wrf_monthly, wrf_daily] = load_wrf_and_behr_data(plot_loc, date_in, varargin)
            p = inputParser;
            p.addParameter('load_gridded', false);
            p.addParameter('max_frac_nans', 1);
            p.parse(varargin{:});
            pout = p.Results;
            
            load_gridded = pout.load_gridded;
            max_frac_nans = pout.max_frac_nans;
            
            if load_gridded
                [Native, Data] = load_behr_file(date_in, 'daily', 'us');
            else
                Data = load_behr_file(date_in, 'daily', 'us');
                Native = Data; % solely for the "is location in swath" check
            end
            
            wrf_int_mode = 'box';
            plot_radius = 10; %TODO: make this resolution agnostic. Will probably need to modify find_loc_indices
            
            wrf_daily.lon = {};
            wrf_daily.lat = {};
            wrf_daily.no2_vcds = {};
            
            behr_daily.lon = {};
            behr_daily.lat = {};
            behr_daily.no2_vcds = {};
            behr_daily.no2_scds = {};
            if load_gridded
                behr_daily.areaweights = {};
            end
            
            for a=1:numel(Data)
                % Loop over every swath. If the location is in the box
                % defined by the four corner pixel centers of the swath
                % plot it
                corner_x = [Native(a).Longitude(1,1), Native(a).Longitude(1,end), Native(a).Longitude(end,end), Native(a).Longitude(end,1)];
                corner_y = [Native(a).Latitude(1,1), Native(a).Latitude(1,end), Native(a).Latitude(end,end), Native(a).Latitude(end,1)];
                if ~inpolygon(plot_loc.Longitude, plot_loc.Latitude, corner_x, corner_y);
                    continue
                end
                
                badpix = mod(Data(a).BEHRQualityFlags, 2) ~= 0;
                Data(a).BEHRColumnAmountNO2Trop(badpix) = NaN;
                [xx_sat, yy_sat] = misc_behr_v3_validation.find_loc_indices(plot_loc, Data(a).Longitude, Data(a).Latitude, plot_radius);
                
                these_scds = Data(a).BEHRColumnAmountNO2Trop(xx_sat,yy_sat) .* Data(a).BEHRAMFTrop(xx_sat, yy_sat);
                if fracnan(these_scds) > max_frac_nans
                    continue
                end
                
                behr_daily.no2_vcds{end+1} = Data(a).BEHRColumnAmountNO2Trop(xx_sat,yy_sat);
                behr_daily.no2_scds{end+1} = these_scds;
                
                if load_gridded
                    behr_daily.lon{end+1} = Data(a).Longitude(xx_sat,yy_sat);
                    behr_daily.lat{end+1} = Data(a).Latitude(xx_sat,yy_sat);
                else
                    % If not loaded gridded data, using the pixel corners
                    % if usually better with pcolor() because of how it
                    % matches colors and coordinates.
                    behr_daily.lon{end+1} = squeeze(Data(a).FoV75CornerLongitude(1,xx_sat,yy_sat));
                    behr_daily.lat{end+1} = squeeze(Data(a).FoV75CornerLatitude(1,xx_sat,yy_sat));
                end
                
                if load_gridded
                    behr_daily.areaweights{end+1} = Data(a).Areaweight;
                end
                
                % Get the WRF file name, but just retrieve the date b/c
                % files produced on the cluster will have different paths
                % and use the subset files, which aren't stored locally.
                [~, daily_wrf_file_tmp] = fileparts(Data(a).BEHRWRFFile);
                wrf_date = date_from_wrf_filenames(daily_wrf_file_tmp);
                daily_file = fullfile(find_wrf_path('us','daily',wrf_date), sprintf('wrfout_d01_%s', datestr(wrf_date, 'yyyy-mm-dd_HH-MM-SS')));
                daily_wrf_no2_vcds = compute_wrf_trop_columns(daily_file, wrf_int_mode, 200);
                daily_wrf_lon = ncread(daily_file, 'XLONG');
                daily_wrf_lat = ncread(daily_file, 'XLAT');
                
                [xx_daily, yy_daily] = misc_behr_v3_validation.find_loc_indices(plot_loc, daily_wrf_lon, daily_wrf_lat, plot_radius);
                
                wrf_daily.lon{end+1} = daily_wrf_lon(xx_daily, yy_daily);
                wrf_daily.lat{end+1} = daily_wrf_lat(xx_daily, yy_daily);
                wrf_daily.no2_vcds{end+1} = daily_wrf_no2_vcds(xx_daily, yy_daily);
            end
            
            % Load a new monthly file, if needed
            if numel(behr_daily.lon) > 0
                wrf_monthly_file = fullfile(find_wrf_path('us','monthly',date_in), sprintf('WRF_BEHR_monthly_%02d.nc', month(date_in)));
                monthly_wrf_no2_vcds_tmp = compute_wrf_trop_columns(wrf_monthly_file, wrf_int_mode, 200);
                monthly_wrf_lon_tmp = ncread(wrf_monthly_file, 'XLONG');
                monthly_wrf_lat_tmp = ncread(wrf_monthly_file, 'XLAT');
                
                [xx_monthly, yy_monthly] = misc_behr_v3_validation.find_loc_indices(plot_loc, monthly_wrf_lon_tmp, monthly_wrf_lat_tmp, plot_radius);
                
                wrf_monthly.lon = monthly_wrf_lon_tmp(xx_monthly, yy_monthly);
                wrf_monthly.lat = monthly_wrf_lat_tmp(xx_monthly, yy_monthly);
                wrf_monthly.no2_vcds = monthly_wrf_no2_vcds_tmp(xx_monthly, yy_monthly);
            else
                wrf_monthly.lon = [];
                wrf_monthly.lat = [];
                wrf_monthly.no2_vcds = [];
            end
        end
        
        function [box_vals, box_pres, box_quartiles] = make_boxplot_bins(pres, no2)
            [box_vals, box_pres, box_quartiles] = bin_omisp_pressure(pres, no2, 'median');
            return
            
            %             % old way, would need for true box plots
            %             [binned_no2, binned_pres] = bin_omisp_pressure(pres, no2, 'binonly');
            %             % Because BOXPLOT is pretty simpleminded, we need to make the
            %             % values and pressures vectors equal lengths
            %             box_vals = veccat(binned_no2{:});
            %             box_pres = nan(size(box_vals));
            %             i = 1;
            %             for a=1:numel(binned_no2)
            %                 j = i + numel(binned_no2{a}) - 1;
            %                 box_pres(i:j) = binned_pres(a);
            %                 i = j + 1;
            %             end
            %
            %             box_nans = isnan(box_vals);
            %             box_vals(box_nans) = [];
            %             box_pres(box_nans) = [];
        end
        
        function [aircraft_comp, pandora_comp] = match_aircraft_and_pandora_sites(aircraft_comp, pandora_comp, varargin)
            p = inputParser;
            p.addParameter('match_time', false);
            p.parse(varargin{:});
            pout = p.Results;
            match_time = pout.match_time;
            
            aircraft_comp.details = match_verify_struct_details(aircraft_comp, aircraft_comp.details);
            
            aircraft_coords = [aircraft_comp.profile_lon(:), aircraft_comp.profile_lat(:)];
            pandora_coords = [pandora_comp.profile_lon(:), pandora_comp.profile_lat(:)];
            
            xx_air = false(size(aircraft_coords,1),1);
            xx_pandora = false(size(pandora_coords,1),1);
            distance_crit = 0.05;
            for i_pan = 1:size(pandora_coords, 1)
                distances = sqrt((pandora_coords(i_pan, 1) - aircraft_coords(:,1)).^2 + (pandora_coords(i_pan, 2) - aircraft_coords(:,2)).^2);
                xx_air = xx_air | distances < distance_crit;
                xx_pandora(i_pan) = any(distances < distance_crit);
            end
            
            if match_time
                air_dates = cellfun(@(x) datenum(x, 'yyyy-mm-dd'), aircraft_comp.profile_dates);
                pandora_dates = floor(pandora_comp.omi_time);
                xx_air = xx_air & ismember(air_dates, pandora_dates);
                xx_pandora = xx_pandora & ismember(pandora_dates, air_dates);
            end
            
            fns_air = fieldnames(aircraft_comp);
            for i_air = 1:numel(fns_air)
                if ~ischar(aircraft_comp.(fns_air{i_air}))
                    aircraft_comp.(fns_air{i_air}) = aircraft_comp.(fns_air{i_air})(xx_air);
                end
            end
            
            fns_pandora = fieldnames(pandora_comp);
            for i_pan = 1:numel(fns_pandora)
                if ~ischar(pandora_comp.(fns_pandora{i_pan}))
                    pandora_comp.(fns_pandora{i_pan}) = pandora_comp.(fns_pandora{i_pan})(xx_pandora);
                end
            end
        end
        
        function save_scd_results(filename, results, eval_opts)
            
            h5_filename = strrep(filename, 'mat', 'h5');
            loc_names = {results.loc_name};
            valid_results = ~cellfun(@isempty, loc_names);
            results = results(valid_results);
            save(filename, 'results');
            
            for i_res = 1:numel(results)
                if isempty(results(i_res).loc_name)
                    continue
                end
                
                sanitized_locname = regexprep(results(i_res).loc_name, '\W', '');
                this_date = results(i_res).date;
                misc_behr_v3_validation.scd_subgroup(h5_filename, sprintf('/%s/%s', sanitized_locname, this_date), results(i_res));
            end
            
            description = ['This file contains the locations randomly selected to compare OMI SCDs to WRF VCDs\n',...
                'to determine how well WRF is capturing the wind fields. The hierarchy is SiteName/Date/<datasets>.\n',...
                'Each day will have a user value and user confidence field. Confidence ranges from 1 (low) to 3 (high).\n',...
                'User value will have one of the following numeric values, with its corresponsing meaning:\n\n%s\n\n',...
                'The BEHR SCDs and WRF monthly and daily VCDs will be included as subgroups with lat/lon coordinates.'];
            opts_description = sprintfmulti('    %d: %s', num2cell(1:numel(eval_opts)), eval_opts);
            h5writeatt(h5_filename, '/', 'Description', sprintf(description, strjoin(opts_description, '\n')));
        end
        
        function [data_structs, opts] = load_comparison_data(varargin)
            p = inputParser;
            p.addParameter('comp_file','');
            p.addParameter('data_source', '');
            p.addParameter('extend_method', '');
            p.addParameter('x_var', '');
            p.addParameter('y_var', '');
            p.addParameter('prof_mode', '');
            p.addParameter('campaigns', {});
            p.addParameter('time_range', '');
            p.addParameter('version', '');
            
            % This allows us to just pass all parameters from any calling
            % function
            p.KeepUnmatched = true;
            
            p.parse(varargin{:});
            pout = p.Results;
            
            comp_file = pout.comp_file;
            data_source = pout.data_source;
            prof_extend_method = pout.extend_method;
            x_var = pout.x_var;
            y_var = pout.y_var;
            prof_mode = pout.prof_mode;
            campaigns = pout.campaigns;
            time_range = pout.time_range;
            product_version = pout.version;
            
            if ~isempty(comp_file)
                comp_struct = load(comp_file);
            else
                data_source = opt_ask_multichoice('Which data source to use?', {'aircraft', 'pandora'}, data_source, '"data_source"', 'list', true);
                if strcmpi(data_source, 'pandora')
                    comp_struct = load(misc_behr_v3_validation.pandora_comp_file);
                else
                    opts.extend_method = misc_behr_v3_validation.get_profile_extend_method(prof_extend_method);
                    comp_struct = load(misc_behr_v3_validation.profile_comp_file(opts.extend_method));
                end
            end
            
            [opts.data_source, allowed_vars, opts.labels] = misc_behr_v3_validation.comp_struct_type(comp_struct);
            
            opts.x_var = opt_ask_multichoice('Which variable to plot on the x-axis?', allowed_vars, x_var, '"x_var"', 'list', true);
            opts.y_var = opt_ask_multichoice('Which variable to plot on the y-axis?', allowed_vars, y_var, '"y_var"', 'list', true);
            
            opts.labels.x = opts.labels.(opts.x_var);
            opts.labels.y = opts.labels.(opts.y_var);
            
            allowed_prof_modes = {'monthly', 'daily', 'both', 'md-swap'};
            if regcmpi(opts.x_var, 'behr') || regcmpi(opts.y_var, 'behr')
                opts.prof_mode = opt_ask_multichoice('Which profile mode to use for v3 data?', allowed_prof_modes, prof_mode, '"prof_mode"', 'list', true);
            else
                % There's no difference in the SP data between our monthly
                % and daily profile product
                opts.prof_mode = 'monthly';
            end
            
            allowed_campaigns = fieldnames(comp_struct.v2.us_monthly);
            if isempty(campaigns)
                opts.campaigns = ask_multiselect('Which campaign(s) to plot?', allowed_campaigns);
            else
                if ischar(campaigns)
                    campaigns = {campaigns};
                elseif ~iscellstr(campaigns)
                    E.badinput('CAMPAIGN must be a string or cell array of strings')
                end
                if any(~ismember(campaigns, allowed_campaigns))
                    E.badinput('All strings in CAMPAIGN must be one of: %s', strjoin(allowed_campaigns, ', '));
                end
                opts.campaigns = campaigns;
            end
            
            allowed_times = fieldnames(comp_struct.v2.us_monthly.(allowed_campaigns{1}));
            opts.time_range = opt_ask_multichoice('Which time range to use?', allowed_times, time_range, '"time_range"', 'list', true);
            
            % Form the data into a common structure so that the calling
            % function doesn't necessarily need to know any of the options
            if regcmpi(opts.x_var, 'behr') || regcmpi(opts.y_var, 'behr')
                v2_string = 'v2.1C';
                v3M_string = 'v3.0B (M)';
                v3D_string = 'v3.0B (D)';
                opts.title_product_str = 'BEHR';
            else
                v2_string = 'v2.1';
                v3M_string = 'v3.0';
                v3D_string = 'v3.0';
                opts.title_product_str = 'NASA';
            end
            
            opts.version = opt_ask_multichoice('Which version of the product?', {'v2', 'v3', 'both'}, product_version, '"version"', 'list', true);
            
            swap_md_air = false;
            switch lower(opts.prof_mode)
                case 'monthly'
                    tmp_data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_monthly};
                    opts.legend_strings = {v2_string, v3M_string};
                    xx_v3 = [false, true];
                case 'daily'
                    tmp_data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_daily};
                    opts.legend_strings = {v2_string, v3D_string};
                    xx_v3 = [false, true];
                case 'md-swap'
                    tmp_data_structs = {comp_struct.v3.us_monthly, comp_struct.v3.us_daily};
                    opts.legend_strings = {v3M_string, v3D_string};
                    xx_v3 = [true, true];
                    swap_md_air = true;
                    if ~strcmpi(opts.version, 'v3')
                        warning('For "md-swap" data source, version must be set to "v3", so setting version to "v3" now');
                    end
                    opts.version = 'v3';
                case 'both'
                    tmp_data_structs = {comp_struct.v2.us_monthly, comp_struct.v3.us_monthly, comp_struct.v3.us_daily};
                    opts.legend_strings = {v2_string, v3M_string, v3D_string};
                    xx_v3 = [false, true, true];
                otherwise
                    E.notimplemented('prof_mode == %s', opts.prof_mode);
            end
            
            switch lower(opts.version)
                case 'v2'
                    tmp_data_structs = tmp_data_structs(~xx_v3);
                    opts.legend_strings = opts.legend_strings(~xx_v3);
                case 'v3'
                    tmp_data_structs = tmp_data_structs(xx_v3);
                    opts.legend_strings = opts.legend_strings(xx_v3);
                case 'both'
                    % Don't need to do anything
                otherwise
                    E.notimplemented('version == %s', opts.version);
            end
            
            data_structs = cell(size(tmp_data_structs));
            for i_dat = 1:numel(data_structs)
                this_comp = tmp_data_structs{i_dat};
                for i_cam = 1:numel(opts.campaigns)
                    data_substruct = this_comp.(opts.campaigns{i_cam}).(opts.time_range);
                    if ~isempty(data_substruct)
                        % Special case: if we're swapping out the monthly
                        % and daily profiles' aircraft VCDs (since the
                        % integration of those VCDs uses limits derived
                        % from the profiles), then we need to override
                        % where the x-variable comes from. We assume that
                        % tmp_data_structs{1} is monthly and
                        % tmp_data_structs{2} is daily profiles.
                        if swap_md_air
                            if i_dat == 1
                                i_dat_swap = 2;
                            elseif i_dat == 2
                                i_dat_swap = 1;
                            else
                                E.notimplemented('Cannot swap monthly/daily aircraft VCDs if more than two data structs');
                            end
                            x_substruct = tmp_data_structs{i_dat_swap}.(opts.campaigns{i_cam}).(opts.time_range);
                            if ~isempty(x_substruct)
                                data_substruct.x = x_substruct.(opts.x_var);
                            else
                                warning('No daily prof data for %s, so not swapping x variable', opts.campaigns{i_cam});
                                data_substruct.x = data_substruct.(opts.x_var);
                            end
                        else
                            data_substruct.x = data_substruct.(opts.x_var);
                        end
                            
                        data_substruct.y = data_substruct.(opts.y_var);
                    end
                    data_substruct.campaign = opts.campaigns{i_cam};
                    data_structs{i_dat}(i_cam) = data_substruct;
                end
            end
            
            dum=1;
        end
        
        function [data_struct, air_opts, pan_opts] = load_combined_aircraft_pandora(varargin)
            p = advInputParser;
            p.addParameter('match', false);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pout = p.Results;
            
            match_aircraft_pandora = pout.match;
            
            [data_struct, air_opts] = misc_behr_v3_validation.load_comparison_data(varargin{:}, 'data_source', 'aircraft');
            
            x_var = air_opts.x_var;
            if regcmp(x_var, '^air')
                x_var = 'pandora_no2';
            end
            y_var = air_opts.y_var;
            if regcmp(y_var, '^air')
                y_var = 'pandora_no2';
            end
            
            [pan_data, pan_opts] = misc_behr_v3_validation.load_comparison_data(air_opts, 'data_source', 'pandora', 'x_var', x_var, 'y_var', y_var, 'time_range', 't1230_1430');
            
            for i_dat = 1:numel(data_struct)
                if ~isfield(data_struct{i_dat}, 'x')
                    continue
                end
                
                if match_aircraft_pandora
                    [data_struct{i_dat}, pan_data{i_dat}] = misc_behr_v3_validation.match_aircraft_and_pandora_sites(data_struct{i_dat}, pan_data{i_dat}, 'match_time', true);
                end
                
                
                data_struct{i_dat}.x = veccat(data_struct{i_dat}.x, pan_data{i_dat}.x);
                data_struct{i_dat}.y = veccat(data_struct{i_dat}.y, pan_data{i_dat}.y);
            end
        end
        
        function plot_sites(ax, varargin)
            p = inputParser;
            p.addParameter('site_type', 'all')
            
            p.parse(varargin{:});
            pout = p.Results;
            
            site_type = pout.site_type;
            
            locs = misc_behr_v3_validation.read_locs_file;
            if ~strcmpi(site_type, 'all')
                xx = strcmpi({locs.SiteType}, site_type);
                locs = locs(xx);
            end
            
            x = [locs.Longitude];
            y = [locs.Latitude];
            names = {locs.ShortName};
            
            % Add a few extra cities
            x = veccat(x, -86.8025);
            y = veccat(y, 33.5207);
            names = veccat(names, {'Birmingham'});
            
            % Control which ones should be aligned right instead to avoid
            % overlapping names
            sites_align_right = {'Birmingham'};
            
            x_lims = get(ax, 'XLim');
            y_lims = get(ax, 'YLim');
            
            % Text will get printed outside the axes if we don't cut it
            % down
            in_lims = x > x_lims(1) & x < x_lims(2) - 0.5 & y > y_lims(1) + 0.25 & y < y_lims(2) - 0.25;
            x = x(in_lims);
            y = y(in_lims);
            names = names(in_lims);
            
            % Anything too close to the right side of the map will need
            % it's text put on the left instead of the right
            %align_right = x > x_lims(2) - diff(x_lims) * 0.2;
            align_right = ismember(names, sites_align_right);
            line(x, y, 'marker', 'p', 'color', 'k', 'linestyle', 'none');
            text(x(~align_right)+0.25, y(~align_right), names(~align_right), 'BackgroundColor', 'w');
            text(x(align_right)-0.25, y(align_right), names(align_right),  'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
        end
        
        function slopes = gather_slopes(fit_info, campaigns_list)
            % Gather the slopes into a structure for each campaign in
            % discover_campaigns. Prefer ones with V > 0 only, if available.
            discover_campaigns = misc_behr_v3_validation.discover_campaigns;
            for i_cam = 1:numel(discover_campaigns)
                this_campaign = upper(strrep(discover_campaigns{i_cam}, '_', '-'));
                xx = strcmp(campaigns_list(:,1), sprintf('%s (V > 0)', this_campaign));
                if sum(xx) == 0
                    xx = strcmp(campaigns_list(:,1), this_campaign);
                end
                
                slopes.(discover_campaigns{i_cam}).sp_no2 = fit_info.sp_no2{xx};
                slopes.(discover_campaigns{i_cam}).behr_no2 = fit_info.behr_no2{xx};
            end
        end
        
        function [unwt_avg, unwt_std, wt_avg, wt_std, wt_eff_sample_size] = average_aircraft_pandora_slopes(air_data, pandora_data)
            air = [air_data.P]';
            % Concatenating P should alternate slope, intercept, slope,
            % intercepts, etc. We just want the slopes.
            air = air(1:2:end);
            air_n= [air_data.num_pts]';
            air_std = [air_data.StdDevM]';
            pandora = [pandora_data.P]';
            pandora = pandora(1:2:end);
            pandora_n = [pandora_data.num_pts]';
            pandora_std = [pandora_data.StdDevM]';
            
            unwt_avg = (air + pandora)/2;
            unwt_std = sqrt( (0.5 .* air_std).^2 + (0.5 .* pandora_std).^2 );
            
            sum_n = air_n + pandora_n;
            air_wt = air_n / sum_n;
            pandora_wt = pandora_n / sum_n;
            wt_avg = air .* air_wt + pandora .* pandora_wt;
            wt_std = sqrt( ( air_wt .* air_std ).^2 + ( pandora_wt .* pandora_std ).^2 );
            both_wts = [air_wt, pandora_wt];
            
            % I'm taking the definition of effective sample size from http://www.analyticalgroup.com/download/Alternative%20Approaches.pdf
            wt_eff_sample_size = nansum(both_wts).^2/nansum(both_wts.^2);
        end
        
        function xx_not_outlier = is_not_outlier(x,y)
            xx_not_outlier = ~isoutlier(x) & ~isoutlier(y);
        end
        
        function [fit_data, fit_line, varargout] = calculate_fit(x_no2, y_no2, varargin)
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('remove_outliers',nan);
            p.addParameter('remove_neg_sat', nan);
            p.addParameter('force_origin', nan);
            
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pout = p.Results;
            
            do_remove_outliers = pout.remove_outliers;
            do_remove_neg_sat = pout.remove_neg_sat;
            if pout.force_origin
                regression_type = 'orth-origin';
            else
                regression_type = 'rma';
            end
            
            % Make these parameters to be consistent with the other
            % methods, so that parameters can just be all passed to this
            % function.
            if any(structfun(@isnan, pout))
                E.badinput('All parameters (%s) must be given', strjoin(fieldnames(pout), ', '));
            end
            
            if strcmpi(x_no2, 'load')
                [data, opts] = misc_behr_v3_validation.load_comparison_data(y_no2{:});
                varargout{1} = opts;
            elseif strcmpi(x_no2, 'load-combined')
                [data, air_opts, pan_opts] = misc_behr_v3_validation.load_combined_aircraft_pandora(y_no2{:});
                varargout = {air_opts, pan_opts};
            else
                data{1}.x = x_no2;
                data{1}.y = y_no2;
            end
            
            for i_dat = 1:numel(data)
                if ~isfield(data{i_dat}, 'x')
                    continue
                end
                
                x = data{i_dat}.x;
                y = data{i_dat}.y;
                
                xx_keep = true(size(x));
                
                xx_keep = xx_keep & (~isnan(x) & ~isnan(y));
                
                % Temporary until I fix the imaginary column issue
                not_imag = difftol(imag(x), zeros(size(x))) & difftol(imag(y), zeros(size(y)));
                xx_keep = xx_keep & not_imag;
                
                if do_remove_neg_sat
                    xx_keep = xx_keep & x >= 0 & y >= 0;
                end
                
                
                if do_remove_outliers
                    % Must remove the points we're not considering for
                    % other reasons before removing outliers. Trying to
                    % keep xx_keep the length of the original data vectors
                    % so that it can be returned to cut down other vectors
                    % if necessary.
                    xx_keep(xx_keep) = xx_keep(xx_keep) & misc_behr_v3_validation.is_not_outlier(x(xx_keep), y(xx_keep));
                end
                [fit_line(i_dat).x, fit_line(i_dat).y, fit_line(i_dat).legend, fit_data(i_dat)] = calc_fit_line(x(xx_keep), y(xx_keep), 'regression', regression_type, 'xcoord', [-1e17, 1e17], 'significance', true);
                
            end
            
            function [x, y] = check_data(data)
                if isscalar(data)
                    x = data.x;
                    y = data.y;
                else
                    E.notimplemented('getting multiple data structures simultaneously');
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions %
        %%%%%%%%%%%%%%%%%%%%%%
        
        function plot_one_gcas_comparison(varargin)
            allowed_vars = {'GCAS_NO2vcd', 'ColumnAmountNO2Trop', 'BEHRColumnAmountNO2Trop'};
            
            labels = struct('GCAS_NO2vcd', 'GCAS NO_2 VCD (molec. cm^{-2})',...
                'ColumnAmountNO2Trop', 'NASA NO_2 VCD (molec. cm^{-2})',...
                'BEHRColumnAmountNO2Trop', 'BEHR NO_2 VCD (molec. cm^{-2})');
            
            misc_behr_v3_validation.plot_one_vcd_comparison(misc_behr_v3_validation.gcas_vec_comp_file, allowed_vars, labels, varargin{:});
        end
        
        function [values, column_names, row_names, section_end_rows, fit_data ] = tabulate_insitu_comparisons(varargin)
            % Create a table of slopes, intercepts, and R2 values for sat
            % vs. aircraft or pandora comparisons.
            %
            % The follow parameters allow you to bypass the interactive
            % questions:
            %
            %   'data_source' - either 'aircraft' or 'pandora'; indicates
            %   which data set to compare sat columns against
            %
            %   'campaigns' - cell array of strings indicating which
            %   campaigns to include.
            %
            %   'remove_neg_vcd' - whether or not to remove negative
            %   satellite VCDs. Can either be a scalar logical (in which
            %   case it will apply to all campaigns) or a logical array the
            %   same size as 'campaigns', in which case you can specify
            %   which campaigns to remove negative VCDs for.
            %
            %   'extend_method' - 'wrf', 'geos', or 'extrap'; only useful
            %   if 'data_source' is 'aircraft' as it controls how the
            %   aircraft profiles are extended to surface and tropopause.
            E = JLLErrors;
            p = inputParser;
            p.addParameter('data_source','')
            p.addParameter('campaigns',{});
            p.addParameter('remove_neg_vcd', []);
            p.addParameter('extend_method', '');
            p.addParameter('force_origin', nan);
            p.addParameter('swap_md', nan)
            p.addParameter('return_as_table', nan);
            p.addParameter('remove_outliers', nan);
            p.addParameter('fit_quantities', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            data_source = pout.data_source;
            allowed_data_sources = {'aircraft', 'pandora'};
            if isempty(data_source)
                data_source = ask_multichoice('Which data source to use?', allowed_data_sources, 'list', true);
            elseif ~ismember(data_source, allowed_data_sources)
                E.badinput('DATA_SOURCE must be one of: %s', strjoin(allowed_data_sources, ', '));
            end
            
            
            campaigns = pout.campaigns;
            extend_method = opt_ask_multichoice('Which profile extension method to use?', misc_behr_v3_validation.profile_extend_methods, pout.extend_method, 'extend_method', 'list', true);
            do_remove_outliers = opt_ask_yn('Remove outliers?', pout.remove_outliers, '"remove_outliers"');
            force_origin = opt_ask_yn('Force fits through origin?', pout.force_origin, '"force_origin"');
            do_swap_monthly_daily = opt_ask_yn('Swap monthly/daily prof aircraft VCDs?', pout.swap_md, '"swap_md"');
            return_as_table = opt_ask_yn('Return as table (not array)?', pout.return_as_table, '"return_as_table"');
            fit_quantities = opt_ask_multichoice('Include slopes and std dev. or slope, int, R2?', {'std. dev.', 'int+R2'}, pout.fit_quantities, '"fit_quantities"', 'list', true);
            do_intercept_r2 = strcmpi(fit_quantities, 'int+R2');
            
            remove_neg_vcd = opt_ask_yn('Remove negative VCDs?', pout.remove_neg_vcd, '"remove_neg_vcd"',...
                'test_fxn', @(x) islogical(x) && (isscalar(x) || numel(x) == numel(campaigns)), 'ask_condition', @isempty,...
                'test_msg', '%s must be a logical array, either scalar or the same size as "campaigns"');
            
            switch lower(data_source)
                case 'aircraft'
                    allowed_vars = {'air_no2_nasa', 'air_no2_behr';...
                        'sp_no2', 'behr_no2'};
                    comp_file = misc_behr_v3_validation.profile_comp_file(extend_method);
                    time_range = 't1200_1500';
                case 'pandora'
                    allowed_vars = {'pandora_no2', 'pandora_no2';...
                        'sp_no2', 'behr_no2'};
                    comp_file = misc_behr_v3_validation.pandora_comp_file();
                    time_range = 't1230_1430';
                otherwise
                    E.badinput('data_source == "%s" not supported, may be "aircraft" or "pandora"', data_source);
            end
            
            
            
            if isempty(campaigns)
                % Need this to get the available campaigns
                comp_struct = load(comp_file);
                campaigns = fieldnames(comp_struct.v2.us_monthly);
            end
            
            if isscalar(remove_neg_vcd)
                remove_neg_vcd = repmat(remove_neg_vcd, size(campaigns));
            end
            
            if do_swap_monthly_daily
                prof_mode = 'md-swap';
                version = 'v3';
            else
                prof_mode = 'both';
                version = 'both';
            end
            
            common_opts = {'data_source', data_source, 'extend_method', extend_method, 'prof_mode', prof_mode, 'version', version, 'time_range', time_range, 'remove_outliers', do_remove_outliers, 'plot_type', 'scatter',...
                'color_by', 'none', 'match_pandora_aircraft', false, 'force_origin', force_origin};
            
            for i_var = 1:2
                x_var = allowed_vars{1, i_var};
                y_var = allowed_vars{2, i_var};
                for i_campaign = 1:numel(campaigns)
                    [fig, fit_substruct] = misc_behr_v3_validation.plot_one_vcd_comparison('campaigns', campaigns{i_campaign}, 'x_var', x_var, 'y_var', y_var, 'remove_neg_sat', remove_neg_vcd(i_campaign), common_opts{:});
                    close(fig)
                    fit_data.(y_var){i_campaign} = fit_substruct;
                end
            end
            
            % Ultimately want the table to look like
            %
            %                              | Slope | Intercept | R2 or p-value |
            %            | SP   | v2       |       |           |               |
            %   campaign |      | v3       |       |           |               |
            %            | BEHR | v2       |       |           |               |
            %            |      | v3 (M)   |       |           |               |
            %            |      | v3 (D)   |       |           |               |
            row_names = {};
            values = [];
            section_end_rows = [];
            products = {'sp_no2','behr_no2'};
            product_table_names = {'SP','BEHR'};
            for i_campaign = 1:numel(campaigns)
                if remove_neg_vcd(i_campaign)
                    vcd_str = ' ($V > 0$)';
                else
                    vcd_str = '';
                end
                for i_prod = 1:numel(products)
                    substruct = fit_data.(products{i_prod}){i_campaign};
                    for i_sub = 1:numel(substruct)
                        this_row_names = {sprintf('%s%s', upper(strrep(campaigns{i_campaign},'_', '-')), vcd_str), sprintf('%s %s', product_table_names{i_prod}, substruct(i_sub).prof_type)};
                        if do_intercept_r2
                            this_row_values = [substruct(i_sub).P(1), substruct(i_sub).P(2), substruct(i_sub).R2];
                        else
                            this_row_values = [substruct(i_sub).P(1), substruct(i_sub).StdDevM];
                        end
                        
                        row_names = cat(1, row_names, this_row_names);
                        values = cat(1, values, this_row_values);
                    end
                end
                section_end_rows = veccat(section_end_rows, size(values,1));
            end
            
            if do_intercept_r2
                column_names = {'Campaign', 'Product', 'Slope','Intercept','$R^2$'};
            else
                column_names = {'Campaign', 'Product', 'Slope'};
            end
            
            if return_as_table
                column_names{end} = 'R2';
                table_row_names = cell(size(values,1),1);
                for i_row = 1:numel(table_row_names)
                    table_row_names{i_row} = strjoin(row_names(i_row, :), ': ');
                end
                values = array2table(values, 'RowNames', table_row_names, 'VariableNames', column_names(3:end));
            end
        end
        
        function [values, colnames, rownames, samples, extra_hlines] = make_aircraft_pandora_combo(varargin)
            % This one regresses aircraft and pandora together, rather than
            % averaging their slopes after the fact.
            p = advInputParser;
            p.addParameter('return_as_table', true);
            p.addParameter('products', {});
            p.addParameter('fit_quantities', '');
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pout = p.Results;
            
            % These have an obvious default if running interactively
            do_return_as_table = pout.return_as_table;
            products_to_keep = pout.products;
            
            % These do not, so ask
            fit_quantities = opt_ask_multichoice('Include slopes and std dev. or slope, int, R2?', {'std. dev.', 'int+R2'}, pout.fit_quantities, '"fit_quantities"', 'list', true);
            do_intercept_r2 = strcmpi(fit_quantities, 'int+R2');
            
            campaigns = misc_behr_v3_validation.discover_campaigns;
            variables = {'air_no2_nasa', 'sp_no2';...
                'air_no2_behr', 'behr_no2'};
            values = [];
            rownames = {};
            samples = [];
            extra_hlines = logical([]);
            for i_cam = 1:numel(campaigns)
                this_campaign = campaigns{i_cam};
                for i_var = 1:size(variables)
                    air_var = variables{i_var, 1};
                    sat_var = variables{i_var, 2};
                    data_opts = {'data_source','aircraft','extend_method','geos','x_var',air_var,'y_var',sat_var,...
                        'prof_mode','both','campaigns',{this_campaign},'time_range','t1200_1500','version','both'};
                    [fit_data_match, ~, air_opts] = misc_behr_v3_validation.calculate_fit('load-combined', [data_opts, 'match', true], 'remove_outliers', true, 'remove_neg_sat', true, 'force_origin', false);
                    fit_data_all = misc_behr_v3_validation.calculate_fit('load-combined', [data_opts, 'match', false], 'remove_outliers', true, 'remove_neg_sat', true, 'force_origin', false);
                    
                    if numel(fit_data_match) ~= numel(fit_data_all)
                        E.notimplemented('The matched and all fits returned different numbers of structures');
                    end
                    
                    for i_fit = 1:numel(fit_data_match)
                        this_rowname = {upper(strrep(this_campaign, '_', '-')), sprintf('%s %s', air_opts.title_product_str, air_opts.legend_strings{i_fit})};
                        % Easier way to homogenize this with the average table
                        this_rowname = strrep(this_rowname, 'NASA', 'SP');
                        if ~isempty(products_to_keep) && ~ismember(this_rowname{2}, products_to_keep)
                            continue
                        end
                        if do_intercept_r2
                            these_vals = [fit_data_match(i_fit).P(1), fit_data_match(i_fit).P(2), fit_data_match(i_fit).R2, fit_data_all(i_fit).P(1), fit_data_all(i_fit).P(2), fit_data_all(i_fit).R2];
                        else
                            these_vals = [fit_data_match(i_fit).P(1), fit_data_match(i_fit).StdDevM, fit_data_all(i_fit).P(1), fit_data_all(i_fit).StdDevM];
                        end
                        values = cat(1, values, these_vals);
                        samples = cat(1, samples, [fit_data_match(i_fit).num_pts, fit_data_all(i_fit).num_pts]);
                        rownames = cat(1, rownames, this_rowname);
                        extra_hlines = cat(1, extra_hlines, false);
                    end
                end
                % Mark the last line of each campaign to have a horizontal
                % line after it.
                extra_hlines(end) = true;
            end
            
            
            
            if do_return_as_table
                if do_intercept_r2
                    colnames = {'MatchedData','MatchedDataSigma','AllData','AllDataSigma'};
                else
                    colnames = {'MatchedSlope', 'MatchedInt', 'MatchedR2', 'AllSlope', 'AllInt', 'AllR2'};
                end
                
                for i_row = 1:size(rownames,1)
                    table_rownames{i_row} = strjoin(rownames(i_row, :), ': ');
                end
                values = array2table(values, 'RowNames', table_rownames, 'VariableNames', colnames);
            else
                if do_intercept_r2
                    colnames = {'', '', 'Matched Data', 'Matched Data', 'Matched Data', 'All Data', 'All Data', 'All Data';...
                        'Campaign', 'Product', 'Slope', 'Intercept', 'R2', 'Slope', 'Intercept', 'R2'};
                else
                    colnames = {'Campaign', 'Product', 'Slope (Matched)', 'Slope (All)'};
                end
                
            end
        end
        
        function varargout = make_aircraft_pandora_avg(varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('extend_method', '');
            p.addParameter('plot_or_table', '');
            p.addParameter('products', {});
            p.addParameter('remove_outliers', nan);
            p.addParameter('force_origin', nan);
            p.addParameter('swap_md', nan)
            p.addParameter('return_as_table', nan);
            p.parse(varargin{:});
            pout = p.Results;
            
            requested_products = pout.products;
            
            product_helper = struct('name', {'SP v2.1', 'SP v3.0', 'BEHR v2.1', 'BEHR v3.0 (M)', 'BEHR v3.0 (D)'},...
                'get_fxn', {@(S, campaign) S.(campaign).sp_no2(1),...
                @(S, campaign) S.(campaign).sp_no2(2),...
                @(S, campaign) S.(campaign).behr_no2(1),...
                @(S, campaign) S.(campaign).behr_no2(2),...
                @(S, campaign) S.(campaign).behr_no2(3)});
            product_formats = struct('marker', {'o','*','d','v','^'}, 'color', {'r',[1 0.5 0],'b',[0.5 0 0.5],[0 0.5 0]},...
                'markersize',11,'linewidth',2);
            products = {product_helper.name};
            
            extend_method = opt_ask_multichoice('Which profile extend method to use for the aircraft data', misc_behr_v3_validation.profile_extend_methods, pout.extend_method, '"extend_method"', 'list', true);
            plot_or_table = opt_ask_multichoice('Make a plot of slopes or return information for a table?', {'plot', 'table'}, pout.plot_or_table, '"plot_or_table"', 'list', true);
            is_plot = strcmpi(plot_or_table, 'plot');
            is_table = strcmpi(plot_or_table, 'table');
            
            if isempty(requested_products)
                product_inds = ask_multiselect('Which products to include?', products, 'returnindex', true);
            else
                if any(~ismember(requested_products, products))
                    E.badinput('One of the products given is not valid. Allowed values are: %s', strjoin(products, ', '));
                end
                product_inds = find(ismember(products, requested_products));
            end
            
            do_remove_outliers = opt_ask_yn('Remove outliers?', pout.remove_outliers, '"remove_outliers"');
            force_origin = opt_ask_yn('Force fits through origin?', pout.force_origin, '"force_origin"');
            do_swap_monthly_daily = opt_ask_yn('Swap monthly/daily prof aircraft VCDs?', pout.swap_md, '"swap_md"');
            
            product_helper = product_helper(product_inds);
            product_formats = product_formats(product_inds);
            
            
            
            if is_table
                return_as_table = opt_ask_yn('Return as table (not array)?', pout.return_as_table, '"return_as_table"');
            end
            
            discover_campaigns = misc_behr_v3_validation.discover_campaigns;
            
            common_opts = {'campaigns', discover_campaigns, 'remove_neg_vcd', [false, false, false, true], 'extend_method', extend_method, 'force_origin', force_origin,...
                'swap_md', do_swap_monthly_daily, 'return_as_table', false, 'remove_outliers', do_remove_outliers};
            [~, ~, rownames_aircraft, ~, fit_data_aircraft] = misc_behr_v3_validation.tabulate_insitu_comparisons(common_opts{:}, 'data_source', 'aircraft');
            [~, ~, rownames_pandora, ~, fit_data_pandora] = misc_behr_v3_validation.tabulate_insitu_comparisons(common_opts{:}, 'data_source', 'pandora');
            
            air_slopes = misc_behr_v3_validation.gather_slopes(fit_data_aircraft, unique(rownames_aircraft(:, 1), 'stable'));
            pandora_slopes = misc_behr_v3_validation.gather_slopes(fit_data_pandora, unique(rownames_pandora(:, 1), 'stable'));
            
            
            
            if is_table
                table_values = [];
                sample_sizes = [];
                rownames = {};
            elseif is_plot
                plot_empty_array = nan(numel(discover_campaigns), numel(product_helper));
                plot_unwt_values = plot_empty_array;
                plot_unwt_err = plot_empty_array;
                plot_wt_values = plot_empty_array;
                plot_wt_err = plot_empty_array;
            else
                E.notimplemented('initialization for neither plot nor table')
            end
            
            for i_cam = 1:numel(discover_campaigns)
                for i_prod = 1:numel(product_helper)
                    try
                        air_behr = product_helper(i_prod).get_fxn(air_slopes, discover_campaigns{i_cam});
                        pandora_behr = product_helper(i_prod).get_fxn(pandora_slopes, discover_campaigns{i_cam});
                    catch err
                        if strcmpi(err.identifier, 'MATLAB:badsubscript')
                            % This happens if we try to get a BEHR daily
                            % result for a campaign that doesn't have one,
                            % so just move on to the next campaign.
                            continue
                        else
                            rethrow(err)
                        end
                    end
                    
                    [unwt_avg, unwt_std, wt_avg, wt_std, wt_eff_sample_size] = misc_behr_v3_validation.average_aircraft_pandora_slopes(air_behr, pandora_behr);
                    
                    
                    if is_table
                        campaign_name = strrep(upper(discover_campaigns{i_cam}), '_', '-');
                        these_rownames = {campaign_name, product_helper(i_prod).name};
                        rownames = cat(1, rownames, these_rownames);
                        
                        table_values = cat(1, table_values, [unwt_avg, unwt_std, wt_avg, wt_std]);
                        sample_sizes = cat(1, sample_sizes, [2, wt_eff_sample_size]);
                    elseif is_plot
                        plot_unwt_values(i_cam, i_prod) = unwt_avg;
                        plot_unwt_err(i_cam, i_prod) = unwt_std;
                        plot_wt_values(i_cam, i_prod) = wt_avg;
                        plot_wt_err(i_cam, i_prod) = wt_std;
                    else
                        E.notimplemented('data concatenation for neither plot nor table')
                    end
                end
            end
            
            if strcmpi(plot_or_table, 'table')
                if return_as_table
                    colnames = {'AvgSlope', 'AvgSlopeSigma', 'WtAvgSlope', 'WtAvgSlopeSigma'};
                    t_rownames = cell(size(rownames,1),1);
                    for i_row = 1:size(rownames,1)
                        t_rownames{i_row} = strjoin(rownames(i_row,:), ': ');
                    end
                    varargout{1} = array2table(table_values, 'VariableNames', colnames, 'RowNames', t_rownames);
                    varargout{2} = sample_sizes;
                else
                    colnames = {'Campaign','Product','Avg. slope','Weighted avg. slope'};
                    varargout = {table_values, colnames, rownames, sample_sizes};
                end
                
            elseif strcmpi(plot_or_table, 'plot')
                xtick_labels = cellfun(@(x) strrep(upper(x), '_', '-'), discover_campaigns, 'uniformoutput', false);
                varargout{1} = figure;
                l = scatter_grouped(plot_unwt_values, plot_unwt_err, 'group_fmt', product_formats, 'error_bar_fmt', rmfield(product_formats, {'marker','markersize'}));
                
                line([0, 5], [1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 2);
                legend(l, {product_helper.name});
                set(gca, 'xticklabel', xtick_labels,'xticklabelrotation', 30, 'fontsize', 16);
                title('Unweighted average');
                
                varargout{2} = figure;
                l2 = scatter_grouped(plot_wt_values, plot_wt_err, 'group_fmt', product_formats, 'error_bar_fmt', rmfield(product_formats, {'marker','markersize'}));
                line([0, 5], [1 1], 'color', 'k', 'linestyle', '--', 'linewidth', 2);
                legend(l2, {product_helper.name});
                set(gca, 'xticklabel', xtick_labels,'xticklabelrotation', 30, 'fontsize', 16);
                title('Weighted average');
            else
                E.notimplemented('No action defined for plot_or_table == "%s"', plot_or_table);
            end
            
        end
        
        function varargout = plot_one_vcd_comparison(varargin)
            % This should be called from another method in this class that
            % provides the right comp_file, allowed_vars, and labels.
            %
            % Plots one v2-v3 comparison against in situ aircraft data
            % derived columns. prof_mode may be 'monthly', 'daily' or
            % 'both'. campaign and time_range must match the fields in the
            % misc_behr_v3_validation.profile_comp_file. If omitted these
            % inputs will be asked interactively
            E = JLLErrors;
            p = inputParser;
            p.addParameter('remove_outliers',[]);
            p.addParameter('plot_type','');
            p.addParameter('color_by','');
            p.addParameter('match_pandora_aircraft', []);
            p.addParameter('remove_neg_sat', nan);
            p.addParameter('map_quantity', '');
            p.addParameter('force_origin', nan);
            
            p.KeepUnmatched = true;
            
            p.parse(varargin{:});
            pout = p.Results;
            
            do_remove_outliers = pout.remove_outliers;
            plot_type = pout.plot_type;
            color_by = pout.color_by;
            match_pandora_aircraft = pout.match_pandora_aircraft;
            remove_neg_sat = pout.remove_neg_sat;
            map_quantity = pout.map_quantity;
            
            [data_structs, opts] = misc_behr_v3_validation.load_comparison_data(varargin{:});
            x_var = opts.x_var;
            y_var = opts.y_var;
            labels = opts.labels;
            campaigns = opts.campaigns;
            legend_strings = opts.legend_strings;
            
            allowed_plot_types = {'scatter','map','zonal','pres top'};
            if isempty(plot_type)
                plot_type = ask_multichoice('Which type of plot?', allowed_plot_types, 'list', true);
            elseif ~ismember(plot_type, allowed_plot_types)
                E.badinput('plot_type must be one of: %s', strjoin(allowed_plot_types, ', '));
            end
            
            allowed_color_bys = {'prof_max', 'prof_num', 'date', 'none'};
            allowed_map_quantities = {'abs_diff', 'percent_diff', 'yx_ratio'};
            if strcmpi(plot_type, 'scatter')
                if isempty(color_by)
                    color_by = ask_multichoice('Color the scatter plot by what?', allowed_color_bys, 'list', true);
                elseif ~ismember(color_by, allowed_color_bys)
                    E.badinput('color_by must be one of: %s', strjoin(allowed_color_bys), ', ');
                end
            elseif strcmpi(plot_type, 'map')
                map_quantity = opt_ask_multichoice('Which quantity to plot on the map?', allowed_map_quantities, map_quantity, '"map_quantity"', 'list', true);
            end
            
            if isempty(do_remove_outliers)
                do_remove_outliers = ask_yn('Remove outliers?');
            elseif ~isscalar(do_remove_outliers) || ~islogical(do_remove_outliers)
                E.badinput('do_remove_outliers must be a scalar logical')
            end
            
            if isempty(match_pandora_aircraft)
                match_pandora_aircraft = ask_yn('Only use profiles/sites with both aircraft and Pandora data?');
            elseif ~isscalar(match_pandora_aircraft) || ~islogical(match_pandora_aircraft)
                E.badinput('match_pandora_aircraft must be a scalar logical')
            end
            
            remove_neg_sat = opt_ask_yn('Remove negative satellite VCDs?', remove_neg_sat, '"remove_neg_sat"');
            force_through_origin = opt_ask_yn('Force fit through origin?', pout.force_origin, '"force_origin"');
            if force_through_origin
                regression_type = 'orth-origin';
            else
                regression_type = 'RMA';
            end
            
            data_lines = gobjects(size(data_structs));
            fit_lines = gobjects(numel(data_structs),1);
            fit_legends = cell(size(fit_lines));
            line_fmts = struct('color', {'k','b','r', [0, 0.5, 0]}, 'marker', {'s','x','^','o'});
            color_labels = struct('prof_max', 'Max [NO_2] (mixing ratio)', 'prof_num', 'Profile number', 'date', '', 'none', '');
            fit_data = struct('P',[],'R2',[],'StdDevM',[],'StdDevB',[],'p_value',[], 'num_pts', [], 'is_significant', [], 'x_var', labels.x, 'y_var', labels.y, 'prof_type',legend_strings);
            keep_fit_data = false(size(fit_data));
            
            is_aircraft = strcmpi(opts.data_source, 'aircraft');
            if match_pandora_aircraft
                % Is what we were given pandora or aircraft data? If the
                % time range is 12:30 to 14:30, then it's pandora,
                % otherwise it's aircraft.
                if is_aircraft
                    PTmp = load(misc_behr_v3_validation.pandora_comp_file);
                    pandora = PTmp.v2.us_monthly;
                else
                    ATmp = load(misc_behr_v3_validation.profile_comp_file('geos'));
                    aircraft = ATmp.v2.us_monthly;
                end
            end
            
            if strcmpi(plot_type,'scatter') || strcmpi(plot_type,'zonal') || strcmpi(plot_type, 'pres top')
                fig = figure;
            else
                fig = gobjects(numel(data_structs),1);
            end
            limits = [Inf, -Inf];
            for a=1:numel(data_structs)
                lon = [];
                lat = [];
                meas_pres_top = [];
                x_no2 = [];
                y_no2 = [];
                color_val = [];
                % TODO: handle if request daily profiles for a campaign
                % that doesn't have them
                for b=1:numel(data_structs{a})
                    this_struct = data_structs{a}(b);
                    
                    if numel(fieldnames(this_struct)) > 1  % data structures with no data will just have the field "campaign"
                        if match_pandora_aircraft
                            if is_aircraft
                                this_struct = misc_behr_v3_validation.match_aircraft_and_pandora_sites(this_struct, pandora.(this_struct.campaign).t1230_1430, 'match_time', true);
                            else
                                [~, this_struct] = misc_behr_v3_validation.match_aircraft_and_pandora_sites(aircraft.(this_struct.campaign).t1200_1500, this_struct, 'match_time', true);
                            end
                        else
                            % The matching process already also matches the
                            % details (in order to be able to cut them down
                            % with the other fields)
                            try
                                this_details = match_verify_struct_details(this_struct, this_struct.details);
                            catch err
                                % This probably happens because we passed a
                                % pandora struct that doesn't need to have its
                                % detail field matched up. As long as the
                                % details field has the right number of
                                % entries, and "pandora_no2" is a field, just
                                % keep the details field as is.
                                if strcmpi(err.identifier, 'MATLAB:nonExistentField') && isfield(this_struct, 'pandora_no2') && numel(this_struct.details) == numel(this_struct.pandora_no2)
                                    this_details = this_struct.details;
                                else
                                    rethrow(err)
                                end
                            end
                        end
                        
                        if ~strcmpi(plot_type, 'scatter')
                            lon = veccat(lon, this_struct.profile_lon);
                            lat = veccat(lat, this_struct.profile_lat);
                        else
                            switch lower(color_by)
                                case 'prof_max'
                                    color_val = veccat(color_val, max(cat(1,this_details.binned_profile), [], 2));
                                case 'prof_num'
                                    profile_numbers = cellfun(@str2double, this_struct.profile_ids);
                                    color_val = veccat(color_val, profile_numbers);
                                case 'date'
                                    if is_aircraft
                                        color_val = cat(1, color_val, datenum(this_struct.profile_dates));
                                    else
                                        color_val = cat(1, color_val, this_struct.omi_time);
                                    end
                                case 'none'
                                    % do nothing
                                otherwise
                                    E.notimplemented('No method to color by %s', color_by);
                            end
                        end
                        if strcmpi(plot_type, 'pres top')
                            for i_prof = 1:numel(this_details)
                                i_top = find(this_details(i_prof).is_appended_or_interpolated == 0, 1, 'last');
                                this_pres_top = this_details(i_prof).binned_pressure(i_top);
                                meas_pres_top = veccat(meas_pres_top, this_pres_top);
                            end
                        end
                        
                        x_no2 = veccat(x_no2, this_struct.x);
                        y_no2 = veccat(y_no2, this_struct.y);
                        
                    end
                end
                
                not_nans = ~isnan(x_no2) & ~isnan(y_no2);
                if ~strcmpi(plot_type, 'scatter')
                    lon = lon(not_nans);
                    lat = lat(not_nans);
                end
                x_no2 = x_no2(not_nans);
                y_no2 = y_no2(not_nans);
                
                % Temporary until I fix the imaginary column issue
                not_imag = difftol(imag(x_no2), zeros(size(x_no2)));
                x_no2 = x_no2(not_imag);
                y_no2 = y_no2(not_imag);
                
                xx_keep = true(size(x_no2));
                
                
                if remove_neg_sat
                    if any(strcmpi(x_var, {'behr_no2', 'sp_no2'}))
                        xx_keep = xx_keep & x_no2 >= 0;
                    end
                    if any(strcmpi(y_var, {'behr_no2', 'sp_no2'}))
                        xx_keep = xx_keep & y_no2 >= 0;
                    end
                end
                
                
                if do_remove_outliers
                    xx_keep = xx_keep & misc_behr_v3_validation.is_not_outlier(x_no2, y_no2);
                end
                
                if ~isempty(x_no2) && ~isempty(y_no2)
                    limits(1) = min([min(x_no2), min(y_no2)]);
                    limits(2) = max([max(x_no2), max(y_no2)]);
                    
                    if strcmpi(plot_type, 'scatter')
                        if strcmpi(color_by, 'none')
                            data_lines(a) = line(x_no2(xx_keep), y_no2(xx_keep), 'linestyle', 'none', 'color', line_fmts(a).color, 'marker', line_fmts(a).marker);
                        else
                            data_lines(a) = scatter(x_no2(xx_keep), y_no2(xx_keep), [], color_val(xx_keep), line_fmts(a).marker);
                            if a == 1
                                cb = colorbar;
                                cb.Label.String = color_labels.(color_by);
                                if strcmpi(color_by, 'date')
                                    cbdatetick(cb);
                                    cb.TicksMode = 'manual';
                                    cb.FontSize = 8;
                                end
                            end
                            hold on
                        end
                        [fit_x, fit_y, fit_legends{a}, fit_data_tmp] = calc_fit_line(x_no2(xx_keep), y_no2(xx_keep), 'regression', regression_type, 'xcoord', [-1e17, 1e17], 'significance', true);
                        fit_lines(a) = line(fit_x, fit_y, 'color', line_fmts(a).color, 'linestyle', '--');
                        
                        fit_data(a) = copy_structure_fields(fit_data_tmp, fit_data(a), fieldnames(fit_data_tmp));
                        keep_fit_data(a) = true;
                    elseif strcmpi(plot_type, 'zonal')
                        data_lines(a) = line(lon(xx_keep), y_no2(xx_keep) - x_no2(xx_keep), 'linestyle', 'none', 'marker', line_fmts(a).marker, 'color', line_fmts(a).color);
                    elseif strcmpi(plot_type, 'map')
                        fig(a) = figure;
                        
                        used_log = false;
                        switch lower(map_quantity)
                            case 'abs_diff'
                                diff_quantity = y_no2(xx_keep) - x_no2(xx_keep);
                                clabel_str = '\Delta VCD (molec. cm^{-2})';
                            case 'percent_diff'
                                diff_quantity = reldiff(y_no2(xx_keep), x_no2(xx_keep))*100;
                                clabel_str = '%\Delta VCD (molec. cm^{-2})';
                            case 'yx_ratio'
                                diff_quantity = y_no2(xx_keep) ./ x_no2(xx_keep);
                                if all(diff_quantity > 0)
                                    diff_quantity = log10(diff_quantity);
                                    used_log = true;
                                else
                                    warning('Negative Y/X ratios, cannot use log scale');
                                end
                                clabel_str = 'Y/X';
                            otherwise
                                E.notimplemented('No method for map_quantity = %s', map_quantity);
                        end
                        
                        scatter(lon(xx_keep), lat(xx_keep), [], diff_quantity);
                        %set(gca,'xlimmode','manual');
                        title(sprintf('%s -\n %s: %s', labels.(y_var), labels.(x_var), legend_strings{a}));
                        cb = colorbar;
                        cb.Label.String = clabel_str;
                        state_outlines('k','not','ak','hi');
                        set(gca,'fontsize',16);
                        caxis(calc_plot_limits(diff_quantity, 'diff'));
                        if ~strcmpi(diff_quantity, 'yx_ratio') || used_log
                            colormap(blue_red_only_cmap)
                        end
                        if used_log
                            cb.TickLabels = sprintfmulti('%.2f', 10.^(cb.Ticks));
                        end
                        
                        misc_behr_v3_validation.plot_sites(gca, 'site_type', 'Cities');
                    elseif strcmpi(plot_type, 'pres top')
                        data_lines(a) = line(meas_pres_top(xx_keep), y_no2(xx_keep) - x_no2(xx_keep), 'linestyle', 'none', 'marker', line_fmts(a).marker, 'color', line_fmts(a).color);
                    else
                        E.notimplemented('Do not know how to make plot type %s', plot_type);
                    end
                end
            end
            
            if strcmpi(plot_type, 'scatter')
                % Scatter puts all the data structs on the same plot, so we
                % need to handle this after all the series have been
                % plotted.
                xylims(calc_plot_limits(limits, 'pow10'));
                
                xlabel(labels.x);
                ylabel(labels.y);
                
                all_lines = gobjects(numel(data_lines)*2, 1);
                all_lines(1:2:end) = data_lines;
                all_lines(2:2:end) = fit_lines;
                
                all_legends = cell(1, numel(data_structs)*2);
                all_legends(1:2:end) = legend_strings;
                all_legends(2:2:end) = fit_legends;
                
                xx_valid = ishandle(all_lines);
                
                lgnd = legend(all_lines(xx_valid), all_legends(xx_valid),'Location','best');
                title_campaigns = strrep(campaigns, '_', '\_');
                title(strjoin(upper(title_campaigns), ', '));
                ax = gca;
                ax.FontSize = 16;
                lgnd.FontSize = 10;
            elseif strcmpi(plot_type,'zonal')
                x_edges = get(gca,'xlim');
                line(x_edges, [0 0], 'linestyle', '--', 'linewidth', 2, 'color', [0.5 0.5 0.5]);
                legend(data_lines, legend_strings);
                xlabel('Longitude');
                ylabel(sprintf('%s - %s', labels.(y_var), labels.(x_var)));
                set(gca,'fontsize',16);
            elseif strcmpi(plot_type, 'pres top')
                x_edges = get(gca,'xlim');
                line(x_edges, [0 0], 'linestyle', '--', 'linewidth', 2, 'color', [0.5 0.5 0.5]);
                legend(data_lines, legend_strings);
                xlabel('Pressure at top measurement');
                ylabel(sprintf('%s -\n %s', labels.(y_var), labels.(x_var)));
                set(gca,'fontsize',16);
            end
            
            if nargout > 0
                varargout = {fig, fit_data(keep_fit_data)};
            end
        end
        
        function [pandora_far_coords, pandora_close_coords, aircraft_coords, cluster_coords] = plot_aircraft_vs_pandora_locations(varargin)
            %
            p = inputParser;
            p.addParameter('campaigns', {});
            p.addParameter('time_range', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            campaigns = pout.campaigns;
            time_range = pout.time_range;
            
            pandoras = load(misc_behr_v3_validation.pandora_comp_file());
            aircraft = load(misc_behr_v3_validation.profile_comp_file('geos'));
            
            allowed_campaigns = fieldnames(aircraft.v2.us_monthly);
            if isempty(campaigns)
                campaigns = ask_multiselect('Which campaign(s) to plot?', allowed_campaigns);
            else
                if ischar(campaigns)
                    campaigns = {campaigns};
                elseif ~iscellstr(campaigns)
                    E.badinput('CAMPAIGN must be a string or cell array of strings')
                end
                if any(~ismember(campaigns, allowed_campaigns))
                    E.badinput('All strings in CAMPAIGN must be one of: %s', strjoin(allowed_campaigns, ', '));
                end
            end
            
            allowed_times = fieldnames(aircraft.v2.us_monthly.(allowed_campaigns{1}));
            if isempty(time_range)
                time_range = ask_multichoice('Which time range to use?', allowed_times, 'list', true);
            elseif ~ischar(time_range) || ~ismember(time_range, allowed_times)
                E.badinput('TIME_RANGE must be one of the strings: %s', strjoin(allowed_times, ', '));
            end
            
            aircraft_coords = [];
            cluster_coords = [];
            pandora_close_coords = [];
            pandora_far_coords = [];
            
            
            for i_campaign = 1:numel(campaigns)
                pandora_comp = pandoras.v3.us_monthly.(campaigns{i_campaign}).t1230_1430;
                aircraft_comp = aircraft.v3.us_monthly.(campaigns{i_campaign}).(time_range);
                
                [~, pandora_close_comp] = misc_behr_v3_validation.match_aircraft_and_pandora_sites(aircraft_comp, pandora_comp);
                xx_far = ~ismember(pandora_comp.profile_lon, pandora_close_comp.profile_lon) & ~ismember(pandora_comp.profile_lat, pandora_close_comp.profile_lat);
                
                this_aircraft_coords = [aircraft_comp.profile_lon(:), aircraft_comp.profile_lat(:)];
                aircraft_coords = cat(1, aircraft_coords, this_aircraft_coords);
                platlon_far_temp = unique([pandora_comp.profile_lon(xx_far), pandora_comp.profile_lat(xx_far)],'rows');
                pandora_far_coords = cat(1, pandora_far_coords, platlon_far_temp);
                platlon_close_temp = unique([pandora_close_comp.profile_lon(:), pandora_close_comp.profile_lat(:)],'rows');
                pandora_close_coords = cat(1, pandora_close_coords, platlon_close_temp);
                
                % If we're doing a discover campaign, try to cluster the
                % aircraft profiles using kmeans so that we get one point
                % per profile.
                if regcmpi(campaigns{i_campaign}, '^discover')
                    profnums = cellfun(@(x) floor(str2double(x)/1000), aircraft_comp.profile_ids);
                    n_profs = numel(unique(profnums));
                    [~, this_cluster_means] = kmeans(this_aircraft_coords, n_profs);
                    cluster_coords = cat(1, cluster_coords, this_cluster_means);
                end
            end
            leg_strs = {'Pandora too far', 'Pandora nearby', 'Aircraft'};
            figure;
            l(1) = line(pandora_far_coords(:,1), pandora_far_coords(:,2), 'marker', 'x', 'linestyle', 'none', 'color', 'r');
            l(2) = line(pandora_close_coords(:,1), pandora_close_coords(:,2), 'marker', 'x', 'linestyle', 'none', 'color', [0 0.5 0]);
            l(3) = line(aircraft_coords(:,1), aircraft_coords(:,2), 'marker', 's', 'linestyle', 'none', 'color', 'b');
            if ~isempty(cluster_coords)
                l(4) = line(cluster_coords(:,1), cluster_coords(:,2), 'marker', '^', 'linestyle', 'none', 'color', 'c');
                leg_strs{end+1} = 'Aircraft cluster';
            end
            state_outlines('k');
            legend(l(:), leg_strs);
        end
        
        function plot_measured_no2(varargin)
            E = JLLErrors;
            p = inputParser;
            p.addParameter('campaign', '');
            p.addParameter('prof_type', '');
            p.addParameter('bottom_pres', []);
            p.addParameter('top_pres', []);
            p.addParameter('plot_quantity', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            campaign = opt_ask_multichoice('Which campaign?', allowed_campaigns, pout.campaign , '"campaign"', 'list', true);
            
            avail_prof_types = fieldnames(wrfcomp);
            xx_profs = false(size(avail_prof_types));
            for i_prof = 1:numel(avail_prof_types)
                xx_profs(i_prof) = isfield(wrfcomp.(avail_prof_types{i_prof}), campaign);
            end
            allowed_prof_types = avail_prof_types(xx_profs);
            prof_type = opt_ask_multichoice('Which profile type?', allowed_prof_types, pout.prof_type , '"campaign"', 'list', true);
            
            bottom_pres = opt_ask_number('Give the bottom pressure to include, in hPa', pout.bottom_pres, '"bottom_pres"', 'testfxn', @(x) isscalar(x) && isnumeric(x) && x > 0 , 'testmsg', 'Bottom pressure must be a positive scalar');
            top_pres = opt_ask_number('Give the top pressure to include, in hPa', pout.top_pres, '"top_pres"', 'testfxn', @(x) isscalar(x) && isnumeric(x) && x >= 0 && x < bottom_pres,...
                'testmsg', sprintf('Top pressure must be a positive scalar less than the bottom pressure (%.1f)', bottom_pres));
            plot_quantity = opt_ask_multichoice('Which quantity to plot?', {'Aircraft NO2', 'WRF NO2', 'Aircraft - WRF NO2'}, pout.plot_quantity, '"plot_quantity"', 'list', true);
            
            comp_struct = wrfcomp.(prof_type).(campaign).All.match;
            
            air_pres = comp_struct.data.pres;
            xx = air_pres <= bottom_pres & air_pres >= top_pres;
            
            air_no2 = comp_struct.data.no2(xx);
            wrf_no2 = comp_struct.wrf.no2(xx);
            lon = comp_struct.data.lon(xx);
            lat = comp_struct.data.lat(xx);
            
            switch lower(plot_quantity)
                case 'aircraft no2'
                    plot_vals = air_no2;
                    clabel_str = 'Aircraft NO2 between %.1f and %.1f hPa';
                case 'wrf no2'
                    plot_vals = wrf_no2;
                    clabel_str = 'WRF NO2 between %.1f and %.1f hPa';
                case 'aircraft - wrf no2'
                    plot_vals = air_no2 - wrf_no2;
                    clabel_str = 'Aircraft - WRF NO2 between %.1f and %.1f hPa';
                otherwise
                    E.notimplemented('No case for plot_quantity = %s', plot_quantity);
            end
            
            figure;
            scatter(lon, lat, [], plot_vals);
            cb = colorbar;
            cb.Label.String = sprintf(clabel_str, bottom_pres, top_pres);
            state_outlines('k');
        end
        
        function varargout = plot_one_wrf_comparison(prof_types, campaign, prof_number, uncert_type, varargin)
            
            p = inputParser;
            p.addParameter('bin_mode','');
            p.parse(varargin{:});
            pout = p.Results;
            
            bin_mode = pout.bin_mode;
            
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            allowed_prof_nums = fieldnames(wrfcomp.monthly.(campaign));
            if ~exist('prof_number', 'var')
                prof_number = ask_multichoice('Which profile?', allowed_prof_nums, 'list', true);
            else
                if isnumeric(prof_number) && isscalar(prof_number)
                    prof_number = sprintf('p%d', prof_number);
                elseif ~ischar(prof_number)
                    E.badinput('PROF_NUMBER must be a string or a scalar number')
                end
                
                if ~ismember(prof_number, allowed_prof_nums)
                    E.badinput('PROF_NUMBER must be one of these strings or the numeric component of it: %s', strjoin(allowed_prof_nums));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_types', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif any(~ismember(prof_types, allowed_prof_types))
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            allowed_bin_modes = {'mean', 'median'};
            if isempty(bin_mode)
                bin_mode = ask_multichoice('How to bin the profiles?', allowed_bin_modes, 'default', 'mean', 'list', true);
            elseif ~ismember(bin_mode, allowed_bin_modes)
                E.badinput('BIN_MODE must be one of: %s', strjoin(allowed_bin_modes, ', '));
            end
            
            % Do the actual plotting. First get each profile type's raw
            % data, and bin it to give the average shape. Convert from ppm
            % to ppt (hence the 1e6)
            for a=1:numel(prof_types)
                if a == 1
                    % The aircraft data should be the same for all profile
                    % types
                    no2.aircraft = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.data.no2*1e6;
                    pres.aircraft = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.data.pres;
                    [binned_no2.aircraft, binned_pres.aircraft, binned_no2_std.aircraft] = bin_omisp_pressure(pres.aircraft, no2.aircraft, bin_mode);
                end
                
                no2.(prof_types{a}) = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.wrf.no2*1e6;
                pres.(prof_types{a}) = wrfcomp.(prof_types{a}).(campaign).(prof_number).match.wrf.pres;
                [binned_no2.(prof_types{a}), binned_pres.(prof_types{a}), binned_no2_std.(prof_types{a})] = bin_omisp_pressure(pres.(prof_types{a}), no2.(prof_types{a}), bin_mode);
            end
            
            % We'll always use the same colors for the various profile
            % types. Use lighter versions of the colors for the raw data
            plt_cols = misc_behr_v3_validation.plot_colors;
            
            fig=figure;
            fns = fieldnames(no2);
            l = gobjects(numel(fns),1);
            legstr = cell(1, numel(fns));
            for a=1:numel(fns)
                if any(strcmpi(uncert_type, {'raw', 'both'}))
                    line(no2.(fns{a}), pres.(fns{a}), 'linestyle', 'none', 'marker', '.', 'color', plt_cols.(fns{a}).raw);
                end
                if any(strcmpi(uncert_type, {'std', 'both'}))
                    if strcmpi(bin_mode, 'mean')
                        scatter_errorbars(binned_no2.(fns{a}), binned_pres.(fns{a}), binned_no2_std.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'direction', 'x');
                    else
                        scatter_errorbars(binned_no2.(fns{a}), binned_pres.(fns{a}), binned_no2_std.(fns{a})(1,:), binned_no2_std.(fns{a})(2,:), 'color', plt_cols.(fns{a}).avg, 'direction', 'x');
                    end
                end
                
                l(a) = line(binned_no2.(fns{a}), binned_pres.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'linewidth', 2);
                legstr{a} = capitalize_words(fns{a});
            end
            
            legend(l, legstr);
            set(gca,'ydir','reverse','fontsize',14);
            xlabel('[NO_2] (pptv)');
            ylabel('Pressure (hPa)');
            title(prof_number);
            
            if nargout > 0
                varargout{1} = fig;
            end
        end
        
        function plot_aircraft_on_wrf(prof_type, campaign, prof_number, wrf_level)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            allowed_prof_nums = fieldnames(wrfcomp.monthly.(campaign));
            if ~exist('prof_number', 'var')
                prof_number = ask_multichoice('Which profile?', allowed_prof_nums, 'list', true);
            else
                if isnumeric(prof_number) && isscalar(prof_number)
                    prof_number = sprintf('p%d', prof_number);
                elseif ~ischar(prof_number)
                    E.badinput('PROF_NUMBER must be a string or a scalar number')
                end
                
                if ~ismember(prof_number, allowed_prof_nums)
                    E.badinput('PROF_NUMBER must be one of these strings or the numeric component of it: %s', strjoin(allowed_prof_nums));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_type', 'var')
                prof_type = ask_multichoice('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_type, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            if ~exist('wrf_level', 'var')
                wrf_level = ask_number('Which WRF level to plot (0-29)? 0 will plot the lowest level that matched the aircraft.', 'testfxn', @(x) isscalar(x) && x >= 0 && x <= 29, 'testmsg', 'Must be between 0 and 29');
            end
            
            % Okay, this one's a little different. We need to look at the
            % date used for the WRF file (which is just a month number for
            % v2 and monthly) and load that WRF file, then extract the
            % subset of NO2 data to plot with PCOLOR() behind the aircraft
            % data to plot with SCATTER() on top of it.
            match = wrfcomp.(prof_type).(campaign).(prof_number).match;
            buffer = 3;
            we_inds = (min(match.indicies.west_east)-buffer):(max(match.indicies.west_east)+buffer);
            sn_inds = (min(match.indicies.south_north)-buffer):(max(match.indicies.south_north)+buffer);
            if wrf_level == 0
                wrf_level = min(match.indicies.bottom_top);
            end
            
            % Convert from ppm to ppt
            conv = 1e6;
            
            [wrf_no2, wrf_lon, wrf_lat] = misc_behr_v3_validation.load_wrf_no2(prof_type, match.wrf.time, we_inds, sn_inds, wrf_level);
            wrf_no2 = wrf_no2 * conv;
            
            % To make pcolor match up with the actual positioning of the
            % WRF grid cells better, we'll use the bottom left corner,
            % because pcolor would plot wrf_no2(i,j) with corners
            % lon(i,j):lon(i+1,j+1) and likewise for lat.
            [wrf_loncorn, wrf_latcorn] = wrf_grid_corners(wrf_lon, wrf_lat);
            wrf_loncorn = squeeze(wrf_loncorn(1,:,:));
            wrf_latcorn = squeeze(wrf_latcorn(1,:,:));
            
            % Do 10-second averaging to reduce the number of points to plot
            n_sec_per_avg = 10;
            air_lon = avg_n_elements(match.data.lon, n_sec_per_avg, 'op', 'nanmean');
            air_lat = avg_n_elements(match.data.lat, n_sec_per_avg, 'op', 'nanmean');
            air_no2 = avg_n_elements(match.data.no2, n_sec_per_avg, 'op', 'nanmean');
            air_pres = avg_n_elements(match.data.pres, n_sec_per_avg, 'op', 'nanmean');
            
            air_no2 = air_no2 * conv;
            
            % Sort everything so that the highest altitudes (lowest
            % pressure) are plotted first, so that the lower altitudes are
            % plotted on top
            [air_pres, xx_sort] = sort(air_pres);
            air_lon = air_lon(xx_sort);
            air_lat = air_lat(xx_sort);
            air_no2 = air_no2(xx_sort);
            
            % Scale pressure to make the scatter a reasonable size.
            size_from_pres = scale_to_range(air_pres, [50 200]);
            
            figure;
            pcolor(wrf_loncorn, wrf_latcorn, wrf_no2);
            hold on;
            scatter(air_lon, air_lat, size_from_pres, air_no2, 'filled', 'markeredgecolor','k','linewidth',0.5);
            colorbar;
            set(gca,'fontsize',14);
            title(sprintf('%s: Profile %s vs. %s WRF', upper(strrep(campaign, '_', '-')), strrep(prof_number, 'p', ''), prof_type));
        end
        
        function plot_one_site_comparison(prof_types, campaign, site_number, uncert_type)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            prof_fns = fieldnames(wrfcomp.monthly.(campaign));
            tmp = strrep(prof_fns, 'p', '');
            site_numbers = floor(str2double(tmp)/1000); % DISCOVER profile numbers are either snnn or sssnnn where s or sss is the site number and nnn the profile number at that location
            allowed_site_numbers_strs = cellfun(@num2str, num2cell(unique(site_numbers)), 'uniformoutput', false);
            if ~exist('site_number', 'var')
                site_number = ask_multichoice('Which profile?', allowed_site_numbers_strs, 'list', true);
                site_number = str2double(site_number);
            else
                if ~isnumeric(site_number) || ~isscalar(site_number)
                    E.badinput('SITE_NUMBER must be a scalar number')
                end
                
                if isnan(site_number)
                    return
                elseif ~ismember(site_number, site_numbers)
                    E.badinput('SITE_NUMBER must be one of: %s', strjoin(allowed_site_numbers_strs));
                end
            end
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_types', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_types, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            %%% PLOTTING %%%
            conv = 1e6; % convert NO2 from ppm to ppt
            
            % First make a list of all profiles for a given site then
            % concatenate all of the matched data. We will then bin that
            % and plot it.
            xx = site_numbers == site_number;
            prof_fns = prof_fns(xx);
            for a=1:numel(prof_types)
                if a == 1
                    no2.aircraft = [];
                    pres.aircraft = [];
                end
                no2.(prof_types{a}) = [];
                pres.(prof_types{a}) = [];
                
                for b=1:numel(prof_fns)
                    if a==1
                        no2.aircraft = veccat(no2.aircraft, wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.data.no2*conv, 'column');
                        pres.aircraft = veccat(pres.aircraft, wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.data.pres, 'column');
                    end
                    no2.(prof_types{a}) = veccat(no2.(prof_types{a}), wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.wrf.no2*conv, 'column');
                    pres.(prof_types{a}) = veccat(pres.(prof_types{a}), wrfcomp.(prof_types{a}).(campaign).(prof_fns{b}).match.wrf.pres, 'column');
                end
                
                % Bin
                if a==1
                    [binned_no2.aircraft, binned_pres.aircraft, binned_no2_std.aircraft] = bin_omisp_pressure(pres.aircraft, no2.aircraft, 'mean');
                end
                [binned_no2.(prof_types{a}), binned_pres.(prof_types{a}), binned_no2_std.(prof_types{a})] = bin_omisp_pressure(pres.(prof_types{a}), no2.(prof_types{a}), 'mean');
            end
            
            plt_cols = misc_behr_v3_validation.plot_colors;
            
            fns = fieldnames(no2);
            l = gobjects(numel(fns),1);
            legstr = cell(1, numel(fns));
            figure;
            for a=1:numel(fns)
                l(a) = line(binned_no2.(fns{a}), binned_pres.(fns{a}), 'color', plt_cols.(fns{a}).avg, 'linewidth', 2);
                legstr{a} = capitalize_words(fns{a});
            end
            set(gca,'ydir','reverse','fontsize',14);
            legend(l, legstr);
            title(sprintf('%s: Site %d', upper(strrep(campaign,'_','-')), site_number));
        end
        
        function plot_all_sites_comparison(prof_types, campaign, uncert_type)
            E = JLLErrors;
            wrfcomp = load(misc_behr_v3_validation.wrf_comp_file);
            
            allowed_campaigns = fieldnames(wrfcomp.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('For which campaign?', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns));
            end
            
            prof_fns = fieldnames(wrfcomp.monthly.(campaign));
            tmp = strrep(prof_fns, 'p', '');
            site_numbers = unique(floor(str2double(tmp)/1000)); % DISCOVER profile numbers are either snnn or sssnnn where s or sss is the site number and nnn the profile number at that location
            
            allowed_prof_types = {'v2', 'monthly'};
            % Not all campaigns will have daily profiles, so only offer it
            % if it is available.
            if isfield(wrfcomp.daily, campaign)
                allowed_prof_types{end+1} = 'daily';
            end
            if ~exist('prof_type', 'var')
                prof_types = ask_multiselect('Which profile types to include?', allowed_prof_types);
            elseif ~ismember(prof_types, allowed_prof_types)
                E.badinput('PROF_TYPE must be one of: %s', strjoin(allowed_prof_types));
            end
            
            allowed_uncert_types = {'raw','std', 'both', 'none'};
            if ~exist('uncert_type', 'var')
                uncert_type = ask_multichoice('How to show the variability in the profiles?', allowed_uncert_types, 'default', 'std');
            elseif ~ismember(uncert_type, allowed_uncert_types)
                E.badinput('UNCERT_TYPE must be one of: %s', strjoin(allowed_uncert_types, ', '));
            end
            
            for a=1:numel(site_numbers)
                misc_behr_v3_validation.plot_one_site_comparison(prof_types, campaign, site_numbers(a), uncert_type);
            end
        end
        
        function varargout = plot_scd_vs_wrf_columns(varargin)
            % Plots OMI tropospheric SCDs, WRF monthly, and WRF daily
            % columns for a given date range to show whether WRF is
            % capturing the wind direction correctly. If daily BEHR files
            % do not exist for a date in the range given, that day will
            % just be skipped.
            
            E = JLLErrors;
            p = advInputParser;
            p.addOptional('location_name', '', @(x) ischar(x) || isstruct(x));
            p.addOptional('start_date', '', @(x) ischar(x) || isnumeric(x));
            p.addOptional('end_date', '', @(x) ischar(x) || isnumeric(x));
            p.addParameter('max_frac_nans', 0);
            p.addParameter('titles', true);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            location_name = pout.location_name;
            start_date = pout.start_date;
            end_date = pout.end_date;
            max_frac_nans = pout.max_frac_nans;
            include_titles = pout.titles;
            
            if ischar(location_name)
                locs = misc_behr_v3_validation.read_locs_file();
                loc_names = {locs.ShortName};
                if isempty(location_name)
                    loc_ind = ask_multichoice('Choose a location', loc_names, 'list', true, 'index', true);
                elseif ~ismember(location_name, loc_names)
                    E.badinput('LOCATION_NAME must be one of: %s', strjoin(loc_names, ', '));
                else
                    loc_ind = strcmpi(location_name, loc_names);
                end
                
                plot_loc = locs(loc_ind);
                
                if isempty(start_date)
                    start_date = datenum(ask_date('Enter the start date'));
                else
                    start_date = validate_date(start_date);
                end
                
                if isempty(end_date)
                    end_date = datenum(ask_date('Enter the end date'));
                else
                    end_date = validate_date(end_date);
                end
                
                dvec = start_date:end_date;
                load_data = true;
            elseif isstruct(location_name)
                % if given one input, it must be a structure of results
                if ~isstruct(location_name)
                    E.badinput('With one input, it must be a structure of results')
                end
                results = location_name;
                dvec = 1:numel(results);
                load_data = false;
            end
            
            for d=1:numel(dvec)
                % Load the daily BEHR file, if it exists
                if load_data
                    try
                        [behr_daily, wrf_monthly, wrf_daily] = misc_behr_v3_validation.load_wrf_and_behr_data(plot_loc, dvec(d), 'max_frac_nans', max_frac_nans);
                    catch err
                        if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
                            fprintf('No BEHR daily profile file available for %s\n', datestr(dvec(d)));
                            % Create empty structs in case we need to
                            % return after this date
                            behr_daily = struct('lon', {{}},'lat', {{}}, 'no2_scds', {{}}, 'no2_vcds', {{}});
                            wrf_monthly = struct('lon', {{}},'lat', {{}}, 'no2_vcds', {{}});
                            wrf_daily = struct('lon', {{}},'lat', {{}}, 'no2_vcds', {{}});
                            continue
                        else
                            rethrow(err)
                        end
                    end
                else
                    behr_daily = results(d).behr_data;
                    wrf_monthly = results(d).wrf_monthly;
                    wrf_daily = results(d).wrf_daily;
                    plot_loc.Longitude = results(d).loc_longitude;
                    plot_loc.Latitude = results(d).loc_latitude;
                end
                
                fig = [];
                for a=1:numel(behr_daily.lon)
                    
                    % Make the plot, 3 side-by-side figures
                    fig=figure;
                    fig.Position(3) = fig.Position(3)*2;
                    for p=1:3
                        subplot(1,3,p);
                        if p==1
                            % Convert VCD back to SCD
                            pcolor(behr_daily.lon{a}, behr_daily.lat{a}, behr_daily.no2_scds{a});
                            if include_titles
                                title('OMI Trop. SCD');
                            end
                        elseif p==2
                            pcolor(wrf_monthly.lon, wrf_monthly.lat, wrf_monthly.no2_vcds);
                            if include_titles
                                title('Monthly WRF');
                            end
                        elseif p==3
                            pcolor(wrf_daily.lon{a}, wrf_daily.lat{a}, wrf_daily.no2_vcds{a});
                            if include_titles
                                title('Daily WRF');
                            end
                        else
                            E.notimplemented('>3 plots')
                        end
                        
                        %caxis([0 1e16]);
                        cb=colorbar;
                        cb.Label.String = 'molec. cm^{-2}';
                        line(plot_loc.Longitude, plot_loc.Latitude, 'linestyle', 'none', 'marker', 'p', 'color', 'w', 'markersize',16,'linewidth',2);
                        set(gca,'fontsize',14);
                    end
                end
                
                
            end
            
            if nargout > 0
                if numel(dvec) > 1
                    E.notimplemented('Returning data when >1 day requested');
                end
                varargout = {behr_daily, wrf_monthly, wrf_daily, fig};
            end
        end
        
        function plot_campaign_profile_boxplot(campaign, include_error_bars)
            E = JLLErrors;
            
            profs = load(misc_behr_v3_validation.wrf_comp_file);
            allowed_campaigns = fieldnames(profs.monthly);
            if ~exist('campaign', 'var')
                campaign = ask_multichoice('Choose a campaign to plot', allowed_campaigns, 'list', true);
            elseif ~ismember(campaign, allowed_campaigns)
                E.badinput('CAMPAIGN must be one of: %s', strjoin(allowed_campaigns, ', '));
            end
            
            if ~exist('include_error_bars', 'var')
                include_error_bars = ask_yn('Plot with error bars?');
            elseif ~isscalar(include_error_bars) || ~islogical(include_error_bars)
                E.badinput('INCLUDE_ERROR_BARS must be a scalar logical')
            end
            
            have_daily = isfield(profs.daily, campaign);
            
            % Bin monthly, v2, and (if available) version 3 data
            [bins.air.no2, bins.air.pres, bins.air.quartiles] = misc_behr_v3_validation.make_boxplot_bins(profs.monthly.(campaign).All.match.data.pres, profs.monthly.(campaign).All.match.data.no2);
            [bins.monthly_wrf.no2, bins.monthly_wrf.pres, bins.monthly_wrf.quartiles] = misc_behr_v3_validation.make_boxplot_bins(profs.monthly.(campaign).All.match.wrf.pres, profs.monthly.(campaign).All.match.wrf.no2);
            [bins.v2_wrf.no2, bins.v2_wrf.pres, bins.v2_wrf.quartiles] = misc_behr_v3_validation.make_boxplot_bins(profs.v2.(campaign).All.match.wrf.pres, profs.v2.(campaign).All.match.wrf.no2);
            if have_daily
                [bins.daily_wrf.no2, bins.daily_wrf.pres, bins.daily_wrf.quartiles] = misc_behr_v3_validation.make_boxplot_bins(profs.daily.(campaign).All.match.wrf.pres, profs.daily.(campaign).All.match.wrf.no2);
            end
            
            % Set up offsets and legend strings
            if have_daily
                offsets.air = 0.25;
                offsets.daily_wrf = 0.12;
                offsets.monthly_wrf = -0.12;
                offsets.v2_wrf = -0.25;
            else
                offsets.air = 0.25;
                offsets.monthly_wrf = 0;
                offsets.v2_wrf = -0.25;
            end
            
            legend_strs = struct('air', 'Aircraft', 'daily_wrf', 'Daily WRF', 'monthly_wrf', 'Monthly WRF', 'v2_wrf', 'BEHR v2 WRF');
            
            if include_error_bars
                hPa_offset_factor = 10;
                line_opts = {'marker', 'o', 'linestyle', 'none', 'linewidth', 2, 'markersize', 8};
            else
                hPa_offset_factor = 0;
                line_opts = {'linestyle', '-', 'linewidth', 2, 'markersize', 8};
            end
            
            % Plot styles
            styles.air = struct('color', misc_behr_v3_validation.plot_colors.aircraft.avg);
            styles.daily_wrf = struct('color', misc_behr_v3_validation.plot_colors.daily.avg);
            styles.monthly_wrf = struct('color', misc_behr_v3_validation.plot_colors.monthly.avg);
            styles.v2_wrf = struct('color', misc_behr_v3_validation.plot_colors.v2.avg);
            
            fns = fieldnames(bins);
            l = gobjects(numel(fns),1);
            legstr = cell(1,numel(fns));
            
            unit_conv = 1e3; % convert ppm to ppb
            
            figure;
            for f=1:numel(fns)
                l(f) = line(bins.(fns{f}).no2 * unit_conv, bins.(fns{f}).pres + offsets.(fns{f}) * hPa_offset_factor, 'color', styles.(fns{f}).color, line_opts{:});
                if include_error_bars
                    scatter_errorbars(bins.(fns{f}).no2 * unit_conv, bins.(fns{f}).pres + offsets.(fns{f}) * hPa_offset_factor, bins.(fns{f}).quartiles(1,:), bins.(fns{f}).quartiles(2,:), 'direction', 'x', 'color', styles.(fns{f}).color);
                end
                legstr{f} = legend_strs.(fns{f});
            end
            ax = gca;
            ax.FontSize = 16;
            ax.YDir = 'reverse';
            ax.XLim(1) = 0;
            xlabel('[NO_2] (ppbv)');
            ylabel('Pressure (hPa)');
            legend(l,legstr);
            title(upper(strrep(campaign,'_','-')));
        end
        
        function plot_behr_wrf_surfpres_diffs()
            % Generated by generate_behr_wrf_surfpres_comparison
            Pres = load(misc_behr_v3_validation.pres_comp_file);
            figure;
            rdel = reldiff(Pres.all_behr_pres, Pres.all_wrf_pres)*100;
            [sp_y, sp_x] = square_subplot_dims(size(rdel,3));
            labels = cell(1, size(rdel,3));
            for m = 1:size(rdel, 3)
                subplot(sp_y, sp_x, m);
                pcolor(Pres.lon, Pres.lat, rdel(:,:,m));
                shading flat;
                state_outlines('k','not','ak','hi');
                cb = colorbar;
                cb.Label.String = '%\Delta Surf. P. (hPa)';
                caxis(calc_plot_limits(rdel(:,:,m), 'pow10', 'diff'));
                labels{m} = datestr(datenum(2012, m, 1), 'mmm');
                title(labels{m});
            end
            colormap(blue_red_cmap);
            
            figure;
            boxplot(reshape(rdel,[],size(rdel,3)));
            set(gca,'XTickLabel',labels,'ygrid','on');
            ylabel('%\Delta Surf. P. (hPa, BEHR - WRF)');
        end
        
        function fig = plot_uncertainty_estimate(varargin)
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('titles', true);
            p.addParameter('region', 'us');
            p.addParameter('change_field', '')
            p.addParameter('normalize_by', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            include_titles = pout.titles;
            region = pout.region;
            change_field = opt_ask_multichoice('Which field to plot?', {'PercentChangeNO2', 'PercentChangeNO2Vis', 'PercentChangeAMF', 'PercentChangeAMFVis'}, pout.change_field, '"change_field"', 'list', true);
            normalize_by = opt_ask_multichoice('Normalize the uncertainty in the bar plot by:', {'None', 'Percent difference', 'Absolute difference'}, pout.normalize_by, '"normalize_by"', 'list', true);
            
            % First we need to get what uncertainty parameters are
            % available and for which months. The months should be the same
            % for all parameters. We'll load the uncertainty at the same
            % time
            
            F = dir(behr_paths.BEHRUncertSubdir(region));
            uncert_params = {F([F.isdir]).name};
            % Need to remove ., .., and BaseCase
            uncert_params(regcmp(uncert_params, '\.+')) = [];
            uncert_params(strcmp(uncert_params, 'BaseCase')) = [];
            
            % Assume for now that there's four month of uncertainty
            uncertainties = make_empty_struct_from_cell(uncert_params);
            
            plot_param_in_bar_graph = true(size(uncert_params));
            switch lower(normalize_by)
                case 'none'
                    input_diff_fxn = @(A,B,denom) ones(size(A));
                    plot_in_bar_if_cant_norm = true;
                    bar_y_label = '% uncertainty';
                case 'percent difference'
                    input_diff_fxn = @(A,B,denom) reldiff(A,B)*100./denom;
                    plot_in_bar_if_cant_norm = false;
                    bar_y_label = '%\Delta VCD / %\Delta input';
                case 'absolute difference'
                    input_diff_fxn = @(A,B,denom) (A - B)./denom;
                    plot_in_bar_if_cant_norm = false;
                    bar_y_label = '%\Delta VCD / \Delta input';
                otherwise
                    E.notimplemented('No calculation implemented for normalize_by = "%s"', normalize_by)
            end
            
            for i_param = 1:numel(uncert_params)
                fprintf('Loading %s files: ', uncert_params{i_param});
                
                uncert_files_tmp = dirff(fullfile(behr_paths.BEHRUncertSubdir(region), uncert_params{i_param}, 'BEHR*.mat'));
                % Order the files DJF, MAM, JJA, SON. The default
                % alphabetical order is very strange to see in the plot.
                if numel(uncert_files_tmp) ~= 4
                    E.notimplemented('Not 4 uncertainty files')
                end
                regexes = {'DJF', 'MAM', 'JJA', 'SON'};
                uncert_names = {uncert_files_tmp.name};
                for i_re = 1:numel(regexes)
                    uncert_files(i_re) = uncert_files_tmp(regcmp(uncert_names, regexes{i_re}));
                end
                
                file_dates = cellfun(@(x) regexp(x, '(DJF|MAM|JJA|SON)', 'match', 'once'), {uncert_files.name}, 'UniformOutput', false);
                if i_param == 1
                    check_file_dates = file_dates;
                elseif ~isequal(file_dates, check_file_dates)
                    E.notimplemented('Different parameters produced for different months');
                end
                
                substruct = make_empty_struct_from_cell({'date', 'percent_diff', 'norm_factor'}, cell(size(uncert_files)));
                
                for i_file = 1:numel(uncert_files)
                    UAvg = load(uncert_files(i_file).name);
                    ErrorAvg = UAvg.ErrorAvg;
                    
                    fprintf('%d ', i_file);
                    
                    error_fields = fieldnames(ErrorAvg);
                    xx_base_fn = regcmp(error_fields, 'Base$');
                    if sum(xx_base_fn) == 1
                        base_field = error_fields{xx_base_fn};
                        perturbed_field = regexprep(base_field, 'Base$', '');
                        this_input_diff_fxn = input_diff_fxn;
                        % If we can normalize this parameter, always plot
                        % it in the bar graph
                    elseif sum(xx_base_fn) == 0
                        base_field = change_field;
                        perturbed_field = change_field;
                        this_input_diff_fxn = @(A,B,denom) ones(size(A));
                        % If we can't normalize this parameter, whether we
                        % plot it depends on if we're trying to normalize
                        % or not
                        plot_param_in_bar_graph(i_param) = plot_in_bar_if_cant_norm;
                    else
                        E.callError('multiple_base_fields','Multiple fields found matching the regular expression "Base$"');
                    end
                    
                    substruct(i_file).date = file_dates(i_file);
                    if numel(ErrorAvg) == 1
                        perdiff = ErrorAvg.(change_field);
                        input_diff = this_input_diff_fxn(ErrorAvg.(perturbed_field), ErrorAvg.(base_field), 1);
                    elseif numel(ErrorAvg) == 2
                        % If we calculated a percent difference by raising
                        % and lowering the perturbed value, this will
                        % effectively calculate a mean percent difference.
                        % I do it this way rather than summing and dividing
                        % by two or summing the absolute value and dividing
                        % by two to keep the sign. This assumes implicitly
                        % that the change w.r.t. the input parameter is
                        % monotonic; if, e.g. increasing or decreasing the
                        % input parameter both increased the NO2 VCDs, then
                        % this method will calculate a reduced uncertainty
                        % than if we summed or just took one side. However,
                        % I argue that is both unlikely and, even if it
                        % happens, correct because the range over which the
                        % VCDs varies is smaller in that case than if it
                        % was monotonically increasing or decreasing with
                        % the input parameter.
                        perdiff = (ErrorAvg(2).(change_field) - ErrorAvg(1).(change_field))./2;
                        input_diff = this_input_diff_fxn(ErrorAvg(2).(perturbed_field), ErrorAvg(1).(perturbed_field), 2);
                    else
                        E.notimplemented('numel(ErrorAvg) > 2')
                    end
                    % Sometime the perturbation produces an absurdly high
                    % uncertainty. Restrict to 3 sigma.
                    %perdiff(abs(perdiff) > 3*nanstd(perdiff)) = nan;
                    substruct(i_file).percent_diff = perdiff;
                    substruct(i_file).norm_factor = input_diff;
                end
                uncertainties.(uncert_params{i_param}) = substruct;
                fprintf('\n');
            end
            
            % Whether we include the total in the bar graph depends on
            % whether we're including parameters than can't be normalized.
            plot_param_in_bar_graph(end+1) = plot_in_bar_if_cant_norm;
            
            % Now we want to plot the quadrature sum of all contributions
            % to the uncertainty for each month, and a plot of the average
            % contribution of each component. For plotting, will have to
            % assume that lat/lon is the same as the normal BEHR grid
            [~, OMI] = load_behr_file('2012-01-01', 'monthly', region);
            lon = OMI(1).Longitude;
            lat = OMI(1).Latitude;
            
            n_months = numel(check_file_dates);
            n_params = numel(uncert_params);
            
            median_contributions = nan(n_months, n_params+1);
            contribution_quantiles = nan([2,size(median_contributions)]);
            mean_contributions = nan(size(median_contributions));
            contribution_sigmas = nan(size(median_contributions));
            
            fig = figure;
            fig.Position(3:4) = [2 3] .* fig.Position(3:4);
            
            if n_months ~= 4
                E.notimplemented('# of months ~= 4');
            end
            for i_month = 1:n_months
                pAvg = RunningAverage();
                for i_param = 1:n_params
                    this_struct = uncertainties.(uncert_params{i_param})(i_month);
                    if ~strcmp(this_struct.date, check_file_dates{i_month})
                        E.callError('wrong_date', 'Dates in substruct out of order compared to what was expected')
                    end
                    pAvg.addData((this_struct.percent_diff).^2);
                    this_norm_perdiff = this_struct.percent_diff(:) ./ this_struct.norm_factor(:);
                    mean_contributions(i_month, i_param) = abs(nanmean(this_norm_perdiff));
                    contribution_sigmas(i_month, i_param) = nanstd(this_norm_perdiff);
                    median_contributions(i_month, i_param) = nanmedian(abs(this_norm_perdiff));
                    contribution_quantiles(:, i_month, i_param) = quantile(abs(this_norm_perdiff),[0.05;0.95]);
                end
                subplot(3,2,i_month);
                perdiff = sqrt(pAvg.values());
                mean_contributions(i_month, end) = nanmean(abs(perdiff(:)));
                contribution_sigmas(i_month, end) = nanstd(abs(perdiff(:)));
                median_contributions(i_month, end) = nanmedian(abs(perdiff(:)));
                contribution_quantiles(:, i_month, end) = quantile(abs(perdiff(:)), [0.05; 0.95]);
                
                % This would make the plots look nicer
                %perdiff(isnan(perdiff)) = 0;
                pcolor(lon, lat, perdiff);
                shading flat
                caxis(calc_plot_limits(perdiff, 'zero', 'max', [0 200]));
                cb = colorbar;
                cb.Label.String = 'Total % Uncertainty';
                set(gca, 'fontsize', 14);
                state_outlines('w');
                if include_titles
                    title(check_file_dates{i_month});
                end
                label_subfigs(gcf, 'xshift', 0.2)
            end
            
            ax = subplot(3,2,n_months+1);
            %bar(median_contributions(:, plot_param_in_bar_graph));
            %bar_errors(median_contributions(:, plot_param_in_bar_graph), squeeze(contribution_quantiles(1,:,plot_param_in_bar_graph)), squeeze(contribution_quantiles(2,:,plot_param_in_bar_graph)));
            
            bar(mean_contributions(:, plot_param_in_bar_graph));
            %bar_errors(mean_contributions(:, plot_param_in_bar_graph), contribution_sigmas(:,plot_param_in_bar_graph));
            
            set(gca, 'ygrid', 'on', 'xticklabel', check_file_dates);
            legend_parms = veccat(uncert_params, 'Total');
            legend(legend_parms{plot_param_in_bar_graph}, 'Location', 'eastoutside');
            % Expand this plot to fill the width of the figure with some
            % room for the legend
            ax.Position(3) = 1.8*ax.Position(3);
            %ch = get(gcf,'children');
            %center_axes(ax, ch(4), ch(6));
            set(gca,'fontsize',14);
            ylabel(bar_y_label)
        end
        
        function fig = plot_prof_regression(varargin)
            % This function will plot the R2 value at each model level for
            % either the profile concentration or shape factor between WRF
            % and aircraft profiles
            %
            %   'campaigns' - which campaigns to include, as a cell array
            %
            %   'versions' - which versions to include (v2, monthly, daily)
            %   as a cell array of strings
            %
            %   'title' - boolean, whether or not to include titles (which
            %   will be the campaign name). Default is true.
            %
            %   'with_diffs' - boolean, whether or not to include bar
            %   graphs of the differences in R2 if exactly two versions
            %   requested. Default is true.
            %
            %   'num_pts' - boolean, whether to include the number of
            %   comparisons that went into each R2 calculation on the
            %   difference plot. Has no effect if numel(versions) ~= 2.
            %   Default is true.
            
            p = advInputParser;
            p.addParameter('campaigns', '');
            p.addParameter('versions', {});
            p.addParameter('title',true);
            p.addParameter('with_diffs',true);
            p.addParameter('num_pts',true);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            wrf_comp = load(misc_behr_v3_validation.wrf_comp_file);
            available_campaigns = fieldnames(wrf_comp.v2);
            
            include_title = pout.title;
            include_diffs = pout.with_diffs;
            include_num_pts = pout.num_pts;
            
            campaigns = opt_ask_multiselect('Which campaigns to plot?', available_campaigns, pout.campaigns, '"campaign"');
            
            available_versions = {'v2','monthly'};
            if isfield(wrf_comp.daily, campaigns)
                available_versions = veccat(available_versions, 'daily');
            end
            versions = opt_ask_multiselect('Which versions to include?', available_versions, pout.versions, '"versions"');
            
            include_diffs = include_diffs && numel(versions) == 2;
            
            plot_colors = misc_behr_v3_validation.plot_colors;
            plot_markers = misc_behr_v3_validation.plot_markers;
            plot_leg_str = misc_behr_v3_validation.versions2legend(versions);
            
            fig = figure;
            if include_diffs
                subplot_x = 2;
            else
                subplot_x = 1;
            end
            subplot_y = numel(campaigns);
            
            % We only want to label the first axes in each row, not the
            % differences.
            axes_to_label = gobjects(subplot_y, 1);
            
            for i_cam = 1:numel(campaigns)
                this_campaign = campaigns{i_cam};
                
                if include_diffs
                    subplot_ind = (i_cam-1)*2 + 1;
                else
                    subplot_ind = i_cam;
                end
                
                left_ax = subplot(subplot_y, subplot_x, subplot_ind);
                axes_to_label(i_cam) = left_ax;
                
                r2 = [];
                num_pts = [];
                for i_version = 1:numel(versions)
                    this_version = versions{i_version};
                    [r2(:,i_version), pres_bins, num_pts(:,i_version)] = misc_behr_v3_validation.calculate_profile_regressions(wrf_comp.(this_version).(this_campaign));
                end
                
                % Plot the individual levels' R2 values, only including
                % levels where all versions have an R2 value. Some may not
                % have an R2 value in some levels due to differences in
                % pressure for each WRF point between the versions.
                notnans = all(~isnan(r2),2);
                bin_number = (1:size(r2,1))';
                plot_ylim = [find(notnans,1,'first')-1, find(notnans,1,'last')+1];
                i_lastval = find(notnans,1,'last');
                l = gobjects(numel(versions),1);
                for i_version = 1:numel(versions)
                    this_version = versions{i_version};
                    if plot_markers.(this_version).filled
                        extra_args = {'markerfacecolor', plot_colors.(this_version).avg};
                    else
                        extra_args = {};
                    end
                    % Plot using model level as the y-coordinate to avoid
                    % squishing up the bottom levels that are closer in
                    % pressure.
                    l(i_version) = line(r2(notnans,i_version), bin_number(notnans), 'color', plot_colors.(this_version).avg, 'marker', plot_markers.(this_version).avg, 'linestyle','none','markersize',10,'linewidth',2, extra_args{:});
                end
                
                % Only put the legend on the first plot.
                if i_cam == 1
                    legend(plot_leg_str)
                end
                % Since we plotted with model level as the y-coordinate,
                % redo the ticks to be pressure. I just randomly picked
                % every 4 levels to give a decent spacing.
                set(gca,'ytick',bin_number(1:4:i_lastval),'yticklabel',pres_bins(1:4:i_lastval), 'fontsize', 16, 'ygrid', 'on');
                xlabel('R^2');
                ylabel('Pressure (hPa)');
                if include_title
                    title(upper(strrep(this_campaign,'_','-')))
                end
                
                if include_num_pts
                    combined_num_pts_strings = cell(size(num_pts,1),1);
                    for i_level = 1:size(num_pts,1)
                        combined_num_pts_strings{i_level} = sprintf('(%d, %d)', num_pts(i_level,1), num_pts(i_level,2));
                    end
                end
                
                % If we're plotting exactly 2 versions and we want to
                % include diffs, add them now.
                if numel(versions) == 2 && include_diffs
                    set(gca,'ylim', plot_ylim);
                    right_ax = subplot(subplot_y, subplot_x, subplot_ind+1);
                    r2_diffs = r2(:,2) - r2(:,1);
                    pos_diffs = r2_diffs > 0;
                    positive_bars = barh(bin_number(pos_diffs), r2_diffs(pos_diffs), plot_colors.(versions{2}).avg);
                    hold on
                    negative_bars = barh(bin_number(~pos_diffs), r2_diffs(~pos_diffs), plot_colors.(versions{1}).avg);
                    
                    if include_num_pts
                        text(positive_bars.YData*1.1, positive_bars.XData, combined_num_pts_strings(positive_bars.XData));
                        text(negative_bars.YData*1.1, negative_bars.XData, combined_num_pts_strings(negative_bars.XData), 'HorizontalAlignment', 'right');
                        plot_xlims = calc_plot_limits(r2_diffs*1.5, 0.2);
                    else
                        plot_xlims = calc_plot_limits(r2_diffs, 0.2);
                    end
                    
                    set(gca, 'fontsize', 16, 'ylim', plot_ylim, 'ytick', [], 'xlim', plot_xlims);
                    xlabel(sprintf('R^2 %s - R^2 %s', plot_leg_str{2}, plot_leg_str{1}));
                    
                    % Shift the axes so that the profiles have more space than
                    % the differences, and make sure their size in the
                    % y-dimension is the same
                    shift = 0.1;
                    right_ax.Position = right_ax.Position + [shift, 0, -shift, 0];
                    right_ax.Position([2 4]) = left_ax.Position([2 4]);
                    left_ax.Position(3) = left_ax.Position(3) + shift;
                end
            end
            
            % Give the plots a little more size
            fig.Position(3:4) = fig.Position(3:4) .* [1.5 subplot_y*1.2];
            
            if ~include_title
                label_subfigs(gcf, 'xshift', 0.2, 'ax', axes_to_label);
            end
        end
        
        function [r2_vals, pres_bins, num_pts] = calculate_profile_regressions(prof_struct, varargin)
            % prof_struct should be one campaign's structure from the wrf
            % comparison file.
            
            p = advInputParser;
            % y-resid and rma seem to return the same R^2 values
            p.addParameter('reg_type', 'y-resid');
            % If I need to add shape factors for a reviewer (because they
            % seem to love shape factors) probably the simplest way would
            % be to calculate a partial VCD just for the extent that the
            % aircraft covers
            
            p.parse(varargin{:});
            pout = p.Results;
            
            reg_type = pout.reg_type;
            
            E = JLLErrors;
            % First, get the list of individual profiles. If there are
            % none, throw and error b/c this only works for DISCOVER-type
            % campaigns
            fns = fieldnames(prof_struct);
            fns = fns(regcmp(fns, 'p\d+'));
            if isempty(fns)
                E.badinput('Must pass a profile structure containing at least one field named p#### or p######')
            end
            n_profs = numel(fns);
            
            % Now we need to loop through each profile field and bin the
            % WRF and aircraft data for that profile into pressure bins. If
            % this is the first profile, use the size of the bin vectors to
            % initialize the full array for the profiles.
            for i_prof = 1:n_profs
                wrf_pres = prof_struct.(fns{i_prof}).match.wrf.pres;
                wrf_no2 = prof_struct.(fns{i_prof}).match.wrf.no2;
                air_pres = prof_struct.(fns{i_prof}).match.data.pres;
                air_no2 = prof_struct.(fns{i_prof}).match.data.no2;
                
                [wrf_no2_binned, pres_bins] = bin_omisp_pressure(wrf_pres, wrf_no2, 'mean');
                air_no2_binned = bin_omisp_pressure(air_pres, air_no2, 'mean');
                
                if ~isequal(size(wrf_no2_binned), size(air_no2_binned))
                    E.notimplemented('wrf and aircraft binned NO2 are different sizes')
                end
                
                if i_prof == 1
                    all_wrf_no2_binned = nan(numel(wrf_no2_binned), n_profs);
                    all_air_no2_binned = nan(numel(air_no2_binned), n_profs);
                end
                
                all_wrf_no2_binned(:, i_prof) = wrf_no2_binned(:);
                all_air_no2_binned(:, i_prof) = air_no2_binned(:);
            end
            
            % Now, finally, calculate the correlation for each level
            n_levels = size(all_wrf_no2_binned,1);
            r2_vals = nan(n_levels, 1);
            num_pts = nan(n_levels, 1);
            warning('off', 'fitting:nans_removed')
            for i_level = 1:n_levels
                [~,~,~,line_data] = calc_fit_line(all_air_no2_binned(i_level, :), all_wrf_no2_binned(i_level, :), 'regression', reg_type);
                r2_vals(i_level) = line_data.R2;
                num_pts(i_level) = line_data.num_pts;
            end
            warning('on', 'fitting:nans_removed')
        end
        
        function figs = plot_vcd_prof_ensembles(varargin)
            p = advInputParser;
            p.addParameter('extend_methods', {});
            p.addParameter('prof_type', '');
            p.addParameter('campaigns', {});
            p.addParameter('title', true);
            p.parse(varargin{:});
            pout = p.Results;
            
            extend_methods = opt_ask_multiselect('Which extend method(s) to include?', misc_behr_v3_validation.profile_extend_methods, pout.extend_methods, '"extend_methods"');
            prof_type = opt_ask_multichoice('Which a priori profile type to use?', {'monthly','daily'}, pout.prof_type, '"prof_type"', 'list', true, 'default', 'monthly');
            
            % We need to load at least one comparison structure to figure
            % out which campaigns are available
            comp_structs = make_empty_struct_from_cell(extend_methods);
            for i_extend = 1:numel(extend_methods)
                comp_structs.(extend_methods{i_extend}) = load(misc_behr_v3_validation.profile_comp_file(extend_methods{i_extend})); %#ok<LOAD>
            end
            
            available_campaigns = misc_behr_v3_validation.list_available_campaigns(comp_structs.(extend_methods{1}), prof_type);
            campaigns = opt_ask_multiselect('Which campaign(s) to include?', available_campaigns, pout.campaigns, '"campaigns"');
            
            include_title = pout.title;
            
            prof_fn = sprintf('us_%s', prof_type);
            figs = gobjects(size(campaigns));
            extend_plot_fmts = struct('wrf', struct('linecolor', 'r', 'envcolor', [1 0.5 0], 'legend_str', 'WRF-Chem'),...
                'geos', struct('linecolor', 'b', 'envcolor', 'c', 'legend_str', 'GEOS-Chem'),...
                'extrap', struct('linecolor', [0 0.5 0], 'envcolor', 'g', 'legend_str', 'Extrapolation'));
            
            for i_cam = 1:numel(campaigns)
                this_campaign = campaigns{i_cam};
                figs(i_cam) = figure;
                
                l = gobjects(numel(extend_methods),1);
                legend_str = cell(1,numel(extend_methods));
                for i_extend = 1:numel(extend_methods)
                    this_extend = extend_methods{i_extend};
                    plot_fmt = extend_plot_fmts.(this_extend);
                    prof_color = plot_fmt.linecolor;
                    env_color = plot_fmt.envcolor;
                    details_struct = comp_structs.(this_extend).v3.(prof_fn).(this_campaign).t1200_1500.details;
                    prof_ens = cat(1, details_struct.binned_profile)'*1e9;
                    % These should be all at the same level, but this will
                    % inherently deal with empty pres values
                    prof_pres = nanmean(cat(1, details_struct.binned_pressure)',2);
                    
                    [prof_mean, ~, ~, l(i_extend)] = plot_profile_envelope(prof_ens, prof_pres, 'quantile', [0.1 0.9], 'linecolor', prof_color, 'envcolor', env_color, 'envalpha', 0.5);
                    legend_str{i_extend} = plot_fmt.legend_str;
                    was_app_mat = cat(1, details_struct.is_appended_or_interpolated)';
                    xx_appended = sum(was_app_mat == 1,2)/size(was_app_mat, 2) > 0.5;
                    prof_mean(~xx_appended) = nan;
                    line(prof_mean, prof_pres, 'color', prof_color, 'linestyle', 'none', 'marker', 'o', 'linewidth', 1);
                    %legend_str{by_2_ind+1} = sprintf('%s (extension)', this_extend);
                end
                legend(l, legend_str);
                xlabel('[NO_2] (ppbv)');
                ylabel('Pressure (hPa)');
                set(gca,'ydir','reverse','fontsize',16);
                xlim(max(get(gca,'xlim'),0));
                ylim([0 1050]);
                if include_title
                    title(this_campaign);
                end
            end
        end
        
        function plot_campaign_boxplots(varargin)
            p = advInputParser;
            p.addParameter('campaigns', {});
            p.addParameter('extend_method', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            aircraft = load(misc_behr_v3_validation.profile_comp_file(pout.extend_method));
            pandora = load(misc_behr_v3_validation.pandora_comp_file);
            
            campaigns = opt_ask_multiselect('Which campaigns?', fieldnames(pandora.v3.us_monthly), pout.campaigns, '"campaigns"');
            all_versions = {'v2','v3'};
            
            air_labels = {'Air SP %s','SP %s-air', 'Air BEHR %s', 'BEHR %s-air'};
            daily_labels = {'Air BEHR %s (D)', 'BEHR %s (D)-air'};
            pandora_labels = {'Pandora-%s', 'SP %s-Pan', 'BEHR %s-Pan'};
            daily_pan_labels = {'BEHR %s (D)-Pan'};
            for i_cam = 1:numel(campaigns)
                this_campaign = campaigns{i_cam};
                for i_vers = 1:numel(all_versions)
                    version = all_versions{i_vers};
                    this_air = aircraft.(version).us_monthly.(this_campaign).t1200_1500;
                    data = {this_air.air_no2_nasa, this_air.sp_no2, this_air.air_no2_behr, this_air.behr_no2};
                    labels = air_labels;
                    
                    is_daily = isfield(aircraft.(version), 'us_daily') && ~isempty(aircraft.(version).us_daily.(this_campaign).t1200_1500);
                    if is_daily
                        this_daily_air = aircraft.(version).us_daily.(this_campaign).t1200_1500;
                        data = veccat(data, {this_daily_air.air_no2_behr, this_daily_air.behr_no2});
                        labels = veccat(labels, daily_labels);
                    end
                    
                    this_pan = pandora.(version).us_monthly.(this_campaign).t1230_1430;
                    data = veccat(data, {this_pan.pandora_no2,  this_pan.sp_no2, this_pan.behr_no2});
                    labels = veccat(labels, pandora_labels);
                    if is_daily
                        this_daily_pan = pandora.(version).us_daily.(this_campaign).t1230_1430;
                        data = veccat(data, {this_daily_pan.behr_no2});
                        labels = veccat(labels, daily_pan_labels);
                    end
                    
                    labels = cprintf(labels, version);
                    misc_behr_v3_validation.boxplot_helper(data, labels);
                    title(upper(strrep(this_campaign, '_', '-')));
                end
            end
        end
        
        function fig = boxplot_helper(x_data, labels)
            xall = [];
            gall = [];
            for i_dat = 1:numel(x_data)
                x = x_data{i_dat}(:);
                xall = cat(1, xall, x);
                gall = cat(1, gall, i_dat * ones(size(x)));
            end
            
            fig = figure;
            boxplot(xall,gall,'labels',labels);
            ylabel('VCD (molec cm^{-2})');
            set(gca,'ygrid','on','fontsize',14,'XTickLabelRotation',30)
        end
        
        
        function [unwt_table, wt_table] = are_avg_slopes_different(varargin)
            E = JLLErrors;
            
            p = advInputParser;
            p.addParameter('avg_or_combo', '');
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            pout = p.Results;
            
            avg_or_combo = opt_ask_multichoice('Use the average or combined aircraft and Pandora slopes?', {'avg', 'combo'}, pout.avg_or_combo, '"avg_or_combo"', 'list', true);
            
            switch lower(avg_or_combo)
                case 'avg'
                    [slopes_and_stds, ~, rownames, samples] = misc_behr_v3_validation.make_aircraft_pandora_avg(varargin{:}, 'plot_or_table', 'table', 'return_as_table', false);
                case 'combo'
                    [slopes_and_stds, ~, rownames, samples] = misc_behr_v3_validation.make_aircraft_pandora_combo(varargin{:}, 'return_as_table', false);
                otherwise
                    E.notimplemented('No action for case = %s', avg_or_combo);
            end
            slopes = slopes_and_stds(:, 1:2:end);
            stddevs = slopes_and_stds(:,2:2:end);
            
            comp_mat = build_comp_mat(rownames);
            labels = cell(size(rownames,1), 1);
            for i=1:numel(labels)
                labels{i} = strjoin(rownames(i,:), ' ');
                labels{i} = regexprep(capitalize_words(labels{i}), '\W', '');
            end
            
            unwt_table = misc_behr_v3_validation.display_slope_sig_diff(slopes(:,1), stddevs(:,1), samples(:,1), 'labels', labels, 'do_comp', comp_mat);
            wt_table = misc_behr_v3_validation.display_slope_sig_diff(slopes(:,2), stddevs(:,2), samples(:,2), 'labels', labels, 'do_comp', comp_mat);
            
            function mat = build_comp_mat(rows)
                % We want to compare three/five pairs within each campaign:
                % SP v3.0 <-> BEHR v3.0 (M/D), BEHR v2.1C <-> BEHR v3.0
                % (M/D), and BEHR v3.0 (M) <-> BEHR v3.0 (D)
                mat = false(size(rows,1));
                pairs = {'SP v3\.0', 'BEHR v3\.0. \(M\)';...
                    'SP v3\.0', 'BEHR v3\.0. \(D\)';...
                    'BEHR v2\.1.', 'BEHR v3\.0. \(M\)';...
                    'BEHR v2\.1.', 'BEHR v3\.0. \(D\)';...
                    'BEHR v3\.0. \(M\)', 'BEHR v3\.0. \(D\)'};
                u_cams = unique(rows(:,1));
                for i_cam = 1:numel(u_cams)
                    for i_pair = 1:size(pairs,1)
                        [pair_i, pair_j] = find_pair_inds(rows, u_cams{i_cam}, pairs{i_pair, 1}, pairs{i_pair, 2});
                        if ~isempty(pair_i) && ~isempty(pair_j)
                            mat(pair_i, pair_j) = true;
                        end
                    end
                end
            end
            
            function [i,j] = find_pair_inds(rows, campaign, first_name, second_name)
                xx_cam = strcmp(rows(:,1), campaign);
                xx_1 = regcmp(rows(:,2), first_name);
                xx_2 = regcmp(rows(:,2), second_name);
                i = find(xx_cam & xx_1);
                j = find(xx_cam & xx_2);
                
                if numel(i) > 1 || numel(j) > 1
                    E.callError('pair_matching:multiple_found', 'More than one match found for a pair of products');
                end
            end
        end
        
        function out = display_slope_sig_diff(slopes, std_devs, n_pts, varargin)
            % Given a list 
            
            p = advInputParser;
            p.addParameter('do_comp', []);
            p.addParameter('labels', {});
            p.addParameter('out_format', 'table');
            p.parse(varargin{:});
            pout = p.Results;
            
            do_comp_matrix = pout.do_comp;
            labels = pout.labels;
            out_format = pout.out_format;
            
            if isempty(do_comp_matrix)
                do_comp_matrix = true(numel(slopes));
            end
            if isempty(labels)
                labels = num2cell(1:numel(slopes));
                labels = cellfun(@num2str, labels, 'uniformoutput', false);
            end
            
            allowed_out_fmts = {'table', 'array', 'plot'};
            if ~ismember(out_format, allowed_out_fmts)
                E.badinput('"out_format" can only be one of: %s', strjoin(allowed_out_fmts, ', '));
            end
            
            is_diff = nan(numel(slopes));
            
            for i=1:numel(slopes)
                for j=i+1:numel(slopes)
                    if do_comp_matrix(i,j) || do_comp_matrix(j,i)
                        is_diff(i,j) = slope_significant_difference( slopes(i), std_devs(i), n_pts(i), slopes(j), std_devs(j), n_pts(j));
                    end
                end
            end
            
            switch lower(out_format)
                case 'table'
                    str_cell = cell(size(is_diff));
                    for i_cell = 1:numel(is_diff)
                        if is_diff(i_cell) == 0
                            str_cell{i_cell} = 'n';
                        elseif is_diff(i_cell) == 1
                            str_cell{i_cell} = 'y';
                        end
                    end
                    out = cell2table(str_cell, 'RowNames', labels, 'VariableNames', labels);
                case 'array'
                    out = is_diff;
                case 'plot'
                    E.notimplemented('plotting')
                otherwise
                    E.notimplemented(out_format)
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions %
        %%%%%%%%%%%%%%%%%%%%%%
        
        function results = record_scd_comparison()
            nsites = 10;
            ndays = 5;
            dvec = datenum('2012-01-01'):datenum('2012-12-31');
            locs = misc_behr_v3_validation.read_locs_file;
            
            allowed_site_types = {'Cities', 'PowerPlants'};
            %site_types = ask_multiselect('Which site types to use?', allowed_site_types);
            
            % Remove the undesired site types
            %xx = ismember({locs.SiteType}, site_types);
            %locs = locs(xx);
            xx = ismember({locs.ShortName}, {'Atlanta','Chicago','Las Vegas','Los Angeles','New York','Four Corners'});
            %xx = ismember({locs.ShortName}, {'Four Corners'});
            locs = locs(xx);
            
            nsites = min(nsites, numel(locs));
            
            results = struct('loc_name', '', 'date', '', 'user_value', [], 'user_confidence', [], 'user_note', '');
            results = repmat(results, nsites*ndays, 1);
            
            eval_opts = {'Daily - good agreement with SCD', 'Daily - bad agreement with SCD', 'Similar to monthly - good agreement with SCD', 'Similar to monthly - bad agreement with SCD', 'Not enough data'};
            ned_ind = numel(eval_opts); % index for "Not enough data", if this is the response, we don't want to store the result. Should always be the last one.
            for a=1:nsites
                % Make a copy so that we can remove days as we test them
                isite = randi(numel(locs), 1);
                site_dvec = dvec;
                b = 1;
                while b <= ndays && ~isempty(site_dvec)
                    idate = randi(numel(site_dvec), 1);
                    this_date = site_dvec(idate);
                    
                    [behr, wrf_monthly, wrf_daily] = misc_behr_v3_validation.plot_scd_vs_wrf_columns(locs(isite).ShortName, site_dvec(idate), site_dvec(idate), 'max_frac_nans', 0.1);
                    
                    if numel(behr.lon) > 0
                        try
                            user_ans = ask_multichoice('Evaluate this day', eval_opts, 'list', true, 'index', true);
                        catch err
                            if strcmp(err.identifier, 'ask_multichoice:user_cancel')
                                % If we cancel the work, save the results
                                % completed so far.
                                misc_behr_v3_validation.save_scd_results(misc_behr_v3_validation.scd_comp_file(true), results, eval_opts(1:end-1));
                                return
                            else
                                rethrow(err)
                            end
                        end
                    else
                        fprintf('Skipping %s %s because not enough SCD data\n', locs(isite).ShortName, datestr(this_date));
                    end
                    % Whether or not we use this day, we don't want to
                    % repeat it, and we want to close the figures
                    close all
                    
                    site_dvec(idate) = [];
                    if numel(behr.lon) == 0 || user_ans == ned_ind
                        % If the user responded "Not enough data", then
                        % don't store any result.
                        continue
                    end
                    
                    % Ask the user to rate their confidence
                    confidence = ask_number('Rate confidence (1 low to 3 high)', 'testfxn', @(x) isscalar(x) && x >= 1 && x <= 3, 'testmsg', 'Enter 1-3');
                    
                    iresult = sub2ind([ndays, nsites], b, a);
                    results(iresult).loc_name = locs(isite).ShortName;
                    results(iresult).date = datestr(this_date);
                    results(iresult).loc_longitude = locs(isite).Longitude;
                    results(iresult).loc_latitude = locs(isite).Latitude;
                    results(iresult).user_value = user_ans;
                    results(iresult).user_note = input('Enter a note if you wish: ', 's');
                    results(iresult).user_confidence = confidence;
                    results(iresult).behr_data = behr;
                    results(iresult).wrf_monthly = wrf_monthly;
                    results(iresult).wrf_daily = wrf_daily;
                    b = b + 1;
                    fprintf('%d of %d completed...\n', iresult, nsites*ndays);
                end
                
                locs(isite) = [];
                
            end
            
            misc_behr_v3_validation.save_scd_results(misc_behr_v3_validation.scd_comp_file, results, eval_opts(1:end-1));
        end
        
        function calculate_wrf_behr_correlation(varargin)
            % The idea of this function is to compute the correlation
            % between BEHR SCDs and WRF monthly and daily VCDs, to see if
            % the daily profiles provide better correlation than the
            % monthly profiles. This will (if it works) be a more
            % quantitative metric than the eyeballed agreement of
            % record_scd_comparison().
            E = JLLErrors;
            
            p = inputParser;
            p.addParameter('locations', []);
            p.addParameter('start', '');
            p.addParameter('end', '');
            
            p.parse(varargin{:});
            pout = p.Results;
            
            locations = pout.locations;
            start_date = pout.start;
            end_date = pout.end;
            
            if isempty(locations)
                locations = misc_behr_v3_validation.read_locs_file();
                loc_inds = ask_multiselect('Select which city(ies) to test', {locations.ShortName}, 'returnindex', true);
                locations = locations(loc_inds);
            elseif ~isstruct(locations) || any(~isfield('locations',{'ShortName', 'Latitude', 'Longitude', 'Radius'}))
                E.badinput('LOCATIONS must be the struct returned by misc_behr_v3_validation.read_locs_file, or a subset of it');
            end
            
            if isempty(start_date)
                start_date = datenum(ask_date('Enter the first date to test'));
            else
                start_date = validate_date(start_date);
            end
            
            if isempty(end_date)
                end_date = datenum(ask_date('Enter the last date to test'));
            else
                end_date = validate_date(end_date);
            end
            
            % Load each day/location's data and interpolate the BEHR pixels
            % to the WRF grid
            dvec = start_date:end_date;
            for a=1:numel(locations)
                for b=1:numel(dvec)
                    % Load the gridded BEHR data so that we can average
                    % over multiple orbits, if present
                    [behr_daily, wrf_monthly, wrf_daily] = misc_behr_v3_validation.load_wrf_and_behr_data(locations(a), dvec(b), 'load_gridded', true);
                    
                end
            end
            
        end
        
        
        
        function [is_significant, t] = calculate_slope_significant_difference(varargin)
            p = inputParser;
            p.addParameter('slope1_args', {});
            p.addParameter('slope2_args', {});
            p.addParameter('match_pandora_aircraft', nan);
            p.addParameter('remove_outliers', nan);
            p.addParameter('combine_air_pan', nan);
            p.addParameter('plot', true);
            
            p.parse(varargin{:});
            pout = p.Results;
            
            slope1_args = pout.slope1_args;
            slope2_args = pout.slope2_args;
            remove_outliers = pout.remove_outliers;
            do_plot = pout.plot;
            do_combine_aircraft_pandora = opt_ask_yn('Combine aircraft and Pandora observations?', pout.combine_air_pan, '"combine_air_pan"');
            if do_combine_aircraft_pandora
                match_pandora_aircraft = opt_ask_yn('Subset the data so that only sites/times with Pandora and aircraft data are used?', pout.match_pandora_aircraft, 'match_pandora_aircraft');
                load_fxn = @misc_behr_v3_validation.load_comparison_data;
            else
                load_fxn = @misc_behr_v3_validation.load_combined_aircraft_pandora;
            end
            
            
            fprintf('Loading data for slope 1...\n');
            [data_structs_1, opts1] = load_fxn(slope1_args{:}, 'match', match_pandora_aircraft);
            fprintf('Loading data for slope 2...\n');
            [data_structs_2, opts2] = load_fxn(slope2_args{:},  'match', match_pandora_aircraft);
            
            if numel(data_structs_1) > 1 || numel(data_structs_2) > 1
                E.notimplemented('comparing multiple slopes at once');
            end
            
            data1 = data_structs_1{1};
            data2 = data_structs_2{1};
            
            data_types = {opts1.data_source, opts2.data_source};
            if ~do_combine_aircraft_pandora && all(ismember({'pandora', 'aircraft'}, data_types))
                % If we have both a pandora and aircraft measurement, offer
                % the option to match them
                if opt_ask_yn('Subset the data so that only sites/times with Pandora and aircraft data are used?', match_pandora_aircraft, 'match_pandora_aircraft')
                    if strcmpi(opts1.data_source, 'aircraft')
                        [data1, data2] = misc_behr_v3_validation.match_aircraft_and_pandora_sites(data1, data2, 'match_time', true);
                    else
                        [data2, data1] = misc_behr_v3_validation.match_aircraft_and_pandora_sites(data2, data1, 'match_time', true);
                    end
                end
            end
            
            if opt_ask_yn('Remove outliers?', remove_outliers, 'remove_outliers')
                not_out1 = misc_behr_v3_validation.is_not_outlier(data1.x, data1.y);
                not_out2 = misc_behr_v3_validation.is_not_outlier(data2.x, data2.y);
            else
                not_out1 = true(size(data1.x));
                not_out2 = true(size(data2.x));
            end
            
            [line1x,line1y,line1_legstr,fit1_info] = calc_fit_line(data1.x(not_out1), data1.y(not_out1), 'regression', 'rma', 'xcoord', [-1e17 1e17]);
            [line2x,line2y,line2_legstr,fit2_info] = calc_fit_line(data2.x(not_out2), data2.y(not_out2), 'regression', 'rma', 'xcoord', [-1e17 1e17]);
            
            if do_plot
                figure;
                line(data1.x(not_out1), data1.y(not_out1), 'color', 'b', 'marker', 's', 'linestyle', 'none');
                l(1) = line(line1x, line1y, 'color', 'b');
                line(data2.x(not_out2), data2.y(not_out2), 'color', 'r', 'marker', 'x', 'linestyle', 'none');
                l(2) = line(line2x, line2y, 'color', 'r');
                legend(l', {line1_legstr, line2_legstr});
            end
            
            [is_significant, t] = slope_significant_difference( fit1_info.P, fit1_info.StdDevM, numel(data1.x), fit2_info.P, fit2_info.StdDevM, numel(data2.x) );
        end
        
        
    end
    
    methods(Static = true, Access = private)
        function leg_strings = versions2legend(versions)
            % Convert a cell array of version names into a cell array of
            % legend-formatted names
            leg_strings = cell(size(versions));
            for i=1:numel(versions)
                leg_strings{i} = misc_behr_v3_validation.plot_wrf_prof_legend_names.(versions{i});
            end
        end
        
        function scd_subgroup(h5_filename, group_name, results)
            E = JLLErrors;
            % Make sure group name does NOT end in a slash to avoid double
            % slashes in the name when we join it with the dataset
            if ~isscalar(results) || ~isstruct(results)
                E.badinput('RESULTS must be a scalar struct')
            end
            group_name = regexprep(group_name, '/$', '');
            deferred_fields = {};
            fns = fieldnames(results);
            n_fields = numel(fns);
            for i_fn = 1:n_fields
                this_field = fns{i_fn};
                this_data = results.(this_field);
                dset_name = strjoin({group_name, this_field}, '/');
                if isnumeric(this_data) || iscell(this_data)
                    if iscell(this_data)
                        this_data = cat_cell_contents(this_data);
                    end
                    sz = size(this_data);
                    h5create(h5_filename, dset_name, sz);
                    h5write(h5_filename, dset_name, this_data);
                elseif ischar(this_data)
                    deferred_fields{end+1} = this_field; %#ok<AGROW>
                elseif isstruct(this_data)
                    if ~isstruct(this_data)
                        E.notimplemented('Non-scalar substructure')
                    end
                    misc_behr_v3_validation.scd_subgroup(h5_filename, dset_name, this_data);
                end
            end
            
            for i_fn = 1:numel(deferred_fields)
                this_field = deferred_fields{i_fn};
                this_attr = results.(this_field);
                h5writeatt(h5_filename, group_name, this_field, this_attr);
            end
        end
        
        function [struct_type, plot_vars, labels] = comp_struct_type(comp_struct)
            if isstruct(comp_struct)
                % Dig down through the structure until we find the time range
                % field names, which are e.g. t1200_1500
                fns = fieldnames(comp_struct);
                while ~regcmp(fns{1}, 't\d{4}_\d{4}')
                    comp_struct = comp_struct.(fns{1});
                    fns = fieldnames(comp_struct);
                end
            elseif ischar(comp_struct)
                fns{1} = comp_struct;
            end
            
            % As of 4 May 2018, Pandora structs have time range 1230-1430
            % as their only range.
            if any(strcmpi(fns{1}, {'t1230_1430','pandora'}))
                struct_type = 'pandora';
                plot_vars = {'pandora_no2', 'sp_no2', 'behr_no2'};
                labels = struct('pandora_no2', 'Pandora NO_2 VCD (molec. cm^{-2})',...
                    'sp_no2', 'NASA NO_2 VCD (molec. cm^{-2})',...
                    'behr_no2', 'BEHR NO_2 VCD (molec. cm^{-2})');
            else
                struct_type = 'aircraft';
                plot_vars = {'air_no2_nasa', 'air_no2_behr', 'sp_no2', 'behr_no2'};
                labels = struct('air_no2_nasa', 'Aircraft NO_2 VCD (molec. cm^{-2})',...
                    'air_no2_behr', 'Aircraft NO_2 VCD (molec. cm^{-2})',...
                    'sp_no2', 'NASA NO_2 VCD (molec. cm^{-2})',...
                    'behr_no2', 'BEHR NO_2 VCD (molec. cm^{-2})');
            end
        end
        
        function campaigns = list_available_campaigns(comp_struct, prof_mode)
            prof_fn = sprintf('us_%s', prof_mode);
            substruct = comp_struct.v3.(prof_fn);
            campaigns = fieldnames(substruct);
            keep_campaigns = true(size(campaigns));
            for i_cam = 1:numel(campaigns)
                keep_campaigns(i_cam) = ~isempty(substruct.(campaigns{i_cam}).t1200_1500);
            end
            campaigns = campaigns(keep_campaigns);
        end
    end
    
    
end
