function [] = behr_validation_plots()
%BEHR_VALIDATION_PLOTS Generate plots and tables for the BEHR validation paper

mydir = fileparts(mfilename('fullpath'));
images_dir = fullfile(mydir,'..','Images');
paper_tex_file = fullfile(mydir, '..', 'BEHR3-validation.tex');
supp_tex_file = fullfile(mydir, '..', 'Supplement', 'behr3-validation-supplement.tex');

trop_test=false;

discover_campaigns = {'discover_md','discover_ca','discover_tx','discover_co'};
remove_neg_vcds_discover = false(size(discover_campaigns));
remove_neg_vcds_discover(strcmp(discover_campaigns, 'discover_co')) = true;

ltng_validation_campaigns = {'soas', 'seac4rs'};
remove_neg_vcds_ltng = false(size(ltng_validation_campaigns));

all_vcd_val_campaigns = veccat(discover_campaigns, ltng_validation_campaigns);
remove_neg_vcds_all = veccat(remove_neg_vcds_discover, remove_neg_vcds_ltng);

profile_campaigns = [discover_campaigns, {'dc3', 'seac4rs'}];

mat2latex_slope_sigma = {'u10', 2};
mat2latex_intR2 = {'%.3g'};

if ask_yn('Regenerate intermediate data? This will overwrite any existing files.')
    generate_comparison_data();
end

plot_surfp_diff()

senex_seac4rs_caption = 'Slopes and 1$\sigma$ uncertainties for RMA regression of satellite VCDs against in situ calculated VCDs. Both methods of extending the profiles (using GEOS-Chem modeled profiles or extrapolating the top/bottom ten points) are included. Outliers are removed before calculating these parameters.';
insert_vcd_insitu_table(ltng_validation_campaigns, false, 'aircraft',  senex_seac4rs_caption, 'tab:vcd-insitu-stats', paper_tex_file, 'SENEX-SEAC4RS-TABLE', 'fit_vals', 'std. dev.');
 
aircraft_caption = 'Slopes, intercepts, and $R^2$ values for RMA regression of satellite VCDs against in situ calculated VCDs. Outliers are removed before calculating these parameters; negative VCDs are retained unless noted.';
insert_vcd_insitu_table(all_vcd_val_campaigns, remove_neg_vcds_all, 'aircraft',  aircraft_caption, 'tab:supp:vcd-insitu-stats', supp_tex_file, 'VCDTABLE', 'env', 'sidewaystable', 'center', true, 'fit_vals', 'int+R2');
pandora_caption = 'Slopes, intercepts, and $R^2$ values for RMA regression of satellite VCDs against Pandora VCDs. Outliers are removed before calculating these parameters.';
insert_vcd_insitu_table(discover_campaigns, remove_neg_vcds_discover, 'pandora', pandora_caption, 'tab:supp:vcd-pandora-stats', supp_tex_file, 'PANDORATABLE', 'center', true, 'fit_vals', 'int+R2');
 
main_avg_caption = 'Slopes and 1$\sigma$ uncertainties of BEHR vs. combined aircraft (extended with GEOS-Chem profiles) and Pandora VCDs. Matched slopes use only Pandora data approximately coincident with aircraft profiles to get similar sampling; all uses all valid Pandora data. Outliers and negative VCDs are removed before computing slopes.';
make_avg_table(paper_tex_file, 'tab:vcd-avg-slopes', main_avg_caption, 'AVGTABLE', 'fit_vals', 'std. dev.');
make_avg_table(supp_tex_file, 'tab:supp:vcd-avg-slopes-intR2', main_avg_caption, 'AVGTABLE-SUPP', 'fit_vals', 'int+R2', 'center', true, 'env', sidewaystable);
 
make_individual_scatter_plots();
plot_scd_comparisons()
wrf_emis_change();
profile_mean_bias_table();
plot_average_wrf_profiles();
prof_r2_plots()
make_special_lightning_scatter();
dc3_vs_soas_flights();
uncertainty_plot();

%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting subfunctions %
%%%%%%%%%%%%%%%%%%%%%%%%%
    function profile_mean_bias_table()
        campaigns = get_available_campaigns_wrf_profs(profile_campaigns);
        xx_seac4rs = strcmp(campaigns, 'seac4rs');
        no_seac4rs = campaigns(~xx_seac4rs);
        
        bl_pres = [1100 775];
        ut_pres = [775 0];
        mab_bl = misc_behr_v3_validation.tabulate_profile_mab(campaigns, 'bias_output', 'total', 'pres_range', bl_pres);
        mab_bl_nosec = misc_behr_v3_validation.tabulate_profile_mab(no_seac4rs, 'bias_output', 'total', 'pres_range', bl_pres);
        mab_ut = misc_behr_v3_validation.tabulate_profile_mab(campaigns, 'bias_output', 'total', 'pres_range', ut_pres);
        mab_ut_nosec = misc_behr_v3_validation.tabulate_profile_mab(no_seac4rs, 'bias_output', 'total', 'pres_range', ut_pres);
        
        A = [struct2array(mab_bl)', struct2array(mab_bl_nosec)', struct2array(mab_ut)', struct2array(mab_ut_nosec)'];
        fprintf('Here''s the table--will just put in paper manually:\n');
        T = array2table(A, 'RowNames', fieldnames(mab_bl), 'VariableNames', {'BL','BL_NoSEAC4RS','UT','UT_NoSEAC4RS'})
        
        
    end

    function wrf_emis_change()
        fig = misc_behr_v3_validation.behr_v2_vs_v3_emis();
        save_all_the_formats(fig, 'v2-v3-no-emis-change', true);
    end

    function plot_average_wrf_profiles()
       [campaigns, has_daily] = get_available_campaigns_wrf_profs(profile_campaigns);
       figs = gobjects(numel(campaigns),1);
       for i_campaign = 1:numel(campaigns)
           if has_daily(i_campaign)
               prof_types = {'monthly', 'v2','daily'};
           else
               prof_types = {'monthly', 'v2'};
           end
           [figs(i_campaign), profiles] = misc_behr_v3_validation.plot_one_wrf_comparison(prof_types, campaigns{i_campaign}, 'All', 'std', 'bin_mode', 'mean');
           %title(strrep(upper(campaigns{i_campaign}),'_',' '));
           title('');
           ax = gca;
           ax.XLim(1) = 0; % sometimes a negative xlimit occurs, avoid that
           % Some campaigns have extreme values in the v2.1 profile due to
           % extrapolation, restrict their upper limits to ensure the other
           % profiles are visible
           switch lower(campaigns{i_campaign})
               case 'soas'
                   ax.XLim(2) = 2000;
               case 'seac4rs'
                   ax.XLim(2) = 2500;
           end
           
           ax.YLim(1) = 200;
           ax.YLim(2) = min(ax.YLim(2), 1020);
       end
       
       combo_fig = combine_plots(figs, 'dims', [0, 2], 'scale', 1);
       
       label_subfigs(combo_fig, 'xshift', 0.2);
       if mod(numel(figs),2) ~= 0
            % if we had an odd number of figures, then center the bottom
            % one. Assume that the children are in last-in-first-out order.
            ch = findall(combo_fig.Children, 'type', 'axes');
            center_axes(ch(1), ch(3), ch(2));
       end
       close(figs);
       %move_profile_legend(combo_fig);
       
       save_all_the_formats(combo_fig, 'WRF-Prof-Comparison');
    end

    function insert_vcd_insitu_table(campaigns, remove_neg_vcds, data_source, tex_caption, tex_label, tex_file, tex_marker, varargin)
        % TODO: include daily with "N/A" if not present, reverse order of
        % products (v3.0D, v3.0M, v2.1C, SP)
        p = advInputParser;
        p.addParameter('env', 'table');
        p.addParameter('center', false);
        p.addParameter('fit_vals', 'std. dev.');
        p.parse(varargin{:});
        pout = p.Results;
        
        environment = pout.env;
        center = pout.center;
        fit_vals = pout.fit_vals;
        
        using_stddev = strcmp(fit_vals, 'std. dev.');
        
        if using_stddev
            m2l = mat2latex_slope_sigma;
        else
            m2l = mat2latex_intR2;
        end
        
        common_opts = {'campaigns', campaigns, 'data_source', data_source, 'remove_neg_vcd', remove_neg_vcds, 'remove_outliers', true, 'force_origin', false, 'swap_md', false, 'return_as_table', false, 'fit_quantities', fit_vals};
        [values, colnames, rownames, extra_hlines] = misc_behr_v3_validation.tabulate_insitu_comparisons(common_opts{:}, 'extend_method', 'geos');
        [values, rownames, extra_hlines] = remove_spv2(values, rownames, extra_hlines);
        if strcmpi(data_source, 'aircraft')
            
            [values_ex, colnames_ex, rownames_ex] = misc_behr_v3_validation.tabulate_insitu_comparisons(common_opts{:}, 'extend_method', 'extrap');
            values_ex = remove_spv2(values_ex, rownames_ex);
            values = cat(2, values, values_ex);
            
            if using_stddev
                % Using std. deviation, the only quantity is slope +/-
                % sigma, so we can directly insert the extension method
                colnames = strrep(colnames, 'Slope', 'Slope (GEOS-Chem)');
                colnames_ex = strrep(colnames_ex, 'Slope', 'Slope (Extrap.)');
                ex_method = [];
            else
                % Two of the columns names are for the row headers.
                ex_method = [repmat({''}, 1, size(rownames,2)), repmat({'Extended with GEOS-Chem'}, 1, numel(colnames)-2), repmat({'Extended by extrapolation'}, 1, numel(colnames_ex)-2)];
            end
            
            colnames_all = veccat(colnames, colnames_ex(size(rownames,2)+1:end));
            % If we add in the second extend method, we need to overwrite
            % the column names
            colnames = cat(1, ex_method, colnames_all);
        end
        
        make_latex_table(values, 'colnames', colnames, 'rownames', rownames, 'extra_hlines', extra_hlines, 'insert', true, 'marker', tex_marker, 'file', tex_file,...
           'caption', tex_caption, 'label', tex_label, 'lines', {'\tophline','\middlehline','\bottomhline'}, 'm2l', m2l, 'environment', environment, 'center', center);
    end

    function make_avg_table(tex_file, tex_label, tex_caption, insert_marker, varargin)
        % TODO: include daily with "N/A" if not present, reverse order of
        % products (v3.0D, v3.0M, v2.1C, SP)
        p = advInputParser;
        p.addParameter('env', 'table');
        p.addParameter('center', false);
        p.addParameter('fit_vals', 'std. dev.');
        p.parse(varargin{:});
        pout = p.Results;
        
        environment = pout.env;
        center = pout.center;
        fit_vals = pout.fit_vals;

        if strcmp(fit_vals, 'std. dev.')
            m2l = mat2latex_slope_sigma;
        else
            m2l = mat2latex_intR2;
        end
        
        [values, colnames, rownames, ~, extra_hlines] = misc_behr_v3_validation.make_aircraft_pandora_combo('return_as_table', false, 'fit_quantities', fit_vals, 'products', {'SP v3.0', 'BEHR v2.1C', 'BEHR v3.0B (M)', 'BEHR v3.0B (D)'});

        make_latex_table(values, 'colnames', colnames, 'rownames', rownames, 'insert', true, 'marker', insert_marker, 'file', tex_file,...
           'caption', tex_caption, 'label', tex_label, 'lines', {'\tophline','\middlehline','\bottomhline'}, 'extra_hlines', extra_hlines,  'm2l', m2l, 'environment', environment, 'center', center);
    end

    function plot_surfp_diff()
        data_opts = {'data_source','aircraft','extend_method','geos','x_var','air_no2_behr','y_var','behr_no2',...
            'prof_mode','daily','campaigns',{'discover_co'},'time_range','t1200_1500','version','v3'};
        scale_ht_dir = fullfile(misc_behr_v3_validation.validation_root_dir, 'VCD-Comparison','WRF-Tropopause');
        [fit_data_hypso, ~, ~, ~, x_hypso, y_hypso] = misc_behr_v3_validation.calculate_fit('load-combined', [data_opts, 'match', true, 'alt_dir', {''}], 'remove_outliers', true, 'remove_neg_sat', true, 'force_origin', false);
        [fit_data_scaleht, ~, ~, ~, x_scaleht, y_scaleht] = misc_behr_v3_validation.calculate_fit('load-combined', [data_opts, 'match', true, 'alt_dir', scale_ht_dir], 'remove_outliers', true, 'remove_neg_sat', true, 'force_origin', false);
        
        lims = [0 2e16];
        fig = figure;
        l = gobjects(4,1);
        l(1) = line(x_scaleht{1}, y_scaleht{1}, 'marker', 'o', 'markersize', 10, 'linestyle', 'none', 'linewidth', 1, 'color', 'b');
        l(2) = line(lims, fit_data_scaleht.P(1) .* lims + fit_data_scaleht.P(2), 'linestyle', '--', 'linewidth', 2, 'color', 'b');
        l(3) = line(x_hypso{1}, y_hypso{1}, 'marker', '^', 'markersize', 10, 'linestyle', 'none', 'linewidth', 1, 'color', 'r');
        l(4) = line(lims, fit_data_hypso.P(1) .* lims + fit_data_hypso.P(2), 'linestyle', '--', 'linewidth', 2, 'color', 'r');
        line(lims, lims, 'linestyle', ':', 'linewidth', 2, 'color', 'k');
        
        xylims(lims);
        scaleht_str = sprintf('Fit (%.2fx + %.2g)', fit_data_scaleht.P(1), fit_data_scaleht.P(2));
        hypso_str = sprintf('Fit (%.2fx + %.2g)', fit_data_hypso.P(1), fit_data_hypso.P(2));
        legend(l, {'Scale height', scaleht_str, 'Hypsometric eqn.', hypso_str})
        xlabel('Aircraft/Pandora NO_2 VCDs (molec. cm^{-2})');
        ylabel('BEHR v3.0 (D) NO_2 VCDs (molec. cm^{-2})');
        set(gca, 'fontsize', 12, 'xtick', lims(1):5e15:lims(2), 'ytick', lims(1):5e15:lims(2));
        reset_fig_size(fig);
        
        save_all_the_formats(fig, 'Hypsometric-SurfP-Comparison', true);
    end
    
    function plot_scd_comparisons()
        results_file = fullfile(misc_behr_v3_validation.scd_comp_dir, 'scd_daily_2018-05-11_15-17-00.mat');
        R = load(results_file);
        
        [~,~,~,good_agreement_fig] = misc_behr_v3_validation.plot_scd_vs_wrf_columns(R.results(24), 'titles', false);
        [~,~,~,bad_agreement_fig] = misc_behr_v3_validation.plot_scd_vs_wrf_columns(R.results(5), 'titles', false);
        
        label_subfigs(good_agreement_fig, 'xshift', 0.2);
        save_all_the_formats(good_agreement_fig, 'SCD-Comparison-Good', false);
        label_subfigs(bad_agreement_fig, 'xshift', 0.2);
        save_all_the_formats(bad_agreement_fig, 'SCD-Comparison-Bad', false);
        
%         % Note that combine_plots reverses the order of plots if >1 is
%         % present in a figure being merged, so the daily WRF will be on the
%         % left and SCDs on the right.
%         combo_fig = combine_plots(figs, 'scale', 1);
%         label_subfigs(combo_fig, 'xshift', 0.2);
%         save_all_the_formats(combo_fig, 'SCD-Comparison', false);
    end

    function make_individual_scatter_plots()
        campaigns = all_vcd_val_campaigns;
        for i_cam = 1:numel(campaigns)
            subfigs = gobjects(0);
            common_air_opts = {'data_source', 'aircraft', 'prof_mode', 'both', 'version', 'both', 'campaigns', campaigns{i_cam}, 'time_range', 't1200_1500',...
                'remove_outliers', true, 'plot_type', 'scatter', 'color_by', 'none', 'remove_neg_sat', false, 'match_pandora_aircraft', false};
            common_pandora_opts = update_params(common_air_opts, 'data_source', 'pandora', 'time_range', 't1230_1430');
            subfigs(1) = misc_behr_v3_validation.plot_one_vcd_comparison('extend_method', 'geos', 'x_var', 'air_no2_nasa', 'y_var', 'sp_no2', common_air_opts{:});
            title('');
            subfigs(2) = misc_behr_v3_validation.plot_one_vcd_comparison('extend_method', 'geos', 'x_var', 'air_no2_behr', 'y_var', 'behr_no2', common_air_opts{:});
            title('');
            subfigs(3) = misc_behr_v3_validation.plot_one_vcd_comparison('extend_method', 'extrap', 'x_var', 'air_no2_nasa', 'y_var', 'sp_no2', common_air_opts{:});
            title('');
            subfigs(4) = misc_behr_v3_validation.plot_one_vcd_comparison('extend_method', 'extrap', 'x_var', 'air_no2_behr', 'y_var', 'behr_no2', common_air_opts{:});
            title('');
            
            % Pandora is only available for the DISCOVER campaigns
            if ismember(campaigns{i_cam}, discover_campaigns)
                subfigs(5) = misc_behr_v3_validation.plot_one_vcd_comparison('x_var', 'pandora_no2', 'y_var', 'sp_no2', common_pandora_opts{:});
                title('');
                subfigs(6) = misc_behr_v3_validation.plot_one_vcd_comparison('x_var', 'pandora_no2', 'y_var', 'behr_no2', common_pandora_opts{:});
                title('');
            end
            combo_fig = combine_plots(subfigs, 'dims', [0 2], 'scale', 1.25);
            close(subfigs);
            % Need to adjust the axes position a little bit to keep the
            % x-label on the plot
            ax = findobj(combo_fig,'type','axes');
            for i_ax = 1:numel(ax)
                ax(i_ax).Position([2,4]) = [0.02, -0.02] + ax(i_ax).Position([2,4]);
            end
            label_subfigs(combo_fig,'xshift',0.2);
            save_all_the_formats(combo_fig, sprintf('vcd-scatter-%s', strrep(campaigns{i_cam},'_','-')), true);
        end
    end

    function make_special_lightning_scatter()
        campaigns = {'soas','seac4rs'};
        fig = gobjects(3, numel(campaigns));
        for i_cam = 1:numel(campaigns)
            campaign = campaigns{i_cam};
            
            fig(1, i_cam) = misc_behr_v3_validation.plot_vcd_prof_ensembles('extend_methods', {'wrf','geos'}, 'prof_type', 'daily', 'campaigns', {campaign}, 'title', false);
            
            common_air_opts = {'data_source', 'aircraft', 'prof_mode', 'daily', 'version', 'v3', 'campaigns', campaign, 'time_range', 't1200_1500',...
                'remove_outliers', true, 'plot_type', 'scatter', 'color_by', 'none', 'remove_neg_sat', false, 'match_pandora_aircraft', false,...
                'x_var', 'air_no2_behr', 'y_var', 'behr_no2'};
            [geos_ex, geos_opts] = misc_behr_v3_validation.load_comparison_data('extend_method', 'geos', common_air_opts{:});
            wrf_ex = misc_behr_v3_validation.load_comparison_data('extend_method', 'wrf', common_air_opts{:});
            
            geos_ex = geos_ex{1};
            wrf_ex = wrf_ex{1};
            
            fig(2, i_cam) = figure;
            geos_x = veccat(geos_ex.x);
            wrf_x = veccat(wrf_ex.x);
            geos_y = veccat(geos_ex.y);
            wrf_y = veccat(wrf_ex.y);
            lon = [geos_ex.profile_lon];
            lat = [geos_ex.profile_lat];
            l = plot_changes([geos_x(:), wrf_x(:)], [geos_y(:), wrf_y(:)], 'group_fmts', struct('marker', {'^','*'}, 'color', {[0 0.5 0],'m'}));
            line([0 1e16], [0 1e16], 'color', 'r', 'linestyle', '--', 'linewidth', 2);
            
            set(gca,'fontsize', 14);
            xlabel(geos_opts.labels.x);
            ylabel(geos_opts.labels.y);
            legend(l, {'GEOS-Chem extension', 'WRF-Chem extension'},'location','best');
            xlim([0 1e16]);
            ylim([0 4e15]);
            
            %map_params = update_params(common_air_opts, 'plot_type', 'map');
            %fig(2) = misc_behr_v3_validation.plot_one_vcd_comparison('extend_method', 'wrf', map_params{:});
            %title('');
            % Move the scatter plot up to the top of the render stack
            %uistack(fig(2).Children(end).Children(end), 'top');
            fig(3,i_cam) = figure;
            delta_vcd = wrf_x - geos_x;
            scatter(lon, lat, 48, delta_vcd,'filled');
            state_outlines('k');
            cb = colorbar;
            cb.Label.String = '\Delta VCD (WRF - GEOS)';
            caxis(calc_plot_limits(delta_vcd,'diff'));
            set(gca,'fontsize',16);
        end
        
        combo_fig = combine_plots(fig', 'dims', [0 numel(campaigns)], 'scale', 1);
        label_subfigs(combo_fig, 'xshift', 0.2);
        colormap(blue_red_only_cmap);
        
        % Need to adjust the axes position a little bit to keep the
        % x-label on the plot
        ax = findobj(combo_fig,'type','axes');
        for i_ax = 1:numel(ax)
            ax(i_ax).Position([2,4]) = [0.02, -0.02] + ax(i_ax).Position([2,4]);
        end
        
        close(fig(:));
        
        save_all_the_formats(combo_fig, 'Lightning-GEOSvWRF-scatter', false);
    end

    function dc3_vs_soas_flights()
        wrf = load(misc_behr_v3_validation.wrf_comp_file);
        fig = figure;
        plot(wrf.daily.soas.All.match.data.lon, wrf.daily.soas.All.match.data.lat,...
            wrf.daily.dc3.All.match.data.lon, wrf.daily.dc3.All.match.data.lat,...
            wrf.daily.seac4rs.All.match.data.lon, wrf.daily.seac4rs.All.match.data.lat);
        state_outlines('k');
        set(gca,'fontsize',14);
        legend('SOAS','DC3','SEAC4RS');
        
        save_all_the_formats(fig, 'DC3-SOAS-Flightpaths', true);
    end

    function uncertainty_plot()
        fig = misc_behr_v3_validation.plot_uncertainty_estimate('change_field', 'PercentChangeNO2','normalize_by','None','titles',false);
        save_all_the_formats(fig, 'Uncertainty', true);
    end

    function prof_r2_plots
        campaigns = {'discover_ca','discover_tx','discover_co'};
        fig = misc_behr_v3_validation.plot_prof_regression('campaigns', campaigns, 'versions', {'monthly', 'daily'}, 'title', false, 'num_pts', false);
        save_all_the_formats(fig, 'Profile-R2', false);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generative subfunctions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function generate_comparison_data()
        if trop_test
            extra_opts = {'v3_dir', '/Volumes/share-sat/SAT/BEHR/BEHR_Files_SPandTrop_Test/'};
        else
            extra_opts = {};
        end

        extend_methods = {'geos','wrf','extrap'};
        for i_method = 1:numel(extend_methods)
            misc_behr_v3_validation.generate_vcd_comparison_structs('data_source', 'aircraft', 'extend_method', extend_methods{i_method}, 'do_save', true, 'overwrite', true, extra_opts{:});
        end
        misc_behr_v3_validation.generate_vcd_comparison_structs('data_source', 'pandora', 'do_save', true, 'overwrite', true, extra_opts{:});
        misc_behr_v3_validation.generate_profile_comparison_struct('do_save', true);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Helper nested functions %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function save_all_the_formats(hfig, filename, is_supplement)
        if ~exist('is_supplement','var')
            is_supplement = false;
        end
        mydir = fileparts(mfilename('fullpath'));
        if ~is_supplement
            images_dir = fullfile(mydir,'..','Images');
        else
            images_dir = fullfile(mydir, '..', 'Supplement', 'Images');
        end
        set(hfig, 'paperpositionmode', 'auto');
        print(hfig, '-depsc2', '-loose', fullfile(images_dir, filename));
        
        fig_name = [fullfile(images_dir,filename), '.fig'];
        savefig(hfig, fig_name);
        
        saveas(hfig, fullfile(images_dir, filename), 'png');
    end

    

end



% Helper subfunctions
function [campaigns, with_daily] = get_available_campaigns_wrf_profs(desired_campaigns)
% Return a cellstr array of campaign names available for WRF profile
% comparison and a logical array that is true for campaigns with daily
% profiles available.
S = load(misc_behr_v3_validation.wrf_comp_file);
campaigns = fieldnames(S.monthly);
campaigns(~ismember(campaigns, desired_campaigns)) = [];
with_daily = isfield(S.daily, campaigns);
end

function move_profile_legend(combo_fig)

% Find the legends in the figure children
is_legend = false(size(combo_fig.Children));
for i=1:numel(is_legend)
    is_legend(i) = strcmp(combo_fig.Children(i).Type, 'legend');
end

all_legends = combo_fig.Children(is_legend);
all_axes = combo_fig.Children(~is_legend);

if mod(numel(all_axes),2) ~= 1
    warning('Even number of plots - the legend will cover one of them')
end

% Move the first legend that includes "Daily" profiles; delete the rest
for i=1:numel(all_legends)
    if ismember('Daily', all_legends(i).String)
        the_legend = all_legends(i);
        all_legends(i) = [];
        delete(all_legends);
        break
    end
end

% Assume that the first axes are the bottom left and the second are in the
% row above but on the right (plotting is FILO in Matlab apparently). Use
% that to position the legend in the center of the empty bottom right
% space. Make the legend text bigger first so that the position reflects
% its final size.

the_legend.FontSize = 20;
legend_size = the_legend.Position(3:4);
vertical_center = all_axes(1).Position(2) + 0.5*all_axes(1).Position(4);
horiz_center = all_axes(2).Position(1) + 0.5*all_axes(2).Position(3);

legend_pos(1) = horiz_center - 0.5*legend_size(2);
legend_pos(2) = vertical_center - 0.5*legend_size(1);


the_legend.Position(1:2) = legend_pos;

end

function [vals, rows, varargout] = remove_spv2(vals, rows, hlines)
% Quick fix since my main tabulate function always includes all products
xx = ~strcmpi(rows(:,2), 'SP v2.1');
if nargin > 2
    % Must do this before cutting down the rows
    hlines_bool = false(size(rows,1),1);
    hlines_bool(hlines) = true;
    varargout{1} = hlines_bool(xx);
end
vals = vals(xx, :);
rows = rows(xx,:);
end
