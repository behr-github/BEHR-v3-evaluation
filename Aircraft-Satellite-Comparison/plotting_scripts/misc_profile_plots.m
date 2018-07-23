function [ varargout ] = misc_profile_plots( plttype, varargin )
%MISC_PROFILE_PLOTS A collection of plotting function needed infrequently
%   This function will contain various miscellaneous plotting functions for
%   working with aircraft profiles. Collecting them here allows me to reuse
%   them without dozens of .m files taking up space in my directories.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Jul 2015
E = JLLErrors;
plttype = lower(plttype);
switch plttype
    case 'range_shape'
        plot_range_shape(varargin{1:3});
    case 'range_overlap'
        plot_range_overlap(varargin{:});
    case 'range_prof'
        plot_profile_no2_aerosol(varargin{1:3});
    case 'cat_comparisons'
        varargout{1} = concatenate_comparisons(varargin{:});
    case 'dist_curtain'
        plot_distance_curtain(varargin{:});
    otherwise
        fprintf('Could not recognize plot type\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OTHER FUNCTIONS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
    function O = concatenate_comparisons(varargin)
        fns = fieldnames(varargin{1});
        O(numel(varargin{1})) = make_empty_struct_from_cell(fns);
        xx = ~iscellcontents(regexp(fns,'db_iall'),'isempty');
        fns(xx) = []; % db_iall needs handled separately
        fns2 = fieldnames(varargin{1}(1).db_iall);
        for s=1:numel(varargin{1})
            O(s).db_iall = make_empty_struct_from_cell(fns2);
        end
        for a=1:numel(varargin)
            for s=1:numel(varargin{1})
                for b=1:numel(fns)
                    O(s).(fns{b}) = cat(1, O(s).(fns{b}), varargin{a}(s).(fns{b}));
                end
                for b=1:numel(fns2)
                    O(s).db_iall.(fns2{b}) = cat(2, O(s).db_iall.(fns2{b}), varargin{a}(s).db_iall.(fns2{b}));
                end
            end
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PLOTTING FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_range_shape(campaign_name, date, utcrange)
        % Given the name of a campaign, a date, and a 2-element vector with
        % the starting and ending times for the range, this will plot the
        % altitude and lat/lon shape of that part of the flight.
        [Names, ~, directory] = merge_field_names(campaign_name);
        file_datestr = datestr(date, 'yyyy_mm_dd');
        Files = dir(fullfile(directory, '*.mat'));
        file_names = {Files(:).name};
        find_file = regexp(file_names, file_datestr);
        xx = ~iscellcontents(find_file,'isempty');
        if sum(xx) > 1
            E.toomanyfiles(sprintf('%s %s merge', file_datestr, campaign_name));
        elseif sum(xx) < 1
            E.filenotfound(sprintf('%s %s merge', file_datestr, campaign_name));
        else
            load(fullfile(directory, file_names{xx}),'Merge'); % will add the Merge file to the workspace
        end
        
        utc = remove_merge_fills(Merge, 'UTC');
        lon = remove_merge_fills(Merge, Names.longitude);
        lat = remove_merge_fills(Merge, Names.latitude);
        alt = remove_merge_fills(Merge, Names.gps_alt);
        
        figure;
        subplot(1,2,1);
        
        tt = utc >= utcrange(1) & utc < utcrange(2);
        line(utc(tt), alt(tt));
        xlabel('UTC (s)');
        ylabel(sprintf('Alt (%s)', Merge.Data.(Names.gps_alt).Unit));
        
        subplot(1,2,2);
        scatter(lon(tt), lat(tt), 8, alt(tt));
        cb = colorbar;
        cb.Label.String = sprintf('Alt (%s)', Merge.Data.(Names.gps_alt).Unit);
    end

    function [  ] = plot_range_overlap( varargin )
        %PLOT_RANGE_OVERLAP Plots the UTC ranges in any number of Range structures
        %   A Range structure contains many lists of UTC time ranges that are used
        %   to identify portions of an aircraft flight that have interesting data.
        %   Often, I will set up multiple Range structures for a single campaign
        %   identifying different flight patterns, e.g. profiles vs. wall vs.
        %   porpoising. This will plot any number of range files from the same
        %   campaign. It's original intention is to be able to check that I've not
        %   tagged a particular maneuver as two different flight patterns, but I'm
        %   sure there will be other uses.
        %
        %   Josh Laughner <joshlaugh5@gmail.com> 10 Jul 2015
        
        E = JLLErrors;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% INPUT CHECKS %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % Check that the only thing input are Range type data structures and that
        % they are all the same.
        
        if any(~iscellcontents(varargin,'isstruct'))
            E.badinput('All inputs must be Range structures')
        end
        
        for a=1:numel(varargin)
            if ~isfield(varargin{a},'Date') || ~isfield(varargin{a},'Ranges')
                E.badinput('Input %d does not have the expected fields ''Date'' and ''Ranges''', a)
            end
            
            if a==1
                n = numel(varargin{a});
                dates = datenum({varargin{a}(:).Date});
            else
                if numel(varargin{a}) ~= n
                    E.badinput('Input %d has a different number of days than the first input', a)
                end
                newdates = datenum({varargin{a}(:).Date});
                if any(dates ~= newdates)
                    E.badinput('Input %d has different dates than the first input', a)
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%% PLOT MAKING %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Calculate the vertical offset between different Range structures
        if numel(varargin) > 1
            offsets = linspace(-0.2,0.2,numel(varargin));
        else
            offsets = 0;
        end
        
        % Use this to plot each Range set as a different color
        colord = get(groot,'defaultAxesColorOrder');
        markerlist = {'o','x','^','v','s','h','d'};
        
        figure;
        for a=1:numel(varargin)
            mark = markerlist{mod(a,7)};
            col = colord(mod(a,7),:);
            for b=1:n
                y = datenum(varargin{a}(b).Date) + offsets(a);
                for c=1:size(varargin{a}(b).Ranges,1)
                    x = varargin{a}(b).Ranges(c,:);
                    line(x, [y y], 'color', col, 'marker', mark);
                end
            end
        end
        
        set(gca, 'ytick',min(dates):max(dates));
        set(gca, 'yticklabels', cellstr(datestr(min(dates):max(dates),29)));
        set(gca, 'ygrid', 'on');
        
    end

    function [  ] = plot_profile_no2_aerosol(campaign_name, date, utcrange)
        % Given the name of a campaign, a date, and a 2-element vector with
        % the starting and ending times for the range, this will plot the
        % NO2 and aerosol extinction profiles for that part of the flight.
        [Names, ~, directory] = merge_field_names(campaign_name);
        file_datestr = datestr(date, 'yyyy_mm_dd');
        Files = dir(fullfile(directory, '*.mat'));
        file_names = {Files(:).name};
        find_file = regexp(file_names, file_datestr);
        xx = ~iscellcontents(find_file,'isempty');
        if sum(xx) > 1
            E.toomanyfiles(sprintf('%s %s merge', file_datestr, campaign_name));
        elseif sum(xx) < 1
            E.filenotfound(sprintf('%s %s merge', file_datestr, campaign_name));
        else
            load(fullfile(directory, file_names{xx}),'Merge'); % will add the Merge file to the workspace
        end
        
        utc = remove_merge_fills(Merge, 'UTC');
        no2 = remove_merge_fills(Merge, Names.no2_lif);
        aer = remove_merge_fills(Merge, Names.aerosol_extinction);
        alt = remove_merge_fills(Merge, Names.radar_alt);
        
        tt = utc >= utcrange(1) & utc < utcrange(2);
        [no2bins, no2binalt] = bin_vertical_profile(alt(tt), no2(tt), 0.25);
        [aerbins, aeraltbins] = bin_vertical_profile(alt(tt), aer(tt), 0.25);
        
        figure;
        ax = plotxx(no2bins, no2binalt, aerbins, aeraltbins, {'[NO_2]', 'Aerosol ext.'}, {'Alt (km)', 'Alt (km)'});
        
        % make sure the y limits are the same for both
        y1 = get(ax(1),'ylim');
        y2 = get(ax(2),'ylim');
        ymax = max([y1, y2]);
        set(ax(1),'ylim',[0 ymax]);
        set(ax(2),'ylim',[0 ymax]);
    end

    function [  ] = plot_distance_curtain(data, lon, lat, alt, center_lon, center_lat)
        if ndims(data) ~= ndims(lon) || any(size(data) ~= size(lon))
            E.badinput('data and lon must be the same size')
        elseif ndims(data) ~= ndims(lat) || any(size(data) ~= size(lat))
            E.badinput('data and lat must be the same size')
        elseif ndims(data) ~= ndims(alt) || any(size(data) ~= size(alt))
            E.badinput('data and alt must be the same size')
        end
        
        if ~isscalar(center_lon) || ~isnumeric(center_lon)
            E.badinput('center_lon must be scalar numeric');
        elseif ~isscalar(center_lat) || ~isnumeric(center_lat)
            E.badinput('center_lat must be scalar numeric')
        end
        
        dist = ( (lon - center_lon).^2 + (lat - center_lat).^2 ).^0.5;
        dist_lims = [min(min(dist(:)),0), max(max(dist(:)),0)];
        
        alt_lims = [min(alt(:)), max(alt(:))];
        
        % Ensure the max and min will always be a bin limit, but try to
        % divide distance into bins of ~0.05 deg and altitude into bins of
        % ~0.25 km.
        dist_bins = linspace(dist_lims(1), dist_lims(2), ceil((dist_lims(2) - dist_lims(1))/0.05));
        alt_bins = linspace(alt_lims(1), alt_lims(2), ceil((alt_lims(2) - alt_lims(1))/0.25));
        
        bins = nan(numel(alt_bins)-1, numel(dist_bins)-1);
        
        for a=1:(numel(alt_bins)-1)
            for b=1:(numel(dist_bins)-1)
                xx = dist >= dist_bins(b) & dist < dist_bins(b+1) & alt >= alt_bins(a) & alt < alt_bins(a+1);
                bins(a,b) = nanmean(data(xx));
            end
        end
        
        [X,Y] = meshgrid(dist_bins(1:end-1), alt_bins(1:end-1));
        
        figure;
        pcolor(X,Y,bins);
        xlabel('Distance (degrees)')
        ylabel('Altitude (km)')
        colorbar;
    end
end

