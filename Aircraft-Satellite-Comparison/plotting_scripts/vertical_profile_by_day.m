% vertical_profile_by_day
%	A script that will plot vertical profiles (1 or 2) using the functions
%	plot_vertical_profile and add_vertical_profile.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 May 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%	USER INPUT VAR	  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The location to plot.  Must match the file name location.
location = 'Texas';

% Data Type: Aircraft or Sondes
dtype = 'Aircraft';

% Start and end dates.  Leave as empty strings to plot all days from a location.
start_date = '';
end_date = '';

% Set to 1 to mark boundary layer heights on the plots
blh_bool = 0;

% Set to 1 to automatically save the figures
save_bool = 0;
save_dir = '/Users/Josh/Dropbox/Berkeley/Research/NO2 Profile Project/Preliminary Figs/NO2 vs Theta/';
save_format = 'png'; %Leave blank for fig

% Set to 1 to close each figure after generating it
close_bool = 0;

% Fields to plot vertical profiles for.  Leave field2 as an empty str if no second profile is desired.
field = 'NO2_LIF';
field2 = '';

% Bin width in kilometers.  An empty binwidth2 will be set equal to binwidth
binwidth1 = 0.3;
binwidth2 = [];

% Only plot profiles within a certain range of profile numbers or site flags.
% To specify a range, use [min max].  Leave empty (i.e. []) to disregard that parameter.
profnums = [101000,101999];
siteflags = [];

% The location where the .mat files from reading in .ict merge files can be found
mat_file_dir = '/Volumes/share/GROUP/DISCOVER-AQ/Matlab Files';

DEBUG_LEVEL = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%	INPUT PARSING	  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(location);
    error('vert_prof_by_day:location','Location cannot be an empty string');
end

if ~any(strcmp(dtype,{'Aircraft','Sondes'}));
    error('vert_prof_by_day:data_type','Data type is not valid.')
end

if ~exist(fullfile(mat_file_dir,dtype), 'dir');
    error('vert_prof_by_day:mat_files','.mat file directory not valid');
end

if ~exist(save_dir, 'dir') && ~isempty(save_dir);
    error('vert_prof_by_day:mat_files','Save file directory not valid');
end

if close_bool && ~save_bool
    error('vert_prof_by_day:closefig','Figures are being automatically closed without saving.');
end

% Set wide ranges for profnums or site flags if they are empty -
% this will include all values in the vertical profile
if isempty(profnums); profnums = [-1e10 1e10]; end
if isempty(siteflags); siteflags = [-1000 1000]; end

% Set binwidths equal to no value given for the second one.
if isempty(binwidth2); binwidth2 = binwidth1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%		LOAD FILES	  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_pattern = [location, '*.mat'];
mat_files = dir(fullfile(mat_file_dir,dtype, file_pattern));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%		MAIN LOOP	  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all
blh_vec = zeros(1,numel(mat_files));
for a=1:numel(mat_files);
    % Handle date ranges
    fname = mat_files(a).name;
    [s, e] = regexp(fname,'\d\d\d\d.*\d\d.*\d\d'); %Find the date in the file name, assuming it is of the form yyyymmdd with any or no single character separator
    fdate = datenum(fname(s:e));
    if ~isempty(start_date) && ~isempty(end_date) %If the user specified start and end dates, check that the file is within the date range
        % If the date of the file falls outside the range give, skip this file.
        if fdate < datenum(start_date) || fdate > datenum(end_date)
            if DEBUG_LEVEL > 0; fprintf('Skipping %s\n',fname(s:e)); end
            continue
        end
    end
    
    if DEBUG_LEVEL > 0; fprintf('Plotting %s\n',fname(s:e)); end
    load(fullfile(mat_file_dir,dtype, mat_files(a).name), 'Merge');
    
    if strcmp(dtype,'Aircraft')   
        [bin_vals, bin_alt] = plot_vertical_profile(Merge, field, 'binwidth',binwidth1,'profnum',profnums,'siteflag',siteflags);
    elseif strcmp(dtype,'Sondes')
        [bin_vals, bin_alt] = plot_sonde_vertical_profile(Merge, field, 'binwidth',binwidth1);
    end
    bl_height = find_bdy_layer_height(bin_vals, bin_alt, field);
    blh_vec(a) = bl_height; 
    % Get the current line for formatting, and legend
    l1 = findall(gca,'Type','Line');
    
    if blh_bool; line([0 max(bin_vals)],[bl_height, bl_height],'linewidth',2,'linestyle','--','color','k'); end
    
    if ~isempty(field2)
        if strcmp(dtype,'Aircraft') 
            add_vertical_profile(Merge, field2, 'binwidth',binwidth2,'profnum',profnums','siteflag',siteflags);
        elseif strcmp(dtype,'Sondes')
            add_sonde_vertical_profile(Merge, field2, 'binwidth', binwidth2);
        end
        l2 = findall(gca,'Type','Line');
        legend([l1;l2],{field,field2},'Interpreter','none');
    end
    
    
    
    if save_bool
        if ~isempty(field2); savename = sprintf('%s_vertical_profile_%s',field,fname(s:e));
        else savename = sprintf('%s_and_%s_vertical_profiels_%s',field,field2,fname(s:e));
        end
        
        if strcmp(save_format,'fig') || isempty(save_format);
            savefig(savename);
        else
            saveas(gcf,savename,save_format)
        end
    end
    
    if close_bool;
        close all
    end
    
end