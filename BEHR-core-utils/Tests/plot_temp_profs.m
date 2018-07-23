function [  ] = plot_temp_profs( month_in )
%PLOT_TEMP_PROFS Plot the NASA temperature profiles for one month
%   PLOT_TEMP_PROFS( MONTH_IN ) Given the numeric MONTH_IN, plots the NASA
%   temperature profiles using PLOT_SLICE_GUI() (general utils repo). This
%   was used to verify that the temperature profiles were being assigned
%   wrong in v2.1C.

E = JLLErrors;
Grid = GlobeGrid(5,2);
if ~exist('month_in', 'var')
    month_in = ask_number('Enter the month', 'testfxn', @(x) isscalar(x) && x >= 1 && x <= 12);
elseif ~isscalar(month_in) || ~isnumeric(month_in)
    E.badinput('MONTH_IN must be a scalar number');
end

month_in = repmat(month_in, size(Grid.GridLon));

fileTmp = fullfile(behr_paths.amf_tools_dir,'nmcTmpYr.txt');
temp = rNmcTmp2(fileTmp, behr_pres_levels, Grid.GridLon, Grid.GridLat, month_in);
temp = permute(temp, [2 3 1]);
M = load('coast');
plot_slice_gui(temp, Grid.GridLon, Grid.GridLat, M.long, M.lat);
end

