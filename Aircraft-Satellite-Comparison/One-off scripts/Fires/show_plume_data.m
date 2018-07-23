function [ output_args ] = show_plume_data( Fires, campaign_name, data_field, flight_date, avg_bool, linespec )
%SHOW_PLUME_DATA Function to highlight plumes in a plot to examine
%   This function is specifically designed to let me study the plumes
%   identified by find_campaign_fires. It will plot a time series of the
%   requested data field for the requested date, highlighting the time
%   periods that the Fires structure identified as fire plumes. The data
%   field can be one recognized by merge_field_names or in Merge.Data
%   itself.  If you want to plot multiple species on a subplot, pass them
%   as a cell array. The date should be a string.
%
%   avg_bool is optional (and defaults to false). If set to true, it will
%   average the data using the same averaging period as the Fires
%   structure.
%
%   linespec is also optional, and can be any combination of marker shape,
%   line style, and 1 letter color spec that will be passed along to the
%   timeseries (will be passed to plot).  This can also be passed as a cell
%   array if data_field is; that will let you set the linespec for each
%   subplot.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 June 2015

E = JLLErrors;

if ~isstruct(Fires) 
    E.badinput('Fires needs to be a structure output from find_campaign_fires')
end

if ~ischar(campaign_name)
    E.badinput('campaign_name should be a string')
end

if ~ischar(data_field) && (~iscell(data_field) || any(~iscellcontents(data_field,'ischar')))
    E.badinput('data_field should be a string or a cell array of strings')
end

if ~ischar(flight_date)
    E.badinput('flight_date should be a string')
end

% When parsing the optional arguments, there are four cases to handle:
%   1) Neither is passed and must be set to defaults
%   2) Only avg_bool is passed (linespec must be set to default)
%   3) Only linespec is passed: it will be read in as avg_bool and so the
%   arguments will need to be flipped and avg_bool set to default
%   4) Both are passed in
if nargin < 5
    avg_bool = false;
    linespec = 'b-';
elseif nargin == 5
    if ischar(avg_bool) || iscell(avg_bool) % only linespec passed
        linespec = avg_bool;
        avg_bool = false;
    elseif ~isscalar(avg_bool) || (~isnumeric(avg_bool) && ~islogical(avg_bool))
        E.badinput('avg_bool must be scalar and convert to a logical value. If you were trying to pass linespec, it must be a string or cell array of strings to be identified as such.')
    else
        linespec = 'b-';
    end
elseif narargin > 5;
    if ~isscalar(avg_bool) || (~isnumeric(avg_bool) && ~islogical(avg_bool))
        E.badinput('avg_bool must be scalar and convert to a logical value')
    elseif ~ischar(linespec) && (~iscell(linespec) || any(~iscellcontents(linespec, 'ischar')))
        E.badinput('linespec must be a string or cell array of strings');
    end
end

if iscell(linespec) && ~iscell(data_field)
    E.badinput('linespec can only be a cell array if data_field is')
elseif iscell(linespec) && numel(linespec) ~= numel(data_field)
    E.badinput('linespec and data_field must have the same number of elements if both are cell arrays');
end

[Names, dates, merge_dir] = merge_field_names(campaign_name);

% Check the dates
fire_datestrings = glob(fieldnames(Fires), 'Flight*');
fire_dates = nan(size(fire_datestrings));
for a=1:numel(fire_datestrings)
    f = regexp(fire_datestrings{a},'\d\d\d\d\d\d\d\d');    
    fire_dates(a) = datenum(fire_datestrings{a}(f:f+7),'yyyymmdd');
end

flight_datenum = datenum(flight_date);

if any(fire_dates < datenum(dates{1})) || any(fire_dates > datenum(dates{2}))
    E.callError('date_mismatch', 'Dates in the fire structure don''t match dates expected for the campaign %s', campaign_name);
elseif flight_datenum < datenum(dates{1}) || flight_datenum > datenum(dates{2})
    E.callError('date_mismatch', 'Date given is outside the range of dates for the campaign')
elseif ~ismember(flight_datenum, fire_dates)
    E.callError('date_mismatch', 'Date given doesn''t correspond to a day with fire data')
end

% Load the Merge file
file_datestr = sprintf('*%d_%02d_%02d.mat',year(flight_datenum), month(flight_datenum), day(flight_datenum));
F = dir(fullfile(merge_dir, file_datestr));
if numel(F) == 1
    M = load(fullfile(merge_dir, F(1).name)); % loads the Merge variable
    Merge = M.Merge;
    clear('M');
elseif numel(F) < 1
    E.filenotfound(sprintf('Merge file fitting pattern %s', file_datestr));
else
    E.toomanyfiles(sprintf('Merge file fitting pattern %s', file_datestr));
end



% Time series plot with highlighting
figure;
if ~iscell(data_field)
    data_field = {data_field};
else
    [nx,ny] = design_subplot(numel(data_field));
    subplot(ny, nx, 1);
end

% Correctly convert linespec into a cell array that matches data_field. It
% should only replicate if it is not a cell array already, but it always
% needs to be a cell array to be passed correctly to the plot drawing
% subfunction.
if ~iscell(linespec)
    linespec = {linespec};
    linespec = repmat(linespec, size(data_field));
end

for a=1:numel(data_field)
    
    % Get the requested data & UTC (it is a time series after all)
    if isfield(Names, data_field{a})
        [data,utc] = remove_merge_fills(Merge, Names.(data_field{a}));
        data_unit = Merge.Data.(Names.(data_field{a})).Unit;
    else
        [data,utc] = remove_merge_fills(Merge, data_field{a});
        data_unit = Merge.Data.(data_field{a}).Unit;
    end
    
    
    bckgnd = calc_day_background(Merge,data_field{a},campaign_name);
    
    if avg_bool
        data = avg_n_elements(data, Fires.run_params.n_sec_avg, 'op', 'nanmean');
        utc = avg_n_elements(utc, Fires.run_params.n_sec_avg, 'op', 'nanmean');
    end
    
    % a will only be greater than 1 if multiple data fields were passed,
    % this prevents us from creating a subplot when we don't need to.  The
    % first subplot is created outside the for loop.
    if a>1; subplot(ny,nx,a); end
    
    draw_plot(utc, data, data_field{a}, linespec{a});
end

    function draw_plot(utc, data, data_fieldname, linespec)
        plot(utc, data, linespec);
        line([min(utc), max(utc)], [bckgnd, bckgnd], 'linewidth', 3, 'color', 'r');
        fire_fieldname = sprintf('Flight_%s',datestr(flight_datenum,'yyyymmdd'));
        highlight_plot(Fires.(fire_fieldname).utc_ranges);
        set(gca,'fontsize',16);
        xlabel('UTC');
        ylabel(sprintf('%s (%s)', data_fieldname, data_unit));
        if ~avg_bool
            title(sprintf('Timeseries of %s on %s\nwith "fires" highlighed',regexprep(data_fieldname,'_',' '),flight_date));
        else
            title(sprintf('Timeseries of %s on %s\nwith "fires" highlighed\nAveraged in %d sec periods',data_fieldname,flight_date,Fires.run_params.n_sec_avg));
        end
    end

    

end

function [nx, ny] = design_subplot(n_plots)
ny = round(sqrt(n_plots));
nx = ceil(n_plots / ny);
end
