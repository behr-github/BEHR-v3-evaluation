function [ Matched_Datas ] = run_gcas_verification( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

noinput = 'noinput';

E = JLLErrors;
p = inputParser;
p.addParameter('campaign', noinput);
p.addParameter('behr_dir', noinput);
p.addParameter('vectorize', noinput);

% Inputs passed onto gcas_sao_verification
p.addParameter('cloud_prod', noinput);
p.addParameter('cloud_frac_max', noinput);
p.addParameter('row_anomaly', noinput);
p.addParameter('sat_fields', noinput);
p.addParameter('time_window', noinput);
p.addParameter('DEBUG_LEVEL', noinput);


p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = set_input('DEBUG_LEVEL', 2);

% Choose the campaign. It must be recognized by the subfunction
% PATH_BY_CAMPAIGN(). Get a structure representing dates where there is
% GCAS data and the list of files for each date.

campaign = set_input('campaign', 'discover-tx');
campaign_path = path_by_campaign(campaign);
campaign_files = list_gcas_files_by_date(campaign_path);

% Should the output be reduced to vectors of matched data, rather than
% day-by-day 2D subsets?
do_vectorize = set_input('vectorize', false);

% Filtering arguments passed through to gcas_sao_verification
cloud_prod = set_input('cloud_prod', 'omi');
cloud_frac_max = set_input('cloud_frac_max', 0.2);
row_anomaly = set_input('row_anomaly', 'XTrackFlags');
sat_fields = set_input('sat_fields', {'BEHRColumnAmountNO2Trop', 'ColumnAmountNO2Trop'});
time_window = set_input('time_window', 1.5);

pass_through_args = {'cloud_prod', cloud_prod, 'cloud_frac_max', cloud_frac_max,...
                     'row_anomaly', row_anomaly, 'sat_fields', sat_fields,...
                     'time_window', time_window, 'DEBUG_LEVEL', DEBUG_LEVEL};

% Loop through each day, load the right BEHR file, and run
% gcas_sao_verification to compare the GCAS files with the BEHR one
behr_dir = set_input('behr_dir', behr_paths.BEHRMatSubdir('us','monthly'));
Matched_Datas = [];
for a=1:numel(campaign_files)
    Data = load_behr_file_for_gcas(behr_dir, campaign_files(a).date);
    this_matched_data = gcas_sao_verification(campaign_files(a).files, Data, pass_through_args{:});
    % Remove any orbits that don't have any matched data
    xx = true(size(this_matched_data));
    for b=1:numel(this_matched_data)
        xx(b) = ~isempty(this_matched_data(b).Matches);
    end
    
    Matched_Datas = veccat(Matched_Datas, this_matched_data(xx), 'column');
end

% If requested, convert Matched_Datas from a structure of separate days
% where the 2D shape of the pixels is retained to a scalar structure of
% vectors of matched data that can be easily scatter plotted.
if do_vectorize
    Matched_Datas = vectorize_gcas_matches(Matched_Datas);
end


    function inpt = set_input(parameter, default_value)
        if isequal(pout.(parameter), noinput)
            inpt = default_value;
        else
            inpt = pout.(parameter);
        end
    end

end

function files = list_gcas_files_by_date(gcas_path)
E = JLLErrors;
F = dirff(fullfile(gcas_path, '*.h5'));
% Identify the files by their dates
datelist = [];
filecell = {}; % will be a cell array of cell arrays, a cell array of file names for each day
for a=1:numel(F)
    % Extract just the file name itself, in case a date like number shows
    % up somewhere else in the path.
    [~, this_file] = fileparts(F(a).name);
    this_date = datenum(regexp(this_file, '(?<=_)\d\d\d\d\d\d\d\d(?=_)','match','once'), 'yyyymmdd');
    if isempty(this_date)
        % As of R2014b, datenum('') returns an empty matrix
        E.callError('date_not_found', 'Could not find date in file %s', F(a).name);
    end
    
    xx = datelist == this_date;
    if any(xx)
        filecell{xx}{end+1} = F(a).name; %#ok<*AGROW>
    else
        datelist(end+1) = this_date;
        filecell{end+1} = {F(a).name};
    end
end

files = struct('date', num2cell(datelist), 'files', filecell);
end

function Data = load_behr_file_for_gcas(behr_path, date_in)
E = JLLErrors;
file_pattern = sprintf('*_%s.mat', datestr(date_in, 'yyyymmdd'));
F = dir(fullfile(behr_path, file_pattern));
if numel(F) == 1
    D = load(fullfile(behr_path, F(1).name), 'Data');
    Data = D.Data;
else
    E.filenotfound('BEHR file for %s in %s', datestr(date_in), behr_path);
end
end

function p = path_by_campaign(campaign)
E = JLLErrors;
root_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw/';
if ~isempty(regexpi(campaign, 'discover[_\-]?tx'))
    p = fullfile(root_dir, 'DISCOVER-AQ_TX', 'B200', 'GCAS-SAO');
else
    E.badinput('Path for campaign "%s" not defined', campaign);
end
end