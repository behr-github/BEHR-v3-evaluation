function matches = run_pandora_verification(campaign, varargin)
%RUN_PANDORA_VERIFICATION Run verify_sat_vs_pandora for a whole campaign.
%   MATCHES = RUN_PANDORA_VERIFICATION( CAMPAIGN ) Runs verification of OMI
%   columns against Pandora columns for all Pandora sites available during
%   CAMPAIGN. CAMPAIGN must be one of 'discover_md', 'discover_ca',
%   'discover_tx', or 'discover_co'. MATCHES will be a structure with
%   fields 'sp_no2', 'behr_no2', 'pandora_no2', 'omi_time', and 'details'.
%   'sp_no2', 'behr_no2', and 'pandora_no2' are the total VCDs from the
%   NASA SP and BEHR products, and the Pandora spectrometers, respectively.
%   'omi_time' is a datenum giving the overpass time for the pixels matched
%   to the Pandora location. 'details' is another structure that contains
%   the BEHR Data structure cut down to just the matched pixels.
%
%   All parameters are passed through to verify_sat_vs_pandora. The
%   following parameters are recognized by this function:
%
%   Of the following four parameters either 'behr_prof_mode' or
%   'behr_prefix', 'behr_version', and 'behr_dir' must be given. If all are
%   given, the latter three take precedence.
%
%       'behr_prof_mode' - must be 'daily' or 'monthly'. Determines which
%       version of BEHR to compare against the Pandoras.
%
%       'behr_prefix' - the prefix of the BEHR file (the part before the
%       version string) to load.
%
%       'behr_version' - which BEHR version to load.
%
%       'behr_dir' - which directory to load BEHR files from.
%
%   If only 'behr_prof_mode' is given, then BEHR files are loaded with:
%
%       Data = load_behr_file(<date>, behr_prof_mode, 'us')
%
%   where <date> is detemined automatically from the Pandora data. If the
%   other three are given, then the file given by:
%
%       file_name = behr_filename(<date>, 'prefix', behr_prefix, '.mat',
%       behr_version)
%
%   is loaded from the 'behr_dir' directory.

pandora_dir = verify_pandora_dir(campaign);
pandora_mats = dirff(fullfile(pandora_dir, '*.mat'));
matches = [];
for i_file = 1:numel(pandora_mats)
    fprintf('Loading Pandora merge file %s...\n', pandora_mats(i_file).name);
    P = load(pandora_mats(i_file).name);
    Merge = P.Merge;
    
    % Iterate over the dates in the Merge structure - Pandora files place
    % an entire campaign's worth of data in one file per site. Try to load
    % the proper BEHR file for each day, if one cannot be found, try the
    % next.
    pandora_dates = unique(floor(pandora_gmt_to_datenum(remove_merge_fills(Merge,'DateGMT'))));
    for i_date = 1:numel(pandora_dates)
        fprintf('  Loading BEHR data for %s...\n', datestr(pandora_dates(i_date)));
        try
            Data = load_behr_for_comparison(pandora_dates(i_date), varargin{:});
        catch err
            if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
                fprintf('    No BEHR files for %s\n', datestr(pandora_dates(i_date)));
                continue
            else
                rethrow(err)
            end
        end
        fprintf('    Matching Pandora to OMI...\n');
        this_match = verify_sat_vs_pandora(Merge, Data, varargin{:});
        if isempty(matches)
            matches = this_match;
        else
            matches = cat_fields(matches, this_match);
        end
    end
end

end

function pandora_dir = verify_pandora_dir(campaign)
E = JLLErrors;
pandora_base_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignInstrMats';
campaign = strrep(upper(campaign),'DISCOVER','DISCOVER-AQ');
pandora_dir = fullfile(pandora_base_dir, campaign, 'Pandora', 'Filtered');
if ~exist(pandora_dir, 'dir')
    E.badinput('Could not find Pandora directory for "%s" campaign (looking for %s)', campaign, pandora_dir);
end
end

function Data = load_behr_for_comparison(behr_date, varargin)
E = JLLErrors;

p = inputParser;
p.addParameter('behr_prof_mode', '');
p.addParameter('behr_prefix', '');
p.addParameter('behr_version','');
p.addParameter('behr_dir', '');
p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

behr_prof_mode = pout.behr_prof_mode;
behr_prefix = pout.behr_prefix;
behr_vers = pout.behr_version;
behr_dir = pout.behr_dir;

% Check that all or none of prefix, version, and dir are given.
input_check = [isempty(behr_prefix), isempty(behr_vers), isempty(behr_dir)];
if any(input_check) && ~all(input_check)
    E.badinput('All or none of the parameters ''behr_prefix'', ''behr_vers'', ''behr_dir'' must be given');
elseif ~all(input_check) && isempty(behr_prof_mode)
    E.badinput('''behr_prof_mode'' must be given in ''behr_prefix'', ''behr_vers'', and ''behr_dir'' are not');
end

if all(input_check)
    % Load the most recent version
    Data = load_behr_file(behr_date, behr_prof_mode, 'us');
else
    % Otherwise load the data from the specified location
    file_name = behr_filename(behr_date, 'prefix', behr_prefix, '.mat', behr_vers);
    D = load(fullfile(behr_dir, file_name), 'Data');
    Data = D.Data;
end


end