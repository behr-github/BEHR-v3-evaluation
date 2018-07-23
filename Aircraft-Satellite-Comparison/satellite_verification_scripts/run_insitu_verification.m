function [all_profs_final, all_profs_struct] = run_insitu_verification(campaign, behr_prof_mode, varargin)
%RUN_INSITU_VERIFICATION Driver function for verify_sat_vs_aircraft_vcds()
%   [ ALL_PROFS_FINAL, ALL_PROFS_STRUCT ] = RUN_INSITU_VERIFICATION(
%   CAMPAIGN, BEHR_PROF_MODE ) Internally, this will load each Merge file
%   for the given CAMPAIGN (which must be a campaign name that
%   merge_field_names() recognizes) and the corresponding BEHR file and
%   pass them to verify_sat_vs_aircraft. ALL_PROFS_FINAL will be a scalar
%   struct with the same fields as the ALL_PROFS output from
%   verify_sat_vs_aircraft(), but with each field being the concatenation
%   of the outputs from all the calls to verify_sat_vs_aircraft().
%   Similarly, ALL_PROFS_STRUCT is a vector structure that is the
%   concatenation of all the PROFS_STRUCT outputs from
%   verify_sat_vs_aircraft(). In its case, it is the structures themselves
%   that are concatenated, not their fields.
%
%   All parameter arguments are passed down into verify_sat_vs_aircraft(),
%   so see its help text for a description of the available options. This
%   does have several parameters that it accepts:
%
%       'behr_dir' - the directory to look for BEHR files in. If not given,
%       or given as an empty string, the directory is automatically
%       determined from behr_paths. If given, 'behr_prefix' and
%       'behr_version' must be given as well.
%
%       'behr_prefix' - the prefix string for the BEHR filenames, i.e. what
%       you would pass as the third argument to BEHR_FILENAME( date,
%       'prefix' ) in order to generate a file name. If not given or given
%       as an empty string, then the filename is generated for the current
%       BEHR version for the apropriate profile mode and the 'us' region.
%       If given, 'behr_dir' and 'behr_version' must be given as well.
%
%       'behr_version' - the version string found in the BEHR files to
%       load, e.g. 'v3-0A' or 'v2-1C'. If not given or given as an empty
%       string, then the current BEHR version is used. If given, then
%       'behr_dir' and 'behr_prefix' must be given also.

[merge_files, behr_files] = list_files_to_load(campaign, behr_prof_mode, varargin{:});

for i_day = 1:numel(merge_files)
    [profs_final, profs_struct] = call_verify(campaign, merge_files{i_day}, behr_files{i_day}, varargin{:});
    % Join together the final output so that all the VCDs and related data
    % are concatenated into vectors (easy to plot) and the detailed
    % structures are concatenated into a single structure.
    if i_day == 1
        all_profs_final = profs_final;
        all_profs_struct = profs_struct;
    else
        all_profs_final = cat_fields(all_profs_final, profs_final);
        all_profs_struct = veccat(all_profs_struct, profs_struct);
    end
end

end

function [merge_files, behr_files] = list_files_to_load(campaign, behr_prof_mode, varargin)

E = JLLErrors;

p = inputParser;
p.addParameter('behr_dir', '');
p.addParameter('behr_version','');
p.addParameter('behr_prefix','');
p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;
behr_dir = pout.behr_dir;
behr_vers = pout.behr_version;
behr_prefix = pout.behr_prefix;

behr_file_params_check = [isempty(behr_dir), isempty(behr_vers), isempty(behr_prefix)];
if (any(behr_file_params_check) && ~all(behr_file_params_check)) || (any(~behr_file_params_check) && ~all(~behr_file_params_check))
    E.badinput('All or none of ''behr_dir'', ''behr_version'', and ''behr_prefix'' must be specified');
end

[~, ~, merge_dir] = merge_field_names(campaign);

% Make a list of the available aircraft merge files
merge_files = dirff(fullfile(merge_dir, '*.mat'));
merge_files = {merge_files.name};
% Make a list of the corresponding BEHR files
behr_files = cell(size(merge_files));
if isempty(behr_dir)
    behr_dir = behr_paths.BEHRMatSubdir('us', behr_prof_mode);
end
for i_merge = 1:numel(merge_files)
    % Extract the date of the Merge file, which is always at the end of the
    % file name right before the .mat extension
    merge_date = regexp(merge_files{i_merge}, '\d\d\d\d_\d\d_\d\d(?=\.mat$)','once','match');
    if isempty(behr_prefix)
        behr_file = behr_filename(datenum(merge_date), behr_prof_mode, 'us');
    else
        behr_file = behr_filename(datenum(merge_date), 'prefix', behr_prefix, 'mat', behr_vers);
    end
    behr_files{i_merge} = fullfile(behr_dir, behr_file);
end

end

function [profs_final, profs_struct] = call_verify(campaign, merge_file, behr_file, varargin)
E = JLLErrors;
AIR = load(merge_file,'Merge');
try
    BEHR = load(behr_file,'Data');
catch err
    % If we can't load the BEHR file, give a clearer error message
    if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
        prof_mode = regexp(behr_file, '(daily|monthly)', 'match', 'once');
        behr_date = datestr(datenum(regexp(behr_file, '\d\d\d\d\d\d\d\d', 'match', 'once'), 'yyyymmdd'));
        E.filenotfound('No %s BEHR file available for %s\n', prof_mode, behr_date);
    else
        rethrow(err)
    end
end

[profs_final, profs_struct] = verify_sat_vs_aircraft_vcds(BEHR.Data, AIR.Merge, campaign, varargin{:});
end