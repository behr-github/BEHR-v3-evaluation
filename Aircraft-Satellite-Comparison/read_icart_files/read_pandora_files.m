function read_pandora_files(input_base_dir, output_base_dir, campaign, varargin)
%READ_PANDORA_FILES Read a series of Pandora ICARTT files into .mat files
%   READ_PANDORA_FILES( INPUT_BASE_DIR, OUTPUT_BASE_DIR, CAMPAIGN ) Will
%   read all *.ict files in fullfile(INPUT_BASE_DIR, CAMPAIGN, 'Pandora',
%   'Filtered') with read_icartt_file() and save each resulting Merge
%   structure in fullfile(OUTPUT_BASE_DIR, CAMPAIGN, 'Pandora',
%   'Filtered'). Each merge file is given the name
%   'Pandora_<location>_P<n>.mat' where <location> is the lower case
%   location name in the icartt file name and <n> is the Pandora number. If
%   no Pandora number is given in the file name, <n> is set to 'X'.
%
%   Parameters:
%       'field_mismatch' - controls what read_icartt_file() does if the
%       list of variable at the top of the file does not match the data
%       table header. This is passed through as the second argument to
%       read_icartt_file, see it's documentation for possible values.

E = JLLErrors;

if nargin == 0
    input_base_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignRaw';
    output_base_dir = '/Volumes/share2/USERS/LaughnerJ/CampaignInstrMats';
    campaign = 'DISCOVER-AQ_CO';
elseif nargin < 2
    E.badinput('Must give either 0 or >= 2 inputs');
end

if ~exist(input_base_dir, 'dir')
    E.badinput('INPUT_DIR (%s) is not a valid directory', input_base_dir)
elseif ~exist(output_base_dir, 'dir')
    E.badinput('OUTPUT_DIR (%s) is not a valid directory', output_base_dir)
end

p = inputParser;
p.addParameter('field_mismatch','');
p.parse(varargin{:});
pout = p.Results;

field_mismatch_preference = pout.field_mismatch;

% List the icartt files in the input directory, assuming that the directory
% structure is <base_dir>/<campaign>/Pandora/Filtered
icartt_files = dirff(fullfile(input_base_dir, campaign, 'Pandora', 'Filtered', '*.ict'));

full_out_path = fullfile(output_base_dir, campaign, 'Pandora', 'Filtered');

for i_file = 1:numel(icartt_files)
    [Merge, header, DataTable] = read_icartt_file(icartt_files(i_file).name, field_mismatch_preference);
    
    % Figure out the save name. Pandora files give one file per pandora for
    % the whole campaign, so rather than the file date, we need to extract
    % the location from the file name. Just to make things even more fun,
    % in some campaigns there were multiple pandoras at a given location,
    % so we also want the pandora number, but Maryland doesn't give the
    % pandora number in the file name. So for Maryland we'll just plug in
    % PX for the pandora number instead of e.g. P6 or P21.
    %
    % An example Pandora file name is:
    % discoveraq-pandora-NO2-P31_ground-Rocky-Flats_20140622_R0_thru20140814-Filtered.ict
    % Note that the location ("Boulder") always comes after ground- and is
    % followed by an underscore; inside the name, words are separated by
    % dashes. The capitalization of "ground" varies from campaign to
    % campaign, so we have to use case insensitive matching.
    [~,file_name] = fileparts(icartt_files(i_file).name);
    location_name = regexpi(file_name, '(?<=ground\-)[a-zA-Z0-9\-]+(?=_)', 'match', 'once');
    
    % Try to get the pandora number, Pn or Pnn which comes before "ground".
    % If no match, that's okay, just insert a dummy number.
    pandora_number = regexpi(file_name, 'P\d{1,3}(?=_ground)','match','once');
    if isempty(pandora_number)
        pandora_number = 'PX';
    end
    
    save_name = sprintf('Pandora_%s_%s.mat', lower(location_name), upper(pandora_number));
    full_save_name = fullfile(full_out_path, save_name);
    fprintf('Saving as %s\n', full_save_name);
    save(full_save_name, 'Merge', 'header', 'DataTable');
end

end

