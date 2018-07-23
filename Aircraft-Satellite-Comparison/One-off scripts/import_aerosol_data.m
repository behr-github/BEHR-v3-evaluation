function [  ] = import_aerosol_data( campaign_name, xls_file )
%IMPORT_AEROSOL_DATA Import aerosol data from Lee Thornhill
%   Lee Thornhill, Andreas Byersdorf, etc. were kind enough to send me blue
%   wavelength extinction for DC3, SEAC4RS, and the DISCOVER campaigns.
%   Blue extinction data will be better to use b/c NO2 absorbs in the
%   violet/blue wavelengths, so this will more accurately reflect the
%   aerosol effects on NO2 retrievals.
%
%   This one-off script will read the excel sheet Lee sent and match the
%   data to existing Merge structures. It will save the new .mat file in a
%   temporary location so that you can confirm that everything imported
%   correctly.
%
%   This function will probably not be needed again, it exists more for
%   documentation of how the import was handled in case of future concerns.
%
%   Josh Laughner <joshlaugh5@gmail.com> 13 Aug 2015

E = JLLErrors;

% Check input quickly, both should be strings, and xls_file must be a
% regular file

if ~ischar(campaign_name) || ~ischar(xls_file)
    E.badinput('Both inputs must be strings');
elseif ~exist(xls_file, 'file')
    E.badinput('%s does not exist', xls_file);
end

% File paths

% Location where the .mat files should be saved, temporaily until the user
% is satisfied that everything is appended correctly. Set to '' (an empty
% string) to directly overwrite the existing Merge files.
temp_dir = '~/Documents/MATLAB/NO2 Profiles/Temp merges';

% Main function

[~, ~, merge_dir] = merge_field_names(campaign_name);
F = dir(fullfile(merge_dir,'*.mat'));
merge_filenames = {F.name};

% The sheets are named by campaign, SEAC4RS, DC3, DAQ-SS where SS is the
% state. We can convert the campaign name.

sheet_name = regexprep(upper(campaign_name), 'DISCOVER','DAQ');
if ~ismember(sheet_name, {'SEAC4RS','DC3','DAQ-MD','DAQ-CA','DAQ-TX','DAQ-CO'});
    E.badinput('The campaign name cannot be converted to one of the assumed sheet names in the file')
end

% Now load the sheet as a cell array (this will take a bit). Each flight is
% 7 columns, counting the column describing the flight. Empty cells will be
% imported as NaNs. The second cell from the top in the first column of
% each will have the month and day specified, which we can use to load the
% merge file. The next 6 columns are, in order, time in sec after midnight
% UTC, ambient RH, fRH_amb, ambient scattering, dry scattering, and dry
% absorbance.
%
% In order we need to 1) the flights data as vectors (Colorado will need
% special handling b/c of multiple flights per day) 2) find the
% corresponding day's Merge file and load it, and 3) match the data up by
% UTC time.

fprintf('Reading the excel file (be patient)... ');
[~,~,sheet] = xlsread(xls_file, sheet_name);
fprintf('Done.\n');

num_col = size(sheet,2);

a=1;
while a <= num_col
    flight_date = regexprep(sheet{2,a},'[/-]','');
    
    fprintf('Reading data from %s\n',flight_date); 
    flight_utc = cell2mat(sheet(2:end,a+1));
    flight_RH_amb = cell2mat(sheet(2:end,a+2));
    flight_fRH = cell2mat(sheet(2:end,a+3));
    flight_sc_amb = cell2mat(sheet(2:end,a+4));
    flight_sc_dry = cell2mat(sheet(2:end,a+5));
    flight_abs_dry = cell2mat(sheet(2:end,a+6));
    
    % I manually added "L1" and "L2" in the cell below the date for the
    % Colorado days split into two legs. This will check if that's there,
    % and if so, it knows it needs to concatenate the two legs.
    if strcmpi(sheet_name, 'DAQ-CO') && ischar(sheet{3,a}) && ~isempty(regexp(sheet{3,a},'L\d', 'once'))
        fprintf('\t(reading second leg)\n');
        flight_utc = cat(1, flight_utc, cell2mat(sheet(2:end,a+8)));
        flight_RH_amb = cat(1, flight_RH_amb, cell2mat(sheet(2:end,a+9)));
        flight_fRH = cat(1, flight_fRH, cell2mat(sheet(2:end,a+10)));
        flight_sc_amb = cat(1, flight_sc_amb, cell2mat(sheet(2:end,a+11)));
        flight_sc_dry = cat(1, flight_sc_dry, cell2mat(sheet(2:end,a+12)));
        flight_abs_dry = cat(1, flight_abs_dry, cell2mat(sheet(2:end,a+13)));
        
        a=a+14;
    else
        a=a+7;
    end
    
    % If the UTC is a NaN, that probably means that this is a cell that
    % exists simply to keep the cell array square, so remove it.
    nans = isnan(flight_utc);
    flight_utc(nans) = [];
    flight_RH_amb(nans) = [];
    flight_fRH(nans) = [];
    flight_sc_amb(nans) = [];
    flight_sc_dry(nans) = [];
    flight_abs_dry(nans) = [];
    
    % Load the Merge .mat file, which will contain the Merge structure,
    % DataTable matrix, and header cell array.
    
    this_merge_date = sprintf('%s_%s.mat',flight_date(1:2),flight_date(3:4));
    xx = ~iscellcontents(strfind(merge_filenames, this_merge_date), 'isempty');
    if sum(xx) > 1
        E.toomanyfiles(sprintf('Merge file for %s: %s', campaign_name, flight_date));
    elseif sum(xx) < 1
        E.filenotfound(sprintf('Merge file for %s: %s', campaign_name, flight_date));
    end
    
    fprintf('Loading Merge file...\n');
    load(fullfile(merge_dir, F(xx).name));
    
    % Make the new Merge fields and DataTable columns by matching
    % measurements up to UTC times
    merge_utc = Merge.Data.UTC.Values;
    
    merge_RH_amb = nan(size(merge_utc));
    merge_fRH = nan(size(merge_utc));
    merge_sc_amb = nan(size(merge_utc));
    merge_sc_dry = nan(size(merge_utc));
    merge_abs_dry = nan(size(merge_utc));
    
    fprintf('Appending data, matching UTCs\n');
    for b=1:numel(merge_utc)
        uu = flight_utc == merge_utc(b);
        if sum(uu) == 0
            continue
        elseif sum(uu) > 1
            E.badvar('uu',sprintf('More than one flight UTC matched the merge_utc %d',merge_utc(b)));
        end
        
        merge_RH_amb(b) = flight_RH_amb(uu);
        merge_fRH(b) = flight_fRH(uu);
        merge_sc_amb(b) = flight_sc_amb(uu);
        merge_sc_dry(b) = flight_sc_dry(uu);
        merge_abs_dry(b) = flight_abs_dry(uu);
    end
    
    Merge.Data.Lee_RH_amb.Unit = '%';
    Merge.Data.Lee_RH_amb.Fill = -9999;
    merge_RH_amb(isnan(merge_RH_amb)) = -9999;
    Merge.Data.Lee_RH_amb.Values = merge_RH_amb;
    
    Merge.Data.Lee_fRH.Unit = 'unitless';
    Merge.Data.Lee_fRH.Fill = -9999;
    merge_fRH(isnan(merge_fRH)) = -9999;
    Merge.Data.Lee_fRH.Values = merge_fRH;
    
    Merge.Data.Lee_sc450nm_amb.Unit = 'Mm^-1';
    Merge.Data.Lee_sc450nm_amb.Fill = -9999;
    merge_sc_amb(isnan(merge_sc_amb)) = -9999;
    Merge.Data.Lee_sc450nm_amb.Values = merge_sc_amb;
    
    Merge.Data.Lee_sc450nm_dry.Unit = 'Mm^-1';
    Merge.Data.Lee_sc450nm_dry.Fill = -9999;
    merge_sc_dry(isnan(merge_sc_dry)) = -9999;
    Merge.Data.Lee_sc450nm_dry.Values = merge_sc_dry;
    
    Merge.Data.Lee_abs450nm_dry.Unit = 'Mm^-1';
    Merge.Data.Lee_abs450nm_dry.Fill = -9999;
    merge_abs_dry(isnan(merge_abs_dry)) = -9999;
    Merge.Data.Lee_abs450nm_dry.Values = merge_abs_dry;
    
    % Calculate extinction as scattering + absorption. Lee indicated that
    % there's no way to get truly ambient absorption, so just to use dry
    % absorption. Usually scattering dominates anyway, so it shouldn't
    % matter much.
    Merge.Data.Lee_ext450nm_amb.Unit = 'Mm^-1';
    Merge.Data.Lee_ext450nm_amb.Fill = -9999;
    Merge.Data.Lee_ext450nm_amb.Values = merge_sc_amb + merge_abs_dry;
    fills = merge_sc_amb < -9000 | merge_abs_dry < -9000;
    Merge.Data.Lee_ext450nm_amb.Values(fills) = -9999;
    
    % Add these to the DataTable and include the names in the header. The
    % Values vectors in the structure are row vectors, so transpose them
    % first.
    DataTable = cat(2, DataTable, merge_RH_amb', merge_fRH', merge_sc_amb', merge_sc_dry', merge_abs_dry', merge_abs_dry' + merge_sc_amb');
    header = cat(2, header, {'Lee_RH_amb','Lee_fRH','Lee_sc450nm_amb','Lee_sc450nm_dry','Lee_abs450nm_dry','Lee_ext450nm_amb'});
    
    if isempty(temp_dir)
        save_name = fullfile(merge_dir, F(xx).name);
    else
        save_name = fullfile(temp_dir, F(xx).name);
    end
    
    fprintf('Saving as %s\n',save_name);
    save(save_name, 'Merge', 'DataTable', 'header');
end

end

