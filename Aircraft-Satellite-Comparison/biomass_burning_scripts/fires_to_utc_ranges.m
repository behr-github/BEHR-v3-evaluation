function [ Ranges ] = fires_to_utc_ranges( Fires )
%FIRES_TO_UTC_RANGES Convert the output from find_campaign_fires to a Ranges structure.
%   My function find_campaign_fires outputs selections of campaign data in
%   a structure that's useful for looking at what parts of the data are
%   considered biomass burning (based on CO, ACN, and HCN) but it's not in
%   quite the right format to be fed into any of my existing profile
%   analysis functions as a set of UTC ranges. This will handle the
%   necessary conversion and output the Ranges variable (you'll need to
%   save it manually). Note that the output variable should always be named
%   Ranges for the profile code to work properly. Since the only
%   requirement for the Ranges structure is that different days are
%   represented by different elements in the structure array (i.e.
%   Ranges(1) is the first day, Ranges(2) the second, etc.) and the fields
%   "Date" and "Ranges" it's okay to include the run parameters from
%   find_campaign_fires so that we'll have a record of how these ranges
%   were arrived at.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 June 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(Fires)
    E.badinput('Fires must be a structure');
elseif ~isfield(Fires, 'run_params');
    E.badinput('Fires doesn''t look like a structure output from find_campaign_fires (there''s no "run_params" field)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Fires structures store each day's data in a field named "Flight_yyyymmdd"
% so we'll convert the date in it to the format used in Ranges (mm/dd/yyyy)
fire_fns = glob(fieldnames(Fires),'Flight*');
date_names = cell(size(fire_fns));
for a=1:numel(fire_fns)
    [ind1, ind2] = regexp(fire_fns{a}, '\d\d\d\d\d\d\d\d');
    date_names{a} = datestr(datenum(fire_fns{a}(ind1:ind2),'yyyymmdd'),'mm/dd/yyyy');
end

Ranges = struct('Date', cell(size(date_names)), 'Ranges', cell(size(date_names)), 'Fire_Run_Parameters', cell(size(date_names)));

for a=1:numel(Ranges)
    Ranges(a).Date = date_names{a};
    Ranges(a).Ranges = Fires.(fire_fns{a}).utc_ranges; % the utc_ranges field in Fires already has the right format, we just needed to handle the structure
    Ranges(a).Fire_Run_Parameters = Fires.run_params;
end

end

