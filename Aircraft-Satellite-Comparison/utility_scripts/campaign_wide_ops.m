function [ Out ] = campaign_wide_ops( campaign_name, species_in, op, varargin )
%CAMPAIGN_WIDE_OPS Operations performed on a campaign's worth of data.
%   CAMPAIGN_WIDE_OPS( CAMPAIGN_NAME, SPECIES_IN, OP) 
%   When it is desirable to look at data for an entire campaign, it is very
%   inconvinient that each day is its own Merge file. This function can run
%   several different operations on all the data in the campaign. Takes at
%   least 3 arguments:
%       campaign_name - the name of the campaign to operate on. Must be
%       recognized by merge_field_names.
%
%       species - Essentially what field to operate on in the Merge file.
%       This will try first to find this in the Names structure returned by
%       merge_field_names (so no2_lif will get you the appropriate NO2 LIF
%       field for the given campaign), then try to find that field in the
%       Merge files for the campaign.
%
%       op - the operation to be carried out. These will be described in
%       more detail below. Some operations may require additional arguments
%       to be passed; see their descriptions below for this information.
%
%   This function will return one or more outputs, again this depends on
%   the operation.
%
%   Possible operations are:
%
%       'cat' - simply concatenate the designated species into one long
%       vector. Will return this vector, plus a vector of UTC times and a
%       cell array of dates describing the date/time of each measurement.
%       Requires no additional arguments.
%
%       'dayavg' - provides one average values per day.
%
%       'bin' - bins the data by GPS altitude. Requires 1 additional
%       argument, the bin width in km. Outputs the bin median values, bin
%       midpoint altitudes, and upper and lower quartiles for each bin.
%
%       'bin_rolling' - bins the data with rolling averaging windows (using
%       bin_rolling_vertical_profile). Need the bin width and bin spacing
%       (in km) passed as additional arguments. Outputs the bin median
%       values, bin midpoint altitudes, and upper and lower quartiles for
%       each bin.
%
%       'bin_pres' - bins by pressure, into OMI standard pressure bins.
%       Requires no additional arguments. Returns the bin median values,
%       bin pressure levels, and upper and lower quartiles.
%
%   Parameters:
%
%       datefmt - changes how the date is given. Can either be a format
%       number or string recognized by datestr, or the string 'datenum',
%       which will output the dates as a vector of date numbers, rather
%       than a cell array of date strings.
%
%       binmode - can be 'median' or 'mean', changes how bins are
%       calculated if a binning operation is chosen. Default is 'median'.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(3,Inf);
p = inputParser;
p.addOptional('oparg1',[]);
p.addOptional('oparg2',[]);
p.addParameter('datefmt','yyyy-mm-dd');
p.addParameter('binmode','median');
p.addParameter('debug',1);
p.parse(varargin{:});
pout = p.Results;

datefmt = pout.datefmt;
binmode = pout.binmode;
DEBUG_LEVEL = pout.debug;
if ~ismember(binmode, {'mean','median'})
    E.badinput('BINMODE can only be ''mean'' or ''median''')
end

% Confirm that the operation requested is allowed, then check that enough
% arguments were passed.

allowed_ops =   {'cat','dayavg','bin','bin_rolling','bin_pres'};
req_args =      [0,     0,       1,    2,            0];
opargs = {pout.oparg1, pout.oparg2};
opargs(iscellcontents(opargs,'isempty')) = [];

op = lower(op);
if ~ismember(op, allowed_ops)
    E.badinput('op %s is not one of the expected values: %s', op, strjoin(allowed_ops, ', '));
end
xx = strcmp(op, allowed_ops);
if numel(opargs) ~= req_args(xx)
    E.badinput('op %s requires exactly %d additional arguments', op, req_args(xx));
end

% We'll check that species is a valid input later, we need a Merge file
% loaded to do that. For now, we need to convert it to a cell array if it
% isn't.
if ischar(species_in)
    species_in = {species_in};
elseif ~iscellstr(species_in)
    E.badinput('SPECIES_IN must be a string or cell array of strings')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

[Names, ~, merge_dir] = merge_field_names(campaign_name);

F = dir(fullfile(merge_dir,'*.mat'));
Ms(numel(F)) = struct('Merge',[]);
merge_fields = cell(1,numel(species_in));
for a=1:numel(F)
    if DEBUG_LEVEL > 0; fprintf('Loading %s\n', F(a).name); end
    load(fullfile(merge_dir, F(a).name),'Merge'); % brings the variable Merge into the workspace
    % Now we'll handle checking species against Names and Merge fields. On
    % successive loops, ensure that the field name is still defined for the
    % new Merge. Do this for each species requested
    
    for s = 1:numel(species_in)
        species = species_in{s};
        if a == 1
            if isfield(Names, species)
                merge_fields{s} = Names.(species);
            elseif isfield(Merge.Data, species)
                merge_fields{s} = species;
            else
                E.badinput('species %s is not a defined field in Names or Merge.Data for the campaign %s', species, campaign_name);
            end
        elseif ~isfield(Merge.Data, merge_fields{s})
            E.callError('inconsistent_merge','Later Merge.Data does not have the field %s', merge_fields{s});
        end
    end
    
    Ms(a).Merge = Merge;
end

% At this point, the operation specified must be called.
switch op
    case 'cat'
        Out = concatenate_output_merges(Ms, merge_fields, species_in, datefmt);
    case 'dayavg'
        Out = day_avg_merges(Ms, merge_fields, species_in, datefmt);
    case 'bin'
        Out = bin_merges(Ms, merge_fields, Names, species_in, 'alt', binmode, opargs{1});
    case 'bin_rolling'
        Out = bin_merges(Ms, merge_fields, Names, species_in, 'rolling', binmode, opargs{1:2});
    case 'bin_pres'
        Out = bin_merges(Ms, merge_fields, Names, species_in, 'pres', binmode);
    otherwise
        E.badinput('Operation %s not recognized',op);
end


end

function varargout = concatenate_merges(Ms, merge_fields, datefmt)
catted_species = cell(1,numel(merge_fields));
utcs = [];
if strcmpi(datefmt,'datenum')
    dates = [];
else
    dates = {};
end
% Loop through all the merges, appending the requested species, times, and
% dates

for a=1:numel(Ms)
    for b=1:numel(merge_fields)
        data = remove_merge_fills(Ms(a).Merge, merge_fields{b});
        catted_species{b} = cat(2, catted_species{b}, data);
    end
    if isfield(Ms(1).Merge.Data.UTC,'Values')
        utcs = cat(2, utcs, Ms(a).Merge.Data.UTC.Values);
    end
    
    if strcmpi(datefmt,'datenum')
        merge_date = datenum(Ms(a).Merge.metadata.date,'yyyy-mm-dd');
    else
        merge_date = {datestr(Ms(a).Merge.metadata.date, datefmt)};
    end
    this_dates = repmat(merge_date,1,numel(data));
    dates = cat(2, dates, this_dates);
    
end

varargout = [catted_species, {utcs}, {dates}]; %#ok<VARARG> each concatenated element is a cell
end

function Out = day_avg_merges(Ms, merge_fields, species_in, datefmt)
all_species = cell(1,numel(merge_fields));
[all_species{:}, utcs, dates] = concatenate_merges(Ms, merge_fields, datefmt);
Out.data = make_empty_struct_from_cell(species_in);

dnums = datenum(dates);
udnums = unique(dnums);
for a=1:numel(merge_fields)
    val = nan(size(udnums));
    for b=1:numel(udnums)
        dd = dnums == udnums(b);
        val(b) = nanmean(all_species{a}(dd));
    end
    Out.data.(species_in{a}) = val;
end
Out.utcs = utcs;
Out.dates = cellstr(datestr(udnums,datefmt))';

end

function Out = concatenate_output_merges(Ms, merge_fields, species_in, datefmt)
% Concatenates the fields given in MERGE_FIELDS across the merges contained
% in MS. Output structure uses the names in SPECIES_IN so that the user
% doesn't have to translate the field names in the output structure.
all_species = cell(1,numel(merge_fields));
[all_species{:}, utcs, dates] = concatenate_merges(Ms, merge_fields, datefmt);
Out.data = make_empty_struct_from_cell(species_in);
for a=1:numel(merge_fields)
    Out.data.(species_in{a}) = all_species{a};
end
Out.utcs = utcs;
Out.dates = dates;
end

function Out = bin_merges(Ms, merge_fields, Names, species_in, bin_type, bin_mode, bin_width, bin_spacing )
E = JLLErrors;

all_species = cell(1, numel(merge_fields));
bin_vals = cell(1, numel(merge_fields));
bin_quarts = cell(1, numel(merge_fields));
[all_species{:}] = concatenate_merges(Ms, merge_fields, 29); % 29 is the yyyy-mm-dd date format
if strcmpi(bin_type,'pres')
    all_pres = concatenate_merges(Ms, {Names.pressure}, 29);
else
    all_alts = concatenate_merges(Ms, {Names.gps_alt}, 29);
end
for a=1:numel(all_species)
    switch lower(bin_type)
        case 'alt'
            [bin_vals{a}, Out.bin_alts, bin_quarts{a}] = bin_vertical_profile(all_alts, all_species{a}, bin_width, 'binmode', bin_mode);
        case 'rolling'
            [bin_vals{a}, Out.bin_alts, bin_quarts{a}] = bin_rolling_vertical_profile(all_alts, all_species{a}, bin_width, bin_spacing, 'binmode', bin_mode);
        case 'pres'
            [bin_vals{a}, Out.bin_pres, bin_quarts{a}] = bin_omisp_pressure(all_pres, all_species{a}, 'binmode', bin_mode);
        otherwise
            E.notimplemented('Binning type %s not recognized', bin_type);
    end
end

Out.data = make_empty_struct_from_cell(species_in);
Out.quartiles = make_empty_struct_from_cell(species_in);
for a=1:numel(species_in)
    Out.data.(species_in{a}) = bin_vals{a};
    Out.quartiles.(species_in{a}) = bin_quarts{a};
end

end
