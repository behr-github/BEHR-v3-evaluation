function [ varargout ] = convert_units( data, unit_in, unit_out, varargin )
%CONVERT_UNITS Convert data between defined units
%   Takes in a matrix or vector of data and two strings defining the unit
%   the data is in and the unit you want the data to be in.  Will convert
%   between the units as long as they are defined and in the same category.
%   You can also pass 'ls', 'list', or 'listunits' as the only input to see
%   what units it knows and their categories.
%
%   Parameters:
%
%       'case' - logical, default is true, controls whether unit
%       abbreviations are matched with case sensitive comparison or not.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 June 2015

E = JLLErrors;

p = inputParser;
p.addParameter('case', true);

p.parse(varargin{:});
pout = p.Results;

abbrev_case = pout.case;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% UNIT DEFINITION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add units to a category (or make a new category here). The code will
% require that both the in and out unit be in the same category.
wstate = warning('off','unitid:dup_abbrev');
MixingRatios = UnitConversion;
MixingRatios.abbrev_case_sensitive = abbrev_case;
MixingRatios.add_unit('parts-per-trillion', 'ppt', 1e-12);
MixingRatios.add_unit('parts-per-trillion-volume', 'pptv', 1e-12);
MixingRatios.add_unit('parts-per-billion', 'ppb', 1e-9);
MixingRatios.add_unit('parts-per-billion-volume', 'ppbv', 1e-9);
MixingRatios.add_unit('parts-per-million', 'ppm', 1e-6);
MixingRatios.add_unit('parts-per-million-volume', 'ppmv', 1e-6);
MixingRatios.add_unit('parts-per-part', 'ppp', 1);
MixingRatios.add_unit('parts-per-part-volume', 'pppv', 1);
UnitCategories.mixing_ratios = MixingRatios;

Pressures = UnitConversion;
Pressures.abbrev_case_sensitive = abbrev_case;
Pressures.add_unit('pascal', 'Pa', 1);
Pressures.add_unit('bar', 'b', 1e5);
Pressures.add_unit('atmosphere', 'atm', 101325);
Pressures.add_unit('Torr', 'torr', 101325/760);
UnitCategories.pressure = Pressures;

Lengths = UnitConversion;
Lengths.abbrev_case_sensitive = abbrev_case;
Lengths.add_unit('meter', 'm', 1);
Lengths.add_unit('metre', 'm', 1);
UnitCategories.lengths = Lengths;
warning(wstate);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(data) && ismember(data, {'ls','list','listunits'})
    list_units(UnitCategories)
    return
elseif ~isnumeric(data) 
    E.badinput('data is expected to be numeric')
elseif ~ischar(unit_in) || ~ischar(unit_out)
    E.badinput('The units are expected to be given as strings')
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONVERTING %%%%%
%%%%%%%%%%%%%%%%%%%%%%

fns = fieldnames(UnitCategories);
for a = 1:numel(fns)
    in_cat = [UnitCategories.(fns{a}).is_unit_defined(unit_in), UnitCategories.(fns{a}).is_unit_defined(unit_out)];
    if all(in_cat)
        varargout{1} = data * UnitCategories.(fns{a}).get_conversion(unit_in, unit_out);
        return
    end
end

E.callError('unit_not_found', 'Could not find a conversion from %s to %s', unit_in, unit_out);

end

function list_units(UnitCategories)
fns = fieldnames(UnitCategories);
for a=1:numel(fns)
    fprintf('%s:\n', fns{a});
    abbrev = UnitCategories.(fns{a}).list_abbreviations();
    names = UnitCategories.(fns{a}).list_long_names();
    for b=1:numel(names)
        fprintf('     %s (%s)\n', names{b}, abbrev{b});
    end
end
end