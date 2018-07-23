classdef UnitConversion < handle
    %UnitConversion - Class for defining unit conversions
    %   This class is used to define relationships between units of the
    %   same physical quantity. The user adds individual units with the
    %   ADD_UNIT() method, then retrieves a conversion factor between two
    %   units with the GET_CONVERSION() method. Metric prefixes (from
    %   "atto" to "exa") are defined internally, so one does not need to
    %   define separate units for different prefixes.
    %
    %   Public properties:
    %       abbrev_case_sensitive (default true) controls whether the
    %       instance compares unit abbreviations with or without case
    %       sensitivity (true = with, i.e. ppbV is different from ppbv).
    
    properties
        abbrev_case_sensitive = true;
    end
    
    properties(SetAccess = protected)
        unit_names = {};
    end
    
    properties(Constant = true)
        short_prefixes = struct('a', 1e-18, 'f', 1e-15, 'p', 1e-12, 'n', 1e-9,...
            'u', 1e-6, 'm', 1e-3, 'c', 1e-2, 'd', 1e-1, 'da', 1e1, 'h', 1e2,...
            'k', 1e3, 'M', 1e6, 'G', 1e9, 'T', 1e12, 'P', 1e15, 'E', 1e18);
        long_prefixes = struct('atto', 1e-18, 'femto', 1e-15, 'pico', 1e-12, ...
            'nano', 1e-9, 'micro', 1e-6, 'milli', 1e-3, 'centi', 1e-2, 'deci', 1e-1,...
            'deca', 1e1, 'hecto', 1e2, 'kilo', 1e3', 'mega', 1e6, 'giga', 1e9,...
            'tera', 1e12, 'peta', 1e15, 'exa', 1e18);
    end
    
    methods
        function add_unit(obj, long_name, abbreviation, conv_factor)
            % U.add_unit( LONG_NAME, ABBREVIATION, CONV_FACTOR)
            %
            %   Adds a unit to this instance that it knows how to convert.
            %   LONG_NAME is the full name of the unit (e.g. "meter",
            %   "pascal"). It should be singular. 
            %
            %   ABBREVIATION is the standard abbreviation for the unit
            %   (e.g. "m", "Pa"). 
            %
            %   CONV_FACTOR is the numerical conversion factor that you
            %   would multiply a value by to convert this unit to one with
            %   a conversion factor of 1. For example, if we were defining
            %   conversion among pressure units, and gave Pa a conv_factor
            %   of 1, then an atm would have a conv_factor of 101325.
            %
            %   Note: when trying to convert units, this class ignores case
            %   for long names and allows long names to be pluralized with
            %   either "s" or "es" on the end. However, abbreviations are
            %   matched including case (i.e. "Pa" != "pa"). Metric prefixes
            %   are also defined internally, so you only need to define the
            %   base unit.
            
            if ~ischar(long_name)
                error('unitid:bad_input', 'LONG_NAME must be a string')
            elseif ~ischar(abbreviation)
                error('unitid:bad_input', 'ABBREVIATION must be a string')
            elseif length(long_name) < length(abbreviation)
                warning('unitid:name_length','LONG_NAME should not be shorter than ABBREVIATION')
            elseif ~isnumeric(conv_factor) || ~isscalar(conv_factor)
                error('unitid:bad_input', 'CONV_FACTOR must be a scalar number')
            elseif conv_factor <= 0
                warning('unitid:neg_conv_factor','CONV_FACTOR should not be <= 0')
            end
            
            long_name = lower(long_name);
            
            if ismember(long_name, obj.list_long_names())
                warning('unitid:dup_name','Duplicate long name: "%s"', long_name);
            elseif ismember(abbreviation, obj.list_abbreviations())
                warning('unitid:dup_abbrev','Duplicate abbreviation: "%s"', long_name);
            end
            
            obj.unit_names{end+1} = struct('name', long_name, 'abbr', abbreviation, 'conv', conv_factor);
        end
        
        function conv = get_conversion(obj, old_unit, new_unit)
            % U.get_conversion(old_unit, new_unit)
            %
            %   Get the multiplicative factor that will convert OLD_UNIT to
            %   NEW_UNIT. This will match both OLD_UNIT and NEW_UNIT
            %   against the units defined internally, accounting for metric
            %   prefixes. To be specific, the number returned if
            %   F_old/F_new, the conversion factors for the old and new
            %   unit (multiplied by metric prefixes) respectively. 
            %
            %   Note: this will expect long metric prefixes with long
            %   names, or short metric prefixes with abbreviations. (So
            %   "mm" or "millimeter" are fine, but "mmeter" or "millim" are
            %   not.) Also, "u" is used as the short prefix for "micro".
            old_conv = obj.get_conv_factor(old_unit);
            new_conv = obj.get_conv_factor(new_unit);
            conv = old_conv/new_conv;
        end
        
        function b = is_unit_defined(obj, unit)
            % U.is_unit_defined(UNIT)
            %
            %   Returns true if UNIT can be identified by this instance,
            %   false otherwise.
            b = true;
            try
                obj.get_conv_factor(unit);
            catch err
                if ~strcmp(err.identifier, 'unit_conv:no_match')
                    rethrow(err);
                else
                    b = false;
                end
            end
        end
        
        function abbrevs = list_abbreviations(obj)
            % U.list_abbreviations()
            %
            %   Returns a cell array of unit abbreviations known to this
            %   instance.
            abbrevs = cell(size(obj.unit_names));
            for a=1:numel(abbrevs)
                abbrevs{a} = obj.unit_names{a}.abbr;
            end
        end
        
        function names = list_long_names(obj)
            % U.list_long_names()
            %
            %   Returns a cell array of unit long names known to this
            %   instance.
            names = cell(size(obj.unit_names));
            for a=1:numel(names)
                names{a} = obj.unit_names{a}.name;
            end
        end
        
        function convs = list_conversions(obj)
            % U.list_conversions()
            %
            %   Returns a vector of unit conversion factors known to this
            %   instance.
            convs = nan(size(obj.unit_names));
            for a=1:numel(convs)
                convs(a) = obj.unit_names{a}.conv;
            end
        end
    end
    
    methods(Access = private)
        function conv = get_conv_factor(obj, unit)
            % Get the conversion factor for the specified unit. Try
            % abbreviations first, then long names. Abbreviations are case
            % sensitive, long names are not. Long names may be pluralized
            % with 's' or 'es'
            
            abbr = obj.list_abbreviations();
            names = obj.list_long_names();
            unit_convs = obj.list_conversions();
            [shortex, longex] = obj.prefix_regex();
            
            if obj.abbrev_case_sensitive
                regex_fxn = @(u, r) regexp(u, r, 'match', 'once');
            else
                regex_fxn = @(u, r) regexpi(u, r, 'match', 'once');
            end
            
            unit_conv = 0;
            for a=1:numel(names)
                ex = sprintf('^(%s)?%s$', shortex, abbr{a});
                if ~isempty(regex_fxn(unit, ex));
                    if unit_conv ~= 0
                        error('unit_conv:multiple_matches', 'Multiple matches found for unit %s', unit)
                    end
                    p_ex = sprintf('^(%s)(?=%s$)', shortex, abbr{a});
                    prefix = regex_fxn(unit, p_ex);
                    if ~isempty(prefix)
                        conv = obj.short_prefixes.(prefix);
                    else
                        conv = 1;
                    end
                    conv = conv * unit_convs(a);
                    return
                end

                ex = sprintf('^(%s)?%se?s?', longex, names{a});
                if ~isempty(regexpi(unit, ex, 'once'));
                    if unit_conv ~= 0
                        error('unit_conv:multiple_matches', 'Multiple matches found for unit %s', unit)
                    end
                    p_ex = sprintf('^(%s)(?=%s)', longex, names{a});
                    prefix = lower(regexpi(unit, p_ex, 'match', 'once'));
                    if ~isempty(prefix)
                        conv = obj.long_prefixes.(prefix);
                    else
                        conv = 1;
                    end
                    conv = conv * unit_convs(a);
                    return
                end
            end
            

            error('unit_conv:no_match', 'No match found for unit %s', unit)

        end
        
        function [shortex, longex] = prefix_regex(obj)
            short_prefix_names = fieldnames(obj.short_prefixes);
            long_prefix_names = fieldnames(obj.long_prefixes);
            
            shortex = strjoin(sprintfmulti('(%s)', short_prefix_names), '|');
            longex = strjoin(sprintfmulti('(%s)', long_prefix_names), '|');
        end
    end
    
end

