function attr_val = ncreadatt_default(filename, variable, attribute, varargin)
% NCREADATT_DEFAULT Read a netCDF attribute with a default fallback if not found
%   ATTR_VAL = NCREADATT( FILENAME, VARIABLE, ATTRIBUTE, DEFAULT ) Reads
%   the ATTRIBUTE of VARIABLE in file given by FILENAME. If the attribute
%   is not found, returns the value DEFAULT instead. (Internally, this
%   calls NCREADATT(FILENAME, VARIABLE, ATTRIBUTE) and returns DEFAULT if
%   an error occurs.) By default, leading and trailing whitespace is
%   trimmed from the attribute value.
%
%   NCREADATT_DEFAULT( FILENAME, VARIABLE, ATTRIBUTE ) This behaves as
%   NCREADATT and will throw an error if ATTRIBUTE cannot be found.
%   Whitespace is still trimmed. This is a convenience to allow this
%   function to replace all uses of NCREADATT even where a default value is
%   not desired.
%
%   NCREADATT_DEFAULT( ___, 'keep_whitespace' ) will prevent the function
%   from trimming leading and trailing whitespace.
%
%   Additional parameters:
%       'fatal_if_missing' - scalar logical value (default is false). If
%       set to true, this will override the default behavior and throw an
%       error if the attribute is not found whether or not DEFAULT is
%       given. This gives you an easier way to switch off default values in
%       a broader netCDF reading function without needing to wrap every
%       call to this in an if statement to not pass the DEFAULT value if
%       errors are desired.
%
%       'DEBUG_LEVEL' - level of verbosity in printing to the console. Must
%       be a scalar number, default is 2. 0 is no printing to the console.

E = JLLErrors;

NO_DEFAULT_VALUE = 'NO_DEFAULT_VALUE_ERROR_IF_MISSING';

p = advInputParser;
p.addOptional('default', NO_DEFAULT_VALUE, @ischar);
p.addParameter('fatal_if_missing', false);
p.addFlag('keep_whitespace');
p.addParameter('DEBUG_LEVEL',2);

p.parse(varargin{:});
pout = p.AdvResults;

default = pout.default;
fatal_if_missing = pout.fatal_if_missing;
keep_whitespace = pout.keep_whitespace;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

%%%%%%%%%%%%%%%%%%
% Input checking %
%%%%%%%%%%%%%%%%%%

% If no default given, an error must be thrown if the attribute is not
% found. Otherwise, put no limits on the default value. 
if ischar(default) && strcmp(default, NO_DEFAULT_VALUE)
    fatal_if_missing = true;
end

if ~isscalar(fatal_if_missing) || (~islogical(fatal_if_missing) && ~isnumeric(fatal_if_missing))
    E.badinput('The parameter "fatal_if_missing" must be a scalar logical or numeric value')
end

if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('The parameter "DEBUG_LEVEL" must be a scalar number');
end

%%%%%%%%%%%%%%%%%
% Main function %
%%%%%%%%%%%%%%%%%

try
    % Read the value. Remove leading and trailing whitespace because CMAQ
    % likes to add a bunch of that.
    attr_val = ncreadatt(filename, variable, attribute);
catch attr_err
    % If we cannot find the attribute in the file fall back on the default
    if strcmp(attr_err.identifier, 'MATLAB:imagesci:netcdf:libraryFailure') && ~fatal_if_missing
        if DEBUG_LEVEL > 0
            % Convert default to a string to insert in the warning message.
            % Handle common types.
            if isnumeric(default)
                default_str = num2str(default);
            else
                default_str = default;
            end
        
            fprintf('\tAssuming "%s" for attribute "%s" not found in variable "%s" of file "%s".\n', default_str, attribute, variable, filename); 
        end
        attr_val = default;
    else
        rethrow(attr_err);
    end
end

if ~keep_whitespace
    attr_val = strtrim(attr_val);
end

end