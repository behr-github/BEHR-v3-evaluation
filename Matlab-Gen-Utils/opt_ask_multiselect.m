function given_value = opt_ask_multiselect(prompt, allowed_values, given_value, value_name, varargin)
%OPT_ASK_MULTISELECT Function that asks the user to choose from a list of choices if a value is not given
%   GIVEN_VAL = OPT_ASK_MULTISELECT( PROMPT, ALLOWED_VALUES, GIVEN_VAL,
%   VALUE_NAME ) Will check if the given value (GIVEN_VAL) is empty and if
%   so, ask the user to choose a valid option with ASK_MULTISELECT( PROMPT,
%   ALLOWED_VALUES, varargin{:}). If not empty, checks that GIVEN_VAL is a
%   char array and a member of ALLOWED_VALUES and if not, throws an error,
%   using VALUE_NAME to describe the value.
%
%   This default behavior assumes that you will ingest options using the
%   input parser and make each value empty if it is not given.
%
%   Parameters:
%
%       'ask_condition' - a function handle that accepts one argument
%       (GIVEN_VAL) and should return true if the user should be asked to
%       choose a value. By default this is "@isempty"
%
%       'test_fxn' - a function handle that accepts two arguments
%       (GIVEN_VAL and ALLOWED_VALUES) and returns true if GIVEN_VAL is a
%       valid option. By default this is "@(val, allowed) iscell(val) &&
%       all(ismember(val, allowed))".

E = JLLErrors;
p = inputParser;
p.addParameter('ask_condition', @isempty);
p.addParameter('test_fxn', @(val, allowed) iscell(val) && all(ismember(val, allowed)));
p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;
ask_condition = pout.ask_condition;
test_fxn = pout.test_fxn;

my_params = fieldnames(pout);
remaining_params = update_params('remove', varargin, my_params{:});

if ask_condition(given_value)
    given_value = ask_multiselect(prompt, allowed_values, remaining_params{:});
elseif ~test_fxn(given_value, allowed_values)
    msg = sprintf('%s must be a cell array containing one or more of: "%s"', value_name, strjoin(allowed_values, '", "'));
    stack = dbstack(1,'-completenames');
    errstruct = struct('identifier','opt_ask_multichoice:bad_input','message',msg,'stack',stack);
    error(errstruct);
end

end

