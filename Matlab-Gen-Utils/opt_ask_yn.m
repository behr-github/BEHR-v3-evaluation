function given_value = opt_ask_yn(prompt, given_value, value_name, varargin)
% OPT_ASK_YN Ask the user a yes/no question if a value isn't given
%   GIVEN_VALUE = OPT_ASK_YN(PROMPT, GIVEN_VALUE, VALUE_NAME) Takes
%   the current value, GIVEN_VALUE, and if it is a NaN, asks the user
%   a yes/no question given by PROMPT. If GIVEN_VALUE is not a NaN,
%   then this function tests if it is a scalar logical and throws an
%   error in the calling function if not, using the name of the value
%   VALUE_NAME in the error message.
%
%   Parameters:
%
%       'ask_condition' - function handle that takes one input (the given
%       value) and must return scalar true or false. If true, the prompt
%       will be asked interactively. Default is "@isnan".
%
%       'test_fxn' - function handle that accepts one argument (the given
%       value) and returns true if it is a valid value. Default is 
%       "@(val) islogical(val) && isscalar(val)".
%
%       'test_msg' - the error message that will be printed if the test_fxn
%       returns false. Is given value_name as its only format value. Default
%       is "%s must be a scalar logical".

E = JLLErrors;
p = inputParser;
p.addParameter('ask_condition', @isnan);
p.addParameter('test_fxn', @(val) islogical(val) && isscalar(val));
p.addParameter('test_msg', '%s must be a scalar logical');
p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;
ask_condition = pout.ask_condition;
test_fxn = pout.test_fxn;
test_msg = pout.test_msg;

my_params = fieldnames(pout);
remaining_params = update_params('remove', varargin, my_params{:});

if ask_condition(given_value)
    given_value = ask_yn(prompt, remaining_params{:});
elseif ~test_fxn(given_value)
    msg = sprintf(test_msg, value_name);
    stack = dbstack(1,'-completenames');
    errstruct = struct('identifier','opt_ask_yn:bad_input','message',msg,'stack',stack);
    error(errstruct);
end

end

