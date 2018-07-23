function [given_val] = opt_ask_number(prompt, given_val, value_name, varargin)
%OPT_ASK_NUMBER Function that asks the user to choose from a list of choices if a value is not given
%   GIVEN_VAL = OPT_ASK_NUMBER( PROMPT, GIVEN_VAL, VALUE_NAME ) Will check
%   if the given value (GIVEN_VAL) is empty and if so, ask the user to
%   choose a valid option with ASK_NUMBER( PROMPT, varargin{:}). If not
%   empty, checks that GIVEN_VAL returns true from the function given by
%   the 'testfxn' parameter, which is also passed through to ask_number if
%   necessary. The 'testmsg' parameter is used here and also by ask_number.
%   This default behavior assumes that you will ingest options using the
%   input parser and make each value empty if it is not given.
%
%   Parameters:
%
%       'ask_condition' - a function handle that accepts one argument
%       (GIVEN_VAL) and should return true if the user should be asked to
%       choose a value. By default this is "@isempty"
%
%       'testfxn' and 'testmsg' - a function handle that given in GIVEN_VAL
%       should return true if it is a valid value. If not, an error is
%       thrown with the message given by testmsg. By default, 'testfxn' is
%       set to "@(x) isnumeric(x) && isscalar(x)" and 'testmsg' is will say
%       that the VALUE_NAME must be numeric and scalar.

p = inputParser;
p.addParameter('ask_condition', @isempty);
p.addParameter('testfxn', @(x) isnumeric(x) && isscalar(x));
p.addParameter('testmsg', sprintf('%s must be a scalar number', value_name));

p.KeepUnmatched = true;

p.parse(varargin{:});
pout = p.Results;

ask_condition = pout.ask_condition;
given_test_fxn = pout.testfxn;
given_test_msg = pout.testmsg;

if ask_condition(given_val)
    given_val = ask_number(prompt, varargin{:});
elseif ~given_test_fxn(given_val)
    stack = dbstack(1,'-completenames');
    errstruct = struct('identifier','opt_ask_multichoice:bad_input','message',given_test_msg,'stack',stack);
    error(errstruct);
end
end

