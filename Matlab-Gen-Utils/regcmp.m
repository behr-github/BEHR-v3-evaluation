function [ b ] = regcmp( str, expression, varargin )
%REGCMP Compares a string against a regular expression
%   B = REGCMP( STR, EXPRESSION, ... ) Will compare STR against the case
%   sensitive regular expression EXPRESSION and return B = true if a match
%   is found.

E = JLLErrors;

if ~ischar(expression)
    E.badinput('EXPRESSION must be a string')
end

if ischar(str)
    b = ~isempty(regexp(str, expression, varargin{:}, 'once'));
elseif iscellstr(str)
    tmp = regexp(str, expression, varargin{:}, 'once');
    b = ~cellfun(@isempty, tmp);
else
    E.badinput('STR must be a string or cell array of strings');
end

end

