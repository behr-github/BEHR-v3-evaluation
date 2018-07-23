function [ strout ] = tex_in_printf( strin, n_printf )
%TEX_IN_PRINTF Escapes backslashes not part of special characters 
%   printf functions (sprintf, fprintf) read special characters like \n
%   (new line) or \b (backspace) which is very useful for special
%   formatting but makes using Tex characters (e.g. \in, \Delta) difficult
%   because the backslashes there need to be escaped. This is especially
%   annoying if the string will be passed through multiple printf
%   functions, as the number of escapes needed increases geometrically.
%
%   This function will simplify this somewhat by escaping all bashslashes
%   not part of an identifiable special characters (listed at
%   http://www.mathworks.com/help/matlab/ref/sprintf.html#inputs). The
%   limiting factor is that in most cases the special character must be
%   followed by a whitespace or another backslash to be identified as such,
%   and not part of a Tex string. Consider two strings:
%       'Wind velocity \nu = 8 m/s'
%       'Frequency \nu = 8 MHz'
%   In the first, the user intended the \n to represent a newline, with 'u
%   = 8 m/s' on the second line, whereas in the second, '\nu' was intended
%   to be interpreted as the Greek letter nu. Since listing every possible
%   Tex string that starts with the same letter as the special characters
%   would be painful and prone to failure, this function assumes that any
%   special characters should be followed by a space or a backslash. So
%   'Long first line \n\n Even longer second line' will correctly interpret
%   the '\n\n ' as two new lines.
%
%   The second argument is optional (it defaults to 1). This handles the
%   case where you will pass the string through several printf functions,
%   which requires the backslashes to be escaped multiple times. Example:
%
%       mystr = sprintf('Values of \Delta E in cells 1-%d: %s', ncells, repmat('%f ', 1, ncells))
%       newstr = sprintf(mystr, cells{:})
%
%   would require \Delta E to be written as \\\\Delta E, since the first
%   sprintf will read \\\\Delta E to \\Delta E, leaving the backslash
%   escaped for the second sprintf.
%
%   Josh Laughner <joshlaugh5@gmail.com> 9 Jul 2015

% The regular expression that will identify all instances of a backslash
% not part of the special characters \a, \b, \f, \n, \r, \t, \v, \xN, \N
% (where N is some number) followed by a space or end of string
%
% The '(?! )' notation is a look-ahead assertion that requires the
% backslash not be followed by any of the patterns listed but does not
% return those patterns as part of the string to be replaced.

regex_pat = '\\(?!a\s|b\s|f\s|n\s|r\s|t\s|v\s|x\d+\s|\d+\s|a$|b$|f$|n$|r$|t$|v$|x\d+$|\d+$)';
%regex_pat = '\\(?!a\s|b\s|f\s|n\s|r\s|t\s|v\s)';

% This will figure out how many times the backslashs must be escaped. We
% need to add one because the regexprep itself needs them escaped
if nargin < 2
    n_printf = 1;
else
    if ~isscalar(n_printf) || ~isnumeric(n_printf)
        error('tex_in_printf:bad_input','The optional input n_printf must be a scalar number');
    end
end
n = 2^(n_printf+1);
rep_pat = repmat('\',1,n);

% Do the replacement
strout = regexprep(strin, regex_pat, rep_pat);

end

