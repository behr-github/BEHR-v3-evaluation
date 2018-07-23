function [ str ] = capitalize_words( str )
%CAPITALIZE_WORDS Capitalize each word in a string
%   STR = CAPITALIZE_WORDS( STR ) Capitalzed each word in STR, i.e. any
%   letter preceeded by at least one whitespace is capitalized. The first
%   character in the string is also capitalized.

% Adapted from https://www.mathworks.com/matlabcentral/answers/107307-function-to-capitalize-first-letter-in-each-word-in-string-but-forces-all-other-letters-to-be-lowerc

% Find all instances of a non-whitespace character preceeded by whitespace.
% Pad the beginning of the string with a space to include the first word,
% then decrement the indices returned to adjust for that.
idx = regexp([' ' str],'(?<=\s+)\S','start')-1;
str(idx) = upper(str(idx));

end

