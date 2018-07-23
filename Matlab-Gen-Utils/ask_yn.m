function [ user_bool ] = ask_yn( prompt )
%ASK_YN Ask a yes or no question
%   USER_BOOL = ASK_YN( PROMPT ) Ask the question PROMPT and accept a yes
%   or no answer. Returns true for yes, false for no.

user_bool = strcmpi(ask_multichoice(prompt, {'y','n'}),'y');

end

