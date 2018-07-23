function [ user_ans ] = ask_multichoice( prompt, allowed_options, varargin )
%ASK_MULTICHOICE Will ask the user to select from the allowed options
%   Asks the user the question given in PROMPT to choose from the cell
%   array ALLOWED_OPTIONS. Two parameters exist:
%
%   'default' - choose a default value for the response.
%
%   'list' - boolean, if false (default), will require the user to type
%   their response. If true, will present a list of options and ask the
%   user to enter the number corresponding to that option.
%
%   'index' - boolean, if false (default), will return the actual value as
%   a string. If true, returns the index of the choice instead.
%
%   'softquit' - boolean, if false (default), will exit if the user enters
%   'q' at any point by throwing an error. If true, this will cause the
%   function to return a 0 (the number not the string).
%
%   This function always returns a string (unless softquit or index is
%   true) so you'll need to parse the string back into a number if you want
%   it so. Answers are also always returned in lower case.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 Jan 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if ~ischar(prompt)
    E.badinput('prompt must be a string')
end
if ~iscellstr(allowed_options)
    E.badinput('allowed_options must be a cell array of strings')
end

p = inputParser;
p.addParameter('default',0);
p.addParameter('list',false);
p.addParameter('index',false);
p.addParameter('softquit',false,@(x) (islogical(x) && isscalar(x)));
p.parse(varargin{:});
pout = p.Results;

default = pout.default;
if default == 0
    use_default = false;
else
    if ~ischar(default) || ~ismember(default, allowed_options)
        E.badinput('default must be a string that is one of those in allowed_options')
    end
    use_default = true;
end
softquit = pout.softquit;
print_list = pout.list;
return_ind = pout.index;

if ~isscalar(print_list) || (~islogical(print_list) && ~isnumeric(print_list))
    E.badinput('The parameter ''list'' must be a scalar logical or number')
end

if print_list
    fprintf('%s:\n',prompt);
    for a=1:numel(allowed_options)
        fprintf('\t%d - %s',a,allowed_options{a});
        if use_default && strcmpi(allowed_options{a},default)
            fprintf(' (default)');
        end
        fprintf('\n');
    end
    fprintf('Enter 1 to %d or q to quit: ',numel(allowed_options));
else
    if use_default
        fprintf('%s (%s - %s is default): ', prompt, strjoin(allowed_options, ', '), default);
    else
        fprintf('%s (%s): ', prompt, strjoin(allowed_options, ', '));
    end
end


while true
    user_ans = lower(input('', 's'));
    if print_list
        user_ind = str2double(user_ans);
    end
    if use_default && isempty(user_ans)
        if return_ind
            user_ans = find(strcmp(default, allowed_options));
        else
            user_ans = default;
        end
        return
    elseif print_list && user_ind >= 1 && user_ind <= numel(allowed_options)
        if return_ind
            user_ans = user_ind;
        else
            user_ans = allowed_options{user_ind};
        end
        return
    elseif ~print_list && ismember(user_ans, lower(allowed_options))
        if return_ind
            user_ans = find(strcmp(user_ans, allowed_options));
        end
        return
    elseif strcmpi(user_ans, 'q')
        if softquit
            user_ans = 0;
        else
            E.userCancel()
        end
    else
        if print_list
            fprintf('You must enter 1 to %d or q to quit: ',numel(allowed_options));
        else
            fprintf('You must choose one of %s, or enter q to quit: ', strjoin(allowed_options, ', '));
        end
    end
    
end

end

