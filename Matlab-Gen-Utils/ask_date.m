function [ user_ans ] = ask_date( prompt, varargin )
%ASK_DATE Will ask the user to provide a date.
%   Asks the user the question given in PROMPT to provide a date input. If
%   the input date string isn't parsable by Matlab, the user will be
%   prompted to reenter it. Two parameters exist:
%
%   'default' - choose a default value for the response.
%
%   'softquit' - boolean, if false (default), will exit if the user enters
%   'q' at any point by throwing an error. If true, this will cause the
%   function to return a 0 (the number not the string).
%
%   This function always returns the string given by the user (unless
%   softquit is true).
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 Jan 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if ~ischar(prompt)
    E.badinput('prompt must be a string')
end

p = inputParser;
p.addParameter('default',0);
p.addParameter('softquit',false,@(x) (islogical(x) && isscalar(x)));
p.parse(varargin{:});
pout = p.Results;

default = pout.default;
if default == 0
    use_default = false;
else
    use_default = true;
    try
        datenum(default);
        datestr(default);
    catch err
        if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString') || strcmp(err.identifier, 'MATLAB:datestr:ConvertDateNumber')
            E.badinput('default must be a valid date string or date number')
        else 
            rethrow(err);
        end
    end
end
softquit = pout.softquit;

if use_default
    fprintf('%s (%s is default): ', prompt, default);
else
    fprintf('%s: ', prompt);
end

while true
    user_ans = lower(input('', 's'));
    if use_default && isempty(user_ans)
        user_ans = default;
        return
    elseif strcmpi(user_ans, 'q')
        if softquit
            user_ans = 0;
            return
        else
            E.userCancel()
        end
    else
        try
            datenum(user_ans);
            return
        catch err
            if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
                fprintf('\tFormat not recognized. Try yyyy-mm-dd or \n\tyyyy-mm-dd HH:MM:SS if a time is required, or q to quit: ');
            else
                rethrow(err);
            end
        end
    end
    
end

end

