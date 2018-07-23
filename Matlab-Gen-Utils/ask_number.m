function [ user_ans ] = ask_number( prompt, varargin )
%ASK_NUMBER Will ask the user to input a numeric value
%   Asks the user the question given in PROMPT. Four parameters exist:
%
%   'default' - choose a default value for the response.
%
%   'softquit' - boolean, if false (default), will exit if the user enters
%   'q' at any point by throwing an error. Unlike other ask functions, if
%   true, this will cause the function to return a NaN.
%
%   'testfxn' - a function handle to a function that will test the value
%   given by the user. 
%
%   'testmsg' - a string that will be printed if the test given by
%   'testfxn' fails. If this is not passed, it will simply say that the
%   number must fulfill the function given.
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
p.addParameter('default',[]);
p.addParameter('softquit',false,@(x) (islogical(x) && isscalar(x)));
p.addParameter('testfxn',@(x) true);
p.addParameter('testmsg','');
p.parse(varargin{:});
pout = p.Results;

testfxn = pout.testfxn;
try
    testfxn(0);
catch err
    E.badinput('Evaluation of testfxn(0) produced the error: "%s" - the function must accept a scalar number as input.',err.message)
end
if ~isa(testfxn,'function_handle') || ~isscalar(testfxn(0)) || ~islogical(testfxn(0))
    E.badinput('testfxn must be a handle to a function that returns a scalar logical value');
end
testmsg = pout.testmsg;
if ~ischar(testmsg)
    E.badinput('testmsg must be a string')
end

default = pout.default;
if isempty(default)
    use_default = false;
else
    use_default = true;
    if ~testfxn(default)
        E.badinput('default fails the given test function (%s)', func2str(testfxn))
    end
end
softquit = pout.softquit;


if use_default
    fprintf('%s (%s is default): ', prompt, sprintf_ranges(default, 'value_sep', ', '));
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
            user_ans = NaN;
            return
        else
            E.userCancel()
        end
    else
        user_ans = parse_list(strtrim(user_ans));
        try
            if any(isnan(user_ans))
                fprintf('\tNumber not recognized by str2double (q to quit): ');
            elseif ~testfxn(user_ans)
                if isempty(testmsg)
                    fprintf('\tValue must cause the function %s to return true. Try again, or q to quit: ', functiontostring(testfxn));
                else
                    fprintf('\t%s (q to quit): ', testmsg);
                end
            else
                return
            end
        catch err
            % If the test function is written expecting a scalar value,
            % then it might cause issues during the call to testfxn. This
            % will most commonly happen with the use of && or || with non
            % scalar logicals, so if that happens, print a special message.
            if strcmp(err.identifier,'MATLAB:nonLogicalConditional')
                fprintf('\tA scalar value seems to be expected. If this is not true, check the use of testfxn in ask_number\n');
                fprintf('\tTry again or q to quit: ');
            end
        end
    end
    
end

end

function list_out = parse_list(string_in)
args = strsplit(string_in, ':');

if numel(args) == 1
    list_out = str2double(strsplit(args{1}));
elseif numel(args) <= 3
    args = num2cell(str2double(args));
    list_out = colon(args{:});
else
    error('parse_input:bad_input', 'Cannot have more than 2 colons')
end

end