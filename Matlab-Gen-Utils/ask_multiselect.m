function [ user_ans ] = ask_multiselect( prompt, options, varargin )
%ASK_MULTISELECT Allows selection of one or more options from a list
%   USER_ANS = ASK_MULTISELECT( PROMPT, OPTIONS ) Displays the string
%   PROMPT followed by the list of options given in the cell array OPTIONS.
%   This returns the values chosen as a cell array.  The user will be given
%   the option to quit, if that option is taken, an error will be thrown to
%   halt execution. Options selected are returned in order and only once,
%   i.e. if the user enters 5 1 3 3 as the options, the options returned
%   will be the first, third, and fifth in that order.
%
%   The behavior can be modified by the following parameters:
%
%   'softquit' - boolean, if false (default), will exit if the user enters
%   'q' at any point by throwing an error. If true, this will cause the
%   function to return a 0 (the number not the string).
%
%   'returnindex' - boolean, if true, this function returns the indices of
%   the selections rather than the values themselves.

E = JLLErrors;
if ~ischar(prompt)
    E.badinput('prompt must be a string')
end
if ~iscellstr(options)
    E.badinput('options must be a cell array of strings')
end

p = inputParser;
p.addParameter('softquit',false,@(x) (islogical(x) && isscalar(x)));
p.addParameter('returnindex',false)
p.parse(varargin{:});
pout = p.Results;

soft_quit = pout.softquit;
return_ind = pout.returnindex;

if ~isscalar(soft_quit) || ~islogical(soft_quit)
    E.badinput('The parameter ''softquit'' must be a scalar logical value');
end
if ~isscalar(return_ind) || ~islogical(return_ind) && ~isnumeric(return_ind)
    E.badinput('The parameter ''returnindex'' must be a scalar logical value');
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s:\n',prompt);
fprintf('(You may select a range with colon notation, e.g. 1:5)\n')
for a=1:numel(options)
    fprintf('\t%d - %s\n',a,options{a});
end

testfxn = @(x) all(x >= 1 & x <= numel(options));
testmsg = sprintf('Only values 1-%d permitted', numel(options));
try
    sel_inds = ask_number(sprintf('Select options (1-%d)', numel(options)), 'softquit', soft_quit, 'testfxn', testfxn, 'testmsg', testmsg);
catch err
    if strcmp(err.identifier,'ask_number:user_cancel')
        E.userCancel();
    else
        rethrow(err);
    end
end

% Return options in order
sel_inds = sort(unique(sel_inds));
if isnan(sel_inds)
    user_ans = 0;
    return
elseif return_ind
    user_ans = sel_inds;
else
    user_ans = options(sel_inds);
end
end

