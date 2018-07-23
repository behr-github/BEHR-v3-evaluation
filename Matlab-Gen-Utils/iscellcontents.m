function [ chk ] = iscellcontents( cell_in, test_fxn )
%iscellcontents Checks if cell contents are of the given type
%   Tired of being unable to check what a cell contains automatically? This
%   is the function for you! Takes two arguments: a cell array and a string
%   of the test function name you wish to apply to each cell.  For example,
%   passing 'ischar' as the second argument will check if each cell
%   contains characters.  This function returns a logical matrix the same
%   size as the cell input.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,2);
if ~iscell(cell_in)
    E.badinput('"cell_in" must be a cell array');
elseif ~ischar(test_fxn) && ~isa(test_fxn,'function_handle')
    E.badinput('"test_fxn" must be a string or function handle');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

chk = false(size(cell_in));
if ischar(test_fxn)
    test_hndl = str2func(test_fxn);
else
    test_hndl = test_fxn;
end

for a=1:numel(cell_in)
    try
        this_chk = test_hndl(cell_in{a});
    catch err
        % Handle the case where the user passed an invalid function name
        % specially
        if strcmp(err.identifier,'MATLAB:UndefinedFunction')
            E.badinput('The function "%s" does not appear to be a valid function.  Check the spelling. If it is a build-in function, check your toolboxes. If it is a custom function, check that it is on the search path.',test_fxn);
        else
            rethrow(err);
        end
    end
    
    % Check that this_chk makes sense as what a logical test should return:
    % it should be a scalar and have a value of 1 or 0.
    if isscalar(this_chk)
        if this_chk ~= 1 && this_chk ~= 0
            E.callError('logical_test_failure',sprintf('The output to %s should be 1 or 0 - only logical test functions should be used',test_fxn));
        end
    else
        E.callError('logical_test_failure',sprintf('The output to the test function should be scalar - iscellcontents is meant to test the type of variable in each cell array'));
    end
    
    % If the output of the test function makes sense, save it to the output
    % variable
    chk(a) = this_chk;
end
    
end

