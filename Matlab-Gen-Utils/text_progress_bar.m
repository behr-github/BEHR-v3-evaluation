function [  ] = text_progress_bar( i, n )
%TEXT_PROGRESS_BAR Print a text-based progress bar
%   TEXT_PROGRESS_BAR( I, N ) Typical use is in a for loop from 1:n to call
%   TEXT_PROGRESS_BAR(I, N) where I is the for loop index. This will print
%   a progress bar like
%
%   [====>      ]
%
%   where the number of spaces between the [] is N+1 and the length of the
%   arrow increases with I. Note that this will behave strangely is any
%   other text is printed between calls to this function.

E = JLLErrors;

if ~isnumeric(n) || ~isscalar(n)
    E.badinput('N must be a scalar number');
end
if ~isnumeric(i) || ~isscalar(i) || i < 1 || i > n
    E.badinput('I must be a scalar number that is between 1 and N');
end

progress_str = [repmat('=',1,i), '>', repmat(' ',1,n-i), ']'];
delete_str = repmat('\b',1,numel(progress_str));

if i > 1
    fprintf(delete_str);
else
    % Printing this here and then only printing the progress bar after it
    % fixed an issue where it would erase a previous line. It may be that I
    % just accidentally had the number of backspaces wrong before, but this
    % works, so not going to mess with it.
    fprintf('[');
end
fprintf(progress_str);
if i == n
    fprintf('\n');
end


end

