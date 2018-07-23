function [ carray ] = wrap_string( str, width )
%CARRAY = WRAP_STRING( STR ) Wraps STR to 80 character wide character array
%
%CARRAY = WRAP_STRING( STR, WIDTH ) Wraps STR to a WIDTH character wide
%character array.

E = JLLErrors;

if ~ischar(str)
    E.badinput('STR must be a string')
end

if ~exist('width','var')
    width = 80;
elseif ~isscalar(width) || ~isnumeric(width) || width < 1
    E.badinput('WIDTH must be a scalar number > 0');
end


% Make a character array into a single string; joining rows with spaces.
if ~isvector(str)
    str = strjoin(cellstr(str),' ');
end
l = length(str);
carray = [];

a = 1;
while a+width < l+1
    % Get the next WIDTH characters and find the last space in the line. If
    % no space in the line, will break in the middle of the word.
    line = str(a:a+width-1);
    xx = find(line == ' ',1,'last');
    if ~isempty(xx)
        line = strtrim(line(1:xx));
        if ~isempty(line) % this test should handle cases where the line's only space is at the beginning.
            spacer = repmat(' ',1,width-length(line));
            line = [line, spacer]; %#ok<AGROW>
            carray = cat(1, carray, line);
        end
        a = a+xx;
    else
        carray = cat(1, carray, line);
        a = a+width;
    end
end
% The last line with most likely need a spacer.
line = str(a:end);
spacer = repmat(' ',1,width-length(line));
line = [line, spacer]; 
carray = cat(1, carray, line);

% for a=1:width:l
%     b = a+width-1; 
%     if b <= l
%         carray = cat(1, carray, str(a:b));
%     else
%         line = str(a:end);
%         spacer = repmat(' ',1,width - length(line));
%         line = [line, spacer];
%         carray = cat(1, carray, line);
%     end
%         
% end


end

