function [  ] = classhelp( class_in, method )
%CLASSHELP Prints help on a class or methods therein
%   The matlab help functions don't seem to handle methods in classes very
%   well. This function will print the help comments from a class and list
%   every method in it. Additionally you can specify a method name for
%   more help on that method. Pass 'allmethods' as the method name to list
%   every non-private method in the class with its help.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(class_in)
    E.badinput('Must pass class name as a string')
elseif nargin > 1 && ~ischar(method)
    E.badinput('Must pass the method name as a string')
elseif nargin < 2
    method = 'classonly';
end

C = which(class_in);
if isempty(C)
    fprintf('%s is not on your MATLAB path\n', class_in);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(method,'classonly')
    % Print the class help, then each public method but not the help
    % for any method
    help(class_in)
    methods(class_in)
    return
end

fid = fopen(C);
tline = fgetl(fid);

cleanobj = onCleanup(@() fclose(fid));

public_method_block = false;

while ischar(tline)
    if ~isempty(strfind(tline,'methods')) && isempty(regexpi(tline,'Access\o{40}*=\o{40}*private')) && isempty(regexpi(tline,'Access\o{40}*=\o{40}*protected'))
        public_method_block = true;
    elseif ~isempty(strfind(tline,'methods'))
        public_method_block = false;
    end
    
    if public_method_block
        if ~isempty(strfind(tline,'function'))
            [s,e] = regexp(tline,'(?<=function\o{40}.*\o{40})\w*(?=\()');
            fxn_name = tline(s:e);
            if ~isempty(s) && (strcmpi(method,'allmethods') || strcmpi(fxn_name, method))
                tline = strtrim(tline);
                i = find(tline ~= '%', 1); % this will be used to remove the % symbols beginning the line
                tline = tline(i:end);
                fprintf('%s\n',tline);
                help(sprintf('%s>%s',class_in,fxn_name))
            end
        end
    end
    tline = fgetl(fid);
end

%     function print_help
%         tline = fgets(fid);
%         while ischar(tline)
%             tline = strtrim(tline);
%             if ~strcmp(tline(1),'%')
%                 return
%             end
%             i = find(tline ~= '%', 1); % this will be used to remove the % symbols beginning the line
%             tline = tline(i:end);
%             fprintf(tline);
%             
%             tline = fgets(fid);
%         end
%     end

end
