function ret = savefigs(allflag)
% SAVEFIGS Quickly save all open figures
%
% This function allows you to quickly save all currently open figures with
% a custom filename for each in multiple formats.  To use the function
% simply call savefigs with no arguments, then follow the prompts
%
% Upon execution this function will first ask what type of file to save by
% extension. Accepted extensions are:
%       .fig
%       .cmp
%       .emf
%       .eps
%       .jpg
%       .pbm
%       .pcx
%       .pdf
%       .pgm
%       .png
%       .ppm
%       .tif
%
% It will then one-by-one bring each currently open figure to the foreground
% and supply a text prompt in the main console window asking you for a
% filename. If you enter a blank string, it will save the figure with its
% title (after some reformatting to make it path safe, see below) and
% todays date.
%
% Some characters will be removed or changed to make titles path safe:
%       / replaced by -
%       : replaced by -
%       \ removed
%       % replaced by "percent"
%
% 
% SAVEFIGS ALL will ask for an extension, but then will save every open
% figure in the current directory using its title and todays date as the
% filename.
%
% Copyright 2010 Matthew Guidry 
% matt.guidry ATT gmail DOTT com  (Email reformatted for anti-spam)
% Modified 2016 Josh Laughner


if exist('allflag','var') && strcmpi(allflag,'all')
    all_bool = true;
else
    all_bool = false;
end

extension = input('Enter extension to use (blank for .fig). Include the dot.  ','s');

if strcmp(extension,'')
    extension = '.fig';
end

%Check that the input is a valid extension to save figures
if (strcmp(extension,'.cmp') || strcmp(extension,'.emf') || strcmp(extension,'.eps') || strcmp(extension,'.fig') || strcmp(extension,'.jpg') || strcmp(extension,'.pbm') || strcmp(extension,'.pcx') || strcmp(extension,'.pdf') || strcmp(extension,'.pgm') || strcmp(extension,'.png') || strcmp(extension,'.ppm') || strcmp(extension,'.tif'))
fprintf('Extension %s accepted.\n',extension);
pause
hfigs = get(0, 'children');                          %Get list of figures

if all_bool
    extra_name = input('If you wish, enter something to append to each filename: ', 's');
end

for m = 1:length(hfigs)
    figure(hfigs(m));                                %Bring Figure to foreground
    if ~all_bool
        filename = input('Filename? (0 to skip or blank to use title with today''s date)  ', 's');%Prompt user
    end    
    if ~all_bool && strcmp(filename, '0')                        %Skip figure when user types 0
        continue
    elseif all_bool || strcmp(filename,'')
        ax = findall(gcf,'Type','axes');
        for i=1:numel(ax)
            htitle = get(ax(i),'Title');
            filename = strtrim(get(htitle,'String'));
            if iscell(filename); filename = cat_str_in_cell(filename); end
            if all_bool && ~isempty(extra_name)
                filename = [filename, ' ',extra_name];
            end
            if ~isempty(filename); break; end
        end
        filename = strrep(filename,'/','-');
        filename = strrep(filename,':','-');
        filename = strrep(filename,'\','');
        % Replace % signs with the word percent with reasonable
        % capitalization
        filename = regexprep(filename,'(?<=\w\s+)\%','percent'); % preceded by a letter before a whitespace - lowercase percent
        filename = regexprep(filename, '%', 'Percent'); % otherwise capital
        filename = strtrim(filename);
        filename = sprintf('%s - %s',datestr(today,29),filename);
        fprintf('Saving figure %u as %s.\n',m,filename);
        saveas(hfigs(m), [filename, extension]);
    else
        %saveas(hfigs(m), [filename '.fig']) %Matlab .FIG file
        %saveas(hfigs(m), [filename '.emf']) %Windows Enhanced Meta-File (best for powerpoints)
        saveas(hfigs(m), [filename, extension]); %Standard PNG graphics file (best for web)
        %eval(['print -depsc2 ' filename])   %Enhanced Postscript (Level 2 color) (Best for LaTeX documents)
    end
end
else
    error('extension:invalid','Extension %s not recognized. Be sure to include the ''.''.', extension);
end
