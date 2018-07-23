function [  ] = print_eps( fig, filename )
%PRINT_EPS Saves a specified figures as an EPS image
%   PRINT_EPS( FIG, FILENAME ) Save the figure given by the handle FIG to
%   FILENAME; FILENAME should not include the extension.

E = JLLErrors;

if ~isgraphics(fig, 'figure')
    E.badinput('FIG must be a handle to a figure');
end

if ~ischar(filename)
    E.badinput('FILENAME must be a string')
end
% Test if an extension was given, and if the path exists
[p,f] = fileparts(filename);
eps_fname = fullfile(p,f);
if ~strcmp(eps_fname, filename)
    E.badinput('FILENAME must not include an extension');
elseif ~exist(p, 'dir') && ~isempty(p) % if p is empty, then we are saving in the current directory
    E.badinput('Directory %s does not exist', p);
end

% This ensures that any resizing of the figure is represented in the EPS
% file.
set(fig, 'paperpositionmode', 'auto');
print(fig, '-depsc2', '-loose', filename);

end

