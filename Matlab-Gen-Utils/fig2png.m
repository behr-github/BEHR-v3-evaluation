function [] = fig2png(file_in);
% Saves the file located at "file_in" as a .png.  Useful for remote conversion of figures.

if ~exist(file_in, 'file')
	error('fig2png:fileDNE','File %s does not exist.',file_in);
end

h = hgload(file_in);
savename = regexprep(file_in,'.png','');
saveas(h, savename, '.png');
fprintf('File saved at %s',savename);

end
