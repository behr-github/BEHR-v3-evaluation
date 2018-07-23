function [  ] = save_individual_geos_vars( struct_in, savepath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

struct_name = inputname(1);

for a=1:numel(struct_in)
    vname = sprintf('%s_%s', struct_name, regexprep(struct_in(a).fullName,'[-_$\ ]','')); % remove special characters
    vname = regexprep(vname,'\(\w*\)',''); % remove anything in parentheses
    eval(sprintf('%s = struct_in(a);',vname));
    fname = sprintf('%s-%s.mat', struct_name, regexprep(struct_in(a).fullName,'-\$',''));% remove special characters
    fname = regexprep(fname,'\(.*\)',''); % remove anything in parentheses
    fname = regexprep(fname,' *.mat','.mat'); % trim space before the file extension
    fprintf('Saving %s\n',fname);
    save(fullfile(savepath,fname),vname);
    clear(vname);
end

end

