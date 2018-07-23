function [  ] = fix_soas_no2_units(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Names, ~, merge_dir] = merge_field_names('soas');
my_dir = fileparts(mfilename('fullpath'));
save_dir = fullfile(my_dir, '..', 'Workspaces', 'Scratch', 'SOAS unit fix');
if ~exist(save_dir, 'dir')
    if ask_yn(sprintf('%s doesn''t exist, create it?', save_dir));
        mkdir(save_dir);
    end
end
F = dir(fullfile(merge_dir, '*.mat'));
for a=1:numel(F)
    M = load(fullfile(merge_dir, F(a).name));
    M.Merge.Data.(Names.no2_ncar).Unit = 'ppbv';
    save(fullfile(save_dir, F(a).name), '-struct', 'M');
end
end

