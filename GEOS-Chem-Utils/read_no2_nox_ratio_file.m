function [ dc3_pres, dc3_ratio, gc_pres, gc_ratio ] = read_no2_nox_ratio_file( ratio_file )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

fid = fopen(ratio_file, 'r');
if fid < 0
    E.badinput('Cannot open %s', ratio_file);
end

tline = fgetl(fid); % header
tline = fgetl(fid); % first time through gets the first data line
dc3_ratio = [];
dc3_pres = [];
gc_ratio = [];
gc_pres = [];
while ischar(tline)
    data = strsplit(tline);
    data = data(~iscellcontents(data, 'isempty'));
    if numel(data) == 2
        gc_ratio = [gc_ratio; str2double(data{2})];
        gc_pres = [gc_pres; str2double(data{1})];
    elseif numel(data) == 4
        dc3_pres = [dc3_pres; str2double(data{1})];
        dc3_ratio = [dc3_ratio; str2double(data{4})];
        gc_pres = [gc_pres; str2double(data{2})];
        gc_ratio = [gc_ratio; str2double(data{3})];
    else
        fprintf('Skipping line, not 2 or 4 entries\n');
    end
    tline = fgetl(fid);
end

end

