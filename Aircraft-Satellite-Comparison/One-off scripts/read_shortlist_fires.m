%Read_short_fire_list
fname = '';
fid = fopen(fname,'r');
line = 1;
i=1;
while line~=-1;
    line = fgetl(fid);
    line = regexprep(line,'<\td>',',');
    line = regexprep(line,'<[^>]*>','');
    readin = textscan(line,'%s %s %s %d %d','Delimiter',',');
end