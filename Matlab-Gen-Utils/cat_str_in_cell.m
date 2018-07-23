function [ str_out ] = cat_str_in_cell( cell_in )
%cat_str_in_cell: Concatenates string parts in a cell array.

n = numel(cell_in);
str_out = '';

for a=1:n;
    str_out = strcat(str_out, sprintf(' %s',cell_in{a}'));
end


end

