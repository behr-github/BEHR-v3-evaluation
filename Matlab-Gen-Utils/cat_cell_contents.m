function array_out = cat_cell_contents(cell_in)
% CAT_CELL_CONTENTS Concatenate numeric arrays from a cell array
%
%   A = CAT_CELL_CONTENTS( C ) Given a cell array, C, of all numeric 
%   arrays of differing sizes, will concatenate them along a new 
%   dimension, padding as necessary to make them all the same size.
%   That is, if the largest array is 3x4x5, then all arrays in C will
%   be concatenated along the fourth dimension, and padded with NaNs
%   to be 3x4x5.

all_numeric = cellfun(@isnumeric, cell_in);


if ~all(all_numeric)
    E.badinput('All arrays in the input cell array must be numeric')
end

% Get how big the concatenated array will have to be
n_cells = numel(cell_in);
max_ndims = max(cellfun(@ndims, cell_in));
final_size = zeros(1, max_ndims);
for i_arr = 1:n_cells
    sz_i = padded_size(cell_in{i_arr}, max_ndims);
    final_size = max(final_size, sz_i);
end

array_out = nan([n_cells, final_size]);

% Pad the individual arrays to concatenate them, then insert in the output
% array
for i_arr = 1:n_cells
    array_i = cell_in{i_arr};
    pad_size = final_size - padded_size(array_i, max_ndims);
    % Insert the padded array into the final array. Need to reshape to a
    % vector during the insertion, but since it and the final array slice
    % are the same shape (after padding), it's shape will be retained.
    array_out(i_arr, :) = reshape(padarray(array_i, pad_size, nan, 'post'),[],1);
end

perm_vec = [2:(max_ndims+1), 1];
array_out = permute(array_out, perm_vec);

end

function sz = padded_size(array, req_ndims)
sz = size(array);
sz = padarray(sz, [0, req_ndims - ndims(array)], 'post');
end
