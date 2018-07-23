function ind = index_by_n(loop_ind, n_in_loop, n_per_loop)
%INDEX_BY_N Calculate an index that needs to jump by some multiple between loops
%   IND = INDEX_BY_N( LOOP_IND, N_IN_LOOP, N_PER_LOOP ) Assuming you have a
%   for loop with index i, and for each iteration you want to save three
%   values in a vector, the first element would need to go at indices 1, 4,
%   7, etc., the second at 2, 5, 8, and so on. This is a convenience
%   function to help calculate those values. LOOP_IND is the loop index,
%   N_IN_LOOP is the "inner" index, so 1, 2, or 3 in this example.
%   N_PER_LOOP would be 3 in this example, i.e. how many elements will be
%   created within the loop.

ind = (loop_ind - 1) * n_per_loop + n_in_loop;

end

