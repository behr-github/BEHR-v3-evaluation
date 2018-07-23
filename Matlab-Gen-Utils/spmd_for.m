function [ svec ] = spmd_for( for_vec, li, nwork )
%SPMD_FOR Return a vector for a for loop in an SPMD block
%   SVEC = SPMD_FOR( FOR_VEC, LI, NWORK ) Given a vector FOR_VEC that is
%   the vector to iterate over, the lab index LI, and the number of
%   labs/workers NWORK, this will break up FOR_VEC into NWORK pieces and
%   return the piece that lab LI should do as the vector SVEC. If NWORK is
%   0, FOR_VEC is returned unaltered.

if nwork == 0
    svec = for_vec;
    return;
elseif li > nwork
    error('spmd_for:bad_index','LI (labindex) cannot be greater than NWORK (the number of workers)');
elseif li < 0
    error('spmd_for:bad_index','LI (labindex) cannot be < 0');
end

n_per_worker = ceil(numel(for_vec) / nwork);
start_ind = (li-1) * n_per_worker + 1;
end_ind = min(li * n_per_worker, numel(for_vec));
svec = for_vec(start_ind:end_ind);

end

