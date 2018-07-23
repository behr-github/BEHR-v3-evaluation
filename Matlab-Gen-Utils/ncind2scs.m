function [start, count, stride] = ncind2scs(varargin)
%NCIND2SCD Convert Matlab indices into start, count, stride arguments
%   [ START, COUNT, STRIDE ] = NCIND2SCS( DIM1_IND, DIM2_IND, ... ) Given
%   logical or numeric indices along each dimension, will generate the
%   START, COUNT, and STRIDE argument that NCREAD understands. Note that
%   logical and numeric index vectors may be mixed.
%
%   [ ___ ] = NCIND2SCS( ___, 'min_ndim', NDIM ) will make the START,
%   COUNT, and STRIDE vectors have NDIM values. For dimensions that did not
%   have indices specified, they will be set to read in all values along
%   that dimension.

E = JLLErrors;
p = advInputParser;

p.addParameter('min_ndim', 0);
p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

indices = p.Unmatched;
n_dims = pout.min_ndim;
if n_dims == 0
    n_dims = numel(indices);
elseif n_dims < numel(indices)
    E.badinput('Requested fewer dimensions that gave indicies for');
end

start = ones(1, n_dims);
count = inf(1, n_dims);
stride = ones(1, n_dims);

for i = 1:numel(indices)
    this_index = indices{i};
    if ischar(this_index) || isstring(this_index)
        if strcmpi(this_index, 'all')
            % leave the values (1, inf, 1) - which will cause ncread() to
            % include all values along that dimension
            continue
        else
            E.badinput('The only string permitted for an index is "all"')
        end
    elseif islogical(this_index) && isvector(this_index)
        % We can just reuse the logic once we convert this logical index to
        % a numeric one
        this_index = find(this_index);
    end
    
    if isnumeric(this_index) && isvector(this_index)
        start(i) = this_index(1);
        step_size = diff(this_index);
        % diff() of a scalar produces an empty array, in which case we know
        % that count will be 1, so stride doesn't matter
        if ~isempty(step_size)
            if all(step_size == step_size(1))
                stride(i) = step_size(1);
            else
                warning('ncind2scs:stride', 'The index vector for dim %d has inconsistent spacing, stride will be set to 1', i);
            end
        end
        count(i) = (this_index(end) - this_index(1) + stride(i))/stride(i);
        if mod(count(i), 1) ~= 0
            warning('ncind2scs:stride', 'The stride for dim %d does not divide the span of indices evenly, an element may be omitted from the end', i);
            count(i) = floor(count(i));
        end
    else
        E.badinput('All INDEX inputs must be logical or numeric vectors')
    end
end

end

