function S = cat_gc_struct(varargin)
%CAT_GC_STRUCT Concatenates structures output from read_geos_output

E=JLLErrors;
if any(~iscellcontents(varargin,'isstruct'))
    E.badinput('All inputs must be structures')
end

S = varargin{1};
n = numel(varargin);
for a=2:n
    for b=1:numel(S)
        % Some variables have two space and one time dimension, some
        % have three space and one time dimension. The time dimensions is
        % always the last, however this will break if there's only one
        % entry in time.  It will try to check for that, but cannot be sure.
        % We'll check by seeing if there's 47 elements along the time 
        % dimension, since that's the usual number of vertical levels.
        time_dim = ndims(varargin{a}(b).dataBlock);
        if size(varargin{a}(b).dataBlock,time_dim) == 47 && a == 1
            warning('May have concatenated along the vertical dimension')
        end
        S(b).dataBlock = cat(time_dim, S(b).dataBlock, varargin{a}(b).dataBlock);
        S(b).tVec = cat(1, S(b).tVec, varargin{a}(b).tVec);
        S(b).tEdge = cat(1, S(b).tEdge, varargin{a}(b).tEdge(2:end));
    end
end

end
