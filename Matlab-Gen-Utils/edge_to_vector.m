function vec = edge_to_vector(M, varargin)
%EDGE_TO_VECTOR Extract the edges of a matrix into a vector
%
%   VEC = EDGE_TO_VECTOR( M ) Returns the values along the edge of M in a
%   vector, VEC, starting in the upper left (1,1) corner of M and going
%   clockwise.
%
%   VEC = EDGE_TO_VECTOR( M, 'close' ) By default, the (1,1) value is not
%   repeated at the end of the vector. Including the 'close' flag causes it
%   to be repeated.

p = advInputParser;
p.addFlag('close');
p.parse(varargin{:});

pout = p.Results;
do_close = pout.close;

E = JLLErrors;
if ~ismatrix(M)
    E.badinput('M must be a matrix (cannot have >2 dimensions)')
end

vec = veccat(M(1,1:end-1), M(1:end-1,end), fliplr(M(end,2:end)), flipud(M(2:end,1)), 'column');
if do_close
    vec = veccat(vec, M(1,1));
end

end

