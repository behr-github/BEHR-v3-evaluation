function [ varargout ] = line_seg( X, Y, varargin )
%LINE_SEG(X,Y,...) Plot separate line segments
%   Given two matrices, X and Y, this will plot each row as a line segment,
%   by concatenating the rows into one long vector separated by NaNs. The
%   additional arguments that this function can take are exactly the same
%   as those by the normal line() function, as they are passed on
%   unaltered.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismatrix(X) || ~ismatrix(Y)
    E.badinput('X and Y must be 2-D at most')
end

if any(size(X) ~= size(Y))
    E.badinput('X and Y must be the same size')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

x = [];
y = [];

for a=1:size(X,1)
    x = cat(2,x,X(a,:));
    y = cat(2,y,Y(a,:));
    if a < size(X,1)
        x = cat(2,x,nan);
        y = cat(2,y,nan);
    end
end

l = line(x,y,varargin{:});

% Return the line handle only if requested.
if nargout > 0
    varargout{1} = l;
end

end

