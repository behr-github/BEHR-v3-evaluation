function [ bins, bin_centers, bin_edges ] = bin_data( X, Y, binning )
%BIN_DATA Bins values in Y by the corresponding value in X. Returns cell array.
%   Returns a cell array of values from Y binned based on the corresponding
%   values in X. Requires input of X and Y. The third argument is optional
%   and can take on several forms:
%
%       - If it is a scalar, it is assumed to be the number of bins to
%       generate. If unspecified, it is assumed that this number is 10.
%
%       - If it is a vector, V, with length n, then n-1 bins will be
%       generated assuming each value is the edge between bins, i.e. the
%       first bin will be values of X s.t. V(1) <= X < V(2), the second bin
%       will be V(2) <= X < V(3) and so on.
%
%       - If it is an n-by-2 matrix, then n bins will be generated using
%       the first column as the bottom of the bins and the second column as
%       the top.
%
%   This will return a cell array of the bins, plus a vector of bin centers
%   and an n-by-2 array of bin edges.
%
%   Note that inputting X and Y as a matrix vs. a vector makes no
%   difference, they will be binned on an element-by-element basis.
%
%   Josh Laughner <joshlaugh5@gmail.com> 12 June 2015

E = JLLErrors;
DEBUG_LEVEL = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(X) ~= numel(Y)
    E.badinput('X and Y must have the same number of elements');
end

if nargin < 3
    binning = 10;
elseif ~isnumeric(binning)
    E.badinput('binning must be numeric, be is scalar, vector, or matrix.')
elseif isscalar(binning) && binning < 1
    E.badinput('Number of bins must be positive')
elseif ismatrix(binning) && ~isvector(binning)
    if size(binning,2) ~= 2
        E.badinput('If the optional "binning" input is a matrix, it is expected that the second dimension be 2, that is, each row specifies the edges of a bin.')
    elseif any(binning(:,1) >= binning(:,2))
        E.badinput('If the optional "binning" input is a matrix, it is expected that the left column be the lower bound for each bin.')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PREPARATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Make into vectors to simplify some operations later
X = X(:);
Y = Y(:);

% Figure out the binning. 3 options:
%   1) Make evenly spaced bins based on the limits of X
if isscalar(binning)
    nbins = binning;
    % Makes the bin edges - there will be 1 more edge than bins. Making it
    % a column vector makes it easier to make the actual edges matrix.
    bin_edge_vec = (linspace(min(X), max(X), nbins + 1))';
    bin_edges = bin_edge_vec2mat(bin_edge_vec);
elseif isvector(binning)
    bin_edges = bin_edge_vec2mat(binning);
elseif ismatrix(binning)
    bin_edges = binning;
else
    E.unknownError('binning is apparently not a scalar, vector, or matrix.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BINNING AND OUTPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bins = cell(size(bin_edges,1),1);

% Bin Y based on X
for a=1:numel(X)
    xx = X(a) >= bin_edges(:,1) & X(a) < bin_edges(:,2);
    if sum(xx) > 0 %Some pressure values fall between bins
        bins{xx} = cat(2,bins{xx},Y(a));
    end
end

% Make the bin centers and call it a day
bin_centers = mean(bin_edges,2);


end

function edges = bin_edge_vec2mat(vec)
    edges(:,1) = vec(1:end-1);
    edges(:,2) = vec(2:end);
end

