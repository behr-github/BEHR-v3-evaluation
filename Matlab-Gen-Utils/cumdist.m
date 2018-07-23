function [ varargout ] = cumdist( Data, nbins, varargin )
%CUMDIST( DATA ) Plot a cumulative distribution of the values of Data
%   By default, bins the data into 10 bins and plots the normalized
%   cumulative distribution. Because it calls HIST to bin the data, if Data
%   is a matrix, it will plot a cumulative distribution for each column of
%   Data.
%
%   CUMDIST( DATA, NBINS ) will bin the data into NBINS bins instead of the
%   default 10.
%
%   CUMDIST( ___ , 'nonorm' ) will plot the cumulative sum of counts rather
%   than a normalized fraction.
%
%   CUMDIST( ___ , 'reverse' ) will plot the cumulative sum from greatest
%   to least bin rather than least to greatest.
%
%   CUMDIST( ___ , linespec options ) allows you to manipulate the line
%   specification. Any valid input to PLOT will work. Note that if you
%   override color and input a matrix for Data, the series will not be
%   distinguished.
%
%   C = CUMDIST( ___ ) like HIST, when a single output is requested,
%   this will return the cumulative values (normalized or not, as specified
%   by the presence or absence of the 'nonorm' parameter. They will not be
%   plotted.
%
%   [C, BIN_CENTERS] = CUMDIST( ___ ) also returns the bin centers.

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isnumeric(Data)
    E.badinput('Data must be numeric');
end
if ~exist('nbins','var')
    nbins = 10;
elseif ischar(nbins)
    % Handles if a parameter is input without a number of bins
    varargin{end+1} = nbins;
    nbins = 10;
elseif ~isnumeric(nbins) || ~isscalar(nbins) || nbins < 1
    E.badinput('nbins must be a positive scalar number')
end

norm_bool = true;
rev_bool = false;
keep_log = true(size(varargin));
for a = 1:numel(varargin)
    if strcmpi('nonorm',varargin{a})
        norm_bool = false;
        keep_log(a) = false;
    elseif strcmpi('reverse',varargin{a})
        rev_bool = true;
        keep_log(a) = false;
    end
end
varargin = varargin(keep_log);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% bin the data
[counts, bin_centers] = hist(Data, nbins);
if isrow(counts); 
    counts = counts'; 
end
if rev_bool
    counts = flipud(counts);
    c = cumsum(counts,1);
    c = flipud(c);
else
    c = cumsum(counts,1);
end
if norm_bool
    cmax = max(c,[],1);
    cmax = repmat(cmax,size(c,1),1);
    c = c ./ cmax;
end

if nargout == 0
    bc = repmat(bin_centers,1,size(c,2));
    plot(bc,c,varargin{:});
else
    varargout{1} = c;
    varargout{2} = bin_centers;
end

