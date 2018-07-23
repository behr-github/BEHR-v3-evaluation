function [  ] = histsemilog( Y, varargin )
%HISTSEMILOG Plots a histogram with a log-scale Y axis
%   HISTSEMILOG( Y ) Bins data in Y into 10 bins and plots with a log-scale
%   y-axis.
%
%   HISTSEMILOG( Y, M ) or HISTSEMILOG( Y, X ) behaves as would HIST( Y, M )
%   or HIST( Y, X )

E = JLLErrors;

if isa(Y,'matlab.graphics.axis.Axes')
    ax = Y;
    if numel(varargin) >= 1
        Y = varargin{1};
    else
        E.badinput('If giving an axis as the first input, there must be at least one more input specifying the data to plot');
    end
else
    ax = gca;
end

[n, xout] = hist(Y, varargin{:});
bar(ax, xout, n, 'barwidth', 1, 'basevalue', 1);
set(ax, 'yscale', 'log');

end

