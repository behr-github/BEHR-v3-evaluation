function [  ] = split_subplots( fig )
%SPLIT_SUBPLOTS Split subplots into individual figures
%   SPLIT_SUBPLOTS( FIG ) Place each child object in the figure
%   specified by the figure handle FIG into its own figure.

ch = get(fig,'children');
ch = ch(isgraphics(ch,'axes'));

default_pos = [0.1300, 0.1100, 0.7750, 0.8150];

for a=1:numel(ch)
    nf = figure;
    copyobj(ch(a), nf);
    nch = get(nf, 'children');
    nch.Position = default_pos;
end


end

