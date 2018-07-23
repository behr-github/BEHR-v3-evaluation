function [  ] = label_subfigs( fig, varargin )
%LABEL_SUBFIGS Label subfigures (child axes) with letters
%   LABEL_SUBFIGS( FIG, ... ) Labels the axes within figure FIG with
%   letters as e.g. "(a)", "(b)", etc. going horizontally first then
%   vertically. Additional arguments accepted by label_axis_with_letter can
%   be passed.

p = advInputParser;
p.addParameter('ax', []);
p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

ax = pout.ax;

% Remove arguments for this function
varargin = update_params('remove', varargin, pout);

if isempty(ax)
    xx = isgraphics(fig.Children, 'axes');
    figax = fig.Children(xx);
else
    figax = ax;
end
pos = cat(1, figax.Position);
[~,ord_ind] = sortrows(pos,[-2 1]); % sort by vertical position in reverse order, then by horizontal position in forward order

for a=1:numel(ord_ind)
    % add the subfigure letters. 
    label_axis_with_letter(a, 'ax', figax(ord_ind(a)), varargin{:});
end

end

