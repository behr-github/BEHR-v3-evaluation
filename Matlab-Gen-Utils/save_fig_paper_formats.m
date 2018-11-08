function save_fig_paper_formats(fig, save_name)
% SAVE_PAPER_FORMATS Save .fig, .png. and .eps versions of figure for papers
%   SAVE_PAPER_FORMATS( FIG, SAVE_NAME ) Saves the figure given by the
%   handle FIG to as SAVE_NAME in .fig, .png, and .eps formats.

set(fig,'paperpositionmode','auto');
savefig(fig, [save_name, '.fig'], 'compact');
print(fig, save_name, '-dpng','-r0');
print_eps(fig, save_name);
end