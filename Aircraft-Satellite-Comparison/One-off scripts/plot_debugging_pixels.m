% plot_debugging_pixels
%
%   Help plot pixels in a debugging structure from "Run Spiral
%   Verification"

debg = db;

n=numel(debg);
for a=1:n
    s = size(debg(a).db.loncorn);
    if s(2) > 0;
        for b=1:s(2)
            m_line([debg(a).db.loncorn(:,b); debg(a).db.loncorn(1,b)], [debg(a).db.latcorn(:,b); debg(a).db.latcorn(1,b)],'color', 'k', 'linewidth',3);
        end
    end
end