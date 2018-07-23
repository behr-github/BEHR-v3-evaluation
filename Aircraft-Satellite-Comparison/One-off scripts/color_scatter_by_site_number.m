% color by site number
% assumes that you have run Run_Spiral_Verification and have the variables
% airno2_iall, behrno2_iall, and db_iall in the workspace.


pns = cell2mat(db_iall.profnums);
sitenums = floor(mod(pns,100000)/1000);
usitenums = unique(sitenums(~isnan(sitenums)));

figure; 

colors = {'b','r',[0 0.5 0],'m','k','c','y',[0.5 0 0.5]};

for a=1:numel(usitenums)
    xx = sitenums == usitenums(a);
    b=a+2;
    l(b) = line(airno2_iall(xx), behrno2_iall(xx), 'color', colors{a}, 'linestyle', 'none', 'marker', 'o');
    leg{b} = sprintf('Site %d',a);
end

axes_equal();
[xline,yline,str] = calc_fit_line(airno2_iall,behrno2_iall,'regression','rma');
l(1) = line(xline,yline,'color','k','linestyle','--','linewidth',2);
leg{1} = str;

l(2) = line(xline,xline,'color','r','linestyle',':','linewidth',2);
leg{2} = '1:1';

legend(l', leg);
set(gca,'fontsize',16)
xlabel('Aircraft NO_2 columns (molec./cm^2)')
ylabel('BEHR NO_2 columns (molec./cm^2)')