function [  ] = plot_hourly_profiles( BinStruct )
%PLOT_HOURLY_PROFILES Plots the output from bin_profile_by_start_time
%   Generates plots from the output structure of bin_profile_by_start_time.
%   If you pass it the whole structure, it will make subplots for each
%   field; if you pass it one of the first fields (e.g. BinStruct.AllSites
%   or BinStruct.Site01) it will just plot that one.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 May 2015

if isfield(BinStruct,'AllSites') % indicates the whole structure was passed
    fns = fieldnames(BinStruct);
    n = numel(fns);
    x = floor(sqrt(n));
    y = ceil(n/x);
    figure; subplot(x,y,1);
    for a=1:n
        subplot(x,y,a);
        make_individual_plot(BinStruct.(fns{a}))
        title(fns{a})
    end
else
    figure;
    make_individual_plot(BinStruct)
end


end

function make_individual_plot(SiteStruct)
n = size(SiteStruct.final_no2_bins,2);
l = gobjects(n,1);
lstr = cell(1,n);
co = get(groot,'defaultAxesColorOrder');
max_co = size(co,1);
line_styles = {'-','--',':'};
for a=1:n
    if a > max_co; 
        c = mod(a,max_co)+1;
    else
        c = a;
    end
    s = ceil(a/max_co);
    l(a) = line(SiteStruct.final_no2_bins(:,a), SiteStruct.final_pres_bins(:,a), 'color', co(c,:),'linewidth',2,'linestyle',line_styles{s});
    lstr{a} = sprintf('%02d00',SiteStruct.bin_start_hours(a));
end

set(gca,'ydir','reverse');
set(gca,'fontsize',16);
xlabel('[NO_2] (pptv)');
ylabel('Pressure (hPa)');
legend(l,lstr);

end