matdir = '/Volumes/share/GROUP/INTEX-B/Matlab files';
matpat = 'INTEXB*.mat';

f = dir(fullfile(matdir,matpat));
for n=1:numel(f)
    load(fullfile(matdir,f(n).name));
    temp = remove_merge_fills(Merge,'TEMP_STAT_C');
    pres = remove_merge_fills(Merge,'STAT_PRESSURE');
    
    % The wikipedia article on potential temperature gives R/c_p as 0.286
    theta = temp .* (1013./pres).^0.286;
    theta(isnan(theta))=-9999;
    Merge.Data.THETA.Unit = 'K';
    Merge.Data.THETA.Fill = -9999;
    Merge.Data.THETA.Values = theta;
    
    save(fullfile(matdir,'THETA',f(n).name),'Merge');
end