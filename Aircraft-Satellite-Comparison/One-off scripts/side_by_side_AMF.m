function [  ] = side_by_side_AMF(table1, table2)
%side_by_side_AMFs Plot profiles side by side with WRF vs. in-situ AMFs below

profiles_dir_1 = '/Users/Josh/Documents/MATLAB/Figures/Sat Verification/DISCOVER California/Scaled profile shapes generated 15 Aug/';
profiles_dir_2 = '/Users/Josh/Documents/MATLAB/Figures/Sat Verification/DISCOVER California/Unscaled profile shapes generated 15 Aug/';

% Iterate over all AMFs and profiles, plotting profiles side by side with
% the AMFs listed below each plot

files1 = dir(fullfile(profiles_dir_1,'*.fig'));
files2 = dir(fullfile(profiles_dir_2,'*.fig'));

if numel(files1) ~= numel(files2)
    warning('Number of figures is not equal');
end

dates = datestr(table2array(table1.table(:,1)));
profnums = table2array(table1.table(:,2));

wrfamfs1 = table2array(table1.table(:,3));
insituamfs1 = table2array(table1.table(:,4));
wrfamfs2 = table2array(table2.table(:,3));
insituamfs2 = table2array(table2.table(:,4));



for a=1:numel(files1)
    
    filename1 = fullfile(profiles_dir_1, files1(a).name);
    filename2 = fullfile(profiles_dir_2, files2(a).name);
    p1 = openfig(filename1); p2 = openfig(filename2);
    f1 = montagefigures([p1,p2],2,1);
    close([p1,p2]);
    
    %Get the non-legend axes so we can set the title and xlabel
    ax = findall(f1,'type','axes');
    for b=1:numel(ax)
        if ~isempty(get(ax(b),'tag'));
            ax(b) = nan;
        end
    end
    ax(isnan(ax)) = [];
    
    sp1_title = sprintf('%s #%d: %s',dates(a,:),profnums(a),table1.top_extrap);
    sp2_title = sprintf('%s #%d: %s',dates(a,:),profnums(a),table2.top_extrap);
    xname1 = sprintf('WRF AMF: %.4f; In-Situ AMFS: %.4f',wrfamfs1(a),insituamfs1(a));
    xname2 = sprintf('WRF AMF: %.4f; In-Situ AMFS: %.4f',wrfamfs2(a),insituamfs2(a));
    
    
    %figure out which axis is right and which is left
    pos1 = get(ax(1),'position'); pos2 = get(ax(2),'position');
    if pos1(1) < pos2(1);
        set(get(ax(1),'title'),'string',sp1_title,'fontsize',20); set(get(ax(1),'xlabel'),'string',xname1,'fontsize',18);
        set(get(ax(2),'title'),'string',sp2_title,'fontsize',20); set(get(ax(2),'xlabel'),'string',xname2,'fontsize',18);
    else
        set(get(ax(2),'title'),'string',sp1_title,'fontsize',20); set(get(ax(2),'xlabel'),'string',xname1,'fontsize',18);
        set(get(ax(1),'title'),'string',sp2_title,'fontsize',20); set(get(ax(1),'xlabel'),'string',xname2,'fontsize',18);
    end
    
    
    fprintf('Paused\n');
    pause;
    close(f1);
    
end

end

