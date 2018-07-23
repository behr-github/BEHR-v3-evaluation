% Albedo_avg - calculate the average albedo over a point using mcd43c3 files
%
% Josh Laughner <joshlaugh5@gmail.com> 30 Apr 2014

start_date = '1-Jun-2013';
end_date = '31-Jul-2013';

center_lat = 32.90328;
center_lon = -87.24994;

mcd43_dir = '/Volumes/share/GROUP/SAT/MODIS/MCD43C3';

if ~exist(mcd43_dir,'dir'); error('albavg:mcd43dir','Given mcd43 directory does not exist.'); end

dates = datenum(start_date):datenum(end_date);
alb_values_latlon = -99 .* ones(3,numel(dates));
for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    day = num2str(modis_date_to_day(curr_date));
    
    filename = fullfile(mcd43_dir,year,['MCD43C3.A',year,day,'*.hdf']);
    files = dir(filename);
    if isempty(files)
        fprintf('%s not found. \n',filename);
    elseif length(files)==1
        fprintf('%s found. Proceeding... \n',files(1).name);
        [alb,lons,lats] = read_modis_mcd43c3_2014(fullfile(mcd43_dir,year,files(1).name));
        
        lon_diff = abs(lons - center_lon);
        lat_diff = abs(lats - center_lat);
        
        col = find(lon_diff(1,:) == min(lon_diff(1,:)));
        row = find(lat_diff(:,1) == min(lat_diff(:,1)));
        
        alb_values_latlon(1,d) = alb(row,col);
        alb_values_latlon(2,d) = lats(row,col);
        alb_values_latlon(3,d) = lons(row,col);
    else
        error('albavg:mcd43files','Too many files found for comfort. Aborting');
    end
end

zz = find(alb_values_latlon(1,:) == -99);
alb_values_latlon(:,zz) = [];
avg_alb = nanmean(alb_values_latlon(1,:))