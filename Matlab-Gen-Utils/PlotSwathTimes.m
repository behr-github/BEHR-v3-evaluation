%PlotSwathTimes - plots the edges of swaths chosen by the same logic as
%"read_omno2_v_aug2012 for OMI and MODIS Cloud.  The goal is to see how the
%swaths overlap.

lonmin = -125;    lonmax = -65;
latmin = 25;    latmax = 50;

date_start='2007/02/01';
date_end='2007/02/21';

satellite='OMI';
retrieval='SP';
modis_myd06_dir = '/Volumes/share/GROUP/SAT/MODIS/MYD06_L2_Collection5';
he5_dir = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';

Latitude=zeros(60,2000);
Longitude=zeros(60,300);

file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);

last_date=datestr(datenum(date_start)-1,26);

nfig=1;

total_days=datenum(date_end)-datenum(last_date)+1;
for j=1:total_days;
    %Read the desired year, month, and day
    R=addtodate(datenum(last_date), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    Data=struct('Longitude',0,'Latitude',0);
    
    short_filename=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5'];
    file_dir = fullfile(he5_dir,year,month); %Used both here to find all he5 files and in the swath for loop to identify each file.
    file=fullfile(file_dir,short_filename);
    sp_files = dir(file);
    n = length(sp_files);
    E=0;
    if isempty(sp_files);
        disp(['No Data Available For ',month,' ',day,' ',year])
    else
        for e=1:n
            %Read in each file, saving the hierarchy as 'hinfo'
            filename= sp_files(e).name;
            hinfo = h5info(fullfile(file_dir,filename));
            
            %Read in the full latitude data set; this will be used to determine
            %which pixels to read in later.
            Latitude = h5read(fullfile(file_dir,filename), h5dsetname(hinfo,1,2,1,2,'Latitude')); %h5dsetname takes 1) the object returned by h5info, 2) The indicies of the group tree 3) The last argument may be the index or name of the dataset of interest
            Longitude = h5read(fullfile(file_dir,filename), h5dsetname(hinfo,1,2,1,2,'Longitude'));
            Row = 1:size(Latitude,2); Row=repmat(Row,60,1);
            
            %Deletes any point that falls outside of the boundaries specified.
            Lat=Latitude; Lon=Longitude; %LatTest = Latitude; LonTest = Longitude;
            x=(Lon>lonmax | Lon<lonmin); xcol=find(sum(x,1));
            y=(Lat>latmax | Lat<latmin); ycol=find(sum(y,1));
            
            Lon(x) = NaN; Lon(y) = NaN; Lon(isnan(Lon))=[];
            Lat(x) = NaN; Lat(y) = NaN; Lat(isnan(Lat))=[];
            Row(x) = NaN; Row(y) = NaN; Row(isnan(Row))=[];
            
            %Lon(:,~xcol)=NaN;     Lon(:,~ycol)=NaN;     Lon(isnan(Lon))=[];
            %Lat(:,~xcol)=NaN;     Lat(:,~ycol)=NaN;     Lat(isnan(Lat))=[];
            
            if isempty(Lon)==1 || isempty(Lat)==1 || length(Lat)==1;
                fprintf('Swath %u is empty\n',e);
                continue
            else
                fprintf('Swath %u has points\n',e)
                %omi1 = find(Row==min(Row));%,1,'first'); 
                %omi2 = find(Row==min(Row));%,1,'last');
                %omi3 = find(Row==max(Row));%,1,'last');
                %omi4 = find(Row==max(Row));%,1,'first');
                dr = max(Row)-min(Row)+1;
                rows = min(Row):max(Row);
                omi = zeros(1,2*dr);
                for r=1:dr
                    omi(r) = find(Row==rows(r),1,'first');
                    omi(end-(r-1)) = find(Row==rows(r),1,'last');
                end
                
                OMILats = Lat(omi);
                OMILons = Lon(omi);
                %OMILats = [Lat(1,1), Lat(1,end), Lat(end,end), Lat(end,1), Lat(1,1)];
                %OMILons = [Lon(1,1), Lon(1,end), Lon(end,end), Lon(end,1), Lon(1,1)];
                OMITime = {['OMI t=',sp_files(e).name(29:30),':',sp_files(e).name(31:32)]};
                
                %disp(OMILats); disp(OMILons);
  
                %These were debugging plots, not needed now
%                 figure(20);
%                 m_proj('Albers Equal-Area Conic','lon',[min(Lon(:)), max(Lon(:))],'lat',[min(Lat(:)), max(Lat(:))]); m_coast('color','k'); %m_grid;
%                 m_line(Lon,Lat,'linestyle','none','marker','.');
%                 m_line(OMILons,OMILats,'linewidth',5);
                
                %MODIS Clouds
                d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
                x=numel(num2str(d2));
                if x==1;
                    julian_day=(['00',num2str(d2)]);
                elseif x==2;
                    julian_day=(['0',num2str(d2)]);
                elseif x==3;
                    julian_day=num2str(d2);
                end
                
                modis_file=(['MYD06_L2.A',year,julian_day,'*.hdf']);
                modis_files=dir(fullfile(modis_myd06_dir,year,modis_file));
                
                modis_lats = zeros(length(modis_files),5);
                modis_lons = zeros(length(modis_files),5);
                modis_times = cell(length(modis_files),1);
                
                for ii=1:length(modis_files);
                    mod_filename=modis_files(ii).name;
                    if str2double(mod_filename(19:22))<str2double(sp_files(e).name(29:32));
                        continue
                    elseif str2double(mod_filename(19:22))>str2double(sp_files(e+1).name(29:32));
                        continue
                    else
                        %For each file that fits the criteria mentioned
                        %above, import its latitude, longitude, and cloud
                        %fraction.
                        modis_times{ii} = ['MOD t=',mod_filename(19:20),':',mod_filename(21:22)];
                        mod_filename=fullfile(modis_myd06_dir,year,modis_files(ii).name); %Redefine the filename to have the full path to the file
                        mod_fileinfo=hdfinfo(mod_filename);
                        Latitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Latitude')); Latitude=double(Latitude); %Latitude=Latitude(:);
                        Longitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Longitude')); Longitude=double(Longitude); %Longitude=Longitude(:);
                        CloudFraction=hdfread(hdf_dsetID(mod_fileinfo,1,2,'Cloud_Fraction')); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:);
                        CloudFraction(CloudFraction==127)=100; CloudFraction=CloudFraction*0.009999999776482582;
                        
                        x=find(Longitude>lonmax | Longitude<lonmin);
                        y=find(Latitude>latmax | Latitude<latmin);
                                               
                        %westIndMODIS = [1,1];%find(Longitude == min(Longitude));
                        %eastIndMODIS = [end,end];%find(Longitude == max(Longitude));
                        %southIndMODIS = [end,1];%find(Latitude == min(Latitude));
                        %northIndMODIS = [1,end];%find(Latitude == max(Latitude));
                        
                        modis_lats(ii,:) = [Latitude(1,1), Latitude(1,end), Latitude(end,end), Latitude(end,1), Latitude(1,1)];
                        modis_lons(ii,:) = [Longitude(1,1), Longitude(1,end), Longitude(end,end), Longitude(end,1), Longitude(1,1)];
                        
                        
                        
                    end
                end
                
                %Do the map... do the monster map
                %westIndOMI = find(Data(E).Longitude == min(Data(E).Longitude));
                %eastIndOMI = find(Data(E).Longitude == max(Data(E).Longitude));
                %southIndOMI = find(Data(E).Latitude == min(Data(E).Latitude));
                %northIndOMI = find(Data(E).Latitude == max(Data(E).Latitude));
                
                %OMILats = [Data(E).Latitude(westIndOMI), Data(E).Latitude(northIndOMI), Data(E).Latitude(eastIndOMI), Data(E).Latitude(southIndOMI), Data(E).Latitude(westIndOMI)];
                %OMILons = [Data(E).Longitude(westIndOMI), Data(E).Longitude(northIndOMI), Data(E).Longitude(eastIndOMI), Data(E).Longitude(southIndOMI), Data(E).Longitude(westIndOMI)];
                %OMITime = ['OMI t=',sp_files(e).name(29:30),':',sp_files(e).name(31:32)];
                figure(nfig)
                title(sprintf('%s: %s  ',date,OMITime{1}),'fontsize',16);
                m_proj('Albers Equal-Area Conic','lon',[-150 -50],'lat',[0,75]);
                m_coast('color','k'); m_states('k');
                m_grid('linestyle','none');
                
                
                disp(OMITime);
                m_line(OMILons,OMILats,'linewidth',5);
                %m_text(mean(OMILons(1:4)),mean(OMILats(1:4)),OMITime(1))
                
                for m = 1:size(modis_lats,1)
                    fprintf('Adding MODIS swath %u\n',m);
                    m_line(modis_lons(m,:),modis_lats(m,:),'color','r')
                    m_text(mean(modis_lons(m,:)),mean(modis_lats(m,:)),modis_times(m));
                    
                end
                nfig=nfig+1;
                %pause
            end
        end
    end
    
end