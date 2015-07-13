%% Details of the function:
% NAME:
% readMERRA(files2process,plotting)
%----------------------------------
% PURPOSE:
%   read inst6_3d_ana_Np MERRA product
%   Short name:MAI6NPANA 
%   Long name:MERRA DAS 3d analyzed state on pressure, 
%   http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/inst6_3d_ana_Np.shtml
%   http://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset.pl?LOOKUPID_List=MAI3CPASM
%   files2process is a switch to process 1 single file (0) or the whole
%   directory (1); plotting is false or true to show plots
%   
%
% CALLING SEQUENCE:
%  out = readMERRA(files2process,plotting)
%
% INPUT:
%  - MERRA300.prod.assim.inst6_3d_ana_Np.yyyymmdd.SUB.hdf file
%  - files under: F:\ARISE\ArcticCRF\METdata\Sep2014\ARISEdomain\
% 
% 
% OUTPUT:
%  - out struct with all fields:
% Variable Name	Dims	Description	            Units 
% SLP 	        2D	    Sea-level   pressure 	Pa 
% PS	        3D	    Surface     pressure	Pa
% H	            3D	    Geopotential    height	m
% T 	        3D	    Air temperature         K 
% U 	        3D	    Eastward wind component m s-1
% V 	        3D	    Northward wind componentm s-1
% QV 	        3D	    Specific humidity 	    kg kg-1
% O3 	        3D	    Ozone mixing ratio 	    kg kg-1
%
%
% DEPENDENCIES:
%  getfullfilename__
%
% NEEDED FILES/INPUT:
% 
%
% EXAMPLE:
%  - out = readMERRA(1,'true')
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Jan-21-2015
% -------------------------------------------------------------------------
%% function routine
function out = readMERRA

%% choose input directory or file
% ------------------------------ %
% assign default values
% if varargin==0
  files2process = 1;
  plotting      = 'false';
% end

infile = getfullname__('MERRA300.prod.assim.inst6_3d_ana_Np.*.SUB.hdf','F:','Select a MERRA* file');
[pname, fname, ext] = fileparts(infile);
pname = [pname, filesep];

if files2process==0
    fileinfo = hdfinfo(infile);
    dm    = strsplit(fname,'.');
    dates = dm(4);
    allfiles.name = [fname,ext];
else
    allfiles = dir([pname,'MERRA300.prod.assim.inst6_3d_ana_Np.*.SUB.hdf']);
    for i=1:length(allfiles)
        dm{i} = strsplit_ms(allfiles(i).name,'.');
        dates(i) = dm{i}(4);
        fileinfo = hdfinfo([pname,allfiles(i).name]);
    end
end
% some conversion constants %
r    = 6378*1000;       % Earth radius
Mair = 28.97;           % Molecular weight of dry air
Mh2o = 18.01;           % Molecular weight of water
Mo3  = 48.00;           % Molecular weight of ozone
Mo2  = 32.00;           % Molecular weight of oxygen
Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)
%---------------------------%
% read file data
%----------------
for i=1:length(dates)-1
    
    datest = strcat('m',dates{:,i},'_'); 
    out.(datest).time      = double(hdfread([pname,allfiles(i).name],'/time'));
    out.(datest).longitude = double(hdfread([pname,allfiles(i).name],'/longitude')); p =length(out.(datest).longitude);
    out.(datest).latitude  = double(hdfread([pname,allfiles(i).name],'/latitude'));  q =length(out.(datest).latitude);
    out.(datest).sealevelp = double(hdfread([pname,allfiles(i).name],'/slp'));            %[Pa-2D]
    out.(datest).surfpress = double(hdfread([pname,allfiles(i).name],'/ps'));             %[Pa-3D]
    out.(datest).plevels   = double(hdfread([pname,allfiles(i).name],'/levels'));    s =length(out.(datest).plevels);%[hPa];
    out.(datest).geopothgt = double(hdfread([pname,allfiles(i).name],'/h'));              %[m-3D]
    out.(datest).airT      = double(hdfread([pname,allfiles(i).name],'/t'));              %[K-3D]
    out.(datest).eastwind  = double(hdfread([pname,allfiles(i).name],'/u'));              %[m/s-3D]
    out.(datest).northwind = double(hdfread([pname,allfiles(i).name],'/v'));              %[m/s-3D]
    out.(datest).watermmr  = double(hdfread([pname,allfiles(i).name],'/qv'));             %[kg/kg-3D] specific humidity
    out.(datest).ozonemmr  = double(hdfread([pname,allfiles(i).name],'/o3'));             %[kg/kg-3D] ozone mixing ratio
    
                                                           
% choose one time
%      out.(strcat('m',dates{:,i},'_')).sealevelp = squeeze(out.(strcat('m',dates{:,i},'_')).sealevelp(end,:,:));  % take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).surfpress = squeeze(out.(strcat('m',dates{:,i},'_')).surfpress(end,:,:));  % take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).geopothgt = squeeze(out.(strcat('m',dates{:,i},'_')).geopothgt(end,:,:,:));% take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).airT      = squeeze(out.(strcat('m',dates{:,i},'_')).airT(end,:,:,:));     % take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).eastwind  = squeeze(out.(strcat('m',dates{:,i},'_')).eastwind(end,:,:,:)); % take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).northwind = squeeze(out.(strcat('m',dates{:,i},'_')).northwind(end,:,:,:));% take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).watermmr  = squeeze(out.(strcat('m',dates{:,i},'_')).watermmr(end,:,:,:)); % take only 18:00UTC
%      out.(strcat('m',dates{:,i},'_')).ozonemmr  = squeeze(out.(strcat('m',dates{:,i},'_')).ozonemmr(end,:,:,:)); % take only 18:00UTC
     
    % replace missiing values with NaN
    out.(datest).sealevelp ( out.(datest).sealevelp >9E14) = NaN;
    out.(datest).surfpress ( out.(datest).surfpress >9E14) = NaN;
    out.(datest).geopothgt ( out.(datest).geopothgt >9E14) = NaN;
    out.(datest).airT      ( out.(datest).airT      >9E14) = NaN;
    out.(datest).eastwind  ( out.(datest).eastwind  >9E14) = NaN;
    out.(datest).northwind ( out.(datest).northwind >9E14) = NaN;
    out.(datest).watermmr  ( out.(datest).watermmr  >9E14) = NaN;
    out.(datest).ozonemmr  ( out.(datest).ozonemmr  >9E14) = NaN;
    
    % create additional parameters
    out.(datest).watervsat = 611*exp(17.67*(out.(datest).airT - 273.16)./...
                                                               (out.(datest).airT - 29.65));     % Clausius-Clapeyron;sat vapor pressure [Pa]
    out.(datest).potenT    = (out.(datest).airT ).*(repmat(1000 ,[s,q,p])./...
                                                  repmat(out.(datest).plevels' ,[1,q,p])).^0.286;% this is potential temperature
    out.(datest).st700T    =  squeeze(out.(datest).potenT(13,:,:))  - ...
                                                  squeeze(out.(datest).potenT(1,:,: ));          % lower tropospheric static stability (LTS)
    out.(datest).st850T    =  squeeze(out.(datest).potenT(7,:,:))  - ...
                                                  squeeze(out.(datest).potenT(1,:,: ));          % lower tropospheric static stability (at 850 mb)
    out.(datest).st925T    =  squeeze(out.(datest).potenT(4,:,:))  - ...
                                                  squeeze(out.(datest).potenT(1,:,:));           % low level static stability(LLSS)
    out.(datest).geomhgt   = (out.(datest).geopothgt*r)./...
                                                 (r - out.(datest).geopothgt);                   % conversion from geopotential height to geometrical height
    out.(datest).H2Oppm    =  out.(datest).watermmr * (Mair/Mh2o) *1e6;      % conversion from MMR to ppmv
    out.(datest).H2Opart   =  out.(datest).H2Oppm.*...
                                                  repmat(out.(datest).plevels' ,[1,q,p])*1e-6;   % conversion from ppmv to hPa
    out.(datest).H2Odens   =  out.(datest).H2Opart *100 ./...
                                                  (R*out.(datest).airT);                         % row=P/RT: conversion from hPa to kg/m3
    out.(datest).H2Omolec  =  out.(datest).H2Odens * (1/Mh2o) *...
                                                  Na *(1000/(100)^3);                                                % conversion from kg/m3 to #/cm3
    
    out.(datest).RH        =  ((out.(datest).H2Opart*100)./out.(datest).watervsat)*100;% RH
    
    out.(datest).O3ppm     =  out.(datest).ozonemmr * (Mair/Mo3) *1e6;       % conversion from MMR to ppmv
    out.(datest).O3part    =  out.(datest).O3ppm.*...
                                                  repmat(out.(datest).plevels' ,[1,q,p])*1e-6;   % conversion from ppmv to hPa
    out.(datest).O3dens    =  out.(datest).O3part * 100./...
                                                  (R*out.(datest).airT);                         % row=P/RT: conversion from hPa to kg/m3
    out.(datest).O3molec   =  out.(datest).O3dens * (1/Mo3) *...
                                                  Na *(1000/(100)^3);                                                % conversion from kg/m3 to #/cm3
    out.(datest).airdens   =  repmat(out.(datest).plevels' ,[1,q,p]) * 100./...
                                                  (R*out.(datest).airT);                         % row=P/RT: conversion from hPa to kg/m3
    out.(datest).airmolec  =  out.(datest).airdens * (1/Mair) *...
                                                  Na *(1000/(100)^3);                                                % conversion from kg/m3 to #/cm3
    out.(datest).O2dens    =  repmat(out.(datest).plevels' ,[1,q,p])*0.21*100./...
                                                  (R*out.(datest).airT);                         % row=P/RT: conversion from hPa to kg/m3
    out.(datest).O2molec   =  out.(datest).O2dens * (1/Mo2) *...
                                                  Na *(1000/(100)^3);                                                % conversion from kg/m3 to #/cm3
                                              
    % close current file
    hdfml('closeall');
    
    % read AMSR sea-ice line
    idate_ = dates{:,i}; idate = idate_(1:end);
    tmp = readAMSR2dataV2(idate);
    lab = strcat('amsr',idate);
    amsr.(datest).lon = tmp.(lab).longrid;
    amsr.(datest).lat = tmp.(lab).latgrid;
    amsr.(datest).gradline = tmp.(lab).Ith;
    clear tmp;
end
% repmat(out.(strcat('m',dates{:,i},'_')).plevels' ,[1,q,p])creates 3D
% level matrix for the whole domain
%%

%% create daily average over the selected domain
ColorSet = varycolor(length(dates));
legendall ={};
legendall1={};
for i=1:length(dates)-1
    datest = strcat('m',dates{:,i},'_'); 
    out.(datest).avgLon       = mean(out.(datest).longitude);
    out.(datest).stdLon       = std (out.(datest).longitude);
    out.(datest).avgLat       = mean(out.(datest).latitude);
    out.(datest).stdLat       = std (out.(datest).latitude);
    out.(datest).avgsealevelp = nanmean(nanmean(out.(datest).sealevelp));
    out.(datest).stdsealevelp = nanstd (nanstd (out.(datest).sealevelp));
    out.(datest).avgsurfpress = nanmean(nanmean(out.(datest).surfpress));
    out.(datest).stdsurfpress = nanstd (nanstd (out.(datest).surfpress));
    out.(datest).avgst700T    = nanmean(nanmean(out.(datest).st700T));
    out.(datest).stdst700T    = nanstd (nanstd (out.(datest).st700T));
    out.(datest).avgst850T    = nanmean(nanmean(out.(datest).st850T));
    out.(datest).stdst850T    = nanstd (nanstd (out.(datest).st850T));
    out.(datest).avgst925T    = nanmean(nanmean(out.(datest).st925T));
    out.(datest).stdst925T    = nanstd (nanstd (out.(datest).st925T));
    out.(datest).avgairT      = nanmean(nanmean(out.(datest).airT,3),2);
    out.(datest).stdairT      = nanstd (nanstd (out.(datest).airT,[],3),[],2);
    out.(datest).avghgt       = nanmean(nanmean(out.(datest).geomhgt,3),2);
    out.(datest).stdhgt       = nanstd (nanstd (out.(datest).geomhgt,[],3),[],2);
    out.(datest).avgMMRw      = nanmean(nanmean(out.(datest).watermmr,3),2);
    out.(datest).stdMMRw      = nanstd (nanstd (out.(datest).watermmr,[],3),[],2);
    out.(datest).avgRH        = nanmean(nanmean(out.(datest).RH,3),2);
    out.(datest).stdRH        = nanstd (nanstd (out.(datest).RH,[],3),[],2);
    
    % center Lat/Lon location to save in file name and plot text
    avgLon = round(out.(datest).avgLon*100)/100;
    avgLat = round(out.(datest).avgLat*100)/100;
    stdLon = round(out.(datest).stdLon*100)/100;
    stdLat = round(out.(datest).stdLat*100)/100;
    % mean surface/sea pressure to save in plot text (from Pa to hPa)
    avgsealevelp = round(out.(datest).avgsealevelp/100*1000)/1000;
    stdsealevelp = round(out.(datest).stdsealevelp/100*1000)/1000;
    avgsurfpress = round(out.(datest).avgsurfpress/100*1000)/1000;
    stdsurfpress = round(out.(datest).stdsurfpress/100*1000)/1000;
    
    % convert pressure levels to altitude levels
    % z=pres2alt(t0(K),p0(Pa),plevels(Pa),hgt(m),t(K))
    out.(datest).altLevels = pres2alt(out.(datest).avgairT(1),...
                                                          avgsealevelp*100,...
                                                          out.(datest).plevels*100,...
                                                          out.(datest).avghgt,...
                                                          out.(datest).avgairT);
    if strcmp(plotting,'true')
        
        % plot mean domain temp profiles per day
        figure(i);
        %xerrorbar2(axestype,xmin, xmax, ymin, ymax, x, y, l,u,symbol,teewidth,colorsym)
        %xerrorbar2('loglog',180,400,min(out.(strcat('m',dates{:,i},'_')).plevels),max(out.(strcat('m',dates{:,i},'_')).plevels),...
                   %(out.(strcat('m',dates{:,i},'_')).avgairT),(out.(strcat('m',dates{:,i},'_')).plevels'),...
                   % out.(strcat('m',dates{:,i},'_')).stdairT,out.(strcat('m',dates{:,i},'_')).stdairT,'o-',0.25,ColorSet(i,:));
        plot(out.(datest).avgairT,out.(datest).plevels','-o',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        plot(out.(datest).avgairT + out.(datest).stdairT,...
             out.(datest).plevels',':','color',ColorSet(i,:),'linewidth',2);hold on;
        plot(out.(datest).avgairT - out.(datest).stdairT,...
             out.(datest).plevels',':','color',ColorSet(i,:),'linewidth',2);hold off;
        set(gca,'Ydir','reverse');
        axis([200 300 min(out.(datest).plevels) max(out.(datest).plevels)]);
        t(i) = text(240,50,dates{:,i});
        set(t(i),'FontSize',12,'color',ColorSet(i,:));
        ti(i)=text(220,150,['avgLon=',num2str(avgLon),'±',num2str(stdLon), ' avgLat=',num2str(avgLat),'±',num2str(stdLat)]);
        set(ti(i),'FontSize',12,'color','k');
        ts(i)=text(200,250,['avgSeaSurfP=',num2str(avgsealevelp),'±',num2str(stdsealevelp),...
                                ' avgSurfP=',num2str(avgsurfpress),'±',num2str(stdsurfpress)]);
        set(ts(i),'FontSize',12,'color','k');
        set(gca,'fontsize',12);grid on;
        fi=[strcat(pname, dates{:,i}, 'avgTdomain')];
        save_fig(i,fi,false);
        
        % plot parameters versus pressure levels
        figure(1000);
        ax(1)=subplot(1,3,1);
        plot(out.(datest).avgairT,out.(datest).plevels','-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[300:100:1000]);set(gca,'YTickLabel',{'300','400','500','600','700','800','900','1000'});
        set(gca,'Ydir','reverse');
        title('Temperature','FontSize',11);
        xlabel('K','FontSize',11);ylabel('pressure levels [hPa]','FontSize',11);
        
        set(gca,'Fontsize',11);
        axis([210 280 300 1000]);
        legendall = [legendall;dates{:,i}];
        grid on;
        ax(2)=subplot(1,3,2);
        plot(out.(datest).avgMMRw,out.(datest).plevels','-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[300:100:1000]);set(gca,'YTickLabel',{''});
        %set(gca,'XTick',[1e-4:1e-4:4e-4]);set(gca,'XTickLabel',{'1E^{-4}','2E^{-4}','3E^{-4}','4E^{-4}'});
        set(gca,'Ydir','reverse');
        title('Specific Humidity','FontSize',11);
        xlabel('kg kg^{-1}','FontSize',11);
        grid on;
        set(gca,'Fontsize',11);
        axis([0 max(out.(strcat('m',dates{:,1},'_')).avgMMRw) min(out.(datest).plevels) max(out.(datest).plevels)]);
        ax(3)=subplot(1,3,3);
        plot(out.(datest).avgRH,out.(datest).plevels','-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[300:100:1000]);set(gca,'YTickLabel',{''});
        set(gca,'XTick',[0:20:100]);set(gca,'XTickLabel',{'0','20','40','60','80','100'});
        set(gca,'Ydir','reverse');
        title('Relative Humidity','FontSize',11);
        xlabel('%','FontSize',11);
        set(gca,'Fontsize',11);
        grid on;
        axis([0 100 300 1000]);
        
        % plot parameters versus altitude
        figure(1001);
        ax1(1)=subplot(1,3,1);
        plot(out.(datest).avgairT,out.(datest).altLevels/1000,'-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[1:9]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9'});
        title('Temperature','FontSize',11);
        xlabel('K','FontSize',11);ylabel('MERRA [km]','FontSize',11);
        
        set(gca,'Fontsize',11);
        axis([210 280 min(out.(datest).altLevels/1000) 9]);
        legendall1 = [legendall1;dates{:,i}];
        grid on;
        ax1(2)=subplot(1,3,2);
        plot(out.(datest).avgMMRw,out.(datest).altLevels/1000,'-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[1:9]);set(gca,'YTickLabel',{''});
        title('Specific Humidity','FontSize',11);
        xlabel('kg kg^{-1}','FontSize',11);
        grid on;
        set(gca,'Fontsize',11);
        axis([0 max(out.(strcat('m',dates{:,1},'_')).avgMMRw) min(out.(datest).altLevels/1000) 9]);
        ax1(3)=subplot(1,3,3);
        plot(out.(datest).avgRH,out.(datest).altLevels/1000,'-',...
             'color',ColorSet(i,:),'markerfacecolor',ColorSet(i,:), 'linewidth',2);hold on;
        set(gca,'YTick',[1:9]);set(gca,'YTickLabel',{''});
        set(gca,'XTick',[0:20:100]);set(gca,'XTickLabel',{'0','20','40','60','80','100'});
        title('Relative Humidity','FontSize',11);
        xlabel('%','FontSize',11);
        set(gca,'Fontsize',11);
        grid on;
        axis([0 100 min(out.(datest).altLevels/1000) 9]);
        
    end
end
if strcmp(plotting,'true')
    legend(ax(3),legendall,'location',  'EastOutside');
    legend(ax1(3),legendall1,'location','EastOutside');
    %lg=legend('location','EastOutside');
    linkaxes(ax,'y');
    linkaxes(ax1,'y');
    % handle legend example for multiple plots:
    % a = get(gca,'Children');
    % b = get(gca,'Children');
    % h = [b;a];legend(h,'decreasing','increasing','sine2','sine1','location','eastoutside')

    % get subplot locations:
    % [left, bottom, width, height].
     p1 = get(ax(1), 'position');
     p2 = get(ax(2), 'position');
     p3 = get(ax(3), 'position');
     % adjust locations
     p2(1) = p1(1) + p1(3) +p1(3)/8;
     p3(1) = p2(1) + p2(3) +p2(3)/8;
     set(ax(2), 'position', p2)
     set(ax(3), 'position', p3);
     fi=[strcat(pname, 'Sep2014', 'avgTdomain3subplotVsPst')];
     save_fig(1000,fi,true);
     % get subplot locations:
     % [left, bottom, width, height].
     p1a = get(ax1(1), 'position');
     p2a = get(ax1(2), 'position');
     p3a = get(ax1(3), 'position');
     % adjust locations
     p2a(1) = p1a(1) + p1a(3) +p1a(3)/8;
     p3a(1) = p2a(1) + p2a(3) +p2a(3)/8;
     set(ax1(2), 'position', p2a)
     set(ax1(3), 'position', p3a);
     fi=[strcat(pname, 'Sep2014', 'avgTdomain3subplotVsAlt')];
     save_fig(1001,fi,true);
end
%%

%% cluster averaged daily profiles
% kcluster =3;   % 4 stability regimes (by Barton et al., 2012)
% % cluster from surface to 300 mbar (1:21)
% % create data to cluster:
% dat = [];
% for i=1:length(dates)
%         datest = strcat('m',dates{:,i},'_'); 
%         vec = [out.(datest).avgairT(1),...
%                out.(datest).avgRH(1),...
%                out.(datest).avgst925T,...
%                out.(datest).avgst700T];
%         dat = [dat; vec];
% end
% [IDX C sumd] = kmeans(dat,kcluster,'replicates',50,'start','cluster','MaxIter',200);
% IDXcolor = [0.1 0.8 0.2;0.8 0.2 0.2;0.5 0.3 0.1];
% % plot temperature profile with clustered data
% for i=1:length(dates)
%         datest = strcat('m',dates{:,i},'_'); 
%         figure(1001);
%         if IDX(i)==1
%             plot(out.(datest).avgairT,out.(datest).plevels','-',...
%              'color',IDXcolor(1,:),'markerfacecolor',ColorSet(1,:), 'linewidth',2);hold on;
%         elseif IDX(i)==2
%             plot(out.(datest).avgairT,out.(datest).plevels','-',...
%              'color',IDXcolor(2,:),'markerfacecolor',ColorSet(2,:), 'linewidth',2);hold on;   
%         elseif IDX(i)==3
%              plot(out.(datest).avgairT,out.(datest).plevels','-',...
%              'color',IDXcolor(3,:),'markerfacecolor',ColorSet(3,:), 'linewidth',2);hold on;  
%         end
%         set(gca,'YTick',[300:100:1000]);set(gca,'YTickLabel',{'300','400','500','600','700','800','900','1000'});
%         set(gca,'Ydir','reverse');
%         title('Temperature','FontSize',11);
%         xlabel('K','FontSize',11);ylabel('pressure levels [hPa]','FontSize',11);
%         
%         set(gca,'Fontsize',11);
%         axis([210 280 300 1000]);
%         legendall = [legendall;dates{:,i}];
%         grid on;
% end
%% create temporal averages over each grid-cell

out.Lon = out.(strcat('m',dates{:,1},'_')).longitude;
out.Lat = out.(strcat('m',dates{:,1},'_')).latitude;

for i=1:q     % latitude
    for j=1:p % longitude
        out.airTmat   = [];
        out.potenTmat = [];
        out.st700Tmat = [];
        out.st850Tmat = [];
        out.st925Tmat = [];
        for k=1:length(dates)-1
            out.airTmat   = [out.airTmat   , squeeze(out.(strcat('m',dates{:,k},'_')).airT(  :,i,j))];
            out.potenTmat = [out.potenTmat , squeeze(out.(strcat('m',dates{:,k},'_')).potenT(:,i,j))];
            out.st700Tmat = [out.st700Tmat ,         out.(strcat('m',dates{:,k},'_')).st700T(  i,j)];
            out.st850Tmat = [out.st850Tmat ,         out.(strcat('m',dates{:,k},'_')).st850T(  i,j)];
            out.st925Tmat = [out.st925Tmat ,         out.(strcat('m',dates{:,k},'_')).st925T(  i,j)];
             if k==length(dates)-1
                 out.airTavg(1:s,i,j)   = nanmean(out.airTmat,2);
                 out.airTstd(1:s,i,j)   = nanstd (out.airTmat,[],2);
                 out.potenTavg(1:s,i,j) = nanmean(out.potenTmat,2);
                 out.potenTstd(1:s,i,j) = nanstd (out.potenTmat,[],2);
                 out.st700Tdif(1:k-1,i,j) = diff (out.st700Tmat');
                 out.st850Tdif(1:k-1,i,j) = diff (out.st850Tmat');
                 out.st925Tdif(1:k-1,i,j) = diff (out.st925Tmat');
                 out.st700Tdif2(1:k-2,i,j)= diff(diff   (out.st700Tmat'));
                 out.st850Tdif2(1:k-2,i,j)= diff(diff   (out.st850Tmat'));
                 out.st925Tdif2(1:k-2,i,j)= diff(diff   (out.st925Tmat'));
                 out.st700Tfull(1:k,i,j) =        (out.st700Tmat');
                 out.st850Tfull(1:k,i,j) =        (out.st850Tmat');
                 out.st925Tfull(1:k,i,j) =        (out.st925Tmat');
             end
        end
    end
end

% plot stability trends with time
if strcmp(plotting,'true')
        figure(110)
        for i=1:q
            for j=1:p
            plot(1:30,out.st700Tdif(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{700} time series 1st derivative');
            end
        end
        figure(111)
        for i=1:q
            for j=1:p
            plot(1:31,out.st700Tfull(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{700} time series');
            end
        end
        figure(112)
        for i=1:q
            for j=1:p
            plot(1:29,out.st700Tdif2(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{700} time series 2nd derivative');
            end
        end

        figure(210)
        for i=1:q
            for j=1:p
            plot(1:30,out.st925Tdif(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{925} time series 1st derivative');
            end
        end
        figure(211)
        for i=1:q
            for j=1:p
            plot(1:31,out.st925Tfull(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{925} time series');
            end
        end
        figure(212)
        for i=1:q
            for j=1:p
            plot(1:29,out.st925Tdif2(:,i,j),'-','color',[0.5 0.5 0.5]);hold on;
            ylabel('\Theta_{925} time series 2nd derivative');
            end
        end
end

%% process sea ice concentration

%------------------------------%
% record stability temp movie
if strcmp(plotting,'true')
    % load C130 profiles
    % load star struct
    star = load('F:\ARISE\C-130data\Met\SeaIceProfiles\ARISEairprocessed_with_insitu_woRH_w_anacompare_w_consolidatedclouds20150318.mat');
    nFields    = sum(~structfun(@isempty,star));
    starfieldNames  = fieldnames(star);
        figure;set(gcf, 'Color','white');
        axis tight manual
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        hold all;
        % set(gca, 'nextplot','replacechildren', 'Visible','off');
        loops = length(dates)-1;
        mov(loops) = struct('cdata',[],'colormap',[]);
                for j = 1:loops
                   % plot profile locations 
                   
                   try
                        sname = strcat('met',dates{:,j},'_');
                        aname = strcat('m',  dates{:,j},'_');
                            if length(star.(sname))>=1
                                names  = fieldnames(star.(sname));
                                pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
                                nProfiles = length(pnames);%
                                    if nProfiles > 0
                                            for kk=1:nProfiles
                                                prof = strcat('profnum',num2str(kk));
                                                plot3(star.(sname).(prof).meanlon,star.(sname).(prof).meanlat,150,'s','markersize',12,'markerFaceColor',[1 1 1]);hold on;
                                            end
                                    end
                            end
                    catch
                   end
                    
                    % temp field
                    
                    surf(out.Lon,out.Lat,squeeze(out.st925Tfull(j,:,:)),'edgecolor','none');hold on;  
                                       
                    % sea-ice edge field
                    contour3((amsr.(aname).lon),(amsr.(aname).lat),100+(amsr.(aname).gradline),'k');hold on;
                    
                    % wind fields
                    U = squeeze(out.(aname).eastwind(4,:,:));% 4 for 925mb;7 for 850 mb
                    V = squeeze(out.(aname).northwind(4,:,:));
                    [Lonmat,Latmat] = meshgrid(out.Lon,out.Lat);
                    Unan = isnan(U);
                    Z    = 50*ones(length(out.Lat),length(out.Lon));
                    Z(Unan==1) = NaN;
                    Zval = zeros(length(out.Lat),length(out.Lon));
                    Zval(Unan==1) = NaN;
                    quiver3(Lonmat(~isnan(U)),Latmat(~isnan(U)),Z(~isnan(U)),U(~isnan(U)),V(~isnan(U)),Zval(~isnan(U)),'color',[0.5 0.5 0.5],'LineWidth',0.1);
                    
                   
                    %cb=colorbarlabeled('LTS [K]');
                    cb=colorbarlabeled('LLSS [K]');
                    caxis([0 15]);
                    %caxis([10 30]);
                    %axis([-165 -115 62.5 82.5]);
                    axis([-165 -115 70 80 0 150]);
                    set(gca,'Ytick',70:1:80);
                    set(gca,'YTickLabel',{'70','71','72','73','74','75',...
                                          '76','77','78','79','80'});
                    set(gca,'Xtick',-165:10:-115);
                    set(gca,'XTickLabel',{'-165','-155','-145',...
                                          '-135','-125','-115'});
                    view(2);
                    
                    % add sea-ice line
                    % Plot  
                    hold off;
%                     plot(-66.75,224,'k.','MarkerSize',30); 
                    % add dates
                    % Title: Days 1:31 inclusive. 20140901 to 20141001
                    date = datenum(2014, 09, 01) + j-1;   % Convert t into serial numbers
                    str = datestr(date, 'dd mmm yyyy'); % Show in the format
                    max_ws = sqrt(max(max(U)^2 + max(max(V)^2)));
                    %max_ws = max(max(max(U)),max(max(V)));
                    
                    set(gca,'fontsize',12);
                    title([str,' Max wind speed [m/s] = ',num2str(round(max_ws*100)/100)]);
                    xlabel('Longitude','fontsize',12);
                    ylabel('Latitude','fontsize',12);
                    drawnow
%                     fi = strcat('Theta925_wProfileData_',dates{:,j});
%                     save_fig(11,fi,'true');
                    mov(j) = getframe(gcf);
                    pause(2);
                end
        %  movie(M,n,fps);
        %# save as AVI file, and open it using system video player
        movie2avi(mov, 'LLSS925_wprofilesv2.avi', 'compression','None', 'fps',1);
        winopen('LLSS925_wprofilesv2.avi');
%         movie2avi(mov, 'LTS700.avi', 'compression','None', 'fps',1);
%         winopen('LTS700.avi');
%         movie2avi(mov, 'LLSS925.avi', 'compression','None', 'fps',1);
%         winopen('LLSS925.avi');
end
%movie(F,1,1);
% figure
% Z = peaks;
% surf(Z)
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% 
% loops = 40;
% F(loops) = struct('cdata',[],'colormap',[]);
% for j = 1:loops
%     X = sin(j*pi/10)*Z;
%     surf(X,Z)
%     drawnow
%     F(j) = getframe;
% end
%movie(M,n,fps);
%movie(F,2,1);
%%
out.dates = dates;
%% save processed struct
si = [pname,'MERRAreanalysis_',dates{:,1},'_',dates{:,end},'_','avgLon',num2str(avgLon),'_','avgLat',num2str(avgLat),'.mat'];
disp(['saving to ' si]);
save(si,'-struct','out');
return;

%% function to convert pressure levels to altitude levels
function z=pres2alt(t0,p0,plevels,hgt,t)
% p in Pa
% t,t0 in K
% L - lapse rate in K/m
% p0,t0, pressure and temp at zero altitude
% par = ((plevels/p0)^-LR/g - 1)
% const = t0/L
% z = const*par
R     = 287.053;           % Joule/kg/Kelvin
g     = 9.80665;           % m/s2
L_    = diff(t(1:22))./diff(hgt(1:22));% use surface:10 km
L = mean(L_);
par   = (((plevels/p0).^(-L*R/g))-1);
const = t0/L;
z     = const*par;
return;
