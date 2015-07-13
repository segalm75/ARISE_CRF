% createSS_WFpdf
%
% this routine creates static stability
% and wind fields pdf plots for ARISE
% BAMS paper
%
% created by Michal Segal, 2015-05-21, NASA Ames
%-----------------------------------------------

%% load data tables
arisedir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';

% aircraft
% filename = 'ARISEairprocessed_with_insitu_withWVparams20150318';
% processed aircraft data
filename_orig = 'ARISEairprocessed_with_insitu_withWVparams20150318';
filename      = strcat(arisedir,filename_orig,'_w_anacompare_w_consolidatedclouds20150318.mat');
air           =  load(filename);

% load MERRA
reana = load('F:\ARISE\ArcticCRF\METdata\Sep2014\ARISEdomain\MERRAreanalysis_20140901_20141001_avgLon-140_avgLat72.5.mat');

%% agregate stability and wind fields over water/ice (>15% ice)

% construct domain around water/ice profiles
% use mean and std of averaged profiles over a day

% find number of days analyzed
nFields    = sum(~structfun(@isempty,air));
starfieldNames  = fieldnames(air);
k=0;
kk = 0;
iind = [];
wind = [];
Theta925i = [];
Theta925w = [];
Theta850i = [];
Theta850w = [];

ws925i = [];
ws925w = [];
ws850i = [];
ws850w = [];

wd925i = [];
wd925w = [];
wd850i = [];
wd850w = [];

for i=1:nFields
    % Extract only the "profnum" fields
    names  = fieldnames(air.(starfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);
    
    legendall ={};
    if nProfiles>0
            kk = kk+1;
            iceLat_   = [];
            iceLon_   = [];  
            waterLat_ = [];
            waterLon_ = [];
            
            % calculate wind field parameters
            % 925 mb
            % ws = wind speed; wd = wind direction
            % u positive wind towards the east (from west)
            % v positive is wind towards the north (from south)
            ws925 = sqrt(squeeze(reana.(anastr).eastwind(4,:,:)).^2 + squeeze(reana.(anastr).northwind(4,:,:)).^2);
            %wd925 = atan2d(squeeze(reana.(anastr).northwind(4,:,:)),squeeze(reana.(anastr).eastwind(4,:,:)));
            %wd925new = atand(squeeze(reana.(anastr).northwind(4,:,:))./squeeze(reana.(anastr).eastwind(4,:,:)));
            u925 = squeeze(reana.(anastr).eastwind(4,:,:));
            v925 = squeeze(reana.(anastr).northwind(4,:,:));
            wd925 = uv2wd(u925,v925);
            
            % 850 mb
            ws850 = sqrt(squeeze(reana.(anastr).eastwind(7,:,:)).^2 + squeeze(reana.(anastr).northwind(7,:,:)).^2);
            %wd850 = atan2d(squeeze(reana.(anastr).northwind(7,:,:)),squeeze(reana.(anastr).eastwind(7,:,:)));
            %wd850p = wd850;
            %wd850p(wd850<0) = wd850(wd850<0) + 360;
            u850 = squeeze(reana.(anastr).eastwind(7,:,:));
            v850 = squeeze(reana.(anastr).northwind(7,:,:));
            wd850 = uv2wd(u850,v850);
            
         for j=1:nProfiles
             profstr = pnames{j,:};
             if air.(starfieldNames{i,:}).(profstr).iceconc > 15 % i.e. ice
                 iceLat_   = [iceLat_ air.(starfieldNames{i,:}).(profstr).meanlat];
                 iceLon_   = [iceLon_ air.(starfieldNames{i,:}).(profstr).meanlon];
             else
                 waterLat_ = [waterLat_ air.(starfieldNames{i,:}).(profstr).meanlat];
                 waterLon_ = [waterLon_ air.(starfieldNames{i,:}).(profstr).meanlon];
             end

         end
        
         % find boundaries for each of the profile days
         % ice
         anastr = strcat('m',starfieldNames{i,:}(4:end));
         
         max_iLat = max(iceLat_); min_iLat = min(iceLat_);
         max_iLon = max(iceLon_); min_iLon = min(iceLon_);
         if ~isempty(max_iLat|max_iLon|min_iLat|min_iLon)
             
             iLat_in = find(reana.(anastr).latitude<=max_iLat & reana.(anastr).latitude>=min_iLat);
             iLon_in = find(reana.(anastr).longitude<=max_iLon & reana.(anastr).longitude>=min_iLon);
              
                     if isempty(iLat_in)
                        iLat_in = find(reana.(anastr).latitude<=max_iLat +...
                                       0.25 & reana.(anastr).latitude>=min_iLat - 0.25);
                     end
                    
                     if isempty(iLon_in)
                        iLon_in = find(reana.(anastr).longitude<=max_iLon +...
                                       0.25 & reana.(anastr).longitude>=min_iLon - 0.25);
                     end
                     
                     % aggregate stability values
                     % 925 mb
                      Theta925i_ = reana.(anastr).st925T(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      Theta925i = [Theta925i;Theta925i_(:)];
                      % 850 mb
                      Theta850i_ = reana.(anastr).st850T(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      Theta850i = [Theta850i;Theta850i_(:)];
                      
                      if ~isempty(Theta925i_)
                            iind = [iind;i*ones(length(Theta925i_(:)),1)];
                      end
                      
                     % aggregate wind field values over ice
                     % 925 mb
                      
                      ws925i_ = ws925(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      ws925i = [ws925i;ws925i_(:)];
                      
                      wd925i_ = wd925(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      wd925i = [wd925i;wd925i_(:)];
                      
                      % 850 mb
                      ws850i_ = ws850(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      ws850i = [ws850i;ws850i_(:)];
                      
                      wd850i_ = wd850(iLat_in(1):iLat_in(end),...
                                                         iLon_in(1):iLon_in(end));
                      wd850i = [wd850i;wd850i_(:)];
                     
                     
          end
         % water
         max_wLat = max(waterLat_); min_wLat = min(waterLat_);
         max_wLon = max(waterLon_); min_wLon = min(waterLon_);
         if ~isempty(max_wLat|max_wLon|min_wLat|min_wLon)
             
                    wLat_in = find(reana.(anastr).latitude<=max_wLat & reana.(anastr).latitude>=min_wLat);
                    wLon_in = find(reana.(anastr).longitude<=max_wLon & reana.(anastr).longitude>=min_wLon);

                     if isempty(wLat_in)
                        wLat_in = find(reana.(anastr).latitude<=max_wLat +...
                                       0.25 & reana.(anastr).latitude>=min_wLat - 0.25);
                     end

                     if isempty(wLon_in)
                        wLon_in = find(reana.(anastr).longitude<=max_wLon +...
                                       0.25 & reana.(anastr).longitude>=min_wLon - 0.25);
                     end

                     % aggregate ss values
                     % 925 mb
                      
                      Theta925w_ = reana.(anastr).st925T(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      Theta925w = [Theta925w;Theta925w_(:)];
                      
                      % 850 mb
                      
                      Theta850w_ = reana.(anastr).st850T(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      Theta850w = [Theta850w;Theta850w_(:)];
                      
                      if ~isempty(Theta925w_)
                            wind = [wind;i*ones(length(Theta925w_(:)),1)];
                      end
                      
                     % aggregate wind field values over water
                     % 925 mb
                      
                      ws925w_ = ws925(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      ws925w = [ws925w;ws925w_(:)];
                      
                      wd925w_ = wd925(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      wd925w = [wd925w;wd925w_(:)];
                      
                      % 850 mb
                            
                      ws850w_ = ws850(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      ws850w = [ws850w;ws850w_(:)];
                      
                      wd850w_ = wd850(wLat_in(1):wLat_in(end),...
                                                         wLon_in(1):wLon_in(end));
                      wd850w = [wd850w;wd850w_(:)];
        


         end
         
    end
end

%% plot Static stability PDF's
% perform histogram bin count
bincounts925w = histc(Theta925w,-1:16);
bincounts925i = histc(Theta925i,-1:16);
bincounts925w_0919 = histc(Theta925w(wind==14),-1:16);
bincounts925i_0919 = histc(Theta925i(iind==14),-1:16);
bincounts850w = histc(Theta850w,-1:16);
bincounts850i = histc(Theta850i,-1:16);
bincounts850w_0919 = histc(Theta850w(wind==14),-1:16);
bincounts850i_0919 = histc(Theta850i(iind==14),-1:16);

stability_w = [100*bincounts925w/sum(bincounts925w) 100*bincounts850w/sum(bincounts925w)...
               100*bincounts925w_0919/sum(bincounts925w) 100*bincounts850w_0919/sum(bincounts925w)];
stability_i = [100*bincounts925i/sum(bincounts925i) 100*bincounts850i/sum(bincounts925i)...
               100*bincounts925i_0919/sum(bincounts925i) 100*bincounts850i_0919/sum(bincounts925i)];

figure(101);
ax(1)=subplot(211);
h1 = bar(stability_w,1.5,'hist');
set(h1(1),'Facecolor',[0.2 0.2 0.6],'EdgeColor','none');
set(h1(2),'Facecolor',[0.2 0.2 0.9],'EdgeColor','none');
set(h1(3),'Facecolor',[0.6 0.2 0.6],'EdgeColor','k');
set(h1(4),'Facecolor',[0.9 0.2 0.9],'EdgeColor','k');
legend('SS (925 mb)','SS (850 mb)','SS (925 mb) -Sep-19','SS (850 mb) -Sep-19','Location','Best');
axis([0 18 0 40]);ylabel('normalized frequency (%)');
set(ax(1),'Xtick',[0:18]);
set(ax(1),'Xticklabel',{''});
ax(2)=subplot(212);
h2 = bar(stability_i,1.5,'hist');
set(h2(1),'Facecolor',[0.2 0.6 0.6],'EdgeColor','none');
set(h2(2),'Facecolor',[0.2 0.9 0.9],'EdgeColor','none');
set(h2(3),'Facecolor',[0.4 0.6 0.6],'EdgeColor','k');
set(h2(4),'Facecolor',[0.4 0.9 0.9],'EdgeColor','k');
legend('SS (925 mb)','SS (850 mb)','SS (925 mb) -Sep-19','SS (850 mb) -Sep-19','Location','NorthWest');
set(ax(2),'Xtick',[0:18]);
set(ax(2),'Xticklabel',{'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18'});
axis([0 18 0 50]);ylabel('normalized frequency (%)');
legend('boxoff');
xlabel('Temperature [degC]');
set(gcf,'color','white');
% get subplot locations:
    % [left, bottom, width, height].
     p1 = get(ax(1), 'position');p1(4)=0.38;p1(2)=0.50;
     p2 = get(ax(2), 'position');p1(4)=0.38;
     set(ax(1),'position',p1);
     set(ax(2),'position',p2);
fid = strcat(arisedir,'StaticStabilityMERRA_arise_stats');
save_fig(101,fid,'true');

%% plot wind roses

% #Set options in a structure:
Options.AngleNorth     = 0;
Options.AngleEast      = 90;
Options.Labels         = {'N (0)','S (180)','E (90)','W (270)'};
Options.FreqLabelAngle = 45;
Options.lablegend      = '';
Options.nDirections    = 36;%8; %(N, NE, E, SE, S, SW, W and NW)
Options.CenteredIn0    = 'true'; % bin not centered in 0 direction
% Options.nSpeeds        = 8;%6;
Options.vWinds       = [0:2:14];% set intensity range once you
% know it
Options.MaxFrequency   = 30;%50;%60;% upper limit for normalized frequency (60 was used for overall)
Options.nFreq          = 6;%5;%6;% number of frequency lines on circle (5 was used for overall plot)
Options.min_radius     = 0;% no hole in the middle

% create 4 subplots - overall campaign
figure(1021);
% ha = tight_subplot2(Nh, Nw, gap, marg_h, marg_w);
ha = tight_subplot2(2, 2, 0.1, 0.1, 0.1);
set(gcf,'color','white');
% axes1 = subplot(2,2,1);
% axes2 = subplot(2,2,2);
% axes3 = subplot(2,2,3);
% axes4 = subplot(2,2,4);
% set(gcf,'units','normalized','position',[0 0 1 1]);

% water - 850 mb
Options.axes = ha(1);%axes1;
Options.TitleString = {'850 mb winds over water'};
dir = wd850w;
spd = ws850w;
[ax1,count1,speeds1,directions1,Table1] = WindRose((270-dir),spd,Options);

% ice - 850 mb
Options.LegendType = 0;
Options.axes = ha(2);%axes2;
Options.TitleString = {'850 mb winds over ice'};
dir = wd850i;
spd = ws850i;
[ax2,count2,speeds2,directions2,Table2] = WindRose((270-dir),spd,Options);

% water - 925 mb
Options.axes = ha(3);%axes3;
Options.TitleString = {'925 mb winds over water'};
dir = wd925w;
spd = ws925w;
[ax3,count3,speeds3,directions3,Table3] = WindRose((270-dir),spd,Options);

% ice - 925 mb
Options.LegendType = 2;
Options.axes = ha(4);%axes4;
Options.TitleString = {'925 mb winds over ice'};
dir = wd925i;
spd = ws925i;
[ax4,count4,speeds4,directions4,Table4] = WindRose((270-dir),spd,Options);

fid = strcat(arisedir,'WindRoseMERRA_arise_stats_6freqbins_correct_wd_20150527_highres');
save_fig(1021,fid,'true');


% create 4 subplots - 20140919
figure(1031);
% ha = tight_subplot2(Nh, Nw, gap, marg_h, marg_w);
ha = tight_subplot2(2, 2, 0.1, 0.1, 0.1);
set(gcf,'color','white');
% axes1 = subplot(2,2,1);
% axes2 = subplot(2,2,2);
% axes3 = subplot(2,2,3);
% axes4 = subplot(2,2,4);
% set(gcf,'units','normalized','position',[0 0 1 1]);
% Options.vWinds       = [0 2.4 4.8 7.2 9.6 14];% set intensity range once you
Options.MaxFrequency   = 100;% upper limit for normalized frequency (60 was used for overall)
Options.nFreq          = 5;% number of frequency lines on circle (5 was used for overall plot)
% water - 850 mb
Options.axes = ha(1);%axes1;
Options.TitleString = {'850 mb winds over water'};
dir = wd850w(wind==14);
spd = ws850w(wind==14);
[ax1,si.count1,si.speeds1,si.directions1,si.Table1] = WindRose((270-dir),spd,Options);

% ice - 850 mb
Options.LegendType = 0;
Options.axes = ha(2);%axes2;
Options.TitleString = {'850 mb winds over ice'};
dir = wd850i(iind==14);
spd = ws850i(iind==14);
[ax2,si.count2,si.speeds2,si.directions2,si.Table2] = WindRose((270-dir),spd,Options);

% water - 925 mb
Options.axes = ha(3);%axes3;
Options.TitleString = {'925 mb winds over water'};
dir = wd925w(iind==14);
spd = ws925w(iind==14);
[ax3,si.count3,si.speeds3,si.directions3,si.Table3] = WindRose((270-dir),spd,Options);

% ice - 925 mb
Options.LegendType = 2;
Options.axes = ha(4);%axes4;
Options.TitleString = {'925 mb winds over ice'};
dir = wd925i(iind==14);
spd = ws925i(iind==14);
[ax4,si.count4,si.speeds4,si.directions4,si.Table4] = WindRose((270-dir),spd,Options);

fid = strcat(arisedir,'WindRoseMERRA_arise_20140919_correct_wd_20150527highres');
save_fig(1031,fid,'true');

%% test plot with wind_rose function

% figure(102);
% % ha = tight_subplot2(Nh, Nw, gap, marg_h, marg_w);
% ax=tight_subplot2(2, 2, 0.1, 0.1, 0.1);
% set(gcf,'color','white');
diWinds = [0:2:14];
pcircles= [10:5:30];
% water - 850 mb
figure(102);
set(gcf,'color','white');
dir = 90- wd850w;
spd = ws850w;
%[HANDLES1,DATA1] = wind_rose(dir,spd,'dtype','meteo','parent',subplot(2,2,1),'di',diWinds,'ci',pcircles,'quad',3);
[HANDLES1,DATA1] = wind_rose(dir,spd,'dtype','meteo','di',diWinds,'ci',pcircles,'quad',3);

% ice - 850 mb
figure(103);
set(gcf,'color','white');

dir = 90- wd850i;
spd = ws850i;
[HANDLES2,DATA2] = wind_rose(dir,spd,'dtype','meteo','di',diWinds,'ci',pcircles,'quad',3);

%[HANDLES2,DATA2] = wind_rose(dir,spd,'dtype','meteo','parent',subplot(2,2,2),'di',diWinds,'ci',pcircles,'quad',3);

% water - 925 mb
figure(104);
set(gcf,'color','white');
dir = 90- wd925w;
spd = ws925w;
[HANDLES3,DATA3] = wind_rose(dir,spd,'dtype','meteo','di',diWinds,'ci',pcircles,'quad',3);

%[HANDLES3,DATA3] = wind_rose(dir,spd,'dtype','meteo','parent',subplot(2,2,3),'di',diWinds,'ci',pcircles,'quad',3);

% ice - 925 mb
figure(105);
set(gcf,'color','white');
dir = 90- wd925i;
spd = ws925i;
[HANDLES4,DATA4] = wind_rose(dir,spd,'dtype','meteo','di',diWinds,'ci',pcircles,'quad',3);

%[HANDLES4,DATA4] = wind_rose(dir,spd,'dtype','meteo','parent',subplot(2,2,4),'di',diWinds,'ci',pcircles,'quad',3);

% set axis before saving
% [left, bottom, width, height].
%      p1 = get(ax(1), 'position');p1(3) = 1;p1(4)=0.5;%p1(4)=0.38;p1(2)=0.50;
%      p2 = get(ax(2), 'position');p2(3) = 1;p2(4)=0.5;%p1(4)=0.38;
%      p3 = get(ax(3), 'position');p3(3) = 1;p3(4)=0.5;%p1(4)=0.38;p1(2)=0.50;
%      p4 = get(ax(4), 'position');p4(3) = 1;p4(4)=0.5;%p1(4)=0.38;
%      set(ax(1),'position',p1);
%      set(ax(2),'position',p2);
%      set(ax(3),'position',p3);
%      set(ax(4),'position',p4);

%% create nice looking subplot

%parameters for figure and panel size
plotheight=15;
plotwidth=15;
subplotsx=2;
subplotsy=2;   
leftedge=0.5;
rightedge=0.5;   
topedge=1;
bottomedge=1;
spacex=0.2;
spacey=0.2;
fontsize=12;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

% loop to create axis
for i=1:subplotsx
    for ii=1:subplotsy

        axx=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        
    end
end

fid = strcat(arisedir,'WindRoseMERRA_arise_wind_rose_correct_wd_20150527');
save_fig(102,fid,'true');

% wind rose plot direction examples
%  intensity = 30*rand(10000,1);
%  direction =   rand(10000,1);
%  figure;
%  wind_rose(direction,intensity);title('East is 0');
%  
%  figure;
%  d = 90 - direction;
%  wind_rose(d,intensity);title('North is 0');
%  
%  figure;
%  d = 90 - direction;
%  wind_rose(direction,intensity,'dtype','meteo');title('meteo direction type');

% these are some tests (see ARISEexploratory analysis presentation from
% 20150527)
% upos_ = logical(u925>0);
%             upos  = zeros(size(u925));
%             uneg_ = logical(u925<0);
%             uneg  = zeros(size(u925));
%             vnew = zeros(size(v925));
%             upos(upos_)  = u925(upos_==1);
%             uneg(uneg_)  = u925(uneg_==1);
%             wstest = sqrt(upos.^2 + vnew.^2);
%             wdtest = uv2wd(upos,vnew);
%             figure;
%             [HANDLES,DATA] = wind_rose((90-wdtest),wstest,'dtype','meteo');title(['u+ only coming from West - meteorological rep ' anastr]);
%             figure;
%             [HANDLES,DATA] = wind_rose((wdtest),wstest);title(['u+ only blowing towards East ' anastr]);
%             
%             figure;
%              [figure_handle,count,speeds,directions,Table] = WindRose(270-wdtest,wstest,'anglenorth',0,'angleeast',90,'labels',{'N (0)','S (180)','E (90)','W (270)'},'freqlabelangle',45);
%    
%             wdtest2 = uv2wd(uneg,vnew);
%             figure;
%             [HANDLES,DATA] = wind_rose((90-wdtest2),wstest,'dtype','meteo');title(['u- only coming from East - meteorological rep ' anastr]);
%             figure;
%             [HANDLES,DATA] = wind_rose((wdtest2),wstest);title(['u- only blowing towards West ' anastr]);
%             
%             figure;
%              [figure_handle,count,speeds,directions,Table] = WindRose(270-wdtest2,wstest,'anglenorth',0,'angleeast',90,'labels',{'N (0)','S (180)','E (90)','W (270)'},'freqlabelangle',45);
%    
%             figure;hist(wd925(:));title(['wd uv2wd distribution ' anastr]);
%             % "simple" wind rose plot
%           figure;
%          [HANDLES1,DATA1] = wind_rose(wd925,ws925,'dtype','meteo');title(['uv2wd function ' anastr]);
%          figure;
%          [HANDLES2,DATA2] = wind_rose(wd925new,ws925,'dtype','meteo');title(['atand function ' anastr]);
         
%             figure;
%             quiver(squeeze(reana.(anastr).northwind(4,:,:)),squeeze(reana.(anastr).eastwind(4,:,:)));
%             figure;
%             [x,y] = meshgrid(reana.(anastr).longitude,reana.(anastr).latitude);
%             quiver(x,y,squeeze(reana.(anastr).northwind(4,:,:)),squeeze(reana.(anastr).eastwind(4,:,:)));
%             wd925p = wd925;
%             wd925p=mod(-90-wd925,360); 
%             figure;
%             subplot(121);
%             hist(ws925(:),[0:14]);title(['wind speed ', anastr]);
%             subplot(122);
%             hist(wd925p(:),[0:360]);title(['wind direction ', anastr]);
            %wd925p(wd925<0) = wd925(wd925<0) + 360;









