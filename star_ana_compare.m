%% Details of the function:
%  NAME:
% star_ana_compare
%----------------------------------
% PURPOSE:
%  - create atmospheric profiles from ARISE flights-vertical legs
%  - add diagnostics and cloud profiles
%  - plot individual profiles and average them, compare differences
%
% CALLING SEQUENCE:
%  out = star_ana_compare
%
% INPUT:
%  - air.mat structure of C-130 SeaIce profile summary
%  - out.mat structure of reanalysis data
% 
% 
% OUTPUT:
%  saved out struct and plots; daily averages and std's plus stats
%  creates combined profiles for libRadTran
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%  - genCloudProf.m
% 
%
% NEEDED FILES/INPUT:
% 
% - e.g.:F:\ARISE\C-130data\Met\ARISEairprocessed.mat
% - e.g. F:\ARISE\ArcticCRF\METdata\Sep2014\ARISEdomain\...
%           MERRAreanalysis_20140901_20141001_avgLon-140_avgLat72.5.mat
%
% EXAMPLE:
%  - out = star_ana_compare
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Feb-02-2015
% MS, 2015-02-24, added calculations of cloud profiles
% MS, 2015-02-25, added calculation of atm profiles
% -------------------------------------------------------------------------
%% function routine
function out = star_ana_compare
clear all; close all;
startup_plotting; 

%% load data
% load 4star
% this is only summary of meteorological data
% next should add cloud data (flight, satellite etc. in separate routine) and analyze here
% star =  load('F:\ARISE\C-130data\Met\SeaIceProfiles\ARISEairprocessed_wTemp_noRH_noCloud.mat');
arisedir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
filename = 'ARISEairprocessed_with_insitu_woRH';
star     =  load([arisedir filename '.mat']);

% load MERRA
reana = load('F:\ARISE\ArcticCRF\METdata\Sep2014\ARISEdomain\MERRAreanalysis_20140901_20141001_avgLon-140_avgLat72.5.mat');

% prepare for plotting
ColorSet = varycolor(14);% max num of profiles per flight
%

% find number of days analyzed
nFields    = sum(~structfun(@isempty,star));
starfieldNames  = fieldnames(star);
k=0;
for i=1:nFields
    % Extract only the "profnum" fields
    names  = fieldnames(star.(starfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);%max(star.(starfieldNames{i,:}).prof.nprof);
    % load ice-conc data
    daystr=starfieldNames{i,:};
    ice = readAMSR2data({daystr(4:11)});
    legendall ={};
    if nProfiles>0
        dtice       = [];
        dtocean     = [];
        tsurf_ice   = [];
        tsurf_ocean = [];
        for j=1:nProfiles
            profstr    = strcat('profnum',num2str(j));
            %tmp        = starfieldNames{j,:};
            reanadate  = strcat('m',daystr(4:12));
            
            % find closest reanalysis grid
            latintrp= interp1(reana.Lat,[1:length(reana.Lat)],...
                              star.(starfieldNames{i,:}).(profstr).meanlat, 'nearest');
            lonintrp= interp1(reana.Lon,[1:length(reana.Lon)],...
                              star.(starfieldNames{i,:}).(profstr).meanlon, 'nearest');
           
           % save reanalysis products into 4star profiles
           star.(starfieldNames{i,:}).(profstr).ana.alt      = reana.(reanadate).altLevels;
           star.(starfieldNames{i,:}).(profstr).ana.pst      = reana.(reanadate).plevels;
           star.(starfieldNames{i,:}).(profstr).ana.tmp      = reana.(reanadate).airT(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.rh       = reana.(reanadate).RH(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.theta    = reana.(reanadate).potenT(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.h2o      = reana.(reanadate).H2Omolec(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.o3       = reana.(reanadate).O3molec(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.air      = reana.(reanadate).airmolec(:,latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.theta700 = reana.(reanadate).st700T(latintrp,lonintrp);
           star.(starfieldNames{i,:}).(profstr).ana.theta925 = reana.(reanadate).st925T(latintrp,lonintrp);
           
           
           %% generate atmospheric profiles for RT simulations
           %  only std and reanalysis
            [star.(starfieldNames{i,:}).(profstr).ana.atmfile] = gen_atmos_rt(star.(starfieldNames{i,:}).(profstr)...
                                                                             ,star.(starfieldNames{i,:}).(profstr).ana,0,daystr(4:11),profstr);
%            %  std, reanalysis and c130
%            [star.(starfieldNames{i,:}).(profstr).atmfile]     = gen_atmos_rt(star.(starfieldNames{i,:}).(profstr)...
%                                                                             ,star.(starfieldNames{i,:}).(profstr).ana,1,daystr(4:11),profstr);
           %% find closest iceconc in 4STAR grid
           % create interpolated grid
           F = TriScatteredInterp(ice.(strcat('amsr',daystr(4:11))).longrid(:),...
                                  ice.(strcat('amsr',daystr(4:11))).latgrid(:),...
                                  ice.(strcat('amsr',daystr(4:11))).iceconc(:),'nearest');
           % evaluate at 4star locations
           star.(starfieldNames{i,:}).(profstr).iceconc = ...
                  F(star.(starfieldNames{i,:}).(profstr).meanlon,...
                    star.(starfieldNames{i,:}).(profstr).meanlat);
            
           % compare temp
           figure(10);
           plot(flipud(star.(starfieldNames{i,:}).(profstr).Static_AirT+273),...
                flipud(star.(starfieldNames{i,:}).(profstr).z/1000),':','color',ColorSet(j,:),'linewidth',2);hold on;
           plot(reana.(reanadate).airT(:,latintrp,lonintrp),...
                reana.(reanadate).altLevels/1000,'-','color',ColorSet(j,:),'linewidth',2);hold on;
           xlabel('Temperature [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([230 280 0 7]);
           ha=text(231,6.5,daystr(4:11)) ;set(ha,'fontsize',12,'color','k');
           hb=text(270,6.3,'.. C130')    ;set(hb,'fontsize',12,'color','k');
           hc=text(270,5.8,'- MERRA')    ;set(hc,'fontsize',12,'color','k');
           %legend('C130','MERRA');
           set(gca,'Fontsize',12);
           clear ha hb hc
           %% generate cloud fields for each profile
           doy = datestr2doy(daystr(4:11),'yyyymmdd');
           dprof = j;
           k = k+1;
           allprof = k;
           [star.(starfieldNames{i,:}).(profstr).cld] = genCloudProf(star.(starfieldNames{i,:}).(profstr),dprof,allprof,daystr(4:11),num2str(doy));
%            
           %% compare temp with cloud flags
           figure(101);
           plot(flipud(star.(starfieldNames{i,:}).(profstr).Static_AirT+273),...
                flipud(star.(starfieldNames{i,:}).(profstr).z/1000),':','color',ColorSet(j,:),'linewidth',2);hold on;
           plot(reana.(reanadate).airT(:,latintrp,lonintrp),...
                reana.(reanadate).altLevels/1000,'-','color',ColorSet(j,:),'linewidth',2);hold on;
           if (star.(starfieldNames{i,:}).(profstr).cld.cldnum>0)
           plot(flipud(star.(starfieldNames{i,:}).(profstr).Static_AirT(star.(starfieldNames{i,:}).(profstr).cld.cldflag==1)+273),...
                flipud(star.(starfieldNames{i,:}).(profstr).z(star.(starfieldNames{i,:}).(profstr).cld.cldflag==1)/1000),...
                'o','color',ColorSet(j,:),'markerFaceColor',ColorSet(j,:),'markersize',6);hold on;
           end
           xlabel('Temperature [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[0:0.5:3]);set(gca,'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5','3.0'});
           axis([255 275 0 3]);
           ha=text(256,2.8,daystr(4:11))   ;set(ha,'fontsize',12,'color','k');
           hb=text(271,2.7,'.. C130')      ;set(hb,'fontsize',12,'color','k');
           hc=text(271,2.5,'- MERRA')      ;set(hc,'fontsize',12,'color','k');
           hd=text(271,2.3,'o Cloud flag') ;set(hd,'fontsize',12,'color','k');
           %legend('C130','MERRA');
           set(gca,'Fontsize',12);
           
           %% compare RH/MMR
           
           %% smooth, then compare difference in Temp over ice/ocean
           starx = smooth(star.(starfieldNames{i,:}).(profstr).z/1000,0.05);%lowess used
           stary = smooth(starx,star.(starfieldNames{i,:}).(profstr).Static_AirT+273,0.15);%
          %starp = smooth(starx,star.(starfieldNames{i,:}).(profstr).StaticP,0.15);%
           if length(unique(stary)) < length(stary) % then bin
               [npointsx,starx] = bin_profiles_simple(star.(starfieldNames{i,:}).(profstr).z/1000,...
                                           star.(starfieldNames{i,:}).(profstr).z/1000);            
                                           starx=starx(~isnan(starx));
               [npointsy,stary] = bin_profiles_simple(star.(starfieldNames{i,:}).(profstr).z/1000,...
                                           star.(starfieldNames{i,:}).(profstr).Static_AirT+273);
                                           stary=stary(~isnan(stary));
           end
           if length(starx)>2
               star.(starfieldNames{i,:}).(profstr).t700ind = interp1(starx,[1:length(starx)],2.500,'nearest'); % this is for ~2500 m at 700hPa  
               star.(starfieldNames{i,:}).(profstr).t925ind = interp1(starx,[1:length(starx)],0.700 ,'nearest');% this is for ~700  m at 925hPa    
               % figure;plot(stary,starx,'-r');
               starTinterp=interp1(starx,stary,...
                                   reana.(reanadate).altLevels/1000);
                dt = reana.(reanadate).airT(:,latintrp,lonintrp) - starTinterp';
           else
                dt = NaN(length(reana.(reanadate).altLevels),1);
           end
           if ~isnan(star.(starfieldNames{i,:}).(profstr).iceconc) && star.(starfieldNames{i,:}).(profstr).iceconc > 0
               an = [0.1 0.9 0.9];% sea-ice is cyan
               dtice = [dtice;dt'];
               tsurf_ice = [tsurf_ice;stary(1)];
           else
               an = [0.1 0.1 0.9];% open ocean is blue
               dtocean = [dtocean;dt'];
               tsurf_ocean = [tsurf_ocean;stary(1)];
           end
           
           %% plot Temp differences
           figure(11);
           plot((dt),(reana.(reanadate).altLevels/1000),'-o','color',an,...
                'markerfacecolor',an,'linewidth',0.5,'markersize',4);hold on;
           xlabel('\DeltaT (MERRA-C130) [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([-10 10 0 7]);
           %legend('blue is over open ocean','cyan is over sea ice');
           set(gca,'Fontsize',12);
           
           %% plot RH/MMR differences over ice/ocean
           
           
           %% create simple auxilliary profile dat file 
           % (date profnum, timeat bottom,timeattop)
           if strcmp(star.(starfieldNames{i,:}).(profstr).direction,'descend')
               % profile bottom
               t1 = datevec(star.(starfieldNames{i,:}).(profstr).time(end)/24);
               time1 = strcat(num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)));
               % profile top
               t2 = datevec(star.(starfieldNames{i,:}).(profstr).time(1)/24);
               time2 = strcat(num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)));
           else
               t1 = datevec(star.(starfieldNames{i,:}).(profstr).time(1)/24);
               time1 = strcat(num2str(t1(4)),':',num2str(t1(5)),':',num2str(t1(6)));
               t2 = datevec(star.(starfieldNames{i,:}).(profstr).time(end)/24);
               time2 = strcat(num2str(t2(4)),':',num2str(t2(5)),':',num2str(t2(6)));
           end
           
           filen1=['F:\ARISE\ArcticCRF\METdata\profiles_bottom.txt'];
           disp( ['Writing to file: ' filen1]);
           line1=[{'date ' daystr(4:11) ' profilenum ' num2str(j) ' DOY ' num2str(doy) ' UTtime ' time1}];
           dlmwrite(filen1,line1,'-append','delimiter','');
           clear line1;
           
           filen2=['F:\ARISE\ArcticCRF\METdata\profiles_top.txt'];
           disp( ['Writing to file: ' filen2]);
           line2=[{'date ' daystr(4:11) ' profilenum ' num2str(j) ' DOY ' num2str(doy) ' UTtime ' time2}];
           dlmwrite(filen2,line2,'-append','delimiter','');
           clear line;
           
        end% num of profiles
        %% save figure 10
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'tempCompare')];
        save_fig(10,fi,false);
        close(10);
        
        %% save figure 101
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'tempCompare_w_cldflag')];
        save_fig(101,fi,false);
        close(101);
        %% continue and save figure 11
        %-----------------------------%
            figure(11)
            hold on;
            if size(dtice,1)>1
                boxplot(dtice,  'boxstyle','filled','orientation','horizontal',...
                                'positions',reana.(reanadate).altLevels/1000,'colors',[0.1 0.9 0.9]); hold on;
                set(gca,'YTickLabel',{' '})   % Erase ylabels  
            end
            if size(dtocean,1)>1
                boxplot(dtocean,'boxstyle','filled','orientation','horizontal',...
                                'positions',reana.(reanadate).altLevels/1000,'colors',[0.1 0.1 0.9]);
                %set(gca,'XTickLabel',{' '})  % Erase xlabels   
                set(gca,'YTickLabel',{' '})   % Erase ylabels 
            end
            ytix = {'1','2','3','4','5','6','7','8','9','10'};   %  labels
            ytixloc = [1:10];      %  label locations
            set(gca,'YTickMode','auto','YTickLabel',ytix,'YTick',ytixloc);
            axis([-10 10 0 7]);
            h1=text(-9.5,6.5,daystr(4:11))       ;set(h1,'fontsize',12,'color','k');
            h2=text(-9.5,5.8,'over sea-ice')     ;set(h2,'fontsize',12,'color',[0.1 0.9 0.9]);
            h3=text(-9.5,5.5,'over open ocean');set(h3,'fontsize',12,'color',[0.1 0.1 0.9]);
        fi1=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'deltatemp_ana_4star')];
        save_fig(11,fi1,false);
        close(11);
    %% save derived quantities
    % save profile difference
    star.(starfieldNames{i,:}).dtice  = dtice;
    star.(starfieldNames{i,:}).docean = dtocean;
    % save surface temperatures
    star.(starfieldNames{i,:}).tsurf_ice_mean    = nanmean(tsurf_ice);
    star.(starfieldNames{i,:}).tsurf_ice_std     = nanstd( tsurf_ice);
    star.(starfieldNames{i,:}).tsurf_ocean_mean  = nanmean(tsurf_ocean);
    star.(starfieldNames{i,:}).tsurf_ocean_std   = nanstd( tsurf_ocean);

    end
clear ice;   

end% number of days

%% calculate LTS and LLSS per profile
           % LTS  - lower tropospheric stability (theta700-theta_surface)
           % LLSS - lower level static stability (theta925-theta_surface)
    LTSice_mean    = [];
    LTSice_std     = [];
    LLSSice_mean   = [];
    LTSocean_mean  = [];
    LTSocean_std   = [];
    LLSSocean_mean = [];
    LLSSocean_std  = [];
    anaLTSice_mean    = [];
    anaLTSice_std     = [];
    anaLLSSice_mean   = [];
    anaLTSocean_mean  = [];
    anaLTSocean_std   = [];
    anaLLSSocean_mean = [];
    anaLLSSocean_std  = [];
    label_str = {};
for i=1:nFields
    nProfiles = max(star.(starfieldNames{i,:}).prof.nprof);
    daystr=starfieldNames{i,:};
    
    if nProfiles>0
        LTSice = [];
        LLSSice= [];
        LTSocean = [];
        LLSSocean= [];
        anaLTSice= [];
        anaLLSSice=[];
        anaLTSocean=[];
        anaLLSSocean=[];
        for j=1:nProfiles
            profstr    = strcat('profnum',num2str(j));
            if star.(starfieldNames{i,:}).(profstr).iceconc > 0
                if ~isnan(star.(starfieldNames{i,:}).(profstr).t700ind)
                star.(starfieldNames{i,:}).(profstr).LTS =...
                      star.(starfieldNames{i,:}).(profstr).theta(star.(starfieldNames{i,:}).(profstr).t700ind) - ...
                      star.(starfieldNames{i,:}).tsurf_ice_mean;% these numbers need to be calculated as nearest when staticP is available
                LTSice = [LTSice ;star.(starfieldNames{i,:}).(profstr).LTS ];
                end
                if ~isnan(star.(starfieldNames{i,:}).(profstr).t925ind)
                star.(starfieldNames{i,:}).(profstr).LLSS=...
                      star.(starfieldNames{i,:}).(profstr).theta(star.(starfieldNames{i,:}).(profstr).t925ind) - ...
                      star.(starfieldNames{i,:}).tsurf_ice_mean;% these numbers need to be calculated as nearest when staticP is available
                
                LLSSice= [LLSSice;star.(starfieldNames{i,:}).(profstr).LLSS];
                end
                anaLTSice  = [anaLTSice;  star.(starfieldNames{i,:}).(profstr).ana.theta700];
                anaLLSSice = [anaLLSSice; star.(starfieldNames{i,:}).(profstr).ana.theta925];
                % store reanalysis data
                % ana_theta700
            else
                if ~isnan(star.(starfieldNames{i,:}).(profstr).t700ind)
                star.(starfieldNames{i,:}).(profstr).LTS =...
                      star.(starfieldNames{i,:}).(profstr).theta(star.(starfieldNames{i,:}).(profstr).t700ind) - ...
                      star.(starfieldNames{i,:}).tsurf_ocean_mean;% these numbers need to be calculated as nearest when staticP is available
                      LTSocean = [LTSocean ;star.(starfieldNames{i,:}).(profstr).LTS ];
                end
                if ~isnan(star.(starfieldNames{i,:}).(profstr).t925ind)
                star.(starfieldNames{i,:}).(profstr).LLSS=...
                      star.(starfieldNames{i,:}).(profstr).theta(star.(starfieldNames{i,:}).(profstr).t925ind) - ...
                      star.(starfieldNames{i,:}).tsurf_ocean_mean;% these numbers need to be calculated as nearest when staticP is available
                
                LLSSocean= [LLSSocean;star.(starfieldNames{i,:}).(profstr).LLSS];
                end
                anaLTSocean  = [anaLTSocean;  star.(starfieldNames{i,:}).(profstr).ana.theta700];
                anaLLSSocean = [anaLLSSocean; star.(starfieldNames{i,:}).(profstr).ana.theta925];
            end 
        
        
            
        end% profiles
         % calculate daily mean/std for C130
            LTSice_mean    = [LTSice_mean;nanmean(LTSice)];
            LTSice_std     = [LTSice_std;nanstd(LTSice)];
            LLSSice_mean   = [LLSSice_mean;nanmean(LLSSice)];
            LLSSice_std    = [LLSSocean_std;nanstd(LLSSocean)];
            LTSocean_mean  = [LTSocean_mean;nanmean(LTSocean)];
            LTSocean_std   = [LTSocean_std;nanstd(LTSocean)];
            LLSSocean_mean = [LLSSocean_mean;nanmean(LLSSocean)];
            LLSSocean_std  = [LLSSocean_std;nanstd(LLSSocean)];
          
         % calculate daily mean/std for MERRA
            anaLTSice_mean    = [anaLTSice_mean;nanmean(anaLTSice)];
            anaLTSice_std     = [anaLTSice_std;nanstd(anaLTSice)];
            anaLLSSice_mean   = [anaLLSSice_mean;nanmean(anaLLSSice)];
            anaLLSSice_std    = [anaLLSSocean_std;nanstd(anaLLSSocean)];
            anaLTSocean_mean  = [anaLTSocean_mean;nanmean(anaLTSocean)];
            anaLTSocean_std   = [anaLTSocean_std;nanstd(anaLTSocean)];
            anaLLSSocean_mean = [anaLLSSocean_mean;nanmean(anaLLSSocean)];
            anaLLSSocean_std  = [anaLLSSocean_std;nanstd(anaLLSSocean)];
          
            label_str = [label_str,{daystr(8:11)}];%  label_str{:,3}     
            
    end% if any profiles
end% days

%% create and save figure 13
    % plot comparison for each day
    % LTS
    figure(13);
    % plot C130
    errorbar(1:length(label_str),LTSice_mean  ,LTSice_std,  'p','color'  ,[0.1 0.9 0.8],'markerfacecolor',[0.1 0.9 0.8],'markersize',10);  hold on;
    errorbar(1:length(label_str),LTSocean_mean,LTSocean_std,'p','color'  ,[0.1 0.1 0.8],'markerfacecolor',[0.1 0.1 0.8],'markersize',10);  hold on;
    % plot MERRA
    errorbar(1:length(label_str),anaLTSice_mean  ,anaLTSice_std,  'p','color'  ,[0.1 0.9 0.5],'markerfacecolor',[0.1 0.9 0.5],'markersize',6);  hold on;
    errorbar(1:length(label_str),anaLTSocean_mean,anaLTSocean_std,'p','color'  ,[0.1 0.8 0.1],'markerfacecolor',[0.1 0.8 0.1],'markersize',6);  hold on;
    set(gca,'XTickLabel',{' '})  % Erase xlabels 
    set(gca,'XTick',[1:length(label_str)]);  % Erase xlabels 
    set(gca,'XTickLabel',{label_str{:,1},label_str{:,2},label_str{:,3},label_str{:,4},label_str{:,5},...
                           label_str{:,6},label_str{:,7},label_str{:,8},label_str{:,9},label_str{:,10},...
                           label_str{:,11},label_str{:,12},label_str{:,13}});
    axis([0 14 -10 70]);
    th=rotateticklabel(gca,90);
    % add trend lines
    s1=spline(1:length(label_str),LTSice_mean,1:length(label_str));
    s2=spline(1:length(label_str),LTSocean_mean,1:length(label_str));
    s3=spline(1:length(label_str),anaLTSice_mean,1:length(label_str));
    s4=spline(1:length(label_str),anaLTSocean_mean,1:length(label_str));
    plot(1:length(label_str),s1,':','color'  ,[0.1 0.9 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s2,':','color'  ,[0.1 0.1 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s3,':','color'  ,[0.1 0.9 0.5],'linewidth',1.5);hold on;
    plot(1:length(label_str),s4,':','color'  ,[0.1 0.8 0.1],'linewidth',1.5);hold on;
    
    ylabel('\theta_{700} - LML [K]','fontsize',12);
    h11=text(8,65,'C130 over sea-ice')  ; set(h11,'fontsize',12,'color',[0.1 0.9 0.8]);
    h22=text(8,61,'C130over open ocean'); set(h22,'fontsize',12,'color',[0.1 0.1 0.8]);
    h33=text(8,57,'MERRA over sea-ice')  ;set(h33,'fontsize',12,'color',[0.1 0.9 0.5]);
    h44=text(8,53,'MERRA open ocean');    set(h44,'fontsize',12,'color',[0.1 0.8 0.1]);
    
    fi2=[strcat('F:\ARISE\ArcticCRF\figures\',label_str{:,1},'_',label_str{:,end},'LTScompare')];
    save_fig(13,fi2,false);
    close(13);
    
%% create and save figure 14
    % plot comparison for each day
    % LLSS
    figure(14);
    % plot C130
    errorbar(1:length(label_str),LLSSice_mean  ,LLSSice_std,  'p','color'  ,[0.1 0.9 0.8],'markerfacecolor',[0.1 0.9 0.8],'markersize',10);  hold on;
    errorbar(1:length(label_str),LLSSocean_mean,LLSSocean_std,'p','color'  ,[0.1 0.1 0.8],'markerfacecolor',[0.1 0.1 0.8],'markersize',10);  hold on;
    % plot MERRA
    errorbar(1:length(label_str),anaLLSSice_mean  ,anaLLSSice_std,  'p','color'  ,[0.1 0.9 0.5],'markerfacecolor',[0.1 0.9 0.5],'markersize',6);  hold on;
    errorbar(1:length(label_str),anaLLSSocean_mean,anaLLSSocean_std,'p','color'  ,[0.1 0.8 0.1],'markerfacecolor',[0.1 0.8 0.1],'markersize',6);  hold on;
    set(gca,'XTickLabel',{' '})  % Erase xlabels 
    set(gca,'XTick',[1:length(label_str)]);  % Erase xlabels 
    set(gca,'XTickLabel',{label_str{:,1},label_str{:,2},label_str{:,3},label_str{:,4},label_str{:,5},...
                           label_str{:,6},label_str{:,7},label_str{:,8},label_str{:,9},label_str{:,10},...
                           label_str{:,11},label_str{:,12},label_str{:,13}});
    axis([0 14 -10 70]);
    th=rotateticklabel(gca,90);
    % add trend lines
    s1=spline(1:length(label_str),LLSSice_mean,1:length(label_str));
    s2=spline(1:length(label_str),LLSSocean_mean,1:length(label_str));
    s3=spline(1:length(label_str),anaLLSSice_mean,1:length(label_str));
    s4=spline(1:length(label_str),anaLLSSocean_mean,1:length(label_str));
    plot(1:length(label_str),s1,':','color'  ,[0.1 0.9 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s2,':','color'  ,[0.1 0.1 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s3,':','color'  ,[0.1 0.9 0.5],'linewidth',1.5);hold on;
    plot(1:length(label_str),s4,':','color'  ,[0.1 0.8 0.1],'linewidth',1.5);hold on;
    
    ylabel('\theta_{925} - LML [K]','fontsize',12);
    h11=text(8,65,'C130 over sea-ice')  ; set(h11,'fontsize',12,'color',[0.1 0.9 0.8]);
    h22=text(8,61,'C130over open ocean'); set(h22,'fontsize',12,'color',[0.1 0.1 0.8]);
    h33=text(8,57,'MERRA over sea-ice')  ;set(h33,'fontsize',12,'color',[0.1 0.9 0.5]);
    h44=text(8,53,'MERRA open ocean');    set(h44,'fontsize',12,'color',[0.1 0.8 0.1]);
    
    fi2=[strcat('F:\ARISE\ArcticCRF\figures\',label_str{:,1},'_',label_str{:,end},'LLSScompare')];
    save_fig(14,fi2,false);
    close(14);
%% save processed star struct
file2save=[strcat(arisedir,filename,'_w_anacompare_w_cldflag_andProperties.mat')];
save(file2save,'-struct','star');
%% function to convert temperature to potential temperature
function theta=temp2theta(temp,pst)
% pst is pressure in hPa/mbar
% P0 is 1000mbar/1000hPa
% temp,in K
% R - gas constant
% Cp, specific heat
% const=R/Cp = 0.286 for air
% theta=temp*(P0/pst)^const
% calculate
const = 0.286;           %
P0    =1000;             % hPa
theta = temp.*(P0/pst).^const;
return;