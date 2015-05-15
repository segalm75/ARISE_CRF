% readCRFout
% subroutine to plot CRF output
% from arise libRadtran runs
%------------------------------------------
%% Details of the function:
% NAME:
%   plotariseCRFout
% 
% PURPOSE:
%   read e.g. cases_MERRABASE__MERRABASE_noFog__MERRABASE_wAlb_test.mat
%   plot surface/aircraft/TOA levels for SW/LW/total and other plots
% 
%
% CALLING SEQUENCE:
%  function [d] = plotariseCRFout
%
% INPUT:
%  - 
% 
% 
% OUTPUT:
%  - d processed struct: clear/cloudt TOA and surface quantities
%
% DEPENDENCIES: 
%
% NEEDED FILES:
%  - .mat structure processed by readariseCRFout.m
%
% EXAMPLE:
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer, NASA Ames, MAr, 16, 2015
%
% -------------------------------------------------------------------------
function d = plotariseCRFout
clear all;close all;
%% read in libradtran output struct file

filein = 'F:\ARISE\ArcticCRF\libRadTran_output\arise\cases_MERRABASE__MERRABASE_noFog__MERRABASE_wAlb__MERRABASE_noFogwAlb.mat';

cases     = {'MERRABASE','MERRABASE_noFog','MERRABASE_wAlb','MERRABASE_noFogwAlb'};% base case; no fog cloud below; weighted albedo
groupcase = {'MERRA','MERRA_noFog','MERRA_wAlb','MERRA_noFogwAlb'};
wln     = {'lw','sw'};
sky     = {'allsky','clear'};

d = load(filein);

d.header = {'lambda_nm','lambda_nm','sza','zout','direct','global','diff_down','diff_up','net_down','sum_up'};
isza = 3;
izout= 4;
idir =5;
iglobal=6;
idiffup=8;
inetdn=9;

% calculate CRE
% for each SW/LW, get net = global (dir+diffdn) - diffup
% then subtract net(allsky) - net(clear)
%
dfieldNames  = fieldnames(d.(wln{:,1}).(sky{:,1}).(cases{:,1}));

% this is already inetdn but just a test
% for i=1:length(cases)
%     for j=1:length(wln)
%         for k=1:length(sky)
%             %% calculate net flux at each level
%                 for kk=1:length(dfieldNames)
%                     name  = dfieldNames{kk,:};
%                     l = size(d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name),2);
%                     if ~isempty(d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name))
%                             %for l = 1:size(d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name),1)
%                             d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name)(:,l+1) = ...
%                             d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name)(:,iglobal) - ...
%                             d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name)(:,idiffup);
% %                           d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name).(d.header{isza}).(l) = ...
% %                           d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name)(l,isza);
% %                           d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name).(d.header{izout}).(l) = ...
% %                           d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(name)(l,izout);
%                             %end
%                     end
%                 end
%         end
%     end
% end

%% calculate CRE for each star profile
for i=1:length(cases)
    for j=1:length(wln)
        
            %% calculate net CRE (allsky-clear)
                for kk=1:length(dfieldNames)
                    name  = dfieldNames{kk,:};
                    l = size(d.(wln{:,j}).allsky.(cases{:,i}).(name),2);
                    if ~isempty(d.(wln{:,j}).allsky.(cases{:,i}).(name))
                            d.(wln{:,j}).allsky.(cases{:,i}).(name)(:,l+1) = ...
                            d.(wln{:,j}).allsky.(cases{:,i}).(name)(:,inetdn) - ...
                            d.(wln{:,j}).clear.(cases{:,i}).(name)(:,inetdn);
                            % add normalized CRE
                            
                    end
                end
        
    end
end
% load star struct
s = load('F:\ARISE\C-130data\Met\SeaIceProfiles\ARISEairprocessed_with_insitu_woRH_w_anacompare_w_consolidatedclouds20150318.mat');

%% save CRE to s
nFields    = sum(~structfun(@isempty,s));
sfieldNames  = fieldnames(s);

% initalize fields

swCREsur = [];swCREair = [];swCREtoa = [];
lwCREsur = [];lwCREair = [];lwCREtoa = [];
totCREsur= [];totCREair= [];totCREtoa= [];
swNETsur = [];swNETair = [];swNETtoa = [];
lwNETsur = [];lwNETair = [];lwNETtoa = [];
totNETsur= [];totNETair= [];totNETtoa= [];
zkmair = []; sza = []; iceconc  = []; 
airLTS = []; airLLSS = []; 
anaLTS = []; anaLLSS = [];
ncld = []; cldbot = []; cldtop = [];
intwc = []; intthick = []; int_lwp = [];
allcld_ncld = [];
date = {}; profile = {};
allcld_date = {}; allcld_profile = {};
swgroup = {}; lwgroup = {}; totgroup = {};
iceconcgroup = {};
allcld_iceconcgroup = {};
intwc_group = {};
intthick_group = {};
cldtop_group = {};
cldbot_group = {};
anaLTS_group = {};
intlwp_group = {};
allcld_cldbot = [];
allcld_cldtop = [];
allcld_cldref = [];
allcld_cldwc  = [];
allcld_cldwp  = [];
allcld_cldthick = [];
allcld_iceconc = [];
anaTprof = []; anaMMRprof = []; 
airTprof = []; airMMRprof = [];
cldProf = [];
cldabove = [];

for i=1:nFields
    % Extract only the "profnum" fields
    names  = fieldnames(s.(sfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);
    daystr=sfieldNames{i,:};
    
    if nProfiles>0
        for j=1:nProfiles
            profstr    = strcat('profnum',num2str(j));
            for kk=1:length(cases)
                if length(profstr)==8
                    pname = strcat('day',daystr(4:11),'profile',profstr(end));
                else
                    pname = strcat('day',daystr(4:11),'profile',profstr(8:end));
                end
                % SW CRE
                if ~isempty(d.sw.allsky.(cases{:,kk}).(pname))
                        % sw net down
                         s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swNet = ...
                         d.sw.allsky.(cases{:,kk}).(pname)(:,inetdn);
                     
                        swNETsur  = [swNETsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swNet(1)];
                        swNETair  = [swNETair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swNet(2)];
                        swNETtoa  = [swNETtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swNet(3)];
                        
                        % cre
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swCRF  = ...
                        d.sw.allsky.(cases{:,kk}).(pname)(:,end);
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).szaCRF = ...
                        d.sw.allsky.(cases{:,kk}).(pname)(:,isza);
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).zkmCRF   = ...
                        d.sw.allsky.(cases{:,kk}).(pname)(:,izout);
                    
                        swCREsur = [swCREsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swCRF(1)];
                        swCREair = [swCREair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swCRF(2)];
                        swCREtoa = [swCREtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swCRF(3)];
                        
                        
                        zkmair   = [zkmair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).zkmCRF(2)];
                        sza      = [sza   ;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).szaCRF(2)];
                        swgroup  = [swgroup ;groupcase{:,kk}];
                end
                % LW CRE
                if ~isempty(d.lw.allsky.(cases{:,kk}).(pname))
                    
                         % lw net down
                         s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwNet = ...
                         d.lw.allsky.(cases{:,kk}).(pname)(:,inetdn);
                     
                        lwNETsur  = [lwNETsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwNet(1)];
                        lwNETair  = [lwNETair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwNet(2)];
                        lwNETtoa  = [lwNETtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwNet(3)];
                        
                        % cre
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwCRF  = ...
                        d.lw.allsky.(cases{:,kk}).(pname)(:,end);
                    
                        lwCREsur = [lwCREsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwCRF(1)];
                        lwCREair = [lwCREair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwCRF(2)];
                        lwCREtoa = [lwCREtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwCRF(3)];
                        
                        lwgroup  = [lwgroup ;groupcase{:,kk}];
                end
                % total CRE
                if ~isempty(d.sw.allsky.(cases{:,kk}).(pname))&& ~isempty(d.lw.allsky.(cases{:,kk}).(pname))
                    
                    
                        % total net down
                         s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).Net = ...
                         s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swNet + ...
                         s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwNet;
                     
                        totNETsur  = [totNETsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).Net(1)];
                        totNETair  = [totNETair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).Net(2)];
                        totNETtoa  = [totNETtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).Net(3)];
                     
                        % cre
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).totCRF  = ...
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).swCRF  + ...
                        s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).lwCRF;
                    
                        totCREsur = [totCREsur;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).totCRF(1)];
                        totCREair = [totCREair;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).totCRF(2)];
                        totCREtoa = [totCREtoa;s.(sfieldNames{i,:}).(profstr).(cases{:,kk}).totCRF(3)];
                        
                        totgroup  = [totgroup ;groupcase{:,kk}];
                end
                
            end
            
                %% save other relevant parameters
                iceconc = [iceconc; s.(sfieldNames{i,:}).(profstr).iceconc];
                anaLTS  = [anaLTS;  s.(sfieldNames{i,:}).(profstr).ana.theta700 ];
                anaLLSS = [anaLLSS; s.(sfieldNames{i,:}).(profstr).ana.theta925 ];
                ncld    = [ncld;    s.(sfieldNames{i,:}).(profstr).cld.cldnum   ];
                cldabove= [cldabove;s.(sfieldNames{i,:}).(profstr).cld.cldabove ];
                intwc_tmp = sum(s.(sfieldNames{i,:}).(profstr).cld.cldwc(1:end));
                intwc   = [intwc;   intwc_tmp];
                intthick_tmp = sum(s.(sfieldNames{i,:}).(profstr).cld.cldthick(1:end));
                intthick= [intthick;intthick_tmp];
                int_lwp_tmp = intwc_tmp*intthick_tmp;
                int_lwp  = [int_lwp;int_lwp_tmp];
                if ~isempty((s.(sfieldNames{i,:}).(profstr).cld.cldbot))
                    cldbot_tmp =  min(s.(sfieldNames{i,:}).(profstr).cld.cldbot);
                    cldbot  = [cldbot;  min(s.(sfieldNames{i,:}).(profstr).cld.cldbot) ];
                else
                    cldbot_tmp = NaN;
                    cldbot  = [cldbot;  NaN ];
                end
                if ~isempty((s.(sfieldNames{i,:}).(profstr).cld.cldtop))
                    cldtop_tmp = max(s.(sfieldNames{i,:}).(profstr).cld.cldtop);
                    cldtop  = [cldtop;  max(s.(sfieldNames{i,:}).(profstr).cld.cldtop) ];
                else
                    cldtop = NaN;
                    cldtop  = [cldtop;  NaN ];
                end
                
                % iceconc group name
                if s.(sfieldNames{i,:}).(profstr).iceconc>=0 && s.(sfieldNames{i,:}).(profstr).iceconc <15
                        iceconcgroup  = [iceconcgroup ;'0-15'];
                elseif s.(sfieldNames{i,:}).(profstr).iceconc>15 && s.(sfieldNames{i,:}).(profstr).iceconc <30
                        iceconcgroup  = [iceconcgroup ;'15-30'];
                elseif s.(sfieldNames{i,:}).(profstr).iceconc>=30 && s.(sfieldNames{i,:}).(profstr).iceconc <50
                         iceconcgroup  = [iceconcgroup ;'30-50'];
                elseif s.(sfieldNames{i,:}).(profstr).iceconc>=50 && s.(sfieldNames{i,:}).(profstr).iceconc <70
                        iceconcgroup  = [iceconcgroup ;'50-70'];
                elseif s.(sfieldNames{i,:}).(profstr).iceconc>=70 && s.(sfieldNames{i,:}).(profstr).iceconc <85
                         iceconcgroup  = [iceconcgroup ;'70-85'];
                elseif s.(sfieldNames{i,:}).(profstr).iceconc>=85 && s.(sfieldNames{i,:}).(profstr).iceconc <=100
                        iceconcgroup  = [iceconcgroup ;'85-100'];
                end
                
                 % integrated wc group name

                if intwc_tmp>=0 && intwc_tmp <0.1
                        intwc_group  = [intwc_group ;'0-0.1'];
                elseif intwc_tmp>0.1 && intwc_tmp <0.4
                        intwc_group  = [intwc_group ;'0.1-0.4'];
                elseif intwc_tmp>=0.4 && intwc_tmp <1
                         intwc_group  = [intwc_group ;'0.4-1.0'];
                elseif intwc_tmp>=1 && intwc_tmp <4
                        intwc_group  = [intwc_group ;'1.0-4.0'];
                elseif intwc_tmp>=4
                         intwc_group  = [intwc_group ;'4.0-'];
                end
                
                % integrated cloud thickness group name
                
                 if intthick_tmp>=0 && intthick_tmp <100
                        intthick_group  = [intthick_group ;'0-100'];
                elseif intthick_tmp>100 && intthick_tmp <300
                        intthick_group  = [intthick_group ;'100-300'];
                elseif intthick_tmp>=300 && intthick_tmp <600
                         intthick_group  = [intthick_group ;'300-600'];
                elseif intthick_tmp>=600 && intthick_tmp <1000
                        intthick_group  = [intthick_group ;'600-1000'];
                elseif intthick_tmp>=1000
                         intthick_group  = [intthick_group ;'1000-'];
                 end
                
                  % max cloud top group name
                
                 if cldtop_tmp>=0 && cldtop_tmp <200
                        cldtop_group  = [cldtop_group ;'0-200'];
                elseif cldtop_tmp>200 && cldtop_tmp <600
                        cldtop_group  = [cldtop_group ;'200-600'];
                elseif cldtop_tmp>=600 && cldtop_tmp <1000
                         cldtop_group  = [cldtop_group ;'600-1000'];
                elseif cldtop_tmp>=1000 && cldtop_tmp <3000
                        cldtop_group  = [cldtop_group ;'1000-3000'];
                elseif cldtop_tmp>=3000
                         cldtop_group  = [cldtop_group ;'3000-'];
                 end
                
                   % min cloud bot group name
                
                 if cldbot_tmp>=0 && cldbot_tmp <100
                        cldbot_group  = [cldbot_group ;'0-100'];
                elseif cldbot_tmp>100 && cldbot_tmp <200
                        cldbot_group  = [cldbot_group ;'100-200'];
                elseif cldbot_tmp>=200 && cldbot_tmp <500
                         cldbot_group  = [cldbot_group ;'200-500'];
                elseif cldbot_tmp>=500
                         cldbot_group  = [cldbot_group ;'500-'];
                 end
                
                 
                  % Lower tropospheric stability anaLTS
                
                 if s.(sfieldNames{i,:}).(profstr).ana.theta700 < 12
                        anaLTS_group  = [anaLTS_group ;'0-12'];
                 elseif  s.(sfieldNames{i,:}).(profstr).ana.theta700 > 12 &&...
                        s.(sfieldNames{i,:}).(profstr).ana.theta700 < 20
                        anaLTS_group  = [anaLTS_group ;'12-20'];
                elseif s.(sfieldNames{i,:}).(profstr).ana.theta700 >=20
                         anaLTS_group  = [anaLTS_group ;'20-30'];
                 end
                
                 
                   % LWP [g/m2]
                
                 if int_lwp_tmp < 40
                        intlwp_group  = [intlwp_group ;'0-40'];
                 elseif  int_lwp_tmp >= 40 &&...
                        int_lwp_tmp < 100
                        intlwp_group  = [intlwp_group ;'40-100'];
                elseif int_lwp_tmp >= 100 &&...
                        int_lwp_tmp < 200
                         intlwp_group  = [intlwp_group ;'100-200'];
                elseif int_lwp_tmp >= 200 &&...
                        int_lwp_tmp < 400
                        intlwp_group  = [intlwp_group ;'200-400'];
                end
                
                date    = [date;    sfieldNames{i,:}(4:11)];
                profile = [profile; profstr];
%                 anaTprof= [anaTprof;   s.(sfieldNames{i,:}).(profstr).ana.tmp];
%                 anaMMRprof= [anaMMRprof;   s.(sfieldNames{i,:}).(profstr).ana.mmr];
%                 airTprof= [airTprof;   s.(sfieldNames{i,:}).(profstr).airTinterp];
%                 airMMRprof= [airMMRprof;   s.(sfieldNames{i,:}).(profstr).airMMRinterp];
                try
                    airLTS  = [airLTS;  s.(sfieldNames{i,:}).(profstr).LTS    ];
                catch
                    airLTS  = [airLTS;  NaN    ];
                end
                try
                    airLLSS = [airLLSS  s.(sfieldNames{i,:}).(profstr).LLSS   ];
                catch
                    airLLSS = [airLLSS  NaN   ];
                end
                
                % create array for all cloud properties
                for ll= 1:s.(sfieldNames{i,:}).(profstr).cld.cldnum 
                     if s.(sfieldNames{i,:}).(profstr).iceconc>=0 && s.(sfieldNames{i,:}).(profstr).iceconc <15
                        allcld_iceconcgroup  = [allcld_iceconcgroup ;'0-15'];
                    elseif s.(sfieldNames{i,:}).(profstr).iceconc>15 && s.(sfieldNames{i,:}).(profstr).iceconc <30
                            allcld_iceconcgroup  = [allcld_iceconcgroup ;'15-30'];
                    elseif s.(sfieldNames{i,:}).(profstr).iceconc>=30 && s.(sfieldNames{i,:}).(profstr).iceconc <50
                             allcld_iceconcgroup  = [allcld_iceconcgroup ;'30-50'];
                    elseif s.(sfieldNames{i,:}).(profstr).iceconc>=50 && s.(sfieldNames{i,:}).(profstr).iceconc <70
                            allcld_iceconcgroup  = [allcld_iceconcgroup ;'50-70'];
                    elseif s.(sfieldNames{i,:}).(profstr).iceconc>=70 && s.(sfieldNames{i,:}).(profstr).iceconc <85
                             allcld_iceconcgroup  = [allcld_iceconcgroup ;'70-85'];
                    elseif s.(sfieldNames{i,:}).(profstr).iceconc>=85 && s.(sfieldNames{i,:}).(profstr).iceconc <=100
                            allcld_iceconcgroup  = [allcld_iceconcgroup ;'85-100'];
                    end
                    
                    if ~isempty((s.(sfieldNames{i,:}).(profstr).cld.cldbot(ll)))
                        allcld_cldbot  = [allcld_cldbot;  s.(sfieldNames{i,:}).(profstr).cld.cldbot(ll) ];
                    else
                        allcld_cldbot  = [allcld_cldbot;  NaN ];
                    end
                    
                    if ~isempty((s.(sfieldNames{i,:}).(profstr).cld.cldtop(ll)))
                        allcld_cldtop  = [allcld_cldtop;  s.(sfieldNames{i,:}).(profstr).cld.cldtop(ll) ];
                    else
                        allcld_cldtop  = [allcld_cldtop;  NaN ];
                    end
                    
                    allcld_iceconc   = [allcld_iceconc; s.(sfieldNames{i,:}).(profstr).iceconc];
                    allcld_cldref    = [allcld_cldref; s.(sfieldNames{i,:}).(profstr).cld.cldreff(ll)];
                    allcld_cldwc     = [allcld_cldwc; s.(sfieldNames{i,:}).(profstr).cld.cldwc(ll)];
                    allcld_cldthick  = [allcld_cldthick; s.(sfieldNames{i,:}).(profstr).cld.cldthick(ll)];
                    allcld_cldwp     = [allcld_cldwp;s.(sfieldNames{i,:}).(profstr).cld.cldwc(ll)*s.(sfieldNames{i,:}).(profstr).cld.cldthick(ll)];
                    allcld_date      = [allcld_date;    sfieldNames{i,:}(4:11)];
                    allcld_profile   = [allcld_profile; profstr];
                    allcld_ncld      = [allcld_ncld; s.(sfieldNames{i,:}).(profstr).cld.cldnum];
                    
                end
            
        end
    end
end

%% plot potential temperature for each day
%xerrorbar2(axestype,xmin, xmax, ymin, ymax, x, y, l,u,symbol,teelength,colorsym)
ColorSet = varycolor(length(sfieldNames));% max num of profiles per flight
legendall={};
for i=1:nFields
    names  = fieldnames(s.(sfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);
        if nProfiles>0
            figure(10);
            xerrorbar2('linlin',250,300, 0, 10, s.(sfieldNames{i,:}).prof.thetamean,s.(sfieldNames{i,:}).prof.zmean/1000,...
                                       s.(sfieldNames{i,:}).prof.thetastd,s.(sfieldNames{i,:}).prof.thetastd,'-',0.01,ColorSet(i,:));hold on;
            legendall = [legendall;sfieldNames{i,:}(4:11)];
        end

end
legend(legendall);
xlabel('\Theta [K]');ylabel('Altitude [km]');
%% plot SW/LW/tot with level
figure(333);
range = {'Surface CRE [W/m^{2}]','C-130 CRE [W/m^{2}]','TOA CRE [W/m^{2}]'};
sw = [swCREsur swCREair swCREtoa];
lw = [lwCREsur lwCREair lwCREtoa];
tot= [totCREsur totCREair totCREtoa];

for i=1:3
    ax(i)=subplot(length(range),1,i);
     boxplot(sw(:,i),  swgroup, 'colors',[0.1 0.3 0.9],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h = findobj(gca,'Tag','Box');
     for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),[0.1 0.3 0.9],'FaceAlpha',.5);
     end
     boxplot(lw(:,i),  lwgroup, 'colors',[0.8 0.3 0.2],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h1 = findobj(gca,'Tag','Box','-and',{'Color',[0.8 0.3 0.2]});
     for j=1:length(h1)
        patch(get(h1(j),'XData'),get(h1(j),'YData'),[0.8 0.3 0.2],'FaceAlpha',.8);
     end
     boxplot(tot(:,i), totgroup,'colors',[0.2 0.8 0.2],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h2 = findobj(gca,'Tag','Box','-and',{'Color',[0.2 0.8 0.2]});
     for j=1:length(h2)
        patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.2 0.8 0.2],'FaceAlpha',.5);
     end
     hold on;
     line(0:5,[0,0,0,0,0,0],'Color','k','linewidth',2,'LineStyle','--')
     set(gca,'fontsize',10);
     ylabel(range{:,i},'fontsize',10);
     ylim([-300 200]);
     % get subplot locations:
     % [left, bottom, width, height].
     p(i,:) = get(ax(i), 'position');
     if i==1
         ht = text(2.00,150,'SW');
         set(ht,'fontsize',14,'color',[0.1 0.3 0.9]);
         ht1 = text(2.3,150,'LW');
         set(ht1,'fontsize',14,'color',[0.8 0.3 0.2]);
         ht2 = text(2.6,150,'Total');
         set(ht2,'fontsize',14,'color',[0.2 0.8 0.2]);
     end
end

% tighten subplots
% adjust locations
p(1,2) = 0.70;
p(1,4) = 0.26;
%p(2,2) = p(2,2) + 0.05; %p1(3) +p1(3)/8;
p(2,2) = 0.40;
p(2,4) = 0.265;
%p(3,2) = p(3,2) + 0.1; %p1(3) +p1(3)/8;
p(3,2) = 0.10;
p(3,4) = 0.275;
set(ax(1), 'position', p(1,:))
set(ax(2), 'position', p(2,:));
set(ax(3), 'position', p(3,:));

% ytix = {'1','2','3','4','5','6','7','8','9','10'};   %  labels
xtixloc = [1:4];      %  label locations
set(gca,'XTickMode','auto','XTickLabel',groupcase,'XTick',xtixloc);
% 

%% save figure 333
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'CREbyLevel_MERRA20150320')];
        save_fig(333,fi,true);
        
%%       
figure(3331);
range = {'Surface Flux [W/m^{2}]','C-130 Flux [W/m^{2}]','TOA Flux [W/m^{2}]'};
swNet = [swNETsur swNETair swNETtoa];
lwNet = [lwNETsur lwNETair lwNETtoa];
totNet= [totNETsur totNETair totNETtoa];

for i=1:3
    ax(i)=subplot(length(range),1,i);
     boxplot(swNet(:,i),  swgroup, 'colors',[0.1 0.3 0.9],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h = findobj(gca,'Tag','Box');
     for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),[0.1 0.3 0.9],'FaceAlpha',.5);
     end
     boxplot(lwNet(:,i),  lwgroup, 'colors',[0.8 0.3 0.2],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h1 = findobj(gca,'Tag','Box','-and',{'Color',[0.8 0.3 0.2]});
     for j=1:length(h1)
        patch(get(h1(j),'XData'),get(h1(j),'YData'),[0.8 0.3 0.2],'FaceAlpha',.8);
     end
     boxplot(totNet(:,i), totgroup,'colors',[0.2 0.8 0.2],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
     h2 = findobj(gca,'Tag','Box','-and',{'Color',[0.2 0.8 0.2]});
     for j=1:length(h2)
        patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.2 0.8 0.2],'FaceAlpha',.5);
     end
     hold on;
     line(0:5,[0,0,0,0,0,0],'Color','k','linewidth',2,'LineStyle','--')
     set(gca,'fontsize',10);
     ylabel(range{:,i},'fontsize',10);
     ylim([-300 300]);
     % get subplot locations:
     % [left, bottom, width, height].
     p(i,:) = get(ax(i), 'position');
     if i==1
         ht = text(2.00,-200,'SW');
         set(ht,'fontsize',14,'color',[0.1 0.3 0.9]);
         ht1 = text(2.3,-200,'LW');
         set(ht1,'fontsize',14,'color',[0.8 0.3 0.2]);
         ht2 = text(2.6,-200,'Total');
         set(ht2,'fontsize',14,'color',[0.2 0.8 0.2]);
     end
end

% tighten subplots
% adjust locations
p(1,2) = 0.70;
p(1,4) = 0.26;
%p(2,2) = p(2,2) + 0.05; %p1(3) +p1(3)/8;
p(2,2) = 0.40;
p(2,4) = 0.265;
%p(3,2) = p(3,2) + 0.1; %p1(3) +p1(3)/8;
p(3,2) = 0.10;
p(3,4) = 0.275;
set(ax(1), 'position', p(1,:))
set(ax(2), 'position', p(2,:));
set(ax(3), 'position', p(3,:));

% ytix = {'1','2','3','4','5','6','7','8','9','10'};   %  labels
xtixloc = [1:4];      %  label locations
set(gca,'XTickMode','auto','XTickLabel',groupcase,'XTick',xtixloc);
% 

%% save figure 3331
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'FluxbyLevel_MERRA20150320')];
        save_fig(3331,fi,true);
        
%% plot SW/LW by various params 
% (surface)
swMERRAsur   = swCREsur(strcmp(swgroup,'MERRA'));
lwMERRAsur   = lwCREsur(strcmp(swgroup,'MERRA'));
totMERRAsur  = totCREsur(strcmp(swgroup,'MERRA'));

swMERRAnfsur   = swCREsur(strcmp(swgroup,'MERRA_noFog'));
lwMERRAnfsur   = lwCREsur(strcmp(swgroup,'MERRA_noFog'));
totMERRAnfsur  = totCREsur(strcmp(swgroup,'MERRA_noFog'));

swMERRAwasur   = swCREsur(strcmp(swgroup,'MERRA_wAlb'));
lwMERRAwasur   = lwCREsur(strcmp(swgroup,'MERRA_wAlb'));
totMERRAwasur  = totCREsur(strcmp(swgroup,'MERRA_wAlb'));

swMERRAnfwasur   = swCREsur(strcmp(swgroup,'MERRA_noFogwAlb'));
lwMERRAnfwasur   = lwCREsur(strcmp(swgroup,'MERRA_noFogwAlb'));
totMERRAnfwasur  = totCREsur(strcmp(swgroup,'MERRA_noFogwAlb'));

% (air)
swMERRAair   = swCREair(strcmp(swgroup,'MERRA'));
lwMERRAair   = lwCREair(strcmp(swgroup,'MERRA'));
totMERRAair  = totCREair(strcmp(swgroup,'MERRA'));

swMERRAnfair   = swCREair(strcmp(swgroup,'MERRA_noFog'));
lwMERRAnfair   = lwCREair(strcmp(swgroup,'MERRA_noFog'));
totMERRAnfair  = totCREair(strcmp(swgroup,'MERRA_noFog'));

swMERRAwaair   = swCREair(strcmp(swgroup,'MERRA_wAlb'));
lwMERRAwaair   = lwCREair(strcmp(swgroup,'MERRA_wAlb'));
totMERRAwaair  = totCREair(strcmp(swgroup,'MERRA_wAlb'));

swMERRAnfwaair   = swCREair(strcmp(swgroup,'MERRA_noFogwAlb'));
lwMERRAnfwaair   = lwCREair(strcmp(swgroup,'MERRA_noFogwAlb'));
totMERRAnfwaair  = totCREair(strcmp(swgroup,'MERRA_noFogwAlb'));

% (toa)
swMERRAtoa   = swCREtoa(strcmp(swgroup,'MERRA'));
lwMERRAtoa   = lwCREtoa(strcmp(swgroup,'MERRA'));
totMERRAtoa  = totCREtoa(strcmp(swgroup,'MERRA'));

swMERRAnftoa   = swCREtoa(strcmp(swgroup,'MERRA_noFog'));
lwMERRAnftoa   = lwCREtoa(strcmp(swgroup,'MERRA_noFog'));
totMERRAnftoa  = totCREtoa(strcmp(swgroup,'MERRA_noFog'));

swMERRAwatoa   = swCREtoa(strcmp(swgroup,'MERRA_wAlb'));
lwMERRAwatoa   = lwCREtoa(strcmp(swgroup,'MERRA_wAlb'));
totMERRAwatoa  = totCREtoa(strcmp(swgroup,'MERRA_wAlb'));

swMERRAnfwatoa   = swCREtoa(strcmp(swgroup,'MERRA_noFogwAlb'));
lwMERRAnfwatoa   = lwCREtoa(strcmp(swgroup,'MERRA_noFogwAlb'));
totMERRAnfwatoa  = totCREtoa(strcmp(swgroup,'MERRA_noFogwAlb'));

%% plot relative change from MERRA base case
% noFog
dswMERRAnfsur = swMERRAnfsur - swMERRAsur;
dlwMERRAnfsur = lwMERRAnfsur - lwMERRAsur;
dtotMERRAnfsur = totMERRAnfsur - totMERRAsur;

dswMERRAnfair = swMERRAnfair - swMERRAair;
dlwMERRAnfair = lwMERRAnfair - lwMERRAair;
dtotMERRAnfair = totMERRAnfair - totMERRAair;

dswMERRAnftoa = swMERRAnftoa - swMERRAtoa;
dlwMERRAnftoa = lwMERRAnftoa - lwMERRAtoa;
dtotMERRAnftoa = totMERRAnftoa - totMERRAtoa;

dsw_nofog = [dswMERRAnfsur dswMERRAnfair dswMERRAnftoa];
dlw_nofog = [dlwMERRAnfsur dlwMERRAnfair dlwMERRAnftoa];
dtot_nofog = [dtotMERRAnfsur dtotMERRAnfair dtotMERRAnftoa];

%wAlb
dswMERRAwasur = swMERRAwasur - swMERRAsur;
dlwMERRAwasur = lwMERRAwasur - lwMERRAsur;
dtotMERRAwasur = totMERRAwasur - totMERRAsur;

dswMERRAwaair = swMERRAwaair - swMERRAair;
dlwMERRAwaair = lwMERRAwaair - lwMERRAair;
dtotMERRAwaair = totMERRAwaair - totMERRAair;

dswMERRAwatoa = swMERRAwatoa - swMERRAtoa;
dlwMERRAwatoa = lwMERRAwatoa - lwMERRAtoa;
dtotMERRAwatoa = totMERRAwatoa - totMERRAtoa;

dsw_walb = [dswMERRAwasur dswMERRAwaair dswMERRAwatoa];
dlw_walb = [dlwMERRAwasur dlwMERRAwaair dlwMERRAwatoa];
dtot_walb = [dtotMERRAwasur dtotMERRAwaair dtotMERRAwatoa];

% delta_walb - merrabase by ice conc
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    dsw_wa_byice.(lab)  = dsw_walb(strcmp(iceconcgroup,ice(i)),:);
    dlw_wa_byice.(lab)  = dlw_walb(strcmp(iceconcgroup,ice(i)),:);
    dtot_wa_byice.(lab) = dtot_walb(strcmp(iceconcgroup,ice(i)),:);
end

% sw_walb_byice = {dsw_wa_byice.ice0_15;dsw_wa_byice.ice15_30;dsw_wa_byice.ice30_50;dsw_wa_byice.ice50_70;dsw_wa_byice.ice70_85;dsw_wa_byice.ice85_100}; % Create a cell array with the data for each group
% lw_walb_byice = {dlw_wa_byice.ice0_15;dlw_wa_byice.ice15_30;dlw_wa_byice.ice30_50;dlw_wa_byice.ice50_70;dlw_wa_byice.ice70_85;dlw_wa_byice.ice85_100}; 
% tot_walb_byice = {dtot_wa_byice.ice0_15;dtot_wa_byice.ice15_30;dtot_wa_byice.ice30_50;dtot_wa_byice.ice50_70;dtot_wa_byice.ice70_85;dtot_wa_byice.ice85_100}; 

% remove cases that are the same (i.e. 0)
sw_walb_byice = {dsw_wa_byice.ice0_15(dsw_wa_byice.ice0_15(:,1)~=0,:);dsw_wa_byice.ice15_30(dsw_wa_byice.ice15_30(:,1)~=0,:);dsw_wa_byice.ice30_50(dsw_wa_byice.ice30_50(:,1)~=0,:);...
                 dsw_wa_byice.ice50_70(dsw_wa_byice.ice50_70(:,1)~=0,:);dsw_wa_byice.ice70_85(dsw_wa_byice.ice70_85(:,1)~=0,:);dsw_wa_byice.ice85_100(dsw_wa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_walb_byice = {dlw_wa_byice.ice0_15(dlw_wa_byice.ice0_15(:,1)~=0,:);dlw_wa_byice.ice15_30(dlw_wa_byice.ice15_30(:,1)~=0,:);dlw_wa_byice.ice30_50(dlw_wa_byice.ice30_50(:,1)~=0,:);...
                 dlw_wa_byice.ice50_70(dlw_wa_byice.ice50_70(:,1)~=0,:);dlw_wa_byice.ice70_85(dlw_wa_byice.ice70_85(:,1)~=0,:);dlw_wa_byice.ice85_100(dlw_wa_byice.ice85_100(:,1)~=0,:)}; 
tot_walb_byice = {dtot_wa_byice.ice0_15(dtot_wa_byice.ice0_15(:,1)~=0,:);dtot_wa_byice.ice15_30(dtot_wa_byice.ice15_30(:,1)~=0,:);dtot_wa_byice.ice30_50(dtot_wa_byice.ice30_50(:,1)~=0,:);...
                  dtot_wa_byice.ice50_70(dtot_wa_byice.ice50_70(:,1)~=0,:);dtot_wa_byice.ice70_85(dtot_wa_byice.ice70_85(:,1)~=0,:);dtot_wa_byice.ice85_100(dtot_wa_byice.ice85_100(:,1)~=0,:)}; 
%lw_walb_byice = {dlw_wa_byice.ice0_15;dlw_wa_byice.ice15_50;dlw_wa_byice.ice50_85;dlw_wa_byice.ice85_100}; % Create a cell array with the data for each group
%tot_walb_byice= {dtot_wa_byice.ice0_15;dtot_wa_byice.ice15_50;dtot_wa_byice.ice50_85;dtot_wa_byice.ice85_100}; % Create a cell array with the data for each group

%% delta_walb - merrabase by intwc/intThick/top/bot
% sw/lw/tot
param = {'intwc_group','intthick_group','cldtop_group','cldbot_group','anaLTS_group','intlwp_group'};
pre.par1 = unique(intwc_group);
pre.par2_ = unique(intthick_group);
pre.par2 = pre.par2_; pre.par2{3,:} = pre.par2_{4,:};
pre.par2{4,:} = pre.par2_{5,:};pre.par2{5,:} = pre.par2_{3,:};
pre.par3_ = unique(cldtop_group);
pre.par3 = pre.par3_; pre.par3{2,:} = pre.par3_{3,:};pre.par3{3,:} = pre.par3_{5,:};
pre.par3{4,:} = pre.par3_{2,:};pre.par3{5,:} = pre.par3_{4,:};
pre.par4 = unique(cldbot_group);
pre.par5 = unique(anaLTS_group);
pre.par6_ = unique(intlwp_group);
pre.par6 = pre.par6_; pre.par6{2,:} = pre.par6_{4,:};pre.par6{3,:} = pre.par6_{2,:};pre.par6{4,:} = pre.par6_{3,:};

sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};
% tmp_ = 
% 
%     [18x3 double]
%     [30x3 double]
%      
% tmp=[tmp_(:,end);a]
% 
% tmp = 
% 
%     [18x3 double]
%     [30x3 double]
%     [ 5x3 double]
    
for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             dsw_wa_bywc.(lab)  = dsw_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             dlw_wa_bywc.(lab)  = dlw_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             dtot_wa_bywc.(lab) = dtot_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_wa_bywc.(lab)(dsw_wa_bywc.(lab)(:,1)~=0,:)};
                 %sw_plot={sw_plot;dsw_wa_bywc.(lab)(dsw_wa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_wa_bywc.(lab)(dlw_wa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_wa_bywc.(lab)(dtot_wa_bywc.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_wa_bywc.(lab)(dsw_wa_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_wa_bywc.(lab)(dlw_wa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_wa_bywc.(lab)(dtot_wa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 %sw_plot={sw_plot;dsw_wa_bywc.(lab)(dsw_wa_bywc.(lab)(:,1)~=0,:)};
                 sw_plot__={dsw_wa_bywc.(lab)(dsw_wa_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_wa_bywc.(lab)(dlw_wa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_wa_bywc.(lab)(dtot_wa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             dsw_wa_bythick.(lab)  = dsw_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             dlw_wa_bythick.(lab)  = dlw_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             dtot_wa_bythick.(lab) = dtot_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_wa_bythick.(lab)(dsw_wa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_wa_bythick.(lab)(dlw_wa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_wa_bythick.(lab)(dtot_wa_bythick.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_wa_bythick.(lab)(dsw_wa_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_wa_bythick.(lab)(dlw_wa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_wa_bythick.(lab)(dtot_wa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_wa_bythick.(lab)(dsw_wa_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_wa_bythick.(lab)(dlw_wa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_wa_bythick.(lab)(dtot_wa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             dsw_wa_byctop.(lab)  = dsw_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dlw_wa_byctop.(lab)  = dlw_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dtot_wa_byctop.(lab) = dtot_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={dsw_wa_byctop.(lab)(dsw_wa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_wa_byctop.(lab)(dlw_wa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_wa_byctop.(lab)(dtot_wa_byctop.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_wa_byctop.(lab)(dsw_wa_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_wa_byctop.(lab)(dlw_wa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_wa_byctop.(lab)(dtot_wa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_wa_byctop.(lab)(dsw_wa_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_wa_byctop.(lab)(dlw_wa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_wa_byctop.(lab)(dtot_wa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             dsw_wa_bycbot.(lab)  = dsw_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dlw_wa_bycbot.(lab)  = dlw_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dtot_wa_bycbot.(lab) = dtot_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_wa_bycbot.(lab)(dsw_wa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_wa_bycbot.(lab)(dlw_wa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_wa_bycbot.(lab)(dtot_wa_bycbot.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_wa_bycbot.(lab)(dsw_wa_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_wa_bycbot.(lab)(dlw_wa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_wa_bycbot.(lab)(dtot_wa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_wa_bycbot.(lab)(dsw_wa_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_wa_bycbot.(lab)(dlw_wa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_wa_bycbot.(lab)(dtot_wa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
         elseif kkk==5
             dsw_wa_byLTS.(lab)  = dsw_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dlw_wa_byLTS.(lab)  = dlw_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dtot_wa_byLTS.(lab) = dtot_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_wa_byLTS.(lab)(dsw_wa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_wa_byLTS.(lab)(dlw_wa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_wa_byLTS.(lab)(dtot_wa_byLTS.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_wa_byLTS.(lab)(dsw_wa_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_wa_byLTS.(lab)(dlw_wa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_wa_byLTS.(lab)(dtot_wa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_wa_byLTS.(lab)(dsw_wa_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_wa_byLTS.(lab)(dlw_wa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_wa_byLTS.(lab)(dtot_wa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        end
    end
    
%% plot boxplot per param
%     sw_plot=sw_plot(2:end,:);
%     lw_plot=lw_plot(2:end,:);
%     tot_plot=tot_plot(2:end,:);
    
    figure(111+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
     elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-100 100]);
    set(gca,'XTickLabel',{' '});title('\Delta wAlb - MERRABASE');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CTB 0-100 [m]','CTB 100-200 [m]','CTB 200-500 [m]','CTB > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-100 100]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-100 100]);  
    
end

%noFogwAlb
dswMERRAnfwasur = swMERRAnfwasur - swMERRAsur;
dlwMERRAnfwasur = lwMERRAnfwasur - lwMERRAsur;
dtotMERRAnfwasur = totMERRAnfwasur - totMERRAsur;

dswMERRAnfwaair = swMERRAnfwaair - swMERRAsur;
dlwMERRAnfwaair = lwMERRAnfwaair - lwMERRAsur;
dtotMERRAnfwaair = totMERRAnfwaair - totMERRAsur;

dswMERRAnfwatoa = swMERRAnfwatoa - swMERRAsur;
dlwMERRAnfwatoa = lwMERRAnfwatoa - lwMERRAsur;
dtotMERRAnfwatoa = totMERRAnfwatoa - totMERRAsur;

dsw_nofog_walb = [dswMERRAnfwasur dswMERRAnfwaair dswMERRAnfwatoa];
dlw_nofog_walb = [dlwMERRAnfwasur dlwMERRAnfwaair dlwMERRAnfwatoa];
dtot_nofog_walb = [dtotMERRAnfwasur dtotMERRAnfwaair dtotMERRAnfwatoa];

%noFogwAlb - MERRAwAlb
dswMERRAfogsur = swMERRAwasur - swMERRAnfwasur;%swMERRAnfwasur - swMERRAwasur;
dlwMERRAfogsur = lwMERRAwasur - lwMERRAnfwasur;%lwMERRAnfwasur - lwMERRAwasur;
dtotMERRAfogsur = totMERRAwasur - totMERRAnfwasur;%totMERRAnfwasur - totMERRAwasur;

dswMERRAfogair = swMERRAwaair - swMERRAnfwaair;%swMERRAnfwaair - swMERRAwaair;
dlwMERRAfogair = lwMERRAwaair - lwMERRAnfwaair;%lwMERRAnfwaair - lwMERRAwaair;
dtotMERRAfogair = totMERRAwaair - totMERRAnfwaair;%totMERRAnfwaair - totMERRAwaair;

dswMERRAfogtoa = swMERRAwatoa - swMERRAnfwatoa;% swMERRAnfwatoa - swMERRAwatoa;
dlwMERRAfogtoa = lwMERRAwatoa - lwMERRAnfwatoa;%lwMERRAnfwatoa - lwMERRAwatoa;
dtotMERRAfogtoa = totMERRAwatoa - totMERRAnfwatoa;%totMERRAnfwatoa - totMERRAwatoa;

dsw_fog = [dswMERRAfogsur dswMERRAfogair dswMERRAfogtoa];
dlw_fog = [dlwMERRAfogsur dlwMERRAfogair dlwMERRAfogtoa];
dtot_fog = [dtotMERRAfogsur dtotMERRAfogair dtotMERRAfogtoa];

% MERRABASE  no delta

swmerra = [swMERRAsur swMERRAair swMERRAtoa];
lwmerra= [lwMERRAsur lwMERRAair lwMERRAtoa];
totmerra = [totMERRAsur totMERRAair totMERRAtoa];

% MERRA wAlb no delta

swmerra_wAlb = [swMERRAwasur swMERRAwaair swMERRAwatoa];
lwmerra_wAlb = [lwMERRAwasur lwMERRAwaair lwMERRAwatoa];
totmerra_wAlb = [totMERRAwasur totMERRAwaair totMERRAwatoa];

%% plot surface SW by date

swMERRAwasur_water = swMERRAwasur(iceconc<15);date_water = date(iceconc<15);dw = unique(date_water);
swMERRAwasur_ice   = swMERRAwasur(iceconc>15);date_ice   = date(iceconc>15);di = unique(date_ice);
lwMERRAwasur_water = lwMERRAwasur(iceconc<15);
lwMERRAwasur_ice   = lwMERRAwasur(iceconc>15);
totMERRAwasur_water = totMERRAwasur(iceconc<15);
totMERRAwasur_ice   = totMERRAwasur(iceconc>15);

cloudnum_water = ncld(iceconc<15);
cloudnum_ice   = ncld(iceconc>15);

figure(3333);
% water
%ax(1)=subplot(211)
boxplot(swMERRAwasur_water,  date_water, 'colors',[0.1 0.3 0.9],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),[0.1 0.3 0.9],'FaceAlpha',.5);
end
hold on;
boxplot(lwMERRAwasur_water,  date_water, 'colors',[0.9 0.3 0.1],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h1 = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h1(j),'XData'),get(h1(j),'YData'),[0.9 0.3 0.1],'FaceAlpha',.5);
end
hold on;
boxplot(totMERRAwasur_water,  date_water, 'colors',[0.1 0.9 0.1],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h2 = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.1 0.9 0.1],'FaceAlpha',.2);
end
ylim([-300 100]);
set(gca,'XTick',[1:length(dw)]);
set(gca,'XTickLabel',{dw{1,:},dw{2,:},dw{3,:},dw{4,:},dw{5,:},...
                           dw{6,:},dw{7,:},dw{8,:},dw{9,:}},'fontsize',12);
set(gca,'linewidth',3);hold on;
line(0:14,zeros(15,1),'Color','k','linewidth',2,'LineStyle','--');
ylabel('Surface CRE over Water [W/m^{2}]','fontsize',12);
         ht = text(8.5,-200,'SW');
         set(ht,'fontsize',14,'color',[0.1 0.3 0.9]);
         ht1 = text(8.5,-225,'LW');
         set(ht1,'fontsize',14,'color',[0.9 0.3 0.1]);
         ht2 = text(8.5,-250,'Total');
         set(ht2,'fontsize',14,'color',[0.3 0.7 0.3]);
%ax(2)=subplot(212)
% ice
%ax(1)=subplot(211)
figure(33333)
boxplot(swMERRAwasur_ice,  date_ice, 'colors',[0.1 0.7 0.7],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),[0.1 0.7 0.7],'FaceAlpha',.5);
end
hold on;
boxplot(lwMERRAwasur_ice,  date_ice, 'colors',[0.8 0.5 0.1],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h1 = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h1(j),'XData'),get(h1(j),'YData'),[0.8 0.5 0.1],'FaceAlpha',.7);
end
hold on;
boxplot(totMERRAwasur_ice,  date_ice, 'colors',[0.1 0.6 0.1],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h2 = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h2(j),'XData'),get(h2(j),'YData'),[0.1 0.6 0.1],'FaceAlpha',.7);
end
ylim([-300 100]);
set(gca,'XTick',[1:length(di)]);
set(gca,'XTickLabel',{di{1,:},di{2,:},di{3,:},di{4,:},di{5,:},...
                           di{6,:},di{7,:},di{8,:},di{9,:},di{10,:},di{11,:}},'fontsize',12);
set(gca,'linewidth',3);hold on;
line(0:14,zeros(15,1),'Color','k','linewidth',2,'LineStyle','--')
ylabel('Surface CRE over Sea-Ice [W/m^{2}]','fontsize',12);
ht = text(10.5,-200,'SW');
         set(ht,'fontsize',14,'color',[0.1 0.7 0.7]);
         ht1 = text(10.5,-225,'LW');
         set(ht1,'fontsize',14,'color',[0.8 0.5 0.1]);
         ht2 = text(10.5,-250,'Total');
         set(ht2,'fontsize',14,'color',[0.1 0.6 0.1]);
         
         
dateice = unique(date);         
figure(2222);
% ice conc with time
boxplot(iceconc,  date, 'colors',[0.1 0.7 0.9],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h4 = findobj(gca,'Tag','Box');
for j=1:length(h4)
   patch(get(h4(j),'XData'),get(h4(j),'YData'),[0.1 0.7 0.9],'FaceAlpha',.8);
end

ylim([-5 105]);
set(gca,'XTick',[1:length(dateice)]);
set(gca,'XTickLabel',{dateice{1,:},dateice{2,:},dateice{3,:},dateice{4,:},dateice{5,:},...
                           dateice{6,:},dateice{7,:},dateice{8,:},dateice{9,:},dateice{10,:},...
                           dateice{11,:},dateice{12,:},dateice{13,:}},'fontsize',8);
set(gca,'linewidth',3);hold on;
xlabel('Flight Date','fontsize',14);
ylabel('% sea ice','fontsize',14);

%% cloud num with water/ice
figure(4444);
% water
ax(1)=subplot(211)
boxplot(cloudnum_water,  date_water, 'colors',[0.1 0.3 0.9],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),[0.1 0.3 0.9],'FaceAlpha',.8);
end
ylabel('# clouds above open water','fontsize',10);
ylim([-1 7]);
set(gca,'XTick',[1:length(dw)]);
set(gca,'XTickLabel',{dw{1,:},dw{2,:},dw{3,:},dw{4,:},dw{5,:},...
                           dw{6,:},dw{7,:},dw{8,:},dw{9,:}},'fontsize',12);% water
ax(2)=subplot(212)
boxplot(cloudnum_ice,  date_ice, 'colors',[0.1 0.9 0.7],'outliersize',4); hold on;set(gca,'XTickLabel',{' '})   % Erase ylabels  
h3 = findobj(gca,'Tag','Box');
for j=1:length(h3)
   patch(get(h3(j),'XData'),get(h3(j),'YData'),[0.1 0.9 0.7],'FaceAlpha',.8);
end
ylabel('# clouds above sea-ice','fontsize',10);
ylim([-1 7]);
set(gca,'XTick',[1:length(di)]);
set(gca,'XTickLabel',{di{1,:},di{2,:},di{3,:},di{4,:},di{5,:},...
                           di{6,:},di{7,:},di{8,:},di{9,:},di{10,:},di{11,:}},'fontsize',12);

%% wAlb by iceconc
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    sw_wa_byice.(lab)  = swmerra_wAlb(strcmp(iceconcgroup,ice(i)),:);
    lw_wa_byice.(lab)= lwmerra_wAlb(strcmp(iceconcgroup,ice(i)),:);
    tot_wa_byice.(lab) = totmerra_wAlb(strcmp(iceconcgroup,ice(i)),:);
end

%
sw_byice_wAlb = {sw_wa_byice.ice0_15(sw_wa_byice.ice0_15(:,1)~=0,:);sw_wa_byice.ice15_30(sw_wa_byice.ice15_30(:,1)~=0,:);sw_wa_byice.ice30_50(sw_wa_byice.ice30_50(:,1)~=0,:);...
                   sw_wa_byice.ice50_70(sw_wa_byice.ice50_70(:,1)~=0,:);sw_wa_byice.ice70_85(sw_wa_byice.ice70_85(:,1)~=0,:);sw_wa_byice.ice85_100(sw_wa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_byice_wAlb = {lw_wa_byice.ice0_15(lw_wa_byice.ice0_15(:,1)~=0,:);lw_wa_byice.ice15_30(lw_wa_byice.ice15_30(:,1)~=0,:);lw_wa_byice.ice30_50(lw_wa_byice.ice30_50(:,1)~=0,:);...
                   lw_wa_byice.ice50_70(lw_wa_byice.ice50_70(:,1)~=0,:);lw_wa_byice.ice70_85(lw_wa_byice.ice70_85(:,1)~=0,:);lw_wa_byice.ice85_100(lw_wa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
tot_byice_wAlb= {tot_wa_byice.ice0_15(tot_wa_byice.ice0_15(:,1)~=0,:);tot_wa_byice.ice15_30(tot_wa_byice.ice15_30(:,1)~=0,:);tot_wa_byice.ice30_50(tot_wa_byice.ice30_50(:,1)~=0,:);...
                   tot_wa_byice.ice50_70(tot_wa_byice.ice50_70(:,1)~=0,:);tot_wa_byice.ice70_85(tot_wa_byice.ice70_85(:,1)~=0,:);tot_wa_byice.ice85_100(tot_wa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group


% wAlb by anaLTS no delta
% sw/lw/tot
for i = 1:length(unique(anaLTS_group))
    LTS = unique(anaLTS_group); cLTS = LTS{i,1};
    LTSlab = strrep(cLTS, '-', '_');
    lab = strcat('anaLTS',LTSlab);
    sw_wa_byLTS.(lab)  = swmerra_wAlb(strcmp(anaLTS_group,LTS(i)),:);
    lw_wa_byLTS.(lab)= lwmerra_wAlb(strcmp(anaLTS_group,LTS(i)),:);
    tot_wa_byLTS.(lab) = totmerra_wAlb(strcmp(anaLTS_group,LTS(i)),:);
end

%
sw_byLTS_wAlb = {sw_wa_byLTS.anaLTS0_12;sw_wa_byLTS.anaLTS12_20;sw_wa_byLTS.anaLTS20_30}; % Create a cell array with the data for each group
lw_byLTS_wAlb = {lw_wa_byLTS.anaLTS0_12;lw_wa_byLTS.anaLTS12_20;lw_wa_byLTS.anaLTS20_30}; % Create a cell array with the data for each group
tot_byLTS_wAlb= {tot_wa_byLTS.anaLTS0_12;tot_wa_byLTS.anaLTS12_20;tot_wa_byLTS.anaLTS20_30}; % Create a cell array with the data for each group
               
               

% delta_nofog_walb - merrabase by ice conc
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    dsw_nfwa_byice.(lab)  = dsw_nofog_walb(strcmp(iceconcgroup,ice(i)),:);
    dlw_nfwa_byice.(lab)  = dlw_nofog_walb(strcmp(iceconcgroup,ice(i)),:);
    dtot_nfwa_byice.(lab) = dtot_nofog_walb(strcmp(iceconcgroup,ice(i)),:);
end

% sw_nfwalb_byice = {dsw_nfwa_byice.ice0_15;dsw_nfwa_byice.ice15_30;dsw_nfwa_byice.ice30_50;dsw_nfwa_byice.ice50_70;dsw_nfwa_byice.ice70_85;dsw_nfwa_byice.ice85_100}; % Create a cell array with the data for each group
% lw_nfwalb_byice = {dlw_nfwa_byice.ice0_15;dlw_nfwa_byice.ice15_30;dlw_nfwa_byice.ice30_50;dlw_nfwa_byice.ice50_70;dlw_nfwa_byice.ice70_85;dlw_nfwa_byice.ice85_100}; % Create a cell array with the data for each group
% tot_nfwalb_byice= {dtot_nfwa_byice.ice0_15;dtot_nfwa_byice.ice15_30;dtot_nfwa_byice.ice30_50;dtot_nfwa_byice.ice50_70;dtot_nfwa_byice.ice70_85;dtot_nfwa_byice.ice85_100}; % Create a cell array with the data for each group

% remove cases that are zero
sw_nfwalb_byice = {dsw_nfwa_byice.ice0_15(dsw_nfwa_byice.ice0_15(:,1)~=0,:);dsw_nfwa_byice.ice15_30(dsw_nfwa_byice.ice15_30(:,1)~=0,:);dsw_nfwa_byice.ice30_50(dsw_nfwa_byice.ice30_50(:,1)~=0,:);...
                   dsw_nfwa_byice.ice50_70(dsw_nfwa_byice.ice50_70(:,1)~=0,:);dsw_nfwa_byice.ice70_85(dsw_nfwa_byice.ice70_85(:,1)~=0,:);dsw_nfwa_byice.ice85_100(dsw_nfwa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_nfwalb_byice = {dlw_nfwa_byice.ice0_15(dlw_nfwa_byice.ice0_15(:,1)~=0,:);dlw_nfwa_byice.ice15_30(dlw_nfwa_byice.ice15_30(:,1)~=0,:);dlw_nfwa_byice.ice30_50(dlw_nfwa_byice.ice30_50(:,1)~=0,:);...
                   dlw_nfwa_byice.ice50_70(dlw_nfwa_byice.ice50_70(:,1)~=0,:);dlw_nfwa_byice.ice70_85(dlw_nfwa_byice.ice70_85(:,1)~=0,:);dlw_nfwa_byice.ice85_100(dlw_nfwa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
tot_nfwalb_byice= {dtot_nfwa_byice.ice0_15(dtot_nfwa_byice.ice0_15(:,1)~=0,:);dtot_nfwa_byice.ice15_30(dtot_nfwa_byice.ice15_30(:,1)~=0,:);dtot_nfwa_byice.ice30_50(dtot_nfwa_byice.ice30_50(:,1)~=0,:);...
                   dtot_nfwa_byice.ice50_70(dtot_nfwa_byice.ice50_70(:,1)~=0,:);dtot_nfwa_byice.ice70_85(dtot_nfwa_byice.ice70_85(:,1)~=0,:);dtot_nfwa_byice.ice85_100(dtot_nfwa_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group

               
% delta_nofog_walb - merrabase_wAlb by ice conc
% this is merrabase_wAlb - nofog_wAlb
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    dsw_fog_byice.(lab)  = dsw_fog(strcmp(iceconcgroup,ice(i)),:);
    dlw_fog_byice.(lab)  = dlw_fog(strcmp(iceconcgroup,ice(i)),:);
    dtot_fog_byice.(lab) = dtot_fog(strcmp(iceconcgroup,ice(i)),:);
end

% sw_nfwalb_byice = {dsw_nfwa_byice.ice0_15;dsw_nfwa_byice.ice15_30;dsw_nfwa_byice.ice30_50;dsw_nfwa_byice.ice50_70;dsw_nfwa_byice.ice70_85;dsw_nfwa_byice.ice85_100}; % Create a cell array with the data for each group
% lw_nfwalb_byice = {dlw_nfwa_byice.ice0_15;dlw_nfwa_byice.ice15_30;dlw_nfwa_byice.ice30_50;dlw_nfwa_byice.ice50_70;dlw_nfwa_byice.ice70_85;dlw_nfwa_byice.ice85_100}; % Create a cell array with the data for each group
% tot_nfwalb_byice= {dtot_nfwa_byice.ice0_15;dtot_nfwa_byice.ice15_30;dtot_nfwa_byice.ice30_50;dtot_nfwa_byice.ice50_70;dtot_nfwa_byice.ice70_85;dtot_nfwa_byice.ice85_100}; % Create a cell array with the data for each group

% remove cases that are zero (i.e. had no fog layers to begin with)
sw_fog_byice = {dsw_fog_byice.ice0_15(dsw_fog_byice.ice0_15(:,1)~=0,:);dsw_fog_byice.ice15_30(dsw_fog_byice.ice15_30(:,1)~=0,:);dsw_fog_byice.ice30_50(dsw_fog_byice.ice30_50(:,1)~=0,:);...
                   dsw_fog_byice.ice50_70(dsw_fog_byice.ice50_70(:,1)~=0,:);dsw_fog_byice.ice70_85(dsw_fog_byice.ice70_85(:,1)~=0,:);dsw_fog_byice.ice85_100(dsw_fog_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_fog_byice = {dlw_fog_byice.ice0_15(dlw_fog_byice.ice0_15(:,1)~=0,:);dlw_fog_byice.ice15_30(dlw_fog_byice.ice15_30(:,1)~=0,:);dlw_fog_byice.ice30_50(dlw_fog_byice.ice30_50(:,1)~=0,:);...
                   dlw_fog_byice.ice50_70(dlw_fog_byice.ice50_70(:,1)~=0,:);dlw_fog_byice.ice70_85(dlw_fog_byice.ice70_85(:,1)~=0,:);dlw_fog_byice.ice85_100(dlw_fog_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
tot_fog_byice= {dtot_fog_byice.ice0_15(dtot_fog_byice.ice0_15(:,1)~=0,:);dtot_fog_byice.ice15_30(dtot_fog_byice.ice15_30(:,1)~=0,:);dtot_fog_byice.ice30_50(dtot_fog_byice.ice30_50(:,1)~=0,:);...
                   dtot_fog_byice.ice50_70(dtot_fog_byice.ice50_70(:,1)~=0,:);dtot_fog_byice.ice70_85(dtot_fog_byice.ice70_85(:,1)~=0,:);dtot_fog_byice.ice85_100(dtot_fog_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group

% this is merrabase_wAlb - nofog_wAlb
% does not include numcloud==1
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    dsw_fogonly_byice.(lab)  = dsw_fog(strcmp(iceconcgroup,ice(i))&ncld>1,:);
    dlw_fogonly_byice.(lab)  = dlw_fog(strcmp(iceconcgroup,ice(i))&ncld>1,:);
    dtot_fogonly_byice.(lab) = dtot_fog(strcmp(iceconcgroup,ice(i))&ncld>1,:);
end

% remove cases that are zero (i.e. had no fog layers to begin with) and
% cases transition from fog to clear
sw_fogonly_byice = {dsw_fogonly_byice.ice0_15(dsw_fogonly_byice.ice0_15(:,1)~=0,:);dsw_fogonly_byice.ice30_50(dsw_fogonly_byice.ice30_50(:,1)~=0,:);...
                   dsw_fogonly_byice.ice70_85(dsw_fogonly_byice.ice70_85(:,1)~=0,:);dsw_fogonly_byice.ice85_100(dsw_fogonly_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_fogonly_byice = {dlw_fogonly_byice.ice0_15(dlw_fogonly_byice.ice0_15(:,1)~=0,:);dlw_fogonly_byice.ice30_50(dlw_fogonly_byice.ice30_50(:,1)~=0,:);...
                   dlw_fogonly_byice.ice70_85(dlw_fogonly_byice.ice70_85(:,1)~=0,:);dlw_fogonly_byice.ice85_100(dlw_fogonly_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
tot_fogonly_byice= {dtot_fogonly_byice.ice0_15(dtot_fogonly_byice.ice0_15(:,1)~=0,:);dtot_fogonly_byice.ice30_50(dtot_fogonly_byice.ice30_50(:,1)~=0,:);...
                   dtot_fogonly_byice.ice70_85(dtot_fogonly_byice.ice70_85(:,1)~=0,:);dtot_fogonly_byice.ice85_100(dtot_fogonly_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group

%% wAlb by iceconc for ncld>1
% sw/lw/tot
for i = 1:length(unique(iceconcgroup))
    ice = unique(iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    sw_wa_mcld_byice.(lab)  = swmerra_wAlb(strcmp(iceconcgroup,ice(i))&ncld>1,:);
    lw_wa_mcld_byice.(lab)= lwmerra_wAlb(strcmp(iceconcgroup,ice(i))&ncld>1,:);
    tot_wa_mcld_byice.(lab) = totmerra_wAlb(strcmp(iceconcgroup,ice(i))&ncld>1,:);
end

%
sw_byice_wAlb_mcld = {sw_wa_mcld_byice.ice0_15(sw_wa_mcld_byice.ice0_15(:,1)~=0,:);sw_wa_mcld_byice.ice30_50(sw_wa_mcld_byice.ice30_50(:,1)~=0,:);...
                   sw_wa_mcld_byice.ice70_85(sw_wa_mcld_byice.ice70_85(:,1)~=0,:);sw_wa_mcld_byice.ice85_100(sw_wa_mcld_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
lw_byice_wAlb_mcld = {lw_wa_mcld_byice.ice0_15(lw_wa_mcld_byice.ice0_15(:,1)~=0,:);lw_wa_mcld_byice.ice30_50(lw_wa_mcld_byice.ice30_50(:,1)~=0,:);...
                  lw_wa_mcld_byice.ice70_85(lw_wa_mcld_byice.ice70_85(:,1)~=0,:);lw_wa_mcld_byice.ice85_100(lw_wa_mcld_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group
tot_byice_wAlb_mcld= {tot_wa_mcld_byice.ice0_15(tot_wa_mcld_byice.ice0_15(:,1)~=0,:);tot_wa_mcld_byice.ice30_50(tot_wa_mcld_byice.ice30_50(:,1)~=0,:);...
                   tot_wa_mcld_byice.ice70_85(tot_wa_mcld_byice.ice70_85(:,1)~=0,:);tot_wa_mcld_byice.ice85_100(tot_wa_mcld_byice.ice85_100(:,1)~=0,:)}; % Create a cell array with the data for each group

               
%% delta_walbnoFog - merrabase by intwc/intThick/top/bot
% sw/lw/tot
sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};

for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             dsw_nfwa_bywc.(lab)  = dsw_nofog_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             dlw_nfwa_bywc.(lab)  = dlw_nofog_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             dtot_nfwa_bywc.(lab) = dtot_nofog_walb(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_nfwa_bywc.(lab)(dsw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nfwa_bywc.(lab)(dlw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nfwa_bywc.(lab)(dtot_nfwa_bywc.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nfwa_bywc.(lab)(dsw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nfwa_bywc.(lab)(dlw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nfwa_bywc.(lab)(dtot_nfwa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nfwa_bywc.(lab)(dsw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nfwa_bywc.(lab)(dlw_nfwa_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nfwa_bywc.(lab)(dtot_nfwa_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             dsw_nfwa_bythick.(lab)  = dsw_nofog_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             dlw_nfwa_bythick.(lab)  = dlw_nofog_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             dtot_nfwa_bythick.(lab) = dtot_nofog_walb(strcmp(intthick_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_nfwa_bythick.(lab)(dsw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nfwa_bythick.(lab)(dlw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nfwa_bythick.(lab)(dtot_nfwa_bythick.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nfwa_bythick.(lab)(dsw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nfwa_bythick.(lab)(dlw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nfwa_bythick.(lab)(dtot_nfwa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nfwa_bythick.(lab)(dsw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nfwa_bythick.(lab)(dlw_nfwa_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nfwa_bythick.(lab)(dtot_nfwa_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             dsw_nfwa_byctop.(lab)  = dsw_nofog_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dlw_nfwa_byctop.(lab)  = dlw_nofog_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dtot_nfwa_byctop.(lab) = dtot_nofog_walb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={dsw_nfwa_byctop.(lab)(dsw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nfwa_byctop.(lab)(dlw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nfwa_byctop.(lab)(dtot_nfwa_byctop.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nfwa_byctop.(lab)(dsw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nfwa_byctop.(lab)(dlw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nfwa_byctop.(lab)(dtot_nfwa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nfwa_byctop.(lab)(dsw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nfwa_byctop.(lab)(dlw_nfwa_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nfwa_byctop.(lab)(dtot_nfwa_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             dsw_nfwa_bycbot.(lab)  = dsw_nofog_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dlw_nfwa_bycbot.(lab)  = dlw_nofog_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dtot_nfwa_bycbot.(lab) = dtot_nofog_walb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_nfwa_bycbot.(lab)(dsw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nfwa_bycbot.(lab)(dlw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nfwa_bycbot.(lab)(dtot_nfwa_bycbot.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nfwa_bycbot.(lab)(dsw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nfwa_bycbot.(lab)(dlw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nfwa_bycbot.(lab)(dtot_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nfwa_bycbot.(lab)(dsw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nfwa_bycbot.(lab)(dlw_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nfwa_bycbot.(lab)(dtot_nfwa_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==5
             dsw_nfwa_byLTS.(lab)  = dsw_nofog_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dlw_nfwa_byLTS.(lab)  = dlw_nofog_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dtot_nfwa_byLTS.(lab) = dtot_nofog_walb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_nfwa_byLTS.(lab)(dsw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nfwa_byLTS.(lab)(dlw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nfwa_byLTS.(lab)(dtot_nfwa_byLTS.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nfwa_byLTS.(lab)(dsw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nfwa_byLTS.(lab)(dlw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nfwa_byLTS.(lab)(dtot_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nfwa_byLTS.(lab)(dsw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nfwa_byLTS.(lab)(dlw_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nfwa_byLTS.(lab)(dtot_nfwa_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        end
    end
    
%% plot boxplot per param
    
    figure(222+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4 
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-200 100]);
    set(gca,'XTickLabel',{' '});title('\Delta wAlbnoFog - MERRABASE');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-100 100]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-200 100]);  
    
end
               
%% delta_noFog - merrabase by intwc/intThick/top/bot
% sw/lw/tot
sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};

for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             dsw_nf_bywc.(lab)  = dsw_nofog(strcmp(intwc_group,pre.(parlab)(i)),:);
             dlw_nf_bywc.(lab)  = dlw_nofog(strcmp(intwc_group,pre.(parlab)(i)),:);
             dtot_nf_bywc.(lab) = dtot_nofog(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_nf_bywc.(lab)(dsw_nf_bywc.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nf_bywc.(lab)(dlw_nf_bywc.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nf_bywc.(lab)(dtot_nf_bywc.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nf_bywc.(lab)(dsw_nf_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nf_bywc.(lab)(dlw_nf_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nf_bywc.(lab)(dtot_nf_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nf_bywc.(lab)(dsw_nf_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nf_bywc.(lab)(dlw_nf_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nf_bywc.(lab)(dtot_nf_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             dsw_nf_bythick.(lab)  = dsw_nofog(strcmp(intthick_group,pre.(parlab)(i)),:);
             dlw_nf_bythick.(lab)  = dlw_nofog(strcmp(intthick_group,pre.(parlab)(i)),:);
             dtot_nf_bythick.(lab) = dtot_nofog(strcmp(intthick_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)) sw_plot_ = {0*zeros(1,3)};end
                 lw_plot_={dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)) lw_plot_ = {0*zeros(1,3)};end
                 tot_plot_={dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)) tot_plot_ = {0*zeros(1,3)};end
             elseif  i==2
                 sw_plot={dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)) sw_plot = {0*zeros(1,3)};end
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)) lw_plot = {0*zeros(1,3)};end
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)) tot_plot = {0*zeros(1,3)};end
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_nf_bythick.(lab)(dsw_nf_bythick.(lab)(:,1)~=0,:)) sw_plot__ = {0*zeros(1,3)};end
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_nf_bythick.(lab)(dlw_nf_bythick.(lab)(:,1)~=0,:)) lw_plot__ = {0*zeros(1,3)};end
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_nf_bythick.(lab)(dtot_nf_bythick.(lab)(:,1)~=0,:)) tot_plot__ = {0*zeros(1,3)};end
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             dsw_nf_byctop.(lab)  = dsw_nofog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dlw_nf_byctop.(lab)  = dlw_nofog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dtot_nf_byctop.(lab) = dtot_nofog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={dsw_nf_byctop.(lab)(dsw_nf_byctop.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nf_byctop.(lab)(dlw_nf_byctop.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nf_byctop.(lab)(dtot_nf_byctop.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nf_byctop.(lab)(dsw_nf_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nf_byctop.(lab)(dlw_nf_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nf_byctop.(lab)(dtot_nf_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nf_byctop.(lab)(dsw_nf_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nf_byctop.(lab)(dlw_nf_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nf_byctop.(lab)(dtot_nf_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             dsw_nf_bycbot.(lab)  = dsw_nofog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dlw_nf_bycbot.(lab)  = dlw_nofog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dtot_nf_bycbot.(lab) = dtot_nofog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_nf_bycbot.(lab)(dsw_nf_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nf_bycbot.(lab)(dlw_nf_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nf_bycbot.(lab)(dtot_nf_bycbot.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nf_bycbot.(lab)(dsw_nf_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nf_bycbot.(lab)(dlw_nf_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nf_bycbot.(lab)(dtot_nf_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nf_bycbot.(lab)(dsw_nf_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nf_bycbot.(lab)(dlw_nf_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nf_bycbot.(lab)(dtot_nf_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==5
             dsw_nf_byLTS.(lab)  = dsw_nofog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dlw_nf_byLTS.(lab)  = dlw_nofog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dtot_nf_byLTS.(lab) = dtot_nofog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_nf_byLTS.(lab)(dsw_nf_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_nf_byLTS.(lab)(dlw_nf_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_nf_byLTS.(lab)(dtot_nf_byLTS.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_nf_byLTS.(lab)(dsw_nf_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_nf_byLTS.(lab)(dlw_nf_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_nf_byLTS.(lab)(dtot_nf_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_nf_byLTS.(lab)(dsw_nf_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_nf_byLTS.(lab)(dlw_nf_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_nf_byLTS.(lab)(dtot_nf_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        end
    end
    
%% plot boxplot per param
    
    figure(333+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
       legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-150 100]);
    set(gca,'XTickLabel',{' '});title('\Delta noFog - MERRABASE');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-150 100]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
       legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-150 100]);  
    
end            

%% delta_walbnoFog - merrabasewAlb by intwc/intThick/top/bot
% sw/lw/tot
sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};

for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             dsw_fog_bywc.(lab)  = dsw_fog(strcmp(intwc_group,pre.(parlab)(i)),:);
             dlw_fog_bywc.(lab)  = dlw_fog(strcmp(intwc_group,pre.(parlab)(i)),:);
             dtot_fog_bywc.(lab) = dtot_fog(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_fog_bywc.(lab)(dsw_fog_bywc.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_fog_bywc.(lab)(dlw_fog_bywc.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_fog_bywc.(lab)(dtot_fog_bywc.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_fog_bywc.(lab)(dsw_fog_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_fog_bywc.(lab)(dlw_fog_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_fog_bywc.(lab)(dtot_fog_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_fog_bywc.(lab)(dsw_fog_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_fog_bywc.(lab)(dlw_fog_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_fog_bywc.(lab)(dtot_fog_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             dsw_fog_bythick.(lab)  = dsw_fog(strcmp(intthick_group,pre.(parlab)(i)),:);
             dlw_fog_bythick.(lab)  = dlw_fog(strcmp(intthick_group,pre.(parlab)(i)),:);
             dtot_fog_bythick.(lab) = dtot_fog(strcmp(intthick_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)) sw_plot_ = {0*zeros(1,3)};end
                 lw_plot_={dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)) lw_plot_= {0*zeros(1,3)};end
                 tot_plot_={dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)) tot_plot_ = {0*zeros(1,3)};end
             elseif  i==2
                 sw_plot={dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)) sw_plot = {0*zeros(1,3)};end
                 sw_plot=[sw_plot;sw_plot_];
                 lw_plot={dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)) lw_plot = {0*zeros(1,3)};end
                 lw_plot=[lw_plot;lw_plot_];
                 tot_plot={dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)) tot_plot = {0*zeros(1,3)};end
                 tot_plot=[tot_plot;tot_plot_];
             else
                 sw_plot__={dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dsw_fog_bythick.(lab)(dsw_fog_bythick.(lab)(:,1)~=0,:)) sw_plot__ = {0*zeros(1,3)};end
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dlw_fog_bythick.(lab)(dlw_fog_bythick.(lab)(:,1)~=0,:)) lw_plot__ = {0*zeros(1,3)};end
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)};
                 if isempty(dtot_fog_bythick.(lab)(dtot_fog_bythick.(lab)(:,1)~=0,:)) tot_plot__ = {0*zeros(1,3)};end
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             dsw_fog_byctop.(lab)  = dsw_fog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dlw_fog_byctop.(lab)  = dlw_fog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             dtot_fog_byctop.(lab) = dtot_fog(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={dsw_fog_byctop.(lab)(dsw_fog_byctop.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_fog_byctop.(lab)(dlw_fog_byctop.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_fog_byctop.(lab)(dtot_fog_byctop.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_fog_byctop.(lab)(dsw_fog_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_fog_byctop.(lab)(dlw_fog_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_fog_byctop.(lab)(dtot_fog_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_fog_byctop.(lab)(dsw_fog_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_fog_byctop.(lab)(dlw_fog_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_fog_byctop.(lab)(dtot_fog_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             dsw_fog_bycbot.(lab)  = dsw_fog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dlw_fog_bycbot.(lab)  = dlw_fog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             dtot_fog_bycbot.(lab) = dtot_fog(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_fog_bycbot.(lab)(dsw_fog_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_fog_bycbot.(lab)(dlw_fog_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_fog_bycbot.(lab)(dtot_fog_bycbot.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_fog_bycbot.(lab)(dsw_fog_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_fog_bycbot.(lab)(dlw_fog_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_fog_bycbot.(lab)(dtot_fog_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_fog_bycbot.(lab)(dsw_fog_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_fog_bycbot.(lab)(dlw_fog_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_fog_bycbot.(lab)(dtot_fog_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
             
        elseif kkk==5
             dsw_fog_byLTS.(lab)  = dsw_fog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dlw_fog_byLTS.(lab)  = dlw_fog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             dtot_fog_byLTS.(lab) = dtot_fog(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={dsw_fog_byLTS.(lab)(dsw_fog_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot_={dlw_fog_byLTS.(lab)(dlw_fog_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot_={dtot_fog_byLTS.(lab)(dtot_fog_byLTS.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={dsw_fog_byLTS.(lab)(dsw_fog_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={dlw_fog_byLTS.(lab)(dlw_fog_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={dtot_fog_byLTS.(lab)(dtot_fog_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={dsw_fog_byLTS.(lab)(dsw_fog_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={dlw_fog_byLTS.(lab)(dlw_fog_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={dtot_fog_byLTS.(lab)(dtot_fog_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        
        end
    end
    
%% plot boxplot per param
    
    figure(444+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
       legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-70 70]);
    set(gca,'XTickLabel',{' '});title('\Delta noFog_wAlb - MERRABASEwAlb');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
       legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-70 70]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-70 70]);  
    
end    

%% plot swmerra_wAlb instead of ds_fog

% sw/lw/tot
sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};

for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             sw_wAlb_bywc.(lab)  = swmerra_wAlb(strcmp(intwc_group,pre.(parlab)(i)),:);
             lw_wAlb_bywc.(lab)  = lwmerra_wAlb(strcmp(intwc_group,pre.(parlab)(i)),:);
             tot_wAlb_bywc.(lab) = totmerra_wAlb(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_bywc.(lab)(sw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_bywc.(lab)(lw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_bywc.(lab)(tot_wAlb_bywc.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_bywc.(lab)(sw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_wAlb_bywc.(lab)(lw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_wAlb_bywc.(lab)(tot_wAlb_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_wAlb_bywc.(lab)(sw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_bywc.(lab)(lw_wAlb_bywc.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_bywc.(lab)(tot_wAlb_bywc.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             sw_wAlb_bythick.(lab)  = swmerra_wAlb(strcmp(intthick_group,pre.(parlab)(i)),:);
             lw_wAlb_bythick.(lab)  = lwmerra_wAlb(strcmp(intthick_group,pre.(parlab)(i)),:);
             tot_wAlb_bythick.(lab) = totmerra_wAlb(strcmp(intthick_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_bythick.(lab)(sw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_bythick.(lab)(lw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_bythick.(lab)(tot_wAlb_bythick.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_bythick.(lab)(sw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_wAlb_bythick.(lab)(lw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_wAlb_bythick.(lab)(tot_wAlb_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_wAlb_bythick.(lab)(sw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_bythick.(lab)(lw_wAlb_bythick.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_bythick.(lab)(tot_wAlb_bythick.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             sw_wAlb_byctop.(lab)  = swmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             lw_wAlb_byctop.(lab)  = lwmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             tot_wAlb_byctop.(lab) = totmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_byctop.(lab)(sw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_byctop.(lab)(lw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_byctop.(lab)(tot_wAlb_byctop.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_byctop.(lab)(sw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot_];
                 lw_plot={lw_wAlb_byctop.(lab)(lw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot_];
                 tot_plot={tot_wAlb_byctop.(lab)(tot_wAlb_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot_];
             else
                 sw_plot__={sw_wAlb_byctop.(lab)(sw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_byctop.(lab)(lw_wAlb_byctop.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_byctop.(lab)(tot_wAlb_byctop.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             sw_wAlb_bycbot.(lab)  = swmerra_wAlb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             lw_wAlb_bycbot.(lab)  = lwmerra_wAlb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             tot_wAlb_bycbot.(lab) = totmerra_wAlb(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_bycbot.(lab)(sw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_bycbot.(lab)(lw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_bycbot.(lab)(tot_wAlb_bycbot.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_bycbot.(lab)(sw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_wAlb_bycbot.(lab)(lw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_wAlb_bycbot.(lab)(tot_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_wAlb_bycbot.(lab)(sw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_bycbot.(lab)(lw_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_bycbot.(lab)(tot_wAlb_bycbot.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
             
        elseif kkk==5
             sw_wAlb_byLTS.(lab)  = swmerra_wAlb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             lw_wAlb_byLTS.(lab)  = lwmerra_wAlb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             tot_wAlb_byLTS.(lab) = totmerra_wAlb(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_byLTS.(lab)(sw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_byLTS.(lab)(lw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_byLTS.(lab)(tot_wAlb_byLTS.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_byLTS.(lab)(sw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_wAlb_byLTS.(lab)(lw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_wAlb_byLTS.(lab)(tot_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_wAlb_byLTS.(lab)(sw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_byLTS.(lab)(lw_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_byLTS.(lab)(tot_wAlb_byLTS.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
             
         elseif kkk==6
             sw_wAlb_byLWP.(lab)  = swmerra_wAlb(strcmp(intlwp_group,pre.(parlab)(i)),:);
             lw_wAlb_byLWP.(lab)  = lwmerra_wAlb(strcmp(intlwp_group,pre.(parlab)(i)),:);
             tot_wAlb_byLWP.(lab) = totmerra_wAlb(strcmp(intlwp_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={sw_wAlb_byLWP.(lab)(sw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 lw_plot_={lw_wAlb_byLWP.(lab)(lw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 tot_plot_={tot_wAlb_byLWP.(lab)(tot_wAlb_byLWP.(lab)(:,1)~=0,:)};
             elseif  i==2
                 sw_plot={sw_wAlb_byLWP.(lab)(sw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_wAlb_byLWP.(lab)(lw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_wAlb_byLWP.(lab)(tot_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_wAlb_byLWP.(lab)(sw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_wAlb_byLWP.(lab)(lw_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_wAlb_byLWP.(lab)(tot_wAlb_byLWP.(lab)(:,1)~=0,:)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        
        end
    end
    
%% plot boxplot per param
    
    figure(555+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{3}]','WC 0.1-0.4 [g/m^{3}]]','WC 0.4-1.0 [g/m^{3}]','WC 1.0-4.0 [g/m^{3}]','WC > 4.0 [g/m^{3}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    elseif kkk==6
        legend('LWP 0-40 [g/m^{2}]','LWP 40-100 [g/m^{2}]','LWP 100-200 [g/m^{2}]','LWP 200-400 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    end
    ylim([-300 150]);
    set(gca,'XTickLabel',{' '});title('MERRABASEwAlb');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
     elseif kkk==6
        legend('LWP 0-40 [g/m^{2}]','LWP 40-100 [g/m^{2}]','LWP 100-200 [g/m^{2}]','LWP 200-400 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    end
    ylim([0 100]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
     elseif kkk==6
        legend('LWP 0-40 [g/m^{2}]','LWP 40-100 [g/m^{2}]','LWP 100-200 [g/m^{2}]','LWP 200-400 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    end
    ylim([-300 150]);  
    hold on; line(0:5,zeros(6,1),'Color','k','linewidth',1,'LineStyle',':')
end    
%% plot surface CRE vs. LWP/cldbot, with sea ice conc/sza
figure(123)
sza4prof = sza(1:4:end);
scatter3(cldbot,swmerra_wAlb(:,1),iceconc,24,iceconc,'filled');
axis([0 1000 -200 50 -5 100]);


figure(1231)
sza4prof = sza(1:4:end);
scatter3(swmerra_wAlb(:,1),anaLTS,iceconc,24,iceconc,'filled');
axis([-200 50 5 30 -5 100]);
xlabel('SW CRE [W/m^{2}]');ylabel('LTS (MERRA)');zlabel('Ice conc (%)');

figure(1232)
sza4prof = sza(1:4:end);
scatter3(lwmerra_wAlb(:,1),anaLTS,iceconc,24,iceconc,'filled');
axis([40 80 5 30 -5 100]);
xlabel('LW CRE [W/m^{2}]');ylabel('LTS (MERRA)');zlabel('Ice conc (%)');


figure(1234)
scatter3(cldbot,lwmerra_wAlb(:,1),sza4prof,24,sza4prof,'filled');
axis([-5 800 40 85 60 80]);
xlabel('CBH [m]');ylabel('LW CRE [W/m^{2}]');zlabel('SZA');
figure(12341)
scatter3(cldbot,swmerra_wAlb(:,1),sza4prof,24,sza4prof,'filled');
axis([-5 800 -200 50 60 80]);
xlabel('CBH [m]');ylabel('SW CRE [W/m^{2}]');zlabel('SZA');
%% plot merra

% sw/lw/tot
sw_plot={};sw_plot_={};sw_plot__={};
lw_plot={};lw_plot_={};lw_plot__={};
tot_plot={};tot_plot_={};tot_plot__={};

for kkk=1:length(param)
    parlab = strcat('par',num2str(kkk));
    for i = 1:length(pre.(parlab))
        pt = pre.(parlab){i,1};
        tlab_ = strrep(pt, '-', '_');tlab = strrep(tlab_, '.', '_');
        lab = strcat('par',tlab);
        if kkk==1
             sw_bywc.(lab)  = swmerra(strcmp(intwc_group,pre.(parlab)(i)),:);
             lw_bywc.(lab)  = lwmerra(strcmp(intwc_group,pre.(parlab)(i)),:);
             tot_bywc.(lab) = totmerra(strcmp(intwc_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={sw_bywc.(lab)};
                 lw_plot_={lw_bywc.(lab)};
                 tot_plot_={tot_bywc.(lab)};
             elseif  i==2
                 sw_plot={sw_bywc.(lab)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_bywc.(lab)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_bywc.(lab)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_bywc.(lab)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_bywc.(lab)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_bywc.(lab)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==2
             sw_bythick.(lab)  = swmerra(strcmp(intthick_group,pre.(parlab)(i)),:);
             lw_bythick.(lab)  = lwmerra(strcmp(intthick_group,pre.(parlab)(i)),:);
             tot_bythick.(lab) = totmerra(strcmp(intthick_group,pre.(parlab)(i)),:);
             % store for plotting
             if i==1
                 sw_plot_={sw_bythick.(lab)};
                 lw_plot_={lw_bythick.(lab)};
                 tot_plot_={tot_bythick.(lab)};
             elseif  i==2
                 sw_plot={sw_bythick.(lab)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_bythick.(lab)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_bythick.(lab)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_bythick.(lab)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_bythick.(lab)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_bythick.(lab)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==3
             sw_byctop.(lab)  = swmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             lw_byctop.(lab)  = lwmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             tot_byctop.(lab) = totmerra_wAlb(strcmp(cldtop_group,pre.(parlab)(i)),:);
             %
              % store for plotting
             if i==1
                 sw_plot_={sw_byctop.(lab)};
                 lw_plot_={lw_byctop.(lab)};
                 tot_plot_={tot_byctop.(lab)};
             elseif  i==2
                 sw_plot={sw_byctop.(lab)};
                 sw_plot=[sw_plot;sw_plot_];
                 lw_plot={lw_byctop.(lab)};
                 lw_plot=[lw_plot;lw_plot_];
                 tot_plot={tot_byctop.(lab)};
                 tot_plot=[tot_plot;tot_plot_];
             else
                 sw_plot__={sw_byctop.(lab)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_byctop.(lab)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_byctop.(lab)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        elseif kkk==4
             sw_bycbot.(lab)  = swmerra(strcmp(cldbot_group,pre.(parlab)(i)),:);
             lw_bycbot.(lab)  = lwmerra(strcmp(cldbot_group,pre.(parlab)(i)),:);
             tot_bycbot.(lab) = totmerra(strcmp(cldbot_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={sw_bycbot.(lab)};
                 lw_plot_={lw_bycbot.(lab)};
                 tot_plot_={tot_bycbot.(lab)};
             elseif  i==2
                 sw_plot={sw_bycbot.(lab)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_bycbot.(lab)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_bycbot.(lab)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_bycbot.(lab)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_bycbot.(lab)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_bycbot.(lab)};
                 tot_plot=[tot_plot;tot_plot__];
             end
             
        elseif kkk==5
             sw_byLTS.(lab)  = swmerra(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             lw_byLTS.(lab)  = lwmerra(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             tot_byLTS.(lab) = totmerra(strcmp(anaLTS_group,pre.(parlab)(i)),:);
             %
             % store for plotting
             if i==1
                 sw_plot_={sw_byLTS.(lab)};
                 lw_plot_={lw_byLTS.(lab)};
                 tot_plot_={tot_byLTS.(lab)};
             elseif  i==2
                 sw_plot={sw_byLTS.(lab)};
                 sw_plot=[sw_plot_;sw_plot];
                 lw_plot={lw_byLTS.(lab)};
                 lw_plot=[lw_plot_;lw_plot];
                 tot_plot={tot_byLTS.(lab)};
                 tot_plot=[tot_plot_;tot_plot];
             else
                 sw_plot__={sw_byLTS.(lab)};
                 sw_plot=[sw_plot;sw_plot__];
                 lw_plot__={lw_byLTS.(lab)};
                 lw_plot=[lw_plot;lw_plot__];
                 tot_plot__={tot_byLTS.(lab)};
                 tot_plot=[tot_plot;tot_plot__];
             end
        
        end
    end
    
%% plot boxplot per param
    
    figure(666+kkk);
    ax2(1)=subplot(311);
    aboxplot(sw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
    ylabel('SW CRE [W/m^{2}]'); % Set the X-axis label
    if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-300 150]);
    set(gca,'XTickLabel',{' '});title('MERRABASE');
    % LW delta by ice for wAlb case
    ax2(2)=subplot(312);
    aboxplot(lw_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
    ylabel('LW CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
         legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-50 150]);
    set(gca,'XTickLabel',{' '})
    % Total delta by ice for wAlb case
    ax2(3)=subplot(313);
    aboxplot(tot_plot,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
    ylabel('CRE [W/m^{2}]'); % Set the X-axis label
     if kkk==1
        legend('WC 0-0.1 [g/m^{2}]','WC 0.1-0.4 [g/m^{2}]]','WC 0.4-1.0 [g/m^{2}]','WC 1.0-4.0 [g/m^{2}]','WC > 4.0 [g/m^{2}]','Orientation','horizontal'); % Add a legend
    elseif kkk==2
        legend('CTCK 0-100 [m]','CTCK 100-300 [m]','CTCK 300-600 [m]','CTCK 600-1000 [m]','CTCK > 1000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==3
         legend('CTH 0-200 [m]','CTH 200-600 [m]','CTH 600-1000 [m]','CTH 1000-3000 [m]','CTH > 3000 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==4
        legend('CBH 0-100 [m]','CBH 100-200 [m]','CBH 200-500 [m]','CBH > 500 [m]','Orientation','horizontal'); % Add a legend
    elseif kkk==5
        legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
    end
    ylim([-300 150]);  
    
end    
              
% concatanate all for figure 555 surface
dsw_sur = [dswMERRAnfsur dswMERRAwasur dswMERRAnfwasur];
dlw_sur = [dlwMERRAnfsur dlwMERRAwasur dlwMERRAnfwasur];
dtot_sur= [dtotMERRAnfsur dtotMERRAwasur dtotMERRAnfwasur];

% air
dsw_air = [dswMERRAnfair dswMERRAwaair dswMERRAnfwaair];
dlw_air = [dlwMERRAnfair dlwMERRAwaair dlwMERRAnfwaair];
dtot_air= [dtotMERRAnfair dtotMERRAwaair dtotMERRAnfwaair];

% toa
dsw_toa = [dswMERRAnftoa dswMERRAwatoa dswMERRAnfwatoa];
dlw_toa = [dlwMERRAnftoa dlwMERRAwatoa dlwMERRAnfwatoa];
dtot_toa= [dtotMERRAnftoa dtotMERRAwatoa dtotMERRAnfwatoa];

% Concatenate each group in a  3 x n x 3 matrix
hh1 = cat(1, reshape(dsw_sur,[1 size(dsw_sur)]), reshape(dlw_sur,[1 size(dlw_sur)]), reshape(dtot_sur,[1 size(dtot_sur)]));
hh2 = cat(1, reshape(dsw_air,[1 size(dsw_air)]), reshape(dlw_air,[1 size(dlw_air)]), reshape(dtot_air,[1 size(dtot_air)]));
hh3 = cat(1, reshape(dsw_toa,[1 size(dsw_toa)]), reshape(dlw_toa,[1 size(dlw_toa)]), reshape(dtot_toa,[1 size(dtot_toa)]));

%% plot difference from base case

figure(555);
% surface CRE delta
ax1(1)=subplot(311);
aboxplot(hh1,'Colormap',[0.1 0.3 0.9;0.8 0.3 0.2;0.2 0.8 0.2]); % Advanced box plot
ylabel('\Delta CRE Surface [W/m^{2}]'); % Set the X-axis label
legend('SW','LW','Total','orientation','horizontal'); % Add a legend
ylim([-10 10]);
set(gca,'XTickLabel',{' '})
% aircraft CRE delta
ax1(2)=subplot(312);
aboxplot(hh2,'Colormap',[0.1 0.3 0.9;0.8 0.3 0.2;0.2 0.8 0.2]); % Advanced box plot
ylabel('\Delta CRE C-130 [W/m^{2}]'); % Set the X-axis label
%legend('SW','LW','Total'); % Add a legend
ylim([-20 20]);
set(gca,'XTickLabel',{' '})
% TOA CRE delta
ax1(3)=subplot(313);
aboxplot(hh3,'labels',groupcase(:,2:end),'Colormap',[0.1 0.3 0.9;0.8 0.3 0.2;0.2 0.8 0.2]); % Advanced box plot
ylabel('\Delta CRE TOA [W/m^{2}]'); % Set the X-axis label
%legend('SW','LW','Total'); % Add a legend
ylim([-100 100]);
%% save figure 555
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'deltaCREbyLevel_MERRA20150325')];
        save_fig(555,fi,true);

%% figure(666)
% differences by ice concentration for 2 cases
% 1 - wAlb-Merrabase
figure(666);
% SW CRE delta by ice for wAlb case
ax2(1)=subplot(311);
aboxplot(sw_walb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);
set(gca,'XTickLabel',{' '});title('\Delta wAlb - MERRABASE');
% LW delta by ice for wAlb case
ax2(2)=subplot(312);
aboxplot(lw_walb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for wAlb case
ax2(3)=subplot(313);
aboxplot(tot_walb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);

% tight subplots
% p(i,:) = get(ax(i), 'position');
%% save figure 666
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'deltaCREbyice_wAlb_nozerocases_20150325')];
        save_fig(666,fi,true);
%% noFogwAlb-Merrabase
figure(777);
% SW CRE delta by ice for wAlb case
ax2(1)=subplot(311);
aboxplot(sw_nfwalb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);title('\Delta noFogwAlb - MERRABASE');
set(gca,'XTickLabel',{' '})
% LW delta by ice for wAlb case
ax2(2)=subplot(312);
aboxplot(lw_nfwalb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for wAlb case
ax2(3)=subplot(313);
aboxplot(tot_nfwalb_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 100]);
set(gcf,'color','white');

% tight subplots
% p(i,:) = get(ax(i), 'position');
%% save figure 777
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'deltaCREbyice_noFogwAlb_20150320')];
        save_fig(777,fi,true);
        
        
%% noFogwAlb-MerrabasewAlb to extract fog layer influence above weighted ice
% this includes all fog layers, even when the fog is the only cloud layer
% (i.e. change from cloud to clear)
figure(888);
% SW CRE delta by ice for wAlbnoFog versus wAlb case
ax2(1)=subplot(311);
aboxplot(sw_fog_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-200 100]);title('\Delta MERRABASEwAlb - noFogwAlb');
set(gca,'XTickLabel',{' '})
% LW delta by ice for fog case
ax2(2)=subplot(312);
aboxplot(lw_fog_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-50 150]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for fog case
ax2(3)=subplot(313);
aboxplot(tot_fog_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-200 100]);
set(gcf,'color','white');
linkaxes(ax2,'x');
% tight subplots
% p(i,:) = get(ax(i), 'position');
%% save figure 888
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'deltaCREbyice_fog_20150406')];
        save_fig(888,fi,true);

%% noFogwAlb-MerrabasewAlb to extract fog layer influence above weighted ice
% this includes only fog layers that are part of a greater cloud system (i.e. more than 1)

figure(8881);
% SW CRE delta by ice for wAlbnoFog versus wAlb case
ax2(1)=subplot(311);
aboxplot(sw_fogonly_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 50]);title('\Delta MERRABASEwAlb - noFogwAlb');
set(gca,'XTickLabel',{' '})
% LW delta by ice for fog case
ax2(2)=subplot(312);
aboxplot(lw_fogonly_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-50 100]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for fog case
ax2(3)=subplot(313);
aboxplot(tot_fogonly_byice,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-100 50]);
set(gcf,'color','white');
linkaxes(ax2,'x');
% tight subplots
% p(i,:) = get(ax(i), 'position');

%% save figure 8881
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'deltaCREbyice_fog_mcld_20150512')];
        save_fig(8881,fi,true);    
        
%% MERAA_wAlb for ncld>1
figure(8882);
% SW CRE delta by ice for wAlbnoFog versus wAlb case
ax2(1)=subplot(311);
aboxplot(sw_byice_wAlb_mcld,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-250 100]);title('MERRABASEwAlb');
set(gca,'XTickLabel',{' '})
% LW delta by ice for fog case
ax2(2)=subplot(312);
aboxplot(lw_byice_wAlb_mcld,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-50 150]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for fog case
ax2(3)=subplot(313);
aboxplot(tot_byice_wAlb_mcld,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 30-50 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-200 100]);
set(gcf,'color','white');
linkaxes(ax2,'x');
% tight subplots
% p(i,:) = get(ax(i), 'position');

%% save figure 8881
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'CREbyice_wAlb_fog_mcld_20150512')];
        save_fig(8882,fi,true);        
%% MerrabasewAlb by ice
figure(999);
% SW CRE  by ice for wAlb case
ax2(1)=subplot(311);
aboxplot(sw_byice_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('SW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-200 100]);title('MERRABASEwAlb');
set(gca,'XTickLabel',{' '})
% LW CRE  by ice for wAlb case
ax2(2)=subplot(312);
aboxplot(lw_byice_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('LW CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-20 100]);
set(gca,'XTickLabel',{' '})
% Total CRE  by ice for wAlb case
ax2(3)=subplot(313);
aboxplot(tot_byice_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('Total CRE [W/m^{2}]'); % Set the X-axis label
legend('ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]','Orientation','horizontal'); % Add a legend
ylim([-200 100]);
set(gcf,'color','white');
 hold on; line(0:5,zeros(6,1),'Color','k','linewidth',1,'LineStyle',':')
% tight subplots
% p(i,:) = get(ax(i), 'position');
%% save figure 999
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'CREbyice_wAlb_20150406')];
        save_fig(999,fi,true);
        
%% MERRAwAlb by LTS
figure(1000);
% SW CRE delta by ice for wAlbnoFog versus wAlb case
ax2(1)=subplot(311);
aboxplot(sw_byLTS_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','blue_down'); % Advanced box plot
ylabel('\Delta SW CRE [W/m^{2}]'); % Set the X-axis label
legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
ylim([-300 200]);title('MERRABASEwAlb');
set(gca,'XTickLabel',{' '})
% LW delta by ice for fog case
ax2(2)=subplot(312);
aboxplot(lw_byLTS_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','red_down'); % Advanced box plot
ylabel('\Delta LW CRE [W/m^{2}]'); % Set the X-axis label
legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
ylim([-50 150]);
set(gca,'XTickLabel',{' '})
% Total delta by ice for fog case
ax2(3)=subplot(313);
aboxplot(tot_byLTS_wAlb,'labels',{'Surface', 'C-130','TOA'},'colorgrad','green_down'); % Advanced box plot
ylabel('\Delta Total CRE [W/m^{2}]'); % Set the X-axis label
legend('LTS 0-12 [k]','LTS 12-20 [K]','LTS 20-30 [K]','Orientation','horizontal'); % Add a legend
ylim([-300 200]);
set(gcf,'color','white');

% tight subplots
% p(i,:) = get(ax(i), 'position');
%% save figure 1000
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'CREbyLTS_wAlb_20150406')];
        save_fig(1000,fi,true);        


%% plot cloud properties versus ice concentration
% create cloud property sets by ice concentration

for i = 1:length(unique(allcld_iceconcgroup))
    ice = unique(allcld_iceconcgroup); cice = ice{i,1};
    icelab = strrep(cice, '-', '_');
    lab = strcat('ice',icelab);
    reff.(lab)      = allcld_cldref(strcmp(allcld_iceconcgroup,ice(i)),:);
    wc.(lab)        = allcld_cldwc(strcmp(allcld_iceconcgroup,ice(i)),:);
    wp.(lab)        = allcld_cldwp(strcmp(allcld_iceconcgroup,ice(i)),:);
    cldthick.(lab)  = allcld_cldthick(strcmp(allcld_iceconcgroup,ice(i)),:);
    hcldtop.(lab)   = allcld_cldtop(strcmp(allcld_iceconcgroup,ice(i)),:);
    hcldbot.(lab)   = allcld_cldbot(strcmp(allcld_iceconcgroup,ice(i)),:);
    hncld.(lab)     = allcld_ncld(strcmp(allcld_iceconcgroup,ice(i)),:);
end
hreff = {reff.ice0_15;reff.ice15_30;reff.ice30_50;reff.ice50_70;reff.ice70_85;reff.ice85_100}; % Create a cell array with the data for each group
hwc   = {wc.ice0_15;wc.ice15_30;wc.ice30_50;wc.ice50_70;wc.ice70_85;wc.ice85_100};
hwp   = {wp.ice0_15;wp.ice15_30;wp.ice30_50;wp.ice50_70;wp.ice70_85;wp.ice85_100};
hcldthick = {cldthick.ice0_15;cldthick.ice15_30;cldthick.ice30_50;cldthick.ice50_70;cldthick.ice70_85;cldthick.ice85_100};
hhcldtop = {hcldtop.ice0_15;hcldtop.ice15_30;hcldtop.ice30_50;hcldtop.ice50_70;hcldtop.ice70_85;hcldtop.ice85_100};
hhcldbot = {hcldbot.ice0_15;hcldbot.ice15_30;hcldbot.ice30_50;hcldbot.ice50_70;hcldbot.ice70_85;hcldbot.ice85_100};
hhncld   = {hncld.ice0_15;hncld.ice15_30;hncld.ice30_50;hncld.ice50_70;hncld.ice70_85;hncld.ice85_100};

%% figures od cloud properties
% plot ncld with ice concentration
figure(88)
aboxplot(hhncld,'colorgrad','blue_down','colorrev','true'); % Advanced box plot
ylabel('# of clouds in examined profiles','fontsize',12); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'linewidth',2)
xtixloc = [0.70:0.12:1.35];      %  label locations
%xtixloc = [0.72 0.82 0.94 1.06 1.16 1.28];
set(gca,'XTickMode','auto','XTickLabel',{'ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]'},'XTick',xtixloc);
xlim([0.6 1.4]);set(gca,'TickLength',[0,0]);
ylim([0 7.5]);xlim([0.6 1.4]);
%% save figure 88
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'cldnum_byiceconc_20150323')];
        save_fig(88,fi,true);
% plot cloud numbers by proportion
maxcld = 7;
figure(887)
for i=1:length(hhncld)
    for ii=1:maxcld
        msize = sum(hhncld{i,:}==ii);
        if msize>0
             plot(i,ii,'o','color','k','markerfacecolor',[0.5 0.5 0.5],'markersize',...
             sum(hhncld{i,:}==ii));hold on;
        end
    end
end
xtixloc = [1:6];      %  label locations
%xtixloc = [0.72 0.82 0.94 1.06 1.16 1.28];
set(gca,'XTickMode','auto','XTickLabel',{'ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]'},'XTick',xtixloc);
ylabel('# of clouds in examined profiles','fontsize',12); 
xlim([0 7]);
ylim([0 7.5]);
%% save figure 87
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'cldnum_byiceconc_scalebymarkersizenormalizedtototalnumber20150323')];
        save_fig(887,fi,true);
%% plot cloud properties by ice concentration
figure(888)
ax3(1)=subplot(511);
aboxplot(hreff,'colorgrad','red_up'); % Advanced box plot
ylabel('R_{eff} [\mum]'); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'XTick',[])
ylim([5 40]);xlim([0.6 1.4]);
ax3(2)=subplot(512);
aboxplot(hwp,'colorgrad','blue_up');
               %'colorgrad','red_down','colorrev','true'); % Advanced box plot
ylabel('CWP [g/m^{2}]]'); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'XTick',[])
ylim([0 100]);%ylim([0 2.5]);this is for wc
xlim([0.6 1.4]);
ax3(3)=subplot(513);
aboxplot(hcldthick,'colorgrad','green_up');
               %'colorgrad','red_down','colorrev','true'); % Advanced box plot
ylabel('CTCK [m]'); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'XTick',[])
ylim([0 1000]);xlim([0.6 1.4]);
ax3(4)=subplot(514);
aboxplot(hhcldtop,'colorgrad','orange_up');
               %'colorgrad','red_down','colorrev','true'); % Advanced box plot
ylabel('CTH [m]'); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'XTick',[])
ylim([0 2500]);xlim([0.6 1.4]);
ax3(5)=subplot(515);
aboxplot(hhcldbot,'colorgrad','orange_down');
               %'colorgrad','red_down','colorrev','true'); % Advanced box plot
ylabel('CBH [m]'); % Set the X-axis label
set(gca,'XTickLabel',{' '});set(gca,'XTick',[])
ylim([0 1500]);xlim([0.6 1.4]);
xl = get(gca,'xlim');
xtixloc = [0.70:0.12:1.35];      %  label locations
%xtixloc = [0.72 0.82 0.94 1.06 1.16 1.28];
set(gca,'XTickMode','auto','XTickLabel',{'ice 0-15 [%]','ice 15-30 [%]','ice 30-50 [%]','ice 50-70 [%]','ice 70-85 [%]','ice 85-100 [%]'},'XTick',xtixloc);
xlim([0.6 1.4]);set(gca,'TickLength',[0,0])
%% save figure 888
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', 'cld_properties_by_iceconc_20150512')];
        save_fig(888,fi,true);


%% create profile cloud properties per ice conc
figure(2222);
ax(1)=subplot(511);
plot(iceconc,anaLTS,   'p','color',[0.2 0.8 0.2],'markerfacecolor',[0.2 0.8 0.2]);hold on;
plot(iceconc,anaLLSS,'p','color',[0.2 0.8 0.2]);hold on;
ax(2)=subplot(512);
plot(iceconc,intthick, 'p','color',[0.5 0.2 0.2],'markerfacecolor',[0.5 0.2 0.2]);hold on;
ax(3)=subplot(513);
plot(iceconc,int_lwp,  'p','color',[0.2 0.2 0.8],'markerfacecolor',[0.2 0.2 0.8]);hold on;
ax(4)=subplot(514);
plot(iceconc,cldbot,   's','color',[0.2 0.8 0.8],'markerfacecolor',[0.2 0.8 0.8]);hold on;
ax(5)=subplot(515);
plot(iceconc,cldbot,   's','color',[0.8 0.8 0.2],'markerfacecolor',[0.8 0.8 0.2]);hold on;


%plot(iceconc,airLTS,'o','color',[0.2 0.2 0.8],'markerfacecolor',[0.2 0.8 0.2]);hold on;
%plot(iceconc,anaLTS,'o','color',[0.2 0.2 0.8]);hold on;

%% create all cloud properties per ice conc
figure(3333);
ax(1)=subplot(311);
plot(allcld_iceconc,allcld_cldref,   'p','color',[0.8 0.5 0.2],'markerfacecolor',[0.8 0.5 0.2]);hold on;
axis([0 100 0 40]);
ax(2)=subplot(312);
plot(allcld_iceconc,allcld_cldthick, 'p','color',[0.5 0.2 0.2],'markerfacecolor',[0.5 0.2 0.2]);hold on;
axis([0 100 0 1500]);
ax(3)=subplot(313);
plot(allcld_iceconc(allcld_cldwp~=40),allcld_cldwp(allcld_cldwp~=40),  'p','color',[0.2 0.2 0.8],'markerfacecolor',[0.2 0.2 0.8]);hold on;

% 30-50% ice
indall30_50 = find(allcld_iceconc>=30 & allcld_iceconc<=50);
dateall_30_50 = allcld_date(indall30_50);
profall_30_50 = allcld_profile(indall30_50);
crefall_30_50 = allcld_cldref(indall30_50);

ind30_50 = find(iceconc>=30 & iceconc<=50);
date_30_50 = date(ind30_50);
prof_30_50 = profile(ind30_50);
lwp_30_50  = int_lwp(ind30_50);
analts_30_50  = anaLTS(ind30_50);
airlts_30_50  = airLTS(ind30_50);
anallss_30_50  = anaLLSS(ind30_50);
airllss30_50  = airLLSS(ind30_50);

% 50-70% ice
indall50_70 = find(allcld_iceconc>=50 & allcld_iceconc<=70);
dateall_50_70 = allcld_date(indall50_70);
profall_50_70 = allcld_profile(indall50_70);
crefall_50_70 = allcld_cldref(indall50_70);

ind50_70 = find(iceconc>=50 & iceconc<=70);
date_50_70 = date(ind50_70);
prof_50_70 = profile(ind50_70);
lwp_50_70  = int_lwp(ind50_70);
analts_50_70  = anaLTS(ind50_70);
airlts_50_70  = airLTS(ind50_70);
anallss_50_70  = anaLLSS(ind50_70);
airllss50_70  = airLLSS(ind50_70);

%% plot versus various params
figure(444);
axx(1)=subplot(311);
plot(iceconc,swMERRAsur      ,'o','color',[0.1 0.3 0.9],'markerfacecolor',[0.1 0.3 0.9]);hold on;
plot(iceconc,lwMERRAsur      ,'o','color',[0.8 0.3 0.2],'markerfacecolor',[0.8 0.3 0.2]);hold on;
plot(iceconc,totMERRAsur     ,'o','color',[0.2 0.8 0.2],'markerfacecolor',[0.2 0.8 0.2]);hold on;
axx(2)=subplot(312);
plot(iceconc,swMERRAnfsur      ,'o','color',[0.1 0.3 0.9],'markerfacecolor',[0.1 0.3 0.9]);hold on;
plot(iceconc,lwMERRAnfsur      ,'o','color',[0.8 0.3 0.2],'markerfacecolor',[0.8 0.3 0.2]);hold on;
plot(iceconc,totMERRAnfsur     ,'o','color',[0.2 0.8 0.2],'markerfacecolor',[0.2 0.8 0.2]);hold on;
axx(3)=subplot(313);
plot(iceconc,swMERRAwasur      ,'o','color',[0.1 0.3 0.9],'markerfacecolor',[0.1 0.3 0.9]);hold on;
plot(iceconc,lwMERRAwasur      ,'o','color',[0.8 0.3 0.2],'markerfacecolor',[0.8 0.3 0.2]);hold on;
plot(iceconc,totMERRAwasur     ,'o','color',[0.2 0.8 0.2],'markerfacecolor',[0.2 0.8 0.2]);hold on;
linkaxes(axx,'xy');


figure(555);
plot(swMERRAwasur,lwMERRAwasur      ,'o','color',[0.8 0.8 0.8],'markerfacecolor',[0.8 0.8 0.8]);hold on;
plot(swMERRAwasur(ncld<=1),lwMERRAwasur(ncld<=1),'.c');hold on;
plot(swMERRAwasur(ncld>1),lwMERRAwasur(ncld>1),'.b');hold on;
plot(swMERRAwasur(cldtop<=500),lwMERRAwasur(cldtop<=500),'.g');hold on;
plot(swMERRAwasur(cldtop>500) ,lwMERRAwasur(cldtop>500) ,'.r');hold on;
plot(swMERRAwasur(anaLTS<18),lwMERRAwasur(anaLTS<18),'.g');hold on;
plot(swMERRAwasur(anaLTS>=18) ,lwMERRAwasur(anaLTS>=18) ,'.r');hold on;
plot(swMERRAwasur(anaLLSS<6),lwMERRAwasur(anaLLSS<6),'.g');hold on;
plot(swMERRAwasur(anaLLSS>=6) ,lwMERRAwasur(anaLLSS>=6) ,'.r');hold on;
plot(swMERRAwasur(cldtop<=1000),lwMERRAwasur(cldtop<=1000),'.g');hold on;
plot(swMERRAwasur(cldtop>1000) ,lwMERRAwasur(cldtop>1000) ,'.r');hold on;
xlabel('SW');ylabel('LW');

plot(iceconc(ncld<=1),swMERRAsur(ncld<=1),'.c');hold on;
plot(iceconc(ncld> 1),swMERRAsur(ncld>1) ,'.b');hold on;
plot(iceconc(cldtop<=500),swMERRAsur(cldtop<=500),'.g');hold on;
plot(iceconc(cldtop>500) ,swMERRAsur(cldtop>500) ,'.r');hold on;
plot(iceconc(anaLTS<20),swMERRAsur(anaLTS<20),'.g');hold on;
plot(iceconc(anaLTS>=20) ,swMERRAsur(anaLTS>=20) ,'.r');hold on;
plot(iceconc(airLLSS<10),swMERRAsur(airLLSS<10),'.g');hold on;
plot(iceconc(airLLSS>=10) ,swMERRAsur(airLLSS>=10) ,'.r');hold on;

plot(iceconc,swCREsur(MERRAnfind==1)    ,'o','color',[0.5 0.5 0.9],'markerfacecolor',[0.5 0.5 0.9]);hold on;
plot(iceconc,MERRAwa    ,'o','color',[0.1 0.1 0.9],'markerfacecolor',[0.1 0.1 0.9]);hold on;
%% save profile data in .dat file for clustering
X = [iceconc ncld cldtop cldbot anaLTS anaLLSS];
X = X(~isnan(cldtop),:);
%Z = linkage(X,method,metric);
Z = linkage(X,'average','euclidean');
figure;[H,T] = dendrogram(Z,'colorthreshold',1200);
set(H,'LineWidth',2);title('average euclidean');
c = clusterdata(X,'linkage','average','distance','euclidean',...
    'maxclust',4);
% plot SW/LW with 4 clusters
figure(555);
SW = swMERRAwasur(~isnan(cldtop));
LW = lwMERRAwasur(~isnan(cldtop));
plot(SW,LW      ,'o','color',[0.8 0.8 0.8],'markerfacecolor',[0.8 0.8 0.8]);hold on;
plot(SW(c==1),LW(c==1),'.','color',[0.2 0.2 0.9]);hold on;
plot(SW(c==2),LW(c==2),'.','color',[0.2 0.9 0.3]);hold on;
plot(SW(c==3),LW(c==3),'.','color',[0 0 0]);hold on;
plot(SW(c==4),LW(c==4),'.','color',[0.9 0.2 0.3]);hold on;
xlabel('SW');ylabel('LW');

Z1 = linkage(X,'complete','euclidean');
figure;[H1,T1] = dendrogram(Z1,'colorthreshold',1200);
set(H1,'LineWidth',2);title('complete euclidean');
c1 = clusterdata(X,'linkage','complete','distance','euclidean',...
    'maxclust',5);
% plot SW/LW with 5 clusters
figure(556);
SW = swMERRAwasur(~isnan(cldtop));
LW = lwMERRAwasur(~isnan(cldtop));
plot(SW,LW      ,'o','color',[0.8 0.8 0.8],'markerfacecolor',[0.8 0.8 0.8]);hold on;
plot(SW(c1==1),LW(c1==1),'.','color',[0.8 0.2 0.8]);hold on;
plot(SW(c1==2),LW(c1==2),'.','color','c');hold on;
plot(SW(c1==3),LW(c1==3),'.','color',[0 0 0]);hold on;
plot(SW(c1==4),LW(c1==4),'.','color',[0.2 0.9 0.3]);hold on;
plot(SW(c1==5),LW(c1==5),'.','color',[0.9 0.2 0.3]);hold on;
xlabel('SW');ylabel('LW');

%% save new s struct

casestr = [];
for i=1:length(cases)
    casestr = strcat(casestr,'_',cases{:,i},'_');
end
fiout = strcat('F:\ARISE\ArcticCRF\libRadTran_output\arise\','cases',casestr,'test.mat');
save(fiout,'-struct','d');


return;