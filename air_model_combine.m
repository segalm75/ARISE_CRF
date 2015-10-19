%% Details of the function:
%  NAME:
% air_model_combine
%----------------------------------
% PURPOSE:
%  - combine vertical profiles for ARISE
%    from aircraft and model reanalysis
%    and generate profiles for RT simulations
%
% CALLING SEQUENCE:
%  air_model_combine
%
% INPUT:
%  - air.mat structure of C-130 SeaIce profile summary
% 
% 
% OUTPUT:
%  - native and interpolated profiles from aircraft (.dat files)
%  - native model relevant profiles (.dat file)
%  - cloud profile .dat fiel for RT simulations
%  - atmospheric state profiles for RT simulations
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%  - genCloudProf.m
%  - gen_atmos_rt.m
% 
%
% NEEDED FILES/INPUT:
% 
% - e.g.:F:\ARISE\C-130data\Met\ARISEairprocessed.mat
% - e.g. F:\ARISE\ArcticCRF\METdata\GEOS_FP\ariseFPdata\...
%           GEOS.fp.asm.inst3_3d_asm_Np.20140904_0000.V01.nc4
%
% EXAMPLE:
%  - air_model_combine
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Sep-25-2015
% modified from star_ana_compare to use only relevant profiles
% when reading model files
% -------------------------------------------------------------------------
%% function routine
function air_model_combine
clear all; close all;
startup_plotting; 

% prepare for plotting
ColorSet_ = varycolor(90);% max num of profiles per flight
ColorSet  = ColorSet_(1:6:90,:);
%% load air data
%  .mat file of selected profiles

arisedir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
filename = 'ARISEairprocessed_with_insitu_withWVparams20150921';

reanadir = 'F:\ARISE\ArcticCRF\METdata\GEOS_FP\ariseFPdata\';

%  load air data
air     =  load([arisedir filename '.mat']);

% find number of days analyzed
nFields    = sum(~structfun(@isempty,air));
starfieldNames  = fieldnames(air);
k=0;

        % initialize vars for saving air data for plotting

        yearm       = [];
        monthm      = [];
        daym        = [];
        profilem    = [];
        iceconcm    = [];
        iceconc_stdm= [];
        lonm        = [];
        latm        = [];
        zm          = [];
        Tm          = [];
        Thetam      = [];
        RHm         = [];
        Qm          = [];
        windSm      = [];
        windDm      = [];
        staticPm    = [];
        SatVPwaterm = [];
        SatVPicem   = [];
        cldflagm    = [];
        cldwcm      = [];
        cldicm      = [];

        % initialize vars for saving model data for plotting
        date = [];
        nprof= [];
        profile = [];
        plev = [];
        galt_mean = []; galt_std = [];
        zalt_mean = []; zalt_std = [];
        o3mmr_mean= []; o3mmr_std= [];
        omega_mean= []; omega_std= [];
        qice_mean = []; qice_std = [];
        qliq_mean = []; qliq_std = [];
        qwv_mean  = []; qwv_std  = [];
        rh_mean   = []; rh_std   = [];
        t_mean    = []; t_std    = [];
        theta_mean= []; theta_std= [];
        winds_mean= []; winds_std= [];
        windd_mean= []; windd_std= [];
        cflag_ana = [];
        ctype_ana = [];
        cth_ana   = [];
        cbh_ana   = [];
        ctk_ana   = [];
        cwc_ana   = [];
        lwc_ana   = [];
        iwc_ana   = [];
        iceconc_mean = [];
        iceconc_std  = [];
        lon_mean = []; lon_std = [];
        lat_mean = []; lat_std = [];
        airTinterp = []; airMMRinterp = []; airRHinterp = []; airThetainterp = [];
        s700_air = []; s850_air = []; s925_air = [];
        s700_ana = []; s850_ana = []; s925_ana = [];
        cflag_air  = []; cwc_air = []; cic_air = [];

% loop over all days
for i=1:nFields
    
    % Extract only the "profnum" fields
    names  = fieldnames(air.(starfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);%max(star.(starfieldNames{i,:}).prof.nprof);
    
    % load ice-conc data
    daystr=starfieldNames{i,:};
    %ice = readAMSR2data({daystr(4:11)});% this is with gradient calculated
    ice = readAMSR2dataV2({daystr(4:11)});%
    
    % get reanalysis files for each day
    day        = daystr(4:11);  
    file       = strcat('GEOS.fp.asm.inst3_3d_asm_Np.',day,'*','nc4');
    anaDirInfo = dir([reanadir file]);
    
    % extract forecast time fields
        for kk=1:length(anaDirInfo)
             dm{kk}       = strsplit_ms(anaDirInfo(kk).name,'.');
             daytimes(kk) = dm{kk}(4);
             stime        = (strsplit_ms(daytimes(kk),'_')); ttmp = cellstr(stime{:});
             % convert to decimal time
             dectimes(kk) = str2double(ttmp{:}(1:2)) + ...
                            str2double(ttmp{:}(3:4))/60;
        end
    
    legendall ={};
    
    if nProfiles>0
        
        dtice       = [];
        dtocean     = [];
        dmmrice     = [];
        dmmrocean   = [];
        drhice      = [];
        drhocean    = [];
        tsurf_ice   = [];
        tsurf_ocean = [];
        dtsurf_ice  = [];
        dtsurf_ocean= [];
        
        s700ice_a   = [];
        s850ice_a   = [];
        s925ice_a   = [];
        s700ocean_a = [];
        s850ocean_a = [];
        s925ocean_a = [];
        
        s700ice_m   = [];
        s850ice_m   = [];
        s925ice_m   = [];
        s700ocean_m = [];
        s850ocean_m = [];
        s925ocean_m = [];
        
        for j=1:nProfiles
            
            profstr    = strcat('profnum',num2str(j));
            %tmp       = starfieldNames{j,:};
            
            % get profile time
            
            proftime = nanmean(air.(daystr).(profstr).time);
            
            % choose closest file in time
            [imin,ind] = min(abs(dectimes-proftime));
            cfile      = strcat('GEOS.fp.asm.inst3_3d_asm_Np.',daytimes{ind},'*','nc4');
            fileInfo   = dir([reanadir cfile]);
            fileinfoNp = ncinfo([reanadir,fileInfo.name]);
            
            % read essential vars
            lon = ncread([reanadir,fileInfo.name], 'lon');
            lat = ncread([reanadir,fileInfo.name], 'lat');
            
            % find closest reanalysis grids
            lonintrp= interp1(lon,[1:length(lon)],...
                              min(air.(starfieldNames{i,:}).(profstr).lon):max(air.(starfieldNames{i,:}).(profstr).lon),...
                              'nearest'); if length(lonintrp)>1 lonintrp=lonintrp(1):lonintrp(2); else lonintrp = lonintrp(1); end
            latintrp= interp1(lat,[1:length(lat)],...
                              min(air.(starfieldNames{i,:}).(profstr).lat):max(air.(starfieldNames{i,:}).(profstr).lat),...
                              'nearest'); if length(latintrp)>1 latintrp=latintrp(1):latintrp(2); else latintrp = latintrp(1); end
            
            nmodelProf = length(lonintrp)*length(latintrp);
                          
            % store vars
           
             vars = {'lev','H','O3','OMEGA','QI','QL','QV','RH','T','U','V'};
             % units: [hPa], [m], [kg kg-1], [pa/s],[kg kg-1],[kg kg-1],[kg
             % kg-1],[fraction],[K],[m/s],[m/s]
             
                 for l = 1:length(vars)

                     % read vars
                     vartmp = ncread([reanadir,fileInfo.name], vars{l});

                     % read in relevant profiles, make stats
                     if l>1 % all vars except pressure levels
                         
                         
                            if (length(lonintrp)==1 && length(latintrp)==1)
                                     tmp_avg = squeeze((vartmp([lonintrp],[latintrp],:)));  
                                     tmp_std = squeeze((vartmp([lonintrp], [latintrp],:)));
                                     % replace missiing values with NaN
                                     tmp_avg ( tmp_avg >1E14) = NaN;
                                     tmp_std ( tmp_avg >1E14) = NaN;

                                     % save into air struct
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_mean')) = tmp_avg;
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_std'))  = tmp_std;
                                     
                            elseif (length(lonintrp)==1 && length(latintrp)> 1) || ...
                               (length(lonintrp)> 1 && length(latintrp)==1)
                                     tmp_avg = squeeze(nanmean(vartmp([lonintrp],[latintrp],:)));  
                                     tmp_std = squeeze(nanstd(vartmp([lonintrp], [latintrp],:)));
                                     % replace missiing values with NaN
                                     tmp_avg ( tmp_avg >1E14) = NaN;
                                     tmp_std ( tmp_avg >1E14) = NaN;

                                     % save into air struct
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_mean')) = tmp_avg;
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_std'))  = tmp_std;
                            else
                                     tmp_avg = squeeze(nanmean(nanmean(vartmp([lonintrp],[latintrp],:))));  
                                     tmp_std = squeeze(nanstd(nanstd(vartmp([lonintrp], [latintrp],:))));
                                     % replace missiing values with NaN
                                     tmp_avg ( tmp_avg >1E14) = NaN;
                                     tmp_std ( tmp_avg >1E14) = NaN;

                                     % save into air struct
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_mean')) = tmp_avg;
                                     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{l},'_std'))  = tmp_std;
                            end

                     else
                             % save pressure levels
                             air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})) = vartmp;
                     end

                 end
            
           % generate calculated fields
           
           air.(starfieldNames{i,:}).(profstr).ana.zalt_mean = ...
                                     geo2z(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{2},'_mean')));
                                 
           air.(starfieldNames{i,:}).(profstr).ana.theta_mean = ...
                                     temp2theta(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean')), ...
                                                           air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})));
           air.(starfieldNames{i,:}).(profstr).ana.h2o        = mmr2molec(air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})),...
                                                                          air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{7},'_mean')),...
                                                                          'h2o',air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean')));
           air.(starfieldNames{i,:}).(profstr).ana.o3         = mmr2molec(air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})),...
                                                                          air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{3},'_mean')),...
                                                                          'o3',air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean')));
           air.(starfieldNames{i,:}).(profstr).ana.air        = air2molec(air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})),...
                                                                          'air',air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean')));
           air.(starfieldNames{i,:}).(profstr).ana.o2         = air2molec(air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1})),...
                                                                          'o2',air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean')));
                                                                      
           air.(starfieldNames{i,:}).(profstr).ana.theta700   = air.(starfieldNames{i,:}).(profstr).ana.theta_mean(13) - ...
                                                                air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1);%13 is 700 mbar
           air.(starfieldNames{i,:}).(profstr).ana.theta850   = air.(starfieldNames{i,:}).(profstr).ana.theta_mean(7) - ...
                                                                air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1);%7  is 850 mbar
           air.(starfieldNames{i,:}).(profstr).ana.theta925   = air.(starfieldNames{i,:}).(profstr).ana.theta_mean(4) - ...
                                                                air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1);%4  is 925 mbar
                                                            
                                                            
           %% find closest iceconc in aircraft grid
           % create interpolated grid
           F = TriScatteredInterp(ice.(strcat('amsr',daystr(4:11))).longrid(:),...
                                  ice.(strcat('amsr',daystr(4:11))).latgrid(:),...
                                  ice.(strcat('amsr',daystr(4:11))).iceconc(:),'nearest');
           
           % evaluate at aircraft locations
           
           air.(starfieldNames{i,:}).(profstr).iceconc_mean = ...
                  nanmean(F(air.(starfieldNames{i,:}).(profstr).lon,...
                    air.(starfieldNames{i,:}).(profstr).lat));
           air.(starfieldNames{i,:}).(profstr).iceconc_std = ...
                  nanstd(F(air.(starfieldNames{i,:}).(profstr).lon,...
                    air.(starfieldNames{i,:}).(profstr).lat));
           %---------------------------------------------------
           
                                                            
           %% generate model cloud flags, type and properties
           %-------------------------------------------------
           
           ml = 42; %model levels = 42
           % cflag=0; no cloud, cflag=1; cloud
           % ctype=0, liquid, ctype=1, ice, ctype=2, mixed
           
           qliq = air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{5},'_mean'));
           qice = air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{6},'_mean'));
           p    = air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1}));
           t    = air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean'));
           
           air.(starfieldNames{i,:}).(profstr).ana.lwc   = qliq.*air2dens(p,t)*1000;
           air.(starfieldNames{i,:}).(profstr).ana.iwc   = qice.*air2dens(p,t)*1000;
           
           lwc = air.(starfieldNames{i,:}).(profstr).ana.lwc;
           iwc = air.(starfieldNames{i,:}).(profstr).ana.iwc;
           air.(starfieldNames{i,:}).(profstr).ana.cwc  = lwc + iwc;
           cwc = lwc + iwc;
           
           for c=1:ml
               
               if       lwc(c)>0.0001 && ...
                        iwc(c)<0.0001

                                air.(starfieldNames{i,:}).(profstr).ana.cflag(c) = 1;
                                air.(starfieldNames{i,:}).(profstr).ana.ctype(c) = 0;
                        
               elseif  lwc(c)<0.0001 && ...
                       iwc(c)>0.0001
                   
                               air.(starfieldNames{i,:}).(profstr).ana.cflag(c) = 1;
                               air.(starfieldNames{i,:}).(profstr).ana.ctype(c) = 1;

               elseif   lwc(c)>0.0001 || ...
                        iwc(c)>0.0001
                        
                               air.(starfieldNames{i,:}).(profstr).ana.cflag(c) = 1;
                               air.(starfieldNames{i,:}).(profstr).ana.ctype(c) = 2;
                               
               elseif   lwc(c)<0.0001 && ...
                        iwc(c)<0.0001
                        
                               air.(starfieldNames{i,:}).(profstr).ana.cflag(c) = 0;
                               air.(starfieldNames{i,:}).(profstr).ana.ctype(c) = NaN;
               end
           end
           
           air.(starfieldNames{i,:}).(profstr).ana.cflag = air.(starfieldNames{i,:}).(profstr).ana.cflag';
           air.(starfieldNames{i,:}).(profstr).ana.ctype = air.(starfieldNames{i,:}).(profstr).ana.ctype';
           
           % store only up to 7km (plev=19)
             p = 19;
             
%              if strcmp(starfieldNames{i,:},'met20140907_') && strcmp(profstr,'profnum2')
%                  p = 17;
%              end
             
           if sum(air.(starfieldNames{i,:}).(profstr).ana.cflag)>0
               if  min(find(air.(starfieldNames{i,:}).(profstr).ana.cflag(1:p)==1))==1
                   % CTH and CBH are taken at center of pressure level
                   % plus/minus 100m (which is about half layer)
                   air.(starfieldNames{i,:}).(profstr).ana.cbh   = 0;
                   air.(starfieldNames{i,:}).(profstr).ana.cth   = ...
                   air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(max(find(air.(starfieldNames{i,:}).(profstr).ana.cflag(1:p)==1)))+100;
                   air.(starfieldNames{i,:}).(profstr).ana.ctk   = air.(starfieldNames{i,:}).(profstr).ana.cth - ...
                                                                   air.(starfieldNames{i,:}).(profstr).ana.cbh;
               else
                   air.(starfieldNames{i,:}).(profstr).ana.cbh   = ...
                   air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(min(find(air.(starfieldNames{i,:}).(profstr).ana.cflag(1:p)==1)))-100;
                   air.(starfieldNames{i,:}).(profstr).ana.cth   = ...
                   air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(max(find(air.(starfieldNames{i,:}).(profstr).ana.cflag(1:p)==1)))+100;
                   
                   air.(starfieldNames{i,:}).(profstr).ana.ctk   = air.(starfieldNames{i,:}).(profstr).ana.cth - ...
                                                                   air.(starfieldNames{i,:}).(profstr).ana.cbh;
               end
           end
           
           %-----------------------------------------------------------------------------------------------------------------
                                                                      
           % store data for further plotting
           % store only up to 7km (plev=19)
           % p = 19;
           
            date = [date; repmat(daystr(4:11),p,1)];
            profile   = [profile; repmat(j,p,1)];
            nprof= [nprof; repmat(nmodelProf,p,1)];
            plev = [plev; air.(starfieldNames{i,:}).(profstr).ana.(strcat('p',vars{1}))(1:p)];
            galt_mean = [galt_mean;  air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{2},'_mean'))(1:p)]; 
            galt_std  = [galt_std;   air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{2},'_std'))(1:p)];
            zalt_mean = [zalt_mean;  air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)]; 
            o3mmr_mean= [o3mmr_mean; air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{3},'_mean'))(1:p)]; 
            o3mmr_std=  [o3mmr_std;  air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{3},'_std'))(1:p)];
            omega_mean= [omega_mean; air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{4},'_mean'))(1:p)]; 
            omega_std=  [omega_std;  air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{4},'_std'))(1:p)];
            qice_mean = [qice_mean;  air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{5},'_mean'))(1:p)]; 
            qice_std =  [qice_std;   air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{5},'_std'))(1:p)];
            qliq_mean = [qliq_mean;  air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{6},'_mean'))(1:p)]; 
            qliq_std =  [qliq_std;   air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{6},'_std'))(1:p)];
            qwv_mean  = [qwv_mean;   air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{7},'_mean'))(1:p)]; 
            qwv_std  =  [qwv_std;    air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{7},'_std'))(1:p)];
            rh_mean   = [rh_mean;    air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{8},'_mean'))(1:p)]; 
            rh_std   =  [rh_std;     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{8},'_std'))(1:p)];
            t_mean    = [t_mean;     air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean'))(1:p)]; 
            t_std    =  [t_std;      air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_std'))(1:p)];
            theta_mean = [theta_mean; air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1:p)];
            cflag_ana  = [cflag_ana;  air.(starfieldNames{i,:}).(profstr).ana.cflag(1:p)]; 
            ctype_ana  = [ctype_ana;  air.(starfieldNames{i,:}).(profstr).ana.ctype(1:p)]; 
            cth_ana    = [cth_ana;    repmat(air.(starfieldNames{i,:}).(profstr).ana.cth,p,1)]; 
            cbh_ana    = [cbh_ana;    repmat(air.(starfieldNames{i,:}).(profstr).ana.cbh,p,1)]; 
            ctk_ana    = [ctk_ana;    repmat(air.(starfieldNames{i,:}).(profstr).ana.ctk,p,1)]; 
            cwc_ana    = [cwc_ana;    air.(starfieldNames{i,:}).(profstr).ana.cwc(1:p)]; 
            lwc_ana    = [lwc_ana;    air.(starfieldNames{i,:}).(profstr).ana.lwc(1:p)]; 
            iwc_ana    = [iwc_ana;    air.(starfieldNames{i,:}).(profstr).ana.iwc(1:p)]; 
            iceconc_mean = [iceconc_mean; repmat(air.(starfieldNames{i,:}).(profstr).iceconc_mean,p,1)];
            iceconc_std  = [iceconc_std ; repmat(air.(starfieldNames{i,:}).(profstr).iceconc_std ,p,1)];
            lon_mean     = [lon_mean;     repmat(air.(starfieldNames{i,:}).(profstr).meanlon ,p,1)];
            lon_std      = [lon_std;      repmat(nanstd(air.(starfieldNames{i,:}).(profstr).lon) ,p,1)];
            lat_mean     = [lat_mean;     repmat(air.(starfieldNames{i,:}).(profstr).meanlat ,p,1)];
            lat_std      = [lat_std;      repmat(nanstd(air.(starfieldNames{i,:}).(profstr).lat) ,p,1)];
            
            
            % add wind field direction and speed
            % ws = wind speed; wd = wind direction
            % u positive wind towards the east (from west)
            % v positive is wind towards the north (from south)
            u  = air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{10},'_mean'))(1:p);
            v  = air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{11},'_mean'))(1:p);
            ws = sqrt(u.^2 + ...
                      v.^2);
            wd = uv2wd(u,v);
            winds_mean= [winds_mean; ws];
            windd_mean= [windd_mean; wd];
            
            
            
           %% generate atmospheric profiles for RT simulations
           %  this is for model/reanalysis (using std atmosphere and model)
           %  mode=0 is for MERRA, 0.1 is for GEOS-FP
            [air.(starfieldNames{i,:}).(profstr).ana.atmfile] = gen_atmos_prof4rt(air.(starfieldNames{i,:}).(profstr)...
                                                                             ,air.(starfieldNames{i,:}).(profstr).ana,0.1,daystr(4:11),profstr);
           %  this is for aircraft (using std atm, model for above aircrat, and aircraft)
           %  this uses interpolated aircrft fields
%            [air.(starfieldNames{i,:}).(profstr).airAtmNative]= gen_atmos_prof4rt(air.(starfieldNames{i,:}).(profstr)...
%                                                                             ,air.(starfieldNames{i,:}).(profstr).ana,1,daystr(4:11),profstr);
          
           %% compare temp
           figure(10);
           plot(flipud(air.(starfieldNames{i,:}).(profstr).Static_AirT+273),...
                flipud(air.(starfieldNames{i,:}).(profstr).z/1000),':','color',ColorSet(j,:),'linewidth',2);hold on;
           plot(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean'))(1:p),...
                air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                '-','color',ColorSet(j,:),'linewidth',2);hold on;
           xlabel('Temperature [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([230 280 0 7]);
           ha=text(231,6.5,daystr(4:11)) ;set(ha,'fontsize',12,'color','k');
           hb=text(270,6.3,'.. C130')    ;set(hb,'fontsize',12,'color','k');
           hc=text(270,5.8,'- GEOS-5-FP')    ;set(hc,'fontsize',12,'color','k');
           %legend('C130','MERRA');
           set(gca,'Fontsize',12);
           clear ha hb hc
           %% generate cloud fields for each aircraft profile
           doy = datestr2doy(daystr(4:11),'yyyymmdd');
           dprof = j;
           k = k+1;
           allprof = k;
           [air.(starfieldNames{i,:}).(profstr).cld] = genAirCloudProf(air.(starfieldNames{i,:}).(profstr),dprof,allprof,daystr(4:11),num2str(doy));        
           %---------------------------------------------------------------------------------------------------------------------------------------%
           
           %% compare temp with cloud flags
           figure(101);
           plot(flipud(air.(starfieldNames{i,:}).(profstr).Static_AirT+273),...
                flipud(air.(starfieldNames{i,:}).(profstr).z/1000),...
                ':','color',ColorSet(j,:),'linewidth',2);hold on;
           plot(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean'))(1:p),...
                air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                '-','color',ColorSet(j,:),'linewidth',2);hold on;
           % cloud from aircraft 
           if (air.(starfieldNames{i,:}).(profstr).cld.cldnum>0)
                   plot(flipud(air.(starfieldNames{i,:}).(profstr).Static_AirT(air.(starfieldNames{i,:}).(profstr).cld.cldflag==1)+273),...
                        flipud(air.(starfieldNames{i,:}).(profstr).z(air.(starfieldNames{i,:}).(profstr).cld.cldflag==1)/1000),...
                        'o','color',ColorSet(j,:),'markerFaceColor',ColorSet(j,:),'markersize',6);hold on;
                    if air.(starfieldNames{i,:}).(profstr).cld.cldbelow==1
                           %ib = find(min(air.(starfieldNames{i,:}).(profstr).z));% == min(air.(starfieldNames{i,:}).(profstr).cld.cldbot));
                           [mb,ib] = min(air.(starfieldNames{i,:}).(profstr).z);
                           plot((air.(starfieldNames{i,:}).(profstr).Static_AirT(ib)+273),...
                                (min(air.(starfieldNames{i,:}).(profstr).z(ib))/1000),...
                                's','color',ColorSet(j,:),'markersize',12);hold on;
                    end
                    if air.(starfieldNames{i,:}).(profstr).cld.cldabove>0
                           [mt,it] = max(air.(starfieldNames{i,:}).(profstr).z);
                           if air.(starfieldNames{i,:}).(profstr).cld.cldabove==1
                               % thin cloud; probably cirrus
                               plot((air.(starfieldNames{i,:}).(profstr).Static_AirT(it)+273),...
                                    ((air.(starfieldNames{i,:}).(profstr).z(it))/1000),...
                                    'd','color',ColorSet(j,:),'markersize',12);hold on;
                           else
                               % liquid/thicker cloud (cldabove==2)
                               plot((air.(starfieldNames{i,:}).(profstr).Static_AirT(it)+273),...
                                   ((air.(starfieldNames{i,:}).(profstr).z(it))/1000),...
                                    'o','color',ColorSet(j,:),'markersize',12);hold on;
                           end
                    end
           end
           % cloud from model
           if sum(air.(starfieldNames{i,:}).(profstr).ana.cflag)>0
                   plot(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{9},'_mean'))(air.(starfieldNames{i,:}).(profstr).ana.cflag==1),...
                        air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(air.(starfieldNames{i,:}).(profstr).ana.cflag==1)/1000,...
                        'p','color',ColorSet(j,:),'markerFaceColor',ColorSet(j,:),'markersize',6);hold on;
           end
           
           xlabel('Temperature [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[0:0.5:4]);set(gca,'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'});
           axis([245 280 0 4]);
           ha=text(246,3.8,daystr(4:11))   ;set(ha,'fontsize',12,'color','k');
           hb=text(265,3.7,'.. C130')      ;set(hb,'fontsize',12,'color','k');
           hc=text(265,3.5,'- GEOS-5-FP')      ;set(hc,'fontsize',12,'color','k');
           hd=text(265,3.3,'o aircraft Cloud flag') ;set(hd,'fontsize',12,'color','k');
           he=text(265,3.1,'p GEOS-5-FP Cloud flag')    ;set(he,'fontsize',12,'color','k');
           hd=text(265,2.9,'s cloud below aircraft');set(hd,'fontsize',12,'color','k');
           he=text(265,2.7,'d cloud above aircraft');set(he,'fontsize',12,'color','k');
           %legend('C130','MERRA');
           set(gca,'Fontsize',12);
           %----------------------
           
           %% compare RH
           
           figure(11);
           plot(flipud(air.(starfieldNames{i,:}).(profstr).RH/100),...
                flipud(air.(starfieldNames{i,:}).(profstr).z/1000),':','color',ColorSet(j,:),'linewidth',2);hold on;
           plot(air.(starfieldNames{i,:}).(profstr).ana.(strcat(vars{8},'_mean'))(1:p),...
                air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                '-','color',ColorSet(j,:),'linewidth',2);hold on;
           xlabel('RH [fraction]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([0.1 1 0 7]);
           ha=text(0.31,6.5,daystr(4:11)) ;set(ha,'fontsize',12,'color','k');
           hb=text(0.85,6.3,'.. C130')    ;set(hb,'fontsize',12,'color','k');
           hc=text(0.85,5.8,'- GEOS-5-FP')    ;set(hc,'fontsize',12,'color','k');
           %legend('C130','MERRA');
           set(gca,'Fontsize',12);
           clear ha hb hc
           %----------------------
           
           %% smooth, then compare difference in Temp/RH over ice/ocean
           starx = smooth(air.(starfieldNames{i,:}).(profstr).z/1000,0.05);%lowess used
           stary = smooth(starx,air.(starfieldNames{i,:}).(profstr).Static_AirT+273,0.15);%
           starz = smooth(starx,air.(starfieldNames{i,:}).(profstr).MMR);%
           starw = smooth(starx,air.(starfieldNames{i,:}).(profstr).RH/100);%
           starv = smooth(starx,air.(starfieldNames{i,:}).(profstr).theta);%
          
           if (length(unique(stary)) < length(stary)) ||  (length(unique(starz)) < length(starz))% then bin
               [npointsx,starx] = bin_profiles_simple(air.(starfieldNames{i,:}).(profstr).z/1000,...
                                           air.(starfieldNames{i,:}).(profstr).z/1000);            
                                           starx=starx(~isnan(starx));
               [npointsy,stary] = bin_profiles_simple(air.(starfieldNames{i,:}).(profstr).z/1000,...
                                           air.(starfieldNames{i,:}).(profstr).Static_AirT+273);
                                           stary=stary(~isnan(stary));
               [npointsz,starz] = bin_profiles_simple(air.(starfieldNames{i,:}).(profstr).z/1000,...
                                           air.(starfieldNames{i,:}).(profstr).MMR);
                                           starz=starz(~isnan(starz));
               [npointsw,starw] = bin_profiles_simple(air.(starfieldNames{i,:}).(profstr).z/1000,...
                                           air.(starfieldNames{i,:}).(profstr).RH/100);
                                           starw=starw(~isnan(starw));
               [npointsv,starv] = bin_profiles_simple(air.(starfieldNames{i,:}).(profstr).z/1000,...
                                           air.(starfieldNames{i,:}).(profstr).theta);
                                           starv=starv(~isnan(starv));
           end
           
           if length(starx)>2
               air.(starfieldNames{i,:}).(profstr).t700ind = interp1(starx,[1:length(starx)],2.500,'nearest'); % this is for ~2500 m at 700hPa 
               air.(starfieldNames{i,:}).(profstr).t850ind = interp1(starx,[1:length(starx)],1.400,'nearest'); % this is for ~1400 m at 850hPa 
               air.(starfieldNames{i,:}).(profstr).t925ind = interp1(starx,[1:length(starx)],0.700 ,'nearest');% this is for ~700  m at 925hPa    
                % figure;plot(stary,starx,'-r');
                % t [K]
                air.(starfieldNames{i,:}).(profstr).airTinterp=interp1(starx,stary,...
                                   air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000);
                               
                air.(starfieldNames{i,:}).(profstr).airThetainterp=interp1(starx,starv,...
                                   air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000);
                               
                dt = air.(starfieldNames{i,:}).(profstr).ana.T_mean - air.(starfieldNames{i,:}).(profstr).airTinterp;
                
                % MMR [g/kg]
                if length(starz)>2
                    air.(starfieldNames{i,:}).(profstr).airMMRinterp=interp1(starx,starz,...
                                       air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000);
                    dmmr = air.(starfieldNames{i,:}).(profstr).ana.QV_mean*1000 - air.(starfieldNames{i,:}).(profstr).airMMRinterp;
                else
                    air.(starfieldNames{i,:}).(profstr).airMMRinterp = [];
                    dmmr = NaN(length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean),1);
                end
                
                % RH [fraction]
                if length(starz)>2
                    air.(starfieldNames{i,:}).(profstr).airRHinterp=interp1(starx,starw,...
                                       air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000);
                    drh = air.(starfieldNames{i,:}).(profstr).ana.RH_mean - air.(starfieldNames{i,:}).(profstr).airRHinterp;
                else
                    air.(starfieldNames{i,:}).(profstr).airRHinterp = [];
                    drh = NaN(length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean),1);
                end
                
           else
                % t [K]
                air.(starfieldNames{i,:}).(profstr).airTinterp = [];
                dt = NaN(length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean),1);
                % MMR [g/kg]
                air.(starfieldNames{i,:}).(profstr).airMMRinterp = [];
                dmmr = NaN(length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean),1);
                % RH [fraction]
                air.(starfieldNames{i,:}).(profstr).airRHinterp = [];
                drh = NaN(length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean),1);
           end
           
           % interp airT/airMMR/airRH to reana levels, including LML
           if sum(~isNaN(air.(starfieldNames{i,:}).(profstr).airTinterp))>2
                air.(starfieldNames{i,:}).(profstr).airTinterp(1) = ...
                       interp1(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(2:4),air.(starfieldNames{i,:}).(profstr).airTinterp(2:4),...
                       air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1),'linear','extrap');
                       air.(starfieldNames{i,:}).(profstr).airTinterp(p+1:end) = NaN;
                       
                       
                air.(starfieldNames{i,:}).(profstr).airThetainterp(1) = ...
                       interp1(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(2:4),air.(starfieldNames{i,:}).(profstr).airThetainterp(2:4),...
                       air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1),'linear','extrap');
                       air.(starfieldNames{i,:}).(profstr).airThetainterp(p+1:end) = NaN;
                       
                dt = air.(starfieldNames{i,:}).(profstr).ana.T_mean - air.(starfieldNames{i,:}).(profstr).airTinterp;
                
           else
               air.(starfieldNames{i,:}).(profstr).airTinterp     = (NaN(1,length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean)))';
               air.(starfieldNames{i,:}).(profstr).airThetainterp = (NaN(1,length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean)))';
           end
           
           airTinterp     = [airTinterp;     air.(starfieldNames{i,:}).(profstr).airTinterp(1:p)];
           airThetainterp = [airThetainterp; air.(starfieldNames{i,:}).(profstr).airThetainterp(1:p)];
           
           % create stability fields for saving
           % aircraft
            s700_air = [s700_air;repmat(air.(starfieldNames{i,:}).(profstr).airThetainterp(13) - ... 
                                        air.(starfieldNames{i,:}).(profstr).airThetainterp(1),p,1)];
            s850_air = [s850_air;repmat(air.(starfieldNames{i,:}).(profstr).airThetainterp(7) - ... 
                                        air.(starfieldNames{i,:}).(profstr).airThetainterp(1),p,1)];
            s925_air = [s925_air;repmat(air.(starfieldNames{i,:}).(profstr).airThetainterp(4) - ... 
                                        air.(starfieldNames{i,:}).(profstr).airThetainterp(1),p,1)];
                                    
           % model/analysis
            s700_ana = [s700_ana;repmat(air.(starfieldNames{i,:}).(profstr).ana.theta_mean(13) - ... 
                                        air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1),p,1)];
            s850_ana = [s850_ana;repmat(air.(starfieldNames{i,:}).(profstr).ana.theta_mean(7) - ... 
                                        air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1),p,1)];
            s925_ana = [s925_ana;repmat(air.(starfieldNames{i,:}).(profstr).ana.theta_mean(4) - ... 
                                        air.(starfieldNames{i,:}).(profstr).ana.theta_mean(1),p,1)];
                                         
           
           % MMR
           if sum(~isNaN(air.(starfieldNames{i,:}).(profstr).airMMRinterp))>2
               air.(starfieldNames{i,:}).(profstr).airMMRinterp(1) = ...
                   interp1(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(2:4),air.(starfieldNames{i,:}).(profstr).airMMRinterp(2:4),...
                           air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1),'linear','extrap');
               air.(starfieldNames{i,:}).(profstr).airMMRinterp(p+1:end) = NaN;
               
               dmmr = air.(starfieldNames{i,:}).(profstr).ana.QV_mean*1000 - air.(starfieldNames{i,:}).(profstr).airMMRinterp;
               
           else
               air.(starfieldNames{i,:}).(profstr).airMMRinterp = (NaN(1,length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean)))';
           end
           
           airMMRinterp = [airMMRinterp; air.(starfieldNames{i,:}).(profstr).airMMRinterp(1:p)];
           
           % RH
           if sum(~isNaN(air.(starfieldNames{i,:}).(profstr).airRHinterp))>2
               air.(starfieldNames{i,:}).(profstr).airRHinterp(1) = ...
                   interp1(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(2:4),air.(starfieldNames{i,:}).(profstr).airRHinterp(2:4),...
                           air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1),'linear','extrap');
               air.(starfieldNames{i,:}).(profstr).airRHinterp(19:end) = NaN;
               
               drh = air.(starfieldNames{i,:}).(profstr).ana.RH_mean - air.(starfieldNames{i,:}).(profstr).airRHinterp;
               
           else
               air.(starfieldNames{i,:}).(profstr).airRHinterp = (NaN(1,length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean)))';
           end
           
           airRHinterp = [airRHinterp; air.(starfieldNames{i,:}).(profstr).airRHinterp(1:p)];
           
           if ~isnan(air.(starfieldNames{i,:}).(profstr).iceconc_mean) && air.(starfieldNames{i,:}).(profstr).iceconc_mean > 15
               an = [0.1 0.9 0.9];% sea-ice is cyan
               dtice = [dtice;(dt(1:p))'];
               tsurf_ice = [tsurf_ice;
                            repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               dtsurf_ice= [dtsurf_ice; repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(1) - ... 
                                        air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               dmmrice = [dmmrice;dmmr(1:p)'];
               drhice  = [drhice; drh(1:p)'];
               
               s700ice_a = [s700ice_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(13) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               s850ice_a = [s850ice_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(7) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               s925ice_a = [s925ice_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(4) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
                                         
               s700ice_m = [s700ice_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(13) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
               s850ice_m = [s850ice_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(7) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
               s925ice_m = [s925ice_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(4) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
                                         
           else
               an = [0.1 0.1 0.9];% open ocean is blue
               dtocean = [dtocean;dt(1:p)'];
               tsurf_ocean = [tsurf_ocean;
                              repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               dtsurf_ocean= [dtsurf_ocean;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(1) - ... 
                                           air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               dmmrocean =   [dmmrocean;dmmr(1:p)'];
               drhocean =    [drhocean;drh(1:p)'];
               
               s700ocean_a = [s700ocean_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(13) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               s850ocean_a = [s850ocean_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(7) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
               s925ocean_a = [s925ocean_a;repmat(air.(starfieldNames{i,:}).(profstr).airTinterp(4) - ... 
                                             air.(starfieldNames{i,:}).(profstr).airTinterp(1),p,1)];
                                         
               s700ocean_m = [s700ocean_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(13) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
               s850ocean_m = [s850ocean_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(7) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
               s925ocean_m = [s925ocean_m;repmat(air.(starfieldNames{i,:}).(profstr).ana.T_mean(4) - ... 
                                             air.(starfieldNames{i,:}).(profstr).ana.T_mean(1),p,1)];
                                         
           end
           
           %% compile cloud profile into reanalysis levels
            [airCldprof,npointscld,~] = bin_profiles_ana(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean,...
                                                         air.(starfieldNames{i,:}).(profstr).z/1000,...
                                                         air.(starfieldNames{i,:}).(profstr).cld.cldflag);
            [ibot, ~] = min(find(airCldprof==1));
            itop      = interp1(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000,...
                                [1:length(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000)],...
                                 max(air.(starfieldNames{i,:}).(profstr).z/1000),'nearest');
            % add cloud below aircraft instances
            if air.(starfieldNames{i,:}).(profstr).cld.cldbelow==1  airCldprof (1:ibot)=1; end  
            % eliminate non corresponding altitudes
            airCldprof(itop+1:end) = NaN;
            
            cflag_air = [cflag_air; airCldprof(1:p)'];
            
            %% compile cloud wc into reanalysis levels
            [airWCprof,npointswc,airWCmean] = bin_profiles_ana(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean,...
                                              air.(starfieldNames{i,:}).(profstr).z/1000,...
                                              air.(starfieldNames{i,:}).(profstr).LWC1_gm3);
            airWCmean(airWCmean<0) = NaN;
            % eliminate non corresponding altitudes
            airWCmean(itop+1:end) = NaN;
            cwc_air = [cwc_air; airWCmean(1:p)'];
            
            %% compile cloud ic into reanalysis levels
            [airICprof,npointsic,airICmean] = bin_profiles_ana(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean,...
                                              air.(starfieldNames{i,:}).(profstr).z/1000,...
                                              (air.(starfieldNames{i,:}).(profstr).TWC_gm3 - ...
                                              abs(air.(starfieldNames{i,:}).(profstr).LWC1_gm3)));
                                                   
            airICmean(airICmean<0) = NaN; 
             % eliminate non corresponding altitudes
            airICmean(itop+1:end) = NaN;
            cic_air = [cic_air; airICmean(1:p)'];
            
            
           %% generate atmospheric profiles for RT simulations
           %  this is for aircraft (using std atm, model for above aircrat, and aircraft)
           %  this uses interpolated aircrft fields
           [air.(starfieldNames{i,:}).(profstr).airAtmInterp]     = gen_atmos_prof4rt(air.(starfieldNames{i,:}).(profstr)...
                                                                            ,air.(starfieldNames{i,:}).(profstr).ana,2,daystr(4:11),profstr);
          
           %% plot Temp differences
           figure(100);
           plot((dt),(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000),'-o','color',an,...
                'markerfacecolor',an,'linewidth',0.5,'markersize',4);hold on;
           xlabel('\DeltaT (GEOS-5-FP-C130) [K]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([-10 10 0 7]);
           %legend('blue is over open ocean','cyan is over sea ice');
           set(gca,'Fontsize',12);
           
           %% plot RH differences over ice/ocean
           figure(111);
           plot((drh),(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000),'-o','color',an,...
                'markerfacecolor',an,'linewidth',0.5,'markersize',4);hold on;
           xlabel('\DeltaRH (GEOS-5-FP-C130) [g/kg]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           axis([-0.5 0.5 0 7]);
           %legend('blue is over open ocean','cyan is over sea ice');
           set(gca,'Fontsize',12);
           
           %% plot MMR differences over ice/ocean
           figure(121);
           plot((dmmr),(air.(starfieldNames{i,:}).(profstr).ana.zalt_mean/1000),'-o','color',an,...
                'markerfacecolor',an,'linewidth',0.5,'markersize',4);hold on;
           xlabel('\DeltaMMR (GEOS-5-FP-C130) [g/kg]','FontSize',12);
           ylabel('Altitude [km]','fontSize',12);
           set(gca,'YTick',[1:10]);set(gca,'YTickLabel',{'1','2','3','4','5','6','7','8','9','10'});
           set(gca,'XTick',[-5:2.5:5]);set(gca,'YTickLabel',{'-5.0','-2.5','0.0','2.5','5.0'});
           axis([-5 5 0 7]);
           %legend('blue is over open ocean','cyan is over sea ice');
           set(gca,'Fontsize',12);
           
             %% store parameters from air struct profiles to .txt file for plotting

            year = str2double(daystr(4:7));         yearm      = [yearm;   repmat(year,             length(air.(starfieldNames{i,:}).(profstr).z),1)];
            month= str2double(daystr(8:9));         monthm     = [monthm;  repmat(month,            length(air.(starfieldNames{i,:}).(profstr).z),1)];
            day  = str2double(daystr(10:11));       daym       = [daym;    repmat(day,              length(air.(starfieldNames{i,:}).(profstr).z),1)];
            profilenumber = j;                      profilem   = [profilem;repmat(profilenumber,    length(air.(starfieldNames{i,:}).(profstr).z),1)];
            
            iceconcm      = [iceconcm; repmat(air.(starfieldNames{i,:}).(profstr).iceconc_mean,...
                                                  length(air.(starfieldNames{i,:}).(profstr).z),1)];
            iceconc_stdm  = [iceconc_stdm; repmat(air.(starfieldNames{i,:}).(profstr).iceconc_std,...
                                                  length(air.(starfieldNames{i,:}).(profstr).z),1)];
            Tm            = [Tm;         (air.(starfieldNames{i,:}).(profstr).Static_AirT+273)];
            lonm          = [lonm;        air.(starfieldNames{i,:}).(profstr).lon];
            latm          = [latm;        air.(starfieldNames{i,:}).(profstr).lat];
            zm            = [zm;          air.(starfieldNames{i,:}).(profstr).z/1000];
            staticPm      = [staticPm;    air.(starfieldNames{i,:}).(profstr).staticP];
            RHm           = [RHm;         air.(starfieldNames{i,:}).(profstr).RH];
            Qm            = [Qm;          air.(starfieldNames{i,:}).(profstr).MMR];%C130 is g/Kg
            SatVPwaterm   = [SatVPwaterm; air.(starfieldNames{i,:}).(profstr).SatVPwater];
            SatVPicem     = [SatVPicem;   air.(starfieldNames{i,:}).(profstr).SatVPice];
            windSm        = [windSm;      air.(starfieldNames{i,:}).(profstr).WindS];
            windDm        = [windDm;      air.(starfieldNames{i,:}).(profstr).WindD];
            cldflagm      = [cldflagm;    air.(starfieldNames{i,:}).(profstr).cld.cldflag];
            cldwcm        = [cldwcm;      air.(starfieldNames{i,:}).(profstr).LWC1_gm3];
            cldicm        = [cldicm;      air.(starfieldNames{i,:}).(profstr).TWC_gm3 - ...
                                          abs(air.(starfieldNames{i,:}).(profstr).LWC1_gm3)];

            % add theta
            theta=temp2theta((air.(starfieldNames{i,:}).(profstr).Static_AirT+273),...
                                    air.(starfieldNames{i,:}).(profstr).staticP);

            Thetam = [Thetam; theta];

           
           
        end % end of profiles
        
        %% save figure 10
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'tempCompare_newprofiles20150930')];
        save_fig(10,fi,false);
        close(10);
        
        %% save figure 101
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'tempCompare_w_cldflag_air_model20150930')];
        save_fig(101,fi,false);
        close(101);
        
        %% save figure 11
        fi=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'RHCompare_newprofiles20150930')];
        save_fig(11,fi,false);
        close(11);
        
        %% continue and save figure 100
        %------------------------------%
            figure(100)
            hold on;
            if size(dtice,1)>1
                boxplot(dtice,  'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.9 0.9]); hold on;
                set(gca,'YTickLabel',{' '})   % Erase ylabels  
            end
            if size(dtocean,1)>1
                boxplot(dtocean,'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.1 0.9]);
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
            fi1=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'deltatemp_ana_air')];
            save_fig(100,fi1,false);
            close(100);
            
        %% continue and save figure 111
        %-------------------------------%
            figure(111)
            hold on;
            if size(drhice,1)>1
                boxplot(drhice,  'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.9 0.9]); hold on;
                set(gca,'YTickLabel',{' '})   % Erase ylabels  
            end
            if size(drhocean,1)>1
                boxplot(drhocean,'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.1 0.9]);
                %set(gca,'XTickLabel',{' '})  % Erase xlabels   
                set(gca,'YTickLabel',{' '})   % Erase ylabels 
            end
            ytix = {'1','2','3','4','5','6','7','8','9','10'};   %  labels
            ytixloc = [1:10];      %  label locations
            set(gca,'YTickMode','auto','YTickLabel',ytix,'YTick',ytixloc);
            axis([-0.5 0.5 0 7]);
            h1=text(-9.5,6.5,daystr(4:11))       ;set(h1,'fontsize',12,'color','k');
            h2=text(-9.5,5.8,'over sea-ice')     ;set(h2,'fontsize',12,'color',[0.1 0.9 0.9]);
            h3=text(-9.5,5.5,'over open ocean');set(h3,'fontsize',12,'color',[0.1 0.1 0.9]);
            fi1=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'deltaRH_ana_air')];
            save_fig(111,fi1,false);
            close(111);
            
        %% continue and save figure 121
        %-------------------------------%
            figure(121)
            hold on;
            if size(dmmrice,1)>1
                boxplot(dmmrice,  'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.9 0.9]); hold on;
                set(gca,'YTickLabel',{' '})   % Erase ylabels  
            end
            if size(dmmrocean,1)>1
                boxplot(dmmrocean,'boxstyle','filled','orientation','horizontal',...
                                'positions',air.(starfieldNames{i,:}).(profstr).ana.zalt_mean(1:p)/1000,...
                                'colors',[0.1 0.1 0.9]);
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
            fi1=[strcat('F:\ARISE\ArcticCRF\figures\', daystr(4:11),'deltaMMR_ana_air')];
            save_fig(121,fi1,false);
            close(121);
            
            %% save derived quantities
            % save profile difference
            air.(starfieldNames{i,:}).dtice  = dtice;
            air.(starfieldNames{i,:}).dtocean = dtocean;
            air.(starfieldNames{i,:}).dmmrice  = dmmrice;
            air.(starfieldNames{i,:}).dmmrocean = dmmrocean;
            % save surface temperatures
            air.(starfieldNames{i,:}).tsurf_ice_mean    = nanmean(tsurf_ice);
            air.(starfieldNames{i,:}).tsurf_ice_std     = nanstd( tsurf_ice);
            air.(starfieldNames{i,:}).tsurf_ocean_mean  = nanmean(tsurf_ocean);
            air.(starfieldNames{i,:}).tsurf_ocean_std   = nanstd( tsurf_ocean);
            
            air.(starfieldNames{i,:}).dtsurf_ice_mean    = nanmean(dtsurf_ice);
            air.(starfieldNames{i,:}).dtsurf_ice_std     = nanstd( dtsurf_ice);
            air.(starfieldNames{i,:}).dtsurf_ocean_mean  = nanmean(dtsurf_ocean);
            air.(starfieldNames{i,:}).dtsurf_ocean_std   = nanstd( dtsurf_ocean);

    end     % end if profnum>0
    
    clear ice;
     
end         % end number of days

        %% save air data for plotting profiles
        %  create dat array to save
        cldwcm (cldwcm<0) = 0;
        cldicm (cldicm<0) = 0;
        dat2save = [yearm monthm daym profilem lonm latm zm staticPm iceconcm iceconc_stdm windSm windDm Tm Thetam RHm Qm SatVPwaterm SatVPicem cldflagm cldwcm cldicm];
        fi2sav   = ['C:\Users\msegalro.NDC\Documents\R\ArcticCRE\data\air_atmProfiles_4plotting_created_on_' datestr(now,'yyyy-mm-dd') '.txt'];
        save(fi2sav, '-ASCII','dat2save');
               
        %% save ana and interpolated air data for plotting profiles
        %  create dat array to save

        dat2save = [str2num(date) profile nprof lon_mean lon_std lat_mean lat_std iceconc_mean iceconc_std plev galt_mean galt_std zalt_mean o3mmr_mean o3mmr_std ...
                    omega_mean omega_std qice_mean qice_std qliq_mean qliq_std qwv_mean qwv_std rh_mean rh_std t_mean t_std theta_mean ...
                    winds_mean windd_mean cflag_ana ctype_ana cth_ana cbh_ana ctk_ana cwc_ana lwc_ana iwc_ana cflag_air cwc_air cic_air airTinterp airMMRinterp airRHinterp ...
                    airThetainterp s700_ana s850_ana s925_ana s700_air s850_air s925_air];
        fi2sav   = ['C:\Users\msegalro.NDC\Documents\R\ArcticCRE\data\ana_and_interp_air_atmProfiles_4plotting_created_on_' datestr(now,'yyyy-mm-dd') '.txt'];
        save(fi2sav, '-ASCII','dat2save');
        
        
%% add plots and save entire struct for further processing

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
    nProfiles = max(air.(starfieldNames{i,:}).prof.nprof);
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
            if air.(starfieldNames{i,:}).(profstr).iceconc_mean > 15
                if ~isnan(air.(starfieldNames{i,:}).(profstr).t850ind)
                air.(starfieldNames{i,:}).(profstr).LTS =...
                      air.(starfieldNames{i,:}).(profstr).theta(air.(starfieldNames{i,:}).(profstr).t850ind) - ...
                      air.(starfieldNames{i,:}).tsurf_ice_mean;% these numbers need to be calculated as nearest when staticP is available
                LTSice = [LTSice ;air.(starfieldNames{i,:}).(profstr).LTS ];
                end
                if ~isnan(air.(starfieldNames{i,:}).(profstr).t925ind)
                air.(starfieldNames{i,:}).(profstr).LLSS=...
                      air.(starfieldNames{i,:}).(profstr).theta(air.(starfieldNames{i,:}).(profstr).t925ind) - ...
                      air.(starfieldNames{i,:}).tsurf_ice_mean;% these numbers need to be calculated as nearest when staticP is available
                
                LLSSice= [LLSSice;air.(starfieldNames{i,:}).(profstr).LLSS];
                end
                anaLTSice  = [anaLTSice;  air.(starfieldNames{i,:}).(profstr).ana.theta850];
                anaLLSSice = [anaLLSSice; air.(starfieldNames{i,:}).(profstr).ana.theta925];
                % store reanalysis data
                % ana_theta850
            else
                if ~isnan(air.(starfieldNames{i,:}).(profstr).t850ind)
                air.(starfieldNames{i,:}).(profstr).LTS =...
                      air.(starfieldNames{i,:}).(profstr).theta(air.(starfieldNames{i,:}).(profstr).t850ind) - ...
                      air.(starfieldNames{i,:}).tsurf_ocean_mean;% these numbers need to be calculated as nearest when staticP is available
                      LTSocean = [LTSocean ;air.(starfieldNames{i,:}).(profstr).LTS ];
                end
                if ~isnan(air.(starfieldNames{i,:}).(profstr).t925ind)
                air.(starfieldNames{i,:}).(profstr).LLSS=...
                      air.(starfieldNames{i,:}).(profstr).theta(air.(starfieldNames{i,:}).(profstr).t925ind) - ...
                      air.(starfieldNames{i,:}).tsurf_ocean_mean;% these numbers need to be calculated as nearest when staticP is available
                
                LLSSocean= [LLSSocean;air.(starfieldNames{i,:}).(profstr).LLSS];
                end
                anaLTSocean  = [anaLTSocean;  air.(starfieldNames{i,:}).(profstr).ana.theta850];
                anaLLSSocean = [anaLLSSocean; air.(starfieldNames{i,:}).(profstr).ana.theta925];
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
    errorbar(1:length(label_str),LTSice_mean  ,LTSice_std,  'p','color'  ,[0.8 0.1 0.9],'markerfacecolor',[0.8 0.1 0.9],'markersize',10);  hold on;
    errorbar(1:length(label_str),LTSocean_mean,LTSocean_std,'p','color'  ,[0.1 0.1 0.8],'markerfacecolor',[0.1 0.1 0.8],'markersize',10);  hold on;
    % plot MERRA
    errorbar(1:length(label_str),anaLTSice_mean  ,anaLTSice_std,  'p','color'  ,[0.1 0.9 0.5],'markerfacecolor',[0.1 0.9 0.5],'markersize',6);  hold on;
    errorbar(1:length(label_str),anaLTSocean_mean,anaLTSocean_std,'p','color'  ,[0.8 0.8 0.1],'markerfacecolor',[0.8 0.8 0.1],'markersize',6);  hold on;
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
    plot(1:length(label_str),s1,':','color'  ,[0.8 0.1 0.9],'linewidth',1.5);hold on;
    plot(1:length(label_str),s2,':','color'  ,[0.1 0.1 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s3,':','color'  ,[0.1 0.9 0.5],'linewidth',1.5);hold on;
    plot(1:length(label_str),s4,':','color'  ,[0.8 0.8 0.1],'linewidth',1.5);hold on;
    
    ylabel('\theta_{850} - LML [K]','fontsize',12);
    h11=text(8,65,'C130 over sea-ice')  ; set(h11,'fontsize',12,'color',[0.8 0.1 0.9]);
    h22=text(8,61,'C130over open ocean'); set(h22,'fontsize',12,'color',[0.1 0.1 0.8]);
    h33=text(8,57,'GEOS-FP over sea-ice')  ;set(h33,'fontsize',12,'color',[0.1 0.9 0.5]);
    h44=text(8,53,'GEOS-FP open ocean');    set(h44,'fontsize',12,'color',[0.8 0.8 0.1]);
    
    fi2=[strcat('F:\ARISE\ArcticCRF\figures\',label_str{:,1},'_',label_str{:,end},'LTScompareGEOS_FP_iceconc15')];
    save_fig(13,fi2,false);
    close(13);
    
%% create and save figure 14
    % plot comparison for each day
    % LLSS
    figure(14);
    % plot C130
    errorbar(1:length(label_str),LLSSice_mean  ,LLSSice_std,  'p','color'  ,[0.8 0.1 0.9],'markerfacecolor',[0.8 0.1 0.9],'markersize',10);  hold on;
    errorbar(1:length(label_str),LLSSocean_mean,LLSSocean_std,'p','color'  ,[0.1 0.1 0.8],'markerfacecolor',[0.1 0.1 0.8],'markersize',10);  hold on;
    % plot MERRA
    errorbar(1:length(label_str),anaLLSSice_mean  ,anaLLSSice_std,  'p','color'  ,[0.1 0.9 0.5],'markerfacecolor',[0.1 0.9 0.5],'markersize',6);  hold on;
    errorbar(1:length(label_str),anaLLSSocean_mean,anaLLSSocean_std,'p','color'  ,[0.8 0.8 0.1],'markerfacecolor',[0.8 0.8 0.1],'markersize',6);  hold on;
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
    plot(1:length(label_str),s1,':','color'  ,[0.8 0.1 0.9],'linewidth',1.5);hold on;
    plot(1:length(label_str),s2,':','color'  ,[0.1 0.1 0.8],'linewidth',1.5);hold on;
    plot(1:length(label_str),s3,':','color'  ,[0.1 0.9 0.5],'linewidth',1.5);hold on;
    plot(1:length(label_str),s4,':','color'  ,[0.8 0.8 0.1],'linewidth',1.5);hold on;
    
    ylabel('\theta_{925} - LML [K]','fontsize',12);
    h11=text(8,65,'C130 over sea-ice')  ; set(h11,'fontsize',12,'color',[0.8 0.1 0.9]);
    h22=text(8,61,'C130over open ocean'); set(h22,'fontsize',12,'color',[0.1 0.1 0.8]);
    h33=text(8,57,'GEOS-FP over sea-ice')  ;set(h33,'fontsize',12,'color',[0.1 0.9 0.5]);
    h44=text(8,53,'GEOS-FP open ocean');    set(h44,'fontsize',12,'color',[0.8 0.8 0.1]);
    
    fi2=[strcat('F:\ARISE\ArcticCRF\figures\',label_str{:,1},'_',label_str{:,end},'LLSScompare_GEOS_FP_iceconc15')];
    save_fig(14,fi2,false);
    close(14);
    
%% save processed air struct
file2save=[strcat(arisedir,filename,'_w_anacompare_w_consolidatedcloudsAir_andModel', datestr(now,'yyyy-mm-dd'), '.mat')];
save(file2save,'-struct','air');
        

%% function to convert geopotential height to geometrical height [m]
function z=geo2z(geoH)
    r    = 6378*1000;       % Earth radius
    z = (geoH*r)./(r - geoH); 
return;

%% function to convert mass mixing ratio [kg/kg] to molec/cm3
function molec=mmr2molec(plev,mmr,gas,temp)

% mmr is in Kg/Kg
% gas is gas name string
% temp is Temperature in K

Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)
mwair= 28.97;           % Molecular weight of dry air

if strcmp(gas,'h2o')
    mw = 18.01;           % Molecular weight of water
elseif strcmp(gas,'o3')
    mw = 48.00;           % Molecular weight of ozone
end

% calculate
ppmv = mmr*(mwair/mw)*1e6;              % conversion from mmr to ppm
partp= ppmv.*plev*1e-6;                 % conversion from ppmv to hPa
dens = partp*100./(R*temp);             % row=P/RT: conversion from hPa to kg/m3
molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3

return;      

%% function to convert air levels to number of air/o2 molec
function molec=air2molec(plev,gas,temp)

% mmr is in Kg/Kg
% gas is gas name string
% temp is Temperature in K

Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)

if strcmp(gas,'air')
    mw = 28.97;                             % Molecular weight of dry air
    dens = plev*100 ./(R*temp);             % row=P/RT: conversion from hPa to kg/m3
    molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3
    
elseif strcmp(gas,'o2')
    mw = 32.00;                             % Molecular weight of oxygen
    dens = plev*0.21*100 ./(R*temp);        % row=P/RT: conversion from hPa to kg/m3
    molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3
end

return;  

%% function to convert air levels to number of air/o2 molec
function dens=air2dens(plev,temp)

% temp is Temperature in K
% plev is pressure level in hPa

R    = 287.05;                          % gas constant J/(kgK)

dens = plev*100 ./(R*temp);             % row=P/RT: conversion from hPa to kg/m3

return;  
  