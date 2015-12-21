%% Details of the function:
%  NAME:
% genModCloudProf
%----------------------------------
% PURPOSE:
%  - generates modeled cloud profiles from model inputs/interpolated
%    aircraft input
%    and cloud summary file for RT runs input
%    
%
% CALLING SEQUENCE:
%   [c] = genModCloudProf(model,prof,dprof,allprof,daystr,doy,datsource,atmsource,ctype,cReff);
%
% INPUT:
%  - model is the model source (e.g., 'MERRA2','GEOS-FP' etc.)
%  - prof is a struct containing vertical profile data (from either model/aircraft)
%  - dprof is profile number within each day - integer
%  - allprof is consecutive profile number (from total) - integer
%  - daystr is date ('yyyy-mm-dd') - string
%  - doy is day of yesr (string)
%  - datsource is data source ('mod'/'air');
%  - atmsource is atmospheric profile source ('mod/'air')
%  - ctype is cloud type 'wc'/'ic'/'mx' or ice as water (iawc),
%    mxflp(ice/liq is flipped)
%  - cReff is Reff 10,15,20 for 'mx' both ice and liquid are similar
%  - 0 is original Reff from in-situ
% 
% 
% OUTPUT:
%  creates cloud data according to selected profile
%  generates cloud profile for RT runs
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%
% NEEDED FILES/INPUT:
%
% EXAMPLE:
%  - [c] = genModCloudProf(model,prof,dprof,allprof,daystr,doy,datsource,atmsource,ctype,cReff)
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Oct-12,2015  
% Based on genAirCloudProf
% Modified: 2015-10-26, MS: added function capabilities to accept more
%                           options
% Modified: 2015-11-11, MS: changed cloud bot to 0 at LML
% Modified: 2015-11-13, MS: chnaged Reff definition to 10/15 etc.
% Modified: 2015-11-23, MS: changed Reff definition back to 10/10 etc.
% ---------------------------------------------------------------------------
%% function routine
function [c] = genModCloudProf(model,prof,dprof,allprof,daystr,doy,datsource,atmsource,ctype,cReff)

version_set('1.1');
startup_plotting;

%% cloud flag and cloud water content assignments

% constrain model vertical data to correspond with aircraft

   hl = find(isnan(prof.cflagInterp));

if strcmp(datsource,'air') && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'mx'))
    cloudflag = (prof.cflagInterp)';
    cloudwc   = (prof.cwcInterp)';
    cloudic   = (prof.cicInterp)';
elseif strcmp (datsource,'mod') && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'mx'))
    cloudflag = prof.ana.cflag; cloudflag(hl(1):end) = 0;
    cloudwc   = prof.ana.lwc;   cloudwc(hl(1):end)   = 0;
    cloudic   = prof.ana.iwc;   cloudic(hl(1):end)   = 0;
elseif strcmp (datsource,'mod') && (strcmp(ctype,'iawc') || strcmp(ctype,'mxflp'))
    cloudflag = prof.ana.cflag; cloudflag(hl(1):end) = 0;
    cloudwc   = prof.ana.iwc;   cloudwc(hl(1):end)   = 0;
    cloudic   = prof.ana.lwc;   cloudic(hl(1):end)   = 0;
end

cdiff = diff(double(cloudflag));
if cloudflag(1)  ==1 cdiff(1)  = 1; end;
if cloudflag(end)==1 cdiff(end)=-1; end;
cind     = find(cdiff~=0&~isnan(cdiff));


% define cloud boundaries
if mod(length(cind),2)==0
    cstart = NaN(length(cind)/2,1);
    cend   = NaN(length(cind)/2,1);
elseif length(cind)==1
    cstart = NaN(1,1);
    cend   = NaN(1,1);
else
    cstart = NaN(floor(length(cind)/2),1);
    cend   = NaN(floor(length(cind)/2),1);
end

k=0;

if length(cind) > 1

        for i=1:2:(length(cind)-1)

            k=k+1;

            % cloud begin
            if cind(i) < length(cloudflag)
                if  cdiff(cind(i))==1 && cind(i) > 1
                    cstart(k) = cind(i)+1;
                elseif cdiff(cind(i))==1 && cind(i) == 1
                    cstart(k) = cind(i);
                end
            end
            % cloud end
            if cind(i+1) >0
                if cdiff(cind(i+1))==-1
                    cend(k)   = cind(i+1);
                end
            end
        end
        
        
elseif length(cind)==1
         k = k+1;
         cstart(k) = cind(1);
         cend(k)   = cind(1)+1;

end

%  assign cloud start and end indices
c.cldnum   = sum(~isnan(cstart)&~isnan(cend));
if c.cldnum>0
    c.cldstart = cstart(~isnan(cstart)&~isnan(cend));
    c.cldend   = cend(  ~isnan(cstart)&~isnan(cend));
else
    c.cldstart = NaN;
    c.cldend   = NaN;
end

if c.cldnum==0
    c.cldstart = [];
    c.cldend   = [];
    c.cldthick = [];
    c.cldtop   = [];
    c.cldbot   = [];
    c.cldphase = [];
    c.cldreff  = [];
    c.cldreffWC  = [];
    c.cldreffIC  = [];
    c.cldwc    = [];
    c.cldic    = [];
    c.cldflag  = zeros(length(prof.time),1);
end

%% cloud geometry
if c.cldnum>0
    %  calculate cloud thickness [m]
    c.cldthick = NaN(c.cldnum,1);
    c.cldtop   = NaN(c.cldnum,1);
    c.cldbot   = NaN(c.cldnum,1);

    c.cldthick = prof.ana.zalt_mean(c.cldend) - prof.ana.zalt_mean(c.cldstart);
    c.cldtop   = double(prof.ana.zalt_mean(c.cldend));
    c.cldbot   = double(prof.ana.zalt_mean(c.cldstart));
    
    % change cloud bot from LML to 0
    
    if c.cldstart(1)==1
        c.cldbot(1) = 0;
    end
    
    % check if ctop & cbot are in same layer
    
    levels = double(diff(prof.ana.zalt_mean));
    
    for kk=1:c.cldnum
        
        if c.cldtop(kk)==c.cldbot(kk)
            % compute again accoding to layer thickness
            c.cldbot(kk)   = double(c.cldbot(kk));
            c.cldtop(kk)   = double(c.cldbot(kk) + levels(c.cldstart(kk)-1));
            c.cldthick(kk) = double(c.cldtop(kk) - c.cldbot(kk)); % or: levels(c.cldstart(kk)-1)
        end
        
    end


%% cloud Reff

    c.cldreff   = NaN(c.cldnum,1);
    
    if     cReff>=0  && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'iawc'))
                    c.cldreff = NaN(c.cldnum,1);
    elseif cReff>=10 && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
                    c.cldreffWC = NaN(c.cldnum,1);
                    c.cldreffIC = NaN(c.cldnum,1);
    end
    
    % this is for 10/15/20 for both ice and liquid
    
        for i=1:c.cldnum
            if     cReff==0  && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'iawc'))
                        c.cldreff(i) = nanmean(prof.cld.cldreff);
            elseif cReff>=10 && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'iawc'))
                        c.cldreff(i) = cReff;
            elseif cReff==0  && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
                        c.cldreffWC(i) = nanmean(prof.cld.cldreff);
                        c.cldreffIC(i) = nanmean(prof.cld.cldreff);
            elseif cReff>=10 && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
                        c.cldreffWC(i) = cReff;
                        c.cldreffIC(i) = cReff;
            end
        end
    
    % this is for 10/15, 15/20, 20/30 Reff combinations for liquid/ice
    
%     for i=1:c.cldnum
%         if     cReff==0  && (strcmp(ctype,'wc') || strcmp(ctype,'ic') || strcmp(ctype,'iawc'))
%                     c.cldreff(i) = nanmean(prof.cld.cldreff);
%         elseif cReff==0  && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
%                     c.cldreffWC(i) = nanmean(prof.cld.cldreff);
%                     c.cldreffIC(i) = nanmean(prof.cld.cldreff);
%         elseif cReff==10 && (strcmp(ctype,'wc') || strcmp(ctype,'iawc'))
%                     c.cldreff(i) = cReff;
%         elseif cReff==10 && (strcmp(ctype,'ic'))
%                     c.cldreff(i) = 15;
%         elseif cReff==10 && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
%                     c.cldreffWC(i) = cReff;
%                     c.cldreffIC(i) = 15;
%         elseif cReff==15 && (strcmp(ctype,'wc') || strcmp(ctype,'iawc'))
%                     c.cldreff(i) = cReff;
%         elseif cReff==15 && (strcmp(ctype,'ic'))
%                     c.cldreff(i) = 20;
%         elseif cReff==15 && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
%                     c.cldreffWC(i) = cReff;
%                     c.cldreffIC(i) = 20;
%         elseif cReff==20 && (strcmp(ctype,'wc') || strcmp(ctype,'iawc'))
%                     c.cldreff(i) = cReff;
%         elseif cReff==20 && (strcmp(ctype,'ic'))
%                     c.cldreff(i) = 30;
%         elseif cReff==20 && (strcmp(ctype,'mx') || strcmp(ctype,'mxflp'))
%                     c.cldreffWC(i) = cReff;
%                     c.cldreffIC(i) = 30;
%         end
%     end



%% cloud phase and cloud water content
% !!! this is only by pre-assignment to test sensitivity (w/i/m)

    c.cldphase = zeros(c.cldnum,1);
    c.cldwc    = zeros(c.cldnum,1);
    c.cldic    = zeros(c.cldnum,1);

    % 0 is water (or supercooled), 1 is ice, 2 is mixed

    for i=1:c.cldnum

        if strcmp(ctype,'wc') || strcmp(ctype,'iawc')
            c.cldphase(i) = 0;
            if c.cldstart(i)>2
                c.cldwc(i)    = nanmean(cloudwc(c.cldstart(i)-2:c.cldstart(i)+2));
            else
                c.cldwc(i)    = nanmean(cloudwc(c.cldstart(i):c.cldstart(i)+2));
            end
        elseif strcmp(ctype,'ic')
            c.cldphase(i) = 1;
            if c.cldstart(i)>2
                c.cldic(i)    = nanmean(cloudic(c.cldstart(i)-2:c.cldstart(i)+2));
            else
                c.cldic(i)    = nanmean(cloudic(c.cldstart(i):c.cldstart(i)+2));
            end
        else
            c.cldphase(i) =2;
            if c.cldstart(i)>2
                c.cldwc(i)    = nanmean(cloudwc(c.cldstart(i)-2:c.cldstart(i)+2));
                c.cldic(i)    = nanmean(cloudic(c.cldstart(i)-2:c.cldstart(i)+2));
            else
                c.cldwc(i)    = nanmean(cloudwc(c.cldstart(i):c.cldstart(i)+2));
                c.cldic(i)    = nanmean(cloudic(c.cldstart(i):c.cldstart(i)+2));
            end
        end

        if (c.cldwc(i)<=0 || isnan(c.cldwc(i))) c.cldwc(i) = 0.001; end;
        if (c.cldic(i)<=0 || isnan(c.cldic(i))) c.cldic(i) = 0.001; end; 

    end

end% if cloud

%% cloud above

   c.cldabove = prof.cld.cldabove;
   
%% cloud below

   c.cldbelow = prof.cld.cldbelow;
   
%% save platform lowest altitude

    if strcmp(prof.direction,'descend')
        c.c130sur_alt = prof.z(end)/1000;
    else
        c.c130sur_alt = prof.z(1)/1000;
    end
    
%% create data file as input to libRadTran routine
           % calculate sza at bottom of each profile
           day = datevec(daystr,'yyyymmdd');
           if strcmp(prof.direction,'descend')
               v = datevec(prof.time(end)/24);
               lat = prof.lat(end);
               lon = prof.lon(end);
           else
               v = datevec(prof.time(1)/24);
               lat = prof.lat(1);
               lon = prof.lon(1);
           end
           day(4:6) = v(4:6);
           
           [sza, az, soldst, ha, dec, el, am] = sunae(lat, lon, datenum(day));
           
          
           % assign params
           cldtop   = flipud(c.cldtop);
           cldbot   = flipud(c.cldbot);
           cldthick = flipud(c.cldthick);
           cldphase = flipud(c.cldphase);
           cldreff  = flipud(c.cldreff);
           
           %% write cloud params to file for RT
           %-----------------------------------
           % write data to file
           % this is for RT simulation input
           dirsave = 'C:\cygwin\home\msegalro\libradtran\input\CRF\arise\';
           filen=[dirsave, 'cloud4RT_datsource_', datsource, '_atmsource_', atmsource, '_', ctype, '_', num2str(cReff),'_for',model, '.txt'];
           disp( ['Writing to file: ' filen]);
           prof.iceconc = prof.iceconc_mean;
           if strcmp(atmsource,'air')
                [filepath filename ext] = fileparts(prof.airAtmInterp);
           else
                [filepath filename ext] = fileparts(prof.ana.atmfile);
           end
           
           atmosfile = strcat(filename, ext);
           
           if strcmp(ctype,'mx') || strcmp(ctype,'mxflp')
                       if prof.iceconc >15
                           s_albedo_file = 'albedo_ice.dat';%strcat('F:\ARISE\ArcticCRF\METdata\albedo\albedo_ice.dat');
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreffwc ' num2str(c.cldreffWC')  ' cloudreffic ' num2str(c.cldreffIC') ...
                                  ' cloudphase ' num2str(cldphase') ' cloudwc '   num2str(c.cldwc') ' cloudic '   num2str(c.cldic') '  '}];
                       else
                           s_albedo_file = 'albedo_water.dat';
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreffwc ' num2str(c.cldreffWC') ' cloudreffic ' num2str(c.cldreffIC') ...
                                  ' cloudphase ' num2str(cldphase') ' cloudwc '   num2str(c.cldwc') ' cloudic '   num2str(c.cldic') '  '}];
                       end
           elseif strcmp(ctype,'ic')
               
                        if prof.iceconc >15
                           s_albedo_file = 'albedo_ice.dat';%strcat('F:\ARISE\ArcticCRF\METdata\albedo\albedo_ice.dat');
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff')  ' cloudphase ' num2str(cldphase')...
                                  ' cloudic '   num2str(c.cldic') '  '}];
                       else
                           s_albedo_file = 'albedo_water.dat';
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff') ' cloudphase ' num2str(cldphase')...
                                  ' cloudic '   num2str(c.cldic') '  '}];
                       end
                       
           elseif strcmp(ctype,'wc') || strcmp(ctype,'iawc')
               
                        if prof.iceconc >15
                           s_albedo_file = 'albedo_ice.dat';%strcat('F:\ARISE\ArcticCRF\METdata\albedo\albedo_ice.dat');
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff')  ' cloudphase ' num2str(cldphase')...
                                  ' cloudwc '   num2str(c.cldwc') '  '}];
                       else
                           s_albedo_file = 'albedo_water.dat';
                           w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                           line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                                  ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                                  ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                                  ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff') ' cloudphase ' num2str(cldphase')...
                                  ' cloudwc '   num2str(c.cldwc') '  '}];
                       end
               
           end
               
               dlmwrite(filen,line,'-append','delimiter','');

               clear line;
           

return;

