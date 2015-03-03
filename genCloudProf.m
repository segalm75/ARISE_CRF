%% Details of the function:
%  NAME:
% genCloudProf
%----------------------------------
% PURPOSE:
%  - combines 4star/ssfr/in-situ data to generate cloud field for a
%    selected vertical profile
%
% CALLING SEQUENCE:
%   [c] = genCloudProf(prof,dprof,allprof,daystr,doy);
%
% INPUT:
%  - prof is a struct containing vertical profile data (this is input from
%    createARISEAtmProf.m)
%  - dprof is profile number within each day
%  - allprof is general profile number
%  - daystr is date
% 
% 
% OUTPUT:
%  creates cloud data according to selected profile
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%
% NEEDED FILES/INPUT:
%
% EXAMPLE:
%  - [c] = genCloudProf(prof,dprof,allprof,daystr,doy)
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Feb-23,2015  
% MS, 2015-03-03, added profin variable to account for serial profile number
%     adjusted dprof, allprof for saving purposes and data extraction
% -------------------------------------------------------------------
%% function routine
function [c] = genCloudProf(prof,dprof,allprof,daystr,doy)

version_set('1.0');
startup_plotting;

%% cloud flag

% plot nCDP to determine threshold
% figure(22);
% plot(prof.time,...
%      prof.nCDP_cm3,'.b');hold on;
% plot(prof.time,...
%      nanmean(prof.nCDP_cm3),'--g');hold on;
% figure(222);hist(prof.nCDP_cm3,[1:2:100]);
% nCDPq = quantile(prof.nCDP_cm3,[.025 .25 .50 .75 .975]); % Summary of nCDP
nCDPt = 2;% set constant for now; nCDPq(4);

% generate cloud flag
cloudflag = logical(prof.nCDP_cm3>nCDPt);
% figure(2222);
% plot(prof.time(cloudflag==0),cloudflag(cloudflag==0),'.','color','b');hold on;
% plot(prof.time(cloudflag==1),cloudflag(cloudflag==1),'.','color',[0.5 0.5 0.5]);hold on;
% legend('no cloud','cloud');

cdiff = diff(double(cloudflag));
if cloudflag(1)  ==1 cdiff(1)  = 1; end;
if cloudflag(end)==1 cdiff(end)=-1; end;
cind     = find(cdiff~=0);
% define cloud boundaries
if mod(length(cind),2)==0
    cstart = NaN(length(cind)/2,1);
    cend   = NaN(length(cind)/2,1);
else
    cstart = NaN(floor(length(cind)/2),1);
    cend   = NaN(floor(length(cind)/2),1);
end
k=0;
for i=1:2:(length(cind)-1)
    k=k+1;
    % check for consistancy of at least 5 seconds
    % cloud begin
    if cind(i)+5 < length(cloudflag)
        if (cdiff(cind(i))==1 && sum(cloudflag([cind(i)+1:cind(i)+5]))==5)
            cstart(k) = cind(i)+1;
        end
    end
    % cloud end
    if cind(i+1)-4 >0
        if cdiff(cind(i+1))==-1 && sum(cloudflag([cind(i+1)-4:cind(i+1)]))==5
            cend(k)   = cind(i+1);
        end
    end
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
    c.cldwc    = [];
    c.cldflag  = zeros(length(prof.time),1);
end

%% cloud geometry
if c.cldnum>0
    %  calculate cloud thickness [m]
    c.cldthick = NaN(c.cldnum,1);
    c.cldtop   = NaN(c.cldnum,1);
    c.cldbot   = NaN(c.cldnum,1);
    if strcmp(prof.direction,'descend')
        c.cldthick = prof.z(c.cldstart) - prof.z(c.cldend);
        c.cldtop   = prof.z(c.cldstart);
        c.cldbot   = prof.z(c.cldend);
    else
        c.cldthick = prof.z(c.cldend) - prof.z(c.cldstart);
        c.cldtop   = prof.z(c.cldend);
        c.cldbot   = prof.z(c.cldstart);
    end

%% refine cloud flag
cind       = logical(c.cldthick>0);
c.cldnum   = sum(cind);
c.cldstart = c.cldstart(cind==1);
c.cldend   = c.cldend(  cind==1);
c.cldthick = c.cldthick(cind==1);
c.cldtop   = c.cldtop(  cind==1);
c.cldbot   = c.cldbot(  cind==1);
c.cldflag  = zeros(length(prof.time),1);
for i=1:c.cldnum
        c.cldflag(c.cldstart(i):c.cldend(i)) = 1;
end
%% cloud Reff
reff = {'05','10','20','30','40','50'}; % in um
reffum = [5 10 20 30 40 50];
c.cldreff = zeros(c.cldnum,1);
partsum   = zeros(c.cldnum,length(reff));

% sum and avg normalized particle number per cloud
for i=1:c.cldnum
    for j=1:length(reff)
        pstr = strcat('CDP',reff{:,j},'um');
        partsum(i,j) = sum(prof.(pstr)(c.cldstart(i):c.cldend(i)));
    end
    c.cldreff(i) = round(sum((partsum(i,:).*reffum))/sum(partsum(i,:)));
end
clear partsum;
%% cloud phase
% !!! this is only by Reff/temp; next stage to add 4star data/RH over ice
c.cldphase = zeros(c.cldnum,1);
% 0 is water (or supercooled), 1 is ice, 2 is mixed
for i=1:c.cldnum
    if c.cldreff(i) <= 30 && nanmean(prof.Static_AirT(c.cldstart(i):c.cldend(i)))> -20
        c.cldphase(i) = 0;
    elseif c.cldreff(i) > 30 && nanmean(prof.Static_AirT(c.cldstart(i):c.cldend(i)))< -20
        c.cldphase(i) = 1;
    elseif c.cldreff(i) >= 30 && nanmean(prof.Static_AirT(c.cldstart(i):c.cldend(i)))> -20 &&...
                                 nanmean(prof.Static_AirT(c.cldstart(i):c.cldend(i)))<  0
        c.cldphase(i) = 2;
    else
        c.cldphase(i) = 0;
    end
end
%% cloud water content (ice/water)
% the density of liquid water as 1 g/cm3 = (1 g/cm3)(0.1 cm/mm) = 0.1 g/(cm2-mm)
cnvrt = 0.1;
c.cldwc = zeros(c.cldnum,1);
for i=1:c.cldnum
    if strcmp(prof.direction,'descend')
        wp = nanmean(prof.LWP_mm(c.cldend(i)-5:c.cldend(i)));    % in mm (from GVR)
    else
        wp = nanmean(prof.LWP_mm(c.cldstart(i):c.cldstart(i)+5));% in mm (from GVR)
    end
    if wp==0 wp=0.04; end;
    c.cldwc(i)= (wp*cnvrt*100^2)/c.cldthick(i);                    % this is water path [g/m3]
end

end% if cloud

%% cloud above
% 4star/ssfr flags
if ~isnan(c.cldstart)
    if strcmp(prof.direction,'descend')
        [mini,I] = min(c.cldstart);
        ind = 1:mini;
    else
        [mini,I] = max(c.cldend);
        ind = mini:length(prof.z);
    end
else
     if strcmp(prof.direction,'descend')
         if length(prof.z)>=30
            ind = 1:30;% 30 sec during descend
         else
            ind = 1:length(prof.z);
         end
     else
         if length(prof.z)>=30
            ind = [length(prof.z)-30:length(prof.z)];
         else
            ind = [length(prof.z)-length(prof.z)+1:length(prof.z)];
         end
    end
end
    if      mean(prof.starStr(ind)) <= 1 &&...
            mean(prof.starMd(ind))  <=1
            c.cldabove = 0;% clear skies above highest cloud
    elseif mean(prof.starStr(ind)) >  1 &&...
           mean(prof.starStr(ind)) <= 2 &&...
           mean(prof.starMd(ind))  >  1 &&...
           mean(prof.starMd(ind))  <  7
           c.cldabove = 0; % clear sky; this is aerosol scan
    elseif mean(prof.starStr(ind)) <= 1 &&...
           mean(prof.starMd(ind)) > 1
           c.cldabove = 1; % thin cloud above (mode 7 exist)
    elseif mean(prof.starStr(ind)) >  1 &&...
           mean(prof.starStr(ind)) <= 2 &&...
           mean(prof.starMd(ind))  >  7 &&...
           mean(prof.starMd(ind))  <= 10
           c.cldabove = 2; % cloud above (zenith/cloud scan)
    else
           c.cldabove = 0;
    end
    
%% save platform lowest altitude
if strcmp(prof.direction,'descend')
    c.c130sur_alt = prof.z(end)/1000;
else
    c.c130sur_alt = prof.z(1)/1000;
end
%% clouds below
%LVIS/ssfr data/flags/nadir video

% open below (by video inspection) cloud flag below
% bfile = 'F:\ARISE\ArcticCRF\METdata\cldbelow_nadircam_20150303.dat';
bfile = 'F:\ARISE\ArcticCRF\METdata\clouds_surface_below20150303.txt';
bflag = importdata(bfile);
c.cldbelow = bflag.data(allprof,2);
% add below cloud to cloud layers
if c.cldbelow==1
    c.cldnum = c.cldnum + 1;
    if strcmp(prof.direction,'descend')
        c.cldstart = [c.cldstart ;length(prof.z)-5];
        c.cldend   = [c.cldend   ;length(prof.z)  ];
        c.cldthick = [c.cldthick ;prof.z(c.cldstart(end)) - 0];%bottom of aircraft to ground
        c.cldtop   = [c.cldtop   ;prof.z(end)];
        c.cldbot   = [c.cldbot   ;0];
    else
        c.cldstart = [c.cldstart ;1];
        c.cldend   = [c.cldend   ;5];
        c.cldthick = [c.cldthick ;prof.z(c.cldend(end)) - 0]; %bottom of aircraft to ground
        c.cldtop   = [c.cldtop   ;prof.z(1)];
        c.cldbot   = [c.cldbot   ;0];
    end
    
    % reff
    reff = {'05','10','20','30','40','50'}; % in um
    reffum = [5 10 20 30 40 50];
     for j=1:length(reff)
        pstr = strcat('CDP',reff{:,j},'um');
        partsum(j) = sum(prof.(pstr)(c.cldstart(end):c.cldend(end)));
     end
     if sum(partsum)>0
        c.cldreff = [c.cldreff ;round(sum((partsum.*reffum))/sum(partsum))];
     else
        c.cldreff = [c.cldreff ;10];% default value
     end
     
     % phase
     if c.cldreff(end) <= 30 && nanmean(prof.Static_AirT(c.cldstart(end):c.cldend(end)))> -20
        c.cldphase = [c.cldphase ;0];
     elseif c.cldreff(end) > 30 && nanmean(prof.Static_AirT(c.cldstart(end):c.cldend(end)))< -20
        c.cldphase = [c.cldphase ;1];
     elseif c.cldreff(end) >= 30 && nanmean(prof.Static_AirT(c.cldstart(end):c.cldend(end)))> -20 &&...
                                 nanmean(prof.Static_AirT(c.cldstart(end):c.cldend(end)))< 0
        c.cldphase = [c.cldphase ;2];
     else
        c.cldphase = [c.cldphase ;0];
     end
     
     % wc
      if strcmp(prof.direction,'descend')
        wp = nanmean(prof.LWP_mm(c.cldend(end)-5:c.cldend(end)));    % in mm (from GVR)
      else
        wp = nanmean(prof.LWP_mm(c.cldstart(end):c.cldstart(end)+5));% in mm (from GVR)
      end
      if wp==0 wp=0.04; end;
      c.cldwc = [c.cldwc ;(wp*0.1*100^2)/c.cldthick(end)];
      
      % cloudflag
      if strcmp(prof.direction,'descend')
            c.cldflag(end-5:end)=1;
      else
            c.cldflag(1:5)=1;
      end
     
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
           
           % write data to file
           filen=['C:\cygwin\home\msegalro\libradtran\input\CRF\arise\MERRABASE\MERRABASEdatin.txt'];
           disp( ['Writing to file: ' filen]);
           [filepath filename ext] = fileparts(prof.ana.atmfile);
           atmosfile = strcat(filename, ext);
           for i=1:c.cldnum
               if prof.iceconc >15
                   s_albedo_file = 'albedo_ice.dat';%strcat('F:\ARISE\ArcticCRF\METdata\albedo\albedo_ice.dat');
                   w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                   line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                          ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                          ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(c.cldtop(i)/1000)...
                          ' cloudbot '  num2str(c.cldbot(i)/1000) ' cloudthick ' num2str(c.cldthick(i)/1000) ' cloudreff ' num2str(c.cldreff(i)) ' cloudphase ' num2str(c.cldphase(i))...
                          ' cloudwc '   num2str(c.cldwc(i)) ' c130_sur_alt ' num2str(c.c130sur_alt)}];
               else
                   s_albedo_file = 'albedo_water.dat';
                   w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                   line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                          ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                          ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(c.cldtop(i)/1000)...
                          ' cloudbot '  num2str(c.cldbot(i)/1000) ' cloudthick ' num2str(c.cldthick(i)/1000) ' cloudreff ' num2str(c.cldreff(i)) ' cloudphase ' num2str(c.cldphase(i))...
                          ' cloudwc '   num2str(c.cldwc(i)) ' c130_sur_alt ' num2str(c.c130sur_alt)}];
               end
               
               dlmwrite(filen,line,'-append','delimiter','');

               clear line;
           end

return;
