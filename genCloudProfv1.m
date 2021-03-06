%% Details of the function:
%  NAME:
% genCloudProfv1
% this is the non-succesful version - don't use!
%-----------------------------------------------
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
%     changed writing format into data file (from line for each cloud to
%     line for each profile
% MS, 2015-03-13, changed nCDPt threshold from 2 to 0.1
%                 added 2 options to consolidate cloud layers
% ---------------------------------------------------------------------------
%% function routine
function [c] = genCloudProf(prof,dprof,allprof,daystr,doy)

version_set('1.0');
startup_plotting;

%% cloud flag

% plot nCDP to determine threshold
%   figure(22);
%   plot(prof.time,...
%        prof.nCDP_cm3,'.b');hold on;
%   plot(prof.time,...
%        nanmean(prof.nCDP_cm3),'--g');hold on;
%   figure;plot(prof.nCDP_cm3,'.g');
% figure(222);hist(prof.nCDP_cm3,[1:2:100]);
% nCDPq = quantile(prof.nCDP_cm3,[.025 .25 .50 .75 .975]); % Summary of nCDP
nCDPt = 0.1;% set constant for now; nCDPq(4);

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

% check for cloud layer by alt/time derivative
zdiff  = diff(smooth(double(prof.z)));   zdiff = [zdiff;zdiff(end)];
zdiff2 = diff(smooth(double(zdiff)));    zdiff2 = [zdiff2;zdiff2(end)];
tdiff = diff(smooth(double(prof.time)));tdiff = [tdiff;tdiff(end)];
dzdt  = zdiff./tdiff;

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
    c.cldwc(i)= (wp*cnvrt*100^2)/c.cldthick(i);                  % this is water path [g/m3]
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
%c.cldbelow = 0; % this is to test without cloud below
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
           
          
           %% consolidate cloud layers
           % sort clouds by high to low
           [cldtop ix]   = sort(c.cldtop,'descend');
           %[cldstart ix] = sort(c.cldstart,'descend');
            cldbot     = c.cldbot(ix);
            cldthick   = c.cldthick(ix);
            cldreff    = c.cldreff(ix);
            cldphase   = c.cldphase(ix);
            cldwc      = c.cldwc(ix);
            cldstart   = c.cldstart(ix);
            cldend     = c.cldend(ix);
            
           % consolidate cloud layers
            logic = ones(length(cldtop),1);
            cloud = [cldbot,cldtop,cldthick,cldstart,cldend,cldphase,cldreff,cldwc,logic];
            %clear cldtop cldbot cldstart cldend cldthick cldphase cldwc cldreff
%             for k=length(cloud)-1:-1:2
%                 if  (cloud(k+1,2) >= cloud(k,1))
%                     % top of lower cloud higher than bottom
%                     % of above cloud, then consolidate
%                     cloud(k,1)     = cloud(k+1,1);    % replace cloud bot
%                     cloud(k,4)     = cloud(k+1,4);    % replace all above line
%                     cloud(k+1,end)   = 0;             % don't keep
%                     %cloud(k-1,end) = 1;              % keep for next iteration
%                     
%                 end
%             end
           % check for thickness
           if ~isempty(cloud)
               for k=size(cloud,1):-1:1
                   if cloud(k,3)<50 && cloud(k,1)~=0
                       cloud(k,end)   = 0;               % don't keep
                   end
               end
               % consolidate I
               %cloudnew = cloud(cloud(:,end)==1,:);
               cflagind = find(cloudflag==1);
           
           if     2==1% test(cloud(end,1)==0)&& (size(cloud,1) >1 && (sum(cloud(:,end))==1 || sum(cloud(:,end))==0))
                  % choose next cloud from cloudflag
                  if strcmp(prof.direction,'descend')
                      cldstart = cflagind(end);
                      cldend   = cflagind(1);
                      cldbot   = prof.z(cldstart);%prof.z(cflagind(end));
                      cldtop   = prof.z(cldend);  %prof.z(cflagind(1));
                      cldthick = cldtop - cldbot;
                  else
                      cldstart = cflagind(1);
                      cldend   = cflagind(end);
                      cldbot   = prof.z(cldstart);     %prof.z(cflagind(1));
                      cldtop   = prof.z(cldend);       %prof.z(cflagind(end));
                      cldthick = cldtop - cldbot;
                  end
                   

                   % sum and avg normalized particle number per cloud
                   partsum   = zeros(1,length(reff));
                   for j=1:length(reff)
                       pstr = strcat('CDP',reff{:,j},'um');
                       partsum(j) = sum(prof.(pstr)(cldstart:cldend));
                   end
                   cldreff = round(sum((partsum.*reffum))/sum(partsum));
                   clear partsum; 
                   
                   % cloudphase
                    if cldreff <= 30 && nanmean(prof.Static_AirT(cldstart:c.cldend))> -20
                        cldphase = 0;
                    elseif cldreff > 30 && nanmean(prof.Static_AirT(cldstart:cldend))< -20
                        cldphase = 1;
                    elseif cldreff >= 30 && nanmean(prof.Static_AirT(cldstart:cldend))> -20 &&...
                                                 nanmean(prof.Static_AirT(cldstart:cldend))<  0
                        cldphase = 2;
                    else
                        cldphase = 0;
                    end
                    
                    % wc
                    if strcmp(prof.direction,'descend')
                        wp = nanmean(prof.LWP_mm(cldend-5:cldend));    % in mm (from GVR)
                    else
                        wp = nanmean(prof.LWP_mm(cldstart:cldstart+5));% in mm (from GVR)
                    end
                    if wp==0 wp=0.04; end;
                    cldwc = (wp*0.1*100^2)/cldthick;                   % this is water path [g/m3]
                    % construct consolidated cloud struct
                    clear c.cldnum c.cldstart c.cldend c.cldthick c.cldreff c.cldwc c.cldphase c.cldtop c.cldbot
                    if cloud(end,1)==0
                        c.cldstart = [cloud(end,4);cldstart];
                        c.cldend   = [cloud(end,5);cldend];
                        c.cldthick = [cloud(end,3);cldthick]; cldthick = c.cldthick;
                        c.cldbot   = [cloud(end,1);cldbot];   cldbot   = c.cldbot;
                        c.cldtop   = [cloud(end,2);cldtop];   cldtop   = c.cldtop;
                        c.cldphase = [cloud(end,6);cldphase]; cldphase = c.cldphase;
                        c.cldreff  = [cloud(end,7);cldreff];  cldreff  = c.cldreff;
                        c.cldwc    = [cloud(end,8);cldwc];    cldwc    = c.cldwc;
                        c.cldnum   = length(c.cldstart);
                    else
                        c.cldstart = [cldstart];
                        c.cldend   = [cldend];
                        c.cldthick = [cldthick];
                        c.cldbot   = [cldbot];
                        c.cldtop   = [cldtop];
                        c.cldphase = [cldphase];
                        c.cldreff  = [cldreff];
                        c.cldwc    = [cldwc];
                        c.cldnum   = length(c.cldstart);
                    end
           elseif size(cloud,1) > 1 && sum(cloud(:,end)) < (size(cloud,1))%size(cloud,1) > 1 && sum(cloud(:,end)) == (size(cloud,1)-1)
                  % consolidate II
                   for k=2:size(cloud,1)
                       dis(1) = 10000;
                       dis(k) = abs(cloud(k-1,5)-cloud(k,4));
                   end
                   cl = cloud(:,end);
                   [cmin I] = min(dis);
                   % assign new values
                   % consolidate cloud=0 with min cloud distance
                   if cl(I-1)==0
                           cloud(I,2) = cloud(I-1,2);% cldtop
                           cloud(I,4) = cloud(I-1,4);% cldstart
                           cloud(I,3) = cloud(I,2) - cloud(I,1);% cldthick

                           % reff
                           partsum   = zeros(1,length(reff));
                           for j=1:length(reff)
                               pstr = strcat('CDP',reff{:,j},'um');
                               partsum(j) = sum(prof.(pstr)(cloud(I,4):cloud(I,5)));
                           end
                           cloud(I,7) = round(sum((partsum.*reffum))/sum(partsum));
                           clear partsum; 

                           % cloudphase
                            if cloud(I,7) <= 30 && nanmean(prof.Static_AirT(cloud(I,4):cloud(I,5)))> -20
                                cloud(I,6) = 0;
                            elseif cloud(I,7) > 30 && nanmean(prof.Static_AirT(cloud(I,4):cloud(I,5)))< -20
                                cloud(I,6) = 1;
                            elseif cloud(I,7) >= 30 && nanmean(prof.Static_AirT(cloud(I,4):cloud(I,5)))> -20 &&...
                                                         nanmean(prof.Static_AirT(cloud(I,4):cloud(I,5)))<  0
                                cloud(I,6) = 2;
                            else
                                cloud(I,6) = 0;
                            end

                            % wc
                            if strcmp(prof.direction,'descend')
                                wp = nanmean(prof.LWP_mm(cloud(I,5)-5:cloud(I,5)));    % in mm (from GVR)
                            else
                                wp = nanmean(prof.LWP_mm(cloud(I,4):cloud(I,4)+5));    % in mm (from GVR)
                            end
                            if wp==0 wp=0.04; end;
                            cloud(I,8) = (wp*0.1*100^2)/cloud(I,3);                    % this is water path [g/m3]

                            % construct consolidated cloud struct
                            clear c.cldnum c.cldstart c.cldend c.cldthick c.cldreff c.cldwc c.cldphase c.cldtop c.cldbot
                            if sum(cl) == size(cloud,1) -1
                                cloud = cloud(cl==1,:);
                                c.cldnum   = sum(cl==1);
                            else
                                cloudtmp = cloud(1:I-2,:);cloud = [cloudtmp;cloud(I:end,:)];
                                c.cldnum = size(cloud,1);
                            end
                             c.cldstart = cloud(:,4);
                             c.cldend   = cloud(:,5);
                             c.cldthick = cloud(:,3); cldthick = c.cldthick;
                             c.cldbot   = cloud(:,1); cldbot   = c.cldbot;
                             c.cldtop   = cloud(:,2); cldtop   = c.cldtop;
                             c.cldphase = cloud(:,6); cldphase = c.cldphase;
                             c.cldreff  = cloud(:,7); cldreff  = c.cldreff;
                             c.cldwc    = cloud(:,8); cldwc    = c.cldwc;
                       else

                             % use original c.cld struct
                            cldbot     = cloud(:,1);
                            cldtop     = cloud(:,2);
                            cldthick   = cloud(:,3);
                            cldreff    = cloud(:,7);
                            cldphase   = cloud(:,6);
                            cldwc      = cloud(:,8);
                            cldstart   = cloud(:,4);
                            cldend     = cloud(:,5);
                   end
           else
                    % use original c.cld struct
                    cldbot     = cloud(:,1);
                    cldtop     = cloud(:,2);
                    cldthick   = cloud(:,3);
                    cldreff    = cloud(:,7);
                    cldphase   = cloud(:,6);
                    cldwc      = cloud(:,8);
                    cldstart   = cloud(:,4);
                    cldend     = cloud(:,5);
                   
           end
           end
           
          
           %% save cloud params to file
           % write data to file
           filen=['C:\cygwin\home\msegalro\libradtran\input\CRF\arise\MERRABASE\MERRABASE_consolidateclouds.txt'];
           disp( ['Writing to file: ' filen]);
           [filepath filename ext] = fileparts(prof.ana.atmfile);
           atmosfile = strcat(filename, ext);
           %for i=1:c.cldnum
               if prof.iceconc >15
                   s_albedo_file = 'albedo_ice.dat';%strcat('F:\ARISE\ArcticCRF\METdata\albedo\albedo_ice.dat');
                   w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                   line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                          ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                          ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                          ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff')  ' cloudphase ' num2str(cldphase')...
                          ' cloudwc '   num2str(cldwc') '  '}];
               else
                   s_albedo_file = 'albedo_water.dat';
                   w_albedo_file = strcat('albedo_ice_' ,num2str(round(prof.iceconc)),'.dat');
                   line=[{'date ' daystr ' profilenum ' num2str(dprof) ' DOY ' doy ' sza ' num2str(sza) ' lat ' num2str(lat) ' lon ' num2str(lon) ...
                          ' single albedo file ' s_albedo_file ' ice_conc ' num2str(round(prof.iceconc)) ' weighted albedo file ' w_albedo_file ' atmos file ' atmosfile...
                          ' platform_low_alt ' num2str(c.c130sur_alt) ' numclouds ' num2str(c.cldnum) ' cloudabove ' num2str(c.cldabove) ' cloudbelow ' num2str(c.cldbelow) ' cloudtop ' num2str(cldtop'/1000)...
                          ' cloudbot '  num2str(cldbot'/1000) ' cloudthick ' num2str(cldthick'/1000) ' cloudreff ' num2str(cldreff') ' cloudphase ' num2str(cldphase')...
                          ' cloudwc '   num2str(cldwc') '  '}];
               end
               
               dlmwrite(filen,line,'-append','delimiter','');

               clear line;
           %end

return;

