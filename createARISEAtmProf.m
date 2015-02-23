%% Details of the function:
%  NAME:
% createARISEAtmProf
%----------------------------------
% PURPOSE:
%  - create profiles from ARISE flights-vertical legs
%  - iclude Temp/Rh and in-situ/cloud data
%  - plot individual profiles and average them
%  - these profiles are used for compare with MERRA
%  - for creating atm profile for RT, and for analysis
%  - of cloud profile effect
%
% CALLING SEQUENCE:
%  out = createARISEAtmProf
%
% INPUT:
%  - arise-C130-Hskping_c130_20140830_RA_Preliminary.ict air data
% 
% 
% OUTPUT:
%  saved out struct and plots; daily averages and std's
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%  - bin_profiles.m
%
% NEEDED FILES/INPUT:
% 
% - arise-C130-Hskping_c130_20140830_RA_Preliminary.ict air data
%
% EXAMPLE:
%  - out = createARISEAtmProf
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Jan-27-2015
% MS: 2015-02-23, added version set tracking (v1.0)
%                 created repo, test
%                 added probe .ice file read
% -------------------------------------------------------------------------
%% function routine
function out = createARISEAtmProf
version_set('1.0');
startup_plotting;
%compare = 'false'; % defualt is false
%% choose input directory or file
% ------------------------------ %
infile = getfullname__('arise-C130-Hskping_c130_*.ict','F:','Select a met ict* files directory');
[dirn, fname, ext] = fileparts(infile);
 dirn = [dirn, filesep];

allfiles = dir([dirn,'arise-C130-Hskping_c130_*.ict']);
% read all data
    for i=1:length(allfiles)
        dm{i} = strsplit(allfiles(i).name,'_');
        dates(i) = dm{:,i}(2);
        tmp=ictread([dirn allfiles(i).name]);
        air.(strcat('met',dates{:,i},'_')).UTCsec        = tmp.Start_UTC;
        air.(strcat('met',dates{:,i},'_')).UTChr         =(air.(strcat('met',dates{:,i},'_')).UTCsec/86400)*24;
        air.(strcat('met',dates{:,i},'_')).longitude     = tmp.Longitude;
        air.(strcat('met',dates{:,i},'_')).latitude      = tmp.Latitude;
        air.(strcat('met',dates{:,i},'_')).GPS_Altitude  = tmp.GPS_Altitude;%[m]
        air.(strcat('met',dates{:,i},'_')).Pst_Altitude  = tmp.Pressure_Altitude;%[ft]
        air.(strcat('met',dates{:,i},'_')).Static_AirT   = tmp.Static_Air_Temp;%[C]
        air.(strcat('met',dates{:,i},'_')).poten_Temp    = tmp.Potential_Temp;%[K]
        air.(strcat('met',dates{:,i},'_')).DewPoint      = tmp.Dew_Point_3Stage;%[C]
        air.(strcat('met',dates{:,i},'_')).TotalT        = tmp.Total_Air_Temp;%[C]
        air.(strcat('met',dates{:,i},'_')).IRsurfT       = tmp.IR_Surf_Temp;%[C]
        air.(strcat('met',dates{:,i},'_')).StaticP       = tmp.Static_Pressure;%[mbar]
        air.(strcat('met',dates{:,i},'_')).WindS         = tmp.Wind_Speed;%[m/s]
        air.(strcat('met',dates{:,i},'_')).WindD         = tmp.Wind_Direction;%deg[0-360]
        air.(strcat('met',dates{:,i},'_')).SZA           = tmp.Solar_Zenith_Angle;
        air.(strcat('met',dates{:,i},'_')).MMR           = tmp.Mixing_Ratio;%[g/kg]
        air.(strcat('met',dates{:,i},'_')).PartPwater    = tmp.Part_Press_Water_Vapor;%[mbar]
        air.(strcat('met',dates{:,i},'_')).SatVPwater    = tmp.Sat_Vapor_Press_H2O;%[mbar]
        air.(strcat('met',dates{:,i},'_')).SatVPice      = tmp.Sat_Vapor_Press_Ice;%[mbar]
        air.(strcat('met',dates{:,i},'_')).RH            = tmp.Relative_Humidity;  %[%] 
        clear tmp;
    end
    
%% create vertical profiles
Tempmov(length(dates)) = struct('cdata',[],'colormap',[]);
RHmov(length(dates))   = struct('cdata',[],'colormap',[]);
for i=1:length(allfiles)
    
    % plot flight altitude
    figure(998);plot(air.(strcat('met',dates{:,i},'_')).latitude,...
                air.(strcat('met',dates{:,i},'_')).GPS_Altitude,'-b');
           xlabel('latitude');ylabel('Altitude [m]');title(dates{:,i});
    figure(999);plot(air.(strcat('met',dates{:,i},'_')).UTChr,...
                air.(strcat('met',dates{:,i},'_')).latitude,'-b');
           xlabel('time [UThr]');ylabel('latitude');title(dates{:,i});
    figure(1000);plot(air.(strcat('met',dates{:,i},'_')).UTChr,...
                air.(strcat('met',dates{:,i},'_')).GPS_Altitude,'-b');
           xlabel('time [UThr]');ylabel('Altitude [m]');title(dates{:,i});
    % form separate profiles
    profnum = input('input number of profiles');
    for j=1:profnum
        disp('choose bottom to top points')
        prof = ginput(2);
        if prof(2,1)-prof(1,1)>0
            air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).direction = 'ascend';
            ut1 = prof(1,1);
            ut2 = prof(2,1);
        else
            air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).direction = 'descend';
            ut1 = prof(2,1);
            ut2 = prof(1,1);
        end
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).lon    =...
             air.(strcat('met',dates{:,i},'_')).longitude(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).lat    =...
             air.(strcat('met',dates{:,i},'_')).latitude(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).z    =...
             air.(strcat('met',dates{:,i},'_')).GPS_Altitude(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).time =...
             air.(strcat('met',dates{:,i},'_')).UTChr(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).Static_AirT =...
             air.(strcat('met',dates{:,i},'_')).Static_AirT(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).theta =...
             air.(strcat('met',dates{:,i},'_')).poten_Temp(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).DewPoint =...
             air.(strcat('met',dates{:,i},'_')).DewPoint(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).staticP =...
             air.(strcat('met',dates{:,i},'_')).StaticP(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).MMR =...
             air.(strcat('met',dates{:,i},'_')).MMR(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).RH =...
             air.(strcat('met',dates{:,i},'_')).RH(air.(strcat('met',dates{:,i},'_')).UTChr<=ut2 &...
                                                             air.(strcat('met',dates{:,i},'_')).UTChr>=ut1);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).meanlon = ...
                             nanmean(air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).lon);
        air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).meanlat = ...
                             nanmean(air.(strcat('met',dates{:,i},'_')).(strcat('profnum',num2str(j))).lat);
    
    
    % compare with reanalysis
    
%         if strcmp(compare,'true')
%             reana = readMERRA(0,'false');
%             % choose closest profile
%             %nm_340 = interp1(vis.nm,[1:length(vis.nm)],340.0, 'nearest');
%             air.(strcat('met',dates{:,i},'_')).(strcat('reanum',num2str(j))) = reana;
%         end
    
    end
    % average profiles
    [binnedprof] = bin_profiles(air.(strcat('met',dates{:,i},'_')),profnum);
    air.(strcat('met',dates{:,i},'_')).prof = binnedprof;
    % plot avg and std profiles for each day
    fh = figure(i);
    % pressure
    ax(1) = subplot(2,2,1);
    plot(air.(strcat('met',dates{:,i},'_')).prof.Pmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-b',...
         air.(strcat('met',dates{:,i},'_')).prof.Pmean + air.(strcat('met',dates{:,i},'_')).prof.Pstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--b',...
         air.(strcat('met',dates{:,i},'_')).prof.Pmean - air.(strcat('met',dates{:,i},'_')).prof.Pstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--b','linewidth',2);
    xlabel('Pressure [mbar]');
    ylabel('Altitude [km]');
    text(410,1,dates{:,i});
    axis([400 1000 0 7]);
    legend('mean','std');grid on;
    % temperature
    ax(2) = subplot(2,2,2);
    plot(air.(strcat('met',dates{:,i},'_')).prof.Tmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-r',...
         air.(strcat('met',dates{:,i},'_')).prof.Tmean + air.(strcat('met',dates{:,i},'_')).prof.Tstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--r',...
         air.(strcat('met',dates{:,i},'_')).prof.Tmean - air.(strcat('met',dates{:,i},'_')).prof.Tstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--r','linewidth',2);
    xlabel('Temperature [C]');
    set(ax(2),'YTickLabel',{''})
    axis([-40 5 0 7]);
    legend('mean','std');grid on;
    % mixing ratio
    ax(3) = subplot(2,2,3);
    plot(air.(strcat('met',dates{:,i},'_')).prof.MMRmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-c',...
         air.(strcat('met',dates{:,i},'_')).prof.MMRmean + air.(strcat('met',dates{:,i},'_')).prof.MMRstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--c',...
         air.(strcat('met',dates{:,i},'_')).prof.MMRmean - air.(strcat('met',dates{:,i},'_')).prof.MMRstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--c','linewidth',2);
    xlabel('water vapor mixing ratio [g/kg]');
    ylabel('Altitude [km]');
    axis([0 0.4 0 7]);
    legend('mean','std');grid on;
    % relative humidity
    ax(4) = subplot(2,2,4);
    plot(air.(strcat('met',dates{:,i},'_')).prof.RHmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-g',...
         air.(strcat('met',dates{:,i},'_')).prof.RHmean + air.(strcat('met',dates{:,i},'_')).prof.RHstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--g',...
         air.(strcat('met',dates{:,i},'_')).prof.RHmean - air.(strcat('met',dates{:,i},'_')).prof.RHstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--g','linewidth',2);
    xlabel('relative humidity [%]');
    set(ax(4),'YTickLabel',{''});
    set(ax(4),'XTick',[0:20:100]);
    set(ax(4),'XTickLabel',{'0','20','40','60','80','100'});
    axis([0 100 0 7]);
    legend('mean','std');grid on;
    linkaxes(ax,'y');
    % get subplot locations:
    % [left, bottom, width, height].
     p1 = get(ax(1), 'position');
     p2 = get(ax(2), 'position');
     p3 = get(ax(3), 'position');
     p4 = get(ax(4), 'position');
     % adjust locations
     p2(1) = p1(1) + p1(3) +p1(3)/8;
     p4(1) = p3(1) + p3(3) +p3(3)/8;
     set(ax(2), 'position', p2)
     set(ax(4), 'position', p4);
     fi=[strcat(dirn, dates{:,i}, 'dailymeanProfiles4subplot')];
     save_fig(i,fi,false);
     close 1000;
     close(fh);
   % create temperature profile movie
                    h1 = figure;set(h1, 'Color','white');
                    bx=plot(air.(strcat('met',dates{:,i},'_')).prof.Tmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-k',...
                    air.(strcat('met',dates{:,i},'_')).prof.Tmean + air.(strcat('met',dates{:,i},'_')).prof.Tstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--k',...
                    air.(strcat('met',dates{:,i},'_')).prof.Tmean - air.(strcat('met',dates{:,i},'_')).prof.Tstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--k','linewidth',2);
                    axis tight manual
                    %bx = gca;
                    bx.NextPlot = 'replaceChildren';
                    xlabel('Temperature [C]','fontSize',12);
                    ylabel('Altitude [km]','FontSize',12);
                    set(gca,'fontsize',12);
                    title(dates{:,i});
                    legend('mean','std');
                    axis([-20 30 0 7]);
                    date = datenum(str2double(dates{:,i}(1:4)), str2double(dates{:,i}(5:6)), str2double(dates{:,i}(7:8)));   % Convert date into serial numbers
                    str = datestr(date, 'dd mmm yyyy'); % Show in the format
                    drawnow
                    Tempmov(i) = getframe(h1);
    % create RH profile movie
                    h2 = figure;set(h2, 'Color','white');
                    cx=plot(air.(strcat('met',dates{:,i},'_')).prof.RHmean,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'-k',...
                    air.(strcat('met',dates{:,i},'_')).prof.RHmean + air.(strcat('met',dates{:,i},'_')).prof.RHstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--k',...
                    air.(strcat('met',dates{:,i},'_')).prof.RHmean - air.(strcat('met',dates{:,i},'_')).prof.RHstd,air.(strcat('met',dates{:,i},'_')).prof.zmean/1000,'--k','linewidth',2);
                    axis tight manual
                    %cx = gca;
                    cx.NextPlot = 'replaceChildren';    
                    xlabel('Relative Humidity [%]','fontSize',12);
                    ylabel('Altitude [km]','FontSize',12);
                    set(gca,'fontsize',12);
                    title(dates{:,i});
                    legend('mean','std');
                    axis([0 100 0 7]);
                    date = datenum(str2double(dates{:,i}(1:4)), str2double(dates{:,i}(5:6)), str2double(dates{:,i}(7:8)));   % Convert date into serial numbers
                    str = datestr(date, 'dd mmm yyyy'); % Show in the format
                    drawnow
                    RHmov(i) = getframe(h2);
    % make binned profiles for each day (all flight)
    
    
    
end

%% save movies
movie2avi(Tempmov, 'ARISEmeanTempprof.avi', 'compression','None', 'fps',1);
winopen('ARISEmeanTempprof.avi');
movie2avi(RHmov, 'ARISEmeanRHprof.avi', 'compression','None', 'fps',1);
winopen('ARISEmeanRHprof.avi');

%% save processed air struct for all flights
disp(strcat('saving',dirn,'ARISEairprocessed.mat'));
save([dirn,'ARISEairprocessed.mat'],'-struct','air');
return;



