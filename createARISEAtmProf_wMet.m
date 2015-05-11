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
%  out = createARISEAtmProf_wMet
%
% INPUT:
%  - arise-C130-Hskping_c130_20140830_RA_Preliminary.ict air data
%  - Arise-Waterparams-13Sept2014.csv complimentary water vapor air data
%  - 20140901Probes.csv probes .csv files
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
% - Arise-Waterparams-13Sept2014.csv complimentary water vapor air data
% - 20140901Probes.csv probes .csv files
%
%
% EXAMPLE:
%  - out = createARISEAtmProf_wMET
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Jan-27-2015
% MS: 2015-02-23, added version set tracking (v1.0)
%                 created repo, test
%                 added probe .ice file read
%                 test
%                 change version set to v1.1 due to probe data addition
% -------------------------------------------------------------------------
%% function routine
function out = createARISEAtmProf_wMet
version_set('1.1');
startup_plotting;
%compare = 'false'; % defualt is false
%% choose input directory or file
% ------------------------------ %
% air data
% under: F:\ARISE\C-130data\Met
airinfile = getfullname__('arise-C130-Hskping_c130_*.ict','F:','Select a met ict* files directory');
[airdir, fname, ext] = fileparts(airinfile);
 airdir = [airdir, filesep];
 air_files = dir([airdir,'arise-C130-Hskping_c130_*.ict']);
 
 % complimentary water vapor air data in .csv format
 % under F:\ARISE\C-130data\Met\NewArchiveData
 wvinfile = getfullname__('Arise_Waterparams*.csv','F:','Select water vapor .csv files directory');
[wvdir, fname, ext] = fileparts(wvinfile);
 wvdir = [wvdir, filesep];
 wv_files = dir([wvdir,'Arise-Waterparams*.csv']);
 
 % in-situ probe data
 % under: F:\ARISE\C-130data\Probes
 isinfile = getfullname__('*Probes.csv','F:','Select probes *.csv files directory');
[isdir, fname, ext] = fileparts(isinfile);
 isdir = [isdir, filesep];
 is_files = dir([isdir,'*Probes.csv']);
 
 % 4star data
 % under F:\ARISE\starmat_ARISE
 starinfile = getfullname__('*star.mat','F:','Select star.mat files directory');
[stardir, fname, ext] = fileparts(starinfile);
 stardir = [stardir, filesep];
 star_files = dir([stardir,'*star.mat']);

% read all data
    for i=1:length(air_files)
        dm{i} = strsplit(air_files(i).name,'_');
        dates(i) = dm{:,i}(2);
        dates_str = strcat('met',dates{:,i},'_');
        tmp =ictread([airdir air_files(i).name]);        % air data    
        air.(dates_str).UTCsec        = tmp.Start_UTC;
        air.(dates_str).UTChr         =(air.(dates_str).UTCsec/86400)*24;
        air.(dates_str).longitude     = tmp.Longitude;
        air.(dates_str).latitude      = tmp.Latitude;
        air.(dates_str).GPS_Altitude  = tmp.GPS_Altitude;%[m]
        air.(dates_str).Pst_Altitude  = tmp.Pressure_Altitude;%[ft]
        air.(dates_str).Static_AirT   = tmp.Static_Air_Temp;%[C]
        air.(dates_str).poten_Temp    = tmp.Potential_Temp;%[K]
        air.(dates_str).DewPoint      = tmp.Dew_Point_3Stage;%[C]
        air.(dates_str).TotalT        = tmp.Total_Air_Temp;%[C]
        air.(dates_str).IRsurfT       = tmp.IR_Surf_Temp;%[C]
        air.(dates_str).StaticP       = tmp.Static_Pressure;%[mbar]
        air.(dates_str).WindS         = tmp.Wind_Speed;%[m/s]
        air.(dates_str).WindD         = tmp.Wind_Direction;%deg[0-360]
        air.(dates_str).SZA           = tmp.Solar_Zenith_Angle;
        air.(dates_str).MMR           = tmp.Mixing_Ratio;%[g/kg]
        air.(dates_str).PartPwater    = tmp.Part_Press_Water_Vapor;%[mbar]
        air.(dates_str).SatVPwater    = tmp.Sat_Vapor_Press_H2O;%[mbar]
        air.(dates_str).SatVPice      = tmp.Sat_Vapor_Press_Ice;%[mbar]
        air.(dates_str).RH            = tmp.Relative_Humidity;  %[%] 
        clear tmp;
        
        % replace ict met data with existing .csv water vapor data
        wvfname = strcat(wvdir,'Arise_Waterparams_',dates(i),'.csv');
        if exist(wvfname{:},'file')
            
            tmp =importdata(wvfname{:});   % water vapor data 
            timeutc = (tmp.data(:,1)/86400)*24;
            pressure= tmp.data(:,5);
            sat     = tmp.data(:,6);
            dewpoint= tmp.data(:,7);
            partpres= tmp.data(:,8);
            mmr     = tmp.data(:,9);
            satPwat = tmp.data(:,10);
            satPice = tmp.data(:,11);
            rh      = tmp.data(:,12);
            
            if strcmp(dates(i),'20140924')
                timeutc(2) = timeutc(1) + 0.0001;
            end
            
            air.(dates_str).StaticP      = interp1(timeutc,pressure, air.(dates_str).UTChr);
            air.(dates_str).SatVPwater   = interp1(timeutc,satPwat , air.(dates_str).UTChr);
            air.(dates_str).SatVPice     = interp1(timeutc,satPice , air.(dates_str).UTChr);
            air.(dates_str).DewPoint     = interp1(timeutc,dewpoint, air.(dates_str).UTChr);
            air.(dates_str).MMR          = interp1(timeutc,mmr     , air.(dates_str).UTChr);
            air.(dates_str).PartPwater   = interp1(timeutc,partpres, air.(dates_str).UTChr);
            air.(dates_str).RH           = interp1(timeutc,rh      , air.(dates_str).UTChr);
            
            clear tmp;
        end
        
        tmp=csvread([isdir  is_files(i).name]);          % probes data
        % interpolate params to air time
        is_params = [5:9 14 18 23 30 35 40 45];%{'TWC_gm3','LWC1_gm3','LWC2_gm3','PWV_cm','LWP_mm','nCDP_cm3','CDP03_dNdlogD','CDP08_dNdlogD','CDP15_dNdlogD','CDP20_dNdlogD','CDP25_dNdlogD','CDP30_dNdlogD'};
        is_labels = {'TWC_gm3','LWC1_gm3','LWC2_gm3','PWV_cm','LWP_mm','nCDP_cm3','CDP05um',...
                     'CDP10um','CDP20um','CDP30um','CDP40um','CDP50um'};
        %CDP03_dNdlogD% cm-3; 5um
        %CDP08_dNdlogD% cm-3; 10um
        %CDP15_dNdlogD% cm-3; 20 um
        %CDP20_dNdlogD% cm-3; 30um
        %CDP25_dNdlogD% cm-3; 40um
        %CDP30_dNdlogD% cm-3; 50um
        for kk=1:length(is_labels)
                air.(dates_str).(is_labels{:,kk}) = ...
                     interp1(24*(tmp(:,1)/86400),tmp(:,is_params(kk)),...
                     air.(dates_str).UTChr,'nearest');
        end
        clear tmp;
        % load 4star
        s=load([stardir  star_files(i).name]);           % 4star data
        utc_range = [air.(dates_str).UTChr(1),...
                     air.(dates_str).UTChr(end)];
        star=combine_star_params(s,utc_range);
        clear s;
        star_params = {'Str','Md'}; star_labels = {'starStr','starMd'};
        if length(star.utc)~=length(unique(star.utc))
            star.utc = star.utc + (([1:length(star.utc)])*1e-12)';
        end
        for kk=1:length(star_params)
                air.(dates_str).(star_labels{:,kk}) = ...
                     interp1(star.utc,star.(star_params{:,kk}),...
                     air.(dates_str).UTChr,'nearest');
        end
    end
    
%% create vertical profiles
Tempmov(length(dates)) = struct('cdata',[],'colormap',[]);
RHmov(length(dates))   = struct('cdata',[],'colormap',[]);
for i=1:length(air_files)
    dates_str = strcat('met',dates{:,i},'_');
    % plot flight altitude
    figure(998);plot(air.(dates_str).latitude,...
                air.(dates_str).GPS_Altitude,'-b');
           xlabel('latitude');ylabel('Altitude [m]');title(dates{:,i});
    figure(999);plot(air.(dates_str).UTChr,...
                air.(dates_str).latitude,'-b');
           xlabel('time [UThr]');ylabel('latitude');title(dates{:,i});
    figure(1000);plot(air.(dates_str).UTChr,...
                air.(dates_str).GPS_Altitude,'-b');
           xlabel('time [UThr]');ylabel('Altitude [m]');title(dates{:,i});
    % form separate profiles
    profnum = input('input number of profiles');
    for j=1:profnum
                prof_str = strcat('profnum',num2str(j));
                disp('choose bottom to top points')
                prof = ginput(2);
                if prof(2,1)-prof(1,1)>0
                    air.(dates_str).(prof_str).direction = 'ascend';
                    ut1 = prof(1,1);
                    ut2 = prof(2,1);
                else
                    air.(dates_str).(prof_str).direction = 'descend';
                    ut1 = prof(2,1);
                    ut2 = prof(1,1);
                end
                air.(dates_str).(prof_str).lon    =...
                     air.(dates_str).longitude(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).lat    =...
                     air.(dates_str).latitude(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).z    =...
                     air.(dates_str).GPS_Altitude(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).time =...
                     air.(dates_str).UTChr(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).Static_AirT =...
                     air.(dates_str).Static_AirT(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).theta =...
                     air.(dates_str).poten_Temp(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).DewPoint =...
                     air.(dates_str).DewPoint(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).staticP =...
                     air.(dates_str).StaticP(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).MMR =...
                     air.(dates_str).MMR(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).RH =...
                     air.(dates_str).RH(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).SatVPwater =...
                     air.(dates_str).SatVPwater(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).SatVPice =...
                     air.(dates_str).SatVPice(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).starStr =...
                     air.(dates_str).starStr(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                air.(dates_str).(prof_str).starMd =...
                     air.(dates_str).starMd(air.(dates_str).UTChr<=ut2 &...
                                                                     air.(dates_str).UTChr>=ut1);
                                                                 
                % assimilate in-situ data into profile
                for kk=1:length(is_labels)
                        air.(dates_str).(prof_str).(is_labels{:,kk}) = ...
                        air.(dates_str).(is_labels{:,kk})(air.(dates_str).UTChr<=ut2 &...
                                                          air.(dates_str).UTChr>=ut1);
                end
                                                                 
                % mean profile location
                air.(dates_str).(prof_str).meanlon = ...
                nanmean(air.(dates_str).(prof_str).lon);
                air.(dates_str).(strcat('profnum',num2str(j))).meanlat = ...
                nanmean(air.(dates_str).(prof_str).lat);
    
    end
    % average profiles
    [binnedprof] = bin_profiles(air.(dates_str),profnum);
    air.(dates_str).prof = binnedprof;
    % plot avg and std profiles for each day
    fh = figure(i);
    % pressure
    ax(1) = subplot(2,2,1);
    plot(air.(dates_str).prof.Pmean,air.(dates_str).prof.zmean/1000,'-b',...
         air.(dates_str).prof.Pmean + air.(dates_str).prof.Pstd,air.(dates_str).prof.zmean/1000,'--b',...
         air.(dates_str).prof.Pmean - air.(dates_str).prof.Pstd,air.(dates_str).prof.zmean/1000,'--b','linewidth',2);
    xlabel('Pressure [mbar]');
    ylabel('Altitude [km]');
    text(410,1,dates{:,i});
    axis([400 1000 0 7]);
    legend('mean','std');grid on;
    % temperature
    ax(2) = subplot(2,2,2);
    plot(air.(dates_str).prof.Tmean,air.(dates_str).prof.zmean/1000,'-r',...
         air.(dates_str).prof.Tmean + air.(dates_str).prof.Tstd,air.(dates_str).prof.zmean/1000,'--r',...
         air.(dates_str).prof.Tmean - air.(dates_str).prof.Tstd,air.(dates_str).prof.zmean/1000,'--r','linewidth',2);
    xlabel('Temperature [C]');
    set(ax(2),'YTickLabel',{''})
    axis([-40 5 0 7]);
    legend('mean','std');grid on;
    % mixing ratio
    ax(3) = subplot(2,2,3);
    plot(air.(dates_str).prof.MMRmean,air.(dates_str).prof.zmean/1000,'-c',...
         air.(dates_str).prof.MMRmean + air.(dates_str).prof.MMRstd,air.(dates_str).prof.zmean/1000,'--c',...
         air.(dates_str).prof.MMRmean - air.(dates_str).prof.MMRstd,air.(dates_str).prof.zmean/1000,'--c','linewidth',2);
    xlabel('water vapor mixing ratio [g/kg]');
    ylabel('Altitude [km]');
    axis([0 0.4 0 7]);
    legend('mean','std');grid on;
    % relative humidity
    ax(4) = subplot(2,2,4);
    plot(air.(dates_str).prof.RHmean,air.(dates_str).prof.zmean/1000,'-g',...
         air.(dates_str).prof.RHmean + air.(dates_str).prof.RHstd,air.(dates_str).prof.zmean/1000,'--g',...
         air.(dates_str).prof.RHmean - air.(dates_str).prof.RHstd,air.(dates_str).prof.zmean/1000,'--g','linewidth',2);
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
     fi=[strcat(airdir, dates{:,i}, 'dailymeanProfiles4subplot')];
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
disp(strcat('saving',airdir,'ARISEairprocessed_with_insitu_withWVparams.mat'));
save([airdir,'ARISEairprocessed_with_insitu_withWVparams.mat'],'-struct','air');
return;



