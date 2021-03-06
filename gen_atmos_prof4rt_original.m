%% PURPOSE:
%   Generate the atmosphere profile from a standard atmospheric, reanalysis
%   profile and the temperatures, pressures, and relative humidity measured
%   onboard the NASA c130
%   this is for specific profiles defined by star_ana_compare
%
% CALLING SEQUENCE:
%   filen = gen_atmos_prof4rt(air,reana,mode,daystr,profstr)
%
% INPUT:
%   air - structure of atmospheric parameter profiles from airborne
%         platform (e.g. c130)
%   reana - struct of atmospheric profile parameters from reanalysis (e.g.
%           MERRA)
%   mode  - 0 combine std atm and reanalysis only, str=MERRA/MERRA2
%           0.1 same as 0, just for FP
%           1 combine, std atm, reanalysis and airborne, str=C130
%           2 same as 1, str = C130interp
%  daystr is date yyyymmdd
%  profstr is the profile number for each day (e.g. 'profnum1')
% 
% OUTPUT:
%  - plots of profiles of temperature, pressure, and relative humidity
%  - file of atmosheric profile for libRadTran simulations
%  - filein is name of atmospheric file generated from  given profile
%    yyyymmdd_profnumx_MERRA (mode 0), yyyymmdd_profnumx_C130 (mode 1)
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%  - t2utch.m
%  - magnus.m ; for ideal gas law 
%  - version_set.m
%
% NEEDED FILES:
%  - standard atmosphere file ('afglss.dat')
%
% EXAMPLE:
%  
%
% MODIFICATION HISTORY:
% Written: Michal Segal, NASA Ames, 2015-02-25
% MS, 2015-03-13, fixed bug in interpolation to 0 altitude in mode 0
% MS, based on gen_atmos_rt. this version is for latest FP reanalysis data
% MS, 2015-10-12, added mode==2 to diff between native and interpolated
%                 C130 fields, added mode=0.5, for FP product in filename
% -------------------------------------------------------------------------

%% Start of function
function [filen]=gen_atmos_prof4rt(air,ana,mode,daystr,profstr)
startup_plotting;
version_set('1.0');
plotting = false;
%% handle inputs

% airborne
if mode==1
    % original airborne data
    air.t = air.Static_AirT+273.15;
    air.p = air.staticP;
    air.rh= air.RH;
    air.z = air.z;
    air.h2ond=magnus(100*air.rh,air.t);% rh needs to be in %
    air.airnd=TP2numd(air.t,air.p,'air');
elseif mode==2
    % interpolated airborne data
    air.t = air.airTinterp;
    air.p = ana.plev;
    air.rh= air.airRHinterp;
    air.z = ana.zalt_mean;
    air.h2ond=magnus(100*air.rh,air.t);% rh needs to be in %
    air.airnd=TP2numd(air.t,air.p,'air');
end

% standard
dir = 'F:\ARISE\ArcticCRF\METdata\atm_profiles\';
std=importdata([dir 'afglss.dat']);
std.z=std.data(:,1).*1000.;
std.p=std.data(:,2);
std.t=std.data(:,3);%-273.15;
std.h2ond=std.data(:,7);
std.airnd=std.data(:,4);
std.o3nd=std.data(:,5);
std.o2nd=std.data(:,6);
std.co2nd=std.data(:,8);
std.no2nd=std.data(:,9);

% reanalysis/model
ana.z = flipud(ana.zalt_mean);
ana.p = flipud(ana.plev);          
ana.t = flipud(ana.T_mean);         
ana.rh= flipud(ana.RH_mean);
ana.h2ond=flipud(ana.h2o);
ana.o3nd = flipud(ana.o3);
ana.airnd= flipud(ana.air);

%% prepare for plotting
if plotting

figure(1000);

% Pressure
subplot(2,2,1);
plot(air.p,air.z,'ob',ana.p,ana.z,'r.',std.p,std.z,'go');
xlabel('Pressure [mb]');
ylabel('Altitude [m]');
title('Pressure');
legend('C130','reanalysis','Standard');

%Temperature
subplot(2,2,2);
plot(air.t,air.z,'ob',ana.t,ana.z,'r.',std.t,std.z,'go');
xlabel('Temperature [�C]');
ylabel('Altitude [m]');
title('Temperature');
legend('C130','reanalysis','Standard');

%Relative Humidity
subplot(2,2,3);
% plot(air.rh,air.z,'b.',ana.rh,ana.z,'r.');%,std.RH,std.z,'go');
% xlabel('Relative Humidity [%]');
% ylabel('Altitude [m]');
% title('Relative Humidity');
% legend('C130','reanalysis');
%air density
plot(air.airnd,air.z,'ob',ana.airnd,ana.z,'r.',std.airnd,std.z,'go');
xlabel('air number density [#/cm^{-3}]');
ylabel('Altitude [m]');
title('Relative Humidity');
legend('C130','reanalysis','standard');
title('Air density');

%water vapor mixture
subplot(2,2,4);
plot(air.h2ond,air.z,'ob',ana.h2ond,ana.z,'r.',std.h2ond,std.z,'go');
xlabel('Water vapor number density [#/cm^{-3}]');
ylabel('Altitude [m]');
title('Water vapor');
legend('C130','reanalysis','Standard');

save_fig(1000,[dir,'atm_profiles_' daystr '_' profstr],true);
 
end

%% Select the correct values, first from airborne, then from reanalysis, then from standard

mode_ = mode;
mode = round(mode_);

if mode==0
    % assign altitude according to reanalysis
    % only up to 30km - for Tinterp
    % only up to 15km - for T and P interp
    intind = interp1(ana.z,[1:length(ana.z)],15000,'nearest');
    z1=ana.z(intind:end);%ana.z(1):-100:0;%1800:-200:1000;
    zind = interp1(std.z,[1:length(std.z)],max(z1),'nearest');
    if z1(end)~=0 z1 = [z1; 0.0]; end
    Z=[std.z(1:zind-1);z1];
    % interp reanalysis
    [ana.zz,ana.i]=sort(ana.z);
    [ana.zz,ana.ia]=unique(ana.zz);
    ana.i=ana.i(ana.ia);

    ana.Tz=interp1(ana.z,ana.t,Z,'linear','extrap');
    ana.Pz=interp1([ana.z;120000],[ana.p; 2e-10],Z,'linear','extrap');
    ana.Hz=interp1(ana.z,ana.h2ond,Z,'linear','extrap');
    ana.Az=interp1(ana.z,ana.airnd,Z,'linear','extrap');
    ana.Oz=interp1(ana.z,ana.o3nd,Z,'linear','extrap');
    % ana.Oz=interp1(ana.zz,ana.o3nd(ana.i),Z,'linear','extrap');
    
     if     mode_==0
         modestr = 'MERRA2';
    elseif mode_==0.1
         modestr = 'GEOS_FP';
    end
   
    
elseif mode==1 % C130 original
    
    % assign altitude according to airborne, analysis and std atm
    
        %z1=max(air.z):-20:min(air.z);
        %if z1(end)~=min(air.z) z1=[z1 min(air.z)]; end
        %z1 = [z1 0];
        
        % sort air levels to avoid duplicates
        [air.zz,air.i]=sort(air.z);
        [air.zz,air.ia]=unique(air.zz);
        air.i=air.i(air.ia);
        air.zz = flipud(air.zz);
        air.i  = flipud(air.i);
        
        % select unique values
        t     = air.t(air.i);     
        p     = air.p(air.i);
        h2ond = air.h2ond(air.i);
        airnd = air.airnd(air.i);
        
        % fill-in missing air values with model
        
        if sum(p==0)>10
            p_     = interp1(ana.z,ana.p,air.zz);     p=p_;
            h2ond_ = interp1(ana.z,ana.h2ond,air.zz); h2ond=h2ond_;
            airnd_ = interp1(ana.z,ana.airnd,air.zz); airnd=airnd_;
        end
        
        % o3
        
            o3nd   = interp1(ana.z,ana.o3nd,air.zz);
        
        % complete missing values with model values
        zind = interp1(ana.z,[1:length(ana.z)],max(air.zz),'nearest');
        Z1= [ana.z(1:zind-1)    ;air.zz];
        Tz1=[ana.t(1:zind-1)    ;t];
        Pz1=[ana.p(1:zind-1)    ;p];
        Hz1=[ana.h2ond(1:zind-1);h2ond];
        Az1=[ana.airnd(1:zind-1);airnd];
        Oz1=[ana.o3nd(1:zind-1); o3nd];
        
        % complete rest with std atm
        zind = interp1(std.z,[1:length(std.z)],max(Z1),'nearest');
        Z= [std.z(1:zind-1)    ;Z1];
        Tz=[std.t(1:zind-1)    ;Tz1];
        Pz=[std.p(1:zind-1)    ;Pz1];
        Hz=[std.h2ond(1:zind-1);Hz1];
        Az=[std.airnd(1:zind-1);Az1];
        Oz=[std.o3nd(1:zind-1) ;Oz1];
        
        % interp values to zero altitude
        Z  = [Z;0];     
        Tz = [Tz; 0];     Tz(end)=interp1(Z(end-1:-1:end-4),Tz(end-1:-1:end-4),Z(end),'linear','extrap');
        Pz = [Pz; 0];     Pz(end)=interp1(Z(end-1:-1:end-4),Pz(end-1:-1:end-4),Z(end),'linear','extrap');
        Hz = [Hz; 0];     Hz(end)=interp1(Z(end-1:-1:end-4),Hz(end-1:-1:end-4),Z(end),'linear','extrap');
        Az = [Az; 0];     Az(end)=interp1(Z(end-1:-1:end-4),Az(end-1:-1:end-4),Z(end),'linear','extrap');
        Oz = [Oz; 0];     Oz(end)=interp1(Z(end-1:-1:end-4),Oz(end-1:-1:end-4),Z(end),'linear','extrap');

        % complete missing variables from std atm
        
        
        
    %         ana.Tz=interp1(ana.z,ana.t,Z);
    %         ana.Pz=interp1(ana.z,ana.p,Z);
    %         ana.Hz=interp1(ana.z,ana.h2ond,Z);
    %         ana.Az=interp1(ana.z,ana.airnd,Z);
    %         ana.Oz=interp1(ana.z,ana.o3nd,Z);
    
        modestr = 'C130';
   
    
    % test figure;
    %     figure;
    %     plot(ana.z,ana.airnd,'ob');hold on;
    %     plot(Z,ana.Az,'.r')       ;hold on;
    %     plot(std.z,std.airnd,'.g');hold on;
    %     legend('ana-original','ana-interp','std');

elseif mode==2 %C130-interp
    
        % assign altitude according to airborne
        z1=max(air.z):-50:0;%1800:-200:1000;
        zind = interp1(std.z,[1:length(std.z)],max(air.z),'nearest');
        Z=[std.z(1:zind-1);z1'];
        % interp air
        [air.zz,air.i]=sort(air.z);
        [air.zz,air.ia]=unique(air.zz);
        air.i=air.i(air.ia);

        % need to interp the last (lowest pressure by ana/air)
        air.Tz=interp1(air.zz,air.t(air.i),Z);     air.Tz(end)=interp1(air.zz(1:3),air.t(1:3),Z(end),    'linear','extrap');
        air.Pz=interp1(air.zz,air.p(air.i),Z);     air.Pz(end)=interp1(air.zz(1:3),air.p(1:3),Z(end),    'linear','extrap');
        air.Hz=interp1(air.zz,air.h2ond(air.i),Z); air.Hz(end)=interp1(air.zz(1:3),air.h2ond(1:3),Z(end),'linear','extrap');
        air.Az=interp1(air.zz,air.airnd(air.i),Z); air.Az(end)=interp1(air.zz(1:3),air.airnd(1:3),Z(end),'linear','extrap');

        % interp reanalysis according to airborne
        ana.Tz=interp1(ana.z,ana.t,Z);
        ana.Pz=interp1(ana.z,ana.p,Z);
        ana.Hz=interp1(ana.z,ana.h2ond,Z);
        ana.Az=interp1(ana.z,ana.airnd,Z);
        ana.Oz=interp1(ana.z,ana.o3nd,Z);

        
        modestr = 'C130interp';
    
end

% interp std
std.Tz=interp1(std.z,std.t,Z);
std.Pz=interp1(std.z,std.p,Z);
std.Hz=interp1(std.z,std.h2ond,Z);
std.Az=interp1(std.z,std.airnd,Z);
std.Oz=interp1(std.z,std.o3nd,Z);

% assign temperature,preesure, water vapor, air, o3 density
if mode==0
    
    % temp
    atm.T=ana.Tz;
%     it2=find(isnan(atm.T));
%     atm.T(it2)=std.Tz(it2);
    atm.T(1:zind)=std.Tz(1:zind);
    
    % pres
    atm.P=ana.Pz;
%     ip2=find(isnan(atm.P));
%     atm.P(ip2)=std.Pz(ip2);
    atm.P(1:zind)=std.Pz(1:zind);
%     if sum(diff(atm.P)<0)>1
%         atm.P = smooth(atm.P);
%     end
    
    % water vapor
    atm.H=ana.Hz;
%     ih2=find(isnan(atm.H));
%     atm.H(ih2)=std.Hz(ih2);
    atm.H(1:zind)=std.Hz(1:zind);
    
    % air density
    atm.air=ana.Az;
%     ia2=find(isnan(atm.air));
%     atm.air(ia2)=std.Az(ia2);
    atm.air(1:zind)=std.Az(1:zind);
    
    % o3 density
    atm.o3=ana.Oz;
%     io2=find(isnan(atm.o3));
%     atm.o3(io2)=std.Oz(io2);
    atm.o3(1:zind)=std.Oz(1:zind);
    
else
    
    % temp
    atm.T=air.Tz; 
    it1=find(isnan(atm.T));
    atm.T(it1)=ana.Tz(it1);
    it2=find(isnan(atm.T));
    atm.T(it2)=std.Tz(it2);
    
    % pres
    atm.P=air.Pz; 
    ip1=find(isnan(atm.P));
    atm.P(ip1)=ana.Pz(ip1);
    ip2=find(isnan(atm.P));
    atm.P(ip2)=std.Pz(ip2);
    
    % water vapor
    atm.H=air.Hz; 
    ih1=find(isnan(atm.H));
    atm.H(ih1)=ana.Hz(ih1);
    ih2=find(isnan(atm.H));
    atm.H(ih2)=std.Hz(ih2);
    
    % air density
    atm.air=air.Az; 
    ia1=find(isnan(atm.air));
    atm.air(ia1)=ana.Az(ia1);
    ia2=find(isnan(atm.air));
    atm.air(ia2)=std.Az(ia2);
    
    % o3 density
    atm.o3=ana.Oz; 
    io2=find(isnan(atm.o3));
    atm.o3(io2)=std.Oz(io2);
    
end

atm.Z=Z./1000.0;

%complete the rest
%atm.ratio=atm.P./std.Pz;
atm.air=gas_law(atm.P,atm.T); %interp1(std.z,std.air,Z).*atm.ratio;
atm.air_int=interp1(std.z,std.airnd,Z);
atm.ratio=atm.air./atm.air_int;
%atm.ratio2=atm.air2./atm.air_int;
%atm.o3=interp1(std.z,std.o3,Z).*atm.ratio;
atm.o2=interp1(std.z,std.o2nd,Z).*atm.ratio;
atm.co2=interp1(std.z,std.co2nd,Z).*atm.ratio;
atm.no2=interp1(std.z,std.no2nd,Z).*atm.ratio;

if plotting
 figure; plot(atm.P,Z); title('Pressure');
 figure; plot(atm.T,Z,'.b'); title('Temperature');
 figure; plot(atm.H,Z); title('Water vapor');
 figure; plot(atm.air,Z); title('air density');
 figure; plot(atm.o3,Z); title('Ozone');
 figure; plot(atm.o2,Z); title('Oxygen');
 figure; plot(atm.co2,Z); title('Carbon Dioxide');
 figure; plot(atm.no2,Z); title('Nitrogen Dioxide');
end

filen=[dir 'atmos_' daystr '_' profstr '_' modestr '.dat'];
disp(['Writing to file: ' filen]);

head1=[{'# new atmosphere file created from atmosphere based on profiles of ' modestr ' on: ' daystr ' profile number ' profstr(8)}];
head2=[{'# z(km) p(mb)   T(K) '    'air(cm^-3)  o3(cm^-3) '   ' o2(cm^-3) '   ' h2o(cm^-3) '  ' co2(cm^-3) '  ' no2(cm^-3)'}];
dlmwrite(filen,head1,'');
dlmwrite(filen,head2,'-append','delimiter','');
dlmwrite(filen,[atm.Z,atm.P,atm.T,atm.air,atm.o3,atm.o2,atm.H,atm.co2,atm.no2],'-append','delimiter','\t','precision',7);


return;

%% function to convert T,P to # of molecules
function d=TP2numd(temp,pst,gas)
% pst is pressure in hPa/mbar
% temp,in K
% R - gas constant
% Na - avogadro number
% Mw - gas molecular weight
% dens    = P/RT
% numdens = (dens/Mw)*Na
% constants
if strcmp(gas,'air')
    Mw = 28.97; 
elseif strcmp(gas,'h2o')
    Mw = 18.01; 
elseif strcmp(gas,'o3')
    Mw = 48.00;
elseif strcmp(gas,'o2')
     Mw = 32.00;
end
Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)

% calculate
dens    = (pst*100)./(R*temp);              % row=P/RT: conversion from hPa to kg/m3
numdens = dens*(1/Mw)*Na *(1000/(100)^3);   % conversion from kg/m3 to #/cm3
d = numdens;
return;

%% function to calculate the number density from pressure and temperature
function nd=gas_law(p,t)
% p in mb
% T in Kelvin
R = 8.314;           % Joule/mol/Kelvin
Na= 6.02297.*10^23.; % molecules/mol

nd=p./t./R.*Na./10000.0;
return;
