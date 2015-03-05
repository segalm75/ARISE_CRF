%% function that examines various spectral albedo and constructs ice/water albedo
% to process in ARISE CRF libRadTran simulations
% MS, Mar-04-2015
%--------------------------------------------------------------------------------
function prepare_albedo
startup_plotting;

% load 4star/ssfr data
albedo_ice_star     = load('C:\cygwin\home\msegalro\libradtran\data\albedo\albedo_ice.dat');
albedo_snowice_star = load('C:\cygwin\home\msegalro\libradtran\data\albedo\albedo_snow_ice.dat');
albedo_20140919_2093= load('C:\cygwin\home\msegalro\libradtran\data\albedo\albedo_20140919_2093.dat');
albedo_20140919_2158= load('C:\cygwin\home\msegalro\libradtran\data\albedo\albedo_20140919_2158.dat');

% load John's Hopkins University lab data
albedo_jhu_snow24umgrain  = importdata('F:\ARISE\ArcticCRF\METdata\albedo\jhu.becknic.water.snow.fine.24um.fine.spectrum.txt');
albedo_jhu_snow82umgrain  = importdata('F:\ARISE\ArcticCRF\METdata\albedo\jhu.becknic.water.snow.granular.82um.medium.spectrum.txt');
albedo_jhu_snow178umgrain = importdata('F:\ARISE\ArcticCRF\METdata\albedo\jhu.becknic.water.snow.granular.178um.coarse.spectrum.txt');
albedo_jhu_seawater       = load      ('F:\ARISE\ArcticCRF\METdata\albedo\jhu.becknic.seawater.spectrum_dataonly.txt');
albedo_jhu_seawater       = flipud(albedo_jhu_seawater);
% plot 4star/ssfr data
figure(1);
plot(albedo_ice_star(:,1)    ,albedo_ice_star(:,2)    ,'--c'); hold on;
plot(albedo_snowice_star(:,1),albedo_snowice_star(:,2),'--b'); hold on;
plot(albedo_20140919_2093(:,1),albedo_20140919_2093(:,2),'--g'); hold on;
plot(albedo_20140919_2158(:,1),albedo_20140919_2158(:,2),'--','color',[0.5 0.5 0.5]); hold on;
xlabel('wavelength [nm]');ylabel('Surface Albedo'); legend('ice','snow-ice','20140919-2093-SSFR','20140919-2158-SSFR');
fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\','SSFR_arise_snow_ice_albedo')];
save_fig(1,fi,true);

% plot JHU lab data
figure(2);
plot(albedo_jhu_snow24umgrain.data(:,1)  ,albedo_jhu_snow24umgrain.data(:,2),'color',[0.7 0.5 0.6]);hold on;
plot(albedo_jhu_snow82umgrain.data(:,1)  ,albedo_jhu_snow82umgrain.data(:,2),'color',[0.9 0.5 0.6]);hold on;
plot(albedo_jhu_snow178umgrain.data(:,1),albedo_jhu_snow178umgrain.data(:,2),'color',[0.7 0.2 0.2]);hold on;
plot(albedo_jhu_seawater(:,1)           ,albedo_jhu_seawater(:,2)           ,'-b');hold on;
xlabel('wavelength [\mum]');ylabel('Reflectance [%]');legend('JHU-snow-24um','JHU-snow-82um','JHU-snow-178um','JHU-SeaWater');
fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\','JHU_snow_reflectance_wSeaWaterZoom')];
save_fig(2,fi,true);

% compare both for UV-VIS-SWIR region
figure(3)
plot(albedo_ice_star(:,1)/1000    ,100*albedo_ice_star(:,2)    ,'--c'); hold on;
plot(albedo_snowice_star(:,1)/1000,100*albedo_snowice_star(:,2),'--b'); hold on;
plot(albedo_jhu_snow24umgrain.data(:,1)  ,albedo_jhu_snow24umgrain.data(:,2),'color',[0.7 0.5 0.6]);hold on;
plot(albedo_jhu_snow82umgrain.data(:,1)  ,albedo_jhu_snow82umgrain.data(:,2),'color',[0.9 0.5 0.6]);hold on;
plot(albedo_jhu_snow178umgrain.data(:,1),albedo_jhu_snow178umgrain.data(:,2),'color',[0.7 0.2 0.2]);hold on;
xlabel('wavelength [\mum]');ylabel('Reflectance [%]');legend('ice-SSFR','snow-ice-SSFR','JHU-snow-24um','JHU-snow-82um','JHU-snow-178um');
axis([0.3 4 0 100]);set(gca,'Xtick',[0.4:0.2:4]);set(gca,'XtickLabel',{'0.4','0.6','0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2','3.4','3.6','3.8','4.0'});
fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\','JHU_SSFR_compared_reflectance')];
save_fig(3,fi,true);


% prepare water albedo for SW runs
savefolder = 'C:\cygwin\home\msegalro\libradtran\input\CRF\arise\MERRABASE\';
sw_water_file = 'sw_albedo_water.dat';

% set 7% to 300-1800 nm
wln = 300:1:35000;

[ix,iy] = unique(albedo_jhu_seawater(:,1)*1000);
cat_wln_water = [(300:1800)';albedo_jhu_seawater(iy,1)*1000];
cat_alb_water = [(ones(length(300:1800),1)*0.07);albedo_jhu_seawater(iy,2)/100]; % constant 7% albedo for water
albedo_water  = interp1(cat_wln_water,cat_alb_water,wln,'nearest','extrap');
sw_wln          = wln(1:3701)';
sw_water_albedo = albedo_water(1:3701)';
sw_water = [sw_wln,sw_water_albedo];
filen=[savefolder sw_water_file];
disp( ['Writing to file: ' filen]);
save(filen,'-ASCII','sw_water');

figure(111);
subplot(211);plot(wln,albedo_water,'-b');hold on;

% prepare ice albedo for SW runs
sw_ice_file   = 'sw_albedo_ice.dat';
[jx,jy]       = unique(albedo_jhu_snow178umgrain.data(:,1)*1000);
cat_wln_ice   = [albedo_jhu_snow178umgrain.data(jy,1)*1000];
cat_alb_ice   = [albedo_jhu_snow178umgrain.data(jy,2)/100];
albedo_ice    = interp1(cat_wln_ice,cat_alb_ice,wln,'nearest','extrap');
sw_wln        = wln(1:3701)';
sw_ice_albedo = albedo_ice(1:3701)';
sw_ice = [sw_wln,sw_ice_albedo];
filen=[savefolder sw_ice_file];
disp( ['Writing to file: ' filen]);
save(filen,'-ASCII','sw_ice');

figure(111)
subplot(211);plot(wln,albedo_ice,'-c');hold on;xlabel('wavelength [nm]');ylabel('SW Surface Albedo');
axis([300 4000 0 1]); legend('water','ice');


% prepare water albedo for LW runs
lw_water_file = 'lw_albedo_water.dat';
lw_wln          = wln(3702:end)';
lw_water_albedo = albedo_water(3702:end)';
lw_water = [lw_wln,lw_water_albedo];
filen=[savefolder lw_water_file];
disp( ['Writing to file: ' filen]);
save(filen,'-ASCII','lw_water');

figure(111);
subplot(212);plot(wln,albedo_water,'-r');hold on;

% prepare ice albedo for LW runs
lw_ice_file = 'lw_albedo_ice.dat';
lw_ice_albedo = albedo_ice(3702:end)';
lw_ice = [lw_wln,lw_ice_albedo];
filen=[savefolder lw_ice_file];
disp( ['Writing to file: ' filen]);
save(filen,'-ASCII','lw_ice');

figure(111)
subplot(212);plot(wln,albedo_ice,'-','color',[0.9 0.5 0.5]);hold on;xlabel('wavelength [nm]');ylabel('LW Surface Albedo');
axis([4000 35000 0 1]); legend('water','ice');

%% save figure
fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\','SW_LW_water_ice_albedo')];
save_fig(111,fi,true);

