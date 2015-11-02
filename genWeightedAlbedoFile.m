%% Details of the function:
%  NAME:
% genWeightedAlbedoFile
%----------------------------------
% PURPOSE:
%  - generates weighted albedo file for libRadTran RT runs
%
% CALLING SEQUENCE:
%   [filen] = genWeightedAlbedoFile
%
% INPUT:
%  - none at command line
% 
% 
% OUTPUT:
%  creates weighted albedo .dat file with dim nx2 [wavelength, weighted albedo]
%
%
% DEPENDENCIES:
%  - startup_plotting.m
%  - save_fig.m
%
% NEEDED FILES/INPUT:
% -
%
% EXAMPLE:
%  - [filen] = genWeightedAlbedoFile
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Mar-06,2015  
% MS, 2015-10-27, added new profile run with GEOS_FP
%                 changed iceconc to iceconc_mean
% ---------------------------------------------------------------------------
%% function routine
function [filen] = genWeightedAlbedoFile

version_set('1.0');
startup_plotting;


%% load arise struct .mat that contains ice_conc per profile
arisedir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
filename = 'ARISEairprocessed_with_insitu_woRH';
%file2load=[strcat(arisedir,filename,'_w_anacompare_w_cldflag_noCloudBelow_andProperties20150306.mat')];
file2load=[arisedir,'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-10-26.mat'];
s = load(file2load);

%% load "pure" (water/ice) albedo files
sw_albedo_water = load('F:\ARISE\ArcticCRF\METdata\albedo\sw_albedo_water.dat');
sw_albedo_ice   = load('F:\ARISE\ArcticCRF\METdata\albedo\sw_albedo_ice.dat');
lw_albedo_water = load('F:\ARISE\ArcticCRF\METdata\albedo\lw_albedo_water.dat');
lw_albedo_ice   = load('F:\ARISE\ArcticCRF\METdata\albedo\lw_albedo_ice.dat');

%% read in each profile and generate albedo file
% find number of days analyzed
nFields    = sum(~structfun(@isempty,s));
starfieldNames  = fieldnames(s);

for i=1:nFields
    % Extract only the "profnum" fields
    names  = fieldnames(s.(starfieldNames{i,:}));
    pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
    nProfiles = length(pnames);%max(star.(starfieldNames{i,:}).prof.nprof);
    % load ice-conc data
    daystr=starfieldNames{i,:};
    
    if nProfiles>0
      
        for j=1:nProfiles
            profstr    = strcat('profnum',num2str(j));
            iceconc    = (s.(starfieldNames{i,:}).(profstr).iceconc_mean)/100;
            
            % sw albedo
            sw_wln        = sw_albedo_water(:,1);%in nm
            w_sw_albedo   = iceconc*sw_albedo_ice(:,2) + (1-iceconc)*sw_albedo_water(:,2);
            sw_w_albedo_file = strcat('sw_albedo_ice_' ,num2str(round(iceconc*100)),'.dat');
            
            % plot original and weighted
            figure(333);
            plot(sw_albedo_water(:,1),sw_albedo_water(:,2),'-b');hold on;
            plot(sw_albedo_ice(:,1)  ,sw_albedo_ice(:,2)  ,'-c');hold on;
            plot(sw_albedo_water(:,1),w_sw_albedo         ,'--','color',[0.5 0.5 0.5]);hold on;
            xlabel('wavelength');ylabel('SW Surface Albedo');title([daystr(4:11) ' ' profstr ' iceconc= ' num2str(round(iceconc*100))]);
            legend('water','ice','weighted');axis([200 4700 0 1]);
            fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\SW\sw_', daystr(4:11), profstr,'iceconc',num2str(round(iceconc*100)))];
            save_fig(333,fi,false);
            close(333);
            
            % save SW albedo file
            sw_filen = ['F:\ARISE\ArcticCRF\METdata\albedo\SW\',sw_w_albedo_file];
            disp( ['Writing to file: ' sw_filen]);
            savedat = [sw_wln,w_sw_albedo];
            save(sw_filen,'-ASCII','savedat');
            
            % lw albedo
            lw_wln        = lw_albedo_water(:,1);%in nm
            w_lw_albedo   = iceconc*lw_albedo_ice(:,2) + (1-iceconc)*lw_albedo_water(:,2);
            lw_w_albedo_file = strcat('lw_albedo_ice_' ,num2str(round(iceconc*100)),'.dat');
            
            % plot original and weighted
            figure(334);
            plot(lw_albedo_water(:,1),lw_albedo_water(:,2),'-b');hold on;
            plot(lw_albedo_ice(:,1)  ,lw_albedo_ice(:,2)  ,'-c');hold on;
            plot(lw_albedo_water(:,1),w_lw_albedo         ,'--','color',[0.5 0.5 0.5]);hold on;
            xlabel('wavelength');ylabel('LW Surface Albedo');title([daystr(4:11) ' ' profstr ' iceconc= ' num2str(round(iceconc*100))]);
            legend('water','ice','weighted');axis([4500 35000 0 0.1]);
            fi=[strcat('F:\ARISE\ArcticCRF\METdata\albedo\LW\lw_', daystr(4:11), profstr,'iceconc',num2str(round(iceconc*100)))];
            save_fig(334,fi,false);
            close(334);
            
            % save LW albedo file
            lw_filen = ['F:\ARISE\ArcticCRF\METdata\albedo\LW\',lw_w_albedo_file];
            disp( ['Writing to file: ' lw_filen]);
            savedat = [lw_wln,w_lw_albedo];
            save(lw_filen,'-ASCII','savedat');
            
        end
    end
end


