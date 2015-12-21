%% Details of the function:
%  NAME:
% makeCloudFile4RT
%----------------------------------
% PURPOSE:
%  - generates cloud summary file for RT runs input
%  - calls genModCloudProf.m for multiple options
%    
%
% CALLING SEQUENCE:
%   makeCloudFile4RT
%
% INPUT:
%  - none at the moment
% 
% 
% OUTPUT:
% 
% - generates cloud profile file summary for RT runs
%
%
% DEPENDENCIES:
% - version_set
%
% NEEDED FILES/INPUT:
% - model/aircraft profile summary file: e.g.,:
% - ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-10-13.mat
%   found in dir: F:\ARISE\C-130data\Met\SeaIceProfiles\
%
% EXAMPLE:
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer (MS), NASA Ames,Oct-13,2015  
% Modified: MS, 2015-11-04: added more options based on DOE excel table
% 
% ---------------------------------------------------------------------------
%% routine
%  makeCloudFile4RT

   clear all;close all;
   
%% load data base

% choose model
    
    model = 'MERRA2';%'GEOS-FP'% make sure to make this correspond with input .mat file

    wdir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
    % wfile= 'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-10-26.mat';% this is for GOES-FP
    % ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-11-03_GEOS.mat
    % ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-11-02_MERRA2.mat
    % wfile = 'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-11-02_MERRA2.mat';
    wfile = 'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-11-12_MERRA2.mat';
    dat    = load([wdir wfile]);
    nFields= sum(~structfun(@isempty,dat));
    airfieldNames= fieldnames(dat);
    
    % case names from excel DOE table
    
    cname = {'BASE-AIR-INT', 'BASE-MOD-ATM-INT', 'BASE-AIR-INT-MX', 'BASE-MOD-ATM-INT-MX', 'MOD-AIR-ATM-IAL10', 'MOD-AIR-ATM-IAL15','MOD-AIR-ATM-IAL20',...
             'MOD-AIR-ATM-MX1FLP','MOD-AIR-ATM-MX2FLP','MOD-AIR-ATM-MX3FLP','MOD-I10','MOD-I15','MOD-I20','MOD-MX1',...
             'MOD-MX2','MOD-MX3','MOD-IAL10','MOD-IAL15','MOD-IAL20','MOD-MX1FLP','MOD-MX2FLP','MOD-MX3FLP',...
             'AIR-L10','AIR-L15','AIR-L20','AIR-MX1','AIR-MX2','AIR-MX3',...
             'AIR-MOD-ATM-L10','AIR-MOD-ATM-L15','AIR-MOD-ATM-L20','AIR-MOD-ATM-MX1','AIR-MOD-ATM-MX2','AIR-MOD-ATM-MX3'};
         
    % case list for processing: x_y_z_w
    % x - measurement platform (air/mod)
    % y - atmospheric state profile (air/mod) - all interpolated
    % z - cloud phase (wc/ic/iawc/mx/mxflp); iawc - ice amount but as
    % liquid, mx - both phases, mxflp - ice amount as liquid 
    % and liquid amount as ice
    % w - Reff (10,15,20 or 0, i.e., aircraft in-situ)
    
    clist = {'air_air_wc_0','air_mod_wc_0','air_air_mx_0','air_mod_mx_0','mod_air_iawc_10','mod_air_iawc_15','mod_air_iawc_20',...
             'mod_air_mxflp_10','mod_air_mxflp_15','mod_air_mxflp_20','mod_mod_ic_10','mod_mod_ic_15',...
             'mod_mod_ic_20','mod_mod_mx_10','mod_mod_mx_15','mod_mod_mx_20','mod_mod_iawc_10','mod_mod_iawc_15',...
             'mod_mod_iawc_20','mod_mod_mxflp_10','mod_mod_mxflp_15','mod_mod_mxflp_20',...
             'air_air_wc_10','air_air_wc_15','air_air_wc_20','air_air_mx_10','air_air_mx_15','air_air_mx_20',...
             'air_mod_wc_10','air_mod_wc_15','air_mod_wc_20','air_mod_mx_10','air_mod_mx_15','air_mod_mx_20'};
         
         
%% loop over all options

% atm profiles does not exist for first half of campaign
% for option 2, need to interpolate Reff into model level

    for i=[5:length(clist)];%[1:2 11:22 29:34];%1:length(clist)
            
            % find case options
            
            op = strsplit(clist{:,i},'_');%op{:,1}
            
            
            %% loop over data file for each option
                    k = 0;
                    for ii=1:nFields
                        
                            daystr=airfieldNames{ii,:};
                            k = k+1;
                            doy = datestr2doy(daystr(4:11),'yyyymmdd');
                            % Extract only the "profnum" fields
                            names  = fieldnames(dat.(airfieldNames{ii,:}));
                            pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
                            nProfiles = length(pnames);
                            
                            
                            if nProfiles>0
                                
                                    for j=1:nProfiles
            
                                            profstr    = strcat('profnum',num2str(j));
                                            
                                            % generate cloud profile
                                            
                                              % genModCloudProf(model,prof,dprof,allprof,daystr,doy,datsource,atmsource,ctype,cReff);
                                                genModCloudProf(model,dat.(daystr).(profstr),j,k,daystr(4:11),num2str(doy),op{1},op{2},op{3},str2double(op{4}));
                                            
                                    end % nProfiles
                                
                            end % nProfiles>0
                            
                    end % number of fields/days
            
            
            
            
        
    end
         
         
