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
% 
% ---------------------------------------------------------------------------
%% routine
%  makeCloudFile4RT
   clear all;close all;
   
%% load data base

    wdir = 'F:\ARISE\C-130data\Met\SeaIceProfiles\';
    wfile= 'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-10-26.mat';
    dat    = load([wdir wfile]);
    nFields= sum(~structfun(@isempty,dat));
    airfieldNames= fieldnames(dat);
    
    clist = {'air_wc_10','air_wc_15','air_wc_20','air_ic_15','air_ic_20','air_ic_30',...
             'mod_wc_10','mod_wc_15','mod_wc_20','mod_ic_15','mod_ic_20','mod_ic_30',...
             'air_mx_1' ,'air_mx_2', 'air_mx_3',...
             'mod_mx_1' ,'mod_mx_2', 'mod_mx_3'};
         
         
%% loop over all options

    for i=1:length(clist)
            
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
                                            
                                              % genModCloudProf(prof,dprof,allprof,daystr,doy,datsource,ctype,cReff);
                                                genModCloudProf(dat.(daystr).(profstr),j,k,daystr(4:11),num2str(doy),op{1},op{2},str2double(op{3}));
                                            
                                    end % nProfiles
                                
                            end % nProfiles>0
                            
                    end % number of fields/days
            
            
            
            
        
    end
         
         
