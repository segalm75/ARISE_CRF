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
    wfile= 'ARISEairprocessed_with_insitu_withWVparams20150921_w_anacompare_w_consolidatedcloudsAir_andModel2015-10-13.mat';
    air    = load([wdir wfile]);
    nFields= sum(~structfun(@isempty,air));
    airfieldNames= fieldnames(air);
    
    clist = {'air_lwc_10','air_lwc_15','air_lwc_20','air_iwc_15','air_iwc_20','air_iwc_30',...
             'ana_lwc_10','ana_lwc_15','ana_lwc_20','ana_iwc_15','ana_iwc_20','ana_iwc_30'};
         
         
%% loop over all options

    for i=1:length(clist)
            
            % find case options
            
            op = strsplit(clist{:,1},'_');%op{:,1}
            
            
            %% loop over data file for each option
            
                    for ii=1:nFields
    
                            % Extract only the "profnum" fields
                            names  = fieldnames(air.(airfieldNames{ii,:}));
                            pnames = names(~cellfun('isempty', strfind(names, 'profnum')));
                            nProfiles = length(pnames);
                            daystr=airfieldNames{ii,:};
                            
                            if nProfiles>0
                                
                                    for j=1:nProfiles
            
                                            profstr    = strcat('profnum',num2str(j));
                                            
                                            % generate cloud profile
                                            
                                              genModCloudProf(prof,dprof,allprof,daystr,doy,ctype,cReff);
                                            
                                    end % nProfiles
                                
                            end % nProfiles>0
                            
                    end % number of fields/days
            
            
            
            
        
    end
         
         
