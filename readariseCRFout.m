% readCRFout
% subroutine to read CRF output
% files from arise libRadtran runs
%------------------------------------------
%% Details of the function:
% NAME:
%   readariseCRFout
% 
% PURPOSE:
%   read relevant parameters from libRadtran .out files
%   save relevant parameters whether they are clear/cloudy sky
% 
%
% CALLING SEQUENCE:
%  function [d] = readariseCRFout
%
% INPUT:
%  - .out libradtran files
% 
% 
% OUTPUT:
%  - d struct: clear/cloudt TOA and surface quantities
%
% DEPENDENCIES: 
%
% NEEDED FILES:
%  - *allsky*.out/*clear*.out CRF libradtran files
%
% EXAMPLE:
%
%
% MODIFICATION HISTORY:
% Written: Michal Segal-Rozenhaimer, NASA Ames, Nov, 3, 2014
% MS, Modified, Mar-09, 2015: modified to read arise CRF simulations
%
% -------------------------------------------------------------------------
function d = readariseCRFout

%% read in output files
folders = {'MERRABASE','MERRABASE_noCloudBelow','MERRABASE_wICE','MERRABASE_noFogwIceAlb'};
cases   = {'MERRABASE','MERRABASE_noFog','MERRABASE_wAlb','MERRABASE_noFogwAlb'};% base case; no fog cloud below; weighted albedo
wln     = {'lw','sw'};
sky     = {'allsky','clear'};

d.header = {'lambda_nm','lambda_nm','sza','zout','direct','global','diff_down','diff_up','net_down','sum_up'};

for i=1:length(folders)
    outfolder = strcat('F:\ARISE\ArcticCRF\libRadTran_output\arise\',folders{:,i},'\');
    for j=1:length(wln)
        for k=1:length(sky)
            list = ls(strcat(outfolder,'\','*',wln{:,j},'*',sky{:,k},'*'));
            %% read data and store variable
                for kk=1:length(list)
                    tmp = load(strcat(outfolder,list(kk,:)));
                    if strcmp(sky{:,k},'allsky')
                        if strcmp(list(kk,34),'_')
                            dateprofstr = strcat('day',list(kk,17:24),list(kk,26:33));
                        else
                            dateprofstr = strcat('day',list(kk,17:24),list(kk,26:34));
                        end
                    else
                        if strcmp(list(kk,33),'_')
                            dateprofstr = strcat('day',list(kk,16:23),list(kk,25:32));
                        else
                            dateprofstr = strcat('day',list(kk,16:23),list(kk,25:33));
                        end
                    end
                    d.(wln{:,j}).(sky{:,k}).(cases{:,i}).(dateprofstr) = tmp;
                    clear tmp;
                end
        end
    end
end


%% save to new struct and data files
casestr = [];
for i=1:length(cases)
    casestr = strcat(casestr,'_',cases{:,i},'_');
end
fiout = strcat('F:\ARISE\ArcticCRF\libRadTran_output\arise\','cases',casestr,'.mat');
save(fiout,'-struct','d');


