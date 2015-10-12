%% PURPOSE:
%   binning of atmospheric profile parameters from flight data
%   this is only for one profile per run according to reanalysis
%   altitude levels (can be onto any other altitude grid as well)
%
% CALLING SEQUENCE:
%   [npoints,varmean,varstat] = bin_profiles_simple(z,alt,var)
%   z is altitude array from reanalysis (m)
%   alt and var need to have same size
%
% INPUT:
%   alt is atitude (input in km), var is variable to bin (can be alt, temp, wc etc.)
%   reanalysis z is in meters, hence need to be converted
% 
% OUTPUT:
%  - binned atmospheric parameters
%  - varflag is 0/1 flag for each of the altitude bins
%  - varsum is number of points found in each bin
%  - varmean is the avg value of the var in each bin (non relevant for
%    logical variables such as cloud flag)
%
% DEPENDENCIES:
%
%
% NEEDED FILES:
% 
%
% EXAMPLE:
%  [varflag,varsum,varmean] = bin_profiles_ana(ana.z,alt,temp)
%
% MODIFICATION HISTORY:
% Written: Michal Segal, NASA Ames, Mar, 18, 2015
% MS, 2015-09-30, added varstat for mean values in bin
% 
% -------------------------------------------------------------------------

%% Start of function

function [varflag,varsum,varmean] = bin_profiles_ana(z,alt,var)

% set common altitude range [km]
nlayer = length(z);
zlow=[0;z(1:end-1)];zlow  = zlow/1000;
zhigh=[z(1:end)];   zhigh = zhigh/1000;

% create dummy variables

    varmean = NaN(1,1);   % variable mean within profile bins 
    varsum  = NaN(1,1);   % number of points averaged within each bin
    varflag = NaN(1,1);   % 0/1 flag (1 = the quantity exists in the layer)
    

% bin and average
for iz = 1:nlayer
   
        jzuse           = find(alt>=zlow(iz) & alt<=zhigh(iz));% find all indices in each layer in each profile
        varmean(iz)     = nanmean(var(jzuse));                          
        varsum(iz)      = sum(logical(var(jzuse)~=0));
        varflag(iz)     = logical(varsum(iz)>0);
end

return;