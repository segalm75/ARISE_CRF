%% PURPOSE:
%   binning of atmospheric profile parameters from flight data
%   this is only for one profile per run according to reanalysis
%   altitude levels
%
% CALLING SEQUENCE:
%   [npoints,varmean] = bin_profiles_simple(z,alt,var)
%   z is altitude array from reanalysis
%   alt and var need to have same size
%
% INPUT:
%   alt is atitude (in km), var is variable to bin (can be alt, temp etc.)
% 
% OUTPUT:
%  - binned atmospheric parameters
%
% DEPENDENCIES:
%
%
% NEEDED FILES:
% 
%
% EXAMPLE:
%  [npoints,varmean] = bin_profiles_ana(ana.z,alt,temp)
%
% MODIFICATION HISTORY:
% Written: Michal Segal, NASA Ames, Mar, 18, 2015
% 
% -------------------------------------------------------------------------

%% Start of function

function [npoints,varsum] = bin_profiles_ana(z,alt,var)

% set common altitude range [km]
nlayer = length(z);
zlow=[0,z(1:end-1)];zlow = zlow/1000;
zhigh=[z(1:end)];zhigh = zhigh/1000;

% create dummy variables

    varmean = NaN(1,1);   % variable mean within profile bins          
    npoints = NaN(1,1);   % number of points averaged within each bin

% bin and average
for iz = 1:nlayer
   
        jzuse        = find(alt>=zlow(iz) & alt<=zhigh(iz));% find all indices in each layer in each profile
        %varmean      = nanmean(var(jzuse));                          
        varsum(iz)      = sum(var(jzuse));
        npoints(iz)     = logical(varsum(iz)>0);
end

return;