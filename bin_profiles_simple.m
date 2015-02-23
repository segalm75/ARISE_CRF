%% PURPOSE:
%   binning of atmospheric profile parameters from flight data
%   this is only for one profile per run
%
% CALLING SEQUENCE:
%   [npoints,varmean] = bin_profiles_simple(alt,var)
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
%  [npoints,varmean] = bin_profiles(alt,temp)
%
% MODIFICATION HISTORY:
% Written: Michal Segal, NASA Ames, Feb 10, 2015
% 
% -------------------------------------------------------------------------

%% Start of function

function [npoints,varmean] = bin_profiles_simple(alt,var)

% set common altitude range [m]
ztop = 12000/1000;
zbot = 0;
deltaz = 50/1000;
nlayer=(ztop-zbot)/deltaz;
zlow=[zbot:deltaz:ztop-deltaz];
zhigh=[zbot+deltaz:deltaz:ztop];

% create dummy variables

    varmean = NaN(1,1);   % variable mean within profile bins          
    npoints = NaN(1,1);   % number of points averaged within each bin

% bin and average
for iz = 1:nlayer
   
        jzuse        = find(alt>=zlow(iz) & alt<=zhigh(iz));% find all indices in each layer in each profile
        %varmean      = nanmean(var(jzuse));                          
        varmean(iz)     = nanmean(var(jzuse));
        npoints(iz)     = sum(~isNaN(varmean));
end

return;