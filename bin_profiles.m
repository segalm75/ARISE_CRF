%% PURPOSE:
%   binning of atmospheric profile parameters from flight data
%
% CALLING SEQUENCE:
%   d = bin_profiles
%
% INPUT:
%   prof is a struct array of dates, profile numbers and parameters to bin
%   in the form of: prof.(strcat('met',date,'_')).param
%   dates is a date string in the form yyyymmdd
%   profnum is number of profiles in each flight
% 
% OUTPUT:
%  - binned atmospheric parameters according to pre-defined flight profiles
%
% DEPENDENCIES:
%
%
% NEEDED FILES:
% 
%
% EXAMPLE:
%  [d] = bin_profiles(air.(strcat('met',dates{:,i},'_')),profnum)
%
% MODIFICATION HISTORY:
% Written: Michal Segal, NASA Ames, Jan 28, 2015
% 
% -------------------------------------------------------------------------

%% Start of function

function [d] = bin_profiles(prof,profnum)
startup_plotting;

% set common altitude range [m]
ztop = 12000;
zbot = 0;
deltaz = 50;
nlayer=(ztop-zbot)/deltaz;
zlow=[zbot:deltaz:ztop-deltaz];
zhigh=[zbot+deltaz:deltaz:ztop];

if nlayer==length(zlow)
    nlayer=nlayer;
else
    nlayer = nlayer-1;
end

% create dummy variables

    d.zmean(    nlayer,:) = NaN(1,1);             
    d.Tmean(    nlayer,:) = NaN(1,1);
    d.Tstd(     nlayer,:) = NaN(1,1);
    d.thetamean(nlayer,:) = NaN(1,1);
    d.thetastd( nlayer,:) = NaN(1,1);
    d.dewpmean( nlayer,:) = NaN(1,1);
    d.dewpstd(  nlayer,:) = NaN(1,1);
    d.Pmean(    nlayer,:) = NaN(1,1);
    d.Pstd(     nlayer,:) = NaN(1,1);
    d.MMRmean(  nlayer,:) = NaN(1,1);
    d.MMRstd(   nlayer,:) = NaN(1,1);
    d.RHmean(   nlayer,:) = NaN(1,1);
    d.RHstd(    nlayer,:) = NaN(1,1);
    d.nprof(    nlayer,:) = NaN(1,1);

% bin and average

for iz = 1:nlayer
    zmean = [];         
    Tmean = [];
    Tstd = [];
    thetamean = [];
    thetastd = [];
    dewpmean = [];
    dewpstd = [];
    Pmean = [];
    Pstd = [];
    MMRmean = [];
    MMRstd = [];
    RHmean = [];
    RHstd = [];
    nprof = [];
    
    % daystr    = strcat('met',date,'_');
    
    for ipr = 1:profnum
        profstr      = strcat('profnum',num2str(ipr));
        jzuse        = find(prof.(profstr).z>=zlow(iz) & prof.(profstr).z<=zhigh(iz));% find all indices in each layer in each profile
        zmean_       = nanmean(prof.(profstr).z(jzuse));                              % average one profile at a time and store
        zmean        = [zmean zmean_];                                                % array with avg for each bin each profile
        Tmean_       = nanmean(prof.(profstr).Static_AirT(jzuse));                    % average one profile at a time and store
        Tmean        = [Tmean Tmean_];                                                % array with avg for each bin each profile
        Tstd_        = nanstd(prof.(profstr).Static_AirT(jzuse));                                  
        Tstd         = [Tstd Tstd_];
        thetamean_   = nanmean(prof.(profstr).theta(jzuse));                                                   
        thetamean    = [thetamean thetamean_];                                                                                      
        thetastd_    = nanstd(prof.(profstr).theta(jzuse));                                  
        thetastd     = [thetastd thetastd_]; 
        dewpmean_    = nanmean(prof.(profstr).DewPoint(jzuse));                                                   
        dewpmean     = [dewpmean dewpmean_];                                                                                      
        dewpstd_     = nanstd(prof.(profstr).DewPoint(jzuse));                                  
        dewpstd      = [dewpstd dewpstd_];   
        Pmean_       = nanmean(prof.(profstr).staticP(jzuse));                                                   
        Pmean        = [Pmean Pmean_];                                                                                      
        Pstd_        = nanstd(prof.(profstr).staticP(jzuse));                                  
        Pstd         = [Pstd Pstd_];  
        MMRmean_     = nanmean(prof.(profstr).MMR(jzuse));                                                   
        MMRmean      = [MMRmean MMRmean_];                                                                                      
        MMRstd_      = nanstd(prof.(profstr).MMR(jzuse));                                  
        MMRstd       = [MMRstd MMRstd_];  
        RHmean_      = nanmean(prof.(profstr).RH(jzuse));                                                   
        RHmean       = [RHmean RHmean_];                                                                                      
        RHstd_       = nanstd(prof.(profstr).RH(jzuse));                                  
        RHstd        = [RHstd RHstd_];  
        
        
    end
        d.zmean(iz)     = nanmean(zmean);
        d.Tmean(iz)     = nanmean(Tmean);
        d.Tstd(iz)      = nanmean(Tstd);
        d.thetamean(iz) = nanmean(thetamean);
        d.thetastd(iz)  = nanmean(thetastd);
        d.dewpmean(iz)  = nanmean(dewpmean);
        d.dewpstd(iz)   = nanmean(dewpstd);
        d.Pmean(iz)     = nanmean(Pmean);
        d.Pstd(iz)      = nanmean(Pstd);
        d.MMRmean(iz)   = nanmean(MMRmean);
        d.MMRstd(iz)    = nanmean(MMRstd);
        d.RHmean(iz)    = nanmean(RHmean);
        d.RHstd(iz)     = nanmean(RHstd);
        % # of profiles that are stored in each bin
        d.nprof(iz)     = sum(~isNaN(zmean));
end

return;