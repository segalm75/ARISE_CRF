%% function to convert mass mixing ratio [kg/kg] to molec/cm3
function molec=mmr2molec(plev,mmr,gas,temp)

% mmr is in Kg/Kg
% gas is gas name string
% temp is Temperature in K

Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)
mwair= 28.97;           % Molecular weight of dry air

if strcmp(gas,'h2o')
    mw = 18.01;           % Molecular weight of water
elseif strcmp(gas,'o3')
    mw = 48.00;           % Molecular weight of ozone
end

% calculate
ppmv = mmr*(mwair/mw)*1e6;              % conversion from mmr to ppm
partp= ppmv.*plev*1e-6;                 % conversion from ppmv to hPa
dens = partp*100./(R*temp);             % row=P/RT: conversion from hPa to kg/m3
molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3

return;      