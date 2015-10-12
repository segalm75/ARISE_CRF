%% function to convert air levels to number of air/o2 molec
function molec=air2molec(plev,gas,temp)

% mmr is in Kg/Kg
% gas is gas name string
% temp is Temperature in K

Na   = 6.0221413e23;    % Avogadro number
R    = 287.05;          % gas constant J/(kgK)

if strcmp(gas,'air')
    mw = 28.97;                             % Molecular weight of dry air
    dens = plev*100 ./(R*temp);             % row=P/RT: conversion from hPa to kg/m3
    molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3
    
elseif strcmp(gas,'o2')
    mw = 32.00;                             % Molecular weight of oxygen
    dens = plev*0.21*100 ./(R*temp);        % row=P/RT: conversion from hPa to kg/m3
    molec= dens*(1/mw) *Na *(1000/(100)^3); % conversion from kg/m3 to #/cm3
end

return;      