%% function to convert geopotential height to geometrical height [m]
function z=geo2z(geoH)
    r    = 6378*1000;       % Earth radius
    z = (geoH*r)./(r - geoH); 
return;