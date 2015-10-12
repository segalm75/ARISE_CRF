%% function to convert temperature to potential temperature
function theta=temp2theta(temp,pst)
% pst is pressure in hPa/mbar
% P0 is 1000mbar/1000hPa
% temp,in K
% R - gas constant
% Cp, specific heat
% const=R/Cp = 0.286 for air
% theta=temp*(P0/pst)^const
% calculate
const = 0.286;           %
P0    =1000;             % hPa
P0mat =repmat(P0,length(pst),1);
theta = temp.*(P0mat./pst).^const;
return;