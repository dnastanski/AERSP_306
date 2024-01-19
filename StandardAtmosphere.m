







%%

function [T,P,rho]=StandardAtmosphere(h)

%Calculates atmosphere properties for ISA on a standard day.

%Inputs
%h = 0; %altitude (in meters)

%outputs
%T temperature (in K)
%P pressure (in Pa)
%rho: density (in kg/m3)

g = 9.80665;
R = 287.05;

if and(h>0,h<=11000) % we're in the troposphere

    h0 = 0;
    a = -6.5/1000;
    T0 = 288.15;
    P0 = 101325;
    rho0 = 1.225;
  
    T = T0 + a(h-h0); 
    P = P0*(T/T0)^(-g/(a*R)); 
    rho = rho0*(T/T0)^((-g/(a*R))-1); 

    fprintf('troposphere!!\n'); % you can uncomment this line if you wish
    return;
elseif and(h>11000,h<=20000) % we're in the tropopause
    
    h0 = 11000;
    T0 = 216.65;
    P0 = 22632;
    rho0 = 0.3639;

    T = T0;
    P = P0*exp((-g/(R*T))*(h-h0));
    rho = rho0*exp((-g/(R*T))*(h-h0));

    fprintf('tropopause!!\n');
    return;
elseif and(h>20000,h<=32000) % we're in the lower stratosphere
    
    h0 = 20000;
    a = 1/1000;
    T0 = 216.65;
    P0 = 5475;
    rho0 = 0.08803;
  
    T = T0 + a(h-h0); 
    P = P0*(T/T0)^(-g/(a*R)); 
    rho = rho0*(T/T0)^((-g/(a*R))-1);

    fprintf('stratosphere!!\n'); % you can uncomment this line if you wish
    return;
elseif and(h>32000,h<=47000) % we're in the upper stratosphere
    
    h0 = 32000;
    a = 1/1000;
    T0 = 228.65;
    P0 = 868.019;
    rho0 = 0.013225;
  
    T = T0 + a(h-h0); 
    P = P0*(T/T0)^(-g/(a*R)); 
    rho = rho0*(T/T0)^((-g/(a*R))-1);
    
    return;
elseif and(h>47000,h<=51000) % we're in the stratopause
% lapse rate is 0.
        
    h0 = 47000;
    T0 = 270.65;
    P0 = 110.906;
    rho0 = 0.00142753;

    T = T0;
    P = P0*exp((-g/(R*T))*(h-h0));
    rho = rho0*exp((-g/(R*T))*(h-h0));

    return;
elseif and(h>51000,h<=71000) % we're in the lower mesosphere
% lapse rate is -2.8K/km
        
    h0 = 51000;
    a = -2.8/1000;
    T0 = 270.65;
    P0 = 66.9389;
    rho0 = 0.0000172341;
  
    T = T0 + a(h-h0); 
    P = P0*(T/T0)^(-g/(a*R)); 
    rho = rho0*(T/T0)^((-g/(a*R))-1);
    
    return;
elseif and(h>71000,h<=85000) % we're in the upper mesosphere
% lapse rate is -2.0K/km
        
    h0 = 71000;
    a = -2.0;
    T0 = 214.65;
    P0 = 3.95642;
    rho0 = 0.0000642110;
  
    T = T0 + a(h-h0); 
    P = P0*(T/t)^(-g/(a*R)); 
    rho = rho0*(T/T0)^((-g/(a*R))-1);
    
    return;
elseif and(h>85000,h<=90000) % we're in the mesopause. The standard atmosphere definition ends at h=86km.
% lapse rate is zero
        
    h0 = ;
    T0 = ;
    P0 = ;
    rho0 = ;

    T = T0;
    P = P0*exp((-g/(R*T))*(h-h0));
    rho = rho0*exp((-g/(R*T))*(h-h0));

    return;
else % too high or too low! Return not-a-numbers.
    fprintf('Altitude out of bounds!!\n');
    T=NaN;
    P=NaN;
    rho=NaN;
    return;

end


end

