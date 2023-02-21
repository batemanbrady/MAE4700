% Brady Bateman
% MAE 4720/7720
% Homework 1
% Due 2023/02/24

clear
syms V gamma h s % state variables
syms n % control
syms g S C_D0 rho_SL H_s  m K % constants

assume([V gamma h s n g S C_D0 rho_SL H_s  m K],'real')

rho = rho_SL*exp(-h/H_s);
C_L = (n*g)/(rho*V^2*S/(2*m)) ;
C_D = C_D0 + K*C_L^2;
D = 0.5*rho*V^2*S*C_D;


Vdot = (-D/m) - g*sin(gamma);       %F1
gammadot = (g/V)*(n-cos(gamma));    %F2
hdot = V*sin(gamma);                %F3
sdot = V*cos(gamma);                %F4


dF1dV = diff(Vdot,V); dF1dgamma = diff(Vdot,gamma);dF1dh = diff(Vdot,h);dF1ds = diff(Vdot,s);
dF2dV = diff(gammadot,V); d21dgamma = diff(gammadot,gamma);dF2dh = diff(gammadot,h);dF2ds = diff(gammadot,s);
dF3dV = diff(hdot,V); dF3dgamma = diff(hdot,gamma);dF3dh = diff(hdot,h);dF3ds = diff(hdot,s);
dF4dV = diff(sdot,V); dF4dgamma = diff(sdot,gamma);dF4dh = diff(sdot,h);dF4ds = diff(sdot,s);

dF1dn = diff(Vdot,n); dF2dn = diff(gammadot,n); dF3dn = diff(hdot,n); dF4dn = diff(sdot,n);

A=[dF1dV dF1dgamma dF1dh dF1ds ;
dF2dV  d21dgamma dF2dh dF2ds ;
dF3dV  dF3dgamma dF3dh dF3ds ;
dF4dV  dF4dgamma dF4dh dF4ds];

B = [dF1dn dF2dn dF3dn dF4dn]'; 
