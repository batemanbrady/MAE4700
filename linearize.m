function [A,B] = linearize(state,n)

g=32.174;           % Earth gravitational acceleration, ft/s^2
rho_SL=0.002377;    % Sea level density, slugs/ft^3
H_s=31500;          % Scale height, ft
m=4200;             % RLV mass, slugs
S=2370;             % RLV wing area, slugs
C_D0=0.048;         % Zero-lift drag coefficient, unitless
K=0.218;            % Induced drag coefficient, unitless

V=state(1);
gamma=state(2) ;
h=state(3);
s=state(4);

a11=(8*K*g^2*m*n^2*exp(h/H_s))/(S*V^3*rho_SL) - (S*V*rho_SL*exp(-h/H_s)*(C_D0 + (4*K*g^2*m^2*n^2*exp((2*h)/H_s))/(S^2*V^4*rho_SL^2)))/m;
a12=-g*cos(gamma);
a13=(S*V^2*rho_SL*exp(-h/H_s)*(C_D0 + (4*K*g^2*m^2*n^2*exp((2*h)/H_s))/(S^2*V^4*rho_SL^2)))/(2*H_s*m) - (4*K*g^2*m*n^2*exp(h/H_s))/(H_s*S*V^2*rho_SL);
a14=0;
a21=-(g*(n - cos(gamma)))/V^2;
a22=(g*sin(gamma))/V;
a23=0;
a24=0;
a31=sin(gamma);
a32=V*cos(gamma);
a33=0;
a34=0;
a41=cos(gamma);
a42=-V*sin(gamma);
a43=0;
a44=0;

b1=-(4*K*g^2*m*n*exp(h/H_s))/(S*V^2*rho_SL);
b2=g/V; b3=0; b4=0;

A=[a11 a12 a13 a14;
    a21 a22 a23 a24;
    a31 a32 a33 a34;
    a41 a42 a43 a44];

B=[b1 b2 b3 b4]';
end