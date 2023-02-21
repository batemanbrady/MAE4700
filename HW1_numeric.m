% Brady Bateman
% MAE 4720/7720
% Homework 1
% Due 2023/02/24
clear;clf
%syms g S C_D0 rho_SL H_s  m K % constants
g=32.174;           % Earth gravitational acceleration, ft/2^2
rho_SL=0.002377;    % Sea level density, slugs/ft^3
H_s=31500;          % Scale height, ft
m=4200;             % RLV mass, slugs
S=2370;             % RLV wing area, slugs
C_D0=0.048;         % Zero-lift drag coefficient, unitless
K=0.218;            % Induced drag coefficient, unitless

Ref=readtable("RLV_ref_traj.xls"); %Import data from an excel file
dx=readtable("RLV_dx.xls");
Ref{:,3}=Ref{:,3}*pi/180; %convert degrees to radians
dx{:,3}=dx{:,3}*pi/180;



deltaX0=[0 -0.2*pi/180 200 0]';
deltan0=0;
deltaX=zeros(size(Ref{:,2:5}));
X=zeros(size(Ref{:,2:5}));

for t=1:length(Ref{:,1})
    Ref_St=Ref{t,2:5};
    Ref_C=Ref{t,6};
    [A,B]=linearize(Ref_St,Ref_C);


    lambda=eig(A);

    STM=expm(A*10);
    deltaX(t,:)=STM*deltaX0;
    deltaX0=deltaX(t,:)';
    X(t,:)=Ref_St+deltaX(t,:);

end
str = {'\deltaV',' \delta\gamma', '\deltah', '\deltas'};
xlabstr = {'ft/s', 'rad', 'ft', 'ft'};
for i=1:4
    figure(i)
    hold on
    %plot(Ref{:,1},X(:,i))
    %plot(Ref{:,1},Ref{:,i+1})
   plot(Ref{:,1},deltaX(:,i))
    plot(dx{:,1},dx{:,i+1})
     title(str{i})
    ylabel(xlabstr{i})
    xlabel('time (s)')
    legend('Linear results','Numerical Integration')
    hold off
    
end

