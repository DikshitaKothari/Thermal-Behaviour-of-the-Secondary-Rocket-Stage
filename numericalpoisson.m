clc 
clear 
close all 

onedim = importdata("onedim.mat")
b = 13.8;
nz = 257;
WB = 0;
T = flip(onedim);

maximum = max(T);
minimum = min(T);

[z, W] = poisson( b, nz, WB , T);

WW = max(W);

figure 
plot(z,W,'g+')
xlabel 'z in meters'
ylabel 'Temperature in Kelvin'
title 'One Dimensional Temperature Profile'

%%

% Calculate grid spacing
dz = b/(nz-1);

% Set equation parameters 
E = 7800; % MPa

% Empty Matrix for solution
    Wz = zeros(1, nz);

% Loop over points to calculate stress
    for j = 1:nz-1
         Wz(j) = E*(W(j+1)-W(j))/dz;
    end

% For the Neumann boundary
    Wz(nz) = 0;

% Stress equation S=E*Wz in MPa
    %S = E*Wz; 

% Plot the stress distribution
figure 
plot(z,Wz,'bo')
hold on 
title 'Stress Distribution'
xlabel 'z in meters'
ylabel 'Stress in MPa'

S = min(Wz);






