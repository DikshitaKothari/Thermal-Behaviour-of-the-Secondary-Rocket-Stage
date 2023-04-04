function Tss = laplaceanalytical()

% File name: laplaceanalytical
%
% Description: This file is the analytical solution to the Laplace equation describing the temperature profile for the rocket.
%
% Input parameters:
% T0       - concentration at the bottom boundary
% T1       - concentration at the top boundary
%
%               T1
%          _____________
%         |             |
%         |             |
%     0   |             |  0
%         |             |
%         |_____________|
%
%               T0
%

% Spatial parameters
L = 13.8;                                     % Length of rocket (z dimension)
D = 3.7;                                      % Diameter of rocket (x dimension)
dL = 0.01;                                    % Spatial step size in z
dD = 0.01;                                    % Spatial step size in x
nL = round((L/dL)+1);                         % Number of grid points in z
nD = round((D/dD)+1);                         % Number of grid points in x
z = 0:dL:L;                                   % Vector of grid points in z
x = 0:dD:D;                                   % Vector of grid points in x
nfs = 60;                                     % Number of Fourier terms

% Temperatures
T0 = 212;                                   % Temperature at bottom boundary
T1 = 3382;                                  % Temperature at top boundary

% Calculate Fourier Coefficients
B0 = zeros(1,nfs);                          % Bottom side
B1 = zeros(1,nfs);                          % Top side
for m=1:nfs
    B0(m) = (2*T0)/(m*pi)*(1-cos(m*pi));    % Bottom side 
    B1(m) = (2*T1)/(m*pi)*(1-cos(m*pi));    % Top side
end

% Initialise solution array
Tss = zeros(nD,nL);

% Analytical solution
for i=1:nD
   for j=1:nL
      for m=1:nfs
         % Fill in the analytical solution here
         Tss(i,j) = Tss(i,j) + B1(m)*(sinh((m*pi*z(j))/D)/sinh((m*pi*L/D))*sin(m*pi*x(i)/D)) + B0(m)*(sinh((m*pi*(z(j)-L))/D)/sinh((m*pi*-L/D))*sin(m*pi*x(i)/D));
      end
   end
end

% Display the steady-state result
pcolor(x,z,flip(Tss',1)),shading interp,title('Temperature (Steady State)'),xlabel('x in meters'),ylabel('z in meters'),colorbar;

end
