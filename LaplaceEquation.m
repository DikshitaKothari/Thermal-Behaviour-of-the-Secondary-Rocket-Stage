function [X, Y, T, count] = LaplaceEquation( a, b, nx, ny, TL, TR, TT, TB )
% LaplaceEquation - Solves the 2D Laplace equation for domain 0 < x < a and
% 0 < y < b with nx and ny points in the x and y dimensions
% Input :
%    a - maximum value of the x domain [m]
%    b - maximum value of the y domain [m]
%    nx - number of points in the x domain
%    ny - number of points in the y domain
%    TL, TR, TT, TB - temperature at left, right, top, and bottom
%         boundaries [deg C]
% Output :
%    X - matrix of x-coordinate values for each point in the domain [m]
%    Y - matrix of y-coordinate values for each point in the domain [m]
%    T - steady state temperature distribution [deg C]

% Create domain
x = linspace(0, a, nx);
y = linspace(0, b, ny);
[X, Y] = meshgrid(x,y);

% Calculate grid spacing and maximum stable timestep
dx = a/(nx-1);
dy = b/(ny-1);
dt = min(dx,dy)^2/4;

% Initialise solution array for timestep n+1
% Note indexing: rows = y, columns = x to agree with meshgrid output
Tnp1 = zeros(ny, nx);

% Set boundary conditions
Tnp1(:,1) = TL;   % Left boundary
Tnp1(:,end) = TR; % Right boundary
Tnp1(1,:) = TB;   % Bottom boundary
Tnp1(end,:) = TT; % Top boundary

% Initialise error and set tolerance for convergence
err = 1;
tol = 1e-10;
count = 0;
alpha = 0.5*dy^2/(dx^2+dy^2);
gamma = 0.5*dx^2/(dx^2+dy^2);

while  err >= tol
    % Update solution array for this timestep
    Tn = Tnp1;
    
    count = count + 1;

    % Loop over internal points
    for i = 2:nx-1
        for j = 2:ny-1
            Tnp1(j,i) = (gamma*(Tn(j+1,i) + Tnp1(j-1,i)) + alpha*(Tn(j,i+1) + Tnp1(j,i-1)));
        end
    end

    % Compute error as maximum change in domain (absolute value)
    err = norm((Tnp1(:) - Tn(:)),1)*dx*dy;
end

% Set output variable
T = flip(Tnp1,1);
end