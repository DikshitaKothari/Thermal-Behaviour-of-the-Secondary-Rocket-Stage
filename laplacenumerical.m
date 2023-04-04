clc 
clear all 

[X, Y, T, count] = LaplaceEquation( 3.7, 13.8, 257, 257, 0, 0, 3382, 212 );

% Display the steady-state result
pcolor(X,Y,T),shading interp,title('Numerical Steady State Temperature (n = 8)'),xlabel('x in meters'),ylabel('z in meters'),colorbar;
%%
a = 3.7;
b = 13.8; 
nx = 257;
ny = 257; 
x = linspace(0, a, nx);
y = linspace(0, b, ny);
[X, Y] = meshgrid(x,y);
onedim = zeros(ny,2)
for i = 1:ny 
    onedim(i,:) = mean(T(i,:));
end
size(onedim)

% figure
% plot(onedim, y, 'b+');
% hold on 
% ylabel 'z in centimeters'
% xlabel 'Temperature in Kelvin'
% title 'One Dimensional Temperature Profile'

figure
pcolor(1:2, y, onedim),shading interp,title('One Dimensional Temperature Profile'),ylabel('z in centimeters'),colorbar;
