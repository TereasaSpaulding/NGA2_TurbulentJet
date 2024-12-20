% To test method for creating locally refined mesh before implementation in FORTRAN

clear all
clc

% Problem Specifications
Ly = 0.0610*4;
ny = 64;
djet = 0.0061;
njet = 20;
dy_min = djet/njet;


% Function b*tanh(a*x-njet/(2d)) + c varies -1 to 1 on x in [0, ~ a/d + b]

% Function Parameters
b = 1;              % Determines max value of dy -> altered later by scale_factor
c = dy_min + b;     % Determines min value of dy
d = 0.85;           % Determined through preliminary sensitivity analysis for Ly ~ 10D

% Loop over upper domain 
for j = ny/2 + 2: ny +1
    jmask = (j - 0.5*ny - 1);
    dy(j-1) = b*tanh(d*jmask - (njet/2)/d) + c;        
end
    
% Mirror domain
for j = 1:ny/2
    dy(j) = dy(ny - j +1);
end

% Scale b so that the sum of dy = Ly
scale_factor = (Ly -  ny*dy_min)/sum(dy);
dy = (dy-dy_min)*scale_factor + dy_min;

% Plot upper domain - debugging
% plot(dy(ny/2 + 1: ny)*10^4, "-o")

% Build y grid
y(1) = -0.5*Ly;
for j = 2:ny+1
    y(j) = y(j-1) + dy(j-1);
end


% For debugging: 
% - check max change in dy 
% - number of elements in jet domain
% - sum(dy) = Ly

change_dy = dy(2:end) - dy(1:end-1);
maxchange_dy = max(change_dy);
num_el = 0;

for j = 2:ny+1
    if (y(j) >= -0.5*djet & y(j) <= 0.5*djet)
        num_el = num_el +1;
    end
end


figure
plot(dy, "-o")
fprintf("Maximum change in dy: %0.5f \n", maxchange_dy)
fprintf("Number of elements : %0.0f \n", num_el)

fprintf("Max y : %0.4f \n", max(y))
fprintf("Min y : %0.4f \n", min(y))
fprintf("Sum of dy : %0.4f \n", sum(dy))



