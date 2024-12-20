% To analyze importance of values for b and d parameters in the number of elements in the refined region
% Seeks to identify appropriate values for b and d for LocallyRefinedGeometry_TestScript.m

clear all
clc

% Problem Specifications
Ly = 0.0610*10;
ny = 64;
djet = 0.0061;
njet = 32;
dy_min = djet/njet;

b = 1;              % Determines max value of dy
c = dy_min + b;     % Determines min value of dy
d = 0.95;
b_idx = 0;


for b = 0.1:0.01:1
    b_idx = b_idx + 1;
    c = dy_min + b;     % Determines min value of dy
    d_idx = 0;

    for d = 0.4:0.01:1
        d_idx = d_idx + 1;
    
        % Loop over upper domain 
        for j = ny/2 + 2: ny +1
            jmask = (j - 0.5*ny - 1);
            dy(j-1) = b*tanh(d*jmask - (njet/2)/d) + c;        
        end
            
        % Mirror domain
        for j = 1:ny/2
            dy(j) = dy(ny - j);
        end
        
        % Scale b so that the sum of dy = Ly
        scale_factor = (Ly -  ny*dy_min)/sum(dy);
        dy = (dy-dy_min)*scale_factor + dy_min;
                
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
        maxchange_dy(b_idx,d_idx) = max(change_dy);
        num_el(b_idx,d_idx) = 0;
        
        for j = 2:ny+1
            if (y(j) >= -0.5*djet & y(j) <= 0.5*djet)
                num_el(b_idx,d_idx) = num_el(b_idx,d_idx) +1;
            end
        end


    end
end

dspan = 0.4:0.01:1; %0.25:0.01:1;
bspan = 0.1:0.01:1; %0.01:0.01:1;
legend_array = {};
figure
hold on
% Plot number of elements in domain for every value of b
for n = 1:length(bspan)
    % Plot number of elements in domain for every value of d
    plot(dspan, num_el(n,:))
    legend_array{n} = "b = " + num2str(bspan(n));
end

legend(legend_array)

% fprintf("Maximum change in dy: %0.5f \n", maxchange_dy)
% fprintf("Number of elements : %0.0f \n", num_el)
% 
% fprintf("Max y : %0.4f \n", max(y))
% fprintf("Min y : %0.4f \n", min(y))
% fprintf("Sum of dy : %0.4f \n", sum(dy))



