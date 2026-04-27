% =========================================================================
% 3.m
%
% This script creates a specific scenario with three non-intersecting circles
% and calls a solver to find the tangent solution circle.
% =========================================================================
clc;
clear;
close all;

% --- Step 1: define the properties of the three circles ---
Rs = 30; % define the radius for all circles.
x_loc = [50, 150, 100];
y_loc = [50, 60, 180];
C1_center = [x_loc(1), y_loc(1)];
C2_center = [x_loc(2), y_loc(2)];
C3_center = [x_loc(3), y_loc(3)];
triplet_ids = [1, 2, 3];

% --- Step 2: print the (x, y) coordinates of each circle's center ---
disp('--- Initial Circle Setup ---');
fprintf('Circle 1 (S1) Center: (%.2f, %.2f)\n', C1_center(1), C1_center(2));
fprintf('Circle 2 (S2) Center: (%.2f, %.2f)\n', C2_center(1), C2_center(2));
fprintf('Circle 3 (S3) Center: (%.2f, %.2f)\n', C3_center(1), C3_center(2));
fprintf('All circles have a radius Rs = %.2f\n\n', Rs);

% --- Step 3: plot the initial three circles ---
figure;
hold on;
grid on;
axis equal;
title('Apollonius Problem: Three Circles and Solution');
xlabel('X-axis');
ylabel('Y-axis');

plot_circle(C1_center, Rs, 'k', '-', 1.5, 'S_1');
plot_circle(C2_center, Rs, 'k', '-', 1.5, 'S_2');
plot_circle(C3_center, Rs, 'k', '-', 1.5, 'S_3');

xlim([min(x_loc)-2*Rs, max(x_loc)+2*Rs]);
ylim([min(y_loc)-2*Rs, max(y_loc)+2*Rs]);

% --- Step 4: run the solver ---
disp('--- Solving the Apollonius Problem ---');
% This assumes your final, working solveApolloniusForTriplet.m is available
[optimal_center, solution_radius] = solveApolloniusForTriplet(triplet_ids, x_loc, y_loc, Rs);

% --- Step 5: display the results and plot the solution circle ---
if ~isempty(optimal_center)
    fprintf('Solution Found!\n');
    fprintf('  Optimal Center: (%.2f, %.2f)\n', optimal_center(1), optimal_center(2));
    fprintf('  Solution Radius: %.2f\n', solution_radius);
    
    disp('Plotting the solution circle...');
    
    % Define points for the solution circle
    theta = linspace(0, 2*pi, 100);
    circle_x = optimal_center(1) + solution_radius * cos(theta);
    circle_y = optimal_center(2) + solution_radius * sin(theta);
    
    % light blue fill
fill(circle_x, circle_y, [0.8 0.9 1.0], 'FaceAlpha', 0.7, 'LineStyle', 'none');
25    
    %   blue dashed outline
    plot(circle_x, circle_y, 'b--', 'LineWidth', 2);
    
    % mark center with a filled red square
    plot(optimal_center(1), optimal_center(2), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
else
    fprintf('No solution could be found for the given three circles.\n');
end

hold off;
legend('Initial Circles', '', '', '', 'Solution Outline', 'Solution Center', 'Location', 'best');

function plot_circle(center, radius, color, style, width, label)
    theta = linspace(0, 2*pi, 100);
    circle_x = center(1) + radius * cos(theta);
    circle_y = center(2) + radius * sin(theta);
    % plot outline
    plot(circle_x, circle_y, 'Color', color, 'LineStyle', style, 'LineWidth', width);
    % plot center point
    plot(center(1), center(2), 'k.', 'MarkerSize', 15);
end
