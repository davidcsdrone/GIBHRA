% =========================================================================
% SCRIPT: run_hora_simulation.m
% DESCRIPTION: Main script to set up a WSN and run the HORA algorithm.
% =========================================================================
function run_hora_simulation()
    % --- Main Simulation Setup ---
    clear; clc; close all;
    rng('default'); % For repeatable results

    % --- Simulation Parameters ---
    deployment_area = [1000, 1000];
    sensing_radius = 50;
    communication_radius = 2 * sensing_radius;
    num_nodes = 300; % You can change this
    
    % --- Network Generation ---
    [x_loc, y_loc] = generate_network(deployment_area, num_nodes);
    
    % --- Plot Initial Network State ---
    figure('Name', 'HORA Simulation');
    subplot(1, 2, 1);
    plot_network_state(x_loc, y_loc, sensing_radius, 'Initial Deployment');
    
    % --- Run the HORA Restoration Algorithm ---
    fprintf('--- Starting HORA Hole Restoration ---\n');
    [x_final, y_final] = restore_network_with_HORA(x_loc, y_loc, sensing_radius, communication_radius, deployment_area);
    fprintf('--- HORA Restoration Complete ---\n');
    
    % --- Plot Final Network State ---
    subplot(1, 2, 2);
    plot_network_state(x_final, y_final, sensing_radius, 'After HORA Restoration');

end

% --- Helper function to generate the network ---
function [x, y] = generate_network(area, num_nodes)
    x = area(1) * rand(1, num_nodes);
    y = area(2) * rand(1, num_nodes);
end

% --- Helper function for plotting ---
function plot_network_state(x, y, R, title_str)
    hold off;
    % Plot sensing circles
    theta = linspace(0, 2*pi, 100);
    for i = 1:length(x)
        fill(x(i) + R*cos(theta), y(i) + R*sin(theta), [0, 0.8, 0.2], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        hold on;
    end
    % Plot node centers
    plot(x, y, 'k.', 'MarkerSize', 12);
    axis equal;
    grid on;
    box on;
    title(title_str, 'FontSize', 14);
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
end