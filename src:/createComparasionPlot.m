% =========================================================================
% --- Algorithm Performance Comparison Plot ---
% =========================================================================
clc; clear; close all;

% --- Paste Data From Your Simulation Runs Here ---

% Data from your NOVEL Apollonius Method (Aggressive Adaptive)
% Format: [iteration_number, nodes_added, coverage_gained]
apollonius_data = [
    1, 24, (58.80 - 31.18);
    2, 49, (82.91 - 58.80);
    3, 25, (85.88 - 82.91);
    4, 3,  (89.56 - 85.88);
    5, 6,  (91.51 - 89.56);
    6, 1,  (91.57 - 91.51)
];

% Data from the BASELINE Circumcenter Method
% Format: [iteration_number, nodes_added, coverage_gained]
circumcenter_data = [
    1, 40, (64.81 - 31.18);
    2, 87, (88.55 - 64.81);
    3, 22, (90.37 - 88.55);
    4, 3,  (91.36 - 90.37);
    5, 1,  (92.25 - 91.36);
    6, 2,  (92.85 - 92.25)
];


% --- Calculate Per-Node Efficiency ---
% Efficiency = Coverage Gained / Nodes Added

% Apollonius Efficiency
apollonius_iterations = apollonius_data(:, 1);
apollonius_efficiency = apollonius_data(:, 3) ./ apollonius_data(:, 2);
% Handle any case where 0 nodes were added to avoid division by zero
apollonius_efficiency(apollonius_data(:, 2) == 0) = 0;

% Circumcenter Efficiency
circumcenter_iterations = circumcenter_data(:, 1);
circumcenter_efficiency = circumcenter_data(:, 3) ./ circumcenter_data(:, 2);
circumcenter_efficiency(circumcenter_data(:, 2) == 0) = 0;


% --- Create the Professional Plot ---
figure;
hold on;
grid on;

plot(apollonius_iterations, apollonius_efficiency, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Novel Apollonius Method');
plot(circumcenter_iterations, circumcenter_efficiency, '--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Baseline Circumcenter Method');

% --- Add Labels and Title ---
title('Algorithm Efficiency Comparison: Coverage Gained per Node', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Restoration Iteration Number', 'FontSize', 12);
ylabel('Coverage Gained per Deployed Node (%)', 'FontSize', 12);
legend('show', 'Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11); % Set axis tick font size

hold off;
