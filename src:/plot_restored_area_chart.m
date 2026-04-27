% =========================================================================
% FUNCTION: plot_restored_area_chart
% DESCRIPTION: Creates a professional, publication-quality bar chart.
% =========================================================================
function plot_restored_area_chart(x_data, y_data, legend_names)
    
    figure('Name', 'Algorithm Performance: Restored Area');
    
    % --- Create the Grouped Bar Chart ---
    b = bar(x_data, y_data, 'grouped');
    
    % --- Customize Appearance for Publication ---
    
    % Define professional colors (example: a Red-Blue-Green theme)
    color_our_method = [0.8, 0, 0];    % A strong red
    color_hcha = [0, 0.4, 0.8];        % A deep blue
    color_benchmark = [0.1, 0.5, 0.1];  % A dark green
    
    b(1).FaceColor = color_our_method;
    b(2).FaceColor = color_hcha;
    b(3).FaceColor = color_benchmark;
    
    % Set bar width to be slightly less cluttered
    b(1).BarWidth = 0.8;
    
    % --- Add Labels, Title, and Legend ---
    xlabel('Number of Deployed Nodes', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Avg. Restored Area per Node (m^2)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Effectiveness of Hole Restoration Algorithms', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Create the legend using the provided names
    legend(legend_names, 'Location', 'northwest', 'FontSize', 11, 'Interpreter', 'none');
    
    % --- Final Touches ---
    ax = gca; % Get the current axes handle
    ax.FontSize = 11; % Set font size for axis numbers
    ax.YGrid = 'on';  % Turn on the horizontal grid lines
    ax.XGrid = 'off'; % Keep vertical grid lines off
    ax.Box = 'on';    % Draw a box around the whole plot
    
    % Set Y-axis to start at 0
    ylim([0, max(y_data(:)) * 1.15]);
    
    % Optional: Add the numerical values on top of each bar
    % This can be cluttered, but is useful for precise data presentation.
    for i = 1:size(y_data, 2) % Loop through each algorithm (each set of bars)
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(ytips, 1)); % Round to one decimal place
        text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
    
    % --- Save the Figure ---
    % High-quality EPS is best for LaTeX papers
    print('restored_area_comparison_chart', '-depsc', '-r300');
    fprintf('Chart saved as restored_area_comparison_chart.eps\n');
end