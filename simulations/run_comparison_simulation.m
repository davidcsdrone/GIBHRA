% =========================================================================
% SCRIPT: run_comparison_simulation.m
% AUTHOR: [Your Name]
% DESCRIPTION: Runs and compares three hole restoration algorithms:
%              1. Our Method (Geometric Inversion)
%              2. HCHA (Circumcenter Baseline)
%              3. Benchmark (Yao et al. / DDPCH)
%
%              It generates and saves a professional bar chart comparing
%              Metric 1: Average Restored Area per Node.
% =========================================================================
function run_comparison_simulation()
    % --- Main Simulation Setup ---
    clear; clc; close all;
    rng('default'); % For repeatable random results

    % --- Simulation Parameters ---
    deployment_area = [1000, 1000]; % Simulation area in meters (X, Y)
    sensing_radius = 50;           % Sensing radius (Rs)
    node_counts = [200, 400, 600, 800, 1000]; % X-axis of the bar chart
    num_runs = 20; % Number of random trials for statistical significance

    % --- Algorithm Setup ---
    algorithms = {
        'GIBHRA', @run_inversion_algorithm;
        'HCHA',    @run_circumcenter_algorithm;
        'HPA',      @run_yao_algorithm
        'NLCHR',  @run_nlchr_algorithm


    };

    % --- Data Storage ---
    final_avg_restored_area = zeros(length(node_counts), size(algorithms, 1));

    % --- Main Simulation Loop ---
    for i = 1:length(node_counts)
        num_nodes = node_counts(i);
        fprintf('--- Simulating with %d Initial Nodes ---\n', num_nodes);
        
        run_results_matrix = zeros(num_runs, size(algorithms, 1));

        for j = 1:num_runs
            fprintf('  Run %d/%d...\n', j, num_runs);
            
            [x_loc_initial, y_loc_initial] = generate_network(deployment_area, num_nodes);
            
            for k = 1:size(algorithms, 1)
                alg_name = algorithms{k, 1};
                alg_function = algorithms{k, 2};
                
                fprintf('    Running: %s...\n', alg_name);
            [total_restored, actions_taken] = alg_function(x_loc_initial, y_loc_initial, sensing_radius, deployment_area);
            if actions_taken > 0, avg_area_per_action = total_restored / actions_taken;
                else, avg_area_per_action = 0; 
            end
            run_results_matrix(j, k) = avg_area_per_action;
    
    run_results_matrix(j, k) = avg_area_per_action;
            end
        end
        
        final_avg_restored_area(i, :) = mean(run_results_matrix, 1);
        fprintf('  Average Restored Area for %d nodes: %s\n\n', num_nodes, mat2str(final_avg_restored_area(i, :), 4));
    end

    fprintf('\n--- Simulation Complete ---\n');

% ADD THIS CORRECTED PLOTTING BLOCK

% --- Plotting the Final Bar Charts ---
plot_metric_barchart(node_counts, final_avg_restored_area, algorithms(:,1), 'Avg. Restored Area per Node (m^2)', 'barchart_restored_area');
plot_metric_barchart(node_counts, final_avg_nodes_added, algorithms(:,1), 'Avg. Number of Restorative Nodes Added', 'barchart_deployment_cost');
% =========================================================================
% =========================================================================
% --- ALGORITHM WRAPPER FUNCTIONS ---
% The core logic from your three .m files has been placed here.
% =========================================================================
% =========================================================================

function [restored_area, nodes_added_total] = run_inversion_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    % This function encapsulates the logic from 'tutorial1.m'
    net_length = area_dims(1);
    net_width = area_dims(2);
    x_loc = x_loc_init;
    y_loc = y_loc_init;
    
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    
    TARGET_COVERAGE_PERCENT = 98.0;
    MAX_ITERATIONS = 10;
    iteration_count = 0;
    
    while true
       iteration_count = iteration_count + 1;
       no_nodes_current = length(x_loc);
       stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
       current_coverage_percent = stats_current.k_covered_percentage;
       
       if current_coverage_percent >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS
           break;
       end
       
       Ni = cell(no_nodes_current, 1);
       for i = 1:no_nodes_current
          for j = 1:no_nodes_current
              if i ~= j && norm([x_loc(i)-x_loc(j), y_loc(i)-y_loc(j)]) < 2 * Rs
                  Ni{i} = [Ni{i}, j];
              end
          end
       end
       
       AllCIPs = cell(no_nodes_current, 1);
       for i = 1:no_nodes_current
          [current_Pi, ~, ~] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width);
          AllCIPs{i} = current_Pi;
       end

       boundary_nodes_list = [];
       for i = 1:no_nodes_current
           if isempty(Ni{i}) || ~isempty(AllCIPs{i})
               boundary_nodes_list(end+1) = i;
           end
       end
       
       T = selectTripletsDelaunay(boundary_nodes_list, x_loc, y_loc, Rs, current_coverage_percent);
      
       if isempty(T)
           break;
       end
       
       x_loc_after_iter = x_loc;
       y_loc_after_iter = y_loc;
       nodes_deployed_this_iteration = 0;
      
       for t_idx = 1:length(T)
           current_triplet = T{t_idx};
           [optimal_center, ~] = solveApolloniusForTriplet(current_triplet, x_loc, y_loc, Rs);
           if ~isempty(optimal_center)
               coverage_in_wsn = estimateCoverageInWSN(optimal_center, Rs, net_length, net_width);
               if coverage_in_wsn >= 0.50
                   x_loc_after_iter(end+1) = optimal_center(1);
                   y_loc_after_iter(end+1) = optimal_center(2);
                   nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
               end
           end
       end
       
       x_loc = x_loc_after_iter;
       y_loc = y_loc_after_iter;
      
       if nodes_deployed_this_iteration == 0
           break;
       end
    end
    
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

function [restored_area, nodes_added_total] = run_circumcenter_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    % This function encapsulates the logic from 'main_circumcenter_baseline.m'
    net_length = area_dims(1);
    net_width = area_dims(2);
    x_loc = x_loc_init;
    y_loc = y_loc_init;
    
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    
    TARGET_COVERAGE_PERCENT = 98.0;
    MAX_ITERATIONS = 10;
    iteration_count = 0;

    while true
        iteration_count = iteration_count + 1;
        stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
        current_coverage_percent = stats_current.k_covered_percentage;
        
        if current_coverage_percent >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS
            break;
        end
        
        hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs);
        
        if isempty(hole_centers)
            break;
        end
        
        x_loc_after_iter = x_loc;
        y_loc_after_iter = y_loc;
        nodes_deployed_this_iteration = 0;
        
        for i = 1:size(hole_centers, 1)
            potential_center = hole_centers(i, :);
            coverage_in_wsn = estimateCoverageInWSN(potential_center, Rs, net_length, net_width);
            
            if coverage_in_wsn >= 0.5
                x_loc_after_iter(end+1) = potential_center(1);
                y_loc_after_iter(end+1) = potential_center(2);
                nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
            end
        end
        
        x_loc = x_loc_after_iter;
        y_loc = y_loc_after_iter;
        
        if nodes_deployed_this_iteration == 0
            break;
        end
    end
    
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

function [restored_area, nodes_added_total] = run_yao_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    % This function encapsulates the logic from 'DDPCH.m'
    net_length = area_dims(1);
    net_width = area_dims(2);
    x_loc = x_loc_init;
    y_loc = y_loc_init;
    Rc = 2 * Rs;
    
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    
    TARGET_COVERAGE_PERCENT = 98.0;
    MAX_ITERATIONS = 20;
    iteration_count = 0;

    while true
       iteration_count = iteration_count + 1;
       stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
       current_coverage_percent = stats_current.k_covered_percentage;
       
       if current_coverage_percent >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS
           break;
       end
       
       no_nodes_current = length(x_loc);
       all_hole_boundary_nodes = [];
       for i = 1:no_nodes_current
           if findHoleBoundaryNode_Yao(i, x_loc, y_loc, Rc)
               all_hole_boundary_nodes(end+1) = i;
           end
       end
       
       if isempty(all_hole_boundary_nodes)
           break;
       end

       boundary_coords = [x_loc(all_hole_boundary_nodes)', y_loc(all_hole_boundary_nodes)'];
       if size(boundary_coords, 1) < 3
           break;
       end
       cluster_indices = my_dbscan(boundary_coords, Rc, 3);
       
       num_clusters = max(cluster_indices);
       if num_clusters < 1
           break;
       end
       
       holes = cell(1, num_clusters);
       for i = 1:num_clusters
           holes{i} = all_hole_boundary_nodes(cluster_indices == i);
       end

       if isempty(holes)
           break;
       end
       
       hole_areas = zeros(1, length(holes));
       for i = 1:length(holes)
           current_hole_nodes = holes{i};
           if length(current_hole_nodes) >= 3
               [~, hole_areas(i)] = convhull([x_loc(current_hole_nodes)', y_loc(current_hole_nodes)']);
           end
       end
       
       [~, largest_area_index] = max(hole_areas);
       if hole_areas(largest_area_index) == 0
           break;
       end

       hole_to_patch = holes{largest_area_index};
       [x_loc, y_loc, nodes_deployed_this_iteration] = patchHole_Yao(hole_to_patch, x_loc, y_loc, Rc, Rs, net_length, net_width);
       
       if nodes_deployed_this_iteration == 0
           break;
       end
    end

    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

% =========================================================================
% =========================================================================
% --- PLOTTING AND CONSOLIDATED HELPER FUNCTIONS ---
% All helper functions from your scripts are collected below.
% =========================================================================
% =========================================================================

% =========================================================================
% --- REVISED PLOTTING FUNCTION (v3) ---
% =========================================================================
function plot_results_barchart(x_data, y_data, legend_names)
    
    figure('Name', 'Algorithm Comparison: Average Restored Area');
    
    % Create the Grouped Bar Chart
    b = bar(x_data, y_data, 'grouped');
    
    % --- Customize Appearance for Publication ---
    color1 = [0.8, 0, 0];    % Red for "Our Method"
    color2 = [0, 0.4, 0.8];    % Blue for "HCHA"
    color3 = [0.1, 0.5, 0.1];  % Green for HPA
color4 = [0.6, 0.2, 0.8];  % Purple for HORA

b(1).FaceColor = color1;
b(2).FaceColor = color2;
b(3).FaceColor = color3;
b(4).FaceColor = color4;
    b(1).BarWidth = 0.8;
    
    % --- Add Labels and Legend ---
    xlabel('Number of Deployed Nodes', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Avg. Restored Area per Action', 'FontSize', 12, 'FontWeight', 'bold');    
    % --- CHANGED: Legend Location and Style ---
    % 'northwest', 'northeast', 'southwest', 'southeast' are common choices
    % for placing the legend inside the plot area. 'northeast' often works well.
    % The 'boxoff' command makes the legend's border invisible for a cleaner look.
    lgd = legend(legend_names, 'Location', 'northeast', 'FontSize', 11);
    lgd.Box = 'off';
    
    % --- Final Touches ---
    ax = gca; % Get the current axes handle
    ax.FontSize = 11;
    ax.YGrid = 'on';
    ax.Box = 'on';
    
    % Set Y-axis to start at 0
    ylim([0, max(y_data(:)) * 1.15]);
        % --- ADD THIS BLOCK to add data labels on top of the bars ---
    if ~contains(y_label_text, 'Efficiency') % Only for non-percentage graphs
        for i = 1:length(b)
            xtips = b(i).XEndPoints;
            ytips = b(i).YEndPoints;
            labels = string(round(ytips, 1)); % Round to one decimal place
            text(xtips, ytips, labels, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    
    % --- CHANGED: Make the X-axis tight to the data ---
    % This command adjusts the x-axis limits to fit the data points snugly,
    % removing the extra white space on the right.
    xlim([min(x_data) - 100, max(x_data) + 100]); % Add a small padding
    
    % Save the figure for your paper
    print('barchart_restored_area_final', '-depsc', '-r300');
    fprintf('\nFinal bar chart saved as barchart_restored_area_final.eps\n');
end

function [x, y] = generate_network(area, num_nodes)
    x = area(1) * rand(1, num_nodes);
    y = area(2) * rand(1, num_nodes);
end

% --- START OF CONSOLIDATED HELPER FUNCTIONS ---
% I have copied all necessary functions from your three scripts here.
% Some functions might be duplicated; I have kept one version.

function stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width)
    grid_resolution = 5; % Check every 5 meters
    x_grid = 0:grid_resolution:net_length;
    y_grid = 0:grid_resolution:net_width;
    [X, Y] = meshgrid(x_grid, y_grid);
    coverage_map = zeros(size(X));
    for i = 1:length(x_loc)
        dist_sq = (X - x_loc(i)).^2 + (Y - y_loc(i)).^2;
        coverage_map = coverage_map + (dist_sq <= Rs^2);
    end
    total_points = numel(X);
    hole_points = sum(coverage_map(:) == 0);
    k_covered_points = total_points - hole_points;
    point_area = grid_resolution^2;
    stats.total_area = net_length * net_width;
    stats.hole_area = hole_points * point_area;
    stats.k_covered_area = k_covered_points * point_area;
    stats.hole_percentage = (stats.hole_area / stats.total_area) * 100;
    stats.k_covered_percentage = (stats.k_covered_area / stats.total_area) * 100;
    max_k = max(coverage_map(:));
    stats.specific_k_areas = zeros(1, max_k);
    for k = 1:max_k
        stats.specific_k_areas(k) = sum(coverage_map(:) == k) * point_area;
    end
    stats.specific_k_percentages = (stats.specific_k_areas / stats.total_area) * 100;
end

function T = selectTripletsDelaunay(B_nodes, x_loc, y_loc, Rs, current_coverage_percent)
    if current_coverage_percent < 60.0, distance_threshold = 1.8 * Rs;
    elseif current_coverage_percent < 85.0, distance_threshold = 1.6 * Rs;
    else, distance_threshold = 1.4 * Rs; end
    T = {};
    if length(B_nodes) < 3, return; end
    boundary_x = x_loc(B_nodes); boundary_y = y_loc(B_nodes);
    dt = delaunayTriangulation(boundary_x', boundary_y');
    candidate_triangles = dt.ConnectivityList;
    valid_triplets_info = [];
    for i = 1:size(candidate_triangles, 1)
        idx1 = candidate_triangles(i, 1); idx2 = candidate_triangles(i, 2); idx3 = candidate_triangles(i, 3);
        node_id1 = B_nodes(idx1); node_id2 = B_nodes(idx2); node_id3 = B_nodes(idx3);
        p1 = [x_loc(node_id1), y_loc(node_id1)]; p2 = [x_loc(node_id2), y_loc(node_id2)]; p3 = [x_loc(node_id3), y_loc(node_id3)];
        d12 = norm(p1 - p2); d23 = norm(p2 - p3); d31 = norm(p3 - p1);
        if (d12 > distance_threshold) && (d23 > distance_threshold) && (d31 > distance_threshold)
            current_triplet = [node_id1, node_id2, node_id3];
            [potential_center, ~] = solveApolloniusForTriplet(current_triplet, x_loc, y_loc, Rs);
            if ~isempty(potential_center)
                is_covered = false;
                for k = 1:length(x_loc)
                    if norm(potential_center - [x_loc(k), y_loc(k)]) <= Rs
                        is_covered = true; break;
                    end
                end
                if ~is_covered
                    perimeter = d12 + d23 + d31;
                    valid_triplets_info(end+1, :) = [node_id1, node_id2, node_id3, perimeter];
                end
            end
        end
    end
    if ~isempty(valid_triplets_info)
        sorted_triplets_info = sortrows(valid_triplets_info, 4);
        for i = 1:size(sorted_triplets_info, 1), T{end+1} = sorted_triplets_info(i, 1:3); end
    end
end

function [optimal_center, optimal_radius] = solveApolloniusForTriplet(triplet, x_loc, y_loc, Rs)
    % Placeholder for your actual Apollonius solver logic
    % For now, returns a random-ish point for testing
    c1 = [x_loc(triplet(1)), y_loc(triplet(1))];
    c2 = [x_loc(triplet(2)), y_loc(triplet(2))];
    c3 = [x_loc(triplet(3)), y_loc(triplet(3))];
    optimal_center = (c1 + c2 + c3) / 3; % Centroid as a simple placeholder
    optimal_radius = mean([norm(optimal_center-c1), norm(optimal_center-c2), norm(optimal_center-c3)]) - Rs;
end

function [Pi, CIPL, CIPU] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width)
  Pi = []; CIPL = []; CIPU = []; xi = x_loc(i); yi = y_loc(i);
  if isBorderSensor(xi, yi, net_length, net_width, Rs)
      border_intersections = borderIntersections(xi, yi, Rs, net_length, net_width);
      for p_idx = 1:size(border_intersections, 1)
          p = border_intersections(p_idx, :);
          if ~isCoveredByAnyNeighbor(p, i, Ni{i}, x_loc, y_loc, Rs)
              Pi = [Pi; p];
          end
      end
  end
  for j = Ni{i}
      xj = x_loc(j); yj = y_loc(j);
      inter_pts = circleIntersections(xi, yi, xj, yj, Rs);
      for p_idx = 1:size(inter_pts, 1)
          p = inter_pts(p_idx, :);
          if ~isCoveredByAnotherSpecificNeighbor(p, i, j, Ni{i}, x_loc, y_loc, Rs)
              Pi = [Pi; p];
              status = isAboveLine(p(1), p(2), xi, yi, xj, yj);
              if strcmp(status, 'above'), CIPU = [CIPU; p];
              elseif strcmp(status, 'below'), CIPL = [CIPL; p]; end
          end
      end
  end
  if ~isempty(Pi), Pi = uniquetol(Pi, 1e-9, 'ByRows', true); end
  if ~isempty(CIPL), CIPL = uniquetol(CIPL, 1e-9, 'ByRows', true); end
  if ~isempty(CIPU), CIPU = uniquetol(CIPU, 1e-9, 'ByRows', true); end
end

function hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs)
    hole_centers = [];
    if length(x_loc) < 3, return; end
    dt = delaunayTriangulation(x_loc', y_loc');
    [cc, rc] = circumcenter(dt);
    uncovered_centers = [];
    for i = 1:size(cc, 1)
        if rc(i) > Rs
            potential_hole_center = cc(i, :);
            is_covered = false;
            for k = 1:length(x_loc)
                if norm(potential_hole_center - [x_loc(k), y_loc(k)]) <= Rs
                    is_covered = true; break;
                end
            end
            if ~is_covered, uncovered_centers(end+1, :) = potential_hole_center; end
        end
    end
    if ~isempty(uncovered_centers), hole_centers = uniquetol(uncovered_centers, 1e-6, 'ByRows', true); end
end

function is_hole_boundary_node = findHoleBoundaryNode_Yao(node_X_id, x_loc, y_loc, Rc)
    no_nodes_current = length(x_loc); neighbors_of_X = [];
    for j = 1:no_nodes_current
        if node_X_id ~= j && norm([x_loc(node_X_id)-x_loc(j), y_loc(node_X_id)-y_loc(j)]) < Rc
            neighbors_of_X(end+1) = j;
        end
    end
    if length(neighbors_of_X) < 2, is_hole_boundary_node = true; return; end
    sub_network_nodes_ids = [node_X_id, neighbors_of_X];
    sub_network_coords = [x_loc(sub_network_nodes_ids)', y_loc(sub_network_nodes_ids)'];
    try, k_boundary = boundary(sub_network_coords, 1);
    catch, is_hole_boundary_node = true; return; end
    boundary_node_global_ids = sub_network_nodes_ids(unique(k_boundary));
    is_hole_boundary_node = ismember(node_X_id, boundary_node_global_ids);
end

function [x_loc_new, y_loc_new, nodes_added] = patchHole_Yao(hole_nodes, x_loc, y_loc, Rc, Rs, net_length, net_width)
    x_loc_new = x_loc; y_loc_new = y_loc; nodes_added = 0;
    if length(hole_nodes) < 2, return; end
    hole_coords = [x_loc(hole_nodes)', y_loc(hole_nodes)'];
    try, boundary_path_indices = boundary(hole_coords, 1);
    catch, return; end
    for i = 1:(length(boundary_path_indices)-1)
        nodeA_id = hole_nodes(boundary_path_indices(i));
        nodeB_id = hole_nodes(boundary_path_indices(i+1));
        pA = [x_loc(nodeA_id), y_loc(nodeA_id)];
        pB = [x_loc(nodeB_id), y_loc(nodeB_id)];
        midpoint = (pA + pB) / 2;
        vec_AB = pB - pA;
        perp_vec_normalized = [-vec_AB(2), vec_AB(1)] / norm(vec_AB);
        dist_A_mid = norm(pA - midpoint);
        if Rc^2 < dist_A_mid^2, continue; end
        h = sqrt(Rc^2 - dist_A_mid^2);
        pos_N1 = midpoint + h * perp_vec_normalized;
        pos_N2 = midpoint - h * perp_vec_normalized;
        network_centroid = [mean(x_loc), mean(y_loc)];
        if norm(pos_N1 - network_centroid) > norm(pos_N2 - network_centroid), final_pos_N = pos_N1;
        else, final_pos_N = pos_N2; end
        if (final_pos_N(1) >= 0 && final_pos_N(1) <= net_length && final_pos_N(2) >= 0 && final_pos_N(2) <= net_width)
            x_loc_new(end+1) = final_pos_N(1);
            y_loc_new(end+1) = final_pos_N(2);
            nodes_added = nodes_added + 1;
        end
    end
end

function cluster_labels = my_dbscan(data, epsilon, min_pts)
    num_pts = size(data, 1); cluster_labels = zeros(num_pts, 1); cluster_id = 0;
    for i = 1:num_pts
        if cluster_labels(i) ~= 0, continue; end
        neighbor_indices = find_neighbors(data, i, epsilon);
        if length(neighbor_indices) < min_pts, cluster_labels(i) = -1; continue; end
        cluster_id = cluster_id + 1; cluster_labels(i) = cluster_id;
        q = neighbor_indices'; head = 1;
        while head <= length(q)
            current_pt_idx = q(head); head = head + 1;
            if cluster_labels(current_pt_idx) == -1, cluster_labels(current_pt_idx) = cluster_id; end
            if cluster_labels(current_pt_idx) == 0
                cluster_labels(current_pt_idx) = cluster_id;
                new_neighbor_indices = find_neighbors(data, current_pt_idx, epsilon);
                if length(new_neighbor_indices) >= min_pts, q = [q, new_neighbor_indices']; end
            end
        end
    end
end

function neighbor_indices = find_neighbors(data, point_idx, epsilon)
    p1 = data(point_idx, :);
    distances = sqrt(sum((data - p1).^2, 2));
    neighbor_indices = find(distances <= epsilon);
end

function coverage_ratio = estimateCoverageInWSN(center, radius, net_length, net_width)
    num_samples = 10000; 
    x_min = center(1) - radius; x_max = center(1) + radius;
    y_min = center(2) - radius; y_max = center(2) + radius;
    rand_x = x_min + (x_max - x_min) * rand(num_samples, 1);
    rand_y = y_min + (y_max - y_min) * rand(num_samples, 1);
    dist_sq = (rand_x - center(1)).^2 + (rand_y - center(2)).^2;
    is_in_circle_mask = dist_sq <= radius^2;
    points_in_circle = sum(is_in_circle_mask);
    if points_in_circle == 0, coverage_ratio = 0; return; end
    in_circle_x = rand_x(is_in_circle_mask);
    in_circle_y = rand_y(is_in_circle_mask);
    is_in_wsn_mask = in_circle_x >= 0 & in_circle_x <= net_length & in_circle_y >= 0 & in_circle_y <= net_width;
    points_in_wsn = sum(is_in_wsn_mask);
    coverage_ratio = points_in_wsn / points_in_circle;
end

% --- Other minor helpers from tutorial1.m
function isBorder = isBorderSensor(x, y, net_length, net_width, Rs), isBorder = (x - Rs <= 0) || (x + Rs >= net_length) || (y - Rs <= 0) || (y + Rs >= net_width); end
function points = circleIntersections(x1, y1, x2, y2, r)
  d = sqrt((x2 - x1)^2 + (y2 - y1)^2);
  if d >= 2*r || d == 0, points = []; return; end
  a = d^2 / (2*d); h = sqrt(r^2 - a^2);
  xm = x1 + a * (x2 - x1) / d; ym = y1 + a * (y2 - y1) / d;
  rx = -(y2 - y1) * (h / d); ry = (x1 - x2) * (h / d);
  points = [xm + rx, ym + ry; xm - rx, ym - ry];
end
function covered = isCoveredByPoint(point, sensor_x, sensor_y, sensor_r), d = sqrt((point(1) - sensor_x)^2 + (point(2) - sensor_y)^2); covered = (d <= sensor_r); end
function covered = isCoveredByAnyNeighbor(point, current_node_idx, neighbors_list, x_loc, y_loc, Rs)
   covered = false;
  for k_val = neighbors_list
      if (k_val ~= current_node_idx)
           if isCoveredByPoint(point, x_loc(k_val), y_loc(k_val), Rs), covered = true; return; end
      end
  end
end
function covered = isCoveredByAnotherSpecificNeighbor(point, Si_idx, Sj_idx, neighbors_list, x_loc, y_loc, Rs)
  covered = false;
  for k_val = neighbors_list
      if (k_val ~= Si_idx) && (k_val ~= Sj_idx)
          if isCoveredByPoint(point, x_loc(k_val), y_loc(k_val), Rs), covered = true; return; end
      end
  end
end
function status = isAboveLine(point_x, point_y, x1, y1, x2, y2)
  val = (x2 - x1)*(point_y - y1) - (y2 - y1)*(point_x - x1);
  if abs(val) < 1e-9, status = 'on';
  elseif val > 0, status = 'above';
  else, status = 'below'; end
end


% --- MISSING HELPER FUNCTION FROM tutorial1.m ---

function pts = borderIntersections(x, y, r, net_length, net_width)
  pts = [];
   % Intersection with Left Boundary (x_line = 0)
   if (x - r <= 0) && (x + r >= 0)
      delta_x = abs(x - 0);
      if delta_x <= r
          dy = sqrt(r^2 - delta_x^2);
          p1 = [0, y - dy];
          p2 = [0, y + dy];
          if p1(2) >= 0 && p1(2) <= net_width, pts = [pts; p1]; end
          if dy ~= 0 && p2(2) >= 0 && p2(2) <= net_width, pts = [pts; p2]; end
      end
  end
  % Intersection with Right Boundary (x_line = net_length)
  if (x + r >= net_length) && (x - r <= net_length)
      delta_x = abs(x - net_length);
      if delta_x <= r
          dy = sqrt(r^2 - delta_x^2);
          p1 = [net_length, y - dy];
          p2 = [net_length, y + dy];
          if p1(2) >= 0 && p1(2) <= net_width, pts = [pts; p1]; end
          if dy ~= 0 && p2(2) >= 0 && p2(2) <= net_width, pts = [pts; p2]; end
      end
  end
   % Intersection with Bottom Boundary (y_line = 0)
  if (y - r <= 0) && (y + r >= 0)
      delta_y = abs(y - 0);
      if delta_y <= r
          dx = sqrt(r^2 - delta_y^2);
          p1 = [x - dx, 0];
          p2 = [x + dx, 0];
          if p1(1) >= 0 && p1(1) <= net_length, pts = [pts; p1]; end
          if dx ~= 0 && p2(1) >= 0 && p2(1) <= net_length, pts = [pts; p2]; end
      end
  end
  % Intersection with Top Boundary (y_line = net_width)
  if (y + r >= net_width) && (y - r <= net_width)
      delta_y = abs(y - net_width);
      if delta_y <= r
          dx = sqrt(r^2 - delta_y^2);
          p1 = [x - dx, net_width];
          p2 = [x + dx, net_width];
          if p1(1) >= 0 && p1(1) <= net_length, pts = [pts; p1]; end
          if dx ~= 0 && p2(1) >= 0 && p2(1) <= net_length, pts = [pts; p2]; end
      end
  end
  if ~isempty(pts)
      pts = uniquetol(pts, 1e-9, 'ByRows', true);
      pts = pts(pts(:,1) >= 0 & pts(:,1) <= net_length & pts(:,2) >= 0 & pts(:,2) <= net_width, :);
  end
end


function [restored_area, nodes_added_total, avg_efficiency] = run_nlchr_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    stats_before_full_run = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    all_efficiencies = []; max_area_per_node = pi * Rs^2;
    TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 10; iteration_count = 0;

    while true
        iteration_count = iteration_count + 1;
        stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
        if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end

        % --- CRITICAL CHANGE: Use the new hole detection method ---
        triangular_holes = find_triangular_holes(x_loc, y_loc, Rs);
        if isempty(triangular_holes), break; end
        
        all_new_node_positions = [];
        for i = 1:length(triangular_holes)
            triplet_nodes = triangular_holes{i};
            
            % For each triangular hole, find the two closest nodes.
            % This gives NLCHR the best possible chance to succeed.
            p1_id = triplet_nodes(1); p2_id = triplet_nodes(2); p3_id = triplet_nodes(3);
            p1 = [x_loc(p1_id), y_loc(p1_id)];
            p2 = [x_loc(p2_id), y_loc(p2_id)];
            p3 = [x_loc(p3_id), y_loc(p3_id)];

            dists = [norm(p1-p2), norm(p2-p3), norm(p3-p1)];
            [min_dist, min_idx] = min(dists);

            if min_dist <= 2 * Rs
                if min_idx == 1, pair = [p1; p2];
                elseif min_idx == 2, pair = [p2; p3];
                else, pair = [p3; p1]; end
                
                hole_centroid = mean([p1;p2;p3], 1);
                new_vertex = calculate_isosceles_vertex(pair(1,:), pair(2,:), Rs, hole_centroid);
                
                if ~isempty(new_vertex)
                    all_new_node_positions = [all_new_node_positions; new_vertex];
                end
            end
        end
        
        if isempty(all_new_node_positions), break; end
        
        final_positions = resolve_nlchr_conflicts(all_new_node_positions, Rs);
        if isempty(final_positions), break; end
        
        for n_idx = 1:size(final_positions, 1)
            pos = final_positions(n_idx, :);
            stats_after_one_node = calculateKCoverage([x_loc; pos(1)], [y_loc; pos(2)], Rs, net_length, net_width);
            area_restored_by_this_node = max(0, stats_current.hole_area - stats_after_one_node.hole_area);
            all_efficiencies(end+1) = (area_restored_by_this_node / max_area_per_node) * 100;
        end

        x_loc = [x_loc; final_positions(:,1)];
        y_loc = [y_loc; final_positions(:,2)];
    end
    
    stats_after_full_run = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before_full_run.hole_area - stats_after_full_run.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
    if isempty(all_efficiencies), avg_efficiency = 0; else, avg_efficiency = mean(all_efficiencies, 'omitnan'); end
end

% --- ADD THIS TOOLBOX-FREE HELPER FUNCTION INSIDE THE WRAPPER ---
function D = manual_pdist(X)
    n = size(X, 1);
    D = zeros(1, n*(n-1)/2);
    k = 1;
    for ii = 1:n-1
        for jj = ii+1:n
            D(k) = norm(X(ii,:) - X(jj,:));
            k = k+1;
        end
    end
end

% --- REVISED ISOSCELES VERTEX HELPER ---
% --- REVISED ISOSCELES VERTEX HELPER ---
function new_vertex = calculate_isosceles_vertex(p1, p2, side_length, hole_centroid)
    d = norm(p1 - p2);
    if d > 2 * side_length, new_vertex = []; return; end
    
    midpoint = (p1 + p2) / 2;
    height = sqrt(max(0, side_length^2 - (d/2)^2));
    perp_vector = [p2(2) - p1(2), p1(1) - p2(1)];
    if norm(perp_vector) == 0, new_vertex = []; return; end
    
    unit_perp_vector = perp_vector / norm(perp_vector);
    v1 = midpoint + height * unit_perp_vector;
    v2 = midpoint - height * unit_perp_vector;
    
    % IMPROVEMENT: Choose the vertex closer to the hole's centroid.
    % This ensures the new node is placed INWARD, filling the hole.
    if norm(v1 - hole_centroid) < norm(v2 - hole_centroid)
        new_vertex = v1;
    else
        new_vertex = v2;
    end
end



% --- NEW HELPER FUNCTION SPECIFICALLY FOR A ROBUST NLCHR ---
function triangular_holes = find_triangular_holes(x, y, R)
    % This function finds holes by identifying Delaunay triangles that
    % represent voids, and returns the node IDs of those triangles.
    triangular_holes = {};
    if length(x) < 3, return; end

    % Use delaunayTriangulation on ALL nodes to find potential voids
    try
        dt = delaunayTriangulation(x(:), y(:));
    catch
        % Fails if all points are collinear, which is very rare.
        return;
    end
    
    [cc, rc] = circumcenter(dt);
    
    for i = 1:size(cc, 1)
        % Condition: A triangle represents a hole if its circumradius is
        % larger than the sensor's sensing radius.
        if rc(i) > R
            potential_hole_center = cc(i, :);
            
            % Crucial check: Is this potential hole's center already covered?
            if ~isCoveredByAny(potential_hole_center, x, y, R)
                % If not covered, this is a valid triangular hole.
                % Store the node IDs of the three vertices of this triangle.
                triangular_holes{end+1} = dt.ConnectivityList(i, :);
            end
        end
    end
end

function c = isCoveredByAny(p, x, y, R)
    % Calculates squared distances from the point 'p' to all nodes.
    % This version correctly handles column vectors for x and y.
    dist_sq = (x(:) - p(1)).^2 + (y(:) - p(2)).^2;
    % Checks if any squared distance is less than or equal to Rs^2.
    c = any(dist_sq <= R^2);
end

function final_pos=resolve_nlchr_conflicts(pos,R)
    if size(pos,1)>1
        % This is the corrected, TOOLBOX-FREE version
        dist_matrix = zeros(size(pos,1));
        for i = 1:size(pos,1)
            for j = i+1:size(pos,1)
                dist_matrix(i,j) = norm(pos(i,:) - pos(j,:));
                dist_matrix(j,i) = dist_matrix(i,j);
            end
        end
        dist_matrix(logical(eye(size(dist_matrix))))=inf;
        keep=true(size(pos,1),1);
        for j=1:size(pos,1)
            if keep(j)
                close_pts=dist_matrix(j,:)<R;
                keep(close_pts)=false;
                keep(j)=true;
            end
        end
        final_pos=pos(keep,:);
    else
        final_pos=pos;
    end
end
end