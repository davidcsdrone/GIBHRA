classdef NLCHR_Algorithm
    methods (Static)
        % =========================================================================
        % --- Main Wrapper for the "Restored Area" Simulation ---
        % =========================================================================
       % --- REVISED NLCHR WRAPPER ---
function [restored_area, nodes_added_total] = run_nlchr_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 10; iteration_count = 0;

    while true
        iteration_count = iteration_count + 1;
        stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
        if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end

        all_holes = findAndClusterHoles_NLCHR(x_loc, y_loc, Rs, area_dims);
        if isempty(all_holes), break; end
        
        all_new_node_positions = [];
        for i = 1:length(all_holes)
            boundary_nodes = all_holes{i};
            if length(boundary_nodes) < 2, continue; end
            
            % --- IMPROVEMENT: Find the CLOSEST pair of nodes in the hole boundary ---
            boundary_coords = [x_loc(boundary_nodes), y_loc(boundary_nodes)];
            dist_matrix = squareform(pdist(boundary_coords));
            dist_matrix(logical(eye(size(dist_matrix)))) = inf;
            [min_dist, linear_idx] = min(dist_matrix(:));
            
            if min_dist <= 2 * Rs
                [row_idx, col_idx] = ind2sub(size(dist_matrix), linear_idx);
                node1_id = boundary_nodes(row_idx);
                node2_id = boundary_nodes(col_idx);
                
                p1 = [x_loc(node1_id), y_loc(node1_id)];
                p2 = [x_loc(node2_id), y_loc(node2_id)];
                
                % IMPROVEMENT: Pass the hole's centroid to choose the inward-facing vertex
                hole_centroid = mean(boundary_coords, 1);
                new_vertex = calculate_isosceles_vertex(p1, p2, Rs, hole_centroid);
                
                if ~isempty(new_vertex), all_new_node_positions = [all_new_node_positions; new_vertex]; end
            end
        end
        if isempty(all_new_node_positions), break; end
        
        final_positions = resolve_nlchr_conflicts(all_new_node_positions, Rs);
        if isempty(final_positions), break; end
        
        x_loc = [x_loc; final_positions(:,1)];
        y_loc = [y_loc; final_positions(:,2)];
    end
    
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end
    end
end

% =========================================================================
function all_holes = findAndClusterHoles_NLCHR(x, y, Rs, area_dims)
    num_nodes = length(x);
    net_length = area_dims(1); net_width = area_dims(2);
    Rc = 2 * Rs;
    Ni = cell(num_nodes, 1);
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j && norm([x(i)-x(j), y(i)-y(j)]) < 2 * Rs
                Ni{i} = [Ni{i}, j];
            end
        end
    end
    boundary_nodes = [];
    for i = 1:num_nodes
        [Pi, ~, ~] = computeCP_and_CIP(i, x, y, Ni, Rs, net_length, net_width);
        if isempty(Ni{i}) || ~isempty(Pi)
            boundary_nodes(end+1) = i;
        end
    end
    if isempty(boundary_nodes), all_holes = {}; return; end
    boundary_coords = [x(boundary_nodes)', y(boundary_nodes)'];
    cluster_indices = cluster_nodes_by_distance(boundary_coords, Rc);
    num_clusters = max(cluster_indices);
    if isempty(num_clusters) || num_clusters == 0, all_holes = {}; return; end
    all_holes = cell(1, num_clusters);
    for i = 1:num_clusters
        all_holes{i} = boundary_nodes(cluster_indices == i);
    end
end