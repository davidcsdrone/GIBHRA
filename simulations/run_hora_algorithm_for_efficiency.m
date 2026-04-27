% =========================================================================
% FUNCTION: run_hora_algorithm_for_efficiency.m
% DESCRIPTION: Implements the HORA placement strategy for the efficiency comparison.
%              This file contains all necessary local helper functions for HORA.
% =========================================================================
function avg_efficiency = run_hora_algorithm_for_efficiency(x_loc_init, y_loc_init, Rs, area_dims)
    % --- Setup ---
    net_length = area_dims(1); 
    net_width = area_dims(2);
    Rc = 2 * Rs;
    max_area_per_node = pi * Rs^2;
    
    % --- Get initial network state ---
    stats_before_iter = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    
    % --- HORA's core logic for finding a target location ---
    % 1. Find and cluster all holes in the network
    all_holes = findAndClusterHoles_HORA(x_loc_init, y_loc_init, Rs, Rc);
    
    if isempty(all_holes)
        avg_efficiency = NaN;
        return;
    end
    
    % 2. Find the largest hole to provide a representative test
    largest_hole_area = 0;
    largest_hole_idx = 0;
    for i = 1:length(all_holes)
        boundary_nodes = all_holes{i};
        if length(boundary_nodes) > 2
            try
                pgon = polyshape(x_loc_init(boundary_nodes), y_loc_init(boundary_nodes));
                if pgon.area > largest_hole_area
                    largest_hole_area = pgon.area;
                    largest_hole_idx = i;
                end
            catch
                continue;
            end
        end
    end
    
    if largest_hole_idx == 0
        avg_efficiency = NaN;
        return;
    end
    
    boundary_nodes_of_largest_hole = all_holes{largest_hole_idx};

    % 3. HORA's Strategy: The target is the centroid of the boundary polygon.
    pgon = polyshape(x_loc_init(boundary_nodes_of_largest_hole), y_loc_init(boundary_nodes_of_largest_hole));
    [target_x, target_y] = centroid(pgon);
    
    % --- Calculate efficiency for this single placement ---
    x_loc_temp = [x_loc_init, target_x];
    y_loc_temp = [y_loc_init, target_y];
    stats_after_one_node = calculateKCoverage(x_loc_temp, y_loc_temp, Rs, net_length, net_width);
    
    area_restored_by_this_node = max(0, stats_before_iter.hole_area - stats_after_one_node.hole_area);
    
    efficiency = (area_restored_by_this_node / max_area_per_node) * 100;
    
    avg_efficiency = efficiency;
end

% =========================================================================
% --- LOCAL HELPER FUNCTIONS FOR HORA ---
% =========================================================================

function all_holes = findAndClusterHoles_HORA(x, y, Rs, Rc)
    num_nodes = length(x);
    neighbors = cell(1, num_nodes);
    for i = 1:num_nodes
        dists = sqrt((x - x(i)).^2 + (y - y(i)).^2);
        neighbors{i} = find(dists < 2*Rs & dists > 0);
    end
    
    boundary_nodes = [];
    for i = 1:num_nodes
        if isempty(neighbors{i}), boundary_nodes(end+1) = i; continue; end
        
        neighbor_coords = [x(neighbors{i})', y(neighbors{i})'];
        if size(neighbor_coords, 1) < 3, boundary_nodes(end+1) = i; continue; end
        
        try
            hull_indices = convhull(neighbor_coords);
            if length(unique(hull_indices)) < length(neighbors{i})
                 boundary_nodes(end+1) = i;
            end
        catch
            boundary_nodes(end+1) = i;
        end
    end
    
    if isempty(boundary_nodes), all_holes = {}; return; end
    
    boundary_coords = [x(boundary_nodes)', y(boundary_nodes)'];
    cluster_indices = dbscan(boundary_coords, Rc, 3); 
    
    num_clusters = max(cluster_indices);
    all_holes = cell(1, num_clusters);
    for i = 1:num_clusters
        all_holes{i} = boundary_nodes(cluster_indices == i);
    end
end

function mobile_nodes = identifyMobileNodes_HORA(x, y, Rc)
    % This is not strictly needed for the placement strategy, 
    % but we keep it for completeness of the HORA concept.
    % In the efficiency script, we assume any node can be placed.
    num_nodes = length(x);
    neighbors = cell(1, num_nodes);
    for i = 1:num_nodes
        dists = sqrt((x - x(i)).^2 + (y - y(i)).^2);
        neighbors{i} = find(dists < Rc & dists > 0);
    end
    
    mobile_nodes = [];
    for i = 1:num_nodes
        if length(neighbors{i}) < 3, continue; end
        combos = nchoosek(neighbors{i}, 3);
        for j = 1:size(combos, 1)
            group = [i, combos(j,:)];
            is_complete = true;
            group_pairs = nchoosek(group, 2);
            for k = 1:size(group_pairs, 1)
                p1 = group_pairs(k,1);
                p2 = group_pairs(k,2);
                if ~ismember(p2, neighbors{p1})
                    is_complete = false;
                    break;
                end
            end
            if is_complete
                mobile_nodes(end+1) = i;
                break;
            end
        end
    end
    mobile_nodes = unique(mobile_nodes);
end