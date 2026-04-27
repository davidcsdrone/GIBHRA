% =========================================================================
% FUNCTION: restore_network_with_HORA.m
% DESCRIPTION: Implements the main iterative loop of the HORA algorithm.
% =========================================================================
function [x_final, y_final] = restore_network_with_HORA(x_loc, y_loc, Rs, Rc, area_dims)

    x_current = x_loc;
    y_current = y_loc;
    
    max_rounds = 5; % Limit the number of healing rounds
    for round = 1:max_rounds
        fprintf('  Round %d of %d...\n', round, max_rounds);
        
        % --- Phase 1: Identify Mobile Nodes (CT/HCT) ---
        mobile_node_ids = identifyMobileNodes_HORA(x_current, y_current, Rc);
        if isempty(mobile_node_ids)
            fprintf('    No mobile (redundant) nodes found. Stopping.\n');
            break;
        end
        
        % --- Phase 2: Find Holes and their Boundary Nodes ---
        all_holes = findAndClusterHoles_HORA(x_current, y_current, Rs, Rc);
        if isempty(all_holes)
            fprintf('    No coverage holes detected. Stopping.\n');
            break;
        end
        
        fprintf('    Found %d mobile nodes and %d holes.\n', length(mobile_node_ids), length(all_holes));
        
        nodes_moved_this_round = 0;
        
        % --- Phase 3 & 4: Select Best Mover and Calculate Target ---
        for i = 1:length(all_holes)
            boundary_nodes = all_holes{i};
            
            % Find the geometric center of the hole
            hole_center = mean([x_current(boundary_nodes)', y_current(boundary_nodes)'], 1);
            
            % Find the mobile node closest to this hole
            dists = sqrt(sum(( [x_current(mobile_node_ids)', y_current(mobile_node_ids)'] - hole_center ).^2, 2));
            [~, min_idx] = min(dists);
            selected_mobile_node_id = mobile_node_ids(min_idx);
            
            % Calculate the target destination (centroid of the hole polygon)
            pgon = polyshape(x_current(boundary_nodes), y_current(boundary_nodes));
            [target_x, target_y] = centroid(pgon);
            
            % --- Phase 5: Verification Check ---
            % Calculate Kh-value at current and potential new location
            kh_current = calculateKhValue_HORA(selected_mobile_node_id, x_current, y_current, Rs);
            
            % Temporarily move the node to calculate its potential new Kh-value
            x_temp = x_current;
            y_temp = y_current;
            x_temp(selected_mobile_node_id) = target_x;
            y_temp(selected_mobile_node_id) = target_y;
            kh_new = calculateKhValue_HORA(selected_mobile_node_id, x_temp, y_temp, Rs);
            
            fprintf('    Mobile node %d: Current Kh=%.2f, Potential New Kh=%.2f\n', selected_mobile_node_id, kh_current, kh_new);
            
            % Only move if the new location is less dense
            if kh_new < kh_current
                fprintf('    --> Moving node %d to patch hole %d.\n', selected_mobile_node_id, i);
                x_current(selected_mobile_node_id) = target_x;
                y_current(selected_mobile_node_id) = target_y;
                nodes_moved_this_round = nodes_moved_this_round + 1;
                
                % Remove the moved node from the list of available mobile nodes for this round
                mobile_node_ids(min_idx) = [];
                if isempty(mobile_node_ids), break; end
            else
                fprintf('    --> Move for node %d rejected (new location is not less dense).\n', selected_mobile_node_id);
            end
        end
        
        if nodes_moved_this_round == 0
            fprintf('  No beneficial moves found this round. Halting.\n');
            break;
        end
    end
    
    x_final = x_current;
    y_final = y_current;
end

function mobile_nodes = identifyMobileNodes_HORA(x, y, Rc)
    num_nodes = length(x);
    neighbors = cell(1, num_nodes);
    for i = 1:num_nodes
        dists = sqrt((x - x(i)).^2 + (y - y(i)).^2);
        neighbors{i} = find(dists < Rc & dists > 0);
    end
    
    mobile_nodes = [];
    for i = 1:num_nodes
        if length(neighbors{i}) < 3
            continue; % Cannot form a 4-node group
        end
        
        % Get all combinations of 3 neighbors
        combos = nchoosek(neighbors{i}, 3);
        
        for j = 1:size(combos, 1)
            group = [i, combos(j,:)]; % The 4-node group
            
            % Check if it's a complete graph (all are neighbors of each other)
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
                break; % Found a CT/HCT, so node 'i' is mobile. Move to next node.
            end
        end
    end
    mobile_nodes = unique(mobile_nodes);
end


% =========================================================================
% FUNCTION: calculateKhValue_HORA.m
% DESCRIPTION: Calculates the Kh-value for a given node. This is a complex
%              geometric task. This implementation provides a robust approximation.
% =========================================================================
function kh_value = calculateKhValue_HORA(node_id, x, y, Rs)
    
    node_pos = [x(node_id), y(node_id)];
    
    % Find all neighbors of the node at its current/potential position
    dists = sqrt((x - node_pos(1)).^2 + (y - node_pos(2)).^2);
    neighbor_ids = find(dists < 2*Rs & dists > 0);
    
    if isempty(neighbor_ids)
        kh_value = 1; % The node covers itself
        return;
    end
    
    % The relevant nodes for checking overlap are the node itself and its neighbors
    group_ids = [node_id, neighbor_ids];
    
    % --- Kh-Value Approximation Strategy ---
    % 1. Find all intersection points between all pairs of circles in the group.
    % 2. For each intersection point, count how many circles in the group contain it.
    % 3. The Kh-value is the maximum count found.
    
    all_intersection_points = [];
    circle_pairs = nchoosek(group_ids, 2);
    
    for i = 1:size(circle_pairs, 1)
        id1 = circle_pairs(i,1);
        id2 = circle_pairs(i,2);
        
        d = norm([x(id1)-x(id2), y(id1)-y(id2)]);
        if d < 2*Rs && d > 0 % Circles intersect
             a = (Rs^2 - Rs^2 + d^2) / (2*d);
             h = sqrt(max(0, Rs^2 - a^2));
             x_mid = x(id1) + a*(x(id2)-x(id1))/d;
             y_mid = y(id1) + a*(y(id2)-y(id1))/d;
             
             p1 = [x_mid + h*(y(id2)-y(id1))/d, y_mid - h*(x(id2)-x(id1))/d];
             p2 = [x_mid - h*(y(id2)-y(id1))/d, y_mid + h*(x(id2)-x(id1))/d];
             
             all_intersection_points = [all_intersection_points; p1; p2];
        end
    end
    
    if isempty(all_intersection_points)
        % If no intersections, the max overlap is likely at the node's center
        % Count how many neighbors cover the node's own center
        kh_value = 1 + sum(sqrt((x(neighbor_ids) - node_pos(1)).^2 + (y(neighbor_ids) - node_pos(2)).^2) < Rs);
        return;
    end
    
    % Check each intersection point
    max_coverage_count = 0;
    for i = 1:size(all_intersection_points, 1)
        point = all_intersection_points(i,:);
        
        % Count how many circles in the group cover this point
        coverage_count = sum(sqrt((x(group_ids) - point(1)).^2 + (y(group_ids) - point(2)).^2) <= Rs);
        if coverage_count > max_coverage_count
            max_coverage_count = coverage_count;
        end
    end
    
    kh_value = max_coverage_count;
end

function all_holes = findAndClusterHoles_HORA(x, y, Rs, Rc)
    num_nodes = length(x);
    neighbors = cell(1, num_nodes);
    for i = 1:num_nodes
        dists = sqrt((x - x(i)).^2 + (y - y(i)).^2);
        neighbors{i} = find(dists < 2*Rs & dists > 0);
    end
    
    % --- Find all boundary nodes using a simplified method ---
    % A node is a boundary node if its neighbors don't form a closed loop around it
    boundary_nodes = [];
    for i = 1:num_nodes
        if isempty(neighbors{i}), boundary_nodes(end+1) = i; continue; end
        
        neighbor_coords = [x(neighbors{i})', y(neighbors{i})'];
        try
            hull_indices = convhull(neighbor_coords);
            if length(hull_indices) < length(neighbors{i}) + 1
                 boundary_nodes(end+1) = i;
            end
        catch
            boundary_nodes(end+1) = i;
        end
    end
    
    if isempty(boundary_nodes)
        all_holes = {};
        return;
    end
    
    % --- Cluster boundary nodes to identify distinct holes ---
    boundary_coords = [x(boundary_nodes)', y(boundary_nodes)'];
    % Use DBSCAN for clustering. e.g., points within a comms radius belong to the same hole.
    cluster_indices = dbscan(boundary_coords, Rc, 2); % (data, epsilon, min_pts)
    
    num_clusters = max(cluster_indices);
    all_holes = cell(1, num_clusters);
    for i = 1:num_clusters
        all_holes{i} = boundary_nodes(cluster_indices == i);
    end
end