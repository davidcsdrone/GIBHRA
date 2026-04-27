clc; clear; close all;
rng('default');

% --- Initial Setup ---
no_nodes = input('Enter the number of nodes: ');
net_length = input('Enter the Length of the network: ');
net_width = input('Enter the width of the network: ');

% --- Simulation Parameters ---
TARGET_COVERAGE_PERCENT = 98.0;
MAX_ITERATIONS = 20; % Increased slightly to allow for more complex scenarios
Rs = 20; % Sensing range radius
Rc = 2 * Rs; % Communication range

% --- Sanity Check for Yao et al. Assumption ---
if Rs < Rc / sqrt(3)
    warning('Rs < Rc/sqrt(3) condition from Yao et al. paper is NOT met. Results may be inaccurate.');
end

% --- Initialization ---
x_loc = zeros(1, no_nodes);
y_loc = zeros(1, no_nodes);

% --- Plotting Setup ---
figure;
hold on;
grid on;
axis equal;
xlim([0 net_length]);
ylim([0 net_width]);
xlabel('Network Length');
ylabel('Network Width');
plot([0 net_length net_length 0 0], [0 0 net_width net_width 0], 'k-', 'LineWidth', 1.5);

% --- Initial Node Deployment ---
for i = 1:no_nodes
  x_loc(i) = net_length * rand;
  y_loc(i) = net_width * rand;
  theta = linspace(0, 2*pi, 100);
  circle_x = x_loc(i) + Rs * cos(theta);
  circle_y = y_loc(i) + Rs * sin(theta);
  fill(circle_x, circle_y, 'g', 'FaceAlpha', 0.1, 'EdgeColor','b', 'LineWidth',1);
  plot(x_loc(i), y_loc(i), 'k.', 'MarkerSize', 15);
  pause(0.01);
end

disp('-----Calculating k-coverage before hole restoration-----');
stats_before = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Initial Hole Area: %.2fm^2 (%.2f%%)\n', stats_before.hole_area, stats_before.hole_percentage);
fprintf('Initial Covered Area (k>=1): %.2f m^2 (%.2f%%)\n', stats_before.k_covered_area, stats_before.k_covered_percentage);

% =========================================================================
% --- START OF THE ITERATIVE RESTORATION LOOP ---
% =========================================================================
iteration_count = 0;
while true
   iteration_count = iteration_count + 1;
   fprintf('\n\n====================================================\n');
   fprintf('--- Starting Restoration Iteration #%d ---\n', iteration_count);
  
   no_nodes_current = length(x_loc);
   fprintf('Current node count: %d\n', no_nodes_current);
   
   stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
   current_coverage_percent = stats_current.k_covered_percentage;
   fprintf('Current Coverage (k>=1): %.2f%%\n', current_coverage_percent);
   
   if current_coverage_percent >= TARGET_COVERAGE_PERCENT
       fprintf('\nTarget coverage of %.2f%% reached. Stopping iterations.\n', TARGET_COVERAGE_PERCENT);
       break;
   end
  
   if iteration_count > MAX_ITERATIONS
       fprintf('\nMaximum number of iterations (%d) reached. Stopping.\n', MAX_ITERATIONS);
       break;
   end
  
   disp('--- Implementing Yao et al. Hole Detection ---');
  
   all_hole_boundary_nodes = [];
   for i = 1:no_nodes_current
       if findHoleBoundaryNode_Yao(i, x_loc, y_loc, Rc)
           all_hole_boundary_nodes(end+1) = i;
       end
   end
   
   if isempty(all_hole_boundary_nodes)
       fprintf('No more hole-boundary-nodes detected. Stopping iterations.\n');
       break;
   end
   
   fprintf('Found %d hole-boundary-nodes.\n', length(all_hole_boundary_nodes));
   plot(x_loc(all_hole_boundary_nodes), y_loc(all_hole_boundary_nodes), 'yo', 'MarkerSize', 10, 'LineWidth', 2);
   
   fprintf('Grouping %d boundary nodes into proximity-based clusters...\n', length(all_hole_boundary_nodes));
   boundary_coords = [x_loc(all_hole_boundary_nodes)', y_loc(all_hole_boundary_nodes)'];
   epsilon = Rc; 
   minpts = 3;
   
   cluster_indices = my_dbscan(boundary_coords, epsilon, minpts);
   
   num_clusters = max(cluster_indices);
   if num_clusters < 1
       fprintf('Boundary nodes are too sparse to form any clusters. Stopping.\n');
       break;
   end
   
   holes = cell(1, num_clusters);
   for i = 1:num_clusters
       nodes_in_this_cluster = all_hole_boundary_nodes(cluster_indices == i);
       holes{i} = nodes_in_this_cluster;
   end
   fprintf('Grouped boundary nodes into %d distinct hole cluster(s).\n', num_clusters);

   if isempty(holes)
       fprintf('No valid hole clusters found. Stopping.\n');
       break;
   end
   
   hole_areas = zeros(1, length(holes));
   for i = 1:length(holes)
       current_hole_nodes = holes{i};
       if length(current_hole_nodes) >= 3
           hole_coords = [x_loc(current_hole_nodes)', y_loc(current_hole_nodes)'];
           [~, current_area] = convhull(hole_coords);
           hole_areas(i) = current_area;
       else
           hole_areas(i) = 0;
       end
   end
   
   [~, largest_area_index] = max(hole_areas);
   
   if hole_areas(largest_area_index) == 0
       fprintf('No valid holes with a measurable area found to patch. Stopping.\n');
       break;
   end

   hole_to_patch = holes{largest_area_index};
   fprintf('Attempting to patch the largest area hole (Index %d) with area %.2f.\n', ...
           largest_area_index, hole_areas(largest_area_index));
           
   [x_loc, y_loc, nodes_deployed_this_iteration] = patchHole_Yao(hole_to_patch, x_loc, y_loc, Rc, Rs, net_length, net_width);
   
   if nodes_deployed_this_iteration == 0
       fprintf('Could not deploy any new valid nodes for the selected hole. Stopping.\n');
       break;
   end
   
   drawnow;
end

% --- FINAL SUMMARY REPORT ---
fprintf('\n\n====================================================\n');
fprintf('--- FINAL RESTORATION SUMMARY ---\n');
final_iteration_count = max(0, iteration_count - 1);
fprintf('Algorithm stopped after %d restoration iteration(s).\n', final_iteration_count);
total_nodes_added = length(x_loc) - no_nodes;
fprintf('Total initial nodes: %d\n', no_nodes);
fprintf('Total new nodes added: %d\n', total_nodes_added);
fprintf('Final total nodes: %d\n', length(x_loc));
final_stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Final Coverage (k>=1): %.2f%%\n', final_stats.k_covered_percentage);
fprintf('====================================================\n');
hold off;

% =========================================================================
% --- ALGORITHM HELPER FUNCTIONS ---
% =========================================================================

function is_hole_boundary_node = findHoleBoundaryNode_Yao(node_X_id, x_loc, y_loc, Rc)
    no_nodes_current = length(x_loc);
    neighbors_of_X = [];
    for j = 1:no_nodes_current
        if node_X_id ~= j && norm([x_loc(node_X_id)-x_loc(j), y_loc(node_X_id)-y_loc(j)]) < Rc
            neighbors_of_X(end+1) = j;
        end
    end
    if length(neighbors_of_X) < 2, is_hole_boundary_node = true; return; end
    sub_network_nodes_ids = [node_X_id, neighbors_of_X];
    sub_network_coords = [x_loc(sub_network_nodes_ids)', y_loc(sub_network_nodes_ids)'];
    try
        k_boundary = boundary(sub_network_coords, 1);
    catch
        is_hole_boundary_node = true; return;
    end
    boundary_node_global_ids = sub_network_nodes_ids(unique(k_boundary));
    is_hole_boundary_node = ismember(node_X_id, boundary_node_global_ids);
end

function [x_loc_new, y_loc_new, nodes_added] = patchHole_Yao(hole_nodes, x_loc, y_loc, Rc, Rs, net_length, net_width)
    x_loc_new = x_loc;
    y_loc_new = y_loc;
    nodes_added = 0;
    if length(hole_nodes) < 2, return; end
    hole_coords = [x_loc(hole_nodes)', y_loc(hole_nodes)'];
    try
        boundary_path_indices = boundary(hole_coords, 1);
    catch
        return;
    end
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
        if norm(pos_N1 - network_centroid) > norm(pos_N2 - network_centroid)
            final_pos_N = pos_N1;
        else
            final_pos_N = pos_N2;
        end
        if (final_pos_N(1) >= 0 && final_pos_N(1) <= net_length && ...
            final_pos_N(2) >= 0 && final_pos_N(2) <= net_width)
            x_loc_new(end+1) = final_pos_N(1);
            y_loc_new(end+1) = final_pos_N(2);
            nodes_added = nodes_added + 1;
            plot(final_pos_N(1), final_pos_N(2), 'bs', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
            theta = linspace(0, 2*pi, 100);
            fill(final_pos_N(1) + Rs*cos(theta), final_pos_N(2) + Rs*sin(theta), 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'b', 'LineStyle', '--');
        else
            fprintf('Skipping node deployment at (%.2f, %.2f): Outside WSN boundary.\n', final_pos_N(1), final_pos_N(2));
        end
    end
end

function cluster_labels = my_dbscan(data, epsilon, min_pts)
    num_pts = size(data, 1);
    cluster_labels = zeros(num_pts, 1);
    cluster_id = 0;
    for i = 1:num_pts
        if cluster_labels(i) ~= 0, continue; end
        neighbor_indices = find_neighbors(data, i, epsilon);
        if length(neighbor_indices) < min_pts
            cluster_labels(i) = -1; % Noise
            continue;
        end
        cluster_id = cluster_id + 1;
        cluster_labels(i) = cluster_id;
        
        q = neighbor_indices'; % *** INITIALIZE QUEUE AS ROW VECTOR ***
        
        head = 1;
        while head <= length(q)
            current_pt_idx = q(head);
            head = head + 1;
            
            if cluster_labels(current_pt_idx) == -1
                cluster_labels(current_pt_idx) = cluster_id;
            end
            
            if cluster_labels(current_pt_idx) == 0
                cluster_labels(current_pt_idx) = cluster_id;
                new_neighbor_indices = find_neighbors(data, current_pt_idx, epsilon);
                if length(new_neighbor_indices) >= min_pts
                    q = [q, new_neighbor_indices']; % *** APPEND ROW VECTOR TO ROW VECTOR ***
                end
            end
        end
    end
end

function neighbor_indices = find_neighbors(data, point_idx, epsilon)
    p1 = data(point_idx, :);
    distances = sqrt(sum((data - p1).^2, 2));
    neighbor_indices = find(distances <= epsilon);
end