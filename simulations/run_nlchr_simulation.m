% =========================================================================
% SCRIPT: run_nlchr_iterative_simulation.m
% DESCRIPTION: An interactive script to visualize the round-by-round 
%              healing process of the NLCHR algorithm.
% =========================================================================

% --- Main Simulation Setup ---
clear; clc; close all;
rng('default'); 

% --- User Input for Simulation Parameters ---
num_nodes_initial = input('Enter the number of initial nodes: ');
net_length = input('Enter the Length of the network: ');
net_width = input('Enter the width of the network: ');
sensing_radius = 50; 
deployment_area = [net_length, net_width];

% --- Network Generation ---
fprintf('Generating initial network with %d nodes...\n', num_nodes_initial);
[x_current, y_current] = generate_network(deployment_area, num_nodes_initial);

% --- Plot Initial Network State ---
figure('Name', 'NLCHR Iterative Healing', 'Position', [100, 100, 1000, 800]);
subplot(2, 2, 1);
plot_network_state(x_current, y_current, sensing_radius, ...
    sprintf('Initial Deployment (%d Nodes)', length(x_current)), deployment_area);

% --- Iterative Healing and Plotting ---
fprintf('--- Starting NLCHR Iterative Hole Restoration ---\n');
max_rounds = 3; % We will show the first 2 rounds + the final state
nodes_added = false;

for round = 1:max_rounds
    fprintf('  Round %d of healing...\n', round);
    
    % Perform ONE round of healing
    [x_new, y_new, num_added] = perform_one_round_of_NLCHR(x_current, y_current, sensing_radius, deployment_area);
    
    if num_added == 0
        fprintf('    No more nodes could be added. Halting.\n');
        break; % Exit if no progress was made
    end
    
    nodes_added = true;
    x_current = x_new;
    y_current = y_new;
    
    % Plot the intermediate state if it's one of the first rounds
    if round <= 2
        subplot(2, 2, round + 1);
        plot_network_state(x_current, y_current, sensing_radius, ...
            sprintf('After Round %d (%d Nodes)', round, length(x_current)), deployment_area);
    end
end

% --- Plot Final Network State ---
fprintf('--- NLCHR Restoration Complete ---\n');
subplot(2, 2, 4);
if ~nodes_added
    % If no nodes were ever added, just copy the initial plot
    plot_network_state(x_current, y_current, sensing_radius, ...
        sprintf('Final State (No Changes, %d Nodes)', length(x_current)), deployment_area);
else
    plot_network_state(x_current, y_current, sensing_radius, ...
        sprintf('Final State (%d Nodes)', length(x_current)), deployment_area);
end


% =========================================================================
% --- CORE NLCHR LOGIC (LOCAL FUNCTIONS) ---
% =========================================================================
function [x_new, y_new, num_added] = perform_one_round_of_NLCHR(x_loc, y_loc, Rs, area_dims)
    x_current = x_loc;
    y_current = y_loc;
    
    all_holes = findAndClusterHoles_NLCHR(x_current, y_current, Rs, area_dims);
    
    if isempty(all_holes)
        x_new = x_current;
        y_new = y_current;
        num_added = 0;
        return;
    end
    
    fprintf('    Found %d hole(s) to patch this round.\n', length(all_holes));
    
    all_new_node_positions = [];
    
    for i = 1:length(all_holes)
        boundary_nodes = all_holes{i};
        if length(boundary_nodes) < 2, continue; end
        
        boundary_coords = [x_current(boundary_nodes)', y_current(boundary_nodes)'];
        hull_indices = convhull(boundary_coords);
        ordered_nodes = boundary_nodes(hull_indices);
        
        heal_node_values = mod(1:length(ordered_nodes), 2);
        heal_base_node_indices = find(heal_node_values == 1);
        
        if isempty(heal_base_node_indices)
            continue;
        end
        
        for j = 1:length(heal_base_node_indices)
            base_node_idx_in_list = heal_base_node_indices(j);
            node1_id = ordered_nodes(base_node_idx_in_list);
            next_node_idx_in_list = mod(base_node_idx_in_list, length(ordered_nodes)) + 1;
            node2_id = ordered_nodes(next_node_idx_in_list);
            
            p1 = [x_current(node1_id), y_current(node1_id)];
            p2 = [x_current(node2_id), y_current(node2_id)];
            
            new_vertex = calculate_isosceles_vertex(p1, p2, Rs);
            
            if ~isempty(new_vertex)
                all_new_node_positions = [all_new_node_positions; new_vertex];
            end
        end
    end
    
    if isempty(all_new_node_positions)
        x_new = x_current;
        y_new = y_current;
        num_added = 0;
        return;
    end
    
    if size(all_new_node_positions, 1) > 1
        num_new_pos = size(all_new_node_positions, 1);
        dist_matrix = zeros(num_new_pos, num_new_pos);
        for p_idx1 = 1:num_new_pos
            for p_idx2 = p_idx1 + 1:num_new_pos
                dist = norm(all_new_node_positions(p_idx1,:) - all_new_node_positions(p_idx2,:));
                dist_matrix(p_idx1, p_idx2) = dist;
                dist_matrix(p_idx2, p_idx1) = dist;
            end
        end
        
        dist_matrix(logical(eye(size(dist_matrix)))) = inf;
        keep_indices = true(size(all_new_node_positions, 1), 1);
        for j = 1:size(all_new_node_positions, 1)
            if keep_indices(j)
                close_points = dist_matrix(j,:) < Rs;
                keep_indices(close_points) = false;
                keep_indices(j) = true;
            end
        end
        final_positions = all_new_node_positions(keep_indices, :);
    else
        final_positions = all_new_node_positions;
    end
    
    if isempty(final_positions)
        x_new = x_current;
        y_new = y_current;
        num_added = 0;
        return;
    end
    
    x_new = [x_current, final_positions(:,1)'];
    y_new = [y_current, final_positions(:,2)'];
    num_added = size(final_positions, 1);
    fprintf('    Deployed %d new nodes.\n', num_added);
end

% =========================================================================
% --- OTHER LOCAL HELPER FUNCTIONS ---
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

function new_vertex = calculate_isosceles_vertex(p1, p2, side_length)
    d = norm(p1 - p2);
    if d > 2 * side_length, new_vertex = []; return; end
    midpoint = (p1 + p2) / 2;
    height = sqrt(max(0, side_length^2 - (d/2)^2));
    perp_vector = [p2(2) - p1(2), p1(1) - p2(1)];
    if norm(perp_vector) == 0, new_vertex = []; return; end
    unit_perp_vector = perp_vector / norm(perp_vector);
    v1 = midpoint + height * unit_perp_vector;
    v2 = midpoint - height * unit_perp_vector;
    hole_center_approx = mean([p1;p2],1);
    if norm(v1 - hole_center_approx) > norm(v2 - hole_center_approx), new_vertex = v1; else, new_vertex = v2; end
end

function cluster_indices = cluster_nodes_by_distance(coords, threshold)
    num_points = size(coords, 1);
    adj_matrix = false(num_points, num_points);
    for i = 1:num_points
        for j = i + 1:num_points
            dist = sqrt(sum((coords(i,:) - coords(j,:)).^2));
            if dist < threshold, adj_matrix(i, j) = true; adj_matrix(j, i) = true; end
        end
    end
    G = graph(adj_matrix, 'omitselfloops');
    cluster_indices = conncomp(G)';
end

function [x, y] = generate_network(area, num_nodes)
    x = area(1) * rand(1, num_nodes);
    y = area(2) * rand(1, num_nodes);
end

function plot_network_state(x, y, R, title_str, area)
    hold off;
    theta = linspace(0, 2*pi, 100);
    for i = 1:length(x)
        fill(x(i) + R*cos(theta), y(i) + R*sin(theta), [0, 0.8, 0.2], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        hold on;
    end
    plot(x, y, 'k.', 'MarkerSize', 12);
    axis equal;
    grid on;
    box on;
    title(title_str, 'FontSize', 14);
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    xlim([0, area(1)]); ylim([0, area(2)]);
end

function [Pi, CIPL, CIPU] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width)
  Pi = []; CIPL = []; CIPU = [];
  xi = x_loc(i); yi = y_loc(i);
  if isBorderSensor(xi, yi, net_length, net_width, Rs)
      border_intersections = borderIntersections(xi, yi, Rs, net_length, net_width);
      for p_idx = 1:size(border_intersections, 1)
          p = border_intersections(p_idx, :);
          if ~isCoveredByAnyNeighbor(p, i, Ni{i}, x_loc, y_loc, Rs), Pi = [Pi; p]; end
      end
  end
  for j = Ni{i}
      inter_pts = circleIntersections(xi, yi, x_loc(j), y_loc(j), Rs);
      for p_idx = 1:size(inter_pts, 1)
          p = inter_pts(p_idx, :);
          if ~isCoveredByAnotherSpecificNeighbor(p, i, j, Ni{i}, x_loc, y_loc, Rs), Pi = [Pi; p]; end
      end
  end
  if ~isempty(Pi), Pi = uniquetol(Pi, 1e-9, 'ByRows', true); end
end

function isB = isBorderSensor(x, y, net_length, net_width, Rs)
  isB = (x - Rs <= 0) || (x + Rs >= net_length) || (y - Rs <= 0) || (y + Rs >= net_width);
end

function points = circleIntersections(x1, y1, x2, y2, r)
  d = norm([x1-x2,y1-y2]);
  if d >= 2*r || d == 0, points = []; return; end
  a = d/2; h = sqrt(max(0, r^2 - a^2));
  xm = x1 + a * (x2-x1)/d; ym = y1 + a * (y2-y1)/d;
  rx = -(y2-y1)*(h/d); ry = (x2-x1)*(h/d);
  points = [xm + rx, ym + ry; xm - rx, ym - ry];
end

function pts = borderIntersections(x, y, r, L, W)
  pts = [];
  if (x - r <= 0), dy=sqrt(max(0,r^2-x^2)); pts=[pts;[0,y-dy];[0,y+dy]]; end
  if (x + r >= L), dy=sqrt(max(0,r^2-(L-x)^2)); pts=[pts;[L,y-dy];[L,y+dy]]; end
  if (y - r <= 0), dx=sqrt(max(0,r^2-y^2)); pts=[pts;[x-dx,0];[x+dx,0]]; end
  if (y + r >= W), dx=sqrt(max(0,r^2-(W-y)^2)); pts=[pts;[x-dx,W];[x+dx,W]]; end
  if ~isempty(pts), pts = pts(all(pts >= 0 & pts <= [L, W], 2),:); pts = uniquetol(pts, 1e-9, 'ByRows', true); end
end

function covered = isCoveredByPoint(point, sx, sy, sr)
  d = sqrt((point(1)-sx)^2 + (point(2)-sy)^2);
  covered = (d <= sr + 1e-9);
end

function covered = isCoveredByAnyNeighbor(p, i, N_i, x, y, Rs)
   covered = false;
   for k = N_i
      if (k ~= i) && isCoveredByPoint(p, x(k), y(k), Rs), covered = true; return; end
  end
end

function covered = isCoveredByAnotherSpecificNeighbor(p, i, j, N_i, x, y, Rs)
  covered = false;
  for k = N_i
      if (k ~= i) && (k ~= j) && isCoveredByPoint(p, x(k), y(k), Rs), covered = true; return; end
  end
end