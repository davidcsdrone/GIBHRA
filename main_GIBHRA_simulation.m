clc; clear; close all;
no_nodes = input('Enter the number of nodes: ');
net_length = input('Enter the Length of the network: ');
net_width = input('Enter the width of the network: ');
% --- Iteration Control Parameters ---
TARGET_COVERAGE_PERCENT = 98.0;
MAX_ITERATIONS = 25; % A safe limit to prevent infinite loops
iteration_count = 0;
nodes_added_this_run = 0; % To track nodes added in the current run
Rs = 30; % Sensing range radius
Rc = 2 * Rs; % communication range (currently unused for plotting, but good to define)
Ni = cell(no_nodes, 1); % Each cell holds sensing neighbors
x_loc = zeros(1, no_nodes);
y_loc = zeros(1, no_nodes);  
node_id = zeros(1, no_nodes);
figure;
hold on;
axis equal; 
xlim([0 net_length]); % Set x-axis limits
ylim([0 net_width]);  % Set y-axis limits
xlabel('Network Length');
ylabel('Network Width');
% Explicitly draw the network boundary
plot([0 net_length net_length 0 0], [0 0 net_width net_width 0], 'k-', 'LineWidth', 1.5); % Black border
for i = 1:no_nodes
  x_loc(i) = net_length * rand;
  y_loc(i) = net_width * rand;
  node_id(i) = i;
 
  theta = linspace(0, 2*pi, 100); % For a smooth circle
  circle_x = x_loc(i) + Rs * cos(theta);
  circle_y = y_loc(i) + Rs * sin(theta);
% A light grey fill with a solid, darker grey outline
    fill(circle_x, circle_y, [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.5); 
  % text(x_loc(i) + 0.02 * net_length, y_loc(i) + 0.02 * net_width, ...
  %     ['S_{', num2str(i), '}'], 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 8);
  %
end
disp('-----Calculating k-coverage before our hole restoration approach-----');
stats_before = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Initial Hole Area: %.2fm^2(%.2f%%\n', stats_before.hole_area, stats_before.hole_percentage);
fprintf('Initial Covered Area (k>=1): %.2f m^2 (%.2f%%)\n', stats_before.k_covered_area, stats_before.k_covered_percentage);
disp('--- Detailed Coverage Distribution (Before) ---');
for k = 1:length(stats_before.specific_k_areas)
   fprintf('Area with exactly k=%d coverage: %.2f m^2 (%.2f%%)\n', ...
           k, stats_before.specific_k_areas(k), stats_before.specific_k_percentages(k));
end
source = randi(no_nodes);
destination = randi(no_nodes);
for i = 1:no_nodes
  for j = 1:no_nodes
      if i ~= j
          d = sqrt((x_loc(i) - x_loc(j))^2 + (y_loc(i) - y_loc(j))^2);
          if d < 2 * Rs % This defines sensing neighbors (used for Ni)
              Ni{i} = [Ni{i}, j];
          end
      end
  end
end
AllCIPs = cell(no_nodes, 1); % To store CIPs for each node
AllCIPLs = cell(no_nodes, 1); % To store CIPLs for each node
AllCIPUs = cell(no_nodes, 1); % To store CIPUs for each node
nodeHoleInfo = cell(no_nodes, 1);
allBoundedHolesToPlot = {}; % To collect bounded hole polygons for plotting later
% Compute Critical Intersection Points (CIPs) and Covered Points (CPs) and Plot them
for i = 1:no_nodes
  % Get the classified CIPs and CPs for the current node
  [current_Pi, current_CIPL, current_CIPU] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width);
   % Store for plotting and table (combining all CIP types into AllCIPs for general red plotting)
  AllCIPs{i} = uniquetol([current_Pi; current_CIPL; current_CIPU], 1e-9, 'ByRows', true);
  AllCIPLs{i} = current_CIPL; % Keep original CIPL/CIPU for hole detection
  AllCIPUs{i} = current_CIPU; % Keep original CIPL/CIPU for hole detection
  % Plot CIPs (red squares as in research paper)
  if ~isempty(AllCIPs{i})
scatter(AllCIPs{i}(:,1), AllCIPs{i}(:,2), 30, 'r', '.', 'filled');   end
   % this is hole detection for Node i
  hole_type_CIPL = Hole_Type(AllCIPLs{i}, x_loc(i), y_loc(i));
  hole_type_CIPU = Hole_Type(AllCIPUs{i}, x_loc(i), y_loc(i));
   nodeHoleInfo{i}.CIPL_type = hole_type_CIPL;
  nodeHoleInfo{i}.CIPU_type = hole_type_CIPU;

  if ~strcmp(hole_type_CIPL, 'None')
      fprintf('Node S%d CIPL set detected: %s hole\n', i, hole_type_CIPL);
      if strcmp(hole_type_CIPL, 'Bounded')
          sorted_CIPL_points = sortPointsClockwise(AllCIPLs{i}, x_loc(i), y_loc(i));
          allBoundedHolesToPlot{end+1} = [sorted_CIPL_points; sorted_CIPL_points(1,:)]; % Close the loop for plotting
      end
  end
  if ~strcmp(hole_type_CIPU, 'None')
      fprintf('Node S%d CIPU set detected: %s hole\n', i, hole_type_CIPU);
      if strcmp(hole_type_CIPU, 'Bounded')
          sorted_CIPU_points = sortPointsClockwise(AllCIPUs{i}, x_loc(i), y_loc(i));
          allBoundedHolesToPlot{end+1} = [sorted_CIPU_points; sorted_CIPU_points(1,:)]; % Close the loop for plotting
      end
  end
end % End of the main plotting loop

while true % Loop will be broken from the inside
   iteration_count = iteration_count + 1;
   fprintf('\n\n====================================================\n');
   fprintf('--- Starting Restoration Iteration #%d ---\n', iteration_count);
  
   no_nodes_current = length(x_loc); % Use the CURRENT number of nodes
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
  
   disp(' ');
   disp('--- Identifying Final Boundary Node List ---');
  

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
      % We assume computeCP_and_CIP exists and works correctly
      [current_Pi, ~, ~] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width);
      AllCIPs{i} = current_Pi;
   end
   boundary_nodes_list = [];
   for i = 1:no_nodes_current
       if isempty(Ni{i}) || ~isempty(AllCIPs{i})
           boundary_nodes_list(end+1) = i;
       end
   end
   fprintf('Found %d boundary nodes.\n', length(boundary_nodes_list));
  
   % Triplets code
    T = selectTripletsDelaunay(boundary_nodes_list, x_loc, y_loc, Rs, current_coverage_percent);
  
   if isempty(T)
       fprintf('No more valid triplets found to restore. Stopping iterations.\n');
       break;
   end
  
   % deploying new nodes for this iteration
   x_loc_after = x_loc;
   y_loc_after = y_loc;
   nodes_deployed_this_iteration = 0;
  
   for t_idx = 1:length(T)
       current_triplet = T{t_idx};
      
       % note to self: we pass x_loc, not x_loc_after, to the solver to ensure that
       % all triplet decisions in this pass are based on the same initial state.
       [optimal_center, ~] = solveApolloniusForTriplet(current_triplet, x_loc, y_loc, Rs);
      
       if ~isempty(optimal_center)
           coverage_in_wsn = estimateCoverageInWSN(optimal_center, Rs, net_length, net_width);
           if coverage_in_wsn >= 0.50
               x_loc_after(end+1) = optimal_center(1);
               y_loc_after(end+1) = optimal_center(2);
               nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
              
               % this will plot the new node immediately
               new_node_id = length(x_loc_after);
% A light blue fill with a solid, darker blue outline
fill(optimal_center(1) + Rs*cos(theta), optimal_center(2) + Rs*sin(theta), [0.6 0.7 0.9], 'FaceAlpha', 0.6, 'EdgeColor', [0.2 0.4 0.6], 'LineWidth', 1.5);
% we also plot a center dot to show the new node's location
           end
       end
   end
  
   fprintf('Deployed %d new nodes in this iteration.\n', nodes_deployed_this_iteration);
   nodes_added_this_run = nodes_added_this_run + nodes_deployed_this_iteration;
  
   % Update the main node list for the next pass
   x_loc = x_loc_after;
   y_loc = y_loc_after;
  
   drawnow; % Update the plot to show the new nodes
   if nodes_deployed_this_iteration == 0
       fprintf('Could not deploy any new valid nodes. Stopping iterations.\n');
       break;
   end
end 
% --- our final summary report! ---
% =========================================================================
fprintf('\n\n====================================================\n');
fprintf('--- FINAL RESTORATION SUMMARY ---\n');
% We subtract 1 from iteration_count because the last pass only checks coverage
% and doesn't deploy new nodes.
fprintf('Algorithm stopped after %d restoration iteration(s).\n', iteration_count - 1);
% calculate how many nodes were initial vs. added
no_nodes_initial = no_nodes; % Assuming 'no_nodes' is your initial input
total_nodes_added = length(x_loc) - no_nodes_initial;
fprintf('Total initial nodes: %d\n', no_nodes_initial);
fprintf('Total new nodes added: %d\n', total_nodes_added);
fprintf('Final total nodes: %d\n', length(x_loc));
%final coverage calculation using the final state of x_loc and y_loc
final_stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Final Coverage (k>=1): %.2f%%\n', final_stats.k_covered_percentage);
fprintf('Final Hole Area (k=0): %.2f%%\n', 100 - final_stats.k_covered_percentage);
fprintf('====================================================\n');
fprintf('\n\n--- COVERAGE COMPARISON ---\n');
fprintf('                           BEFORE       AFTER\n');
fprintf('--------------------------------------------------\n');
fprintf('Total Nodes:               %-12d %d\n', no_nodes_initial, length(x_loc));
fprintf('Coverage (k>=1):          %-12.2f%% %.2f%%\n', stats_before.k_covered_percentage, final_stats.k_covered_percentage);
fprintf('Hole Area (k=0):           %-12.2f%% %.2f%%\n', stats_before.hole_percentage, final_stats.hole_percentage);
fprintf('====================================================\n');
% You can place any other final analysis or plotting code after this.
% For example, you might want to remove the CIP markers for a cleaner final plot.
% findobj(gca, 'Type', 'scatter', 'Marker', 's').delete(); % Example to remove red squares
% =========================================================================
% --- END OF THE ITERATIVE RESTORATION LOOP ---
% =========================================================================
% --- FINAL STEP: Recalculate K-Coverage AFTER Restoration ---
% disp(' ');
% disp('-----Calculating k-coverage AFTER our hole restoration approach-----');
% % fprintf('Final Covered Area (k>=1): %.2f m^2 (%.2f%%)\n', stats_after.k_covered_area, stats_after.k_covered_percentage);
% disp('--- Detailed Coverage Distribution (After) ---');
% % for k = 1:length(stats_after.specific_k_areas)
%    fprintf('Area with exactly k=%d coverage: %.2f m^2 (%.2f%%)\n', ...
%            k, stats_after.specific_k_areas(k), stats_after.specific_k_percentages(k));
% end
% Compare before and after
% area_restored = stats_before.hole_area - stats_after.hole_area;
% fprintf('\n--- SUMMARY ---\n');
% fprintf('Total area restored: %.2f m^2\n', area_restored);
% fprintf('Coverage hole reduced from %.2f%% to %.2f%%.\n', stats_before.hole_percentage, stats_after.hole_percentage);
% --- Plot communication links (dashed lines) ---
for i = 1:no_nodes
  xi = x_loc(i);
  yi = y_loc(i);
  for j = Ni{i} % Iterate through neighbors of node i
      xj = x_loc(j);
      yj = y_loc(j);
    
      % Only draw each link once (e.g., draw i-j but not j-i)
      if i < j
          % plot([xi, xj], [yi, yj], 'k--', 'LineWidth', 0.5); % Black dashed line
      end
  end
end
% --- START OF CODE FOR LABELLING ALL INTERSECTION POINTS ---
% 1. Collect all CIPs and CPs into a single list
all_intersection_points = [];
for i = 1:no_nodes
  % Add CIPs from current node
  if ~isempty(AllCIPs{i})
      all_intersection_points = [all_intersection_points; AllCIPs{i}];
  end
end % Correct end of this loop
% 2. Remove duplicate points
if ~isempty(all_intersection_points)
  all_unique_points = uniquetol(all_intersection_points, 1e-9, 'ByRows', true); % Using uniquetol for floating point comparisons
  % 3. Iterate through unique points and add sequential labels
  point_counter = 0;
  for p_idx = 1:size(all_unique_points, 1)
      point_counter = point_counter + 1;
      px = all_unique_points(p_idx, 1);
      py = all_unique_points(p_idx, 2);
      % Place the label slightly offset from the point for better visibility
      % text(px + 0.01 * net_length, py + 0.01 * net_width, ...
      %      ['P_{', num2str(point_counter), '}'], ... % Label format P_n using TeX interpreter
      %      'Interpreter', 'tex', ...
      %      'Color', 'k', ... % Black color for the label
      %      'FontSize', 8);
  end
end
% Initialize cell arrays to store data for the table
node_names_for_table = cell(no_nodes, 1);
neighbors_for_table = cell(no_nodes, 1);
cip_points_for_table = cell(no_nodes, 1);
disp(' '); % Add a blank line for readability
disp('--- Table of Nodes, Neighbors, and Critical Intersection Points (CIPs) ---');
disp(' ');
for i = 1:no_nodes
  % 1. Node Name (e.g., 'S1', 'S2')
  node_names_for_table{i} = ['S', num2str(i)];
  % 2. Neighbors (e.g., 'S2, S10, S11')
  if isempty(Ni{i})
      neighbors_for_table{i} = 'None (Isolated)';
  else
      neighbor_str = '';
      for k = 1:length(Ni{i})
          neighbor_str = [neighbor_str, 'S', num2str(Ni{i}(k))];
          if k < length(Ni{i})
              neighbor_str = [neighbor_str, ', '];
          end
      end
      neighbors_for_table{i} = neighbor_str;
  end
  % 3. Corresponding CIPs (Formatted as a string, e.g., '[x1 y1; x2 y2; ...]')
  if isempty(AllCIPs{i}) % <--- CHANGED: Use AllCIPs
      cip_points_for_table{i} = 'None';
  else
      % Format the CIP matrix into a readable string
      cip_str = '[';
      for row = 1:size(AllCIPs{i}, 1) % <--- CHANGED: Use AllCIPs
          cip_str = [cip_str, sprintf('%.2f %.2f', AllCIPs{i}(row, 1), AllCIPs{i}(row, 2))]; % <--- CHANGED: Use AllCIPs
          if row < size(AllCIPs{i}, 1) % <--- CHANGED: Use AllCIPs
              cip_str = [cip_str, '; ']; % Add semicolon to separate rows in matrix representation
          end
      end
      cip_str = [cip_str, ']'];
      cip_points_for_table{i} = cip_str;
  end
 end
% Create the MATLAB table
T = table(node_names_for_table, neighbors_for_table, cip_points_for_table, ...
        'VariableNames', {'Node', 'Neighbors', 'CorrespondingCIP'});
% Add this right before `hold off;` at the end of the main script, after the CIP labeling.
% New cell array to store uncovered arcs for plotting
allUncoveredArcsToPlot = cell(no_nodes, 1);
% % --- START OF ARC CONSTRUCTION AND CLASSIFICATION ---
% disp(' ');
% disp('--- Constructing Covered and Uncovered Arcs ---');
% for i = 1:no_nodes
%    S1_x = x_loc(i);
%    S1_y = y_loc(i);
%    current_node_neighbors = Ni{i}; % Neighbors of S1
%
%    current_node_cip_sets = {AllCIPLs{i}, AllCIPUs{i}};
%
%    uncovered_arcs_for_S1 = {}; % To store [start_x, start_y, end_x, end_y] for S1
%    for set_idx = 1:length(current_node_cip_sets)
%        current_cip_set = current_node_cip_sets{set_idx};
%        if size(current_cip_set, 1) < 2
%            continue; % Need at least two points to form an arc
%        end
%
%         sorted_cip_set = sortPointsClockwise(current_cip_set, S1_x, S1_y);
%
%         for p_k = 1:size(sorted_cip_set, 1)
%            p_start = sorted_cip_set(p_k, :);
%
%             p_end = sorted_cip_set(mod(p_k, size(sorted_cip_set, 1)) + 1, :);
%
%             neighbor_indices = Ni{i};
%             neighbor_coords = []; % Initialize as empty for nodes with no neighbors
%             if ~isempty(neighbor_indices)
%      % Create an Mx2 matrix of neighbor coordinates
%              neighbor_coords = [x_loc(neighbor_indices)', y_loc(neighbor_indices)'];
%             end
%
% % Call the new, corrected function with the correct arguments
%
% is_arc_covered_flag = isArcCovered(x_loc(i), y_loc(i), p_start, p_end, Rs, Ni{i}, x_loc, y_loc, i);            if ~is_arc_covered_flag
%                 uncovered_arcs_for_S1{end+1} = struct('S1_idx', i, ...
%                                                    'p_start', p_start, ...
%                                                    'p_end', p_end);
%                fprintf('Node S%d: Uncovered arc detected between (%.2f,%.2f) and (%.2f,%.2f)\n', ...
%                        i, p_start(1), p_start(2), p_end(1), p_end(2));
%             end
%        end
%    end
%    allUncoveredArcsToPlot{i} = uncovered_arcs_for_S1;
% end
% % --- END OF ARC CONSTRUCTION AND CLASSIFICATION ---
% --- Plotting Uncovered Arcs ---
hold on; % Ensure we plot on the existing figure
for i = 1:no_nodes
  uncovered_arcs_S1 = allUncoveredArcsToPlot{i};
  S1_x = x_loc(i);
  S1_y = y_loc(i);
  for k = 1:length(uncovered_arcs_S1)
      arc_info = uncovered_arcs_S1{k};
      p_start = arc_info.p_start;
      p_end = arc_info.p_end;
      % Get angles of start and end points relative to S1
      angle_start = atan2(p_start(2) - S1_y, p_start(1) - S1_x);
      angle_end = atan2(p_end(2) - S1_y, p_end(1) - S1_x);
      % Adjust angles to sweep clockwise for plotting (consistent with sortPointsClockwise)
      if angle_end > angle_start
          angle_end = angle_end - 2*pi;
      end
    
      % Generate points for the arc
      num_arc_points = 50; % For a smooth arc
      arc_theta = linspace(angle_start, angle_end, num_arc_points);
      arc_x = S1_x + Rs * cos(arc_theta);
      arc_y = S1_y + Rs * sin(arc_theta);
      plot(arc_x, arc_y, 'r-', 'LineWidth', 2); % Plot uncovered arc in red
  end
end
% --- Plotting Bounded Holes (filled green) ---
hold on;
if ~isempty(allBoundedHolesToPlot)
  disp(' ');
  disp('--- Plotting Bounded Coverage Holes (Green) ---');
  hole_count = 0;
  for k = 1:length(allBoundedHolesToPlot)
      hole_polygon = allBoundedHolesToPlot{k};
      if ~isempty(hole_polygon)
          hole_count = hole_count + 1;
          fill(hole_polygon(:,1), hole_polygon(:,2), 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Green fill with some transparency
        
          % Optionally, label the holes (similar to H1, H2 in the image)
          % Calculate centroid for label placement
          if size(hole_polygon, 1) > 2
              centroid_x = mean(hole_polygon(1:end-1, 1)); % Exclude last point (duplicate of first)
              centroid_y = mean(hole_polygon(1:end-1, 2));
              text(centroid_x, centroid_y, ['H_{', num2str(hole_count), '}'], ...
                   'Interpreter', 'tex', 'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold', ...
                   'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
          end
      end
  end
  if hole_count == 0
      disp('No bounded holes detected by the current algorithm.');
  end
else
  disp('No bounded holes detected by the current algorithm.');
end
hold on;
h = zeros(2, 1);
h(1) = fill(NaN, NaN, [0.8 0.8 0.8], 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1.5);
h(2) = fill(NaN, NaN, [0.7 0.8 0.9], 'EdgeColor', [0.2 0.4 0.6], 'LineWidth', 1.5);
legend(h, {'Initial Sensor', 'Restorative Sensor'}, 'Location', 'northwest', 'FontSize', 12, 'Box', 'on');

hold off;
hold off;
function sorted_points = sortPointsClockwise(points, center_x, center_y)
  % Sorts a set of 2D points (Nx2 matrix) in clockwise order around a given center.
  % If center is not provided, it uses the centroid of the points.
  if nargin < 2
      % Calculate centroid if no center is provided
      center_x = mean(points(:,1));
      center_y = mean(points(:,2));
  end
  % Calculate angles relative to the center
  angles = atan2(points(:,2) - center_y, points(:,1) - center_x);
   angles_positive = mod(angles, 2*pi); % Convert to [0, 2*pi]
  [~, sort_idx] = sort(angles_positive, 'descend'); % Sort descending for clockwise
  sorted_points = points(sort_idx, :);
end
function is_closed = isLoopClosed(points_sequence, tolerance)
   if isempty(points_sequence) || size(points_sequence, 1) < 2
      is_closed = false;
      return;
  end
   if nargin < 2
      tolerance = 1e-9;
  end
  first_point = points_sequence(1, :);
  last_point = points_sequence(end, :);
  distance = sqrt((first_point(1) - last_point(1))^2 + (first_point(2) - last_point(2))^2);
  is_closed = distance < tolerance;
end
function hole_type = Hole_Type(points_set, node_x, node_y)
  % Implements Algorithm 4, lines 31-35.
  % points_set: CIPL or CIPU set for a specific node (Nx2 matrix)
  % node_x, node_y: Coordinates of the current sensor node (Si)
  hole_type = 'None'; % Default: No hole detected for this set
  if isempty(points_set)
      return;
  end
  % Arrange points in clockwise order around the current node as a reference
  % The paper implies ordering based on the "line connecting to the location of any two neighboring sensors"
  % A simple approach for ordering points that define a 'hole' is around their centroid.
  % However, the paper's description suggests it's from the perspective of Si.
  % Let's use the node's location as a reference point for sorting.
  sorted_points = sortPointsClockwise(points_set, node_x, node_y);
  % Check if the initial and terminal points are the same (forms a loop)
  if isLoopClosed(sorted_points)
      hole_type = 'Bounded';
  else
      hole_type = 'Unbounded';
  end
end
% --- Functions below ---
function isBorder = isBorderSensor(x, y, net_length, net_width, Rs)
  isBorder = (x - Rs <= 0) || (x + Rs >= net_length) || (y - Rs <= 0) || (y + Rs >= net_width);
end
function points = circleIntersections(x1, y1, x2, y2, r)
  d = sqrt((x2 - x1)^2 + (y2 - y1)^2);
  if d >= 2*r || d == 0 || d <= abs(r - r)
      points = [];
      return;
  end
  a = (r^2 - r^2 + d^2) / (2*d);
  h = sqrt(r^2 - a^2);
  xm = x1 + a * (x2 - x1) / d;
  ym = y1 + a * (y2 - y1) / d;
  rx = -(y2 - y1) * (h / d);
  ry = -(x2 - x1) * (h / d);
  points = [
      xm + rx, ym - ry;
      xm - rx, ym + ry
  ];
end
function [Pi, CIPL, CIPU] = computeCP_and_CIP(i, x_loc, y_loc, Ni, Rs, net_length, net_width)
  Pi = [];   % Initialize Pi (Critical Intersection Points) for node i
  CIPL = []; % Initialize CIPL (CIPs Below Line) for node i
  CIPU = []; % Initialize CIPU (CIPs Upper Line) for node i
  xi = x_loc(i);
  yi = y_loc(i);
  % --- Step 1-7: Check for Border Sensor Intersection Points ---
  if isBorderSensor(xi, yi, net_length, net_width, Rs)
      border_intersections = borderIntersections(xi, yi, Rs, net_length, net_width);
      for p_idx = 1:size(border_intersections, 1)
          p = border_intersections(p_idx, :);
        
          if ~isCoveredByAnyNeighbor(p, i, Ni{i}, x_loc, y_loc, Rs)
              Pi = [Pi; p]; % p is a CIP
              % A border intersection point is considered a CIP (Pi) if it's
              % not covered by any *other* neighbor.
              % Algorithm 3 (Localization) does NOT apply to border CIPs,
              % only to CIPs formed by intersections of two sensor circles.
          end
      end
  end
  % --- Step 8-20: Scan all members of Ni for intersection points with node i ---
  for j = Ni{i} % Iterate through each neighbor 'j' of node 'i'
      xj = x_loc(j);
      yj = y_loc(j);
    
      % Find intersection points of node i's sensing circle with node j's sensing circle
      inter_pts = circleIntersections(xi, yi, xj, yj, Rs);
    
      for p_idx = 1:size(inter_pts, 1)
          p = inter_pts(p_idx, :);
        
          % Check if the intersection point 'p' is covered by any other neighbor 'Sk' (k != i, k != j)
          if ~isCoveredByAnotherSpecificNeighbor(p, i, j, Ni{i}, x_loc, y_loc, Rs)
              Pi = [Pi; p]; % p is a Critical Intersection Point (CIP)
            
              % --- Algorithm 3: Localization of CIP ---
              % Check the location of this CIP with respect to the line Lij (connecting Si and Sj)
              status = isAboveLine(p(1), p(2), xi, yi, xj, yj); % Pass point coordinates and line endpoints
            
              if strcmp(status, 'above')
                  CIPU = [CIPU; p]; % Assign CIP to CIPU (Upper)
              elseif strcmp(status, 'below')
                  CIPL = [CIPL; p]; % Assign CIP to CIPL (Lower)
            
              end
              % --- End Algorithm 3 Localization ---
          end
      end
  end
  % Post-processing: Remove any duplicate CIPs, CPs, CIPLs, and CIPUs.
  % Use 'uniquetol' for robust comparison of floating-point coordinates.
  if ~isempty(Pi)
      Pi = uniquetol(Pi, 1e-9, 'ByRows', true);
  end
 
  if ~isempty(CIPL)
      CIPL = uniquetol(CIPL, 1e-9, 'ByRows', true);
  end
  if ~isempty(CIPU)
      CIPU = uniquetol(CIPU, 1e-9, 'ByRows', true);
  end
end
function pts = borderIntersections(x, y, r, net_length, net_width)
  pts = [];
   if (x - r <= 0) && (x + r >= 0) % Check if circle potentially crosses or touches x=0
      delta_x = abs(x - 0);
      if delta_x <= r % If it's not just tangent at (0,y) or outside
          dy = sqrt(r^2 - delta_x^2);
          p1 = [0, y - dy];
          p2 = [0, y + dy];
          % Add points if they are within the y-bounds of the network
          if p1(2) >= 0 && p1(2) <= net_width
              pts = [pts; p1];
          end
          if dy ~= 0 && p2(2) >= 0 && p2(2) <= net_width
              pts = [pts; p2];
          end
      end
  end
  % Intersection with Right Boundary (x_line = net_length)
  if (x + r >= net_length) && (x - r <= net_length) % Check if circle potentially crosses or touches x=net_length
      delta_x = abs(x - net_length);
      if delta_x <= r
          dy = sqrt(r^2 - delta_x^2);
          p1 = [net_length, y - dy];
          p2 = [net_length, y + dy];
          if p1(2) >= 0 && p1(2) <= net_width
              pts = [pts; p1];
          end
          if dy ~= 0 && p2(2) >= 0 && p2(2) <= net_width
              pts = [pts; p2];
          end
      end
  end
   % Intersection with Bottom Boundary (y_line = 0)
  if (y - r <= 0) && (y + r >= 0) % Check if circle potentially crosses or touches y=0
      delta_y = abs(y - 0);
      if delta_y <= r
          dx = sqrt(r^2 - delta_y^2);
          p1 = [x - dx, 0];
          p2 = [x + dx, 0];
          if p1(1) >= 0 && p1(1) <= net_length
              pts = [pts; p1];
          end
          if dx ~= 0 && p2(1) >= 0 && p2(1) <= net_length
              pts = [pts; p2];
          end
      end
  end
  % Intersection with Top Boundary (y_line = net_width)
  if (y + r >= net_width) && (y - r <= net_width) % Check if circle potentially crosses or touches y=net_width
      delta_y = abs(y - net_width);
      if delta_y <= r
          dx = sqrt(r^2 - delta_y^2);
          p1 = [x - dx, net_width];
          p2 = [x + dx, net_width];
          if p1(1) >= 0 && p1(1) <= net_length
              pts = [pts; p1];
          end
          if dx ~= 0 && p2(1) >= 0 && p2(1) <= net_length
              pts = [pts; p2];
          end
      end
  end
  % Remove duplicate points and points slightly outside the boundary due to floating point
  if ~isempty(pts)
      pts = uniquetol(pts, 1e-9, 'ByRows', true); % Use uniquetol for robustness
      % Further filter to ensure points are strictly within or on the boundary
      pts = pts(pts(:,1) >= 0 & pts(:,1) <= net_length & pts(:,2) >= 0 & pts(:,2) <= net_width, :);
  end
end
% NEW FUNCTION: isCoveredByPoint
function covered = isCoveredByPoint(point, sensor_x, sensor_y, sensor_r)
  % Checks if a given point is within the sensing range of a sensor.
  % point: [px, py]
  % sensor_x, sensor_y: coordinates of the covering sensor
  % sensor_r: sensing radius of the covering sensor
  d = sqrt((point(1) - sensor_x)^2 + (point(2) - sensor_y)^2);
  covered = (d <= sensor_r);
end
function covered = isCoveredByAnyNeighbor(point, current_node_idx, neighbors_list, x_loc, y_loc, Rs)
   covered = false;
  for k_val = neighbors_list % Iterate through actual neighbor IDs
   
      if (k_val ~= current_node_idx) % Ensure we're not checking the current node itself (which is Si)
           if isCoveredByPoint(point, x_loc(k_val), y_loc(k_val), Rs)
              covered = true;
              return;
          end
      end
  end
end
% NEW FUNCTION: isCoveredByAnotherSpecificNeighbor (for circle-circle points)
function covered = isCoveredByAnotherSpecificNeighbor(point, Si_idx, Sj_idx, neighbors_list, x_loc, y_loc, Rs)
  % Checks if 'point' (an intersection of Si and Sj) is covered by a neighbor Sk.
  % Sk must be a member of neighbors_list, but not Si_idx or Sj_idx.
  covered = false;
  for k_val = neighbors_list % Iterate through actual neighbor IDs
      % k must be distinct from both the current node (Si) and the neighbor (Sj) that formed the intersection.
      if (k_val ~= Si_idx) && (k_val ~= Sj_idx)
          if isCoveredByPoint(point, x_loc(k_val), y_loc(k_val), Rs)
              covered = true;
              return;
          end
      end
  end
end
function status = isAboveLine(point_x, point_y, x1, y1, x2, y2)
  tolerance = 1e-9;
  % Case 1: Vertical Line (x1 == x2)
  if abs(x2 - x1) < tolerance
      if abs(point_x - x1) < tolerance % Point is on the vertical line
          status = 'on';
      elseif point_x > x1 % Point is to the right of the vertical line
          status = 'above'; % Convention for vertical lines in this context
      else % Point is to the left of the vertical line
          status = 'below'; % Convention for vertical lines in this context
      end
  % Case 2: Horizontal Line (y1 == y2)
  elseif abs(y2 - y1) < tolerance
      if abs(point_y - y1) < tolerance % Point is on the horizontal line
          status = 'on';
      elseif point_y > y1
          status = 'above';
      else
          status = 'below';
      end
  % Case 3: Sloped Line
  else
      % Calculate the y-coordinate on the line at 'point_x'
      % y = m*x + b, where m = (y2-y1)/(x2-x1) and b = y1 - m*x1
      y_on_line = y1 + (point_x - x1) * (y2 - y1) / (x2 - x1);
      if point_y > y_on_line + tolerance
          status = 'above';
      elseif point_y < y_on_line - tolerance
          status = 'below';
      else
          status = 'on';
      end
  end
end
% --- NEW FUNCTION: calculateAngle ---
function angle_rad = calculateAngle(S1_x, S1_y, p1_x, p1_y, p2_x, p2_y)
  % Calculates the angle between two vectors S1->p1 and S1->p2.
  % S1: center point (sensor node)
  % p1, p2: end points of the vectors (CIPs)
   vec_S1p1 = [p1_x - S1_x, p1_y - S1_y];
  vec_S1p2 = [p2_x - S1_x, p2_y - S1_y];
   dot_product = dot(vec_S1p1, vec_S1p2);
  magnitude_S1p1 = norm(vec_S1p1);
  magnitude_S1p2 = norm(vec_S1p2);
   % Avoid division by zero or NaN if points are coincident
  if magnitude_S1p1 == 0 || magnitude_S1p2 == 0
      angle_rad = 0; % Or handle as an error, depending on desired behavior
      return;
  end
   cosine_theta = dot_product / (magnitude_S1p1 * magnitude_S1p2);
   % Ensure cosine_theta is within the valid range [-1, 1] to prevent issues with acos due to floating point inaccuracies
  cosine_theta = max(-1, min(1, cosine_theta));
   angle_rad = acos(cosine_theta);
end

function covered = isArcCovered(S1_center, p_start, p_end, Rs, neighbor_centers)
% Determines if an arc is "covered" by checking if its midpoint is covered
% by any neighboring node's sensing disc.
   covered = false; % Assume the arc is UNCOVERED by default.
   % If there are no neighbors, no arcs can be covered by them.
   if isempty(neighbor_centers)
       return;
   end
  
   % --- Step 1: Calculate the geometric midpoint of the arc ---
   % Create vectors from the center to the start and end points of the arc
   vec_start = p_start - S1_center;
   vec_end = p_end - S1_center;
  
   % The vector to the midpoint is the sum of these two vectors.
   vec_mid = vec_start + vec_end;
  
   % Handle the rare case where the arc is a perfect 180 degrees.
   if norm(vec_mid) < 1e-9
       vec_mid = [-vec_start(2), vec_start(1)]; % Use a perpendicular vector
   end
   midpoint_on_circle = S1_center + Rs * (vec_mid / norm(vec_mid));
  
   % --- Step 2: Check if this single midpoint is covered by any neighbor ---
   for k = 1:size(neighbor_centers, 1)
       neighbor_center = neighbor_centers(k, :);
      
       dist_to_neighbor = sqrt((midpoint_on_circle(1) - neighbor_center(1))^2 + (midpoint_on_circle(2) - neighbor_center(2))^2);
      
       if dist_to_neighbor <= Rs + 1e-9 % Add tolerance for floating point
           covered = true;
           return; % The arc is COVERED. Exit immediately.
       end
   end
end
function T = selectTriplets(B_nodes, x_loc, y_loc, Rs)
% Implements the O(m^2) greedy algorithm to select non-overlapping triplets.
% INPUTS:
%   B_nodes: A vector of node IDs that are boundary nodes.
%   x_loc:   The global array of all node x-coordinates.
%   y_loc:   The global array of all node y-coordinates.
%   Rs:      The uniform sensing radius.
%
% OUTPUT:
%   T:       A cell array, where each cell contains a 1x3 vector of node IDs
%            representing a selected triplet.
   fprintf('Starting triplet selection...\n');
  
   T = {}; % Initialize the list of final triplets
   used_nodes = false(1, max(B_nodes)); % A fast way to track used nodes
   % The pseudocode sorts by node_id. Since B_nodes is just a list of IDs,
   % sorting it numerically achieves this.
   B_sorted = sort(B_nodes);
   m = length(B_sorted);
   % --- Main Loop (FOR EACH node Si in B) ---
   for i = 1:m
       Si_id = B_sorted(i);
      
       % Skip if this node has already been used in another triplet
       if used_nodes(Si_id)
           continue;
       end
      
       % --- Build the Eligible_Partners list for Si ---
       eligible_partners = [];
       for p = 1:m
           Sp_id = B_sorted(p);
          
           % A partner must be a different, unused node
           if Sp_id ~= Si_id && ~used_nodes(Sp_id)
              
               % Check for non-overlapping condition
               dist = sqrt((x_loc(Si_id) - x_loc(Sp_id))^2 + (y_loc(Si_id) - y_loc(Sp_id))^2);
               if dist >= 2 * Rs
                   % Store the partner's ID and its distance to Si
                   eligible_partners(end+1, :) = [Sp_id, dist];
               end
           end
       end
      
       % --- Find the two closest partners ---
       if size(eligible_partners, 1) < 2
           continue; % Not enough partners to form a triplet, move to next Si
       end
      
       % Sort partners by distance (the second column) to find the closest
       sorted_partners = sortrows(eligible_partners, 2);
      
       Sj_id = sorted_partners(1, 1); % The closest partner
       Sk_id = sorted_partners(2, 1); % The second-closest partner
      
       % --- Final Validation Check ---
       % Check if the two closest partners overlap with each other
       dist_jk = sqrt((x_loc(Sj_id) - x_loc(Sk_id))^2 + (y_loc(Sj_id) - y_loc(Sk_id))^2);
      
       if dist_jk >= 2 * Rs
           % Success! We found a valid triplet.
           new_triplet = [Si_id, Sj_id, Sk_id];
           T{end+1} = new_triplet;
          
           % Mark these three nodes as used
           used_nodes(Si_id) = true;
           used_nodes(Sj_id) = true;
           used_nodes(Sk_id) = true;
          
           fprintf('Triplet found: S%d, S%d, S%d\n', Si_id, Sj_id, Sk_id);
       end
   end
  
   fprintf('Triplet selection complete. Found %d triplets.\n', length(T));
end

function coverage_ratio = estimateCoverageInWSN(center, radius, net_length, net_width)
% estimates the percentage of a circle's area that is inside the WSN
% boundaries using a Monte Carlo method.
   num_samples = 10000; % Increase for more accuracy, decrease for more speed.
  
   % define the bounding box for the circle
   x_min = center(1) - radius;
   x_max = center(1) + radius;
   y_min = center(2) - radius;
   y_max = center(2) + radius;
  
   % generate random points within the circle's bounding box
   rand_x = x_min + (x_max - x_min) * rand(num_samples, 1);
   rand_y = y_min + (y_max - y_min) * rand(num_samples, 1);
  
   % calculate the squared distance from the center to each random point
   dist_sq = (rand_x - center(1)).^2 + (rand_y - center(2)).^2;
  
   % find out which points are actually inside the circle
   is_in_circle_mask = dist_sq <= radius^2;
   points_in_circle = sum(is_in_circle_mask);
  
   if points_in_circle == 0
       coverage_ratio = 0;
       return;
   end
  
   % of the points inside the circle, find which are also inside the WSN
   in_circle_x = rand_x(is_in_circle_mask);
   in_circle_y = rand_y(is_in_circle_mask);
  
   is_in_wsn_mask = in_circle_x >= 0 & in_circle_x <= net_length & ...
                      in_circle_y >= 0 & in_circle_y <= net_width;
                     
   points_in_wsn = sum(is_in_wsn_mask);
  
   % calculate the ratio
   coverage_ratio = points_in_wsn / points_in_circle;
end

