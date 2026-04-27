clc; clear; close all;
rng('default');
no_nodes = input('Enter the number of nodes: ');
net_length = input('Enter the Length of the network: ');
net_width = input('Enter the width of the network: ');
Rs = 20; % Sensing range radius
Rc = 2 * Rs; % communication range (currently unused for plotting, but good to define)
Ni = cell(no_nodes, 1); % Each cell holds sensing neighbors
x_loc = zeros(1, no_nodes);
y_loc = zeros(1, no_nodes);
node_id = zeros(1, no_nodes);
figure;
hold on;
grid on;
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
 
  % Plot sensing range as a blue circle (replacing the rectangle)
  theta = linspace(0, 2*pi, 100); % For a smooth circle
  circle_x = x_loc(i) + Rs * cos(theta);
  circle_y = y_loc(i) + Rs * sin(theta);
  fill(circle_x, circle_y, 'g', 'FaceAlpha', 0.1, 'EdgeColor','b', 'LineWidth',1);
 
  plot(x_loc(i), y_loc(i), 'k.', 'MarkerSize', 15);
  text(x_loc(i) + 0.02 * net_length, y_loc(i) + 0.02 * net_width, ...
      ['S_{', num2str(i), '}'], 'Interpreter', 'tex', 'Color', 'k', 'FontSize', 8);
 
  pause(0.1);
end

% =========================================================================
% --- START OF THE ITERATIVE RESTORATION LOOP (CIRCUMCENTER METHOD) ---
% =========================================================================




% =========================================================================
% --- SETUP FOR ITERATIVE RESTORATION ---
% =========================================================================
% Define the control parameters for the iterative loop
TARGET_COVERAGE_PERCENT = 98.0;
MAX_ITERATIONS = 10; % A safe limit to prevent infinite loops
% --- Data Collection for Plotting ---
iteration_data = [];
coverage_per_iteration = [];
nodes_added_per_iteration = [];
iteration_count = 0;

% Store the initial number of nodes for the final report
no_nodes_initial = length(x_loc); 

% =========================================================================
% --- START OF THE ITERATIVE RESTORATION LOOP (CIRCUMCENTER METHOD) ---
% =========================================================================
while true % Loop will be broken from the inside
    iteration_count = iteration_count + 1;
    fprintf('\n\n====================================================\n');
    fprintf('--- Starting Restoration Iteration #%d ---\n', iteration_count);
    
    no_nodes_current = length(x_loc); % Use the CURRENT number of nodes
    fprintf('Current node count: %d\n', no_nodes_current);

    % --- Analyze Current Network State ---
    stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    current_coverage_percent = stats_current.k_covered_percentage;
    fprintf('Current Coverage (k>=1): %.2f%%\n', current_coverage_percent);

    % --- Check Stopping Conditions ---
    if current_coverage_percent >= TARGET_COVERAGE_PERCENT
        fprintf('\nTarget coverage of %.2f%% reached. Stopping iterations.\n', TARGET_COVERAGE_PERCENT);
        break;
    end
    
    if iteration_count > MAX_ITERATIONS
        fprintf('\nMaximum number of iterations (%d) reached. Stopping.\n', MAX_ITERATIONS);
        break;
    end
    
    % --- HOLE RESTORATION STAGE (CIRCUMCENTER METHOD) ---
    
    % Step 1: Find all uncovered holes using the paper's method.
    hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs);
    
    if isempty(hole_centers)
        fprintf('No more valid holes found to restore. Stopping iterations.\n');
        break;
    end
    
    % Step 2: Attempt to deploy a new node at each identified hole center.
    x_loc_after = x_loc;
    y_loc_after = y_loc;
    nodes_deployed_this_iteration = 0;
    
    for i = 1:size(hole_centers, 1)
        potential_center = hole_centers(i, :);
        
        % Use the same 50% boundary check for fair comparison
        coverage_in_wsn = estimateCoverageInWSN(potential_center, Rs, net_length, net_width);
        
        if coverage_in_wsn >= 0.5
            % Add the new node's location to a temporary list
            x_loc_after(end+1) = potential_center(1);
            y_loc_after(end+1) = potential_center(2);
            nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
            
            % Plot the new node immediately using a distinct professional style
            % For the baseline, we'll use magenta squares to distinguish them
            % from your Apollonius algorithm's blue triangles.
            plot(potential_center(1), potential_center(2), 'ms', 'MarkerSize', 6, 'MarkerFaceColor', 'm');
            theta = linspace(0, 2*pi, 100);
            fill(potential_center(1) + Rs*cos(theta), potential_center(2) + Rs*sin(theta), [1.0 0.8 0.9], 'FaceAlpha', 0.7, 'EdgeColor', 'm', 'LineWidth', 1.5);
        end
    end
    
    fprintf('Deployed %d new nodes in this iteration.\n', nodes_deployed_this_iteration);
    
    % Update the main node list for the next pass
    x_loc = x_loc_after;
    y_loc = y_loc_after;
    % --- Record data for this iteration's plot ---
    stats_after_pass = calculateKCoverage(x_loc_after, y_loc_after, Rs, net_length, net_width);
    coverage_after_pass = stats_after_pass.k_covered_percentage;
    
    coverage_gained_this_iteration = coverage_after_pass - current_coverage_percent;
    
    iteration_data(end+1, :) = [iteration_count, nodes_deployed_this_iteration, coverage_gained_this_iteration];
    drawnow; % Update the plot to show the new nodes

    if nodes_deployed_this_iteration == 0
        fprintf('Could not deploy any new valid nodes. Stopping iterations.\n');
        break;
    end
end % This is the end for the 'while true' loop

% =========================================================================
% --- FINAL SUMMARY REPORT ---
% =========================================================================
fprintf('\n\n====================================================\n');
fprintf('--- FINAL RESTORATION SUMMARY (Circumcenter Method) ---\n');
fprintf('Algorithm stopped after %d restoration iteration(s).\n', iteration_count - 1);
total_nodes_added = length(x_loc) - no_nodes_initial;
fprintf('Total initial nodes: %d\n', no_nodes_initial);
fprintf('Total new nodes added: %d\n', total_nodes_added);
fprintf('Final total nodes: %d\n', length(x_loc));
final_stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Final Coverage (k>=1): %.2f%%\n', final_stats.k_covered_percentage);
fprintf('Final Hole Area (k=0): %.2f%%\n', 100 - final_stats.k_covered_percentage);
fprintf('====================================================\n');

hold off;






% =========================================================================
% --- FINAL SUMMARY REPORT ---
% =========================================================================
fprintf('\n\n====================================================\n');
fprintf('--- FINAL RESTORATION SUMMARY (Circumcenter Method) ---\n');
fprintf('Algorithm stopped after %d restoration iteration(s).\n', iteration_count - 1);
no_nodes_initial = no_nodes;
total_nodes_added = length(x_loc) - no_nodes_initial;
fprintf('Total initial nodes: %d\n', no_nodes_initial);
fprintf('Total new nodes added: %d\n', total_nodes_added);
fprintf('Final total nodes: %d\n', length(x_loc));
final_stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
fprintf('Final Coverage (k>=1): %.2f%%\n', final_stats.k_covered_percentage);
fprintf('Final Hole Area (k=0): %.2f%%\n', 100 - final_stats.k_covered_percentage);
fprintf('====================================================\n');

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
% --- REVISED AND CORRECTED isArcCovered FUNCTION ---
% --- FINAL, CORRECT isArcCovered FUNCTION ---
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
% =========================================================================
% --- HELPER FUNCTION FOR THE BASELINE ALGORITHM ---
% =========================================================================
function hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs)
% Implements the hole detection method from the research paper.
% Identifies coverage holes by finding Delaunay triangles whose circumcircles
% are larger than the sensing range.

    hole_centers = [];
    if length(x_loc) < 3
        return;
    end

    fprintf('Finding holes using Delaunay circumcenter method...\n');
    
    % Perform Delaunay triangulation on ALL current nodes
    dt = delaunayTriangulation(x_loc', y_loc');
    
    % Get the circumcenter and circumradius for every triangle
    % The circumcenter function is built into MATLAB.
    [cc, rc] = circumcenter(dt);
    
    uncovered_centers = [];
    
    % Iterate through every triangle to check for holes
    for i = 1:size(cc, 1)
        % Condition from paper: A potential hole exists if the circumradius
        % is larger than the sensor's sensing radius.
        if rc(i) > Rs
            potential_hole_center = cc(i, :);
            
            % Now, we must check if this potential hole is already covered
            % by a nearby node (a crucial step for iterative algorithms).
            is_covered = false;
            for k = 1:length(x_loc)
                % If the distance from the potential hole center to any
                % existing node is less than Rs, it's already covered.
                if norm(potential_hole_center - [x_loc(k), y_loc(k)]) <= Rs
                    is_covered = true;
                    break; % This spot is covered, no need to check other nodes
                end
            end
            
            % If, after checking all nodes, the spot is still uncovered,
            % then it is a true hole.
            if ~is_covered
                uncovered_centers(end+1, :) = potential_hole_center;
            end
        end
    end
    
    % Remove duplicate hole centers that might be found
    if ~isempty(uncovered_centers)
        hole_centers = uniquetol(uncovered_centers, 1e-6, 'ByRows', true);
    end
    
    fprintf('Found %d potential holes to heal.\n', size(hole_centers, 1));
end

% =========================================================================
% --- HELPER FUNCTION FOR AREA-BASED BOUNDARY CHECK ---
% =========================================================================
function coverage_ratio = estimateCoverageInWSN(center, radius, net_length, net_width)
% Estimates the percentage of a circle's area that is inside the WSN
% boundaries using a Monte Carlo method.

    num_samples = 10000; % Increase for more accuracy, decrease for more speed.
    
    % Define the bounding box for the circle
    x_min = center(1) - radius;
    x_max = center(1) + radius;
    y_min = center(2) - radius;
    y_max = center(2) + radius;
    
    % Generate random points within the circle's bounding box
    rand_x = x_min + (x_max - x_min) * rand(num_samples, 1);
    rand_y = y_min + (y_max - y_min) * rand(num_samples, 1);
    
    % Calculate the squared distance from the center to each random point
    dist_sq = (rand_x - center(1)).^2 + (rand_y - center(2)).^2;
    
    % Find which points are actually inside the circle
    is_in_circle_mask = dist_sq <= radius^2;
    points_in_circle = sum(is_in_circle_mask);
    
    if points_in_circle == 0
        coverage_ratio = 0;
        return;
    end
    
    % Of the points inside the circle, find which are also inside the WSN
    in_circle_x = rand_x(is_in_circle_mask);
    in_circle_y = rand_y(is_in_circle_mask);
    
    is_in_wsn_mask = in_circle_x >= 0 & in_circle_x <= net_length & ...
                       in_circle_y >= 0 & in_circle_y <= net_width;
                       
    points_in_wsn = sum(is_in_wsn_mask);
    
    % Calculate the ratio
    coverage_ratio = points_in_wsn / points_in_circle;
end
