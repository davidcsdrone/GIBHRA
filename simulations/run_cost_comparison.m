% =========================================================================
% SCRIPT: run_final_simulation.m (Final Corrected and Complete Version)
% DESCRIPTION: Runs and compares FIVE hole restoration algorithms.
% =========================================================================

% --- Main Simulation Setup ---
clear; clc; close all;
rng('default'); % For repeatable random results

% --- Simulation Parameters ---
deployment_area = [1000, 1000];
sensing_radius = 50;
node_counts = [200, 400, 600];
num_runs = 20;

% --- Algorithm Setup ---
algorithms = {
    'GIBHRA', @run_inversion_algorithm;
    'HCHA',   @run_circumcenter_algorithm;
    'HPA',    @run_yao_algorithm;
    'NLCHR',  @run_nlchr_algorithm
};

% --- Data Storage ---
final_avg_restored_area = zeros(length(node_counts), size(algorithms, 1));
final_avg_nodes_added = zeros(length(node_counts), size(algorithms, 1));

% --- Main Simulation Loop ---
for i = 1:length(node_counts)
    num_nodes = node_counts(i);
    fprintf('--- Simulating with %d Initial Nodes ---\n', num_nodes);
    
    run_results_area = zeros(num_runs, size(algorithms, 1));
    run_results_nodes = zeros(num_runs, size(algorithms, 1));

    for j = 1:num_runs
        fprintf('  Run %d/%d...\n', j, num_runs);
        
        [x_loc_initial, y_loc_initial] = generate_network(deployment_area, num_nodes);
        
        for k = 1:size(algorithms, 1)
            alg_name = algorithms{k, 1};
            alg_function = algorithms{k, 2};
            
            fprintf('    Testing: %s...\n', alg_name);
            [total_restored, nodes_added] = alg_function(x_loc_initial, y_loc_initial, sensing_radius, deployment_area);
            
            if nodes_added > 0
                run_results_area(j, k) = total_restored / nodes_added;
            else
                run_results_area(j, k) = 0;
            end
            run_results_nodes(j, k) = nodes_added;
        end
    end
    
    final_avg_restored_area(i, :) = mean(run_results_area, 1);
    final_avg_nodes_added(i, :) = mean(run_results_nodes, 1);
    fprintf('  Avg Restored Area per Node: %s\n', mat2str(final_avg_restored_area(i, :), 4));
    fprintf('  Avg Nodes Added: %s\n\n', mat2str(final_avg_nodes_added(i, :), 3));
end

fprintf('\n--- Simulation Complete ---\n');

% --- Plotting the Final Bar Charts ---
plot_metric_barchart(node_counts, final_avg_restored_area, algorithms(:,1), 'Avg. Restored Area per Node (m^2)', 'barchart_restored_area');
plot_metric_barchart(node_counts, final_avg_nodes_added, algorithms(:,1), 'Avg. Number of Restorative Nodes Added', 'barchart_deployment_cost');


% =========================================================================
% --- UNIVERSAL PLOTTING FUNCTION (Updated for 5 algorithms) ---
% =========================================================================
function plot_metric_barchart(x_data, y_data, legend_names, y_label_text, filename)
    figure('Name', y_label_text);
    b = bar(x_data, y_data, 'grouped');
    
    colors = {[0.8,0,0], [0,0.4,0.8], [0.1,0.5,0.1], [0.6,0.2,0.8], [1,0.6,0]};
    for i = 1:length(b)
        b(i).FaceColor = colors{i};
    end
    
    xlabel('Number of Deployed Nodes', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(y_label_text, 'FontSize', 12, 'FontWeight', 'bold');
    
    legend(legend_names, 'Location', 'northwest', 'FontSize', 11);
    
    ax = gca; ax.FontSize = 11; ax.YGrid = 'on'; ax.Box = 'on';
    ylim([0, max([1, y_data(:)']) * 1.15]);
    
    print(filename, '-depsc', '-r300');
    fprintf('\nChart saved as %s.eps\n', filename);
end

% =========================================================================
% --- ALGORITHM WRAPPER FUNCTIONS (All 5 included) ---
% =========================================================================
function [restored_area, nodes_added_total] = run_inversion_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:); % Ensure column vectors
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 10; iteration_count = 0;
    while true
       iteration_count = iteration_count + 1;
       stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
       if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end
       
       boundary_nodes_list = findBoundaryNodes(x_loc, y_loc, Rs, net_length, net_width);
       T = selectTripletsDelaunay(boundary_nodes_list, x_loc, y_loc, Rs, stats_current.k_covered_percentage);
       if isempty(T), break; end
       
       nodes_deployed_this_iteration = 0;
       for t_idx = 1:length(T)
           [optimal_center, ~] = solveApolloniusForTriplet(T{t_idx}, x_loc, y_loc, Rs);
           if ~isempty(optimal_center) && estimateCoverageInWSN(optimal_center, Rs, net_length, net_width) >= 0.5
               x_loc(end+1,1) = optimal_center(1); 
               y_loc(end+1,1) = optimal_center(2);
               nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
           end
       end
       if nodes_deployed_this_iteration == 0, break; end
    end
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

function [restored_area, nodes_added_total] = run_circumcenter_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 10; iteration_count = 0;
    while true
        iteration_count = iteration_count + 1;
        stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
        if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end
        
        hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs);
        if isempty(hole_centers), break; end
        
        nodes_deployed_this_iteration = 0;
        for i = 1:size(hole_centers, 1)
            if estimateCoverageInWSN(hole_centers(i, :), Rs, net_length, net_width) >= 0.5
                x_loc(end+1,1) = hole_centers(i, 1); 
                y_loc(end+1,1) = hole_centers(i, 2);
                nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
            end
        end
        if nodes_deployed_this_iteration == 0, break; end
    end
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

function [restored_area, nodes_added_total] = run_yao_algorithm(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    Rc = 2 * Rs;
    stats_before = calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
    TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 20; iteration_count = 0;
    while true
       iteration_count = iteration_count + 1;
       stats_current = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
       if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end
       
       all_hbn = [];
       for i = 1:length(x_loc), if findHoleBoundaryNode_Yao(i, x_loc, y_loc, Rc), all_hbn(end+1) = i; end, end
       if isempty(all_hbn), break; end
       
       boundary_coords = [x_loc(all_hbn), y_loc(all_hbn)];
       if size(boundary_coords, 1) < 3, break; end
       
       cluster_indices = cluster_nodes_by_distance(boundary_coords, Rc);
       num_clusters = max(cluster_indices);
       if isempty(num_clusters) || num_clusters < 1, break; end
       
       holes = cell(1, num_clusters);
       for i = 1:num_clusters, holes{i} = all_hbn(cluster_indices == i); end
       
       hole_areas = zeros(1, num_clusters);
for i = 1:num_clusters
    current_hole_nodes = holes{i};
    if length(current_hole_nodes) >= 3
        % --- START OF THE ROBUSTNESS CHECK ---
        hole_coords = [x_loc(current_hole_nodes), y_loc(current_hole_nodes)];
        
        % First, get only the unique points using a small tolerance
        unique_hole_coords = uniquetol(hole_coords, 1e-9, 'ByRows', true);
        
        % Only call convhull if there are at least 3 unique points
        if size(unique_hole_coords, 1) >= 3
            try
                % We still wrap in a try-catch in case of collinearity
                [~, hole_areas(i)] = convhull(unique_hole_coords);
            catch
                % If convhull still fails, the points are collinear.
                % Assign zero area and safely continue.
                hole_areas(i) = 0;
            end
        else
            % Not enough unique points to form a polygon.
            hole_areas(i) = 0;
        end
        % --- END OF THE ROBUSTNESS CHECK ---
    end
end
       
       [max_area, idx] = max(hole_areas);
       if max_area == 0, break; end
       
       [x_loc, y_loc, nodes_added_iter] = patchHole_Yao(holes{idx}, x_loc, y_loc, Rc, Rs, net_length, net_width);
       if nodes_added_iter == 0, break; end
    end
    stats_after = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    restored_area = stats_before.hole_area - stats_after.hole_area;
    nodes_added_total = length(x_loc) - length(x_loc_init);
end

% --- REVISED NLCHR WRAPPER ---
% --- REVISED NLCHR WRAPPER ---
% --- REVISED AND IMPROVED NLCHR WRAPPER ---
% --- FINAL, CORRECTED NLCHR WRAPPER ---
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


% =========================================================================
% --- CONSOLIDATED HELPER FUNCTIONS ---
% =========================================================================
function [x, y] = generate_network(area, num_nodes)
    x = area(1) * rand(num_nodes, 1); % Creates an N x 1 column vector
    y = area(2) * rand(num_nodes, 1); % Creates an N x 1 column vector
end
function stats = calculateKCoverage(x,y,R,L,W), res=10;xg=0:res:L;yg=0:res:W;[X,Y]=meshgrid(xg,yg);cm=zeros(size(X));for i=1:length(x),cm=cm+(((X-x(i)).^2+(Y-y(i)).^2)<=R^2);end;pa=res^2;ta=L*W;stats.hole_area=sum(cm(:)==0)*pa;stats.k_covered_area=ta-stats.hole_area;stats.k_covered_percentage=(stats.k_covered_area/ta)*100;end
function boundary_nodes = findBoundaryNodes(x,y,R,L,W)
    num_nodes = length(x); Ni = cell(num_nodes,1);
    for i=1:num_nodes, for j=1:num_nodes, if i~=j && norm([x(i)-x(j),y(i)-y(j)])<2*R, Ni{i}=[Ni{i},j]; end,end,end
    boundary_nodes = [];
    for i=1:num_nodes, [Pi,~,~]=computeCP_and_CIP(i,x,y,Ni,R,L,W); if isempty(Ni{i})||~isempty(Pi), boundary_nodes(end+1)=i; end, end
end
function T = selectTripletsDelaunay(B_nodes, x_loc, y_loc, Rs, current_coverage_percent)
    if current_coverage_percent < 60.0, distance_threshold = 1.8 * Rs;
    elseif current_coverage_percent < 85.0, distance_threshold = 1.6 * Rs;
    else, distance_threshold = 1.4 * Rs; end
    T = {};
    if length(B_nodes) < 3, return; end
    boundary_x = x_loc(B_nodes); boundary_y = y_loc(B_nodes);
    dt = delaunayTriangulation(boundary_x, boundary_y);
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
function [pc, pr] = solveApolloniusForTriplet(t,x,y,r),c(1)=struct('x',x(t(1)),'y',y(t(1)),'r',r);c(2)=struct('x',x(t(2)),'y',y(t(2)),'r',r);c(3)=struct('x',x(t(3)),'y',y(t(3)),'r',r);p1=struct('x',c(1).x,'y',c(1).y);c2r=struct('x',c(2).x,'y',c(2).y,'r',0);c3r=struct('x',c(3).x,'y',c(3).y,'r',0);sols=solve_pcc_by_inversion(p1,c2r,c3r);if isempty(sols),pc=[];pr=[];return;end;[~,idx]=max([sols.r]);sol=sols(idx);pc=[sol.x,sol.y];pr=sol.r+r;end
function sols=solve_pcc_by_inversion(p,c1,c2),ic=p;ir=1;c1i=invert_circle(c1,ic,ir);c2i=invert_circle(c2,ic,ir);tl=find_common_tangents(c1i,c2i);sols=struct('x',{},'y',{},'r',{});for i=1:length(tl),if isstruct(tl(i)),sols(end+1)=invert_line(tl(i),ic,ir);end,end,end
function ic=invert_circle(c,center,k),d2=(c.x-center.x)^2+(c.y-center.y)^2;s=k^2/(d2-c.r^2);ic=struct('x',center.x+s*(c.x-center.x),'y',center.y+s*(c.y-center.y),'r',abs(s)*c.r);end
function ic=invert_line(l,c,k),d=2*(l.c+l.a*c.x+l.b*c.y);ic=struct('x',c.x-k^2*l.a/d,'y',c.y-k^2*l.b/d,'r',abs(k^2/d));end
function t=find_common_tangents(c1,c2),t=cell(1,4);vx=c2.x-c1.x;vy=c2.y-c1.y;d2=vx^2+vy^2;dr=c2.r-c1.r;if d2>=dr^2,d=sqrt(d2);rt=sqrt(d2-dr^2);a=(vx*dr+vy*rt)/d2;b=(vy*dr-vx*rt)/d2;t{1}=struct('a',a,'b',b,'c',c1.r-(a*c1.x+b*c1.y));a=(vx*dr-vy*rt)/d2;b=(vy*dr+vx*rt)/d2;t{2}=struct('a',a,'b',b,'c',c1.r-(a*c1.x+b*c1.y));end;sr=c1.r+c2.r;if d2>=sr^2,d=sqrt(d2);rt=sqrt(d2-sr^2);a=(vx*sr+vy*rt)/d2;b=(vy*sr-vx*rt)/d2;t{3}=struct('a',a,'b',b,'c',-c1.r-(a*c1.x+b*c1.y));a=(vx*sr-vy*rt)/d2;b=(vy*sr+vx*rt)/d2;t{4}=struct('a',a,'b',b,'c',-c1.r-(a*c1.x+b*c1.y));end;t=[t{:}];end
% --- REPLACE this entire function ---
function c = isCoveredByAny(p, x, y, R)
    % Calculates squared distances from the point 'p' to all nodes.
    % This version correctly handles column vectors for x and y.
    dist_sq = (x(:) - p(1)).^2 + (y(:) - p(2)).^2;
    % Checks if any squared distance is less than or equal to Rs^2.
    c = any(dist_sq <= R^2);
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
function hc=findHoles_Circumcenter(x,y,R),hc=[];if length(x)<3,return;end;[cc,rc]=circumcenter(delaunayTriangulation(x(:),y(:)));uc=[];for i=1:size(cc,1),if rc(i)>R&&~isCoveredByAny(cc(i,:),x,y,R),uc(end+1,:)=cc(i,:);end,end;if~isempty(uc),hc=uniquetol(uc,1e-6,'ByRows',true);end,end
function is_hbn=findHoleBoundaryNode_Yao(id,x,y,Rc),N_X=find(sqrt(sum(( [x,y] - [x(id),y(id)] ).^2,2))<Rc&(1:length(x))'~=id);if length(N_X)<2,is_hbn=true;return;end;sub_ids=[id;N_X];try,k_b=boundary(x(sub_ids),y(sub_ids),1);catch,is_hbn=true;return;end;is_hbn=ismember(id,sub_ids(unique(k_b)));end
function [x_new,y_new,added]=patchHole_Yao(h_nodes,x,y,Rc,~,L,W),x_new=x;y_new=y;added=0;if length(h_nodes)<2,return;end;try,b_path=boundary(x(h_nodes),y(h_nodes),1);catch,return;end;for i=1:(length(b_path)-1),pA=[x(h_nodes(b_path(i))),y(h_nodes(b_path(i)))];pB=[x(h_nodes(b_path(i+1))),y(h_nodes(b_path(i+1)))];mid=(pA+pB)/2;v_AB=pB-pA;pv_norm=[-v_AB(2),v_AB(1)]/norm(v_AB);d_A_mid=norm(pA-mid);if Rc^2<d_A_mid^2,continue;end;h=sqrt(Rc^2-d_A_mid^2);p1=mid+h*pv_norm;p2=mid-h*pv_norm;if norm(p1-mean([x,y]))>norm(p2-mean([x,y])),f_pos=p1;else,f_pos=p2;end;if all(f_pos>=0 & f_pos<=[L,W]),x_new(end+1)=f_pos(1);y_new(end+1)=f_pos(2);added=added+1;end,end,end
function C=cluster_nodes_by_distance(coords,th),N=size(coords,1);adj=false(N);for i=1:N,for j=i+1:N,if norm(coords(i,:)-coords(j,:))<th,adj(i,j)=true;adj(j,i)=true;end,end,end;C=conncomp(graph(adj,'omitselfloops'))';end
function cr=estimateCoverageInWSN(c,r,L,W),N=10000;p=rand(N,2).*[2*r,2*r]+(c-r);ip=p(sum((p-c).^2,2)<=r^2,:);if isempty(ip),cr=0;return;end;cr=sum(all(ip>=0&ip<=[L,W],2))/size(ip,1);end
function pts=borderIntersections(x,y,r,L,W),pts=[];if(x-r<=0),dy=sqrt(max(0,r^2-x^2));pts=[pts;[0,y-dy];[0,y+dy]];end;if(x+r>=L),dy=sqrt(max(0,r^2-(L-x)^2));pts=[pts;[L,y-dy];[L,y+dy]];end;if(y-r<=0),dx=sqrt(max(0,r^2-y^2));pts=[pts;[x-dx,0];[x+dx,0]];end;if(y+r>=W),dx=sqrt(max(0,r^2-(W-y)^2));pts=[pts;[x-dx,W];[x+dx,W]];end;if~isempty(pts),pts=pts(all(pts>=0&pts<=[L,W],2),:);pts=uniquetol(pts,1e-9,'ByRows',true);end,end
function isB=isBorderSensor(x,y,L,W,R),isB=(x-R<=0)||(x+R>=L)||(y-R<=0)||(y+R>=W);end
function c=isCoveredByAnyNeighbor(p,i,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
function c=isCoveredByAnotherSpecificNeighbor(p,i,j,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&k~=j&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
function all_holes=findAndClusterHoles_HORA(x,y,~,Rc)
    num_nodes=length(x);boundary_nodes=[];
    for i=1:num_nodes
        neighbors=find(sqrt((x-x(i)).^2+(y-y(i)).^2)<Rc & (1:num_nodes)'~=i);
        if isempty(neighbors)||length(neighbors)<2, boundary_nodes(end+1)=i; continue; end
        try
            neighbor_coords=[x(neighbors),y(neighbors)];
            if size(uniquetol(neighbor_coords,1e-9,'ByRows'),1)<3, boundary_nodes(end+1)=i; continue; end
            hull_indices=convhull(neighbor_coords);
            if length(unique(hull_indices))<length(neighbors), boundary_nodes(end+1)=i; end
        catch
            boundary_nodes(end+1)=i;
        end
    end
    if isempty(boundary_nodes),all_holes={};return;end
    cluster_indices=cluster_nodes_by_distance([x(boundary_nodes)',y(boundary_nodes)'],Rc);
    num_clusters=max(cluster_indices);
    if isempty(num_clusters)||num_clusters==0,all_holes={};return;end
    all_holes=cell(1,num_clusters);
    for i=1:num_clusters,all_holes{i}=boundary_nodes(cluster_indices==i);end
end
function all_holes=findAndClusterHoles_NLCHR(x,y,R,dims),all_holes={};boundary_nodes=findBoundaryNodes(x,y,R,dims(1),dims(2));if isempty(boundary_nodes),return;end;cluster_indices=cluster_nodes_by_distance([x(boundary_nodes)',y(boundary_nodes)'],2*R);num_clusters=max(cluster_indices);if isempty(num_clusters)||num_clusters==0,return;end;all_holes=cell(1,num_clusters);for i=1:num_clusters,all_holes{i}=boundary_nodes(cluster_indices==i);end,end
% --- REVISED ISOSCELES VERTEX HELPER ---

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
function s=isAboveLine(px,py,x1,y1,x2,y2),v=(x2-x1)*(py-y1)-(y2-y1)*(px-x1);if abs(v)<1e-9,s='on';elseif v>0,s='above';else,s='below';end,end
function pts=circleIntersections(x1,y1,x2,y2,r),d=norm([x1-x2,y1-y2]);if d>=2*r||d==0,pts=[];return;end;a=d/2;h=sqrt(r^2-a^2);xm=x1+a*(x2-x1)/d;ym=y1+a*(y2-y1)/d;pts=[xm-h*(y2-y1)/d,ym+h*(x2-x1)/d;xm+h*(y2-y1)/d,ym-h*(x2-x1)/d];end



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