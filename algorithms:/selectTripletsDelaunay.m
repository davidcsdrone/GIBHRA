function T = selectTripletsDelaunay(B_nodes, x_loc, y_loc, Rs, current_coverage_percent)
    % determines the distance threshold based on current network coverage
    if current_coverage_percent < 60.0
        distance_threshold = 1.8 * Rs;
    elseif current_coverage_percent < 85.0
        distance_threshold = 1.6 * Rs;
    else
        distance_threshold = 1.4 * Rs;
    end
    
    T = {}; % Initialize the output list of triplets
    if length(B_nodes) < 3, return; end
    
    % Perform Delaunay triangulation on only the boundary nodes
    boundary_x = x_loc(B_nodes); 
    boundary_y = y_loc(B_nodes);
    dt = delaunayTriangulation(boundary_x', boundary_y');
    candidate_triangles = dt.ConnectivityList;
    
    valid_triplets_info = []; % To store [node1, node2, node3, perimeter]
    
    % check each triangle from the triangulation
    for i = 1:size(candidate_triangles, 1)
        % get the global IDs of the nodes in the current triangle
        node_ids = B_nodes(candidate_triangles(i, :));
        
        % get coordinates and distances for the three nodes
        p1 = [x_loc(node_ids(1)), y_loc(node_ids(1))]; 
        p2 = [x_loc(node_ids(2)), y_loc(node_ids(2))]; 
        p3 = [x_loc(node_ids(3)), y_loc(node_ids(3))];
        
        d12 = norm(p1 - p2); 
        d23 = norm(p2 - p3); 
        d31 = norm(p3 - p1);
        
        % this is the heuristic non-overlap check from your original algorithm
        if (d12 > distance_threshold) && (d23 > distance_threshold) && (d31 > distance_threshold)
            
            [potential_center, ~] = solveApolloniusForTriplet(node_ids, x_loc, y_loc, Rs);
            
            % check if the potential placement is not already covered
            if ~isempty(potential_center) && ~isCoveredByAny(potential_center, x_loc, y_loc, Rs)
                perimeter = d12 + d23 + d31;
                valid_triplets_info(end+1, :) = [node_ids, perimeter];
            end
        end
    end
    
    if ~isempty(valid_triplets_info)
        % sort by the largest perimeter first to prioritize the biggest holes
        sorted_triplets_info = sortrows(valid_triplets_info, 4, 'descend'); 
        for i = 1:size(sorted_triplets_info, 1)
            T{end+1} = sorted_triplets_info(i, 1:3); 
        end
    end
end