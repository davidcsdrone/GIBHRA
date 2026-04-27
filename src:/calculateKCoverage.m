function stats = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width)
    grid_resolution = 10; % using a slightly coarser grid also speeds things up
    x_grid = 0:grid_resolution:net_length;
    y_grid = 0:grid_resolution:net_width;
    [X, Y] = meshgrid(x_grid, y_grid);
    coverage_map = zeros(size(X));

    for i = 1:length(x_loc)
        % we can optimize: create a smaller bounding box for each node
        x_min_idx = find(x_grid >= x_loc(i) - Rs, 1, 'first');
        x_max_idx = find(x_grid <= x_loc(i) + Rs, 1, 'last');
        y_min_idx = find(y_grid >= y_loc(i) - Rs, 1, 'first');
        y_max_idx = find(y_grid <= y_loc(i) + Rs, 1, 'last');
        
        if isempty(x_min_idx) || isempty(x_max_idx) || isempty(y_min_idx) || isempty(y_max_idx)
            continue;
        end
        
        % create a small sub-grid for this node's potential coverage area
        [sub_X, sub_Y] = meshgrid(x_grid(x_min_idx:x_max_idx), y_grid(y_min_idx:y_max_idx));
        
        % perform the expensive calculation ONLY on the small sub-grid
        dist_sq = (sub_X - x_loc(i)).^2 + (sub_Y - y_loc(i)).^2;
        
        % add the result to the corresponding part of the main coverage map
        coverage_map(y_min_idx:y_max_idx, x_min_idx:x_max_idx) = ...
            coverage_map(y_min_idx:y_max_idx, x_min_idx:x_max_idx) + (dist_sq <= Rs^2);
    end
    
    total_grid_points = numel(X); % Total number of grid points
    point_area = grid_resolution^2;
    stats.total_area = net_length * net_width;

    % calculate hole area (k=0 coverage)
    stats.hole_area = sum(coverage_map(:) == 0) * point_area;
    stats.hole_percentage = (stats.hole_area / stats.total_area) * 100; % Corrected percentage calculation

    % calculate k>=1 covered area
    stats.k_covered_area = sum(coverage_map(:) >= 1) * point_area;
    stats.k_covered_percentage = (stats.k_covered_area / stats.total_area) * 100;

    % --- this will calculate specific k-coverage areas and percentages ---
    max_k = max(coverage_map(:));
    stats.specific_k_areas = zeros(1, max_k);
    stats.specific_k_percentages = zeros(1, max_k);

    for k_val = 1:max_k
        stats.specific_k_areas(k_val) = sum(coverage_map(:) == k_val) * point_area;
        stats.specific_k_percentages(k_val) = (stats.specific_k_areas(k_val) / stats.total_area) * 100;
    end
end