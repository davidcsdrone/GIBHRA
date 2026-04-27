% FILE: solveApolloniusForTriplet.m (Corrected and better more robust version)

function [optimal_center, solution_radius] = solveApolloniusForTriplet(triplet_ids, x_loc, y_loc, Rs)
% solves the Apollonius problem for the small circle tangent to three other
% circles of the same radius Rs.

    % --- Step 1: Get Circle Data ---
    c1 = struct('x', x_loc(triplet_ids(1)), 'y', y_loc(triplet_ids(1)), 'r', Rs);
    c2 = struct('x', x_loc(triplet_ids(2)), 'y', y_loc(triplet_ids(2)), 'r', Rs);
    c3 = struct('x', x_loc(triplet_ids(3)), 'y', y_loc(triplet_ids(3)), 'r', Rs);
    
  
    A = 2 * [c2.x - c1.x, c2.y - c1.y;
             c3.x - c1.x, c3.y - c1.y];
         
    B = [c2.x^2 + c2.y^2 - c1.x^2 - c1.y^2;
         c3.x^2 + c3.y^2 - c1.x^2 - c1.y^2];
    
    % --- Step 3: solve for the center of the Solution Circle ---
    % Check if the matrix A is invertible. If not, the three points are collinear.
    if rank(A) < 2
        optimal_center = [];
        solution_radius = [];
        % fprintf('Debug: Triplet nodes are collinear. Cannot solve.\n');
        return;
    end
    
    % Solve the linear system A * center = B
    center = A \ B;
    optimal_center = center'; % Transpose to get a 1x2 vector [x, y]
    
    % --- Step 4: Solve for the Radius of the Solution Circle ---
    % Now that we have the center, we can find the radius.
    % The distance from the optimal center to any of the three circle centers
    % should be equal to (Rs + solution_radius).
    dist_to_c1 = norm(optimal_center - [c1.x, c1.y]);
    
    % The radius of the new node is the distance minus the existing node's radius.
    solution_radius = dist_to_c1 - Rs;
    
    % --- Step 5: Final Sanity Check ---
    % The solution radius must be positive.
    if solution_radius <= 0
        optimal_center = [];
        solution_radius = [];
        % fprintf('Debug: Solver resulted in a non-positive radius.\n');
    end
end