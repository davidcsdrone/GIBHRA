% FILE NAME: isCoveredByAny.m

function covered = isCoveredByAny(point, x_loc, y_loc, Rs)


    covered = false; % Assume not covered by default

    for i = 1:length(x_loc)
        % calculate the squared distance from the point to the center of node i
        dist_sq = (point(1) - x_loc(i))^2 + (point(2) - y_loc(i))^2;

        % if the squared distance is less than or equal to the radius squared, it's covered.
        % we use squared values to avoid a costly square root operation inside the loop.
        if dist_sq <= Rs^2
            covered = true;
            return; % exit the function early since we found a covering node
        end
    end
end
