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