% =========================================================================
% CLASS: HORA_Algorithm.m (FINAL CORRECTED VERSION - TOOLBOX FREE)
% =========================================================================
classdef HORA_Algorithm
    methods (Static)
        function [restored_area, nodes_added_total] = run_hora_iterative(x_loc_init, y_loc_init, Rs, area_dims)
            net_length = area_dims(1); net_width = area_dims(2);
            x_loc = x_loc_init(:); y_loc = y_loc_init(:);
            Rc = 2 * Rs;
            
            stats_before = HORA_Algorithm.calculateKCoverage(x_loc_init, y_loc_init, Rs, net_length, net_width);
            
            TARGET_COVERAGE_PERCENT = 98.0; MAX_ITERATIONS = 15; iteration_count = 0;
            
            while true
                iteration_count = iteration_count + 1;
                stats_current = HORA_Algorithm.calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
                if stats_current.k_covered_percentage >= TARGET_COVERAGE_PERCENT || iteration_count > MAX_ITERATIONS, break; end
                
                all_holes = HORA_Algorithm.findAndClusterHoles_HORA(x_loc, y_loc, Rs, Rc, area_dims);
                if isempty(all_holes), break; end
                
                nodes_deployed_this_iteration = 0;
                for i = 1:length(all_holes)
                    boundary_nodes = all_holes{i};
                    if length(boundary_nodes) > 2
                        try
                            bnodes_coords = [x_loc(boundary_nodes), y_loc(boundary_nodes)];
                            if size(uniquetol(bnodes_coords, 1e-9, 'ByRows'), 1) < 3, continue; end
                            
                            hull_indices = convhull(bnodes_coords);
                            pgon = polyshape(bnodes_coords(hull_indices,:));
                            [target_x, target_y] = centroid(pgon);
                            
                            if ~HORA_Algorithm.isCoveredByAny([target_x, target_y], x_loc, y_loc, Rs)
                                x_loc(end+1) = target_x;
                                y_loc(end+1) = target_y;
                                nodes_deployed_this_iteration = nodes_deployed_this_iteration + 1;
                            end
                        catch
                            continue;
                        end
                    end
                end
                if nodes_deployed_this_iteration == 0, break; end
            end
            
            stats_after = HORA_Algorithm.calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
            restored_area = stats_before.hole_area - stats_after.hole_area;
            nodes_added_total = length(x_loc) - length(x_loc_init);
        end

        function all_holes = findAndClusterHoles_HORA(x, y, Rs, Rc, area_dims)
            num_nodes = length(x);
            Ni = cell(num_nodes,1);
            for i=1:num_nodes, for j=1:num_nodes, if i~=j && norm([x(i)-x(j),y(i)-y(j)])<2*Rs, Ni{i}=[Ni{i},j]; end,end,end
            boundary_nodes = [];
            for i = 1:num_nodes
                Pi = HORA_Algorithm.computeCP_and_CIP(i, x, y, Ni, Rs, area_dims(1), area_dims(2));
                if isempty(Ni{i}) || ~isempty(Pi), boundary_nodes(end+1) = i; end
            end
            if isempty(boundary_nodes), all_holes = {}; return; end
            boundary_coords = [x(boundary_nodes), y(boundary_nodes)];
            
            % THIS NOW CALLS THE CORRECT TOOLBOX-FREE FUNCTION
            cluster_indices = HORA_Algorithm.cluster_nodes_by_distance(boundary_coords, Rc);
            
            num_clusters = max(cluster_indices);
            if isempty(num_clusters) || num_clusters < 1, all_holes = {}; return; end
            all_holes = cell(1, num_clusters);
            for i = 1:num_clusters, all_holes{i} = boundary_nodes(cluster_indices == i); end
        end
        
        % --- THIS IS YOUR ORIGINAL, CORRECTED, TOOLBOX-FREE FUNCTION ---
        function cluster_indices = cluster_nodes_by_distance(coords, threshold)
            num_points = size(coords, 1);
            adj_matrix = false(num_points, num_points);
            for i = 1:num_points
                for j = i + 1:num_points
                    dist = sqrt(sum((coords(i,:) - coords(j,:)).^2));
                    if dist < threshold
                        adj_matrix(i, j) = true;
                        adj_matrix(j, i) = true;
                    end
                end
            end
            G = graph(adj_matrix, 'omitselfloops');
            cluster_indices = conncomp(G)';
        end

        % --- ALL OTHER HELPER FUNCTIONS ARE NOW INSIDE THE CLASS ---
        function stats=calculateKCoverage(x,y,R,L,W), res=10;xg=0:res:L;yg=0:res:W;[X,Y]=meshgrid(xg,yg);cm=zeros(size(X));for i=1:length(x),cm=cm+(((X-x(i)).^2+(Y-y(i)).^2)<=R^2);end;pa=res^2;ta=L*W;stats.hole_area=sum(cm(:)==0)*pa;stats.k_covered_area=ta-stats.hole_area;stats.k_covered_percentage=(stats.k_covered_area/ta)*100;end
        function c=isCoveredByAny(p,x,y,R), c=any(sqrt(sum((p-[x,y]).^2,2))<=R);end
        function Pi=computeCP_and_CIP(i,x,y,Ni,R,L,W),Pi=[];xi=x(i);yi=y(i);if HORA_Algorithm.isBorderSensor(xi,yi,L,W,R),b_pts=HORA_Algorithm.borderIntersections(xi,yi,R,L,W);for p_idx=1:size(b_pts,1),if~HORA_Algorithm.isCoveredByAnyNeighbor(b_pts(p_idx,:),i,Ni{i},x,y,R),Pi=[Pi;b_pts(p_idx,:)];end,end,end;for j=Ni{i},i_pts=HORA_Algorithm.circleIntersections(xi,yi,x(j),y(j),R);for p_idx=1:size(i_pts,1),p=i_pts(p_idx,:);if~HORA_Algorithm.isCoveredByAnotherSpecificNeighbor(p,i,j,Ni{i},x,y,R),Pi=[Pi;p];end,end,end;if~isempty(Pi),Pi=uniquetol(Pi,1e-9,'ByRows',true);end,end
        function isB=isBorderSensor(x,y,L,W,R),isB=(x-R<=0)||(x+R>=L)||(y-R<=0)||(y+R>=W);end
        function pts=borderIntersections(x,y,r,L,W),pts=[];if(x-r<=0),dy=sqrt(max(0,r^2-x^2));pts=[pts;[0,y-dy];[0,y+dy]];end;if(x+r>=L),dy=sqrt(max(0,r^2-(L-x)^2));pts=[pts;[L,y-dy];[L,y+dy]];end;if(y-r<=0),dx=sqrt(max(0,r^2-y^2));pts=[pts;[x-dx,0];[x+dx,0]];end;if(y+r>=W),dx=sqrt(max(0,r^2-(W-y)^2));pts=[pts;[x-dx,W];[x+dx,W]];end;if~isempty(pts),pts=pts(all(pts>=0&pts<=[L,W],2),:);pts=uniquetol(pts,1e-9,'ByRows',true);end,end
        function c=isCoveredByAnyNeighbor(p,i,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
        function c=isCoveredByAnotherSpecificNeighbor(p,i,j,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&k~=j&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
        function pts=circleIntersections(x1,y1,x2,y2,r),d=norm([x1-x2,y1-y2]);if d>=2*r||d==0,pts=[];return;end;a=d/2;h=sqrt(r^2-a^2);xm=x1+a*(x2-x1)/d;ym=y1+a*(y2-y1)/d;pts=[xm-h*(y2-y1)/d,ym+h*(x2-x1)/d;xm+h*(y2-y1)/d,ym-h*(x2-x1)/d];end
    end
end