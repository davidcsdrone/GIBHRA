% =========================================================================
% SCRIPT: run_efficiency_simulation.m (Final Corrected and Complete Version)
% DESCRIPTION: Simulates and compares the placement quality of four hole
%              restoration algorithms.
% METRIC:      Coverage Efficiency (%)
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
    'GIBHRA', @run_inversion_algorithm_for_efficiency;
    'HCHA',   @run_circumcenter_algorithm_for_efficiency;
    'HPA',    @run_yao_algorithm_for_efficiency;
    'NLCHR',  @run_nlchr_algorithm_for_efficiency
};

% --- Data Storage ---
final_avg_efficiency = zeros(length(node_counts), size(algorithms, 1));

% --- Main Simulation Loop ---
for i = 1:length(node_counts)
    num_nodes = node_counts(i);
    fprintf('--- Simulating Efficiency for %d Initial Nodes ---\n', num_nodes);
    
    run_results_matrix = zeros(num_runs, size(algorithms, 1));

    for j = 1:num_runs
        fprintf('  Run %d/%d...\n', j, num_runs);
        
        [x_loc_initial, y_loc_initial] = generate_network(deployment_area, num_nodes);
        
        for k = 1:size(algorithms, 1)
            alg_name = algorithms{k, 1};
            alg_function = algorithms{k, 2};
            
            fprintf('    Testing: %s...\n', alg_name);
            avg_efficiency_for_run = alg_function(x_loc_initial, y_loc_initial, sensing_radius, deployment_area);
            
            run_results_matrix(j, k) = avg_efficiency_for_run;
        end
    end
    
    final_avg_efficiency(i, :) = mean(run_results_matrix, 1, 'omitnan');
    fprintf('  Average Coverage Efficiency for %d nodes: %s %%\n\n', num_nodes, mat2str(final_avg_efficiency(i, :), 3));
end

fprintf('\n--- Simulation Complete ---\n');

% --- Plotting the Final Bar Chart ---
plot_efficiency_barchart(node_counts, final_avg_efficiency, algorithms(:,1));


% =========================================================================
% --- UNIVERSAL PLOTTING FUNCTION ---
% =========================================================================
function plot_efficiency_barchart(x_data, y_data, legend_names)
    figure('Name', 'Algorithm Comparison: Coverage Efficiency');
    b = bar(x_data, y_data, 'grouped');
    
    colors = {[0.8,0,0], [0,0.4,0.8], [0.1,0.5,0.1], [1,0.6,0]}; % Red, Blue, Green, Orange
    for i = 1:length(b)
        b(i).FaceColor = colors{i};
    end
    
    xlabel('Number of Deployed Nodes', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Coverage Efficiency (%)', 'FontSize', 12, 'FontWeight', 'bold');
    
    legend(legend_names, 'Location', 'northeast', 'FontSize', 11);
    
    ax = gca; ax.FontSize = 11; ax.YGrid = 'on'; ax.Box = 'on';
    ylim([0, 100]);
    
    print('barchart_coverage_efficiency', '-depsc', '-r300');
    fprintf('\nCoverage efficiency chart saved as barchart_coverage_efficiency.eps\n');
end

% =========================================================================
% --- ALGORITHM WRAPPERS FOR EFFICIENCY METRIC (Corrected) ---
% =========================================================================
function avg_efficiency = run_inversion_algorithm_for_efficiency(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    all_efficiencies = []; max_area_per_node = pi * Rs^2;

    stats_before_iter = calculateKCoverage(x_loc, y_loc, Rs, net_length, net_width);
    
    boundary_nodes_list = findBoundaryNodes(x_loc, y_loc, Rs, net_length, net_width);
    T = selectTripletsDelaunay(boundary_nodes_list, x_loc, y_loc, Rs, stats_before_iter.k_covered_percentage);
    
    if isempty(T), avg_efficiency = 0; return; end
    
    x_loc_current_iter = x_loc; y_loc_current_iter = y_loc;
    
    for t_idx = 1:length(T)
        [optimal_center, ~] = solveApolloniusForTriplet(T{t_idx}, x_loc, y_loc, Rs);
        if ~isempty(optimal_center) && estimateCoverageInWSN(optimal_center, Rs, net_length, net_width) >= 0.5
            if isCoveredByAny(optimal_center, x_loc_current_iter, y_loc_current_iter, Rs), continue; end

            stats_before_this_node = calculateKCoverage(x_loc_current_iter, y_loc_current_iter, Rs, net_length, net_width);
            stats_after_this_node = calculateKCoverage([x_loc_current_iter; optimal_center(1)], [y_loc_current_iter; optimal_center(2)], Rs, net_length, net_width);
            
            area_restored = max(0, stats_before_this_node.hole_area - stats_after_this_node.hole_area);
            efficiency = (area_restored / max_area_per_node) * 100;
            all_efficiencies(end+1) = efficiency;
            
            x_loc_current_iter(end+1) = optimal_center(1);
            y_loc_current_iter(end+1) = optimal_center(2);
        end
    end
    if isempty(all_efficiencies), avg_efficiency = 0; else, avg_efficiency = mean(all_efficiencies, 'omitnan'); end
end

function avg_efficiency = run_circumcenter_algorithm_for_efficiency(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    all_efficiencies = []; max_area_per_node = pi * Rs^2;

    hole_centers = findHoles_Circumcenter(x_loc, y_loc, Rs);
    if isempty(hole_centers), avg_efficiency = 0; return; end
    
    x_loc_current_iter = x_loc;
    y_loc_current_iter = y_loc;

    for i = 1:size(hole_centers, 1)
        potential_center = hole_centers(i, :);
        if estimateCoverageInWSN(potential_center, Rs, net_length, net_width) >= 0.5
            if isCoveredByAny(potential_center, x_loc_current_iter, y_loc_current_iter, Rs), continue; end

            stats_before_this_node = calculateKCoverage(x_loc_current_iter, y_loc_current_iter, Rs, net_length, net_width);
            stats_after_one_node = calculateKCoverage([x_loc_current_iter; potential_center(1)], [y_loc_current_iter; potential_center(2)], Rs, net_length, net_width);
            area_restored = max(0, stats_before_this_node.hole_area - stats_after_one_node.hole_area);
            efficiency = (area_restored / max_area_per_node) * 100;
            all_efficiencies(end+1) = efficiency;
            
            x_loc_current_iter(end+1) = potential_center(1);
            y_loc_current_iter(end+1) = potential_center(2);
        end
    end
    if isempty(all_efficiencies), avg_efficiency = 0; else, avg_efficiency = mean(all_efficiencies, 'omitnan'); end
end

function avg_efficiency = run_yao_algorithm_for_efficiency(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    Rc = 2 * Rs;
    all_efficiencies = []; max_area_per_node = pi * Rs^2;
    
    all_hbn = [];
    for i = 1:length(x_loc), if findHoleBoundaryNode_Yao(i, x_loc, y_loc, Rc), all_hbn(end+1) = i; end, end
    if isempty(all_hbn), avg_efficiency = 0; return; end
    
    boundary_coords = [x_loc(all_hbn), y_loc(all_hbn)];
    if size(boundary_coords, 1) < 3, avg_efficiency = 0; return; end
    
    cluster_indices = cluster_nodes_by_distance(boundary_coords, Rc);
    num_clusters = max(cluster_indices);
    if isempty(num_clusters) || num_clusters < 1, avg_efficiency = 0; return; end
    
    holes = cell(1, num_clusters);
    for i = 1:num_clusters, holes{i} = all_hbn(cluster_indices == i); end
    
    x_loc_current_iter = x_loc;
    y_loc_current_iter = y_loc;

    for i=1:length(holes)
        stats_before_this_hole = calculateKCoverage(x_loc_current_iter, y_loc_current_iter, Rs, net_length, net_width);
        [~, ~, ~, new_node_coords] = patchHole_Yao_for_efficiency(holes{i}, x_loc_current_iter, y_loc_current_iter, Rc, Rs, net_length, net_width);
        if ~isempty(new_node_coords)
            for n_idx = 1:size(new_node_coords, 1)
                new_pos = new_node_coords(n_idx, :);
                if isCoveredByAny(new_pos, x_loc_current_iter, y_loc_current_iter, Rs), continue; end

                stats_after_one_node = calculateKCoverage([x_loc_current_iter; new_pos(1)], [y_loc_current_iter; new_pos(2)], Rs, net_length, net_width);
                area_restored = max(0, stats_before_this_hole.hole_area - stats_after_one_node.hole_area);
                efficiency = (area_restored / max_area_per_node) * 100;
                all_efficiencies(end+1) = efficiency;

                x_loc_current_iter(end+1) = new_pos(1);
                y_loc_current_iter(end+1) = new_pos(2);
            end
        end
    end
    if isempty(all_efficiencies), avg_efficiency = 0; else, avg_efficiency = mean(all_efficiencies, 'omitnan'); end
end

function avg_efficiency = run_nlchr_algorithm_for_efficiency(x_loc_init, y_loc_init, Rs, area_dims)
    net_length = area_dims(1); net_width = area_dims(2);
    x_loc = x_loc_init(:); y_loc = y_loc_init(:);
    all_efficiencies = []; max_area_per_node = pi * Rs^2;

    potential_holes = find_triangular_holes(x_loc, y_loc, Rs);
    if isempty(potential_holes), avg_efficiency = 0; return; end
    
    all_new_node_positions = [];
    for i = 1:length(potential_holes)
        triplet_nodes = potential_holes{i};
        p1=[x_loc(triplet_nodes(1)),y_loc(triplet_nodes(1))]; p2=[x_loc(triplet_nodes(2)),y_loc(triplet_nodes(2))]; p3=[x_loc(triplet_nodes(3)),y_loc(triplet_nodes(3))];
        dists=[norm(p1-p2),norm(p2-p3),norm(p3-p1)]; [min_dist,~]=min(dists);
        if min_dist<=2*Rs
            [~,min_idx]=min(dists);
            if min_idx==1,pair=[p1;p2];elseif min_idx==2,pair=[p2;p3];else,pair=[p3;p1];end
            new_vertex=calculate_isosceles_vertex(pair(1,:),pair(2,:),Rs,mean([p1;p2;p3],1));
            if~isempty(new_vertex),all_new_node_positions=[all_new_node_positions;new_vertex];end
        end
    end
    
    if isempty(all_new_node_positions), avg_efficiency = 0; return; end
    final_positions = resolve_nlchr_conflicts(all_new_node_positions, Rs);
    if isempty(final_positions), avg_efficiency = 0; return; end
    
    x_loc_current_iter = x_loc;
    y_loc_current_iter = y_loc;
    for n_idx = 1:size(final_positions, 1)
        pos = final_positions(n_idx, :);
        if isCoveredByAny(pos, x_loc_current_iter, y_loc_current_iter, Rs), continue; end
        
        stats_before_this_node = calculateKCoverage(x_loc_current_iter, y_loc_current_iter, Rs, net_length, net_width);
        stats_after_one_node = calculateKCoverage([x_loc_current_iter; pos(1)], [y_loc_current_iter; pos(2)], Rs, net_length, net_width);
        area_restored = max(0, stats_before_this_node.hole_area - stats_after_one_node.hole_area);
        all_efficiencies(end+1) = (area_restored / max_area_per_node) * 100;

        x_loc_current_iter(end+1) = pos(1);
        y_loc_current_iter(end+1) = pos(2);
    end
    if isempty(all_efficiencies), avg_efficiency = 0; else, avg_efficiency = mean(all_efficiencies, 'omitnan'); end
end

% =========================================================================
% --- CONSOLIDATED HELPER FUNCTIONS (ALL OF THEM) ---
% =========================================================================
function [x,y]=generate_network(a,n),x=a(1)*rand(n,1);y=a(2)*rand(n,1);end
function stats=calculateKCoverage(x,y,R,L,W),res=10;xg=0:res:L;yg=0:res:W;[X,Y]=meshgrid(xg,yg);cm=zeros(size(X));for i=1:length(x),cm=cm+(((X-x(i)).^2+(Y-y(i)).^2)<=R^2);end;pa=res^2;ta=L*W;stats.hole_area=sum(cm(:)==0)*pa;stats.k_covered_area=ta-stats.hole_area;stats.k_covered_percentage=(stats.k_covered_area/ta)*100;end
function boundary_nodes=findBoundaryNodes(x,y,R,L,W),n=length(x);Ni=cell(n,1);for i=1:n,for j=1:n,if i~=j&&norm([x(i)-x(j),y(i)-y(j)])<2*R,Ni{i}=[Ni{i},j];end,end,end;bn=[];for i=1:n,[Pi,~,~]=computeCP_and_CIP(i,x,y,Ni,R,L,W);if isempty(Ni{i})||~isempty(Pi),bn(end+1)=i;end,end,boundary_nodes=bn;end
function T=selectTripletsDelaunay(B,x,y,R,p),if p<60,th=1.8*R;elseif p<85,th=1.6*R;else,th=1.4*R;end;T={};if length(B)<3,return;end;dt=delaunayTriangulation(x(B),y(B));C=dt.ConnectivityList;V=[];for i=1:size(C,1),n=B(C(i,:));p1=[x(n(1)),y(n(1))];p2=[x(n(2)),y(n(2))];p3=[x(n(3)),y(n(3))];d12=norm(p1-p2);d23=norm(p2-p3);d31=norm(p3-p1);if d12>th&&d23>th&&d31>th,[pc,~]=solveApolloniusForTriplet(n,x,y,R);if~isempty(pc)&&~isCoveredByAny(pc,x,y,R),V(end+1,:)=[n,d12+d23+d31];end;end;end;if~isempty(V),S=sortrows(V,4,'descend');for i=1:size(S,1),T{end+1}=S(i,1:3);end;end,end
function [pc,pr]=solveApolloniusForTriplet(t,x,y,r),c(1)=struct('x',x(t(1)),'y',y(t(1)),'r',r);c(2)=struct('x',x(t(2)),'y',y(t(2)),'r',r);c(3)=struct('x',x(t(3)),'y',y(t(3)),'r',r);p1=struct('x',c(1).x,'y',c(1).y);c2r=struct('x',c(2).x,'y',c(2).y,'r',0);c3r=struct('x',c(3).x,'y',c(3).y,'r',0);sols=solve_pcc_by_inversion(p1,c2r,c3r);if isempty(sols),pc=[];pr=[];return;end;[~,idx]=max([sols.r]);sol=sols(idx);pc=[sol.x,sol.y];pr=sol.r+r;end
function sols=solve_pcc_by_inversion(p,c1,c2),ic=p;ir=1;c1i=invert_circle(c1,ic,ir);c2i=invert_circle(c2,ic,ir);tl=find_common_tangents(c1i,c2i);sols=struct('x',{},'y',{},'r',{});for i=1:length(tl),if isstruct(tl(i)),sols(end+1)=invert_line(tl(i),ic,ir);end,end,end
function ic=invert_circle(c,center,k),d2=(c.x-center.x)^2+(c.y-center.y)^2;s=k^2/(d2-c.r^2);ic=struct('x',center.x+s*(c.x-center.x),'y',center.y+s*(c.y-center.y),'r',abs(s)*c.r);end
function ic=invert_line(l,c,k),d=2*(l.c+l.a*c.x+l.b*c.y);ic=struct('x',c.x-k^2*l.a/d,'y',c.y-k^2*l.b/d,'r',abs(k^2/d));end
function t=find_common_tangents(c1,c2),t=cell(1,4);vx=c2.x-c1.x;vy=c2.y-c1.y;d2=vx^2+vy^2;dr=c2.r-c1.r;if d2>=dr^2,d=sqrt(d2);rt=sqrt(d2-dr^2);a=(vx*dr+vy*rt)/d2;b=(vy*dr-vx*rt)/d2;t{1}=struct('a',a,'b',b,'c',c1.r-(a*c1.x+b*c1.y));a=(vx*dr-vy*rt)/d2;b=(vy*dr+vx*rt)/d2;t{2}=struct('a',a,'b',b,'c',c1.r-(a*c1.x+b*c1.y));end;sr=c1.r+c2.r;if d2>=sr^2,d=sqrt(d2);rt=sqrt(d2-sr^2);a=(vx*sr+vy*rt)/d2;b=(vy*sr-vx*rt)/d2;t{3}=struct('a',a,'b',b,'c',-c1.r-(a*c1.x+b*c1.y));a=(vx*sr-vy*rt)/d2;b=(vy*sr+vx*rt)/d2;t{4}=struct('a',a,'b',b,'c',-c1.r-(a*c1.x+b*c1.y));end;t=[t{:}];end
function c=isCoveredByAny(p,x,y,R),c=any(sqrt(sum((p-[x(:),y(:)]).^2,2))<=R+1e-9);end
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
function hc=findHoles_Circumcenter(x,y,R),hc=[];if length(x)<3,return;end;try,dt=delaunayTriangulation(x(:),y(:));catch,return;end;[cc,rc]=circumcenter(dt);uc=[];for i=1:size(cc,1),if rc(i)>R&&~isCoveredByAny(cc(i,:),x,y,R),uc(end+1,:)=cc(i,:);end,end;if~isempty(uc),hc=uniquetol(uc,1e-6,'ByRows',true);end,end
function is_hbn=findHoleBoundaryNode_Yao(id,x,y,Rc),N_X=find(sqrt(sum(( [x,y] - [x(id),y(id)] ).^2,2))<Rc&(1:length(x))~=id);if length(N_X)<2,is_hbn=true;return;end;sub_ids=[id;N_X];try,k_b=boundary(x(sub_ids),y(sub_ids),1);catch,is_hbn=true;return;end;is_hbn=ismember(id,sub_ids(unique(k_b)));end
function [x_new,y_new,added,new_coords]=patchHole_Yao_for_efficiency(h_nodes,x,y,Rc,~,L,W),x_new=x;y_new=y;added=0;new_coords=[];if length(h_nodes)<2,return;end;try,b_path=boundary(x(h_nodes),y(h_nodes),1);catch,return;end;for i=1:(length(b_path)-1),pA=[x(h_nodes(b_path(i))),y(h_nodes(b_path(i)))];pB=[x(h_nodes(b_path(i+1))),y(h_nodes(b_path(i+1)))];mid=(pA+pB)/2;v_AB=pB-pA;pv_norm=[-v_AB(2),v_AB(1)]/norm(v_AB);d_A_mid=norm(pA-mid);if Rc^2<d_A_mid^2,continue;end;h=sqrt(Rc^2-d_A_mid^2);p1=mid+h*pv_norm;p2=mid-h*pv_norm;if norm(p1-mean([x,y]))>norm(p2-mean([x,y])),f_pos=p1;else,f_pos=p2;end;if all(f_pos>=0 & f_pos<=[L,W]),x_new(end+1)=f_pos(1);y_new(end+1)=f_pos(2);added=added+1;new_coords(end+1,:)=f_pos;end,end,end
function C=cluster_nodes_by_distance(coords,th),N=size(coords,1);adj=false(N);for i=1:N,for j=i+1:N,if norm(coords(i,:)-coords(j,:))<th,adj(i,j)=true;adj(j,i)=true;end,end,end;C=conncomp(graph(adj,'omitselfloops'))';end
function cr=estimateCoverageInWSN(c,r,L,W),N=10000;p=rand(N,2).*[2*r,2*r]+(c-r);ip=p(sum((p-c).^2,2)<=r^2,:);if isempty(ip),cr=0;return;end;cr=sum(all(ip>=0&ip<=[L,W],2))/size(ip,1);end
function pts=borderIntersections(x,y,r,L,W),pts=[];if(x-r<=0),dy=sqrt(max(0,r^2-x^2));pts=[pts;[0,y-dy];[0,y+dy]];end;if(x+r>=L),dy=sqrt(max(0,r^2-(L-x)^2));pts=[pts;[L,y-dy];[L,y+dy]];end;if(y-r<=0),dx=sqrt(max(0,r^2-y^2));pts=[pts;[x-dx,0];[x+dx,0]];end;if(y+r>=W),dx=sqrt(max(0,r^2-(W-y)^2));pts=[pts;[x-dx,W];[x+dx,W]];end;if~isempty(pts),pts=pts(all(pts>=0&pts<=[L,W],2),:);pts=uniquetol(pts,1e-9,'ByRows',true);end,end
function isB=isBorderSensor(x,y,L,W,R),isB=(x-R<=0)||(x+R>=L)||(y-R<=0)||(y+R>=W);end
function c=isCoveredByAnyNeighbor(p,i,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
function c=isCoveredByAnotherSpecificNeighbor(p,i,j,Ni_i,x,y,R),c=false;for k=Ni_i,if k~=i&&k~=j&&norm(p-[x(k),y(k)])<=R,c=true;return;end,end,end
function pts=circleIntersections(x1,y1,x2,y2,r),d=norm([x1-x2,y1-y2]);if d>=2*r||d==0,pts=[];return;end;a=d/2;h=sqrt(r^2-a^2);xm=x1+a*(x2-x1)/d;ym=y1+a*(y2-y1)/d;pts=[xm-h*(y2-y1)/d,ym+h*(x2-x1)/d;xm+h*(y2-y1)/d,ym-h*(x2-x1)/d];end
function v=calculate_isosceles_vertex(p1,p2,s,hc),d=norm(p1-p2);if d>2*s,v=[];return;end;mid=(p1+p2)/2;h=sqrt(max(0,s^2-(d/2)^2));pv=[p2(2)-p1(2),p1(1)-p2(1)];if norm(pv)==0,v=[];return;end;upv=pv/norm(pv);v1=mid+h*upv;v2=mid-h*upv;if norm(v1-hc)<norm(v2-hc),v=v1;else,v=v2;end,end
function final_pos=resolve_nlchr_conflicts(pos,R)
    if size(pos,1)>1
        dist_matrix=zeros(size(pos,1)); for i=1:size(pos,1),for j=i+1:size(pos,1),dist_matrix(i,j)=norm(pos(i,:)-pos(j,:));dist_matrix(j,i)=dist_matrix(i,j);end,end
        dist_matrix(logical(eye(size(dist_matrix))))=inf;
        keep=true(size(pos,1),1);
        for j=1:size(pos,1)
            if keep(j),close_pts=dist_matrix(j,:)<R;keep(close_pts)=false;keep(j)=true;end
        end
        final_pos=pos(keep,:);
    else,final_pos=pos;
    end
end
function s=isAboveLine(px,py,x1,y1,x2,y2),v=(x2-x1)*(py-y1)-(y2-y1)*(px-x1);if abs(v)<1e-9,s='on';elseif v>0,s='above';else,s='below';end,end
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