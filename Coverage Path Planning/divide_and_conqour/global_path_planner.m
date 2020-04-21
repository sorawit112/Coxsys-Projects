%% FIND PATH OF AREA TO CLEAN

%% find all possible paths

k = 1;
current_graph = start_area;
candidate = cell(1,0);
max_iter = 5;
break_loop = false;
while(~break_loop)
    while (k<=max_iter)
        temp_graph = [];
        temp_edges = [];
        for i = 1:size(current_graph,2)
            neighbor_area = area_graph.neighbors(current_graph(end,i));
            temp_graph = [temp_graph [repmat(current_graph(:,i),[1,numel(neighbor_area)]);neighbor_area']];
        end
        current_graph = temp_graph;
        SG = current_graph(:,current_graph(end,:)==goal_area);
        for i = 1:size(SG,2)
            if all(ismember(isSelected,SG(:,i)))
                candidate{1,end+1} = SG(:,i);
            end
        end
        k = k+1;
    end
    
    if isempty(candidate)
        max_iter = max_iter + 1;
    else
        break_loop = true;
    end
    
end

% find minimum path
dist_path = zeros(1,numel(candidate));
for i = 1:numel(candidate)
    cand_i = candidate{i};
    path_edge = zeros(1,numel(cand_i)-1);
    for j = 1:numel(cand_i)-1
        path_edge(j) = intersect(area_graph.Nodes.Edges{cand_i(j)},area_graph.Nodes.Edges{cand_i(j+1)});
    end
    [~,idx]=ismember(path_edge,mid_edge_graph.Nodes.EdgeNumber);
    mid_edge_path = mid_edge_graph.Nodes.EdgeNumber(idx);
    idx_a =  mid_edge_path(1:end-1);
    idx_b =  mid_edge_path(2:end);
    dist_path(i) = sum(sum((wall_graph.Edges.MidEdge(idx_a,:)-wall_graph.Edges.MidEdge(idx_b,:)).^2));
end

[min_dist,idx] = min(dist_path);
path_nodes = candidate{idx};
cleanList = isSelected;
path_clean = zeros(numel(path_nodes),1);
for i = 1:numel(path_nodes)
    [isClean,idx] = ismember(path_nodes(i),cleanList);
    if isClean
        path_clean(i) = 1;
        cleanList(idx) = [];
    end
end

% area to travel to & status whether to clean
path_edge = zeros(1,numel(path_nodes)-1);
for i = 1:numel(path_nodes)-1
    path_edge(i) = intersect(area_graph.Nodes.Edges{path_nodes(i)},area_graph.Nodes.Edges{path_nodes(i+1)});
end

[~,idx]=ismember(path_edge,mid_edge_graph.Nodes.EdgeNumber);