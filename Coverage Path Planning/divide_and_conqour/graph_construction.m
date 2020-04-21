
% Vertices (Selected Points)
% x = [0, 0, 4, 4.5, 5, 5, 7, 7, 5, 7, 8, 8];
% y = [0, 3, 3,   0, 3, 7, 7, 0, 0, 1, 1, 0];
% p = [x;y];
% 
% s =      [1, 1, 2, 3, 3, 4, 5, 5, 6,  7, 8,  8,  8, 10, 11];
% t =      [2, 4, 3, 4, 6, 5, 6, 9, 7, 10, 9, 10, 12, 11, 12];
% isWall = [1, 1, 1, 0, 1, 1, 0, 1, 1,  1, 1,  0,  1,  1,  1];
x = [-75 -75 -75+110 -75+110 -75+110 -75+110 -75+110  -75+110+45 -75+110+45 -75+110+45 -75+110+45 -75+110+50 -75+110+50 -75+110+50 -75+100+45+22 -75+100+45+22 -75+110+50+25 -75+110+50+25+45 -75+110+50+80 -75+110+192-40 -75+110+192-40 -75+110+192-40 -75+110+192 -75+110+192 -75+110+192 -75+110+192 -75+110+192+550 -75+110+192+550];
y = [10-40 10 10-40 10 10-40+100 10-40+100+50 10-40+100+50+562 10-40+100+50 10-40+100+50+45 10-40+100+50+45+50 10-40+100+50+562 10-40 10-40+60 10-40+100 10-40+100+50+45 10-40+100+50+45+50 10-40 10-40-45 10-40+60-50 10-40-45 10-40+60-50 10-40+100 10-40-45 10-40 10-40+100 10-40+150 10-40-45 10-40];
p  = [x;y]*4/90;

% Connectivity (drag and draw)

s =      [1 1 2 3 3  4 5 5  6 6 7  8 8  9  9  10 10 12 12 13 13 14 15 17 17 18 19 20 20 21 22 23 23 24 24 25 27];
t =      [2 3 4 4 12 5 6 14 7 8 11 9 26 10 15 11 16 13 17 14 19 22 16 18 19 20 21 21 23 22 25 24 27 25 28 26 28];
isWall = [1 1 1 0 1  1 1 0  1 0 1  1 1  0  1  1  1  0  1  1  1  1  1  1  0  1  1  0  1  1  0  0  1  1  1  1  1]; % select when draw


name_Nodes = 1:size(p,2);
name_Edges = 1:size(s,2);
mid_edge = zeros(1,numel(isWall));
EdgeTable = table([s' t'],name_Edges',isWall',mid_edge','VariableNames',{'EndNodes','Number','isWall','MidEdge'});
NodeTable = table(name_Nodes',p','VariableNames',{'Number','Coordinates'});
wall_graph = graph(EdgeTable,NodeTable);

point_a = wall_graph.Nodes.Coordinates(wall_graph.Edges.EndNodes(:,1),:);
point_b = wall_graph.Nodes.Coordinates(wall_graph.Edges.EndNodes(:,2),:);
wall_graph.Edges.MidEdge = (point_a+point_b)/2;
% Conncet area
% select consecutive edges only

% edges = {[1, 2, 4, 3],[4, 6, 7, 5],[7, 8, 11, 12, 10, 9],[12, 13, 15, 14]};
% s = [1 2 3];
% t = [2 3 4];

edges = {[2 4 3 1],[5 18 20 8 6 4],[8 22 31 36 13 10 7],[10 12 14 16 11 9],[15 23 17 14],[19 25 21 18],[24 26 28 27 25],[29 32 34 31 30 28],[33 37 35 32]};
% Connectivity (drag and draw)
s = [1 2 2 3 3 4 6 7 8];
t = [2 3 6 4 8 5 7 8 9];


% automatically detect vertices of an area based on edge
num_area = numel(edges);

vertices = cell(1,num_area);
actual_vertices = cell(1,num_area);
% belongToArea = cell(1,wall_graph.numedges); % tell which edge belong to which area
for i = 1:num_area
    wall_i = wall_graph.Edges.EndNodes(edges{i},:);
    vertices_i = zeros(1,numel(edges{i}));
    edge_i = edges{i};
    for j = 1:numel(edges{i})
        if j == 1
            vertices_i(j) = intersect(wall_i(numel(edges{i}),:),wall_i(j,:));
        else
            vertices_i(j) = intersect(wall_i(j-1,:),wall_i(j,:));
        end
        vertices{i} = vertices_i;
        %         belongToArea{edge_i(j)} = [belongToArea{edge_i(j)} i];
    end
end
name_Nodes = 1:num_area;

num_edges = size(s,2);
center = zeros(2,num_area);  % area--nodes
area = zeros(1,num_area);
rho_offset = cell(1,num_area);
theta = cell(1,num_area);
alpha = cell(1,num_area);
inner_vertices = cell(1,num_area);
corner_vertices = cell(1,num_area);
ver_opp = cell(1,num_area);

distance = zeros(1,num_edges);
EdgeTable = table([s' t'],distance','VariableNames',{'EndNodes','Distance'});
NodeTable = table(name_Nodes',center',area',vertices',actual_vertices',edges','VariableNames',{'Number','Centroid','Area','Vertices','ActualVertices','Edges'});
area_graph = graph(EdgeTable,NodeTable);

% find centroid and actual vertices

for i = area_graph.Nodes.Number'
    
    coordinates = wall_graph.Nodes.Coordinates(area_graph.Nodes.Vertices{i},:)';
    area_graph.Nodes.Area(i) = polyarea(coordinates(1,:),coordinates(2,:));
    num_ver = numel(area_graph.Nodes.Vertices{i});
    area_graph.Nodes.Centroid(i,:) = centroid(coordinates');
    
    ver_shift = coordinates'-area_graph.Nodes.Centroid(i,:);
    theta_i = zeros(1,num_ver);
       
    for j = 1:num_ver
        point_a = ver_shift(j,:);
        if j~=num_ver
            point_b = ver_shift(j+1,:);
        else
            point_b = ver_shift(1,:);
        end
        theta_i(j) = mod(atan2(point_b(2)-point_a(2),point_b(1)-point_a(1)),2*pi);
    end
    % must compute the above regard less of perpendicularity
    
    actual_theta = [];
    actual_vertices = [];
    ver_list = area_graph.Nodes.Vertices{i};
    for j = 1:num_ver
        if isempty(actual_theta)
            actual_theta(end+1) = theta_i(j);
            actual_vertices(end+1) = ver_list(j);
        else
            if abs(theta_i(j-1)-theta_i(j))>0.000001
                actual_theta(end+1) = theta_i(j);
                actual_vertices(end+1) = ver_list(j);
            end
        end
    end
    
    area_graph.Nodes.ActualVertices{i} = actual_vertices;
    
    actual_coordinates = wall_graph.Nodes.Coordinates(area_graph.Nodes.ActualVertices{i},:)';
    num_actual_ver = numel(area_graph.Nodes.ActualVertices{i});
    area_graph.Nodes.Centroid(i,:) = (sum(actual_coordinates,2)/num_actual_ver)';
    
end

%% construct mid_edge_graph for traversal

connectivity_mid_edge = [];

for i = area_graph.Nodes.Number'
   edge_i = area_graph.Nodes.Edges{i};
   isConnected = ~wall_graph.Edges.isWall(edge_i);
    if sum(isConnected)>1
       temp_list = edge_i(isConnected);
       for j = 1:numel(temp_list)
           for k = j+1:numel(temp_list)
               connectivity_mid_edge = [connectivity_mid_edge;[temp_list(j) temp_list(k)]];
           end
       end
   end
end

sorted_list = sort(connectivity_mid_edge,2,'ascend'); 
connectivity_mid_edge = sortrows(sorted_list);

point_a = wall_graph.Edges.MidEdge(connectivity_mid_edge(:,1),:);
point_b = wall_graph.Edges.MidEdge(connectivity_mid_edge(:,2),:);
weight = sum((point_a-point_b).^2,2);

node_mid_edge = unique(connectivity_mid_edge);
[~,mid_edge_idx] = ismember(connectivity_mid_edge,node_mid_edge);

EdgeTable = table(mid_edge_idx,weight,'VariableNames',{'EndNodes','Weight'});
NodeTable = table(node_mid_edge,'VariableNames',{'EdgeNumber'});
% 
mid_edge_graph = graph(EdgeTable,NodeTable);
