% Vertices (Selected Points)
x = [-8 -4 0 2 0 -5 5 10 3 4 6 11 17 13 21];
y = [5 0 9 5 -4 -10 -1 2 -8 -15 -6 -8 4 -15 -5];
p = [x;y];

% Connectivity (drag and draw)
s = [1 1 2 2 2 3 4 4 5 5 5 6 7 7 8 8 9 9 10 10 11 12 12 13];
t = [2 3 4 5 6  4 7 8 6 7 9 10 8 11 12 13 10 11 12 14 12 14 15 15];
isWall = [1 1 0  0 1 1 0 1 0 1 1 1 0 1 0 1 0 1 0 1 0 1 1 1]; % select when draw

name_Nodes = 1:size(p,2);
name_Edges = 1:size(s,2);
EdgeTable = table([s' t'],name_Edges',isWall','VariableNames',{'EndNodes','Number','isWall'});
NodeTable = table(name_Nodes',p','VariableNames',{'Number','Coordinates'});
wall_graph = graph(EdgeTable,NodeTable);

%%
% select consecutive edges only

edges = {[1 3 6 2],[3 4 10 7],[4 5 9],[7 13 8],[9 12 17 11],[13 14 21 15],[17 19 21 18],[19 20 22],[15 23 24 16]};
% Connectivity (drag and draw)
s = [1 2 2 3 4 5 6 6 7];
t = [2 3 4 5 6 7 7 9 8];


% automatically detect vertices of an area based on edge
num_area = numel(edges);

vertices = cell(1,num_area);
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
centroid = zeros(2,num_area);  % area--nodes
area = zeros(1,num_area);
rho_offset = cell(1,num_area);
theta = cell(1,num_area);
alpha = cell(1,num_area);
inner_vertices = cell(1,num_area);
corner_vertices = cell(1,num_area);
ver_opp = cell(1,num_area);

distance = zeros(1,num_edges);
EdgeTable = table([s' t'],distance','VariableNames',{'EndNodes','Distance'});
NodeTable = table(name_Nodes',centroid',area',vertices',inner_vertices',corner_vertices',edges',theta',alpha',ver_opp',rho_offset','VariableNames',{'Number','Centroid','Area','Vertices','InnerVertices','Corner','Edges','Theta','Alpha','VerOpp','RhoOffset'});
area_graph = graph(EdgeTable,NodeTable);

%% find offset

width = 0.75;

for i = area_graph.Nodes.Number'
    
    coordinates = wall_graph.Nodes.Coordinates(area_graph.Nodes.Vertices{i},:)';
    area_graph.Nodes.Area(i) = polyarea(coordinates(1,:),coordinates(2,:));
    num_ver = numel(area_graph.Nodes.Vertices{i});
    area_graph.Nodes.Centroid(i,:) = (sum(coordinates,2)/num_ver)';
    
    ver_shift = coordinates'-area_graph.Nodes.Centroid(i,:);
    theta_i = zeros(1,num_ver);
    alpha_i = zeros(1,num_ver);
    rho = zeros(1,num_ver);
    rho_offset_i = zeros(1,num_ver);
    rho_temp = zeros(1,num_ver);
    p = zeros(2,num_ver); % inner vertices
    p_temp = zeros(2,num_ver);
    corner_vertices_i = zeros(1,num_ver);
    
    for j = 1:num_ver
        point_a = ver_shift(j,:);
        if j~=num_ver
            point_b = ver_shift(j+1,:);
        else
            point_b = ver_shift(1,:);
        end
        theta_i(j) = mod(atan2(point_b(2)-point_a(2),point_b(1)-point_a(1)),2*pi);
        rho(j) = -cos(theta_i(j))*point_b(2)+sin(theta_i(j))*point_b(1);
        rho_offset_i(j) = rho(j)-width*0.9; %***************
        rho_temp(j) = rho(j)-width/(sqrt(2)+0.1); % outter
    end
    
    for j = 1:num_ver
        theta_a = theta_i(j);
        rho_a = rho_offset_i(j);
        rho_c = rho_temp(j);
        if j~=num_ver
            idx = j+1;
        else
            idx = 1;
        end
        theta_b = theta_i(idx);
        rho_b = rho_offset_i(idx);
        rho_d = rho_temp(idx);
        theta_b = theta_b+(theta_b<theta_a)*2*pi;
        alpha_i(idx) = (theta_a+theta_b)/2+pi/2;
        p(:,idx) = [sin(theta_a) -cos(theta_a); sin(theta_b) -cos(theta_b)]\[rho_a;rho_b];
        p_temp(:,idx) = [sin(theta_a) -cos(theta_a); sin(theta_b) -cos(theta_b)]\[rho_c;rho_d];
        p_x = p_temp(1,idx)+area_graph.Nodes.Centroid(i,1);
        p_y = p_temp(2,idx)+area_graph.Nodes.Centroid(i,2);

    end
    area_graph.Nodes.Theta{i} = theta_i;
    area_graph.Nodes.Alpha{i} = alpha_i;
    area_graph.Nodes.RhoOffset{i} = rho_offset_i;
    area_graph.Nodes.InnerVertices{i} = p+area_graph.Nodes.Centroid(i,:)';
    area_graph.Nodes.Corner{i} = p_temp+area_graph.Nodes.Centroid(i,:)';
end

for i = 1:area_graph.numedges
    cm_ab = area_graph.Nodes.Centroid(area_graph.Edges.EndNodes(i,:),:);
    area_graph.Edges.Distance(i) = norm(cm_ab(1,:)-cm_ab(2,:));
end

%% find opposite vertices


for i = area_graph.Nodes.Number'
    ver_list = area_graph.Nodes.Vertices{i};
    
    num_ver = numel(ver_list);
    ver_opp_i = zeros(1,num_ver);
    for j = 1:num_ver
        ver_start = ver_list(j);
        wall_sub = wall_graph.subgraph(area_graph.Nodes.Vertices{i});
        alpha_i = area_graph.Nodes.Alpha{i};
        c = cos(-alpha_i(j));
        s = sin(-alpha_i(j));
        R = [c -s; s c];
        offset = wall_sub.Nodes.Coordinates(j,:)';
        new_area = R*(wall_sub.Nodes.Coordinates'-offset);
        [~,idx] = max(new_area(1,:));
        ver_opp_i(j) = ver_list(idx);
    end
    area_graph.Nodes.VerOpp{i} = ver_opp_i;
end
