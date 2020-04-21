clear
radius = 3;

graph_construction;
visualize_graph;
drawnow;

isSelected = [2 3 4 5  8 9];
start_pose = [-0.9 -0.4];
goal_pose = [-0.9 -0.4]; % select goal location

% detect area for start pose and goal pose
i = 1;
while i <= area_graph.numnodes
    ver_list = area_graph.Nodes.ActualVertices{i};
    coor = wall_graph.Nodes.Coordinates(ver_list,:);
    if inpolygon(start_pose(1),start_pose(2),coor(:,1),coor(:,2))
        start_area = i;
        break;
    end
    i = i +1;
end

i = 1;
while i <= area_graph.numnodes
    ver_list = area_graph.Nodes.ActualVertices{i};
    coor = wall_graph.Nodes.Coordinates(ver_list,:);
    if inpolygon(goal_pose(1),goal_pose(2),coor(:,1),coor(:,2))
        goal_area = i;
        break;
    end
    i = i +1;
end



% [path_nodes path_clean]

tic;
global_path_planner;
toc

num_node_path = numel(path_nodes);


current_pose = start_pose;

isPlot = false;
time = 0;
for i = 1:num_node_path
    if path_clean(i)
        path_idx = path_nodes(i);
        if i~=num_node_path
            end_pose = wall_graph.Edges.MidEdge(path_edge(i),:);
        else
            end_pose = area_graph.Nodes.Centroid(path_nodes(i),:);
        end

        [via_points,region_radius] = viaPointGenerator(wall_graph.Nodes.Coordinates(area_graph.Nodes.ActualVertices{path_idx},:),radius,current_pose,end_pose,isPlot);
        time = time + calculateTimeKill(region_radius);
        current_pose = via_points(end,:);
        for j = 1:size(via_points,1)
            h = plot(via_points(j,1),via_points(j,2),'mo','linewidth',3);
            pause(0.5);
        end
    else
        if i ~= num_node_path
            loc = wall_graph.Edges.MidEdge(path_edge(i),:);
            h = plot(loc(1),loc(2),'go','linewidth',3);
            pause(0.5);
            delete(h)
        else
            loc = area_graph.Nodes.Centroid(path_nodes(i),:);
            h = plot(loc(1),loc(2),'go','linewidth',3);
        end
    end
end



hours = floor(time/3600);
minutes = floor(mod(time,3600)/60);
seconds = floor(time-3600*hours-60*minutes);
area = sum(area_graph.Nodes.Area(path_nodes).*path_clean);

performance = time/area; % second/m^2

fprintf('Time Spent : %d Hours %d Minutes %d Seconds\n',hours,minutes,seconds);
fprintf('Performance : %4.4f s/m^2\n',performance);