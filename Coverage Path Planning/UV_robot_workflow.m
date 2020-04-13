
construct_graph;
visualize_graph;
isSelected = [2 3 5 6 7 8];
start_area = 7; % select start location
current_pose = area_graph.Nodes.Centroid(start_area,:)';
goal_atea = 1; % select goal location

findPathArea;

[path_nodes path_clean]
theta = 0:0.01:(2*pi);
for k = 1:size(path_nodes,1)
    area = path_nodes(k);
    if path_clean(k)
        [~,j] = min(sum((area_graph.Nodes.Corner{area}-current_pose).^2));
        i = area;
        via_point_generate;
        current_pose = via_point_area(:,end);
    else
        robot_pos = area_graph.Nodes.Centroid(area,:)';
        h_plot = fill(robot_pos(1)+width*cos(theta), robot_pos(2)+width*sin(theta), 'r','LineStyle','none');
        pause(1.5)
        if k ~= size(path_nodes,1)
            delete(h_plot)
        end
        current_pose = robot_pos;
    end
end
%
% i = 2; % automatically select area in sequence
% j = 3;  % index will automatically select based on nearest
