energy = 67;   % to kill [J/m^2]

width_inner = 0.3; 
power_inner = 0.76;

width_outer = 1.3; % 1.20;
power_outer = 0.53;  % 0.33; => %51.12

offset_ratio = 1.5; % relative to width

% width_inner = 0.3;
% power_inner = 0.30;
% 
% width_outer = 1.3;
% power_outer = 0.30;
% 
% offset_ratio = 1.5; % relative to width
%{
result;
hours/room = 0.6224
room = 40.9444
%}

% width_inner = 0.1;
% power_inner = 0.20;
% 
% width_outer = 1.9;
% power_outer = 0.20;
% 
% offset_ratio = 1.5; % relative to width

%{
result;
hours/room = 0.3287
room = 40.9444
%}


grid_width = (width_outer-width_inner);
edge_width = (width_outer-width_inner);
power_actual = min(power_inner,power_outer);

duration = energy/power_actual;%duration_test*(width_outer/width_test)^2;
minute = duration/60;

%%

construct_graph_modified;
visualize_graph;

%%
isSelected = [2 3 4 5  8 9];
start_area = 1; % select start location
current_pose = area_graph.Nodes.Centroid(start_area,:)';
goal_atea = 1; % select goal location

findPathArea;

[path_nodes path_clean]
%%
path_edge = zeros(1,numel(path_nodes)-1);
for k = 1:numel(path_nodes)-1
    path_edge(k) = intersect(area_graph.Nodes.Edges{path_nodes(k)},area_graph.Nodes.Edges{path_nodes(k+1)});
end
   
%%
theta = 0:0.01:(2*pi);
num_via_edge = 0;
num_via_area = 0;

for k = 1:size(path_nodes,1)
    area = path_nodes(k);
    if path_clean(k)
        [~,j] = min(sum((area_graph.Nodes.Corner{area}-current_pose).^2));
        i = area;
        via_point_generate;
        if ~isempty(via_point_area)
            current_pose = via_point_area(:,end);
        else
            current_pose = via_point_edge(:,end);
        end
    else
        robot_pos = area_graph.Nodes.Centroid(area,:)';
        h_plot = fill(robot_pos(1)+width_outer*cos(theta), robot_pos(2)+width_outer*sin(theta), 'r','LineStyle','none');
        pause(1.5)
        if k ~= size(path_nodes,1)
            delete(h_plot)
        end
        current_pose = robot_pos;
    end
end

%%

(num_via_edge+num_via_area)*minute/(60*numel(isSelected))

total_area = sum(area_graph.Nodes.Area(isSelected));
average_area = total_area/numel(isSelected)

%
% i = 2; % automatically select area in sequence
% j = 3;  % index will automatically select based on nearest
