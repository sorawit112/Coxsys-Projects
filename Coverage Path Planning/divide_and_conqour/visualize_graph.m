%% Visualize Map
hold on;
for i = wall_graph.Edges.Number'
    wall_i = wall_graph.Edges.EndNodes(i,:);
    points = wall_graph.Nodes.Coordinates(wall_i',:)';
    if wall_graph.Edges.isWall(i)
        c = 'k';
    else
        c = 'r';
    end
    plot(points(1,:),points(2,:),c);
end

for i = area_graph.Nodes.Number'
    text(area_graph.Nodes.Centroid(i,1),area_graph.Nodes.Centroid(i,2),sprintf('%d',area_graph.Nodes.Number(i)));
end

min_x = min(wall_graph.Nodes.Coordinates(:,1));
min_y = min(wall_graph.Nodes.Coordinates(:,2));
max_x = max(wall_graph.Nodes.Coordinates(:,1));
max_y = max(wall_graph.Nodes.Coordinates(:,2));
range_x = max_x-min_x;
range_y = max_y-min_y;

axis equal
axis([min_x-0.1*range_x max_x+0.1*range_x min_y-0.1*range_y max_y+0.1*range_y])
grid on;
