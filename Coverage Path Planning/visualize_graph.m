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
    num_ver = numel(area_graph.Nodes.Vertices{i});
    corner_i = area_graph.Nodes.Corner{i};
    alpha_i = area_graph.Nodes.Alpha{i};
    for j = 1:num_ver
        
        plot(corner_i(1,j),corner_i(2,j),'mo');
        L = 1;
        plot([corner_i(1,j) corner_i(1,j)+L*cos(alpha_i(j))],[corner_i(2,j) corner_i(2,j)+L*sin(alpha_i(j))],'b')
        
    end
end
axis equal
grid on;
