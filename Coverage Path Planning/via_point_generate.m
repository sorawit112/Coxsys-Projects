%% find all points in area for given start vertex

%%
%(5,2)

ver_list = area_graph.Nodes.Vertices{i};
num_ver = numel(ver_list);


ver_start = ver_list(j);

inner_ver = area_graph.Nodes.InnerVertices{i};
opp_ver_i = area_graph.Nodes.VerOpp{i};
alpha_i = area_graph.Nodes.Alpha{i};
c = cos(-alpha_i(j));
s = sin(-alpha_i(j));
R = [c -s; s c];
offset = inner_ver(:,j);
new_area = R*(inner_ver-offset);

min_y = min(new_area(2,:));
max_y = max(new_area(2,:));
max_x = max(new_area(1,:));

grid_width = width*(sqrt(2)+0.2);

num_y = ceil(max_y/grid_width)-floor(min_y/grid_width)+1;
y_vector = linspace(floor(min_y/grid_width)*grid_width,ceil(max_y/grid_width)*grid_width,num_y);
num_x = ceil(max_x/(sqrt(3)*grid_width/2))+1;
x_vector = linspace(0,ceil(max_x/grid_width)*grid_width,num_x);

[X,Y] = meshgrid(x_vector,y_vector);
Y(:,1:2:end)=Y(:,1:2:end)-grid_width/2;
[in,on] = inpolygon(X,Y,new_area(1,:),new_area(2,:));
occ_grid = or(in,on);
shifted_grid = [permute(X,[3 1 2]) ; permute(Y,[3 1 2])];
actual_grid = zeros(size(shifted_grid));
for k = 1:size(occ_grid,2)
    actual_grid(:,:,k) = R'*shifted_grid(:,:,k)+offset;
end
actual_grid = permute(actual_grid,[2 3 1]);

num_point = 0;
via_point_area = [];
dir = 1;
for idx = 1:size(occ_grid,2)
    via_temp = [];
    for k = 1:size(occ_grid,1)
        
        if occ_grid(k,idx)
            if dir>0
                via_temp = [via_temp permute(actual_grid(k,idx,:),[3 1 2])];
            else
                via_temp = [permute(actual_grid(k,idx,:),[3 1 2]) via_temp];
            end
            num_point = num_point+1;
        end
    end
    via_point_area = [via_point_area via_temp];
    dir = -dir;
end

alpha_i = area_graph.Nodes.Alpha{i};
vertices_i = area_graph.Nodes.Vertices{i};
corner_i = area_graph.Nodes.Corner{i};
theta_i = area_graph.Nodes.Theta{i};

edge_width = width*(sqrt(2)-0.15);

via_point_edge = [];


for k = circshift(1:num_ver,-j+1)
    corner_a = corner_i(:,k);
    if k~= num_ver
        idx = k+1;
    else
        idx = 1;
    end
    corner_b = corner_i(:,idx);
    dcorner = norm(corner_a-wall_graph.Nodes.Coordinates(vertices_i(k),:)')-width;
    if dcorner>0
        corner_width = width+0.1;
        via_extra = corner_a+(dcorner+0.1)*[cos(alpha_i(k)+pi);sin(alpha_i(k)+pi)];
    else
        via_extra = [];
    end
    via_edge_shifted = 0:edge_width:norm(corner_a-corner_b);
    c = cos(theta_i(k));
    s = sin(theta_i(k));
    R = [c -s; s c];
    p = [via_extra corner_a+[via_edge_shifted*cos(theta_i(k));via_edge_shifted*sin(theta_i(k))]];
    via_point_edge = [via_point_edge p];
end


theta = 0:0.01:(2*pi);

hold on;
plot(via_point_area(1,:),via_point_area(2,:),'go');
plot(via_point_edge(1,:),via_point_edge(2,:),'kx');

% for k = 1:size(via_point_edge,2)
%     plot(via_point_edge(1,k)+width*cos(theta),via_point_edge(2,k)+width*sin(theta),'k')
% end
% for k = 1:size(via_point_area,2)
%     plot(via_point_area(1,k)+width*cos(theta),via_point_area(2,k)+width*sin(theta),'b')
% end

%% animation
c_e = [0.5 0.5 1];
c_a = [0.5 1 0.5];

for k = 1:size(via_point_edge,2)
    fill(via_point_edge(1,k)+width*cos(theta), via_point_edge(2,k)+width*sin(theta), c_e,'LineStyle','none')
    pause(0.1)
end
for k = 1:size(via_point_area,2)
    fill(via_point_area(1,k)+width*cos(theta), via_point_area(2,k)+width*sin(theta), c_a,'LineStyle','none')
    pause(0.1)
end

grid on;
axis equal;


% size(via_point_edge,2)+size(via_point_area,2)
