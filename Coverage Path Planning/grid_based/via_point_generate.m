%% find all points in area for given start vertex

%%
%(5,2)

ver_list = area_graph.Nodes.ActualVertices{i};
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

% rectanuglar
num_y = ceil(max_y/grid_width)-floor(min_y/grid_width)+1;
y_vector = linspace(floor(min_y/grid_width)*grid_width,ceil(max_y/grid_width)*grid_width,num_y);
num_x = ceil(max_x/grid_width)+1;
x_vector = linspace(0,ceil(max_x/grid_width)*grid_width,num_x);

[X,Y] = meshgrid(x_vector,y_vector);
[in,on] = inpolygon(X,Y,new_area(1,:),new_area(2,:));
occ_grid = or(in,on);
shifted_grid = [permute(X,[3 1 2]) ; permute(Y,[3 1 2])];
actual_grid = zeros(size(shifted_grid));
for k = 1:size(occ_grid,2)
    actual_grid(:,:,k) = R'*shifted_grid(:,:,k)+offset;
end
actual_grid = permute(actual_grid,[2 3 1]);

% %triangle
% num_y = ceil(max_y/grid_width)-floor(min_y/grid_width)+1;
% y_vector = linspace(floor(min_y/grid_width)*grid_width,ceil(max_y/grid_width)*grid_width,num_y);
% num_x = ceil(max_x/(sqrt(3)*grid_width/2))+1;
% x_vector = linspace(0,ceil(max_x/grid_width)*grid_width,num_x);
% 
% [X,Y] = meshgrid(x_vector,y_vector);
% Y(:,1:2:end)=Y(:,1:2:end)-grid_width/2;
% [in,on] = inpolygon(X,Y,new_area(1,:),new_area(2,:));
% occ_grid = or(in,on);
% shifted_grid = [permute(X,[3 1 2]) ; permute(Y,[3 1 2])];
% actual_grid = zeros(size(shifted_grid));
% for k = 1:size(occ_grid,2)
%     actual_grid(:,:,k) = R'*shifted_grid(:,:,k)+offset;
% end
% actual_grid = permute(actual_grid,[2 3 1]);

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
vertices_i = area_graph.Nodes.ActualVertices{i};
corner_i = area_graph.Nodes.Corner{i};
theta_i = area_graph.Nodes.Theta{i};


via_point_edge = [];


for k = circshift(1:num_ver,-j+1)
    corner_a = corner_i(:,k);
    if k~= num_ver
        idx = k+1;
    else
        idx = 1;
    end
    corner_b = corner_i(:,idx);
    dcorner = norm(corner_a-wall_graph.Nodes.Coordinates(vertices_i(k),:)')-width_outer;
    if dcorner>0
        corner_width = width_outer+0.1;
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
% if ~isempty(via_point_area)
%     plot(via_point_area(1,:),via_point_area(2,:),'go');
%     plot(via_point_edge(1,:),via_point_edge(2,:),'kx');
%     
% end
via_point_edge_temp = via_point_edge;

for k = 1:size(via_point_edge_temp)
    for kk = (k+1):size(via_point_edge_temp)
        if norm(via_point_edge_temp(:,k)-via_point_edge_temp(:,kk))<width_outer*2
            via_point_edge_temp(:,kk) = [];
        end
    end
end

via_point_edge = via_point_edge_temp
%% animation
c_e = [0.5 0.5 1];
c_a = [0.5 1 0.5];

dur = 0.01;
for k = 1:size(via_point_edge,2)
    outer_x = via_point_edge(1,k)+width_outer*cos(theta);
    inner_x = via_point_edge(1,k)+width_inner*cos(fliplr(theta));
    outer_y = via_point_edge(2,k)+width_outer*sin(theta);
    inner_y = via_point_edge(2,k)+width_inner*sin(fliplr(theta));
    
    fill([outer_x inner_x], [outer_y inner_y], c_e,'LineStyle','none')
    pause(dur)
end
for k = 1:size(via_point_area,2)
    outer_x = via_point_area(1,k)+width_outer*cos(theta);
    inner_x = via_point_area(1,k)+width_inner*cos(fliplr(theta));
    outer_y = via_point_area(2,k)+width_outer*sin(theta);
    inner_y = via_point_area(2,k)+width_inner*sin(fliplr(theta));
    
    fill([outer_x inner_x], [outer_y inner_y], c_a,'LineStyle','none')
    
    pause(dur)
end

grid on;
axis equal;

%<<<<<<< Updated upstream

% size(via_point_edge,2)+size(via_point_area,2)
%=======
num_via_edge = num_via_edge+size(via_point_edge,2);
num_via_area = num_via_area+size(via_point_area,2);
%>>>>>>> Stashed changes
