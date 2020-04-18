function sub_poly = cutPolyMinor(vertices)
vertices = vertices';
centroid = sum(vertices,2)/size(vertices,2);

dp = (vertices-centroid);
[xxyy] = sum(dp.^2,2);

xy = sum(dp(1,:).*dp(2,:),2);
I = [xxyy(1) xy;xy xxyy(2)];

[V,~] = eig(I);

theta_cut = -atan2(V(1,1),V(2,1));
rho_cut = cos(theta_cut)*centroid(1)+sin(theta_cut)*centroid(2);


edge = [1 2 3 4 ; 2 3 4 1];

EdgeTable = table(edge',(1:size(edge,2))','VariableNames',{'EndNodes','Label'});
NodeTable = table((1:size(vertices,2))',vertices','VariableName',{'Label','Coordinates'});
sub_g = graph(EdgeTable,NodeTable);


num_ver = sub_g.numnodes;

start_ver = 1;
next_idx = start_ver;
sub_poly = cell(1,2);
for i = 1:2
    current_ver = next_idx;
    
    loop_break = false;
    coordinates = sub_g.Nodes.Coordinates;
    ver_list_minus = [coordinates(current_ver,:)];
    while ~loop_break
        if ~(current_ver == 1)
            idx = current_ver-1;
        else
            idx = num_ver;
        end
        current = coordinates(current_ver,:);
        look_ahead = coordinates(idx,:);
        theta_look_ahead = -atan2(look_ahead(1)-current(1),look_ahead(2)-current(2));
        rho_look_ahead = cos(theta_look_ahead)*current(1)+sin(theta_look_ahead)*current(2);
        
        A = [cos(theta_look_ahead) sin(theta_look_ahead);
            cos(theta_cut) sin(theta_cut)];
        b = [rho_look_ahead;rho_cut];
        cut_intersect = A\b;
        cut_direction = -atan2(cut_intersect(1)-current(1),cut_intersect(2)-current(2));
        if abs(cut_direction-theta_look_ahead)>0.001
            ver_list_minus = [ver_list_minus; look_ahead];
            current_ver = idx;
        else
            distance_look_ahead = sum((look_ahead-current).^2);
            distance_cut = sum((cut_intersect'-current).^2);
            if distance_cut <= distance_look_ahead
                ver_list_minus = [ver_list_minus; cut_intersect'];
                loop_break = true;
            else
                ver_list_minus = [ver_list_minus; look_ahead];
                current_ver = idx;
            end
            
        end
    end
    
    loop_break = false;
    current_ver = next_idx;
    coordinates = sub_g.Nodes.Coordinates;
    ver_list_plus = [coordinates(current_ver,:)];
    while ~loop_break
        if ~(current_ver == num_ver)
            idx = current_ver+1;
        else
            idx = 1;
        end
        current = coordinates(current_ver,:);
        look_ahead = coordinates(idx,:);
        theta_look_ahead = -atan2(look_ahead(1)-current(1),look_ahead(2)-current(2));
        rho_look_ahead = cos(theta_look_ahead)*current(1)+sin(theta_look_ahead)*current(2);
        
        A = [cos(theta_look_ahead) sin(theta_look_ahead);
            cos(theta_cut) sin(theta_cut)];
        b = [rho_look_ahead;rho_cut];
        cut_intersect = A\b;
        cut_direction = -atan2(cut_intersect(1)-current(1),cut_intersect(2)-current(2));
        if abs(cut_direction-theta_look_ahead)>0.001
            ver_list_plus = [ver_list_plus; look_ahead];
            current_ver = idx;
        else
            distance_look_ahead = sum((look_ahead-current).^2);
            distance_cut = sum((cut_intersect'-current).^2);
            if distance_cut <= distance_look_ahead
                ver_list_plus = [ver_list_plus; cut_intersect'];
                loop_break = true;
                next_idx = idx;
            else
                ver_list_plus = [ver_list_plus; look_ahead];
                current_ver = idx;
            end
            
        end
    end
    sub_poly{i} = [flipud(ver_list_minus(2:end,:));ver_list_plus];
end

end


