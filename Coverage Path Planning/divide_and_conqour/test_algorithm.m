clear
current_pose = [2 -1];
end_pose = [10 3];
radius = 2;
vertices = [0 11 7 1; 0 0 6 5]';
% vertices = [1 5 6 6 5 1 0; 0 0 1 5 6 6 5]'; % cut exactly at vertex ** special case
% vertices = [1 5 6 5 1 0 0; 0 0 1 6 6 5 1]'; 
theta = pi/6;

% vertices = ([cos(theta) -sin(theta); sin(theta) cos(theta)]*[1 5 6 6 5 1 0 0; 0 0 1 5 6 6 5 1])'; % symmetric octagon
% plot(vertices([1:end 1],1),vertices([1:end 1],2))

% vertices = [3 5 6 6 3; 0 0 1 3 3]'; % cut exactly at vertex ** special case
% vertices = [3 6 6 5 3; 3 3 5 6 6]';
% vertices = [1 3 3  0 0;0 0 3 3 1]';

isPlot = true;
[via_points,radius] = viaPointGenerator(vertices,radius,current_pose,end_pose,isPlot);