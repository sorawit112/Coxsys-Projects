function [via_points,path_radius] = viaPointGenerator(vertices,radius,current_pose,end_pose,varargin)
% for given vertices, find path inside
if ~isempty(varargin)
    isPlot = varargin{1};
end

[polygons,region_cm,region_radius] = coverage_points(vertices,radius);
d_current = region_cm-current_pose;
d_end = region_cm-end_pose;

[~,idx_start] = min(hypot(d_current(:,1),d_current(:,2)));
[~,idx_end] = min(hypot(d_end(:,1),d_end(:,2)));

del = delaunay(region_cm); % calculate triangular mesh
[edge_delaunay,~,~] = unique(sort([del(:,1:2);del(:,2:3);del(:,[3 1])],2,'ascend'),'rows');
weight = sqrt(sum((region_cm(edge_delaunay(:,1),:)-region_cm(edge_delaunay(:,2),:)).^2,2));

% remove edge that's too long
% new_idx = abs(weight-median(weight))>std(weight);
% edge_delaunay(new_idx,:) = [];
% weight(new_idx,:) = [];

EdgeTable = table(edge_delaunay,weight,'VariableNames',{'EndNodes','Weight'});
NodeTable = table((1:numel(polygons))',polygons',region_cm,region_radius,'VariableNames',{'Label','Vertices','Centroid','Radius'});

delaunay_graph = graph(EdgeTable,NodeTable); % graph of sub-area

%
nStops = delaunay_graph.numnodes;
idxs = delaunay_graph.Edges.EndNodes;
dist = delaunay_graph.Edges.Weight;
lendist = delaunay_graph.numedges;

% hGraph = plot(delaunay_graph,'XData',region_cm(:,1),'YData',region_cm(:,2),'LineStyle','none','NodeLabel',{});
% hold on

Aeq = spalloc(nStops,length(idxs),nStops*(nStops-1)); % Allocate a sparse matrix

for ii = 1:nStops
    if ii == idx_start
        whichIdxs = (idxs == ii);
        whichIdxs = sparse(sum(whichIdxs,2));
        Aeq(ii,:) = whichIdxs';
    elseif ii == idx_end
        whichIdxs = (idxs == ii);
        whichIdxs = sparse(sum(whichIdxs,2));
        Aeq(ii,:) = whichIdxs';
    else
        whichIdxs = (idxs == ii); % Find the trips that include stop ii
        whichIdxs = sparse(sum(whichIdxs,2)); % Include trips where ii is at either end
        Aeq(ii,:) = whichIdxs'; % Include in the constraint matrix
    end
end
beq = 2*ones(nStops,1);

if idx_start ~= idx_end
    beq(idx_start) = 1;
    beq(idx_end) = 1;
end

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = 2*ones(lendist,1);

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

number_x_tsp = round(x_tsp);
x_tsp = logical(number_x_tsp);
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));
% hold on
% highlight(hGraph,Gsol,'LineStyle','-')
% title('Solution with Subtours')

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
% fprintf('# of subtours: %d\n',numtours);

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,nStops)]; % A guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1) + 1; % Counter for indexing
        subTourIdx = find(tourIdxs == ii); % Extract the current subtour
        %         The next lines find all of the variables associated with the
        %         particular subtour, then add an inequality constraint to prohibit
        %         that subtour and all subtours that use those stops.
        s_graph = delaunay_graph.subgraph(subTourIdx);
        variations = s_graph.Edges.EndNodes;
        for jj = 1:size(variations,1)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
    end
    
    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
    
    x_tsp = logical(round(x_tsp));
    Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));
    % Visualize result
    %     hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    %     highlight(hGraph,Gsol,'LineStyle','-')
    %     drawnow
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    
end

edge_local_path = Gsol.Edges.EndNodes;
current_idx = idx_start;

path_idx = zeros(size(edge_local_path,1)+1,1);
path_idx(1) = current_idx;

for i = 1:size(edge_local_path,1)
    [idx_i,idx_j] = find(edge_local_path==current_idx);
    if numel(idx_i)>1
        idx_i = idx_i(1);
        idx_j = idx_j(1);
    end
    current_idx = edge_local_path(idx_i,-idx_j+3);
    edge_local_path(idx_i,:) = [];
    path_idx(i+1)=current_idx;
end


via_points = region_cm(path_idx,:);
path_radius = region_radius(path_idx);

%visualization
if isPlot
    hold on;
    grid on;
    axis equal
    for i = 1:numel(polygons)
        ver_i = polygons{i};
        cm = region_cm(i,:);
        %         plot(ver_i([1:end 1],1),ver_i([1:end 1],2),'k','linewidth',2)
        plot(ver_i([1:end 1],1),ver_i([1:end 1],2),'linewidth',2)
        %         pause(1)
        plot(cm(1),cm(2),'ro','linewidth',4)
        theta = 0:0.01:(2*pi);
        plot(region_radius(i)*cos(theta)+cm(1),region_radius(i)*sin(theta)+cm(2),'g');
        
    end
    for i = 1:size(edge_delaunay,1)
        point_a = region_cm(edge_delaunay(i,1),:);
        point_b = region_cm(edge_delaunay(i,2),:);
        
        plot([point_a(1) point_b(1)],[point_a(2) point_b(2)],'r--')
    end
    plot(vertices([1:end 1],1),vertices([1:end 1],2),'b','linewidth',3);
    
    for i = 1:Gsol.numedges
        point_a = region_cm(Gsol.Edges.EndNodes(i,1),:);
        point_b = region_cm(Gsol.Edges.EndNodes(i,2),:);
        
        plot([point_a(1) point_b(1)],[point_a(2) point_b(2)],'g','linewidth',4)
    end
    
end