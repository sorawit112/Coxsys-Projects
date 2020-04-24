function cm = centroid(vertices_i)
vertices_i = vertices_i([1:end 1],:);
A = vertices_i(1:end-1,1).*vertices_i(2:end,2)-vertices_i(2:end,1).*vertices_i(1:end-1,2);
As = sum(A)/2;
x_bar = (sum((vertices_i(2:end,1)+vertices_i(1:end-1,1)).*A)*1/6)/As;
y_bar = (sum((vertices_i(2:end,2)+vertices_i(1:end-1,2)).*A)*1/6)/As;
cm = [x_bar y_bar];
end
