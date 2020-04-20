vertices = [0 11 7 1; 0 0 6 5]';
sub_poly = cutPolyMinor(vertices);
c = ['r','b'];
grid on;
hold on;
axis equal
for i = 1:2
    ver = sub_poly{i};
    plot(ver([1:end 1],1),ver([1:end 1],2),c(i))
end