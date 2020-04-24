function [polygons,region_cm,region_radius] = coverage_points(vertices_i,radius)
cm = centroid(vertices_i);
if any(sum((vertices_i-cm).^2,2)>(radius^2))
    sub_poly = cutPolyMinor(vertices_i);
    hold on;
    grid on;
    axis equal;
    s_1 = sub_poly{1};
    s_2 = sub_poly{2};
    [polygon_1,region_cm_1,region_radius_1] = coverage_points(sub_poly{1},radius);
    [polygon_2,region_cm_2,region_radius_2] = coverage_points(sub_poly{2},radius);
    cm_1 = centroid(sub_poly{1});
    cm_2 = centroid(sub_poly{2});
    if cm_1(1)<cm_2(1)
        polygons = cat(2,polygon_1,polygon_2);
        region_cm = [region_cm_1;region_cm_2];
        region_radius = [region_radius_1;region_radius_2];
    else
        polygons = cat(2,polygon_2,polygon_1);
        region_cm = [region_cm_2;region_cm_1];
        region_radius = [region_radius_2;region_radius_1];
    end
else
    actual_radii = sqrt((sum((vertices_i-cm).^2,2)));
    [~,idx] = max(actual_radii);
    actual_radius = max(actual_radii(idx));
    polygons = {vertices_i};
    region_cm = cm;
    region_radius = actual_radius;
end

end

