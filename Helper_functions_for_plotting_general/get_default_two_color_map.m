function two_color_map = get_default_two_color_map(length_1, length_2, RGB_1, RGB_2)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

map_1(:,1) = linspace(RGB_1(1),1,length_1)';
map_1(:,2) = linspace(RGB_1(2),1,length_1)';
map_1(:,3) = linspace(RGB_1(3),1,length_1)';

map_2(:,1) = linspace(1,RGB_2(1),length_2)';
map_2(:,2) = linspace(1,RGB_2(2),length_2)';
map_2(:,3) = linspace(1,RGB_2(3),length_2)';

two_color_map = [map_1; map_2(2:end,:)];

end