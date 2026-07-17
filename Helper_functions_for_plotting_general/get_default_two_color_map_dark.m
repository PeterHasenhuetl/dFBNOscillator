function two_color_map = get_default_two_color_map_dark(length_1, length_2, ...
    percentage_1, percentage_2, RGB_1, RGB_2)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

% defines the first colormap
map_1_dark(:,1) = linspace(0,RGB_1(1),length_1)';    
map_1_dark(:,2) = linspace(0,RGB_1(2),length_1)';
map_1_dark(:,3) = linspace(0,RGB_1(3),length_1)'; 

map_1_white(:,1) = linspace(RGB_1(1),1,length_1)';    
map_1_white(:,2) = linspace(RGB_1(2),1,length_1)'; 
map_1_white(:,3) = linspace(RGB_1(3),1,length_1)';

% constructs the first colormap and concatenates it according to 'percentage 1'
map_1 = [map_1_dark; map_1_white(2:end,:)];
idx_1 = length(map_1)-(round(length(map_1)*percentage_1)-1);
map_1 = map_1(idx_1:end,:);

% defines the second colormap
map_2_white(:,1) = linspace(1,RGB_2(1),length_2)';    
map_2_white(:,2) = linspace(1,RGB_2(2),length_2)'; 
map_2_white(:,3) = linspace(1,RGB_2(3),length_2)';

map_2_dark(:,1) = linspace(RGB_2(1),0,length_2)';    
map_2_dark(:,2) = linspace(RGB_2(2),0,length_2)';
map_2_dark(:,3) = linspace(RGB_2(3),0,length_2)'; 

% constructs the second colormap and concatenates it according to 'percentage 2'
map_2 = [map_2_white; map_2_dark(2:end,:)];
idx_2 = round(length(map_2)*percentage_2);
map_2 = map_2(1:idx_2,:);

% concatenates both colormaps
two_color_map = [map_1; map_2(2:end,:)];

end