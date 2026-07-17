function color_map = get_default_color_map(length_1, ...
    percentage_1, RGB_1, cond_dark)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

% defines the colormap
map_1_white(:,1) = linspace(1,RGB_1(1),length_1)';    
map_1_white(:,2) = linspace(1,RGB_1(2),length_1)'; 
map_1_white(:,3) = linspace(1,RGB_1(3),length_1)';

map_1_dark(:,1) = linspace(RGB_1(1),0,length_1)';    
map_1_dark(:,2) = linspace(RGB_1(2),0,length_1)';
map_1_dark(:,3) = linspace(RGB_1(3),0,length_1)'; 

if cond_dark == "normal"
    map_1 = map_1_white;
elseif cond_dark == "dark"
    % constructs the colormap and concatenates it according to 'percentage'
    map_1 = [map_1_white; map_1_dark(2:end,:)];
    idx_2 = round(length(map_1)*percentage_1);
    map_1 = map_1(1:idx_2,:);
end

color_map =  map_1;

end