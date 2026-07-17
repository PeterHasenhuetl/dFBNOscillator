function get_default_colorbar(curr_ax, x_pos, y_pos, sz_1, ht_1, lim_1, min_max, label_1, orientation_1)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

label_rotation = 0;
ac_1 = colorbar(curr_ax,string(orientation_1));
ac_1.Units = 'centimeters';    
ac_1.Position = [x_pos, y_pos, sz_1, ht_1];
ylabel(ac_1,string(label_1))
ac_1.FontSize = get_default_font_size;
ac_1.Box = 'on';
ac_1.Color = 'k';
ac_1.TickLength = 0;
   
% Tests if colormap cuts off some values (e.g. if data are very skewed, or if there are outliers). 
% If so, it will indicate by "≤ " and "≥ ".
if min_max(1) < lim_1(1)
    curr_string1 = string(['\leq ', num2str(lim_1(1))]);
else
    curr_string1 = string(lim_1(1));
end
    
if min_max(2) > lim_1(2)    
    curr_string2 = char(['\geq ', num2str(lim_1(2))]);
else    
    curr_string2 = string(lim_1(2));
end

% First, tests if colormap includes negative AND positive values.
if lim_1(1) < 0 && lim_1(2) > 0
    ac_1.Ticks = [lim_1(1), 0, lim_1(2)];
    ac_1.TickLabels = [curr_string1, 0, curr_string2];  
else
    ac_1.Ticks = [lim_1(1), lim_1(2)];
    ac_1.TickLabels = [curr_string1, curr_string2];  
end
ac_1.Limits = lim_1;
ac_1.FontAngle = 'normal';
ac_1.Ruler.TickLabelRotation = label_rotation;

end






