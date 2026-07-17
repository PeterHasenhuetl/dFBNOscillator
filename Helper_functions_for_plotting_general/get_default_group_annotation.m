function get_default_group_annotation(x_pos, y_pos, group_name, symbol_color, text_color, annotation_style)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

sz_1 = 0.25;
xlm_1 = [0.5, 1.5];

curr_ax = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos+0.03, sz_1, 0.2]);
if annotation_style == "dot"  
    scatter(curr_ax, 1, 1, 10, 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', symbol_color, 'MarkerFaceAlpha', 1);
    x_dist = 0.05;
elseif annotation_style == "line"
    plot(curr_ax, xlm_1, [1, 1], 'Color', symbol_color, 'LineWidth', 1.5);
    x_dist = 0.05;
elseif annotation_style == "square"
    scatter(curr_ax, 1, 1, 10, 'Marker', "square", 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', symbol_color, 'MarkerFaceAlpha', 1);
    x_dist = 0.05;
end
curr_ax.Color = 'none';
curr_ax.XAxis.Color = 'none';
curr_ax.YAxis.Color = 'none';
xlim(xlm_1);    
ylim([0.5, 1.5])
        

annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos+sz_1+x_dist, y_pos, 10, 0.2], ...
    'string', char(group_name), 'EdgeColor', 'none',...
    'FontSize', get_default_font_size, 'FontWeight', 'normal', 'Color', text_color,...
    'FontAngle', 'italic', 'Margin', 0,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
    

end