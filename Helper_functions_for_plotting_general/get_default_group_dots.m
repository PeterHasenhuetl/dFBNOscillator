function get_default_group_dots(dot_matrix, x_pos, y_pos, sz_1, text_array)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

dot_matrix = dot_matrix';
n_col1 = size(dot_matrix,2);
n_row1 = size(dot_matrix,1);


curr_y_pos = y_pos-(0.2*n_row1);

curr_y_dist = 0;
for loop_idx = 1:n_row1
    
    dot_data = (1:n_col1).*dot_matrix(loop_idx,:);
    dot_data(dot_data == 0) = [];
    curr_ax = axes('Units', 'Centimeters', 'Position', [x_pos, curr_y_pos+curr_y_dist, sz_1, 0.2]);
    scatter(curr_ax, dot_data, ones(length(dot_data),1), 10, 'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 1);
    
    curr_ax.Color = 'none';    
    curr_ax.XAxis.Color = 'none';
    curr_ax.YAxis.Color = 'none';
    xlim([0.5, n_col1+0.5]);
    ylim([0.5, 1.5])
    
    curr_text = char(text_array(loop_idx));
        
    annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos, curr_y_pos+curr_y_dist-0.05, 0, 0.2], ...
        'string', curr_text, 'EdgeColor', 'none',...
        'FontSize', get_default_font_size, 'FontWeight', 'normal', 'Color', 'k',...
        'FontAngle', 'italic', 'Margin', 0,...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
    
    curr_y_dist = curr_y_dist+0.2;

end


end