function get_default_scatter_n_numbers_annotation(x_pos, y_pos, scatter_x_size, n_numbers, group_text, text_color, font_angle)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

y_dist = 0.4;
norm_fac = scatter_x_size/length(n_numbers);
group_dist = 0.5:(length(n_numbers)-0.5);
norm_group_dist = group_dist*norm_fac;
num_groups = length(n_numbers);
font_size = 5;

for loop_idx = 1:num_groups
    
    if isnan(n_numbers(loop_idx)) == 0
        annotation('textbox', 'Units','centimeters','Position',...
            [x_pos+norm_group_dist(loop_idx), y_pos-y_dist, 0, 0.25], ...
            'string', string(n_numbers(loop_idx)), 'EdgeColor','none',...
            'FontSize',font_size,'FontWeight','normal','Color',text_color,...
            'Margin',0,'HorizontalAlignment','center')
    end

end

% Condition if it writes "n flies", for example, or simply "n"
if isempty(group_text) == 0
    
    annotation('textbox', 'Units','centimeters','Position',[x_pos-10, y_pos-y_dist, 10, 0.25], 'string', ...
        ['{\it n} ', char(group_text)], 'EdgeColor','none','Margin',0,...
        'FontSize',font_size,'FontWeight','normal','Color',[0,0,0],'FontAngle',font_angle,'HorizontalAlignment','right')

elseif isempty(group_text) == 1
    
    annotation('textbox', 'Units','centimeters','Position',[x_pos-10, y_pos-y_dist, 10, 0.25], 'string', ...
        '{\it n} ', 'EdgeColor','none','Margin',0,...
        'FontSize',font_size,'FontWeight','normal','Color',[0,0,0],'FontAngle',font_angle,'HorizontalAlignment','right')
    
end


end