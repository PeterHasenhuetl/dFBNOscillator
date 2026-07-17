function  get_default_img_trace_plot(curr_pnl1, curr_trace, frame_rate, ...
    x_sb_length, y_sb_length, plotting_color_1, line_width_plot, cond_y_axis, y_limit, y_tick_start_end, ytick_space,...
    xlabel, x_label_offset, y_label, y_label_offset)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

x_pos = curr_pnl1.Position(1);
y_pos = curr_pnl1.Position(2);
sz_1 = curr_pnl1.Position(3);
ht_1 = curr_pnl1.Position(4);
x_scale_bar = [length(curr_trace)-(frame_rate*x_sb_length), length(curr_trace)];
x_limit = [1, length(curr_trace)];

% Plots the trace
plot(curr_pnl1, curr_trace, 'Color', plotting_color_1, 'LineWidth', line_width_plot, 'LineStyle', '-')


% Plots y axis, y scale bar, or nothing, based on "cond_y_axis"
if cond_y_axis == "with_y_axis"
    
    get_default_ax(curr_pnl1, x_limit(1), x_limit(2), x_limit(1), x_limit(2), y_limit(1), y_limit(2),...
        y_tick_start_end(1), y_tick_start_end(2), [], ytick_space, "linear", "linear", [], char(y_label),...
        'none','k','none','k', sz_1, ht_1)

elseif cond_y_axis == "with_y_scale_bar"
    
    get_default_ax(curr_pnl1, x_limit(1), x_limit(2), x_limit(1), x_limit(2), y_limit(1), y_limit(2),...
        y_tick_start_end(1), y_tick_start_end(2), [], [], "linear", "linear", [], [],...
        'none', 'none', 'none', 'none', sz_1, ht_1)
    
    % Plots y axis scale bar and annotation/label
    curr_pnl_y_sb = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos, 0.1, ht_1]);
    plot(curr_pnl_y_sb, [1, 1], [y_limit(1), y_limit(1)+y_sb_length], 'k', 'LineWidth', get_default_scale_bar_width)
    get_default_ax(curr_pnl_y_sb, 1, 2, 1, 2, y_limit(1), y_limit(2),...
        [], [], [], [], "linear", "linear", [], [],...
        'none', 'none', 'none', 'none', 0.1, ht_1)
    get_default_annotation(x_pos+sz_1+0.2, y_pos+y_label_offset, ...
        char(y_label), 'k', 'normal', "left")
    
elseif cond_y_axis == "without_y_scale_bar"
    
    get_default_ax(curr_pnl1, x_limit(1), x_limit(2), x_limit(1), x_limit(2), y_limit(1), y_limit(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1)
    
end

% Plots the x scale bar and annotation/label
curr_pnl_x_sb = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-0.1, sz_1, 0.1]);
if x_sb_length > 0
    plot(curr_pnl_x_sb, x_scale_bar, [1, 1], 'k', 'LineWidth', get_default_scale_bar_width)
    get_default_annotation(x_pos+sz_1+x_label_offset, y_pos-0.2, ...
    char(xlabel), 'k', 'normal', "right")
end
get_default_ax(curr_pnl_x_sb, x_limit(1), x_limit(2), x_limit(1), x_limit(2), 1, 2,...
        1, 2, [], [], "linear", "linear", [], [],...
        'none', 'none', 'none', 'none', sz_1, 0.1)

end


