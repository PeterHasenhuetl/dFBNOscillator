function get_default_scatter_group(curr_ax, input_data, x_axis_coordinate, scatter_color)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

error_bar_width = 0.75;
mean_width = 1.5;
x_1 = -0.15;
x_2 = 0.15;
scatter_jitter = x_1 + (x_2-x_1)*(rand([length(input_data),1]));
scatter_alpha = 1;

scatter_dot_size = get_default_scatter_dot_size;
scatter(curr_ax,x_axis_coordinate*ones(length(input_data),1)+scatter_jitter,input_data,scatter_dot_size,'MarkerEdgeColor','none',...
    'MarkerFaceColor',scatter_color,'MarkerFaceAlpha',scatter_alpha);
hold on
% plotting the mean
plot(curr_ax,[-0.25,0.25]+x_axis_coordinate,[mean(input_data),mean(input_data)],'k','LineWidth',mean_width)

% plotting the SEM error bars
plot(curr_ax,[-0.125,0.125]+x_axis_coordinate,...
    [mean(input_data)+(std(input_data)/sqrt(length(input_data))),mean(input_data)+(std(input_data)/sqrt(length(input_data)))],...
    'k','LineWidth',error_bar_width)
plot(curr_ax,[-0.125,0.125]+x_axis_coordinate,...
    [mean(input_data)-(std(input_data)/sqrt(length(input_data))),mean(input_data)-(std(input_data)/sqrt(length(input_data)))],...
    'k','LineWidth',error_bar_width)
plot(curr_ax,[x_axis_coordinate,x_axis_coordinate],...
    [mean(input_data)-(std(input_data)/sqrt(length(input_data))),mean(input_data)+(std(input_data)/sqrt(length(input_data)))],...
    'k','LineWidth',error_bar_width)

end