function get_default_scatter_and_whisker(curr_pnl, x_axis_coordinate, curr_feature, plotting_color)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

if sum(isnan(curr_feature)) > 0 % If any of the datapoints contains a "NaN", throws a warning
    warning('!!!DATA CONTAIN NaN!!!')
    pause

elseif sum(isinf(curr_feature)) > 0 % If any of the datapoints contains a "inf", throws a warning
    warning('!!!DATA CONTAIN Inf!!!')
    pause

end


% Plots mean as scatter dot
scatter(curr_pnl,x_axis_coordinate,...
        mean(curr_feature),get_default_dot_size,'MarkerEdgeColor','none',...
    'MarkerFaceColor',plotting_color,'MarkerFaceAlpha',1);   
hold on

% Plots the error bars as ± SEM
plot(curr_pnl,[x_axis_coordinate,x_axis_coordinate],...
        [mean(curr_feature)-(std(curr_feature)/sqrt(length(curr_feature))),...
        mean(curr_feature)+(std(curr_feature)/sqrt(length(curr_feature)))],...
        'Color',plotting_color,'LineWidth',get_default_error_bar_width)

scatter(curr_pnl,x_axis_coordinate,...
        mean(curr_feature)-(std(curr_feature)/sqrt(length(curr_feature))),get_default_whisker_width,'Marker','_',...
        'MarkerEdgeColor',plotting_color,...
    'MarkerFaceColor',plotting_color,'MarkerFaceAlpha',1,'LineWidth',get_default_error_bar_width);
    
scatter(curr_pnl,x_axis_coordinate,...
        mean(curr_feature)+(std(curr_feature)/sqrt(length(curr_feature))),get_default_whisker_width,'Marker','_',...
        'MarkerEdgeColor',plotting_color,...
    'MarkerFaceColor',plotting_color,'MarkerFaceAlpha',1,'LineWidth',get_default_error_bar_width);

end