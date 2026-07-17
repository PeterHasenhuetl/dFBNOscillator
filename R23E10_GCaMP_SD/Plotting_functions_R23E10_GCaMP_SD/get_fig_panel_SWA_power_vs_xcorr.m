function corrcoef_delta = get_fig_panel_SWA_power_vs_xcorr(cross_corr_table, x_pos, y_pos, color, mean_cond, source_data_details)
% Plotting mean dFBN SWA power vs. interhemispheric crosscorrelations.
% Code written by Peter Hasenhuetl.


sz_1 = 1.25;
sz_2 = 1.25;
ht_1 = 1.8;
ht_auto = 1;
ht_cross = 1.6;
y_dist_traces = 0.4;
xlm_scatter = [0, 0.6];
ylm_scatter = [-0.6, 0.15];
plotting_color_xcorr = color.yellow;
scatter_dot_color = color.dark_gray;
scatter_dot_size = get_default_scatter_dot_size;
scatter_dot_transparency = 0.5;
line_color = [0, 0, 0];
xlm_xcorr = [-5, 5];
ylm_xcorr = [-0.225, 0.15];
ylm_autocorr = [-0.3, 0.2];

major_x_ticks_scatter = 0.2;
major_y_ticks_scatter = 0.25;

% Empty arrays for loop below
autocorr_dendrite_left = [];
autocorr_dendrite_right = [];
cross_corr_delta = [];
cross_corr_delta_over_flies = [];
corrcoef_delta = [];
delta_left = [];
delta_right = [];
autocorr_both_hemi = [];
autocorr_both_hemi_over_flies = [];

% Loops through flies (rows of cross_corr_table) and generates summary data
% for plotting
for loop_idx = 1:size(cross_corr_table,1)

    autocorr_dendrite_left = [autocorr_dendrite_left, cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_left];
    autocorr_dendrite_right = [autocorr_dendrite_right, cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_right];
    curr_ac_both_hemi = (cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_left +...
        cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_right)/2;
    autocorr_both_hemi = [autocorr_both_hemi, curr_ac_both_hemi];
    autocorr_both_hemi_over_flies = [autocorr_both_hemi_over_flies, mean(curr_ac_both_hemi,2)];
    
    cross_corr_delta = [cross_corr_delta, cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta];
    cross_corr_delta_over_flies = [cross_corr_delta_over_flies, mean(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2)];
    corrcoef_delta = [corrcoef_delta, cross_corr_table.ihc_stuff(loop_idx,1).corrcoef_delta];
    delta_left = [delta_left, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_left];
    delta_right = [delta_right, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_right];
    
end

autocorr_both_hemi = (autocorr_dendrite_left + autocorr_dendrite_right)/2;

SWA_power = mean([delta_left; delta_right],1);
time_lag = cross_corr_table.ihc_stuff(1,1).time_lag;
lin_mdl = fitlm(SWA_power,corrcoef_delta);

quartile_4 = prctile(mean([delta_left;delta_right],1),75);
quartile_1 = prctile(mean([delta_left;delta_right],1),25);
median_1 = prctile(mean([delta_left;delta_right],1),50);

idx_q1 = find(mean([delta_left;delta_right],1) <= quartile_1);
idx_q2 = find(mean([delta_left;delta_right],1) > quartile_1 & mean([delta_left;delta_right],1) <= median_1);
idx_q3 = find(mean([delta_left;delta_right],1) > median_1 & mean([delta_left;delta_right],1) <= quartile_4);
idx_q4 = find(mean([delta_left;delta_right],1) > quartile_4); 
    
% Plots the mean integrated delta power vs. interhemi. corr. coef.
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(pnl_1, xlm_scatter, [0, 0], ':k', 'LineWidth', 0.5)
scatter(pnl_1, SWA_power, corrcoef_delta, scatter_dot_size, 'MarkerFaceColor', scatter_dot_color,...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', scatter_dot_transparency);
get_default_axis_limits_warning(SWA_power, xlm_scatter)
get_default_axis_limits_warning(corrcoef_delta, ylm_scatter)
plot(pnl_1, SWA_power, lin_mdl.Fitted, 'Color', line_color, 'LineWidth', 0.5)
r = corrcoef(SWA_power, corrcoef_delta);
get_default_separated_ax(pnl_1, xlm_scatter(1), xlm_scatter(2), xlm_scatter(1), xlm_scatter(2), ylm_scatter(1), ylm_scatter(2),...
    ylm_scatter(1), ylm_scatter(2), major_x_ticks_scatter, major_y_ticks_scatter, "linear", "linear",...
    'Mean SWA power', 'Correlation coefficient', 'k', 'k', 'k', 'k', sz_1, ht_1)

% Displays correlation coefficient of x vs. y axes in plot
norm_dist1 = ht_1+0.15;
get_default_annotation(x_pos+sz_1, y_pos+norm_dist1, ['r = ', num2str(round(r(2),2))], 'k', 'normal', "right")

% Plots the cross correlation
get_fig_panel_xcorr_plot(cross_corr_delta, cross_corr_delta_over_flies, ...
    autocorr_both_hemi, autocorr_both_hemi_over_flies, time_lag, ...
    x_pos, y_pos+ht_1+y_dist_traces, sz_2, ht_auto, ht_cross, plotting_color_xcorr,...
    xlm_xcorr, ylm_autocorr, ylm_xcorr, idx_q1, idx_q2, idx_q3, idx_q4, mean_cond, source_data_details)

%% Saves the source data
if isempty(source_data_details) == 0
power_vs_cross_corr = table;
power_vs_cross_corr.x_values_mean_SWA_power = SWA_power';
power_vs_cross_corr.interhemi_corrcoef = corrcoef_delta';
power_vs_cross_corr.linear_regression = lin_mdl.Fitted;
writetable(power_vs_cross_corr, [source_data_details.data_path, source_data_details.file_name],'Sheet','power_vs_cross_corr')
end

end


function get_fig_panel_xcorr_plot(crosscorr_trace, cross_corr_delta_over_flies, ...
    autocorr_both_hemi, autocorr_both_hemi_over_flies, time_lag, ...
    curr_x_pos, curr_y_pos, sz_1,ht_auto, ht_cross, plotting_color_contra,...
    xlm_xcorr, ylm_autocorr, ylm_xcorr, idx_q1, idx_q2, idx_q3, idx_q4, mean_cond, source_data_details)


y_dist = 0.1;

% Generates time vector (xcorr time lag divided by frame-rate))
t_v = time_lag/14.56;

% Plotting colors for sorting according to different quartiles
quartile_color1 = [0.3, 0.3, 0.3];
quartile_color2 = quartile_color1+0.1;
quartile_color3 = quartile_color2+0.1;
quartile_color4 = quartile_color3+0.1;

% Plots mean xcorr either as mean over flies or mean over ROIs
if mean_cond == "over ROIs"
    mean_crosscorr_trace = unique(crosscorr_trace','rows')';
    mean_autocorr_trace = unique(autocorr_both_hemi','rows')';
elseif mean_cond == "over flies"
    mean_crosscorr_trace = cross_corr_delta_over_flies;
    mean_autocorr_trace = autocorr_both_hemi_over_flies;
end

arrow_idx = [-2.129, 0];

marker_size = 5;
triangle_color = [0.3, 0.3, 0.3];
pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos+ht_cross+y_dist, sz_1, ht_auto+0.15]);
scatter(pnl_1, arrow_idx, ones(1, length(arrow_idx)), marker_size, 'Marker', 'v', 'MarkerFaceColor', triangle_color,...
    'MarkerEdgeColor', 'none')
get_default_ax(pnl_1, xlm_xcorr(1), xlm_xcorr(2), [], [], 0, 1,...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_auto)

pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos, sz_1, ht_cross]);
scatter(pnl_1, arrow_idx, [0, 0].*ones(1, length(arrow_idx)), marker_size, 'Marker', '^', 'MarkerFaceColor', triangle_color, ...
    'MarkerEdgeColor', 'none')
get_default_ax(pnl_1, xlm_xcorr(1), xlm_xcorr(2), xlm_xcorr(1), xlm_xcorr(2), 0, 1,...
    0, 1, [], 0.1, "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_cross)

% Plots the data
pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos+ht_cross+y_dist, sz_1, ht_auto]);
hold on
plot(pnl_1, xlm_xcorr, [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, autocorr_both_hemi(:,idx_q1), t_v, quartile_color1)
get_default_SEM_area_plot(pnl_1, autocorr_both_hemi(:,idx_q2), t_v, quartile_color2)
get_default_SEM_area_plot(pnl_1, autocorr_both_hemi(:,idx_q3), t_v, quartile_color3)
get_default_SEM_area_plot(pnl_1, autocorr_both_hemi(:,idx_q4), t_v, quartile_color4)
get_default_SEM_area_plot(pnl_1, mean_autocorr_trace, t_v, plotting_color_contra)
get_default_ax(pnl_1, xlm_xcorr(1), xlm_xcorr(2), [], [], ylm_autocorr(1), ylm_autocorr(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_cross)
y_scale_bar = [-0.1, 0]-0.1;
scale_bar_x_dist = 1;
plot(pnl_1, [xlm_xcorr(1)+scale_bar_x_dist, xlm_xcorr(1)+scale_bar_x_dist], y_scale_bar,...
    'LineWidth', get_default_scale_bar_width, 'Color', 'k')
get_default_annotation(curr_x_pos-0.4, curr_y_pos+ht_cross+y_dist+(ht_auto*0.3)+0.1, ...
    {'0.1 r'}, 'k', 'normal', "left")

% Plots the data
pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos, sz_1, ht_cross]);
hold on
plot(pnl_1, xlm_xcorr, [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, crosscorr_trace(:,idx_q1), t_v, quartile_color1)
get_default_SEM_area_plot(pnl_1, crosscorr_trace(:,idx_q2), t_v, quartile_color2)
get_default_SEM_area_plot(pnl_1, crosscorr_trace(:,idx_q3), t_v, quartile_color3)
get_default_SEM_area_plot(pnl_1, crosscorr_trace(:,idx_q4), t_v, quartile_color4)
get_default_SEM_area_plot(pnl_1, mean_crosscorr_trace, t_v, plotting_color_contra)
y_scale_bar = [-0.1, 0]-0.05;
scale_bar_x_dist = 1;
plot(pnl_1,[xlm_xcorr(1)+scale_bar_x_dist, xlm_xcorr(1)+scale_bar_x_dist],y_scale_bar,...
    'LineWidth', get_default_scale_bar_width, 'Color', 'k')
get_default_annotation(curr_x_pos-0.4, curr_y_pos+(ht_cross*0.4), ...
    {'0.1 r'}, 'k', 'normal', "left")
get_default_ax(pnl_1, xlm_xcorr(1), xlm_xcorr(2), [], [], ylm_xcorr(1), ylm_xcorr(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none','none','none','none', sz_1, ht_cross)

% Scale bar plus annotation
scale_bar = [xlm_xcorr(2)-2, xlm_xcorr(2)];
plot(pnl_1, scale_bar,...
    [ylm_xcorr(1), ylm_xcorr(1)], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
get_default_annotation(curr_x_pos+sz_1, curr_y_pos+0.2, ...
    '2 s', 'k', 'normal', "right")

%% Saves the source data

if isempty(source_data_details) == 0

    auto_corr_traces = table;
    auto_corr_traces.t_v = t_v.';
    auto_corr_traces.mean_traces = mean_autocorr_trace;
    auto_corr_traces.blank_space1 = NaN(size(mean_autocorr_trace,1),1);
    auto_corr_traces.traces_quart1 = autocorr_both_hemi(:,idx_q1);
    auto_corr_traces.blank_space2 = NaN(size(mean_autocorr_trace,1),1);
    auto_corr_traces.traces_quart2 = autocorr_both_hemi(:,idx_q2);
    auto_corr_traces.blank_space3 = NaN(size(mean_autocorr_trace,1),1);
    auto_corr_traces.traces_quart3 = autocorr_both_hemi(:,idx_q3);
    auto_corr_traces.blank_space4 = NaN(size(mean_autocorr_trace,1),1);
    auto_corr_traces.traces_quart4 = autocorr_both_hemi(:,idx_q4);
    
    xcorr_traces = table;
    xcorr_traces.t_v = t_v.';
    xcorr_traces.mean_traces = mean_crosscorr_trace;
    xcorr_traces.blank_space1 = NaN(size(mean_crosscorr_trace,1),1);
    xcorr_traces.traces_quart1 = crosscorr_trace(:,idx_q1);
    xcorr_traces.blank_space2 = NaN(size(mean_crosscorr_trace,1),1);
    xcorr_traces.traces_quart2 = crosscorr_trace(:,idx_q2);
    xcorr_traces.blank_space3 = NaN(size(mean_crosscorr_trace,1),1);
    xcorr_traces.traces_quart3 = crosscorr_trace(:,idx_q3);
    xcorr_traces.blank_space4 = NaN(size(mean_crosscorr_trace,1),1);
    xcorr_traces.traces_quart4 = crosscorr_trace(:,idx_q4);
    
    writetable(auto_corr_traces, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'auto_corr_traces')
    writetable(xcorr_traces, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'xcorr_traces')

end

end

