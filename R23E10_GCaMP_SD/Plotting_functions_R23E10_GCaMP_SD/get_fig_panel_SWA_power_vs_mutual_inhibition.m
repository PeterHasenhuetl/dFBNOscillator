function contra_mean = get_fig_panel_SWA_power_vs_mutual_inhibition(cross_corr_table, ...
    x_pos, y_pos, color, mean_cond, cond_n_flies, source_data_details)
% Plotting mean dFBN SWA power vs. interhemispheric inhibition.
% Code written by Peter Hasenhuetl.

sz_1 = 1.25;
ht_1 = 1.8;
ht_ipsi = 1.25;
ht_contra = 1.25;
xlm_scatter = [0, 0.6];
ylm_scatter = [-0.1, 0.02];
plotting_color_ipsi = color.yellow;
plotting_color_contra = color.yellow;
scatter_dot_color = color.dark_gray;
scatter_dot_size = 2;
scatter_dot_transparency = 0.5;
line_color = [0, 0, 0];
ylm_ipsi = [-0.05, 0.3];
ylm_contra = [-0.048, 0.015];
y_dist_traces = 0.4;

major_x_ticks_scatter = 0.2;
major_y_ticks_scatter = 0.05;


% empty arrays for loop below
GCaMP_ipsi = [];
GCaMP_contra = [];
delta_ipsi = [];
delta_contra = [];
ipsi_over_flies = [];
contra_over_flies = [];



% loops through flies (rows of cross_corr_table) and generates summary data
% for plotting
for loop_idx = 1:size(cross_corr_table,1)
    
    % gets GCaMP transients (hence called "ipsi" for ipsilateral)
    GCaMP_ipsi = [GCaMP_ipsi, ...
        cross_corr_table.ihc_stuff(loop_idx,1).ipsi_trace_left, cross_corr_table.ihc_stuff(loop_idx,1).ipsi_trace_right];
    curr_ipsi_left = mean(unique(cross_corr_table.ihc_stuff(loop_idx,1).ipsi_trace_left','rows')',2);
    curr_ipsi_right = mean(unique(cross_corr_table.ihc_stuff(loop_idx,1).ipsi_trace_right','rows')',2);
    ipsi_over_flies = [ipsi_over_flies,  mean([curr_ipsi_left, curr_ipsi_right],2)];

    % gets GCaMP traces on other hemisphere (hence called "contra" for contralateral)
    GCaMP_contra = [GCaMP_contra, ...
        cross_corr_table.ihc_stuff(loop_idx,1).contra_trace_left, cross_corr_table.ihc_stuff(loop_idx,1).contra_trace_right];
    curr_contra_left = mean(unique(cross_corr_table.ihc_stuff(loop_idx,1).contra_trace_left','rows')',2);
    curr_contra_right = mean(unique(cross_corr_table.ihc_stuff(loop_idx,1).contra_trace_right','rows')',2);
    contra_over_flies = [contra_over_flies,  mean([curr_contra_left, curr_contra_right],2)];
    
    % delta power of hemisphere where transients occured ("ipsi")
    delta_ipsi = [delta_ipsi, ...
        cross_corr_table.ihc_stuff(loop_idx,1).delta_power_left, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_right];
    % delta power of hemisphere opposite of where transients occured ("contra")
    delta_contra = [delta_contra, ...
        cross_corr_table.ihc_stuff(loop_idx,1).delta_power_right, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_left];
    
end

wind_img = 8:20;
SWA_power = mean([delta_ipsi; delta_contra],1);
contra_mean = mean(GCaMP_contra(wind_img,:));

% Fits linear regression to data
lin_mdl = fitlm(SWA_power,contra_mean);

% Plots the scatter plot of mean integrated delta power vs. mean
% contralateral dF/F
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(pnl_1, xlm_scatter, [0, 0], ':k', 'LineWidth', 0.5)
contra_mean_capped = contra_mean;
contra_mean_capped(contra_mean < ylm_scatter(1)) = ylm_scatter(1);
contra_mean_capped(contra_mean < ylm_scatter(1)) = [];
sz_capped = size(contra_mean,2) - size(contra_mean_capped,2);
contra_mean_triangles = ones(sz_capped).*ylm_scatter(1);
SWA_power_capped = SWA_power;
SWA_power_capped(contra_mean < ylm_scatter(1)) = [];
triangle_x_values = SWA_power(contra_mean < ylm_scatter(1));
get_default_axis_limits_warning(SWA_power, xlm_scatter)
get_default_axis_limits_warning(contra_mean, ylm_scatter)
warning('!!!Capped some points to be within y-axis limits!!!')
scatter(pnl_1, SWA_power_capped, contra_mean_capped, scatter_dot_size, 'MarkerFaceColor', scatter_dot_color,...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', scatter_dot_transparency);
scatter(pnl_1, triangle_x_values, contra_mean_triangles, scatter_dot_size, 'Marker', '^', 'MarkerFaceColor', scatter_dot_color,...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', scatter_dot_transparency);
plot(pnl_1, SWA_power, lin_mdl.Fitted, 'Color', line_color, 'LineWidth', 0.5)

get_default_separated_ax(pnl_1, xlm_scatter(1), xlm_scatter(2), xlm_scatter(1), xlm_scatter(2), ylm_scatter(1), ylm_scatter(2),...
    ylm_scatter(1), ylm_scatter(2), major_x_ticks_scatter, major_y_ticks_scatter, "linear", "linear", 'Mean SWA power',...
    {'Mean'; 'contralateral \DeltaF/F'}, 'k', 'k', 'k', 'k', sz_1, ht_1)

% Computes correlation coefficient of data on x and y axes in scatter plot
r = corrcoef(SWA_power,contra_mean);
norm_dist1 = ht_1+0.15;
get_default_annotation(x_pos+sz_1, y_pos+norm_dist1, ['r = ', num2str(round(r(2),2))], 'k', 'normal', "right")

quartile_4 = prctile(mean([delta_ipsi;delta_contra],1),75);
quartile_1 = prctile(mean([delta_ipsi;delta_contra],1),25);
median_1 = prctile(mean([delta_ipsi;delta_contra],1),50);

idx_q1 = find(mean([delta_ipsi;delta_contra],1) <= quartile_1);
idx_q2 = find(mean([delta_ipsi;delta_contra],1) > quartile_1 & mean([delta_ipsi;delta_contra],1) <= median_1);
idx_q3 = find(mean([delta_ipsi;delta_contra],1) > median_1 & mean([delta_ipsi;delta_contra],1) <= quartile_4);
idx_q4 = find(mean([delta_ipsi;delta_contra],1) > quartile_4); 
   
get_fig_panel_transient_aligned_contra(GCaMP_ipsi, GCaMP_contra, ipsi_over_flies, contra_over_flies, ...
    x_pos, y_pos+ht_1+y_dist_traces, sz_1, ht_ipsi, ht_contra, plotting_color_ipsi, plotting_color_contra,...
    ylm_ipsi, ylm_contra, idx_q1, idx_q2, idx_q3, idx_q4, mean_cond, source_data_details)

if cond_n_flies == 1
    n_flies = loop_idx;
    curr_annotation_y_pos = y_pos+ht_1+y_dist_traces+ht_ipsi+ht_contra+0.4;
    get_default_annotation(x_pos, curr_annotation_y_pos, [num2str(n_flies), ' flies'], 'k', 'normal', "left")
end


%% Saves the source data
if isempty(source_data_details) == 0
    power_vs_contra_dF = table;
    power_vs_contra_dF.x_values_mean_SWA_power = SWA_power';
    power_vs_contra_dF.contra_mean_DF = contra_mean';
    power_vs_contra_dF.linear_regression = lin_mdl.Fitted;
    writetable(power_vs_contra_dF, [source_data_details.data_path, source_data_details.file_name],'Sheet','power_vs_contra_dF')
end

end


function get_fig_panel_transient_aligned_contra(ipsi_trace, contra_trace, ipsi_over_flies, contra_over_flies, ...
    curr_x_pos, curr_y_pos, sz_1, ht_ipsi, ht_contra, plotting_color_ipsi, plotting_color_contra,...
    ylm_ipsi, ylm_contra, idx_q1, idx_q2, idx_q3, idx_q4, mean_cond, source_data_details)


y_dist = 0.3;
xlm_1 = [1, length(mean(ipsi_trace,2))];
t_v = xlm_1(1):xlm_1(2);

scale_bar = [xlm_1(2)-(0.5*14.56), xlm_1(2)];

quartile_color1 = [0.3, 0.3, 0.3];
quartile_color2 = quartile_color1+0.1;
quartile_color3 = quartile_color2+0.1;
quartile_color4 = quartile_color3+0.1;

if mean_cond == "over ROIs"
    mean_ipsi_trace = unique(ipsi_trace','rows')';
    mean_contra_trace = unique(contra_trace','rows')';
elseif mean_cond == "over flies"
    mean_ipsi_trace = ipsi_over_flies;
    mean_contra_trace = contra_over_flies;
    size(mean_contra_trace)
end


pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos+ht_contra+y_dist, sz_1, ht_ipsi]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, ipsi_trace(:,idx_q1), t_v, quartile_color1)
get_default_SEM_area_plot(pnl_1, ipsi_trace(:,idx_q2), t_v, quartile_color2)
get_default_SEM_area_plot(pnl_1, ipsi_trace(:,idx_q3), t_v, quartile_color3)
get_default_SEM_area_plot(pnl_1, ipsi_trace(:,idx_q4), t_v, quartile_color4)
get_default_SEM_area_plot(pnl_1, mean_ipsi_trace, t_v, plotting_color_ipsi)
y_scale_bar = [0, 0.1]+0.1;
scale_bar_x_dist = 3;
plot(pnl_1, [xlm_1(1)+scale_bar_x_dist, xlm_1(1)+scale_bar_x_dist], y_scale_bar,...
    'LineWidth', get_default_scale_bar_width, 'Color', 'k')
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_ipsi(1), ylm_ipsi(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_ipsi)    
get_default_annotation(curr_x_pos-0.55, curr_y_pos+ht_contra+y_dist+(ht_ipsi*0.7), ...
    {'0.1'; '\DeltaF/F'}, 'k', 'normal', "left")

pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos, curr_y_pos, sz_1, ht_contra]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, contra_trace(:,idx_q1), t_v, quartile_color1)
get_default_SEM_area_plot(pnl_1, contra_trace(:,idx_q2), t_v, quartile_color2)
get_default_SEM_area_plot(pnl_1, contra_trace(:,idx_q3), t_v, quartile_color3)
get_default_SEM_area_plot(pnl_1, contra_trace(:,idx_q4), t_v, quartile_color4)
get_default_SEM_area_plot(pnl_1, mean_contra_trace, t_v, plotting_color_contra)
y_scale_bar = [-0.02, 0]-0.01;
scale_bar_x_dist = 3;
plot(pnl_1,[xlm_1(1)+scale_bar_x_dist, xlm_1(1)+scale_bar_x_dist],y_scale_bar,...
    'LineWidth',get_default_scale_bar_width,'Color','k')
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_contra(1), ylm_contra(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_contra)
plot(pnl_1, scale_bar, [ylm_contra(1), ylm_contra(1)], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')

% some annotations
get_default_annotation(curr_x_pos-0.55, curr_y_pos+(ht_ipsi*0.55), ...
    {'0.02'; '\DeltaF/F'}, 'k', 'normal', "left")
get_default_annotation(curr_x_pos+sz_1, curr_y_pos+0.2, ...
    '0.5 s', 'k', 'normal', "right")
annot_dist = 0.1;
get_default_annotation(curr_x_pos+sz_1, curr_y_pos+ht_contra+y_dist+ht_ipsi+annot_dist, ...
    'Ipsi', 'k', 'normal', "right")
get_default_annotation(curr_x_pos+sz_1, curr_y_pos+ht_contra+annot_dist, ...
    'Contra', 'k', 'normal', "right")

%% Saves the source data
if isempty(source_data_details) == 0
    traces_ipsi = table;
    traces_ipsi.mean_traces = mean_ipsi_trace;
    traces_ipsi.blank_space1 = NaN(size(mean_ipsi_trace,1),1);
    traces_ipsi.traces_quart1 = ipsi_trace(:,idx_q1);
    traces_ipsi.blank_space2 = NaN(size(mean_ipsi_trace,1),1);
    traces_ipsi.traces_quart2 = ipsi_trace(:,idx_q2);
    traces_ipsi.blank_space3 = NaN(size(mean_ipsi_trace,1),1);
    traces_ipsi.traces_quart3 = ipsi_trace(:,idx_q3);
    traces_ipsi.blank_space4 = NaN(size(mean_ipsi_trace,1),1);
    traces_ipsi.traces_quart4 = ipsi_trace(:,idx_q4);
    
    
    traces_contra = table;
    traces_contra.mean_traces = mean_contra_trace;
    traces_contra.blank_space1 = NaN(size(mean_contra_trace,1),1);
    traces_contra.traces_quart1 = contra_trace(:,idx_q1);
    traces_contra.blank_space2 = NaN(size(mean_contra_trace,1),1);
    traces_contra.traces_quart2 = contra_trace(:,idx_q2);
    traces_contra.blank_space3 = NaN(size(mean_contra_trace,1),1);
    traces_contra.traces_quart3 = contra_trace(:,idx_q3);
    traces_contra.blank_space4 = NaN(size(mean_contra_trace,1),1);
    traces_contra.traces_quart4 = contra_trace(:,idx_q4);
    
    writetable(traces_ipsi, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_ipsi')
    writetable(traces_contra, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_contra')
end

end
