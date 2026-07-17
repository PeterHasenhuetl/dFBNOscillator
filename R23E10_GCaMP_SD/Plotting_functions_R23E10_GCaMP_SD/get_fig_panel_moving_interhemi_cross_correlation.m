function get_fig_panel_moving_interhemi_cross_correlation(cross_corr_table, curr_fly_name, dend_idx, ...
    x_pos, y_pos, r_lim, img_lim, lim_1, x_axis_cond, hemi_annotation,curr_blue, curr_red,...
    source_data_details)
% Plotting the time-varying crosscorrelation of left and right hemispheres
% of an example recording (dFBN>GCaMP).
% Code written by Peter Hasenhuetl.

color_map_1 = get_default_two_color_map(10000, 10000, curr_blue, curr_red);

sz_1 = 8;
ht_1 = 1.5;
ht_2 = 0.5;
ht_3 = 0.75;
y_dist_1 = 0.3;
y_dist_2 = 0.3;

clim_1 = [-1, 1];
lag_lim = [-3, 3];
major_lag = 1.5;
major_time = 5000;
major_time_axis = 2;

sliding_window = 300;


curr_table = cross_corr_table(cross_corr_table.fly_name == curr_fly_name,:);
curr_trace_1 = (curr_table.ihc_stuff.dendrite_1(:,dend_idx));
curr_trace_2 = (curr_table.ihc_stuff.dendrite_2(:,dend_idx));



curr_trace1_plot = zscore(curr_trace_1);
curr_trace2_plot = zscore(curr_trace_2);

curr_trace_1 = [zscore(diff(curr_trace_1)); 0];
curr_trace_2 = [zscore(diff(curr_trace_2)); 0];

length_trace = length(curr_trace_1);

% If sliding window is even, adds 1 (for symmetric sliding window)
is_even = ((sliding_window/2) == round(sliding_window/2));
if is_even == 1
    sliding_window = sliding_window+1;
end

% Adds a padding window before and after the raw trace
size_padding_window = (sliding_window-1)/2;
padding_window = zeros(size_padding_window,1);
curr_trace_1 = [padding_window; curr_trace_1; padding_window];
curr_trace_2 = [padding_window; curr_trace_2; padding_window];
curr_trace1_plot = [padding_window; curr_trace1_plot; padding_window];
curr_trace2_plot = [padding_window; curr_trace2_plot; padding_window];

% Loops through individual points of trace to create time-varying percentile
mov_x_corr = zeros(length_trace,sliding_window);
window_center_idx = size_padding_window+1;

for loop_idx = 1:length_trace
    
    curr_t1 = curr_trace_1(window_center_idx-size_padding_window:window_center_idx+size_padding_window,1);
    curr_t2 = curr_trace_2(window_center_idx-size_padding_window:window_center_idx+size_padding_window,1);

    [a,b] = xcov(curr_t1,curr_t2,length(padding_window),'coef');
    mov_x_corr(loop_idx,:) = a';

    curr_corr = corrcoef(curr_t1,curr_t2);
    mov_corr(loop_idx,1) = curr_corr(2,1);
    window_center_idx = window_center_idx+1;
end

% Plots the time-varying cross-correlation
t_v_mov_x_corr = 1:size(mov_x_corr,1);
xlm_trace = [1, size(mov_x_corr,1)];
time_lag_xcorr = b./14.56;
min_max = [min(mov_x_corr,[],"all"), max(mov_x_corr,[],"all")];
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
imagesc(pnl_1, t_v_mov_x_corr, time_lag_xcorr,mov_x_corr')
hold on
plot(pnl_1, xlm_trace, [0, 0], 'k:','LineWidth', 1)
colormap(pnl_1, color_map_1)
clim(clim_1)
get_default_ax(pnl_1, xlm_trace(1), xlm_trace(2), [], [], ...
    lag_lim(1), lag_lim(2), lag_lim(1), lag_lim(2),...
    major_time, major_lag, "linear", "linear", 'Time (s)', 'Lag (s)',...
    'none', 'k', 'none', 'k', sz_1, ht_1)
ht_color_bar = ht_1;
get_default_colorbar(pnl_1, x_pos+sz_1+0.1, y_pos, 0.2, ht_color_bar, ...
    clim_1, min_max, 'r', 'eastoutside')


% Plots the time-varying correlation coefficient
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos+ht_1+y_dist_1, sz_1, ht_2]);
plot(pnl_1, t_v_mov_x_corr, mov_corr, 'k', 'LineWidth', 0.5)
hold on
plot(pnl_1, xlm_trace, [0, 0], 'k:', 'LineWidth', 1)
get_default_ax(pnl_1, xlm_trace(1), xlm_trace(2), [], [], ...
    r_lim(1), r_lim(2), r_lim(1), r_lim(2),...
    500, abs(r_lim(1)), "linear", "linear", [], 'r',...
    'none', 'k', 'none', 'k', sz_1, ht_2)


% Plots the imaging traces
lim_2 = lim_1 + round(90*14.56);
curr_trace1_plot = curr_trace1_plot(lim_1:lim_2,1);
curr_trace2_plot = curr_trace2_plot(lim_1:lim_2,1);
xlm_trace_inset = [1, size(curr_trace1_plot,1)];
plotting_color_left = [0, 0, 0];
plotting_color_right = [0.5, 0.5, 0.5];
y_pos_img_trace = y_pos+ht_1+ht_2+y_dist_1+y_dist_2;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos_img_trace, sz_1, ht_3]);
plot(pnl_1, 1:size(curr_trace1_plot,1), curr_trace1_plot, 'Color', plotting_color_left, 'LineWidth', 0.5)
hold on
plot(pnl_1, 1:size(curr_trace1_plot,1), curr_trace2_plot, 'Color', plotting_color_right, 'LineWidth', 0.5)
plot(pnl_1, [size(curr_trace1_plot,1)-14.56*5, size(curr_trace1_plot,1)],...
    [img_lim(1), img_lim(1)], 'Color', 'k', 'LineWidth', get_default_scale_bar_width)
get_default_ax(pnl_1, xlm_trace_inset(1), xlm_trace_inset(2), [], [], ...
    img_lim(1), img_lim(2), [], [],...
    [], [], "linear", "linear", [],[],...
    'none', 'none', 'none', 'none', sz_1, ht_3)


pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos_img_trace, 0, ht_3]);
plot(pnl_1, [xlm_trace_inset(1), xlm_trace_inset(1)], [img_lim(1), img_lim(1)+2], 'k', 'LineWidth', get_default_scale_bar_width)
get_default_ax(pnl_1, xlm_trace_inset(1), xlm_trace_inset(2), [], [], ...
    img_lim(1), img_lim(2), [], [],...
    [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_3)


% Annotation for x- and y-scale bars of imaging traace
get_default_annotation(x_pos+sz_1, y_pos_img_trace-0.12, ...
    '5 s', 'k', 'normal', "right")
get_default_annotation_rotated(x_pos+sz_1+0.2, y_pos_img_trace+0.2, ...
    '2 s.d.', 'k', 'normal', "left")

if hemi_annotation == "with_hemi_annotation"
    get_default_annotation(x_pos, y_pos_img_trace+ht_3+0.2, ...
        'Left hemisphere', plotting_color_left, 'normal', "left")
    get_default_annotation(x_pos+1.75, y_pos_img_trace+ht_3+0.2, ...
        'Right hemisphere', plotting_color_right, 'normal', "left")
end


% Plots the lines of highlighting the inset
plot_height = ht_1+(2*y_dist_1)+ht_2;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, plot_height]);
plot(pnl_1, [1, lim_1, lim_1], [plot_height, plot_height-y_dist_1, 0], 'k', 'LineWidth', 0.2)
hold on
plot(pnl_1, [lim_2, lim_2, xlm_trace(2)], [0, plot_height-y_dist_1, plot_height], 'k', 'LineWidth', 0.2)
if (x_axis_cond == "with_x_axis") == 0
    plot(pnl_1,[xlm_trace(2)-round(120*14.56), xlm_trace(2)], [0, 0],...
        'k', 'LineWidth', get_default_scale_bar_width)
end
get_default_ax(pnl_1, xlm_trace(1), xlm_trace(2), [], [], ...
    0, plot_height, [], [],...
    major_time, [], "linear", "linear", 'Time (s)', 'Lag (s)',...
    'none', 'none', 'none', 'none', sz_1, plot_height)
pnl_1.Color = 'none';


% Plots the x-axis underneath the x-corr heatmap
if x_axis_cond == "with_x_axis"

    t_v_mov_x_corr = ((1:size(mov_x_corr,1))./14.56)./60;
    t_v_mov_x_corr = t_v_mov_x_corr - t_v_mov_x_corr(1);
    xlm_axis = [0, t_v_mov_x_corr(end)];
    pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, 0]);
    get_default_ax(pnl_1, xlm_axis(1), xlm_axis(2), xlm_axis(1), xlm_axis(2), ...
        lag_lim(1), lag_lim(2), [], [],...
        major_time_axis, major_lag, "linear", "linear", 'Time (min)', [],...
        'k', 'none', 'k', 'none', sz_1, ht_1)

else
    get_default_annotation(x_pos+sz_1, y_pos-0.12, ...
    '2 min', 'k', 'normal', "right")
end




if isempty(source_data_details) == 0
    
   

    bilateral_recording_traces = table;
    bilateral_recording_traces.img_trace = [curr_trace1_plot, curr_trace2_plot];
    writetable(bilateral_recording_traces, [source_data_details.data_path, ...
        [char(curr_fly_name), source_data_details.file_name2]],...
        'Sheet', 'bilateral_recordings')


    time_steps = (t_v_mov_x_corr./14.56)./60;
    interhemi_mov_corr_trace = table;
    interhemi_mov_corr_trace.x_values = time_steps';
    interhemi_mov_corr_trace.mov_corr = mov_corr;
    writetable(interhemi_mov_corr_trace, [source_data_details.data_path, ...
        [char(curr_fly_name), source_data_details.file_name2]],...
        'Sheet', 'interhemi_mov_corr_trace')


    time_steps = (t_v_mov_x_corr./14.56)./60;
    interhemi_mov_xcorr = table;
    interhemi_mov_xcorr.x_values = [NaN; time_steps'];
    interhemi_mov_xcorr.mov_x_corr = [time_lag_xcorr; mov_x_corr];
    writetable(interhemi_mov_xcorr, [source_data_details.data_path, ...
        [char(curr_fly_name), source_data_details.file_name2]],...
        'Sheet', 'mov_xcorr')

end


end



