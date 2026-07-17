function get_fig_panel_sleep_depr_img_traces(mp_data_table, x_pos, y_pos, ...
    selected_ZT, sd_color, ctrl_color, curr_start_trace_C, curr_end_trace_C, curr_start_trace_SD, curr_end_trace_SD, source_data_details)
% Plotting dFBN GCaMP traces in a rested and a sleep-deprived fly.
% Code written by Peter Hasenhuetl.


idx_1_SD = 2;
idx_2_SD = 4;
idx_1_C = 48;
idx_2_C = 3;
sz_1 = 2.8;
ht_1 = 2.4;
ylim_1 = [-0.075, 1.1];
line_width_traces = 0.5;
y_scale_bar_length = 0.2;
x_dist_cartoon = 1.2;
x_pos = x_pos + x_dist_cartoon;

ZT_vec = 0:4:28;
curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "no" & mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time < ZT_vec(selected_ZT+1),:);
curr_trace_ctrl = curr_tbl.img_trace(idx_1_C,1).trace_img(:,idx_2_C);

curr_tbl(idx_1_C,:)
curr_trace_ctrl = curr_trace_ctrl(curr_start_trace_C:curr_end_trace_C,1);

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl, curr_trace_ctrl, 14.56, ...
    5, y_scale_bar_length, ctrl_color, line_width_traces, "without_y_scale_bar", ylim_1, ...
    [0, 1], 0.2, "5 s", 0, [y_scale_bar_length;"\DeltaF/F"], 0.25)

curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "yes" & mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time < ZT_vec(selected_ZT+1),:);

curr_trace_SD = curr_tbl.img_trace(idx_1_SD,1).trace_img(:,idx_2_SD);
curr_trace_SD = curr_trace_SD(curr_start_trace_SD:curr_end_trace_SD,1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.2, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl, curr_trace_SD, 14.56, ...
    5, y_scale_bar_length, sd_color, line_width_traces,  "with_y_scale_bar", ...
    ylim_1, [0, 1], 0.2, "5 s", 0, [y_scale_bar_length;"\DeltaF/F"], 0.25)

%% Saves the source data
if isempty(source_data_details) == 0
    example_imaging_trace_SD = array2table(curr_trace_SD);
    writetable(example_imaging_trace_SD, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'example_imaging_trace_SD')
    
    example_imaging_trace_ctrl = array2table(curr_trace_ctrl);
    writetable(example_imaging_trace_ctrl, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'example_imaging_trace_ctrl')
end


end