function get_fig_panel_axon_vs_dendrites_traces(ad_corr_table, ...
    x_pos, y_pos, cv_idx, s_r, r_1, r_2, r_3, r_4, ylm_1, ylm_2, source_data_details)
% Plotting simultaneous GCaMP recordings from dFBN axons and dendrites.
% Code written by Peter Hasenhuetl.

ylm_3 =  ylm_2;
color = get_color;
sz_1 = 4;
sz_hist = 2.2;
ht_dend = 0.4;
ht_axon = 2.7;
ht_axon_inset = 2.3;
y_dist_axon = 0.5;
x_dist = 2;

color_model = color.navy;
axon_color = color.medium_gray;
dendrite_color = color.medium_gray;

col_1(:,1) = linspace(1,0,1000);
col_1(:,2) = linspace(1,0,1000);
col_1(:,3) = linspace(1,0,1000);
curr_colormap = [col_1(:,1), col_1(:,2), col_1(:,3)];

img_tbl = ad_corr_table;
selected_row = img_tbl(s_r,:);

img_trace_dendrites = selected_row.axo_dendritic_corr.dend_pred;
img_trace_dendrites = img_trace_dendrites(:,isnan(img_trace_dendrites(1,:,cv_idx)) == 0,cv_idx);

img_trace_axon = selected_row.axo_dendritic_corr.ax_to_pred;
img_trace_axon = img_trace_axon(:,cv_idx);

img_trace_prediction = selected_row.axo_dendritic_corr.ax_pred;
img_trace_prediction = img_trace_prediction(:,cv_idx);

for i = 1:size(img_trace_dendrites,2)
    curr_y_pos = y_pos+ht_axon_inset+y_dist_axon+((i-1)*ht_dend)+((i-1)*0.1);
    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, curr_y_pos, sz_1, ht_dend]);
    hold on
    if i == 1
        cond_scale_bar = "with_y_scale_bar";
        x_scale_bar_label = "20 s";
        x_scale_bar = 20;
    else
        cond_scale_bar = "without_y_scale_bar";
        x_scale_bar_label = [];
        x_scale_bar = [];
    end
    max(img_trace_dendrites(r_1:r_2,i))
    min(img_trace_dendrites(r_1:r_2,i))
    get_default_img_trace_plot(curr_pnl, img_trace_dendrites(r_1:r_2,i), 14.56, ...
        x_scale_bar, 0.5, dendrite_color,0.5,  cond_scale_bar, ...
        ylm_1, [0, 1], 0.5, x_scale_bar_label, 0, ["0.5"; "\DeltaF/F"], 0.25);

end

s_b = length(img_trace_axon(r_1:r_2,:));
scale_bar = [s_b-(14.56*20), s_b];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist, y_pos+ht_axon_inset+y_dist_axon, sz_1, ht_axon]);
plot(curr_pnl, img_trace_axon(r_1:r_2,:), 'Color', axon_color)
hold on
plot(curr_pnl, img_trace_prediction(r_1:r_2,:), 'Color', color_model)
plot(curr_pnl, [s_b, s_b], [ylm_2(1), ylm_2(1)+0.3], 'k', 'LineWidth', get_default_scale_bar_width);
plot(curr_pnl, [r_3, r_3], ylm_2, 'k:')
plot(curr_pnl, [r_4, r_4], ylm_2, 'k:')
curr_pnl.Color = 'none';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
ylim(ylm_2)
xlim([1, length((img_trace_axon(r_1:r_2,:)))])
plot(curr_pnl, scale_bar, [ylm_2(1), ylm_2(1)], 'k', 'LineWidth', get_default_scale_bar_width);
get_default_annotation(x_pos+(x_dist+sz_1)+sz_1, y_pos+ht_axon_inset+y_dist_axon-0.1, ...
    '20 s', [0, 0, 0], 'normal', "right")


s_b = length(img_trace_axon(r_3:r_4,:));
scale_bar = [s_b-(14.56*5), s_b];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist, y_pos, sz_1, ht_axon_inset]);
plot(curr_pnl, img_trace_axon(r_3:r_4,:), 'Color', axon_color)
hold on
plot(curr_pnl, img_trace_prediction(r_3:r_4,:), 'Color', color_model)
plot(curr_pnl, scale_bar, [ylm_3(1), ylm_3(1)], 'k', 'LineWidth', get_default_scale_bar_width);
plot(curr_pnl, [s_b, s_b], [ylm_3(1), ylm_3(1)+0.3], 'k', 'LineWidth', get_default_scale_bar_width);
curr_pnl.Color = 'none';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
ylim(ylm_3)
xlim([1, length((img_trace_axon(r_3:r_4,:)))])
get_default_annotation(x_pos+(x_dist+sz_1)+sz_1, y_pos-0.1, ...
    '5 s', [0, 0, 0], 'normal', "right")

get_default_annotation(x_pos+(x_dist+sz_1), y_pos-0.1, ...
    'Data', axon_color, 'normal', "left")
get_default_annotation(x_pos+(x_dist+sz_1)+0.75, y_pos-0.1, ...
    'model', color_model, 'normal', "left")
get_default_annotation(x_pos, y_pos+ht_axon_inset+y_dist_axon+ht_axon+0.3, ...
    'Dendrites', [0, 0, 0], 'normal', "left")
get_default_annotation(x_pos+(x_dist+sz_1), y_pos+ht_axon_inset+y_dist_axon+ht_axon+0.3, ...
    'Axons', [0, 0, 0], 'normal', "left")


curr_hist = get_data_vs_model_axon_dendrite(ad_corr_table);
curr_edges = ad_corr_table.axo_dendritic_corr(3,1).x_edges;
xy_vec = curr_edges(2:end)-(mean(diff(curr_edges))/2); 
xlm_1 = [-3, 3];
ylm_1 = [-3, 3];
x_pos_model_vs_data = x_pos+(sz_1/2)-(sz_hist/2);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos_model_vs_data, y_pos+0.1, sz_hist, sz_hist]);
hold on
m_scat_hist = curr_hist.*100;
imagesc(curr_pnl, xy_vec, xy_vec, m_scat_hist.') % Note the transpose
curr_pnl.Color = 'none';
colormap(curr_pnl, curr_colormap)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), 1.5, 1.5, "linear", "linear", 'Data (z-score)', 'Model (z-score)',...
    'k', 'k', 'k', 'k', sz_hist, sz_hist)

get_default_annotation(x_pos_model_vs_data+sz_hist, y_pos+0.2+0.1, [num2str(size(ad_corr_table,1)), ' flies'], 'k', 'normal', "right")


%% Adds colorbar

f_s = get_default_font_size;
clim_1 = [0, max(max(m_scat_hist))];
clrbr_1 = colorbar(curr_pnl,'eastoutside');
clrbr_1.Units = 'centimeters';
clrbr_1.Position = [x_pos_model_vs_data+sz_hist+0.1, y_pos+0.1,  0.2, sz_hist];
ylabel(clrbr_1,'Percent')
clrbr_1.FontSize = f_s;
clrbr_1.Color = 'k';
clrbr_1.TickLength = 0;
clrbr_1.Limits = clim_1;
clrbr_1.YLabel.Visible = 'on';
clrbr_1.YLabel.Color = 'k';


if isempty(source_data_details) == 0

    data_and_model_traces = table;
    data_and_model_traces.img_trace_axon = img_trace_axon;
    data_and_model_traces.img_trace_prediction = img_trace_prediction;
    writetable(data_and_model_traces, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'axon_data_and_model')

    dendritic_recordings = table;
    dendritic_recordings.dendritic_trace = img_trace_dendrites;
    writetable(dendritic_recordings, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'img_traces_dendrites')
  
    data_vs_model = table;
    data_vs_model.data_vec = [NaN; xy_vec'];
    data_vs_model.model_vec = [xy_vec; m_scat_hist];
    writetable(data_vs_model, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'data_vs_model')

end


end


function curr_hist = get_data_vs_model_axon_dendrite(ad_corr_table)

for loop_idx = 1:size(ad_corr_table,1)
    
    curr_hist(:,:,loop_idx) = ad_corr_table.axo_dendritic_corr(loop_idx,1).scatter_pred_hist;
    
end

curr_hist = mean(curr_hist,3);

end
