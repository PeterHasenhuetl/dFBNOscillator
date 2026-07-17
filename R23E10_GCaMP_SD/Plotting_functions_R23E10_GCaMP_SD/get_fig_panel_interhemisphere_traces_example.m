function get_fig_panel_interhemisphere_traces_example(mp_data_table, ...
    x_pos, y_pos, plane_id1, plane_id2, r_1, r_2, plotting_color_1, plotting_color_2,...
    image_1, image_2, color_map_1, source_data_details)
% Plotting an example bilateral recording of dFBN GCaMP dynamics.
% Code written by Peter Hasenhuetl.


plotting_color_left = plotting_color_1;
plotting_color_right = plotting_color_2;

p_wind = round(14.56*1.75);
n_wind = 15;

sz_1 = 1.25;
ht_ipsi = 1.25;
ht_contra = 1.25;
ht_2 = 2;
x_dist_between_images = 1.25;
x_dist_between_traces = 1;
size_cartoon = 1;
x_pos_aligned_traces = x_pos;
x_pos_full_trace = x_pos+size_cartoon;
y_pos_aligned_traces = y_pos;
y_dist_images = 0;
y_dist_full_trace = 1;

sz_image = (2*sz_1)+x_dist_between_traces;
h_t_norm = size(image_1,1)/size(image_1,2);
h_t_image = sz_image*h_t_norm;

y_pos_full_trace = y_pos+ht_ipsi+h_t_image+y_dist_full_trace+y_dist_images;
sz_2 = (2*sz_image)+x_dist_between_images-size_cartoon;

ylm_2 = [-2, 7];
ylm_ipsi = [-0.05, 0.25];
ylm_contra = [-0.05, 0.025];
clim_1 = ([-0.1, 0.5])*1;

curr_trace_left = mp_data_table.img_trace(plane_id1,1).trace_img(:,3);
curr_trace_right = mp_data_table.img_trace(plane_id2,1).trace_img(:,4);

curr_trace_left_raw = mp_data_table.img_raw(plane_id1,1).trace_img(:,3);
curr_trace_right_raw = mp_data_table.img_raw(plane_id2,1).trace_img(:,4);

curr_burst_left = mp_data_table.transient_vec(plane_id1,1).trace_img(:,3);
curr_burst_right = mp_data_table.transient_vec(plane_id2,1).trace_img(:,4);

y = find(curr_burst_left > 0);
F0_wind = round(n_wind);   

n = 1;
for i = 4:length(y)-5
    i_trace_left = (curr_trace_left_raw(y(i)-n_wind:y(i)+p_wind,1));
    ipsi_trace_left(:,n) = (i_trace_left-mean(i_trace_left(1:F0_wind,1)))/mean(i_trace_left(1:F0_wind,1));
    c_trace_right = (curr_trace_right_raw(y(i)-n_wind:y(i)+p_wind,1));
    contra_trace_right(:,n) = (c_trace_right-mean(c_trace_right(1:F0_wind,1)))/mean(c_trace_right(1:F0_wind,1));
    n = n+1;
end

ipsi_trace_left = ((ipsi_trace_left));
contra_trace_right = ((contra_trace_right));
n_transients_1 = n;
  
y = find(curr_burst_right > 0);
n = 1;
for i = 4:length(y)-5  
    i_trace_right = (curr_trace_right_raw(y(i)-n_wind:y(i)+p_wind,1));
    ipsi_trace_right(:,n) = (i_trace_right-mean(i_trace_right(1:F0_wind,1)))/mean(i_trace_right(1:F0_wind,1));
    c_trace_left = (curr_trace_left_raw(y(i)-n_wind:y(i)+p_wind,1));
    contra_trace_left(:,n) = (c_trace_left-mean(c_trace_left(1:F0_wind,1)))/mean(c_trace_left(1:F0_wind,1));
    n = n+1;
end

n_transients_2 = n;
    
%% plots the aligned transients

y_dist_annotation = 0.3;

get_fig_panel_transient_aligned_contra(ipsi_trace_left, contra_trace_right, ...
    x_pos_aligned_traces, y_pos_aligned_traces, sz_1, ht_ipsi, ht_contra, x_dist_between_traces, plotting_color_left, plotting_color_right,...
    ylm_ipsi, ylm_contra, 'none', "left")
get_default_annotation(x_pos_aligned_traces, y_pos_aligned_traces+ht_ipsi+y_dist_annotation+y_dist_images+h_t_image, ...
    [num2str(n_transients_1),' transients (L)'], plotting_color_left, 'normal', "left")
get_default_annotation(x_pos_aligned_traces+sz_image, y_pos_aligned_traces+ht_ipsi+y_dist_annotation+y_dist_images+h_t_image, ...
    'Contra (R)', plotting_color_right, 'normal', "right")


x_pos_second_panel = x_pos_aligned_traces+((sz_image)+x_dist_between_images);
get_fig_panel_transient_aligned_contra(ipsi_trace_right, contra_trace_left, ...
    x_pos_second_panel, y_pos_aligned_traces, sz_1, ht_ipsi, ht_contra, x_dist_between_traces, plotting_color_right, plotting_color_left,...
    ylm_ipsi, ylm_contra, 'none', "right")
get_default_annotation(x_pos_second_panel+sz_image, y_pos_aligned_traces+ht_ipsi+y_dist_annotation+y_dist_images+h_t_image, ...
    [num2str(n_transients_2),' transients (R)'], plotting_color_right, 'normal', "right")
get_default_annotation(x_pos_second_panel, y_pos_aligned_traces+ht_ipsi+y_dist_annotation+y_dist_images+h_t_image, ...
    'Contra (L)', plotting_color_left, 'normal', "left")

%% Plots the images

image_sb = 28;
pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos_aligned_traces, y_pos_aligned_traces+ht_ipsi+y_dist_images, sz_image, h_t_image]);
imagesc(pnl_1, image_1)
hold on
plot(pnl_1, [size(image_1,2)-image_sb, size(image_1,2)],...
    [size(image_1,1), size(image_1,1)], 'k', 'LineWidth', get_default_scale_bar_width)
colormap(color_map_1)
clim(clim_1)
pnl_1.Box = 'off';
pnl_1.Color = 'none';
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
xlim(pnl_1, [1, size(image_1,2)])
ylim(pnl_1, [1, size(image_1,1)])

scale_annotation_cond = 1;
if scale_annotation_cond == 1
    get_default_annotation(x_pos_aligned_traces+sz_image, ...
        y_pos_aligned_traces+ht_ipsi+y_dist_images-0.1, ['30 \mu','m'], 'k', 'normal', "right")
end

pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos_second_panel, y_pos_aligned_traces+ht_ipsi+y_dist_images, sz_image, h_t_image]);
imagesc(pnl_1, image_2)
colormap(color_map_1)
clim(clim_1)
pnl_1.Box = 'off';
pnl_1.Color = 'none';
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
xlim(pnl_1, [1, size(image_2,2)])
ylim(pnl_1, [1, size(image_2,1)])


ht_color_bar = h_t_image;
ac = colorbar(pnl_1,'eastoutside');
ac.Units = 'centimeters';
ac.Position = [x_pos_second_panel+sz_image+0.1, y_pos_aligned_traces+ht_ipsi+y_dist_images, 0.1, ht_color_bar];
ylabel(ac,'\DeltaF/F')
ac.FontSize = 5;
ac.Box = 'off';
ac.Color = 'k';
ac.TickLength = 0;
ac.YLabel.Visible = 'on';
ac.YLabel.Color = 'k';

%% Plots the z-scored traces

curr_trace_left = zscore(curr_trace_left);
curr_trace_right = zscore(curr_trace_right);

curr_trace_left = curr_trace_left(r_1:r_2,:);
curr_trace_right = curr_trace_right(r_1:r_2,:);

scale_bar = [length(curr_trace_left)-(14.56*5), length(curr_trace_left)];
xlm_2 = [1, length(curr_trace_left)];

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos_full_trace, y_pos_full_trace, sz_2, ht_2]);
hold on
plot(pnl_1, curr_trace_right, 'Color', plotting_color_right, 'LineWidth', 0.5, 'LineStyle', '-')
plot(pnl_1, curr_trace_left, 'Color', plotting_color_left, 'LineWidth', 0.5, 'LineStyle', '-')
plot(pnl_1, scale_bar, [-2, -2], 'k', 'LineWidth', get_default_scale_bar_width)
plot(pnl_1, [length(curr_trace_left), length(curr_trace_left)], [ylm_2(end)-3, ylm_2(end)-1],...
    'k', 'LineWidth', get_default_scale_bar_width)
ylim(ylm_2)
get_default_ax(pnl_1, xlm_2(1), xlm_2(2), [], [], ylm_2(1), ylm_2(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_2, ht_2)
get_default_annotation(x_pos_full_trace+sz_2-0.35, y_pos_full_trace-0.125, ...
    '5 s', 'k', 'normal', "left")
get_default_annotation(x_pos_full_trace+sz_2+0.2, y_pos_full_trace+ht_2-0.3, ...
    '2 s.d.', 'k', 'normal', "left")
get_default_annotation(x_pos_full_trace, y_pos_full_trace+ht_2+0.2, ...
    'Left hemisphere', plotting_color_left, 'normal', "left")

get_default_annotation(x_pos_full_trace+1.5, y_pos_full_trace+ht_2+0.2, ...
    'Right hemisphere', plotting_color_right, 'normal', "left")

%% Saves the source data

traces_left_leads_ipsi = table;
traces_left_leads_ipsi.ipsi_trace_left = ipsi_trace_left;
writetable(traces_left_leads_ipsi, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_left_leads_ipsi')

traces_left_leads_contra = table;
traces_left_leads_contra.contra_trace_right = contra_trace_right;
writetable(traces_left_leads_contra, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_left_leads_contra')

traces_right_leads_ipsi = table;
traces_right_leads_ipsi.ipsi_trace_right = ipsi_trace_right;
writetable(traces_right_leads_ipsi, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_right_leads_ipsi')

traces_right_leads_contra = table;
traces_right_leads_contra.contra_trace_left = contra_trace_left;
writetable(traces_right_leads_contra, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_right_leads_contra')

example_imaging_traces = table;
example_imaging_traces.left_hemisphere = curr_trace_left;
example_imaging_traces.right_hemisphere = curr_trace_right;
writetable(example_imaging_traces, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'example_imaging_traces')

end


function get_fig_panel_transient_aligned_contra(ipsi_trace, contra_trace, ...
    curr_x_pos, curr_y_pos, sz_1, ht_ipsi, ht_contra, x_dist_between_traces, plotting_color_ipsi, plotting_color_contra,...
    ylm_ipsi, ylm_contra, y_axis_color, left_right)

dist_1 = abs(ylm_ipsi(1)/(ylm_ipsi(1)-ylm_ipsi(end)));
norm_dist1 = ht_ipsi*dist_1;

dist_2 = abs(ylm_contra(1)/(ylm_contra(1)-ylm_contra(end)));
norm_dist2 = ht_ipsi*dist_2;

xlm_1 = [1, length(mean(ipsi_trace,2))];
t_v = xlm_1(1):xlm_1(2);

scale_bar = [xlm_1(2)-(0.5*14.56),xlm_1(2)];
ht_contra = ht_ipsi;

if left_right == "left"    
    curr_x_pos1 = curr_x_pos+sz_1+x_dist_between_traces;
    curr_x_pos2 = curr_x_pos;
elseif left_right == "right"
    curr_x_pos2 = curr_x_pos+sz_1+x_dist_between_traces;
    curr_x_pos1 = curr_x_pos;
    get_default_annotation(curr_x_pos2+sz_1+0.1, curr_y_pos+ht_contra-0.75, ...
    ['0.1 ','\DeltaF/F'], 'k', 'normal', "left")
    get_default_annotation(curr_x_pos1+sz_1+0.1, curr_y_pos-norm_dist2+0.3, ...
    ['0.02 ','\DeltaF/F'], 'k', 'normal', "left")  
end

pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos1, curr_y_pos-norm_dist2, sz_1, ht_contra]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, contra_trace, t_v, plotting_color_contra)
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_contra(1), ylm_contra(2),...
    0, ylm_contra(2), [], 0.1, "linear", "linear", [], 'z-score',...
    'none', y_axis_color, 'none',y_axis_color, sz_1, ht_contra)
if left_right == "right"
    plot(pnl_1, [xlm_1(2)-(0.25*14.56), xlm_1(2)-(0.25*14.56)],...
        [ylm_contra(1)+0.01, ylm_contra(1)+0.03], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
end


pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos2, curr_y_pos-norm_dist1, sz_1, ht_ipsi]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, ipsi_trace, t_v, plotting_color_ipsi)
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_ipsi(1), ylm_ipsi(2),...
    0, ylm_ipsi(2), [], 0.1, "linear", "linear", [], 'z-score',...
    'none', y_axis_color, 'none', y_axis_color, sz_1, ht_ipsi)
if left_right == "left"
    plot(pnl_1,scale_bar,...
        [ylm_ipsi(1), ylm_ipsi(1)], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
end
if left_right == "right"
    plot(pnl_1, [xlm_1(2)-(0.25*14.56), xlm_1(2)-(0.25*14.56)],...
        [ylm_ipsi(2)-0.2, ylm_ipsi(2)-0.1], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
end

if left_right == "left"
    get_default_annotation(curr_x_pos2+sz_1, curr_y_pos-0.3, ...
        '0.5 s', 'k', 'normal', "right")
end

triangle_color = [0.3, 0.3, 0.3];
pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos1, curr_y_pos-norm_dist1, sz_1, ht_ipsi]);
scatter(pnl_1, 13, (ylm_ipsi(2)*0.75), 5, 'Marker', 'v', 'MarkerFaceColor', triangle_color, 'MarkerEdgeColor', 'none')
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_ipsi(1), ylm_ipsi(2),...
    0, ylm_ipsi(2), [], 0.1, "linear", "linear", [], 'z-score',...
    'none', y_axis_color, 'none', y_axis_color, sz_1, ht_ipsi)

pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos2, curr_y_pos-norm_dist1, sz_1, ht_ipsi]);
scatter(pnl_1, 13, (ylm_ipsi(2)*0.75), 5, 'Marker', 'v', 'MarkerFaceColor', triangle_color, 'MarkerEdgeColor', 'none')
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_ipsi(1), ylm_ipsi(2),...
    0, ylm_ipsi(2), [], 0.1,"linear", "linear", [],'z-score',...
    'none', y_axis_color, 'none', y_axis_color, sz_1, ht_ipsi)

end
