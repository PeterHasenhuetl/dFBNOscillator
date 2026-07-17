function get_fig_panel_intro_Fig1(mp_data_table, x_pos, y_pos, ...
    color, start_trace1, end_trace1, source_data_details)
% Plotting an example recording of dFBN GCaMP dynamics and the pooled power spectrum of all recordings.
% Code written by Peter Hasenhuetl.

curr_fly = mp_data_table(mp_data_table.sleep_deprivation == "no",:);

example_fly_id = 7;
plane_idx_1 = 2;
dendrite_idx_1 = 4;

sz_1 = 7;
ht_1 = 1;
ht_inset = 1;
x_dist_cartoon = 1.2;
x_dist_spectrum = 2.3;
sz_inset = sz_1-x_dist_cartoon-x_dist_spectrum;
plot_color = [0, 0, 0];
y_dist_inset = 0.5;
ylim_1 = [-0.075, 1.3];
line_width_plot = 0.25;
y_axis_scale_bar = 0.3;

curr_tbl = curr_fly(curr_fly.fly_id == example_fly_id,:);
curr_trace = curr_tbl.img_trace(plane_idx_1,1).trace_img(:,dendrite_idx_1);

curr_pnl1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl1, [start_trace1, start_trace1], ylim_1, 'r')
plot(curr_pnl1, [end_trace1, end_trace1], ylim_1, 'r')
get_default_img_trace_plot(curr_pnl1, curr_trace, 14.56, ...
    60, y_axis_scale_bar, plot_color, line_width_plot,  "with_y_scale_bar", ...
    ylim_1, [0, 1], y_axis_scale_bar, "1 min", 0, ["0.3"; "\DeltaF/F"], 0.25);

curr_trace = curr_trace(start_trace1:end_trace1,1);
y_pos_inset = y_pos+ht_1+y_dist_inset;
curr_pnl2 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+x_dist_cartoon, y_pos_inset, sz_inset, ht_inset]);
get_default_img_trace_plot(curr_pnl2, curr_trace, 14.56, ...
    10, y_axis_scale_bar, plot_color, line_width_plot, "with_y_scale_bar", ...
    [min(curr_trace)-0.01, max(curr_trace)+0.01], [0, 1], y_axis_scale_bar, "10 s", 0, [], 0.25)

sz_pwr_panel = 1;
ht_pwr_panel = ht_inset;
get_fig_panel_pwr_full(mp_data_table, x_pos+sz_1-sz_pwr_panel, y_pos_inset, sz_pwr_panel,...
    ht_pwr_panel, source_data_details)

%% Saves the source data
if isempty(source_data_details) == 0
    example_imaging_trace = array2table(curr_trace);
    writetable(example_imaging_trace, ...
        [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'example_imaging_trace')

end

end

function get_fig_panel_pwr_full(mp_data_table, x_pos, y_pos, sz_1, ht_1, source_data_details)


pwr_spec_color = [0, 0, 0];

scale_fac = 10^3;
xlm_f = [0, 2];
ylm_pwr = [0, 0.001]*scale_fac;

ac_rate = mp_data_table.frqu_trc(1,1).trace_img;

curr_tbl = mp_data_table;
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:2000;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
SP_full = [];
IDP_full = [];
IIP_full = [];
cum_size = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);
    [n_pwr1, n_delta_power1, n_infra_power1, curr_size] = get_n_pwr(curr_img_tbl);
    if isempty(n_pwr1) == 0
        SP_full = [SP_full, mean(n_pwr1,2,'omitnan')];
        IDP_full = [IDP_full, mean(n_delta_power1,2,'omitnan')];
        IIP_full = [IIP_full, mean(n_infra_power1,2,'omitnan')];
        cum_size = [cum_size, curr_size];
    end
end

n_ROIs = sum(cum_size);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, SP_full*scale_fac, ac_rate', pwr_spec_color)
get_default_separated_ax(curr_pnl, xlm_f(1), xlm_f(2), xlm_f(1), xlm_f(2), ...
    ylm_pwr(1), ylm_pwr(2), ylm_pwr(1), ylm_pwr(2),...
    1, 0.5, "linear", "linear", 'Frequency (Hz)', {'(\DeltaF/F)^2 x10^{-3}'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)

get_default_annotation(x_pos+sz_1, y_pos+ht_1-0.2, [num2str(n_flies), ' flies'], 'k', 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos+ht_1, [num2str(n_ROIs), ' ROIs'], 'k', 'normal', "right")

    
%% Saves the source data
if isempty(source_data_details) == 0
    power_spectrum_full = table;
    power_spectrum_full.x_values_frequ = ac_rate;
    power_spectrum_full.SP_full = SP_full;
    writetable(power_spectrum_full, ...
        [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'power_spectrum_full')
end
end

function [n_pwr, n_delta_power, n_infra_power, curr_size] = get_n_pwr(data_tbl)
   
curr_table_pwr = data_tbl.power_spectrum;
curr_delta_table = data_tbl.int_delta_power;
curr_infra_table = data_tbl.int_infra_slow_power;
excl_vec = data_tbl.excl_idx;

curr_pwr = [];
curr_delta_power = [];
curr_infra_power = [];
curr_size = [];

for i = 1:size(data_tbl,1)
    loop_cond = 1;
    if excl_vec(i,3) ~= 0 && excl_vec(i,4) ~= 0
        loop_cond = 0;
    elseif excl_vec(i,3) ~= 0 && excl_vec(i,4) == 0
        idx_1 = 4;
    elseif excl_vec(i,4) ~= 0 && excl_vec(i,3) == 0
        idx_1 = 3;
    elseif excl_vec(i,3) == 0 && excl_vec(i,4) == 0
        idx_1 = [3,4];
    end
    if loop_cond == 1
        curr_pwr = [curr_pwr, curr_table_pwr(i,1).trace_img(:,idx_1)];
        curr_delta_power = [curr_delta_power, curr_delta_table(i,idx_1)];
        curr_infra_power = [curr_infra_power, curr_infra_table(i,idx_1)];
        curr_size = [curr_size, size(curr_table_pwr(i,1).trace_img(:,idx_1),2)];
    end
end
curr_size = sum(curr_size);

n_pwr = [];
n_delta_power = [];
n_infra_power = [];
n = 1;
for j = 1:size(curr_pwr,2)
    if sum(curr_pwr(:,j)) > 0
        
    n_pwr(:,n) = (curr_pwr(:,j));
    n_delta_power(1,n) = curr_delta_power(j);
    n_infra_power(1,n) = curr_infra_power(j);
    n = n+1;
    end
    
end

end
