function get_fig_panel_dFB_tuning_curve_SD(mp_data_table, x_pos, y_pos, ...
    source_data_details)
% Plotting dFBN SWA power with varying levels of sleep pressure.
% Code written by Peter Hasenhuetl.

plotting_color = [0, 0, 0];
sz_1_bar = 1.8;
ht_1 = 1.8;
error_bar_width = get_default_error_bar_width;
norm_cond = 0;
if norm_cond == 1
    ylm_1 = [-0.3, 0.9];
    major_y_tick = 0.3;
    y_label = "Fold change in SWA power";
else
    ylm_1 = [0.6, 1.8];
    major_y_tick = 0.4;
    y_label = "Normalized SWA power";
end
xlm_1 = [0, 12];
major_x_tick = 3;
plot_sz = sz_1_bar;
SD_idx = [0, 3, 6, 9, 12];
ZT_threshold = 4;

mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n6",:) = "no";
mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n3",:) = "no";

feature_array_ctrl = NaN(100,length(SD_idx));

n_vector = [];
curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos, plot_sz, ht_1]);
hold on
for loop_idx = 1:length(SD_idx)
    hold on
    
    if SD_idx(loop_idx) == 0
        curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "no" & ...
            mp_data_table.zeitgeber_time < ZT_threshold,:);
    elseif SD_idx(loop_idx) == 3
        curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "y3" & ...
            mp_data_table.zeitgeber_time < ZT_threshold,:);
    elseif SD_idx(loop_idx) == 6
        curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "y6" & ...
            mp_data_table.zeitgeber_time < ZT_threshold,:);
    elseif SD_idx(loop_idx) == 9
        curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "y9" & ...
            mp_data_table.zeitgeber_time < ZT_threshold,:);
    elseif SD_idx(loop_idx) == 12
        curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "yes" & ...
            mp_data_table.zeitgeber_time < ZT_threshold,:);
    end
    
    n_flies = round(size(curr_tbl,1)/4);
    num_array = 1:400;
    curr_fly_ids = curr_tbl.fly_id(:,:);
    id_array = ismember(num_array,curr_fly_ids);

    ind_fly_ids = num_array(id_array == 1);
    SP_ctrl = [];
    curr_feature_ctrl = [];
    IIP_ctrl = [];
    for fly_loop_idx = 1:n_flies
        curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(fly_loop_idx),:);
        [n_pwr1, n_delta_power1, n_infra_power1, ~] = get_n_pwr_vs_SD_dur(curr_img_tbl);
        if isempty(n_pwr1) == 0
            SP_ctrl = [SP_ctrl, mean(n_pwr1,2,'omitnan')];
            curr_feature_ctrl = [curr_feature_ctrl, mean(n_delta_power1,2,'omitnan')];
            IIP_ctrl = [IIP_ctrl, mean(n_infra_power1,2,'omitnan')];
        end
    end
    
    mean_array_ctrl(loop_idx) = mean(curr_feature_ctrl);
    if loop_idx == 1
        mean_array_baseline = mean_array_ctrl(loop_idx);
    end
    
    if norm_cond == 1
        curr_feature_ctrl = (curr_feature_ctrl - mean_array_baseline)./mean_array_baseline;
    else
        curr_feature_ctrl = (curr_feature_ctrl./mean_array_baseline);
    end

    if loop_idx ==2
        curr_feature_ctrl(2) = [];
    end
    
    feature_array_ctrl(1:length(curr_feature_ctrl),loop_idx) = curr_feature_ctrl;

    mean_array_ctrl(loop_idx) = mean(curr_feature_ctrl);
    x_axis_coordinate = SD_idx(loop_idx);
    get_default_scatter_and_whisker(curr_pnl, x_axis_coordinate, curr_feature_ctrl, plotting_color)

    empty_space = NaN;
    n_vector = [n_vector, length(curr_feature_ctrl), empty_space];

end

plot(curr_pnl, SD_idx, mean_array_ctrl, 'Color', plotting_color, 'LineWidth', error_bar_width)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), major_x_tick, major_y_tick, "linear", "linear", 'SD duration (h)', string(y_label),...
    'k', 'k', 'k', 'k', plot_sz, ht_1)

n_cond = 1;
if n_cond == 1
    get_default_scatter_n_numbers_annotation(x_pos, y_pos+ht_1+0.3, plot_sz, ...
        n_vector(1:end-1), [], 'k', 'italic')
end
%% Saves the source data

if isempty(source_data_details) == 0

    x_values = SD_idx(1:loop_idx);
    power_array = table;
    power_array.x_values = x_values';
    power_array.C = feature_array_ctrl';
    writetable(power_array, [source_data_details.data_path, ...
        source_data_details.file_name], 'Sheet', 'SD_tuning_curve')

end

end


function [n_pwr, n_delta_power, n_infra_power, curr_size] = get_n_pwr_vs_SD_dur(data_tbl)
   
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

n = 1;
n_pwr = [];
n_delta_power = [];
n_infra_power = [];
for j = 1:size(curr_pwr,2)
    
    if sum(curr_pwr(:,j)) > 0
        
        n_pwr(:,n) = (curr_pwr(:,j));
        n_delta_power(1,n) = curr_delta_power(j);
        n_infra_power(1,n) = curr_infra_power(j);
        n = n+1;
    
    end
    
end

end

