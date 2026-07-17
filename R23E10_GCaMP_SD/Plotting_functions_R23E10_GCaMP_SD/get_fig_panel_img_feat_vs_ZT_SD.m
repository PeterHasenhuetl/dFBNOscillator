function get_fig_panel_img_feat_vs_ZT_SD(mp_data_table, x_pos, y_pos, ...
    sd_color, ctrl_color, pwr_cond, source_data_details)
% Plotting dFBN SWA power or GCaMP transient amplitudes at different
% zeitgeber times. In rested and sleep-deprived flies.
% Code written by Peter Hasenhuetl.

sz_1_bar = 1.8;
ht_1 = 1.8;
error_bar_width = get_default_error_bar_width;
ht_ZT = 0.1;
y_dist = ht_ZT;
ZT_lim = [0, 24];
ylim_ZT = [1.5, 2.5];

if pwr_cond == "pwr"
    source_data_sheet_name_SD = 'power_vs_ZT_SD';
    source_data_sheet_name_C = 'power_vs_ZT_C';
    ylm_1 = [0, 0.18];
    major_y = 0.06;
    y_label = 'Mean SWA power';
elseif pwr_cond == "transients"
    ylm_1 = [0.15, 0.3];
    major_y = 0.05;
    y_label = 'Mean transient amplitude';
    source_data_sheet_name_SD = 'trans_vs_ZT_SD';
    source_data_sheet_name_C = 'trans_vs_ZT_C';
end


ZT_idx = [2, 6, 10, 14, 18, 24];
ZT_array = [0,4; 4,8; 8,12; 12,16; 16,23; 23,30];

mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n3",:) = "no";
mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n6",:) = "no";

feature_array_SD = NaN(100,length(ZT_idx));
feature_array_ctrl = NaN(100,length(ZT_idx));

n_vector = [];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1_bar, ht_1]);
hold on
for loop_idx = 1:length(ZT_idx)
    hold on
    curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "no"& ...
        mp_data_table.zeitgeber_time >= ZT_array(loop_idx,1) & ...
        mp_data_table.zeitgeber_time < ZT_array(loop_idx,2),:);
    n_flies = round(size(curr_tbl,1)/4);
    num_array = 1:400;
    curr_fly_ids = curr_tbl.fly_id(:,:);
    id_array = ismember(num_array,curr_fly_ids);

    ind_fly_ids = num_array(id_array == 1);
    SP_ctrl = [];
    curr_feature_ctrl = [];
    for fly_loop_idx = 1:n_flies
        curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(fly_loop_idx),:);
        if pwr_cond == "pwr"
            [n_pwr1, n_delta_power1, ~, ~] = get_n_pwr_vs_ZT(curr_img_tbl);
        elseif pwr_cond == "transients"
            [n_pwr1, n_delta_power1] = get_n_img(curr_img_tbl);
        end

        if isempty(n_pwr1) == 0
            SP_ctrl = [SP_ctrl, mean(n_pwr1,2,'omitnan')];
            curr_feature_ctrl = [curr_feature_ctrl, mean(n_delta_power1,2,'omitnan')];
        end
    end
    
    curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "yes" &...
        mp_data_table.zeitgeber_time > ZT_array(loop_idx,1) & ...
        mp_data_table.zeitgeber_time < ZT_array(loop_idx,2),:);
    n_flies = round(size(curr_tbl,1)/4);
    num_array = 1:400;
    curr_fly_ids = curr_tbl.fly_id(:,:);
    id_array = ismember(num_array,curr_fly_ids);
    ind_fly_ids = num_array(id_array == 1);
    SP_SD = [];
    curr_feature_SD = [];
    for fly_loop_idx = 1:n_flies
        curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(fly_loop_idx),:);
        if pwr_cond == "pwr"
            [n_pwr2, n_delta_power2, ~, ~] = get_n_pwr_vs_ZT(curr_img_tbl);
        elseif pwr_cond == "transients"
            [n_pwr2, n_delta_power2] = get_n_img(curr_img_tbl);
        end
        if isempty(n_pwr2) == 0
            SP_SD = [SP_SD, mean(n_pwr2,2,'omitnan')];
            curr_feature_SD = [curr_feature_SD, mean(n_delta_power2,2,'omitnan')];
        end
    end


   

    mean_array_SD(loop_idx) = mean(curr_feature_SD);
    mean_array_ctrl(loop_idx) = mean(curr_feature_ctrl);

    feature_array_SD(1:length(curr_feature_SD),loop_idx) = curr_feature_SD';
    feature_array_ctrl(1:length(curr_feature_ctrl),loop_idx) = curr_feature_ctrl';
 
    x_axis_coordinate = ZT_idx(loop_idx);
    get_default_scatter_and_whisker(curr_pnl, x_axis_coordinate, curr_feature_ctrl, ctrl_color)
    get_default_scatter_and_whisker(curr_pnl, x_axis_coordinate, curr_feature_SD, sd_color)

    empty_space = NaN;
    n_vector = [n_vector, length(curr_feature_ctrl), length(curr_feature_SD), empty_space];

end


plot(curr_pnl, ZT_idx, mean_array_ctrl, 'Color', ctrl_color, 'LineWidth', error_bar_width)
plot(curr_pnl, ZT_idx, mean_array_SD, 'Color', sd_color, 'LineWidth', error_bar_width)

get_default_separated_ax(curr_pnl, ZT_lim(1), ZT_lim(2), [], [], ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), [], major_y, "linear", "linear", [], string(y_label),...
    'none', 'k', 'none', 'k', sz_1_bar, ht_1)

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-y_dist, sz_1_bar, ht_ZT]);
ZT_array = [ones(1200,3); zeros(1200,3); ones(800,3)];
x_vec = linspace(0,24+8,length(ZT_array));
imagesc(curr_pnl,x_vec, 1:3, ZT_array')
colormap(curr_pnl, [0,0,0; 1,1,1])
xlim(ZT_lim)
ylim(ylim_ZT)
curr_pnl.Box = 'on';
curr_pnl.XTick = [];
curr_pnl.YTick = [];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-y_dist, sz_1_bar, ht_ZT]);
get_default_ax(curr_pnl, ZT_lim(1), ZT_lim(2), ZT_lim(1), ZT_lim(2), ylim_ZT(1), ylim_ZT(2),...
    [], [], 12, [], "linear", "linear", 'Zeitgeber time', [],...
    'k', 'k', 'k', 'none', sz_1_bar, ht_ZT)

n_cond = 0;
if n_cond == 1
    get_default_scatter_n_numbers_annotation(x_pos, y_pos+ht_1+0.3, sz_1_bar, n_vector(1:end-1), [], 'k', 'italic')
end
%% Saves the source data

if isempty(source_data_details) == 0

    x_values = ZT_idx(1:loop_idx);
    size(x_values)
    power_vs_ZT_SD = table;
    power_vs_ZT_SD.x_values = x_values';
    power_vs_ZT_SD.power_SD = feature_array_SD';
    writetable(power_vs_ZT_SD, [source_data_details.data_path, ...
        source_data_details.file_name],'Sheet',source_data_sheet_name_SD)
    
    power_vs_ZT_C = table;
    power_vs_ZT_C.x_values = x_values';
    power_vs_ZT_C.C = feature_array_ctrl';
    writetable(power_vs_ZT_C, [source_data_details.data_path, ...
        source_data_details.file_name],'Sheet',source_data_sheet_name_C)

end
end


function [n_pwr, n_delta_power, n_infra_power, curr_size] = get_n_pwr_vs_ZT(data_tbl)
   
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
        idx_1 = [3, 4];
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


function [trans_traces, amp_img] = get_n_img(data_tbl)


n_wind = 15;
p_wind = 40;
cum_idx = 1;
excl_vec = data_tbl.excl_idx;
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
        for yoyo = 1:length(idx_1)
            curr_id = idx_1(yoyo);
          
            curr_onsets_nu = data_tbl.trans_ids(i,1).transient_onsets(20:end-50,curr_id);
            curr_img_trace = data_tbl.img_raw(i,1).trace_img(:,curr_id);

            for loop_idx_trans_1 = 1:length(curr_onsets_nu)
                
                if curr_onsets_nu(loop_idx_trans_1)+p_wind < length(curr_img_trace) && curr_onsets_nu(loop_idx_trans_1)-n_wind > 1
                    trans_trace_raw = (curr_img_trace(curr_onsets_nu(loop_idx_trans_1)-n_wind:curr_onsets_nu(loop_idx_trans_1)+p_wind,1));  
                    F0 = mean(trans_trace_raw(1:n_wind,1));
                    trans_trace_dF = (trans_trace_raw-F0)./F0;
                    
                    trans_traces(:,cum_idx) = trans_trace_dF;
                    cum_idx = cum_idx+1;
                end
            end
            
             
        end
    end
    
end


wind_1 = 15;
for loop_idx1 = 1:size(trans_traces,2)
    amp_img(loop_idx1) = max(trans_traces(wind_1+1:wind_1+15,loop_idx1));
end


trans_traces = mean(trans_traces,2);
amp_img = mean(amp_img);



end

