function get_fig_panel_power_spectrum_SD(mp_data_table, x_pos, y_pos, ...
    selected_ZT, sd_color, ctrl_color, source_data_details)
% Plotting the power spectra of dFBN GCaMP signals from rested and
% sleep-deprived flies.
% Code written by Peter Hasenhuetl.

sz_1 = 1.8;
sz_2 = 0.4;
ht_1 = 1.8;
ht_2 = 1.4;

xlm_delta_scatter = [0.5, 2.5];
ylm_delta_scatter = [0, 0.4];
scale_fac = 10^3;
xlm_f = [0, 2];

cond_SD = "SNAP";
if cond_SD == "SNAP"
    ylm_pwr = [0, 0.001]*scale_fac;
    major_y_ticks_scatter = 0.100;
    major_power_tick = 0.2;
elseif cond_SD == "vortex"
    ylm_pwr = [0, 0.001]*scale_fac;
    major_y_ticks_scatter = 0.100;
    major_power_tick = 0.5;
end

ac_rate = mp_data_table.frqu_trc(1,1).trace_img;
ZT_vec = 0:4:28;

mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n6",:) = "no";
mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n3",:) = "no";
mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n9",:) = "no";


curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "no" & ...
    mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time <= ZT_vec(selected_ZT+1),:);
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:400;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
SP_ctrl = [];
IDP_ctrl = [];
IIP_ctrl = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);
    [n_pwr1, n_delta_power1, n_infra_power1, ~] = get_n_pwr_SD_dFB(curr_img_tbl);
    SP_ctrl = [SP_ctrl, mean(n_pwr1,2,'omitnan')]; %
    IDP_ctrl = [IDP_ctrl, mean(n_delta_power1,2,'omitnan')]; %
    IIP_ctrl = [IIP_ctrl, mean(n_infra_power1,2,'omitnan')]; %
end

curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "yes" & ...
    mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time <= ZT_vec(selected_ZT+1),:);
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:400;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
SP_SD = [];
IDP_SD = [];
IIP_SD = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);    
    [n_pwr2, n_delta_power2, n_infra_power2, ~] = get_n_pwr_SD_dFB(curr_img_tbl);
    SP_SD = [SP_SD, mean(n_pwr2,2,'omitnan')]; %
    IDP_SD = [IDP_SD, mean(n_delta_power2,2,'omitnan')]; %
    IIP_SD = [IIP_SD, mean(n_infra_power2,2,'omitnan')]; %
end

curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, SP_SD*scale_fac, ac_rate', sd_color)
get_default_SEM_area_plot(curr_pnl, SP_ctrl*scale_fac, ac_rate', ctrl_color)
get_default_separated_ax(curr_pnl, xlm_f(1), xlm_f(2), xlm_f(1), xlm_f(2), ...
    ylm_pwr(1), ylm_pwr(2), ylm_pwr(1), ylm_pwr(2),...
    1, major_power_tick,"linear", "linear",'Frequency (Hz)',{'(\DeltaF/F)^2 x10^{-3}'},...
    'k','k','k','k', sz_1, ht_1)

scatter_cond = 1;
if scatter_cond == 1
    curr_pnl = axes('Units','Centimeters','Position',[x_pos+sz_1-sz_2, y_pos+ht_1-ht_2, sz_2, ht_2]);
    hold on
    get_default_scatter_group(curr_pnl, IDP_ctrl, 1, ctrl_color)
    get_default_scatter_group(curr_pnl, IDP_SD, 2, sd_color)
    get_default_axis_limits_warning([IDP_SD, IDP_ctrl], ylm_delta_scatter)
    get_default_ax(curr_pnl, xlm_delta_scatter(1), xlm_delta_scatter(2), ...
    xlm_delta_scatter(1), xlm_delta_scatter(2), ...
    ylm_delta_scatter(1), ylm_delta_scatter(2), ylm_delta_scatter(1), ylm_delta_scatter(2),...
    [], major_y_ticks_scatter,"linear", "linear",[],{'Mean SWA power'},...
    'none','k','k','k', sz_1, ht_1)
    n_vector = [length(IDP_ctrl), length(IDP_SD)];
    get_default_scatter_n_numbers_annotation(x_pos+sz_1-sz_2, y_pos+ht_1-ht_2, sz_2, n_vector, [], 'k', 'italic')
end

[p_ranksum_pwr,~] = ranksum(IDP_ctrl, IDP_SD)
[~,p_ttest_pwr] = ttest2(IDP_ctrl, IDP_SD)


[p_ranksum_infra_pwr,~] = ranksum(IIP_ctrl, IIP_SD)
[~,p_ttest_infra_pwr] = ttest2(IIP_ctrl, IIP_SD)

%% Saves the source data
save_it = 1;
if save_it == 1
power_spectrum_SD = table;
power_spectrum_SD.x_values_frequ = ac_rate;
power_spectrum_SD.power_SD = SP_SD;
writetable(power_spectrum_SD, ...
    [source_data_details.data_path, source_data_details.file_name],'Sheet','power_spectrum_SD')

power_spectrum_ctrl = table;
power_spectrum_ctrl.x_values_frequ = ac_rate;
power_spectrum_ctrl.power_ctrl = SP_ctrl;
writetable(power_spectrum_ctrl, ...
    [source_data_details.data_path, source_data_details.file_name],'Sheet','power_spectrum_ctrl')

mean_power = NaN(100,2);
mean_power(1:length(IDP_SD),1) = IDP_SD';
mean_power(1:length(IDP_ctrl),2) = IDP_ctrl';
mean_power = array2table(mean_power);
mean_power.Properties.VariableNames = {'SD','ctrl'};
writetable(mean_power, ...
    [source_data_details.data_path, source_data_details.file_name],'Sheet','mean_power')
end
end


function [n_pwr, n_delta_power, n_infra_power, curr_size] = get_n_pwr_SD_dFB(data_tbl)
   
curr_table_pwr = data_tbl.power_spectrum;
curr_delta_table = data_tbl.int_delta_power;
curr_infra_table = data_tbl.int_infra_slow_power;
excl_vec = data_tbl.excl_idx;

curr_pwr = [];
curr_delta_power = [];
curr_infra_power = [];
curr_size = [];

for loop_idx1 = 1:size(data_tbl,1)
    loop_cond = 1;
    if excl_vec(loop_idx1,3) ~= 0 && excl_vec(loop_idx1,4) ~= 0
        loop_cond = 0;
    elseif excl_vec(loop_idx1,3) ~= 0 && excl_vec(loop_idx1,4) == 0
        idx_1 = 4;
    elseif excl_vec(loop_idx1,4) ~= 0 && excl_vec(loop_idx1,3) == 0
        idx_1 = 3;
    elseif excl_vec(loop_idx1,3) == 0 && excl_vec(loop_idx1,4) == 0
        idx_1 = [3,4];
    end
    if loop_cond == 1
        curr_pwr = [curr_pwr, curr_table_pwr(loop_idx1,1).trace_img(:,idx_1)];
        curr_delta_power = [curr_delta_power, curr_delta_table(loop_idx1,idx_1)];
        curr_infra_power = [curr_infra_power, curr_infra_table(loop_idx1,idx_1)];
        curr_size = [curr_size, size(curr_table_pwr(loop_idx1,1).trace_img(:,idx_1),2)];
    end
end
curr_size = sum(curr_size);

n = 1;
for loop_idx2 = 1:size(curr_pwr,2)
    if sum(curr_pwr(:,loop_idx2)) > 0
        n_pwr(:,n) = (curr_pwr(:,loop_idx2));
        n_delta_power(1,n) = curr_delta_power(loop_idx2);
        n_infra_power(1,n) = curr_infra_power(loop_idx2);
        n = n+1;
    end
    
end

end

