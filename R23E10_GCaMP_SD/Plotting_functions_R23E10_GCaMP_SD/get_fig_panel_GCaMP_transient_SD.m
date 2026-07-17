function get_fig_panel_GCaMP_transient_SD(mp_data_table, x_pos, y_pos, selected_ZT, plot_cond, sd_color, ctrl_color, source_data_details)
% Plotting dFBN GCaMP transients after sleep deprivation.
% Code written by Peter Hasenhuetl.

ZT_vec = 0:4:28;
n_wind = 15;
xlm1 = [-0.5, 2.5];
xlm2 = [0, 3];
if plot_cond == 1
    
    cond_SD = "SNAP";
    if cond_SD == "SNAP"
        ylm_1 = [-0.05, 0.4];
    elseif cond_SD == "vortex"
        ylm_1 = [-0.075, 0.35];
    end
    
elseif plot_cond == 2
    ylm_1 = [0, 0.4];
end

ht_1 = 2.8;
sz_1 = 1;

mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n6",:) = "no";
mp_data_table.sleep_deprivation(mp_data_table.sleep_deprivation == "n3",:) = "no";


curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "no" & ...
    mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time < ZT_vec(selected_ZT+1),:);
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:1000;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
img_ctrl = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);    
    img_ctrl = [img_ctrl, mean(get_n_img(curr_img_tbl, n_wind, plot_cond),2,'omitnan')];
end

curr_tbl = mp_data_table(mp_data_table.sleep_deprivation == "yes" & ...
    mp_data_table.zeitgeber_time > ZT_vec(selected_ZT) & ...
    mp_data_table.zeitgeber_time < ZT_vec(selected_ZT+1),:);
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:1000;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
img_SD = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);    
    img_SD = [img_SD, mean(get_n_img(curr_img_tbl, n_wind, plot_cond),2,'omitnan')];
end

t_v = (0:size(img_ctrl,1)-1)./14.56;
t_v = t_v-t_v(n_wind);

%% traces

scale_bar = [xlm1(2)-0.5, xlm1(2)];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl, [xlm1(1), xlm1(2)], [0, 0], ':k')
plot(curr_pnl, scale_bar,...
    [ylm_1(1), ylm_1(1)], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
plot(curr_pnl, [xlm1(1), xlm1(1)],...
    [0.1, 0.2], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')
get_default_SEM_area_plot(curr_pnl, img_ctrl, t_v, ctrl_color)
get_default_SEM_area_plot(curr_pnl, img_SD, t_v, sd_color)
get_default_ax(curr_pnl, xlm1(1), xlm1(2), [], [], ylm_1(1), ylm_1(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1)

get_default_annotation(x_pos+sz_1-0.25, y_pos-0.1, '0.5 s', 'k', 'normal', "left")
get_default_annotation(x_pos-0.5, y_pos+0.3*ht_1, {'0.1';'\DeltaF/F'}, 'k','normal', "left")

%% scatter plot with amplitudes with same y-axis as plot of traces

for loop_idx1 = 1:size(img_ctrl,2)
    amp_img_ctrl(loop_idx1) = max(img_ctrl(n_wind+1:n_wind+15,loop_idx1));
end

for loop_idx2 = 1:size(img_SD,2)
    amp_img_SD(loop_idx2) = max(img_SD(n_wind+1:n_wind+15,loop_idx2));
end

x_pos_scatter = x_pos+sz_1;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos_scatter, y_pos, 0.75, ht_1]);
get_default_scatter_group(curr_pnl, amp_img_ctrl, 1, ctrl_color)
hold on
get_default_scatter_group(curr_pnl, amp_img_SD, 2, sd_color)
get_default_ax(curr_pnl, xlm2(1), xlm2(2), [], [], ylm_1(1), ylm_1(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.75, ht_1)
get_default_axis_limits_warning([amp_img_ctrl, amp_img_SD], ylm_1)


[p_ranksum_amp,~] = ranksum(amp_img_ctrl, amp_img_SD)
[~,p_ttest_amp] = ttest2(amp_img_ctrl, amp_img_SD)


%% Saves the source data
save_it = 1;
if save_it == 1
    traces_SD = table;
    traces_SD.img_SD = img_SD;
    writetable(traces_SD, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'imaging_traces_SD')

    traces_ctrl = table;
    traces_ctrl.img_ctrl = img_ctrl;
    writetable(traces_ctrl, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'imaging_traces_ctrl')

    peak_dF = NaN(100,2);
    peak_dF(1:length(amp_img_SD),1) = amp_img_SD;
    peak_dF(1:length(amp_img_ctrl),2) = amp_img_ctrl;
    peak_dF = array2table(peak_dF);
    peak_dF.Properties.VariableNames = {'SD','ctrl'};
    writetable(peak_dF, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'peak_dF')
end

end


function [trans_traces] = get_n_img(data_tbl, n_wind, plot_cond)

excl_vec = data_tbl.excl_idx;

p_wind = 40;
cum_idx = 1;
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
        for yoyo = 1:length(idx_1)
            curr_id = idx_1(yoyo);
          
            curr_onsets_nu = data_tbl.trans_ids(i,1).transient_onsets(20:end-50,curr_id);
            curr_img_trace = data_tbl.img_raw(i,1).trace_img(:,curr_id);

            for loop_idx_trans_1 = 1:length(curr_onsets_nu)
                
                if curr_onsets_nu(loop_idx_trans_1)+p_wind < length(curr_img_trace) & ...
                        curr_onsets_nu(loop_idx_trans_1)-n_wind > 1
                    trans_trace_raw = (curr_img_trace(curr_onsets_nu(loop_idx_trans_1) - ...
                        n_wind:curr_onsets_nu(loop_idx_trans_1)+p_wind,1));  
                    F0 = mean(trans_trace_raw(1:n_wind,1));
                    trans_trace_dF = (trans_trace_raw-F0)./F0;
                    
                    trans_traces(:,cum_idx) = trans_trace_dF;
                    cum_idx = cum_idx+1;
                end
            end
            
             
        end
    end
    
end

trans_traces = mean(trans_traces,2);


end
