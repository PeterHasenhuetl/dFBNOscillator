function [autocorr_all_sorted, example_flies] = get_fig_panel_pooled_xcorr(mp_data_table, ...
    cross_corr_table, x_pos, y_pos, curr_blue, curr_red, source_data_details)
% Plots autocorrelations of all dFBN GCaMP traces.
% Code written by Peter Hasenhuetl.

color_map_1 = get_default_two_color_map_dark(10000, 10000, 0.8, 0.8, curr_blue, curr_red);
color_map_2 = get_default_color_map(10000, 1, curr_blue, "dark");

% Some specs for plotting
sz_1 = 4;
ht_1 = 4;
xlm_lags = [-4, 4];
clim_cross = [-0.3, 0.3];
clim_auto = [-0.3, 0.3];
scale_fac = 10^3;

% Empty arrays for loop below
autocorr_dendrite_left = [];
autocorr_dendrite_right = [];
cross_corr_delta = [];
cross_corr_delta_over_flies = [];
corrcoef_delta = [];
delta_left = [];
delta_right = [];
autocorr_both_hemi = [];
autocorr_both_hemi_over_flies = [];
xcorr_per = [];
cross_corr_table_row = [];
corr_ID = [];
interhemi_coherence = [];
associated_fly_id = [];
associated_fly_name = [];

% Loops through flies (rows of cross_corr_table) and generates summary data for plotting
for loop_idx = 1:size(cross_corr_table,1)

    autocorr_dendrite_left = [autocorr_dendrite_left, cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_left];
    autocorr_dendrite_right = [autocorr_dendrite_right, cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_right];
    curr_ac_both_hemi = (cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_left +...
        cross_corr_table.ihc_stuff(loop_idx,1).auto_corr_dendrite_right)/2;
    autocorr_both_hemi = [autocorr_both_hemi, curr_ac_both_hemi];
    autocorr_both_hemi_over_flies = [autocorr_both_hemi_over_flies, mean(curr_ac_both_hemi,2)];
    cross_corr_delta = [cross_corr_delta, cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta];
    cross_corr_delta_over_flies = [cross_corr_delta_over_flies, mean(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2)];
    corrcoef_delta = [corrcoef_delta, cross_corr_table.ihc_stuff(loop_idx,1).corrcoef_delta];
    delta_left = [delta_left, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_left];
    delta_right = [delta_right, cross_corr_table.ihc_stuff(loop_idx,1).delta_power_right];
    xcorr_per = [xcorr_per, cross_corr_table.ihc_stuff(loop_idx,1).pr_delta];
    cross_corr_table_row = [cross_corr_table_row, loop_idx.*ones(1,size(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2))];
    associated_fly_id = [associated_fly_id, ...
        (cross_corr_table.fly_id(loop_idx,1)).*ones(1,size(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2))];
    curr_fly_name = [];
    for i = 1:size(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2)
        curr_fly_name = [curr_fly_name; cross_corr_table.fly_name(loop_idx,1)];
    end
    associated_fly_name = [associated_fly_name; curr_fly_name];
    corr_ID = [corr_ID, 1:size(cross_corr_table.ihc_stuff(loop_idx,1).cross_corr_delta,2)];
    interhemi_coherence = [interhemi_coherence, cross_corr_table.ihc_stuff(loop_idx,1).w_coh];

end

autocorr_all = [autocorr_dendrite_left, autocorr_dendrite_right];
SWA_power = mean([delta_left; delta_right],1);
time_lag = (cross_corr_table.ihc_stuff(1,1).time_lag)./14.56;

% Sorts cross-correlations according to the mean cross-corr value at around
% ± 300 ms lag (similar to just taking the corr-coef., but less sensitive
% to noisy traces)
idx_1 = (find(time_lag >= -0.3,1, 'first'));
idx_2 = (find(time_lag >= 0.3,1, 'first'))-1;
sort_array = mean(cross_corr_delta(idx_1:idx_2,:),1);
[~, b] = sort(sort_array,'ascend');
cross_corr_sorted = cross_corr_delta(:,b);
interhemi_coherence_sorted = interhemi_coherence(:,b);
cross_corr_table_row = cross_corr_table_row(1,b);
corr_ID = corr_ID(1,b);
associated_fly_id = associated_fly_id(1,b);
associated_fly_name = associated_fly_name(b,:);



example_flies = struct;
example_flies.fly_1.fly_name = "fly20240723_1DFB_SD3";
example_flies.fly_1.corr_ID = 1;

example_flies.fly_2.fly_name = "fly20240907_4DFB_B";
example_flies.fly_2.corr_ID = 5;

example_flies.fly_3.fly_name = "fly20210605_12DFB_C";
example_flies.fly_3.corr_ID = 2;

example_flies.fly_4.fly_name = "fly20240824_1DFB_B";
example_flies.fly_4.corr_ID = 1;

idx_1 = associated_fly_name == example_flies.fly_1.fly_name;
idx_2 = (corr_ID == example_flies.fly_1.corr_ID);
fi_1 = find(idx_1 == 1 & idx_2.' == 1);

idx_1 = associated_fly_name == example_flies.fly_2.fly_name;
idx_2 = (corr_ID == example_flies.fly_2.corr_ID);
fi_2 = find(idx_1 == 1 & idx_2.' == 1);

idx_1 = associated_fly_name == example_flies.fly_3.fly_name;
idx_2 = (corr_ID == example_flies.fly_3.corr_ID);
fi_3 = find(idx_1 == 1 & idx_2.' == 1);

idx_1 = associated_fly_name == example_flies.fly_4.fly_name;
idx_2 = (corr_ID == example_flies.fly_4.corr_ID);
fi_4 = find(idx_1 == 1 & idx_2.' == 1);



%% Plots the autocorrelations

auto_corr_cond = "all_flies";
if auto_corr_cond == "interhemi_only"
    [~, b] = sort([SWA_power, SWA_power],'descend');
    autocorr_all_sorted = autocorr_all(:,b);
    autocorr_all_sorted = (unique(autocorr_all_sorted','rows','stable'));
    autocorr_all_sorted = autocorr_all_sorted';
    ylm_xcorr = [0.5, size(autocorr_all_sorted,2)+0.5];
elseif auto_corr_cond == "all_flies"
    [autocorr_all_sorted, pwr_all_sorted] = get_auto_full(mp_data_table);
    ylm_xcorr = [0.5, size(autocorr_all_sorted,2)+0.5];
end

autocorr_pos = y_pos+ht_1+2.1;
min_max = [min(autocorr_all_sorted,[],"all"), max(autocorr_all_sorted,[],"all")];

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, autocorr_pos, sz_1, ht_1]);
imagesc(pnl_1, time_lag, 1:size(autocorr_all_sorted,2), autocorr_all_sorted') % Note the transpose
colormap(pnl_1, color_map_1)
clim(clim_auto)
get_default_ax(pnl_1, xlm_lags(1), xlm_lags(2), xlm_lags(1), xlm_lags(2), ...
    ylm_xcorr(1), ylm_xcorr(2), 0, ylm_xcorr(2),...
    1, 100, "linear", "linear", 'Lag (s)', 'ROI',...
    'k', 'k', 'k', 'k', sz_1, ht_1)
get_default_colorbar(pnl_1, x_pos, autocorr_pos+ht_1+0.1, sz_1, 0.2, ...
    clim_auto, min_max, 'r', 'northoutside')
get_default_annotation(x_pos+sz_1, autocorr_pos+0.25, [num2str(size(mp_data_table,1)/4), ' flies'], 'k', 'normal', "right")
get_default_annotation(x_pos+0.05, autocorr_pos+0.25, 'Auto-correlations', 'k', 'normal', "left")

power_cond = "with_power_spectrum";
if power_cond == "with_power_spectrum"
    sz_pwr = 1;
    clim_pwr = [0, 0.002].*scale_fac;
    xlm_pwr = [0, 2];
    ylm_pwr = [1, size(pwr_all_sorted,2)];
    x_dist_pwr = 0.3;
    frequency_axis_pwr = mp_data_table.frqu_trc(1,1).trace_img;
    min_max = [min(pwr_all_sorted.*scale_fac,[],"all"), max(pwr_all_sorted.*scale_fac,[],"all")];
    pnl_pwr = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist_pwr, ...
        autocorr_pos, sz_pwr, ht_1]);
    imagesc(pnl_pwr, frequency_axis_pwr, 1:size(pwr_all_sorted,2), (pwr_all_sorted.*scale_fac).') % Note the transpose
    hold on
    plot(pnl_pwr, [1, 1], ylm_pwr, 'k:', 'LineWidth', 0.5)
    colormap(pnl_pwr, color_map_2)
    clim(clim_pwr)
    get_default_ax(pnl_pwr, xlm_pwr(1), xlm_pwr(2), xlm_pwr(1), xlm_pwr(2), ...
        ylm_pwr(1), ylm_pwr(2), [], [],...
        1, [], "linear", "linear", 'Frequ. (Hz)', [],...
        'k', 'none', 'k', 'none', sz_pwr, ht_1)    
    get_default_colorbar(pnl_pwr, x_pos+sz_1+x_dist_pwr, autocorr_pos+ht_1+0.1, ...
        sz_pwr, 0.2, clim_pwr, min_max, {'(\DeltaF/F)^2 x10^{-3}'}, 'northoutside') 
end


%% plots the cross-correlations

crosscorr_pos = y_pos;
ylm_xcorr = [0.5, size(cross_corr_sorted,2)+0.5];
min_max = [min(cross_corr_sorted,[],"all"), max(cross_corr_sorted,[],"all")];
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, crosscorr_pos, sz_1, ht_1]);
imagesc(pnl_1, time_lag, 1:size(cross_corr_sorted,2), cross_corr_sorted') % Note the transpose
hold on
plot(pnl_1,[xlm_lags(2)-0.2,xlm_lags(2)],[fi_1, fi_1],'Color','k','LineWidth',2)
plot(pnl_1,[xlm_lags(2)-0.2,xlm_lags(2)],[fi_2, fi_2],'Color','k','LineWidth',2)
plot(pnl_1,[xlm_lags(2)-0.2,xlm_lags(2)],[fi_3, fi_3],'Color','k','LineWidth',2)
plot(pnl_1,[xlm_lags(2)-0.2,xlm_lags(2)],[fi_4, fi_4],'Color','k','LineWidth',2)
colormap(pnl_1, color_map_1)
clim(clim_cross)
get_default_ax(pnl_1, xlm_lags(1), xlm_lags(2), xlm_lags(1), xlm_lags(2), ...
    ylm_xcorr(1), ylm_xcorr(2), 0, ylm_xcorr(2),...
    1, 100, "linear", "linear", 'Lag (s)', 'Interhemispheric crosscorrelation',...
    'k', 'k', 'k', 'k', sz_1, ht_1)
get_default_colorbar(pnl_1, x_pos, crosscorr_pos+ht_1+0.1, sz_1, 0.2, ...
    clim_cross, min_max, 'r', 'northoutside')
get_default_annotation(x_pos+sz_1, crosscorr_pos+0.25, [num2str(loop_idx), ' flies'], 'k', 'normal', "right")
get_default_annotation(x_pos+0.05, crosscorr_pos+0.25, 'Cross-correlations', 'k', 'normal', "left")

coherence_cond = "with_coherence";
if coherence_cond == "with_coherence"
    sz_coherence = sz_pwr;
    clim_coherence = [0, 0.6];
    xlm_coherence = [0, 2];
    ylm_coherence = [1, size(interhemi_coherence_sorted,2)];
    x_dist_coherence = x_dist_pwr;
    frequency_axis_coherence = cross_corr_table.ihc_stuff(9,1).f_wcor;

    min_max = [min(interhemi_coherence_sorted,[],"all"), max(interhemi_coherence_sorted,[],"all")];
    pnl_2 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist_coherence, ...
        crosscorr_pos, sz_coherence, ht_1]);
    imagesc(pnl_2, frequency_axis_coherence, 1:size(interhemi_coherence_sorted,2), interhemi_coherence_sorted') % Note the transpose
    hold on
    plot(pnl_2, [1, 1], ylm_coherence, 'k:', 'LineWidth', 0.5)
    colormap(pnl_2, color_map_2)
    clim(clim_coherence)
    get_default_ax(pnl_2, xlm_coherence(1), xlm_coherence(2), xlm_coherence(1), xlm_coherence(2), ...
        ylm_coherence(1), ylm_coherence(2), [], [],...
        1, [], "linear", "linear", 'Frequ. (Hz)', [],...
        'k', 'none', 'k', 'none', sz_coherence, ht_1)
    get_default_colorbar(pnl_2, x_pos+sz_1+x_dist_coherence, crosscorr_pos+ht_1+0.1, ...
        sz_coherence, 0.2, clim_coherence, min_max, 'Coherence', 'northoutside')

end


%% Saves the source data

if isempty(source_data_details) == 0

    x_values = time_lag;
    autocorr_full = table;
    autocorr_full.x_values = x_values';
    autocorr_full.autocorr_all = autocorr_all_sorted;
    writetable(autocorr_full, [source_data_details.data_path, source_data_details.file_name1],'Sheet','autocorr_full')

    x_values = frequency_axis_pwr;
    power_full = table;
    power_full.x_values = x_values;
    power_full.power_all = pwr_all_sorted;
    writetable(power_full, [source_data_details.data_path, source_data_details.file_name1],'Sheet','power_full')

    x_values = time_lag;
    crosscorr_full = table;
    crosscorr_full.x_values = x_values';
    crosscorr_full.crosscorr_all = cross_corr_sorted;
    writetable(crosscorr_full, [source_data_details.data_path, source_data_details.file_name1],'Sheet','crosscorr_full')

    x_values = frequency_axis_coherence;
    interhemi_coherence = table;
    interhemi_coherence.x_values = x_values;
    interhemi_coherence.coherence_all = interhemi_coherence_sorted;
    writetable(interhemi_coherence, [source_data_details.data_path, source_data_details.file_name1],'Sheet','coherence_full')

end

end


function [AUTO_full, SP_full] = get_auto_full(mp_data_table)

curr_tbl = mp_data_table;
n_flies = round(size(curr_tbl,1)/4);
num_array = 1:10000;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
SP_full = [];
IDP_full = [];
AUTO_full = [];
cum_size = [];
for loop_idx = 1:n_flies
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);
    [n_pwr1, n_delta_power1, n_auto, curr_size] = get_pooled_measures(curr_img_tbl);
    if isempty(n_pwr1) == 0
        SP_full = [SP_full, n_pwr1];
        IDP_full = [IDP_full, n_delta_power1];
        AUTO_full = [AUTO_full, n_auto];
        cum_size = [cum_size, curr_size];
    end
end

n_ROIs = sum(cum_size);

[~, b] = sort(IDP_full,'descend');

AUTO_full = AUTO_full(:,b);
SP_full = SP_full(:,b);
end

function [n_pwr, n_delta_power, n_auto, curr_size] = get_pooled_measures(data_tbl)
   
curr_table_pwr = data_tbl.power_spectrum;
curr_delta_table = data_tbl.int_delta_power;
curr_table_auto = data_tbl.auto_corr_delta;
excl_vec = data_tbl.excl_idx;

curr_pwr = [];
curr_delta_power = [];
curr_size = [];
curr_auto = [];
for loop_idx1 = 1:size(data_tbl,1)
    loop_cond = 1;
    if excl_vec(loop_idx1,3) ~= 0 && excl_vec(loop_idx1,4) ~= 0
        loop_cond = 0;
    elseif excl_vec(loop_idx1,3) ~= 0 && excl_vec(loop_idx1,4) == 0
        idx_1 = 4;
    elseif excl_vec(loop_idx1,4) ~= 0 && excl_vec(loop_idx1,3) == 0
        idx_1 = 3;
    elseif excl_vec(loop_idx1,3) == 0 && excl_vec(loop_idx1,4) == 0
        idx_1 = [3, 4];
    end
    if loop_cond == 1
        curr_pwr = [curr_pwr, curr_table_pwr(loop_idx1,1).trace_img(:,idx_1)];
        curr_delta_power = [curr_delta_power, curr_delta_table(loop_idx1,idx_1)];
        curr_auto = [curr_auto, curr_table_auto(loop_idx1,1).trace_img(:,idx_1)];
        curr_size = [curr_size, size(curr_table_pwr(loop_idx1,1).trace_img(:,idx_1),2)];
    end
end
curr_size = sum(curr_size);

n_pwr = [];
n_delta_power = [];
n_auto = [];
n = 1;
for loop_idx2 = 1:size(curr_pwr,2)
    if sum(curr_pwr(:,loop_idx2)) > 0
        
        n_pwr(:,n) = (curr_pwr(:,loop_idx2));
        n_delta_power(1,n) = curr_delta_power(loop_idx2);
        n_auto(:,n) = (curr_auto(:,loop_idx2));
        n = n+1;

    end
    
end

end





