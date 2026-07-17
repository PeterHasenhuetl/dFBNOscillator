% Analysis of dFBN GCaMP dynamics under various levels of sleep pressure
% and at different zeitgeber times.
% Code written by Peter Hasenhuetl.

clear all
tic

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'GCaMP_dendritic_example_recording.xlsx';

source_data_SNAP_details.data_path = []; %Add path as character array
source_data_SNAP_details.file_name = 'SNAP_sleep_deprivation_GCaMP_data.xlsx';

source_data_vortex_details.data_path = []; %Add path as character array
source_data_vortex_details.file_name = 'vortex_sleep_deprivation_GCaMP_data.xlsx';

source_data_interhemi_details.data_path = []; %Add path as character array
source_data_interhemi_details.file_name = 'GCaMP_interhemispheric_coordination_data.xlsx';

source_data_axon_vs_dendrite_details.data_path = []; %Add path as character array
source_data_axon_vs_dendrite_details.file_name = 'axon_vs_dendrite_data.xlsx';

source_data_pooled_xcorr_and_SD_details.data_path = []; %Add path as character array
source_data_pooled_xcorr_and_SD_details.file_name1 = 'pooled_xcorr.xlsx';
source_data_pooled_xcorr_and_SD_details.file_name2 = 'example_moving_xcorr.xlsx';

source_data_survival_imaging_flies.data_path = []; %Add path as character array
source_data_survival_imaging_flies.file_name = 'survival_imaging_flies.xlsx';
%
% goes into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
mp_data_table_SNAP = [];
cross_corr_table_SNAP = [];
ad_corr_table_SNAP = [];
curr_filter = get_mp_osc_filt(14.56, 0.005, 0.1, 0.1, 1.5);
do_manual_exlusions = "yes"
if do_manual_exlusions == "yes"
    % Excludes some traces due to motion artifacts (identified blinded to experimental condition).
    manual_exclusion_table = get_dFB_manual_exclusion_table;
else
    manual_exclusion_table = table;
    manual_exclusion_table.fly_id = NaN;
    manual_exclusion_table.plane_id = NaN;
    manual_exclusion_table.ROI_id = NaN;
end
quality_check_table = table;
curr_quality_check_table = table;
% loops through individual flies
for fly_idx = 1:n_flies
    
     cd([]) %Add path as character array
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        loaded_fly_original_number = char(fly_input(fly_idx).name);
        % loops through extracted planes
        [ind_data_table, cross_corr_table_row, ad_corr_table_row] = ...
            get_individual_multiplane_fly(loaded_fly, fly_idx, loaded_fly_original_number, curr_filter, manual_exclusion_table);
        mp_data_table_SNAP = [mp_data_table_SNAP; ind_data_table];
        if istable(cross_corr_table_row) == 1
            cross_corr_table_SNAP = [cross_corr_table_SNAP; cross_corr_table_row];
        end
        
        if istable(ad_corr_table_row) == 1
            ad_corr_table_SNAP = [ad_corr_table_SNAP; ad_corr_table_row];
        end
        
        curr_name = loaded_fly.fly_details.flyname;
        fly_names_SNAP(fly_idx,1:size(curr_name,2)) = curr_name;
        
        curr_quality_check_table.fly_name = string(curr_name);
        curr_quality_check_table.ZT = ind_data_table.zeitgeber_time(1);
        curr_quality_check_table.SD = ind_data_table.sleep_deprivation(1);
        curr_quality_check_table.mon_pos = loaded_fly.fly_details.monitor_position;
        
        quality_check_table = [quality_check_table; curr_quality_check_table];

end

cd([]) %Add path as character array

R2 = get_axo_dendritic_pooled_R2(ad_corr_table_SNAP);

writetable(quality_check_table,[source_data_details.data_path, 'quality_check_table_SNAP2026.xlsx']);

toc

%% Plots Fig. 1 example trace and power spectrum full dataset
color = get_color;

close all
figure('Name','full dataset','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 12],'Resize','off')
inset_lngth = round(14.56*100);
start_idx = 7350;
end_idx = start_idx+inset_lngth;
get_fig_panel_intro_Fig1(mp_data_table, 1, 1, color, start_idx, end_idx, source_data_details)

cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'GCaMP_example_recording_fig.pdf')
cd([]) %Add path as character array

%% Plots sleep deprivation data

color = get_color;
close all
figure('Name','SD Fig','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 12],'Resize','off')
sd_color = color.sleep_deprived;
ctrl_color = color.rested;
idx_c = 8000;
get_fig_panel_sleep_depr_img_traces(mp_data_table_SNAP, -1, 7.2, 1, [0,0,0], [0,0,0],...
    idx_c, idx_c+900, 13100, 14000, source_data_SNAP_details)
%
get_fig_panel_GCaMP_transient_SD(mp_data_table_SNAP, 7.1, 7.2, 1, 1, ...
    sd_color, ctrl_color, source_data_SNAP_details)
get_fig_panel_power_spectrum_SD(mp_data_table_SNAP, 0.8, 4.5, 1, ...
    sd_color, ctrl_color, source_data_SNAP_details)
get_fig_panel_img_feat_vs_ZT_SD(mp_data_table_SNAP, 6.75, 4.5, sd_color, ...
    ctrl_color, 'pwr', source_data_SNAP_details)
get_fig_panel_dFB_tuning_curve_SD(mp_data_table_SNAP, 3.75, 4.5, source_data_SNAP_details)

%
cd(source_data_SNAP_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'sleep_deprivation_GCaMP_fig.pdf')
cd([]) %Add path as character array

%% Plots transients tuning curve

close all
figure('Name','Transients tuning curve','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 12],'Resize','off')
get_fig_panel_dFB_tuning_curve_transients(mp_data_table_SNAP, 2, 2, ...
    source_data_SNAP_details)

cd(source_data_SNAP_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'GCaMP_transient_tuning_curve.pdf')
cd([]) %Add path as character array

%% Plots interhemispheric coordination

color = get_color;
plotting_color_2 = color.medium_gray;
plotting_color_1 = [0, 0, 0];
curr_blue = (color.navy)*255;
curr_red = [200, 85, 97];
two_color_map = get_default_two_color_map(1000, 5000, (curr_blue)./255, (curr_red)./255);

close all
figure('Name','GCaMP_interhemispheric_coordination_fig','Color','white',...
    'Units','centimeters','Position',[10 12 18 10],'Resize','off')

start_trace = 17380;
length_trace = 1456;
load([]) %Add path as character array
[dF, dF_contra, table_idx] = get_fig_panel_interhemisphere_images_example(mp_data_table);
get_fig_panel_interhemisphere_traces_example(mp_data_table_SNAP, ...
    0.125, 2.9, table_idx, table_idx, start_trace, start_trace+length_trace,plotting_color_1, plotting_color_2,...
    dF, dF_contra, two_color_map, source_data_interhemi_details)

contra_mean = get_fig_panel_SWA_power_vs_mutual_inhibition(cross_corr_table, 11, 3, color, ...
    "over flies", 0, source_data_interhemi_details);
corrcoef_delta = get_fig_panel_SWA_power_vs_xcorr(cross_corr_table, 14, 3, color, ...
    "over flies", source_data_interhemi_details);

cd(source_data_interhemi_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'GCaMP_interhemispheric_coordination_fig.pdf')
cd([]) %Add path as character array

%% Extended Data Figs

selected_fly = 14;
axon_idx = 1;
cv_idx = 1;
close all
figure('Name','axon vs dendrite fig','Color','white',...
    'Units','centimeters','Position',[10 12 18 12],'Resize','off')
ylm_1 = [-0.14, 1.2];
ylm_2 = [-0.1, 1];
get_fig_panel_axon_vs_dendrites_traces(ad_corr_table_SNAP, ...
    0.5, 2, cv_idx, selected_fly, 1, 3900, 2000, 3000, ylm_1, ylm_2, source_data_axon_vs_dendrite_details)

cd(source_data_axon_vs_dendrite_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'axon_vs_dendrite_fig.pdf')
cd([]) %Add path as character array

%% Plots pooled crosscorrelations

color = get_color;
sd_color = color.sleep_deprived;
ctrl_color = color.rested;
neg_blue = color.navy;
pos_red = ([200, 85, 97])./255;

close all
figure('Name','axon vs dendrite fig','Color','white',...
    'Units','centimeters','Position',[10 12 18 18+3.75+5],'Resize','off')
[autocorr_all, example_flies] = get_fig_panel_pooled_xcorr(mp_data_table, ...
    cross_corr_table, 0.9, 10.5, neg_blue, pos_red, source_data_pooled_xcorr_and_SD_details);

img_lim = [-2.65, 4.5];


x_pos = 9.25;
curr_fly_name = example_flies.fly_1.fly_name;
dend_idx = example_flies.fly_1.corr_ID;
r_lim = [-0.7, 0.1];
lim_1 = 6000;
get_fig_panel_moving_interhemi_cross_correlation(cross_corr_table, curr_fly_name, dend_idx, x_pos, 18.25+(2*2)-5, ...
    r_lim, img_lim, lim_1, "without_x_axis", "with_hemi_annotation", neg_blue, pos_red,...
    source_data_pooled_xcorr_and_SD_details)


curr_fly_name = example_flies.fly_2.fly_name;
dend_idx = example_flies.fly_2.corr_ID;
r_lim = [-0.7, 0.3];
lim_1 = 10000;
get_fig_panel_moving_interhemi_cross_correlation(cross_corr_table, curr_fly_name, dend_idx, x_pos, 14.75+(2*1.5)-5, ...
    r_lim, img_lim, lim_1, "without_x_axis", "without_hemi_annotation", neg_blue, pos_red,...
    source_data_pooled_xcorr_and_SD_details)


curr_fly_name = example_flies.fly_3.fly_name;
dend_idx = example_flies.fly_3.corr_ID;
r_lim = [-0.5, 0.4];
lim_1 = 3000;
get_fig_panel_moving_interhemi_cross_correlation(cross_corr_table, curr_fly_name, dend_idx, x_pos, 11.25+(2*1)-5, ...
    r_lim, img_lim, lim_1, "without_x_axis", "without_hemi_annotation", neg_blue, pos_red,...
    source_data_pooled_xcorr_and_SD_details)


curr_fly_name = example_flies.fly_4.fly_name;
dend_idx = example_flies.fly_4.corr_ID;
r_lim = [-0.5, 0.5];
lim_1 = 10000;
get_fig_panel_moving_interhemi_cross_correlation(cross_corr_table, curr_fly_name, dend_idx, x_pos, 7.5+(2*0.5)-5, ...
    r_lim, img_lim, lim_1, "with_x_axis", "without_hemi_annotation", neg_blue, pos_red,...
    source_data_pooled_xcorr_and_SD_details)


get_fig_panel_img_feat_vs_ZT_SD(mp_data_table_SNAP, 1, 7.6, sd_color, ctrl_color,...
    "transients", source_data_SNAP_details)
get_fig_panel_power_spectrum_SD(mp_data_table_vortex, 3.75, 7.6, 1, ...
    sd_color, ctrl_color, source_data_vortex_details)
get_fig_panel_GCaMP_transient_SD(mp_data_table_vortex, 6.25, 7.5, 1, 1, ...
    sd_color, ctrl_color, source_data_vortex_details)

if isempty(source_data_pooled_xcorr_and_SD_details) == 0
    cd(source_data_pooled_xcorr_and_SD_details.data_path)
    set(gcf,'renderer','Painters')
    saveas(gcf,'pooled_xcorr_and_SD_details.pdf')
    cd([]) %Add path as character array
end

%% Plots the survival of imaging flies after surgery

close all
figure('Name','survival_plot','Color','white',...
    'Units','centimeters','Position',[10 12 18 12],'Resize','off')
[C, W] = get_fig_panel_survival;
cd(source_data_survival_imaging_flies.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'survival_plot.pdf')
cd([]) %Add path as character array

%%

function [ind_data_table, cross_corr_table_row, ad_corr_table_row] = ...
    get_individual_multiplane_fly(input_fly, fly_idx, loaded_fly_original_number, curr_filter, manual_exclusion_table)

fields_img_sessions = fieldnames(input_fly.img_traces);
    
for plane_idx = 1:length(fields_img_sessions)
    plane_name = fields_img_sessions{plane_idx}; 
    input_data = input_fly.img_traces.(string(plane_name));

    curr_manual_ex_tbl = manual_exclusion_table(manual_exclusion_table.fly_name == string(input_fly.fly_details.flyname),:);
    
    if isempty(curr_manual_ex_tbl) == 0
        curr_manual_ex_tbl = curr_manual_ex_tbl(curr_manual_ex_tbl.plane_id == plane_idx,:);
        if isempty(curr_manual_ex_tbl) == 0
            input_data = get_dFB_manual_exclusions(input_data, curr_manual_ex_tbl.ROI_id);
        end
    end

    input_fly.fly_details.loaded_fly_original_number = loaded_fly_original_number;
    output_img = get_mp_img_measures(input_data, curr_filter);
    ind_data_table(plane_idx,:) = get_mp_sd_data_table_row(output_img, input_fly.fly_details, fly_idx, plane_idx);
end           

ind_data_table = get_ROI_selection(ind_data_table);

inter_hemispheric_correlation = get_mp_inter_hemispheric_correlation(ind_data_table, curr_filter);

if isstruct(inter_hemispheric_correlation) == 1
    cross_corr_table_row = get_mp_inter_corr_table_row(inter_hemispheric_correlation, input_fly.fly_details, fly_idx);
else
    cross_corr_table_row = NaN;
end

axo_dendritic_correlation = get_axo_dendritic_correlation(ind_data_table);

if isstruct(axo_dendritic_correlation) == 1
    ad_corr_table_row = get_mp_ad_corr_table_row(axo_dendritic_correlation, input_fly.fly_details, fly_idx);
else
    ad_corr_table_row = NaN;
end

end

%%

function ind_data_table = get_ROI_selection(curr_table)

for hemi_idx = 3:4  
    curr_dendrite = NaN(19500,4);
    % loops through planes and extracts data for ROIs
    for plane_idx = 1:4
        curr_dendrite(:,plane_idx) = curr_table.img_trace(plane_idx,1).trace_img(:,hemi_idx);
    end
        inter_ROI_corr = corrcoef(curr_dendrite);
        [high_corr_idx, ~] = find(inter_ROI_corr >= 0.85 & inter_ROI_corr < 1);
        high_corr_idx = unique(high_corr_idx,'stable');
        mean_img = mean(curr_dendrite(:,high_corr_idx),1);
        
        leave_idx = find(mean_img == max(mean_img));
        high_corr_idx(leave_idx) = [];
        excl_ROI = high_corr_idx;
        curr_table.excl_idx(excl_ROI, hemi_idx) = 2;
end

[row_idx, col_idx] = (find(curr_table.excl_idx == 2));

for loop_idx = 1:length(row_idx)
    curr_table.trans_ids(row_idx(loop_idx),1).transient_onsets(:,col_idx(loop_idx)) = NaN;
end

ind_data_table = curr_table;

end

%%

function dm_row = get_mp_sd_data_table_row(dm_input, fly_details, fly_idx, plane_idx)

dm_behav_row = table;
dm_behav_row.zeitgeber_time = hours(fly_details.start_time-fly_details.beginning_of_entire_experiment);
dm_behav_row.sleep_deprivation = string(fly_details.sleep_deprivation);
dm_behav_row.fly_id = fly_idx;
dm_behav_row.plane_id = plane_idx;
dm_behav_row.fly_name = string(fly_details.flyname);
dm_behav_row.original_id = string(fly_details.loaded_fly_original_number);
dm_behav_row.details = fly_details;

dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;

dm_img_row = table;

dF_trc = struct;
dF_trc.trace_img = dm_input.dF;
dm_img_row.img_trace = dF_trc;

trc_raw = struct;
trc_raw.trace_img = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;

trc_infra = struct;
trc_infra.trace_img = dm_input.trace_infra;
dm_img_row.trace_infra = trc_infra;

trc_delta = struct;
trc_delta.trace_img = dm_input.trace_delta;
dm_img_row.trace_delta = trc_delta;

transient_vec = struct;
transient_vec.trace_img = dm_input.transient_vec;
dm_img_row.transient_vec = transient_vec;

a_c_delta = struct;
a_c_delta.trace_img = dm_input.auto_corr_delta';
dm_img_row.auto_corr_delta = a_c_delta;

pwr = struct;
pwr.trace_img = dm_input.amp;
dm_img_row.power_spectrum = pwr;

frqu_trc = struct;
frqu_trc.trace_img = dm_input.f;
dm_img_row.frqu_trc = frqu_trc;

trans_trace_raw = struct;
trans_trace_raw.trace_img = dm_input.trans_trace_raw;
dm_img_row.trans_trace_raw = trans_trace_raw;

trans_trace = struct;
trans_trace.trace_img = dm_input.trans_trace;
dm_img_row.trans_trace = trans_trace;

trans_ids = struct;
trans_ids.transient_onsets = dm_input.transient_onsets;
dm_img_row.trans_ids = trans_ids;

dm_img_row.int_delta_power = dm_input.int_delta_power;
dm_img_row.int_infra_slow_power = dm_input.int_infra_slow_power;
dm_img_row.ITI_CV = dm_input.ITI_CV;
dm_img_row.ITI_mean = dm_input.ITI_mean;
dm_img_row.transient_mean = dm_input.transient_mean;

dm_img_row.delta_period = dm_input.pr_delta;
dm_img_row.delta_corr = dm_input.amp_delta;
dm_img_row.background_noise = dm_input.background_noise;
dm_img_row.SNR = dm_input.SNR;

% Concatenates them
dm_row = [dm_behav_row, dm_excl_row, dm_img_row];

end

%%

function interhemi_c = get_mp_inter_hemispheric_correlation(ind_data_table, curr_filter)
% analyses interhemispheric coordination and saves data in struct.

% empty arrays for loop below
dendrite_1 = [];
dendrite_2 = [];
dendrite1_raw = [];
dendrite2_raw = [];
dendrite1_transient_vec = [];
dendrite2_transient_vec = [];
pln_coord1 = [];
pln_coord2 = [];
CA = [];
IA = [];
delta_pwr1 = [];
delta_pwr2 = [];
interhemi_scatter_cond = 0;

% checks if both hemispheres were recorded
if any(ind_data_table.excl_idx(:,3) == 0) == 1 && any(ind_data_table.excl_idx(:,4) == 0) == 1
    
    % loops through planes and extracts data for ROIs of left and right dendrites
    for plane_idx = 1:4 
        if ind_data_table.excl_idx(plane_idx,3) == 0
            dendrite_1 = [dendrite_1, ind_data_table.img_trace(plane_idx,1).trace_img(:,3)];
            dendrite1_raw = [dendrite1_raw, ind_data_table.img_raw(plane_idx,1).trace_img(:,3)];
            dendrite1_transient_vec = [dendrite1_transient_vec, ind_data_table.transient_vec(plane_idx,1).trace_img(:,3)];
            pln_coord1 = [pln_coord1, plane_idx];
            delta_pwr1 = [delta_pwr1, ind_data_table.int_delta_power(plane_idx,3)];      
        end
    
        if ind_data_table.excl_idx(plane_idx,4) == 0
            dendrite_2 = [dendrite_2, ind_data_table.img_trace(plane_idx,1).trace_img(:,4)];
            dendrite2_raw = [dendrite2_raw, ind_data_table.img_raw(plane_idx,1).trace_img(:,4)];
            dendrite2_transient_vec = [dendrite2_transient_vec, ind_data_table.transient_vec(plane_idx,1).trace_img(:,4)];
            pln_coord2 = [pln_coord2, plane_idx];
            delta_pwr2 = [delta_pwr2, ind_data_table.int_delta_power(plane_idx,4)];
        end
    end

    % loops through both dendrites and gets summary measures of
    % interhemispheric coordination
    summary_idx = 1;
    % loops through ROIs of left dendrite
    for left_idx = 1:size(dendrite_1,2)
        
        % loops through ROIs of right dendrite
        for right_idx = 1:size(dendrite_2,2) 

            % Computes Magnitude-squared coherence between hemispheres.
            % Coherence window was chosen to satisfy (14.56/0.15)*5: 
            % --> Sampling rate divided by frequency range of interest
            % (here, 0.15 to cover also slower frequencies of coherence), 
            % times 5 to cover several periods of the rhythm of interest
            coherence_window = 500; 
            coherence_overlap = 250;
            [interhemi_c.w_coh(:,summary_idx), interhemi_c.f_wcor] = ...
                mscohere(zscore(dendrite_1(:,left_idx)), zscore(dendrite_2(:,right_idx)),...
                coherence_window, coherence_overlap, [], 14.56, 'onesided');


            % gets cross correlations
            [interhemi_c.pr_delta(summary_idx), interhemi_c.amp_delta(summary_idx), ...
                interhemi_c.cross_corr_delta(:,summary_idx), ...
                interhemi_c.corrcoef_delta(summary_idx), interhemi_c.corrcoef_infra(summary_idx), ...
                interhemi_c.auto_corr_dendrite_left(:,summary_idx), interhemi_c.auto_corr_dendrite_right(:,summary_idx), time_lag] = ...
                get_mp_cross_corr(dendrite_1(:,left_idx), dendrite_2(:,right_idx), curr_filter);
            % gets plane coordinates of dendrites
            interhemi_c.cross_coord(summary_idx,1:2) = [pln_coord1(left_idx), pln_coord2(right_idx)];
            % gets aligned GCaMP transients raw
            [interhemi_c.al_trc_ipsi, interhemi_c.al_trc_contra] = get_aligned_transients(dendrite_1(:,left_idx), dendrite_2(:,right_idx),...
                dendrite1_transient_vec(:,left_idx), dendrite2_transient_vec(:,right_idx));
            % gets aligned GCaMP transients dF/F
            [interhemi_c.ipsi_trace_left(:,summary_idx), interhemi_c.contra_trace_right(:,summary_idx), ...
                interhemi_c.ipsi_trace_right(:,summary_idx), ...
                interhemi_c.contra_trace_left(:,summary_idx),interhemi_c.n_transients_left(summary_idx), ...
                interhemi_c.n_transients_right(summary_idx), ...
                ipsi_all, contra_all] = get_aligned_dF_transients(dendrite1_raw(:,left_idx), dendrite2_raw(:,right_idx),...
                dendrite1_transient_vec(:,left_idx), dendrite2_transient_vec(:,right_idx));
            % all individual transients as one large array
            CA = [CA, contra_all];
            IA = [IA, ipsi_all];
             
            interhemi_c.delta_power_left(1,summary_idx) = delta_pwr1(1,left_idx);
            interhemi_c.delta_power_right(1,summary_idx) = delta_pwr2(1,right_idx);
            interhemi_c.delta_power_for_corr_corr(1,summary_idx) = mean([delta_pwr1(1,left_idx), delta_pwr2(1,right_idx)]);           
            interhemi_c.dendrite_1(:,summary_idx) = dendrite_1(:,left_idx);
            interhemi_c.dendrite_2(:,summary_idx) = dendrite_2(:,right_idx);
            interhemi_c.time_lag = time_lag;
            
            if interhemi_scatter_cond == 1
                [interhemi_c.scatter_hist(:,:,summary_idx), a_1, a_2] = ...
                    histcounts2(zscore(diff(dendrite_1(:,left_idx))), zscore(diff(dendrite_2(:,right_idx))),...
                    'NumBins', 81, 'XBinLimits', [-3, 3], 'YBinLimits', [-3, 3], 'Normalization', 'probability');
                interhemi_c.scatter_hist_edges_1 = a_1;
                interhemi_c.scatter_hist_edges_2 = a_2;
            end
            
            summary_idx = summary_idx+1;
        end    
    end

    interhemi_c.CA = CA;
    interhemi_c.IA = IA;
    
else
    interhemi_c = NaN;
end





end

%%

function [al_trc_ipsi, al_trc_contra] = get_aligned_transients(dendrite1_img, ...
    dendrite2_img, dendrite1_transient_vec, dendrite2_transient_vec)

trc_1 = [];
trc_2 = [];
p_wind = 40;
n_wind = 15;
aligned_trans_idx = 1;
d_burst = [dendrite1_transient_vec, dendrite2_transient_vec];

for hemi_idx = 1:2
    curr_transient_vec = d_burst(:,hemi_idx);
    transient_onsets = find(curr_transient_vec(:,1) > 0);   
    % Only considers the transient onsets that allow for the full        
    % length of transients to be used (i.e., excludes transients at        
    % the beginning and end of the imaging trace).        
    transient_onsets((transient_onsets+p_wind) > length(dendrite1_img) | ...
        (transient_onsets-n_wind) < 1) =  [];
    
    if hemi_idx == 1
        img_1 = dendrite1_img;
        img_2 = dendrite2_img;
    elseif hemi_idx == 2
        img_1 = dendrite2_img;
        img_2 = dendrite1_img;
    end
    
    for ind_trans_idx = 1:length(transient_onsets)
        trc_1(:,aligned_trans_idx) = (img_2(transient_onsets(ind_trans_idx)-n_wind:transient_onsets(ind_trans_idx)+p_wind,1));
        trc_2(:,aligned_trans_idx) = (img_1(transient_onsets(ind_trans_idx)-n_wind:transient_onsets(ind_trans_idx)+p_wind,1));
        aligned_trans_idx = aligned_trans_idx + 1;
    end
end

al_trc_contra = trc_1;
al_trc_ipsi = trc_2;

end

%%

function [ipsi_trace_left, contra_trace_right, ipsi_trace_right, contra_trace_left,...
    n_transients_left, n_transients_right, ipsi_all, contra_all] = ...
    get_aligned_dF_transients(curr_trace_left_raw, curr_trace_right_raw, curr_transient_vec_left, curr_transient_vec_right)
 
p_wind = 40;
n_wind = 10;
F0_wind = round(n_wind);

transient_onsets = find(curr_transient_vec_left > 0);
% Only considers the transient onsets that allow for the full            
% length of transients to be used (i.e., excludes transients at            
% the beginning and end of the imaging trace).          
transient_onsets((transient_onsets+p_wind) > length(curr_trace_left_raw) | ...
    (transient_onsets-n_wind) < 1) =  [];

ipsi_trace_left = NaN(n_wind+p_wind+1,length(transient_onsets));
contra_trace_right = NaN(n_wind+p_wind+1,length(transient_onsets));
for loop_idx1 = 1:length(transient_onsets)
    i_trace_left = (curr_trace_left_raw(transient_onsets(loop_idx1)-n_wind:transient_onsets(loop_idx1)+p_wind,1));
    ipsi_trace_left(:,loop_idx1) = (i_trace_left-mean(i_trace_left(1:F0_wind,1)))/mean(i_trace_left(1:F0_wind,1));
    c_trace_right = (curr_trace_right_raw(transient_onsets(loop_idx1)-n_wind:transient_onsets(loop_idx1)+p_wind,1));
    contra_trace_right(:,loop_idx1) = (c_trace_right-mean(c_trace_right(1:F0_wind,1)))/mean(c_trace_right(1:F0_wind,1));
end

n_transients_left = loop_idx1;
  
transient_onsets = find(curr_transient_vec_right > 0);
% Only considers the transient onsets that allow for the full            
% length of transients to be used (i.e., excludes transients at            
% the beginning and end of the imaging trace). 
transient_onsets((transient_onsets+p_wind) > length(curr_trace_right_raw) | ...
    (transient_onsets-n_wind) < 1) =  [];

ipsi_trace_right = NaN(n_wind+p_wind+1,length(transient_onsets));
contra_trace_left = NaN(n_wind+p_wind+1,length(transient_onsets));
for loop_idx2 = 1:length(transient_onsets)
    i_trace_right = (curr_trace_right_raw(transient_onsets(loop_idx2)-n_wind:transient_onsets(loop_idx2)+p_wind,1));
    ipsi_trace_right(:,loop_idx2) = (i_trace_right-mean(i_trace_right(1:F0_wind,1)))/mean(i_trace_right(1:F0_wind,1));
    c_trace_left = (curr_trace_left_raw(transient_onsets(loop_idx2)-n_wind:transient_onsets(loop_idx2)+p_wind,1));
    contra_trace_left(:,loop_idx2) = (c_trace_left-mean(c_trace_left(1:F0_wind,1)))/mean(c_trace_left(1:F0_wind,1));
end

n_transients_right = loop_idx2;

ipsi_all = [ipsi_trace_left, ipsi_trace_right];
contra_all = [contra_trace_left, contra_trace_right];
ipsi_trace_left = mean(ipsi_trace_left,2);
contra_trace_right = mean(contra_trace_right,2);
ipsi_trace_right = mean(ipsi_trace_right,2);
contra_trace_left = mean(contra_trace_left,2);

end

%%

function interhemi_dm_row = get_mp_inter_corr_table_row(ihc, fly_details, fly_idx)

dm_behav_row = table;
dm_behav_row.zeitgeber_time = hours(fly_details.start_time-fly_details.beginning_of_entire_experiment);
dm_behav_row.sleep_deprivation = string(fly_details.sleep_deprivation);
dm_behav_row.fly_id = fly_idx;
dm_behav_row.fly_name = string(fly_details.flyname);
dm_behav_row.original_id = string(fly_details.loaded_fly_original_number);


dm_img_row = table;

ihc_stuff = ihc;
dm_img_row.ihc_stuff = ihc_stuff;

cross_corr_coordinates = struct;
cross_corr_coordinates.cross_coord_values = ihc.cross_coord;
dm_img_row.cross_corr_delta = cross_corr_coordinates;

c_c_delta = struct;
c_c_delta.cross_corr_delta_trace = ihc.cross_corr_delta;
dm_img_row.cross_corr_delta = c_c_delta;

al_trc_ipsi = struct;
al_trc_ipsi.ipsi_trace = ihc.al_trc_ipsi;
dm_img_row.al_trc_ipsi = al_trc_ipsi;

al_trc_contra = struct;
al_trc_contra.contra_trace = ihc.al_trc_contra;
dm_img_row.al_trc_contra = al_trc_contra;

tbl_inpt = NaN(1,16);
in_dt = ihc.pr_delta;
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.delta_period = tbl_inpt;

tbl_inpt = NaN(1,16);
in_dt = ihc.amp_delta;
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.delta_corr = tbl_inpt;

tbl_inpt = NaN(1,16);
in_dt = ihc.corrcoef_delta;
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.corrcoef_delta = tbl_inpt;

tbl_inpt = NaN(1,16);
in_dt = ihc.corrcoef_infra;
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.corrcoef_infra = tbl_inpt;

tbl_inpt = NaN(1,16);
in_dt = ihc.cross_coord(:,1)';
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.coordinates_dendrite1 = tbl_inpt;

tbl_inpt = NaN(1,16);
in_dt = ihc.cross_coord(:,2)';
tbl_inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.coordinates_dendrite2 = tbl_inpt;

% concatenates them

interhemi_dm_row = [dm_behav_row, dm_img_row];

end

%%

function ad_dm_row = get_mp_ad_corr_table_row(ad_corr, fly_details, fly_idx)

dm_behav_row = table;
dm_behav_row.zeitgeber_time = hours(fly_details.start_time-fly_details.beginning_of_entire_experiment);
dm_behav_row.sleep_deprivation = string(fly_details.sleep_deprivation);
dm_behav_row.fly_id = fly_idx;

dm_img_row = table;
dm_img_row.axo_dendritic_corr = ad_corr;

ad_dm_row = [dm_behav_row, dm_img_row];

end

%%

function ad_corr = get_axo_dendritic_correlation(ind_data_table)

axon_signal = [];
dendrite_1 = [];
dendrite_2 = [];
dendrite1_transient_vec = [];
dendrite2_transient_vec = [];
pln_coord1 = [];
pln_coord2 = [];

for loop_idx1 = 1:4
    if ind_data_table.excl_idx(loop_idx1,3) == 0
        dendrite_1 = [dendrite_1, ind_data_table.img_trace(loop_idx1,1).trace_img(:,3)];
        dendrite1_transient_vec = [dendrite1_transient_vec, ind_data_table.transient_vec(loop_idx1,1).trace_img(:,3)];
        pln_coord1 = [pln_coord1, loop_idx1];
    end
    
    if ind_data_table.excl_idx(loop_idx1,4) == 0
        dendrite_2 = [dendrite_2, ind_data_table.img_trace(loop_idx1,1).trace_img(:,4)];
        dendrite2_transient_vec = [dendrite2_transient_vec, ind_data_table.transient_vec(loop_idx1,1).trace_img(:,4)];
        pln_coord2 = [pln_coord2, loop_idx1];
    end
    
    if ind_data_table.excl_idx(loop_idx1,1) == 0 && ind_data_table.excl_idx(loop_idx1,2) == 0
        ind_trace = [ind_data_table.img_trace(loop_idx1,1).trace_img(:,1), ind_data_table.img_trace(loop_idx1,1).trace_img(:,2)];
        axon_signal = [axon_signal, mean(ind_trace,2)];  
    end
    
end


axon_signal = mean(axon_signal,2,'omitnan');

ad_corr.dend_pred = NaN(3900,8,5);
ad_corr.r2 = NaN(5,1);
ad_corr.ax_pred = NaN(3900,5);
ad_corr.ax_to_pred = NaN(3900,5);

if isempty(dendrite_1) == 0 && isempty(dendrite_2) == 0 && isempty(axon_signal) == 0
    cv_idx = 1:3900:19500;
    cumulative_idx = 1;
          
    for loop_idx2 = 1:5

        cv_interval = cv_idx(loop_idx2):cv_idx(loop_idx2)+3899;
        X = [dendrite_1, dendrite_2];
        Y = axon_signal;
        X(cv_interval,:) = [];
        Y(cv_interval,:) = [];    
        mdl = fitlm(X, Y);

        X_pred = [dendrite_1(cv_interval,:), dendrite_2(cv_interval,:)];
        Y_to_pred_data = axon_signal(cv_interval,1);
        Y_pred_model = predict(mdl,X_pred);
        r = corrcoef(Y_to_pred_data,Y_pred_model);
        r_2 = r(2)^2;
        ad_corr.r2(loop_idx2,1) = r_2;
        ad_corr.dend_pred(:,1:size(X_pred,2),loop_idx2) = X_pred;
        ad_corr.ax_pred(:,loop_idx2) = Y_pred_model;
        ad_corr.ax_to_pred(:,loop_idx2) = Y_to_pred_data;     
        cumulative_idx = cumulative_idx+1;

        x_data = zscore(Y_to_pred_data);
        y_model = zscore(Y_pred_model);
        
        % For quality check of 2-D histogram orientation (not to mix up x- and y-axes)
        quality_check = "off";
        if quality_check == "on"
            x_data(1:(length(x_data)/10)) = 0;
            y_model(1:(length(x_data)/10)) = 2;
        end
        
        % Computes bivariate histogram of z-scored data vs model.
        [scatter_pred_hist(:,:,loop_idx2), ad_corr.x_edges, ad_corr.y_edges] = ...
        histcounts2(x_data, y_model, 'NumBins', 15,...
        'XBinLimits', [-3, 3], 'YBinLimits', [-3, 3], 'Normalization', 'probability');

    end
    
    ad_corr.scatter_pred_hist = mean(scatter_pred_hist,3);

else
    ad_corr = NaN;
end

end

%%

function R2 = get_axo_dendritic_pooled_R2(ad_corr_table)

R2 = [];
for loop_idx1 = 1:size(ad_corr_table,1)
   curr_r2 = ad_corr_table.axo_dendritic_corr(loop_idx1,1).r2;
   curr_r2 = mean(reshape(curr_r2,size(curr_r2,1)*size(curr_r2,2),1),'omitnan');
   R2(loop_idx1) = curr_r2;
end

mean_pooled_R2 = mean(R2)
std_pooled_R2 = std(R2)

end


