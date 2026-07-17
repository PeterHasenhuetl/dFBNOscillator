% Analysis of dFBN GCaMP dynamics recorded from their somata.
% Code written by Peter Hasenhuetl.

clear all
tic

source_data_cellbody_details.data_path = []; %Add path as character array
source_data_cellbody_details.file_name = 'GCaMP_cellbody_recording.xlsx';

cb_data_table = [];
pooled_SW_proportion = [];
pooled_img_trace_array = [];
pooled_fly_id_array = [];
pooled_hemi_id_array = [];
cell_numbers = [];
pairwise_correllations = [];
pooled_auto_cov = [];
pooled_auto_amp = [];
pooled_power_spectrum = []; 
TTA_full = [];
TRANS_full = [];
CORR_full = [];
R_AUTO = [];
R_PWR = [];
POS_PWR = [];
SW = [];

% loop through individual flies
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
cum_idx = 0;
for fly_idx = 1:n_flies     
    fly_input = dir('*fly*.mat');
    loaded_fly = load(char(fly_input(fly_idx).name));
    curr_fly_name = char(fly_input(fly_idx).name);
    loaded_fly.details = 1; 
    % loop through extracted planes    
    if isempty(loaded_fly.hemi_1) == 0
        curr_h_idx = 1;
        curr_hemi_nu = get_curr_hemi_nu(loaded_fly, 1, curr_fly_name);
        [output_img, ind_data_table, curr_SW_proportion, curr_auto_amp, ...
            img_trace_array, fly_id_array, curr_pairwise_corr] = ...
            get_individual_cellbody_fly(curr_hemi_nu, fly_idx);      
        cb_data_table = [cb_data_table, ind_data_table.power_spectrum(1,1).trace_pwr];
        pooled_auto_cov = [pooled_auto_cov, ind_data_table.auto_cov.trace_auto];  
        pooled_power_spectrum = [pooled_power_spectrum, ind_data_table.power_spectrum.trace_pwr]; 
        pooled_SW_proportion = [pooled_SW_proportion; curr_SW_proportion];
        pooled_auto_amp = [pooled_auto_amp; curr_auto_amp];
        pooled_img_trace_array = [pooled_img_trace_array, img_trace_array];    
        pooled_fly_id_array = [pooled_fly_id_array; fly_id_array];    
        pooled_hemi_id_array = [pooled_hemi_id_array; ones(length(fly_id_array),1)];    
        pairwise_correllations = [pairwise_correllations, curr_pairwise_corr];    
        cell_numbers = [cell_numbers, size(curr_SW_proportion,1)];
        cum_idx = cum_idx+1;



        rhythmic_cell_ids_hemi1 = find((curr_auto_amp>0.05));
        if isempty(rhythmic_cell_ids_hemi1) == 0
            curr_cell = rhythmic_cell_ids_hemi1;
            x = length(curr_cell);
            n_TTAs = (size(img_trace_array,2)-1)*x;

            for loop_idx = 1:length(curr_cell)

                [TTA, TRANS_1, CORR] = get_cellbody_TTA(img_trace_array, curr_cell(loop_idx));

                TTA_full = [TTA_full, TTA];
                TRANS_full = [TRANS_full, TRANS_1];
                CORR_full = [CORR_full, CORR];
            end
            R_AUTO = [R_AUTO, ind_data_table.auto_cov.trace_auto(:,curr_cell)];
            R_PWR = [R_PWR, ind_data_table.power_spectrum.trace_pwr(:,curr_cell)];
            POS_PWR = [POS_PWR, output_img.pos_1(curr_cell)];
            SW = [SW; curr_SW_proportion(curr_cell)];
        end

    end
    
    if isempty(loaded_fly.hemi_2) == 0
        curr_h_idx = 2;
        curr_hemi_nu = get_curr_hemi_nu(loaded_fly, 2, curr_fly_name);
        [output_img, ind_data_table, curr_SW_proportion, curr_auto_amp, img_trace_array, fly_id_array, curr_pairwise_corr] = ...
            get_individual_cellbody_fly(curr_hemi_nu, fly_idx);     
        pairwise_correllations = [pairwise_correllations, curr_pairwise_corr];
        cell_numbers = [cell_numbers, size(curr_SW_proportion,1)];
        cb_data_table = [cb_data_table, ind_data_table.power_spectrum(1,1).trace_pwr];
        pooled_auto_cov = [pooled_auto_cov, ind_data_table.auto_cov.trace_auto];  
        pooled_power_spectrum = [pooled_power_spectrum, ind_data_table.power_spectrum.trace_pwr]; 
        pooled_SW_proportion = [pooled_SW_proportion; curr_SW_proportion];
        pooled_auto_amp = [pooled_auto_amp; curr_auto_amp];
        pooled_hemi_id_array = [pooled_hemi_id_array; 2*ones(length(fly_id_array),1)];
        pooled_img_trace_array = [pooled_img_trace_array, img_trace_array];
        pooled_fly_id_array = [pooled_fly_id_array; fly_id_array];
        cum_idx = cum_idx+1;



        rhythmic_cell_ids_hemi2 = find((curr_auto_amp>0.05));
        if isempty(rhythmic_cell_ids_hemi2) == 0
            curr_cell = rhythmic_cell_ids_hemi2;
            x = length(curr_cell);
            n_TTAs = (size(img_trace_array,2)-1)*x;
            for loop_idx = 1:length(curr_cell)

                [TTA, TRANS_1, CORR] = get_cellbody_TTA(img_trace_array, curr_cell(loop_idx));

                TTA_full = [TTA_full, TTA];
                TRANS_full = [TRANS_full, TRANS_1];
                CORR_full = [CORR_full, CORR];
            end
            R_AUTO = [R_AUTO, ind_data_table.auto_cov.trace_auto(:,curr_cell)];
            R_PWR = [R_PWR, ind_data_table.power_spectrum.trace_pwr(:,curr_cell)];
            POS_PWR = [POS_PWR, output_img.pos_1(curr_cell)];
            SW = [SW; curr_SW_proportion(curr_cell)];
            
        end

    end

end

pooled_SW_proportion = pooled_SW_proportion';
clust_ids_quart = get_clust_ids_quart(pooled_SW_proportion);
master_ids_vec = [pooled_fly_id_array+0.1.*pooled_hemi_id_array, clust_ids_quart];
id_array = unique(master_ids_vec(:,1));

cd([]) %Add path as character array

toc

pooled_auto_amp(isnan(pooled_auto_amp) == 1,1) = 0;
scatter((pooled_SW_proportion),pooled_auto_amp)


%%

close all
figure('Name','cellbody_plot','Color','white',...
    'Units','centimeters','Position',[5 6 8.9 10])
auto_plot = pooled_auto_cov;
curr_x = (output_img.lag_dist);
color = get_color;
fly_plot_id = [10, 1]; % first entry is fly; second entry is hemisphere
get_fig_panel_log_hist_plot(0.1+7-1.9, 2.5, pooled_SW_proportion, pooled_auto_amp, ...
    pooled_auto_cov, curr_x, pooled_fly_id_array, pooled_hemi_id_array, fly_plot_id(1), fly_plot_id(2),...
    color, source_data_cellbody_details)
get_fig_panel_cellbody_traces_v2(0.1, 2, pooled_img_trace_array, ...
    pooled_fly_id_array, pooled_hemi_id_array, fly_plot_id(1), fly_plot_id(2), ...
    pooled_SW_proportion, color, "yes", source_data_cellbody_details)

cd(source_data_cellbody_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'cell_body_fig.pdf')
cd([]) %Add path as character array

%%
x_pos = 2;
y_pos = 2;
close all

get_fig_panel_rhythmic_cell_power_spectrum(x_pos, y_pos, ...
    ind_data_table, R_PWR, POS_PWR, SW, pooled_power_spectrum)

get_fig_panel_cellbody_TTA(x_pos+4, y_pos, TTA_full, TRANS_full, CORR_full)


%% Custom functions

function curr_hemi_nu = get_curr_hemi_nu(loaded_fly, curr_hemi_id, curr_fly_name)

if curr_hemi_id == 1
    curr_hemi = loaded_fly.hemi_1;
    curr_ROI_centroids = loaded_fly.ROI_centroids_h1;
    curr_plane_ids = loaded_fly.plane_IDs_h1;
elseif curr_hemi_id == 2
    curr_hemi = loaded_fly.hemi_2;
    curr_ROI_centroids = loaded_fly.ROI_centroids_h2;
    curr_plane_ids = loaded_fly.plane_IDs_h2;
end
background = curr_hemi(:,end);
curr_hemi = curr_hemi(:,1:end-1);

for loop_idx1 = 1:size(curr_ROI_centroids,1)
    for loop_idx2 = 1:size(curr_ROI_centroids,1)
        centroids_1(loop_idx1,loop_idx2) = ...
            abs(curr_ROI_centroids(loop_idx1,1)-curr_ROI_centroids(loop_idx2,1));
        centroids_2(loop_idx1,loop_idx2) = ...
            abs(curr_ROI_centroids(loop_idx1,2)-curr_ROI_centroids(loop_idx2,2));
        centroids_1(loop_idx1,loop_idx1) = NaN;
        centroids_2(loop_idx1,loop_idx1) = NaN;
    end
end

cent_ids = zeros(size(centroids_1));
cent_ids(centroids_1 < 250 & centroids_2 < 250) = 1;

curr_traces = curr_hemi;
cell_coeffs = corrcoef(detrend(curr_traces));

[a, b] = find(cell_coeffs > 0.5 & cell_coeffs < 1);

cell_hits = unique(a);
cell_hits_orig = unique(a);

for loop_idx3 = 1:length(cell_hits)
    
    if isnan(cell_hits(loop_idx3)) == 0
        curr_selection_ids = [cell_hits(loop_idx3),b(a==cell_hits(loop_idx3))'];
        curr_array = curr_traces(:,curr_selection_ids);
        curr_winner = find(sum(curr_array,1) == max(sum(curr_array,1)));
        curr_eclusion_ids = curr_selection_ids;
        curr_eclusion_ids(curr_winner) = [];
    end
    for loop_idx4 = 1:length(curr_eclusion_ids)
        cell_hits(cell_hits == curr_eclusion_ids(loop_idx4)) = NaN;
    end
    
end

to_exclude = cell_hits_orig(isnan(cell_hits(:,1)) == 1,1);
curr_traces(:,to_exclude) = 0;
curr_traces = get_additonal_excl(curr_traces, curr_fly_name, curr_hemi_id);


plot_cond = 0;
if plot_cond == 1
	close all
    get_cellbody_exclusion_plot(loaded_fly, curr_hemi_id, curr_traces)        
    pause
end

excl_idx = sum(curr_traces,1) == 0;
curr_hemi_nu = [curr_traces(:,excl_idx == 0),background];

end

function [output_img, ind_data_table, curr_SW_proportion, curr_auto_amp, curr_img_trace_array, ...
    curr_fly_id_array, pairwise_corr] = get_individual_cellbody_fly(img_traces, fly_idx)

ipsi_corr = corrcoef(detrend(img_traces(101:end,1:end-1)));
n = size(ipsi_corr,2);
z = ((n*n)-n)/2;

pairwise_corr = [];
for loop_idx = 1:size(ipsi_corr,1)
    pairwise_corr = [pairwise_corr, ipsi_corr(loop_idx,loop_idx+1:end)];
end

output_img = get_cellbody_img_measures(img_traces);    
ind_data_table = get_cellbody_data_table_row(output_img, fly_idx, 1);    
[curr_SW_proportion, curr_auto_amp, curr_img_trace_array, curr_fly_id_array] = ...
    get_cell_body_table(ind_data_table);    

end

function output_img = get_cellbody_img_measures(imaging_traces)

f_r = 14.56; % frame-rate
trace_start = 101;
raw_green = imaging_traces(trace_start:end,1:end-1)-mean(imaging_traces(trace_start:end,end),1);
img_green = raw_green;
length_corr = size(raw_green,1);

for cell_idx = 1:size(img_green,2)
    
    r = img_green(:,cell_idx);
    dF_F = (smoothdata(r,'Gaussian',10));
    dF_F = zscore(dF_F); 
    output_img.dF_F(:,cell_idx) = dF_F;
    output_img.raw_img(:,cell_idx) = img_green(:,cell_idx);
     
    curr_img = detrend(dF_F);    
    diff_vec = diff(curr_img);
    [output_img.curr_autocov(:,cell_idx), lag_dist] = xcov(curr_img,length_corr,'coef');   
    PRD = get_auto_period(output_img.curr_autocov(:,cell_idx), 0.1, 0.1, 0.01, lag_dist/f_r);
    output_img.per_1(cell_idx) = PRD.curr_pr;
    output_img.amp_1(cell_idx) = PRD.curr_amp;
    output_img.prom_1(cell_idx) = PRD.p;
    if PRD.curr_pr > 0 && PRD.curr_pr <= 3
        output_img.per_1(cell_idx) = PRD.curr_pr;
        output_img.amp_1(cell_idx) = PRD.curr_amp;
        output_img.prom_1(cell_idx) = PRD.p;
    else
        output_img.per_1(cell_idx) = NaN;
        output_img.amp_1(cell_idx) = NaN;
        output_img.prom_1(cell_idx) = NaN;
    end
    output_img.curr_autocov_filt(:,cell_idx) = xcov(diff_vec,length_corr,'coef');
    PRD = get_auto_period(output_img.curr_autocov_filt(:,cell_idx), 0.1, 0.1, 0.01, lag_dist/f_r);
    if PRD.curr_pr > 0 && PRD.curr_pr <= 3
        output_img.per_2(cell_idx) = PRD.curr_pr;
        output_img.amp_2(cell_idx) = PRD.curr_amp;
        output_img.prom_2(cell_idx) = PRD.p;
    else
        output_img.per_2(cell_idx) = NaN;
        output_img.amp_2(cell_idx) = NaN;
        output_img.prom_2(cell_idx) = NaN;
    end
    [curr_spectrum, f] = pspectrum(curr_img,f_r,'FrequencyLimits',[0,14.56]);
    output_img.power_spectrum(:,cell_idx) = curr_spectrum;
    d_idx1 = find(f >= 0.2,1,'first');
    d_idx2 = find(f >= 2,1,'first');
    curr_peak = max(findpeaks(curr_spectrum(d_idx1:d_idx2)));
    curr_pos = find(curr_spectrum(d_idx1:d_idx2) == (max(findpeaks(curr_spectrum(d_idx1:d_idx2)))));
    output_img.p_amp_1(cell_idx) = curr_peak;
    output_img.p_amp_2(cell_idx) = max(curr_spectrum(d_idx1:d_idx2));
    output_img.diff_trace(:,cell_idx) = diff_vec;
    output_img.pos_1(cell_idx) = curr_pos;
end

output_img.lag_dist = lag_dist/f_r;
output_img.f = f;
ipsi_corr = corrcoef(output_img.diff_trace);
ipsi_corr(ipsi_corr == 1) = [];
output_img.ipsi_corr_mean = mean(ipsi_corr);
output_img.ipsi_corr = ipsi_corr;

end

function dm_row = get_cellbody_data_table_row(dm_input, fly_idx, plane_idx)

%% First, general variables
dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.plane_id = plane_idx;

%% Then, adds the imaging metrics 
dm_img_row = table;
dm_img_row.int_delta_power = dm_input.p_amp_1;
dm_img_row.per_1 = dm_input.per_1;
dm_img_row.amp_1 = dm_input.amp_1;
dm_img_row.prom_1 = dm_input.prom_1;
dm_img_row.per_2 = dm_input.per_2;
dm_img_row.amp_2 = dm_input.amp_2;
dm_img_row.prom_2 = dm_input.prom_2; 
dm_img_row.p_amp_1 = dm_input.p_amp_1;
dm_img_row.p_amp_2 = dm_input.p_amp_2;
dm_img_row.ipsi_corr_mean = dm_input.ipsi_corr_mean;

dm_img_row.ipsi_corr = dm_input.ipsi_corr;

trc_img = struct;
trc_img.trace_img = dm_input.dF_F;
dm_img_row.img_trace = trc_img;

trc_raw = struct;
trc_raw.trace_img = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;

pwr_spc = struct;
pwr_spc.trace_pwr = dm_input.power_spectrum;
dm_img_row.power_spectrum = pwr_spc;

frqu_trc = struct;
frqu_trc.trace_f = dm_input.f;
dm_img_row.frqu_trc = frqu_trc;

auto_cov = struct;
auto_cov.trace_auto = dm_input.curr_autocov;
dm_img_row.auto_cov = auto_cov;

lag_trc = struct;
lag_trc.lag_dist = dm_input.lag_dist;
dm_img_row.lag_trc = lag_trc;

%% Concatenates them

dm_row = [dm_behav_row, dm_img_row];

end

function [curr_SW_proportion, curr_auto_amp, img_trace_array, fly_id_array] = get_cell_body_table(curr_table)

curr_SW_proportion = curr_table.p_amp_1';
curr_auto_amp = curr_table.prom_1';

img_trace_array = curr_table.img_trace.trace_img;
fly_id_array = curr_table.fly_id(1)*ones(size(img_trace_array,2),1);

end

function PRD = get_auto_period(img_signal, min_dist, min_height, min_prom, curr_lag_dist)

[peak_amp, lag_idx, ~, p] = findpeaks(img_signal, 'MinPeakDistance', min_dist, 'MinPeakProminence', min_prom,...
    'MinPeakHeight', min_height);
lag_time = abs(curr_lag_dist(lag_idx));

if size(lag_idx,1) == 1
    PRD.curr_amp = NaN;
    PRD.curr_pr = NaN;
    PRD.p = NaN;
else
    [nu_l, curr_idx] = sort(lag_time, 'ascend');
    PRD.curr_amp = (peak_amp(curr_idx(2)));
    PRD.p = (p(curr_idx(2)));
    PRD.curr_pr = abs(nu_l(2));
end

end

function clust_ids_quart = get_clust_ids_quart(pooled_SW_proportion)

quart_1 = prctile(pooled_SW_proportion,25);
quint_2 = prctile(pooled_SW_proportion,50);
quint_3 = prctile(pooled_SW_proportion,75);

idx_q1 = find(pooled_SW_proportion <= quart_1);
idx_q2 = find(pooled_SW_proportion > quart_1 & pooled_SW_proportion <= quint_2);
idx_q3 = find(pooled_SW_proportion > quint_2 & pooled_SW_proportion <= quint_3);
idx_q4 = find(pooled_SW_proportion > quint_3);

clust_ids_quart = zeros(length(pooled_SW_proportion),1);
clust_ids_quart(idx_q1,1) = 1;
clust_ids_quart(idx_q2,1) = 2;
clust_ids_quart(idx_q3,1) = 3;
clust_ids_quart(idx_q4,1) = 4;


end

function curr_traces = get_additonal_excl(curr_traces, curr_fly_name, curr_hemi_id)
% additional exclusions based on manual cross-checking of ROIs.

if curr_fly_name == "fly20220425_1DFB_cellbodies.mat" && curr_hemi_id == 1
    curr_traces(:,7) = 0;
elseif curr_fly_name == "fly20230321_4DFB_cellbodies.mat" && curr_hemi_id == 2
    curr_traces(:,[2,4]) = 0;
elseif curr_fly_name == "fly20230322_2DFB_cellbodies.mat" && curr_hemi_id == 1
    curr_traces(:,[4,5]) = 0;
end

end

function [TTA, TRANS_1, CORR] = get_cellbody_TTA(pooled_img_trace_array, curr_cell)

img_trace = pooled_img_trace_array(:,curr_cell);

[~, ~, trans_ids] = get_GCaMP_transients(img_trace);

trans_ids(1) = [];
trans_ids(end) = [];
excl_cell = curr_cell;
cells_to_analyse = 1:size(pooled_img_trace_array,2);
cells_to_analyse(excl_cell) = [];
TTA = [];

pooled_img_trace_array = ((pooled_img_trace_array));
n_wind = round(0.5*14.56);
p_wind = round(2*14.56);

for loop_idx1 = 1:length(cells_to_analyse)
    curr_TTA = [];
    for loop_idx2 = 1:length(trans_ids)
        curr_TTA(:,loop_idx2) = ...
            ((pooled_img_trace_array(trans_ids(loop_idx2)-n_wind:trans_ids(loop_idx2)+p_wind,cells_to_analyse(loop_idx1))));
    end
    TTA = [TTA, mean(curr_TTA,2)];

    [a, ~] = corrcoef(diff((pooled_img_trace_array(:,cells_to_analyse(loop_idx1)))),...
            diff((pooled_img_trace_array(:,excl_cell))));
        curr_corr(1,loop_idx1) = a(2);

end
CORR = curr_corr;

for loop_idx2 = 1:length(trans_ids)
        curr_trans(:,loop_idx2) = ...
            ((pooled_img_trace_array(trans_ids(loop_idx2)-n_wind:trans_ids(loop_idx2)+p_wind,excl_cell)));
        
end

TRANS_1 = mean(curr_trans,2);

end

function [transient_vec, transient_amplitudes, trans_ids] = get_GCaMP_transients(img_trace)

gaussian_sliding_window = 15;
min_peak_height = 0.0001;
trans_signal = smoothdata(diff(img_trace), 'gaussian', gaussian_sliding_window);
[~, trans_ids] = findpeaks(trans_signal, 'MinPeakHeight', min_peak_height);

padded_img = [img_trace; zeros(20,1)];
transient_amplitudes = NaN(1,length(trans_ids));
for loop_idx = 1:length(trans_ids)
    transient_amplitudes(loop_idx) = max(padded_img(trans_ids(loop_idx):trans_ids(loop_idx)+10));   
end

transient_vec = zeros(length(img_trace),1);
transient_vec(trans_ids,1) = 1;

end

%% Plotting functions

function get_fig_panel_cellbody_traces_v2(x_pos, y_pos, pooled_img_trace_array, ...
    pooled_fly_id_array, pooled_hemi_id_array, curr_fly_id, curr_hemi_id, pooled_SW_proportion, ...
    color, cond_labeling, source_data_cellbody_details)

img_idx = find(pooled_fly_id_array == curr_fly_id & pooled_hemi_id_array == curr_hemi_id);
full_trace = pooled_img_trace_array(:,img_idx);
curr_SW = pooled_SW_proportion(1,img_idx);

[~, b] = sort(curr_SW,'ascend');
full_trace = full_trace(:,b);

ht_1 = 1.5;
ht_trace = ht_1/size(full_trace,2);

f_s = get_default_font_size;
sz_1 = 3;
cartoon_dist = 1.2;
sz_inset = sz_1-cartoon_dist;
y_dist_traces = 0;
ht_2 = 1;
y_scale_bar_size = 1;
y_dist_inset = 0.4;

lim_1_inset = 1401;
lim_2_inset = 2100;

for idx_plot = 1:size(full_trace,2)

    curr_trace = detrend(full_trace(:,idx_plot));
    curr_trace = curr_trace-min(curr_trace);
    curr_trace = curr_trace./max(curr_trace);
    
    pnl_1 = axes('Units', 'Centimeters', 'Position',...
        [x_pos, y_pos+((idx_plot-1)*(ht_trace+y_dist_traces)), sz_1, ht_trace]);     
    plot(pnl_1, curr_trace, 'Color', color.dark_gray, 'LineWidth', 0.25)
    hold on
    if idx_plot == size(full_trace,2)
        plot(pnl_1, [lim_1_inset, lim_1_inset], [0,1], 'Color', color.red, 'LineWidth', 0.25)
        plot(pnl_1, [lim_2_inset, lim_2_inset], [0,1], 'Color', color.red, 'LineWidth', 0.25)
    end
    pnl_1.XAxis.Color = 'none';
    pnl_1.YAxis.Color = 'none';
    pnl_1.Color = 'none';  
    ylim([min(curr_trace), max(curr_trace)])
    xlim([1, length(curr_trace)])
    
    if idx_plot == 1
        % imaging scale bar    
        pnl_1 = axes('Units', 'Centimeters', 'Position',...
            [x_pos+sz_1+0.1, y_pos+((idx_plot-1)*(ht_trace+y_dist_traces)), 0.2, ht_trace]);  
        plot(pnl_1, [0, 0], [min(curr_trace), min(curr_trace)+y_scale_bar_size], 'Color', 'k', 'LineWidth', get_default_scale_bar_width)    
        pnl_1.XAxis.Color = 'none';    
        pnl_1.YAxis.Color = 'none';    
        pnl_1.Color = 'none';    
        ylim([min(curr_trace), max(curr_trace)])    
        xlim([0, 1])
        
        pnl_1 = axes('Units', 'Centimeters', 'Position',...
            [x_pos, y_pos+((idx_plot-1)*(ht_trace+y_dist_traces))-0.1, sz_1, 0.1]);  
        plot(pnl_1, [length(curr_trace)-(14.56*10), length(curr_trace)], [0, 0], 'k', 'LineWidth', get_default_scale_bar_width) 
        pnl_1.XAxis.Color = 'none';    
        pnl_1.YAxis.Color = 'none';    
        pnl_1.Color = 'none';    
        ylim([0, 1])
        xlim([1, length(curr_trace)])
        get_default_annotation(x_pos+sz_1, y_pos+((idx_plot-1)*(ht_trace+y_dist_traces))-0.15, ...
        '10 s', 'k', 'normal', "right")
 
    end
    
end

if cond_labeling == "yes"
    annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos+sz_1, y_pos-0.15, 5, 0.5], ...
        'string', {'1 a.U.'}, 'EdgeColor', 'none', 'FontSize', f_s, 'FontWeight', 'normal')
end

y_scale_bar_size_inset = 0.2;
curr_trace = curr_trace(lim_1_inset:lim_2_inset,1);
y_pos_inset = y_pos+((idx_plot-1)*(ht_trace+y_dist_traces))+ht_trace+y_dist_inset;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+cartoon_dist, y_pos_inset, sz_inset, ht_2]);
plot(pnl_1, curr_trace, 'Color', color.dark_gray, 'LineWidth', 0.25)    
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim([min(curr_trace), max(curr_trace)])    
xlim([1, length(curr_trace)])

% imaging scale bar   
pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+cartoon_dist+sz_inset+0.1, y_pos_inset, 0.2, ht_2]);
plot(pnl_1, [0, 0], [min(curr_trace), min(curr_trace)+y_scale_bar_size_inset],...
    'Color', 'k', 'LineWidth', get_default_scale_bar_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim([min(curr_trace), max(curr_trace)])
xlim([0, 1])
        
pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+cartoon_dist, y_pos_inset-0.1, sz_inset, 0.1]);  
plot(pnl_1, [length(curr_trace)-(14.56*5), length(curr_trace)], [0, 0], 'k', 'LineWidth', get_default_scale_bar_width) 
pnl_1.XAxis.Color = 'none';    
pnl_1.YAxis.Color = 'none';    
pnl_1.Color = 'none';    
ylim([0, 1])
xlim([1, length(curr_trace)])
get_default_annotation(x_pos+cartoon_dist+sz_inset, y_pos_inset-0.2, ...
    '5 s', 'k', 'normal', "right")

if cond_labeling == "yes"
    annotation('textbox', 'Units','centimeters','Position',[x_pos+cartoon_dist+sz_inset, y_pos_inset-0.1, 5, 0.5], ...
        'string', {'0.2'; 'a.U.'}, 'EdgeColor','none','FontSize',f_s,'FontWeight','normal')
end

%% Saves the source data

example_traces = table;
example_traces.example_fly = full_trace;
writetable(example_traces, ...
    [source_data_cellbody_details.data_path, source_data_cellbody_details.file_name], 'Sheet', 'example_traces')

end

function get_fig_panel_log_hist_plot(x_pos, y_pos, pooled_SW_proportion, ...
    pooled_auto_amp, pooled_auto_cov, curr_x, pooled_fly_id_array, ...
    pooled_hemi_id_array, curr_fly_id, curr_hemi_id,color, source_data_cellbody_details)

img_idx = find(pooled_fly_id_array == curr_fly_id & pooled_hemi_id_array == curr_hemi_id);
curr_SW = pooled_SW_proportion(1,img_idx);
pooled_auto_cov = pooled_auto_cov(:,img_idx);
[~, b] = sort(curr_SW,'ascend');
curr_auto = pooled_auto_cov(:,b);

sz_1 = 1.8;
ht_1 = 1.2;
xlim_hist = [-3.5, 0.5];
ylm_hist = [0, 60];
ylm_prop = [0, 1];
prop_color = color.medium_gray;
hist_color = [0, 0, 0];

curr_pooled_auto_amp = pooled_auto_amp;

pooled_auto_amp_bin = zeros(length(pooled_auto_amp),3);
pooled_auto_amp_bin(curr_pooled_auto_amp > 0.05,1) = 1;
pooled_auto_amp_bin(curr_pooled_auto_amp > 0.2,2) = 1;
pooled_auto_amp_bin(curr_pooled_auto_amp > 0.3,3) = 1;
curr_SW_prop = log10(pooled_SW_proportion);
[a_1, a_2] = histcounts(curr_SW_prop, 'Normalization', 'count', 'NumBins', 9, 'BinLimits', [-3.5, 0]);


proportion_vec = [];
n_vec = [];
for loop_idx = 1:length(a_2)-1
   idx_1 = find(curr_SW_prop > a_2(loop_idx) & curr_SW_prop <= a_2(loop_idx+1));
   proportion_data = pooled_auto_amp_bin(idx_1,1);
   proportion_vec(loop_idx) = sum(proportion_data)/length(proportion_data);
   n_vec(loop_idx) = length(proportion_data);
   SW_vec(loop_idx) = sum(proportion_data);
end

plot_vec = a_2+mean(diff(a_2))/2;
plot_vec = plot_vec(1:end-1);
n_cells = size(pooled_SW_proportion,2);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
histogram(curr_SW_prop, 'Normalization', 'count', 'NumBins', 9, 'BinLimits', [-3.5, 0],...
    'DisplayStyle', 'stairs', 'EdgeColor', hist_color, 'LineWidth', 0.5);
get_default_separated_ax(curr_pnl, xlim_hist(1), xlim_hist(2), xlim_hist(1), ...
    xlim_hist(2), ylm_hist(1), ylm_hist(2),...
    ylm_hist(1), ylm_hist(2), 2, 20, "linear", "linear", 'Log(SO peak)', 'n cells',...
    'k', 'k', 'k', 'k', sz_1, ht_1)
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
plot(plot_vec, proportion_vec, 'Color', prop_color, 'LineWidth', 1)
curr_pnl.YAxisLocation = 'right';
get_default_separated_ax(curr_pnl, xlim_hist(1), xlim_hist(2), [], ...
    [], ylm_prop(1), ylm_prop(2),...
    ylm_prop(1), ylm_prop(2), [], 0.5, "linear", "linear", [], 'Fraction SO cells',...
    'none', 'k', 'none', prop_color, sz_1, ht_1)


get_default_annotation(x_pos, y_pos+ht_1, [num2str(n_cells), ' cells'], 'k', 'normal', "left")
get_default_annotation(x_pos+sz_1, y_pos+(ht_1*0.35), ...
    [num2str(sum(SW_vec)), ' cells'], prop_color, 'normal', "right")
percent_cond = 0;
if percent_cond == 1
    get_default_annotation(x_pos+sz_1, y_pos+(ht_1*0.35)-0.15, ...
        ['~',num2str(round((sum(SW_vec)/n_cells)*100)), ' %'], prop_color, 'normal', "right")
end

idx_1 = 1;
idx_2 = size(curr_auto,2);
get_fig_panel_auto_corr(curr_auto, ...
    curr_x, idx_1, idx_2, x_pos, y_pos+ht_1+0.4, color, source_data_cellbody_details)


%% saves the source data

cell_body_SO_histogram = table;
cell_body_SO_histogram.SO_peak_bins = plot_vec';
cell_body_SO_histogram.SO_peak_counts = a_1';
cell_body_SO_histogram.fraction_SO_cell = proportion_vec';
writetable(cell_body_SO_histogram, ...
    [source_data_cellbody_details.data_path, source_data_cellbody_details.file_name], 'Sheet', 'cell_body_SO_histogram')

end

function get_fig_panel_auto_corr(curr_auto, ...
    curr_x, idx_1, idx_2, x_pos, y_pos, color, source_data_cellbody_details)

sz_2 = 0.85;
ht_1 = 0.8;
ht_2 = 0.8;
xlim_auto = [-8, 8];
ylim_auto = [-0.5, 1];
auto_color = color.dark_gray;
dotted_line_width = 0.25;
autocorr_line_width = 0.25;

scale_bar = 5;
curr_auto_plot_1 = curr_auto(:,idx_1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_2, ht_2]);
plot(curr_pnl, xlim_auto, [0, 0], 'k:', 'LineWidth', dotted_line_width)
hold on
plot(curr_pnl, curr_x, curr_auto_plot_1, 'LineWidth', autocorr_line_width, 'Color', auto_color)
plot(curr_pnl, [xlim_auto(2)-scale_bar, xlim_auto(2)], [ylim_auto(1), ylim_auto(1)],...
    'k', 'LineWidth', get_default_scale_bar_width)
get_default_separated_ax(curr_pnl, xlim_auto(1), xlim_auto(2), [], ...
    [], ylim_auto(1), ylim_auto(2),...
    ylim_auto(1), ylim_auto(2), [], 0.5, "linear", "linear", [], 'r',...
    'none', 'k', 'none', 'k', sz_2, ht_1)
get_default_annotation(x_pos+sz_2, y_pos+0.2, ...
    '5 s', 'k', 'normal', "right")

curr_auto_plot_2 = curr_auto(:,idx_2);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_2+0.1, y_pos, sz_2, ht_2]);
plot(curr_pnl, xlim_auto, [0,0], 'k:', 'LineWidth', dotted_line_width)
hold on
plot(curr_pnl, curr_x, curr_auto_plot_2, 'LineWidth', autocorr_line_width, 'Color', auto_color)
get_default_separated_ax(curr_pnl, xlim_auto(1), xlim_auto(2), [], ...
   [], ylim_auto(1), ylim_auto(2), [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_2, ht_2)

%% Saves the source data
lag_auto = curr_x.';
example_autocorr = table;
example_autocorr.example_autocorrs = ...
    [lag_auto, curr_auto_plot_1, curr_auto_plot_2];
writetable(example_autocorr, ...
    [source_data_cellbody_details.data_path, source_data_cellbody_details.file_name],...
    'Sheet', 'example_autocorrelations')

end

function get_cellbody_exclusion_plot(loaded_fly, curr_hemi_id, curr_hemi_nu)

if curr_hemi_id == 1
    curr_hemi = loaded_fly.hemi_1;
    curr_ROI_centroids = loaded_fly.ROI_centroids_h1;
    curr_plane_ids = loaded_fly.plane_IDs_h1;
elseif curr_hemi_id == 2
    curr_hemi = loaded_fly.hemi_2;
    curr_ROI_centroids = loaded_fly.ROI_centroids_h2;
    curr_plane_ids = loaded_fly.plane_IDs_h2;
end

y_pos = 0.5;
ht_1 = 1;
sz_1 = 18;
close all
figure('Name', 'cellbody_exclusion_plot', 'Color', 'white',...
    'Units', 'centimeters', 'Position', [5, 6, 19, 27])
for loop_idx = 1:size(curr_hemi,2)
    curr_pnl = axes('Units', 'Centimeters', 'Position', [1, y_pos+(loop_idx-1)*ht_1, sz_1, ht_1]);
    plot(zscore(curr_hemi(:,loop_idx)))
    hold on
    if loop_idx <= size(curr_hemi_nu,2)
        plot(zscore(curr_hemi_nu(:,loop_idx)))
    end
    curr_pnl.XAxis.Color = 'none';
    curr_pnl.YAxis.Color = 'none';
    annotation('textbox', 'Units', 'centimeters', 'Position', [1+sz_1-4, y_pos+(loop_idx-1)*ht_1, 5, 0.5], ...
        'string', [num2str(curr_ROI_centroids(loop_idx,1)), '   ', num2str(curr_ROI_centroids(loop_idx,2))], ...
        'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'normal')
end

end

function get_fig_panel_cellbody_TTA(x_pos, y_pos, TTA_full, TRANS_full, CORR_full)

sz_1 = 1.5;
sz_2 = 1;
ht_1 = 1;
ht_2 = 0.75;
ht_3 = 0.75;
ht_4 = 0.75;
xlm_img = [0, 2.5];
xlm_r = [-0.3, 0.6];
ylm_map = [0.5, size(TTA_full,2)+0.5];
ylm_hist = [0, 0.5];
ylm_img_TRANS = [-1, 1];
ylm_img_TTA = [-0.1, 0.1];
clim_1 = [-0.2, 0.2];
major_tick_cell = 25;
major_tick_time = 1;
major_tick_img_TRANS = 1;
major_tick_img_TTA = 0.1;
y_dist_1 = 0.1;
y_dist_2 = 0.3;
x_dist_1 = 0;
red_3 = [200, 85, 97];
color_1 = [0.5, 0.5, 0.5];
color_map_1 = get_default_two_color_map_dark(5000, 5000, 0.6, 0.6, color_1, (red_3)./255);


[a, b] = sort(CORR_full,'ascend');
a = round(a,2);
TTA_full = TTA_full(:,b);
TTA_full = TTA_full-mean(TTA_full,1);

window_1 = 20;
TTA_full1 = TTA_full(:,1:window_1);
TTA_full2 = TTA_full(:,(window_1+1):(end-window_1));
TTA_full3 = TTA_full(:,end-(window_1-1):end);


t_v = (1:size(TTA_full,1))./14.56;
t_v = t_v-t_v(1);
pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos, y_pos, sz_1, ht_1]);
imagesc(pnl_1, t_v, [1, size(TTA_full,2)],TTA_full')
clim(clim_1)
colormap(pnl_1, color_map_1)
get_default_ax(pnl_1, xlm_img(1), xlm_img(2), xlm_img(1), xlm_img(2), ...
    ylm_map(1), ylm_map(2), 0, ylm_map(2),...
    major_tick_time, major_tick_cell, "linear", "linear", 'Time (s)', 'TAA',...
    'none', 'k', 'none', 'k', sz_1, ht_1)
min_max = [min(TTA_full,[],"all"), max(TTA_full,[],"all")];
get_default_colorbar(pnl_1, x_pos, y_pos-0.2, sz_1/2, 0.1, clim_1, min_max, '\DeltaF/F', 'southoutside')

pnl_scalebar = axes('Units', 'Centimeters', 'Position',...
    [x_pos, y_pos-0.1, sz_1, 0.1]);
plot(pnl_scalebar, [xlm_img(2)-0.5, xlm_img(2)], [0, 0], 'Color', [0, 0, 0],...
    'LineWidth', get_default_scale_bar_width)
get_default_ax(pnl_scalebar, xlm_img(1), xlm_img(2), [], [], ...
    0, 1, [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, 0.1)
get_default_annotation(x_pos+sz_1+x_dist_1, y_pos-0.2, ...
    '0.5 s', 'k', 'normal', "right")


pnl_2 = axes('Units', 'Centimeters', 'Position',...
    [x_pos, y_pos+ht_1+y_dist_1, sz_1, ht_2]);
hold on
get_default_SEM_area_plot(pnl_2, TTA_full1, t_v, [0.5, 0.5, 0.5])
get_default_SEM_area_plot(pnl_2, TTA_full2, t_v, [0.25, 0.25, 0.25])
get_default_SEM_area_plot(pnl_2, TTA_full3, t_v, [0, 0, 0])
get_default_ax(pnl_2, xlm_img(1), xlm_img(2), [], [], ...
    ylm_img_TTA(1), ylm_img_TTA(2), ylm_img_TTA(1), ylm_img_TTA(2),...
    [], major_tick_img_TTA, "linear", "linear", [], '\DeltaF/F',...
    'none', 'k', 'none', 'k', sz_1, ht_2)
get_default_annotation(x_pos+sz_1, y_pos+(ht_1+y_dist_1)+ht_2, ...
    'TAA', 'k', 'normal', "right")

pnl_3 = axes('Units','Centimeters','Position',...
    [x_pos, y_pos+(ht_1+y_dist_1)+(ht_2+y_dist_2), sz_1, ht_3]);
get_default_SEM_area_plot(pnl_3, TRANS_full, t_v, [0, 0, 0])
get_default_ax(pnl_3, xlm_img(1), xlm_img(2), [], [], ...
    ylm_img_TRANS(1), ylm_img_TRANS(2), ylm_img_TRANS(1), ylm_img_TRANS(2),...
    [], major_tick_img_TRANS, "linear", "linear", [], '\DeltaF/F',...
    'none', 'k', 'none', 'k', sz_1, ht_3)

marker_size = 5;
triangle_color = [0.3, 0.3, 0.3];
pnl_4 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+sz_1+x_dist_1, y_pos, sz_2, ht_1]);
hold on
plot(pnl_4, [0, 0], [ylm_map(1), ylm_map(2)], 'k:', 'LineWidth', 0.5)
plot(pnl_4, a,1:length(a), 'k', 'LineWidth', 0.75)
scatter(pnl_4, mean(a), ylm_map(1), marker_size, 'Marker', 'v', 'MarkerFaceColor', triangle_color, 'MarkerEdgeColor', 'none')
plot(pnl_4, [xlm_r(1), xlm_r(2)], [(window_1+0.5), (window_1+0.5)], 'k:')
plot(pnl_4, [xlm_r(1), xlm_r(2)], [(length(a)-window_1-0.5), (length(a)-window_1-0.5)], 'k:')


pnl_4.YDir = "reverse";
get_default_ax(pnl_4, xlm_r(1), xlm_r(2), xlm_r(1), xlm_r(2), ...
    ylm_map(1), ylm_map(2), ylm_map(1), ylm_map(2),...
    [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_2, ht_1)

pnl_4 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+sz_1+x_dist_1, y_pos-0.1, sz_2, 0.1]);
plot(pnl_4, [xlm_r(2)-0.4, xlm_r(2)], [0, 0], 'Color', [0, 0, 0],...
    'LineWidth', get_default_scale_bar_width)
get_default_ax(pnl_4, xlm_r(1), xlm_r(2), xlm_r(1), xlm_r(2), ...
    0, 1, 0, 1, [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_2, 0.1)
get_default_annotation(x_pos+sz_1+x_dist_1+sz_2, y_pos-0.2, ...
    '0.4 r', 'k', 'normal', "right")

[a_1, a_2] = histcounts(a, 'Normalization', 'probability', 'NumBins', 15);
plot_vec = a_2+mean(diff(a_2))/2;
plot_vec = plot_vec(1:end-1);
hist_color = [0, 0, 0];
scale_bar_xpos = 0.25;
scale_bar_ypos = 0.15;
pnl_5 = axes('Units', 'Centimeters', 'Position',...
    [x_pos+sz_1+x_dist_1, y_pos+(ht_1+y_dist_1), sz_2, ht_4]);
bar(pnl_5, plot_vec, a_1, 1, 'FaceColor', hist_color, 'EdgeColor', 'none')
hold on
plot(pnl_5, [scale_bar_xpos, scale_bar_xpos], [scale_bar_ypos, scale_bar_ypos+0.2], 'Color', [0, 0, 0],...
    'LineWidth', get_default_scale_bar_width)
get_default_ax(pnl_5, xlm_r(1), xlm_r(2), xlm_r(1), xlm_r(2), ...
    ylm_hist(1), ylm_hist(2), ylm_hist(1), ylm_hist(2),...
    [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_2, ht_1)
get_default_annotation_rotated(x_pos+sz_1+x_dist_1+0.75, y_pos+(ht_1+y_dist_1)+0.45, ...
    '20%', 'k', 'normal', "left")

end

function get_fig_panel_rhythmic_cell_power_spectrum(x_pos, y_pos, ...
    ind_data_table, R_PWR, POS_PWR, SW, pooled_power_spectrum)


sz_1 = 1.8;
ht_1 = 1.8;
ht_2 = 1;
ylm_map = [0.5, size(R_PWR,2)+0.5];
ylm_img_TTA = [0, 0.06];
xlm_1 = [0, 3];
major_tick_cell = 7;
major_tick_frequency = 1;
major_tick_power = 0.03;
clim_1 = [0, 1];
y_dist_1 = 0.1;
RGB_1 = [0, 0, 0];
col_1(:,1) = linspace(1,RGB_1(1),100);
col_1(:,2) = linspace(1,RGB_1(2),100);
col_1(:,3) = linspace(1,RGB_1(3),100);    
color_map_1 = [col_1(:,1), col_1(:,2), col_1(:,3)];

[~, b] = sort(POS_PWR,'ascend');
norm_PWR = R_PWR./(SW');
curr_pwr = norm_PWR(:,b);

f_v = ind_data_table.frqu_trc.trace_f;
pnl_1 = axes('Units', 'Centimeters', 'Position',...
    [x_pos, y_pos+ht_2+y_dist_1, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(pnl_1, R_PWR, f_v', [0, 0, 0])
get_default_SEM_area_plot(pnl_1, pooled_power_spectrum, f_v', [0.5, 0.5, 0.5])
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ...
    ylm_img_TTA(1), ylm_img_TTA(2), ylm_img_TTA(1), ylm_img_TTA(2),...
    [], major_tick_power, "linear", "linear", [], 'Power',...
    'none', 'k', 'none', 'k', sz_1, ht_1)
get_default_annotation(x_pos+sz_1, y_pos+(ht_2+y_dist_1)+ht_1+0.2, ...
    'Rhythmic dFBNs', 'k', 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos+(ht_2+y_dist_1)+ht_1, ...
    'All dFBNs', [0.5, 0.5, 0.5], 'normal', "right")


pnl_2 = axes('Units', 'Centimeters', 'Position',...
    [x_pos, y_pos, sz_1, ht_2]);
imagesc(pnl_2, ind_data_table.frqu_trc.trace_f, 1:ylm_map(2), curr_pwr')
clim(clim_1)
colormap(pnl_2, color_map_1)
get_default_ax(pnl_2, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ...
    ylm_map(1), ylm_map(2), 0, ylm_map(2),...
    major_tick_frequency, major_tick_cell, "linear", "linear", 'Frequency (Hz)', 'Cell',...
    'k', 'k', 'k', 'k', sz_1, ht_2)
min_max = [min(curr_pwr,[],"all"), max(curr_pwr,[],"all")];
get_default_colorbar(pnl_2, x_pos+sz_1+0.1, y_pos, 0.1, ht_2, clim_1, min_max, 'Norm. Power', 'eastoutside')


end







