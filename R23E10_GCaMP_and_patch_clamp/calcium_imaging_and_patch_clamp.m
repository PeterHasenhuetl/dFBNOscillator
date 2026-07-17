% Analysis of simultaneous calcium imaging and patch-clamp.
% Code written by Peter Hasenhuetl.

tic
clear all
fly_IDs = char('fly01','fly02','fly03','fly04','fly05','fly06','fly07','fly08','fly09','fly10');
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = '2P_patch_data.xlsx';


% go into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);

% loop through individual flies
for fly_idx = 1:n_flies
    
     cd([]) %Add path as character array
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        fly_name = (char(fly_input(fly_idx).name));
        analyzed_flies.(fly_IDs(fly_idx,:)).fly_name = fly_name;
        cd([]) %Add path as character array
        [analyzed_flies.(fly_IDs(fly_idx,:)).align_img, analyzed_flies.(fly_IDs(fly_idx,:)).t_v, ...
            analyzed_flies.(fly_IDs(fly_idx,:)).AP, ...
            analyzed_flies.(fly_IDs(fly_idx,:)).b_w] = ...
            get_img_align(loaded_fly.patch_array,loaded_fly.img_array(:,[1,2]), fly_idx);
        curr_fly = analyzed_flies.(fly_IDs(fly_idx,:));
        analyzed_flies.(fly_IDs(fly_idx,:)).rev_corr = ...
            get_IRF_estimate(curr_fly, "impulse");
        analyzed_flies.(fly_IDs(fly_idx,:)).transient_aligned_PSTH = get_transient_aligned_PSTH(curr_fly);          
        analysis_summary.linear_filter(:,fly_idx) = ...
            analyzed_flies.(fly_IDs(fly_idx,:)).rev_corr.linear_filter;
        
end

norm_fac_lf = mean(max(analysis_summary.linear_filter));
analysis_summary.linear_filter = analysis_summary.linear_filter/norm_fac_lf;

cd([]) %Add path as character array

% Quantifies out-of-sample model performance
fields_flies = fieldnames(analyzed_flies);
analysis_summary.scatter_hist_all = [];
analysis_summary.R2 = [];

for fly_idx = 1:length(fields_flies)
    
    fly_name = fields_flies{fly_idx};
    input_fly = analyzed_flies.(string(fly_name));        

    analyzed_flies.(string(fly_name)).model_prediction = get_rev_corr_pred(input_fly, analysis_summary, fly_idx);
    analysis_summary.scatter_hist_all(:,:,fly_idx) = ...
        analyzed_flies.(string(fly_name)).model_prediction.scatter_pred_hist;
    analysis_summary.R2(fly_idx) = analyzed_flies.(string(fly_name)).model_prediction.R2;
    analysis_summary.cellbody_GCaMP_transients(:,fly_idx) = ...
        mean(analyzed_flies.(string(fly_name)).transient_aligned_PSTH.cellbody_GCaMP_transients,2);
    analysis_summary.aligned_PSTH(:,fly_idx) = ...
        mean(analyzed_flies.(string(fly_name)).transient_aligned_PSTH.aligned_PSTH,2);

end

analysis_summary.curr_diff_img = ([zeros(1,10); diff((analysis_summary.cellbody_GCaMP_transients))]);
mean_img = mean(analysis_summary.curr_diff_img,2);
analysis_summary.curr_diff_img = analysis_summary.curr_diff_img/max(mean_img);
mean_img = mean(analysis_summary.curr_diff_img,2);
trans_onset = find(mean_img >= 0,2,'first');
analysis_summary.aligned_tv = analyzed_flies.fly01.transient_aligned_PSTH.t_v1;
analysis_summary.aligned_tv = analysis_summary.aligned_tv-(mean(diff(analysis_summary.aligned_tv))*trans_onset(2));

toc

%% Plots the data
color = get_color;
rev_corr = analyzed_flies.fly04.rev_corr;
model_prediction = analyzed_flies.fly04.model_prediction;

close all
figure('Name','Fig1_1', 'Color','White','Units','centimeters','Position',[30, 10, 8.9, 24]);
get_fig_panel_reverse_correlation_example(rev_corr, model_prediction, analysis_summary, 0.5, 4, color, source_data_details);
get_fig_panel_data_vs_model(model_prediction, analysis_summary.scatter_hist_all, 1.5, 1.5, source_data_details);
get_fig_panel_GCaMP_aligned_PSTH(analysis_summary, 7.5-1.8, 1.5, color, source_data_details)
%
cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'2P_patch_fig.pdf')

cd([]) %Add path as character array


%% Power spectra of ephys and imaging traces


sliding_window = 500;
all_cells = fieldnames(analyzed_flies);
pooled_f = [];
pooled_v = [];
for cell_idx = 1:size(all_cells,1)

    cell_name = all_cells{cell_idx};
    input_cell = analyzed_flies.(string(cell_name));
    curr_f = zscore((downsample(input_cell.align_img(:,1),40)));
    curr_v = zscore(detrend((downsample(input_cell.align_img(:,2),40))));

    clear pwr_f pwr_v
    for i = 1:(length(curr_f)-sliding_window)
        [pwr_f(:,i), ~] = pspectrum(curr_f(i:(i+sliding_window-1),:),50);
        [pwr_v(:,i), f] = pspectrum(curr_v(i:(i+sliding_window-1),:),50);
    end
   
    pooled_f = [pooled_f, mean(pwr_f,2)];
    pooled_v = [pooled_v, mean(pwr_v,2)];

end


close all
get_fig_panel_power(pooled_v, f.', 1, 1, [0, 0.6], 'Membrane voltage')
get_fig_panel_power(pooled_f, f.', 4.5, 1, [0, 0.6], 'GCaMP')


%% Custom analysis functions used above

function [align_img, t_v, AP, b_w] = get_img_align(input_patch, img_data, fly_idx)
% This function aligns simultaneously-recorded imaging and patch-clamp data
% and prepares them for reverse correlation. Here, the imaging data (which
% were recorded at a substantially lower frame-rate than the patch-clamp
% data) are interpolated to 2 kHz for alignment, before both the imaging and
% patch-clamp are down-sampled later to 10 Hz in a different function
% ('get_reverse_correlation'). The reason for this interpolation step is to
% avoid inaccuracies in the TTL alignemnt and to allow for plotting
% aligned imaging and voltage traces simultaneously for quality checks.

%% Patch-clamp data

patch_data = input_patch(:,[1,2,end]);
b_w = mean(diff(patch_data(:,1)))/1000;

% Finds alignment index in patch-clamp data
patch_data(:,3) = (patch_data(:,3)-min(patch_data(:,3)))/max(patch_data(:,3)-min(patch_data(:,3)));
idx_1 = 1;
for loop_idx = 2:length(patch_data(:,3))
    if patch_data(loop_idx,3) < 0.15 && patch_data(loop_idx-1,3) >= 0.15
        patch_align_idx(idx_1) = loop_idx;
        idx_1 = idx_1+1;
    end
end

patch_align_idx = patch_align_idx(end);

%% Imaging data

% bin size imaging
bs_img = 1/58.26;

img_data(:,1) = img_data(:,1)+10000; % Just adds a large number to TTL signal to avoid negative values
% Finds alignment index in imaging data; for quality check only
img_data(:,1) = (img_data(:,1)-min(img_data(:,1)))/max(img_data(:,1)-min(img_data(:,1)));

idx_1 = 1;
for loop_idx = 2:length(img_data(:,1))
    if img_data(loop_idx,1) < 0.15 && img_data(loop_idx-1,1) >= 0.15
        img_align_idx_original(idx_1) = loop_idx;
        idx_1 = idx_1+1;
    end
end

img_align_idx_original = img_align_idx_original(end); % For quality check only

%% Computes dF/F
% Note: Here, a time-varying baseline fluorescence (F_0) is used. The
% 300-element moving average of the raw fluorescence trace is used for this
% purpose. This is neccessary here because some imaging traces showed slow
% decay of overal fluorescence (either bleaching or slow drift out of
% focus). This is fine in this particular case, because we are not
% interested in the actual dF/F values, but only in the
% time-course/waveform (we z-score the imaging data for model fitting).

F_0 = movmean(img_data(:,2),300);
img_data(:,2) = (img_data(:,2)-F_0)./F_0;

%% Interpolation of imaging data

time1 = [0:bs_img:(bs_img*(length(img_data)-1))]'; % Time vector for original imaging
time_2 = [0:b_w:time1(end,1)]'; % Time vector for patch-clamp and interpolation

% Interpolates imaging data to accomodate patch-clamp frame rate.
interp_img = interp1(time1,img_data,time_2);

%% Finds alignment index in INTERPOLATED imaging data for alignment below

idx_1 = 1;
for loop_idx = 2:length(interp_img(:,1))
    if interp_img(loop_idx,1) < 0.15 && interp_img(loop_idx-1,1) >= 0.15
        align_idx_img(idx_1) = loop_idx;
        idx_1 = idx_1+1;
    end
end

align_idx_img = align_idx_img(end);

%% ALIGNS IMAGING AND PATCH DATA

% Cuts the BEGINNING according to which of the traces is shorter
if align_idx_img > patch_align_idx
    interp_img(1:(align_idx_img-patch_align_idx),:) = [];
elseif align_idx_img < patch_align_idx
    patch_data(1:(patch_align_idx-align_idx_img),:) = [];
end

% cuts the END according to which of the traces is shorter
if length(interp_img) > length(patch_data)
    interp_img = interp_img(1:length(patch_data),:);
elseif length(interp_img) < length(patch_data)
    patch_data = patch_data(1:length(interp_img),:); 
end

% Aligned imaging and patch-clamp data. 
% First column: interpolated and pre-processed imaging data
% Second column: patch-clamp data
% Third column: TTL signal in imaging
% Fourth column: TTL signal in patch-clamp
% For quality check: plot third and fourth column together, they should
% overlay (both are the TTL signals).
align_img = [interp_img(:,2), patch_data(:,2), interp_img(:,1), patch_data(:,3)];

% Two cells had to be visually 'trouble-shooted' for which part of the traces were
% actually usable. Check by plotting the original data!
if fly_idx == 1
    align_img = align_img(1:150000,:);
elseif fly_idx == 2        
    align_img = align_img(25000:105000,:);
end

% Time vector for aligned data
t_v = time_2;
t_v = t_v(1:length(align_img));

%% Finds spike times and creates vector indicating spike times with '1'.
% This was viually checked for its accuracy to detect spikes. Double-check
% by simultaneously plotting the voltage trace of 'align_img' and AP.

[~, spktms_1] = findpeaks(zscore(movmean(diff(align_img(:,2)),10)),'MinPeakHeight',3);
AP = zeros(length(align_img(:,2)),1);
AP(spktms_1,1) = 1;

end

function rev_corr = get_IRF_estimate(curr_fly, filt_cond)
% This function estimates the linear filter / impulse response function
% (IRF) of the aligned and interpolated patch-clamp and imaging data. 
% 'filt_cond' specifies if the model is viewed from the 'filter perspective'
% or 'impulse response function / convolution perspective'. This just means
% two views of the same thing. The former just 'looks into the past' and
% the latter 'into the future' like in convolution.
% When you call this function in the masterscript, try replacing "impulse" 
% with "filter" and see what happens when you plot the figure.

%% Defines some of the variables needed below; just for simplifying the code below. 

AP = curr_fly.AP;
img_signal = curr_fly.align_img(:,1);
voltage_signal = curr_fly.align_img(:,2);
b_w = curr_fly.b_w;
t_v = curr_fly.t_v;

%% Computes PSTH and downsamples imaging trace accordingly.

s_f = 200; % Scaling factor for down-sampling/PSTH
% Computes PSTH
idx_1 = 1;
AP_nu = [AP; zeros(s_f,1)];
for idx_2 = 1:s_f:length(AP)
    PSTH(idx_1,1) = sum(AP_nu(idx_2:idx_2+(s_f-1),1));
    idx_1 = idx_1+1;  
end

% Down-samples imaging trace
calcium_response = movmean(img_signal(:,1),300); % Smoothing imaging trace
calcium_response = downsample(calcium_response,s_f);
calcium_response = calcium_response(1:length(PSTH),1);

%% Creates the predictor matrix

% Length of linear filter / IRF; given the current settings (2 kHz
% down-sampled to 10 Hz using a s_f of 200), a filter_length of 20 means a
% duration of 2 seconds.
filter_length = 20; 

% Pads early bins of PSTH with zeros
padded_PSTH = [zeros(filter_length-1,1); PSTH]; 
padded_PSTH = padded_PSTH(1:length(PSTH),1);

% Creates predictor matrix; 
% number of rows: length of PSTH
% number of columns: filter-length
% Using a 'for-loop', creates a matrix where each column is a 'shifted'
% version of the previous column.
for pm_idx = 1:filter_length
    pred_mat(:,pm_idx) = padded_PSTH;
    padded_PSTH = [padded_PSTH(2:end); 0]; % For each iteration (column of X), 'moves pSX one index up'.
end

% Depending on whether 'filter view' or 'iRF / convolution view' is used,
% predictor matrix either 'looks into past' or 'looks into future',
% respectively.
% To double-check: Plot the  predictor matrix (using the 'imagesc' function ) 
% of a fly when you use "filter" and "impulse" and compare the results.
if filt_cond == "filter"
    predictor_matrix = pred_mat;
elseif filt_cond == "impulse"
    predictor_matrix = flip(pred_mat,2);
end

%% Estimates linear filter / IRF by fitting linear model to data

mdl = fitlm(predictor_matrix,calcium_response); % Fits linear model
% Linear filter / IRF is coefficient weights without the y-intercept
linear_filter = mdl.Coefficients.Estimate(2:end,1);
% Estimates 'in-sample' fit of the imaging data (This is not prediction!)
% The actual out-of-sample model performance for each cell will be assessed
% in a different function ('get_rev_corr_pred').
img_fit = predict(mdl,predictor_matrix);

%% Saves some of the variables as output variables of the function

rev_corr.linear_filter = linear_filter; % Linear filter
rev_corr.img_fit = img_fit;
rev_corr.calcium_response = calcium_response;
rev_corr.inst_rate = PSTH;
rev_corr.voltage_signal = voltage_signal;
rev_corr.predictor_matrix = predictor_matrix;
rev_corr.t_v1 = (0:b_w*s_f:b_w*s_f*(length(img_fit)-1));
rev_corr.t_v = t_v;

% Specifies the time vector of the linear filter / IRF depending on whether
% filter version or IRF version is used.
if filt_cond == "filter"
    rev_corr.t_v2 = (-1*(b_w*s_f*(length(linear_filter)-1)):b_w*s_f:0);
elseif filt_cond == "impulse"
    rev_corr.t_v2 = (0:b_w*s_f:1*(b_w*s_f*(length(linear_filter)-1)));
end


end

function model_prediction = get_rev_corr_pred(input_fly, analysis_summary, fly_idx)
% This function assesses out-of-sample model performance. 
% Briefly, this is done by predicting the imaging trace of a given cell by
% using the mean coefficient weights of all the other cells.


% Constructs the linear filter / IRF for model prediction:
% For a particular cell, the mean of the linear filters / IRFs of all the other
% cells.
pred_filt = analysis_summary.linear_filter;
pred_filt(:,fly_idx) = []; % deletes the filter / IRF of current cell.
pred_filt = mean(pred_filt,2); % computes the mean of all the other cells.

% Predicts imaging trace by computing inner product between predictor matrix
% and linear filter / IRF.
model_prediction.pred_trace = input_fly.rev_corr.predictor_matrix*pred_filt;

% Computes cross-correlation between z-scored data and model and finds
% maximum correlation coefficient. This is done instead of simply
% computing the correlation coefficient to correct for small time-lag
% differences.
[norm_cross_corr, lag_vector_norm_cross_corr] = ...
    xcov(zscore(input_fly.rev_corr.calcium_response),zscore(model_prediction.pred_trace),'coef',60);
[max_corr, lag_idx] = max(norm_cross_corr);

% Finds the lag of the maximum correlation coefficient and check for
% quality control: lags should not be larger than one or two indices!
model_prediction.lags_norm_cross_corr = lag_vector_norm_cross_corr(lag_idx);

% For linear model, computes R2 by squaring correlation coefficient (corrected for lag, see above).
model_prediction.R2 = max_corr^2;

x_data = zscore(input_fly.rev_corr.calcium_response);
y_model = zscore(model_prediction.pred_trace);

% For quality check of 2-D histogram orientation (not to mix up x- and y-axes)
quality_check = "off";
if quality_check == "on"
    x_data(1:(length(x_data)/10)) = 0;
    y_model(1:(length(x_data)/10)) = 2;
end

% Computes bivariate histogram of z-scored data vs model.
[model_prediction.scatter_pred_hist, model_prediction.x_edges, model_prediction.y_edges] = ...
    histcounts2(x_data, y_model, 'NumBins', 15, 'XBinLimits', [-3, 3], 'YBinLimits', [-3, 3],...
    'Normalization', 'probability');

end

function transient_aligned_PSTH = get_transient_aligned_PSTH(curr_fly)
%% Defines some of the variables needed below; just for simplifying the code below. 

AP = curr_fly.AP;
img_signal = curr_fly.align_img(:,1);
b_w = curr_fly.b_w;
t_v = curr_fly.t_v;

%% Computes PSTH and downsamples imaging trace accordingly.

s_f = 20; % Scaling factor for down-sampling/PSTH
% Computes PSTH
idx_1 = 1;
AP_nu = [AP; zeros(s_f,1)];
for idx_2 = 1:s_f:length(AP)
    PSTH(idx_1,1) = sum(AP_nu(idx_2:idx_2+(s_f-1),1));
    idx_1 = idx_1+1;  
end

% Down-samples imaging trace
calcium_response = movmean(img_signal(:,1),500); % Smoothing imaging trace
calcium_response = downsample(calcium_response,s_f);
calcium_response = calcium_response(1:length(PSTH),1);

%% Saves the variables into struct.

transient_aligned_PSTH.AP = PSTH;
transient_aligned_PSTH.img_signal = calcium_response;
[transient_aligned_PSTH.transient_vec, ...
    transient_aligned_PSTH.cellbody_GCaMP_transients, transient_aligned_PSTH.aligned_PSTH] = ...
    get_cellbody_GCaMP_transients(transient_aligned_PSTH.img_signal, transient_aligned_PSTH.AP);
transient_aligned_PSTH.t_v1 = (0:b_w*s_f:b_w*s_f*(length(transient_aligned_PSTH.aligned_PSTH)-1));
transient_aligned_PSTH.t_v = t_v;

end

function [transient_vec, cellbody_GCaMP_transients, aligned_PSTH] = ...
    get_cellbody_GCaMP_transients(curr_img, AP)

gauss_sliding_window = 65;
min_peak_height = 0.00075;
min_dist = 45;

transient_vec = smoothdata(diff(curr_img),'gaussian',gauss_sliding_window);
[~, trans_ids] = findpeaks(transient_vec, 'MinPeakHeight',min_peak_height,...
    'MinPeakDistance',min_dist);


window_1 = 100;
window_2 = 150;
cum_idx = 1;
for trans_loop_idx = 1:length(trans_ids)
    if trans_ids(trans_loop_idx)-window_1 > 1 && ...
            trans_ids(trans_loop_idx)+window_2 <= length(curr_img)
        cellbody_GCaMP_transients(:,cum_idx) = ...
            curr_img(trans_ids(trans_loop_idx)-window_1:trans_ids(trans_loop_idx)+window_2,1);
        aligned_PSTH(:,cum_idx) = ...
            AP(trans_ids(trans_loop_idx)-window_1:trans_ids(trans_loop_idx)+window_2,1);
        cum_idx = cum_idx+1;
    end
end

plot_vec = zeros(length(curr_img),1);
plot_vec(trans_ids,1) = 1;

% hold off
% plot(plot_vec)
% hold on
% plot(curr_img)
% pause

end

%% Custom plotting functions used above

function get_fig_panel_reverse_correlation_example(rev_corr, model_prediction, analysis_summary, x_pos, y_pos, color, source_data_details)

sb_width = get_default_scale_bar_width;
ht_voltage = 0.75;
ht_img = 0.75;
ht_inset = 1;
sz_1 = 7;
zscr_lim = [-3, 3];
t_scale_bar_length = 10;
v_scale_bar_length = 15;
i_scale_bar_length = 2;
y_dist_magn = 0.3;
y_pos = y_pos + ht_img;
annotation_x_dist = 0.2;
x_dist_cartoon = 1.2;
x_dist_IRF = 2.3;
sz_magn = sz_1-x_dist_cartoon-x_dist_IRF;

ylm_voltage = [min(rev_corr.voltage_signal),max(rev_corr.voltage_signal)];
voltage_sb_origin = round(min(rev_corr.voltage_signal));

model_color = color.navy;

if voltage_sb_origin < ylm_voltage(1)
    disp("!!! Warning: y axis cale bar out of axis limits !!!")
end

%% full voltage trace

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_voltage]);
plot(pnl_1, rev_corr.voltage_signal, 'Color', [0, 0, 0], 'LineWidth', 0.25)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(ylm_voltage)
xlim([1, length(rev_corr.voltage_signal)])

% voltage scale bar
% CAUTION: this scale bar is an additional plot with its own axes properties
% and x/y-limits etc.! To use a scale bar like this, one MUST double-check
% that the x/y-limits and axes-dimensions etc. match the plot of the actual
% data!
% In this case, this means: double-check y-limits and size of y-axis.
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos, 0.1, ht_voltage]);
plot(pnl_1, [0, 0], [voltage_sb_origin, voltage_sb_origin+v_scale_bar_length], 'Color', [0, 0, 0],...
    'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(ylm_voltage)
xlim([0, 1])

%% full imaging

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-ht_img-0.1, sz_1, ht_img]);
plot(pnl_1, rev_corr.t_v1, zscore(rev_corr.calcium_response), 'Color', color.medium_gray)
hold on
plot(pnl_1, rev_corr.t_v1, zscore(model_prediction.pred_trace), 'Color', model_color)
plot(pnl_1, [rev_corr.t_v1(end)-t_scale_bar_length, rev_corr.t_v1(end)], [zscr_lim(1), zscr_lim(1)], 'k',...
    'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(zscr_lim)
xlim([0, rev_corr.t_v1(end)])

% imaging scale bar
% voltage scale bar
% CAUTION: this scale bar is an additional plot with its own axes properties
% and x/y-limits etc.! To use a scale bar like this, one MUST double-check
% that the x/y-limits and axes-dimensions etc. match the plot of the actual
% data!
% In this case, this means: double-check y-limits and size of y-axis.
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos-ht_img-0.1, 0.1, ht_img]);
plot(pnl_1, [0, 0], [zscr_lim(1), zscr_lim(1)+i_scale_bar_length], 'Color', model_color, 'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(zscr_lim)
xlim([0, 1])

%% magnified imaging
idx_1 = 85.5;    
idx_2 = 95.5;

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+x_dist_cartoon, y_pos+y_dist_magn+ht_voltage, ...
    sz_magn, ht_inset]);
plot(pnl_1, rev_corr.t_v1, zscore(rev_corr.calcium_response), 'Color', color.medium_gray, 'LineWidth', 1.5)
hold on
plot(pnl_1, rev_corr.t_v1, zscore(model_prediction.pred_trace), 'Color', model_color, 'LineWidth', 1.5)
plot(pnl_1, [rev_corr.t_v1(end)-t_scale_bar_length, rev_corr.t_v1(end)], [zscr_lim(1), zscr_lim(1)], 'k',...
    'LineWidth', 1)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(zscr_lim)
xlim([idx_1, idx_2])
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';

% imaging scale bar
% voltage scale bar
% CAUTION: this scale bar is an additional plot with its own axes properties
% and x/y-limits etc.! To use a scale bar like this, one MUST double-check
% that the x/y-limits and axes-dimensions etc. match the plot of the actual
% data!
% In this case, this means: double-check y-limits and size of y-axis.
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_magn+x_dist_cartoon+0.1, y_pos+y_dist_magn+ht_voltage, ...
    sz_magn, ht_inset]);
plot(pnl_1, [0, 0], [zscr_lim(2)-i_scale_bar_length, zscr_lim(2)], 'Color', model_color,...
    'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(zscr_lim)
xlim([0, 1])

%% magnified voltage
magn_x_sb_origin = idx_2-1;
magn_x_sb_length = 1;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+x_dist_cartoon, y_pos+y_dist_magn+ht_voltage, ...
    sz_magn, ht_inset]);
plot(pnl_1, rev_corr.t_v, rev_corr.voltage_signal, 'Color', [0, 0, 0], 'LineWidth', 0.25)
hold on
plot(pnl_1, [magn_x_sb_origin, magn_x_sb_origin+magn_x_sb_length],...
    [ylm_voltage(1), ylm_voltage(1)], 'k',...
    'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';

ylim(ylm_voltage)
xlim([idx_1, idx_2]) 
    
% voltage scale bar
% voltage scale bar
% CAUTION: this scale bar is an additional plot with its own axes properties
% and x/y-limits etc.! To use a scale bar like this, one MUST double-check
% that the x/y-limits and axes-dimensions etc. match the plot of the actual
% data!
% In this case, this means: double-check y-limits and size of y-axis.
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_magn+x_dist_cartoon+0.1, y_pos+y_dist_magn+ht_voltage, ...
    sz_magn, ht_inset]);
plot(pnl_1, [0, 0], [voltage_sb_origin, voltage_sb_origin+v_scale_bar_length], 'Color', [0, 0, 0],...
    'LineWidth', sb_width)
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.Color = 'none';
ylim(ylm_voltage)
xlim([0, 1])

%% Adds red lines to signify magnified area in full trace

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-ht_img-0.1, sz_1, ht_img+ht_voltage+0.2]);
plot(pnl_1, [idx_1, idx_1], [ht_img+ht_voltage+0.1, 0],...
    'Color', color.salmon, 'LineWidth', 0.5)
hold on
plot(pnl_1, [idx_2, idx_2], [0, ht_img+ht_voltage+0.1],...
    'Color', color.salmon, 'LineWidth', 0.5)
pnl_1.Color = 'none';
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
ylim([0, ht_img+ht_voltage+0.2])
xlim([0, rev_corr.t_v1(end)])

%% Adds annotations
get_default_annotation(x_pos+x_dist_cartoon, y_pos+y_dist_magn+ht_voltage+ht_inset+0.3, 'Membrane potential', [0, 0, 0], 'normal', "left")
get_default_annotation(x_pos+x_dist_cartoon+2.05, y_pos+y_dist_magn+ht_voltage+ht_inset+0.3, 'GCaMP', color.gray, 'normal', "left")
get_default_annotation(x_pos+x_dist_cartoon+sz_magn, y_pos+y_dist_magn+ht_voltage+ht_inset+0.3, 'model', model_color, 'normal', "right")
get_default_annotation(x_pos+sz_1, y_pos-ht_img-0.2, '10 s', 'k', 'normal', "right")
get_default_annotation(x_pos+sz_1+annotation_x_dist, y_pos+0.2, '15 mV', 'k', 'normal', "left")
get_default_annotation(x_pos+sz_1+annotation_x_dist, y_pos-ht_img+0.1, '2 s.d.', 'k', 'normal', "left")
get_default_annotation(x_pos+sz_magn+x_dist_cartoon,...
    y_pos+y_dist_magn+ht_voltage-0.1, '1 s', 'k', 'normal', "right")

%% Adds sub panel with GCaMP IRF

sz_IRF = 1;
get_fig_panel_GCaMP_IRF(rev_corr, analysis_summary.linear_filter, ...
    x_pos+sz_1-sz_IRF, y_pos+y_dist_magn+ht_voltage, sz_IRF, color, source_data_details);

%% Saves the source data

voltage_trace = array2table(rev_corr.voltage_signal);
model_prediction = array2table(zscore(model_prediction.pred_trace));
img_data = array2table(zscore(rev_corr.calcium_response));
writetable(voltage_trace, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'voltage_trace')
writetable(img_data, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'img_data')
writetable(model_prediction, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'model_prediction')

end

function get_fig_panel_GCaMP_IRF(rev_corr, linear_filter, x_pos, y_pos, sz_1, color, source_data_details)

sb_width = get_default_scale_bar_width;
ht_IRF = 1;
ylm_filter = [-0.6, 1.2];
xlm_1 = ([rev_corr.t_v2(1), rev_corr.t_v2(end)]);

model_color = color.navy;

% Plot the mean IRF.
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_IRF]);
hold on
get_default_SEM_area_plot(curr_pnl, linear_filter, rev_corr.t_v2, model_color)
plot(curr_pnl, [rev_corr.t_v2(end-10), rev_corr.t_v2(end-5)], [ylm_filter(1), ylm_filter(1)], 'k', 'LineWidth', sb_width)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), [], [], ylm_filter(1), ylm_filter(2),...
    ylm_filter(1), ylm_filter(2), 0.5, 0.6, "linear", "linear",[],'Norm. weight',...
    'none','k','none','k', sz_1, ht_IRF)

%% Some annotations

get_default_annotation(x_pos+sz_1-0.1, y_pos-0.1, ...
    '500 ms', 'k', 'normal', "right")

%% Saves the source data

GCaMP_IRF = array2table(linear_filter);
writetable(GCaMP_IRF, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'GCaMP_IRF')

end

function get_fig_panel_data_vs_model(model_prediction, scatter_hist_all, x_pos, y_pos, source_data_details)


s_z = 1.8;
f_s = get_default_font_size;

xy_vec = model_prediction.x_edges(2:end)-(mean(diff(model_prediction.x_edges))/2); 
xlm_1 = [-3, 3];
ylm_1 = [-3, 3];

col_1(:,1) = linspace(1,0,10000);
col_1(:,2) = linspace(1,0,10000);
col_1(:,3) = linspace(1,0,10000);
curr_colormap = [col_1(:,1), col_1(:,2), col_1(:,3)];

%% Plots data (using same axis limits etc as above!!!)

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, s_z, s_z]);
hold on
m_scat_hist = mean(scatter_hist_all,3).*100;
imagesc(curr_pnl, xy_vec, xy_vec, m_scat_hist.') % Note the transpose
curr_pnl.Color = 'none';
colormap(curr_pnl,curr_colormap)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), 1.5, 1.5, "linear", "linear", 'Data (z-score)', 'Model (z-score)',...
    'k', 'k', 'k', 'k', s_z, s_z)

%% Adds colorbar

clim1 = [0,max(max(mean(scatter_hist_all,3))).*100];
clrbr = colorbar(curr_pnl,'eastoutside');
clrbr.Units = 'centimeters';
clrbr.Position = [x_pos+s_z+0.1, y_pos,  0.2, s_z];
ylabel(clrbr,'Percent')
clrbr.FontSize = f_s;
clrbr.Color = 'k';
clrbr.TickLength = 0;
clrbr.Limits = clim1;
clrbr.YLabel.Visible = 'on';
clrbr.YLabel.Color = 'k';

%% Adds annotation

get_default_annotation(x_pos+s_z, y_pos+0.2, '10 cells', 'k', 'normal', "right")


%% Saves the source data


data_vs_model = table;
    data_vs_model.data_vec = [NaN; xy_vec'];
    data_vs_model.model_vec = [xy_vec; m_scat_hist];
    writetable(data_vs_model, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'data_vs_model')

end

function get_fig_panel_GCaMP_aligned_PSTH(analysis_summary, x_pos, y_pos, color, source_data_details)

sb_width = get_default_scale_bar_width;
sz_1 = 1.8;
ht_1 = 1.8;
ylm_img = [-0.4, 1.2];
ylm_psth = [0, 30];
xlm_1 = [-0.5, 1.5];

scale_fac = 100;
curr_psth = analysis_summary.aligned_PSTH*scale_fac;
curr_img = analysis_summary.curr_diff_img;
t_v = analysis_summary.aligned_tv;

PETH_color = [0, 0, 0];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, movmean(curr_psth,3), t_v, PETH_color)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_psth(1), ylm_psth(2),...
    ylm_psth(1), ylm_psth(2), 0.5, 10,"linear", "linear", 'Time from transient onset (s)', 'Firing rate (Hz)',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

img_color = color.medium_gray;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, curr_img, t_v, img_color)
scale_bar_cond = 0;
if scale_bar_cond == 1
    plot(curr_pnl, [xlm_1(2), xlm_1(2)],...
        [0, 0.5], 'LineWidth', sb_width, 'Color', img_color)
    plot(curr_pnl, [0, 0], [ylm_img(1), ylm_img(2)], 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':')
    get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_img(1), ylm_img(2),...
        ylm_img(1), ylm_img(2), [], [], "linear", "linear", [], [],...
        'none', 'none', 'none', 'none', sz_1, ht_1)
    get_default_annotation(x_pos+sz_1+0.2, y_pos+(ht_1/2), '0.5 a.U.', 'k', 'normal', "left")

else
    plot(curr_pnl, [0, 0], [ylm_img(1), ylm_img(2)], 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':')
    curr_pnl.YAxisLocation = "right";
    get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_img(1), ylm_img(2),...
        ylm_img(1), ylm_img(2), [], 0.4, "linear", "linear", [],'Normalized dF/dt',...
        'none', 'k', 'none', color.medium_gray, sz_1, ht_1)
    
end



%% Saves the source data

trans_aligned_PETH = array2table(curr_psth(2:end,:));
trans_aligned_imaging = array2table(curr_img(2:end,:)); % cut first row because trace is 'diff';
writetable(trans_aligned_PETH, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'trans_aligned_PETH')
writetable(trans_aligned_imaging, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'trans_aligned_imaging')

end

function get_fig_panel_power(curr_power, f, x_pos, y_pos, ylm_1, curr_annotation)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_1 = [0, 2];

plotting_color = [0, 0, 0];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, curr_power, f, plotting_color)
get_default_separated_ax(curr_pnl, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), 0.5, 0.2, "linear", "linear", 'Frequency (Hz)', 'Power',...
    'k', 'k', 'k', 'k', sz_1, ht_1)
get_default_annotation(x_pos+sz_1, y_pos+ht_1, curr_annotation, 'k', 'normal', "right")

end

