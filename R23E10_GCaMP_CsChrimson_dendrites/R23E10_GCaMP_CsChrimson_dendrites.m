% Analysis of dendritic dFBN GCaMP signals during optogenetic stimulation.
% Code written by Peter Hasenhuetl.

clear all
tic

% Defines the path for saving the source data for figures
source_data_details_GCaMP_Opto.data_path = []; %Add path as character array
source_data_details_GCaMP_Opto.file_name = 'R23E10_GCaMP_dendrites_Opto_data.xlsx';

% Goes into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
GCaMP_OPTO_data_table = [];

% Loops through individual flies
for fly_idx = 1:n_flies
     
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        fly_details.fly_name = char(fly_input(fly_idx).name);
        if char(fly_input(fly_idx).name) == "fly20260116_7DFB_SWAOpto.mat"
            loaded_fly.TTL(3093,1) = 0; % removing an artefact in TTL trace
        end

        GCaMP_OPTO_data_table = [GCaMP_OPTO_data_table; ...
            get_individual_GCaMP_opto_fly(loaded_fly, fly_details, fly_idx)];
end

cd([]) %Add path as character array

GCaMP_summary = get_GCaMP_summary(GCaMP_OPTO_data_table);

toc

%% Plotting the results, comprehensive version

close all
figure('Name','GCaMP Opto figure','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 15],'Resize','off')
color = get_color;
plotting_color_left = [0, 0, 0];
plotting_color_right = color.medium_gray;

get_fig_panel_GCaMP_opto_traces(0.9, 4.5, GCaMP_OPTO_data_table, plotting_color_left, ...
    plotting_color_right)
get_fig_panel_transient_GCaMP_opto(1, 2, GCaMP_summary, ...
    plotting_color_left, plotting_color_right)
get_fig_panel_auto_corr_GCaMP_opto(4, 2, GCaMP_summary)
get_fig_panel_xcorr_GCaMP_opto(7, 2, GCaMP_summary)

cd(source_data_details_GCaMP_Opto.data_path)
set(gcf,'renderer','painters')
saveas(gcf,'R23E10_GCaMP_opto_fig_v1.pdf')
cd([]) %Add path as character array



%% Plotting the results, concise version

close all
figure('Name','GCaMP Opto figure','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 15],'Resize','off')
color = get_color;
plotting_color_left = [0, 0, 0];
plotting_color_right = color.medium_gray;
get_fig_panel_pwr_GCaMP_opto_v2(7.2, 2.1, GCaMP_summary, ...
    plotting_color_left, plotting_color_right)
example_traces = get_fig_panel_GCaMP_opto_traces_v2(0.7, 2, GCaMP_OPTO_data_table, ...
    plotting_color_left, plotting_color_right, 45, 75);

cd(source_data_details_GCaMP_Opto.data_path)
set(gcf,'renderer','painters')
saveas(gcf,'R23E10_GCaMP_opto_fig_v2.pdf')
cd([]) %Add path as character array

writetable(example_traces, [source_data_details_GCaMP_Opto.data_path, source_data_details_GCaMP_Opto.file_name],...
        'Sheet', 'GCaMP traces')

%% Custom functions called in this script

function ind_data_table = get_individual_GCaMP_opto_fly(input_fly, fly_details, fly_idx)

%% Normalizes TTL trace to be between 0 and 1

TTL_input = input_fly.TTL;
nTTL = TTL_input(:,1)+5000;
nTTL = nTTL-min(nTTL);
nTTL = nTTL./max(nTTL);

%% Generates vector with non-zero elements only at indices of TTL onset

TTLidx = zeros(length(nTTL),1);

for loop_idx1 = 1:length(nTTL)
    if nTTL(loop_idx1,1) < 0.3
        binarized_TTL(loop_idx1,1) = 0;
    else
        binarized_TTL(loop_idx1,1) = 1;
    end
end

for loop_idx2 = 2:(numel(binarized_TTL)-1)
    if binarized_TTL(loop_idx2,1) == 1 && binarized_TTL(loop_idx2-1,1) == 0
        TTLidx(loop_idx2,1) = max(nTTL(loop_idx2:loop_idx2+1,1));
    end
end

TTLoutputs_raw = find(TTLidx > 0);
inter_burst_interval = 1.4*get_GCaMP_opto_frame_rate; % 1.4 seconds, instead of 1.5 to have a tolerance window
TTLoutputs = TTLoutputs_raw;
TTLoutputs([inter_burst_interval*2; diff(TTLoutputs)] < 10) = [];


%% Analyses the imaging traces

curr_img = input_fly.img_traces;
img_output = get_GCaMP_opto_img_measures(curr_img, TTL_input, TTLoutputs);
st_img_left = img_output.st_img_left;
st_img_right = img_output.st_img_right;
img_traces = img_output.dF_aligned;
img_traces_full = [img_output.TTL_aligned, img_traces];
st_TTL = img_output.st_TTL;

ind_data_table = get_GCaMP_opto_data_table_row(img_output, st_img_left, st_img_right, ...
    st_TTL, img_traces_full, img_output.pwr_spctrm, img_output.f, fly_details, fly_idx);

end

function img_output = get_GCaMP_opto_img_measures(raw_img, TTL_input, TTLoutputs)

f_r = get_GCaMP_opto_frame_rate; % frame-rate
img_trace = raw_img(:,2:3) - raw_img(:,1); % element-wise subtraction of background fluorescence 
pre_LED_window = round(0.5*(f_r));
post_LED_window = round(1.5*(f_r));
pre_LED_window_full = 200;
post_LED_window_full = round(201.5*get_GCaMP_opto_frame_rate);
rs_1 = TTLoutputs-pre_LED_window;
rs_2 = TTLoutputs+post_LED_window;
rs_full = TTLoutputs(1)-pre_LED_window_full;
TTL_aligned_full = TTLoutputs-rs_full;
img_output = struct;

%% Then, computes dF/F

F0_mov = 7000;
curr_raw = movmean(img_trace(:,1:2),round(f_r/20));
d_F = curr_raw;
img_output.dF = d_F;
d_F = (d_F - movmean(d_F,F0_mov))./(movmean(d_F,F0_mov));
img_output.dF_aligned = d_F(rs_full:rs_full+post_LED_window_full,:);
img_output.TTL_aligned = TTL_input(rs_full:rs_full+post_LED_window_full,1);
img_output.raw_img = img_trace;
img_output.TTL_times_full = TTL_aligned_full;

%% Computes power spectrum

[pwr_spctrm, f] = pspectrum(img_output.dF_aligned, f_r);

%% Computes auto_correlation

[auto_corr1, ~] = xcov(img_output.dF_aligned(:,1),'coef');
[auto_corr2, lag_auto_corr] = xcov(img_output.dF_aligned(:,1),'coef');
img_output.auto_corr = mean([auto_corr1, auto_corr2], 2);
img_output.lag_auto_corr = lag_auto_corr./get_GCaMP_opto_frame_rate;

%% Computes cross_correlation

[img_output.x_corr, lag_xcorr] = xcov(img_output.dF_aligned(:,1), img_output.dF_aligned(:,2), 'coef');
img_output.lag_x_corr = lag_xcorr./get_GCaMP_opto_frame_rate;


%% ...LED-aligned img

for trial_idx = 1:size(rs_1,1)
    F_0_left = mean(img_trace(rs_1(trial_idx,1):TTLoutputs(trial_idx),1),1);
    F_0_right = mean(img_trace(rs_1(trial_idx,1):TTLoutputs(trial_idx),2),1);
    st_img_left(:,trial_idx) = (img_trace(rs_1(trial_idx,1):rs_2(trial_idx,1),1)-F_0_left)/abs(F_0_left);
    st_img_right(:,trial_idx) = (img_trace(rs_1(trial_idx,1):rs_2(trial_idx,1),2)-F_0_right)/abs(F_0_right);
    st_TTL(:,trial_idx) = TTL_input(rs_1(trial_idx,1):rs_2(trial_idx,1),1);
end

st_tv = ((1:size(st_img_right,1))-(pre_LED_window))/f_r;

%% Stores data in structure

img_output.st_img_left = st_img_left;
img_output.st_img_right = st_img_right;
img_output.st_TTL = st_TTL;
img_output.pwr_spctrm = pwr_spctrm;
img_output.f = f.';
img_output.st_tv = st_tv;

end

function dm_row = get_GCaMP_opto_data_table_row(dm_input, st_img_left, st_img_right, st_TTL, ...
    img_traces_full, pwr_spctrm, f, details, fly_idx)

dm_details_row = table;
dm_details_row.fly_id = fly_idx;
det_fl = struct;
det_fl.fly_name = details.fly_name;
dm_details_row.flyname = det_fl;

dm_img_row = table;

TTL_strct = struct;
TTL_strct.TTL_times_full = dm_input.TTL_times_full;
dm_img_row.TTL_times = TTL_strct;

trc_strct = struct;
trc_strct.st_img_left = st_img_left;
dm_img_row.st_img_left = trc_strct;

trc_strct = struct;
trc_strct.st_img_right = st_img_right;
dm_img_row.st_img_right = trc_strct;

trc_strct = struct;
trc_strct.st_tv = dm_input.st_tv;
dm_img_row.st_tv = trc_strct;

st_TTL_strct = struct;
st_TTL_strct.TTL_trace = st_TTL;
dm_img_row.st_TTL = st_TTL_strct;

trc_strct = struct;
trc_strct.traces_aligned = img_traces_full;
dm_img_row.img_traces_full = trc_strct;

trc_strct = struct;
trc_strct.pwr_spctrm = pwr_spctrm;
dm_img_row.pwr_spctrm = trc_strct;

trc_strct = struct;
trc_strct.f = f;
dm_img_row.f = trc_strct;

trc_strct = struct;
trc_strct.x_corr = dm_input.x_corr;
dm_img_row.x_corr = trc_strct;

trc_strct = struct;
trc_strct.lag_x_corr = dm_input.lag_x_corr;
dm_img_row.lag_x_corr = trc_strct;

trc_strct = struct;
trc_strct.auto_corr = dm_input.auto_corr;
dm_img_row.auto_corr = trc_strct;

trc_strct = struct;
trc_strct.lag_auto_corr = dm_input.lag_auto_corr;
dm_img_row.lag_auto_corr = trc_strct;

dm_row = [dm_details_row, dm_img_row];

end

function GCaMP_summary = get_GCaMP_summary(GCaMP_OPTO_data_table)

for loop_idx = 1:size(GCaMP_OPTO_data_table, 1)

    GCaMP_summary.ST_left(:,loop_idx) = mean(GCaMP_OPTO_data_table.st_img_left(loop_idx,1).st_img_left,2);
    GCaMP_summary.ST_right(:,loop_idx) = mean(GCaMP_OPTO_data_table.st_img_right(loop_idx,1).st_img_right,2);
    GCaMP_summary.ST_TTL(:,loop_idx) = mean(GCaMP_OPTO_data_table.st_TTL(loop_idx,1).TTL_trace,2);
    GCaMP_summary.P_left(:,loop_idx) = GCaMP_OPTO_data_table.pwr_spctrm(loop_idx,1).pwr_spctrm(:,1);
    GCaMP_summary.P_right(:,loop_idx) = GCaMP_OPTO_data_table.pwr_spctrm(loop_idx,1).pwr_spctrm(:,2);
    GCaMP_summary.X_CORR(:,loop_idx) = GCaMP_OPTO_data_table.x_corr(loop_idx,1).x_corr;
    GCaMP_summary.AUTO_CORR(:,loop_idx) = GCaMP_OPTO_data_table.auto_corr(loop_idx,1).auto_corr;

end

GCaMP_summary.lag_x_corr = GCaMP_OPTO_data_table.lag_x_corr(1,1).lag_x_corr;
GCaMP_summary.lag_auto_corr = GCaMP_OPTO_data_table.lag_auto_corr(1,1).lag_auto_corr;
GCaMP_summary.f = GCaMP_OPTO_data_table.f(1,1).f;
GCaMP_summary.st_tv = GCaMP_OPTO_data_table.st_tv(1,1).st_tv;

end

function get_fig_panel_transient_GCaMP_opto(x_pos, y_pos, GCaMP_summary, ...
    plotting_color_left, plotting_color_right)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_transient = [-0.5, 1.5];
ylm_transient = [-0.03, 0.3];
color = get_color;
color_stim = color.red;
[~, stim_idx] = findpeaks(mean(GCaMP_summary.ST_TTL,2), MinPeakHeight=10);

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
scatter(curr_pnl, GCaMP_summary.st_tv(stim_idx), ylm_transient(2).*ones(length(stim_idx),1), 2.5, 'Marker', '|',...
    'MarkerEdgeColor', color_stim, 'MarkerFaceColor', color_stim, 'LineWidth', 0.25)
hold on
get_default_SEM_area_plot(curr_pnl, GCaMP_summary.ST_right, GCaMP_summary.st_tv, plotting_color_right)
get_default_SEM_area_plot(curr_pnl, GCaMP_summary.ST_left, GCaMP_summary.st_tv, plotting_color_left)
get_default_separated_ax(curr_pnl, xlm_transient(1), xlm_transient(2), xlm_transient(1), xlm_transient(2), ...
    ylm_transient(1), ylm_transient(2), 0, ylm_transient(2), 0.5, 0.1, "linear", "linear", 'Time (s)', '\DeltaF/F',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

n_flies = size(GCaMP_summary.ST_left,2);
get_default_annotation(x_pos+sz_1, y_pos+ht_1, [num2str(n_flies), ' flies'], 'k', 'normal', "right")


end

function get_fig_panel_xcorr_GCaMP_opto(x_pos, y_pos, GCaMP_summary)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_xcorr = [-3, 3];
ylm_xcorr = [-0.5, 1];
plotting_color_xcorr = [0, 0, 0];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, GCaMP_summary.X_CORR, GCaMP_summary.lag_x_corr, plotting_color_xcorr)
get_default_separated_ax(curr_pnl, xlm_xcorr(1), xlm_xcorr(2), xlm_xcorr(1), xlm_xcorr(2), ...
    ylm_xcorr(1), ylm_xcorr(2), ylm_xcorr(1), ylm_xcorr(2),...
    1, 0.5, "linear", "linear", 'Lag (s)', 'Cross-correlation',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

end

function get_fig_panel_auto_corr_GCaMP_opto(x_pos, y_pos, GCaMP_summary)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_auto_corr = [-3, 3];
ylm_auto_corr = [-0.5, 1];
plotting_color_xcorr = [0, 0, 0];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, GCaMP_summary.AUTO_CORR, GCaMP_summary.lag_auto_corr, plotting_color_xcorr)
get_default_separated_ax(curr_pnl, xlm_auto_corr(1), xlm_auto_corr(2), xlm_auto_corr(1), xlm_auto_corr(2), ...
    ylm_auto_corr(1), ylm_auto_corr(2), ylm_auto_corr(1), ylm_auto_corr(2),...
    1, 0.5, "linear", "linear", 'Lag (s)', 'Auto-correlation',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

end

function get_fig_panel_GCaMP_opto_traces(x_pos, y_pos, GCaMP_OPTO_data_table, ...
    color_left, color_right)

sz_1 = 7.9;
ht_1 = 1.2;
y_dis = 0.2;
f_r = get_GCaMP_opto_frame_rate;
ylim_left = [-0.3, 0.8];
ylim_right = [-0.1, 0.3];
ylm_stim = [0, 1];
line_width_plot = 0.25;
color = get_color;
color_LED_bar = color.red;

curr_fly = 4;
curr_trace_GCaMP_left = GCaMP_OPTO_data_table.img_traces_full(curr_fly,1).traces_aligned(:,2);
curr_trace_GCaMP_right = GCaMP_OPTO_data_table.img_traces_full(curr_fly,1).traces_aligned(:,3);

stim_array = zeros(length(curr_trace_GCaMP_right),1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1+(ht_1+y_dis)]);
hold on
for loop_idx = 1:100
    light_on = GCaMP_OPTO_data_table.TTL_times(curr_fly,1).TTL_times_full(loop_idx);
    light_off = light_on+(f_r*0.5);
    stim_array(light_on:round(light_off),1) = 1;
    fill(curr_pnl,[light_on, light_off, light_off, light_on],...
        [ylm_stim(1), ylm_stim(1), ylm_stim(2), ylm_stim(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
    plot(curr_pnl, [light_on, light_off],...
        [ylm_stim(2), ylm_stim(2)], 'LineWidth', 1, 'Color', color_LED_bar)
end
get_default_separated_ax(curr_pnl, 1, length(curr_trace_GCaMP_right), [], [], ylm_stim(1), ylm_stim(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1+ht_1+y_dis)

curr_pnl1 = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos+ht_1+y_dis, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_GCaMP_left, f_r, ...
    [], [], color_left, line_width_plot,  "with_y_axis", ...
    ylim_left, [0, ylim_left(2)], 0.4, [], 0, "\DeltaF/F", []);

curr_pnl1 = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_GCaMP_right, f_r, ...
    10, [], color_right, line_width_plot,  "with_y_axis", ...
    ylim_right, [0, ylim_right(2)], 0.1, "10 s", 0, "\DeltaF/F", []);

Y_pos_annotation_left = y_pos+ht_1+y_dis+ht_1-0.1;
get_default_annotation(x_pos+sz_1, Y_pos_annotation_left, 'Left hemisphere', color_left, 'normal', "right")
Y_pos_annotation_right = y_pos+ht_1-0.1;
get_default_annotation(x_pos+sz_1, Y_pos_annotation_right, 'Right hemisphere', color_right, 'normal', "right")

end

function f_r = get_GCaMP_opto_frame_rate

f_r = 136.63;

end

function get_fig_panel_pwr_GCaMP_opto_v2(x_pos, y_pos, GCaMP_summary, ...
    plotting_color_left, plotting_color_right)

sz_1 = 1.6;
ht_1 = 1.8;
xlm_fft = [0, 2];
scale_fac = 10^3;
ylm_fft = scale_fac.*[0, 0.0074];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, scale_fac*GCaMP_summary.P_left, GCaMP_summary.f, plotting_color_left)
get_default_SEM_area_plot(curr_pnl, scale_fac*GCaMP_summary.P_right, GCaMP_summary.f, plotting_color_right)
get_default_separated_ax(curr_pnl, xlm_fft(1), xlm_fft(2), xlm_fft(1), xlm_fft(2), ...
    ylm_fft(1), ylm_fft(2), ylm_fft(1), ylm_fft(2),...
    0.5, 3.5, "linear", "linear", 'Frequency (Hz)', {'(\DeltaF/F)^2 x10^{-3}'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)

end

function example_traces = get_fig_panel_GCaMP_opto_traces_v2(x_pos, y_pos, GCaMP_OPTO_data_table, ...
    color_left, color_right, selection_TTL1, selection_TTL2)

sz_1 = 5.5;
ht_1 = 0.9;
y_dis = 0.2;
f_r = get_GCaMP_opto_frame_rate;
ylim_left = [-0.2, 0.3];
ylim_right = [-0.1, 0.1];
ylm_stim = [0, 1];
line_width_plot = 0.25;
color = get_color;
color_LED_bar = color.red;

curr_fly = 4;
curr_TTL = GCaMP_OPTO_data_table.TTL_times(curr_fly,1).TTL_times_full;
curr_range = curr_TTL(selection_TTL1:selection_TTL2);
curr_trace_GCaMP_left = GCaMP_OPTO_data_table.img_traces_full(curr_fly,1).traces_aligned(curr_range(1):curr_range(end),2);
curr_trace_GCaMP_right = GCaMP_OPTO_data_table.img_traces_full(curr_fly,1).traces_aligned(curr_range(1):curr_range(end),3);
curr_range = curr_range-(curr_range(1))+1;


xlim_1 = [1, length(curr_trace_GCaMP_right)];

stim_array = zeros(length(curr_trace_GCaMP_right),1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1+(ht_1+y_dis)]);
hold on
for loop_idx = 1:length(curr_range)-1
    light_on = curr_range(loop_idx);
    light_off = light_on+(f_r*0.5);
    stim_array(light_on:round(light_off),1) = 1;
    fill(curr_pnl,[light_on, light_off, light_off, light_on],...
        [ylm_stim(1), ylm_stim(1), ylm_stim(2), ylm_stim(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
    plot(curr_pnl, [light_on, light_off],...
        [ylm_stim(2), ylm_stim(2)], 'LineWidth', 1, 'Color', color_LED_bar)
end
get_default_separated_ax(curr_pnl, xlim_1(1), xlim_1(2), [], [], ylm_stim(1), ylm_stim(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1+ht_1+y_dis)

curr_pnl1 = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos+ht_1+y_dis, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_GCaMP_left, f_r, ...
    [], [], color_left, line_width_plot,  "with_y_axis", ...
    ylim_left, [0, ylim_left(2)], 0.3, [], 0, "\DeltaF/F", []);

curr_pnl1 = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_GCaMP_right, f_r, ...
    2, [], color_right, line_width_plot,  "with_y_axis", ...
    ylim_right, [ylim_right(1), ylim_right(2)], 0.1, "2 s", 0, "\DeltaF/F", []);

example_traces = table;
example_traces.stim_array = stim_array(1:length(curr_trace_GCaMP_left),1);
example_traces.GCaMP_left = curr_trace_GCaMP_left;
example_traces.GCaMP_right = curr_trace_GCaMP_right;

end
