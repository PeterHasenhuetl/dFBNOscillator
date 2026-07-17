% Analysis of R23E10-GAL4-driven GluSnFR-, GRABAch-, or AchSnFR-signals
% upon optogenetic stimulation of dFBNs.
% Code written by Peter Hasenhuetl.

clear all
tic

% Defines the path for saving the source data for figures
source_data_details_Glu_ACh.data_path = []; %Add path as character array
source_data_details_Glu_ACh.file_name = 'R23E10_Glu_ACh_opto_data.xlsx';

% Goes into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
Ach_OPTO_data_table = [];

% Loops through individual flies
for fly_idx = 1:n_flies
    cd([]) %Add path as character array
    fly_input = dir('*fly*.mat');
    loaded_fly = load(char(fly_input(fly_idx).name));
    fly_details.fly_name = char(fly_input(fly_idx).name);
    Ach_OPTO_data_table = [Ach_OPTO_data_table; ...
        get_individual_Ach_opto_fly(loaded_fly, fly_details, fly_idx)];
end

cd([]) %Add path as character array
Ach_OPTO_data_table(13,:) = []; % excluded because of light artifacts
SNFR_Ach_summary = get_group_summaries(Ach_OPTO_data_table, "AchSnFR");
SNFR_Glu_summary = get_group_summaries(Ach_OPTO_data_table, "GluSnFR");
GRAB_Ach_summary = get_group_summaries(Ach_OPTO_data_table, "GRABAch");

toc
 
%% Plotting the results

close all
figure('Name','NT imaging','Color','white',...
    'Units','centimeters','Position',[10 12 18.3 15],'Resize','off')
color = get_color;
get_fig_panel_Ach_opto_traces(2, 2, Ach_OPTO_data_table, source_data_details_Glu_ACh)
get_fig_panel_aligned_Ach_opto(SNFR_Ach_summary, SNFR_Glu_summary, GRAB_Ach_summary, 7, 2, ...
    color, "with_map", source_data_details_Glu_ACh)

cd(source_data_details_Glu_ACh.data_path)
set(gcf,'renderer','painters')
saveas(gcf,'NT_imaging_fig.pdf')
cd([]) %Add path as character array

%% Custom functions called in this script

function ind_data_table = get_individual_Ach_opto_fly(input_fly, fly_details, fly_idx)

%% Normalizes TTL trace to be between 0 and 1

TTL_input = input_fly.TTL;
nTTL = TTL_input(:,1)+5000;
nTTL = nTTL-min(nTTL);
nTTL = nTTL./max(nTTL);

%% Generates vector with non-zero elements only at indices of TTL onset

TTLidx = zeros(length(nTTL),1);

for loop_idx1 = 1:length(nTTL)
    if nTTL(loop_idx1,1) < 0.2
        binarized_TTL(loop_idx1,1) = 0;
    else
        binarized_TTL(loop_idx1,1) = 1;
    end
end

for loop_idx2 = 2:(numel(binarized_TTL)-1)
    if binarized_TTL(loop_idx2,1) == 1 && binarized_TTL(loop_idx2-1,1) == 0
        TTLidx(loop_idx2,1) = max(nTTL(loop_idx2:loop_idx2+70,1));
    end
end

TTLoutputs = find(TTLidx > 0);

%% Loops through img planes in indivdual sessions

curr_img = input_fly.plane_1;
img_output = get_Ach_opto_img_measures(curr_img, TTL_input, TTLoutputs);
st_img1 = img_output.st_img;
st_amp1 = img_output.img_amplitudes;
st_F0_1 = img_output.F0;
img_trace1 = img_output.dF_aligned;

curr_img = input_fly.plane_2;
img_output = get_Ach_opto_img_measures(curr_img, TTL_input, TTLoutputs);
st_img2 = img_output.st_img;
st_amp2 = img_output.img_amplitudes;
st_F0_2 = img_output.F0;
img_trace2 = img_output.dF_aligned;

st_F0_mean = mean([st_F0_1; st_F0_2],1);
st_amp_mean = mean([st_amp1; st_amp2],1);
img_traces_full = [img_output.TTL_aligned, img_trace1, img_trace2];

st_TTL = img_output.st_TTL;

ind_data_table = get_Ach_opto_data_table_row(img_output, st_img1, st_img2, ...
    st_TTL, st_amp_mean, st_F0_mean, img_traces_full, ...
    fly_details, fly_idx);

end

function img_output = get_Ach_opto_img_measures(current_plane, TTL_input, TTLoutputs)

f_r = get_Ach_opto_frame_rate; % frame-rate
raw_img = current_plane;
img_trace = raw_img(:,2) - raw_img(:,1);
pre_LED_window = round(0.5*(f_r));
post_LED_window = round(1.5*(f_r));
pre_LED_window_full = 200;
post_LED_window_full = 8200;
rs_1 = TTLoutputs-pre_LED_window;
rs_2 = TTLoutputs+post_LED_window;
rs_full = TTLoutputs(1)-pre_LED_window_full;
TTL_aligned_full = TTLoutputs-rs_full;
img_output = struct;

%% First, tests for ROI to be included

if raw_img(:,1) == raw_img(:,2)
    img_output.excl_vec(1,1) = 1;
    img_trace(:,1) = NaN(length(img_trace(:,1)),1);
    warning('!!!SOME DATA NOT INCLUDED!!!')
    pause
else
    img_output.excl_vec(1,1) = 0;
end

%% Then, computes dF

F0_mov = 2000;
curr_raw = img_trace(:,1);
dF = (curr_raw - movmean(curr_raw,F0_mov))./(movmean(curr_raw,F0_mov));
img_output.dF(:,1) = dF;
img_output.dF_aligned(:,1) = dF(rs_full:rs_full+post_LED_window_full,1);
img_output.TTL_aligned(:,1) = TTL_input(rs_full:rs_full+post_LED_window_full,1);
img_output.raw_img(:,1) = img_trace(:,1);
img_output.TTL_times_full = TTL_aligned_full;

%% ...LED-aligned img

for trial_idx = 1:size(rs_1,1)
    F_0 = mean(img_trace(rs_1(trial_idx,1):TTLoutputs(trial_idx),1),1);
    st_img(:,trial_idx) = (img_trace(rs_1(trial_idx,1):rs_2(trial_idx,1),1)-F_0)/abs(F_0);          
    st_TTL(:,trial_idx) = TTL_input(rs_1(trial_idx,1):rs_2(trial_idx,1),1);
    img_amplitudes(1,trial_idx) = mean(st_img(pre_LED_window:pre_LED_window+round(f_r*0.5),trial_idx),1);
    baseline_f(1,trial_idx) = F_0;
end

%% stores data in structure

img_output.st_img = st_img;
img_output.st_TTL = st_TTL;
img_output.img_amplitudes = img_amplitudes;
img_output.F0 = baseline_f;

end

function dm_row = get_Ach_opto_data_table_row(dm_input, st_img1, st_img2, st_TTL, st_amp_mean, st_F0_mean, img_traces_full, ...
    details, fly_idx)

dm_details_row = table;
dm_details_row.fly_id = fly_idx;
det_fl = struct;
det_fl.fly_name = details.fly_name;
dm_details_row.flyname = det_fl;


if contains(char(details.fly_name), "AchSnFR") == 1
    details.exp_group = "AchSnFR";
elseif contains(char(details.fly_name), "GluSnFR") == 1
    details.exp_group = "GluSnFR";    
elseif contains(char(details.fly_name), "GRABACh") == 1
    details.exp_group = "GRABAch";
end

dm_details_row.exp_group = details.exp_group;

dm_img_row = table;
dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;

amp_strct = struct;
amp_strct.img_amp = st_amp_mean;
dm_img_row.img_amplitudes = amp_strct;

F0_strct = struct;
F0_strct.baseline_f = st_F0_mean;
dm_img_row.F0 = F0_strct;

TTL_strct = struct;
TTL_strct.TTL_times_full = dm_input.TTL_times_full;
dm_img_row.TTL_times = TTL_strct;

trc_strct = struct;
trace_pre(:,:,1) = st_img1;
trace_pre(:,:,2) = st_img2;
trc_strct.trace_mean = mean(trace_pre,3);
dm_img_row.st_img_mean = trc_strct;

st_TTL_strct = struct;
st_TTL_strct.TTL_trace = st_TTL;
dm_img_row.st_TTL = st_TTL_strct;

trc_strct = struct;
trc_strct.traces_aligned = img_traces_full;
dm_img_row.img_traces_full = trc_strct;

dm_row = [dm_details_row, dm_excl_row, dm_img_row];

end

function group_summary = get_group_summaries(Ach_OPTO_data_table, exp_group)

group_table = Ach_OPTO_data_table(Ach_OPTO_data_table.exp_group == char(exp_group),:);

group_summary.st_img = [];
group_summary.st_img_full = [];
group_summary.st_amp = [];
group_summary.F_0 = [];
group_summary.st_TTL = [];

for loop_idx = 1:size(group_table,1)
    
    curr_table = group_table(loop_idx,:);
    group_summary.st_img = [group_summary.st_img, mean(curr_table.st_img_mean.trace_mean,2)];
    group_summary.st_img_full = [group_summary.st_img_full, curr_table.st_img_mean.trace_mean];
    group_summary.st_amp = [group_summary.st_amp, mean(curr_table.img_amplitudes.img_amp)];    
    group_summary.F_0 = [group_summary.F_0, mean(curr_table.F0.baseline_f(1,:))];
    group_summary.st_TTL = [group_summary.st_TTL, curr_table.st_TTL.TTL_trace];
    
end

end

function get_fig_panel_aligned_Ach_opto(SNFR_Ach_summary, SNFR_Glu_summary, GRAB_Ach_summary, x_pos, y_pos, ...
    color, map_cond, source_data_details)

sz_1 = 1.5;
ht_1 = 1.8;
ht_2 = 0.5;

if map_cond == "with_map"
    y_dist_1 = 0.2;
    y_dist_2 = 0.3;
    y_pos_trace = y_pos+(3*ht_2)+(2*y_dist_1)+y_dist_2;
    y_pos_map1 = y_pos;
    y_pos_map2 = y_pos_map1+ht_2+y_dist_1;
    y_pos_map3 = y_pos_map2+ht_2+y_dist_1;
    c_lim = [-0.1, 1];
    red_3 = [200, 85, 97];
    curr_colormap = get_default_two_color_map(abs(c_lim(1))*10000, abs(c_lim(2))*10000, ...
        [0.7, 0.7, 0.8], (red_3)./255);    
    ht_color_bar = (3*ht_2)+(2*y_dist_1);
    sz_colorbar = 0.1;
elseif map_cond == "no_map"
    y_pos_trace = y_pos;
end

ylm_img = [-0.05, 0.8];
f_r = get_Ach_opto_frame_rate;
light_on = f_r*0.5;
light_off = (f_r*0.5)+(0.5*f_r);

color_LED_bar = color.red;
plotting_color_AchSnFR = [0.5, 0.5, 0.5];
plotting_color_GluSnFR = color.navy;
plotting_color_GRAB_Ach = [0, 0, 0];

curr_img_AS = SNFR_Ach_summary.st_img;
curr_img_GS = SNFR_Glu_summary.st_img;
curr_img_AG = GRAB_Ach_summary.st_img;
curr_imgAS_full = SNFR_Ach_summary.st_img;
curr_imgGS_full = SNFR_Glu_summary.st_img;
curr_img_AG_full = GRAB_Ach_summary.st_img;

xlm_img = [1, length(mean(curr_img_AS,2))];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos_trace, sz_1, ht_1]);
hold on
fill(curr_pnl, [light_on, light_off, light_off, light_on],...
    [ylm_img(1), ylm_img(1), ylm_img(2), ylm_img(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
plot(curr_pnl, [light_on, light_off],...
    [ylm_img(2), ylm_img(2)], 'LineWidth', get_default_scale_bar_width, 'Color', color_LED_bar)
get_default_SEM_area_plot(curr_pnl, curr_img_AS, 1:size(curr_img_AS,1), plotting_color_AchSnFR)
get_default_SEM_area_plot(curr_pnl, curr_img_GS, 1:size(curr_img_GS,1), plotting_color_GluSnFR)
get_default_SEM_area_plot(curr_pnl, curr_img_AG, 1:size(curr_img_AG,1), plotting_color_GRAB_Ach)
get_default_separated_ax(curr_pnl, xlm_img(1), xlm_img(2), xlm_img(1), xlm_img(2), ylm_img(1), ...
    ylm_img(2), 0, 0.8,...
    [], 0.4, "linear", "linear", [], '\DeltaF/F',...
    'none', 'k', 'none', 'k', sz_1, ht_1)
get_default_annotation(x_pos+(sz_1/2.65), y_pos_trace+ht_1+0.2, '0.5 s', 'k', 'normal', "center")


if map_cond == "with_map"
    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos_map3, sz_1, ht_2]);
    imagesc(curr_pnl, 1:size(curr_imgGS_full,1), 1:size(curr_imgGS_full,2),...
        curr_imgGS_full')    
    colormap(curr_pnl, curr_colormap)
    curr_pnl.XAxis.Color = 'none';
    curr_pnl.YAxis.Color = 'none';
    clim(curr_pnl, c_lim)

    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.05, y_pos_map3, 0.1, ht_2]);
    plot(curr_pnl, [0, 0], [0, 1], 'Color', plotting_color_GluSnFR, 'LineWidth', 1.5)
    get_default_ax(curr_pnl, 0, 1, [], [], 0, 1,...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.1, ht_2)
    
    
    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos_map2, sz_1, ht_2]);
    imagesc(curr_pnl, 1:size(curr_imgAS_full,1), 1:size(curr_imgAS_full,2),...
        curr_imgAS_full')
    colormap(curr_pnl, curr_colormap)
    curr_pnl.XAxis.Color = 'none';
    curr_pnl.YAxis.Color = 'none';
    clim(curr_pnl, c_lim)

    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.05, y_pos_map2, 0.1, ht_2]);
    plot(curr_pnl, [0, 0], [0, 1], 'Color', plotting_color_AchSnFR, 'LineWidth', 1.5)
    get_default_ax(curr_pnl, 0, 1, [], [], 0, 1,...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.1, ht_2)


    curr_pnl = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos_map1, sz_1, ht_2]);
    imagesc(curr_pnl, 1:size(curr_img_AG_full,1), 1:size(curr_img_AG_full,2),...
        curr_img_AG_full')    
    colormap(curr_pnl,curr_colormap)
    curr_pnl.XAxis.Color = 'none';
    curr_pnl.YAxis.Color = 'none';
    clim(curr_pnl, c_lim)

    ac = colorbar(curr_pnl, 'westoutside');
    ac.Units = 'centimeters';
    ac.Position = [x_pos-(sz_colorbar+0.1), y_pos_map1, sz_colorbar, ht_color_bar];
    ylabel(ac, '\DeltaF/F')
    ac.FontSize = get_default_font_size;
    ac.Color = 'k';
    ac.TickLength = 0;
    ac.YLabel.Visible = 'on';
    ac.YLabel.Color = 'k';


    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.05, y_pos_map1, 0.1, ht_2]);
    plot(curr_pnl, [0, 0], [0, 1], 'Color', plotting_color_GRAB_Ach, 'LineWidth', 1.5)
    get_default_ax(curr_pnl, 0, 1, [], [], 0, 1,...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', 0.1, ht_2)

    curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos_map1-0.1, sz_1, 0.1]);
    plot(curr_pnl, [light_on, light_off], [0, 0], 'LineWidth', get_default_scale_bar_width, 'Color', color_LED_bar)
    get_default_ax(curr_pnl, xlm_img(1), xlm_img(2), [], [], 0, 1, ...
    [], [],[], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, 0.1)
    get_default_annotation(x_pos+(sz_1/2.65), y_pos_map1-0.2, '0.5 s', 'k', 'normal', "center")
end

%% Saves the source data
if isempty(source_data_details) == 0

    stimulus_aligned_traces = table;
    stimulus_aligned_traces.iGluSnFR = curr_imgGS_full;
    stimulus_aligned_traces.blank_space1 = NaN(size(curr_imgGS_full,1),1);
    stimulus_aligned_traces.iAchSnFR = curr_imgAS_full;
    stimulus_aligned_traces.blank_space2 = NaN(size(curr_imgGS_full,1),1);
    stimulus_aligned_traces.GRAB_ACh = curr_img_AG_full;
    writetable(stimulus_aligned_traces, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'stimulus_aligned')
  
end

end

function get_fig_panel_Ach_opto_traces(x_pos, y_pos, Ach_OPTO_data_table, source_data_details)

sz_1 = 4;
ht_1 = 1.2;
y_dis = 0.2;
f_r = get_Ach_opto_frame_rate;

ylim_1 = [-0.5, 0.8];
ylm_img = [0, 1];
line_width_plot = 0.25;
y_axis_scale_bar = 0.4;
color = get_color;
plot_color_GluSnFR = color.navy;
plot_color_AchSnFR = [0.5, 0.5, 0.5];
plot_color_GRABAch = [0, 0, 0];
color_LED_bar = color.red;

curr_trace_AchSnFR = Ach_OPTO_data_table.img_traces_full(6,1).traces_aligned(:,3);
curr_trace_GluSnFR = Ach_OPTO_data_table.img_traces_full(10,1).traces_aligned(:,2);
curr_trace_GRAB_Ach = Ach_OPTO_data_table.img_traces_full(20,1).traces_aligned(:,2);

stim_array = zeros(length(curr_trace_GluSnFR),1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1+2*(ht_1+y_dis)]);
hold on
for loop_idx = 1:20

    light_on = Ach_OPTO_data_table.TTL_times(1,1).TTL_times_full(loop_idx);
    light_off = light_on+(f_r*0.5);
    stim_array(light_on:round(light_off),1) = 1;
    fill(curr_pnl,[light_on, light_off, light_off, light_on],...
        [ylm_img(1), ylm_img(1), ylm_img(2), ylm_img(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
    plot(curr_pnl, [light_on, light_off],...
        [ylm_img(2), ylm_img(2)], 'LineWidth', 1, 'Color', color_LED_bar)
end
get_default_separated_ax(curr_pnl, 1, length(curr_trace_GluSnFR), [], [], ylm_img(1), ylm_img(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1+ht_1+y_dis)


curr_pnl1 = axes('Units','Centimeters','Position',[x_pos, y_pos+2*(ht_1+y_dis), sz_1, ht_1]);
hold on
get_default_img_trace_plot(curr_pnl1, curr_trace_GluSnFR, f_r, ...
    [], [], plot_color_GluSnFR, line_width_plot,  "with_y_axis", ...
    ylim_1, [-0.4, ylim_1(2)], y_axis_scale_bar, [], 0, "\DeltaF/F", []);

curr_pnl1 = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos+ht_1+y_dis, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_AchSnFR, f_r, ...
    [], [], plot_color_AchSnFR, line_width_plot,  "with_y_axis", ...
    ylim_1, [-0.4, ylim_1(2)], y_axis_scale_bar, [], 0, "\DeltaF/F", []);

curr_pnl1 = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace_GRAB_Ach, f_r, ...
    10, [], plot_color_GRABAch, line_width_plot,  "with_y_axis", ...
    ylim_1, [-0.4, ylim_1(2)], y_axis_scale_bar, "10 s", 0, "\DeltaF/F", []);

Y_pos_annotation_GluSnFR = y_pos+(2*(ht_1+y_dis))+ht_1+0.2;
get_default_annotation(x_pos+sz_1, Y_pos_annotation_GluSnFR, 'R23E10 > iGluSnFR', plot_color_GluSnFR, 'italic', "right")
Y_pos_annotation_AchSnFR = y_pos+ht_1+y_dis+ht_1;
get_default_annotation(x_pos+sz_1, Y_pos_annotation_AchSnFR, 'R23E10 > iAChSnFR', plot_color_AchSnFR, 'italic', "right")
Y_pos_annotation_GRABAch = y_pos+ht_1;
get_default_annotation(x_pos+sz_1, Y_pos_annotation_GRABAch, 'R23E10 > GRAB(ACh)', plot_color_GRABAch, 'italic', "right")

%% Saves the source data
if isempty(source_data_details) == 0

    example_traces = table;
    example_traces.stim_array = stim_array;
    example_traces.img_trace_iGluSnFR = curr_trace_GluSnFR;
    example_traces.img_trace_iAchSnFR = curr_trace_AchSnFR;
    example_traces.img_trace_GRAB_ACh = curr_trace_GRAB_Ach;
    writetable(example_traces, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'img_data')
  
end

end

function f_r = get_Ach_opto_frame_rate

f_r = 68.32;

end
