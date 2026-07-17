% Analysis of R23E10-GAL4-driven iGluSnFR signals in dFBN axons during
% optogenetic stimulation, in rested vs. sleep-deprived flies.
% Code written by Peter Hasenhuetl.

clear all
tic

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'R23E10_iGluSnfR_opto_sleep_deprivation_data.xlsx';

% go into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
plast_OPTO_data_table = [];

% loop through individual flies
for fly_idx = 1:n_flies
     cd([]) %Add path as character array
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));        
        plast_OPTO_data_table = [plast_OPTO_data_table; ...
            get_individual_plast_opto_fly(loaded_fly, loaded_fly.fly_details, fly_idx)];
end

cd([]) %Add path as character array

SNFR_SD_summary = get_group_summaries(plast_OPTO_data_table, "SD", "SNFR");
SNFR_C_summary = get_group_summaries(plast_OPTO_data_table, "C", "SNFR");

toc

%% Plotting the results

close all
figure('Name', 'R23E10_iGluSnfR_opto_sleep_deprivation', 'Color', 'white',...
    'Units', 'centimeters', 'Position', [10 12 18.3 15], 'Resize', 'off')
color = get_color;
plotting_color_C = color.rested;
plotting_color_SD = color.sleep_deprived;
color_LED_bar = color.red;
color_trace = [0, 0, 0];

get_fig_panel_opto_plast(SNFR_C_summary, SNFR_SD_summary, 2, 2, "no_indiv",...
    color_trace, plotting_color_C, plotting_color_SD, color_LED_bar, source_data_details)

cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'R23E10_iGluSnfR_opto_sleep_deprivation_fig.pdf')
cd([]) %Add path as character array

%% Custom functions called in this script

function ind_data_table = get_individual_plast_opto_fly(input_fly, fly_details, fly_idx)

%% Normalizes TTL trace to be between 0 and 1

TTL_input = input_fly.TTL;
nTTL = TTL_input(:,1)+5000;
nTTL = nTTL-min(nTTL);
nTTL = nTTL./max(nTTL);
  
%% Generates vector with non-zero elements only at indices of TTL onset

TTLidx = zeros(length(nTTL),1);
    
for loop_idx1 = 1:length(nTTL)
    if nTTL(loop_idx1,1) < 0.04
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

%% Extracts light-intensities from the amplitudes of the TTL signals

light_intensities = round(TTLidx(TTLoutputs),2);
for int_idx = 1:numel(light_intensities)
    if light_intensities(int_idx) < 0.06
        light_intensities(int_idx) = 0.05;
    elseif light_intensities(int_idx) >= 0.06 && light_intensities(int_idx) < 0.09
        light_intensities(int_idx) = 0.08;
    elseif light_intensities(int_idx) >= 0.09 && light_intensities(int_idx) < 0.125
        light_intensities(int_idx) = 0.1;
    elseif light_intensities(int_idx) >= 0.125 && light_intensities(int_idx) < 0.6
        light_intensities(int_idx) = 0.5;
    elseif light_intensities(int_idx) >= 0.8
        light_intensities(int_idx) = 1;    
    end  
end

fly_details.light_intensities = light_intensities';

%% Loops through img planes in indivdual sessions

curr_img = input_fly.imaging_opto;
fields_img_planes = fieldnames(curr_img);

for plane_idx = 1:length(fields_img_planes)
    plane_name = fields_img_planes{plane_idx};
    img_output = get_plast_opto_img_measures(curr_img.(string(plane_name)), TTL_input, TTLoutputs);
    
    if plane_idx == 1
        st_img1 = img_output.st_img;
        st_amp1 = img_output.img_amplitudes;
        st_F0_1 = img_output.F0;
    elseif plane_idx == 2
        st_img2 = img_output.st_img;
        st_amp2 = img_output.img_amplitudes;
        st_F0_2 = img_output.F0;
    end
end

st_F0_mean = mean([st_F0_1; st_F0_2],1);
st_amp_mean = mean([st_amp1; st_amp2],1);

ind_data_table = get_plast_opto_data_table_row(img_output, st_img1, st_img2, st_amp_mean, st_F0_mean, ...
    fly_details, fly_idx);

end

function img_output = get_plast_opto_img_measures(imaging, TTL_input, TTLoutputs)

f_r = get_plast_opto_frame_rate; % frame-rate
raw_img = table2array(imaging.img);
img_trace = raw_img(:,2) - raw_img(:,1);
pre_LED_window = round(0.5*(f_r));
post_LED_window = round(1.5*(f_r));
rs_1 = TTLoutputs-pre_LED_window;
rs_2 = TTLoutputs+post_LED_window;

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

curr_raw = img_trace(:,1);    
dF = (curr_raw - mean(curr_raw(1:(pre_LED_window),1),1))./(mean(curr_raw(1:(pre_LED_window),1),1));
img_output.dF(:,1) = dF;    
img_output.raw_img(:,1) = img_trace(:,1);
    
%% ...LED-aligned img
    
for trial_idx = 1:length(rs_1)
    F_0 = mean(img_trace(rs_1(trial_idx):TTLoutputs(trial_idx),1),1);
    st_img(:,trial_idx) = (img_trace(rs_1(trial_idx):rs_2(trial_idx),1)-F_0)/abs(F_0);          
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

function dm_row = get_plast_opto_data_table_row(dm_input, st_img1, st_img2, st_amp_mean, st_F0_mean, ...
    details, fly_idx)

dm_details_row = table;
dm_details_row.fly_id = fly_idx;
det_fl = struct;
det_fl.fly_name = details.flyname;
dm_details_row.flyname = det_fl;
dm_details_row.light_intensities = details.light_intensities;
int_strct = struct;
int_strct.st_int = details.light_intensities;
dm_details_row.light_intensities = int_strct;

cond_SD = contains(char(details.flyname), "SD");

if cond_SD == 1
    details.exp_group = "SD";
else
    details.exp_group = "C";
end

dm_details_row.exp_group = details.exp_group;

cond_genotype = contains(char(details.flyname), "GC");

if cond_genotype == 1
    details.geno_type = "GCAMP";
else
    details.geno_type = "SNFR";
end

dm_details_row.geno_type = details.geno_type;

dm_img_row = table;
dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;

amp_strct = struct;
amp_strct.img_amp = st_amp_mean;
dm_img_row.img_amplitudes = amp_strct;

F0_strct = struct;
F0_strct.baseline_f = st_F0_mean;
dm_img_row.F0 = F0_strct;

trc_strct = struct;
trace_pre(:,:,1) = st_img1;
trace_pre(:,:,2) = st_img2;
trc_strct.trace_mean = mean(trace_pre,3);
dm_img_row.st_img_mean = trc_strct;

dm_row = [dm_details_row, dm_excl_row, dm_img_row];

end

function group_summary = get_group_summaries(plast_OPTO_data_table, exp_group, geno_type)

group_table = plast_OPTO_data_table(plast_OPTO_data_table.exp_group == char(exp_group) & ...
    plast_OPTO_data_table.geno_type == char(geno_type),:);

group_summary.st_img_int1 = [];
group_summary.st_img_int2 = [];
group_summary.st_img_int3 = [];
group_summary.st_img_int4 = [];
group_summary.st_img_int5 = [];

group_summary.st_amp_int1 = [];
group_summary.st_amp_int2 = [];
group_summary.st_amp_int3 = [];
group_summary.st_amp_int4 = [];
group_summary.st_amp_int5 = [];


group_summary.F_0 = [];

for loop_idx = 1:size(group_table,1)
    
    curr_table = group_table(loop_idx,:);
    
    light_intensities = curr_table.light_intensities.st_int;

    group_summary.st_img_int1 = [group_summary.st_img_int1, mean(curr_table.st_img_mean.trace_mean(:,light_intensities == 0.05),2)];
    group_summary.st_img_int2 = [group_summary.st_img_int2, mean(curr_table.st_img_mean.trace_mean(:,light_intensities == 0.08),2)];
    group_summary.st_img_int3 = [group_summary.st_img_int3, mean(curr_table.st_img_mean.trace_mean(:,light_intensities == 0.1),2)];
    group_summary.st_img_int4 = [group_summary.st_img_int4, mean(curr_table.st_img_mean.trace_mean(:,light_intensities == 0.5),2)];
    group_summary.st_img_int5 = [group_summary.st_img_int5, mean(curr_table.st_img_mean.trace_mean(:,light_intensities == 1),2)];
    
    group_summary.st_amp_int1 = [group_summary.st_amp_int1, mean(curr_table.img_amplitudes.img_amp(1,light_intensities == 0.05))];
    group_summary.st_amp_int2 = [group_summary.st_amp_int2, mean(curr_table.img_amplitudes.img_amp(1,light_intensities == 0.08))];
    group_summary.st_amp_int3 = [group_summary.st_amp_int3, mean(curr_table.img_amplitudes.img_amp(1,light_intensities == 0.1))];
    group_summary.st_amp_int4 = [group_summary.st_amp_int4, mean(curr_table.img_amplitudes.img_amp(1,light_intensities == 0.5))];
    group_summary.st_amp_int5 = [group_summary.st_amp_int5, mean(curr_table.img_amplitudes.img_amp(1,light_intensities == 1))];
        
    group_summary.F_0 = [group_summary.F_0, mean(curr_table.F0.baseline_f(1,:))];
    
end

end

function get_fig_panel_opto_plast(SNFR_C_summary, SNFR_SD_summary, x_pos, y_pos, cond_indiv_points,...
    color_trace, plotting_color_C, plotting_color_SD, color_LED_bar, source_data_details)

sz_1 = 0.7;
ht_2 = 1.5;
ht_1 = 1;
x_dist = 0.1;
sz_2 = (2*sz_1)+(x_dist);
ylm_img = [-0.05, 0.8];
ylm_curve = [0, 0.6];
xlm_curve = [0, 1];
dist_y = 0.2;
f_r = get_plast_opto_frame_rate;
light_on = f_r*0.5;
light_off = (f_r*0.5)+(0.5*f_r);
light_intensity_vector = [0.05, 0.08, 0.1, 0.5, 1];

color_diff = 0.175;
plotting_color1 = color_trace;
plotting_color2 = plotting_color1+color_diff;
plotting_color3 = plotting_color2+color_diff;
plotting_color4 = plotting_color3+color_diff;
plotting_color5 = plotting_color4+color_diff;

% first plotting control flies
curr_img1 = SNFR_C_summary.st_img_int1;
curr_img2 = SNFR_C_summary.st_img_int2;
curr_img3 = SNFR_C_summary.st_img_int3;
curr_img4 = SNFR_C_summary.st_img_int4;
curr_img5 = SNFR_C_summary.st_img_int5;

xlm_img = [1,length(mean(curr_img1,2))];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos+ht_2+dist_y, sz_1, ht_1]);
hold on
fill(curr_pnl,[light_on, light_off, light_off, light_on],...
    [ylm_img(1), ylm_img(1), ylm_img(2), ylm_img(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
plot(curr_pnl, [light_on,light_off],...
    [ylm_img(2), ylm_img(2)], 'LineWidth', 1, 'Color', color_LED_bar)
get_default_SEM_area_plot(curr_pnl, curr_img1, 1:size(curr_img1,1), plotting_color1)
get_default_SEM_area_plot(curr_pnl, curr_img2, 1:size(curr_img2,1), plotting_color2)
get_default_SEM_area_plot(curr_pnl, curr_img3, 1:size(curr_img3,1), plotting_color3)
get_default_SEM_area_plot(curr_pnl, curr_img4, 1:size(curr_img4,1), plotting_color4)
get_default_SEM_area_plot(curr_pnl, curr_img5, 1:size(curr_img5,1), plotting_color5)
get_default_ax(curr_pnl, xlm_img(1), xlm_img(2), [], [], ylm_img(1), ylm_img(2), ...
    0, 0.8, [], 0.4, "linear", "linear", [], '\DeltaF/F', 'none', 'k', 'none', 'k', sz_1, ht_1)

get_default_annotation(x_pos+0.1, y_pos+ht_2+ht_1+dist_y+0.2, '0.5 s', 'k', 'normal', "left")

% then plotting sleep-deprived flies
curr_img1 = SNFR_SD_summary.st_img_int1;
curr_img2 = SNFR_SD_summary.st_img_int2;
curr_img3 = SNFR_SD_summary.st_img_int3;
curr_img4 = SNFR_SD_summary.st_img_int4;
curr_img5 = SNFR_SD_summary.st_img_int5;

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+x_dist+sz_1, y_pos+ht_2+dist_y, sz_1, ht_1]);
hold on
fill(curr_pnl,[light_on, light_off, light_off, light_on],...
    [ylm_img(1), ylm_img(1), ylm_img(2), ylm_img(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
plot(curr_pnl, [light_on, light_off],...
    [ylm_img(2), ylm_img(2)], 'LineWidth', 1, 'Color', color_LED_bar)
get_default_SEM_area_plot(curr_pnl, curr_img1, 1:size(curr_img1,1), plotting_color1)
get_default_SEM_area_plot(curr_pnl, curr_img2, 1:size(curr_img2,1), plotting_color2)
get_default_SEM_area_plot(curr_pnl, curr_img3, 1:size(curr_img3,1), plotting_color3)
get_default_SEM_area_plot(curr_pnl, curr_img4, 1:size(curr_img4,1), plotting_color4)
get_default_SEM_area_plot(curr_pnl, curr_img5, 1:size(curr_img5,1), plotting_color5)

curr_pnl.Box = 'off';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
curr_pnl.Color = 'none';
curr_pnl.TickLength = [0, 0];
xlim([1,length(mean(curr_img1,2))])
ylim(ylm_img)

get_default_annotation(x_pos+0.1+x_dist+sz_1, y_pos+ht_2+ht_1+dist_y+0.2, '0.5 s', 'k', 'normal', "left")

%% "dose" response curve

amplitude_plot_C = [SNFR_C_summary.st_amp_int1; SNFR_C_summary.st_amp_int2; ...
    SNFR_C_summary.st_amp_int3; SNFR_C_summary.st_amp_int4; SNFR_C_summary.st_amp_int5];
amplitude_plot_SD = [SNFR_SD_summary.st_amp_int1; SNFR_SD_summary.st_amp_int2; ...
    SNFR_SD_summary.st_amp_int3; SNFR_SD_summary.st_amp_int4; SNFR_SD_summary.st_amp_int5];


curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_2, ht_2]);
hold on

y = mean(amplitude_plot_C,2);
x = [0.05, 0.08, 0.1, 0.5, 1]';
transfer_hyperbola = '(a*x)./(b+x)';
f1 = fit(x, y, transfer_hyperbola);
fitted_curve_C = feval(f1, 0:0.001:1);
plot(0:0.001:1, fitted_curve_C, 'Color', plotting_color_C)

if cond_indiv_points == "show_indiv"
    
    for plot_idx1 = 1:size(amplitude_plot_C,2)
        scatter(curr_pnl, light_intensity_vector', amplitude_plot_C(:,plot_idx1), 4,...
            'MarkerFaceColor', plotting_color_C, 'MarkerEdgeColor', 'none');
    end

    for plot_idx2 = 1:size(amplitude_plot_SD,2)
        scatter(curr_pnl, light_intensity_vector', amplitude_plot_SD(:,plot_idx2), 4,...
            'MarkerFaceColor', plotting_color_SD, 'MarkerEdgeColor', 'none');
    end

end

errorbar_center = mean(amplitude_plot_C,2);
errorbar_extreme = std(amplitude_plot_C,0,2)/sqrt(size(amplitude_plot_C,2));
for error_bar_idx = 1:5
    plot(curr_pnl,[light_intensity_vector(error_bar_idx), light_intensity_vector(error_bar_idx)],...
        [errorbar_center(error_bar_idx)-errorbar_extreme(error_bar_idx),...
        errorbar_center(error_bar_idx)+errorbar_extreme(error_bar_idx)],...
        'Color', plotting_color_C, 'LineWidth', 0.5)
end

scatter(curr_pnl, light_intensity_vector, mean(amplitude_plot_C,2), 4, 'MarkerFaceColor', plotting_color_C,...
   'MarkerEdgeColor', 'none');

y = mean(amplitude_plot_SD,2);
x = [0.05, 0.08, 0.1, 0.5, 1]';
transfer_hyperbola = '(a*x)./(b+x)';
f1 = fit(x, y, transfer_hyperbola);
fitted_curve_SD = feval(f1, 0:0.001:1);
hold on
plot(0:0.001:1, fitted_curve_SD, 'Color', plotting_color_SD)

errorbar_center = mean(amplitude_plot_SD,2);
errorbar_extreme = std(amplitude_plot_SD,0,2)/sqrt(size(amplitude_plot_SD,2));
for error_bar_idx = 1:5
    plot(curr_pnl,[light_intensity_vector(error_bar_idx), light_intensity_vector(error_bar_idx)],...
        [errorbar_center(error_bar_idx)-errorbar_extreme(error_bar_idx),...
        errorbar_center(error_bar_idx)+errorbar_extreme(error_bar_idx)],...
        'Color', plotting_color_SD, 'LineWidth', 0.5)
end
scatter(curr_pnl, light_intensity_vector, mean(amplitude_plot_SD,2), 4, 'MarkerFaceColor', plotting_color_SD,...
    'MarkerEdgeColor', 'none');

get_default_ax(curr_pnl, xlm_curve(1), xlm_curve(2), xlm_curve(1), xlm_curve(2), ylm_curve(1), ...
    ylm_curve(2), ylm_curve(1), ylm_curve(2),...
    0.5, 0.3, "linear", "linear", 'Norm. light intensity', 'Mean \DeltaF/F',...
    'k', 'k', 'k', 'k', sz_2, ht_2)

get_default_annotation(x_pos+0.45, y_pos+ht_1+ht_2+dist_y, 'Ctrl.', plotting_color_C, 'normal', "left")
get_default_annotation(x_pos+0.45+sz_1+x_dist, y_pos+ht_1+ht_2+dist_y, 'S.d.', plotting_color_SD, 'normal', "left")


%% Saves the source data

traces_SD = table;
traces_SD.intensity1_fly = SNFR_SD_summary.st_img_int1;
traces_SD.blank_space1 = NaN(size(SNFR_SD_summary.st_img_int1,1),1);
traces_SD.intensity2_fly = SNFR_SD_summary.st_img_int2;
traces_SD.blank_space2 = NaN(size(SNFR_SD_summary.st_img_int2,1),1);
traces_SD.intensity3_fly = SNFR_SD_summary.st_img_int3;
traces_SD.blank_space3 = NaN(size(SNFR_SD_summary.st_img_int3,1),1);
traces_SD.intensity4_fly = SNFR_SD_summary.st_img_int4;
traces_SD.blank_space4 = NaN(size(SNFR_SD_summary.st_img_int4,1),1);
traces_SD.intensity5_fly = SNFR_SD_summary.st_img_int5;
writetable(traces_SD, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_SD')

traces_ctrl = table;
traces_ctrl.intensity1_fly = SNFR_C_summary.st_img_int1;
traces_ctrl.blank_space1 = NaN(size(SNFR_C_summary.st_img_int1,1),1);
traces_ctrl.intensity2_fly = SNFR_C_summary.st_img_int2;
traces_ctrl.blank_space2 = NaN(size(SNFR_C_summary.st_img_int2,1),1);
traces_ctrl.intensity3_fly = SNFR_C_summary.st_img_int3;
traces_ctrl.blank_space3 = NaN(size(SNFR_C_summary.st_img_int3,1),1);
traces_ctrl.intensity4_fly = SNFR_C_summary.st_img_int4;
traces_ctrl.blank_space4 = NaN(size(SNFR_C_summary.st_img_int4,1),1);
traces_ctrl.intensity5_fly = SNFR_C_summary.st_img_int5;
writetable(traces_ctrl, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_ctrl')

control_fly = amplitude_plot_C;
sd_fly = amplitude_plot_SD;

dose_response_control = table;
dose_response_control.norm_light_intensity = light_intensity_vector';
dose_response_control.control_fly = control_fly;
writetable(dose_response_control, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'dose_response_control')

dose_response_sd = table;
dose_response_sd.norm_light_intensity = light_intensity_vector';
dose_response_sd.sd_fly = sd_fly;
writetable(dose_response_sd, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'dose_response_sd')

x_values = 0:0.001:1;
dose_response_curve = table;
dose_response_curve.x_values = x_values';
dose_response_curve.SD = fitted_curve_SD;
dose_response_curve.C = fitted_curve_C;
writetable(dose_response_curve, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'dose_response_curve')

end

function f_r = get_plast_opto_frame_rate

f_r = 68.32;

end

