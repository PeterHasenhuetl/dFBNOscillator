% Analysis of R84C10AD-R23E10DBD-splitGAL4-driven GCaMP signals ±VGlutRNAi.
% Code written by Peter Hasenhuetl.

clear all
tic

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'R84C10_split_data.xlsx';

% go into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
img_GD = [];
img_ctrl = [];
pwr_GD = [];
pwr_ctrl = [];
int_pwr_GD = [];
int_pwr_ctrl = [];
b_F_GD = [];
b_F_ctrl = [];
average_F_GD = [];
average_F_ctrl = [];
mean_d_F_GD = [];
mean_d_F_ctrl = [];

% loop through individual flies
for fly_idx = 1:n_flies
     
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        fly_details.fly_name = char(fly_input(fly_idx).name);
        output_img = get_curr_R84C10_fly(loaded_fly);
        
        if  contains(char(fly_details.fly_name), "GD_RNAi") == 1
            img_GD = [img_GD, output_img.img_trace];
            pwr_GD = [pwr_GD, output_img.pwr_trace];
            int_pwr_GD = [int_pwr_GD, output_img.int_SWA_power];
            b_F_GD = [b_F_GD, output_img.baseline_F];
            average_F_GD = [average_F_GD, output_img.average_F];
            mean_d_F_GD = [mean_d_F_GD, output_img.mean_d_F];
            
        elseif contains(char(fly_details.fly_name), "GCaMP_control") == 1
            img_ctrl = [img_ctrl, output_img.img_trace];
            pwr_ctrl = [pwr_ctrl, output_img.pwr_trace];
            int_pwr_ctrl = [int_pwr_ctrl, output_img.int_SWA_power];
            b_F_ctrl = [b_F_ctrl, output_img.baseline_F];
            average_F_ctrl = [average_F_ctrl, output_img.average_F];
            mean_d_F_ctrl = [mean_d_F_ctrl, output_img.mean_d_F];
            
        end

end

mean_SWA_pwr_control = mean(int_pwr_ctrl)
mean_SWA_pwr_GD = mean(int_pwr_GD)

mean_baseline_F_control = mean(b_F_ctrl)
mean_baseline_F_GD = mean(b_F_GD)

cd([]) %Add path as character array
toc

%% Plots the data

close all
figure('Name', 'VGlut-RNAi', 'Color', 'white',...
    'Units', 'centimeters', 'Position', [10 12 8.9 20], 'Resize', 'off')
plotting_color_control = ([144, 144, 144])./255;
plotting_color_GD = ([90, 127, 188])./255;
get_fig_panel_R84C10_split_power_spectrum(output_img, pwr_ctrl, pwr_GD, 1, 1,...
    plotting_color_control, plotting_color_GD)
curr_trc_ctrl = img_ctrl(1:2000,9);
ylim_1(1) = min(curr_trc_ctrl)-abs(min(curr_trc_ctrl)*0.1);
ylim_1(2) = max(curr_trc_ctrl)+abs(max(curr_trc_ctrl)*0.1);
get_fig_panel_R84C10_img_trace(curr_trc_ctrl, 0.1, 3.5, plotting_color_control, 1, ylim_1)
curr_trc_GD = img_GD(1:2000,12);
get_fig_panel_R84C10_img_trace(curr_trc_GD, 4.5, 3.5, plotting_color_GD, 1, ylim_1)


cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'R84C10_split_VGlut_RNAi_fig.pdf')
cd([]) %Add path as character array

%% Custom functions called in this script

function output_img = get_curr_R84C10_fly(input_fly)

hemi_fields = fieldnames(input_fly);

img_trace = [];
pwr_trace = [];
curr_int_SWA_power = [];
curr_b_F = [];
curr_average_F = [];
curr_mean_d_F = [];

for loop_idx = 1:length(hemi_fields)

    curr_hemi = input_fly.(char(hemi_fields(loop_idx)));
    curr_output = get_R84C10_split_analysis(curr_hemi);
    img_trace = [img_trace, curr_output.d_F];
    pwr_trace = [pwr_trace, curr_output.pwr_trace];
    curr_int_SWA_power = [curr_int_SWA_power, curr_output.int_delta_power];
    curr_b_F = [curr_b_F, curr_output.baseline_F];
    curr_average_F = [curr_average_F, curr_output.average_F];
    curr_mean_d_F = [curr_mean_d_F, curr_output.mean_d_F];

end

output_img.img_trace = img_trace;
output_img.pwr_trace = mean(pwr_trace, 2);
output_img.f = curr_output.f;
output_img.int_SWA_power = mean(curr_int_SWA_power,2);
output_img.baseline_F = mean(curr_b_F,2);
output_img.average_F = mean(curr_average_F,2);
output_img.mean_d_F = mean(curr_mean_d_F,2);

end

function curr_output = get_R84C10_split_analysis(curr_hemi)

curr_trace = curr_hemi(:,2) - mean(curr_hemi(:,1),1);
curr_trace = curr_trace(501:end,1);
f_r = 30.01; % frame-rate
sliding_window = 1001;
F_0 = get_mov_prctile(curr_trace, sliding_window);
baseline_F = mean(F_0);
average_F = mean(curr_trace);
d_F = (curr_trace - F_0)./(F_0);
mean_d_F = mean(d_F);

[one_sided_pwr_sctrm, f] = pspectrum(d_F, f_r);
curr_output.pwr_trace = one_sided_pwr_sctrm;
d_idx1 = find(f > 0.2,1,'first');
d_idx2 = find(f > 1,1,'first');
curr_output.int_delta_power = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
d_idx1 = find(f > 0.005,1,'first');
d_idx2 = find(f > 0.1,1,'first');

curr_output.int_infra_slow_power = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
curr_output.f = f;
curr_output.d_F = d_F;
curr_output.baseline_F = baseline_F;
curr_output.average_F = average_F;
curr_output.mean_d_F = mean_d_F;


end

%% Custom plotting functions

function get_fig_panel_R84C10_split_power_spectrum(output_img, pwr_ctrl, pwr_GD, x_pos, y_pos,...
    plotting_color_control, plotting_color_GD)

sz_1 = 2;
ht_1 = 1.8;
xlm_f = [0, 2];
scale_fac = 10^3;
ylm_pwr = ([0, 0.006])*scale_fac;
major_x_tick = 1;
major_y_tick = 0.003*scale_fac;
curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, pwr_ctrl*scale_fac, output_img.f.', plotting_color_control)
get_default_SEM_area_plot(curr_pnl, pwr_GD*scale_fac, output_img.f.', plotting_color_GD)
get_default_separated_ax(curr_pnl, xlm_f(1), xlm_f(2), xlm_f(1), xlm_f(2), ...
    ylm_pwr(1), ylm_pwr(2), ylm_pwr(1), ylm_pwr(2),...
    major_x_tick, major_y_tick, "linear", "linear", 'Frequency (Hz)', {'(\DeltaF/F)^2 x10^{-3}'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)

end

function get_fig_panel_R84C10_img_trace(curr_trace, x_pos, y_pos, plot_color, y_axis_scale_bar, ylim_1)

sz_1 = 4;
ht_1 = 1.4;
line_width_plot = 0.3;

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl, curr_trace, 30.01, ...
    5, y_axis_scale_bar, plot_color, line_width_plot,  "with_y_scale_bar", ...
    ylim_1, [0, 1], y_axis_scale_bar, "5 s", 0, [string(y_axis_scale_bar); "\DeltaF/F"], 0.25);



end