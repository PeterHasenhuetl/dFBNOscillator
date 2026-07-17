% Analysis of heat-induced GRAB_DA signals before, during and after
% optogenetic stimulation of dFBNs.
% Code written by Peter Hasenhuetl.

clear all
tic

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'R23E10_GRAB_Chrimson_data.xlsx';

% go into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);

GRAB_OPTO_data_table1 = [];
GRAB_OPTO_data_table2 = [];
GRAB_OPTO_data_table3 = [];
img_CHR_pre = [];
img_CHR_stim = [];
img_CHR_post = [];
amp_CHR_pre = [];
amp_CHR_stim = [];
amp_CHR_post = [];
img_ctrl_pre = [];
img_ctrl_stim = [];
img_ctrl_post = [];
amp_ctrl_pre = [];
amp_ctrl_stim = [];
amp_ctrl_post = [];

% loop through individual flies
for fly_idx = 1:n_flies
    
     cd([]) %Add path as character array
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        
        [align_img, ind_data_table1, ind_data_table2, ind_data_table3] = ...
            get_individual_GRAB_OPTO_fly(loaded_fly, loaded_fly.fly_details, fly_idx);
        GRAB_OPTO_data_table1 = [GRAB_OPTO_data_table1; ind_data_table1];
        GRAB_OPTO_data_table2 = [GRAB_OPTO_data_table2; ind_data_table2];
        GRAB_OPTO_data_table3 = [GRAB_OPTO_data_table3; ind_data_table3];
        
        if ind_data_table1(1,:).CHR == "yes"
            img_CHR_pre = [img_CHR_pre, align_img.img_trace(:,1)];
            img_CHR_stim = [img_CHR_stim, align_img.img_trace(:,2)];
            img_CHR_post = [img_CHR_post, align_img.img_trace(:,3)];
            amp_CHR_pre = [amp_CHR_pre, align_img.amp(1,1)];
            amp_CHR_stim = [amp_CHR_stim, align_img.amp(2,1)];
            amp_CHR_post = [amp_CHR_post, align_img.amp(3,1)];
            
        elseif ind_data_table1(1,:).CHR == "no"
            img_ctrl_pre = [img_ctrl_pre, align_img.img_trace(:,1)];
            img_ctrl_stim = [img_ctrl_stim, align_img.img_trace(:,2)];
            img_ctrl_post = [img_ctrl_post, align_img.img_trace(:,3)];           
            amp_ctrl_pre = [amp_ctrl_pre, align_img.amp(1,1)];
            amp_ctrl_stim = [amp_ctrl_stim, align_img.amp(2,1)];
            amp_ctrl_post = [amp_ctrl_post, align_img.amp(3,1)];
            
        end



end


cd([]) %Add path as character array
toc


img_amplitudes.amp_CHR_pre = amp_CHR_pre;
img_amplitudes.amp_CHR_stim = amp_CHR_stim;
img_amplitudes.amp_CHR_post = amp_CHR_post;
img_amplitudes.amp_ctrl_pre = amp_ctrl_pre;
img_amplitudes.amp_ctrl_stim = amp_ctrl_stim;
img_amplitudes.amp_ctrl_post = amp_ctrl_post;
img_traces.img_CHR_pre = img_CHR_pre;
img_traces.img_CHR_stim = img_CHR_stim;
img_traces.img_CHR_post = img_CHR_post;
img_traces.img_ctrl_pre = img_ctrl_pre;
img_traces.img_ctrl_stim = img_ctrl_stim;
img_traces.img_ctrl_post = img_ctrl_post;

%% Plotting

close all
figure('Name','R23E10_GRAB_Chrimson_fig','Color','white',...
    'Units','centimeters','Position',[10 12 18.3 10],'Resize','off')
color = get_color;
color_no_chrimson1 = color.medium_gray;
color_no_chrimson2 = color.medium_gray;

curr_blue = [102, 153, 204];

color_with_chrimson1 = curr_blue./255;
color_with_chrimson2 = curr_blue./255;

[amplitude_plot1, amplitude_plot2] = get_fig_panel_GRAB_OPTO_punish(img_traces, img_amplitudes, 2, 0.5, ...
    color_no_chrimson1, color_no_chrimson2, color_with_chrimson1, color_with_chrimson2, "right", source_data_details);

cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'R23E10_GRAB_Chrimson_fig.pdf')
cd([]) %Add path as character array

%% Custom functions called in this script

function [align_img, ind_data_table1, ind_data_table2, ind_data_table3] = ...
    get_individual_GRAB_OPTO_fly(input_fly, fly_details, fly_idx)

TTL_input = input_fly.TTL_input;
fly_details.frame_rate_2P = get_GRAB_opto_frame_rate;
fly_details.bl = 10;
fly_details.n_trials_train = 3;
fly_details.trial_length_train = 7000; %trial length for train in ms
fly_details.baseline_jump = round(25*fly_details.frame_rate_2P); %seconds converted in #frames
fly_details.pre_trial_window_img_train = round(5*fly_details.frame_rate_2P); %seconds converted in #frames    
fly_details.post_TTL_lag = round(15.5*fly_details.frame_rate_2P);
fly_details.idx_triallength = round((fly_details.trial_length_train/1000)*fly_details.frame_rate_2P); 

%% Then, imaging
%calculates normalized UNFILTERED TTL signal
for session_idx = 1:3 % loops through the three sessions
    
    nTTL = TTL_input(:,session_idx)+5000;
    nTTL = nTTL-min(nTTL);
    nTTL = nTTL./max(nTTL);

    %logical vector with 'ones' at indices of TTL onset;
    TTLidx = zeros(length(nTTL),1);

    for loop_idx1 = 1:length(nTTL)
        if nTTL(loop_idx1,1) < 0.5 
            binarized_TTL(loop_idx1,1) = 0;
        else 
            binarized_TTL(loop_idx1,1) = 1;
        end
    end

    for loop_idx2 = 2:(numel(binarized_TTL)-1)
        if binarized_TTL(loop_idx2,1) == 1 && binarized_TTL(loop_idx2-1,1) == 0
            TTLidx(loop_idx2,1) = 1;
        end
    end
    
    % creates logical vector which indicates TTL output
    TTLoutputs(session_idx,1) = find(TTLidx == 1,1,'first');

end

%% Loops through img planes in indivdual sessions

curr_img = input_fly.imaging_pre;
ind_data_table1 = get_GRAB_OPTO_session(curr_img, TTL_input(:,1), TTLoutputs(1,1), fly_details, fly_idx);
curr_img = input_fly.imaging_stim;
ind_data_table2 = get_GRAB_OPTO_session(curr_img, TTL_input(:,2), TTLoutputs(2,1), fly_details, fly_idx);
curr_img = input_fly.imaging_post;
ind_data_table3 = get_GRAB_OPTO_session(curr_img, TTL_input(:,3), TTLoutputs(3,1), fly_details, fly_idx);
align_img = get_selected_ROI(ind_data_table1, ind_data_table2, ind_data_table3);

end

function [ind_data_table] = get_GRAB_OPTO_session(curr_img, TTL_input, TTLoutputs, fly_details, fly_idx)

fields_img_planes = fieldnames(curr_img);

for plane_idx = 1:length(fields_img_planes)
    plane_name = fields_img_planes{plane_idx};
    [output_img, ~] = get_GRAB_img_measures(curr_img.(string(plane_name)), TTL_input, TTLoutputs);
    ind_data_table(plane_idx,:) = get_GRAB_data_table_row(output_img, fly_details, fly_idx, plane_idx);
end

end


function align_img = get_selected_ROI(ind_data_table1, ind_data_table2, ind_data_table3)

SNR_array = ind_data_table1.SNR_ROI;
[SNR_row,  SNR_column] = find(SNR_array == max(max(SNR_array)));

align_img.img_trace(:,1) = ind_data_table1.st_img(SNR_row,1).st_trace(:,SNR_column);
align_img.img_trace(:,2) = ind_data_table2.st_img(SNR_row,1).st_trace(:,SNR_column);
align_img.img_trace(:,3) = ind_data_table3.st_img(SNR_row,1).st_trace(:,SNR_column);

align_img.amp(1,1) = ind_data_table1.dF_amplitude(SNR_row,SNR_column);
align_img.amp(2,1) = ind_data_table2.dF_amplitude(SNR_row,SNR_column);
align_img.amp(3,1) = ind_data_table3.dF_amplitude(SNR_row,SNR_column);

end


function [img_output, TTLoutputs_nu] = get_GRAB_img_measures(imaging, TTL_input, TTL_outputs)
% For different imaging sessions, the script calling this function
% loops through all imaging sessions.

f_r = get_GRAB_opto_frame_rate; % frame-rate
raw_green = table2array(imaging.img);
img_green = raw_green(:,2:end) - (raw_green(:,1));
pre_punish_window = round(10*(f_r));
post_punish_window = round(10*(f_r));
cut_off_start = TTL_outputs(1)-pre_punish_window;
cut_off_end = TTL_outputs(1)+post_punish_window;
TTLoutputs_nu = TTL_outputs-cut_off_start;
TTL_input = TTL_input(cut_off_start:cut_off_end,:);

img_output = struct;
for R_idx = 1:size(img_green,2)
    %% First, tests for ROI to be included
    
    if raw_green(:,R_idx+1) == raw_green(:,1)
        img_output.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);     
    elseif mean(img_green(:,R_idx),1) < 10
        img_output.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);
    else
        img_output.excl_vec(1,R_idx) = 0;
    end
    
    %% Then, computes dF
    
    curr_raw = img_green(cut_off_start:cut_off_end,R_idx);
    d_F = (curr_raw - mean(curr_raw(1:(pre_punish_window),1),1))./(mean(curr_raw(1:(pre_punish_window),1),1));    
    img_output.dF(:,R_idx) = d_F;
    img_output.raw_img(:,R_idx) = img_green(:,R_idx);

    %% ...punishment-aligned img

    st_img(:,R_idx) = d_F;
    st_TTL(:,R_idx) = TTL_input(:,1);
    cons_window = 4;
    img_output.SNR_ROI(1,R_idx) = max(movmean(st_img((pre_punish_window):(pre_punish_window)+round(f_r*cons_window),R_idx),10));
    img_output.dF_amplitude(1,R_idx) = mean(movmean(st_img((pre_punish_window):(pre_punish_window)+round(f_r*cons_window),R_idx),3));
            
end
   
%% stores data in structure
img_output.st_img = st_img;
img_output.st_TTL = st_TTL;

end

function dm_row = get_GRAB_data_table_row(dm_input, fly_details, fly_idx, plane_idx)

dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.plane_id = plane_idx;
dm_behav_row.SNR_ROI = dm_input.SNR_ROI;
dm_behav_row.flyname = fly_details.flyname;

cond_chrimson = contains(char(fly_details.flyname), "NO");

if cond_chrimson == 1
    fly_details.chrimson = "no";
else
    fly_details.chrimson = "yes";
end

dm_behav_row.CHR = fly_details.chrimson;

dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;
dm_excl_row.dF_amplitude = dm_input.dF_amplitude;

dm_img_row = table;
curr_trc = struct;
curr_trc.dF_trace = dm_input.dF;
dm_img_row.img_trace = curr_trc;
trc_raw = struct;
trc_raw.raw_trace = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;
st_img = struct;
st_img.st_trace = dm_input.st_img;
dm_img_row.st_img = st_img;
st_TTL = struct;
st_TTL.TTL_trace = dm_input.st_TTL;
dm_img_row.st_TTL = st_TTL;

%% Concatenates them

dm_row = [dm_behav_row, dm_excl_row, dm_img_row];

end

%% Plotting function

function [amplitude_plot1, amplitude_plot2] = get_fig_panel_GRAB_OPTO_punish(img_traces, img_amplitudes, x_pos, y_pos,...
    color_no_chrimson1, color_no_chrimson2, color_with_chrimson1, color_with_chrimson2, x_scale_bar, source_data_details)

f_s = get_default_font_size; % Font size
sz_1 = 0.75; % x-dimension of panel in cm
ht_1 = 2; % y-dimension of panel in cm
ht_2 = 2.2; % y-dimension of panel in cm
dist_y = 0.4;
ht_3 = ht_1+ht_2+dist_y;
f_r = get_GRAB_opto_frame_rate; % Hz
x_dist = 0.1;
sz_2 = (3*sz_1)+(2*x_dist);
ylim_img = [-0.05, 0.6];
ylm_scatter = [-0.05, 0.5];
punish_on = 10;
red_LED_on = punish_on-0.5;
red_LED_off = red_LED_on+3;
xlm_scatter = [0, sz_2];
mov_wind = 5;

%% first, plots the LED stim.

scatter_pun_coordinate1 = (f_r*(punish_on+1)); 
scatter_LED_on_coordinate1 = (f_r*(red_LED_on)); 
scatter_LED_off_coordinate1 = (f_r*(red_LED_off)); 
size_1 = length(img_traces.img_ctrl_pre);
norm_pun_pos = (scatter_pun_coordinate1/size_1)*sz_1;
scatter_pun_pos = [norm_pun_pos, sz_1+(norm_pun_pos)+(x_dist),...
    (2*sz_1)+(norm_pun_pos)+(2*x_dist)];
scatter_LED_pos(1) = sz_1+((scatter_LED_on_coordinate1/size_1)*sz_1)+(x_dist);
scatter_LED_pos(2) = sz_1+((scatter_LED_off_coordinate1/size_1)*sz_1)+(x_dist);

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_2, ht_3]);
fill(curr_pnl, [scatter_LED_pos(1), scatter_LED_pos(2), scatter_LED_pos(2), scatter_LED_pos(1)],...
    [0, 0, 1, 1], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
get_default_ax(curr_pnl, xlm_scatter(1), xlm_scatter(2), [], [], 0, 1,...
    [], [], [], [], "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_2, ht_3)

%% plots traces of punishment responses PRE red LED

curr_img1 = img_traces.img_ctrl_pre;
curr_img2 = img_traces.img_CHR_pre;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos+ht_2+dist_y, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, movmean(curr_img1,mov_wind), 1:size(curr_img1,1), color_no_chrimson1)
get_default_SEM_area_plot(curr_pnl, movmean(curr_img2,mov_wind), 1:size(curr_img2,1), color_with_chrimson1)
plot(curr_pnl, [14.56*punish_on, (14.56*punish_on)+(2*14.56)], [ylim_img(2), ylim_img(2)], 'k', 'LineWidth', 1)
if x_scale_bar == "left"
    plot(curr_pnl, [length(mean(curr_img1,2))-14.56*10, length(mean(curr_img1,2))], [ylim_img(1), ylim_img(1)], 'k', 'LineWidth', 1.5)
    get_default_annotation(x_pos+sz_1, y_pos+ht_2+dist_y-0.1, ...
        '10 s', 'k', 'normal', "right")
end
get_default_ax(curr_pnl, 1, size(curr_img1,1), [], [], ylim_img(1), ylim_img(2),...
    [], [], [], [], "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_1, ht_1)

%% plots traces of punishment responses DURING red LED

curr_img1 = img_traces.img_ctrl_stim;
curr_img2 = img_traces.img_CHR_stim;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+x_dist+sz_1, y_pos+ht_2+dist_y, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, movmean(curr_img1,mov_wind), 1:size(curr_img1,1), color_no_chrimson1)
get_default_SEM_area_plot(curr_pnl, movmean(curr_img2,mov_wind), 1:size(curr_img2,1), color_with_chrimson1)
plot(curr_pnl, [14.56*punish_on, (14.56*punish_on)+(2*14.56)], [ylim_img(2), ylim_img(2)], 'k', 'LineWidth', 1)
plot(curr_pnl, [14.56*9.5, (14.56*9.5)+(3*14.56)], [ylim_img(2)-0.01,ylim_img(2)-0.01], 'r', 'LineWidth', 1)
curr_pnl.Box = 'off';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
curr_pnl.Color = 'none';
xlim([1, length(mean(curr_img1,2))])
ylim(ylim_img)

%% plots traces of punishment responses POST red LED


x_dist_y_sb = 2;

curr_img1 = img_traces.img_ctrl_post;
curr_img2 = img_traces.img_CHR_post;
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos+(2*x_dist)+(2*sz_1), y_pos+ht_2+dist_y, sz_1, ht_1]);
get_default_SEM_area_plot(curr_pnl, movmean(curr_img1,mov_wind), 1:size(curr_img1,1), color_no_chrimson1)
hold on
get_default_SEM_area_plot(curr_pnl, movmean(curr_img2,mov_wind), 1:size(curr_img2,1), color_with_chrimson1)
plot(curr_pnl, [14.56*punish_on, (14.56*punish_on)+(2*14.56)], [ylim_img(2), ylim_img(2)], 'k', 'LineWidth', 1)
plot(curr_pnl, [length(mean(curr_img1,2))-14.56*x_dist_y_sb, length(mean(curr_img1,2))-14.56*x_dist_y_sb], [0.2,0.4], 'k', 'LineWidth', 1.5)
if x_scale_bar == "right"
    plot(curr_pnl, [length(mean(curr_img1,2))-14.56*10, length(mean(curr_img1,2))], [ylim_img(1), ylim_img(1)], 'k', 'LineWidth', 1.5)
    get_default_annotation(x_pos+(2*x_dist)+(3*sz_1), y_pos+ht_2+dist_y-0.1, '10 s', 'k', 'normal', "right")   
end
curr_pnl.Box = 'off';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
curr_pnl.Color = 'none';
xlim([1, length(mean(curr_img1,2))])
ylim(ylim_img)

get_default_annotation(x_pos+(2*x_dist)+(3*sz_1)+0.2, y_pos+ht_2+dist_y+(ht_1*0.6), ...
        {'0.2'; '\DeltaF/F'}, 'k', 'normal', "left")

%% plots the mean dF/F values below as individual points

amplitude_plot1 = [img_amplitudes.amp_ctrl_pre; img_amplitudes.amp_ctrl_stim; img_amplitudes.amp_ctrl_post];
amplitude_plot2 = [img_amplitudes.amp_CHR_pre; img_amplitudes.amp_CHR_stim; img_amplitudes.amp_CHR_post];

curr_pnl = axes('Units', 'Centimeters', 'Position',[x_pos, y_pos, sz_2, ht_2]);
hold on
for plot_idx1 = 1:size(amplitude_plot1,2)
     plot(curr_pnl, scatter_pun_pos, amplitude_plot1(:,plot_idx1), 'Color', color_no_chrimson1, 'LineWidth', 0.2);  
end



for plot_idx2 = 1:size(amplitude_plot2,2)
    plot(curr_pnl, scatter_pun_pos, amplitude_plot2(:,plot_idx2), 'Color', color_with_chrimson1, 'LineWidth', 0.2);
end

bar_center = mean(amplitude_plot1,2);
bar_extreme = std(amplitude_plot1,0,2)/sqrt(length(amplitude_plot1));
for error_bar_idx = 1:3
    plot(curr_pnl, [scatter_pun_pos(error_bar_idx), scatter_pun_pos(error_bar_idx)],...
        [bar_center(error_bar_idx)-bar_extreme(error_bar_idx), bar_center(error_bar_idx)+bar_extreme(error_bar_idx)],...
        'Color', color_no_chrimson2, 'LineWidth', 0.75)
end
plot(curr_pnl, scatter_pun_pos, mean(amplitude_plot1,2), 'Color', color_no_chrimson2, 'LineWidth', 1.5);

bar_center = mean(amplitude_plot2,2);
bar_extreme = std(amplitude_plot2,0,2)/sqrt(length(amplitude_plot2));
for error_bar_idx = 1:3
    plot(curr_pnl, [scatter_pun_pos(error_bar_idx), scatter_pun_pos(error_bar_idx)],...
        [bar_center(error_bar_idx)-bar_extreme(error_bar_idx), bar_center(error_bar_idx)+bar_extreme(error_bar_idx)],...
        'Color', color_with_chrimson2, 'LineWidth', 0.75)
end
plot(curr_pnl, scatter_pun_pos, mean(amplitude_plot2,2), 'Color', color_with_chrimson2, 'LineWidth', 1.5);
plot(curr_pnl, xlm_scatter, [0,0], ':k', 'LineWidth', 0.5);
get_default_ax(curr_pnl, xlm_scatter(1), xlm_scatter(2), [], [], ylm_scatter(1), ylm_scatter(2),...
    0, 0.6, [], 0.25, "linear", "linear", [], 'Mean \DeltaF/F', 'none', 'k', 'none', 'k', sz_2, ht_2)

%% some annotations

n_flies_ctrl = size(curr_img1,2); % computes n numbers of flies
n_flies_CHR = size(curr_img2,2); % computes n numbers of flies

annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos-0.1, y_pos+ht_1+ht_2+dist_y-0.05, 2, 0.5], 'string', ...
    'Pre', 'EdgeColor', 'none', 'FontSize', f_s, 'FontWeight', 'normal', 'Color', [0, 0, 0])
annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos-0.1+(x_dist)+(sz_1), y_pos+ht_1+ht_2+dist_y-0.05, 2, 0.5], 'string', ...
    'Light', 'EdgeColor', 'none', 'FontSize', f_s, 'FontWeight', 'normal', 'Color', [0, 0, 0])
annotation('textbox', 'Units','centimeters','Position',[x_pos-0.1+(2*x_dist)+(2*sz_1), y_pos+ht_1+ht_2+dist_y-0.05, 2, 0.5], 'string', ...
    'Post', 'EdgeColor', 'none', 'FontSize', f_s, 'FontWeight', 'normal', 'Color', [0, 0, 0])

annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos-0, y_pos+ht_1+dist_y+ht_2+0.3, 10, 0.5], 'string', ...
    ['+', 'Chrimson', ' (',num2str(n_flies_CHR), ' flies)'], 'EdgeColor', 'none',...
    'FontSize', f_s, 'FontWeight', 'normal', 'Color', color_with_chrimson1, 'Margin', 0, 'HorizontalAlignment', 'left')
annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos-0, y_pos+ht_1+dist_y+ht_2+0.1, 10, 0.5], 'string', ...
    ['- Chrimson',' (',num2str(n_flies_ctrl), ' flies)'], 'EdgeColor', 'none',...
    'FontSize', f_s, 'FontWeight', 'normal', 'Color', color_no_chrimson1, 'Margin', 0, 'HorizontalAlignment', 'left')


%% Saves the source data

traces_Chr = table;
traces_Chr.PRE = img_traces.img_CHR_pre;
traces_Chr.blank_space1 = NaN(size(img_traces.img_CHR_pre,1),1);
traces_Chr.STIM = img_traces.img_CHR_stim;
traces_Chr.blank_space2 = NaN(size(img_traces.img_CHR_pre,1),1);
traces_Chr.POST = img_traces.img_CHR_post;

traces_ctrl = table;
traces_ctrl.PRE = img_traces.img_ctrl_pre;
traces_ctrl.blank_space = NaN(size(img_traces.img_ctrl_pre,1),1);
traces_ctrl.STIM = img_traces.img_ctrl_stim;
traces_ctrl.blank_space = NaN(size(img_traces.img_ctrl_pre,1),1);
traces_ctrl.POST = img_traces.img_ctrl_post;

writetable(traces_Chr, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_Chrimson')
writetable(traces_ctrl, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'traces_ctrl')

chrimson_flies = array2table(amplitude_plot2');
chrimson_flies.Properties.VariableNames = {'pre', 'stim', 'post'};
writetable(chrimson_flies, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'mean_dF_Chrimson')

control_flies = array2table(amplitude_plot1');
control_flies.Properties.VariableNames = {'pre', 'stim', 'post'};
writetable(control_flies, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'mean_dF_ctrl')

end

function f_r = get_GRAB_opto_frame_rate

f_r = 14.56;

end


