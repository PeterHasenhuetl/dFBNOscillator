% Analysis of dFBN GCaMP dynamics upon optogenetic stimulation of a
% GCaMP-negative, CsChrimson-positive dFBN subset.
% Code written by Peter Hasenhuetl.

tic

clear all

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'FLP_GCaMP_data.xlsx';

% go into folder with extracted traces and get fly names
cd([]) %Add path as character array

ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
aligned_img_withret = [];
aligned_img_NOret = [];
aligned_img_NOchrimson = [];
img_amp_NOchrimson = [];
img_amp_withret = [];
img_amp_NOret = [];
n_cellbodies_withret = [];
n_cellbodies_NOret = [];
n_cellbodies_NOchrimson = [];

img_amp_withret_FULL = [];
img_amp_NOret_FULL = [];
img_amp_NOchrimson_FULL = [];

img_NOchrimson_full = [];
img_withret_full = [];
img_NOret_full = [];

% loop through individual flies
for fly_idx = 1:n_flies
     cd([]) %Add path as character array
        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        curr_group = char(fly_input(fly_idx).name);
        
        if contains(curr_group,'withret') == 1
            [curr_img_trace, curr_img_amplitude, img_amplitudes_full, n_cellbodies, img_traces_full] = ...
                get_individual_conn_opto_fly(loaded_fly);
            aligned_img_withret = [aligned_img_withret, curr_img_trace];
            img_amp_withret = [img_amp_withret, curr_img_amplitude];
            img_amp_withret_FULL = [img_amp_withret_FULL, img_amplitudes_full];
            n_cellbodies_withret = [n_cellbodies_withret; n_cellbodies];
            img_withret_full = [img_withret_full, img_traces_full];
        
        elseif contains(curr_group,'NOret') == 1
            [curr_img_trace, curr_img_amplitude, img_amplitudes_full, n_cellbodies, img_traces_full] = ...
                get_individual_conn_opto_fly(loaded_fly);
            aligned_img_NOret = [aligned_img_NOret, curr_img_trace];
            img_amp_NOret = [img_amp_NOret, curr_img_amplitude];
            img_amp_NOret_FULL = [img_amp_NOret_FULL, img_amplitudes_full];
            n_cellbodies_NOret = [n_cellbodies_NOret; n_cellbodies];
            img_NOret_full = [img_NOret_full, img_traces_full];
            
        elseif contains(curr_group,'NOchrimson') == 1
            [curr_img_trace, curr_img_amplitude, img_amplitudes_full, n_cellbodies, img_traces_full] = ...
                get_individual_conn_opto_fly(loaded_fly);
            aligned_img_NOchrimson = [aligned_img_NOchrimson, curr_img_trace];
            img_amp_NOchrimson = [img_amp_NOchrimson, curr_img_amplitude];
            img_amp_NOchrimson_FULL = [img_amp_NOchrimson_FULL, img_amplitudes_full];
            n_cellbodies_NOchrimson = [n_cellbodies_NOchrimson; n_cellbodies];
            img_NOchrimson_full = [img_NOchrimson_full, img_traces_full];
        end
end

img_amp_NOchrimson_FULL(isnan(img_amp_NOchrimson_FULL)) = [];
img_amp_NOret_FULL(isnan(img_amp_NOret_FULL)) = [];
img_amp_withret_FULL(isnan(img_amp_withret_FULL)) = [];

cd([]) %Add path as character array

toc


histogram(img_amp_NOchrimson_FULL, 'Normalization', 'probability')
hold on
histogram(img_amp_NOret_FULL, 'Normalization', 'probability')
histogram(img_amp_withret_FULL, 'Normalization', 'probability')



%% Plotting the results

color = get_color;
color_img_withret = color.dark_gray;
color_img_NOret = color.medium_gray;
color_img_NOchrimson = color.light_gray;
color_LED_bar = color.red;

close all
figure('Name','mutually-exclusive optogenetics','Color','white',...
    'Units','centimeters','Position',[10 12 9 9],'Resize','off')
get_fig_panel_conn_opto(aligned_img_withret, aligned_img_NOret, aligned_img_NOchrimson, ...
    1, 1, color_img_withret, color_img_NOret, color_img_NOchrimson, color_LED_bar, source_data_details)
get_fig_panel_conn_opto_img_area(img_amp_withret, img_amp_NOret, img_amp_NOchrimson,...
    1+(4.1-1.2), 1, color_img_withret, color_img_NOret, color_img_NOchrimson, source_data_details)
get_fig_panel_conn_opto_confocal_images(1, 3.6, "no scale bar annotation")


cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'FLP_GCaMP_fig.pdf')
cd([]) %Add path as character array


%% Custom functions called in this script

function [aligned_img, img_amplitudes, img_amplitudes_full, n_cellbodies, img_traces_full] = get_individual_conn_opto_fly(input_fly)

if isfield(input_fly,'hemi1') == 1
    curr_trace = input_fly.hemi1;
    img_output_hemi1 = get_conn_opto_img_measures(curr_trace);
    st_img1 = img_output_hemi1.st_img;
    img_amplitudes1 = img_output_hemi1.img_amplitudes;
    img_amplitudes_full1 = img_output_hemi1.img_amplitudes_full;
    n_cellbodies(1) = img_output_hemi1.n_cellbodies;
    traces_full1 = img_output_hemi1.curr_full_trace;
else
    st_img1 = [];
    img_amplitudes1 = [];
    n_cellbodies(1) = NaN;
    img_amplitudes_full1 = NaN;
    traces_full1 = NaN(7900,1);
end  

     
if isfield(input_fly,'hemi2') == 1
    curr_trace = input_fly.hemi2;
    img_output_hemi2 = get_conn_opto_img_measures(curr_trace);
    st_img2 = img_output_hemi2.st_img;
    img_amplitudes2 = img_output_hemi2.img_amplitudes;
    img_amplitudes_full2 = img_output_hemi2.img_amplitudes_full;
    n_cellbodies(2) = img_output_hemi2.n_cellbodies;
    traces_full2 = img_output_hemi2.curr_full_trace;
else
    st_img2 = [];
    img_amplitudes2 = [];
    n_cellbodies(2) = NaN;
    img_amplitudes_full2 = NaN;
    traces_full2 = NaN(7900,1);
end
 
aligned_img = mean([st_img1, st_img2],2);
img_amplitudes = mean([img_amplitudes1, img_amplitudes2],2);
img_amplitudes_full = [img_amplitudes_full1, img_amplitudes_full2];
img_traces_full = [traces_full1, traces_full2];
img_traces_full(:,isnan(img_traces_full(1,:))) = []; 

end

function img_output = get_conn_opto_img_measures(curr_trace)

f_r = 29.13; % frame-rate

% normalizes TTL signal
TTL_input = curr_trace(:,1);
nTTL = TTL_input(:,1)+5000;
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


for loop_idx_2 = 2:(numel(binarized_TTL)-1)
    if binarized_TTL(loop_idx_2,1) == 1 && binarized_TTL(loop_idx_2-1,1) == 0
        TTLidx(loop_idx_2,1) = max(nTTL(loop_idx_2:loop_idx_2+70,1)); 
    end    
end


TTLoutputs = find(TTLidx > 0); % the indices of stimulus onsets
background_img = curr_trace(:,2);
raw_img = curr_trace(:,3:end);
img_traces = raw_img - background_img; % element-wise background subtraction

%%

pre_LED_window = round(5*(f_r));
post_LED_window = round(15*(f_r));
rs_1 = TTLoutputs-pre_LED_window;
rs_2 = TTLoutputs+post_LED_window;

img_output = struct;
    
%% First, tests for ROI to be included

for cell_idx = 1:size(img_traces,2)
    
    curr_cell = img_traces(:,cell_idx);
    
    if mean(curr_cell,1) < 20
        img_output.excl_vec(1,1) = 1;    
        img_traces(:,1) = NaN(length(curr_cell(:,1)),1);
    else
        img_output.excl_vec(1,1) = 0;
    end

    img_output.excl_vec(1,1) = 0;
 
    %% LED-aligned img
    clear st_img st_TTL img_amplitudes baseline_f F_0
    for trial_idx = 1:length(rs_1)
        F_0 = mean(curr_cell(rs_1(trial_idx):TTLoutputs(trial_idx),1),1);
        st_img(:,trial_idx) = (curr_cell(rs_1(trial_idx):rs_2(trial_idx),1)-F_0)/abs(F_0);          
        st_TTL(:,trial_idx) = TTL_input(rs_1(trial_idx):rs_2(trial_idx),1);
        img_amplitudes(1,trial_idx) = ...
            mean(st_img(pre_LED_window:pre_LED_window+round(f_r*5),trial_idx),1);
        baseline_f(1,trial_idx) = F_0;
        
    end
    
    curr_st_img(:,cell_idx) = mean(st_img,2);
    curr_amp(1,cell_idx) = mean(img_amplitudes,2);
    curr_full_trace(:,cell_idx) =  curr_cell(rs_1(1):end,1);
    
end
    
curr_full_trace = curr_full_trace(1:7900,:);

%% stores data in struct

img_output.st_img = mean(curr_st_img,2);
img_output.st_TTL = st_TTL;
img_output.img_amplitudes = mean(curr_amp,2);
img_output.img_amplitudes_full = curr_amp;
img_output.n_cellbodies = cell_idx;
img_output.curr_full_trace = curr_full_trace;
end

function get_fig_panel_conn_opto(aligned_img_withret, aligned_img_NOret, aligned_img_NOchrimson, ...
    x_pos, y_pos, color_img_withret, color_img_NOret, color_img_NOchrimson, color_LED_bar, source_data_details)


sz_1 = 1.5;
ht_1 = 1.8;
ylim_img = [-0.3, 0.1];
beginning_y_scale_bar = 0.05;
y_scale_bar = 0.1;
frame_rate = 29.13;
x_dist_annotation = 1.25;
y_dist_annotation = -0.2;
xlm_1 = [1, length(mean(aligned_img_withret,2))];

%% area plot marking the time of red light stimulation

light_on = frame_rate*5;
light_off = light_on+(frame_rate*5);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
fill(curr_pnl, [light_on, light_off, light_off, light_on], [0, 0, 1, 1], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
plot(curr_pnl, [light_on, light_off],...
    [1, 1], 'LineWidth', get_default_scale_bar_width, 'Color', color_LED_bar)
curr_pnl.Color = 'none';
curr_pnl.Box = 'off';
curr_pnl.XAxis.Color = 'none';
xlim(xlm_1)
ylim([0, 1])
curr_pnl.YAxis.Color = 'none';

%% plotting the actual imaging traces

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl, [1, length(mean(aligned_img_withret,2))], [0, 0], 'k:', 'LineWidth', 0.5)
get_default_SEM_area_plot(curr_pnl, aligned_img_NOret, 1:size(aligned_img_NOret,1), color_img_NOret)
get_default_SEM_area_plot(curr_pnl, aligned_img_NOchrimson, 1:size(aligned_img_NOchrimson,1), color_img_NOchrimson)
get_default_SEM_area_plot(curr_pnl, aligned_img_withret, 1:size(aligned_img_withret,1), color_img_withret)
plot(curr_pnl, [length(mean(aligned_img_withret,2)), length(mean(aligned_img_withret,2))],...
    [ylim_img(1)+beginning_y_scale_bar, ylim_img(1)+beginning_y_scale_bar+y_scale_bar], ...
    'LineWidth', get_default_scale_bar_width, 'Color', 'k')
ylim(ylim_img)
curr_pnl.Color = 'none';
curr_pnl.Box = 'off';
curr_pnl.XAxis.Color = 'none';
xlim([1, length(mean(aligned_img_withret,2))])
curr_pnl.YAxis.Color = 'none';
get_default_annotation(x_pos+sz_1-0.6, y_pos+0.65, {'0.1';'\DeltaF/F'}, 'k', 'normal', "left")

%% some annotations 

get_default_annotation(x_pos+0.425, y_pos+ht_1+0.2, '5 s', 'k','normal',"left")
get_default_group_annotation(x_pos, y_dist_annotation+y_pos+ht_1+0.55, ...
    "With ATR", color_img_withret, [0, 0, 0], "line")
get_default_group_annotation(x_pos+x_dist_annotation, y_dist_annotation+y_pos+ht_1+0.55, ...
    "No ATR", color_img_NOret, [0, 0, 0], "line")
get_default_group_annotation(x_pos+(2*x_dist_annotation), y_dist_annotation+y_pos+ht_1+0.55, ...
    "GCaMP only", color_img_NOchrimson, [0, 0, 0], "line")


%% Saves the source data

aligned_img_NOret = array2table(aligned_img_NOret);
aligned_img_NOchrimson = array2table(aligned_img_NOchrimson);
aligned_img_withret = array2table(aligned_img_withret);
writetable(aligned_img_NOret, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'aligned_img_NOret')
writetable(aligned_img_NOchrimson, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'aligned_img_NOchrimson')
writetable(aligned_img_withret, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'aligned_img_withret')

end

function get_fig_panel_conn_opto_img_area(img_amp_withret, img_amp_NOret, img_amp_NOchrimson,...
    x_pos, y_pos, color_img_withret, color_img_NOret, color_img_NOchrimson, source_data_details)

xlm_1 = [0.5, 3.5];
ylm_1 = [-0.3, 0.1];
sz_1 = 0.4*3;
ht_1 = 1.8;

%% plotting the scatter plot of areas under img curves

input_data = [img_amp_NOchrimson, img_amp_NOret, img_amp_withret];
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(pnl_1, [0.5, 3.5], [0, 0], 'k:', 'LineWidth', 0.5)
get_default_scatter_group(pnl_1, img_amp_NOchrimson, 1, color_img_NOchrimson)
get_default_scatter_group(pnl_1, img_amp_NOret, 2, color_img_NOret)
get_default_scatter_group(pnl_1, img_amp_withret, 3, color_img_withret)
get_default_axis_limits_warning(input_data, ylm_1)
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_1(1), ylm_1(2), -0.3, 0.1,...
    [], 0.1, "linear", "linear", [], 'Mean \DeltaF/F',...
    'none', 'k', 'none', 'k', sz_1, ht_1)

%% annotation with group sizes

n_groups = [length(img_amp_NOchrimson), length(img_amp_NOret), length(img_amp_withret)];
get_default_scatter_n_numbers_annotation(x_pos, y_pos, sz_1, n_groups, [], 'k', 'italic')

%% Saves the source data

mean_dF = NaN(15,3);
mean_dF(1:length(img_amp_NOchrimson),1) = img_amp_NOchrimson;
mean_dF(1:length(img_amp_NOret),2) = img_amp_NOret;
mean_dF(1:length(img_amp_withret),3) = img_amp_withret;
mean_dF = array2table(mean_dF);
mean_dF.Properties.VariableNames = {'NOchrimson','NOret','withret'};
writetable(mean_dF, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'mean_dF')

end

function get_fig_panel_conn_opto_confocal_images(x_pos, y_pos, image_scale_bar_cond)

%% loads images
img1 = imread('full_brain_50micronsSB.png');
img2 = imread('mutually_exclusive_cell_bodies_10micronsSB.png');

%% finds scale bar in image and plots a thicker bar on top of it (thicker in width, not in length!) 

% finds pixels in image that are white (for finding scale bar)
[sb_rows1, sb_columns1] = find(mean(img1,3) == 255);
[sb_rows2, sb_columns2] = find(mean(img2,3) == 255);

%corrects for minor dispalcements of scale bar in y-axis; DOES NOT CHANGE LENGTH OF SCALE BAR
actual_sb_rows1 = [min([sb_rows1(1),sb_rows2(1)]), min([sb_rows1(end),sb_rows2(end)])]; 

%% plots the images

sz_1 = 2.5;
sz_2 = 1.5;
ht_1 = 1.5;
x_dist = 0.1;
annot_y_dist = 0.2;

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
image(curr_pnl, img1)
hold on
plot(curr_pnl, [sb_columns1(1), sb_columns1(end)], [actual_sb_rows1(1), actual_sb_rows1(end)], 'w', 'LineWidth', 1)
curr_pnl.Color = 'none';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
ylim([1, 1200])

if image_scale_bar_cond == "with scale bar annotation"
    get_default_annotation(x_pos+sz_1-0.55, y_pos+0.25, ['50 \mu','m'], 'w', 'normal', "left")
end

curr_pnl = axes('Units', 'Centimeters', 'Position',[x_pos+sz_1+x_dist, y_pos, sz_2, ht_1]);
image(curr_pnl, img2)
hold on
plot(curr_pnl, [sb_columns2(1), sb_columns2(end)], [actual_sb_rows1(1), actual_sb_rows1(end)], 'w', 'LineWidth', 1)
curr_pnl.Color = 'none';
curr_pnl.XAxis.Color = 'none';
curr_pnl.YAxis.Color = 'none';
ylim([1, 1200])

if image_scale_bar_cond == "with scale bar annotation"
    get_default_annotation(x_pos+sz_1+x_dist+sz_2-0.55, y_pos+0.25, ['10 \mu','m'], 'w', 'normal', "left")
end

%% annotation of genetic strategy

% RGB codes taken from actual confocal images
curr_color_GCaMP = num2str([23, 200, 223]/255);
curr_color_Chrimson = num2str([255, 23, 255]/255);

annotation_1 = strcat(['dFBNs >  ','\color[rgb]{',curr_color_GCaMP,'}GCaMP','  \color{black}or  ', ...
    '\color[rgb]{',curr_color_Chrimson,'}CsChrimson::tdTomato']);
get_default_annotation(x_pos, y_pos+ht_1+annot_y_dist, annotation_1, 'k', 'italic', "left")

end
