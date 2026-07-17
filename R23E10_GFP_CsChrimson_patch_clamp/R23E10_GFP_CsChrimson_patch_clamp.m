% Analysis of dFBN patch-clamp recordings during optogenetic stimulation.
% Code written by Peter Hasenhuetl.

clear all
load('Ephys_Opto_Data.mat');

curr_voltage_10Hz = EphysData.Voltage10Hz;
curr_TTL_10Hz = EphysData.TTL10Hz;
[P_spectrum_10Hz, f_10Hz, voltage_trace_10Hz, mean_burst_10Hz, example_burst_10Hz, example_TTL_10Hz] = ...
    get_P_spectrum(curr_voltage_10Hz, curr_TTL_10Hz, 1);

curr_voltage_20Hz = EphysData.Voltage20Hz;
curr_TTL_20Hz = EphysData.TTL20Hz;
[P_spectrum_20Hz, f_20Hz, voltage_trace_20Hz, mean_burst_20Hz, example_burst_20Hz, example_TTL_20Hz] = ...
    get_P_spectrum(curr_voltage_20Hz, curr_TTL_20Hz, 1);


%% Plots the data

color = get_color;

close all
figure('Name','dFBN Chrimson activation','Color','white',...
    'Units','centimeters','Position',[10 12 8.9 12],'Resize','off')
[stim_array_20Hz, curr_trace_20Hz] = get_fig_panel_voltage_opto_traces(0.9, 8.5, (curr_voltage_20Hz(:,1)), example_TTL_20Hz);
get_fig_panel_example_burst_opto(1, 6, example_burst_20Hz)
get_fig_panel_burst_opto(4, 6, mean_burst_20Hz)
get_fig_panel_ephys_power_spectrum(P_spectrum_20Hz, P_spectrum_10Hz, f_20Hz, 7, 6)


[stim_array_10Hz, curr_trace_10Hz] = get_fig_panel_voltage_opto_traces(0.9, 3.5, (curr_voltage_10Hz(:,1)), example_TTL_10Hz);
get_fig_panel_example_burst_opto(1, 1, example_burst_10Hz)
get_fig_panel_burst_opto(4, 1, mean_burst_10Hz)
get_fig_panel_ephys_power_spectrum(P_spectrum_10Hz, P_spectrum_20Hz, f_10Hz, 7, 1)
set(gcf,'renderer','Painters')
saveas(gcf,'Patch_clamp_optogenetics_fig.pdf')


%% Calls the custom functions

function [P_spectrum, f, voltage_trace, mean_burst, example_burst, example_TTL] = ...
    get_P_spectrum(EphysData, curr_TTL, example_idx)

tempData = EphysData;
downsampledData = detrend(downsample(tempData,20)); %imported at 10kHz, downsampled to 500Hz
Fs2 = 500;
[P_spectrum, f] = pspectrum(downsampledData, Fs2); 
voltage_trace = downsampledData;

for loop_idx = 1:size(EphysData,2)
    TTLoutputs = get_TTL(curr_TTL);
    k = 1;
    
    curr_burst = [];
    for TTL_idx = 2:length(TTLoutputs)-1
        curr_burst(:,k) = EphysData(TTLoutputs(TTL_idx)-2000:TTLoutputs(TTL_idx)+10000,loop_idx);
                        
        if loop_idx == example_idx
            example_burst(:,k) = curr_burst(:,k);    
            example_TTL = TTLoutputs;
        end

        curr_burst(:,k) = curr_burst(:,k) - mean(curr_burst(1:2000,k),1);
        
        k = k+1;
    end
    
    mean_burst(:,loop_idx) = mean(curr_burst,2);
    
end


end

function TTLoutputs = get_TTL(TTL_input)
voltage_opto_frame_rate = 20000;

%% Normalizes TTL trace to be between 0 and 1

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
inter_burst_interval = 0.4*voltage_opto_frame_rate; % 0.4 seconds, instead of 0.5 to have a tolerance window
TTLoutputs = TTLoutputs_raw;
TTLoutputs([inter_burst_interval*2; diff(TTLoutputs)] < 4000) = [];

end

function get_fig_panel_ephys_power_spectrum(P_spectrum_main, P_spectrum_comparison, f, x_pos, y_pos)

sz_1 = 1.8;
ht_1 = 1.8;
plotting_color1 = [0, 0, 0];
color = get_color;
plotting_color_comparison = color.light_gray;
xlm_1 = [0, 2];
ylm_1 = [0, 50];
major_x_ticks = 0.5;
major_y_ticks = 25;

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(pnl_1, [1, 1], ylm_1, 'Color', 'r', 'LineStyle', ':', "LineWidth", 0.25)
plot(pnl_1, [10, 10], ylm_1, 'Color', 'r', 'LineStyle', ':', "LineWidth", 0.25)
if ~isempty(P_spectrum_comparison)
    get_default_SEM_area_plot(pnl_1, P_spectrum_comparison, f.', plotting_color_comparison)
end
get_default_SEM_area_plot(pnl_1, P_spectrum_main, f.', plotting_color1)
get_default_separated_ax(pnl_1, xlm_1(1), xlm_1(2), xlm_1(1), xlm_1(2), ylm_1(1), ylm_1(2),...
    ylm_1(1), ylm_1(2), major_x_ticks, major_y_ticks, "linear", "linear", "Frequency (Hz)", "Power",...
    'k', 'k', 'k', 'k', sz_1, ht_1)    

end

function get_fig_panel_burst_opto(x_pos, y_pos, mean_burst)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_transient = [-0.1, 1];
ylm_transient = [-2, 18];
major_y_tick = 6;
st_tv = ((1:size(mean_burst,1))-(2000))/10000;
plotting_color1 = [0, 0, 0];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
get_default_SEM_area_plot(curr_pnl, mean_burst, st_tv, plotting_color1)
get_default_separated_ax(curr_pnl, xlm_transient(1), xlm_transient(2), 0, xlm_transient(2), ...
    ylm_transient(1), ylm_transient(2), 0, ylm_transient(2),...
    0.5, major_y_tick, "linear", "linear", 'Time (s)', '\DeltamV',...
    'k', 'k', 'k', 'k', sz_1, ht_1)


n_cells = size(mean_burst,2);
get_default_annotation(x_pos+sz_1, y_pos+ht_1+0.2, [num2str(n_cells), ' cells'], 'k', 'normal', "right")


end

function get_fig_panel_example_burst_opto(x_pos, y_pos, example_burst)

sz_1 = 1.8;
ht_1 = 1.8;
xlm_transient = [-0.1, 1];
ylm_transient = [-30, -5];
minor_y_ticks = 5;
st_tv = ((1:size(example_burst,1))-(2000))/10000;
plotting_color1 = [0.8, 0.8, 0.8];
plotting_color2 = [0, 0, 0];


curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl, st_tv, example_burst, 'Color', plotting_color1, 'LineWidth', 0.25)
plot(curr_pnl, st_tv, mean(example_burst,2), 'Color', plotting_color2)
get_default_separated_ax(curr_pnl, xlm_transient(1), xlm_transient(2), 0, xlm_transient(2), ...
    ylm_transient(1), ylm_transient(2), ylm_transient(1), ylm_transient(2),...
    0.5, minor_y_ticks, "linear", "linear", 'Time (s)', 'mV',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

get_default_annotation(x_pos+sz_1, y_pos+ht_1+0.2, 'Example cell', 'k', 'normal', "right")


end

function [stim_array, curr_trace] = get_fig_panel_voltage_opto_traces(x_pos, y_pos, curr_voltage, example_TTL)

sz_1 = 7.9;
ht_1 = 1.2;
f_r = 10000;
ylim_voltage = [-30, -5]; % to be set accordingly
major_y_tick = 8; % to be set accordingly
ylm_stim = [0, 1];
line_width_plot = 0.25;
color = get_color;
color_LED_bar = color.red;
plotting_color1 = [0.5, 0.5, 0.5];


limit_trace = 50;
curr_trace = curr_voltage(1:example_TTL(limit_trace+1),1);

stim_array = zeros(length(curr_trace),1);
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
for loop_idx = 1:limit_trace
    light_on = example_TTL(loop_idx);
    light_off = light_on+(f_r*0.5);
    stim_array(light_on:round(light_off),1) = 1;
    fill(curr_pnl,[light_on, light_off, light_off, light_on],...
        [ylm_stim(1), ylm_stim(1), ylm_stim(2), ylm_stim(2)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1)
    plot(curr_pnl, [light_on, light_off],...
        [ylm_stim(2), ylm_stim(2)], 'LineWidth', 1, 'Color', color_LED_bar)
end
get_default_separated_ax(curr_pnl, 1, length(curr_trace), [], [], ylm_stim(1), ylm_stim(2),...
    [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1)


curr_pnl1 = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1]);
get_default_img_trace_plot(curr_pnl1, curr_trace, f_r, ...
    2, [], plotting_color1, line_width_plot,  "with_y_axis", ...
    ylim_voltage, [ylim_voltage(1), ylim_voltage(2)], major_y_tick, "2 s", 0, "mV", []);


end

