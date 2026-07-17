% Analysis of dFBN patch-clamp recordings.
% Code written by Peter Hasenhuetl.

clear all
load('Allcellsloaded.mat')

cell_ids = [8, 15, 5, 17, 9]; % Example cells

V_1 = downsample(V1(:,cell_ids),50);
AP_1 = get_APs_ephys_only(V_1);
analyzed_cells = get_ephys_results(AP_1, 1000);

close all
figure('Name','Patch-clamp figure','Color','white',...
    'Units','centimeters','Position',[10 12 18 20],'Resize','off')
get_fig_panel_R23E10_patch_clamp(V_1, analyzed_cells, 2, 5);


%% Analysis functions

function group_results = get_ephys_results(t_signal_AP, f_r)
% summary measures

group_results.ST = [];
group_results.ISI = NaN(size(t_signal_AP,1),size(t_signal_AP,2));
group_results.spktms = NaN(size(t_signal_AP,1),size(t_signal_AP,2));

for loop_idx1 = 1:size(t_signal_AP,2)
    
    % measures based on inter-spike intervals
    spk_tms = ((find(t_signal_AP(:,loop_idx1) == 1))./f_r).*1000;
    group_results.spktms(1:length(spk_tms),loop_idx1) = spk_tms;
    i_s_i = diff(spk_tms); % inter-spike intervals in ms
    group_results.ISI(1:length(i_s_i),loop_idx1) = i_s_i;
    group_results.ISI_std(loop_idx1) = std(i_s_i); % standard deviation of inter-spike intervals
    group_results.ISI_mean(loop_idx1) = mean(i_s_i); % mean of inter-spike intervals
    group_results.ISI_CV(loop_idx1) = ...
        group_results.ISI_std(loop_idx1)./abs(group_results.ISI_mean(loop_idx1)); % coefficient of variation of inter-spike intervals
    group_results.ST = [group_results.ST; i_s_i];

end

end

function AP = get_APs_ephys_only(input_voltage)

AP = zeros(length(input_voltage),size(input_voltage,2));

for idx_1 = 1:size(input_voltage,2)
    min_peak = 2;
    [~, spktms_1] = findpeaks(zscore(movmean(diff(input_voltage(:,idx_1)),10)),'MinPeakHeight',min_peak);
    AP(spktms_1, idx_1) = 1;
end

end

%% Plotting function

function get_fig_panel_R23E10_patch_clamp(voltage_trace, analysis_summary, x_pos, y_pos)


sz1 = 1.25;
sz2 = 7;
sz3 = 0.5;
ht2 = 0.3;
xlm_hist = [0, 4];
ylm_hist = [0, 1];
major_x_ticks_scatter = 2;
major_y_ticks_scatter = [];
trace_dist = 0.25;
anno_dist = sz3-0.5;
hist_dist = 0.75;
trc_length = round(length(voltage_trace(:,1))./2);

pltrc_1 = voltage_trace(1:trc_length,:);
CV_data = analysis_summary.ISI_CV;
hist_data = analysis_summary.ISI;

x_pos_trace = x_pos+sz1+trace_dist+sz3+hist_dist;
for idx_plot = 1:5

    pnl_1 = axes('Units', 'Centimeters', 'Position',...
        [x_pos_trace, y_pos+((idx_plot-1)*(ht2+0.1)), sz2, ht2]);    
    plot(pnl_1, pltrc_1(:,idx_plot), 'Color', [0, 0, 0], 'LineWidth', 0.25)
    hold on
    if idx_plot == 1
        plot(pnl_1, [length(pltrc_1(:,1))-2000, length(pltrc_1(:,1))], [min(pltrc_1(:,1)), min(pltrc_1(:,1))],...
            'k', 'LineWidth', get_default_scale_bar_width)
        get_default_annotation(x_pos_trace+sz2, y_pos-0.1, '2 s', 'k', 'normal', "right")
    end
    pnl_1.XAxis.Color = 'none';
    pnl_1.YAxis.Color = 'none';
    pnl_1.Color = 'none';
    xlim([1, trc_length])    
    ylim([min(pltrc_1(:,idx_plot)), max(pltrc_1(:,idx_plot))])
    
    % trace scale bar    
    pnl_1 = axes('Units', 'Centimeters', 'Position',...
        [x_pos+sz1+trace_dist+sz2+0.1+sz3+hist_dist, y_pos+((idx_plot-1)*(ht2+0.1)), sz2, ht2]);    
    plot(pnl_1, [0, 0], [min(pltrc_1(:,idx_plot)), min(pltrc_1(:,idx_plot))+10], 'Color', [0, 0, 0],...
        'LineWidth', get_default_scale_bar_width)    
    pnl_1.XAxis.Color = 'none';    
    pnl_1.YAxis.Color = 'none';    
    pnl_1.Color = 'none';    
    ylim([min(pltrc_1(:,idx_plot)), max(pltrc_1(:,idx_plot))])    
    xlim([0, 1])

    pnl_1 = axes('Units', 'Centimeters', 'Position',...
        [x_pos+sz1+hist_dist, y_pos+((idx_plot-1)*(ht2+0.1)), sz3, ht2]);
    hist_plot = hist_data(isnan(hist_data(:,idx_plot))==0,idx_plot)./1000;
    histogram(pnl_1, hist_plot, 'Normalization', 'probability', 'NumBins', 20, 'BinLimits', [0, 4],...
        'FaceColor', [0, 0, 0], 'FaceAlpha', 1, 'EdgeColor', 'none')
    if idx_plot == 1
        x_axis_color = 'k';
    else
        x_axis_color = 'none';
    end
    
    if idx_plot == 5
        hold on
        plot(pnl_1, [0, 0], [0, 1], 'k', 'LineWidth', get_default_scale_bar_width)
    else
        pnl_1.YAxis.Color = 'none';
    end
    
    get_default_separated_ax(pnl_1, xlm_hist(1), xlm_hist(2), xlm_hist(1), xlm_hist(2), ylm_hist(1), ylm_hist(2),...
    ylm_hist(1), ylm_hist(2), major_x_ticks_scatter, major_y_ticks_scatter, "linear", "linear", 'ISI (s)', [],...
    string(x_axis_color), 'none', string(x_axis_color), 'none', sz3, ht2)
    
    annotation('textbox', 'Units','centimeters','Position',...
        [x_pos+sz1+anno_dist+hist_dist, y_pos+((idx_plot-1)*(ht2+0.1))-0.05, 5, 0.5], 'string', ...
        string(round(CV_data(idx_plot),2)), 'EdgeColor', 'none', 'FontSize', 4, 'FontWeight', 'normal')
    
end

get_default_annotation_rotated(x_pos_trace+sz2+0.2, y_pos, '10 mV', 'k', 'normal', "left")
get_default_annotation_rotated(x_pos+sz1+anno_dist+hist_dist-0.2, y_pos, '100%', 'k','normal', "left")
get_default_annotation(x_pos+sz1+anno_dist+hist_dist, y_pos+((idx_plot)*(ht2+0.1))+0.3, 'CV:', 'k','normal', "left")

end
