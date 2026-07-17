% Analysis of dFBN GCaMP dynamics upon dopamine application.
% Code written by Peter Hasenhuetl.


clear all
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'picospritzing_dopamine.xlsx';
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
nu_trace_mean = [];
nu_trace_all = [];
trace_red_aligned = [];
aligned_TTL = [];
baseline_length = round(14.56*40);
for fly_idx = 1:n_flies
    
    fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        left_green = loaded_fly.left_green;
        left_red = loaded_fly.left_red;
        TTL = loaded_fly.TTL;
    trace_green1 = left_green(:,2:end)-left_green(:,1);
    curr_trace_red = left_red(:,2:end)-left_red(:,1);
    curr_trace_green = trace_green1./curr_trace_red;
    
    % Finds alignment index
    TTL = (TTL-min(TTL))/max(TTL-min(TTL));
    idx_1 = 1;
    for loop_idx = 2:length(TTL)
        if TTL(loop_idx,1) < 0.15 && TTL(loop_idx-1,1) >= 0.15
            application_idx(idx_1) = loop_idx;
            idx_1 = idx_1+1;
        end
    end

    application_idx = application_idx(end-2);
    idx_1 = application_idx-baseline_length;
    idx_2 = application_idx+1100;
    F0_range = 1:baseline_length-1;
    curr_trace_green = curr_trace_green(idx_1:idx_2,:);
    curr_trace_green = (curr_trace_green-mean(curr_trace_green(F0_range,:)))./mean(curr_trace_green(F0_range,:));

    nu_trace_all = [nu_trace_all, curr_trace_green];
    nu_trace_mean =  [nu_trace_mean, mean(curr_trace_green,2)];

    curr_trace_red = curr_trace_red(idx_1:idx_2,:);
    curr_trace_red = (curr_trace_red-mean(curr_trace_red(F0_range,:)))./mean(curr_trace_red(F0_range,:));

    trace_red_aligned = [trace_red_aligned, mean(curr_trace_red,2)];
    aligned_TTL = [aligned_TTL,TTL(idx_1:idx_2,1)];

end

%% Plots the data

close all
figure('Name','connectomics analysis','Color','white',...
    'Units','centimeters','Position',[10 12 18 12],'Resize','off')
color = get_color;

sz_1 = 5;
ht_1 = 2.5;
x_pos = 1;
y_pos = 1;
x_sb = round(14.56*20);
y_sb = 0.2;
tdTom_color = color.orange;
ylim_1 = [-0.5, 0.75];
xlim_1 = [1, length(mean(nu_trace_all,2))];
curr_panel = axes('Units','Centimeters','Position',...
    [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_panel,[baseline_length,baseline_length],ylim_1,'r:')
get_default_SEM_area_plot(curr_panel, trace_red_aligned, xlim_1(1):xlim_1(2), tdTom_color)
get_default_SEM_area_plot(curr_panel, nu_trace_mean, xlim_1(1):xlim_1(2), color.dark_gray)
plot(curr_panel,[xlim_1(2)-x_sb,xlim_1(2)],...
    [ylim_1(1),ylim_1(1)],'k','LineWidth',get_default_scale_bar_width)
plot(curr_panel,[xlim_1(2),xlim_1(2)],...
    [ylim_1(1),ylim_1(1)+y_sb],'k','LineWidth',get_default_scale_bar_width)
curr_panel.Box = 'off';
curr_panel.Color = 'none';
curr_panel.XAxis.Color = 'none';
curr_panel.YAxis.Color = 'none';
xlim(xlim_1)
ylim(ylim_1)

scale_annotation_cond = 1;
if scale_annotation_cond == 1
    get_default_annotation(x_pos+sz_1+0.2, ...
        y_pos+0.3, '0.2\DeltaR/R', 'k','normal',"right")
    get_default_annotation(x_pos+sz_1+0.2, ...
        y_pos+0.5, '0.2\DeltaF/F', tdTom_color,'normal',"right")
    get_default_annotation(x_pos+sz_1, ...
        y_pos-0.2, '20 s', 'k','normal',"right")
end


cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'DA_picospritzing_fig.pdf')


if isempty(source_data_details) == 0

    stim_array = zeros(size(nu_trace_mean,1),1);
    stim_array(baseline_length,1) = 1;
    DA_pico_traces = table;
    DA_pico_traces.stim_array = stim_array;
    DA_pico_traces.trace_red_aligned = trace_red_aligned;
    DA_pico_traces.blank_space1 = NaN(size(trace_red_aligned,1),1);
    DA_pico_traces.GCaMP = nu_trace_mean;
    writetable(DA_pico_traces, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'DA_picospritzing_traces')
  
end



cd([]) %Add path as character array




