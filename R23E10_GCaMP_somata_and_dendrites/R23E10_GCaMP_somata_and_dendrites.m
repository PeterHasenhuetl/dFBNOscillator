% Analysis of dFBN GCaMP dynamics, simultaneously recorded from their
% dendrites and somata.
% Code written by Peter Hasenhuetl.

clear all
tic

source_data_cellbody_details.data_path = []; %Add path as character array
source_data_cellbody_details.file_name = 'cellbody_dendrite_recording.xlsx';

pooled_sorted_r2 = [];
data_table = table;

% Loops through individual flies
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
cum_idx = 0;
for fly_idx = 1:n_flies
    fly_input = dir('*fly*.mat');
    loaded_fly = load(char(fly_input(fly_idx).name));
    curr_fly_name = char(fly_input(fly_idx).name);
    loaded_fly.details = 1;

    % Left hemisphere
    if isempty(loaded_fly.dendrites_1) == 0
        
        curr_table = get_individual_recording(curr_fly_name, loaded_fly, ...
            'dendrites_left');
        data_table = [data_table; curr_table];
        curr_size = size(curr_table.img_measures.r2_single_cells);
        curr_padding = NaN(12,curr_size(2));
        curr_padding(1:curr_size(1),:) = curr_table.img_measures.r2_single_cells;
        pooled_sorted_r2 = [pooled_sorted_r2, curr_padding];
        cum_idx = cum_idx+1;

    end
    
    % Right hemisphere
    if isempty(loaded_fly.dendrites_2) == 0

        curr_table = get_individual_recording(curr_fly_name, loaded_fly, ...
            'dendrites_right');
        data_table = [data_table; curr_table];
        curr_size = size(curr_table.img_measures.r2_single_cells);
        curr_padding = NaN(12,curr_size(2));
        curr_padding(1:curr_size(1),:) = curr_table.img_measures.r2_single_cells;
        pooled_sorted_r2 = [pooled_sorted_r2, curr_padding];        
        cum_idx = cum_idx+1;

    end

end

R2_full = [];
for loop_idx = 1:size(data_table,1)
    R2_full = [R2_full, data_table.img_measures(loop_idx,1).r2_full];
end

curr_idx1 = 1;
for loop_idx = 1:size(data_table,1)
    curr_hist =  data_table.img_measures(loop_idx,1).scatter_pred_hist;
    curr_idx2 = curr_idx1+size(curr_hist,3)-1;
    scatter_hist_all(:,:,curr_idx1:curr_idx2) = curr_hist;
    curr_idx1 = curr_idx1+size(curr_hist,3);
end
    
R2_array = [];    
model_size_array = [];
for loop_idx = 1:size(data_table,1)
    curr_R2 = data_table.img_measures(loop_idx,1).r2_full;
    curr_model_size = ones(1,size(curr_R2,2)).*size(data_table.img_measures(loop_idx,1).cell_bodies,2);
    R2_array = [R2_array, mean(curr_R2)];
    model_size_array = [model_size_array, mean(curr_model_size)];
end

%% Plots the data

session_idx = 6;
close all
figure('Name','Fig_cb_dendrites', 'Color','White','Units','centimeters','Position',[30, 10, 18.3, 20])
color = get_color;
get_fig_panel_cb_dendrite_traces(1, 8, data_table, color, session_idx, source_data_cellbody_details)
get_fig_panel_data_vs_model(data_table.img_measures(1,1).x_edges, scatter_hist_all, 1.75, 3.5, n_flies, source_data_cellbody_details)
get_fig_panel_model_performance(7, 0.5, pooled_sorted_r2, color, source_data_cellbody_details)
get_fig_panel_model_size_vs_performance(1, 0.5, model_size_array, R2_array, color, source_data_cellbody_details)
cd(source_data_cellbody_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'cell_body_vs_dendrite_fig.pdf')
cd([]) %Add path as character array

%% Custom analysis functions

function data_table = get_individual_recording(curr_fly_name, loaded_fly, dend_ID)

if dend_ID == "dendrites_left"
    cell_bodies = loaded_fly.cell_bodies_1;
    dendritic_signal = loaded_fly.dendrites_1;
    loop_size = size(loaded_fly.dendrites_1,2);
elseif dend_ID == "dendrites_right"
    cell_bodies = loaded_fly.cell_bodies_2;
    dendritic_signal = loaded_fly.dendrites_2;
    loop_size = size(loaded_fly.dendrites_2,2);
end

cell_bodies = cell_bodies(6:end,:);
dendritic_signal = dendritic_signal(6:end,:);
cell_bodies = (cell_bodies-movmean(cell_bodies,500))./movmean(cell_bodies,500);
cell_bodies = zscore(cell_bodies);
dendritic_signal = (dendritic_signal-movmean(dendritic_signal,500))./movmean(dendritic_signal,500);
dendritic_signal = zscore(dendritic_signal);
cell_bodies = smoothdata(cell_bodies,'Gaussian',10);
dendritic_signal = smoothdata(dendritic_signal,'Gaussian',10);

% Loops through dendritic ROIs
full_model_weights = [];
for loop_idx = 1:loop_size

    [Y_fit_full(:,loop_idx), r2_full(loop_idx), ...
        r2_single_cells(:,loop_idx), r2_single_cells_normalized(:,loop_idx), ...
        sort_idx1(:,loop_idx), sort_idx2(:,loop_idx), mdl_coeff_full,...
        scatter_pred_hist(:,:,loop_idx), x_edges, ~] = ...
        get_img_measures_cb_dendrites(cell_bodies, dendritic_signal(:,loop_idx));

    full_model_weights(1:length(mdl_coeff_full),loop_idx) = mdl_coeff_full;

end

output_img = struct;
output_img.Y_fit_full = Y_fit_full;
output_img.r2_full = r2_full;
output_img.sort_idx1 = sort_idx1;
output_img.r2_single_cells = r2_single_cells;
output_img.sort_idx2 = sort_idx2;
output_img.cell_bodies = cell_bodies;
output_img.dendrites = dendritic_signal;
output_img.full_model_weights = full_model_weights;
output_img.scatter_pred_hist = scatter_pred_hist;
output_img.x_edges = x_edges;

data_table = table;
data_table.fly_name = string(curr_fly_name);
data_table.dend_ID = string(dend_ID);
data_table.img_measures = output_img;

end

function [Y_fit_full, r2_full, r2_single_cells, r2_single_cells_normalized, ...
    sort_idx1, sort_idx2, mdl_coeff_full, scatter_pred_hist, x_edges, y_edges] = ...
    get_img_measures_cb_dendrites(cell_bodies, dendritic_signal)


[~, r2_full, cb_dend_corr] = get_cross_validated_model(cell_bodies, dendritic_signal);

mdl = fitlm(cell_bodies, dendritic_signal);
Y_fit_full = mdl.Fitted;
mdl_coeff_full = mdl.Coefficients.Estimate(2:end,1);
scatter_pred_hist = cb_dend_corr.scatter_pred_hist;
x_edges = cb_dend_corr.x_edges;
y_edges = cb_dend_corr.y_edges;

r2_single_cells = NaN(size(cell_bodies,2),1);
r2_single_cells_normalized = NaN(size(cell_bodies,2),1);

for loop_idx = 1:size(cell_bodies,2)

    curr_X = cell_bodies(:,loop_idx);
    [~, r2_single_cells(loop_idx,1), ~] = get_cross_validated_model(curr_X, dendritic_signal);
    r2_single_cells_normalized(loop_idx,1) = r2_single_cells(loop_idx,1)/r2_full;

end

[r2_single_cells, sort_idx2] = sort(r2_single_cells,'descend');
[r2_single_cells_normalized, sort_idx1] = sort(r2_single_cells_normalized,'descend');

end

function [R2_in_sample, R2_CV, cb_dend_corr] = get_cross_validated_model(cell_bodies, dendritic_signal)

cv_window = 499;
cb_dend_corr.cb_pred = NaN(cv_window,8,5);
cb_dend_corr.r_2 = NaN(5,1);
cb_dend_corr.dend_pred = NaN(cv_window,5);
cb_dend_corr.dend_to_pred = NaN(cv_window,5);
r2_in_sample = NaN(5,1);
cv_idx = 1:cv_window:length(dendritic_signal);  

for loop_idx = 1:5 % Five-fold cross-validation
            
    cv_interval = cv_idx(loop_idx):(cv_idx(loop_idx)+(cv_window-1));
    X = cell_bodies;
    Y = dendritic_signal;

    X(cv_interval,:) = []; % Cuts out current cv-part from trace
    Y(cv_interval,:) = []; % Cuts out current cv-part from trace
    mdl = fitlm(X, Y); % Fits a liner model
    Y_model_pred = predict(mdl,X);    
    r_in_sample = corrcoef(Y,Y_model_pred);    
    r2_in_sample(loop_idx,1) = r_in_sample(2).^2; 

    X_pred = cell_bodies(cv_interval,:); % Trace for testing current model
    Y_data_to_pred = dendritic_signal(cv_interval,1); % Trace for testing current model
    Y_model_pred = predict(mdl,X_pred);    
    r = corrcoef(Y_data_to_pred,Y_model_pred);    
    r_2 = r(2).^2; % R-squared of current iteration of cross-validation
    cb_dend_corr.r_2(loop_idx,1) = r_2;    
    cb_dend_corr.cb_pred(:,1:size(X_pred,2),loop_idx) = X_pred;    
    cb_dend_corr.dend_pred(:,loop_idx) = Y_model_pred;    
    cb_dend_corr.dend_to_pred(:,loop_idx) = Y_data_to_pred;

    x_data = zscore(Y_data_to_pred);
    y_model = zscore(Y_model_pred);
    
    % For quality check of 2-D histogram orientation (not to mix up x- and y-axes)
    quality_check = "off";
    if quality_check == "on"
        x_data(1:(round(length(x_data)/10))) = 0;
        y_model(1:(round(length(x_data)/10))) = 2;
    end

    % Computes bivariate histogram of z-scored data vs model.
    [scatter_pred_hist(:,:,loop_idx), cb_dend_corr.x_edges, cb_dend_corr.y_edges] = ...
        histcounts2(x_data, y_model, 'NumBins', 15,...
        'XBinLimits', [-3, 3], 'YBinLimits', [-3, 3], 'Normalization', 'probability');
    
end

cb_dend_corr.scatter_pred_hist = mean(scatter_pred_hist,3);
R2_in_sample = mean(r2_in_sample);
R2_CV = mean(cb_dend_corr.r_2);

end

%% Plotting functions

function get_fig_panel_model_performance(x_pos, y_pos, pooled_sorted_r2, color, source_data_details)

sz_1 = 4;
ht_1 = 1.8;
indiv_color = color.light_gray;
mean_color = [0, 0, 0];
xlim_1 = [1, 9];
ylm_1 = [0, 1];

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl, pooled_sorted_r2, 'Color', indiv_color, 'LineWidth', 0.5)
plot(curr_pnl, mean(pooled_sorted_r2,2,'omitnan'), 'Color', mean_color, 'LineWidth', 1.5)
get_default_separated_ax(curr_pnl, xlim_1(1), xlim_1(2),xlim_1(1), xlim_1(2), ylm_1(1), ylm_1(2), ...
    ylm_1(1), ylm_1(2), 1, 0.25, "linear", "linear", 'Cell # (sorted)', 'R^2',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

if isempty(source_data_details) == 0

    R2_array = table;
    R2_array.cell_number = (1:length(pooled_sorted_r2)).';
    R2_array.sorted_r2 = pooled_sorted_r2.';
    writetable(R2_array, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'cell_performance')
end

end

function get_fig_panel_cb_dendrite_traces(x_pos, y_pos, data_table, color, session_idx, source_data_details)

dendrite_idx = 1;
curr_table = data_table(session_idx,:);
sort_idx = curr_table.img_measures.sort_idx2(:,dendrite_idx);

sz_1 = 4;
ht_2 = 0.2;
sz_bar = 0.4;
x_pos1 = x_pos+sz_1+2;
h_t1 = length(sort_idx)*ht_2;
y_pos1 = y_pos-h_t1;
h_t3 = h_t1*1;
h_t4 = h_t3*0.75;
lw_trace1 = 0.5;
lw_trace2 = 0.75;
ht_autocorr = 1;
sz_autocorr = 1;
x_dist_bar = (sz_1+sz_bar+0.1)*(-1);
x_dist_autocorr = 0.5;
xlim_auto = [-5, 5];
ylm_auto2 = [-0.25, 1];
idx_1 = 901;
idx_2 = idx_1+700;

contr_model_perf = curr_table.img_measures.r2_single_cells(:,dendrite_idx);
X = curr_table.img_measures.cell_bodies;
X_nu = X(:,sort_idx);
Y_to_pred = curr_table.img_measures.dendrites(:,dendrite_idx);
Y_fit = curr_table.img_measures.Y_fit_full(:,dendrite_idx);

Y_to_pred_expanded = Y_to_pred(idx_1:idx_2,1);
Y_fit_expanded = Y_fit(idx_1:idx_2,1);
cell_1_expanded = X_nu(idx_1:idx_2,1);
cell_2_expanded = X_nu(idx_1:idx_2,2);

autocorr_cell1 = X_nu(:,1);
autocorr_cell2 = X_nu(:,2);

%% Plots somata

ylm_1 = [min([cell_1_expanded;cell_2_expanded])-0.1, ...
    max([cell_1_expanded;cell_2_expanded])+0.1];
get_default_annotation(x_pos, y_pos+0.2, 'Somata', [0, 0, 0], 'normal', "left")
for loop_idx1 = 1:size(X_nu,2)
  
    curr_ylim = [min(X_nu(:,loop_idx1))-(0.1*min(X_nu(:,loop_idx1))),...
        max(X_nu(:,loop_idx1))+(0.1*max(X_nu(:,loop_idx1)))];
    pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos-(loop_idx1*(ht_2)), sz_1, ht_2]);
    plot(X_nu(:,loop_idx1), 'Color', color.gray, 'LineWidth', lw_trace1)    
    if loop_idx1 == size(X_nu,2)
        hold on
        plot(pnl_1, [length(X_nu(:,loop_idx1))-(14.56*20),...
        length(X_nu(:,loop_idx1))],...
        [min(X_nu(:,loop_idx1))-(0.1*min(X_nu(:,loop_idx1))), min(X_nu(:,loop_idx1))-(0.1*min(X_nu(:,loop_idx1)))],...
        'k', 'LineWidth', get_default_scale_bar_width)
        get_default_annotation(x_pos+sz_1, y_pos-(loop_idx1*(ht_2))-0.12, '20 s', 'k', 'normal', "right")
        get_default_annotation_rotated(x_pos+sz_1+0.2, y_pos-(loop_idx1*(ht_2))-0.12, '1 a.u.', 'k','normal', "left")
        plot(pnl_1, [length(X_nu(:,loop_idx1)), length(X_nu(:,loop_idx1))],...
            [curr_ylim(1), 1*(curr_ylim(1)+(curr_ylim(2)-curr_ylim(1)))],...
        'k', 'LineWidth', get_default_scale_bar_width)
    end
    pnl_1.Color = 'none';
    pnl_1.Box = 'off';
    pnl_1.YAxis.Color = 'none';
    pnl_1.XAxis.Color = 'none';
    ylim(curr_ylim)
    xlim([1, length(X_nu(:,1))])
end

%% Plots R2 values

bar_color = color.gray;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist_bar, y_pos-(loop_idx1*(ht_2)),...
    sz_bar, loop_idx1*(ht_2)]);
barh(pnl_1, 1:length(contr_model_perf), contr_model_perf, 0.5, "EdgeColor", "none",...
    "FaceColor", bar_color, "ShowBaseLine", "off", "LineWidth", 0.5);
ylim([0.5, length(contr_model_perf)+0.5])
xlim([0, max(contr_model_perf)])
pnl_1.Box = 'off';
pnl_1.XAxis.Color = 'none';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Direction = 'reverse';
pnl_1.YAxis.Direction = 'reverse';
pnl_1.XAxisLocation = 'top';

%% Plots the markings of inset

y_add = 0.45;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos1, y_pos1-y_add, sz_1, h_t1+y_add]);
hold on
plot(pnl_1, [1, idx_1, idx_1], [0, y_add, h_t1+y_add], 'k:')
plot(pnl_1, [idx_2, idx_2, length(Y_to_pred)], [h_t1+y_add, y_add, 0], 'k:')
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
xlim([1, length(Y_to_pred)])
ylim([0, h_t1+y_add])

%% Plots dendritic trace plus model

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos1, y_pos1, sz_1, h_t1]);
hold on      
plot(Y_to_pred, 'Color', color.gray, 'LineWidth', lw_trace1)        
plot(Y_fit, 'Color', color.navy, 'LineWidth', lw_trace1)        
plot(pnl_1, [length(Y_fit)-(14.56*20), length(Y_fit)], [min(Y_to_pred)-0.1, min(Y_to_pred)-0.1],...
    'k', 'LineWidth', get_default_scale_bar_width)
plot(pnl_1, [length(Y_fit), length(Y_fit)], [min(Y_to_pred), min(Y_to_pred)+2],...
    'k', 'LineWidth', get_default_scale_bar_width)
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
ylim([min(Y_to_pred)-0.1, max(Y_to_pred)+0.1])
xlim([1, length(Y_to_pred)])
get_default_annotation(x_pos1+sz_1, y_pos1-0.12, '20 s', 'k', 'normal', "right")

get_default_annotation(x_pos1, y_pos1+h_t1+0.2, ...
     'Dendrites', [0, 0, 0], 'normal', "left")
get_default_annotation(x_pos1+sz_1-0.75, y_pos1+h_t1+0.2, ...
     'Data', color.gray, 'normal', "right")
get_default_annotation(x_pos1+sz_1, y_pos1+h_t1+0.2, ...
     'Model', color.navy, 'normal', "right")

%% Plots expanded dendritic trace plus model

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos1, y_pos1-h_t3-0.2, sz_1, h_t3]);
hold on      
plot(pnl_1, Y_to_pred_expanded, 'Color', color.gray, 'LineWidth', lw_trace2)        
plot(pnl_1, Y_fit_expanded, 'Color', color.navy, 'LineWidth', lw_trace2)
plot(pnl_1, [length(Y_fit_expanded), length(Y_fit_expanded)], [min(Y_to_pred), min(Y_to_pred)+2],...
    'k', 'LineWidth', get_default_scale_bar_width)
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
ylim([min(Y_to_pred_expanded)-0.1, max(Y_to_pred_expanded)+0.1])
xlim([1, length(Y_to_pred_expanded)])

%% Plots top 2 cell bodies

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos1, y_pos1-h_t3-1.5, sz_1, h_t4]);
hold on      
plot(pnl_1, cell_2_expanded, 'Color', color.light_gray, 'LineWidth', lw_trace2)
plot(pnl_1, cell_1_expanded, 'Color', color.dark_gray, 'LineWidth', lw_trace2) 
plot(pnl_1, [length(Y_to_pred_expanded)-(14.56*5), length(Y_to_pred_expanded)],...
    [ylm_1(1), ylm_1(1)], 'k', 'LineWidth', get_default_scale_bar_width)
plot(pnl_1, [length(Y_fit_expanded), length(Y_fit_expanded)],...
    [min(Y_to_pred), min(Y_to_pred)+2],...
    'k', 'LineWidth', get_default_scale_bar_width)
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
ylim(ylm_1)
xlim([1, length(Y_to_pred_expanded)])
get_default_annotation(x_pos1+sz_1, y_pos1-h_t3-1.5-0.12, '5 s', 'k', 'normal', "right")

get_default_annotation(x_pos1+(sz_1*0.5)-1, y_pos1-h_t3-1.5-0.12, ...
     'Soma 2', color.light_gray, 'normal', "left")
get_default_annotation(x_pos1, y_pos1-h_t3-1.5-0.12, ...
     'Soma 1', color.dark_gray, 'normal', "left")

%% Plots autocorrelation of dendritic signal

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist_autocorr, ...
    y_pos1-h_t3-0.0, sz_autocorr, ht_autocorr]);
hold on
plot(pnl_1, [-10, 10], [0, 0], 'k:')
[a, b] = xcov(Y_to_pred,'coef');
plot(pnl_1, b./14.56, a, 'Color', color.gray)
[a, b] = xcov(Y_fit,'coef');
plot(pnl_1, b./14.56, a, 'Color', color.navy)
plot(pnl_1, [xlim_auto(2)-2, xlim_auto(2)], [ylm_auto2(1), ylm_auto2(1)], 'k', 'LineWidth', get_default_scale_bar_width)
plot(pnl_1, [xlim_auto(1), xlim_auto(1)], [ylm_auto2(2)-0.5, ylm_auto2(2)-0.2], 'k', 'LineWidth', get_default_scale_bar_width)
xlim(xlim_auto)
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
get_default_annotation_rotated(x_pos+sz_1+x_dist_autocorr-0.2, y_pos1-h_t3+0.755, '0.3 r', 'k', 'normal', "left")
get_default_annotation(x_pos+sz_1+x_dist_autocorr+sz_autocorr, y_pos1-h_t3-0.12, '2 s', 'k', 'normal', "right")

%% Plots autocorrelations of somata

pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+x_dist_autocorr, ...
    y_pos1-h_t3-1.5, sz_autocorr, ht_autocorr]);
hold on
plot(pnl_1, [-10, 10], [0, 0], 'k:')
[a, b] = xcov(autocorr_cell2,'coef');
plot(b./14.56, a, 'Color', color.light_gray)
[a, b] = xcov(autocorr_cell1,'coef');
plot(pnl_1, b./14.56, a, 'Color', color.dark_gray)
plot(pnl_1, [xlim_auto(2)-2, xlim_auto(2)], [ylm_auto2(1), ylm_auto2(1)], 'k', 'LineWidth', get_default_scale_bar_width)
plot(pnl_1, [xlim_auto(1), xlim_auto(1)], [ylm_auto2(2)-0.5, ylm_auto2(2)-0.2], 'k', 'LineWidth', get_default_scale_bar_width)
xlim(xlim_auto)
ylim(ylm_auto2)
pnl_1.Color = 'none';    
pnl_1.Box = 'off';
pnl_1.YAxis.Color = 'none';
pnl_1.XAxis.Color = 'none';
get_default_annotation(x_pos+sz_1+x_dist_autocorr+sz_autocorr, y_pos1-h_t3-1.5-0.12, '2 s', 'k', 'normal', "right")
get_default_annotation_rotated(x_pos+sz_1+x_dist_autocorr-0.2, y_pos1-h_t3-1.5+0.755, '0.3 r', 'k', 'normal', "left")

if isempty(source_data_details) == 0

    data_and_model_traces = table;
    data_and_model_traces.img_trace_dendrite = Y_to_pred;
    data_and_model_traces.img_trace_prediction = Y_fit;
    writetable(data_and_model_traces, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'dendrite_data_and_model')

    somatic_recordings = table;
    somatic_recordings.somatic_trace = X_nu;
    writetable(somatic_recordings, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'img_traces_somata')
  
end

end

function get_fig_panel_data_vs_model(x_edges, scatter_hist_all, x_pos, y_pos, n_flies, source_data_details)


s_z = 2.2;
f_s = get_default_font_size;
xy_vec = x_edges(2:end)-(mean(diff(x_edges))/2); 
xlm_1 = [-3, 3];
ylm_1 = [-3, 3];

RGB_1 = [0, 0, 0];
col_1(:,1) = linspace(1,RGB_1(1),10000);
col_1(:,2) = linspace(1,RGB_1(2),10000);
col_1(:,3) = linspace(1,RGB_1(3),10000);
curr_colormap = [col_1(:,1), col_1(:,2), col_1(:,3)];

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

clim_1 = [0, max(max(mean(scatter_hist_all,3))).*100];
clrb_r = colorbar(curr_pnl,'eastoutside');
clrb_r.Units = 'centimeters';
clrb_r.Position = [x_pos+s_z+0.1, y_pos,  0.2, s_z];
ylabel(clrb_r,'Percent')
clrb_r.FontSize = f_s;
clrb_r.Color = 'k';
clrb_r.TickLength = 0;
clrb_r.Limits = clim_1;
clrb_r.YLabel.Visible = 'on';
clrb_r.YLabel.Color = 'k';

%% Adds annotation

get_default_annotation(x_pos+s_z, y_pos+0.2, [num2str(n_flies), ' flies'], 'k', 'normal', "right")

%% Saves the source data

if isempty(source_data_details) == 0

    data_vs_model = table;
    data_vs_model.data_vec = [NaN; xy_vec'];
    data_vs_model.model_vec = [xy_vec; m_scat_hist];
    writetable(data_vs_model, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'data_vs_model')

end

end

function get_fig_panel_model_size_vs_performance(x_pos, y_pos, model_size_array, R2_array, color, source_data_details)

scatter_dot_color = [0, 0, 0];
scatter_dot_size = 7;
scatter_dot_transparency = 0.5;
sz_1 = 4;
ht_1 = 1.8;
line_color = [0, 0, 0];
xlim_1 = [3, 9];
ylm_1 = [0, 1];
y_dist1 = 0.2;

lin_mdl = fitlm(model_size_array, R2_array);
line_slope = lin_mdl.Coefficients.Estimate(2);

curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
hold on
plot(curr_pnl, model_size_array, lin_mdl.Fitted, 'Color', line_color, 'LineWidth', 0.5)
scatter(curr_pnl, model_size_array, R2_array,scatter_dot_size, 'MarkerFaceColor', scatter_dot_color,...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', scatter_dot_transparency)
get_default_separated_ax(curr_pnl, xlim_1(1), xlim_1(2), xlim_1(1), xlim_1(2), ylm_1(1), ylm_1(2), ...
    ylm_1(1), ylm_1(2), 1, 0.25, "linear", "linear", 'Model size ({\itn} cells)', 'R^2',...
    'k', 'k', 'k', 'k', sz_1, ht_1)

[a, ~] = corrcoef(model_size_array, R2_array);
get_default_annotation(x_pos+sz_1, y_pos+y_dist1, ['r = ', num2str(round(a(2),2))], 'k', 'normal', "right")
get_default_annotation(x_pos, y_pos+y_dist1, ['Slope: ', num2str(round(line_slope,2))], 'k', 'normal', "left")

for loop_idx = 1:10
    n_array(loop_idx) = sum(model_size_array == loop_idx);
end

n_array = n_array(3:9);

hist_color = color.medium_gray;
bin_lim = [min(R2_array), max(R2_array)];
[a_1, a_2] = histcounts(R2_array, 'Normalization', 'percentage', 'NumBins', 10, 'BinLimits', bin_lim);
plot_vec = a_2+mean(diff(a_2))/2;
plot_vec = plot_vec(1:end-1);
pnl_hist = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos, 0.4, ht_1]);
bar(pnl_hist, plot_vec, a_1, 1, 'FaceColor', hist_color, 'FaceAlpha', 1,...
    'EdgeColor', hist_color, 'LineWidth', 0.5)
hold on
plot(pnl_hist, [0.9, 0.9], [0, 23], 'k', 'LineWidth', get_default_scale_bar_width)
scatter(pnl_hist, mean(R2_array), 20, scatter_dot_size, 'MarkerFaceColor', color.red,...
    'MarkerEdgeColor', 'none', 'Marker', '<')
xlim([0, 1])
ylim([0, 30])
pnl_hist.XAxis.Color = 'none';
pnl_hist.YAxis.Color = 'none';
set(pnl_hist, 'view', [90 -90]);

get_default_annotation(x_pos+sz_1+0.1, y_pos+ht_1+0.05, ...
    '20%', 'k', 'normal', "left")



if isempty(source_data_details) == 0

    R2_array_vs_n_cells = table;
    R2_array_vs_n_cells.model_size_array = model_size_array.';
    R2_array_vs_n_cells.R_squared = R2_array.';
    writetable(R2_array_vs_n_cells, [source_data_details.data_path, source_data_details.file_name],...
        'Sheet', 'R2_vs_n_cells')

end



end

