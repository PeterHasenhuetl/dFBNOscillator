% Analysis of dFBN GCaMP dynamics during behavior.
% Code written by Peter Hasenhuetl.

clear all
tic

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'R23E10_GCaMP_on_ball_data.xlsx';

% go into folder with extracted traces and get fly names
cd([]) %Add path as character array

ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
behavior_train_data_table = [];
img_train_data_table = [];
cross_corr_table = [];
ad_corr_table = [];
inter_hemi_dm_AP = [];
inter_hemi_dm_AM = [];
ST_IMG_AP = [];
ST_IMG_AM = [];
SWA_trace_AP = [];
SWA_trace_AM = [];
dF_AP = [];
dF_AM = [];
cons_idx1 = 8;
cons_idx2 = 28;

% loop through individual flies
n_1 = 1;
n_2 = 1;
for fly_idx = 1:n_flies
    
     cd([]) %Add path as character array
        
    fly_input = dir('*fly*.mat'); 
    loaded_fly = load(char(fly_input(fly_idx).name));
           
    loaded_fly.fly_details.exp_group = string(loaded_fly.fly_name_imaging(1,end-1:end));
    
    [dm_behav_row, ind_data_table, curr_st_trace] = ...
        get_individual_dFB_fly_behavior(loaded_fly, loaded_fly.fly_details, fly_idx);    
    
    behavior_train_data_table = [behavior_train_data_table; dm_behav_row];
    img_train_data_table = [img_train_data_table; ind_data_table];

    interhemi_corr_data = get_inter_hemispheric_correlation_behav(ind_data_table);
    if isstruct(interhemi_corr_data) == 1 && loaded_fly.fly_details.exp_group == "AP"
        inter_dm_row_AP = get_inter_corr_table_row_behav(interhemi_corr_data, fly_idx, loaded_fly.fly_details);    
        inter_hemi_dm_AP = [inter_hemi_dm_AP; inter_dm_row_AP];
        ipsi_pre_AP(:,n_1) = mean(interhemi_corr_data.al_trc_pre_ipsi,2);
        contra_pre_AP(:,n_1) = mean(interhemi_corr_data.al_trc_pre_contra,2);
        ipsi_post_AP(:,n_1) = mean(interhemi_corr_data.al_trc_post_ipsi,2);
        contra_post_AP(:,n_1) = mean(interhemi_corr_data.al_trc_post_contra,2);
        
        amp_ipsi_pre_AP(n_1,1) = mean(ipsi_pre_AP(cons_idx1:cons_idx2,n_1));
        amp_contra_pre_AP(n_1,1) = mean(contra_pre_AP(cons_idx1:cons_idx2,n_1));
        amp_ipsi_post_AP(n_1,1) = mean(ipsi_post_AP(cons_idx1:cons_idx2,n_1));
        amp_contra_post_AP(n_1,1) = mean(contra_post_AP(cons_idx1:cons_idx2,n_1));
        
        n_1 = n_1+1;
        
    elseif isstruct(interhemi_corr_data) == 1 && loaded_fly.fly_details.exp_group == "AM"
        inter_dm_row_AM = get_inter_corr_table_row_behav(interhemi_corr_data, fly_idx, loaded_fly.fly_details);    
        inter_hemi_dm_AM = [inter_hemi_dm_AM; inter_dm_row_AM];
        ipsi_pre_AM(:,n_2) = mean(interhemi_corr_data.al_trc_pre_ipsi,2);
        contra_pre_AM(:,n_2) = mean(interhemi_corr_data.al_trc_pre_contra,2);
        ipsi_post_AM(:,n_2) = mean(interhemi_corr_data.al_trc_post_ipsi,2);
        contra_post_AM(:,n_2) = mean(interhemi_corr_data.al_trc_post_contra,2);
        
        amp_ipsi_pre_AM(n_2,1) = mean(ipsi_pre_AM(cons_idx1:cons_idx2,n_2));
        amp_contra_pre_AM(n_2,1) = mean(contra_pre_AM(cons_idx1:cons_idx2,n_2));
        amp_ipsi_post_AM(n_2,1) = mean(ipsi_post_AM(cons_idx1:cons_idx2,n_2));
        amp_contra_post_AM(n_2,1) = mean(contra_post_AM(cons_idx1:cons_idx2,n_2));
        
        n_2 = n_2+1;
    end
    
    if loaded_fly.fly_details.exp_group == "AP"
        ST_IMG_AP = [ST_IMG_AP, curr_st_trace];
        mean_dF_trace_full = get_mean_dF_trace_full(ind_data_table);
        dF_AP = [dF_AP, mean_dF_trace_full];
        
    elseif loaded_fly.fly_details.exp_group == "AM"
        ST_IMG_AM = [ST_IMG_AM, curr_st_trace];
        mean_dF_trace_full = get_mean_dF_trace_full(ind_data_table);
        dF_AM = [dF_AM, mean_dF_trace_full];
    end
    
end

cd([]) %Add path as character array
toc


%% plotting

color = get_color;
close all
figure('Name','Fig_GCaMP_ball', 'Color','White','Units','centimeters','Position',[30, 10, 18.3, 20])
y_pos = 2;
x_pos = 2;
color_pre = [0, 0, 0];
color_post = color.medium_gray+0.1;
get_fig_panel_power_spectrum_behavior(img_train_data_table, x_pos, y_pos+3, 'AP', ...
    color_pre, color_post, source_data_details)

coord_1 = 10;
coord_2 = 4;
coord_3 = 3;
get_fig_panel_punish_run(behavior_train_data_table, img_train_data_table, ...
    color, x_pos+1.5, 7, 'AP', coord_1, coord_2, coord_3, dF_AP);

plotting_color_ipsi = color.medium_gray;
plotting_color_contra = color.medium_gray;
get_fig_panel_transient_aligned_contra(ipsi_pre_AP, contra_pre_AP, ...
    x_pos, y_pos, plotting_color_ipsi, plotting_color_contra,...
    [-0.1, 0.45], [-0.12, 0.05], 1)

plotting_color_ipsi = color.medium_gray;
plotting_color_contra = color.medium_gray;
get_fig_panel_transient_aligned_contra(ipsi_post_AP, contra_post_AP, ...
    x_pos+2, y_pos, plotting_color_ipsi, plotting_color_contra,...
    [-0.1, 0.45], [-0.12, 0.05], 0)

cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'R23E10_GCaMP_on_ball_fig.pdf')
cd([]) %Add path as character array

%% custom functions called in this script

function [dm_behav_row, ind_data_table, curr_st_trace] = get_individual_dFB_fly_behavior(input_fly, ...
    fly_details, fly_idx)

input_mov = input_fly.curr_pos;
TTL_input = table2array(input_fly.TTL_pun);

pre_pun_sec = 150;
post_pun_sec = 400;

if fly_details.exp_group == "AP" %Punishment
    input_pwm = input_fly.curr_pwm(:,[1,6]);
    fly_details.b_l = input_pwm(1,2);
elseif fly_details.exp_group == "AM" %Mock
    input_pwm = input_fly.curr_pwm(:,[1,6]);
    fly_details.b_l = input_pwm(1,2);
else 
    fly_details.b_l = 2;
    input_pwm = input_fly.curr_pwm(:,[1,2]);
end

fly_details.frame_rate_img = 14.56;
fly_details.m = 10;
fly_details.test_onset_lag = 40;
fly_details.pre_trial_indices = 9;
fly_details.pre_exclusion_window = 100;
fly_details.pre_running_threshold = 50;
fly_details.running_window = 150;
fly_details.running_threshold = 10;
fly_details.onset_running_threshold = 1;

fly_details.n_trials_train = 3;
fly_details.trial_length_train = 7000; %trial length for train in ms
fly_details.baseline_jump = round(25*fly_details.frame_rate_img); %seconds converted in #frames
fly_details.pre_trial_window_img_train = round(5*fly_details.frame_rate_img); %seconds converted in #frames
fly_details.post_TTL_lag = round(15.5*fly_details.frame_rate_img);
fly_details.idx_triallength = round((fly_details.trial_length_train/1000)*fly_details.frame_rate_img);

%% The actual analysis

% first, behavior
dm_behav_row = get_dFB_behavior_analysis(input_mov, input_pwm, fly_details, fly_idx, pre_pun_sec, post_pun_sec);

% then, imaging
TTLoutputs = get_img_trial_starts(TTL_input, fly_details.n_trials_train);

fields_img_planes = fieldnames(input_fly.imaging_pun);
ind_data_table = [];
for plane_idx = 1:length(fields_img_planes)
    plane_name = fields_img_planes{plane_idx};
    [output_img, TTLoutputs_nu] = ...
        get_dFB_general_img_measures(input_fly.imaging_pun.(string(plane_name)), TTLoutputs, pre_pun_sec, post_pun_sec);
    single_trial_output_img = get_dFB_single_trial_img_responses_behavior(output_img.raw_img, fly_details, TTL_input, TTLoutputs_nu);
    curr_row = get_dFB_behavior_data_table_row(output_img, single_trial_output_img, fly_details, fly_idx, plane_idx);
    ind_data_table = [ind_data_table; curr_row];
end

curr_struct = ind_data_table.st_trace;
n_n = 1;
for loop_idx_st = 1:4
    curr_st_trace(:,n_n:n_n+1,1) = curr_struct(loop_idx_st).trace(:,3:4,1);
    n_n = n_n+2;
end

curr_st_trace(:,isnan(curr_st_trace(1,:)) == 1) = [];
curr_st_trace(:, sum(curr_st_trace(:,:),1) == 0) = [];
curr_st_trace = mean(curr_st_trace,2);

end

function dm_behav_row = get_dFB_behavior_analysis(mov_test1, pwm_test1, ...
       fly_details, fly_idx, pre_pun_sec, post_pun_sec)

curr_mov_test = mov_test1;
curr_pwm_test = pwm_test1;

% sets both time vectors to start from zero
curr_mov_test(:,1) = curr_mov_test(:,1) - curr_mov_test(1,1);
curr_pwm_test(:,1) = curr_pwm_test(:,1) - curr_pwm_test(1,1);

% divides movement trace by 5 (because Matlab-balltracking multiplies it by
% 5 when communicating with LabView)
curr_mov_test(:,2:end) = curr_mov_test(:,2:end)./5;

%Initializes a structure 'output_test' to save the data later on
output_train = struct;

%% calculates ratio between time vectors (mov and pwm) and adjusts them to be the same:

% the factor by which the camera output is higher than the pwm output, should be ~5
s_test = curr_mov_test(end,1)/curr_pwm_test(end,1);

% bin-size of camera in iterations of output from Matlab to Labview
frt_1 = curr_mov_test(end,1)/(length(curr_mov_test)-1);
bin_size_mov = frt_1/s_test; % actual bin-size in milliseconds
curr_mov_test(:,1) = curr_mov_test(:,1)/s_test; % corrects the time-vector of movement to correct time

%% finds starting points of individual trials in pwm-trace and movement-trace

[trialontimes_pwm_test, trialon_indices_position_test, trialontimes_position_test,] = ...
    get_behavior_trial_starts(curr_pwm_test, curr_mov_test, fly_details, fly_details.n_trials_train);

output_train.mov = curr_mov_test;
output_train.PWM = curr_pwm_test;
output_train.trialoff_indices_position = 1;
output_train.trialon_indices_position = trialon_indices_position_test;
output_train.trialontimes_position = trialontimes_position_test;
output_train.trialofftimes_position = 1;
output_train.trialontimes_pwm = trialontimes_pwm_test;

%%

pre_punish_window = round(pre_pun_sec/(bin_size_mov/1000));
post_punish_window = round(post_pun_sec/(bin_size_mov/1000));
cut_off_start = trialon_indices_position_test(1)-pre_punish_window;
cut_off_end = trialon_indices_position_test(1)+post_punish_window;

nu_mov = curr_mov_test(cut_off_start:cut_off_end,:);
nu_mov(:,1) = nu_mov(:,1) - nu_mov(1,1);

%%

dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.exp_group = fly_details.exp_group;
mov = struct;
mov.trace = output_train.mov;
dm_behav_row.mov = mov;

cut_mov = struct;
cut_mov.trace = nu_mov;
dm_behav_row.cut_mov = cut_mov;

P_mov = struct;
P_mov.mov_trace = abs(nu_mov) >= 0.5;
dm_behav_row.P_mov = P_mov;

pwm = struct;
pwm.trace = output_train.PWM;
dm_behav_row.pwm = pwm;
dm_behav_row.bin_size_mov = bin_size_mov;

end

function [trialontimes_pwm, trialon_indices_movement, trialontimes_movement] = get_behavior_trial_starts(pwm_test, ...
    mov_test, fly_details, n_trials)


trial_starts_test = zeros(length(pwm_test),1); 

% creates vector of differentiated pwm signal
dts_test = [0; diff(pwm_test(:,2))];

% Creates a vector that sepcifies where the pwm-signal jumped from baseline
% (bl, usually 2) to the value in that trial (e.g. 10)
tr_test = (pwm_test(:,2)) == (dts_test+fly_details.b_l);


for loop_idx_1 = 1:length(pwm_test)
    
    % ignore time bins before the pre-specified 'test_onset_lag'
    if loop_idx_1 < fly_details.test_onset_lag
        
        trial_starts_test(loop_idx_1) = 0;
        
    else
        
         tt_test(loop_idx_1) = any(dts_test(loop_idx_1-fly_details.pre_trial_indices:loop_idx_1-2,1) ~= 0);
          
         trial_starts_test(loop_idx_1-1) = ...
             (tr_test(loop_idx_1) == 0 && tr_test(loop_idx_1-1) == 1 && tt_test(loop_idx_1) == 0);
         
    end
    
end


% Finds the trial onset times according to the trial onset indices.
trialonindices_pwm_test = find(trial_starts_test == 1);
trialontimes_pwm = pwm_test(find(trial_starts_test == 1),1);

% Takes only the trials that are specified by n_trials.
trialontimes_pwm = trialontimes_pwm(1:n_trials,1);

%% Gets trial on-times in movement-trace

% First, gets the indices in the movement trace that correspond to the 
% on-times in the pwm-trace

for loop_idx_2 = 1:numel(trialontimes_pwm)
    
    trialon_indices_movement(loop_idx_2) = ...
        find(abs(mov_test(:,1) - trialontimes_pwm(loop_idx_2)) < 70,1,'first');
    
end

% Gets a vector of the trial-onset times (in time-vector of movement trace) that
% correspond to trial-on-indices
trialontimes_movement = mov_test(trialon_indices_movement,1);

end

function TTLoutputs = get_img_trial_starts(TTL_input, n_trials)

%calculates normalized TTL signal
nTTL = TTL_input(:,1)+500;
nTTL = nTTL-min(nTTL);
nTTL = nTTL./max(nTTL);

%logical vector with 'ones' at indices of TTL onset;
TTLidx = zeros(length(nTTL),1);

for loop_idz1 = 1:length(nTTL)
    if nTTL(loop_idz1,1) < 0.5 
        curr_TTL(loop_idz1,1) = 0;
    else 
        curr_TTL(loop_idz1,1) = 1;
    end
end

for loop_idx2 = 2:(numel(curr_TTL)-1)
    if curr_TTL(loop_idx2,1) == 1 && curr_TTL(loop_idx2-1,1) == 0
        TTLidx(loop_idx2,1) = 1;
    end
end

% creates logical vector which indicates TTL output
TTLoutputs = find(TTLidx == 1, n_trials, 'first'); 

end

function [output_img, TTLoutputs_nu] = get_dFB_general_img_measures(imaging_input, ...
    TTLoutputs, pre_pun_sec, post_pun_sec)

cumulative_idx1 = 1;
cumulative_idx2 = 1; 
f_r = 14.56; % frame-rate
sliding_window = 501;

raw_green = table2array(imaging_input.img_pun(:,1:5));
raw_red = raw_green;
img_green = raw_green(:,2:end) - mean(raw_green(:,1));

pre_punish_window = round(pre_pun_sec*(f_r));
post_punish_window = round(post_pun_sec*(f_r));
cut_off_start = TTLoutputs(1)-pre_punish_window;
cut_off_end = TTLoutputs(1)+post_punish_window;

TTLoutputs_nu = TTLoutputs-cut_off_start;
img_green = img_green(cut_off_start:cut_off_end,:);

for R_idx = 1:size(img_green,2)
    % First, tests for ROI to be included
    
    if raw_red(:,R_idx+1) == raw_red(:,1)
        output_img.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);
    elseif mean(img_green(:,R_idx),1) < 15
        output_img.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);
    else
        output_img.excl_vec(1,R_idx) = 0;
    end
    
    % then, computes dF/F
    raw_F = img_green(:,R_idx);
    
    if output_img.excl_vec(1,R_idx) == 0   
        F_0 = get_mov_prctile_behavior(raw_F, sliding_window);
        d_F = (raw_F - F_0)./(F_0);
        d_F = smoothdata(d_F,'Gaussian',8);
    else
        d_F = zeros(length(raw_F),1);
    end
    
    output_img.dR(:,R_idx) = d_F;
    output_img.raw_img(:,R_idx) = img_green(:,R_idx);
    
    %% ...(1) COMPUTES THE POWER SPECTRUM OF deltaF/F TRACES

    % full-length    
    [one_sided_pwr_sctrm, f] = pspectrum(d_F,14.56);
    output_img.amp(:,R_idx) = one_sided_pwr_sctrm;    
    d_idx1 = find(f > 0.2,1,'first'); 
    d_idx2 = find(f > 1,1,'first');
    output_img.int_SWA_power(R_idx) = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
    d_idx1 = find(f > 0.005,1,'first'); 
    d_idx2 = find(f > 0.1,1,'first');
    output_img.int_infra_slow_power(R_idx) = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
    output_img.f_full = f';
    
    %% ... (2) ANALYZES Ca TRANSIENTS
    
    [output_img.PSTH(:,R_idx),  output_img.burst_signal(:,R_idx)] = ...
        get_GCaMP_transients_behav(d_F);
 
    %% pre vs post

    idx_pre = 1:2100;
    idx_post = 5500:7600;
    dR_pre = d_F(idx_pre,1);
    dR_post = d_F(idx_post,1);

    if sum(dR_pre) == 0
        
        [amp_temp, f_temp, ~] = pspectrum(zeros(length(dR_pre),1), 14.56);
        output_img.f_pre = f_temp';
        output_img.int_SWA_power_pre(R_idx) = NaN;
        output_img.int_infra_slow_power_pre(R_idx) = NaN;
        output_img.amp_pre(:,R_idx) = NaN(length(amp_temp),1);

        output_img.int_SWA_power_post(R_idx) = NaN;
        output_img.int_infra_slow_power_post(R_idx) = NaN;
        output_img.amp_post(:,R_idx) = NaN(length(output_img.amp_pre(:,R_idx)),1);
        output_img.f_post = f_temp';
    
    else
        [output_img.amp_pre(:,R_idx), f_pre, ~] = pspectrum(dR_pre(:,:), 14.56);
        output_img.f_pre = f_pre';
        d_idx1 = find(output_img.f_pre > 0.2,1,'first'); 
        d_idx2 = find(output_img.f_pre > 1,1,'first');
        output_img.int_SWA_power_pre(R_idx) = sum(output_img.amp_pre(d_idx1:(d_idx2-1),R_idx),'omitnan');
        d_idx1 = find(output_img.f_pre > 0.02,1,'first');
        d_idx2 = find(output_img.f_pre > 0.1,1,'first');
        output_img.int_infra_slow_power_pre(R_idx) = sum(output_img.amp_pre(d_idx1:(d_idx2-1),R_idx),'omitnan');
         
        [output_img.amp_post(:,R_idx), f_post, ~] = pspectrum(dR_post(:,:),14.56);
        output_img.f_post = f_post';
        d_idx1 = find(output_img.f_post > 0.2,1,'first');
        d_idx2 = find(output_img.f_post > 1,1,'first');
        output_img.int_SWA_power_post(R_idx) = sum(output_img.amp_post(d_idx1:(d_idx2-1),R_idx),'omitnan');
        d_idx1 = find(output_img.f_post > 0.02,1,'first'); 
        d_idx2 = find(output_img.f_post > 0.1,1,'first');
        output_img.int_infra_slow_power_post(R_idx) = sum(output_img.amp_post(d_idx1:(d_idx2-1),R_idx),'omitnan');
        
    end
        
    %%
   
    p_window = 40;
    n_window = 15;
    if R_idx > 2
        for u_idx = 1:2
            if u_idx == 1
                img = img_green(:,R_idx);
            else
                img = d_F;
            end
            
            burst_indices = find(output_img.burst_signal(:,R_idx) > 0);
            burst_indices((burst_indices+p_window) > length(img) | ...
                (burst_indices-n_window) < 1) =  [];
            
            if isempty(burst_indices) == 1
                dummy_trace = (img(20-n_window:20+p_window,1));
                output_img.trans_trace = NaN(size(dummy_trace,1),1);
                output_img.trans_trace_raw = NaN(size(dummy_trace,1),1);
            end
            
            if u_idx == 1
                for i = 2:length(burst_indices)-5
                    output_img.trans_trace_raw(:,cumulative_idx1) = ...
                        img(burst_indices(i)-n_window:burst_indices(i)+p_window,1);
                    cumulative_idx1 = cumulative_idx1+1;
                end    
            else
                for i = 2:length(burst_indices)-5
                    output_img.trans_trace(:,cumulative_idx2) = ...
                        img(burst_indices(i)-n_window:burst_indices(i)+p_window,1);
                    cumulative_idx2 = cumulative_idx2+1;
                end    
            end    
        end    
    end
end

end

function single_trial_output_img = get_dFB_single_trial_img_responses_behavior(img_input, ...
    details, TTL_input, TTLoutputs)

n_columns = numel(img_input(1,:));
n_tx = numel(TTLoutputs); %specifies number of trials according to number of TTL outputs
single_trial_onset = details.pre_trial_window_img_train+1;

for t_idx = 1:n_tx %loops through trials
    
    rs_1 = TTLoutputs(t_idx)-round(14.56*10);
    %specifies end of imaging window;
    rs_2 = TTLoutputs(t_idx)+round(14.56*25);
    F0_window = rs_1:TTLoutputs(t_idx);
    
    for column_idx = 1:n_columns %loops through columns
       
        F_0 = prctile(img_input(F0_window,column_idx),5);
        st_img(:,column_idx,t_idx) = (img_input(rs_1:rs_2,column_idx)-F_0)/abs(F_0);
        st_TTL(:,column_idx,t_idx) = TTL_input(rs_1:rs_2,1);
        
    end

end

lag_st_time = (details.frame_rate_img * single_trial_onset)/1000;
singletrialtime = ...
    ((0:details.frame_rate_img:((numel(st_img(:,1,1))*details.frame_rate_img)-1))'/1000)-lag_st_time;

%% stores data in structure

single_trial_output_img = struct;
single_trial_output_img.singletrialtime = singletrialtime;
single_trial_output_img.st_img = st_img;
single_trial_output_img.st_TTL = st_TTL;

end

function mean_dF_trace_full = get_mean_dF_trace_full(input_table)
ROI_idx = [3, 4];

n_idx = 1;
for trc_idx = 1:size(input_table,1)
    curr_input = input_table(trc_idx,:);
    curr_trace = curr_input.img_trace.trace;
    
    for ind_idx = 1:2
        if sum(curr_trace(:,ROI_idx(ind_idx))) == 0
            cond_test = 1;
        elseif sum(curr_trace(:,ROI_idx(ind_idx))) > 0
            mean_dF_trace_full(:,n_idx) = curr_trace(:,ROI_idx(ind_idx));
            n_idx = n_idx + 1;
        end
    end
end

mean_dF_trace_full = mean(mean_dF_trace_full,2);

end

function [amplitude_trace,  event_trace] = get_GCaMP_transients_behav(img_trace)

gaussian_sliding_window = 15;
min_peak_height = 0.01;
trans_signal = (1)*smoothdata(diff(img_trace), 'gaussian', gaussian_sliding_window);
[~, trans_ids] = (findpeaks(trans_signal, 'MinPeakHeight', min_peak_height));

s = [img_trace; zeros(20,1)];
ampli_1 = [];
for loop_idx1 = 1:length(trans_ids)
    ampli_1(loop_idx1) = max(s(trans_ids(loop_idx1):trans_ids(loop_idx1)+10));   
end

amplitude_trace = zeros(length(img_trace),1);
amplitude_trace(trans_ids) = ampli_1;

event_trace = zeros(length(img_trace),1);
event_trace(trans_ids) = 1;


end

function inter_hemi_corr = get_inter_hemispheric_correlation_behav(ind_data_table)

dendrite_1 = [];
dendrite_2 = [];
dendrite1_burst = [];
dendrite2_burst = [];
pln_coord1 = [];
pln_coord2 = [];

for loop_idx1 = 1:4
    if ind_data_table.excl_idx(loop_idx1,3) == 0
        dendrite_1 = [dendrite_1, ind_data_table.img_raw(loop_idx1,1).trace(:,3)];
        dendrite1_burst = [dendrite1_burst, ind_data_table.burst_signal(loop_idx1,1).trace(:,3)];
        pln_coord1 = [pln_coord1, loop_idx1];
    end
    
    if ind_data_table.excl_idx(loop_idx1,4) == 0
        dendrite_2 = [dendrite_2, ind_data_table.img_raw(loop_idx1,1).trace(:,4)];
        dendrite2_burst = [dendrite2_burst, ind_data_table.burst_signal(loop_idx1,1).trace(:,4)];
        pln_coord2 = [pln_coord2, loop_idx1];
    end
    
end

idx_pre = 1:2100;
idx_post = 5500:7600;

if isempty(dendrite_1) == 0 && isempty(dendrite_2) == 0
    n = 1;
    for loop_idx2 = 1:size(dendrite_1,2)
    
        for loop_idx3 = 1:size(dendrite_2,2)
           
            inter_hemi_corr.cross_coord(n,1:2) = [pln_coord1(loop_idx2), pln_coord2(loop_idx3)];
            [inter_hemi_corr.al_trc_ipsi, inter_hemi_corr.al_trc_contra] = get_aligned_transients(dendrite_1(:,loop_idx2), ...
                dendrite_2(:,loop_idx3), dendrite1_burst(:,loop_idx2), dendrite2_burst(:,loop_idx3));
            
            % pre and post heat
            [inter_hemi_corr.al_trc_pre_ipsi, inter_hemi_corr.al_trc_pre_contra] = get_aligned_transients(dendrite_1(idx_pre,loop_idx2), ...
                dendrite_2(idx_pre,loop_idx3), dendrite1_burst(idx_pre,loop_idx2), dendrite2_burst(idx_pre,loop_idx3));
            
            [inter_hemi_corr.al_trc_post_ipsi, inter_hemi_corr.al_trc_post_contra] = get_aligned_transients(dendrite_1(idx_post,loop_idx2), ...
                dendrite_2(idx_post,loop_idx3), dendrite1_burst(idx_post,loop_idx2), dendrite2_burst(idx_post,loop_idx3));

            n = n+1;
        end    
    end
    
else
    inter_hemi_corr = NaN;
end


end

function [ipsi_all, contra_all] = get_aligned_transients(curr_trace_left_raw, ...
    curr_trace_right_raw, curr_transient_vec_left, curr_transient_vec_right)
 
p_wind = 40;
n_wind = 10;
F0_wind = round(n_wind);

transient_onsets = find(curr_transient_vec_left > 0);
% Only considers the transient onsets that allow for the full            
% length of transients to be used (i.e., excludes transients at            
% the beginning and end of the imaging trace).          
transient_onsets((transient_onsets+p_wind) > length(curr_trace_left_raw) | ...
    (transient_onsets-n_wind) < 1) =  [];

ipsi_trace_left = NaN(n_wind+p_wind+1,length(transient_onsets));
contra_trace_right = NaN(n_wind+p_wind+1,length(transient_onsets));
for loop_idx1 = 1:length(transient_onsets)
    i_trace_left = (curr_trace_left_raw(transient_onsets(loop_idx1)-n_wind:transient_onsets(loop_idx1)+p_wind,1));
    ipsi_trace_left(:,loop_idx1) = (i_trace_left-mean(i_trace_left(1:F0_wind,1)))/mean(i_trace_left(1:F0_wind,1));
    c_trace_right = (curr_trace_right_raw(transient_onsets(loop_idx1)-n_wind:transient_onsets(loop_idx1)+p_wind,1));
    contra_trace_right(:,loop_idx1) = (c_trace_right-mean(c_trace_right(1:F0_wind,1)))/mean(c_trace_right(1:F0_wind,1));
end

n_transients_left = loop_idx1;
  
transient_onsets = find(curr_transient_vec_right > 0);
% Only considers the transient onsets that allow for the full            
% length of transients to be used (i.e., excludes transients at            
% the beginning and end of the imaging trace). 
transient_onsets((transient_onsets+p_wind) > length(curr_trace_right_raw) | ...
    (transient_onsets-n_wind) < 1) =  [];

ipsi_trace_right = NaN(n_wind+p_wind+1,length(transient_onsets));
contra_trace_left = NaN(n_wind+p_wind+1,length(transient_onsets));
for loop_idx2 = 1:length(transient_onsets)
    i_trace_right = (curr_trace_right_raw(transient_onsets(loop_idx2)-n_wind:transient_onsets(loop_idx2)+p_wind,1));
    ipsi_trace_right(:,loop_idx2) = (i_trace_right-mean(i_trace_right(1:F0_wind,1)))/mean(i_trace_right(1:F0_wind,1));
    c_trace_left = (curr_trace_left_raw(transient_onsets(loop_idx2)-n_wind:transient_onsets(loop_idx2)+p_wind,1));
    contra_trace_left(:,loop_idx2) = (c_trace_left-mean(c_trace_left(1:F0_wind,1)))/mean(c_trace_left(1:F0_wind,1));
end

n_transients_right = loop_idx2;

ipsi_all = [ipsi_trace_left, ipsi_trace_right];
contra_all = [contra_trace_left, contra_trace_right];
ipsi_trace_left = mean(ipsi_trace_left,2);
contra_trace_right = mean(contra_trace_right,2);
ipsi_trace_right = mean(ipsi_trace_right,2);
contra_trace_left = mean(contra_trace_left,2);

end

function mov_prctile = get_mov_prctile_behavior(raw_trace, sliding_window)

% if sliding window is even, adds 1 (for symmetric sliding window)
is_even = (sliding_window/2) == round(sliding_window/2);
if is_even == 1
    sliding_window = sliding_window + 1;
end

% adds a padding window before and after the raw trace
size_padding_window = (sliding_window-1)/2;
padding_window = NaN(size_padding_window,1);
curr_trace = [padding_window; raw_trace; padding_window];

% loops through individual points of trace to create time_varying
% percentile
mov_prctile = NaN(length(raw_trace),1);
window_center_idx = size_padding_window+1;
for loop_idx = 1:length(raw_trace)
    mov_prctile(loop_idx,1) = prctile(curr_trace(window_center_idx-size_padding_window:window_center_idx+size_padding_window),10);    
    window_center_idx = window_center_idx + 1;
end

end

function dm_row = get_dFB_behavior_data_table_row(dm_input, single_trial_output_img, ...
    details, fly_idx, plane_idx)

dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.exp_group = details.exp_group;
dm_behav_row.plane_id = plane_idx;

dm_img_row = table;
dm_img_row.int_SWA_power = dm_input.int_SWA_power;
dm_img_row.int_infra_slow_power = dm_input.int_infra_slow_power;
dm_img_row.int_SWA_power_pre = dm_input.int_SWA_power_pre;
dm_img_row.int_infra_slow_power_pre = dm_input.int_infra_slow_power_pre;
dm_img_row.int_SWA_power_post = dm_input.int_SWA_power_post;
dm_img_row.int_infra_slow_power_post = dm_input.int_infra_slow_power_post;

dm_img_row.f_pre = dm_input.f_pre;
dm_img_row.f_post = dm_input.f_post;
dm_img_row.f_full = dm_input.f_full;

dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;

trc = struct;
trc.trace = dm_input.dR;
dm_img_row.img_trace = trc;

trc_raw = struct;
trc_raw.trace = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;

brst = struct;
brst.trace = dm_input.PSTH;
dm_img_row.burst_signal = brst;

pwr = struct;
pwr.trace = dm_input.amp;
dm_img_row.power_spectrum = pwr;

pwr_pre = struct;
pwr_pre.trace = dm_input.amp_pre;
dm_img_row.power_spectrum_pre = pwr_pre;

pwr_post = struct;
pwr_post.trace = dm_input.amp_post;
dm_img_row.power_spectrum_post = pwr_post;

trans_trace_raw = struct;
trans_trace_raw.trace = dm_input.trans_trace_raw;
dm_img_row.trans_trace_raw = trans_trace_raw;

st_trace = struct;
st_trace.trace = single_trial_output_img.st_img;
dm_img_row.st_trace = st_trace;

%% Concatenates them

dm_row = [dm_behav_row, dm_excl_row, dm_img_row];

end

function dm_row = get_inter_corr_table_row_behav(ihc, fly_idx, details)
dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.exp_group = details.exp_group;

dm_img_row = table;

inpt = NaN(1,16);
in_dt = ihc.cross_coord(:,1)';
inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.coordinates_dendrite1 = inpt;

inpt = NaN(1,16);
in_dt = ihc.cross_coord(:,2)';
inpt(1,1:size(in_dt,2)) = in_dt;
dm_img_row.coordinates_dendrite2 = inpt;

al_trc_ipsi = struct;
al_trc_ipsi.trace = ihc.al_trc_ipsi;
dm_img_row.al_trc_ipsi = al_trc_ipsi;

al_trc_contra = struct;
al_trc_contra.trace = ihc.al_trc_contra;
dm_img_row.al_trc_contra = al_trc_contra;

al_trc_pre_ipsi = struct;
al_trc_pre_ipsi.trace = ihc.al_trc_pre_ipsi;
dm_img_row.al_trc_pre_ipsi = al_trc_pre_ipsi;

al_trc_post_ipsi = struct;
al_trc_post_ipsi.trace = ihc.al_trc_post_ipsi;
dm_img_row.al_trc_post_ipsi = al_trc_post_ipsi;

al_trc_pre_contra = struct;
al_trc_pre_contra.trace = ihc.al_trc_pre_contra;
dm_img_row.al_trc_pre_contra = al_trc_pre_contra;

al_trc_post_contra = struct;
al_trc_post_contra.trace = ihc.al_trc_post_contra;
dm_img_row.al_trc_post_contra = al_trc_post_contra;


%% Concatenates them

dm_row = [dm_behav_row, dm_img_row];

end

%% plotting functions

function get_fig_panel_power_spectrum_behavior(img_train_data_table, x_pos, y_pos, ...
    exp_group, color_pre, color_post, source_data_details)


sz_1 = 1.8;
ht_1 = 1.8;
ht_2 = sz_1;

scale_fac = 10^3;

xlm_pwr = [0, 2];
ylm_pwr = ([0, 0.004])*scale_fac;

curr_tbl = img_train_data_table(img_train_data_table.exp_group == char(exp_group),:);
n_flies = round(size(curr_tbl,1)/4); % because 4 planes per fly
num_array = 1:100;
curr_fly_ids = curr_tbl.fly_id(:,:);
id_array = ismember(num_array,curr_fly_ids);

ind_fly_ids = num_array(id_array == 1);
SP_pre = [];
SP_post = [];

for loop_idx = 1:n_flies
    
    curr_img_tbl = curr_tbl(curr_tbl.fly_id == ind_fly_ids(loop_idx),:);
    [n_pwr_pre, n_pwr_post] = get_n_pwr(curr_img_tbl);
    SP_pre = [SP_pre, mean(n_pwr_pre,2,'omitnan')];
    SP_post = [SP_post, mean(n_pwr_post,2,'omitnan')];
    
end

ac_rate1 = curr_tbl.f_pre(1,:);
ac_rate2 = curr_tbl.f_post(1,:);


curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_2]);
hold on
get_default_SEM_area_plot(curr_pnl, SP_pre*scale_fac, ac_rate1, color_pre)
get_default_SEM_area_plot(curr_pnl, SP_post*scale_fac, ac_rate2, color_post)
get_default_separated_ax(curr_pnl, xlm_pwr(1), xlm_pwr(2), xlm_pwr(1), xlm_pwr(2), ...
    ylm_pwr(1), ylm_pwr(2), ylm_pwr(1), ylm_pwr(2),...
    1, 2, "linear", "linear", 'Frequency (Hz)', {'(\DeltaF/F)^2 x10^{-3}'},...
    'k', 'k', 'k', 'k', sz_1, ht_1)


%% Saves the source data

power_spectrum_pre = table;
power_spectrum_pre.x_values_frequ = ac_rate1';
power_spectrum_pre.power_pre = SP_pre;
writetable(power_spectrum_pre, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'power_spectrum_pre')

power_spectrum_post = table;
power_spectrum_post.x_values_frequ = ac_rate2';
power_spectrum_post.power_post = SP_post;
writetable(power_spectrum_post, [source_data_details.data_path, source_data_details.file_name], 'Sheet', 'power_spectrum_post')

end

function [n_pwr_pre, n_pwr_post] = get_n_pwr(data_tbl)
   
curr_table_pwr_pre = data_tbl.power_spectrum_pre;
curr_table_pwr_post = data_tbl.power_spectrum_post;

curr_pwr_pre = [];
curr_pwr_post = [];

for loop_idx1 = 1:size(data_tbl,1)
    curr_pwr_pre = [curr_pwr_pre, curr_table_pwr_pre(loop_idx1,1).trace(:,3:4)];
    curr_pwr_post = [curr_pwr_post, curr_table_pwr_post(loop_idx1,1).trace(:,3:4)];
end

n = 1;
for loop_idx2 = 1:size(curr_pwr_pre,2)
    if sum(curr_pwr_pre(:,loop_idx2)) > 0
        n_pwr_pre(:,n) = (curr_pwr_pre(:,loop_idx2));
        n_pwr_post(:,n) = (curr_pwr_post(:,loop_idx2));
        n = n + 1;
    end
end

end

function get_fig_panel_transient_aligned_contra(ipsi_trace, contra_trace, ...
    curr_x_pos, curr_y_pos, plotting_color_ipsi, plotting_color_contra,...
    ylm_ipsi, ylm_contra, annotation_cond)

sz_1 = 1.25;
ht_ipsi = 1.5;
ht_contra = ht_ipsi;
dist_1 = abs(ylm_ipsi(1)/(ylm_ipsi(1)-ylm_ipsi(end)));
norm_dist1 = ht_ipsi*dist_1;
dist_2 = abs(ylm_contra(1)/(ylm_contra(1)-ylm_contra(end)));
norm_dist2 = ht_ipsi*dist_2;
xlm_1 = [1,length(mean(ipsi_trace,2))];
t_v = xlm_1(1):xlm_1(2);

scale_bar = [xlm_1(2)-(0.5*14.56), xlm_1(2)];
curr_x_pos1 = curr_x_pos;    
curr_x_pos2 = curr_x_pos;
  
pnl_1 = axes('Units', 'Centimeters', 'Position', [curr_x_pos1, curr_y_pos-norm_dist2, sz_1, ht_contra]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, contra_trace, t_v, plotting_color_contra)
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_contra(1), ylm_contra(2),...
    [], [], [], [], "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_1, ht_contra)   
plot(pnl_1, scale_bar, [ylm_contra(1), ylm_contra(1)], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')    
plot(pnl_1, [xlm_1(1), xlm_1(1)], [ylm_contra(1)+0.05, ylm_contra(1)+(0.05)+0.04], ...
    'LineWidth', get_default_scale_bar_width, 'Color', 'k')

pnl_1 = axes('Units', 'Centimeters', 'Position',[curr_x_pos2, curr_y_pos-norm_dist1, sz_1, ht_ipsi]);
hold on
plot(pnl_1, [t_v(1), t_v(end)], [0, 0], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 0.5)
get_default_SEM_area_plot(pnl_1, ipsi_trace, t_v, plotting_color_ipsi)
get_default_ax(pnl_1, xlm_1(1), xlm_1(2), [], [], ylm_ipsi(1), ylm_ipsi(2),...
    [], [], [], [], "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_1, ht_ipsi)    
plot(pnl_1, [xlm_1(1), xlm_1(1)], [0.1, 0.2], 'LineWidth', get_default_scale_bar_width, 'Color', 'k')

if annotation_cond == 1
    get_default_annotation(curr_x_pos2-0.2, curr_y_pos+ht_contra-0.75, ...
        {'0.1 ';'\DeltaF/F'}, 'k', 'normal', "left")    
    get_default_annotation(curr_x_pos1-0.2, curr_y_pos-norm_dist2+0.3, ...
        {'0.04 ';'\DeltaF/F'}, 'k', 'normal', "left") 
end
get_default_annotation(curr_x_pos2+sz_1, curr_y_pos-norm_dist2-0.1, ...
        '0.5 s', 'k', 'normal', "right")


end

function get_fig_panel_punish_run(behavior_train_data_table, img_train_data_table, color, ...
    x_pos, y_pos, exp_group, idx_plane, idx_hemisphere, ...
    idx_hemisphere_contra, dF_AP)

img_tbl = img_train_data_table(img_train_data_table.exp_group == char(exp_group),:);

curr_fly = img_train_data_table(idx_plane,:);
curr_fly = curr_fly.fly_id;
curr_behav_table = behavior_train_data_table(behavior_train_data_table.fly_id == curr_fly,:);
t_v = curr_behav_table.cut_mov.trace(:,1)/1000;

sz_1 = 6.5;
ht_1 = 1;
ht_2 = 1;
ht_3 = 1.5;
f_s = get_default_font_size;
y_dist = 0.1;
ylm_img_example = [-0.2, 1.6];
ylm_img_2 = [0, 0.9];
ylm_mov = [0, 2.3];
pun_times = [150, 245, 340];
line_width_traces = 0.5;

curr_data_table = behavior_train_data_table(behavior_train_data_table.exp_group == "AP",:);

for loop_idx = 1:size(curr_data_table,1)    
    movement_traces(:,loop_idx) = behavior_train_data_table.cut_mov(loop_idx,1).trace(1:12676,2);
end

movement_traces = abs(movement_traces);

%%

pun_color = [0, 0, 0];
LineThickness = 0.5;
pnl_2 = axes('Units','Centimeters','Position',[x_pos, y_pos, sz_1, ht_1+ht_2+(2*y_dist)+ht_3]);
hold on
plot(pnl_2, [pun_times(1), pun_times(1)], [0, 2.5], 'LineStyle', ':', 'Color', pun_color, 'LineWidth', LineThickness)
plot(pnl_2, [pun_times(2), pun_times(2)], [0, 2.5], 'LineStyle', ':', 'Color', pun_color, 'LineWidth', LineThickness)
plot(pnl_2, [pun_times(3), pun_times(3)], [0, 2.5], 'LineStyle', ':', 'Color', pun_color, 'LineWidth', LineThickness)
ylim([0, 2.5])
xlim([0, t_v(end)])
pnl_2.FontSize = f_s;
pnl_2.Color = 'none';
pnl_2.Box = 'off';
pnl_2.XAxis.Color = 'none';
pnl_2.YAxis.Color = 'none';
xlabel('time (s)')
pnl_2.TickLength = [0, 0];

%% example trace

curr_img_trace = (img_tbl.img_trace(idx_plane,1).trace(:,idx_hemisphere));
curr_img_panel = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos+ht_2+y_dist+0.1, sz_1, ht_1]);
get_default_img_trace_plot(curr_img_panel, curr_img_trace, 14.56, 0, 0.5, color.medium_gray,...
    line_width_traces, "with_y_scale_bar", ylm_img_example, [0, ylm_img_example(2)], 0.5, [], 0,...
    {'0.5'; '\DeltaF/F'}, 0.25)

%% movement

scale_bar_length = 30*(1/mean(diff(t_v)));
xlim_trc = [1, length(movement_traces)];
curr_pnl = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_2]);
get_default_SEM_area_plot(curr_pnl, movmean(movement_traces,22), 1:length(movement_traces), [0, 0, 0])
x_scale_bar = [xlim_trc(2)-scale_bar_length, xlim_trc(2)];
    plot(curr_pnl, x_scale_bar, [ylm_mov(1), ylm_mov(1)], 'LineWidth', 1.5, 'Color', 'k')
    get_default_annotation(x_pos+sz_1, y_pos-0.1, '30 s', 'k', 'normal', "right")
get_default_ax(curr_pnl, xlim_trc(1), xlim_trc(2), [], [],...
    ylm_mov(1), ylm_mov(2), [], [], [], [], "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_2)


y_label_offset = 0.25;
y_sb_length = 0.5;
curr_pnl3 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos, sz_1, ht_2]);
    plot(curr_pnl3, [1, 1], [0, y_sb_length], 'k', 'LineWidth', 1.5)
    get_default_ax(curr_pnl3, 1, 2, [], [], ylm_mov(1), ylm_mov(2),...
        [], [], [], [], "linear", "linear", [], [],...
        'none', 'none', 'none', 'none', sz_1, ht_1)
    get_default_annotation(x_pos+sz_1+0.2, y_pos+y_label_offset, {'0.5'; '|mm/s|'}, 'k', 'normal', "left")


%% mean imaging

xlim_trc = [1, length(dF_AP)];

curr_pnl = axes('Units','Centimeters','Position',[x_pos, y_pos+2*(y_dist)+ht_1+ht_2, sz_1, ht_3]);
    get_default_SEM_area_plot(curr_pnl, dF_AP, 1:length(dF_AP), [0, 0, 0])
get_default_ax(curr_pnl, xlim_trc(1), xlim_trc(2), xlim_trc(1), xlim_trc(2),...
    ylm_img_2(1), ylm_img_2(2), 0, ylm_img_2(2),...
    1, 0.25, "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_1, ht_2)

y_label_offset = 0.25;
y_sb_length = 0.3;
curr_pnl3 = axes('Units', 'Centimeters', 'Position', [x_pos+sz_1+0.1, y_pos+2*(y_dist)+ht_1+ht_2, sz_1, ht_3]);
plot(curr_pnl3, [1, 1], [0, y_sb_length], 'k', 'LineWidth', 1.5)    
get_default_ax(curr_pnl3, 1, 2, 1, 2, ylm_img_2(1), ylm_img_2(2),...
    ylm_img_2(1), ylm_img_2(2), 0, 0, "linear", "linear", [], [],...
    'none', 'none', 'none', 'none', sz_1, ht_1)
get_default_annotation(x_pos+sz_1+0.2, y_pos+2*(y_dist)+ht_1+ht_2+y_label_offset, ...
    {'0.3'; '\DeltaF/F'}, 'k', 'normal', "left")


%% Inset

sz_inset = 1.25;
ht_inset = 1;
ylm_img_1 = [-0.15, 1.25];
scale_bar_length = 0.5;
pun_onset = [pun_times(1)*14.56, pun_times(1)*14.56];
inset_idx = [2050, 2455];
curr_img_trace_ipsi = (img_tbl.img_trace(idx_plane,1).trace(:,idx_hemisphere));
curr_img_trace_contra = (img_tbl.img_trace(idx_plane,1).trace(:,idx_hemisphere_contra));
curr_inset_trace_ipsi = curr_img_trace_ipsi(:,1);
curr_inset_trace_contra = curr_img_trace_contra(:,1);
xlim_inset = inset_idx;
curr_img_panel = axes('Units', 'Centimeters', 'Position',[x_pos+0, y_pos+0.2, sz_inset, ht_inset]);
hold on
plot(curr_img_panel, pun_onset, ylm_img_1, 'r:')
plot(curr_img_panel, pun_onset+2*14.56, ylm_img_1, 'r:')
plot(curr_img_panel, curr_inset_trace_contra, 'Color', [0, 0, 0], 'LineWidth', 0.25)
plot(curr_img_panel, curr_inset_trace_ipsi, 'Color', color.medium_gray, 'LineWidth', 0.25)
get_default_SEM_area_plot(curr_img_panel, dF_AP, 1:length(dF_AP), color.yellow)
plot(curr_img_panel,[xlim_inset(2), xlim_inset(2)], [ylm_img_1(1), ylm_img_1(1)+scale_bar_length],...
    'k', 'LineWidth', get_default_scale_bar_width)
get_default_ax(curr_img_panel, xlim_inset(1), xlim_inset(2), xlim_inset(1), xlim_inset(2), ylm_img_1(1), ylm_img_1(2),...
    ylm_img_1(1), ylm_img_1(2), [], 2, "linear", "linear", [], [], 'none', 'none', 'none', 'none', sz_inset, ht_inset)
get_default_annotation(x_pos+sz_inset+0.2, y_pos+0.4, ...
    {'0.5'; '\DeltaF/F'}, 'k', 'normal', "left")
get_default_annotation(x_pos+0.2, y_pos+0.2+ht_inset, ...
    {'2 s'}, 'k', 'normal', "left")


end

