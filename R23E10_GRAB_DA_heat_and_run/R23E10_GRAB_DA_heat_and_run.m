% Analysis of R23E10-GAL4-driven GRAB_DA signals during behavior in
% response to a heat shock.
% Code written by Peter Hasenhuetl.

clear all

% define the path for saving the source data for figures
source_data_details.data_path = []; %Add path as character array
source_data_details.file_name = 'R23E10_GRAB_on_ball_data.xlsx';


tic
% go into folder with extracted traces and get fly names
cd([]) %Add path as character array

ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);

behavior_train_data_table = [];
img_train_data_table = [];
behavior_freerun_data_table = [];
locomotion_mat = [];
LTA_mat = [];
randiLTA_mat = [];
pun_img_mat = [];
pun_locomotion_mat = [];


% loop through individual flies
for fly_idx = 1:n_flies
    
    cd([]) %Add path as character array
        
    fly_input = dir('*fly*.mat'); 
    loaded_fly = load(char(fly_input(fly_idx).name));
    
    cd([]) %Add path as character array
    
    [LTA, randiLTA, locomotion, pun_img, pun_locomotion, ...
        ind_data_table1, ind_data_table2, ind_data_table3, ind_data_table4,...
    dm_behav_row1, dm_behav_row2, dm_behav_row3, dm_behav_row4] = ...
    get_individual_GRAB_run_fly(loaded_fly, fly_idx);    
    
    behavior_train_data_table = [behavior_train_data_table; dm_behav_row1];
    behavior_freerun_data_table = [behavior_freerun_data_table; dm_behav_row2];
    LTA_mat = [LTA_mat, LTA];
    randiLTA_mat = [randiLTA_mat, randiLTA];
    locomotion_mat = [locomotion_mat, locomotion];
    
    pun_img_mat = [pun_img_mat, pun_img];
    pun_locomotion_mat = [pun_locomotion_mat, pun_locomotion];
      
end

cd([]) %Add path as character array
toc


%%

close all
get_fig_panel_GRAB_run_pun_resp(pun_img_mat, pun_locomotion_mat, 2, 2, source_data_details)

cd(source_data_details.data_path)
set(gcf,'renderer','Painters')
saveas(gcf,'R23E10_GRAB_on_ball_fig.pdf')
cd([]) %Add path as character array

%%

function [LTA, randiLTA, locomotion, pun_img, pun_locomotion, ...
    ind_data_table1, ind_data_table2, ind_data_table3, ind_data_table4,...
    dm_behav_row1, dm_behav_row2, dm_behav_row3, dm_behav_row4] = ...
    get_individual_GRAB_run_fly(input_fly, fly_idx)

input_fly.pwm_fr1 = [];
input_fly.pwm_fr2 = [];
input_fly.pwm_fr3 = [];

details = input_fly.details;
details.sampling_rate_2P = get_frame_rate_2P;
details.bl = 10;
details.n_trials_train = 3;
details.trial_length_train = 7000; %trial length for train in ms
details.baseline_jump = round(25*details.sampling_rate_2P); %seconds converted in #frames
details.pre_trial_window_img_train = round(5*details.sampling_rate_2P); %seconds converted in #frames    
details.post_TTL_lag = round(15.5*details.sampling_rate_2P);
details.idx_triallength = round((details.trial_length_train/1000)*details.sampling_rate_2P); 

%% imaging

TTL_input = input_fly.TTL_pun;
TTLoutputs_punish = get_TTL_output_punish(table2array(TTL_input));
curr_img = input_fly.imaging_pun;
ind_data_table1 = get_punish_session(curr_img, TTL_input, TTLoutputs_punish, details, fly_idx);
input_pwm = table2array(input_fly.pwm_pun(:,[1, 8]));
input_mov = table2array(input_fly.mov_pun);
[dm_behav_row1, output_behav1] = get_GRAB_behavior_analysis(input_mov, input_pwm, details, fly_idx);

if isempty(input_fly.pwm_fr1) == 0
    TTL_input = input_fly.TTL_fr1;
    TTLoutputs_fr1 = get_TTL_output_punish(table2array(TTL_input));
    curr_img = input_fly.imaging_fr1;
    ind_data_table2 = get_run_session(curr_img, TTL_input, TTLoutputs_fr1, details, fly_idx);
    input_pwm = table2array(input_fly.pwm_fr1(:,[1, 8]));
    input_mov = table2array(input_fly.mov_fr1);
    [dm_behav_row2, output_behav2] = get_GRAB_behavior_analysis(input_mov, input_pwm, details, fly_idx);
    
else
    ind_data_table2 = [];
    dm_behav_row2 = []; 
    output_behav2 = [];
end


if isempty(input_fly.pwm_fr2) == 0
    TTL_input = input_fly.TTL_fr2;
    TTLoutputs_fr2 = get_TTL_output_punish(table2array(TTL_input));
    curr_img = input_fly.imaging_fr2;
    ind_data_table3 = get_run_session(curr_img, TTL_input, TTLoutputs_fr2, details, fly_idx);
    input_pwm = table2array(input_fly.pwm_fr2(:,[1, 8]));
    input_mov = table2array(input_fly.mov_fr2);
    [dm_behav_row3, output_behav3] = get_GRAB_behavior_analysis(input_mov, input_pwm, details, fly_idx);
else
    ind_data_table3 = [];
    dm_behav_row3 = []; 
    output_behav3 = [];
end



if isempty(input_fly.pwm_fr3) == 0
    TTL_input = input_fly.TTL_fr3;
    TTLoutputs_fr3 = get_TTL_output_punish(table2array(TTL_input));
    curr_img = input_fly.imaging_fr3;
    ind_data_table4 = get_run_session(curr_img, TTL_input, TTLoutputs_fr3, details, fly_idx);
    input_pwm = table2array(input_fly.pwm_fr3(:,[1, 8]));
    input_mov = table2array(input_fly.mov_fr3);
    [dm_behav_row4, output_behav4] = get_GRAB_behavior_analysis(input_mov, input_pwm, details, fly_idx);
else
    ind_data_table4 = [];
    dm_behav_row4 = []; 
    output_behav4 = [];
end

[img_trace_pun, img_trace_fr1, img_trace_fr2, img_trace_fr3] = ...
    get_selected_ROI(ind_data_table1, ind_data_table2, ind_data_table3, ind_data_table4);


wind_size_1 = 35;
dm_behav_row1.dm = get_datamatrix(img_trace_pun, TTLoutputs_punish,output_behav1);

pun_idx1 = find(dm_behav_row1.dm(1,1).data_matrix(:,2) == 20,1,'first');
pun_idx1 = pun_idx1+round(2/(mean(diff(dm_behav_row1.dm(1,1).data_matrix(:,1)))/1000));
cons_wind = round(wind_size_1/(mean(diff(dm_behav_row1.dm(1,1).data_matrix(:,1)))/1000));
pun_locomotion = dm_behav_row1.dm(1,1).data_matrix(pun_idx1-cons_wind:pun_idx1+cons_wind,3);
pun_img = dm_behav_row1.dm(1,1).data_matrix(pun_idx1-cons_wind:pun_idx1+cons_wind,6);
pun_img = (pun_img-mean(pun_img(1:cons_wind,1)))/mean(pun_img(1:cons_wind,1));


if isempty(input_fly.pwm_fr1) == 0
    
    
    dm_behav_row2.dm = get_datamatrix(img_trace_fr1, TTLoutputs_fr1,output_behav2);
    dm_behav_row2.movement_initiation = analyze_movement(dm_behav_row2.dm.data_matrix);
    LTA = mean(dm_behav_row2.movement_initiation.LTA,2,'omitnan');
    randiLTA = mean(dm_behav_row2.movement_initiation.randiLTA,2,'omitnan');
    locomotion = mean(dm_behav_row2.movement_initiation.locomotion,2,'omitnan');
    
else
    
    LTA = [];
    randiLTA = [];
    locomotion = [];
    
end

if isempty(input_fly.pwm_fr2) == 0

    dm_behav_row3.dm = get_datamatrix(img_trace_fr2, TTLoutputs_fr2,output_behav3);
    dm_behav_row3.movement_initiation = analyze_movement(dm_behav_row3.dm.data_matrix);    

end

if isempty(input_fly.pwm_fr3) == 0

    dm_behav_row4.dm = get_datamatrix(img_trace_fr3, TTLoutputs_fr3,output_behav4);
    dm_behav_row4.movement_initiation = analyze_movement(dm_behav_row4.dm.data_matrix);

end


end

function [img_trace_pun, img_trace_fr1, img_trace_fr2, img_trace_fr3] = ...
    get_selected_ROI(ind_data_table1, ind_data_table2, ind_data_table3, ind_data_table4)


SNR_array = ind_data_table1.SNR_ROI;
[SNR_row,  SNR_column] = find(SNR_array == max(max(SNR_array)));

img_trace_pun = [ind_data_table1.st_img(SNR_row,1).trace(:,SNR_column), ...
    ind_data_table1.st_TTL(1,1).trace(:,1)];
if isempty(ind_data_table2) == 0
    img_trace_fr1 = [ind_data_table2.st_img(SNR_row,1).trace(:,SNR_column), ...
        ind_data_table2.st_TTL(1,1).trace(:,1)];
else
    img_trace_fr1 = [];
end

if isempty(ind_data_table3) == 0
    img_trace_fr2 = [ind_data_table3.st_img(SNR_row,1).trace(:,SNR_column), ...
        ind_data_table3.st_TTL(1,1).trace(:,1)];
else
    img_trace_fr2 = [];
end

if isempty(ind_data_table4) == 0
    img_trace_fr3 = [ind_data_table4.st_img(SNR_row,1).trace(:,SNR_column), ...
        ind_data_table4.st_TTL(1,1).trace(:,1)];
else
    img_trace_fr3 = [];
end


end

function ind_data_table = get_punish_session(curr_img, TTL_input, TTLoutputs, details, fly_idx)

fields_img_planes = fieldnames(curr_img);

for plane_idx = 1:length(fields_img_planes)
    plane_name = fields_img_planes{plane_idx};
    [output_img, ~] = get_GRAB_img_measures(curr_img.(string(plane_name)), TTL_input, TTLoutputs);
    ind_data_table(plane_idx,:) = get_GRAB_data_table_row(output_img, details, fly_idx, plane_idx);
end

end

function ind_data_table = get_run_session(curr_img, TTL_input, TTLoutputs, details, fly_idx)

fields_img_planes = fieldnames(curr_img);

for plane_idx = 1:length(fields_img_planes)
    plane_name = fields_img_planes{plane_idx};
    [output_img, ~] = get_GRAB_img_measures(curr_img.(string(plane_name)), TTL_input, TTLoutputs);
    ind_data_table(plane_idx,:) = get_GRAB_data_table_row(output_img, details, fly_idx, plane_idx);
end


end

function [img_output, TTLoutputs_nu] = get_GRAB_img_measures(imaging, TTL_input, TTLoutputs)

f_r = get_frame_rate_2P; % frame-rate   
raw_green = table2array(imaging.img);
img_green = raw_green(:,2:end) - mean(raw_green(:,1));

%%

pre_punish_window = round(10*f_r);
post_punish_window = round(10*f_r);
cut_off_start = TTLoutputs(1)-pre_punish_window;
cut_off_end = TTLoutputs(1)+post_punish_window;
TTLoutputs_nu = TTLoutputs-cut_off_start;
TTL_input_nu = TTL_input(cut_off_start:cut_off_end,:);

img_output = struct;
for R_idx = 1:size(img_green,2)
    %% First, tests for ROI to be included
    
    if raw_green(:,R_idx+1) == raw_green(:,1)
        img_output.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);
        
    elseif mean(img_green(:,R_idx),1) < 40
        img_output.excl_vec(1,R_idx) = 1;
        img_green(:,R_idx) = NaN(length(img_green(:,R_idx)),1);
        
    else
        img_output.excl_vec(1,R_idx) = 0;
    end
    
    %% Then, computes dF
    
    r = img_green(cut_off_start:cut_off_end,R_idx);
    dF = (r - mean(r(1:(pre_punish_window),1),1))./(mean(r(1:(pre_punish_window),1),1));
    img_output.dF(:,R_idx) = dF;
    img_output.raw_img(:,R_idx) = img_green(:,R_idx);

    %% ...punishment-aligned img

    st_img(:,R_idx) = dF;
    st_TTL(:,R_idx) = TTL_input_nu(:,1);
    
    cons_window = 4;
    img_output.SNR_ROI(1,R_idx) = ...
        max(movmean(st_img((pre_punish_window):(pre_punish_window)+round(f_r*cons_window),R_idx),10));
    img_output.dF_amplitude(1,R_idx) = ...
        mean(movmean(st_img((pre_punish_window):(pre_punish_window)+round(f_r*cons_window),R_idx),3));
            
end
   
%% stores data in structure
singletrialtime = 1;
img_output.singletrialtime = singletrialtime;
img_output .img_green = img_green;
img_output.st_img = img_green;
img_output.st_TTL = table2array(TTL_input);
img_output.TTLoutputs_nu = TTLoutputs_nu;

end

function TTLoutputs = get_TTL_output_punish(TTL_input)

nTTL = TTL_input+5000;
nTTL = nTTL-min(nTTL);    
nTTL = nTTL./max(nTTL);

    
% Generates vector with 'ones' at indices of TTL onset;    
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

TTLoutputs(1,1) = find(TTLidx == 1,1,'first');


end

function [dm_behav_row, output_behav] = get_GRAB_behavior_analysis(mov_test1, pwm_test1, ...
       details, fly_idx)

mov_test = mov_test1;
pwm_test = pwm_test1;

%sets both time vectors to start from zero
mov_test(:,1) = mov_test(:,1) - mov_test(1,1);
pwm_test(:,1) = pwm_test(:,1) - pwm_test(1,1);

% divides movement trace by 5 (because Matlab-balltracking multiplies is by
% 5 when communicating with LabView)
mov_test(:,2:end) = mov_test(:,2:end)./5;
output_behav = struct;

%% calculates ratio between time vectors (mov and pwm) and adjusts them to be the same:

% the factor by which the camera output is higher than the pwm output, should be ~5
s_test = mov_test(end,1)/pwm_test(end,1);

%bin-size of camera in iterations of output from Matlab to Labview
frt_1 = mov_test(end,1)/(length(mov_test)-1); %bin-size of camera in iterations of output from Matlab to Labview
bin_size_mov = frt_1/s_test; % actual bin-size in milliseconds
mov_test(:,1) = mov_test(:,1)/s_test; % corrects the time-vector of movement to correct time

%% finds starting points of individual trials in pwm-trace and movement-trace

[~, trialonindex_pwm, trialontime_pwm, ~, ...
    trialon_index_position, trialontime_position, trialoff_idx_position, ...
    trialofftime_position] = get_punish_onset_GRAB(pwm_test, mov_test, details, bin_size_mov, details.trial_length_train);

%%

length_st_mov_trace = 668; %Number of indices for length of single-trial movement traces    
window_1 = round(5*(1/(bin_size_mov/1000))); %seconds converted in #frames
window_2 = round(12*(1/(bin_size_mov/1000)));
running.running(:,1) = (mov_test(trialon_index_position(1)-115:trialon_index_position(1)+length_st_mov_trace,2));
running.ptwm = window_1;    
running.ptlm = window_2;
running.pwm(:,1) = (pwm_test(trialonindex_pwm(1)-20:trialonindex_pwm(1)+70,2));

% magnitude of actual movement curve       
running.summed_running(1,1) = (sum(mov_test(trialon_index_position(1):trialoff_idx_position(1)+window_2,2)));

% Area under absolute movement curve    
running.summed_running(1,2) = (sum(abs(mov_test(trialon_index_position(1):trialoff_idx_position(1)+window_2,2))));
 
%%

output_behav.running = running;
output_behav.mov = mov_test;
output_behav.PWM = pwm_test;
output_behav.trialoff_indices_position = trialoff_idx_position;
output_behav.trialon_indices_position = trialon_index_position;
output_behav.trialontimes_position = trialontime_position;
output_behav.trialofftimes_position = trialofftime_position;
output_behav.trialontimes_pwm = trialontime_pwm;

%%

pre_punish_window = round(30/(bin_size_mov/1000));
post_punish_window = round(30/(bin_size_mov/1000));
cut_off_start = trialon_index_position-pre_punish_window;
cut_off_end = trialon_index_position+post_punish_window;

nu_mov = mov_test(cut_off_start:cut_off_end,:);
nu_mov(:,1) = nu_mov(:,1) - nu_mov(1,1);

%%
dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.trialon_index_position = trialon_index_position;

mov = struct;
mov.trace = output_behav.mov;
dm_behav_row.mov = mov;

cut_mov = struct;
cut_mov.trace = nu_mov;
dm_behav_row.cut_mov = cut_mov;

pwm = struct;
pwm.trace = output_behav.PWM;
dm_behav_row.pwm = pwm;
dm_behav_row.bin_size_mov = bin_size_mov;

end

function [trial_starts_test, trialonindex_pwm, trialontime_pwm, trialofftime_pwm, ...
    trialon_index_position, trialontime_position, trialoff_idx_position, ...
    trialofftime_position] = ...
    get_punish_onset_GRAB(pwm_test, mov_test, details, ...
    bin_size_mov, trial_length)


trial_starts_test = zeros(length(pwm_test),1); 
pre_trial_indices = 3;
dts_test = [0; diff(pwm_test(:,2))];

% Creates a vector that sepcifies where the pwm-signal jumped from baseline
% to the value in that trial.
tr_test = (pwm_test(:,2)) == (dts_test+details.bl);

for pwm_idx = 1:length(pwm_test)
    
    % ignore time bins before the pre-specified 'test_onset_lag'
    if pwm_idx < 5
        
        trial_starts_test(pwm_idx) = 0;
        
    else
        
         tt_test(pwm_idx) = any(dts_test(pwm_idx-pre_trial_indices:pwm_idx-2,1) ~= 0); 
         trial_starts_test(pwm_idx-1) = (tr_test(pwm_idx) == 0 && tr_test(pwm_idx-1) == 1 && tt_test(pwm_idx) == 0);
         
    end
    
end


% Finds the trial onset times according to the trial onset indices.
trialonindex_pwm = find(trial_starts_test == 1,1,'first');
trialontime_pwm = pwm_test(trialonindex_pwm,1);

%% Gets trial on- and off-times in movement-trace

% First, gets the indices in the movement trace that correspond to the 
% on-times in the pwm-trace
trialon_index_position = find(abs(mov_test(:,1) - trialontime_pwm) < 50,1,'first');

% Gets a vector of the trial-onset times (in time-vector of movement trace) that
% correspond to trial-on-indices
trialontime_position = mov_test(trialon_index_position,1);
trialoff_idx_position = find(abs(mov_test(:,1) - (trialontime_position+trial_length)) < bin_size_mov,1,'first');
trialoffidx_pwm = find(abs(pwm_test(:,1) - (trialontime_pwm+(trial_length))) < mean(diff(pwm_test(:,1))),1,'first');    
trialofftime_pwm = pwm_test(trialoffidx_pwm,1);
trialofftime_position = mov_test(trialoff_idx_position,1);

end

function dm_row = get_GRAB_data_table_row(dm_input, details, fly_idx, plane_idx)

dm_behav_row = table;
dm_behav_row.fly_id = fly_idx;
dm_behav_row.plane_id = plane_idx;
dm_behav_row.SNR_ROI = dm_input.SNR_ROI;
dm_behav_row.flyname = details.flyname;

%% Adds the imaging metrics

dm_img_row = table;

dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;
dm_excl_row.dF_amplitude = dm_input.dF_amplitude;

trc = struct;
trc.trace = dm_input.dF;
dm_img_row.img_trace = trc;

trc_raw = struct;
trc_raw.trace = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;

st_img = struct;
st_img.trace = dm_input.st_img;
dm_img_row.st_img = st_img;

st_TTL = struct;
st_TTL.trace = dm_input.st_TTL;
dm_img_row.st_TTL = st_TTL;

%% Concatenates them

dm_row = [dm_behav_row, dm_excl_row, dm_img_row];

end

function dm = get_datamatrix(imaging_trace, TTLoutputs_punish, output_test)

cutoff_idx_data_matrix = 1500;
bin_size_mov = mean(diff(output_test.mov(:,1)));
sampling_width_2P = 1/get_frame_rate_2P;

%% First, interpolates imaging data
% Interpolates datapoints to adjust to movement frame-rate

img_time = [0:sampling_width_2P:((length(imaging_trace)-1))/(1/sampling_width_2P)]';
xq_img = ([0:bin_size_mov:1000*(img_time(end,1))]');
xq_img = xq_img/1000;
int_img = interp1(img_time,imaging_trace,xq_img);

%% -----Then interpolates movement/PWM data

xq_test = [0:bin_size_mov:output_test.mov(end,1)]'; %TIME VECTOR FOR INTERPOLATION

% Interpolates
int_pwm = pchip(output_test.PWM(:,1),output_test.PWM(:,2),xq_test);

%% specifies align_mov_trace
align_mov_trace = zeros(length(int_pwm),3);

align_mov_trace(:,1) = xq_test;

% adds pwm-trace
align_mov_trace(:,2) = int_pwm;

% adds all movement traces
align_mov_trace(:,3) = (output_test.mov(:,2));

%% specifies cutoff index for int_img

indx_xq = find((abs(xq_img-img_time(TTLoutputs_punish(1)))) == min(abs(xq_img-img_time(TTLoutputs_punish(1)))),1,'first');
cutoff_img = indx_xq - cutoff_idx_data_matrix; % is the index from which to cut off imaging
indx_mov = find(abs(align_mov_trace(:,1) - output_test.trialontimes_position(1)) == ...
    min(abs(align_mov_trace(:,1) - output_test.trialontimes_position(1))),1,'first');
cutoff_mov = indx_mov - cutoff_idx_data_matrix; % is the index from which to cut off trace 

%% Then, puts them into data matrix

align_img_trace = int_img(cutoff_img:end,:);
align_mov_trace = align_mov_trace(cutoff_mov:end,:);
size_traces1 = size(align_img_trace,2);
size_traces2 = size(align_mov_trace,2);

aligned_traces = NaN(350000,size_traces1+size_traces2);

aligned_traces(1:length(align_mov_trace),1:size(align_mov_trace,2)) = align_mov_trace;

sz_al1 = size(aligned_traces,2);

% then, puts in imaging and TTL (from imaging)
aligned_traces(1:length(align_img_trace),sz_al1+1) = (align_img_trace(:,1));
aligned_traces(1:length(align_img_trace),sz_al1+2) = align_img_trace(:,2);

% cuts data matrix according to which of the traces is shorter
if length(align_mov_trace) > length(align_img_trace)
    aligned_traces = aligned_traces(1:length(align_img_trace),:);
elseif length(align_mov_trace) < length(align_img_trace)
    aligned_traces = aligned_traces(1:length(align_mov_trace),:); 
end

dm = struct;
dm.data_matrix = aligned_traces;

end

function movement_initiation = analyze_movement(curr_matrix)


input_mov = curr_matrix(:,3);
img_trace = curr_matrix(:,6);
run_length = 3; % in seconds
term_length = 20; % in seconds
rest_length = 3; % in seconds
b_t = 0.5; % binarize threshold
run_fac = 80;

b_mov = zeros(size(input_mov,1),size(input_mov,2)); % inititalizing vector for binarized movement;
b_w = mean(diff(curr_matrix(:,1)))/1000;
rest_criterion = round(rest_length/b_w);
run_criterion = round(run_length/b_w); % length of window (as indices) to define running
run_def = round(run_criterion*run_fac); % threshold for proper running
term_criterion = round(term_length/b_w);
cons_wind1 = round(3.5/b_w);

for loop_idx1 = 1:size(input_mov,1)
    %if movement trace is larger than threshold, then it is represented by a '1'   
    if input_mov(loop_idx1,1) > b_t
        b_mov(loop_idx1,1) = 1;   
    end
end

%% finds indices of running-initiation       
% finds the indices where fly was active (i.e. the indices of the '1' in the binarized movement matrix
% and adds a '1' as first index for analyis below

act_idx = [1; find(b_mov == 1)];

%% computes the number of indices between activity bins     

act_window = [1; diff(act_idx)];
term_window = [diff(act_idx); 1];
comb_act = [act_idx, act_window, term_window];

%% finds the indices where movement was initiated and terminated
% ... and then below exludes the ones that don't meet running criteria.

init_idx = comb_act(find(comb_act(:,2) > rest_criterion),1);
term_idx = comb_act(find(comb_act(:,3) > term_criterion),1);

%% Finds criterion for movement vigor


for loop_idx2 = 1:length(init_idx)
    
    if round(init_idx(loop_idx2)+run_criterion) > length(b_mov(:,1))
        
        init_idx(loop_idx2) = NaN;
    
    elseif sum(b_mov(init_idx(loop_idx2):round(init_idx(loop_idx2)+run_criterion),1)) < run_def && ...
            any(ismember([init_idx(loop_idx2):round(init_idx(loop_idx2)+run_criterion)],term_idx)) == 1
    
        init_idx(loop_idx2) = NaN;
                
    elseif any(ismember([init_idx(loop_idx2):round(init_idx(loop_idx2)+run_criterion)],term_idx)) == 1
    
        init_idx(loop_idx2) = NaN;
    
    end
    
end

scale_f = 4;
init_idx(find(isnan(init_idx))) = [];
init_idx(init_idx<(rest_criterion*scale_f)) = [];
randi_idx = randi([min(init_idx),max(init_idx)],[length(init_idx)*3,1]);
randi_idx = sort(randi_idx,'ascend');

if length(init_idx) < 5
    
    movement_initiation.im_corr = [];
    movement_initiation.lags = [];
    movement_initiation.init_idx = [];
    movement_initiation.LTA = [];
    movement_initiation.randiLTA = [];
    movement_initiation.locomotion = [];

else    
    
    %% Computes the locomotion-triggered average
    rest_criterion = rest_criterion*scale_f;
    run_criterion = run_criterion*scale_f;

    time_vector = linspace(-rest_length,run_length,rest_criterion+run_criterion)/1000; % time-vector in seconds

    select_idx1 = rest_criterion-0;
    select_idx2 = rest_criterion+0;

    for loop_idx3 = 1:length(init_idx)-1
        pre_LTA(:,loop_idx3) = img_trace(init_idx(loop_idx3)-rest_criterion:init_idx(loop_idx3)+run_criterion,1);
        max_dF = max(pre_LTA(select_idx1:select_idx2,loop_idx3));
        realignment_idx(loop_idx3) = rest_criterion;
        F0_window = 20;
        LTA(:,loop_idx3) = pre_LTA(realignment_idx(loop_idx3)-cons_wind1:realignment_idx(loop_idx3)+cons_wind1,loop_idx3);
        LTA(:,loop_idx3) = (LTA(:,loop_idx3)-mean(LTA(1:F0_window,loop_idx3),1))./mean(LTA(1:F0_window,loop_idx3),1);
        pre_locomotion(:,loop_idx3) = input_mov(init_idx(loop_idx3)-rest_criterion:init_idx(loop_idx3)+run_criterion,1);
        locomotion(:,loop_idx3) = ...
            pre_locomotion(realignment_idx(loop_idx3)-cons_wind1:realignment_idx(loop_idx3)+cons_wind1,loop_idx3);
    end
    
    
    for loop_idx4 = 1:length(randi_idx)-1
        
        pre_randiLTA(:,loop_idx4) = img_trace(randi_idx(loop_idx4)-rest_criterion:randi_idx(loop_idx4)+run_criterion,1);
        realignment_idx(loop_idx4) = rest_criterion;
        F0_window = 20;
        randiLTA(:,loop_idx4) = pre_randiLTA(realignment_idx(loop_idx4)-cons_wind1:realignment_idx(loop_idx4)+cons_wind1,loop_idx4);
        randiLTA(:,loop_idx4) = (randiLTA(:,loop_idx4)-mean(randiLTA(1:F0_window,loop_idx4),1))./mean(randiLTA(1:F0_window,loop_idx4),1);
    end
    

    LTA_mean = mean(LTA,2,'omitnan');
    std_LTA = std(LTA','omitnan');
    std_LTA = std_LTA';
    locomotion_mean = mean(locomotion,2,'omitnan');
    init_number = 1:size(locomotion(1,1:end),2);

    %% Computes cross-correlation

    nu_LTA = LTA(:,2:end-1);
    nu_locomotion = locomotion(:,2:end-1);
    max_lag = 200;
    lags = NaN(size(nu_LTA,2),(2*max_lag)+1);
    im_corr = NaN((2*max_lag)+1,size(nu_LTA,2));

    for loop_idx5 = 1:size(nu_LTA,2)
        [im_corr(:,loop_idx5), lags(loop_idx5,:)] = xcorr(nu_LTA(:,loop_idx5), nu_locomotion(:,loop_idx5),max_lag);
    end

    %% Stores the results in structure

    movement_initiation.im_corr = im_corr;
    movement_initiation.lags = lags;
    movement_initiation.init_idx = init_idx;
    movement_initiation.LTA = LTA(:,1:end-1);
    movement_initiation.randiLTA = randiLTA(:,1:end-1);
    movement_initiation.locomotion = locomotion(:,1:end-1);

end


end

function frame_rate_2P = get_frame_rate_2P

frame_rate_2P = 14.56;

end
