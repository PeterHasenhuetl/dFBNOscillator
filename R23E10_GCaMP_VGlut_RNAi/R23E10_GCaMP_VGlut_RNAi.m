% Analysis of R23E10-GAL4 driven GCaMP signals ±VGlutRNAi.
% Code written by Peter Hasenhuetl.


clear all
tic

% goes into folder with extracted traces and get fly names
cd([]) %Add path as character array
ind_fly = dir('*fly*');
n_flies = size(ind_fly,1);
mp_data_table_vGlut_RNAi = [];
curr_filter = get_mp_osc_filt(14.56, 0.005, 0.1, 0.025, 1.5);
% loops through individual flies
for fly_idx = 1:n_flies

        fly_input = dir('*fly*.mat');
        loaded_fly = load(char(fly_input(fly_idx).name));
        % loops through extracted planes
        ind_data_table = get_individual_multiplane_fly(loaded_fly, fly_idx, curr_filter);
        mp_data_table_vGlut_RNAi = [mp_data_table_vGlut_RNAi; ind_data_table];
        
        curr_name = loaded_fly.fly_details.flyname;
        fly_names(fly_idx,1:size(curr_name,2)) = curr_name;

end

cd([]) %Add path as character array
save('mp_data_table_vGlut_RNAi','mp_data_table_vGlut_RNAi')
toc


%% Custom functions

function ind_data_table = get_individual_multiplane_fly(input_fly, fly_idx, curr_filter)

fields_img_sessions = fieldnames(input_fly.img_traces);
    
for plane_idx = 1:length(fields_img_sessions)
    plane_name = fields_img_sessions{plane_idx};  
    output_img = get_mp_img_measures(input_fly.img_traces.(string(plane_name)), curr_filter);
    ind_data_table(plane_idx,:) = get_mp_sd_data_table_row(output_img, input_fly.fly_details, fly_idx, plane_idx);
end           

end

%%

function dm_row = get_mp_sd_data_table_row(dm_input, fly_details, fly_idx, plane_idx)

dm_behav_row = table;
dm_behav_row.exp_cond = string(fly_details.flyname(end));
dm_behav_row.fly_id = fly_idx;
dm_behav_row.plane_id = plane_idx;

dm_excl_row = table;
dm_excl_row.excl_idx = dm_input.excl_vec;


dm_img_row = table;

dF_trc = struct;
dF_trc.trace_img = dm_input.dF;
dm_img_row.img_trace = dF_trc;

trc_raw = struct;
trc_raw.trace_img = dm_input.raw_img;
dm_img_row.img_raw = trc_raw;

trc_infra = struct;
trc_infra.trace_img = dm_input.trace_infra;
dm_img_row.trace_infra = trc_infra;

trc_delta = struct;
trc_delta.trace_img = dm_input.trace_delta;
dm_img_row.trace_delta = trc_delta;

transient_vec = struct;
transient_vec.trace_img = dm_input.transient_vec;
dm_img_row.transient_vec = transient_vec;

a_c_delta = struct;
a_c_delta.trace_img = dm_input.auto_corr_delta';
dm_img_row.auto_corr_delta = a_c_delta;

pwr = struct;
pwr.trace_img = dm_input.amp;
dm_img_row.power_spectrum = pwr;

frqu_trc = struct;
frqu_trc.trace_img = dm_input.f;
dm_img_row.frqu_trc = frqu_trc;

trans_trace_raw = struct;
trans_trace_raw.trace_img = dm_input.trans_trace_raw;
dm_img_row.trans_trace_raw = trans_trace_raw;

trans_trace = struct;
trans_trace.trace_img = dm_input.trans_trace;
dm_img_row.trans_trace = trans_trace;

dm_img_row.int_delta_power = dm_input.int_delta_power;
dm_img_row.int_infra_slow_power = dm_input.int_infra_slow_power;
dm_img_row.ITI_CV = dm_input.ITI_CV;
dm_img_row.ITI_mean = dm_input.ITI_mean;
dm_img_row.transient_mean = dm_input.transient_mean;

dm_img_row.delta_period = dm_input.pr_delta;
dm_img_row.delta_corr = dm_input.amp_delta;
dm_img_row.background_noise = dm_input.background_noise;
dm_img_row.SNR = dm_input.SNR;

% Concatenates them
dm_row = [dm_behav_row, dm_excl_row, dm_img_row];

end
