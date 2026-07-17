function output_img = get_mp_img_measures(imaging_data, curr_filter)
% Analysis of dFBN GCaMP dynamics.
% Code written by Peter Hasenhuetl.

cumulative_idx1 = 1;
cumulative_idx2 = 1;
f_r = 14.56; % frame-rate
sliding_window = 501;
threshold_brightness = 14;    
   
cond_dR = 0;
raw_green = table2array(imaging_data.green);
raw_red = table2array(imaging_data.red);
img_green = raw_green(:,2:5) - mean(raw_green(:,1));

img_red = raw_red(:,2:5) - mean(raw_red(:,1));
img_green = img_green(501:end,:);
img_red = img_red(501:end,:);
transient_onsets = NaN(10000,4);

if size(raw_green,2) == 6
    output_img.background_noise = prctile(raw_green(:,end)-mean(raw_green(:,end)),75);
else
    output_img.background_noise = NaN;
end
% loops through individual ROIs
for ROI_idx = 1:4

    % first, tests for ROI to be included
    if raw_red(:,ROI_idx+1) == raw_red(:,1)
        output_img.excl_vec(1,ROI_idx) = 1;
        img_green(:,ROI_idx) = NaN(length(img_green(:,ROI_idx)),1);
        img_red(:,ROI_idx) = NaN(length(img_red(:,ROI_idx)),1);
    elseif mean(img_green(:,ROI_idx),1) < threshold_brightness
        output_img.excl_vec(1,ROI_idx) = 1;
        img_green(:,ROI_idx) = NaN(length(img_green(:,ROI_idx)),1);
        img_red(:,ROI_idx) = NaN(length(img_red(:,ROI_idx)),1);
    else
        output_img.excl_vec(1,ROI_idx) = 0;
    end
    
    % then, computes dF/F
    raw_F = img_green(:,ROI_idx);
    
    if output_img.excl_vec(1,ROI_idx) == 0 
        
        if cond_dR == 1
            F_red = (img_red(:,ROI_idx));
            raw_F = (raw_F)./(F_red);
        end
        F_0 = get_mov_prctile(raw_F, sliding_window);
        d_F = (raw_F - F_0)./(F_0);
        d_F = smoothdata(d_F,'Gaussian',8);
    else
        d_F = zeros(length(raw_F),1);
    end

    output_img.dF(:,ROI_idx) = d_F;
    output_img.raw_img(:,ROI_idx) = img_green(:,ROI_idx);
    output_img.SNR(ROI_idx) = mean(img_green(:,ROI_idx))/output_img.background_noise;
    output_img.fraction_background_noise(ROI_idx) = output_img.background_noise/mean(img_green(:,ROI_idx));

    if isnan(mean(d_F)) == 1
        d_F = zeros(size(d_F,1),size(d_F,2));
    end
    
    [one_sided_pwr_sctrm, f] = pspectrum(d_F,14.56);
    output_img.amp(:,ROI_idx) = one_sided_pwr_sctrm;
    d_idx1 = find(f > 0.2,1,'first');
    d_idx2 = find(f > 1,1,'first');
    output_img.int_delta_power(ROI_idx) = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
    d_idx1 = find(f > 0.005,1,'first');
    d_idx2 = find(f > 0.1,1,'first');
    output_img.int_infra_slow_power(ROI_idx) = sum(one_sided_pwr_sctrm(d_idx1:(d_idx2-1)),'omitnan');
    output_img.f = f;

    % ANALYZES Ca TRANSIENTS
    [output_img.transient_vec(:,ROI_idx), transient_amplitudes, trans_ids] = ...
        get_GCaMP_transients(d_F);

    % some measures based on inter-transient intervals (ITIs)
    curr_ITIs = (diff(find(output_img.transient_vec(:,ROI_idx) > 0))./f_r).*1000;
    output_img.ITI_std(ROI_idx) = std(curr_ITIs);
    output_img.ITI_mean(ROI_idx) = mean(curr_ITIs);
    output_img.ITI_CV(ROI_idx) = ...
        output_img.ITI_std(ROI_idx)./abs(output_img.ITI_mean(ROI_idx));
    output_img.transient_mean(ROI_idx) = mean(transient_amplitudes);

    [output_img.pr_delta(ROI_idx), output_img.amp_delta(ROI_idx),...
        output_img.auto_corr_delta(ROI_idx,:), output_img.trace_infra(:,ROI_idx), ...
        output_img.trace_delta(:,ROI_idx)] = ...
        get_auto_corr(d_F,curr_filter);

    % extracts transients
    p_wind = 40;
    n_wind = 15;
    if ROI_idx > 2
        for u_idx = 1:2 % "u_idx" defines two categories of transients
            if u_idx == 1 % First, raw traces (not of the full dF/F curve)
                curr_img_trace = img_green(:,ROI_idx);
            elseif u_idx == 2 % Second, traces from the full dF/F curve
                curr_img_trace = d_F;
            end

            curr_transient_onsets = trans_ids;
            transient_onsets(1:length(curr_transient_onsets),ROI_idx) = curr_transient_onsets;
            % Only considers the transient onsets that allow for the full
            % length of transients to be used (i.e., excludes transients at
            % the beginning and end of the imaging trace).
            curr_transient_onsets((curr_transient_onsets+p_wind) > length(curr_img_trace) | ...
                (curr_transient_onsets-n_wind) < 1) =  [];
            % Then, extracts the transient traces.
            if isempty(curr_transient_onsets) == 1 
                % If ther are no transients in trace, creates NaNs
                dummy_trace = (curr_img_trace(20-n_wind:20+p_wind,1));
                output_img.trans_trace = NaN(size(dummy_trace,1),1);
                output_img.trans_trace_raw = NaN(size(dummy_trace,1),1);
            end

            if u_idx == 1 % Raw transients (whch can then be dF/F normalized)
                for loop_idx_trans_1 = 1:length(curr_transient_onsets)
                    output_img.trans_trace_raw(:,cumulative_idx1) = ...
                        (curr_img_trace(curr_transient_onsets(loop_idx_trans_1)-...
                        n_wind:curr_transient_onsets(loop_idx_trans_1)+p_wind,1));
                    cumulative_idx1 = cumulative_idx1+1;
                end
            elseif u_idx == 2 % dF/F transients (from full dF/F curve)
                for loop_idx_trans_2 = 1:length(curr_transient_onsets)
                    output_img.trans_trace(:,cumulative_idx2) = ...
                        (curr_img_trace(curr_transient_onsets(loop_idx_trans_2)-...
                        n_wind:curr_transient_onsets(loop_idx_trans_2)+p_wind,1));
                    cumulative_idx2 = cumulative_idx2+1;
                end
            end

        end
    end

end


output_img.transient_onsets = transient_onsets;

end

%%

function [transient_vec, transient_amplitudes, trans_ids] = get_GCaMP_transients(img_trace)

gaussian_sliding_window = 15;
min_peak_height = 0.0001;
trans_signal = (1)*smoothdata(diff(img_trace),'gaussian',gaussian_sliding_window);
[~, trans_ids] = findpeaks(trans_signal,'MinPeakHeight',min_peak_height);

padded_img = [img_trace; zeros(20,1)];
transient_amplitudes = NaN(1,length(trans_ids));
for loop_idx = 1:length(trans_ids)
    transient_amplitudes(loop_idx) = max(padded_img(trans_ids(loop_idx):trans_ids(loop_idx)+10));   
end

transient_vec = zeros(length(img_trace),1);
transient_vec(trans_ids,1) = 1;

end

%%

function [pr_delta, amp_delta, auto_corr_delta, ...
    trace_infra, trace_delta] = get_auto_corr(img_signal, curr_filter)

length_corr = 145; % ~10 seconds

if sum(img_signal) == 0 || isnan(sum(img_signal))
    [auto_corr_delta, ~] = xcorr((img_signal),length_corr,'coef');
    pr_delta = NaN;
    amp_delta = NaN;
    trace_infra = NaN(size(img_signal));
    trace_delta = NaN(size(img_signal));
else
    filt_img = get_filtered_mp_img(curr_filter,img_signal);
    curr_signal = [0; diff(img_signal)];
    [auto_corr_delta, ~] = xcov(curr_signal,length_corr,'coef');
    trace_infra = filt_img.fimg_1';
    trace_delta = filt_img.fimg_2';
    PRD = get_period(auto_corr_delta, 15, 0, 0.05, length_corr);
    pr_delta = PRD.curr_pr;
    amp_delta = PRD.curr_amp;
end

end

%%

function PRD = get_period(img_signal, min_dist, min_height, min_prom, length_corr)

[peak_amp, lag_dist] = findpeaks(img_signal,'MinPeakDistance',min_dist,'MinPeakProminence',min_prom,...
    'MinPeakHeight',min_height);
lag_dist = lag_dist-(length_corr+1);
lag_dist = abs(lag_dist);

if size(lag_dist,1) == 1
    PRD.curr_amp = NaN;
    PRD.curr_pr = NaN;
else
    [nu_l, curr_idx] = sort(lag_dist,'ascend');
    PRD.curr_amp = peak_amp(curr_idx(2));
    PRD.curr_pr = nu_l(2)/14.56;
end

end
