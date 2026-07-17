function [pr_delta, amp_delta, cross_corr_delta, ...
    corrcoef_delta, corrcoef_infra, auto_corr_delta_dendrite_left, auto_corr_delta_dendrite_right, time_lag] = ...
    get_mp_cross_corr(img_signal1, img_signal2, curr_filter)
% Interhemispheric crosscorrelations of dFBN GCaMP traces.
% Code written by Peter Hasenhuetl.

length_corr = 145; % ~10 seconds

filt_img_dendrite1 = get_filtered_mp_img(curr_filter, img_signal1);
filt_img_dendrite2 = get_filtered_mp_img(curr_filter, img_signal2);

curr_dendrite_left = [0; diff(img_signal1)];
curr_dendrite_right = [0; diff(img_signal2)];

[cross_corr_delta, ~] = xcov(curr_dendrite_left,curr_dendrite_right,length_corr,'coef');

[auto_corr_delta_dendrite_left, ~] = xcov(curr_dendrite_left,length_corr,'coef');
[auto_corr_delta_dendrite_right, time_lag] = xcov(curr_dendrite_right,length_corr,'coef');

curr_r = corrcoef(filt_img_dendrite1.fimg_1,filt_img_dendrite2.fimg_1);
corrcoef_infra = curr_r(2);

PRD = get_cross_period(cross_corr_delta, 15, 0, 0.05, length_corr);

pr_delta = PRD.curr_pr;
amp_delta = PRD.curr_amp;
curr_r = corrcoef(curr_dendrite_left,curr_dendrite_right);
corrcoef_delta = curr_r(2);

end

function PRD = get_cross_period(img_signal, min_dist, min_height, min_prom, length_corr)

[peak_amp, lag_dist] = findpeaks(img_signal, 'MinPeakDistance', min_dist, 'MinPeakProminence', min_prom,...
    'MinPeakHeight', min_height);
lag_dist = lag_dist-(length_corr+1);
lag_dist = abs(lag_dist);

if size(lag_dist,1) == 1
    PRD.curr_amp = NaN;
    PRD.curr_pr = NaN;
    
elseif isempty(lag_dist) == 1
    PRD.curr_amp = NaN;
    PRD.curr_pr = NaN;
else
    [nu_l, curr_idx] = sort(lag_dist,'ascend');
    PRD.curr_amp = peak_amp(curr_idx(1));
    PRD.curr_pr = nu_l(1)/14.56;
end

end