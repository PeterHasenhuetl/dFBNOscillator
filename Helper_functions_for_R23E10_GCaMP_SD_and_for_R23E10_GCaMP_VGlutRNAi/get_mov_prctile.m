function mov_prctile = get_mov_prctile(raw_trace, sliding_window)
% Computes the moving 10th percentile of the dFBN GCaMP trace to use as F0.
% Code written by Peter Hasenhuetl.


% if sliding window is even, adds 1 (for symmetric sliding window)
is_even = (sliding_window/2) == round(sliding_window/2);
if is_even == 1
    sliding_window = sliding_window+1;
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
    window_center_idx = window_center_idx+1;
end

end