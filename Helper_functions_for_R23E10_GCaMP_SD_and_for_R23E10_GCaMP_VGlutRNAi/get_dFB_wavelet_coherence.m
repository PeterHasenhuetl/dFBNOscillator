function [w_coh, f] = get_dFB_wavelet_coherence(input1, input2)
% Computing coherence.
% Code written by Peter Hasenhuetl.

x = input1;
y = input2;
[wcoh, ~, f] = wcoherence(x, y, 14.56, 'FrequencyLimits', [0.1, 2]);
w_coh = mean(wcoh,2);

end