function curr_filter = get_mp_osc_filt(s_rt, f1_1, f2_1, f1_2, f2_2)
% Filtering of dFBN GCaMP dynamics.
% Code written by Peter Hasenhuetl.

filter_order = 40;

curr_filter.filt_1 = designfilt('bandpassfir','FilterOrder',filter_order, ...
         'CutoffFrequency1',f1_1,'CutoffFrequency2',f2_1, ...
         'SampleRate',s_rt);
curr_filter.dl_1 = mean(grpdelay(curr_filter.filt_1));

curr_filter.filt_2 = designfilt('bandpassfir','FilterOrder',filter_order, ...
         'CutoffFrequency1',f1_2,'CutoffFrequency2',f2_2, ...
         'SampleRate',s_rt);
curr_filter.dl_2 = mean(grpdelay(curr_filter.filt_2));

end