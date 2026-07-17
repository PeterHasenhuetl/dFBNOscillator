function filt_img = get_filtered_mp_img(curr_filter,curr_img)
% Filtering of dFBN GCaMP dynamics.
% Code written by Peter Hasenhuetl.

fimg_1 = filter(curr_filter.filt_1,[curr_img; zeros(curr_filter.dl_1,1)]); 
filt_img.fimg_1 = fimg_1(curr_filter.dl_1+1:end,1);

fimg_2 = filter(curr_filter.filt_2,[curr_img; zeros(curr_filter.dl_2,1)]); 
filt_img.fimg_2 = fimg_2(curr_filter.dl_2+1:end,1);

end