function [dF, dF_contra, table_idx] = get_fig_panel_interhemisphere_images_example(mp_data_table)
% Plotting an example image of a bilateral recording, for display only.
% For quantification of dF/F traces, see "get_mp_img_measures".
% Code written by Peter Hasenhuetl.

for loop_idx1 = 1:size(mp_data_table,1)
    a(loop_idx1,1) = string(mp_data_table.details(loop_idx1,1).flyname) == 'fly20210601_7DFB_C' && mp_data_table.plane_id(loop_idx1) == 2;
end
table_idx = find(a == 1);
example_table =  mp_data_table(table_idx,:)

curr_burst_1 = mp_data_table.transient_vec(table_idx,1).trace_img(:,3);
curr_burst_2 = mp_data_table.transient_vec(table_idx,1).trace_img(:,4);

curr_trace_1 = mp_data_table.img_trace(table_idx,1).trace_img(:,3);
curr_trace_2 = mp_data_table.img_trace(table_idx,1).trace_img(:,4);

transient_id_cond = "actual transient times";
if transient_id_cond == "actual transient times"
    idx_1 = find(curr_burst_1 > 0);
    idx_1 = idx_1(2:end-1);
    idx_2 = find(curr_burst_2 > 0);
    idx_2 = idx_2(2:end-1);

elseif transient_id_cond == "random transient times" % Shuffled, as quality control
    idx_1 = randperm(19350, length(idx_1));
    idx_1 = idx_1+50;
    idx_2 = randperm(19350, length(idx_2));
    idx_2 = idx_2+50;
end

example_recording = dir('*fly20210601_7DFB_C_plane2_cropped*');
input_file = example_recording.name;
start_idx = 0;
video_idx = 1;
length_video_full = 20000;
for loop_idx2 = 501:length_video_full
    video_full(:,:,video_idx) = imread(input_file,loop_idx2);
    video_idx = video_idx+1;
end
video_full = double(video_full(:,:,:));
background_1 = mean(mean(mean(video_full(1:20,1:20),3,'omitnan')));
video_full = video_full(21:end,:,:) - background_1;
mv_1 = mean(video_full,3);
mv_1(mv_1 < 70) = 0; % Thresholded for display only.


clear curr_dF_video curr_F0_video
for loop_idx3 = 1:length(idx_1)-1
    curr_dF_video(:,:,loop_idx3) = mean(video_full(:,:,idx_1(loop_idx3)+start_idx:idx_1(loop_idx3)+10),3);
    curr_F0_video(:,:,loop_idx3) = mean(video_full(:,:,idx_1(loop_idx3)-15:idx_1(loop_idx3)),3);    
end

curr_dF_video = mean(curr_dF_video,3);
curr_F0_video = mean(curr_F0_video,3);
curr_dF = (curr_dF_video-curr_F0_video)./curr_F0_video;
curr_dF(mv_1 == 0) = 0;    
dF = curr_dF;


clear curr_dF_video curr_F0_video
for loop_idx4 = 1:length(idx_2)-1
    curr_dF_video(:,:,loop_idx4) = mean(video_full(:,:,idx_2(loop_idx4)+start_idx:idx_2(loop_idx4)+10),3);
    curr_F0_video(:,:,loop_idx4) = mean(video_full(:,:,idx_2(loop_idx4)-15:idx_2(loop_idx4)),3);  
end

curr_dF_video = mean(curr_dF_video,3);
curr_F0_video = mean(curr_F0_video,3);
curr_dF = (curr_dF_video-curr_F0_video)./curr_F0_video;  
curr_dF(mv_1 == 0) = 0;    
dF_contra = curr_dF;

clear video_full

end