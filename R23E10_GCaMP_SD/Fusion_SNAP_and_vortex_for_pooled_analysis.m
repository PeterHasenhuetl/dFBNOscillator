% Merging data-tables of SD data obtained with SNAP and  vortex, for
% pooled analysis.
% Code written by Peter Hasenhuetl.


cross_corr_table = [cross_corr_table_SNAP; cross_corr_table_vortex];
mp_data_table = [mp_data_table_SNAP; mp_data_table_vortex];
curr_id = size(fly_names_SNAP,1)+1;
for loop_idx = (size(mp_data_table_SNAP,1)+1):4:size(mp_data_table,1)
    mp_data_table.fly_id(loop_idx:loop_idx+3) = curr_id;
    curr_id = curr_id+1;
end

