function curr_plane = get_dFB_manual_exclusions(curr_plane, curr_idx)
% Implementing the exclusion of traces after inspection of motion artifacts
% (blinded to experimental condition). 
% Code written by Peter Hasenhuetl.

if curr_idx == 0
    curr_idx = [3, 4];
end
curr_idx = curr_idx+1; %(because background is another column!!!)

for i = 1:length(curr_idx)
    curr_idx(i)
    curr_plane.green(:,curr_idx(i)) = curr_plane.green(:,1);
    curr_plane.red(:,curr_idx(i)) = curr_plane.red(:,1);
end

end