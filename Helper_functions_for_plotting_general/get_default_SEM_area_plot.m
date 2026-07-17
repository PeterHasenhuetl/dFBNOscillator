function get_default_SEM_area_plot(curr_ax, curr_traces, curr_tv, curr_color)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

% Input curr_traces as matrix, where each column is one recording
% Input time vector as row vector


sz_1 = sqrt(size(curr_traces,2)); % Square-root of N for computing SEM from STD
STD_1 = std(curr_traces,0,2); % Using the standard normalization for STD (i.e. "N-1")
SEM = STD_1./sz_1;

trace_area = [(mean(curr_traces,2)-SEM)',...
    fliplr((mean(curr_traces,2)+SEM)')];

time_vector = [curr_tv,fliplr(curr_tv)];

hold on
fill(curr_ax, time_vector,trace_area,curr_color,'EdgeColor','none','FaceAlpha',0.5)
plot(curr_ax, curr_tv, mean(curr_traces,2),'Color',curr_color,'LineWidth',0.5,'LineStyle','-')

end