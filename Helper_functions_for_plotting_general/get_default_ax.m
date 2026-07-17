function get_default_ax(ax_1, xlm_1, xlm_2, x_tick_start, x_tick_end, ylm_1, ylm_2,...
    y_tick_start, y_tick_end, xtick_space, ytick_space, x_scale, y_scale, x_label, y_label,...
    x_axis_color, y_axis_color, x_label_color, y_label_color, sz_1, ht_1)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

Minor_tick_cond = "off"; % Condition of whether to display minor ticks

norm_f = max([sz_1, ht_1]); % Factor to normalize tick length to panel size
tick_factor = 1/norm_f;

% Defining the number of ninor ticks, depending on the values between major ticks
multipl1 = [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 1000, 10000]; 
cond_four_minor = [0.2, 0.4, 0.8].*multipl1';
cond_five_minor = [0.125, 0.25, 0.5, 0.75, 1].*multipl1';
cond_three_minor = [0.15, 0.3, 0.6, 0.9].*multipl1';

if ismember(xtick_space, cond_four_minor) == 1
    denom_x = 4;
elseif ismember(xtick_space, cond_five_minor) == 1
    denom_x = 5;    
elseif ismember(xtick_space, cond_three_minor) == 1
    denom_x = 3;    
elseif isempty(xtick_space) == 1
    denom_x = 1;
else
    denom_x = 5;
end


if ismember(ytick_space, cond_four_minor) == 1
    denom_y = 4;
elseif ismember(ytick_space, cond_five_minor) == 1
    denom_y = 5;    
elseif ismember(ytick_space, cond_three_minor) == 1
    denom_y = 3;    
elseif isempty(ytick_space) == 1
    denom_y = 1;
else
    denom_y = 5;
end

ax_1.FontSize = get_default_font_size;
ax_1.FontWeight = 'normal';
ax_1.FontName = 'Helvetica';
ax_1.TickDir = 'out';
ax_1.TickLength = [tick_factor*get_default_tick_length, tick_factor*get_default_tick_length];
ax_1.XLim = [xlm_1, xlm_2];
ax_1.YLim = [ylm_1, ylm_2];
if Minor_tick_cond == "on"
    ax_1.XAxis.MinorTick = 'on';
    ax_1.YAxis.MinorTick = 'on';
    ax_1.XAxis.MinorTickValues = x_tick_start:(xtick_space/denom_x):x_tick_end;
    ax_1.YAxis.MinorTickValues = y_tick_start:(ytick_space/denom_y):y_tick_end;
else
    ax_1.XAxis.MinorTick = 'off';
    ax_1.YAxis.MinorTick = 'off';
end
ax_1.XTick = x_tick_start:xtick_space:x_tick_end;
ax_1.YTick = y_tick_start:ytick_space:y_tick_end;
ax_1.Box = 'off';
ax_1.LineWidth = get_default_axis_width;
ax_1.Color = 'none';    
ax_1.XAxis.Color = x_axis_color;
ax_1.YAxis.Color = y_axis_color;
ax_1.XScale = char(x_scale);
ax_1.YScale = char(y_scale);
ax_1.XLabel.String = string(x_label);
ax_1.YLabel.String = string(y_label);
ax_1.XLabel.Color = x_label_color;
ax_1.YLabel.Color = y_label_color;
ax_1.XLabel.FontSize = get_default_font_size;
ax_1.YLabel.FontSize = get_default_font_size;
xtickangle(ax_1,0)
% ax_1.XLabel.Rotation = 0;
% ax_1.YLabel.Rotation = 0;

end


