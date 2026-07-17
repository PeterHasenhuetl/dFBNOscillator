function get_default_separated_ax(ax_1, xlm_1, xlm_2, x_tick_start, x_tick_end, ylm_1, ylm_2,...
    y_tick_start, y_tick_end, xtick_space, ytick_space, x_scale, y_scale, x_label, y_label,...
    x_axis_color, y_axis_color, x_label_color, y_label_color, sz_1, ht_1)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

% First, defines the original axes limits etc. for plotting (these are the same as
% below, but with invisible axes)
get_default_ax(ax_1, xlm_1, xlm_2, x_tick_start, x_tick_end, ylm_1, ylm_2,...
    y_tick_start, y_tick_end, xtick_space, ytick_space, x_scale, y_scale, x_label, y_label,...
    'none', 'none', 'none', 'none', sz_1, ht_1)

if ax_1.YAxisLocation == "right"
    ax_dist = -0.100;
else
    ax_dist = 0.100;
end


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



%%

curr_pos = ax_1.Position;

ax_X = axes('Units','Centimeters','Position',[curr_pos(1), curr_pos(2)-ax_dist, curr_pos(3), curr_pos(4)]);

ax_X.FontSize = get_default_font_size;
ax_X.FontWeight = 'normal';
ax_X.FontName = 'Helvetica';
ax_X.TickDir = 'out';
ax_X.TickLength = [tick_factor*get_default_tick_length, tick_factor*get_default_tick_length];
ax_X.XLim = [xlm_1, xlm_2];
if Minor_tick_cond == "on"
    ax_X.XAxis.MinorTick = 'on';
    ax_X.XAxis.MinorTickValues = x_tick_start:(xtick_space/denom_x):x_tick_end;
else
    ax_X.XAxis.MinorTick = 'off';
end
ax_X.XTick = x_tick_start:xtick_space:x_tick_end;
ax_X.Box = 'off';
ax_X.LineWidth = get_default_axis_width;
ax_X.Color = 'none';    
ax_X.XAxis.Color = x_axis_color;
ax_X.YAxis.Color = 'none';
ax_X.XScale = char(x_scale);
ax_X.XLabel.String = string(x_label);
ax_X.XLabel.Color = x_label_color;
ax_X.XLabel.FontSize = get_default_font_size;
xtickangle(ax_X,0)

ax_Y = axes('Units','Centimeters','Position',[curr_pos(1)-ax_dist, curr_pos(2), curr_pos(3), curr_pos(4)]);

ax_Y.FontSize = get_default_font_size;
ax_Y.FontWeight = 'normal';
ax_Y.FontName = 'Helvetica';
ax_Y.TickDir = 'out';
ax_Y.TickLength = [tick_factor*get_default_tick_length, tick_factor*get_default_tick_length];
ax_Y.YLim = [ylm_1, ylm_2];
if Minor_tick_cond == "on"
    ax_Y.YAxis.MinorTick = 'on';
    ax_Y.YAxis.MinorTickValues = y_tick_start:(ytick_space/denom_y):y_tick_end;
else
    ax_Y.YAxis.MinorTick = 'off';
end
ax_Y.XTick = x_tick_start:xtick_space:x_tick_end;
ax_Y.YTick = y_tick_start:ytick_space:y_tick_end;
ax_Y.Box = 'off';
ax_Y.LineWidth = get_default_axis_width;
ax_Y.Color = 'none';    
ax_Y.XAxis.Color = 'none';
ax_Y.YAxis.Color = y_axis_color;
ax_Y.YScale = char(y_scale);
ax_Y.YLabel.String = string(y_label);
ax_Y.YLabel.Color = y_label_color;
ax_Y.YLabel.FontSize = get_default_font_size;
xtickangle(ax_Y,0)

if ax_1.YAxisLocation == "right"
    ax_Y.YAxisLocation = "right";
end


end


