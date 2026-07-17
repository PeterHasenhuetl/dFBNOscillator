function [C, W] = get_fig_panel_survival
% Plotting survival of R23E10-GAL4 flies after implantation of chronic imaging windows.
% Code written by Peter Hasenhuetl.


load('n_windows.mat');
C = (n_controls(:,2)./n_controls(:,1)).*100;
W = (n_windows(:,2)./n_windows(:,1)).*100;

x_pos = 5;
y_pos = 5;
color = get_color;
sz_1 = 2;
ht_1 = 1.8;
pnl_1 = axes('Units', 'Centimeters', 'Position', [x_pos, y_pos, sz_1, ht_1]);
bar(pnl_1, [1, 2], [mean(C), mean(W)], 0.5, 'EdgeColor', [0, 0, 0],...
    'FaceColor', color.light_gray, 'ShowBaseLine', 'off')
hold on
get_default_scatter_group(pnl_1, C, 1, [0, 0, 0])
get_default_scatter_group(pnl_1, W, 2, [0, 0, 0])
xlm1 = [0.5, 2.5];
get_default_ax(pnl_1, xlm1(1), xlm1(2), [], [], 0, 100, 0, 100,...
    [], 50, "linear", "linear", [], '% alive',...
    'none', 'k', 'k', 'k', sz_1, ht_1)
get_default_scatter_n_numbers_annotation(x_pos, y_pos-0.05, sz_1, ...
    [sum(n_controls(:,1)), sum(n_windows(:,1))], 'flies', 'k', 'italic')
get_default_scatter_n_numbers_annotation(x_pos, y_pos-0.25, sz_1, ...
    [length(n_controls(:,1)), length(n_windows(:,1))], 'replicates', 'k', 'italic')
get_default_annotation(x_pos, y_pos+ht_1+0.4, 'Imaging window', 'k', 'italic', "right")
get_default_annotation(x_pos+0.35, y_pos+ht_1+0.4, 'no', 'k', 'normal', "left")
get_default_annotation(x_pos+1.35, y_pos+ht_1+0.4, 'yes', 'k', 'normal', "left")


mean(C)
mean(W)

[p, b] = ranksum(C,W)
get_default_p_value_annotation(p,"asterisk", x_pos+0.3, y_pos+ht_1+0.75)

end
