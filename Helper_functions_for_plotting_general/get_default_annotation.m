function get_default_annotation(x_pos, y_pos, my_text, text_color, font_angle, horizontal_alignment)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

if horizontal_alignment == 'right'
    x_pos = x_pos-10;
elseif horizontal_alignment == 'center'
    x_pos = x_pos-5;
end

annotation('textbox', 'Units', 'centimeters', 'Position', [x_pos, y_pos-0.25, 10, 0.25], ...
    'string', string(my_text), 'EdgeColor', 'none',...
    'FontSize', get_default_font_size, 'FontWeight', 'normal', 'Color', text_color,...
    'FontAngle', font_angle, 'Margin', 0,...
    'HorizontalAlignment', horizontal_alignment)

end