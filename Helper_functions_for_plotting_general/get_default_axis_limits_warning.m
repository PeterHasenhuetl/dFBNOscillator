function get_default_axis_limits_warning(input_data, axis_limits)
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.

if min(input_data) <= axis_limits(1)
    exceeded_axis_limit = axis_limits(1)
    at_least_out_of_range = min(input_data)
    warning('!!!SOME DATA MAY BE OUTSIDE OF THE AXIS LIMITS!!!')
    pause
end

if max(input_data) >= axis_limits(2)
    exceeded_axis_limit = axis_limits(2)
    at_least_out_of_range = max(input_data)
    warning('!!!SOME DATA MAY BE OUTSIDE OF THE AXIS LIMITS!!!')
    pause
end

end