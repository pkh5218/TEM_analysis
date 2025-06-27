function [var_int,p_new] = compute_int(var, p_ini, p_start, p_end, N_new_levels)
    [Ny, ~] = size(var);

    p_new = linspace(p_start*100, p_end*100, N_new_levels);
    p_new(1) = 100;

    var_int = zeros(Ny, N_new_levels);
    for i = 1:Ny
        var_int(i, :) = interp1(p_ini, var(i, :), p_new, 'spline', 'extrap');
    end
end