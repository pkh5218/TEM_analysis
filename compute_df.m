function df = compute_df(f)
    % Compute non-uniform dx differences
    m = length(f);
    df = zeros(1, m);

    % Forward difference for the first point
    df(1) = f(2) - f(1);

    % Central difference for the interior points
    for j = 2:m-1
        df(j) = (f(j+1) - f(j-1)) / 2;
    end

    % Backward difference for the last point
    df(m) = f(m) - f(m-1);

end