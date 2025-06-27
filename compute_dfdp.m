function dfdp = compute_dfdp(f,p)
    % Compute non-uniform dx differences
    m = size(f,2);
    dfdp = f;

    % Forward difference for the first point
    dfdp(:,1) = f(:,2) - f(:,1);

    % Central difference for the interior points
    for j = 2:m-1
        dfdp(:,j) = (f(:,j+1) - f(:,j-1)) / 2;
    end

    % Backward difference for the last point
    dfdp(:,m) = f(:,m) - f(:,m-1);
    
    dp = compute_df(p);

    for j = 1:m
        dfdp(:,j) = dfdp(:,j)/dp(j);
    end

end