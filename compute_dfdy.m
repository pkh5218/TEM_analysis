function dfdy = compute_dfdy(f,phi)
    % Compute non-uniform dx differences
    m = size(f,1);
    dfdy = f;

    % Forward difference for the first point
    dfdy(1,:) = f(2,:) - f(1,:);

    % Central difference for the interior points
    for j = 2:m-1
        dfdy(j,:) = (f(j+1,:) - f(j-1,:)) / 2;
    end

    % Backward difference for the last point
    dfdy(m,:) = f(m,:) - f(m-1,:);
    
    dphi = compute_df(phi);

    for j = 1:m
        dfdy(j,:) = dfdy(j,:)/dphi(j);
    end

end