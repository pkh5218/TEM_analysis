function F = solve_PDE_LU_sparse(A, B, C, D, E, H, y, p, F_boundary)
    % Number of latitudes and pressure levels
    N_y = length(y);
    N_p = length(p);
    
    % Compute non-uniform grid spacing
    dy = compute_df(y);
    dp = compute_df(p);
    
    % Initialize sparse coefficient matrix and RHS vector
    coeff_matrix = spalloc(N_y * N_p, N_y * N_p, 5 * N_y * N_p); % Estimate 5 non-zero elements per row
    rhs = zeros(N_y * N_p, 1);
    
    % Fill the coefficient matrix and RHS vector
    for i = 1:N_y
        for j = 1:N_p
            idx = (i - 1) * N_p + j;  % Convert 2D index to 1D index
            if i == 1 || i == N_y || j == N_p
                % Boundary condition: F(1, :) = 0, F(N_y, :) = 0, and F(:, N_p) = 0
                coeff_matrix(idx, idx) = 1;
                rhs(idx) = 0;
            else
            
            % Diagonal element
            coeff_matrix(idx, idx) = -2 * (A(i, j) / dp(j)^2 + C(i, j) / dy(i)^2);
            
            % A*F_pp term
            if j > 1
                coeff_matrix(idx, idx - 1) = A(i, j) / dp(j)^2;
            end
            if j < N_p
                coeff_matrix(idx, idx + 1) = A(i, j) / dp(j)^2;
            end
            
            % C*F_yy term
            if i > 1
                coeff_matrix(idx, idx - N_p) = C(i, j) / dy(i)^2;
            end
            if i < N_y
                coeff_matrix(idx, idx + N_p) = C(i, j) / dy(i)^2;
            end
            
            % 2B*F_py term
            if i > 1 && j > 1
                coeff_matrix(idx, idx - N_p - 1) = coeff_matrix(idx, idx - N_p - 1) + B(i, j) / (4 * dy(i) * dp(j));
            end
            if i < N_y && j < N_p
                coeff_matrix(idx, idx + N_p + 1) = coeff_matrix(idx, idx + N_p + 1) + B(i, j) / (4 * dy(i) * dp(j));
            end
            if i > 1 && j < N_p
                coeff_matrix(idx, idx - N_p + 1) = coeff_matrix(idx, idx - N_p + 1) - B(i, j) / (4 * dy(i) * dp(j));
            end
            if i < N_y && j > 1
                coeff_matrix(idx, idx + N_p - 1) = coeff_matrix(idx, idx + N_p - 1) - B(i, j) / (4 * dy(i) * dp(j));
            end
            
            % D*F_p term
            if j > 1
                coeff_matrix(idx, idx - 1) = coeff_matrix(idx, idx - 1) - D(i, j) / (2 * dp(j));
            end
            if j < N_p
                coeff_matrix(idx, idx + 1) = coeff_matrix(idx, idx + 1) + D(i, j) / (2 * dp(j));
            end
            
            % E*F_y term
            if i > 1
                coeff_matrix(idx, idx - N_p) = coeff_matrix(idx, idx - N_p) - E(i, j) / (2 * dy(i));
            end
            if i < N_y
                coeff_matrix(idx, idx + N_p) = coeff_matrix(idx, idx + N_p) + E(i, j) / (2 * dy(i));
            end
            
            % RHS
            rhs(idx) = H(i, j);
            end
        end
    end
    
    % Apply boundary condition F(:, N_p) = F_boundary
    
    for i = 1:N_y
        idx = (i - 1) * N_p + N_p;
        coeff_matrix(idx, :) = 0;
        coeff_matrix(idx, idx) = 1;
        rhs(idx) = F_boundary(i);
    end
    
    % Solve the linear system using LU decomposition for sparse matrices
    [L, U, P] = lu(coeff_matrix);
    y = L \ (P * rhs);
    F_vector = U \ y;
    
    % Convert solution vector back to matrix form
    F = reshape(F_vector, [N_p, N_y])';
end