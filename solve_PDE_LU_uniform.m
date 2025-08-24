function F = solve_PDE_LU_uniform(A, B, C, D, E, H, phi, p, F_bottom)
% Solve A F_pp + 2B F_py + C F_yy + D F_p + E F_y = H
% on a uniform (φ,p) grid by sparse LU decomposition.
% Author : Pin-Chun Huang
% Date   : 2025-08-23
%
% Purpose:
%   Solve a 2D linear PDE on a uniform grid using a 9-point stencil:
%   - Second derivatives: 5-point (p,y)
%   - Cross derivative  : 4 corners (centered for F_py)
%   Boundary conditions:
%   - Boundary condition: F(1,:) = 0, F(Ny,:) = 0, F(:,1) = 0, F(:,Np) = F_bottom
%
% Inputs:
%   A, B, C, D, E   : coefficient for L operator
%   H               : [Ny x Np] forcing
%   phi, p          : coordinate array (uniform spacing; phi in rad, p in Pa)
%   F_bottom        : [Ny x 1] p = p_max 
%
% Output:
%   F           : [Ny x Np] is the TEM streamfunction response for H

% Number of latitudes and pressure levels
Ny = length(phi);
Np = length(p);

% Compute uniform grid spacing
dphi = phi(2) - phi(1);   % rad
dp   = p(2) - p(1);       % Pa

% Initialize sparse coefficient matrix and RHS vector
N            = Ny * Np;
coeff_matrix = spalloc(N, N, 9 * N); % Estimate 9 non-zero elements per row
rhs          = zeros(N, 1);

%% Fill the coefficient matrix and RHS vector
for i = 1:Ny
    for j = 1:Np
        idx = (i - 1) * Np + j;
        % Boundary condition
        if i == 1 || i == Ny || j == 1
           coeff_matrix(idx, idx) = 1;
           rhs(idx) = 0;
           continue
        elseif j == Np
           coeff_matrix(idx, idx) = 1;
           rhs(idx) = F_bottom(i);
           continue
        end

        % ----- interior stencil -----
        % Diagonal element
        coeff_matrix(idx, idx) = -2 * (A(i,j)/dp^2 + C(i,j)/dphi^2);

        % A * F_pp + D * F_p
        coeff_matrix(idx, idx - 1) = A(i,j)/dp^2 - D(i,j)/(2*dp);
        coeff_matrix(idx, idx + 1) = A(i,j)/dp^2 + D(i,j)/(2*dp);

        % C * F_φφ + E * F_φ
        coeff_matrix(idx, idx - Np) = C(i,j)/dphi^2 - E(i,j)/(2*dphi);
        coeff_matrix(idx, idx + Np) = C(i,j)/dphi^2 + E(i,j)/(2*dphi);

        % 2B * F_pφ
        coeff_matrix(idx, idx - Np - 1) = coeff_matrix(idx, idx - Np - 1) + B(i,j)/(4*dphi*dp);
        coeff_matrix(idx, idx + Np + 1) = coeff_matrix(idx, idx + Np + 1) + B(i,j)/(4*dphi*dp);
        coeff_matrix(idx, idx - Np + 1) = coeff_matrix(idx, idx - Np + 1) - B(i,j)/(4*dphi*dp);
        coeff_matrix(idx, idx + Np - 1) = coeff_matrix(idx, idx + Np - 1) - B(i,j)/(4*dphi*dp);

        % RHS
        rhs(idx) = H(i, j);
    end
end

    % Solve the linear system using LU decomposition for sparse matrices
    [L, U, P] = lu(coeff_matrix);
    y = L \ (P * rhs);
    F_vec = U \ y;
    
    % Convert solution vector back to matrix form
    F = reshape(F_vec, [Np, Ny])';
end
