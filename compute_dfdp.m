function dfdp = compute_dfdp(f, p)
% f: [Ny x Np], p: Pressure level (Pa)
    dp = p(3) - p(2); 
    dfdp = zeros(size(f));
    % Interior (Central Difference)
    dfdp(:,2:end-1) = (f(:,3:end) - f(:,1:end-2)) / (2*dp);
    % Boundaries (Forward/Backward Difference)
    dfdp(:,1) = (f(:,2) - f(:,1)) / dp;
    dfdp(:,end) = (f(:,end) - f(:,end-1)) / dp;
end