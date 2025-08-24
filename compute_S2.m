function S2 = compute_S2(T_bar, p_map)
% Static stability S^2 in pressure coordinates.
% S2(φ,p) = -(Rd/p) * (∂T̄/∂p - (Rd/cp) * T̄/p)
% Inputs:
%   T_bar : [Ny x Np]  zonal-mean temperature (K)
%   p_map : [Ny  x Np] pressure levels (Pa)
% Output:
%   S2    : [Ny x Np]

Rd = 287;
cp = 1004;

dTdp = compute_dfdp(T_bar, p_map(1,:));
S2 = -(Rd ./ p_map) .* (dTdp - (Rd/cp) .* (T_bar ./ p_map));

end