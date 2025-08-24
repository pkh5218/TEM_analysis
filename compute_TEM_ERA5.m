function [eumpsi, tempsi] = compute_TEM_ERA5(v_bar, EHF_dthdp, p, cos_map)
% COMPUTE_TEM_ERA5
% Eulerian/TEM streamfunction:
%   psi_EUM(φ,p) = ∫ (2πa cosφ / g) * v̄  dp'
%   psi_TEM(φ,p) = ∫ (2πa cosφ / g) * v* dp'
%   with v* = v̄ - ∂p(EHF/dthdp)
%
% Inputs:
%   v_bar       : [Ny x Np]
%   EHF_dthdp   : [Ny x Np]    EHF/dthdp
%   p           : [1 x Np]         (Pa)    
%   cos_map     : [Ny x Np]
%
% Outputs:
%   eumpsi : [Ny x Np]        Eulerian-mean streamfunction
%   tempsi : [Ny x Np]        TEM streamfunction

g = 9.8;
a = 6371000;

% v* = v̄ - ∂(EHF/dthdp)/∂p
v_star = v_bar - compute_dfdp(EHF_dthdp, p);

% Metric factor on the (lat, p) grid
M = (2*pi*a/g) * cos_map;                      % [Ny x Np]

% Cumulative trapezoidal integral along pressure (2nd dim)
% p is integrated from top to bottom with top to be 0.
eumpsi = cumtrapz(p, M .* v_bar,  2);          % [Ny x Np]
tempsi = cumtrapz(p, M .* v_star, 2);          % [Ny x Np]
end