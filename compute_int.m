function [var_int, p_new] = compute_int(var, p_old, p_start, p_end, Nlev_new) 
% Interpolate to (nearly) uniform pressure levels. 
% Author : Pin-Chun Huang 
% Date : 2025-08-23 
% 
% Purpose: 
% Interpolates each latitude column of VAR from original pressure 
% levels (Pa) to a nearly uniform grid between p_start and p_end 
% with N_new_levels points. Top level is forced to 100 Pa (1 hPa) 
% to avoid p = 0 issues. Returns p_new in Pa. 
% Uses quadric spline (order=5) so that derivatives up to 3rd order are 
% smooth so derivatives up to 3rd order are smooth later. 
% 
% Inputs: 
% var : [Ny x Nlev_old] % p_old : [1 x Nlev_old], original pressure level in Pa 
% Nlev_new : Number of new levels (set 81 to let a single level to be 12.5 hPa) 
% p_start : Lower Boundary, set 0 % p_end : Upper Boundary, set 1e5 (Pa) 
% 
% Outputs: 
% var_int : [Ny x Nlev_new], Interpolated variable 
% p_new : [1 x Nlev_new], Interpolated new pressure levels（Pa） 
% 
% Notes: 
% - Top level is set to 1 Pa (100 Pa), so the grid is "nearly" uniform. 

[Ny, ~] = size(var); 
% ---- target grid ---- 
p_new = linspace(p_start, p_end, Nlev_new); 
p_new(1) = 100; % avoid 0 Pa 

% ---- allocate ---- 
var_int = zeros(Ny, Nlev_new); 

% ---- interpolate each latitude column ----
%{
for i = 1:Ny 
    yi = var(i,:); 
    sp = spapi(5, p_old.', yi); % quartic spline (order=5, degree=4, C^3 continuity) 
    var_int(i,:) = fnval(sp, p_new);
end
%}
for i = 1:Ny
    var_int(i, :) = interp1(p_old, var(i, :), p_new, 'pchip', 'extrap');
end
end