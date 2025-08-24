function X = compute_X(u_3day, u_prime, v_prime, w_prime, u_bar, v_bar, w_bar, p, lat)
% Zonal-mean Friction X
%   F = (1/(a cos^2 φ)) ∂φ [ (ū v̄ + (u'v')) cos^2 φ ]
%       + ∂p [ ū w̄ + (u'w')̄]  - f v̄ + ∂ū/∂t
%
% Inputs
%   u_3days   : [Nx x Ny x Nlev x 3] = (prev, curr, next) 
%   u_star    : [Nx x Ny x Nlev]  
%   v_star    : [Nx x Ny x Nlev] 
%   w_star    : [Nx x Ny x Nlev] 
%   u_bar     : [Ny x Nlev] 
%   v_bar     : [Ny x Nlev]
%   w_bar     : [Ny x Nlev]
%   p         : [1 x Nlev]           (Pa)
%   lat       : [Ny x 1]             deg）
%
% Output
%   X         : [Ny x Nlev]

    % ----- constants -----
    a     = 6.371e6;     % m
    Omega = 7.292e-5;    % s^-1
    dt    = 86400;       % s/day

    % ----- geometry -----
    phi    = deg2rad(lat(:));    % Ny x 1
    cosphi = cos(phi);
    f      = 2*Omega*sin(phi);

    Ny   = numel(lat);
    Nlev = numel(p);

    % maps
    [~, cos_map] = meshgrid(p, cosphi);

    % ====== time tendency ∂ū/∂t ======
    u_prev = u_3day(:,:,:,1);
    u_curr = u_3day(:,:,:,2);
    u_next = u_3day(:,:,:,3);

    prev_is_halo = all(u_prev(:)==0);
    next_is_halo = all(u_next(:)==0);

    if prev_is_halo && ~next_is_halo
        dudt = (u_next - u_curr)/dt;          % forward at first day
    elseif ~prev_is_halo && next_is_halo
        dudt = (u_curr - u_prev)/dt;          % backward at last day
    else
        dudt = (u_next - u_prev)/(2*dt);      % central for normal days
    end
    dudt_bar = squeeze(mean(dudt,1));         % [Ny x Nlev]

    % ====== eddy fluxes (zonal mean for today) ======
    uv_flux_bar = squeeze(mean(u_prime .* v_prime, 1));   % [Ny x Nlev]
    uw_flux_bar = squeeze(mean(u_prime .* w_prime, 1));   % [Ny x Nlev]

    % ====== terms ======
    % (1) + (2):  (1/(a cos^2φ)) ∂φ [ (ū v̄ + (u'v')̄) cos^2 φ ]
    Aphi   = (u_bar .* v_bar + uv_flux_bar) .* (cos_map.^2);
    dA_dphi= compute_dfdy(Aphi, phi);
    term12 = (1/a) * dA_dphi ./ (cos_map.^2);

    % (3) + (4):  ∂p [ ū w̄ + (u'w')̄ ]
    Bp     = u_bar .* w_bar + uw_flux_bar;
    term34 = compute_dfdp(Bp, p);    

    % (5): - f v̄
    term5  = - f .* v_bar;

    % (6): ∂ū/∂t
    term6  = dudt_bar;

    % total
    X = term12 + term34 + term5 + term6;
end




