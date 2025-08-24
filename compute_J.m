function [J, J_map] = compute_J(u_3day, v_3day, w_3day, T_3day, p, lon, lat, m, day, Nyears, NDays)
% Zonal-mean Diabatic Heating J and the whole map 3-D data J_map
%   F = (1/(a cos^2 φ)) ∂φ [ (ū v̄ + (u'v')) cos^2 φ ]
%       + ∂p [ ū w̄ + (u'w')̄]  - f v̄ + ∂ū/∂t
%
% Inputs
%   u_3days   : [Nx x Ny x Nlev x 3] = (prev, curr, next) 
%   v_3days   : [Nx x Ny x Nlev x 3] = (prev, curr, next)
%   w_3days   : [Nx x Ny x Nlev x 3] = (prev, curr, next)
%   T_3days   : [Nx x Ny x Nlev x 3] = (prev, curr, next)

%   p         : [1 x Nlev]           (Pa)
%   lon       : [Nx x 1]             (deg）
%   lat       : [Ny x 1]             (deg）
%   m         : current running year index
%   day       : current running day index
%   Nyears    : total number of years 
%   Ndays     : total number of days
%
% Output
%   X         : [Ny x Nlev]

    Rd  = 287; 
    cp  = 1004; 
    a   = 6.371e6; 
    dt  = 86400; % s/day
    [Nx,Ny,Nlev,~] = size(T_3day);   % 這裡 Nt 應該就是 3

    phi    = deg2rad(lat);
    cosphi = cos(phi);

    % ---- dT/dt ----
    if (m==1 && day==1)
        dTdt = (T_3day(:,:,:,3) - T_3day(:,:,:,2)) / dt;
    elseif (m==Nyears && day==NDays)
        dTdt = (T_3day(:,:,:,2) - T_3day(:,:,:,1)) / dt;
    else
        dTdt = (T_3day(:,:,:,3) - T_3day(:,:,:,1)) / (2*dt);
    end

    % ---- Space Differentiation ----
    dx = deg2rad(lon(2)-lon(1));
    dy = phi(2)-phi(1);
    dp = p(3)-p(2);

    dTdx = zeros(Nx,Ny,Nlev);
    dTdx(2:Nx-1,:,:) = (T_3day(3:Nx,:,:,2) - T_3day(1:Nx-2,:,:,2)) / (2*dx);
    dTdx(1,:,:)      = (T_3day(2,:,:,2)    - T_3day(Nx,:,:,2))     / (2*dx);
    dTdx(Nx,:,:)     = (T_3day(1,:,:,2)    - T_3day(Nx-1,:,:,2))   / (2*dx);

    dTdy = zeros(Nx,Ny,Nlev);
    dTdy(:,2:Ny-1,:) = (T_3day(:,3:Ny,:,2) - T_3day(:,1:Ny-2,:,2)) / (2*dy);
    dTdy(:,1,:)      = (T_3day(:,2,:,2)    - T_3day(:,1,:,2))      / dy;
    dTdy(:,Ny,:)     = (T_3day(:,Ny,:,2)   - T_3day(:,Ny-1,:,2))   / dy;

    dTdp = zeros(Nx,Ny,Nlev);
    dTdp(:,:,1)   = (T_3day(:,:,2,2) - T_3day(:,:,1,2)) / dp;
    for k = 2:Nlev-1
        dTdp(:,:,k) = (T_3day(:,:,k+1,2) - T_3day(:,:,k-1,2)) / (2*dp);
    end
    dTdp(:,:,Nlev)= (T_3day(:,:,Nlev,2) - T_3day(:,:,Nlev-1,2)) / dp;

    % ---- J_map ----
    J_map = zeros(Nx,Ny,Nlev);
    for j = 1:Ny
        for k = 1:Nlev
        J_map(:,j,k) = dTdt(:,j,k) ...
          + u_3day(:,j,k,2).*dTdx(:,j,k)/(a*cosphi(j)) ...
          + v_3day(:,j,k,2).*dTdy(:,j,k)/a ...
          + w_3day(:,j,k,2).*dTdp(:,j,k) ...
          - w_3day(:,j,k,2).*T_3day(:,j,k,2)/p(k)*(Rd/cp);
        end
    end

    % zonal mean
    J = squeeze(mean(J_map,1));  % [Ny x Nlev]
end