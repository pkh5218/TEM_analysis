function [J,J_map] = compute_J(u, v, w, T, p, lon, lat)
Rd      = 287;
cp      = 1004;
R       = 6371000;

N_p     = length(p);
N_x     = length(lon);
N_y     = length(lat);
N_t     = size(T,4);
phi     = deg2rad(lat);
cosphi  = cos(phi);

dTdt    = T;
dTdx    = T;
dTdy    = T;
dTdp    = T;
J_map   = T;

dt = 86400;
if N_t>=3
dTdt(:,:,:,2:N_t-1) = (T(:,:,:,3:N_t)-T(:,:,:,1:N_t-2))/(2*dt);
dTdt(:,:,:,1) = (T(:,:,:,2)-T(:,:,:,1))/dt;
dTdt(:,:,:,N_t) = (T(:,:,:,N_t)-T(:,:,:,N_t-1))/dt;
else
    if N_t == 2
        dTdt = (T(:,:,:,2)-T(:,:,:,1))/dt;
    end
end

dx = deg2rad(lon(2)-lon(1));
dTdx(2:N_x-1,:,:,:) = (T(3:N_x,:,:,:)-T(1:N_x-2,:,:,:))/(2*dx);
dTdx(1,:,:,:) = (T(2,:,:,:)-T(N_x,:,:,:))/(2*dx);
dTdx(N_x,:,:,:) = (T(1,:,:,:)-T(N_x-1,:,:,:))/(2*dx);

dy = phi(2)-phi(1);
dTdy(:,2:N_y-1,:,:) = (T(:,3:N_y,:,:)-T(:,1:N_y-2,:,:))/(2*dy);
dTdy(:,1,:,:) = (T(:,2,:,:)-T(:,1,:,:))/dy;
dTdy(:,N_y,:,:) = (T(:,N_y,:,:)-T(:,N_y-1,:,:))/dy;

dp = compute_df(p);
dTdp(:,:,1,:) = (T(:,:,2,:) - T(:,:,1,:))/dp(1);
for j = 2:N_p-1
    dTdp(:,:,j,:) = (T(:,:,j+1,:) - T(:,:,j-1,:))/(2*dp(j));
end
dTdp(:,:,N_p,:) = (T(:,:,N_p,:) - T(:,:,N_p-1,:))/dp(N_p);

for j = 1:N_y
    for k = 1:N_p
        J_map(:,j,k,:) = dTdt(:,j,k,:) + u(:,j,k,:).*dTdx(:,j,k,:)/R/cosphi(j) ...
            + v(:,j,k,:).*dTdy(:,j,k,:)/R + w(:,j,k,:).*dTdp(:,j,k,:) ...
            - w(:,j,k,:).*T(:,j,k,:)/p(k)*Rd/cp;
    end
end

J = squeeze(mean(J_map,[1,4]));
