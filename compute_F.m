function F = compute_F(u, u_star, v_star, w_star, u_ave, v_ave, w_ave, p, lat)
R       = 6371000;
Omega   = 7.292e-5;
N_t     = size(u,4);
phi     = deg2rad(lat);
cosphi  = cos(phi);
sinphi  = sin(phi);
f       = 2*Omega*sinphi;

m_flux  = u_star.*v_star;
w_flux  = u_star.*w_star;
m_flux_ave  = squeeze(mean(m_flux,[1,4],'omitnan'));
w_flux_ave  = squeeze(mean(w_flux,[1,4],'omitnan'));

dt = 86400;
if N_t>=3
dudt(:,:,:,2:N_t-1) = (u(:,:,:,3:N_t)-u(:,:,:,1:N_t-2))/(2*dt);
dudt(:,:,:,1) = (u(:,:,:,2)-u(:,:,:,1))/dt;
dudt(:,:,:,N_t) = (u(:,:,:,N_t)-u(:,:,:,N_t-1))/dt;
else
    if N_t == 2
        dudt = (u(:,:,:,2)-u(:,:,:,1))/dt;
    end
end
dudt = squeeze(mean(dudt,[1,4],'omitnan'));

[~,cosphi] = meshgrid(p,cosphi);
term1 = u_ave.*v_ave.*cosphi.*cosphi;
term1 = compute_dfdy(term1,phi);
term1 = term1/R./cosphi./cosphi;

term2 = m_flux_ave.*cosphi.*cosphi;
term2 = compute_dfdy(term2,phi);
term2 = term2/R./cosphi./cosphi;

term3 = compute_dfdp(u_ave.*w_ave,p);

term4 = compute_dfdp(w_flux_ave,p);

term5 = -f.*v_ave;

term6 = dudt;
F = term1 + term2 + term3 + term4 + term5 + term6;

end



