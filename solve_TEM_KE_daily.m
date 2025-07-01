function [FTEM, FJ, FF, FEPphi_QG, FEPphi1, FEPp_QG, FEPp1, FEPp2, ...
    FGphi1, FGphi2, FGp1, FGp2, tempsi, eumpsi, ...
    HJ, HF, ...
    HEPDphi_QG, HEPDphi1, HEPDp_QG, HEPDp1, HEPDp2, ...
    HGDphi1, HGDphi2, HGDp1, HGDp2, ...
    EPDphi_QG, EPDphi1, EPDp_QG, EPDp1, EPDp2, ...
    Fphi_QG, Fphi1, Fp_QG, Fp1, Fp2, ...
    Gphi1, Gphi2, Gp1, Gp2, ...
    J, Fx] = solve_TEM_KE_daily(u, v, w, T, theta, ...
    m, day, N_x, N_y, N_period, ...
    p_ini, p, phi, p_map, f_map, cosphi_map, ...
    p_start, p_end, ...
    BFJ, BFF, BFEPphi_QG, BFEPphi1, BFEPp_QG, BFEPp1, BFEPp2, ...
    BFGphi1, BFGphi2, BFGp1, BFGp2)

% Description:
% Solves Kuo-Eliassen (TEM) equation on a single day using ERA5 data.
% Includes EP flux terms, G-vector, diabatic and friction forcing, and
% outputs TEM streamfunction and intermediate diagnostics.

% Constants
a     = 6371000;
p0    = 100000;
Rd    = 287;
cp    = 1004;
g     = 9.8;
N_lev = length(p);

% Determine day index offset
if m == 1
    dd = 0; % if first year, no previous day
else
    dd = 1;
end

% Interpolate all fields to current pressure grid
u_mid = zeros(N_x, N_y, N_lev);
v_mid = zeros(N_x, N_y, N_lev);
w_mid = zeros(N_x, N_y, N_lev);
T_mid = zeros(N_x, N_y, N_lev);
theta_mid = zeros(N_x, N_y, N_lev);

for j = 1:N_x
    u_mid(j,:,:)     = compute_int(squeeze(u(j,:,:,day+dd)), p_ini, p_start, p_end, N_lev);
    v_mid(j,:,:)     = compute_int(squeeze(v(j,:,:,day+dd)), p_ini, p_start, p_end, N_lev);
    w_mid(j,:,:)     = compute_int(squeeze(w(j,:,:,day+dd)), p_ini, p_start, p_end, N_lev);
    T_mid(j,:,:)     = compute_int(squeeze(T(j,:,:,day+dd)), p_ini, p_start, p_end, N_lev);
    theta_mid(j,:,:) = compute_int(squeeze(theta(j,:,:,day+dd)), p_ini, p_start, p_end, N_lev);
end

% Zonal means and anomalies
u_bar     = squeeze(mean(u_mid,1));
v_bar     = squeeze(mean(v_mid,1));
w_bar     = squeeze(mean(w_mid,1));
theta_bar = squeeze(mean(theta_mid,1));

u_star     = u_mid - reshape(u_bar, [1 N_y N_lev]);
v_star     = v_mid - reshape(v_bar, [1 N_y N_lev]);
w_star     = w_mid - reshape(w_bar, [1 N_y N_lev]);
theta_star = theta_mid - reshape(theta_bar, [1 N_y N_lev]);

% Fluxes and derived quantities
uv_flux = u_star .* v_star;
uw_flux = u_star .* w_star;
wt_flux = w_star .* theta_star;
vt_flux = v_star .* theta_star;

vt_bar  = squeeze(mean(vt_flux,1));
dthdp   = compute_dfdp(theta_bar, p);
dudp    = compute_dfdp(u_bar, p);

vtf     = vt_bar ./ dthdp;
uv_bar  = squeeze(mean(uv_flux,1));
uw_bar  = squeeze(mean(uw_flux,1));
wt_bar  = squeeze(mean(wt_flux,1));

ducosdy = compute_dfdy(u_bar .* cosphi_map, phi);

% EP flux
Fphi_QG = -a * cosphi_map .* uv_bar;
Fphi1   = a * cosphi_map .* dudp .* vtf;
Fp_QG   = a * f_map .* cosphi_map .* vtf;
Fp1     = -ducosdy .* vtf;
Fp2     = -a * cosphi_map .* uw_bar;

% EP divergence
EPDphi_QG = 1./(a*cosphi_map) .* compute_dfdy(Fphi_QG .* cosphi_map, phi);
EPDphi1   = 1./(a*cosphi_map) .* compute_dfdy(Fphi1 .* cosphi_map, phi);
EPDp_QG   = compute_dfdp(Fp_QG, p);
EPDp1     = compute_dfdp(Fp1, p);
EPDp2     = compute_dfdp(Fp2, p);

% G-vector
Gphi1 = -theta_bar .* compute_dfdp(vtf, p);
Gphi2 = -vt_bar;
Gp1   = 1./(a*cosphi_map) .* (theta_bar .* compute_dfdy(vtf .* cosphi_map, phi));
Gp2   = -wt_bar;

% G-vector divergence
GDphi1 = 1./(a*cosphi_map) .* compute_dfdy(Gphi1 .* cosphi_map, phi);
GDphi2 = 1./(a*cosphi_map) .* compute_dfdy(Gphi2 .* cosphi_map, phi);
GDp1   = compute_dfdp(Gp1, p);
GDp2   = compute_dfdp(Gp2, p);

% Streamfunction diagnostics
[eumpsi, tempsi] = compute_TEM_ERA5(v_mid, dthdp, vt_bar, p, cosphi_map);

% Diabatic heating and friction
u3 = compute_var_3days(u, m, day, N_period, p_ini, p_start, p_end, N_lev);
v3 = compute_var_3days(v, m, day, N_period, p_ini, p_start, p_end, N_lev);
w3 = compute_var_3days(w, m, day, N_period, p_ini, p_start, p_end, N_lev);
T3 = compute_var_3days(T, m, day, N_period, p_ini, p_start, p_end, N_lev);

J     = compute_J(u3, v3, w3, T3, p, [], phi);
Fx    = compute_F(u3, u_star, v_star, w_star, u_bar, v_bar, w_bar, p, phi);
dJdy  = compute_dfdy(J, phi);
dFdp  = compute_dfdp(Fx, p);

% Forcings
HJ         = dJdy / a * Rd ./ p_map;
HF         = -f_map .* dFdp;
HEPDphi_QG = -f_map ./ (a*cosphi_map) .* compute_dfdp(EPDphi_QG, p);
HEPDphi1   = -f_map ./ (a*cosphi_map) .* compute_dfdp(EPDphi1, p);
HEPDp_QG   = -f_map ./ (a*cosphi_map) .* compute_dfdp(EPDp_QG, p);
HEPDp1     = -f_map ./ (a*cosphi_map) .* compute_dfdp(EPDp1, p);
HEPDp2     = -f_map ./ (a*cosphi_map) .* compute_dfdp(EPDp2, p);

HGDphi1 = Rd./p_map .* (p_map/p0).^(Rd/cp) .* compute_dfdy(GDphi1, phi) / a;
HGDphi2 = Rd./p_map .* (p_map/p0).^(Rd/cp) .* compute_dfdy(GDphi2, phi) / a;
HGDp1   = Rd./p_map .* (p_map/p0).^(Rd/cp) .* compute_dfdy(GDp1, phi) / a;
HGDp2   = Rd./p_map .* (p_map/p0).^(Rd/cp) .* compute_dfdy(GDp2, phi) / a;

% LU solve
FJ         = solve_PDE_LU_sparse(A, B, C, D, E, HJ, phi, p, BFJ);
FF         = solve_PDE_LU_sparse(A, B, C, D, E, HF, phi, p, BFF);
FEPphi_QG  = solve_PDE_LU_sparse(A, B, C, D, E, HEPDphi_QG, phi, p, BFEPphi_QG);
FEPphi1    = solve_PDE_LU_sparse(A, B, C, D, E, HEPDphi1, phi, p, BFEPphi1);
FEPp_QG    = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp_QG, phi, p, BFEPp_QG);
FEPp1      = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp1, phi, p, BFEPp1);
FEPp2      = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp2, phi, p, BFEPp2);
FGphi1     = solve_PDE_LU_sparse(A, B, C, D, E, HGDphi1, phi, p, BFGphi1);
FGphi2     = solve_PDE_LU_sparse(A, B, C, D, E, HGDphi2, phi, p, BFGphi2);
FGp1       = solve_PDE_LU_sparse(A, B, C, D, E, HGDp1, phi, p, BFGp1);
FGp2       = solve_PDE_LU_sparse(A, B, C, D, E, HGDp2, phi, p, BFGp2);

% Final sum of components
FEPphi = FEPphi_QG + FEPphi1;
FEPp   = FEPp_QG + FEPp1 + FEPp2;
FEP    = FEPphi + FEPp;

FGphi  = FGphi1 + FGphi2;
FGp    = FGp1 + FGp2;
FG     = FGphi + FGp;

FTEM   = FJ + FF + FEP + FG;

end
