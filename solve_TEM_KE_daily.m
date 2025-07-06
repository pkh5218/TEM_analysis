function [FTEM, FJ, FF, FEPphi_QG, FEPphi1, FEPp_QG, FEPp1, FEPp2,...
    FGphi1, FGphi2, FGp1, FGp2, tempsi, eumpsi] = solve_TEM_KE_daily(u, v, w, T, theta,...
    m, day,...
    N_x, N_y, N_period, ...
    p_ini, p, phi, p_map, f_map, cosphi_map, tanphi_map,...
    p_start, p_end,...
    BFJ, BFF, BFEPphi_QG, BFEPphi1, ...
    BFEPp_QG, BFEPp1, BFEPp2,...
    BFGphi1, BFGphi2, BFGp1, BFGp2)

path        = '/jsbach/s0/mup65/ERA5/';
lon         = ncread([path,'UCOMP/u_daily_mean.1979.nc'],'longitude');
lon         = double (lon);
lat         = ncread([path,'UCOMP/u_daily_mean.1979.nc'],'latitude');
lat         = double (lat);
g           = 9.8;
a           = 6371000;
p0          = 1000*100;
Rd          = 287;
cp          = 1004;
N_new_levels = length(p);

if m == 1
   dd = 0;
else
   dd = 1;
end

for j = 1:N_x
    u_mid(j,:,:)    = compute_int(squeeze(u(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    v_mid(j,:,:)    = compute_int(squeeze(v(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    w_mid(j,:,:)    = compute_int(squeeze(w(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    T_mid(j,:,:)    = compute_int(squeeze(T(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
    theta_mid(j,:,:)= compute_int(squeeze(theta(j,:,:,day+dd)),p_ini, p_start, p_end, N_new_levels);
end

% Computation of zonal average, which has a size of 1*lat*lev
u_bar       = mean(u_mid,1);
v_bar       = mean(v_mid,1);
w_bar       = mean(w_mid,1);
T_bar       = mean(T_mid,1);
theta_bar   = mean(theta_mid,1);

% Computation of anomaly, which has a size of lon*lat*lev
u_star      = u_mid;
v_star      = v_mid;
w_star      = w_mid;
T_star      = T_mid;
theta_star  = theta_mid;

for j = 1:N_x
    u_star(j,:,:) = u_mid(j,:,:) - u_bar;
    v_star(j,:,:) = v_mid(j,:,:) - v_bar;
    w_star(j,:,:) = w_mid(j,:,:) - w_bar;
    T_star(j,:,:) = T_mid(j,:,:) - T_bar;
    theta_star(j,:,:) = theta_mid(j,:,:) - theta_bar;
end
    
% Change the array size to lat*lev 
u_bar       = squeeze(mean(u_mid,1));
v_bar       = squeeze(mean(v_mid,1));
w_bar       = squeeze(mean(w_mid,1));
T_bar       = squeeze(mean(T_mid,1));
theta_bar   = squeeze(mean(theta_mid,1));

% Computation of flux terms
uv_flux     = u_star.*v_star;
vt_flux     = v_star.*theta_star;
uw_flux     = u_star.*w_star;
wt_flux     = w_star.*theta_star;

vt_bar      = squeeze(mean(vt_flux,[1,4]));
dthdp       = compute_dfdp(theta_bar,p);
dudp        = compute_dfdp(u_bar,p);

% Zonal mean for flux terms
vt_flux_bar = squeeze(mean(vt_flux,[1,4], "omitnan"));
uv_flux_bar = squeeze(mean(uv_flux,[1,4], "omitnan"));
uw_flux_bar = squeeze(mean(uw_flux,[1,4]));
wt_flux_bar = squeeze(mean(wt_flux,[1,4]));
vtf         = vt_flux_bar./dthdp;
ducosdy     = compute_dfdy(u_bar.*cosphi_map,phi);

% Computation of EP flux
Fphi_QG = -a.*cosphi_map.*uv_flux_bar;
Fphi1   = a.*cosphi_map.*dudp.*vtf;

Fp_QG   = a*f_map.*cosphi_map.*vtf;
Fp1     = -ducosdy.*vtf;
Fp2     = -a.*cosphi_map.*uw_flux_bar;
    
% Computation of EP flux Divergence
EPDphi_QG = 1./(a.*cosphi_map).*compute_dfdy(Fphi_QG.*cosphi_map,phi);
EPDphi1   = 1./(a.*cosphi_map).*compute_dfdy(Fphi1.*cosphi_map,phi);
EPDp_QG   = compute_dfdp(Fp_QG,p);
EPDp1     = compute_dfdp(Fp1,p);
EPDp2     = compute_dfdp(Fp2,p);
    
% Computation of G-vector from Thermodynamic Eqn.
Gphi1   = -theta_bar.*compute_dfdp(vtf,p);
Gphi2   = -vt_flux_bar;
Gp1     = 1./(a.*cosphi_map).*(theta_bar.*compute_dfdy(vtf.*cosphi_map,phi));
Gp2     = -wt_flux_bar;

% Computation of G-vector Divergence
GDphi1  = 1./(a.*cosphi_map).*compute_dfdy(Gphi1.*cosphi_map,phi);
GDphi2  = 1./(a.*cosphi_map).*compute_dfdy(Gphi2.*cosphi_map,phi);
GDp1    = compute_dfdp(Gp1,p);
GDp2    = compute_dfdp(Gp2,p);
    
%% Computation of the TEM streamfunction in the stratosphere
% The real TEM streamfunction
[eumpsi, tempsi]  = compute_TEM_ERA5(v_mid,dthdp,vt_bar,p,cosphi_map);
% The 12-hourly data from 3 points
u_3days = compute_var_3days(u,m,day,N_period,...
    p_ini, p_start, p_end, N_new_levels);
v_3days = compute_var_3days(v,m,day,N_period,...
    p_ini, p_start, p_end, N_new_levels);
w_3days = compute_var_3days(w,m,day,N_period,...
    p_ini, p_start, p_end, N_new_levels);
T_3days = compute_var_3days(T,m,day,N_period,...
    p_ini, p_start, p_end, N_new_levels);

% Computation of diabatic heating and friction
J         = compute_J(u_3days, v_3days, w_3days, T_3days, p, lon, lat);
Fx        = compute_F(u_3days, u_star, v_star, w_star, u_bar, v_bar, w_bar, p, lat);
dJdy      = compute_dfdy(J,phi);
dFdp      = compute_dfdp(Fx,p);

% Compute each Forcing term, H
HJ          = dJdy/a*Rd./p_map;
HF          = -f_map.*dFdp;
HEPDphi_QG  = -f_map./(a.*cosphi_map).*compute_dfdp(EPDphi_QG,p);
HEPDphi1    = -f_map./(a.*cosphi_map).*compute_dfdp(EPDphi1,p);
HEPDp_QG    = -f_map./(a.*cosphi_map).*compute_dfdp(EPDp_QG,p);
HEPDp1      = -f_map./(a.*cosphi_map).*compute_dfdp(EPDp1,p);
HEPDp2      = -f_map./(a.*cosphi_map).*compute_dfdp(EPDp2,p);

HGDphi1     = Rd./p_map.*(p_map./p0).^(Rd/cp).*compute_dfdy(GDphi1,phi)./a;
HGDphi2     = Rd./p_map.*(p_map./p0).^(Rd/cp).*compute_dfdy(GDphi2,phi)./a;
HGDp1       = Rd./p_map.*(p_map./p0).^(Rd/cp).*compute_dfdy(GDp1,phi)./a;
HGDp2       = Rd./p_map.*(p_map./p0).^(Rd/cp).*compute_dfdy(GDp2,phi)./a;

%% Caculation of L terms
% Computation of stability and differential terms
S2          = compute_S2(T_bar,p);
dS2dy       = compute_dfdy(S2,phi);
dTdy        = compute_dfdy(T_bar,phi);
dTdy2       = compute_dfdy(dTdy,phi);
dudp        = compute_dfdp(u_bar,p);
dudp2       = compute_dfdp(dudp,p);
dudy        = compute_dfdy(u_bar,phi);
dudydp      = compute_dfdp(dudy,p); 

Lf          = f_map.^2*g/2/pi/a./cosphi_map;
LS2_1       = g/2/pi/a^3./cosphi_map.*S2;
LS2_2       = g/2/pi/a^3./cosphi_map.*(dS2dy+tanphi_map.*S2);
LT_1        = g/2/pi/a^3./cosphi_map*Rd./p_map.*dTdy;
LT_2        = g/2/pi/a^3./cosphi_map*Rd./p_map.*(dTdy2+tanphi_map.*dTdy);
LU_1        = f_map*g/2/pi/a^2./cosphi_map.*dudp;
LU_2        = f_map*g/2/pi/a^2./cosphi_map.*dudp2;
LU_3        = -f_map*g/2/pi/a^2./cosphi_map.*(dudydp-tanphi_map.*dudp);
LU_4        = -f_map*g/2/pi/a^2./cosphi_map.*(dudy-u_bar.*tanphi_map);

% Given matrices A, B, C, D, E, and H (size lat*lev)
A = Lf + LU_4;
B = (LT_1 + LU_1)/2;
C = LS2_1;
D = LT_2 + LU_3;
E = LS2_2 + LU_2;

if BFEPp_QG == 1
    BFEPp_QG = zeros(N_y,1);
    dw = 1/a./cosphi_map.*compute_dfdy(vtf.*cosphi_map,phi);
    for j = 1:N_y-1
        BFEPp_QG(j+1) = BFEPp_QG(j)-(2*pi*a^2*cosphi_map(j,length(p))/g)*(dw(j,length(p)))*(phi(j+1)-phi(j));
    end
end

% Solve KE equation with LU decomposition
FJ  = solve_PDE_LU_sparse(A, B, C, D, E, HJ, phi, p, BFJ);
FF  = solve_PDE_LU_sparse(A, B, C, D, E, HF, phi, p, BFF);
FEPphi_QG = solve_PDE_LU_sparse(A, B, C, D, E, HEPDphi_QG, phi, p, BFEPphi_QG);
FEPphi1 = solve_PDE_LU_sparse(A, B, C, D, E, HEPDphi1, phi, p, BFEPphi1);
FEPp_QG = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp_QG, phi, p, BFEPp_QG);
FEPp1 = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp1, phi, p, BFEPp1);
FEPp2 = solve_PDE_LU_sparse(A, B, C, D, E, HEPDp2, phi, p, BFEPp2);
FGphi1  = solve_PDE_LU_sparse(A, B, C, D, E, HGDphi1, phi, p, BFGphi1);
FGphi2  = solve_PDE_LU_sparse(A, B, C, D, E, HGDphi2, phi, p, BFGphi2);
FGp1  = solve_PDE_LU_sparse(A, B, C, D, E, HGDp1, phi, p, BFGp1);
FGp2  = solve_PDE_LU_sparse(A, B, C, D, E, HGDp2, phi, p, BFGp2);
    
FEPphi  = FEPphi_QG + FEPphi1;
FEPp    = FEPp_QG + FEPp1 + FEPp2;
FEP     = FEPphi + FEPp;

FGphi   = FGphi1 + FGphi2;
FGp     = FGp1 + FGp2;
FG      = FGphi + FGp;

FTEM    = FJ+FF+FEP+FG;