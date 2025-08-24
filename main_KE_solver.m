%% MAIN SCRIPT: TEM-KE daily
% Author : Pin-Chun Huang
% Date   : 2025-08-23

clear; clc; close all

%% Constants and settings
Omega   = 7.292e-5;
g       = 9.8;
a       = 6371000;
p0      = 1000*100;    % Pa
Rd      = 287;
cp      = 1004;

start_year  = 1979;
end_year    = 2020;
years       = start_year:end_year;
Nyears      = numel(years);

% Uniform pressure grid (Pa): 0~1000 hPa with 41 levels (= 25hPa) or 81 levels (= 12.5 hPa)
Nlev    = 81;
p_start = 0;
p_end   = 1e5;
p       = linspace(p_start, p_end, Nlev);   % increasing top -> bottom
p(1)    = 100;                              % avoid p=0 singularity (set 1 hPa)

%% Load lon/lat and original pressure levels (Pa)
path   = '/jsbach/s0/mup65/ERA5/';
path_U = 'UCOMP/u_daily_mean.';
path_V = 'VCOMP/v_daily_mean.';
path_W = 'OMEGA/w_daily_mean.';
path_T = 'TEMP/T_daily_mean.';
outdir = '/jsbach/s0/pkh5218/TEM_new/';

% Use one file just to read coordinates
lon  = double(ncread([path, path_U, '1979.nc'],'longitude'));
lat  = double(ncread([path, path_U, '1979.nc'],'latitude'));
lev0 = double(ncread([path, path_U, '1979.nc'],'level'))*100; % Pa

phi    = deg2rad(lat);      % rad
sinphi = sin(phi);
cosphi = cos(phi);
tanphi = tan(phi);

% --- Polar guard for cosφ to avoid singularity ---
dphi         = phi(2) - phi(1);
cosphi(1)    = sin(abs(dphi)/2);
cosphi(end)  = sin(abs(dphi)/2);

% recompute tanφ from guarded cosφ
tanphi(1) = sinphi(1) ./ cosphi(1);
tanphi(end) = sinphi(end) ./ cosphi(end);

f        = 2*Omega.*sinphi;     % Coriolis Parameter

Nx       = numel(lon);
Ny       = numel(lat);
Nlev0    = numel(lev0);         % Number of levels for pre-interpolated data

% Pre-compute maps (lat x Nlev)
[~, f_map]       = meshgrid(p, f);
[p_map, cos_map] = meshgrid(p, cosphi);
[~, tan_map]     = meshgrid(p, tanphi);

%% Loop over years
for m = 1:Nyears
    tic
    ycur = years(m);
    fprintf('Currently Running - year: %d\n', ycur);

    % Days in current year (leap year rule)   
    if (mod(ycur,400)==0) || (mod(ycur,4)==0 && mod(ycur,100)~=0)
        Ndays = 366;
    else
        Ndays = 365;
    end

    % Allocate yearly arrays (lat x lev x day)
    FTEM_year        = zeros(Ny, Nlev, Ndays);
    FEUM_year        = zeros(Ny, Nlev, Ndays);
    tempsi_year      = zeros(Ny, Nlev, Ndays);
    eumpsi_year      = zeros(Ny, Nlev, Ndays);

    FJ_year          = zeros(Ny, Nlev, Ndays);
    FF_year          = zeros(Ny, Nlev, Ndays);
    FEPDphi_QG_year  = zeros(Ny, Nlev, Ndays);
    FEPDphi1_year    = zeros(Ny, Nlev, Ndays);
    FEPDp_QG_year    = zeros(Ny, Nlev, Ndays);
    FEPDp1_year      = zeros(Ny, Nlev, Ndays);
    FEPDp2_year      = zeros(Ny, Nlev, Ndays);
    FGDphi1_year     = zeros(Ny, Nlev, Ndays);
    FGDphi2_year     = zeros(Ny, Nlev, Ndays);
    FGDp1_year       = zeros(Ny, Nlev, Ndays);
    FGDp2_year       = zeros(Ny, Nlev, Ndays);

    HJ_year          = zeros(Ny, Nlev, Ndays);
    HX_year          = zeros(Ny, Nlev, Ndays);
    HEPDphi_QG_year  = zeros(Ny, Nlev, Ndays);
    HEPDphi1_year    = zeros(Ny, Nlev, Ndays);
    HEPDp_QG_year    = zeros(Ny, Nlev, Ndays);
    HEPDp1_year      = zeros(Ny, Nlev, Ndays);
    HEPDp2_year      = zeros(Ny, Nlev, Ndays);
    HGDphi1_year     = zeros(Ny, Nlev, Ndays);
    HGDphi2_year     = zeros(Ny, Nlev, Ndays);
    HGDp1_year       = zeros(Ny, Nlev, Ndays);
    HGDp2_year       = zeros(Ny, Nlev, Ndays);

    EPDphi_QG_year   = zeros(Ny, Nlev, Ndays);
    EPDphi1_year     = zeros(Ny, Nlev, Ndays);
    EPDp_QG_year     = zeros(Ny, Nlev, Ndays);
    EPDp1_year       = zeros(Ny, Nlev, Ndays);
    EPDp2_year       = zeros(Ny, Nlev, Ndays);

    EPphi_QG_year    = zeros(Ny, Nlev, Ndays);
    EPphi1_year      = zeros(Ny, Nlev, Ndays);
    EPp_QG_year      = zeros(Ny, Nlev, Ndays);
    EPp1_year        = zeros(Ny, Nlev, Ndays);
    EPp2_year        = zeros(Ny, Nlev, Ndays);

    Gphi1_year       = zeros(Ny, Nlev, Ndays);
    Gphi2_year       = zeros(Ny, Nlev, Ndays);
    Gp1_year         = zeros(Ny, Nlev, Ndays);
    Gp2_year         = zeros(Ny, Nlev, Ndays);

    GDphi1_year      = zeros(Ny, Nlev, Ndays);
    GDphi2_year      = zeros(Ny, Nlev, Ndays);
    GDp1_year        = zeros(Ny, Nlev, Ndays);
    GDp2_year        = zeros(Ny, Nlev, Ndays);

    J_year           = zeros(Ny, Nlev, Ndays);
    X_year           = zeros(Ny, Nlev, Ndays);
    EHF_year         = zeros(Ny, Nlev, Ndays);
    EMF_year         = zeros(Ny, Nlev, Ndays);
    VEHF_year        = zeros(Ny, Nlev, Ndays);
    VEMF_year        = zeros(Ny, Nlev, Ndays);

    % Load daily fields for current year
    if m == 1
    u0 = ncread([path, path_U, num2str(ycur), '.nc'], 'u');
    v0 = ncread([path, path_V, num2str(ycur), '.nc'], 'v');
    w0 = ncread([path, path_W, num2str(ycur), '.nc'], 'w');
    T0 = ncread([path, path_T, num2str(ycur), '.nc'], 't');
    u_b = u0(:,:,:,1); v_b = v0(:,:,:,1); w_b = w0(:,:,:,1); T_b = T0(:,:,:,1);
    u   = u0; v = v0; w = w0; T = T0;
    else
    ub = ncread([path, path_U, num2str(ycur-1), '.nc'], 'u'); u_b = ub(:,:,:,end);
    vb = ncread([path, path_V, num2str(ycur-1), '.nc'], 'v'); v_b = vb(:,:,:,end);
    wb = ncread([path, path_W, num2str(ycur-1), '.nc'], 'w'); w_b = wb(:,:,:,end);
    Tb = ncread([path, path_T, num2str(ycur-1), '.nc'], 't'); T_b = Tb(:,:,:,end);
    u = ncread([path, path_U, num2str(ycur), '.nc'], 'u');
    v = ncread([path, path_V, num2str(ycur), '.nc'], 'v');
    w = ncread([path, path_W, num2str(ycur), '.nc'], 'w');
    T = ncread([path, path_T, num2str(ycur), '.nc'], 't');
    end

    if m == Nyears
    uf = ncread([path, path_U, num2str(ycur), '.nc'], 'u'); u_f = uf(:,:,:,end);
    vf = ncread([path, path_V, num2str(ycur), '.nc'], 'v'); v_f = vf(:,:,:,end);
    wf = ncread([path, path_W, num2str(ycur), '.nc'], 'w'); w_f = wf(:,:,:,end);
    Tf = ncread([path, path_T, num2str(ycur), '.nc'], 't'); T_f = Tf(:,:,:,end);
    else
    uf = ncread([path, path_U, num2str(ycur+1), '.nc'], 'u'); u_f = uf(:,:,:,1);
    vf = ncread([path, path_V, num2str(ycur+1), '.nc'], 'v'); v_f = vf(:,:,:,1);
    wf = ncread([path, path_W, num2str(ycur+1), '.nc'], 'w'); w_f = wf(:,:,:,1);
    Tf = ncread([path, path_T, num2str(ycur+1), '.nc'], 't'); T_f = Tf(:,:,:,1);
    end

    % Concatenate halos → dim4: NDays+2 (back, current, forward)
    u = cat(4, u_b, u, u_f);
    v = cat(4, v_b, v, v_f);
    w = cat(4, w_b, w, w_f);
    T = cat(4, T_b, T, T_f);

    % Initialize
    u_prev = zeros(Nx,Ny,Nlev); v_prev = u_prev; w_prev = u_prev; T_prev = u_prev; th_prev = u_prev;
    u_curr = u_prev;            v_curr = u_prev; w_curr = u_prev; T_curr = u_prev; th_curr = u_prev;
    u_next = u_prev;            v_next = u_prev; w_next = u_prev; T_next = u_prev; th_next = u_next;
 %% Daily loop
    for day = 1:Ndays
    disp(['Day: ', num2str(day),'/',num2str(Ndays)]);
    % ----- 3-Day Sliding-Window Interpolation (prev/curr/next) -----
    % Interpolation
    K = day + 1;
    if day == 1
    for j = 1:Nx
        u_prev(j,:,:)  = compute_int(squeeze(u(j,:,:,K-1)), lev0, p_start, p_end, Nlev);
        v_prev(j,:,:)  = compute_int(squeeze(v(j,:,:,K-1)), lev0, p_start, p_end, Nlev);
        w_prev(j,:,:)  = compute_int(squeeze(w(j,:,:,K-1)), lev0, p_start, p_end, Nlev);
        T_prev(j,:,:)  = compute_int(squeeze(T(j,:,:,K-1)), lev0, p_start, p_end, Nlev);

        u_curr(j,:,:)  = compute_int(squeeze(u(j,:,:,K)),   lev0, p_start, p_end, Nlev);
        v_curr(j,:,:)  = compute_int(squeeze(v(j,:,:,K)),   lev0, p_start, p_end, Nlev);
        w_curr(j,:,:)  = compute_int(squeeze(w(j,:,:,K)),   lev0, p_start, p_end, Nlev);
        T_curr(j,:,:)  = compute_int(squeeze(T(j,:,:,K)),   lev0, p_start, p_end, Nlev);

        u_next(j,:,:)  = compute_int(squeeze(u(j,:,:,K+1)), lev0, p_start, p_end, Nlev);
        v_next(j,:,:)  = compute_int(squeeze(v(j,:,:,K+1)), lev0, p_start, p_end, Nlev);
        w_next(j,:,:)  = compute_int(squeeze(w(j,:,:,K+1)), lev0, p_start, p_end, Nlev);
        T_next(j,:,:)  = compute_int(squeeze(T(j,:,:,K+1)), lev0, p_start, p_end, Nlev);
    end
    end
    

    % Compute Potential Temperature
    th_prev = T_prev; th_curr = T_curr; th_next = T_next; 
    for k = 1:Nlev 
        th_prev(:,:,k,:) = T_prev(:,:,k,:)*((p0/p(k))^(Rd/cp)); 
        th_curr(:,:,k,:) = T_curr(:,:,k,:)*((p0/p(k))^(Rd/cp)); 
        th_next(:,:,k,:) = T_next(:,:,k,:)*((p0/p(k))^(Rd/cp)); 
    end

    % (1) Build 3-day stacks and current day's fields
    u_3day = cat(4, u_prev, u_curr, u_next);
    v_3day = cat(4, v_prev, v_curr, v_next);
    w_3day = cat(4, w_prev, w_curr, w_next);
    T_3day = cat(4, T_prev, T_curr, T_next);

    % (2) Zonal means (lat x lev)
    u_bar = squeeze(mean(u_curr,1));
    v_bar = squeeze(mean(v_curr,1));
    w_bar = squeeze(mean(w_curr,1));
    T_bar = squeeze(mean(T_curr,1));
    th_bar= squeeze(mean(th_curr,1));

    % (3) Eddy components
    u_prime  = u_curr - reshape(u_bar,[1,Ny,Nlev]);
    v_prime  = v_curr - reshape(v_bar,[1,Ny,Nlev]);
    w_prime  = w_curr - reshape(w_bar,[1,Ny,Nlev]);
    th_prime = th_curr - reshape(th_bar,[1,Ny,Nlev]);

    % (4) Zonal Mean Eddy fluxes
    uv_flux = u_prime .* v_prime;
    vt_flux = v_prime .* th_prime;
    uw_flux = u_prime .* w_prime;
    wt_flux = w_prime .* th_prime;

    EMF  = squeeze(mean(uv_flux,1));  % meridional eddy momentum flux
    EHF  = squeeze(mean(vt_flux,1));  % meridional eddy heat flux
    VEMF = squeeze(mean(uw_flux,1));  % vertical eddy momentum flux
    VEHF = squeeze(mean(wt_flux,1));  % vertical eddy heat flux

    % (5) Derivatives on uniform grids
    dthdp   = compute_dfdp(th_bar, p);          % Ny x Nlev
    dudp    = compute_dfdp(u_bar,  p);
    ducosdy = compute_dfdy(u_bar .* cos_map, phi);

    % (6) EP flux components
    EHF_dthdp = EHF ./ dthdp;
    EPphi_QG  = -a .* cos_map .* EMF;
    EPphi1    =  a .* cos_map .* dudp .* EHF_dthdp;
    EPp_QG    =  a .* f_map   .* cos_map .* EHF_dthdp;
    EPp1      = -ducosdy .* EHF_dthdp;
    EPp2      = -a .* cos_map .* VEMF;

    % (7) EP flux divergence
    EPDphi_QG = 1./(a.*cos_map) .* compute_dfdy(EPphi_QG .* cos_map, phi);
    EPDphi1   = 1./(a.*cos_map) .* compute_dfdy(EPphi1   .* cos_map, phi);
    EPDp_QG   = compute_dfdp(EPp_QG, p);
    EPDp1     = compute_dfdp(EPp1,   p);
    EPDp2     = compute_dfdp(EPp2,   p);

    % (8) G-vector terms and divergences
    Gphi1  = -th_bar .* compute_dfdp(EHF_dthdp, p);
    Gphi2  = -EHF;
    Gp1    = 1./(a.*cos_map) .* ( th_bar .* compute_dfdy(EHF_dthdp .* cos_map, phi) );
    Gp2    = -VEHF;

    GDphi1 = 1./(a.*cos_map) .* compute_dfdy(Gphi1 .* cos_map, phi);
    GDphi2 = 1./(a.*cos_map) .* compute_dfdy(Gphi2 .* cos_map, phi);
    GDp1   = compute_dfdp(Gp1, p);
    GDp2   = compute_dfdp(Gp2, p);

    % (9) TEM/EUM streamfunction
    [eumpsi, tempsi] = compute_TEM_ERA5(v_bar, EHF_dthdp, p, cos_map);

    % (10) Diabatic heating & friction (from 3-day fields)
    [J, ~] = compute_J(u_3day, v_3day, w_3day, T_3day, p, lon, lat, m, day, Nyears, Ndays);
    X = compute_X(u_3day, u_prime, v_prime, w_prime, u_bar, v_bar, w_bar, p, lat);

    dJdy = compute_dfdy(J,  phi);
    dXdp = compute_dfdp(X, p);

    % (11) Forcings H
    HJ         = dJdy/a * Rd ./ p_map;
    HX         = -f_map .* dXdp;
    HEPDphi_QG = -f_map./(a.*cos_map) .* compute_dfdp(EPDphi_QG, p);
    HEPDphi1   = -f_map./(a.*cos_map) .* compute_dfdp(EPDphi1,   p);
    HEPDp_QG   = -f_map./(a.*cos_map) .* compute_dfdp(EPDp_QG,   p);
    HEPDp1     = -f_map./(a.*cos_map) .* compute_dfdp(EPDp1,     p);
    HEPDp2     = -f_map./(a.*cos_map) .* compute_dfdp(EPDp2,     p);

    HGDphi1    = Rd./p_map.*(p_map./p0).^(Rd/cp) .* compute_dfdy(GDphi1, phi)./a;
    HGDphi2    = Rd./p_map.*(p_map./p0).^(Rd/cp) .* compute_dfdy(GDphi2, phi)./a;
    HGDp1      = Rd./p_map.*(p_map./p0).^(Rd/cp) .* compute_dfdy(GDp1,   phi)./a;
    HGDp2      = Rd./p_map.*(p_map./p0).^(Rd/cp) .* compute_dfdy(GDp2,   phi)./a;

    % (12) L-operator coefficients
    S2     = compute_S2(T_bar, p);
    dS2dy  = compute_dfdy(S2, phi);
    dTdy   = compute_dfdy(T_bar, phi);
    dTdy2  = compute_dfdy(dTdy,  phi);
    dudp   = compute_dfdp(u_bar, p);
    dudp2  = compute_dfdp(dudp,  p);
    dudy   = compute_dfdy(u_bar, phi);
    dudydp = compute_dfdp(dudy,  p);

    Lf   = f_map.^2 * g/(2*pi*a) ./ cos_map;
    LS2_1= g/(2*pi*a^3) ./ cos_map .* S2;
    LS2_2= g/(2*pi*a^3) ./ cos_map .* (dS2dy + tan_map .* S2);
    LT_1 = g/(2*pi*a^3) ./ cos_map .* Rd./p_map .* dTdy;
    LT_2 = g/(2*pi*a^3) ./ cos_map .* Rd./p_map .* (dTdy2 + tan_map .* dTdy);
    LU_1 = f_map * g/(2*pi*a^2) ./ cos_map .* dudp;
    LU_2 = f_map * g/(2*pi*a^2) ./ cos_map .* dudp2;
    LU_3 = -f_map * g/(2*pi*a^2) ./ cos_map .* (dudydp - tan_map .* dudp);
    LU_4 = -f_map * g/(2*pi*a^2) ./ cos_map .* (dudy - u_bar .* tan_map);

    % Combine A ~ E
    A = Lf + LU_4;
    B = (LT_1 + LU_1)/2;
    C = LS2_1;
    D = LT_2 + LU_3;
    E = LS2_2 + LU_2;

    % --- Boundary Condition --- 
    B0 = zeros(Ny,1);
    % bottom boundary for the FEPp_QG
    omega_star = 1./(a.*cos_map) .* compute_dfdy(EHF_dthdp .* cos_map, phi);
    omega_bot  = omega_star(:, end);                                        
    integrand  = (2*pi*a^2/g) .* cosphi .* omega_bot;
    BFEPp_QG   = -cumtrapz(phi, integrand);                                   % (lat x 1)

    % Solve KE equation with LU decomposition (uniform grid)
    FJ          = solve_PDE_LU_uniform(A,B,C,D,E, HJ,         phi, p, B0);
    FX          = solve_PDE_LU_uniform(A,B,C,D,E, HX,         phi, p, B0);
    FEPDphi_QG  = solve_PDE_LU_uniform(A,B,C,D,E, HEPDphi_QG, phi, p, B0);
    FEPDphi1    = solve_PDE_LU_uniform(A,B,C,D,E, HEPDphi1,   phi, p, B0);
    FEPDp_QG    = solve_PDE_LU_uniform(A,B,C,D,E, HEPDp_QG,   phi, p, BFEPp_QG);
    FEPDp1      = solve_PDE_LU_uniform(A,B,C,D,E, HEPDp1,     phi, p, B0);
    FEPDp2      = solve_PDE_LU_uniform(A,B,C,D,E, HEPDp2,     phi, p, B0);
    FGDphi1     = solve_PDE_LU_uniform(A,B,C,D,E, HGDphi1,    phi, p, B0);
    FGDphi2     = solve_PDE_LU_uniform(A,B,C,D,E, HGDphi2,    phi, p, B0);
    FGDp1       = solve_PDE_LU_uniform(A,B,C,D,E, HGDp1,      phi, p, B0);
    FGDp2       = solve_PDE_LU_uniform(A,B,C,D,E, HGDp2,      phi, p, B0);

    % Sums
    FEPDphi     = FEPDphi_QG + FEPDphi1;
    FEPDp       = FEPDp_QG + FEPDp1 + FEPDp2;
    FEPD        = FEPDphi + FEPDp;
    FGDphi      = FGDphi1 + FGDphi2;
    FGDp        = FGDp1 + FGDp2;
    FGD         = FGDphi + FGDp;
    FTEM        = FJ + FF + FEPD + FGD;
    FEUM        = FJ + FF + FEPDphi_QG + FEPDp2 + FGDphi2 + FGDp2;
    
    % Store
    FJ_year(:,:,day)            = FJ;
    FX_year(:,:,day)            = FX;
    FEPDphi_QG_year(:,:,day)    = FEPDphi_QG;
    FEPDphi1_year(:,:,day)      = FEPDphi1;
    FEPDp_QG_year(:,:,day)      = FEPDp_QG;
    FEPDp1_year(:,:,day)        = FEPDp1;
    FEPDp2_year(:,:,day)        = FEPDp2;
    FGDphi1_year(:,:,day)       = FGDphi1;
    FGDphi2_year(:,:,day)       = FGDphi2;
    FGDp1_year(:,:,day)         = FGDp1;
    FGDp2_year(:,:,day)         = FGDp2;
    FTEM_year(:,:,day)          = FTEM;
    FEUM_year(:,:,day)          = FEUM;
    tempsi_year(:,:,day)        = tempsi;
    eumpsi_year(:,:,day)        = eumpsi;

    HJ_year(:,:,day)            = HJ;
    HX_year(:,:,day)            = HX;
    HEPDphi_QG_year(:,:,day)    = HEPDphi_QG;
    HEPDphi1_year(:,:,day)      = HEPDphi1;
    HEPDp_QG_year(:,:,day)      = HEPDp_QG;
    HEPDp1_year(:,:,day)        = HEPDp1;
    HEPDp2_year(:,:,day)        = HEPDp2;
    HGDphi1_year(:,:,day)       = HGDphi1;
    HGDphi2_year(:,:,day)       = HGDphi2;
    HGDp1_year(:,:,day)         = HGDp1;
    HGDp2_year(:,:,day)         = HGDp2;

    J_year(:,:,day)             = J;
    X_year(:,:,day)             = X;
    EPDphi_QG_year(:,:,day)     = EPDphi_QG;
    EPDphi1_year(:,:,day)       = EPDphi1;
    EPDp_QG_year(:,:,day)       = EPDp_QG;
    EPDp1_year(:,:,day)         = EPDp1;
    EPDp2_year(:,:,day)         = EPDp2;
    GDphi1_year(:,:,day)        = GDphi1;
    GDphi2_year(:,:,day)        = GDphi2;
    GDp1_year(:,:,day)          = GDp1;
    GDp2_year(:,:,day)          = GDp2;

    EPphi_QG_year(:,:,day)      = EPphi_QG;
    EPphi1_year(:,:,day)        = EPphi1;
    EPp_QG_year(:,:,day)        = EPp_QG;
    EPp1_year(:,:,day)          = EPp1;
    EPp2_year(:,:,day)          = EPp2;
    Gphi1_year(:,:,day)         = Gphi1;
    Gphi2_year(:,:,day)         = Gphi2;
    Gp1_year(:,:,day)           = Gp1;
    Gp2_year(:,:,day)           = Gp2;

    EMF_year(:,:,day)           = EMF;
    EHF_year(:,:,day)           = EHF;
    VEMF_year(:,:,day)          = VEMF;
    VEHF_year(:,:,day)          = VEHF;

    % (13) Slide window right by one day: compute new "next"
    if day < Ndays
    K = K + 1; 
    u_prev = u_curr;  v_prev = v_curr;  w_prev = w_curr;  T_prev = T_curr;
    u_curr = u_next;  v_curr = v_next;  w_curr = w_next;  T_curr = T_next;

    for j = 1:Nx
        u_next(j,:,:)  = compute_int(squeeze(u(j,:,:,K+1)),  lev0, p_start, p_end, Nlev);
        v_next(j,:,:)  = compute_int(squeeze(v(j,:,:,K+1)),  lev0, p_start, p_end, Nlev);
        w_next(j,:,:)  = compute_int(squeeze(w(j,:,:,K+1)),  lev0, p_start, p_end, Nlev);
        T_next(j,:,:)  = compute_int(squeeze(T(j,:,:,K+1)),  lev0, p_start, p_end, Nlev);
    end

    end
    end
    %% ====== WRITE NETCDF (once per year) ======
    lev_hPa = p/100;
    
    % 1) TEM_yyyy.nc
    out1 = [outdir, 'TEM_', num2str(ycur), '.nc'];
    write_group(out1, lat, lev_hPa, Ndays, { ...
    'FTEM',        FTEM_year; ...
    'FEUM',        FEUM_year; ...
    'FJ',          FJ_year; ...
    'FX',          FX_year; ...
    'FEPDphi_QG',  FEPDphi_QG_year; ...
    'FEPDphi1',    FEPDphi1_year; ...
    'FEPDp_QG',    FEPDp_QG_year; ...
    'FEPDp1',      FEPDp1_year; ...
    'FEPDp2',      FEPDp2_year; ...
    'FGDphi1',     FGDphi1_year; ...
    'FGDphi2',     FGDphi2_year; ...
    'FGDp1',       FGDp1_year; ...
    'FGDp2',       FGDp2_year; ...
    'tempsi',      tempsi_year; ...
    'eumpsi',      eumpsi_year });

    % 2) RHS_yyyy.nc
    out2 = [outdir, 'RHS_', num2str(ycur), '.nc'];
    write_group(out2, lat, lev_hPa, Ndays, { ...
    'HJ',          HJ_year; ...
    'HX',          HX_year; ...
    'HEPDphi_QG',  HEPDphi_QG_year; ...
    'HEPDphi1',    HEPDphi1_year; ...
    'HEPDp_QG',    HEPDp_QG_year; ...
    'HEPDp1',      HEPDp1_year; ...
    'HEPDp2',      HEPDp2_year; ...
    'HGDphi1',     HGDphi1_year; ...
    'HGDphi2',     HGDphi2_year; ...
    'HGDp1',       HGDp1_year; ...
    'HGDp2',       HGDp2_year });

    % 3) Forcing_yyyy.nc (Including EPD, GD, Diabatic Heating, and Friction)
    out3 = [outdir, 'Forcing_', num2str(ycur), '.nc'];
    write_group(out3, lat, lev_hPa, Ndays, { ...
    'EPDphi_QG',   EPDphi_QG_year; ...
    'EPDphi1',     EPDphi1_year; ...
    'EPDp_QG',     EPDp_QG_year; ...
    'EPDp1',       EPDp1_year; ...
    'EPDp2',       EPDp2_year; ...
    'GDphi1',      GDphi1_year; ...
    'GDphi2',      GDphi2_year; ...
    'GDp1',        GDp1_year; ...
    'GDp2',        GDp2_year; ...
    'J',           J_year; ...
    'X',           X_year });

    % 4) Flux_yyyy.nc (EP-flux, G-vector, EMF/EHF/VEMF/VEHF
    out4 = [outdir, 'Flux_', num2str(ycur), '.nc'];
    write_group(out4, lat, lev_hPa, Ndays, { ...
    'EPphi_QG',    EPphi_QG_year; ...
    'EPphi1',      EPphi1_year; ...
    'EPp_QG',      EPp_QG_year; ...
    'EPp1',        EPp1_year; ...
    'EPp2',        EPp2_year; ...
    'Gphi1',       Gphi1_year; ...
    'Gphi2',       Gphi2_year; ...
    'Gp1',         Gp1_year; ...
    'Gp2',         Gp2_year; ...
    'EMF',         EMF_year; ...
    'EHF',         EHF_year; ...
    'VEMF',        VEMF_year; ...
    'VEHF',        VEHF_year });
    
    toc
end
