clear;
clc;
close all
%% Coefficients
Omega       = 7.292E-5;
p0          = 1000*100;
start_year  = 1979;
end_year    = 2024;
year        = start_year:1:end_year;
N_period    = length(year);
a           = 6371000;
Rd          = 287;
cp          = 1004;

%% Load Data for Longitude, Latitude, and Pressure
path            = '/jsbach/s0/mup65/ERA5/';
lon             = ncread([path,'UCOMP/u_daily_mean.1979.nc'],'longitude');
lon             = double (lon);
lat             = ncread([path,'UCOMP/u_daily_mean.1979.nc'],'latitude');
lat             = double (lat);
lev             = ncread([path,'UCOMP/u_daily_mean.1979.nc'],'level');
lev             = double(lev)*100;
p_ini           = lev;          
phi             = deg2rad(lat);     % Change degree to rad
p1              = linspace(0, 1000*100, 41);
p2              = linspace(0, 850*100, 69);
p3              = linspace(1*100, 50*100, 50);
p1(1)           = 100;
p2(1)           = 100;
lev             = cat(2,p3/100,p2(6:end)/100,p1(36:end)/100);
N_lev           = length(lev);

sinphi          = sin(phi);
cosphi          = cos(phi);
tanphi          = tan(phi);
cosphi(1,:)     = 1e-3;       % Avoid singularity for cosphi near -pi/2 and pi/2.
cosphi(end,:)   = 1e-3;       % Avoid singularity for cosphi near -pi/2 and pi/2.
f               = 2*Omega*sinphi;

% Turn 1D array into 2D array
[~,f_map_1]             = meshgrid(p1,f);
[p_map_1,cosphi_map_1]  = meshgrid(p1,cosphi);
[~,tanphi_map_1]        = meshgrid(p1,tanphi);

[~,f_map_2]             = meshgrid(p2,f);    
[p_map_2,cosphi_map_2]  = meshgrid(p2,cosphi);
[~,tanphi_map_2]        = meshgrid(p2,tanphi);

[~,f_map_3]             = meshgrid(p3,f);    
[p_map_3,cosphi_map_3]  = meshgrid(p3,cosphi);
[~,tanphi_map_3]        = meshgrid(p3,tanphi);

for m = 11:N_period
% Loading data
disp(['Currently Running - year: ', num2str(year(m))]);
if m == 1
    u_back = [];v_back = [];w_back = [];T_back = [];
    u   = ncread([path,'UCOMP/u_daily_mean.',num2str(year(m)),'.nc'], 'u');
    v   = ncread([path,'VCOMP/v_daily_mean.',num2str(year(m)),'.nc'], 'v');
    w   = ncread([path,'OMEGA/w_daily_mean.',num2str(year(m)),'.nc'], 'w');
    T   = ncread([path,'TEMP/T_daily_mean.',num2str(year(m)),'.nc'], 't'); 
end
if m ~= 1
    u_back   = ncread([path,'UCOMP/u_daily_mean.',num2str(year(m)-1),'.nc'], 'u');
    v_back   = ncread([path,'VCOMP/v_daily_mean.',num2str(year(m)-1),'.nc'], 'v');
    w_back   = ncread([path,'OMEGA/w_daily_mean.',num2str(year(m)-1),'.nc'], 'w');
    T_back   = ncread([path,'TEMP/T_daily_mean.',num2str(year(m)-1),'.nc'], 't'); 
    u   = ncread([path,'UCOMP/u_daily_mean.',num2str(year(m)),'.nc'], 'u');
    v   = ncread([path,'VCOMP/v_daily_mean.',num2str(year(m)),'.nc'], 'v');
    w   = ncread([path,'OMEGA/w_daily_mean.',num2str(year(m)),'.nc'], 'w');
    T   = ncread([path,'TEMP/T_daily_mean.',num2str(year(m)),'.nc'], 't'); 
    %{
    u_back   = u;
    v_back   = v;
    w_back   = w;
    T_back   = T; 
    u   = u_forward;
    v   = v_forward;
    w   = w_forward;
    T   = T_forward;
    %}
end
if m ~= N_period
    u_forward   = ncread([path,'UCOMP/u_daily_mean.',num2str(year(m)+1),'.nc'], 'u');
    v_forward   = ncread([path,'VCOMP/v_daily_mean.',num2str(year(m)+1),'.nc'], 'v');
    w_forward   = ncread([path,'OMEGA/w_daily_mean.',num2str(year(m)+1),'.nc'], 'w');
    T_forward   = ncread([path,'TEMP/T_daily_mean.',num2str(year(m)+1),'.nc'], 't');
end
if m == N_period
    u_forward   = [];
    v_forward   = [];
    w_forward   = [];
    T_forward   = [];
end
u = cat(4,u_back(:,:,:,end),u,u_forward(:,:,:,1));
v = cat(4,v_back(:,:,:,end),v,v_forward(:,:,:,1));
w = cat(4,w_back(:,:,:,end),w,w_forward(:,:,:,1));
T = cat(4,T_back(:,:,:,end),T,T_forward(:,:,:,1));

disp('loading data completed')

[N_x, N_y, N_p, ~, ~] = size(u);
if mod(year(m),4) == 0
    N_day = 366;
else
    N_day = 365;
end

FJ_year          = zeros(N_y, N_lev, N_day);
FF_year          = zeros(N_y, N_lev, N_day);
FEPphi_QG_year   = zeros(N_y, N_lev, N_day);
FEPphi1_year     = zeros(N_y, N_lev, N_day);
FEPp_QG_year     = zeros(N_y, N_lev, N_day);
FEPp1_year       = zeros(N_y, N_lev, N_day);
FEPp2_year       = zeros(N_y, N_lev, N_day);
FGphi1_year      = zeros(N_y, N_lev, N_day);
FGphi2_year      = zeros(N_y, N_lev, N_day); 
FGp1_year        = zeros(N_y, N_lev, N_day);
FGp2_year        = zeros(N_y, N_lev, N_day);
FTEM_year        = zeros(N_y, N_lev, N_day);
tempsi_year      = zeros(N_y, N_lev, N_day);
eumpsi_year      = zeros(N_y, N_lev, N_day);

HJ_year          = zeros(N_y, N_lev, N_day);
HF_year          = zeros(N_y, N_lev, N_day);
HEPDphi_QG_year  = zeros(N_y, N_lev, N_day);
HEPDphi1_year    = zeros(N_y, N_lev, N_day);
HEPDp_QG_year    = zeros(N_y, N_lev, N_day);
HEPDp1_year      = zeros(N_y, N_lev, N_day);
HEPDp2_year      = zeros(N_y, N_lev, N_day);
HGDphi1_year     = zeros(N_y, N_lev, N_day);
HGDphi2_year     = zeros(N_y, N_lev, N_day);
HGDp1_year       = zeros(N_y, N_lev, N_day);
HGDp2_year       = zeros(N_y, N_lev, N_day);

EPDphi_QG_year   = zeros(N_y, N_lev, N_day);
EPDphi1_year     = zeros(N_y, N_lev, N_day);
EPDp_QG_year     = zeros(N_y, N_lev, N_day);
EPDp1_year       = zeros(N_y, N_lev, N_day);
EPDp2_year       = zeros(N_y, N_lev, N_day);

Fphi_QG_year     = zeros(N_y, N_lev, N_day);
Fphi1_year       = zeros(N_y, N_lev, N_day);
Fp_QG_year       = zeros(N_y, N_lev, N_day);
Fp1_year         = zeros(N_y, N_lev, N_day);
Fp2_year         = zeros(N_y, N_lev, N_day);

Gphi1_year       = zeros(N_y, N_lev, N_day);
Gphi2_year       = zeros(N_y, N_lev, N_day);
Gp1_year         = zeros(N_y, N_lev, N_day);
Gp2_year         = zeros(N_y, N_lev, N_day);

J_year           = zeros(N_y, N_lev, N_day);
Fx_year          = zeros(N_y, N_lev, N_day);

% Computation of potential temperature
theta = T;
for k = 1:N_p
    theta(:,:,k,:) = T(:,:,k,:)*((p0/p_ini(k))^(Rd/cp));
end

%% Computation starts here
for day = 1:N_day
    tic
    day
    %% Compuation of EP-Flux and Other Forcings
    % Interpolation before Computation
    % Starting from the near-surface layer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [FTEM_1, FJ_1, FF_1, FEPphi_QG_1, FEPphi1_1, FEPp_QG_1, FEPp1_1, FEPp2_1, ...
    FGphi1_1, FGphi2_1, FGp1_1, FGp2_1, tempsi_1, eumpsi_1, ...
    HJ_1, HF_1, ...
    HEPDphi_QG_1, HEPDphi1_1, HEPDp_QG_1, HEPDp1_1, HEPDp2_1, ...
    HGDphi1_1, HGDphi2_1, HGDp1_1, HGDp2_1, ...
    EPDphi_QG_1, EPDphi1_1, EPDp_QG_1, EPDp1_1, EPDp2_1, ...
    Fphi_QG_1, Fphi1_1, Fp_QG_1, Fp1_1, Fp2_1, ...
    Gphi1_1, Gphi2_1, Gp1_1, Gp2_1, ...
    J_1, Fx_1] = solve_TEM_KE_daily...
    (u, v, w, T, theta,...
    m, day,...
    N_x, N_y, N_period, ...
    p_ini, p1, phi, p_map_1, f_map_1, cosphi_map_1, tanphi_map_1,...
    0, 1000,...
    zeros(N_y,1), zeros(N_y,1), zeros(N_y,1), zeros(N_y,1), ...
    1, zeros(N_y,1), zeros(N_y,1),...
    zeros(N_y,1), zeros(N_y,1), zeros(N_y,1), zeros(N_y,1));

    %% Now, compute the troposphere
    [FTEM_2, FJ_2, FF_2, FEPphi_QG_2, FEPphi1_2, FEPp_QG_2, FEPp1_2, FEPp2_2, ...
    FGphi1_2, FGphi2_2, FGp1_2, FGp2_2, tempsi_2, eumpsi_2, ...
    HJ_2, HF_2, ...
    HEPDphi_QG_2, HEPDphi1_2, HEPDp_QG_2, HEPDp1_2, HEPDp2_2, ...
    HGDphi1_2, HGDphi2_2, HGDp1_2, HGDp2_2, ...
    EPDphi_QG_2, EPDphi1_2, EPDp_QG_2, EPDp1_2, EPDp2_2, ...
    Fphi_QG_2, Fphi1_2, Fp_QG_2, Fp1_2, Fp2_2, ...
    Gphi1_2, Gphi2_2, Gp1_2, Gp2_2, ...
    J_2, Fx_2] = solve_TEM_KE_daily...
    (u, v, w, T, theta,...
    m, day,...
    N_x, N_y, N_period, ...
    p_ini, p2, phi, p_map_2, f_map_2, cosphi_map_2, tanphi_map_2,...
    0, 850,...
    FJ_1(:,35), FF_1(:,35), FEPphi_QG_1(:,35), FEPphi1_1(:,35), ...
    FEPp_QG_1(:,35), FEPp1_1(:,35), FEPp2_1(:,35),...
    FGphi1_1(:,35), FGphi2_1(:,35), FGp1_1(:,35), FGp2_1(:,35));
    %% Finally, the startosphere
    [FTEM_3, FJ_3, FF_3, FEPphi_QG_3, FEPphi1_3, FEPp_QG_3, FEPp1_3, FEPp2_3, ...
    FGphi1_3, FGphi2_3, FGp1_3, FGp2_3, tempsi_3, eumpsi_3, ...
    HJ_3, HF_3, ...
    HEPDphi_QG_3, HEPDphi1_3, HEPDp_QG_3, HEPDp1_3, HEPDp2_3, ...
    HGDphi1_3, HGDphi2_3, HGDp1_3, HGDp2_3, ...
    EPDphi_QG_3, EPDphi1_3, EPDp_QG_3, EPDp1_3, EPDp2_3, ...
    Fphi_QG_3, Fphi1_3, Fp_QG_3, Fp1_3, Fp2_3, ...
    Gphi1_3, Gphi2_3, Gp1_3, Gp2_3, ...
    J_3, Fx_3] = solve_TEM_KE_daily...
    (u, v, w, T, theta,...
    m, day,...
    N_x, N_y, N_period, ...
    p_ini, p3, phi, p_map_3, f_map_3, cosphi_map_3, tanphi_map_3,...
    1, 50,...
    FJ_2(:,5), FF_2(:,5), FEPphi_QG_2(:,5), FEPphi1_2(:,5), ...
    FEPp_QG_2(:,5), FEPp1_2(:,5), FEPp2_2(:,5),...
    FGphi1_2(:,5), FGphi2_2(:,5), FGp1_2(:,5), FGp2_2(:,5));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save daily data into array
    FJ_year(:,:,day)          = cat(2,FJ_3,FJ_2(:,6:end),FJ_1(:,36:end));
    FF_year(:,:,day)          = cat(2,FF_3,FF_2(:,6:end),FF_1(:,36:end));
    FEPphi_QG_year(:,:,day)   = cat(2,FEPphi_QG_3,FEPphi_QG_2(:,6:end),FEPphi_QG_1(:,36:end));
    FEPphi1_year(:,:,day)     = cat(2,FEPphi1_3,FEPphi1_2(:,6:end),FEPphi1_1(:,36:end));
    FEPp_QG_year(:,:,day)     = cat(2,FEPp_QG_3,FEPp_QG_2(:,6:end),FEPp_QG_1(:,36:end));
    FEPp1_year(:,:,day)       = cat(2,FEPp1_3,FEPp1_2(:,6:end),FEPp1_1(:,36:end));
    FEPp2_year(:,:,day)       = cat(2,FEPp2_3,FEPp2_2(:,6:end),FEPp2_1(:,36:end));
    FGphi1_year(:,:,day)      = cat(2,FGphi1_3,FGphi1_2(:,6:end),FGphi1_1(:,36:end));
    FGphi2_year(:,:,day)      = cat(2,FGphi2_3,FGphi2_2(:,6:end),FGphi2_1(:,36:end));
    FGp1_year(:,:,day)        = cat(2,FGp1_3,FGp1_2(:,6:end),FGp1_1(:,36:end));
    FGp2_year(:,:,day)        = cat(2,FGp2_3,FGp2_2(:,6:end),FGp2_1(:,36:end));

    J_year(:,:,day)           = cat(2,J_3,J_2(:,6:end),J_1(:,36:end));
    Fx_year(:,:,day)          = cat(2,Fx_3,Fx_2(:,6:end),Fx_1(:,36:end));
    EPDphi_QG_year(:,:,day)    = cat(2,EPDphi_QG_3,EPDphi_QG_2(:,6:end),EPDphi_QG_1(:,36:end));
    EPDphi1_year(:,:,day)      = cat(2,EPDphi1_3,EPDphi1_2(:,6:end),EPDphi1_1(:,36:end));
    EPDp_QG_year(:,:,day)      = cat(2,EPDp_QG_3,EPDp_QG_2(:,6:end),EPDp_QG_1(:,36:end));
    EPDp1_year(:,:,day)        = cat(2,EPDp1_3,EPDp1_2(:,6:end),EPDp1_1(:,36:end));
    EPDp2_year(:,:,day)        = cat(2,EPDp2_3,EPDp2_2(:,6:end),EPDp2_1(:,36:end));

    Fphi_QG_year(:,:,day)    = cat(2,Fphi_QG_3,Fphi_QG_2(:,6:end),Fphi_QG_1(:,36:end));
    Fphi1_year(:,:,day)      = cat(2,Fphi1_3,Fphi1_2(:,6:end),Fphi1_1(:,36:end));
    Fp_QG_year(:,:,day)      = cat(2,Fp_QG_3,Fp_QG_2(:,6:end),Fp_QG_1(:,36:end));
    Fp1_year(:,:,day)        = cat(2,Fp1_3,Fp1_2(:,6:end),Fp1_1(:,36:end));
    Fp2_year(:,:,day)        = cat(2,Fp2_3,Fp2_2(:,6:end),Fp2_1(:,36:end));
    Gphi1_year(:,:,day)       = cat(2,Gphi1_3,Gphi1_2(:,6:end),Gphi1_1(:,36:end));
    Gphi2_year(:,:,day)       = cat(2,Gphi2_3,Gphi2_2(:,6:end),Gphi2_1(:,36:end));
    Gp1_year(:,:,day)         = cat(2,Gp1_3,Gp1_2(:,6:end),Gp1_1(:,36:end));
    Gp2_year(:,:,day)         = cat(2,Gp2_3,Gp2_2(:,6:end),Gp2_1(:,36:end));

    HJ_year(:,:,day)          = cat(2,HJ_3,HJ_2(:,6:end),HJ_1(:,36:end));
    HF_year(:,:,day)          = cat(2,HF_3,HF_2(:,6:end),HF_1(:,36:end));
    HEPDphi_QG_year(:,:,day)   = cat(2,HEPDphi_QG_3,HEPDphi_QG_2(:,6:end),HEPDphi_QG_1(:,36:end));
    HEPDphi1_year(:,:,day)     = cat(2,HEPDphi1_3,HEPDphi1_2(:,6:end),HEPDphi1_1(:,36:end));
    HEPDp_QG_year(:,:,day)     = cat(2,HEPDp_QG_3,HEPDp_QG_2(:,6:end),HEPDp_QG_1(:,36:end));
    HEPDp1_year(:,:,day)       = cat(2,HEPDp1_3,HEPDp1_2(:,6:end),HEPDp1_1(:,36:end));
    HEPDp2_year(:,:,day)       = cat(2,HEPDp2_3,HEPDp2_2(:,6:end),HEPDp2_1(:,36:end));
    HGDphi1_year(:,:,day)      = cat(2,HGDphi1_3,HGDphi1_2(:,6:end),HGDphi1_1(:,36:end));
    HGDphi2_year(:,:,day)      = cat(2,HGDphi2_3,HGDphi2_2(:,6:end),HGDphi2_1(:,36:end));
    HGDp1_year(:,:,day)        = cat(2,HGDp1_3,HGDp1_2(:,6:end),HGDp1_1(:,36:end));
    HGpD2_year(:,:,day)        = cat(2,HGDp2_3,HGDp2_2(:,6:end),HGDp2_1(:,36:end));

    FTEM_year(:,:,day)        = cat(2,FTEM_3,FTEM_2(:,6:end),FTEM_1(:,36:end));
    tempsi_year(:,:,day)      = cat(2,tempsi_3,tempsi_2(:,6:end),tempsi_1(:,36:end));
    eumpsi_year(:,:,day)      = cat(2,eumpsi_3,eumpsi_2(:,6:end),eumpsi_1(:,36:end));
    
    toc
end

%% Save all Kuo-Eliassen TEM results to NetCDF
output_file = ['/jsbach/s0/pkh5218/TEM/daily/TEM_' num2str(year(m)) '.nc'];

% Define dimensions
nccreate(output_file, 'lat', 'Dimensions', {'lat', N_y});
nccreate(output_file, 'lev', 'Dimensions', {'lev', length(lev)});
nccreate(output_file, 'day', 'Dimensions', {'day', N_day});

% Define coordinate variables
ncwrite(output_file, 'lat', lat);
ncwrite(output_file, 'lev', lev);
ncwrite(output_file, 'day', 1:N_day);

% Add attributes
ncwriteatt(output_file, 'lat', 'units', 'degrees');
ncwriteatt(output_file, 'lev', 'units', 'hPa');
ncwriteatt(output_file, 'day', 'units', 'days since Jan 1');

% List of variables to save
var_list = {
    'FTEM', FTEM_year;
    'tempsi', tempsi_year;
    'eumpsi', eumpsi_year;
    'FJ', FJ_year;
    'FF', FF_year;
    'FEPphi_QG', FEPphi_QG_year;
    'FEPphi1', FEPphi1_year;
    'FEPp_QG', FEPp_QG_year;
    'FEPp1', FEPp1_year;
    'FEPp2', FEPp2_year;
    'FGphi1', FGphi1_year;
    'FGphi2', FGphi2_year;
    'FGp1', FGp1_year;
    'FGp2', FGp2_year;
};
% Save each variable
for i = 1:size(var_list,1)
    var_name = var_list{i,1};
    var_data = var_list{i,2};
    nccreate(output_file, var_name, ...
        'Dimensions', {'lat', N_y, 'lev', length(lev), 'day', N_day});
    ncwrite(output_file, var_name, var_data);
end

%% Save all Kuo-Eliassen TEM results to NetCDF
output_file = ['/jsbach/s0/pkh5218/TEM/daily/Forcing_' num2str(year(m)) '.nc'];

% Define dimensions
nccreate(output_file, 'lat', 'Dimensions', {'lat', N_y});
nccreate(output_file, 'lev', 'Dimensions', {'lev', length(lev)});
nccreate(output_file, 'day', 'Dimensions', {'day', N_day});

% Define coordinate variables
ncwrite(output_file, 'lat', lat);
ncwrite(output_file, 'lev', lev);
ncwrite(output_file, 'day', 1:N_day);

% Add attributes
ncwriteatt(output_file, 'lat', 'units', 'degrees');
ncwriteatt(output_file, 'lev', 'units', 'hPa');
ncwriteatt(output_file, 'day', 'units', 'days since Jan 1');

% List of variables to save
var_list = {
    'HJ', HJ_year;
    'HF', HF_year;
    'HEPDphi_QG', HEPDphi_QG_year;
    'HEPDphi1', HEPDphi1_year;
    'HEPDp_QG', HEPDp_QG_year;
    'HEPDp1', HEPDp1_year;
    'HEPDp2', HEPDp2_year;
    'HGDphi1', HGDphi1_year;
    'HGDphi2', HGDphi2_year;
    'HGDp1', HGDp1_year;
    'HGDp2', HGDp2_year;
};
% Save each variable
for i = 1:size(var_list,1)
    var_name = var_list{i,1};
    var_data = var_list{i,2};
    nccreate(output_file, var_name, ...
        'Dimensions', {'lat', N_y, 'lev', length(lev), 'day', N_day});
    ncwrite(output_file, var_name, var_data);
end

%% Save all Kuo-Eliassen TEM results to NetCDF
output_file = ['/jsbach/s0/pkh5218/TEM/daily/Budget_' num2str(year(m)) '.nc'];

% Define dimensions
nccreate(output_file, 'lat', 'Dimensions', {'lat', N_y});
nccreate(output_file, 'lev', 'Dimensions', {'lev', length(lev)});
nccreate(output_file, 'day', 'Dimensions', {'day', N_day});

% Define coordinate variables
ncwrite(output_file, 'lat', lat);
ncwrite(output_file, 'lev', lev);
ncwrite(output_file, 'day', 1:N_day);

% Add attributes
ncwriteatt(output_file, 'lat', 'units', 'degrees');
ncwriteatt(output_file, 'lev', 'units', 'hPa');
ncwriteatt(output_file, 'day', 'units', 'days since Jan 1');

% List of variables to save
var_list = {
    'J', J_year;
    'Fx', Fx_year;
    'EPDphi_QG', EPDphi_QG_year;
    'EPDphi1', EPDphi1_year;
    'EPDp_QG', EPDp_QG_year;
    'EPDp1', EPDp1_year;
    'EPDp2', EPDp2_year;
    'Fphi_QG', Fphi_QG_year;
    'Fphi1', Fphi1_year;
    'Fp_QG', Fp_QG_year;
    'Fp1', Fp1_year;
    'Fp2', Fp2_year;
    'Gphi1', Gphi1_year;
    'Gphi2', Gphi2_year;
    'Gp1', Gp1_year;
    'Gp2', Gp2_year;
};
% Save each variable
for i = 1:size(var_list,1)
    var_name = var_list{i,1};
    var_data = var_list{i,2};
    nccreate(output_file, var_name, ...
        'Dimensions', {'lat', N_y, 'lev', length(lev), 'day', N_day});
    ncwrite(output_file, var_name, var_data);
end


end

load('TEMmap2.mat')
plot_figure_zmlog(lat,lev,mean(FTEM_year,3)*365,mymap,8,[-80 80]);
plot_figure_zmlog(lat,lev,mean(tempsi_year,3)*365,mymap,8,[-80 80]);

A1 = mean(FTEM_year,3);
A2 = mean(tempsi_year,3);
Error = (A1-A2)./A2*100;
plot_figure_zm(lat,lev(2:end),Error(:,2:end),mymap,100,[-80 80]);






