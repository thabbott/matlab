%% rad_cam_wrapper
% MATLAB interface to CAM radiation scheme used in SAM

addpath('..');
SAM_defineConstants;
SOM_defineConstants;

% Load data
dat = SAM_loadVars({'x', 'y', 'z', 'p', 'TABS', 'QV', 'QN', 'QRAD'}, ...
    '/data/thabbott/SimulationData/SAM/RCE_96x96x64_CAM/OUT_3D/perpetual_eqref_TS300_LAT19.45_DAY80.5/perpetual_eqref_TS300_LAT19.45_DAY80.5_96x96x64_1km_8s_6.10.6.CAM_32_0001080000.nc', 1);

% Set up grid
grid = SG_grid(dat.x, dat.y, dat.z, dat.p);
grid.tabs = SG_addVar(dat.tabs, 's');
grid.qv = SG_addVar(dat.qv, 's');
om = SAM_omega(dat.tabs, SOM_T_0n, SOM_T_00n);
grid.qcl = SG_addVar(dat.qn.*om, 's');
grid.qci = SG_addVar(dat.qn.*(1-om), 's');
latitude = 19.45*ones(grid.nx, grid.ny);
longitude = zeros(grid.nx, grid.ny);
sst = 300*ones(grid.nx, grid.ny);
grid.latitude = SG_addVar(latitude, 's2d');
grid.longitude = SG_addVar(longitude, 's2d');
grid.sst = SG_addVar(sst, 's2d');

% Define parameters
params = struct;
params.day = 80.5;
params.day0 = 80.5;
params.dt = 8;
params.nrad = 30;
params.doperpetual = 1;
params.ocean = 1;
params.dolongwave = 1;
params.doshortwave = 1;
params.doseasons = -1;
params.doperpetual = 1;
params.dosolarconstant = 1;
params.doradhomo = 1;
params.zenith_angle = 42.2;
params.solar_constant = 554.0;

% Compare qrad from simulation and offline calculation
grid.qrad_online = SG_addVar(dat.qrad, 's');
grid.qrad_offline = SG_rad_cam(grid, params);

%%

nx = 16;
ny = 16;
nzm = 16;
nz = nzm + 1;

zi = 1:nz;
presi = 1:nz;
day = 120;
day0 = 80.5;
dt = 8;
nrad = 1;
nstep = 6;
doperpetual = -1;
ocean = 1;
dolongwave = 1;
doshortwave = 1;
doseasons = -1;
dosolarconstant = -1;
doradhomo = -1;
zenith_angle = 0.5;
solar_constant = 450;
latitude = 20*ones(nx,ny);
longitude = 10*ones(nx,ny);
sstxy = zeros(nx,ny);
rho = ones(nzm,1);
tabs = 300*ones(nx,ny,nzm);
qv = 0.0001*rand(nx,ny,nzm);
qcl = 0.0001*rand(nx,ny,nzm);
qci = 0.0001*rand(nx,ny,nzm);

qrad = mx_rad_cam_wrapper(...
    zi, ...
    presi, ...
    day, ...
    day0, ...
    dt, ...
    nrad, ...
    nstep, ...
    doperpetual, ...
    ocean, ...
    dolongwave, ...
    doshortwave, ...
    doseasons, ...
    dosolarconstant, ...
    doradhomo, ...
    zenith_angle, ...
    solar_constant, ...
    latitude, ...
    longitude, ...
    sstxy, ...
    rho, ...
    tabs, ...
    qv, ...
    qcl, ...
    qci);