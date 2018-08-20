%% SG_rad_cam
% Call CAM radiation routines
%
% Tristan Abbott // Massachusetts Institute of Technology // 8/20/2018
%
%%% Syntax
%   qrad = SG_rad_cam(grid, params)

%%% Description
% Calls the CAM radiation routines used in SAM version 6.10.6 using data
% stored in a SG_grid object and some additional required parameters.
%
%%% Input Arguments
% *grid - SG_grid object:*
% Model grid and variables returned by a call to SG_grid. Required
% variables are 'tabs' (absolute temperature, K), 'qv' 
% (specific humidity, kg/kg), 'qcl' (cloud water, kg/kg), 'qci'
% (cloud ice, kg/kg), 'latitude' (2d latitude), 'longitude' (2d longitude),
% and 'sst' (2d sea surface temperature, K).
%
% *params - struct of miscellaneous parameters:*
% Additional parameters required by CAM routines but not provided by
% SG_grid. Must include 'day' (Julian day at current time step), 'day0'
% (Julian day at simulation start), 'dt' (model time step), 'nrad' (number
% of time steps per radiation call), 'doperpetual', 'ocean', 'dolongwave',
% 'doshortwave', 'doseasons', 'dosolarconstant', 'doradhomo' (all SAM
% parameters, set to 1 for true and -1 for false), 'zenith_angle' (solar
% zenith angle for perpetual insolation simulations), and 'solar_constant'
% (solar constant for perpetual insolation simulations).
%
% Note that dt and nrad are used only with perpetual insolation for 
% calculating modified solar constants. 
% 
%%% Output Arguments
% *qrad - radiative heating rate:*
% Radiative heating rates at scalar cells, suitable for assignment to an
% SG_grid field.
% 
%%% 
function qrad = SG_rad_cam(grid, params)

    % Prepare fields and parameters
    tabs = check_3d_field(grid, 'tabs');
    qv = check_3d_field(grid, 'qv');
    qcl = check_3d_field(grid, 'qcl');
    qci = check_3d_field(grid, 'qci');
    latitude = check_2d_field(grid, 'latitude');
    longitude = check_2d_field(grid, 'longitude');
    sstxy = check_2d_field(grid, 'sst');
    day = check_param(params, 'day');
    day0 = check_param(params, 'day0');
    dt = check_param(params, 'dt');
    nrad = check_param(params, 'nrad');
    doperpetual = check_param(params, 'doperpetual');
    ocean = check_param(params, 'ocean');
    dolongwave = check_param(params, 'dolongwave');
    doshortwave = check_param(params, 'doshortwave');
    doseasons = check_param(params, 'doseasons');
    dosolarconstant = check_param(params, 'dosolarconstant');
    doradhomo = check_param(params, 'doradhomo');
    zenith_angle = check_param(params, 'zenith_angle');
    solar_constant = check_param(params, 'solar_constant');
    nx = grid.nx;
    ny = grid.ny;
    nz = grid.nz;
    nzm = grid.nzm;
    zi = grid.zi;
    rho = grid.rho;
    presi = grid.pi/1e2;
    nstep = 1;
    
    % Call CAM routines
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
    
    % Return as SAM variable
    qrad = SG_addVar(qrad, 's');

end

function field = check_3d_field(grid, fld)   
    if isfield(grid, fld)
        field = permute(grid.(fld), [3 2 1]);
    else
        error(['SG_rad_cam:no_' fld], ['Field ' fld ' not in grid']);
    end
end

function field = check_2d_field(grid, fld)   
    if isfield(grid, fld)
        field = permute(grid.(fld), [2 1]);
    else
        error(['SG_rad_cam:no_' fld], ['Field ' fld ' not in grid']);
    end
end

function field = check_param(grid, fld)   
    if isfield(grid, fld)
        field = grid.(fld);
    else
        error(['SG_rad_cam:no_' fld], ['Field ' fld ' not in params']);
    end
end
    