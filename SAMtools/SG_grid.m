%% SG_grid
% Generate a data structure suitable for holding variables on the
% C-grid used by SAM.
%
% Tristan Abbott // Massachusetts Institute of Technology // 11/14/2017
%
%%% Syntax
%   grid = SAM_grid(x, y, z, p)

%%% Description
% Generates and returns a data structure containing the C-grid used by SAM.
% SG_addVar can be used to add variables to the grid, and other functions
% are available for performing computations on the grid.
% 
% Most versions of SAM do not output reference density profiles, and this
% function estimates density profiles from height and pressure grids. However,
% more accurate budgets are possible if the reference profiles are saved
% and added to the grid using grid.rho = SG_addVar(..., 'rho') and
% grid.rhoi = SG_addVar(..., 'rhoi').
%
%%% Input Arguments
% *x, y, z - model coordinates:*
% 1d arrays containing model coordinates at tracer cell centers, all in
% units of meters.
%
% *p - reference pressure:*
% 1d array containing the pressure (in mb) at each level in z.
% 
%%% Output Arguments
% *grid - c-grid data structure:*
% A struct containing fields spatial coordinates x, y, and z (all in
% meters) at velocity and tracer points on a C-grid and reference densities
% at interface and scalar heights. Fields are as shown below:
%
% Vertical grid
%
% zi(n+1) -------- w(n+1) = 0; reference density rhoi
%    z(n) ======== stores u, v, scalars; reference density rho
% ....
%    z(2) ======== u(2), v(2), s(2)
%   zi(2) -------- w(2)
%    z(1) ======== u(1), v(1), s(1)
%   zi(1) -------- w(1) = 0; ground level
%
% Horizontal grids
% 
% xi(1)     x(1)        xi(2)       x(2)    ...     x(n)    xi(n+1)
% -----------------------------------------------------------------
% u(1)      v(1)        u(2)        v(2)            v(n)    u(1)
%           w(1)                    w(2)            w(n)
%           s(1)                    s(2)            s(n)
%
% yi(n+1) -------- v(1)
%    y(n) ======== u(n), w(n), s(n)
% ....
%    y(2) ======== u(2), w(2), s(2)
%   yi(2) -------- v(2)
%    y(1) ======== u(1), w(1), s(1)
%   yi(1) -------- v(1)
% 
%%% 
function grid = SG_grid(x, y, z, p)

    global SAM_g;

    % Compute interface locations
    zi = zeros(length(z)+1, 1);
    zi(1) = 0;
    zi(2) = 0.5*(z(1) + z(2));
    for ii = 3:length(zi)-1
        zi(ii) = zi(ii-1) + 0.5*(z(ii) - z(ii-2));
    end
    zi(end) = zi(end-1) + (z(end) - z(end-1));
    xi = zeros(length(x)+1, 1);
    dx = x(2) - x(1);
    xi(1:end-1) = x - 0.5*dx;
    xi(end) = x(end) + 0.5*dx;
    yi = zeros(length(y)+1, 1);
    dy = y(2) - y(1);
    yi(1:end-1) = y - 0.5*dy;
    yi(end) = y(end) + 0.5*dy;

    % Compute scalar level thickness ratios
    dz = 0.5*(z(1)+z(2));
    adz = zeros(length(z),1);
    adz(1) = 1.;
    adz(2:end-1) = 0.5*(z(3:end)-z(1:end-2))/dz;
    adz(end) = (z(end)-z(end-1))/dz;
    
    % Compute density profiles
    rhow = zeros(size(zi));
    for ii = 2:length(z)
        rhow(ii) = 1e2/SAM_g*(p(ii-1) - p(ii))/(z(ii) - z(ii-1));
    end
    rhow(end) = 2*rhow(end-1) - rhow(end-2);
    rhow(1) = 2*rhow(2) - rhow(3);
    rho = 0.5*(rhow(2:end) + rhow(1:end-1));
    
    % Compute interface hydrostatic pressures
    presi = zeros(length(zi),1);
    presi(1) = 1e2*p(1) + SAM_g*rho(1)*(z(1) - zi(1));
    for ii = 2:length(presi)
       presi(ii) = presi(ii-1) - SAM_g*rho(ii-1)*(zi(ii) - zi(ii-1));
    end

    grid = struct;
    grid.nzm = length(z);
    grid.nz = grid.nzm + 1;
    grid.nx = length(x);
    grid.ny = length(y);
    grid.z = z;
    grid.zi = zi;
    grid.adz = adz;
    grid.dz = dz;
    grid.x = x;
    grid.xi = xi;
    grid.dx = dx;
    grid.y = y;
    grid.yi = yi;
    grid.dy = dy;
    grid.rho = rho;
    grid.rhoi = rhow;
    grid.p = p;
    grid.pi = presi;
    grid.pp = repmat(1e2*grid.p, [1 length(y) length(x)]);

end
