%% SOM_condensateIntegral
% Calculate a condensate integral over a reversible moist adiabat based on
% frozen MSE.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%  [C, dpqv, omega] = SOM_condensateIntegral(h_L, q_T, q_p, z, p)
%  [...] = SOM_condensateIntegral(..., 'domain', 'updrafts')
%
%%% Description
%
% Calculate a condensate integral over a moist adiabat defined based on the
% System for Atmospheric Modeling and its single-moment microphysics scheme.
% Specifically, this function numerically approximates
%
% $$
% C = \int_{z_0}^{z_1} \left ( \frac{\partial q_v}{\partial z} \right )_{h_L}
% w(z) \bar{\rho} \mathrm{d}\,z.
% $$
%
% $q_v$ is the saturation specific humidity, and its height derivative is taken
% along a moist adiabat defined by constant frozen MSE:
%
% $$
% h_L = c_p T + gz - L_c q_l - L_s q_i
% $$
%
% where $q_l$ and $q_i$ are specific humidities for liquid and ice water.
% $w$ is the vertical velocity, and $\bar{\rho}$ is a reference density profile.
%
% This definition of a moist adiabat is the one used in the System for
% Atmospheric Modeling, and the partitioning of the water into liquid and ice
% phases is based on the SAM single-moment microphysics scheme. Moist adiabatic
% profile are computed in this function using <html/SOM_moistAdiabat.html 
% SOM_moistAdiabat>.
%
% Additional details about SAM and in single-moment microphysics scheme 
% are given in in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
%
%%% Input Arguments
% *h_L - liquid water/ice moist static energy:*
% Liquid water/ice moist static energy as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be in units of Joules
% per kilogram. If moist adiabats are computed for multiple columns, MSE
% must be provided as a 2D matrix with one value per column.
%
% *q_T - total nonprecipitating water:*
% Total nonprecipitating water specific humidity, as defined in Appendix A 
% of Khairoutdinov and Randall, 2003, in each column.
% Must in units of
% kilograms of water per kilogram of fluid. If moist adiabats are computed
% for multiple columns, q_T must be provided as a 2D matrix with one value
% per column
%
% *q_p - initial total precipitating water:*
% Total precipitating water specific humidity, as defined in Appendix A of 
% Khairoutdinov and Randall, 2003, in each column.
% Must be in units of kilograms
% of water per kilogram of fluid. If moist adiabats are computed for multiple
% columns, q_p must be provided as a 2D matrix with one value per column.
%
% *z - heights:*
% Heights along the condensate integrals. Must be monotonic and in units of
% meters. The
% condensate integral will be evaluated numerically from z(1) to z(end)
% based on integrand values at these points. Increasing the resolution of the 
% grid in z will improve the accuracy of the integral both by increasing the
% number of points in the numeric integral and by increasing the accuracy of
% finite-difference approximations to $\partial q_z$.
%
% This argument
% can be provided in one of two formats. If z is given as a 3D matrix with the
% first and second dimensions equal to the size of h_L, q_T, and q_p, vectors
% along the third dimehsion of z are used as heights for each column, and
% integrals are calculated for each column over the heights in the corresponding
% vector. If z is given as a 3D matrix with singleton first and second
% dimensions, the third dimension is broadcast to each column. (Note that row
% and column vectors are interpreted per the first format, not the second.)
%
% *p - pressures:*
% Pressures at each height given in z. Must in units of Pascals.
%
% *w - vertical velocity profiles:*
% Vertical velocity profiles in each column. Must be in units of m/s. As for
% z, this argument can be provided either as a full 3D matrix (in which case
% condensate integrals are calculated using only a different vertical velocity
% profile for each column) or as a 3D matrix with singleton first and second
% dimensions (in which case condensate integrals will be calculated using the
% same vertical velocity profiles for each column). The size of this argument
% must be equal to either [1, 1, size(z, 3)] or [size(h_L, 1), size(h_L, 2),
% size(z, 3)].
%
%%% Parameters
% *domain - domain criterion for the integral:*
% Specify a criterion that sections of the profiles must meet in order to
% be included in the integral. Options are 'all' (include the entire
% profile) and 'updrafts' (include only sections of the profile where the
% vertical velocity is upwards).
%
% *mask - integral domain mask:* a mask (logical array) that can exclude 
% elements from condensation integrals. Allowable shapes are the same as 
% for z, and elements of integrands that are not masked will be set to 0
% before integrating.
%
%%% Output Arguments
% *C - condensate integral:*
% Value of the condensate integral over each column, in units of kg of water
% per square meter per second.
%
% *dzqv - specific humidity lapse rate:*
% Pressure derivative of water vapor specific humidity, in units of per Pa.
%
% *omega - pressure velocity:*
% Vertical velocities converted to pressure coordinates
%
%%% <../test/html/SOM_condensateIntegral_test.html Tests>
function [C, dzqv, w] = SOM_condensateIntegral(h_L, q_T, q_p, z, p, w,...
    varargin)  
	
    global SAM_g;
    global SAM_R;
    global SAM_R_v;
    eps = SAM_R / SAM_R_v;
    
    % Check for parameters
    domain = 0;
    masked = 0;
    mask = [];
    if numel(varargin) ~= 0
        for i = 2:2:numel(varargin)
            switch lower(varargin{i-1})
                case 'domain'
                    switch lower(varargin{i})
                        case 'updrafts'
                            domain = 1;
                        case 'all'
                        otherwise
                            error(['Unrecognized value for parameter',...
                                ' "domain": %s'], varargin{i});
                    end
                case 'mask'
                    masked = 1;
                    mask = varargin{i};
                otherwise
                    error('Unrecognized parameter: %s', varargin{i-1});
            end
        end
    end

    % Broadcast fields as needed
    h_L = repmat(h_L, [1 1 size(z, 3)]);
	q_T = repmat(q_T, [1 1 size(z, 3)]);
	q_p = repmat(q_p, [1 1 size(z, 3)]);
    if size(z, 1) == 1 && size(z, 2) == 1
        z = repmat(z, [size(h_L,1) size(h_L,2) 1]);
    end
    if size(p, 1) == 1 && size(p, 2) == 1
        p = repmat(p, [size(h_L,1) size(h_L,2) 1]);
    end
    if size(w, 1) == 1 && size(w, 2) == 1
        w = repmat(w, [size(h_L,1) size(h_L,2) 1]);
    end
    if size(mask, 1) == 1 && size(mask, 2) == 1
        mask = repmat(mask, [size(h_L,1) size(h_L,2) 1]);
    end
    
    % Calculate moist adiabatic profiles
    res = SOM_getState(h_L, q_T, q_p, z, p);

    % Compute moisture derivative
    qv = q_T - res.qn;
    dzqv = zeros(size(z));
    dzqv(:,:,1) = (qv(:,:,2) - qv(:,:,1)) ./ (z(:,:,2) - z(:,:,1));
    dzqv(:,:,2:end-1) = (qv(:,:,3:end) - qv(:,:,1:end-2)) ./ ...
        (z(:,:,3:end) - z(:,:,1:end-2));
    dzqv(:,:,end) = (qv(:,:,end) - qv(:,:,end-1)) ./ ...
        (z(:,:,end) - z(:,:,end-1));
    
    % Multiply by vertical velocity
    dzqvw = w.*dzqv;
    
    % Zero part of the integral if needed
    if domain == 1
        dzqvw(logical(w < 0)) = 0;
    end
    if masked == 1
       dzqvw(~mask) = 0; 
    end
        

    % Compute integral with trapezoidal approximation along z
    integrand = (p(:,:,2:end) - p(:,:,1:end-1)) .* ...
        (dzqvw(:,:,2:end) + dzqvw(:,:,1:end-1));
    C = sum(integrand, 3) ./ (2*SAM_g);
    
    % Compute pressure velocity
    res.T = res.T .* (qv + eps) ./ (eps.*(1 + qv));
    w = dz2dp(w, p, res.T);
    
    % Convert lapse rate to pressure coordinates
    dzqv = dp2dz(dzqv, p, res.T);
    
end
