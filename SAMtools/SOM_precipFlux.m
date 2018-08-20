%% SOM_precipFlux
% Calculate precipitation fluxes based on the SAM single-moment 
% microphysics scheme.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/30/2016
%
%%% Syntax
%   [f, f_r, f_s, f_g] = SOM_precipFlux(p, T, q_p)
%
%%% Description
% Calculates downward precipitation fluxes given a pressure, a (virtual)
% temperature, a precipitable water specific humidity. The calculations are
% performed based on the SAM single-moment microphysics scheme.
% Details about the SAM model and microphysics are given in in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
% All inputs can be multidimensional. If multidimensional input
% is given, all inputs must be the same shape and size and the function will
% operate element-wise.
%
%%% Input Arguments
% *p - pressure:*
% Used to compute air density. Must be in units of Pa.
%
% *T - temperature:*
% Used to compute air density. To include the effect of water vapor on air
% density, provide a virtual temperature; otherwise, the density will be
% calculated assuming dry air. Must be in units of Kelvins.
% 
% *q_p - total precipitating water:*
% Total precipitating water specific humidity as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be in units of kilograms
% of water per kilogram of fluid.
%
%%% Output Arguments
% *f - total precipitation flux:*
% Total instantaneous downward precipitation flux (rain, snow, and graupel
% combined). Has units of kg m^2 / s.
%
% *f_r - precipitation flux from rain:*
% Contribution to f from rain. Has units of kg m^2 / s.
%
% *f_s - precipitation flux from snow:*
% Contribution to f from snow. Has units of kg m^2 / s.
%
% *f_g - precipitation flux from graupel:*
% Contribution to f from graupel. Has units of kg m^2 / s.
%
%%% <../test/html/SOM_precipFlux_test.html Tests>
%
%%% Source code
function [f, f_r, f_s, f_g] = SOM_precipFlux(p, T, q_p)

    % Add global variables
    global SOM_a_r;
    global SOM_a_s;
    global SOM_a_g;
    global SOM_b_r;
    global SOM_b_s;
    global SOM_b_g;
    global SOM_N_0r;
    global SOM_N_0s;
    global SOM_N_0g;
    global SOM_T_0p;
    global SOM_T_00p;
    global SOM_T_0g;
    global SOM_T_00g;
    global SOM_rho_r;
    global SOM_rho_s;
    global SOM_rho_g;
    global SOM_rho_0;
    global SAM_R;
    
    % Compute density
    rho = p ./ (SAM_R .* T);
    
    % Compute mixing ratios of different forms of precipitation
    op = SOM_omega(T, SOM_T_0p, SOM_T_00p);
    q_r = op.*q_p;
    og = SOM_omega(T, SOM_T_0g, SOM_T_00g);
    q_s = (1 - op).*(1 - og).*q_p;
    q_g = (1 - op).*og.*q_p;
    
    % Compute fluxes from different forms of precipitation
    f_r = (SOM_a_r .* gamma(4 + SOM_b_r))./6 .* ...
        (pi.*SOM_rho_r.*SOM_N_0r).^(-SOM_b_r/4) .* ...
        (SOM_rho_0 ./ rho).^(0.5) .* ...
        (rho .* q_r).^(1 + SOM_b_r/4);
    f_s = (SOM_a_s .* gamma(4 + SOM_b_s))./6 .* ...
        (pi.*SOM_rho_s.*SOM_N_0s).^(-SOM_b_s/4) .* ...
        (SOM_rho_0 ./ rho).^(0.5) .* ...
        (rho .* q_s).^(1 + SOM_b_s/4);
    f_g = (SOM_a_g .* gamma(4 + SOM_b_g))./6 .* ...
        (pi.*SOM_rho_g.*SOM_N_0g).^(-SOM_b_g/4) .* ...
        (SOM_rho_0 ./ rho).^(0.5) .* ...
        (rho .* q_g).^(1 + SOM_b_g/4);

	% Compute total flux
    f = f_r + f_s + f_g;

end
