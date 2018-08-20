%% SAM_hL
% Liquid water/ice moist static energy from the System for Atmospheric
% Modeling
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   h = SAM_hL(T, z, qn, qp)
%
%%% Description
% Calculates the liquid water/ice moist static energy  of a parcel
% as defined in 
% Appendix A of Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_:
%
% $$
% h_L = c_p T + g_z - L_c \left ( \omega_n q_n + \omega_p q_p \right ) - 
% L_s \left ( (1-\omega_n)q_n + (1 - \omega_p) q_p \right ).
% $$
%
%%% Input Arguments
% *T - absolute temperature (K).*
%
% *qn - total nonprecipitating condensate (nondim).*
% Total nonprecipitating condensate as defined in Appendix A 
% of Khairoutdinov and Randall, 2003. Must be in units of
% kilograms of water per kilogram of fluid.
%
% *qp - total precipitating water (nondim):*
% Total precipitating water as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be in units of kilograms
% of water per kilogram of fluid.
%
% *z - height (m).*
% Height of the parcel. Must be in units of meters.
%
%%% Output Arguments
% *h - liquid water/ice MSE (J/(kg K)).*
% If inputs are non-scalar, computations are performed element-wise.
%
%%% <../test/html/SAM_hL_test.html Tests>
function h = SAM_hL(T, z, qn, qp)

    % Parse input
    p = inputParser;
    check = @(x) assert(isnumeric(x), ...
        'Input must be numeric');
    addRequired(p, 'T', check);
    addRequired(p, 'z', check);
    addRequired(p, 'qn', check);
    addRequired(p, 'qp', check);
    parse(p, T, z, qn, qp);
     

    % Define SAM constants
    % SAM_defineConstants;
    
    % Calculate partitioning factors.
    global SOM_T_0n; global SOM_T_00n;
    global SOM_T_0p; global SOM_T_00p;
    wn = SOM_omega(T, SOM_T_0n, SOM_T_00n);
    wp = SOM_omega(T, SOM_T_0p, SOM_T_00p);
    
    % Calculate MSE
    global SAM_c_p;
    global SAM_g;
    global SAM_L_c;
    global SAM_L_s;
    h = SAM_c_p*T + SAM_g*z - ...
        SAM_L_c * (wn.*qn + wp.*qp) - ...
        SAM_L_s * ((1-wn).*qn + (1-wp).*qp);
end