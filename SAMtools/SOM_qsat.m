%% SOM_qsat
% Saturation specific humidity calculated as a temperature-dependent linear
% combination of values over ice and water.
%
%%% Syntax
%   q = SOM_qsat(T, p)
%
%%% Description
% Calculates saturation specific humidities based on
% polynomial fits and a linear combination of values over ice and water.
% The polynomial is based on Flatau, Walko, and
% Cotton (1992): "Polynomial Fits to Saturation Vapor Pressure" and is
% identical to the one used in the System for Atmospheric Modeling, version
% 6.10.8. The linear combination is described in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
%%% Input Arguments
% *T - temperature (K):*
% May be either scalar or non-scalar. If non-scalar, the output has the same
% size and shape as the input.
%
% *p - pressure (Pa):*
% Ambient pressure. Must have the same dimensions as T.
%
%%% Output Arguments
% *q - saturation specific humidity (nondim):*
% Saturation specific humidity, in units of kg water per total kg.
%

function q = SOM_qsat(T, p)  
    
    global SOM_T_0n; global SOM_T_00n;
    om = SAM_omega(T, SOM_T_0n, SOM_T_00n);
    q = om.*SAM_qsatWater(T, p) + (1-om).*SAM_qsatIce(T, p);

end
