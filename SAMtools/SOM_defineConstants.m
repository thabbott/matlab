%% SOM_defineConstants
% Define constants used in the System for Atmospheric Modeling's
% single-moment microphysics scheme.
%
% Tristan Abbott //
% Massachusetts Institute of Technology //
% 01/23/2017
%
%%% Syntax
%   SOM_defineConstants
%
%%% Description
% Defines the constants used in the System for Atmospheric
% Modeling's single-moment microphysics scheme, 
% as defined in Appendix B of Khairoutdinov and Randall, 2003: 
% "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
% The naming convention is as follows: each variable is prefixed
% with 'SOM_', greek letters are spelled out, subscripts are denoted
% by an underscore, and capitalization follows Appendix B, e.g. SOM_g,
% SOM_C_g, SOM_alpha, SOM_rho_0. Units are as given in Appendix B.
% 
% These constants are declared as global variables. Once you have run this
% script in one workspace, you can access in any constants in it in another
% workspace by declaring them global variables, e.g.
%
%   global SOM_rho_0;
%   m = SOM_rho_0

%%% Constant in fall speed formula for rain
global SOM_a_r; SOM_a_r = 842; % m^{1-b_r}/s
%%% Constant in fall speed formula for snow
global SOM_a_s; SOM_a_s = 4.84; %m^{1-b_s}/s
%%% Constant in fall speed formula for graupel
global SOM_a_g; SOM_a_g = 94.5; %m^{1-b_g}/s
%%% Exponent in fall speed formula for rain
global SOM_b_r; SOM_b_r = 0.8; % nondim
%%% Exponent in fall speed formula for snow
global SOM_b_s; SOM_b_s = 0.25; % nondim
%%% Exponent in fall speed formula for graupel
global SOM_b_g; SOM_b_g = 0.5; % nondim
%%% Intercept parameter for rain
global SOM_N_0r; SOM_N_0r = 8e6; % m^{-4}
%%% Intercept parameter for snow
global SOM_N_0s; SOM_N_0s = 3e6; % m^{-4}
%%% Intercept parameter for graupel
global SOM_N_0g; SOM_N_0g = 4e6; % m^{-4}
%%% Specific gas constant for air
global SOM_R; SOM_R = 287.; % J kg^{-1} K^{-1}
%%% Specific gas constant for water vapor
global SOM_R_v; SOM_R_v = 461.; % J kg^{-1} K^{-1}
%%% Temperature threshold for ice
global SOM_T_0n; SOM_T_0n = 273.16; % K
%%% Temperature threshold for snow/graupel
global SOM_T_0p; SOM_T_0p = 283.16; % K
%%% Temperature threshold for graupel
global SOM_T_0g; SOM_T_0g = 283.16; % K
%%% Temperature threshold for cloud water
global SOM_T_00n; SOM_T_00n = 253.16; % K
%%% Temperature threshold for rain
global SOM_T_00p; SOM_T_00p = 268.16; % K
%%% Temperature threshold for graupel
global SOM_T_00g; SOM_T_00g = 223.16; % K
%%% Reference air density
global SOM_rho_0; SOM_rho_0 = 1.29; % kg/m^3
%%% Density of rain
global SOM_rho_r; SOM_rho_r = 1000; % kg/m^3
%%% Density of snow
global SOM_rho_s; SOM_rho_s = 100; % kg/m^3
%%% Density of graupel
global SOM_rho_g; SOM_rho_g = 400; % kg/m^3
%%% Autoconversion rate
global SOM_alpha; SOM_alpha = 0.001; % 1/s
%%% Ice aggregation rate
global SOM_beta; SOM_beta = 0.001; % 1/s
%%% Threshold cloud water for autoconversion
global SOM_q_co; SOM_q_co = 1e-3; % kg/kg
%%% Threshold ice for aggregation
global SOM_q_io; SOM_q_io = 1e-4; % kg/kg
%%% Collection efficiency of rain for cloud water
global SOM_E_rc; SOM_E_rc = 1.0; % nondim
%%% Collection efficiency of snow for cloud water
global SOM_E_sc; SOM_E_sc = 1.0; % nondim
%%% Collection efficiency of graupel for cloud water
global SOM_E_gc; SOM_E_gc = 1.0; % nondim
%%% Collection efficiency of rain for cloud ice
global SOM_E_ri; SOM_E_ri = 1.0; % nondim
%%% Collection efficiency of snow for cloud ice
global SOM_E_si; SOM_E_si = 0.1; % nondim
%%% Collection efficiency of graupel for cloud ice
global SOM_E_gi; SOM_E_gi = 0.1; % nondim
