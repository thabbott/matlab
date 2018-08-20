%% SAM_defineConstants
% Define constants used in the System for Atmospheric Modeling
%
%%% Syntax
%   SAM_defineConstants
%
%%% Description
% Defines (a subset) of the constants used in the System for Atmospheric
% Modeling, version 6.10.8. 
% The naming convention is as follows: each variable is prefixed
% with 'SAM_', greek letters are spelled out, and subscripts are denoted
% by an underscore, e.g. SAM_g,
% SAM_R, etc. All units are MKS.
% 
% These constants are declared as global variables. Once you have run this
% script in one workspace, you can access in any constants in it in another
% workspace by declaring them global variables, e.g.
%
%   global SAM_R;
%   p = rho * SAM_R * T;
%
% Note: some of these constants are not used in SAM but are included in this
% file because I use them when working with SAM output. These variables are
%
% - SAM_c_pv

%%% Specific heat capacity at constant pressure of dry air
global SAM_c_p; SAM_c_p = 1.004e3; % J kg^{-1} K^{-1}
%%% Specific heat capacity at constant pressure of water vapor
global SAM_c_pv; SAM_c_pv = 1.864e3; % J kg^{-1} K^{-1}
%%% Gravitational acceleration
global SAM_g; SAM_g = 9.81; % m s^{-2}
%%% Latent heat of condensation
global SAM_L_c; SAM_L_c = 2.5104e6; % J kg^{-1}
%%% Latent heat of fusion
global SAM_L_f; SAM_L_f = 0.3336e6; % J kg^{-1}
%%% Latent heat of sublimation
global SAM_L_s; SAM_L_s = 2.8440e6; % J kg^{-1}
%%% Specific gas constant for air
global SAM_R; SAM_R = 287.; % J kg^{-1} K^{-1}
%%% Specific gas constant for water vapor
global SAM_R_v; SAM_R_v = 461.; % J kg^{-1} K^{-1}
