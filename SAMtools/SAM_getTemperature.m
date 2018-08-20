%% SAM_getTemperature
% Solve for the temperature based on thermodynamic prognostic variables.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   T = SAM_getState(h_L, q_T, q_p, z, p)
%
%   solver = SAM_getState(h_L, q_T, q_p, z, p, 'solver')
%
%%% Description
% Calculates the temperature of a parcel of air based on the prognostic 
% thermodynamic variables using in the System for Atmospheric Modeling
% (moist liquid
% water/ice moist static energy, total nonprecipitating water, and total
% precipitating water) and the parcel height.
% Details about the SAM model are given in in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
%
%%% Input Arguments
% *h - liquid water/ice moist static energy:*
% Liquid water/ice moist static energy as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be scalar and in units of Joules
% per kilogram.
%
% *qT - total nonprecipitating water:*
% Total nonprecipitating water specific humidity as defined in Appendix A 
% of Khairoutdinov and Randall, 2003. Must be scalar and in units of
% kilograms of water per kilogram of fluid.
%
% *qp - total precipitating water:*
% Total precipitating water specific humidity as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be scalar and in units of kilograms
% of water per kilogram of fluid.
%
% *z - height:*
% Height of the parcel. Must be scalar and in units of meters.
%
% *p - pressure:*
% Pressure of the parcel. Must be scalar and in units of Pascals.
%
% *format - variable output format (optional):*
% Specify the format of the function output ('vector', 'solver'). 
% By default, the output is a vector. 
% Specifying 'solver' will result in the function returning an
% <NDSecantSolver.html NDSecantSolver> that can be used to solve for the 
% state. *_If a solver is returned, it is returned fully constrained but 
% unsolved._* This output option is intended to
% allow a user to specify their own initial guesses for the state if the
% default initial guesses used in this function do not converge to a
% solution and to avoid the cost of re-instantiating an identical solver
% if SAM_getState is used for multiple computations.
%
%%% Output Arguments
% *T - parcel temperature:*
%
% *solver - an NDSecantSolver for the parcel temperature:*
%
%%% <../test/html/SAM_getTemperature_test.html Tests>
function res = SAM_getTemperature(h_L, q_T, q_p, z, p_tot, varargin)
    
    % Parse input
    p = inputParser;
    checkIn = @(x) assert(isnumeric(x) && isscalar(x), ...
        'SAM_getState() only accepts numeric scalar input');
    validFmt = {'vector', 'solver'};
    checkFmt = @(x) assert(ischar(x) && any(strcmp(x, validFmt)), ...
        'Invalid output format');
    addRequired(p, 'h_L', checkIn);
    addRequired(p, 'q_T', checkIn);
    addRequired(p, 'q_p', checkIn);
    addRequired(p, 'z', checkIn);
    addRequired(p, 'p_tot', checkIn);
    addOptional(p, 'format', validFmt{1}, checkFmt);
    parse(p, h_L, q_T, q_p, z, p_tot, varargin{:});
    fmt = p.Results.format;
    
    % Create solver
    s = NDSecantSolver;
    
    % Add variables
    s.addVar('T');
    
    % Add and set parameters
    s.addParam('h_L');
    s.addParam('q_T');
    s.addParam('q_p');
    s.addParam('z');
    s.addParam('p');
    s.setParams([h_L, q_T, q_p, z, p_tot]);
    
    % Add functions  
    s.setFunctions(@hL);
    
    % If solver requested, return
    if strcmp(fmt, 'solver')
        res = s;
        return;
    end
    
    % Otherwise, insert default guesses
    Tg = 280;
    Tp = 1;
    setVars(s, Tg, 'current');
    setVars(s, Tg-Tp, 'previous');
    
    % Solve
    maxiter = 100;
    try
        solve(s, maxiter);
    catch ME
        switch ME.identifier
            case 'NDSecantSolver:nonConvergence'
                error('NDSecantSolver:nonConvergence', ...
                    ['SAM state solver did not converge with default ',...
                    'guesses within %d iterations. Try returning a ',...
                    'solver by calling ',...
                    'SAM_getState(..., ''solver'') and entering your ',...
                    'own guesses'], maxiter);
            otherwise
                rethrow(ME);
        end
    end
    
    % Return solution
    res = getVars(s, 'current', fmt);

end

function h = hL(x, p)

    global SOM_T_00n;
    global SOM_T_0n;
    global SAM_R;
    global SAM_R_v;
    global SAM_c_p;
    global SAM_g;
    global SOM_T_0p;
    global SOM_T_00p;
    global SAM_L_s;
    global SAM_L_c;

    if x < SOM_T_00n
        q = SAM_R * psatIceFWC(x) / (SAM_R_v * p(5));
    elseif x > SOM_T_0n
        q = SAM_R * psatWaterFWC(x) / (SAM_R_v * p(5));
    else
        q = SAM_omega(x, SOM_T_0n, SOM_T_00n) * ...
            (SAM_R*psatWaterFWC(x))/(SAM_R_v*p(5)) + ...
        (1 - SAM_omega(x, SOM_T_0n, SOM_T_00n)) * ...
            (SAM_R*psatIceFWC(x))/(SAM_R_v*p(5));
    end
    
    h = SAM_c_p * x + SAM_g * p(4) - ... % Temperature and height
        SAM_L_c * ( ...
            SAM_omega(x, SOM_T_0n, SOM_T_00n)*max(0, p(2)-q) + ...
            SAM_omega(x, SOM_T_0p, SOM_T_00p)*p(3) ...
        ) - ... % Liquid water latent heat
        SAM_L_s * ( ...
            (1 - SAM_omega(x, SOM_T_0n, SOM_T_00n))*max(0, p(2)-q) + ...
            (1 - SAM_omega(x, SOM_T_0p, SOM_T_00p))*p(3) ...
        ) - p(1);
    
end