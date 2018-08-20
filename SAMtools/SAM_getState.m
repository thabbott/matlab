%% SAM_getState
% Solve for the thermodynamic state based on prognostic variables.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   vars = SAM_getState(h_L, q_T, q_p, z, p)
%
%   varsWithNames = SAM_getState(h_L, q_T, q_p, z, p, 'struct')
%
%   solver = SAM_getState(h_L, q_T, q_p, z, p, 'solver')
%
%%% Description
% Calculates thermodynamic properties used in the System for Atmospheric
% Modeling for a parcel of air based on the model's prognostic 
% thermodynamic variables (moist liquid
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
% *p_tot - pressure:*
% Pressure of the parcel. Must be scalar and in units of Pascals.
%
% *format - variable output format (optional):*
% Specify the format of the function output ('vector', 'struct', 'solver'). 
% By default, the output is a vector. Specifying 'struct' will result in
% the output being formatted as a struct with variable names as keys, and
% specifying 'solver' will result in the function returning an
% <NDSecantSolver.html NDSecantSolver> that can be used to solve for the 
% state. *_If a solver is returned, it is returned fully constrained but 
% unsolved._* This output option is intended to
% allow a user to specify their own initial guesses for the state if the
% default initial guesses used in this function do not converge to a
% solution and to avoid the cost of re-instantiating an identical solver
% if SAM_getState is used for multiple computations.
%
%%% Output Arguments
% *vars - a vector of model variables:*
% 
% *varsWithNames - a struct with variable names as keys:*
%
% *solver - an NDSecantSolver for the parcel state:*
%
%%% <../test/html/SAM_getState_test.html Tests>
function res = SAM_getState(h_L, q_T, q_p, z, p_tot, varargin)
    
    % Parse input
    p = inputParser;
    checkIn = @(x) assert(isnumeric(x) && isscalar(x), ...
        'SAM_getState() only accepts numeric scalar input');
    validFmt = {'vector', 'struct', 'solver'};
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
    s.addVar('q_v');
    s.addVar('q_n');
    s.addVar('q_sat');
    
    % Add and set parameters
    s.addParam('h_L');
    s.addParam('q_T');
    s.addParam('q_p');
    s.addParam('z');
    s.addParam('p');
    s.setParams([h_L, q_T, q_p, z, p_tot]);
    
    % Add functions
    % Load SAM constants
    % SAM_defineConstants;
    
    % First: definition of moist static energy
    
    % Second: % Definition of saturation specific humidity
        % For now, compute assuming dilute limit: q \approx r,
        % p_a \approx p_total. Could add additional solver constraints
        % to compute explicitly.
        
    % Third: no super-saturation
    
    % Fourth: definition of non-precipitating condensate
    global SAM_c_p;
    global SAM_g;
    global SOM_T_0n;
    global SOM_T_00n;
    global SOM_T_0p;
    global SOM_T_00p;
    global SAM_L_s;
    global SAM_L_c;  
    
    s.setFunctions(@(x, p) ...
        [SAM_c_p * x(1) + SAM_g * p(4) - ... % Temperature and height
        SAM_L_c * ( ...
            SAM_omega(x(1), SOM_T_0n, SOM_T_00n)*x(3) + ...
            SAM_omega(x(1), SOM_T_0p, SOM_T_00p)*p(3) ...
        ) - ... % Liquid water latent heat
        SAM_L_s * ( ...
            (1 - SAM_omega(x(1), SOM_T_0n, SOM_T_00n))*x(3) + ...
            (1 - SAM_omega(x(1), SOM_T_0p, SOM_T_00p))*p(3) ...
        ) ... % Ice latent heat
        - p(1); ... % END DEF'N OF MSE
        qsat(x(1), p(5)) ...
        - x(4); ... % END DEF'N OF SATURATION SPECIFIC HUMIDITY
        max(0, p(2) - x(4)) - x(3); ... % END NO SUPERSATURATION
        p(2) - x(2) - x(3) ... % END DEF'N OF NON-PREC CONDENSATE
        ] ...
    );
    
    % If solver requested, return
    if strcmp(fmt, 'solver')
        res = s;
        return;
    end
    
    % Otherwise, insert default guesses
    Tg = 280;
    qg = 1e-2;
    Tp = 1;
    qp = 1e-3;
    setVars(s, [Tg, qg, qg, qg], 'current');
    setVars(s, [Tg-Tp, qg-qp, qg-qp, qg-qp], 'previous');
    
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

function q = qsat(T, p)
        
    global SOM_T_00n;
    global SOM_T_0n;
    global SAM_R;
    global SAM_R_v;

    if T < SOM_T_00n
        q = SAM_R * psatIceFWC(T) / (SAM_R_v * p);
    elseif T > SOM_T_0n
        q = SAM_R * psatWaterFWC(T) / (SAM_R_v * p);
    else
        q = SAM_omega(T, SOM_T_0n, SOM_T_00n) * ...
            (SAM_R*psatWaterFWC(T))/(SAM_R_v*p) + ...
        (1 - SAM_omega(T, SOM_T_0n, SOM_T_00n)) * ...
            (SAM_R*psatIceFWC(T))/(SAM_R_v*p);
    end
end