%% SAM_moistAdiabat
% Solve for the thermodynamic state along a moist adiabat.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   vars = SAM_moistAdiabat(h_L, q_T, q_p, q_v, z_i, p_i, z_f, p_f)
%
%   vars = SAM_moistAdiabat(..., Name, Value)
%
%%% Description
% Calculates thermodynamic properties along a moist adiabatic using the
% thermodynamics from the System for Atmospheric
% Modeling. In SAM, moist adiabatic processes are those in which the moist
% liquid water/ice static energy is conserved. In this function, properties
% are therefore calculated assuming that MSE remains constant and changes
% in temperature and mixing ratios compensate for changes in height.
% Because these computations contain no sense of passing time, the SAM
% formulations of conversion between non-precipitating to precipitating
% condensate cannot be accounted for in the calculations, and the
% calculations instead assume either that all condensed water is converted
% immediately to precipitation or that none of it is.
%
% Additional details about SAM are given in in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
%
%%% Input Arguments
% *h_L - liquid water/ice moist static energy:*
% Liquid water/ice moist static energy as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be scalar and in units of Joules
% per kilogram.
%
% *qT - initial total nonprecipitating water:*
% Total nonprecipitating water specific humidity, as defined in Appendix A 
% of Khairoutdinov and Randall, 2003, in the initial state.
% Must be scalar and in units of
% kilograms of water per kilogram of fluid.
%
% *q_p - initial total precipitating water:*
% Total precipitating water specific humidity, as defined in Appendix A of 
% Khairoutdinov and Randall, 2003, in the initial state.
% Must be scalar and in units of kilograms
% of water per kilogram of fluid.
%
% *q_v - initial water vapor specific humidity:*
% Water vapor specific humidity, as defined in Appendix A of 
% Khairoutdinov and Randall, 2003, in the initial state.
% Must be scalar and in units of kilograms
% of water per kilogram of fluid. If new condensate is non-precipitable
% (the default, or specified by setting the 'AutoConversion' property to
% 'nonprecipitable'), this argument is unused and can be set to any numeric
% value without affecting the result.
%
% *z_f - final height:*
% Final height of the parcel. Must be scalar and in units of meters.
%
% *p_f - final pressure:*
% Final pressure of the parcel Must be scalar and in units of meters.
%
%%% Name-Value Pair Arguments
%
% *'Format' - variable output format:*
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
% *'AutoConversion' - treatment of new condensation:*
% Specify how to deal with the parcel water budget when the water vapor
% mixing ratio changes. By default, or if 'nonprecipitable' is specified,
% non-precipitable water is conserved and changes in water vapor are offset
% by changes in non-precipitable cloud condensate. If 'precipitable' is
% specified, water vapor plus precipitable water is conserved and changes
% in water vapor are offset by changes in precipitable water. In both
% cases, this function calculates a reversible adiabat (precipitation
% travels with the parcel rather than falling out).
%
%%% Output Arguments
% *vars - a vector of final-state variables:*
% 
% *varsWithNames - a struct with final-state variable name as keys:*
%
% *solver - an NDSecantSolver for the final state:*
%
%%% <../test/html/SAM_moistAdiabat_test.html Tests>
function res = SAM_moistAdiabat(h_L, q_T, q_p, q_v, ...
    z_f, p_f, varargin)  
    
    % Parse input
    p = inputParser;
    checkIn = @(x) assert(isnumeric(x) && isscalar(x), ...
        'SAM_moistAdiabat() only accepts numeric scalar input');
    validFmt = {'vector', 'struct', 'solver'};
    checkFmt = @(x) assert(ischar(x) && any(strcmp(x, validFmt)), ...
        'Invalid output format');
    validAC = {'nonprecipitable', 'precipitable'};
    checkAC = @(x) assert(ischar(x) && any(strcmp(x, validAC)), ...
        'Invalid autoconversion specifier');
    addRequired(p, 'h_L', checkIn);
    addRequired(p, 'q_T', checkIn);
    addRequired(p, 'q_p', checkIn);
    addRequired(p, 'q_v', checkIn);
    addRequired(p, 'z_f', checkIn);
    addRequired(p, 'p_f', checkIn);
    addParameter(p, 'Format', validFmt{1}, checkFmt);
    addParameter(p, 'AutoConversion', validAC{1}, checkAC);
    parse(p, h_L, q_T, q_p, q_v, z_f, p_f, varargin{:});
    fmt = p.Results.Format;
    ac = p.Results.AutoConversion;
    
    % Create solver based on condensate autoconversion option
    switch ac
        
        % Case 1: all new condensate is non-precipitable
        % IMPLEMENT THIS
        case validAC{1}
    
            res = SAM_getState(h_L, q_T, q_p, z_f, p_f, fmt);
        
        % Case 2: all new condensate is precipitable
        case valid
            
            % Declare global constants
            global SAM_c_p; %#ok<*TLEV>
            global SAM_g;
            global SAM_T_0n;
            global SAM_T_00n;
            global SAM_T_0p;
            global SAM_T_00p;
            global SAM_L_s;
            global SAM_L_c;
            
            error('SAM_moistAdiabat:unimplemented', ...
                ['Moist adiabat computations where new condensate is ',...,
                'precipitable have not been implemented.']);
            
        otherwise
            error('SAM_getState:illegalParameter', ...
                '%s is not a valid option for ''AutoConversion''', ...
                ac);
    end

end

function q = qsat(T, p)
        
    global SAM_T_00n;
    global SAM_T_0n;
    global SAM_R;
    global SAM_R_v;

    if T < SAM_T_00n
        q = SAM_R * psatIceFWC(T) / (SAM_R_v * p);
    elseif T > SAM_T_0n
        q = SAM_R * psatWaterFWC(T) / (SAM_R_v * p);
    else
        q = SAM_omega(T, SAM_T_0n, SAM_T_00n) * ...
            (SAM_R*psatWaterFWC(T))/(SAM_R_v*p) + ...
        (1 - SAM_omega(T, SAM_T_0n, SAM_T_00n)) * ...
            (SAM_R*psatIceFWC(T))/(SAM_R_v*p);
    end
end