%% SOM_moistAdiabat
% Solve for the thermodynamic state along a moist adiabat.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   vars = SOM_moistAdiabat(h_L, q_T, q_p, z, p)
%
%%% Description
% Calculates thermodynamic properties along a moist adiabatic using the
% thermodynamics from the System for Atmospheric
% Modeling and the SAM single-moment microphysics scheme. 
% In SAM, moist adiabatic processes are those in which the moist
% liquid water/ice static energy is conserved. In this function, properties
% are therefore calculated assuming that MSE remains constant and changes
% in temperature and mixing ratios compensate for changes in height.
% Because these computations contain no sense of passing time, the SAM
% formulations of conversion between non-precipitating to precipitating
% condensate cannot be accounted for in the calculations, and the
% calculations instead assume that phase changes convert water vapor to and
% from cloud water.
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
% *q_T - initial total nonprecipitating water:*
% Total nonprecipitating water specific humidity, as defined in Appendix A 
% of Khairoutdinov and Randall, 2003, in the initial state.
% Must in units of
% kilograms of water per kilogram of fluid. If moist adiabats are computed
% for multiple columns, q_T must be provided as a 2D matrix with one value
% per column
%
% *q_p - initial total precipitating water:*
% Total precipitating water specific humidity, as defined in Appendix A of 
% Khairoutdinov and Randall, 2003, in the initial state.
% Must be in units of kilograms
% of water per kilogram of fluid. If moist adiabats are computed for multiple
% columns, q_p must be provided as a 2D matrix with one value per column.
%
% *z - heights:*
% Heights along the moist adiabat. Must be in units of meters. This argument
% can be provided in one of two formats. If z is given as a 3D matrix with the
% first and second dimensions equal to the size of h_L, q_T, and q_p, vectors
% along the third dimehsion of z are used as heights for each column, and moist
% adiabats are calculated for each column at the heights in the corresponding
% vector. If z is given as a 3D matrix with singleton first and second
% dimensions, the third dimension is broadcast to each column. (Note that row
% and column vectors are interpreted per the first format, not the second.)
%
% *p - pressures:*
% Pressures at each height given in z. Must in units of Pascals.
%
%%% Output Arguments
% *vars - model variables:*
% Struct with fields qn (total cloud condensate specific humidity) and T
% (absolute temperature) at each point along the moist adiabats. The fields
% are of size [size(h_L, 1), size(h_L, 2), size(z, 3)]. qn is given in units
% of kg water per kg total air and T is given in Kelvins.
%
%%% <../test/html/SOM_moistAdiabat_test.html Tests>
function res = SOM_moistAdiabat(h_L, q_T, q_p, z, p)  
	
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

    res = SOM_getState(h_L, q_T, q_p, z, p);
end
