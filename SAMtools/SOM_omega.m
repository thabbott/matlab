%% SOM_omega
% Water partitioning function from the System for Atmospheric Modeling.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
%
%%% Syntax
%   w = SOM_omega(T, T0, T00)
%
%%% Description
% Calculates a partitioning factor used in the System for Atmospheric
% Modeling's single-moment microphysics scheme, 
% as defined in Equation A13 in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
%
%%% Input Arguments
% *T - temperature.*
% Can be either scalar or non-scalar. If T is non-scalar, the function
% operates element-wise.
%
% *T0, T00 - upper and lower temperature cutoffs:*
% Scalar temperatures used to define the range over which the partitioning
% factor varies.
%
%
%%% Output Arguments
% *w - partitioning factor:*
% The partitioning factor defined based on the temperature cutoffs. w = 0
% for T < T00, w = 1 for T > T0, and w increases linearly from 0 to 1
% between T00 and T0. If T is non-scalar, the output is element-wise.
%
%%% <../test/html/SOM_omega_test.html Tests>
function w = SOM_omega(T, T0, T00)

%     % Parse input
%     p = inputParser;
%     checkT = @(x) assert(isnumeric(x), 'T must be numeric');
%     checkT0 = @(x) assert(isscalar(x) && isnumeric(x), ...
%         'Temperature cutoff must be a numeric scalar');
%     addRequired(p, 'T', checkT);
%     addRequired(p, 'T0', checkT0);
%     addRequired(p, 'T00', checkT0);
%     parse(p, T, T0, T00);
    
    % Assign result
    w = max(0, min(1, (T - T00)/(T0 - T00)));
    
end
