%% SAM_psatIce
% Saturation vapor pressure over ice
%
%%% Syntax
%   e = SAM_psatIce(T)
%
%%% Description
% Calculates the saturation vapor pressure over ice with a
% polynomial fit. The polynomial is based on Flatau, Walko, and
% Cotton (1992): "Polynomial Fits to Saturation Vapor Pressure" and is
% identical to the one used in the System for Atmospheric Modeling, version
% 6.10.8.
%
%%% Input Arguments
% *T - temperature (K):*
% May be either scalar or non-scalar. If non-scalar, the output has the same
% size and shape as the input. 
%
%%% Output Arguments
% *e - saturation vapor pressure (Pa):*
% Saturation vapor pressure in Pascals.
%
%%% <../test/html/SAM_psatIce_test.html Tests>

function e = SAM_psatIce(T)  

    % Define coefficients
    a = [                      ...
           0.252751365e-14, ...
	   0.146898966e-11, ...
	   0.385852041e-9, ...
	   0.602588177e-7, ...
	   0.615021634e-5, ...
	   0.420895665e-3, ...
	   0.188439774e-1, ...
	   0.503160820, ...
	   6.11147274 ...
    ];

    % Calculate polynomial (with additional interpolation)
    e = zeros(size(T));
    
    m1 = logical(T > 273.15);
    e(m1) = SAM_psatWater(T(m1));

    m2 = logical(T > 185 & ~m1);
    Tt = T(m2) - 273.16;
    e(m2) = 100 * polyval(a, Tt);

    m3 = ~(m1 | m2);
    Tt = max(-100, T(m3) - 273.16);
    e(m3) = 100 * (0.00763685 + Tt.*(0.000151069 + Tt*7.48215e-07));

end
