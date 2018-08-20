%% SAM_dpsatIce
% Temperature derivative of saturation vapor pressure over ice
%
%%% Syntax
%   de = SAM_dpsatIce(T)
%
%%% Description
% Calculates the temperature derivative of
% saturation vapor pressure over ice with a
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
% *e - temperature derivative of saturation vapor pressure (Pa/K):*
% Temperature derivative of saturation vapor pressure.
%
%%% <../test/html/SAM_dpsatIce_test.html Tests>

function de = SAM_dpsatIce(T)  

    % Define coefficients
    a = [                      ...
        0.497275778e-16, ...
	0.390204672e-13, ...
	0.132073448e-10, ...
	0.255653718e-8, ...
	0.312668753e-6, ...
	0.249065913e-4, ...
	0.126710138e-2, ...
	 0.377174432e-1, ...
	  0.503223089 ...
    ];

    % Calculate polynomial (with additional interpolation)
    de = zeros(size(T));
    
    m1 = logical(T > 273.15);
    de(m1) = SAM_dpsatWater(T(m1));

    m2 = logical(T > 185 & ~m1);
    Tt = T(m2) - 273.16;
    de(m2) = 100 * polyval(a, Tt);

    m3 = ~(m1 | m2);
    Tt = max(-100, T(m3) - 273.16);
    de(m3) = 100 * (0.0013186 + Tt.*(2.60269e-05 + Tt*1.28676e-07));

end
