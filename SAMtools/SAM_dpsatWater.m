%% SAM_dpsatWater
% Temperature derivative of saturation vapor pressure over liquid water
%
%%% Syntax
%   de = SAM_dpsatWater(T)
%
%%% Description
% Calculates the temperature derivative of
% saturation vapor pressure over liquid water with a
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
%%% <../test/html/SAM_dpsatWater_test.html Tests>

function de = SAM_dpsatWater(T)  

    % Define coefficients
    a = [                      ...
            -0.599634321e-17, ...
	    -0.792933209e-14, ...
	    -0.604119582e-12, ...
	    0.385208005e-9, ...
	    0.103167413e-6, ...
	    0.121167162e-4, ...
	    0.794747212e-3, ...
	    0.285976452e-1, ...
	    0.443956472 ...
    ];

    % Calculate polynomial
    T = max(-80, T - 273.16);
    de = 100 * polyval(a, T);

end
