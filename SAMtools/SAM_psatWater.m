%% SAM_psatWater
% Saturation vapor pressure over liquid water
%
%%% Syntax
%   e = SAM_psatWater(T)
%
%%% Description
% Calculates the saturation vapor pressure over liquid water with a
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
%%% <../test/html/SAM_psatWater_test.html Tests>

function e = SAM_psatWater(T)  

    % Define coefficients
    a = [                      ...
           -0.976195544e-15,   ...
	   -0.952447341e-13,   ...
	    0.640689451e-10,   ...
	     0.206739458e-7, ... 
	     0.302950461e-5, ...
	     0.264847430e-3, ...
	      0.142986287e-1, ...
	       0.443987641, ...
	        6.11239921, ...
    ];

    % Calculate polynomial
    T = max(-80, T - 273.16);
    e = 100 * polyval(a, T);
end
