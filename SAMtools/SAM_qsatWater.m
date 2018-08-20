%% SAM_qsatWater
% Saturation specific humidity over liquid water
%
%%% Syntax
%   q = SAM_qsatWater(T, p)
%
%%% Description
% Calculates the saturation specific humidity over liquid water with a
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
% *p - pressure (Pa):*
% Ambient pressure. Must have the same dimensions as T.
%
%%% Output Arguments
% *q - saturation specific humidity (nondim):*
% Saturation specific humidity, in units of kg water per total kg.
%
%%% <../test/html/SAM_qsatWater_test.html Tests>

function q = SAM_qsatWater(T, p)  

    e = SAM_psatWater(T);
    q = 0.622 * e ./ max(e, p-e);

end
