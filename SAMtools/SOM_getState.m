%% SOM_getState
% Solve for the thermodynamic state based on SAM prognostic variables and the
% SAM single-moment microphysics scheme.
%
% Tristan Abbott // Massachusetts Institute of Technology // 01/16/2016
% Translated from cloud.f90 by Marat Khairoutdinov.
%
%%% Syntax
%   vars = SOM_getState(h_L, q_T, q_p, z, p)
%
%%% Description
% Calculates thermodynamic properties for a parcel of air based on the System
% for Atmospheric Modeling's prognostic 
% thermodynamic variables (moist liquid
% water/ice moist static energy, total nonprecipitating water, and total
% precipitating water) and the SAM single-moment microphysics scheme.
% Details about the SAM model and microphysics are given in in Appendix A of 
% Khairoutdinov and Randall, 2003: "Cloud Resolving
% Modeling of the ARM Summer 1997 IOP: Model Formulation, Results, 
% Uncertainties, and Sensitivities", _Journal of the Atmospheric Sciences_. 
%
% This function computes the same properties as <html/SAM_getState.html
% SAM_getState> but is based on the modified Newton-Raphson iterative method
% used in SAM rather than the <html/NDSecantSolver.html NDSecantSolver> class.
% This function can handle multi-dimensional input and is much, much faster
% than <html/SAM_getState.html SAM_getState>. Unlike the iterative method used
% in SAM, which terminates after 10 iterations and may not resolve limit
% limit cycles, this function falls back on a bisection solver if the
% solution does not converge within 100 iterations.
%
% All inputs can be multidimensional. If multidimensional input
% is given, all inputs must be the same shape and size and the function will
% operate element-wise.
%
% To maximize speed, no input validation is performed.
%
%
%%% Input Arguments
% *h - liquid water/ice moist static energy:*
% Liquid water/ice moist static energy as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be in units of Joules
% per kilogram.
%
% *q_T - total nonprecipitating water:*
% Total nonprecipitating water specific humidity as defined in Appendix A 
% of Khairoutdinov and Randall, 2003. Must be in units of
% kilograms of water per kilogram of fluid.
%
% *q_p - total precipitating water:*
% Total precipitating water specific humidity as defined in Appendix A of 
% Khairoutdinov and Randall, 2003. Must be in units of kilograms
% of water per kilogram of fluid.
%
% *z - height:*
% Height of the parcel. Must be in units of meters.
%
% *p_tot - pressure:*
% Pressure of the parcel. Must be in units of Pascals.
%
%%% Output Arguments
% *vars - model variables:*
% Struct with fields qn (total cloud condensate specific humidity) and T
% (absolute temperature). The fields are the same size and shape as the
% inputs.
%
%%% <../test/html/SOM_getState_test.html Tests>
function [res, qnd, Td] = SOM_getState(h_L, q_T, q_p, z, p_tot)

    % Add global variables
    global SAM_L_c;
    global SAM_L_s;
    global SAM_L_f;
    global SAM_c_p;
    global SAM_g;
    global SOM_T_0n;
    global SOM_T_00n;
    global SOM_T_0p;
    global SOM_T_00p;
    % global SOM_T_0g;
    % global SOM_T_00g;

	% Pre-compute some frequently-used values
	an = 1. / (SOM_T_0n - SOM_T_00n);
	bn = SOM_T_00n * an;
	ap = 1. / (SOM_T_0p - SOM_T_00p);
	bp = SOM_T_00p * ap;
	fac_cond = SAM_L_c / SAM_c_p;
	fac_sub = SAM_L_s / SAM_c_p;
	fac_fus = SAM_L_f / SAM_c_p;
	fac1 = fac_cond + (1 + bp)*fac_fus;
	fac2 = fac_fus*ap;
	% ag = 1. / (SOM_T_0g - SOM_T_00g);

	% Pre-allocate arrays
	qsatt = zeros(size(h_L));
    qn = zeros(size(h_L));

	% Make sure that specific humidities are positive
	q_T = max(0., q_T);

	% Compute initial guesses for temperature
	T = (h_L - SAM_g * z) / SAM_c_p;
	T1 = (T + fac1*q_p) ./ (1 + fac2*q_p);

	% Refine guess for warm cloud
	m1 = logical(T1 > SOM_T_0n);
	T1(m1) = T(m1) + fac_cond*q_p(m1);
	qsatt(m1) = SAM_qsatWater(T1(m1), p_tot(m1));

	% Refine guess for ice cloud
	m2 = logical(T1 < SOM_T_00n);
	T1(m2) = T(m2) + fac_sub*q_p(m2);
	qsatt(m2) = SAM_qsatIce(T1(m2), p_tot(m2));

	% Refine guess for mixed-phase cloud
	m3 = ~(m1 | m2);
	om = an*T1(m3) - bn;
	qsatt(m3) = om.* SAM_qsatWater(T1(m3), p_tot(m3)) + ...
		(1 - om).*SAM_qsatIce(T1(m3), p_tot(m3));

	% Test if condensation is possible
	% Iterate if it is
	m0 = logical(q_T > qsatt);
	Tt = T(m0);
	T1t = T1(m0);
	pt = p_tot(m0);
	qTt = q_T(m0);
	qpt = q_p(m0);
	qst = qsatt(m0);
	dqst = zeros(size(qTt));
	lstn = zeros(size(qTt));
	dlstn = zeros(size(qTt));
	lstp = zeros(size(qTt));
	dlstp = zeros(size(qTt));

	niter = 0;
	maxiter = 10;
	dtabs = 100*ones(size(qTt));
    % Debugging
    qnd = zeros(2*maxiter+1, 1);
    Td = zeros(2*maxiter+1, 1);
	while any(abs(dtabs) > 0.01) && niter < maxiter
		
		% Handle cloud condensate
		m1 = logical(T1t > SOM_T_0n);
		lstn(m1) = fac_cond;
		dlstn(m1) = 0;
		qst(m1) = SAM_qsatWater(T1t(m1), pt(m1));
		dqst(m1) = SAM_dqsatWater(T1t(m1), pt(m1));

		m2 = logical(T1t < SOM_T_00n);
		lstn(m2) = fac_sub;
		dlstn(m2) = 0;
		qst(m2) = SAM_qsatIce(T1t(m2), pt(m2));
		dqst(m2) = SAM_dqsatIce(T1t(m2), pt(m2));

		m3 = ~(m1 | m2);
		om = an*T1t(m3) - bn;
		lstn(m3) = fac_cond + (1-om)*fac_fus;
		dlstn(m3) = an*fac_fus;
		qst(m3) = om.*SAM_qsatWater(T1t(m3), pt(m3)) + ...
			(1-om).*SAM_qsatIce(T1t(m3), pt(m3));
		dqst(m3) = om.*SAM_dqsatWater(T1t(m3), pt(m3)) + ...
			(1-om).*SAM_dqsatIce(T1t(m3), pt(m3));

		% Handle precipitable condensate
		m1 = logical(T1t > SOM_T_0p);
		lstp(m1) = fac_cond;
		dlstp(m1) = 0;

		m2 = logical(T1t < SOM_T_00p);
		lstp(m2) = fac_sub;
		dlstp(m2) = 0;

		m3 = ~(m1 | m2);
		omp = ap*T1t(m3) - bp;
		lstp(m3) = fac_cond + (1-omp)*fac_fus;
		dlstp(m3) = ap*fac_fus;

		% Update
		fff = Tt - T1t + lstn.*(qTt - qst) + lstp.*qpt;
		dfff = dlstn.*(qTt - qst) + dlstp.*qpt - lstn.*dqst - 1;
		dtabs = -fff./dfff;
		niter = niter + 1;
		T1t = T1t + dtabs;
        
        % Debug
        qnd(niter) = max(0, qTt(1) - qst(1));
        Td(niter) = T1t(1);

	end % while any(abs(dtabs) > 0.01 && niter < maxiter)

    % Fall back on bisection if needed
    if any(abs(dtabs)) > 0.01

        % Save previous temperatures
        m00 = logical(abs(dtabs) > 0.01);
        T1told = T1t(m00) - dtabs(m00);
        fffold = fff(m00);
		
        % Recalculate objective
        % Handle cloud condensate
		m1 = logical(T1t > SOM_T_0n);
		lstn(m1) = fac_cond;
		qst(m1) = SAM_qsatWater(T1t(m1), pt(m1));

		m2 = logical(T1t < SOM_T_00n);
		lstn(m2) = fac_sub;
		qst(m2) = SAM_qsatIce(T1t(m2), pt(m2));

		m3 = ~(m1 | m2);
		om = an*T1t(m3) - bn;
		lstn(m3) = fac_cond + (1-om)*fac_fus;
		qst(m3) = om.*SAM_qsatWater(T1t(m3), pt(m3)) + ...
			(1-om).*SAM_qsatIce(T1t(m3), pt(m3));

		% Handle precipitable condensate
		m1 = logical(T1t > SOM_T_0p);
		lstp(m1) = fac_cond;

		m2 = logical(T1t < SOM_T_00p);
		lstp(m2) = fac_sub;

		m3 = ~(m1 | m2);
		omp = ap*T1t(m3) - bp;
		lstp(m3) = fac_cond + (1-omp)*fac_fus;
        
		% Update objective
		fff = Tt - T1t + lstn.*(qTt - qst) + lstp.*qpt;
        
        % Debug
        qnd(niter+1) = max(0, qTt(1) - qst(1));
        Td(niter+1) = T1t(1);
        % END ITERATION

        % Extract unsolved elements
	    Ttt = Tt(m00);
	    T1tt = T1t(m00);
        fff = fff(m00);
	    ptt = pt(m00);
	    qTtt = qTt(m00);
	    qptt = qpt(m00);
	    qstt = qst(m00);
	    lsttn = zeros(size(qTtt));
	    lsttp = zeros(size(qTtt));

        % Perform bisection while temperatures are unconverged
        niter = 0;
        dtabst = dtabs(m00);
        while any(abs(dtabst) > 0.01) && niter < 10*maxiter

            % Calculate midpoint temperature
            T1tmid = (T1tt + T1told) ./ 2;

            % Calculate objective value at midpoint
            % Handle cloud condensate
		    m1 = logical(T1tmid > SOM_T_0n);
		    lsttn(m1) = fac_cond;
		    qstt(m1) = SAM_qsatWater(T1tmid(m1), ptt(m1));

		    m2 = logical(T1tmid < SOM_T_00n);
		    lsttn(m2) = fac_sub;
		    qstt(m2) = SAM_qsatIce(T1tmid(m2), ptt(m2));

		    m3 = ~(m1 | m2);
		    om = an*T1tmid(m3) - bn;
		    lsttn(m3) = fac_cond + (1-om)*fac_fus;
		    qstt(m3) = om.*SAM_qsatWater(T1tmid(m3), ptt(m3)) + ...
		    	(1-om).*SAM_qsatIce(T1tmid(m3), ptt(m3));

		    % Handle precipitable condensate
		    m1 = logical(T1tmid > SOM_T_0p);
		    lsttp(m1) = fac_cond;

		    m2 = logical(T1tmid < SOM_T_00p);
		    lsttp(m2) = fac_sub;

		    m3 = ~(m1 | m2);
		    omp = ap*T1tmid(m3) - bp;
		    lsttp(m3) = fac_cond + (1-omp)*fac_fus;

		    % Update
		    fffmid = Ttt - T1tmid + lsttn.*(qTtt - qstt) + lsttp.*qptt;
		    m1 = logical(fff.*fffmid > 0);
            T1told(~m1) = T1tmid(~m1);
            fffold(~m1) = fffmid(~m1);
            T1tt(m1) = T1tmid(m1);
            fff(m1) = fffmid(m1);
            dtabst = T1tt - T1told;
		    niter = niter + 1;
            
            % Debug
            qnd(maxiter+niter+1) = max(0, qTtt(1) - qstt(1));
            Td(maxiter+niter+1) = T1tmid(1);
        end

        % Merge bisection results
        T1t(m00) = T1tt;
        dtabs(m00) = dtabst;
        qst(m00) = qstt;

    end

    if any(abs(dtabs) > 0.01)
        warning('SOM_getState:nonconvergence',...
            'Failed to converge');
    end
	qst = qst + dqst .* dtabs;

	% Finalize all fields
	qn(m0) = max(0, qTt - qst);
	qn(~m0) = 0;
	T(m0) = T1t;

    % Return temperature and cloud condensate
	res = struct;
	res.T = T;
	res.qn = qn;

end
