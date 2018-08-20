MEX = /usr/local/MATLAB/R2016b/bin/mex
INCLUDE = -I/data/thabbott/matlab/SAMtools/rad_cam_wrapper/netcdff-4.4.4/include
LIBC = -L/data/thabbott/matlab/SAMtools/rad_cam_wrapper/netcdfc-4.6.1/lib
LIBF = -L/data/thabbott/matlab/SAMtools/rad_cam_wrapper/netcdff-4.4.4/lib
NCCLIB = /data/thabbott/matlab/SAMtools/rad_cam_wrapper/netcdfc-4.6.1/lib
NCFLIB = /data/thabbott/matlab/SAMtools/rad_cam_wrapper/netcdff-4.4.4/lib

wrapper:
	ln -s RAD_CAM/*.f90 .
	$(MEX) $(INCLUDE) -c domain.f90
	$(MEX) $(INCLUDE) -c \
		grid.f90 \
		params.f90 \
		vars.f90
	$(MEX) $(INCLUDE) -c \
		shr_kind_mod.f90
	$(MEX) $(INCLUDE) -c \
		rad.f90 \
		ppgrid.f90 \
		shr_orb_mod.f90 \
		volcrad.f90 \
		abortutils.f90 \
		aer_optics.f90
	$(MEX) $(INCLUDE) -c \
		radae.f90 \
		pkg_cldoptics.f90 \
		cam_rad_parameterizations.f90
	$(MEX) $(INCLUDE) \
		$(LIBF) -lnetcdff $(LIBC) -lnetcdf \
		mx_rad_cam_wrapper.F \
		albedo.f90 \
		cldefrint.f90 \
		gffgch.f90 \
		quicksort.f90 \
		radclwmx.f90 \
		radcswmx.f90 \
		raddedmx.f90 \
		rad_full.f90 \
		radini.f90 \
		radinp.f90 \
		radoz2.f90 \
		traceini.f90 \
		trcab.f90 \
		trcabn.f90 \
		trcems.f90 \
		trcplk.f90 \
		trcpth.f90 \
		wrap_nf.f90 \
		zenith.f90 \
		misc.f90 \
		domain.o \
		grid.o \
		params.o \
		vars.o \
		shr_kind_mod.o \
		rad.o \
		ppgrid.o \
		shr_orb_mod.o \
		volcrad.o \
		abortutils.o \
		aer_optics.o \
		radae.o \
		pkg_cldoptics.o \
		cam_rad_parameterizations.o

clean:
	mv grid.f90 grid.f90.tmp
	mv params.f90 params.f90.tmp
	mv vars.f90 vars.f90.tmp
	mv domain.f90 domain.f90.tmp
	mv misc.f90 misc.f90.tmp
	-rm *.f90
	-rm *.o
	-rm *.mod
	-rm *.mexa64
	-rm *.so
	mv grid.f90.tmp grid.f90
	mv params.f90.tmp params.f90
	mv vars.f90.tmp vars.f90
	mv domain.f90.tmp domain.f90
	mv misc.f90.tmp misc.f90
