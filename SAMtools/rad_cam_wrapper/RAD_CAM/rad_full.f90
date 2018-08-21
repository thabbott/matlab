subroutine rad_full()

  ! Interface to the longwave and shortwave radiation code from the
  ! NCAR Community Atmosphere Model (CAM3.0).
  !
  ! Originally written as interface to CCM3 radiation code by Marat
  !     Khairoutdinov
  ! Adapted to CAM3.0 radiation code by Peter Blossey, August 2004.
  !
  use rad
  use ppgrid
  use vars
  use params
  use shr_orb_mod, only: shr_orb_params
  use radae,        only: radaeini, initialize_radbuffer
  use pkg_cldoptics, only: cldefr, cldems
  use aer_optics, only: aer_optics_initialize
!tha  use microphysics, only: Get_reffc, Get_reffi

  implicit none

  integer*4 mexPrintf   !tha
  integer *4 mpf        !tha

  ! Local space:

  real(r4) pmid(pcols,pver)	! Level pressure (Pa)
  real(r4) pint(pcols,pverp)	! Model interface pressure (Pa)
  real(r4) massl(pcols,pver)	! Level mass (g/m2)
  real(r4) pmidrd(pcols,pver)	! Level pressure (dynes/cm2)
  real(r4) pintrd(pcols,pverp)	! Model interface pressure (dynes/cm2)
  real(r4) pmln(pcols,pver)	! Natural Log of pmid
  real(r4) piln(pcols,pverp)	! Natural Log of pint
  real(r4) tlayer(pcols,pver)	! Temperature
  real(r4) qlayer(pcols,pver)	! Specific humidity
  real(r4) cld(pcols,pverp)	! Fractional cloud cover
  real(r4) cliqwp(pcols,pver)	! Cloud liquid water path
  real(r4) cicewp(pcols,pver)	! Cloud ice water path
  real(r4) fice(pcols,pver)	! Fractional ice content within cloud
  real(r4) rel(pcols,pver)	! Liquid effective drop radius (micron)
  real(r4) rei(pcols,pver)	! Ice effective drop size
  real(r4) o3vmr(pcols,pver)	! Ozone volume mixing ratio
  real(r4) o3mmr(pcols,pver)		! Ozone mass mixing ratio

  integer lchnk              ! chunk identifier
  integer ncol               ! number of atmospheric columns
  integer nmxrgn(pcols)      ! Number of maximally overlapped regions

  real(r4) emis(pcols,pver)     ! cloud emissivity (fraction)
  real(r4) landfrac(pcols)      ! Land fraction (seems deprecated)
  real(r4) icefrac(pcols)       ! Ice fraction
  real(r4) psurface(pcols)      ! Surface pressure
  real(r4) player(pcols,pver)   ! Midpoint pressures
  real(r4) landm(pcols)         ! Land fraction
  real(r4) snowh(pcols)         ! snow depth, water equivalent (meters)

  real(r4) pmxrgn(pcols,pverp)  ! Maximum values of pmid for each

  real(r4) qrl(pcols,pver)	! Longwave heating rate (K/s)
  real(r4) qrs(pcols,pver)	! Shortwave heating rate (K/s)

  real(r4) fnl(pcols,pverp)	! Net Longwave Flux at interfaces
  real(r4) fns(pcols,pverp)	! Net Shortwave Flux at interfaces
  real(r4) fcnl(pcols,pverp)	! Net Clearsky Longwave Flux at interfaces
  real(r4) fcns(pcols,pverp)	! Net Clearsky Shortwave Flux at interfaces
  real(r4) flu(pcols,pverp)	! Longwave upward flux
  real(r4) fld(pcols,pverp)	! Longwave downward flux
  real(r4) fsu(pcols,pverp)	! Shortwave upward flux
  real(r4) fsd(pcols,pverp)	! Shortwave downward flux

  !	aerosols:

  real(r4) rh(pcols,pver)		! relative humidity for aerorsol 
  real(r4) aer_mass(pcols,pver,naer_all)  ! aerosol mass mixing ratio
  integer, parameter :: nspint = 19 ! # spctrl intrvls in solar spectrum
  integer, parameter :: naer_groups= 7 ! # aerosol grp for opt diagnostcs
  real(r4) aertau(nspint,naer_groups) ! Aerosol column optical depth
  real(r4) aerssa(nspint,naer_groups) ! Aero col-avg single scattering albedo
  real(r4) aerasm(nspint,naer_groups) ! Aerosol col-avg asymmetry parameter
  real(r4) aerfwd(nspint,naer_groups) ! Aerosol col-avg forward scattering

  !       Diagnostics:

  ! Longwave radiation
  real(r4) flns(pcols)          ! Surface cooling flux
  real(r4) flnt(pcols)          ! Net outgoing flux
  real(r4) flnsc(pcols)         ! Clear sky surface cooing
  real(r4) flntc(pcols)         ! Net clear sky outgoing flux
  real(r4) flwds(pcols)         ! Down longwave flux at surface

  !bloss: New in CAM3.0.
  real(r4) flut(pcols)          ! Upward flux at top of model
  real(r4) flutc(pcols)         ! Upward clear-sky flux at top of model

  ! Shortwave radiation
  real(r4) solin(pcols)        ! Incident solar flux
  real(r4) fsns(pcols)         ! Surface absorbed solar flux
  real(r4) fsnt(pcols)         ! Flux Shortwave Downwelling Top-of-Model
  real(r4) fsntoa(pcols)      ! Total column absorbed solar flux
  real(r4) fsds(pcols)         ! Flux Shortwave Downwelling Surface
  real(r4) fsdsc(pcols)        ! Clearsky Flux Shortwave Downwelling Surface

  real(r4) fsnsc(pcols)        ! Clear sky surface absorbed solar flux
  real(r4) fsntc(pcols)        ! Clear sky total column absorbed solar flx
  real(r4) fsntoac(pcols)      ! Clear sky total column absorbed solar flx
  real(r4) sols(pcols)         ! Direct solar rad incident on surface (< 0.7)
  real(r4) soll(pcols)         ! Direct solar rad incident on surface (>= 0.7)
  real(r4) solsd(pcols)        ! Diffuse solar rad incident on surface (< 0.7)
  real(r4) solld(pcols)        ! Diffuse solar rad incident on surface (>= 0.7)
  real(r4) fsnirtoa(pcols)     ! Near-IR flux absorbed at toa
  real(r4) fsnrtoac(pcols)     ! Clear sky near-IR flux absorbed at toa
  real(r4) fsnrtoaq(pcols)     ! Near-IR flux absorbed at toa >= 0.7 microns

  real(r4) frc_day(pcols)      ! = 1 for daylight, =0 for night columns
  real(r4) coszrs_in(pcols)    ! cosine of solar zenith angle

  real(r4) asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r4) aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r4) asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r4) aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms

  real(r4) lwupsfc(pcols)   ! Longwave up flux in CGS units

  real(r4) qtot
  real(r4) dayy
  integer i,j,k,m,ii,jj,i1,j1,tmp_count,nrad_call
  integer iday, iday0
  real(r4) coef,factor,tmp(1)
  real(8) qradz(nzm),buffer(nzm)
  real perpetual_factor
  real(r4) clat(pcols),clon(pcols)
  real(r4) pii
  real(r4) tmp_ggr, tmp_cp, tmp_eps, tmp_ste, tmp_pst
  
!tha  if(icycle.ne.1) goto 999  ! ugly way to handle the subcycles. add rad heating.

  nrad_call = 3600./dt
  pii = atan2(0.,-1.)

  ncol = 1 ! compute one column of radiation at a time.

  !-------------------------------------------------------
  ! Initialize some stuff
  !


  if(initrad) then

     mpf = mexPrintf("\tBRANCH TAKEN: initrad\n")

     ! check whether this processor's portion of the domain is
     ! evenly divisible into chunks of size ndiv (over which
     ! the absorbtivity/emissivity computations will be performed).
     if(mod(nx,ndiv).ne.0.or.(RUN3D.and.mod(ny,ndiv).ne.0)) then
        if(masterproc) print*,'nx or ny is not divisible by ndiv'
        if(masterproc) print*,'set in RAD_CAM/rad.f90'
        if(masterproc) print*,'Stop.'
        call abort()
     end if

     ! set up size and number of chunks of data for abs/ems computations
     nxdiv=max(1,nx/ndiv)
     nydiv=max(1,ny/ndiv)
     begchunk = 1
     endchunk = nxdiv*nydiv

     !bloss  subroutine initialize_radbuffer
     !bloss  inputs:  none
     !bloss  ouptuts: none (allocates and initializes abs/ems arrays)
     mpf = mexPrintf("\t\tcall initialize_radbuffer\n")
     call initialize_radbuffer()

     !bloss  subroutine shr_orb_params
     !bloss  inputs:  iyear, log_print
     !bloss  ouptuts: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
     mpf = mexPrintf("\t\tcall shr_orb_params\n")
     call shr_orb_params(iyear    , eccen  , obliq , mvelp     ,     &
           &               obliqr   , lambm0 , mvelpp, .false.)

     !bloss  subroutine radaeini
     !bloss  inputs:  pstdx (=1013250 dynes/cm2), mwdry (mwair) and mwco2.
     !bloss  ouptuts: none (sets up lookup tables for abs/ems computat.)
     mpf = mexPrintf("\t\tcall radaeini\n")
     call radaeini( 1.013250e6_r4, mwdry, mwco2 )

     !bloss  subroutine aer_optics_initialize
     !bloss  inputs:  none
     !bloss  ouptuts: none (sets up lookup tables for aerosol properties)
     mpf = mexPrintf("\t\tcall aer_optics_initialize\n")
     call aer_optics_initialize()

     ! sets up initial mixing ratios of trace gases.
     mpf = mexPrintf("\t\tcall tracesini\n")
     call tracesini()

     if(nrestart.eq.0) then

        mpf = mexPrintf("\tBRANCH TAKEN: nrestart .eq. 0\n")

        do k=1,nzm
           do j=1,ny
              do i=1,nx
	         tabs_rad(i,j,k)=0.
	         qv_rad(i,j,k)=0.
	         qc_rad(i,j,k)=0.
	         qi_rad(i,j,k)=0.
	         cld_rad(i,j,k)=0.
	         rel_rad(i,j,k)=25.
	         rei_rad(i,j,k)=25.
	         qrad(i,j,k)=0.
              end do
           end do
        end do
        nradsteps=0	  
        do k=1,nz
           radlwup(k) = 0.
           radlwdn(k) = 0.
           radswup(k) = 0.
           radswdn(k) = 0.
           radqrlw(k) = 0.
           radqrsw(k) = 0.
        end do

     else
        
        mpf = mexPrintf("\tBRANCH NOT TAKEN: nrestart .eq. 0\n")
!tha        call read_rad()

     endif

     if(doperpetual) then

           mpf = mexPrintf("\tBRANCH TAKEN: doperpetual\n")
           ! perpetual sun (no diurnal cycle)
           do j=1,ny
              do i=1,nx
                 p_factor(i,j) = perpetual_factor(day0, latitude(i,j)&
                      &,longitude(i,j))
              end do
           end do
     else
         mpf = mexPrintf("\tBRANCH NOT TAKEN: doperpetual\n")
     end if

  else
     mpf = mexPrintf("\tBRANCH NOT TAKEN: initrad\n")
  endif

  !bloss  subroutine radini
  !bloss  inputs:  ggr, cp, epislo (=0.622), stebol (=5.67e-8), pstd
  !bloss  outputs: none (although it initializes constants, computes
  !                       ozone path lengths).
  tmp_ggr = ggr
  tmp_cp  = cp
  tmp_eps = mwh2o/mwdry
  tmp_ste = 5.67e-8_r4
  tmp_pst = 1.013250e6_r4
  call radini(tmp_ggr, tmp_cp, tmp_eps, tmp_ste, tmp_pst)
  
  !bloss  initialize co2 mass mixing ratio
  co2vmr = 3.550e-4_r4 * nxco2
  co2mmr = co2vmr*rmwco2 ! rmwco2: ratio of mw of co2 to that of dry air
  if ((nstep.eq.1).and.(icycle.eq.1)) then
     if (masterproc) write(*,*) 'CO2 VMR = ', co2vmr
     if (masterproc) write(*,*) 'CO2 MMR = ', co2mmr
  end if

  ! compute pressure levels using pressure at the top
  ! interface and hydrostatic balance below that level.  This
  ! may be important because the radiation code assumes that
  ! the mass of a layer is equal to the difference between the
  ! pressure levels at the upper and lower interfaces.
  pint(:,1)=presi(nz)*100.
  piln(:,1) = log(pint(:,1))
  do k = nzm,1,-1
     m=nz-k
     pint(:,m+1) = pint(:,m) + rho(k)*ggr*(zi(k+1)-zi(k))
     piln(:,m+1) = log(pint(:,m+1))
  end do
  do m=1,nzm
     pmid(:,m) = 0.5*(pint(:,m)+pint(:,m+1))
     pmln(:,m) = log(pmid(:,m))
  end do

  do k=1,nzm
     massl(:,k)=1000.*(pint(:,k+1)-pint(:,k))/ggr
     rh(:,k)=0.
     o3vmr(:,k)=0.6034*o3(k)
     qrl(:,k)=0.
     qrs(:,k)=0.
  end do
  do k=1,nz
     cld(:,k)=0.
     flu(:,k)=0.
     fld(:,k)=0.
     fsu(:,k)=0.
     fsd(:,k)=0.
  end do

  ! Initialize aerosol mass mixing ratio to zero.
  ! TODO: come up with scheme to input aerosol concentrations 
  ! similar to the current scheme for trace gases.
  aer_mass = 0.

  !-----------------------------------------------------------
  ! Check if it is time to compute gas absortion coefficients for
  ! longwave radiation. This coefficients are computed for
  ! horizontally average fields and storred internally, so
  ! they need not to be recalculated on every call of radclw() for
  ! efficiency reasons.

!tha    if(initrad.or.mod(nstep,nrad_call).eq.0) then
    if (initrad) then

     mpf = mexPrintf("\tBRANCH TAKEN: initrad .or. mod(nstep, nrad_call)\n")

     initrad=.false.

     lchnk = 0
     do jj=1,nydiv 
        j1=(jj-1)*(ny/nydiv) 
        do ii=1,nxdiv 
	   i1=(ii-1)*(nx/nxdiv) 
           lchnk = lchnk + 1
	   do k=1,nzm
              tlayer(:,k)=0.
              qlayer(:,k)=0.
              cld(:,k) = 0.
              cliqwp(:,k) = 0.
              fice(:,k) = 0.
              m=nz-k

              tmp_count = 0
              do j=j1+1,j1+(ny/nydiv)
                 do i=i1+1,i1+(nx/nxdiv)	 
                    tlayer(1,k)=tlayer(1,k)+tabs(i,j,m)
                    qlayer(1,k)=qlayer(1,k)+qv(i,j,m)
                    tmp_count = tmp_count+1
                 end do
              end do
              tlayer(1,k)=tlayer(1,k)/float(tmp_count)
              qlayer(1,k)=max(1.e-7,qlayer(1,k)/float(tmp_count))

	   end do

           !bloss  subroutine radinp
           !bloss  inputs:  lchnk, ncol, pmid, pint, o3vmr
           !bloss  outputs: pmidrd,  ! Pressure at mid-levels (dynes/cm*2)
           !                pintrd,  ! Pressure at interfaces (dynes/cm*2)
           !                eccf,    ! Earth-sun distance factor
           !                o3mmr    ! Ozone mass mixing ratio
           call radinp(lchnk,ncol,pmid,pint,o3vmr,pmidrd, &
                pintrd,eccf,o3mmr)
           if ((nstep.eq.1).and.(icycle.eq.1).and.(lchnk.eq.1)) then
              if (masterproc) write(*,*) 'Eccentricity ', eccen
           end if

           lwupsfc(1) = stebol*(sstxy(1,1)+t00)**4 ! CGS units

           !bloss: Set flag for absorbtivity/emissivity computation
           doabsems = .true. 

           ! bloss: Set number of maximally overlapped regions and maximum 
           !        pressure so that only a single region will be computed.
           nmxrgn = 1
           pmxrgn = 1.2e6

           !bloss  subroutine radclwmx
           !bloss  inputs:  lchnk, ncol, nmxrgn, pmxrgn, lwupsfc,
           !                tlayer, qlayer, o3vmr, pmid, pint, pmln, piln,
           !                n2o, ch4, cfc11, cfc12, cld, emis, aer_mass
           !bloss  outputs: qrl,   ! Longwave heating rate
           !                flns,  ! Surface cooling flux
           !                flnt,  ! Net outgoing flux
           !                flut,  ! Upward flux at top of model
           !                flnsc, ! Clear sky surface cooing
           !                flntc, ! Net clear sky outgoing flux
           !                flutc, ! Upward clear-sky flux at top of model
           !                flwds, ! Down longwave flux at surface
           !                fcnl,  ! clear sky net flux at interfaces
           !                fnl    ! net flux at interfaces 
           call radclwmx(lchnk   ,ncol    ,                   &
                lwupsfc ,tlayer  ,qlayer  ,o3vmr   , &
                pmidrd  ,pintrd  ,pmln    ,piln    ,          &
                n2o     ,ch4     ,cfc11   ,cfc12   , &
                cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                flut    ,flutc   , &
                aer_mass,fnl     ,fcnl    ,flu     ,fld)
        end do
     end do

  else
      mpf = mexPrintf("\tBRANCH NOT TAKEN: initrad .or. mod(nstep,nrad_call) .eq. 0\n")
  endif


  !------------------------------------------------------
  !  Accumulate thermodynamical fields over nrad steps 
  !

  do k=1,nzm
     do j=1,ny
        do i=1,nx
           tabs_rad(i,j,k)=tabs_rad(i,j,k)+tabs(i,j,k)
           qv_rad(i,j,k)=qv_rad(i,j,k)+qv(i,j,k)
           qc_rad(i,j,k)=qc_rad(i,j,k)+qcl(i,j,k)
           qi_rad(i,j,k)=qi_rad(i,j,k)+qci(i,j,k)
           if(qcl(i,j,k)+qci(i,j,k).gt.0.) cld_rad(i,j,k) = cld_rad(i,j,k)+1.
        end do
     end do
  end do
  ! Accumulate effective radius by weighting it with mass
!tha  if(compute_reffc) then
!tha    rel_rad(1:nx,1:ny,1:nzm) = rel_rad(1:nx,1:ny,1:nzm) + Get_reffc() * qcl(1:nx,1:ny,1:nzm) 
!tha  end if
!tha  if(compute_reffi) then
!tha    rei_rad(1:nx,1:ny,1:nzm) = rei_rad(1:nx,1:ny,1:nzm) + Get_reffi() * qci(1:nx,1:ny,1:nzm) 
!tha  end if
  nradsteps=nradsteps+1

  !----------------------------------------------------
  ! Update radiation variables if the time is due
  !

  !kzm Oct.14, 03 changed .eq.nrad to .ge.nrad to handle the
  ! case when a smaller nrad is used in restart  
  if(nstep.eq.1.or.nradsteps.ge.nrad) then 

     mpf = mexPrintf("\tBRANCH TAKEN: nstep .eq. 1 .or. nradsteps .ge. nrad\n")
     ! Compute radiation fields for averaged thermodynamic fields


     coef=1./float(nradsteps)

     do k=1,nz
        radlwup(k) = 0.
        radlwdn(k) = 0.
        radswup(k) = 0.
        radswdn(k) = 0.
        radqrlw(k) = 0.
        radqrsw(k) = 0.
     end do

!tha     if(compute_reffc) then
!tha       do k=1,nzm
!tha          do j=1,ny
!tha             do i=1,nx
!tha                 rel_rad(i,j,k) = max(2.5,min(60.,rel_rad(i,j,k)/(1.e-8+qc_rad(i,j,k))))
!tha             end do
!tha          end do
!tha       end do
!tha     end if
!tha     if(compute_reffi) then
!tha       do k=1,nzm
!tha          do j=1,ny
!tha             do i=1,nx
!tha                 rei_rad(i,j,k) = max(5.,min(250.,rei_rad(i,j,k)/(1.e-8+qi_rad(i,j,k))))
!tha             end do
!tha          end do
!tha       end do
!tha     end if
     do k=1,nzm
        do j=1,ny
           do i=1,nx
	      tabs_rad(i,j,k)=tabs_rad(i,j,k)*coef
	      qv_rad(i,j,k)=qv_rad(i,j,k)*coef
	      qc_rad(i,j,k)=qc_rad(i,j,k)*coef
	      qi_rad(i,j,k)=qi_rad(i,j,k)*coef
	      cld_rad(i,j,k)=cld_rad(i,j,k)*coef
           end do
        end do
     end do

     lchnk = 0
     do j=1,ny
        jj = (j-1)/(ny/nydiv)+1
        do i=1,nx
	   ii = (i-1)/(nx/nxdiv)+1
           lchnk = ii + (jj-1)*nxdiv
           do k=1,nzm
              m=nz-k
              tlayer(1,k)=tabs_rad(i,j,m)
              qtot = qc_rad(i,j,m)+qi_rad(i,j,m)
              qlayer(1,k)=max(1.e-7_r4,qv_rad(i,j,m))
              if(qtot.gt.0.) then
	         cliqwp(1,k) = qtot*massl(1,k)
	         fice(1,k) = qi_rad(i,j,m)/qtot
	         cld(1,k) = min(0.99_r4,cld_rad(i,j,m))
              else
	         cld(1,k) = 0.
	         cliqwp(1,k) = 0.
	         fice(1,k) = 0.
              endif
           end do

        ! Default reff computation:
           !bloss  subroutine cldefr
           !bloss  inputs:  lchnk, ncol, landfrac, icefrac, pres0, 
           !                pmid, landm, icrfrac, snowh
           !bloss  outputs: rel, rei (liq/ice effective radii)
           player = pmid
           psurface = 100.*pres0
           icefrac = 0.
           snowh = 0.
           landm = 0.
           landfrac = 0.
           if (.not.OCEAN) then
              if (i+j.eq.2) mpf = mexPrintf("\tBRANCH TAKEN: .not. ocean\n")
              landfrac = 1.
              landm = 1.
           else
              if (i+j.eq.2) mpf = mexPrintf("\tBRANCH NOT TAKEN: .not. ocean\n")
           end if
           call cldefr(lchnk,ncol,landfrac,tlayer,rel,rei, & ! CAM3 interface
              psurface,player,landm,icefrac, snowh)
!tha           if(compute_reffc) then
!tha             rel(1,nzm:1:-1) = rel_rad(i,j,:)
!tha           else
             rel_rad(i,j,:) = rel(1,nzm:1:-1)
!tha           end if
!tha           if(compute_reffi) then
!tha             rei(1,nzm:1:-1) = rei_rad(i,j,:)
!tha           else
             rei_rad(i,j,:) = rei(1,nzm:1:-1)
!tha           end if

           !bloss  subroutine cldems
           !bloss  inputs:  lchnk, ncol, cliqwp, rei, fice
           !bloss  outputs: emis (cloud emissivity)
           call cldems(lchnk,ncol,cliqwp,fice,rei,emis)

           !bloss  subroutine radinp
           !bloss  inputs:  lchnk, ncol, pmid, pint, o3vmr
           !bloss  outputs: pmidrd, pintrd, eccf, o3mmr
           call radinp(lchnk,ncol,pmid,pint,o3vmr,pmidrd, &
                pintrd,eccf,o3mmr)

           lwupsfc(1) = stebol*(sstxy(i,j)+t00)**4 ! CGS units

           if(dolongwave) then

               if (i+j.eq.2) mpf = mexPrintf("\tBRANCH TAKEN: dolongwave\n")

              !bloss: Set flag for absorbtivity/emissivity computation
              doabsems = .false. 

              ! bloss: Set number of maximally overlapped regions 
              !        and maximum pressure so that only a single 
              !        region will be computed.
              nmxrgn = 1
              pmxrgn = 1.2e6

              !bloss  subroutine radclwmx
              !bloss  inputs:  lchnk, ncol, nmxrgn, pmxrgn, lwupsfc,
              !                tlayer, qlayer, o3vmr, pmid, pint, pmln, piln,
              !                n2o, ch4, cfc11, cfc12, cld, emis, aer_mass
              !bloss  outputs: qrl,   ! Longwave heating rate
              !                flns,  ! Surface cooling flux
              !                flnt,  ! Net outgoing flux
              !                flut,  ! Upward flux at top of model
              !                flnsc, ! Clear sky surface cooing
              !                flntc, ! Net clear sky outgoing flux
              !                flutc, ! Upward clear-sky flux at top of model
              !                flwds, ! Down longwave flux at surface
              !                fcnl,  ! clear sky net flux at interfaces
              !                fnl    ! net flux at interfaces 
              call radclwmx(lchnk   ,ncol    ,                   &
                   lwupsfc ,tlayer  ,qlayer  ,o3vmr   , &
                   pmidrd  ,pintrd  ,pmln    ,piln    ,          &
                   n2o     ,ch4     ,cfc11   ,cfc12   , &
                   cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                   flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                   flut    ,flutc   , &
                   aer_mass,fnl     ,fcnl    ,flu    ,fld)
              ! convert radiative heating from units of J/kg/s to K/s
              qrl = qrl/cp
              !
              ! change toa/surface fluxes from cgs to mks units
              !
              flnt     = 1.e-3*flnt
              flntc    = 1.e-3*flntc
              flns     = 1.e-3*flns
              flnsc    = 1.e-3*flnsc
              flwds    = 1.e-3*flwds
              flut     = 1.e-3*flut
              flutc    = 1.e-3*flutc
           else
               qrl = qrl*0.0
               if (i+j.eq.2) mpf = mexPrintf("\tBRANCH NOT TAKEN: dolongwave\n")
           endif

           if(doshortwave) then

              if (i+j.eq.2)mpf = mexPrintf("\tBRANCH TAKEN: doshortwave\n")
              if (doseasons) then
                 if (i+j.eq.2)mpf = mexPrintf("\tBRANCH TAKEN: doseasons\n")
                 ! The diurnal cycle of insolation will vary
                 ! according to time of year of the current day.
                 dayy = day
              else
                 if (i+j.eq.2)mpf = mexPrintf("\tBRANCH NOT TAKEN: doseasons\n")
                 ! The diurnal cycle of insolation from the calendar
                 ! day on which the simulation starts (day0) will be
                 ! repeated throughout the simulation.
                 iday0 = day0
                 iday = day
                 dayy = day-iday
                 dayy = iday0 + dayy
              end if
              if(doperpetual) then
                 if (i+j.eq.2)mpf = mexPrintf("\tBRANCH TAKEN: doperpetual\n")
                 if (dosolarconstant) then
                    if (i+j.eq.2) mpf = mexPrintf("\tBRANCH TAKEN: dosolarconstant\n")
                    ! fix solar constant and zenith angle as specified
                    ! in prm file.
                    coszrs_in(1) = cos(zenith_angle*pii/180.)
                    eccf = solar_constant/(1367.)
                 else
                    if (i+j.eq.2) mpf = mexPrintf("\tBRANCH NOT TAKEN: dosolarconstant\n")
                    ! perpetual sun (no diurnal cycle) - Modeled after Tompkins
                    coszrs_in(1) = 0.637 ! equivalent to zenith angle of 50.5 deg
                    eccf = p_factor(i,j)/coszrs_in(1) ! Adjst solar constant
                 end if
              else
                 if (i+j.eq.2) mpf = mexPrintf("\tBRANCH NOT TAKEN: doperpetual\n")
                 !bloss  subroutine zenith
                 !bloss  inputs:  dayy, latitude, longitude, ncol
                 !bloss  outputs: coszrs  ! Cosine solar zenith angle
                 clat(1) = pie*latitude(i,j)/180.
                 clon(1) = pie*longitude(i,j)/180.
                 coszrs_in(1) = coszrs
                 call zenith(dayy,clat,clon,coszrs_in,ncol)
              end if

	      coszrs = coszrs_in(1) ! needed for the isccp simulator

              !bloss  subroutine albedo
              !bloss  inputs: OCEAN (land/ocean flag), coszrs_in
              !bloss  outputs: 
              !     asdir  ! Srf alb for direct rad   0.2-0.7 micro-ms
              !     aldir  ! Srf alb for direct rad   0.7-5.0 micro-ms
              !     asdif  ! Srf alb for diffuse rad  0.2-0.7 micro-ms
              !     aldif  ! Srf alb for diffuse rad  0.7-5.0 micro-ms
              call albedo(1,1,OCEAN,coszrs_in,sstxy(1,1)+t00,asdir,aldir,asdif,aldif)

              ! bloss: Set number of maximally overlapped regions 
              !        and maximum pressure so that only a single 
              !        region will be computed.
              nmxrgn = 1
              pmxrgn = 1.2e6

              !bloss: compute separate cloud liquid & ice water paths
              cicewp = fice*cliqwp
              cliqwp = cliqwp - cicewp

              !bloss: set up day fraction.
              frc_day(1) = 0.
              if (coszrs_in(1).gt.0.) frc_day(1) = 1.

              !bloss  subroutine radcswmx
              !bloss  inputs:  
              !     lchnk             ! chunk identifier
              !     ncol              ! number of atmospheric columns
              !     pmid     ! Level pressure
              !     pint     ! Interface pressure
              !     qlayer   ! Specific humidity (h2o mass mix ratio)
              !     o3mmr    ! Ozone mass mixing ratio
              !     aer_mass   ! Aerosol mass mixing ratio
              !     rh       ! Relative humidity (fraction)
              !     cld      ! Fractional cloud cover
              !     cicewp   ! in-cloud cloud ice water path
              !     cliqwp   ! in-cloud cloud liquid water path
              !     rel      ! Liquid effective drop size (microns)
              !     rei      ! Ice effective drop size (microns)
              !     eccf     ! Eccentricity factor (1./earth-sun dist^2)
              !     coszrs_in! Cosine solar zenith angle
              !     asdir    ! 0.2-0.7 micro-meter srfc alb: direct rad
              !     aldir    ! 0.7-5.0 micro-meter srfc alb: direct rad
              !     asdif    ! 0.2-0.7 micro-meter srfc alb: diffuse rad
              !     aldif    ! 0.7-5.0 micro-meter srfc alb: diffuse rad
              !     scon     ! solar constant
              !bloss  in/outputs: 
              !     pmxrgn   ! Maximum values of pressure for each
              !              !    maximally overlapped region. 
              !     nmxrgn   ! Number of maximally overlapped regions
              !bloss  outputs: 
              !     solin     ! Incident solar flux
              !     qrs       ! Solar heating rate
              !     fsns      ! Surface absorbed solar flux
              !     fsnt      ! Total column absorbed solar flux
              !     fsntoa    ! Net solar flux at TOA
              !     fsds      ! Flux shortwave downwelling surface
              !     fsnsc     ! Clear sky surface absorbed solar flux
              !     fsdsc     ! Clear sky surface downwelling solar flux
              !     fsntc     ! Clear sky total column absorbed solar flx
              !     fsntoac   ! Clear sky net solar flx at TOA
              !     sols      ! Direct solar rad on surface (< 0.7)
              !     soll      ! Direct solar rad on surface (>= 0.7)
              !     solsd     ! Diffuse solar rad on surface (< 0.7)
              !     solld     ! Diffuse solar rad on surface (>= 0.7)
              !     fsnirtoa  ! Near-IR flux absorbed at toa
              !     fsnrtoac  ! Clear sky near-IR flux absorbed at toa
              !     fsnrtoaq  ! Net near-IR flux at toa >= 0.7 microns
              !     frc_day   ! = 1 for daylight, =0 for night columns
              !     aertau    ! Aerosol column optical depth
              !     aerssa    ! Aerosol column avg. single scattering albedo
              !     aerasm    ! Aerosol column averaged asymmetry parameter
              !     aerfwd    ! Aerosol column averaged forward scattering
              !     fns       ! net flux at interfaces
              !     fcns      ! net clear-sky flux at interfaces
              !     fsu       ! upward shortwave flux at interfaces
              !     fsd       ! downward shortwave flux at interfaces
              call radcswmx(lchnk   ,ncol    ,                   &
                   pintrd  ,pmid    ,qlayer  ,rh      ,o3mmr   , &
                   aer_mass  ,cld     ,cicewp  ,cliqwp  ,rel     , &
                   rei     ,eccf    ,coszrs_in,scon    ,solin   , &
                   asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                   pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                   fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                   fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                   solsd   ,solld   ,frc_day ,                   &
                   aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                   fcns    ,fsu     ,fsd     )
              ! convert radiative heating from units of J/kg/s to K/s
              qrs = qrs/cp
              !
              ! change toa/surface fluxes from cgs to mks units
              !
              fsnt     = 1.e-3*fsnt
              fsntc    = 1.e-3*fsntc
              fsntoa   = 1.e-3*fsntoa
              fsntoac  = 1.e-3*fsntoac
              fsnirtoa = 1.e-3*fsnirtoa
              fsnrtoac = 1.e-3*fsnrtoac
              fsnrtoaq = 1.e-3*fsnrtoaq
              fsns     = 1.e-3*fsns
              fsnsc    = 1.e-3*fsnsc
              fsds     = 1.e-3*fsds
              fsdsc    = 1.e-3*fsdsc

              solin    = 1.e-3*solin
           else
               if (i+j.eq.2) mpf = mexPrintf("\tBRANCH NOT TAKEN: doshortwave\n")
           endif

           do k=1,nzm
              m=nz-k
              qrad(i,j,m)=qrl(1,k)+qrs(1,k)
              radlwup(m)=radlwup(m)+flu(1,k)*1.e-3
              radlwdn(m)=radlwdn(m)+fld(1,k)*1.e-3
              radqrlw(m)=radqrlw(m)+qrl(1,k)
              radswup(m)=radswup(m)+fsu(1,k)*1.e-3
              radswdn(m)=radswdn(m)+fsd(1,k)*1.e-3
              radqrsw(m)=radqrsw(m)+qrs(1,k)
           enddo

           lwnsxy(i,j) = flns(1)
           swnsxy(i,j) = fsns(1)
           lwntxy(i,j) = flut(1)
           lwntmxy(i,j) = flnt(1)
           swntxy(i,j) = fsntoa(1)
           swntmxy(i,j) = fsnt(1)
           lwnscxy(i,j) = flnsc(1)
           swnscxy(i,j) = fsnsc(1)
           lwntcxy(i,j) = flntc(1)
           swntcxy(i,j) = fsntoac(1)
           swdsxy(i,j) = fsds(1)
           lwdsxy(i,j) = flwds(1)
           solinxy(i,j) = solin(1)

        end do
     end do

     tabs_rad(:,:,:)=0.
     qv_rad(:,:,:)=0.
     qc_rad(:,:,:)=0.
     qi_rad(:,:,:)=0.
     cld_rad(:,:,:)=0.
!      if(compute_reffc) then
!        rel_rad(:,:,:) = 0. 
!      end if
!      if(compute_reffi) then
!        rei_rad(:,:,:) = 0. 
!      end if
     nradsteps=0
     
     if(masterproc.and.doshortwave.and..not.doperpetual) &
          print*,'radiation: coszrs=',coszrs_in(1),' solin=',solin(1)
     if(masterproc.and.doshortwave.and.doperpetual) &
          print*,'radiation: perpetual sun, solin=',solin(1)
     if(masterproc.and..not.doshortwave) &
          print*,'longwave radiation is called'

   ! Homogenize radiation:

     if(doradhomo) then    
        mpf = mexPrintf("\tBRANCH TAKEN: doradhomo\n")
        factor = 1./dble(nx*ny)
        do k=1,nzm
          qradz(k) = 0.
           do j=1,ny
             do i=1,nx
              qradz(k) = qradz(k) + qrad(i,j,k)
             end do
           end do
           qradz(k) = qradz(k) * factor
           buffer(k) = qradz(k)
        end do

        factor = 1./float(nsubdomains)

        do k=1,nzm
           qradz(k)=buffer(k)*factor
           do j=1,ny
             do i=1,nx
               qrad(i,j,k) = qradz(k) 
             end do
           end do
        end do

     else
        mpf = mexPrintf("\tBRANCH NOT TAKEN: doradhomo\n")
     end if

  else
     mpf = mexPrintf("\tBRANCH NOT TAKEN: nstep .eq. 1 .or. nradsteps .ge. nrad\n")
  endif ! (nradsteps.eq.nrad) 

!tha !--------------------------------------------------------
!tha ! Prepare statistics:
!tha   
!tha   do j=1,ny
!tha      do i=1,nx
!tha         ! Net surface and toa fluxes
!tha         lwns_xy(i,j) = lwns_xy(i,j) + lwnsxy(i,j) 
!tha         swns_xy(i,j) = swns_xy(i,j) + swnsxy(i,j)
!tha         lwnt_xy(i,j) = lwnt_xy(i,j) + lwntxy(i,j) 
!tha         swnt_xy(i,j) = swnt_xy(i,j) + swntxy(i,j)
!tha         lwnsc_xy(i,j) = lwnsc_xy(i,j) + lwnscxy(i,j)
!tha         swnsc_xy(i,j) = swnsc_xy(i,j) + swnscxy(i,j)
!tha         lwntc_xy(i,j) = lwntc_xy(i,j) + lwntcxy(i,j)
!tha         swntc_xy(i,j) = swntc_xy(i,j) + swntcxy(i,j)
!tha         ! TOA Insolation
!tha         solin_xy(i,j) = solin_xy(i,j) + solinxy(i,j)
!tha      end do
!tha   end do
!tha !----------------------------------------------------------------
!tha   if(dostatisrad) then
!tha 
!tha      do j=1,ny
!tha         do i=1,nx
!tha            s_flns = s_flns + lwnsxy(i,j) 
!tha            s_flnsc = s_flnsc + lwnscxy(i,j) 
!tha            s_flnt = s_flnt + lwntmxy(i,j) 
!tha            s_flntoa = s_flntoa + lwntxy(i,j) 
!tha            s_flntoac = s_flntoac + lwntcxy(i,j) 
!tha            s_fsnt = s_fsnt + swntmxy(i,j) 
!tha            s_fsntoa = s_fsntoa + swntxy(i,j) 
!tha            s_fsntoac = s_fsntoac + swntcxy(i,j) 
!tha            s_fsns = s_fsns + swnsxy(i,j) 
!tha            s_fsnsc = s_fsnsc + swnscxy(i,j) 
!tha            s_fsds = s_fsds + swdsxy(i,j) 
!tha            s_flds = s_flds + lwdsxy(i,j) 
!tha            s_solin = s_solin + solinxy(i,j) 
!tha         end do
!tha      end do
!tha   end if
!tha !----------------------------------------------------------------
!tha !  Write the radiation-restart file:
!tha   
!tha   
!tha   if(mod(nstep,nstat*(1+nrestart_skip)).eq.0.or.nstep.eq.nstop.or.nelapse.eq.0) then
!tha   
!tha      call write_rad() ! write radiation restart file
!tha   
!tha   endif
!tha 
!tha !-------------------------------------------------------
!tha ! Update the temperature field:	
!tha 
!tha 999 continue
!tha 
!tha   do k=1,nzm
!tha      do j=1,ny
!tha         do i=1,nx
!tha 	   t(i,j,k)=t(i,j,k)+qrad(i,j,k)*dtn
!tha         end do
!tha      end do
!tha   end do

end subroutine rad_full


real function perpetual_factor(day, lat, lon)
use shr_kind_mod, only: r4 => shr_kind_r4
use grid, ONLY: dt, nrad
use ppgrid, only: eccf, pie
implicit none

!  estimate the factor to multiply the solar constant
!  so that the sun hanging perpetually right above
!  the head (zenith angle=0) would produce the same
!  total input the TOA as the sun subgect to diurnal cycle.
!  coded by Marat Khairoutdinov, 2004
!
! Input arguments
!
real day             ! Calendar day, without fraction
real lat                ! Current centered latitude (degrees)
real lon          ! Centered longitude (degrees)

! Local:
real(r4) :: tmp
real(r4) :: ttime
real(r4) :: coszrs(1) 
real(r4) :: clat(1), clon(1)
real :: dttime 

ttime = day
dttime = dt*float(nrad)/86400.
tmp = 0.
pie = acos(-1.)

do while (ttime.lt.day+1.)
   clat = pie*lat/180.
   clon = pie*lon/180.
  call zenith(ttime, clat, clon, coszrs, 1)
  tmp = tmp + min(dttime,day+1-ttime)*max(0._r4,eccf*coszrs(1))
  ttime = ttime+dttime
end do

perpetual_factor = tmp

end function perpetual_factor

