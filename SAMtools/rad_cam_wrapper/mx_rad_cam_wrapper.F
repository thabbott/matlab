C=========================================================================
C   mx_rad_cam_wrapper
C   MATLAB wrapper for SAM's CAM-based radiative transfer scheme
C   This file written by Tristan Abbott // MIT EAPS // 8/2/2018
C   All others written by Marat Khairoutdinov and Peter Blossey
C=========================================================================
#include "fintrf.h"

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
        
        use grid
        use vars
        use params
        use rad, only: qrad
        implicit none

C       Arguments
        mwPointer plhs(*), prhs(*)
        integer nlhs, nrhs
        integer ip,ii,jj,kk,ll

C       Function declarations
        mwPointer mxGetPr
        mwPointer mxCreateNumericArray
        mwPointer mxGetNumberOfElements
        integer*4 mexPrintf
        integer*4 mxIsDouble
        integer*4 mxClassIDFromClassName

C       Local variables
        mwPointer tfr
        mwSize dummy_size
        integer, parameter :: sz = kind(dummy_size)
        double precision indp
        double precision, allocatable :: dparr(:)
        mwPointer nel
        mwSize, parameter :: ndim = 3
        mwSize dims(ndim)
        character(20) :: str

C       Check number of arguments
        if (.not. nrhs .eq. 24) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:nrhs",
     +          "24 arguments required")
        endif
        if (.not. nlhs .eq. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:nrhs",
     +          "1 output required")
        endif

C=========================================================================
C=========================================================================
C       Copy inputs and initialize required variables
C       MODULE GRID
        ip = mexPrintf("Transferring inputs\n")
C       zi
        ip = 1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nz) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:zi",
     +          "zi has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nz))
        call mxCopyPtrToReal8(tfr, dparr, nz)
        do ii = 1,nz
            zi(ii) = real(dparr(ii), 4)
        enddo
        deallocate(dparr)
C       presi
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nz) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:presi",
     +          "presi has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nz))
        call mxCopyPtrToReal8(tfr, dparr, nz)
        do ii = 1,nz
            presi(ii) = real(dparr(ii), 4)
        enddo
        deallocate(dparr)
C       pres0
        pres0 = presi(1)
C       day
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:day",
     +          "day has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        day = real(indp, 4)
C       day0
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:day0",
     +          "day0 has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        day0 = real(indp, 4)
C       dt
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:dt",
     +          "dt has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        dt = real(indp, 4)
C       nrad
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:nrad",
     +          "nrad has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        nrad = int(indp)
C       nstep
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:nstep",
     +          "nstep has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        nstep = int(indp)
C       icycle, nrestart, masterproc, dompi
        icycle = 1
        nrestart = 0
        masterproc = .true.
        dompi = .false.

C       MODULE PARAMS
C       doperpetual
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:doperpetual",
     +          "doperpetual has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            doperpetual = .false.
        else
            doperpetual = .true.
        endif
C       ocean
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:ocean",
     +          "ocean has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            ocean = .false.
        else
            ocean = .true.
        endif
C       dolongwave
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:dolongwave",
     +          "dolongwave has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            dolongwave = .false.
        else
            dolongwave = .true.
        endif
C       doshortwave
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:doshortwave",
     +          "doshortwave has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            doshortwave = .false.
        else
            doshortwave = .true.
        endif
C       doseasons
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:doseasons",
     +          "doseasons has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            doseasons = .false.
        else
            doseasons = .true.
        endif
C       dosolarconstant
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:dosolarconstant",
     +          "dosolarconstant has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            dosolarconstant = .false.
        else
            dosolarconstant = .true.
        endif
C       doradhomo
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
        call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:doradhomo",
     +          "doradhomo has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        if (indp < 0) then
            doradhomo = .false.
        else
            doradhomo = .true.
        endif
C       zenith_angle
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:zenith_angle",
     +          "zenith_angle has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        zenith_angle = real(indp,4)
C       solar_constant
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. 1) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:solar_constant",
     +          "solar_constant has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        call mxCopyPtrToReal8(tfr, indp, 1_sz)
        solar_constant = real(indp,4)

C       MODULE VARS
C       latitude
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nx*ny) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:latitude",
     +          "latitude has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny)
        ll = 1
        do jj = 1,ny
            do ii = 1,nx
                latitude(ii,jj) = real(dparr(ll), 4)
                ll = ll+1
            enddo
        enddo
        deallocate(dparr)
C       longitude
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nx*ny) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:longitude",
     +          "longitude has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny)
        ll = 1
        do jj = 1,ny
            do ii = 1,nx
                longitude(ii,jj) = real(dparr(ll), 4)
                ll = ll+1
            enddo
        enddo
        deallocate(dparr)
C       sstxy
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nx*ny) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:sstxy",
     +          "sstxy has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny)
        ll = 1
        do jj = 1,ny
            do ii = 1,nx
                sstxy(ii,jj) = real(dparr(ll), 4)
                ll = ll+1
            enddo
        enddo
        deallocate(dparr)
C       rho
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 .or. nel .ne. nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:rho",
     +          "rho has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nzm))
        call mxCopyPtrToReal8(tfr, dparr, nzm)
        do ii = 1,nzm
            rho(ii) = real(dparr(ii),4)
        enddo
        deallocate(dparr)
C       tabs
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 
     +      .or. nel .ne. nx*ny*nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:tabs",
     +          "tabs has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny*nzm))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny*nzm)
        ll = 1
        do kk = 1,nzm
            do jj = 1,ny
                do ii = 1,nx
                    tabs(ii,jj,kk) = real(dparr(ll),4)
                    ll = ll+1
                enddo
            enddo
        enddo
        deallocate(dparr)
C       qv
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 
     +      .or. nel .ne. nx*ny*nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:qv",
     +          "qv has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny*nzm))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny*nzm)
        ll = 1
        do kk = 1,nzm
            do jj = 1,ny
                do ii = 1,nx
                    qv(ii,jj,kk) = real(dparr(ll),4)
                    ll = ll+1
                enddo
            enddo
        enddo
        deallocate(dparr)
C       qcl
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 
     +      .or. nel .ne. nx*ny*nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:qcl",
     +          "qcl has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny*nzm))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny*nzm)
        ll = 1
        do kk = 1,nzm
            do jj = 1,ny
                do ii = 1,nx
                    qcl(ii,jj,kk) = real(dparr(ll),4)
                    ll = ll+1
                enddo
            enddo
        enddo
        deallocate(dparr)
C       qci
        ip = ip+1
        nel = mxGetNumberOfElements(prhs(ip))
        if (.not. mxIsDouble(prhs(ip)) .eq. 1 
     +       .or. nel .ne. nx*ny*nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:qci",
     +          "qci has wrong type or number of elements")
        endif
        tfr = mxGetPr(prhs(ip))
        allocate(dparr(nx*ny*nzm))
        call mxCopyPtrToReal8(tfr, dparr, nx*ny*nzm)
        ll = 1
        do kk = 1,nzm
            do jj = 1,ny
                do ii = 1,nx
                    qci(ii,jj,kk) = real(dparr(ll),4)
                    ll = ll+1
                enddo
            enddo
        enddo
        deallocate(dparr)

C=========================================================================
C=========================================================================
C       Run radiation scheme
        ip = mexPrintf("Calling rad_full\n")
        call rad_full()

C=========================================================================
C=========================================================================
C       Transfer outputs
        ip = mexPrintf("Transferring outputs\n")
        write(str, *) zi(2)
        ip = mexPrintf("qrad = "//str//"K/day\n")
        write(str, *) zi(1)
        ip = mexPrintf("qrad = "//str//"K/day\n")
C       qrad
        ip = 1
        dims(1) = nx
        dims(2) = ny
        dims(3) = nzm
        plhs(ip) = mxCreateNumericArray(ndim, dims,
     +       mxClassIDFromClassName('double'), 0)
        nel = mxGetNumberOfElements(plhs(ip))
        if (.not. mxIsDouble(plhs(ip)) .eq. 1 
     +       .or. nel .ne. nx*ny*nzm) then
            call mexErrMsgIdAndTxt("mx_rad_cam_wrapper:qrad",
     +          "Allocating qrad failed")
        endif
        allocate(dparr(nx*ny*nzm))
        ll = 1
        do kk = 1,nzm
            do jj = 1,ny
                do ii = 1,nx
                    dparr(ll) = real(qrad(ii,jj,kk), 8)
                    ll = ll+1
                enddo
            enddo
        enddo
        tfr = mxGetPr(plhs(ip))
        call mxCopyReal8ToPtr(dparr, tfr, nx*ny*nzm)
        deallocate(dparr)
        return

      end subroutine mexFunction
