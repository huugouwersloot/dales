!> \file modsurface.f90
!!  Surface parameterization
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

!>
!! Surface routine including a land-surface scheme
!>
!! This module provides an interactive surface parameterization
!!
!! \par Revision list
!! \par Chiel van Heerwaarden
!! \todo documentation
!! \todo implement water reservoir at land surface for dew and interception water
!! \todo add moisture transport between soil layers
!! \deprecated Modsurface replaces the old modsurf.f90
!! \todo: implement fRs[t] based on the tiles.
!
!  Note that rootf should add up to 1 over all layers
!
!Able to handle heterogeneous surfaces using the switch lhetero
!In case of heterogeneity an input file is needed
!EXAMPLE of the old surface.inp.xxx for isurf = 3,4 (use switch loldtable for compatibility):

!#Surface input file - the standard land cover should be listed below. It is marked by typenr = 0
!#typenr    name       z0mav  z0hav    thls   ps    ustin  wtsurf  wqsurf  wsvsurf(01)  wsvsurf(02)  wsvsurf(03)  wsvsurf(04)
!    0   "standard  "  0.035  0.035   300.0  1.0e5   0.1    0.15   0.1e-3          1.0          0.0          0.0       0.0005
!    1   "forest    "  0.500  0.500   300.0  1.0e5   0.1    0.15   0.2e-3          1.0          0.0          0.0       0.0005
!    2   "grass     "  0.035  0.035   300.0  1.0e5   0.1    0.30   0.1e-3          1.0          0.0          0.0       0.0005

!EXAMPLE of surface.interactive.inp.xxx for isurf = 1:

!#Surface input file - the standard land cover should be listed below. It is marked by typenr = 0
!#typenr    name       z0mav  z0hav   ps    albedo  tsoil1 tsoil2 tsoil3 tsoil4 tsoildeep phiw1 phiw2 phiw3 phiw4 rootf1 rootf2 rootf3 rootf4 Cskin  lambdaskin  Qnet  cveg  Wlav   rsmin  LAI     gD    wsvsurf(01)  wsvsurf(02)  wsvsurf(03)  wsvsurf(04)
!    0   "standard  "  0.035  0.035  1.0e5   0.20      290    287    285    283       283   0.3   0.3   0.3   0.3   0.35   0.38   0.23   0.04 20000         5.0   450   0.9   0.0     110  4.0     0.            1.0          0.0          0.0       0.0005
!    1   "forest    "  0.500  0.500  1.0e5   0.15      290    287    285    283       283   0.3   0.3   0.3   0.3   0.35   0.38   0.23   0.04 20000         5.0   450   0.9   0.0     110  6.0     0.            1.0          0.0          0.0       0.0005
!    2   "grass     "  0.035  0.035  1.0e5   0.25      290    287    285    283       283   0.3   0.3   0.3   0.3   0.35   0.38   0.23   0.04 20000         5.0   450   0.9   0.0     110  2.0     0.            1.0          0.0          0.0       0.0005

!EXAMPLE of surface.prescribed.inp.xxx for isurf = 2,3,4:

!#Surface input file - the standard land cover should be listed below. It is marked by typenr = 0
!#typenr    name       z0mav  z0hav    thls   ps    albedo  rsi_s2  ustin  wtsurf  wqsurf  wsvsurf(01)  wsvsurf(02)  wsvsurf(03)  wsvsurf(04)
!    0   "standard  "  0.035  0.035   300.0  1.0e5   0.20     50.0   0.1    0.15   0.1e-3          1.0          0.0          0.0       0.0005
!    1   "forest    "  0.500  0.500   300.0  1.0e5   0.15     50.0   0.1    0.15   0.2e-3          1.0          0.0          0.0       0.0005
!    2   "grass     "  0.035  0.035   300.0  1.0e5   0.25     50.0   0.1    0.30   0.1e-3          1.0          0.0          0.0       0.0005


module modsurface
  use modtimer
  use modsurfdata
  implicit none
  !public  :: initsurface, surface, exitsurface

save

contains
!> Reads the namelists and initialises the soil.
  subroutine initsurface
    use modglobal,  only : i1, j1, i2, j2, itot, jtot,imax,jmax, nsv, ifnamopt, fname_options, ifinput, cexpnr, checknamelisterror, handle_err
    use modraddata, only : iradiation,rad_shortw,irad_par,irad_user,irad_rrtmg,irad_rte_rrtmgp
    use modmpi,     only : myid,  myidx, myidy, comm3d, mpierr, D_MPI_BCAST
    use netcdf

    implicit none

    integer   :: i,j,k, landindex, ierr, defined_landtypes, landtype_0 = -1
    integer   :: tempx,tempy
    integer   :: VARID,STATUS,NCID,timeID
    character(len = nf90_max_name) :: RecordDimName
    character(len=1500) :: readbuffer

    namelist/NAMSURFACE/ & !< Soil related variables
      isurf, tsoilav, tsoildeepav, phiwav, rootfav, &
      ! Land surface related variables
      lmostlocal, lsmoothflux, lneutral, z0mav, z0hav, rsisurf2, Cskinav, lambdaskinav, albedoav, Qnetav, cvegav, Wlav, &
      ! Jarvis-Steward related variables
      rsminav, rssoilminav, LAIav, gDav, &
      ! Prescribed values for isurf 2, 3, 4
      z0, thls, dsst, ps, ustin, wtsurf, wqsurf, wsvsurf, &
      ! Heterogeneous variables
      lhetero, xpatches, ypatches, land_use, loldtable, &
      ! AGS variables
      lrsAgs, lCO2Ags,planttype, &
      ! Delay plant response in Ags
      lrelaxgc, kgc, lrelaxci, kci, &
      ! Soil properties
      phi, phifc, phiwp, R10, T2gm, Q10gm, &
      !2leaf AGS, sunlit/shaded
      lsplitleaf, &
      ! Exponential emission function
      i_expemis, expemis0, expemis1, expemis2, &
      ! heterogeneous tskin
      ltskininp,  min_horv

    call timer_tic('modsurface/initsurface', 0)

    ! 1    -   Initialize soil

    !if (isurf == 1) then

    ! 1.0  -   Read LSM-specific namelist

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSURFACE,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSURFACE')
      write(6 ,NAMSURFACE)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(isurf        , 1       ,  0, comm3d, mpierr)
    call D_MPI_BCAST(tsoilav      , ksoilmax, 0, comm3d, mpierr)
    call D_MPI_BCAST(tsoildeepav  , 1       , 0, comm3d, mpierr)
    call D_MPI_BCAST(phiwav       , ksoilmax, 0, comm3d, mpierr)
    call D_MPI_BCAST(rootfav      , ksoilmax, 0, comm3d, mpierr)

    call D_MPI_BCAST(lmostlocal   , 1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lsmoothflux  , 1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lneutral     , 1,  0, comm3d, mpierr)
    call D_MPI_BCAST(z0mav        , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(z0hav        , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(rsisurf2     , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(Cskinav      , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lambdaskinav , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(albedoav     , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(Qnetav       , 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(rsminav      , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(rssoilminav  , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(cvegav       , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(Wlav         , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(LAIav        , 1, 0, comm3d, mpierr)
    call D_MPI_BCAST(gDav         , 1, 0, comm3d, mpierr)

    call D_MPI_BCAST(z0         ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(ustin      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(wtsurf     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(wqsurf     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(wsvsurf(1:nsv),nsv,0,comm3d,mpierr)
    call D_MPI_BCAST(ps         ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(thls       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dsst       ,1,0,comm3d,mpierr)

    call D_MPI_BCAST(lhetero                    ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(loldtable                  ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lrsAgs                     ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lCO2Ags                    ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(xpatches                   ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(ypatches                   ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(planttype                  ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lrelaxgc                   ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(lrelaxci                   ,            1,  0, comm3d, mpierr)
    call D_MPI_BCAST(kgc                        ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(kci                        ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(phi                        ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(phifc                      ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(phiwp                      ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(R10                        ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(lsplitleaf                 ,            1,  0, comm3d, mpierr)

    call D_MPI_BCAST(land_use(1:mpatch,1:mpatch),mpatch*mpatch,  0, comm3d, mpierr)

    call D_MPI_BCAST(i_expemis                  ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(expemis0                   ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(expemis1                   ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(expemis2                   ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(ltskininp                  ,            1, 0, comm3d, mpierr)
    call D_MPI_BCAST(min_horv                   ,            1, 0, comm3d, mpierr)

    if(lCO2Ags .and. (.not. lrsAgs)) then
      if(myid==0) print *,"WARNING::: You set lCO2Ags to .true., but lrsAgs to .false."
      if(myid==0) print *,"WARNING::: Since AGS does not run, lCO2Ags will be set to .false. as well."
      lCO2Ags = .false.
    endif
    if(lsplitleaf .and. (.not. (rad_shortw .and. ((iradiation.eq.irad_par).or.(iradiation .eq. irad_user) &
                                                                          .or.(iradiation .eq. irad_rrtmg) &
                                                                          .or.(iradiation .eq. irad_rte_rrtmgp))))) then
      if(myid==0) stop "WARNING::: You set lsplitleaf to .true., but that needs direct and diffuse calculations. Make sure you enable rad_shortw"
      if(myid==0) stop "WARNING::: Since there is no direct and diffuse radiation calculated in the atmopshere, we set lsplitleaf to .false."
      lsplitleaf = .false.
    endif

    if(lrsAgs) then
      select case (planttype)
        case (3)              !<  C3 plants based on standard settings; R10, T2gm and Q10gm can be altered in the namelist
        case (4)              !<  C4 plants with standard settings; relevant standard settings for C3 plants are replaced; R10, T2gm and Q10gm cannot be altered in the namelist
          CO2comp298 =    4.3 !<    CO2 compensation concentration
          Q10CO2     =    1.5 !<    Parameter to calculate the CO2 compensation concentration
          gm298      =   17.5 !<    Mesophyll conductance at 298 K
          Q10gm      =    2.0 !<    Parameter to calculate the mesophyll conductance
          T1gm       =  286.0 !<    Reference temperature to calculate the mesophyll conductance
          T2gm       =  309.0 !<    Reference temperature to calculate the mesophyll conductance
          f0         =   0.85 !<    Maximum value Cfrac
          ad         =   0.15 !<    Regression coefficient to calculate Cfrac
          Ammax298   =    1.7 !<    CO2 maximal primary productivity
          Q10am      =    2.0 !<    Parameter to calculate maximal primary productivity
          T1Am       =    286 !<    Reference temperature to calculate maximal primary productivity
          T2Am       =    311 !<    Reference temperature to calculate maximal primary productivity
          alpha0     =  0.014 !<    Initial low light conditions
        case (5)              !<  C4 plants with altered settings; relevant standard settings for C3 plants are replaced, except for T2gm and Q10gm; T2gm and Q10gm should be altered in the namelist; R10 can be altered in the namelist
          CO2comp298 =    4.3 !<    CO2 compensation concentration
          Q10CO2     =    1.5 !<    Parameter to calculate the CO2 compensation concentration
          gm298      =   17.5 !<    Mesophyll conductance at 298 K
          T1gm       =  286.0 !<    Reference temperature to calculate the mesophyll conductance
          f0         =   0.85 !<    Maximum value Cfrac
          ad         =   0.15 !<    Regression coefficient to calculate Cfrac
          Ammax298   =    1.7 !<    CO2 maximal primary productivity
          Q10am      =    2.0 !<    Parameter to calculate maximal primary productivity
          T1Am       =    286 !<    Reference temperature to calculate maximal primary productivity
          T2Am       =    311 !<    Reference temperature to calculate maximal primary productivity
          alpha0     =  0.014 !<    Initial low light conditions
        case default          !<  Plant type is set to a scheme that is not (yet) supported
          if(myid==0) print *,"WARNING::: planttype should be either 3 or 4, corresponding to C3 or C4 plants. It now defaulted to 3."
          planttype = 3
      end select
    endif

    if(lhetero) then

      if(xpatches .gt. mpatch) then
        stop "NAMSURFACE: more xpatches defined than possible (change mpatch in modsurfdata to a higher value)"
      endif
      if(ypatches .gt. mpatch) then
        stop "NAMSURFACE: more ypatches defined than possible (change mpatch in modsurfdata to a higher value)"
      endif
      if (lsmoothflux .eqv. .true.) write(6,*) 'WARNING: You selected to use uniform heat fluxes (lsmoothflux) and ',&
      'heterogeneous surface conditions (lhetero) at the same time'
      if (mod(itot,xpatches) .ne. 0) stop "NAMSURFACE: Not an integer amount of grid points per patch in the x-direction"
      if (mod(jtot,ypatches) .ne. 0) stop "NAMSURFACE: Not an integer amount of grid points per patch in the y-direction"

      allocate(horvpatch(xpatches,ypatches))
      allocate(z0mav_patch(xpatches,ypatches))
      allocate(z0hav_patch(xpatches,ypatches))
      allocate(thls_patch(xpatches,ypatches))
      allocate(qts_patch(xpatches,ypatches))
      allocate(thvs_patch(xpatches,ypatches))
      allocate(ps_patch(xpatches,ypatches))
      allocate(ustin_patch(xpatches,ypatches))
      allocate(wt_patch(xpatches,ypatches))
      allocate(wq_patch(xpatches,ypatches))
      allocate(wsv_patch(100,xpatches,ypatches))
      allocate(rsisurf2_patch(xpatches,ypatches))
      allocate(albedo_patch(xpatches,ypatches))

      allocate(tsoil_patch(ksoilmax,xpatches,ypatches))
      allocate(tsoildeep_patch(xpatches,ypatches))
      allocate(phiw_patch(ksoilmax,xpatches,ypatches))
      allocate(rootf_patch(ksoilmax,xpatches,ypatches))
      allocate(Cskin_patch(xpatches,ypatches))
      allocate(lambdaskin_patch(xpatches,ypatches))
      allocate(Qnet_patch(xpatches,ypatches))
      allocate(cveg_patch(xpatches,ypatches))
      allocate(Wl_patch(xpatches,ypatches))
      allocate(rsmin_patch(xpatches,ypatches))
      allocate(LAI_patch(xpatches,ypatches))
      allocate(gD_patch(xpatches,ypatches))

      allocate(oblpatch(xpatches,ypatches))

      z0mav_patch = -1
      z0hav_patch = -1
      thls_patch  = -1
      qts_patch   = -1
      thvs_patch  = -1
      ps_patch    = -1
      ustin_patch = -1
      wt_patch    = -1
      wq_patch    = -1
      wsv_patch   = -1
      rsisurf2_patch = 0
      albedo_patch= -1

      tsoil_patch      = -1
      tsoildeep_patch  = -1
      phiw_patch       = -1
      rootf_patch      = -1
      Cskin_patch      = -1
      lambdaskin_patch = -1
      Qnet_patch       = -1
      cveg_patch       = -1
      Wl_patch         = -1
      rsmin_patch      = -1
      LAI_patch        = -1
      gD_patch         = -1

      oblpatch         = -0.1

      defined_landtypes = 0
      if(loldtable) then !Old input-file for heterogeneous surfaces: only valid w/o sw-radiation (due to albedo) and isurf = 3,4
        open (ifinput,file='surface.inp.'//cexpnr)
        ierr = 0
        do while (ierr == 0)
          read(ifinput, '(A)', iostat=ierr) readbuffer
          if (ierr == 0) then                               !So no end of file is encountered
            if (readbuffer(1:1)=='#') then
              if (myid == 0)   print *,trim(readbuffer)
            else
              if (myid == 0)   print *,trim(readbuffer)
              defined_landtypes = defined_landtypes + 1
              i = defined_landtypes
              read(readbuffer, *, iostat=ierr) landtype(i), landname(i), z0mav_land(i), z0hav_land(i), thls_land(i), &
                ps_land(i), ustin_land(i), wt_land(i), wq_land(i), wsv_land(1:nsv,i)

              if (ustin_land(i) .lt. 0) then
                if (myid == 0) stop "NAMSURFACE: A ustin value in the surface input file is negative"
              endif
              if(isurf .ne. 3) then
                if(z0mav_land(i) .lt. 0) then
                  if (myid == 0) stop "NAMSURFACE: a z0mav value is not set or negative in the surface input file"
                end if
                if(z0hav_land(i) .lt. 0) then
                  if (myid == 0) stop "NAMSURFACE: a z0hav value is not set or negative in the surface input file"
                end if
              end if

              if (landtype(i) .eq. 0) landtype_0 = i
              do j = 1, (i-1)
                if (landtype(i) .eq. landtype(j)) stop "NAMSURFACE: Two land types have the same type number"
              enddo

            endif
          endif
        enddo
        close(ifinput)
      else
        select case (isurf)
          case (1) ! Interactive land surface
            open (ifinput,file='surface.interactive.inp.'//cexpnr)
            ierr = 0
            do while (ierr == 0)
              read(ifinput, '(A)', iostat=ierr) readbuffer
              if (ierr == 0) then                               !So no end of file is encountered
                if (readbuffer(1:1)=='#') then
                  if (myid == 0)   print *,trim(readbuffer)
                else
                  if (myid == 0)   print *,trim(readbuffer)
                  defined_landtypes = defined_landtypes + 1
                  i = defined_landtypes
                  read(readbuffer, *, iostat=ierr) landtype(i), landname(i), z0mav_land(i), z0hav_land(i), ps_land(i), &
                    albedo_land(i),tsoil_land(1:ksoilmax,i),tsoildeep_land(i),phiw_land(1:ksoilmax,i),rootf_land(1:ksoilmax,i),&
                    Cskin_land(i), lambdaskin_land(i), Qnet_land(i), cveg_land(i), Wl_land(i), rsmin_land(i), LAI_land(i), &
                    gD_land(i), wsv_land(1:nsv,i)

                  if(z0mav_land(i) .lt. 0) then
                    if (myid == 0) stop "NAMSURFACE: a z0mav value is not set or negative in the surface input file"
                  end if
                  if(z0hav_land(i) .lt. 0) then
                    if (myid == 0) stop "NAMSURFACE: a z0hav value is not set or negative in the surface input file"
                  end if
                  if (albedo_land(i) .lt. 0) then
                    if (myid == 0) stop "NAMSURFACE: An albedo value in the surface input file is negative"
                  endif
                  if (albedo_land(i) .gt. 1) then
                    if (myid == 0) stop "NAMSURFACE: An albedo value in the surface input file is greater than 1"
                  endif

                  if (landtype(i) .eq. 0) landtype_0 = i
                  do j = 1, (i-1)
                    if (landtype(i) .eq. landtype(j)) stop "NAMSURFACE: Two land types have the same type number"
                  enddo

                endif
              endif
            enddo
            close(ifinput)

          case default ! Prescribed land surfaces: isurf = 2, 3, 4 (& 10)
            open (ifinput,file='surface.prescribed.inp.'//cexpnr)
            ierr = 0
            do while (ierr == 0)
              read(ifinput, '(A)', iostat=ierr) readbuffer
              if (ierr == -1) then
                if (myid == 0)  print *, "iostat = ",ierr,": No file ",'surface.prescribed.inp.'//cexpnr," found"
              endif
              if (ierr == 0) then                               !So no end of file is encountered
                if (readbuffer(1:1)=='#') then
                  if (myid == 0)   print *,trim(readbuffer)
                else
                  if (myid == 0)   print *,trim(readbuffer)
                  defined_landtypes = defined_landtypes + 1
                  i = defined_landtypes
                  read(readbuffer, *, iostat=ierr) landtype(i), landname(i), z0mav_land(i), z0hav_land(i), thls_land(i), &
                    ps_land(i), albedo_land(i), rsisurf2_land(i), ustin_land(i), wt_land(i), wq_land(i), wsv_land(1:nsv,i)

                  if (ustin_land(i) .lt. 0) then
                    if (myid == 0) stop "NAMSURFACE: A ustin value in the surface input file is negative"
                  endif
                  if (albedo_land(i) .lt. 0) then
                    if (myid == 0) stop "NAMSURFACE: An albedo value in the surface input file is negative"
                  endif
                  if (albedo_land(i) .gt. 1) then
                    if (myid == 0) stop "NAMSURFACE: An albedo value in the surface input file is greater than 1"
                  endif
                  if(isurf .ne. 3) then
                    if(z0mav_land(i) .lt. 0) then
                      if (myid == 0) stop "NAMSURFACE: a z0mav value is not set or negative in the surface input file"
                    end if
                    if(z0hav_land(i) .lt. 0) then
                      if (myid == 0) stop "NAMSURFACE: a z0hav value is not set or negative in the surface input file"
                    end if
                  end if

                  if (landtype(i) .eq. 0) landtype_0 = i
                  do j = 1, (i-1)
                    if (landtype(i) .eq. landtype(j)) stop "NAMSURFACE: Two land types have the same type number"
                  enddo

                endif
              endif
            enddo
            close(ifinput)
        end select
      endif

      if (myid == 0) then
        if (landtype_0 .eq. -1) then
          stop "NAMSURFACE: no standard land type (0) is defined"
        else
          print "(a,i2,a,i2)","There are ",defined_landtypes,&
          " land types defined in the surface input file. The standard land type is defined by line ",landtype_0
        endif
      endif

      select case (isurf)
        case (1) ! Interactive land surface
          tsoilav      = 0
          tsoildeepav  = 0
          phiwav       = 0
          rootfav      = 0
          Cskinav      = 0
          lambdaskinav = 0
          albedoav     = 0
          Qnetav       = 0
          cvegav       = 0
          rsminav      = 0
          LAIav        = 0
          gDav         = 0
          Wlav         = 0

          z0mav        = 0
          z0hav        = 0
          ps           = 0
          thls         = 300

          do i = 1, xpatches
            do j = 1, ypatches
              landindex = landtype_0
              do k = 1, defined_landtypes
                if (landtype(k) .eq. land_use(i,j)) then
                  landindex = k
                endif
              enddo

              tsoil_patch(:,i,j)    = tsoil_land(:,landindex)
              tsoildeep_patch(i,j)  = tsoildeep_land(landindex)
              phiw_patch(:,i,j)     = phiw_land(:,landindex)
              rootf_patch(:,i,j)    = rootf_land(:,landindex)
              Cskin_patch(i,j)      = Cskin_land(landindex)
              lambdaskin_patch(i,j) = lambdaskin_land(landindex)
              albedo_patch(i,j)     = albedo_land(landindex)
              Qnet_patch(i,j)       = Qnet_land(landindex)
              cveg_patch(i,j)       = cveg_land(landindex)
              rsmin_patch(i,j)      = rsmin_land(landindex)
              LAI_patch(i,j)        = LAI_land(landindex)
              gD_patch(i,j)         = gD_land(landindex)
              Wl_patch(i,j)         = Wl_land(landindex)

              z0mav_patch(i,j)     = z0mav_land(landindex)
              z0hav_patch(i,j)     = z0hav_land(landindex)
              ps_patch(i,j)        = ps_land(landindex)

              wsv_patch(1:nsv,i,j) = wsv_land(1:nsv,landindex)

              tsoilav(:)    = tsoilav(:)   +  ( tsoil_patch(:,i,j)    / ( xpatches * ypatches ) )
              tsoildeepav   = tsoildeepav  +  ( tsoildeep_patch(i,j)  / ( xpatches * ypatches ) )
              phiwav(:)     = phiwav(:)    +  ( phiw_patch(:,i,j)     / ( xpatches * ypatches ) )
              rootfav(:)    = rootfav(:)   +  ( rootf_patch(:,i,j)    / ( xpatches * ypatches ) )
              Cskinav       = Cskinav      +  ( Cskin_patch(i,j)      / ( xpatches * ypatches ) )
              lambdaskinav  = lambdaskinav +  ( lambdaskin_patch(i,j) / ( xpatches * ypatches ) )
              albedoav      = albedoav     +  ( albedo_patch(i,j)     / ( xpatches * ypatches ) )
              Qnetav        = Qnetav       +  ( Qnet_patch(i,j)       / ( xpatches * ypatches ) )
              cvegav        = cvegav       +  ( cveg_patch(i,j)       / ( xpatches * ypatches ) )
              rsminav       = rsminav      +  ( rsmin_patch(i,j)      / ( xpatches * ypatches ) )
              LAIav         = LAIav        +  ( LAI_patch(i,j)        / ( xpatches * ypatches ) )
              gDav          = gDav         +  ( gD_patch(i,j)         / ( xpatches * ypatches ) )
              Wlav          = Wlav         +  ( Wl_patch(i,j)         / ( xpatches * ypatches ) )

              z0mav         = z0mav        +  ( z0mav_patch(i,j)      / ( xpatches * ypatches ) )
              z0hav         = z0hav        +  ( z0hav_patch(i,j)      / ( xpatches * ypatches ) )
              ps            = ps           +  ( ps_patch(i,j)         / ( xpatches * ypatches ) )
            enddo
          enddo
        case default ! Prescribed land surfaces: isurf = 2, 3, 4 (& 10)
          thls   = 0
          ps     = 0
          ustin  = 0
          wtsurf = 0
          wqsurf = 0
          wsvsurf(1:nsv) = 0
          if (.not. loldtable) then
            albedoav  = 0
          endif

          z0mav        = 0
          z0hav        = 0

          do i = 1, xpatches
            do j = 1, ypatches
              landindex = landtype_0
              do k = 1, defined_landtypes
                if (landtype(k) .eq. land_use(i,j)) then
                  landindex = k
                endif
              enddo

              z0mav_patch(i,j)     = z0mav_land(landindex)
              z0hav_patch(i,j)     = z0hav_land(landindex)
              thls_patch(i,j)      = thls_land(landindex)
              ps_patch(i,j)        = ps_land(landindex)
              ustin_patch(i,j)     = ustin_land(landindex)
              wt_patch(i,j)        = wt_land(landindex)
              wq_patch(i,j)        = wq_land(landindex)
              wsv_patch(1:nsv,i,j) = wsv_land(1:nsv,landindex)
              if (.not. loldtable) then
                albedo_patch(i,j)  = albedo_land(landindex)
                rsisurf2_patch(i,j)= rsisurf2_land(landindex)
              endif

              thls   = thls   + ( thls_patch(i,j)  / ( xpatches * ypatches ) )
              ps     = ps     + ( ps_patch(i,j)    / ( xpatches * ypatches ) )
              ustin  = ustin  + ( ustin_patch(i,j) / ( xpatches * ypatches ) )
              wtsurf = wtsurf + ( wt_patch(i,j)    / ( xpatches * ypatches ) )
              wqsurf = wqsurf + ( wq_patch(i,j)    / ( xpatches * ypatches ) )
              wsvsurf(1:nsv) = wsvsurf(1:nsv) + ( wsv_patch(1:nsv,i,j) / ( xpatches * ypatches ) )
              if (.not. loldtable) then
                albedoav  = albedoav + ( albedo_patch(i,j) / ( xpatches * ypatches ) )
              endif

              z0mav  = z0mav  + ( z0mav_patch(i,j) / ( xpatches * ypatches ) )
              z0hav  = z0hav  + ( z0hav_patch(i,j) / ( xpatches * ypatches ) )
            enddo
          enddo
      end select
    else  ! not lhetero:
      if((z0mav == -1 .and. z0hav == -1) .and. (z0 .ne. -1)) then
        z0mav = z0
        z0hav = z0
        write(6,*) "WARNING: z0m and z0h not defined, set equal to z0"
      end if

      if(isurf .ne. 3) then
        if(z0mav == -1) then
          stop "NAMSURFACE: z0mav is not set"
        end if
        if(z0hav == -1) then
          stop "NAMSURFACE: z0hav is not set"
        end if
      end if

    endif


    if(isurf == 1) then
      if(tsoilav(1) == -1 .or. tsoilav(2) == -1 .or. tsoilav(3) == -1 .or. tsoilav(4) == -1) then
        stop "NAMSURFACE: tsoil is not set"
      end if
      if(tsoildeepav == -1) then
        stop "NAMSURFACE: tsoildeep is not set"
      end if
      if(phiwav(1) == -1 .or. phiwav(2) == -1 .or. phiwav(3) == -1 .or. phiwav(4) == -1) then
        stop "NAMSURFACE: phiw is not set"
      end if
      if(rootfav(1) == -1 .or. rootfav(2) == -1 .or. rootfav(3) == -1 .or. rootfav(4) == -1) then
        stop "NAMSURFACE: rootf is not set"
      end if
      if(Cskinav == -1) then
        stop "NAMSURFACE: Cskinav is not set"
      end if
      if(lambdaskinav == -1) then
        stop "NAMSURFACE: lambdaskinav is not set"
      end if
      if(albedoav == -1) then
        stop "NAMSURFACE: albedoav is not set"
      end if
      if(Qnetav == -1) then
        stop "NAMSURFACE: Qnetav is not set"
      end if
      if(cvegav == -1) then
        stop "NAMSURFACE: cvegav is not set"
      end if
      if(rsminav == -1) then
        stop "NAMSURFACE: rsminav is not set"
      end if
      if(rssoilminav == -1) then
        print *,"WARNING: RSSOILMINAV is undefined... RSMINAV will be used as a proxy"
        rssoilminav = rsminav
      end if
      if(LAIav == -1) then
        stop "NAMSURFACE: LAIav is not set"
      end if
      if(gDav == -1) then
        stop "NAMSURFACE: gDav is not set"
      end if
      if(Wlav == -1) then
        stop "NAMSURFACE: Wlav is not set"
      end if
    end if

    if (isurf==1) then
      call initlsm
    end if

    allocate(rs(i2,j2))
    allocate(ra(i2,j2))
    if(isurf <= 2) then

      ! CvH set initial values for rs and ra to be able to compute qskin
      ra = 50.
      if(isurf == 1) then
        rs = 100.
      else
        rs = rsisurf2
      end if
    end if

    allocate(albedo(i2,j2))
    allocate(z0m(i2,j2))
    allocate(z0h(i2,j2))
    allocate(obl(i2,j2))
    allocate(tskin(i2,j2))
    allocate(qskin(i2,j2))
    allocate(Cm(i2,j2))
    allocate(Cs(i2,j2))

    if(rad_shortw .and. albedoav == -1) then
      stop "NAMSURFACE: albedoav is not set"
    end if
    if(iradiation == 1) then
      if(albedoav == -1) then
        stop "NAMSURFACE: albedoav is not set"
      end if
      allocate(swdavn(i2,j2,nradtime))
      allocate(swuavn(i2,j2,nradtime))
      allocate(lwdavn(i2,j2,nradtime))
      allocate(lwuavn(i2,j2,nradtime))
      swdavn =  0.
      swuavn =  0.
      lwdavn =  0.
      lwuavn =  0.
    end if

    albedo     = albedoav
    if(lhetero) then
      do j=1,j2
        tempy=patchynr(j)
        do i=1,i2
          tempx=patchxnr(i)
          z0m(i,j)   = z0mav_patch(tempx,tempy)
          z0h(i,j)   = z0hav_patch(tempx,tempy)
          if (.not. loldtable) then
            albedo(i,j) = albedo_patch(tempx,tempy)
            if(isurf .ne. 1) then
              rs(i,j)     = rsisurf2_patch(tempx,tempy)
            endif
          endif
        enddo
      enddo
    else
      z0m        = z0mav
      z0h        = z0hav
    endif

    ! 3. Initialize surface layer
    allocate(ustar   (i2,j2))
    ustar_3D(1:i2, 1:j2, 1:1) => ustar

    allocate(dudz    (i2,j2))
    allocate(dvdz    (i2,j2))

    allocate(thlflux (i2,j2))
    allocate(qtflux  (i2,j2))
    allocate(dqtdz   (i2,j2))
    allocate(dthldz  (i2,j2))
    allocate(svflux  (i2,j2,nsv))
    allocate(svs(nsv))

    allocate(horv(2:i1,2:j1))

    if (lrsAgs) then
      allocate(AnField   (2:i1,2:j1))
      allocate(gcco2Field(2:i1,2:j1))
      allocate(RespField (2:i1,2:j1))
      allocate(wco2Field (2:i1,2:j1))
      allocate(rsco2Field(2:i1,2:j1))
      allocate(fstrField (2:i1,2:j1))
      allocate(gc_old    (2:i1,2:j1))
      allocate(ci_old    (2:i1,2:j1))
      allocate(tauField  (2:i1,2:j1))
      allocate(ciField   (2:i1,2:j1))
      allocate(PARField  (2:i1,2:j1))
      if (lsplitleaf) then
        allocate(PARdirField   (2:i1,2:j1))
        allocate(PARdifField   (2:i1,2:j1))
      endif
    endif

    if(ltskininp) then ! Use tskin.inp.iexpnr.nc to define surface skin temperature
      !--- open tskin.inp.xxx.nc ---
      STATUS = NF90_OPEN('tskin.inp.'//cexpnr//'.nc', nf90_nowrite, NCID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      !--- get time dimensions
      STATUS = NF90_INQ_DIMID(NCID, "time", timeID)
      if (STATUS /= nf90_noerr) call handle_err(status)
      STATUS = nf90_INQUIRE_DIMENSION(NCID, timeID, len=nttskin, name=RecordDimName)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      !--- read time
      allocate(ttskin(nttskin))
      STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, ttskin, start=(/1/), count=(/nttskin/) )
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      !--- read tskin input
      allocate(tskininp(imax,jmax,nttskin))
      STATUS = NF90_INQ_VARID(NCID,'tskin', VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, tskininp, start=(/myidx*imax+1,myidy*jmax+1,1/), &
        & count=(/imax,jmax,nttskin/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_CLOSE(NCID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    endif

    dqtdz = 0 ! need to initialize, otherwise undefined in the first call to thermodynamics, before call surface (cold start)
    ustar = 0 ! need to initialize, otherwise undefined values in the corners in the first exchange

    !$acc enter data copyin(z0m, z0h, obl, tskin, qskin, Cm, Cs, &
    !$acc&                  ustar, dudz, dvdz, thlflux, qtflux, &
    !$acc&                  dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf)

    call timer_toc('modsurface/initsurface')
  end subroutine initsurface

!> Calculates the interaction with the soil, the surface temperature and humidity, and finally the surface fluxes.
  subroutine surface
    use moduser,    only : surf_user
    implicit none

    call timer_tic('modsurface/surface', 0)
    select case (isurf)
      case (1) ! Interactive land surface model
        call calc_mean_wind
        call getobl
        call calc_drag_coefficients
        call calc_aerodynamic_resistance
        call do_lsm
        call calc_friction_velocity
        call calc_surface_flux
        call calc_surface_gradients
      case (2) ! Forced surface temperature, fluxes calculated
        call calc_mean_wind
        call getobl
        call calc_drag_coefficients
        call calc_aerodynamic_resistance
        call presc_skin_temperature
        call qtsurf
        call calc_friction_velocity
        call calc_surface_flux
        call calc_surface_gradients
      case (3) ! Forced fluxes, surface temperature calculated
        call calc_mean_wind
        call calc_drag_coefficients
        call getobl
        call presc_friction_velocity
        call presc_surface_flux
        call calc_surface_gradients
        call calc_surface_scalars
      case (4) ! Forced moisture and heat flux, u_star and surface temperature calculated
        call calc_mean_wind
        call calc_drag_coefficients
        call getobl
        call calc_friction_velocity
        call presc_surface_flux
        call calc_surface_gradients
        call calc_surface_scalars
      case (10) ! User defined surface scheme
        call surf_user
      case (11) ! New LSM, handled by modlsm
        return
      case default
        stop "Invalid option selected for isurf"
    end select

    call timer_toc('modsurface/surface')
  end subroutine surface

  !> Calculates the drag coefficients \f$C_m\f$ and \f$C_s\f$
  subroutine calc_drag_coefficients
    use modglobal, only: i1, j1, fkar, zf
    implicit none

    integer :: i, j
    real :: logz

    !$acc parallel loop collapse(2) default(present)
    do j = 2, j1
      do i = 2, i1
        logz = log(zf(1) / z0m(i,j))
        Cm(i,j) = fkar**2 / (logz - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j)))**2
        Cs(i,j) = fkar**2 / (logz - psim(zf(1) / obl(i,j)) + psim(z0m(i,j) / obl(i,j))) / &
                  (log(zf(1) / z0h(i,j)) - psih(zf(1) / obl(i,j)) + psih(z0h(i,j) / obl(i,j)))
      end do
    end do

  end subroutine calc_drag_coefficients

  !> Calculates the aerodynamic resistance \f$r_a\f$
  subroutine calc_aerodynamic_resistance
    use modglobal, only: i1, j1
    implicit none

    integer :: i, j
    integer :: patchx, patchy

    if (lmostlocal) then
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          ra(i,j) = 1. / (Cs(i,j) * horv(i,j))
        end do
      end do
    else if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j)
          ra(i,j) = 1. / (Cs(i,j) * horvpatch(patchx, patchy))
        end do
      end do
    else
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          ra(i,j) = 1. / (Cs(i,j) * horvav)
        end do
      end do
    end if

  end subroutine calc_aerodynamic_resistance

  !> Prescribes the skin temperature
  subroutine presc_skin_temperature
    use modglobal, only: i1, j1, imax, xsize, dx, pi, rd, cp, pref0
    use modmpi, only: myidx
    implicit none

    integer :: i, j
    real :: x, dthls

    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          tskin(i,j) = thls_patch(patchxnr(i), patchynr(j))
        end do
      end do
    else
      dthls = dsst / (ps/pref0)**(rd/cp) ! convert dsst (temperature) to a thls amplitude
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
         do i = 2, i1
          x = (myidx * imax + i - 2 + 0.5) * dx
          tskin(i,j) = thls - dthls/2 * cos(2*pi*x / xsize) ! sst + optional cos-shaped perturbation (RCEMIP II)
        end do
      end do
    end if

  end subroutine presc_skin_temperature

  !> Calculates the surface temperature and specific humidity
  subroutine calc_surface_scalars
    use modglobal, only: i1, j1, rv, rd, ijtot
    use modfields, only: thl0, qt0
    use modmpi, only: D_MPI_ALLREDUCE, mpi_sum, comm3d, mpierr
    implicit none

    integer :: i, j
    integer :: patchx, patchy
    real :: thlsl, qtsl
    real :: lthls_patch(xpatches, ypatches)
    real :: lqts_patch(xpatches, ypatches)
    integer :: Npatch(xpatches, ypatches), SNpatch(xpatches, ypatches)

    ! TODO: check if splitting these loops speeds things up on the GPU (async)
    !$acc parallel loop collapse(2) default(present) 
    do j = 2, j1
      do i = 2, i1
        tskin(i,j) = min(max(thlflux(i,j) / (Cs(i,j) * horv(i,j)), -10.), 10.) + thl0(i,j,1)
        qskin(i,j) = min(max(qtflux(i,j) / (Cs(i,j) * horv(i,j)), -5.e-2), 5.e-2) + qt0(i,j,1)
      end do
    end do

    thlsl = 0.0
    qtsl = 0.0
    !$acc parallel loop collapse(2) default(present) reduction(+: thlsl, qtsl)
    do j = 2, j1
      do i = 2, i1
        thlsl = thlsl + tskin(i,j)
        qtsl = qtsl + qskin(i,j)
      end do
    end do

    call D_MPI_ALLREDUCE(thlsl, thls, 1, MPI_SUM, comm3d, mpierr)
    call D_MPI_ALLREDUCE(qtsl, qts, 1, MPI_SUM, comm3d, mpierr)

    thls = thls / ijtot
    qts = qts / ijtot
    thvs = thls * (1. + (rv/rd - 1.) * qts)

    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j)
          lthls_patch(patchx, patchy) = lthls_patch(patchx, patchy) + tskin(i,j)
          lqts_patch(patchx, patchy) = lqts_patch(patchx, patchy) + qskin(i,j)
          Npatch(patchx, patchy) = Npatch(patchx, patchy) + 1
        end do
      end do

      call D_MPI_ALLREDUCE(lthls_patch(1:xpatches, 1:ypatches), thls_patch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(lqts_patch(1:xpatches, 1:ypatches), qts_patch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(Npatch(1:xpatches, 1:ypatches), SNpatch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)
      thls_patch = thls_patch / SNpatch
      qts_patch = qts_patch / SNpatch
      thvs_patch = thls_patch * (1. + (rv/rd - 1.) * qts_patch)
    end if

  end subroutine calc_surface_scalars

  !> Calculates the maginitude of the wind vector at the first level
  subroutine calc_mean_wind
    use modglobal, only: i1, j1, cu, cv
    use modfields, only: u0, v0, u0av, v0av
    use modmpi, only: D_MPI_ALLREDUCE, mpi_sum, comm3d, mpierr
    implicit none

    integer :: i, j, patchx, patchy
    real :: upcu, vpcv
    real :: upatch(xpatches, ypatches), vpatch(xpatches, ypatches)
    real :: Supatch(xpatches, ypatches), Svpatch(xpatches, ypatches)
    integer :: Npatch(xpatches, ypatches), SNpatch(xpatches, ypatches)

    !$acc parallel loop collapse(2) default(present) private(upcu, vpcv)
    do j = 2, j1
      do i = 2, i1
        upcu = 0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
        vpcv = 0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
        horv(i,j) = sqrt(upcu**2. + vpcv**2.)
        horv(i,j) = max(horv(i,j), 0.1)
      end do
    end do

    if (lhetero) then
      upatch = 0
      vpatch = 0
      Npatch = 0

      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j)
          upatch(patchx,patchy) = upatch(patchx,patchy) + 0.5 * (u0(i,j,1) + u0(i+1,j,1))
          vpatch(patchx,patchy) = vpatch(patchx,patchy) + 0.5 * (v0(i,j,1) + v0(i,j+1,1))
          Npatch(patchx,patchy) = Npatch(patchx,patchy) + 1
        enddo
      enddo

      call D_MPI_ALLREDUCE(upatch(1:xpatches, 1:ypatches), Supatch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(vpatch(1:xpatches, 1:ypatches), Svpatch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(Npatch(1:xpatches, 1:ypatches), SNpatch(1:xpatches, 1:ypatches), &
                           xpatches*ypatches, MPI_SUM, comm3d, mpierr)

      horvpatch = sqrt(((Supatch/SNpatch) + cu) **2. + ((Svpatch/SNpatch) + cv) ** 2.)
      horvpatch = max(horvpatch, 0.1)
    else
      !$acc update self(u0av(1), v0av(1))
      horvav = sqrt(u0av(1)**2. + v0av(1)**2.)
      horvav = max(horvav, 0.1)
    end if

    !$acc wait

  end subroutine calc_mean_wind

  !> Calculates the friction velocity \f$u_*\f$
  subroutine calc_friction_velocity
    use modglobal, only: i1 ,j1, lopenbc, lboundary, lperiodic
    use modmpi, only: excjs, openboundary_excjs
    implicit none

    integer :: i, j

    if (lmostlocal) then
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          ustar(i,j) = sqrt(Cm(i,j)) * horv(i,j)
          ustar(i,j) = max(ustar(i,j), 1.e-2)
        end do
      end do
    else if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          ustar(i,j) = sqrt(Cm(i,j)) * horvpatch(patchxnr(i), patchynr(j))
          ustar(i,j) = max(ustar(i,j), 1.e-2)
        end do
      end do
    else
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          ustar(i,j) = sqrt(Cm(i,j)) * horvav
          ustar(i,j) = max(ustar(i,j), 1.e-2)
        end do
      end do
    end if
    
    if ( lopenbc ) then
      call openboundary_excjs(ustar_3D, 2,i1,2,j1,1,1,1,1, &
                             (.not.lboundary(1:4)).or.lperiodic(1:4))
    else
       call excjs(ustar_3D,2,i1,2,j1,1,1,1,1)
    endif
  end subroutine calc_friction_velocity

  !> Prescribes the friction velocity \f$u_*\f$
  subroutine presc_friction_velocity
    use modglobal, only: i1, j1, lopenbc, lboundary, lperiodic
    use modmpi, only: excjs, openboundary_excjs
    implicit none

    integer :: i, j

    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          ustar(i,j) = ustin_patch(patchxnr(i), patchynr(j))
          ustar(i,j) = max(ustar(i,j), 1.e-2)
        end do
      end do
    else
      !$acc parallel loop collapse(2) default(present) 
      do j = 2, j1
        do i = 2, i1
          ustar(i,j) = ustin
          ustar(i,j) = max(ustar(i,j), 1.e-2)
        end do
      end do
    end if
    
    if ( lopenbc ) then
      call openboundary_excjs(ustar_3D, 2,i1,2,j1,1,1,1,1, &
                             (.not.lboundary(1:4)).or.lperiodic(1:4))
    else
       call excjs(ustar_3D,2,i1,2,j1,1,1,1,1)
    endif
  end subroutine presc_friction_velocity

  !> Calculates the surfaces fluxes using the scalar values at the surface and
  !> the calculated aerodynamic resistance
  subroutine calc_surface_flux
    use modglobal, only: i1, j1, nsv, ijtot
    use modfields, only: thl0, qt0
    use modmpi, only: D_MPI_ALLREDUCE, mpi_sum, comm3d, mpierr
    implicit none

    integer :: i, j, n
    real :: ust, ustl, wtsurfl, wqsurfl

    ! Uniform sensible and latent heat flux
    if (lsmoothflux) then
      ustl = 0.0
      wtsurfl = 0.0
      wqsurfl = 0.0
      !$acc parallel loop collapse(2) default(present) reduction(+: ustl, wtsurfl, wqsurfl)
      do j = 2, j1
        do i = 2, i1
          ustl = ustl + ustar(i,j)
          wtsurfl = wtsurfl + thlflux(i,j)
          wqsurfl = wqsurfl + qtflux(i,j)
        end do
      end do

      call D_MPI_ALLREDUCE(ustl, ust, 1, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(wtsurfl, wtsurf, 1, MPI_SUM, comm3d, mpierr)
      call D_MPI_ALLREDUCE(wqsurfl, wqsurf, 1, MPI_SUM, comm3d, mpierr)

      ust = ust / ijtot
      wtsurf = wtsurf / ijtot
      wqsurf = wqsurf / ijtot

      call presc_surface_flux
    else
      !$acc parallel loop collapse(2) default(present)
      do j = 2, j1
        do i = 2, i1
          thlflux(i,j) = - (thl0(i,j,1) - tskin(i,j)) / ra(i,j)
          qtflux(i,j) = - (qt0(i,j,1) - qskin(i,j)) / ra(i,j)
        end do
      end do

      ! Passive scalars
      if (nsv > 0) then
        if (lhetero) then
          do n = 1, nsv
            do j = 2, j1
              do i = 2, i1
                svflux(i,j,n) = wsv_patch(n, patchxnr(i), patchynr(j))
              end do
            end do
          end do
        else
          !$acc parallel loop collapse(3) default(present)
          do n = 1, nsv
            do j = 2, j1
              do i = 2, i1
                svflux(i,j,n) = wsvsurf(n)
              end do
            end do
          end do
        end if
      end if

      if (lCO2Ags) then
        !$acc parallel loop collapse(3) default(present)
        do n = 1, indCO2
          do j = 2, j1
            do i = 2, i1
              svflux(i,j,n) = CO2flux(i,j)
            end do
          end do
        end do
      end if
    end if

  end subroutine calc_surface_flux

  !> Prescribes the surface fluxes
  subroutine presc_surface_flux
    use modglobal, only: i1, j1, nsv
    implicit none

    integer :: i, j, n

    if (lhetero) then
      do j = 2, j1
        do i = 2, i1
          thlflux(i,j) = wt_patch(patchxnr(i), patchynr(j))
          qtflux(i,j) = wq_patch(patchxnr(i), patchynr(j))
        end do
      end do

      if (nsv > 0) then
        do n = 1, nsv
          do j = 2, j1
            do i = 2, i1
              svflux(i,j,n) = wsv_patch(n, patchxnr(i), patchynr(j))
            end do
          end do
        end do
      end if
    else
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 2, j1
        do i = 2, i1
          thlflux(i,j) = wtsurf
          qtflux(i,j) = wqsurf
        end do
      end do

      if (nsv > 0) then
        !$acc parallel loop collapse(3) default(present) async(2)
        do n = 1, nsv
          do j = 2, j1
            do i = 2, i1
              svflux(i,j,n) = wsvsurf(n)
            end do
          end do
        end do
      end if
      !$acc wait(1,2)
    end if

  end subroutine presc_surface_flux

  !> Calculates the surface gradients
  subroutine calc_surface_gradients
    use modglobal, only: i1, j1, cu, cv, zf, fkar
    use modfields, only: u0, v0
    implicit none

    integer :: i, j
    real(field_r) :: scaling
    real(field_r) :: upcu, vpcv
    real(field_r) :: phimzf, phihzf

    scaling = 1.0 / (fkar * zf(1))

    !$acc parallel loop collapse(2) default(present)
    do j = 2, j1
      do i = 2, i1
        phimzf = phim(zf(1) / obl(i,j))
        phihzf = phih(zf(1) / obl(i,j))

        ! Momentum fluxes
        upcu = 0.5_field_r * (u0(i,j,1) + u0(i+1,j,1)) + cu
        vpcv = 0.5_field_r * (v0(i,j,1) + v0(i,j+1,1)) + cv
        dudz(i,j) = ustar(i,j) * phimzf * scaling * (upcu / horv(i,j))
        dvdz(i,j) = ustar(i,j) * phimzf * scaling * (vpcv / horv(i,j))

        ! Scalar fluxes
        dthldz(i,j) = - thlflux(i,j) / ustar(i,j) * phihzf * scaling
        dqtdz(i,j) = - qtflux(i,j) / ustar(i,j) * phihzf * scaling
      end do
    end do
  end subroutine calc_surface_gradients

!> Calculate the surface humidity assuming saturation.
  subroutine qtsurf
    use modglobal,   only : tmelt,bt,at,rd,rv,cp,es0,pref0,ijtot,i1,j1
    use modfields,   only : qt0
    !use modsurfdata, only : rs, ra
    use modmpi,      only : mpierr,comm3d,mpi_sum, D_MPI_ALLREDUCE

    implicit none
    real       :: exner, tsurf, qsatsurf, surfwet, es, qtsl
    integer    :: i,j, patchx, patchy
    integer    :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)
    real       :: lqts_patch(xpatches,ypatches)

    patchx = 0
    patchy = 0

    if(isurf <= 2) then
      qtsl = 0.
      !$acc parallel loop collapse(2) default(present) reduction(+: qtsl)
      do j = 2, j1
        do i = 2, i1
          exner      = (ps / pref0)**(rd/cp)
          tsurf      = tskin(i,j) * exner
          es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
          qsatsurf   = rd / rv * es / ps
          surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
          qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * qt0(i,j,1)
          qtsl       = qtsl + qskin(i,j)
        end do
      end do

      call D_MPI_ALLREDUCE(qtsl, qts, 1,  MPI_SUM, comm3d,mpierr)
      qts  = qts / ijtot
      thvs = thls * (1. + (rv/rd - 1.) * qts)

      if (lhetero) then
        lqts_patch = 0.
        Npatch     = 0
        do j = 2, j1
          do i = 2, i1
            patchx     = patchxnr(i)
            patchy     = patchynr(j)
            exner      = (ps_patch(patchx,patchy) / pref0)**(rd/cp)
            tsurf      = tskin(i,j) * exner
            es         = es0 * exp(at*(tsurf-tmelt) / (tsurf-bt))
            qsatsurf   = rd / rv * es / ps_patch(patchx,patchy)
            surfwet    = ra(i,j) / (ra(i,j) + rs(i,j))
            qskin(i,j) = surfwet * qsatsurf + (1. - surfwet) * qt0(i,j,1)

            lqts_patch(patchx,patchy) = lqts_patch(patchx,patchy) + qskin(i,j)
            Npatch(patchx,patchy)     = Npatch(patchx,patchy)     + 1
          enddo
        enddo
        call D_MPI_ALLREDUCE(lqts_patch(1:xpatches,1:ypatches), qts_patch(1:xpatches,1:ypatches),&
        xpatches*ypatches,MPI_SUM, comm3d,mpierr)
        call D_MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)    , SNpatch(1:xpatches,1:ypatches)  ,&
        xpatches*ypatches,MPI_SUM, comm3d,mpierr)
        qts_patch = qts_patch / SNpatch
        thvs_patch = thls_patch * (1. + (rv/rd - 1.) * qts_patch)
      endif
    end if

    return

  end subroutine qtsurf

!> Calculates the Obukhov length iteratively.
  subroutine getobl
    use modglobal, only : zf, rv, rd, grav, i1, j1, i2, j2, cu, cv
    use modfields, only : thl0av, qt0av, u0, v0, thl0, qt0, u0av, v0av
    use modmpi,    only : mpierr,comm3d,mpi_sum,D_MPI_ALLREDUCE, &
                          D_MPI_BCAST
    implicit none

    integer             :: i,j,iter,patchx,patchy
    real                :: thv, thvsl, horv2, oblavl, thvpatch(xpatches,ypatches), horvpatch(xpatches,ypatches)
    real                :: L, Lend, Lstart, Lold
    real                :: Rib, fx, fxdif
    real                :: upcu, vpcv
    real                :: upatch(xpatches,ypatches), vpatch(xpatches,ypatches)
    real                :: Supatch(xpatches,ypatches), Svpatch(xpatches,ypatches)
    integer             :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)
    real                :: lthlpatch(xpatches,ypatches), thlpatch(xpatches,ypatches),&
                           lqpatch(xpatches,ypatches), qpatch(xpatches,ypatches)
    real                :: loblpatch(xpatches,ypatches)


    if (lneutral) then
      !$acc parallel loop collapse(2) default(present)
      do i=1,i2
        do j=1,j2
          obl(i,j) = -1.e10
        end do
      end do
      oblav = -1.e10
      return
    endif

    if(lmostlocal) then

      oblavl = 0.

      do i=2,i1
        do j=2,j1
          thv     =   thl0(i,j,1)  * (1. + (rv/rd - 1.) * qt0(i,j,1))
          thvsl   =   tskin(i,j)   * (1. + (rv/rd - 1.) * qskin(i,j))
          upcu    =   0.5 * (u0(i,j,1) + u0(i+1,j,1)) + cu
          vpcv    =   0.5 * (v0(i,j,1) + v0(i,j+1,1)) + cv
          horv2   =   upcu ** 2. + vpcv ** 2.
          horv2   =   max(horv2, min_horv**2)

          if(lhetero) then
            patchx = patchxnr(i)
            patchy = patchynr(j)
            Rib    = grav / thvs_patch(patchx,patchy) * zf(1) * (thv - thvsl) / horv2
          else
            Rib    = grav / thvs * zf(1) * (thv - thvsl) / horv2
          endif

          if (Rib == 0) then
             ! Rib can be 0 if there is no surface flux
             ! L is capped at 1e6 below, so use the same cap here
             L = 1e6
             write(*,*) 'Obukhov length: Rib = 0 -> setting L=1e6'
          else
             iter = 0
             L = obl(i,j)

             if(Rib * L < 0. .or. abs(L) == 1e5) then
                if(Rib > 0) L = 0.01
                if(Rib < 0) L = -0.01
             end if

             do while (.true.)
                iter    = iter + 1
                Lold    = L
                fx      = Rib - zf(1) / L * (log(zf(1) / z0h(i,j)) - psih(zf(1) / L) + psih(z0h(i,j) / L)) /&
                     (log(zf(1) / z0m(i,j)) - psim(zf(1) / L) + psim(z0m(i,j) / L)) ** 2.
                Lstart  = L - 0.001*L
                Lend    = L + 0.001*L
                fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lstart) + psih(z0h(i,j) / Lstart)) /&
                     (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lstart) + psim(z0m(i,j) / Lstart)) ** 2.) - (-zf(1) / Lend * &
                     (log(zf(1) / z0h(i,j)) - psih(zf(1) / Lend) + psih(z0h(i,j) / Lend)) / (log(zf(1) / z0m(i,j)) - psim(zf(1) / Lend)&
                     + psim(z0m(i,j) / Lend)) ** 2.) ) / (Lstart - Lend)
                L = L - fx / fxdif
                if(Rib * L < 0. .or. abs(L) == 1e5) then
                   if(Rib > 0) L = 0.01
                   if(Rib < 0) L = -0.01
                end if
                if(abs((L - Lold)/L) < 1e-4) exit
                if(iter > 1000) stop 'Obukhov length calculation does not converge!'
             end do

             if (abs(L)>1e6) L = sign(1.0e6,L)
          end if
          obl(i,j) = L

        end do
      end do
    end if

    if(lhetero) then
      upatch    = 0
      vpatch    = 0
      Npatch    = 0
      lthlpatch = 0
      lqpatch   = 0
      loblpatch = 0

      do j = 2, j1
        do i = 2, i1
          patchx = patchxnr(i)
          patchy = patchynr(j)
          upatch(patchx,patchy)    = upatch(patchx,patchy)    + 0.5 * (u0(i,j,1) + u0(i+1,j,1))
          vpatch(patchx,patchy)    = vpatch(patchx,patchy)    + 0.5 * (v0(i,j,1) + v0(i,j+1,1))
          Npatch(patchx,patchy)    = Npatch(patchx,patchy)    + 1
          lthlpatch(patchx,patchy) = lthlpatch(patchx,patchy) + thl0(i,j,1)
          lqpatch(patchx,patchy)   = lqpatch(patchx,patchy)   + qt0(i,j,1)
          loblpatch(patchx,patchy) = loblpatch(patchx,patchy) + obl(i,j)
        enddo
      enddo

      call D_MPI_ALLREDUCE(upatch(1:xpatches,1:ypatches)   ,Supatch(1:xpatches,1:ypatches) ,xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(vpatch(1:xpatches,1:ypatches)   ,Svpatch(1:xpatches,1:ypatches) ,xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)   ,SNpatch(1:xpatches,1:ypatches) ,xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(lthlpatch(1:xpatches,1:ypatches),thlpatch(1:xpatches,1:ypatches),xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(lqpatch(1:xpatches,1:ypatches)  ,qpatch(1:xpatches,1:ypatches)  ,xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(loblpatch(1:xpatches,1:ypatches),oblpatch(1:xpatches,1:ypatches),xpatches*ypatches,&
      MPI_SUM, comm3d,mpierr)

      horvpatch = sqrt(((Supatch/SNpatch) + cu) **2. + ((Svpatch/SNpatch) + cv) ** 2.)
      horvpatch = max(horvpatch, min_horv)

      thlpatch  = thlpatch / SNpatch
      qpatch    = qpatch   / SNpatch
      oblpatch  = oblpatch / SNpatch

      thvpatch  = thlpatch * (1. + (rv/rd - 1.) * qpatch)

      do patchy = 1, ypatches
        do patchx = 1, xpatches
          Rib   = grav / thvs_patch(patchx,patchy) * zf(1) *&
          (thvpatch(patchx,patchy) - thvs_patch(patchx,patchy)) / (horvpatch(patchx,patchy) ** 2.)
          iter = 0
          L = oblpatch(patchx,patchy)

          if(Rib * L < 0. .or. abs(L) == 1e5) then
            if(Rib > 0) L = 0.01
            if(Rib < 0) L = -0.01
          end if

          do while (.true.)
            iter    = iter + 1
            Lold    = L
            fx      = Rib - zf(1) / L * (log(zf(1) / z0hav_patch(patchx,patchy)) - psih(zf(1) / L) +&
            psih(z0hav_patch(patchx,patchy) / L)) / (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / L) &
            + psim(z0mav_patch(patchx,patchy) / L)) ** 2.
            Lstart  = L - 0.001*L
            Lend    = L + 0.001*L
            fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0hav_patch(patchx,patchy)) - psih(zf(1) / Lstart) +&
            psih(z0hav_patch(patchx,patchy) / Lstart)) / (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / Lstart) + &
            psim(z0mav_patch(patchx,patchy) / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0hav_patch(patchx,patchy)) &
            - psih(zf(1) / Lend) + psih(z0hav_patch(patchx,patchy) / Lend)) / &
            (log(zf(1) / z0mav_patch(patchx,patchy)) - psim(zf(1) / Lend) + psim(z0mav_patch(patchx,patchy) / Lend)) ** 2.) )&
            / (Lstart - Lend)
            L       = L - fx / fxdif
            if(Rib * L < 0. .or. abs(L) == 1e5) then
              if(Rib > 0) L = 0.01
              if(Rib < 0) L = -0.01
            end if
            if(abs((L - Lold)/L) < 1e-4) exit
          end do

          if (abs(L)>1e6) L = sign(1.0e6,L)
          oblpatch(patchx,patchy) = L
        enddo
      enddo
      if(.not. lmostlocal) then
        do i=1,i2
          do j=1,j2
            obl(i,j) = oblpatch(patchxnr(i),patchynr(j))
          enddo
        enddo
      endif
    endif ! if lhetero

    !CvH also do a global evaluation if lmostlocal = .true. to get an appropriate local mean
    !$acc update self(thl0av(1), qt0av(1))
    thv    = thl0av(1) * (1. + (rv/rd - 1.) * qt0av(1))

    !$acc update self(u0av(1), v0av(1))
    horv2 = u0av(1)**2. + v0av(1)**2.
    horv2 = max(horv2, min_horv**2)

    Rib   = grav / thvs * zf(1) * (thv - thvs) / horv2
    if (Rib == 0) then
       ! Rib can be 0 if there is no surface flux
       ! L is capped at 1e6 below, so use the same cap here
       L = 1e6
       write(*,*) 'Obukhov length: Rib = 0 -> setting L=1e6 (2nd point)'
    else
       iter = 0
       L = oblav

       if(Rib * L < 0. .or. abs(L) == 1e5) then
          if(Rib > 0) L = 0.01
          if(Rib < 0) L = -0.01
       end if

       do while (.true.)
          iter    = iter + 1
          Lold    = L
          fx      = Rib - zf(1) / L * (log(zf(1) / z0hav) - psih(zf(1) / L) + psih(z0hav / L)) /&
               (log(zf(1) / z0mav) - psim(zf(1) / L) + psim(z0mav / L)) ** 2.
          Lstart  = L - 0.001*L
          Lend    = L + 0.001*L
          fxdif   = ( (- zf(1) / Lstart * (log(zf(1) / z0hav) - psih(zf(1) / Lstart) + psih(z0hav / Lstart)) /&
               (log(zf(1) / z0mav) - psim(zf(1) / Lstart) + psim(z0mav / Lstart)) ** 2.) - (-zf(1) / Lend * (log(zf(1) / z0hav) &
               - psih(zf(1) / Lend) + psih(z0hav / Lend)) / (log(zf(1) / z0mav) - psim(zf(1) / Lend) &
               + psim(z0mav / Lend)) ** 2.) ) / (Lstart - Lend)
          L       = L - fx / fxdif
          if(Rib * L < 0. .or. abs(L) == 1e5) then
             if(Rib > 0) L = 0.01
             if(Rib < 0) L = -0.01
          end if
          if(abs((L - Lold)/L) < 1e-4) exit
          if(iter > 1000) stop 'Obukhov length calculation does not converge!'
       end do

       if (abs(L)>1e6) L = sign(1.0e6,L)
       if(.not. lmostlocal) then
          if(.not. lhetero) then
            !$acc parallel loop collapse(2) default(present)
            do j = 1, j2
              do i = 1, i2
                obl(i,j) = L
              end do
            end do
          endif
       end if
    end if
    oblav = L

    return
  end subroutine getobl

  pure function psim(zeta)
    implicit none

    real             :: psim
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 16. * zeta) ** (0.25)
      psim  = 3.14159265 / 2. - 2. * atan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
      ! CvH use Wilson, 2001 rather than Businger-Dyer for correct free convection limit
      !x     = (1. + 3.6 * abs(zeta) ** (2./3.)) ** (-0.5)
      !psim = 3. * log( (1. + 1. / x) / 2.)
    else
      psim  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - zeta - (10./3.) / 0.35
    end if

    return
  end function psim

  pure function psih(zeta)

    implicit none

    real             :: psih
    real, intent(in) :: zeta
    real             :: x

    if(zeta <= 0) then
      x     = (1. - 16. * zeta) ** (0.25)
      psih  = 2. * log( (1. + x ** 2.) / 2. )
      ! CvH use Wilson, 2001
      !x     = (1. + 7.9 * abs(zeta) ** (2./3.)) ** (-0.5)
      !psih  = 3. * log( (1. + 1. / x) / 2.)
    else
      psih  = -2./3. * (zeta - 5./0.35)*exp(-0.35 * zeta) - (1. + (2./3.) * zeta) ** (1.5) - (10./3.) / 0.35 + 1.
    end if

    return
  end function psih

  ! stability function Phi for momentum.
  ! Many functional forms of Phi have been suggested, see e.g. Optis 2015
  ! Phi and Psi above are related by an integral and should in principle match,
  ! currently they do not.
  ! FJ 2018: For very stable situations, zeta > 1 add cap to phi - the linear expression is valid only for zeta < 1
  function phim(zeta)
    !$acc routine seq
    implicit none
    real             :: phim
    real, intent(in) :: zeta

    if (zeta < 0.) then ! unstable
       phim = (1.-16.*zeta)**(-0.25)
       !phimzf = (1. + 3.6 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
    elseif ( zeta < 1.) then  ! 0 < zeta < 1, stable
       phim = (1.+5.*zeta)
    else
       phim = 6 ! cap phi when z/L > 1
    endif

    return
  end function phim
  
  ! stability function Phi for heat.
  function phih(zeta)
    !$acc routine seq
    implicit none
    real             :: phih
    real, intent(in) :: zeta

    if (zeta < 0.) then ! unstable
       phih = (1.-16.*zeta)**(-0.50)
       !phihzf = (1. + 7.9 * (-zf(1)/obl(i,j))**(2./3.))**(-0.5)
    elseif ( zeta < 1.) then  ! 0 < zf(1) / obl < 1, stable
       phih = (1.+5.*zeta)
    else
     phih = 6  ! cap phi when z/L > 1
    endif

    return
  end function phih


  function E1(x)
  implicit none
    real             :: E1
    real, intent(in) :: x
    real             :: E1sum!, factorial
    integer          :: k

    E1sum = 0.0
    do k=1,99
      !E1sum = E1sum + (-1.0) ** (k + 0.0) * x ** (k + 0.0) / ( (k + 0.0) * factorial(k) )
       E1sum = E1sum + (-1.0 * x) ** k / ( k * factorial(k) )  ! FJ changed this for compilation with cray fortran

    end do
    E1 = -0.57721566490153286060 - log(x) - E1sum

    return
  end function E1

  function factorial(k)
  implicit none
    real                :: factorial
    integer, intent(in) :: k
    integer             :: n
    factorial = 1.0
    do n = 2, k
      factorial = factorial * n
    enddo
  end function factorial

  function patchxnr(xpos)
    use modmpi,     only : myidx
    use modglobal,  only : imax,itot
    implicit none
    integer, intent(in) :: xpos

    integer             :: patchxnr
    integer             :: positionx

    ! Converting the j position to the real j position by taking the processor number into account
    ! First grid point lies at j = 2. Make position = 0 for first grid point

    positionx = xpos + (myidx * imax) - 2

    ! Account for border grid points
    if (positionx .lt. 0)    positionx = positionx + itot
    if (positionx .ge. itot) positionx = positionx - itot

    ! Convert position to patch number
    patchxnr  = 1 + (positionx*xpatches)/itot

    return
  end function

  function patchynr(ypos)
    use modmpi,     only : myidy
    use modglobal,  only : jmax,jtot
    implicit none
    integer, intent(in) :: ypos
    integer             :: patchynr
    integer             :: positiony

    ! Converting the j position to the real j position by taking the processor number into account
    ! First grid point lies at j = 2. Make position = 0 for first gridpoint
    positiony = ypos + (myidy * jmax) - 2

    ! Account for border grid points
    if (positiony .lt. 0)    positiony = positiony + jtot
    if (positiony .ge. jtot) positiony = positiony - jtot

    ! Convert position to patch number
    patchynr  = 1 + (positiony*ypatches)/jtot

    return
  end function

  subroutine exitsurface
    implicit none

    !$acc exit data delete(z0m, z0h, obl, tskin, qskin, Cm, Cs, &
    !$acc&                 ustar, dudz, dvdz, thlflux, qtflux, &
    !$acc&                 dqtdz, dthldz, svflux, svs, horv, ra, rs, wsvsurf)

    return
  end subroutine exitsurface

  subroutine initlsm
    use modglobal, only : i2,j2
    integer :: k,j,i,tempy,tempx

    ! 1.1  -   Allocate arrays
    allocate(zsoil(ksoilmax))
    allocate(zsoilc(ksoilmax))
    allocate(dzsoil(ksoilmax))
    allocate(dzsoilh(ksoilmax))

    allocate(lambda(i2,j2,ksoilmax))
    allocate(lambdah(i2,j2,ksoilmax))
    allocate(lambdas(i2,j2,ksoilmax))
    allocate(lambdash(i2,j2,ksoilmax))
    allocate(gammas(i2,j2,ksoilmax))
    allocate(gammash(i2,j2,ksoilmax))
    allocate(Dh(i2,j2,ksoilmax))
    allocate(phiw(i2,j2,ksoilmax))
    allocate(phiwm(i2,j2,ksoilmax))
    allocate(phifrac(i2,j2,ksoilmax))
    allocate(pCs(i2,j2,ksoilmax))
    allocate(rootf(i2,j2,ksoilmax))
    allocate(tsoil(i2,j2,ksoilmax))
    allocate(tsoilm(i2,j2,ksoilmax))
    allocate(tsoildeep(i2,j2))
    allocate(phitot(i2,j2))

    ! 1.2   -  Initialize arrays
    ! First test, pick ECMWF config
    dzsoil(1) = 0.07
    dzsoil(2) = 0.21
    dzsoil(3) = 0.72
    dzsoil(4) = 1.89

    !! 1.3   -  Calculate vertical layer properties
    zsoil(1)  = dzsoil(1)
    do k = 2, ksoilmax
      zsoil(k) = zsoil(k-1) + dzsoil(k)
    end do
    zsoilc = -(zsoil-0.5*dzsoil)
    do k = 1, ksoilmax-1
      dzsoilh(k) = 0.5 * (dzsoil(k+1) + dzsoil(k))
    end do
    dzsoilh(ksoilmax) = 0.5 * dzsoil(ksoilmax)

    phitot = 0.0
    ! Calculate conductivity saturated soil
    lambdasat = lambdasm ** (1. - phi) * lambdaw ** (phi)

    ! 2.1  -   Allocate arrays
    allocate(Qnet(i2,j2))
    allocate(LE(i2,j2))
    allocate(H(i2,j2))
    allocate(G0(i2,j2))
    allocate(CO2flux(i2,j2))

    allocate(rsveg(i2,j2))
    allocate(rsmin(i2,j2))
    allocate(rssoil(i2,j2))
    allocate(rssoilmin(i2,j2))
    allocate(cveg(i2,j2))
    allocate(cliq(i2,j2))
    allocate(tendskin(i2,j2))
    allocate(tskinm(i2,j2))
    allocate(Cskin(i2,j2))
    allocate(lambdaskin(i2,j2))
    allocate(LAI(i2,j2))
    allocate(gD(i2,j2))
    allocate(Wl(i2,j2))
    allocate(Wlm(i2,j2))

    if(lhetero) then
      do j=1,j2
        tempy=patchynr(j)
        do i=1,i2
          tempx=patchxnr(i)

          phiw( i,j,:)   = phiw_patch( :,  tempx,tempy)
          rootf(i,j,:)   = rootf_patch(:,  tempx,tempy)
          tsoil(i,j,:)   = tsoil_patch(:,  tempx,tempy)
          tsoildeep(i,j) = tsoildeep_patch(tempx,tempy)

          Qnet       (i,j) = Qnet_patch      (tempx,tempy)
          Cskin      (i,j) = Cskin_patch     (tempx,tempy)
          lambdaskin (i,j) = lambdaskin_patch(tempx,tempy)
          rsmin      (i,j) = rsmin_patch     (tempx,tempy)
          rssoilmin  (i,j) = rsmin_patch     (tempx,tempy)
          LAI        (i,j) = LAI_patch       (tempx,tempy)
          gD         (i,j) = gD_patch        (tempx,tempy)
          cveg       (i,j) = cveg_patch      (tempx,tempy)
          Wl         (i,j) = Wl_patch        (tempx,tempy)
        enddo
      enddo

      do k = 1, ksoilmax
        phitot(:,:) = phitot(:,:) + rootf(:,:,k)*phiw(:,:,k)
      end do

      do k = 1, ksoilmax
        phifrac(:,:,k) = rootf(:,:,k)*phiw(:,:,k) / phitot(:,:)
      end do

    else
      ! 1.4   -   Set evaporation related properties
      ! Set water content of soil - constant in this scheme
      phiw(:,:,1) = phiwav(1)
      phiw(:,:,2) = phiwav(2)
      phiw(:,:,3) = phiwav(3)
      phiw(:,:,4) = phiwav(4)

      do k = 1, ksoilmax
        phitot(:,:) = phitot(:,:) + rootfav(k)*phiw(:,:,k)
      end do

      do k = 1, ksoilmax
        phifrac(:,:,k) = rootfav(k)*phiw(:,:,k) / phitot(:,:)
      end do

      ! Set root fraction per layer for short grass
      rootf(:,:,1) = rootfav(1)
      rootf(:,:,2) = rootfav(2)
      rootf(:,:,3) = rootfav(3)
      rootf(:,:,4) = rootfav(4)

      tsoil(:,:,1)   = tsoilav(1)
      tsoil(:,:,2)   = tsoilav(2)
      tsoil(:,:,3)   = tsoilav(3)
      tsoil(:,:,4)   = tsoilav(4)
      tsoildeep(:,:) = tsoildeepav

      ! 2    -   Initialize land surface
      Qnet       = Qnetav
      Cskin      = Cskinav
      lambdaskin = lambdaskinav
      rsmin      = rsminav
      rssoilmin  = rssoilminav
      LAI        = LAIav
      gD         = gDav
      cveg       = cvegav
      Wl         = Wlav
    endif
    cliq       = 0.
  end subroutine initlsm


!> Calculates surface resistance, temperature and moisture using the Land Surface Model
  subroutine do_lsm

    use modglobal, only : pref0,boltz,cp,rd,rhow,rlv,i1,j1,rdt,ijtot,rk3step,nsv,xtime,rtimee,xday,xlat,xlon
    use modfields, only : ql0,qt0,thl0,rhof,presf,svm
    use modraddata,only : iradiation,useMcICA,swd,swu,lwd,lwu,irad_par,swdir,swdif,zenith
    use modmpi, only :comm3d,mpi_sum,mpierr,myid, D_MPI_ALLREDUCE
    use modmicrodata, only : imicro,imicro_bulk

    real     :: f1, f2, f3, f4 ! Correction functions for Jarvis-Stewart
    integer  :: i, j, k, itg
    integer  :: patchx, patchy
    real     :: rk3coef,thlsl

    real     :: swdav, swuav, lwdav, lwuav
    real     :: exner, exnera, tsurfm, Tatm, e,esat, qsat, desatdT, dqsatdT, Acoef, Bcoef
    real     :: fH, fLE, fLEveg, fLEsoil, fLEliq, LEveg, LEsoil, LEliq
    real     :: Wlmx

    real     :: CO2ags, CO2comp, gm, fmin0, fmin, esatsurf, Ds, D0, cfrac, co2abs, ci !Variables for AGS
    real     :: Ammax, betaw, fstr, Am, Rdark, PAR, alphac, tempy, An, AGSa1, Dstar, gcco2 !Variables for AGS
    real     :: rsAgs, rsCO2, fw, Resp, wco2 !Variables for AGS
    real     :: PARdir, PARdif !Variables for 2leaf AGS
    real     :: MW_Air = 28.97
    real     :: MW_CO2 = 44

    real     :: sinbeta, kdrbl, kdf, kdr, ref, ref_dir
    real     :: iLAI, fSL
    real     :: PARdfU, PARdfD, PARdfT, PARdrU, PARdrD, PARdrT, dirPAR, difPAR
    real     :: HdfT, HdrT, dirH, Hshad, Hsun(nr_gauss), Fshad, Fsun, gshad, gsun
    real     :: Hleaf(nr_gauss+1), Fleaf(nr_gauss+1), gleaf(nr_gauss+1), Agl(nr_gauss+1)
    real     :: Fnet(nr_gauss), gnet(nr_gauss)
    real     :: minsinbeta = 1.e-10

    real     :: lthls_patch(xpatches,ypatches)
    integer  :: Npatch(xpatches,ypatches), SNpatch(xpatches,ypatches)

    real     :: local_wco2av
    real     :: local_Anav
    real     :: local_gcco2av
    real     :: local_Respav

    patchx = 0
    patchy = 0

    ! 1.X - Compute water content per layer
    do j = 2,j1
      do i = 2,i1
        phitot(i,j) = 0.0
        do k = 1, ksoilmax
          phitot(i,j) = phitot(i,j) + rootf(i,j,k) * max(phiw(i,j,k),phiwp)
        end do

        do k = 1, ksoilmax
          phifrac(i,j,k) = rootf(i,j,k) * max(phiw(i,j,k),phiwp) / phitot(i,j)
        end do
      end do
    end do

    thlsl = 0.0
    if(lhetero) then
      lthls_patch = 0.0
      Npatch      = 0
    endif

    wco2av       = 0.0
    Anav         = 0.0
    gcco2av      = 0.0
    Respav       = 0.0
    local_wco2av = 0.0
    local_Anav   = 0.0
    local_gcco2av= 0.0
    local_Respav = 0.0

    if (lrsAgs) then
      AnField    = 0.0
      gcco2Field = 0.0
      RespField  = 0.0
      wco2Field  = 0.0
      rsco2Field = 0.0
      fstrField  = 0.0
      ciField    = 0.0
      PARField   = 0.0
      if (lsplitleaf) then
        PARdirField= 0.0
        PARdifField= 0.0
      endif
    endif

    rk3coef = rdt / (4. - dble(rk3step))

    do j = 2, j1
      do i = 2, i1
        if(lhetero) then
          patchx = patchxnr(i)
          patchy = patchynr(j)
        endif

        ! 1.2   -   Calculate the skin temperature as the top boundary conditions for heat transport
        if(iradiation > 0) then
          if(iradiation == 1 .and. useMcICA) then
            if(rk3step == 1) then
              swdavn(i,j,2:nradtime) = swdavn(i,j,1:nradtime-1)
              swuavn(i,j,2:nradtime) = swuavn(i,j,1:nradtime-1)
              lwdavn(i,j,2:nradtime) = lwdavn(i,j,1:nradtime-1)
              lwuavn(i,j,2:nradtime) = lwuavn(i,j,1:nradtime-1)

              swdavn(i,j,1) = swd(i,j,1)
              swuavn(i,j,1) = swu(i,j,1)
              lwdavn(i,j,1) = lwd(i,j,1)
              lwuavn(i,j,1) = lwu(i,j,1)

            end if

            swdav = sum(swdavn(i,j,:)) / nradtime
            swuav = sum(swuavn(i,j,:)) / nradtime
            lwdav = sum(lwdavn(i,j,:)) / nradtime
            lwuav = sum(lwuavn(i,j,:)) / nradtime

            Qnet(i,j) = -(swdav + swuav + lwdav + lwuav)
          elseif(iradiation == irad_par .or. iradiation == 10) then !  Delta-eddington approach (2)  .or. rad_user (10)
            swdav      = -swd(i,j,1)
            Qnet(i,j)  = (swd(i,j,1) - swu(i,j,1) + lwd(i,j,1) - lwu(i,j,1))
          else ! simple radiation scheme and RRTMG
            Qnet(i,j) = -(swd(i,j,1) + swu(i,j,1) + lwd(i,j,1) + lwu(i,j,1))
            swdav     = swd(i,j,1)
          end if
        else
          if(lhetero) then
            Qnet(i,j) = Qnet_patch(patchx,patchy)
          else
            Qnet(i,j) = Qnetav
          endif
        end if

        ! 2.1   -   Calculate the surface resistance
        ! Stomatal opening as a function of incoming short wave radiation
        if (iradiation > 0) then
          f1  = 1. / min(1., (0.004 * max(0.,-swdav) + 0.05) / (0.81 * (0.004 * max(0.,-swdav) + 1.)))
        else
          f1  = 1.
        end if

        ! Soil moisture availability
        f2  = (phifc - phiwp) / (phitot(i,j) - phiwp)
        ! Prevent f2 becoming less than 1
        f2  = max(f2, 1.)
        ! Put upper boundary on f2 for cases with very dry soils
        f2  = min(1.e8, f2)

        ! Response of stomata to vapor deficit of atmosphere
        esat = 0.611e3 * exp(17.2694 * (thl0(i,j,1) - 273.16) / (thl0(i,j,1) - 35.86))
        if(lhetero) then
          e    = qt0(i,j,1) * ps_patch(patchx,patchy) / 0.622
        else
          e    = qt0(i,j,1) * ps / 0.622
        endif

        f3   = 1. / exp(-gD(i,j) * (esat - e) / 100.)

        ! Response to temperature
        exnera  = (presf(1) / pref0) ** (rd/cp)
        Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)
        f4      = 1./ (1. - 0.0016 * (298.0 - Tatm) ** 2.)

        rsveg(i,j)  = rsmin(i,j) / LAI(i,j) * f1 * f2 * f3! * f4 Not considered anymore

        ! 2.1a  - Recalculate vegetation resistance using AGS

        if (lrsAgs) then
          if (.not. linags) then !initialize AGS
            if(nsv .le. 0) then
              if (myid == 0) then
                print *, 'Problem in namoptions NAMSURFACE'
                print *, 'You enabled AGS (with switch lrsAgs), but there are no scalars (and thus no CO2 as needed for AGS) '
                stop 'ERROR: Problem in namoptions NAMSURFACE - AGS'
              endif
            endif
            if(lCHon) then !Chemistry is on
              if (CO2loc .lt. 1) then !No CO2 in chemistry
                if (myid == 0) then
                  print *, 'WARNING ::: There is no CO2 defined in the chemistry scheme'
                  print *, 'WARNING ::: Scalar 1 might be considered to be CO2 '
                  stop 'ERROR: Problem in namoptions NAMSURFACE - AGS'
                endif
                indCO2 = 1
              else  !CO2 present in chemistry
                indCO2 = CO2loc
              endif !Is CO2 present?
            else if (imicro==imicro_bulk) then !chemistry off and bulk_microphysics on
              if (myid==0) then
                print *,'WARNING ::: bulk microphysics and AGS are both ON'
                print *,'WARNING ::: Scalar 1 and 2 are considered qr and Nr,respectively for microphysics scheme'
                print *,'WARNING ::: Scalar 3 is considered CO2 for A-gs'
              endif
              indCO2 = 3
            else !Chemistry and bulk_micro are off
              if (myid == 0) then
                print *, 'WARNING ::: There is no CO2 defined due to the absence of a chemistry scheme'
                print *, 'WARNING ::: Scalar 1 is considered to be CO2 '
              endif
              indCO2 = 1
            endif !Is chemistry or bulk_micro on?
            linags = .true.
          endif !linags

          CO2ags  = svm(i,j,1,indCO2)/1000.0  !From ppb (usual DALES standard) to ppm

          ! Calculate surface resistances using the plant physiological (A-gs) model
          ! Calculate the CO2 compensation concentration
          CO2comp = rhof(1) * CO2comp298 * Q10CO2 ** (0.1 * ( thl0(i,j,1) - 298.0 ) )

          ! Calculate the mesophyll conductance
          gm       = gm298 * Q10gm ** (0.1 * ( thl0(i,j,1) - 298.0) ) / ( (1. + exp(0.3 * ( T1gm - thl0(i,j,1) ))) * (1. + exp(0.3 * (thl0(i,j,1) - T2gm))))
          gm       = gm / 1000   ! conversion from mm s-1 to m s-1

          ! calculate CO2 concentration inside the leaf (ci)
          fmin0    = gmin/nuco2q - (1.0/9.0) * gm
          fmin     = (-fmin0 + ( fmin0 ** 2.0 + 4 * gmin/nuco2q * gm ) ** (0.5)) / (2. * gm)

          esatsurf = 0.611e3 * exp(17.2694 * (tskin(i,j) - 273.16) / (tskin(i,j) - 35.86))
          Ds       = (esatsurf - e) / 1000.0 ! In kPa
          D0       = (f0 - fmin) / ad

          cfrac    = f0 * (1.0 - Ds/D0) + fmin * (Ds/D0)
          co2abs   = CO2ags * (MW_CO2/MW_Air) * rhof(1)

          if (lrelaxci) then
            if (ci_old_set) then
              ci_inf        = cfrac * (co2abs - CO2comp) + CO2comp
              ci            = ci_old(i,j) + min(kci*rk3coef, 1.0) * (ci_inf - ci_old(i,j))
              if (rk3step  == 3) then
                ci_old(i,j) = ci
              endif
            else
              ci            = cfrac * (co2abs - CO2comp) + CO2comp
              ci_old(i,j)   = ci
            endif
          else
            ci              = cfrac * (co2abs - CO2comp) + CO2comp
          endif

          ! Calculate maximal gross primary production in high light conditions (Ag)
          Ammax    = Ammax298 * Q10Am ** ( 0.1 * ( thl0(i,j,1) - 298.0) ) / ( (1.0 + exp(0.3 * ( T1Am - thl0(i,j,1) ))) * (1. + exp(0.3 * (thl0(i,j,1) - T2Am))) )

          ! Calculate the effect of soil moisture stress on gross assimilation rate
          betaw    = max(1.0e-3,min(1.0,(phitot(i,j)-phiwp)/(phifc-phiwp)))

          ! Calculate stress function
          fstr     = betaw

          ! Calculate gross assimilation rate (Am)
          Am       = Ammax * (1 - exp( -(gm * (ci - CO2comp) / Ammax) ) )

          Rdark    = (1.0/9) * Am

          !PAR      = 0.40 * max(0.1,-swdav * cveg(i,j))
          PAR      = 0.50 * max(0.1,abs(swdav)) !Increase PAR to 50 SW
          if (lsplitleaf) then
            PARdir   = 0.50 * max(0.1_field_r,abs(swdir(i,j,1)))
            PARdif   = 0.50 * max(0.1_field_r,abs(swdif(i,j,1)))
          endif

          ! Calculate the light use efficiency
          alphac   = alpha0 * (co2abs  - CO2comp) / (co2abs + 2 * CO2comp)

          if(lsplitleaf) then
            sinbeta  = max(zenith(xtime*3600 + rtimee,xday,xlat,xlon), minsinbeta)
            kdrbl    = 0.5 / sinbeta                                     ! Direct radiation extinction coefficient for black leaves
            kdf      = kdfbl * sqrt(1.0-sigma)
            kdr      = kdrbl * sqrt(1.0-sigma)
            ref      = (1.0 - sqrt(1.0-sigma)) / (1.0 + sqrt(1.0-sigma)) ! Reflection coefficient
            ref_dir  = 2 * ref / (1.0 + 1.6 * sinbeta)

            do itg = 1, nr_gauss ! loop over the different LAI locations
              iLAI   = LAI(i,j) * LAI_g(itg)   ! Integrated LAI between here and canopy top; Gaussian distributed
              fSL    = exp(-kdrbl * iLAI)      ! Fraction of sun-lit leaves

              PARdfD = PARdif * (1.0-ref)     * exp(-kdf * iLAI    )     ! Total downward PAR due to diffuse radiation at canopy top
              PARdrD = PARdir * (1.0-ref_dir) * exp(-kdr * iLAI    )     ! Total downward PAR due to direct radiation at canopy top
              PARdfU = PARdif * (1.0-ref)     * exp(-kdf * LAI(i,j)) * albedo(i,j) * (1.0-ref) * exp(-kdf * (LAI(i,j)-iLAI)) ! Total upward (reflected) PAR that originates as diffuse radiation
              PARdrU = PARdir * (1.0-ref_dir) * exp(-kdr * LAI(i,j)) * albedo(i,j) * (1.0-ref) * exp(-kdf * (LAI(i,j)-iLAI)) ! Total upward (reflected) PAR that originates as direct radiation
              PARdfT = PARdfD + PARdfU                                   ! Total PAR due to diffuse radiation at canopy top
              PARdrT = PARdrD + PARdrU                                   ! Total PAR due to direct radiation at canopy top

              dirPAR = (1.0-sigma) * PARdir * fSL                        ! Purely direct PAR (can only be downward)
              difPAR = PARdfT + PARdrT - dirPAR                          ! Total diffuse radiation

              HdfT   = kdf * PARdfD + kdf * PARdfU
              HdrT   = kdr * PARdrD + kdf * PARdrU
              dirH   = kdrbl * dirPAR
              Hshad  = HdfT + HdrT - dirH

              Hsun   = Hshad + angle_g * (1.0-sigma) * kdrbl * PARdir / sum(angle_g * weight_g)

              Hleaf(1)              = Hshad
              Hleaf(2:(nr_gauss+1)) = Hsun

              Agl    = fstr * (Am + Rdark) * (1 - exp(-alphac*Hleaf/(Am + Rdark)))
              gleaf  = gmin/nuco2q +  Agl/(co2abs-ci)
              !Fleaf  = -(co2abs - ci) / (ra(i,j) + 1.0 / gleaf)
              Fleaf  = Agl - Rdark

              Fshad  = Fleaf(1)
              Fsun   = sum(weight_g * Fleaf(2:(nr_gauss+1)))
              gshad  = gleaf(1)
              gsun   = sum(weight_g * gleaf(2:(nr_gauss+1)))

              Fnet(itg) = Fsun * fSL + Fshad * (1 - fSL)
              gnet(itg) = gsun * fSL + gshad * (1 - fSL)

            end do !itg

            An       = LAI(i,j) * sum(weight_g * Fnet)
            gc_inf   = LAI(i,j) * sum(weight_g * gnet)

          else !lsplitleaf

          ! Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An)
          AGSa1    = 1.0 / (1 - f0)
          Dstar    = D0 / (AGSa1 * (f0 - fmin))

          tempy    = alphac * Kx * PAR / (Am + Rdark)
          An       = (Am + Rdark) * (1 - 1.0 / (Kx * LAI(i,j)) * (E1(tempy * exp(-Kx*LAI(i,j))) - E1(tempy)))
          gc_inf    = LAI(i,j) * (gmin/nuco2q + AGSa1 * fstr * An / ((co2abs - CO2comp) * (1 + Ds / Dstar)))

          endif !lsplitleaf


          if (lrelaxgc) then
            if (gc_old_set) then
              gcco2       = gc_old(i,j) + min(kgc*rk3coef, 1.0) * (gc_inf - gc_old(i,j))
              if (rk3step ==3) then
                gc_old(i,j) = gcco2
              endif
            else
              gcco2 = gc_inf
              gc_old(i,j) = gcco2
            endif
          else
            gcco2 = gc_inf
          endif

          ! Calculate surface resistances for moisture and carbon dioxide
          rsAgs    = 1.0 / (1.6 * gcco2)
          rsCO2    = 1.0 / gcco2

          rsveg(i,j) = rsAgs

          ! Calculate net flux of CO2 into the plant (An)
          An       = - (co2abs - ci) / (ra(i,j) + rsCO2)

          ! CO2 soil respiraion surface flux
          fw       = Cw * wsmax / (phitot(i,j) + wsmin)

!          Resp     = R10 * (1 - fw)*(1+ustar(i,j)) * exp( Eact0 / (283.15 * 8.314) * (1.0 - 283.15 / ( thl0(i,j,1) )))
          Resp     = R10 * (1 - fw)* exp( Eact0 / (283.15 * 8.314) * (1.0 - 283.15 / ( tsoil(i,j,1) )))

          wco2     = (An + Resp) * (MW_Air/MW_CO2) * (1.0/rhof(1)) ! In ppm m/s

          CO2flux(i,j) = wco2 * 1000.0 ! In ppb m/s


          local_wco2av = local_wco2av + wco2
          local_Anav   = local_Anav   + An
          local_gcco2av= local_gcco2av   + gcco2
          local_Respav = local_Respav + Resp

          AnField   (i,j) = An
          gcco2Field(i,j) = gcco2
          RespField (i,j) = Resp
          wco2Field (i,j) = wco2
          rsco2Field(i,j) = rsCO2
          fstrField (i,j) = fstr
          ciField   (i,j) = ci
          PARField  (i,j) = PAR
          if (lsplitleaf)then
            PARdirField(i,j) = PARdir
            PARdifField(i,j) = PARdif
          endif

        endif !lrsAgs

        ! 2.2   - Calculate soil resistance based on ECMWF method

        f2  = (phifc - phiwp) / (phiw(i,j,1) - phiwp)
        f2  = max(f2, 1.)
        f2  = min(1.e8, f2)
        rssoil(i,j) = rssoilmin(i,j) * f2
        ! 1.1   -   Calculate the heat transport properties of the soil
        ! CvH I put it in the init function, as we don't have prognostic soil moisture at this stage

        ! CvH solve the surface temperature implicitly including variations in LWout
        if(rk3step == 1) then
          tskinm(i,j) = tskin(i,j)
          Wlm(i,j)    = Wl(i,j)
        end if

        if(lhetero) then
          exner   = (ps_patch(patchx,patchy) / pref0) ** (rd/cp)
        else
          exner   = (ps / pref0) ** (rd/cp)
        endif
        tsurfm  = tskinm(i,j) * exner

        esat    = 0.611e3 * exp(17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86))
        if(lhetero) then
          qsat    = 0.622 * esat / ps_patch(patchx,patchy)
        else
          qsat    = 0.622 * esat / ps
        endif
        desatdT = esat * (17.2694 / (tsurfm - 35.86) - 17.2694 * (tsurfm - 273.16) / (tsurfm - 35.86)**2.)
        if(lhetero) then
          dqsatdT = 0.622 * desatdT / ps_patch(patchx,patchy)
        else
          dqsatdT = 0.622 * desatdT / ps
        endif

        ! First, remove LWup from Qnet calculation
        Qnet(i,j) = Qnet(i,j) + boltz * tsurfm ** 4.

        ! Calculate coefficients for surface fluxes
        fH      = rhof(1) * cp / ra(i,j)

        ! Allow for dew fall
        if(qsat - qt0(i,j,1) < 0.) then
          rsveg(i,j)  = 0.
          rssoil(i,j) = 0.
        end if

        Wlmx      = LAI(i,j) * Wmax
        Wl(i,j)   = min(Wl(i,j), Wlmx)
        cliq(i,j) = Wl(i,j) / Wlmx

        fLEveg  = (1. - cliq(i,j)) * cveg(i,j) * rhof(1) * rlv / (ra(i,j) + rsveg(i,j))
        fLEsoil = (1. - cveg(i,j))             * rhof(1) * rlv / (ra(i,j) + rssoil(i,j))
        fLEliq  = cliq(i,j) * cveg(i,j)        * rhof(1) * rlv /  ra(i,j)

        fLE     = fLEveg + fLEsoil + fLEliq

        exnera  = (presf(1) / pref0) ** (rd/cp)
        Tatm    = exnera * thl0(i,j,1) + (rlv / cp) * ql0(i,j,1)

        Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. + fH * Tatm + fLE *&
        (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)
!\todo  Acoef   = Qnet(i,j) - boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 4. + fH * Tatm + fLE *&
!       (dqsatdT * tsurfm - qsat + qt0(i,j,1)) + lambdaskin(i,j) * tsoil(i,j,1)- fRs[t]*(1.0 - albedoav(i,j))*swdown
        Bcoef   = 4. * boltz * tsurfm ** 3. + fH + fLE * dqsatdT + lambdaskin(i,j)

        if (Cskin(i,j) == 0.) then
          tskin(i,j) = Acoef * Bcoef ** (-1.) / exner
        else
          tskin(i,j) = (1. + rk3coef / Cskin(i,j) * Bcoef) ** (-1.) * (tsurfm + rk3coef / Cskin(i,j) * Acoef) / exner
        end if

        Qnet(i,j)     = Qnet(i,j) - (boltz * tsurfm ** 4. + 4. * boltz * tsurfm ** 3. * (tskin(i,j) * exner - tsurfm))
        G0(i,j)       = lambdaskin(i,j) * ( tskin(i,j) * exner - tsoil(i,j,1) )
        LE(i,j)       = - fLE * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))

        LEveg         = - fLEveg  * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))
        LEsoil        = - fLEsoil * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))
        LEliq         = - fLEliq  * ( qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat))

        if(LE(i,j) == 0.) then
          rs(i,j)     = 1.e8
        else
          rs(i,j)     = - rhof(1) * rlv * (qt0(i,j,1) - (dqsatdT * (tskin(i,j) * exner - tsurfm) + qsat)) / LE(i,j) - ra(i,j)
        end if

        H(i,j)        = - fH  * ( Tatm - tskin(i,j) * exner )
        tskin(i,j)    = max(min(tskin(i,j),tskinm(i,j)+10.),tskinm(i,j)-10.)
        tendskin(i,j) = Cskin(i,j) * (tskin(i,j) - tskinm(i,j)) * exner / rk3coef

        ! In case of dew formation, allow all water to enter skin reservoir Wl
        if(qsat - qt0(i,j,1) < 0.) then
          Wl(i,j)       =  Wlm(i,j) - rk3coef * ((LEliq + LEsoil + LEveg) / (rhow * rlv))
        else
          Wl(i,j)       =  Wlm(i,j) - rk3coef * (LEliq / (rhow * rlv))
        end if

        thlsl = thlsl + tskin(i,j)
        if (lhetero) then
          lthls_patch(patchx,patchy) = lthls_patch(patchx,patchy) + tskin(i,j)
          Npatch(patchx,patchy)      = Npatch(patchx,patchy)      + 1
        endif

        ! Solve the soil
        if(rk3step == 1) then
          tsoilm(i,j,:) = tsoil(i,j,:)
          phiwm(i,j,:)  = phiw(i,j,:)
        end if

        ! Calculate the soil heat capacity and conductivity based on water content
        do k = 1, ksoilmax
          pCs(i,j,k)    = (1. - phi) * pCm + phiw(i,j,k) * pCw
          Ke            = log10(phiw(i,j,k) / phi) + 1.
          lambda(i,j,k) = Ke * (lambdasat - lambdadry) + lambdadry
        end do

        do k = 1, ksoilmax-1
          lambdah(i,j,k) = (lambda(i,j,k) * dzsoil(k+1) + lambda(i,j,k+1) * dzsoil(k)) / (dzsoil(k+1)+dzsoil(k))
        end do

        lambdah(i,j,ksoilmax) = lambda(i,j,ksoilmax)

        do k = 1, ksoilmax
          gammas(i,j,k)  = gammasat * (phiw(i,j,k) / phi) ** (2. * bc + 3.)
          lambdas(i,j,k) = bc * gammasat * (-1.) * psisat / phi * (phiw(i,j,k) / phi) ** (bc + 2.)
        end do

        do k = 1, ksoilmax-1
          lambdash(i,j,k) = (lambdas(i,j,k) * dzsoil(k+1) + lambdas(i,j,k+1) * dzsoil(k)) / (dzsoil(k+1)+dzsoil(k))
          gammash(i,j,k)  = (gammas(i,j,k)  * dzsoil(k+1) + gammas(i,j,k+1)  * dzsoil(k)) / (dzsoil(k+1)+dzsoil(k))
        end do

        lambdash(i,j,ksoilmax) = lambdas(i,j,ksoilmax)

        ! 1.4   -   Solve the diffusion equation for the heat transport
        tsoil(i,j,1) = tsoilm(i,j,1) + rk3coef / pCs(i,j,1) * ( lambdah(i,j,1) &
        * (tsoil(i,j,2) - tsoil(i,j,1)) / dzsoilh(1) + G0(i,j) ) / dzsoil(1)
        do k = 2, ksoilmax-1
          tsoil(i,j,k) = tsoilm(i,j,k) + rk3coef / pCs(i,j,k) * ( lambdah(i,j,k) * (tsoil(i,j,k+1) - tsoil(i,j,k))&
          / dzsoilh(k) - lambdah(i,j,k-1) * (tsoil(i,j,k) - tsoil(i,j,k-1)) / dzsoilh(k-1) ) / dzsoil(k)
        end do
        tsoil(i,j,ksoilmax) = tsoilm(i,j,ksoilmax) + rk3coef / pCs(i,j,ksoilmax) * ( lambda(i,j,ksoilmax) * &
        (tsoildeep(i,j) - tsoil(i,j,ksoilmax)) / dzsoil(ksoilmax) - lambdah(i,j,ksoilmax-1) * &
        (tsoil(i,j,ksoilmax) - tsoil(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) ) / dzsoil(ksoilmax)
        ! 1.5   -   Solve the diffusion equation for the moisture transport
        phiw(i,j,1) = phiwm(i,j,1) + rk3coef * ( lambdash(i,j,1) * (phiw(i,j,2) - phiw(i,j,1)) / &
        dzsoilh(1) - gammash(i,j,1) - (phifrac(i,j,1) * LEveg + LEsoil) / (rhow*rlv)) / dzsoil(1)
        do k = 2, ksoilmax-1
          phiw(i,j,k) = phiwm(i,j,k) + rk3coef * ( lambdash(i,j,k) * (phiw(i,j,k+1) - phiw(i,j,k)) / &
          dzsoilh(k) - gammash(i,j,k) - lambdash(i,j,k-1) * (phiw(i,j,k) - phiw(i,j,k-1)) / dzsoilh(k-1) + &
          gammash(i,j,k-1) - (phifrac(i,j,k) * LEveg) / (rhow*rlv)) / dzsoil(k)
        end do
        ! closed bottom for now
        phiw(i,j,ksoilmax) = phiwm(i,j,ksoilmax) + rk3coef * (- lambdash(i,j,ksoilmax-1) * &
        (phiw(i,j,ksoilmax) - phiw(i,j,ksoilmax-1)) / dzsoil(ksoilmax-1) + gammash(i,j,ksoilmax-1) &
        - (phifrac(i,j,ksoilmax) * LEveg) / (rhow*rlv) ) / dzsoil(ksoilmax)
      end do
    end do

    if (lrelaxgc .and. (.not. gc_old_set) ) then
      if (rk3step == 3) then
        gc_old_set = .true.
      endif
    endif

    if (lrelaxci .and. (.not. ci_old_set) ) then
      if (rk3step == 3) then
        ci_old_set = .true.
      endif
    endif

    call D_MPI_ALLREDUCE(local_wco2av , wco2av , 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(local_Anav   , Anav   , 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(local_gcco2av, gcco2av, 1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(local_Respav , Respav , 1, MPI_SUM, comm3d,mpierr)

    Anav   = Anav/ijtot
    gcco2av= gcco2av/ijtot
    wco2av = wco2av/ijtot
    Respav = Respav/ijtot


    call D_MPI_ALLREDUCE(thlsl, thls, 1, MPI_SUM, comm3d,mpierr)
    thls = thls / ijtot
    if (lhetero) then
      call D_MPI_ALLREDUCE(lthls_patch(1:xpatches,1:ypatches), thls_patch(1:xpatches,1:ypatches),&
      xpatches*ypatches, MPI_SUM, comm3d,mpierr)
      call D_MPI_ALLREDUCE(Npatch(1:xpatches,1:ypatches)     , SNpatch(1:xpatches,1:ypatches)   ,&
      xpatches*ypatches, MPI_SUM, comm3d,mpierr)
      thls_patch = thls_patch / SNpatch
    endif

    call qtsurf

  end subroutine do_lsm

end module modsurface
