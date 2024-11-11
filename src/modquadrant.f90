!> \file modquadrant.f90
!!  Calculates quadrant-hole statistics


!>
!!  Calculates quadrant-hole statistics
!!  Written to Q#quadrant.<expnr>
!! If netcdf is true, this module also writes to a NetCDF file: quadrant.<expnr>.nc
!!
!!  \author Huug Ouwersloot, MPIC
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
!  Copyright 1993-2014 Delft University of Technology, Wageningen University, Utrecht University, KNMI, Max Planck Institute for Chemistry
!
module modquadrant

use modglobal, only : longint

implicit none
private
PUBLIC :: initquadrant, quadrant, exitquadrant
save
  logical                                    :: lquadrant = .false.
  real                                       :: hole = 0.0
  real                                       :: dtav
  real                                       :: timeav = -1.0
  integer                                    :: iwind = 1 ! 1: u, 2: v, 3: U = sqrt(u^2+v^2)
  integer                                    :: nvar, klow, khigh, knr
  character(80)                              :: fname = 'quadrant.xxx.nc'
  integer                                    :: ncid,nrec = 0
  character(80),allocatable,dimension(:,:,:) :: ncname
  character(80),dimension(1,4)               :: tncname
  integer(kind=longint)                      :: idtav,itimeav,tnext,tnextwrite
  integer, parameter                         :: isamptot=4
  integer                                    :: isamp
  character(30),dimension(4)                 :: samplname,longsamplname

  real, allocatable, dimension(:,:)          :: nrsampl
  real, allocatable, dimension(:,:)          :: uavl,vavl,wavl,utotavl,thlavl,qtavl
  real, allocatable, dimension(:,:)          :: uvarl,vvarl,wvarl,utotvarl,thlvarl,qtvarl
  real, allocatable, dimension(:,:,:)        :: svavl,svvarl
  real, allocatable, dimension(:,:)          :: wuresl,wvresl,wthlresl,wqtresl
  real, allocatable, dimension(:,:)          :: wusubl,wvsubl,wthlsubl,wqtsubl
  real, allocatable, dimension(:,:,:)        :: wsvresl,wsvsubl
  real, allocatable, dimension(:,:)          :: thlqtcovl


contains
!> Initialization routine, reads namelists and inits variables
  subroutine initquadrant
    use modmpi,    only : comm3d,mpierr,myid,D_MPI_BCAST
    use modglobal, only : ladaptive, dtmax,ifnamopt,fname_options,kmax,   &
                           dtav_glob,btime,tres,cexpnr,ifoutput,nsv,lwarmstart,checknamelisterror
    use modstat_nc, only : lnetcdf,define_nc,ncinfo,open_nc,define_nc,ncinfo,nctiminfo,writestat_dims_q_nc
    implicit none

    integer      :: ierr,n
    character(3) :: csvname

    namelist/NAMquadrant/ &
    lquadrant,dtav,timeav,hole,iwind,klow,khigh

    klow   = 2
    khigh  = kmax
    dtav   = dtav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMquadrant,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMquadrant')
      write(6 ,NAMquadrant)
      close(ifnamopt)

      if (timeav .lt. 0.0) timeav = dtav
      if (klow .lt. 2)     klow   = 2
      if (khigh .gt. kmax) khigh  = kmax
    end if ! myid = 0

    call D_MPI_BCAST(lquadrant ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(timeav    ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(hole      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(iwind     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(klow      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(khigh     ,1,0,comm3d,mpierr)

    if(.not. lquadrant) return

    knr              = 1 + khigh - klow

    samplname    (1) = 'Q1'
    samplname    (2) = 'Q2'
    samplname    (3) = 'Q3'
    samplname    (4) = 'Q4'
    longsamplname(1) = 'Quadrant 1: u > 0, w > 0'
    longsamplname(2) = 'Quadrant 2: u < 0, w > 0'
    longsamplname(3) = 'Quadrant 3: u < 0, w < 0'
    longsamplname(4) = 'Quadrant 4: u > 0, w < 0'

    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))
    tnext            = idtav   + btime
    tnextwrite       = itimeav + btime

    if (abs(timeav/dtav-nint(timeav/dtav))>1e-4) then
      stop 'NAMquadrant: timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'NAMquadrant: dtav should be a integer multiple of dtmax'
    end if

    nvar   = 22 + 4 * nsv

    if (myid==0) allocate(ncname(nvar,4,4))

    allocate(nrsampl   (klow:khigh,    4))
    allocate(uavl      (klow:khigh,    4))
    allocate(vavl      (klow:khigh,    4))
    allocate(wavl      (klow:khigh,    4))
    allocate(utotavl   (klow:khigh,    4))
    allocate(thlavl    (klow:khigh,    4))
    allocate(qtavl     (klow:khigh,    4))
    allocate(uvarl     (klow:khigh,    4))
    allocate(vvarl     (klow:khigh,    4))
    allocate(wvarl     (klow:khigh,    4))
    allocate(utotvarl  (klow:khigh,    4))
    allocate(thlvarl   (klow:khigh,    4))
    allocate(qtvarl    (klow:khigh,    4))
    allocate(svavl     (klow:khigh,nsv,4))
    allocate(svvarl    (klow:khigh,nsv,4))
    allocate(wuresl    (klow:khigh,    4))
    allocate(wvresl    (klow:khigh,    4))
    allocate(wthlresl  (klow:khigh,    4))
    allocate(wqtresl   (klow:khigh,    4))
    allocate(wusubl    (klow:khigh,    4))
    allocate(wvsubl    (klow:khigh,    4))
    allocate(wthlsubl  (klow:khigh,    4))
    allocate(wqtsubl   (klow:khigh,    4))
    allocate(wsvresl   (klow:khigh,nsv,4))
    allocate(wsvsubl   (klow:khigh,nsv,4))
    allocate(thlqtcovl (klow:khigh,    4))

    !initialize variables
    nrsampl   = 0.0
    uavl      = 0.0
    vavl      = 0.0
    wavl      = 0.0
    utotavl   = 0.0
    thlavl    = 0.0
    qtavl     = 0.0
    uvarl     = 0.0
    vvarl     = 0.0
    wvarl     = 0.0
    utotvarl  = 0.0
    thlvarl   = 0.0
    qtvarl    = 0.0
    svavl     = 0.0
    svvarl    = 0.0
    wuresl    = 0.0
    wvresl    = 0.0
    wthlresl  = 0.0
    wqtresl   = 0.0
    wusubl    = 0.0
    wvsubl    = 0.0
    wthlsubl  = 0.0
    wqtsubl   = 0.0
    wsvresl   = 0.0
    wsvsubl   = 0.0
    thlqtcovl = 0.0


    if(myid==0 .and. .not. lwarmstart)then
      do isamp = 1,isamptot
        open (ifoutput,file=trim(samplname(isamp))//'quadrant.'//cexpnr,status='replace')
        close (ifoutput)
      enddo
    endif

    if (lnetcdf) then
      if (myid==0) then
        call nctiminfo(tncname(1,:))
        fname(10:12) = cexpnr
        call open_nc(fname,ncid,nrec,nq=knr)
        call define_nc(ncid,1,tncname)
        call writestat_dims_q_nc(ncid,klow,khigh)
        do isamp=1,isamptot
          call ncinfo(ncname(           1,:,isamp),'nrsamp_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'number of points','-','qt')
          call ncinfo(ncname(           2,:,isamp),'uavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'horizontal velocity (u)','m/s','qt')
          call ncinfo(ncname(           3,:,isamp),'vavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'horizontal velocity (v)','m/s','qt')
          call ncinfo(ncname(           4,:,isamp),'wavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'vertical velocity','m/s','qt')
          call ncinfo(ncname(           5,:,isamp),'Umagavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'magnitude horizontal velocity','m/s','qt')
          call ncinfo(ncname(           6,:,isamp),'thlavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'liquid water potential temperature','K','qt')
          call ncinfo(ncname(           7,:,isamp),'qtavg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'specific humidity','kg/kg','qt')
          call ncinfo(ncname(           8,:,isamp),'uvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of horizontal velocity (u)','m^2/s^2','qt')
          call ncinfo(ncname(           9,:,isamp),'vvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of horizontal velocity (v)','m^2/s^2','qt')
          call ncinfo(ncname(          10,:,isamp),'wvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of vertical velocity','m^2/s^2','qt')
          call ncinfo(ncname(          11,:,isamp),'Umagvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of magnitude horizontal velocity','m^2/s^2','qt')
          call ncinfo(ncname(          12,:,isamp),'thlvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of liquid water potential temperature','K^2','qt')
          call ncinfo(ncname(          13,:,isamp),'qtvar_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'variance of specific humidity','kg^2/kg^2','qt')
          do n=1,nsv
            write (csvname(1:3),'(i3.3)') n
            call ncinfo(ncname(      13+n,:,isamp),'sv'//csvname//'avg_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'scalar '//csvname//' mixing ratio','ppb','qt')
            call ncinfo(ncname(  nsv+13+n,:,isamp),'sv'//csvname//'var_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'scalar '//csvname//' variance','ppb^2','qt')
          end do ! nsv
          call ncinfo(ncname(    2*nsv+14,:,isamp),'uwr_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'resolved turb. momentum flux (uw)','m^2/s^2','qt')
          call ncinfo(ncname(    2*nsv+15,:,isamp),'vwr_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'resolved turb. momentum flux (vw)','m^2/s^2','qt')
          call ncinfo(ncname(    2*nsv+16,:,isamp),'wthlr_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'resolved turb. theta_l flux','K m/s','qt')
          call ncinfo(ncname(    2*nsv+17,:,isamp),'wqtr_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'resolved turb. moisture flux','kg/kg m/s','qt')
          call ncinfo(ncname(    2*nsv+18,:,isamp),'uws_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'subgrid turb. momentum flux (uw)','m^2/s^2','qt')
          call ncinfo(ncname(    2*nsv+19,:,isamp),'vws_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'subgrid turb. momentum flux (vw)','m^2/s^2','qt')
          call ncinfo(ncname(    2*nsv+20,:,isamp),'wthls_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'subgrid turb. theta_l flux','K m/s','qt')
          call ncinfo(ncname(    2*nsv+21,:,isamp),'wqts_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'subgrid turb. moisture flux','kg/kg m/s','qt')
          do n=1,nsv
            write (csvname(1:3),'(i3.3)') n
            call ncinfo(ncname(2*nsv+21+n,:,isamp),'wsv'//csvname//'r_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'resolved turb. flux of scalar '//csvname,'ppb m/s','qt')
            call ncinfo(ncname(3*nsv+21+n,:,isamp),'wsv'//csvname//'s_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'subgrid turb. flux of scalar '//csvname,'ppb m/s','qt')
          end do ! nsv
          call ncinfo(ncname(    4*nsv+22,:,isamp),'thlqt_'//samplname(isamp),&
               trim(longsamplname(isamp))//' '//'covariance between theta_l and q','K kg/kg','qt')

          call define_nc( ncid, nvar, ncname(:,:,isamp))
        end do ! isamp
      end if ! myid = 0

    end if ! lnetcdf


  end subroutine initquadrant
!> Cleans up after the run
  subroutine exitquadrant
    use modstat_nc, only : lnetcdf
    use modmpi,     only : myid
    implicit none

    if(.not. lquadrant) return

    deallocate(nrsampl                                  )
    deallocate(uavl,vavl,wavl,utotavl,thlavl,qtavl      )
    deallocate(uvarl,vvarl,wvarl,utotvarl,thlvarl,qtvarl)
    deallocate(svavl,svvarl                             )
    deallocate(wuresl,wvresl,wthlresl,wqtresl           )
    deallocate(wusubl,wvsubl,wthlsubl,wqtsubl           )
    deallocate(wsvresl,wsvsubl                          )
    deallocate(thlqtcovl                                )
    if (lnetcdf .and. myid==0) deallocate(ncname)

  end subroutine exitquadrant

!> General routine, does the timekeeping
  subroutine quadrant
    use modglobal, only : rk3step,timee,dt_lim
#if defined(_OPENACC)
    use modgpu, only: update_host
#endif
    implicit none
    if(.not. lquadrant) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
#if defined(_OPENACC)
      call update_host
#endif
      do isamp = 1,isamptot
        call doquadrant
      end do
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writequadrant
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  return
  end subroutine quadrant
!> Performs the actual quadrant
  subroutine doquadrant
    use modglobal,      only : imax,jmax,i1,j1,i2,j2,dzh,dzf,nsv,ijtot,cu,cv,dx,dy
    use modfields,      only : u0,v0,w0,thl0,qt0,sv0
    use modsubgriddata, only : ekh,ekm
    use modmpi,         only : slabsum
    implicit none

    logical, allocatable, dimension(:,:,:)   :: mask
    real,    allocatable, dimension(:,:,:)   :: uloc, vloc, utot, windfield, uwhole
    real,    allocatable, dimension(:)       :: windfav, uwholeav
    real,    allocatable, dimension(:,:,:)   :: uvar, vvar, wvar, utotvar
    real,    allocatable, dimension(:,:,:)   :: thloc, qtloc, thvar, qtvar
    real,    allocatable, dimension(:,:,:,:) :: svloc, svvar
    real,    allocatable, dimension(:,:,:)   :: wures, wvres, wthlres, wqtres
    real,    allocatable, dimension(:,:,:)   :: wusub, wvsub, wthlsub, wqtsub
    real,    allocatable, dimension(:,:,:)   :: ugrad, vgrad, thgrad, qtgrad
    real,    allocatable, dimension(:,:,:,:) :: wsvres, wsvsub, svgrad
    real,    allocatable, dimension(:,:,:)   :: thlqtcov
    real,    allocatable, dimension(:,:,:)   :: ekhhalf,ekmhalf

    integer :: k,km,n

    allocate(mask     (2:i1,2:j1,klow:khigh    ))
    allocate(uloc     (2:i1,2:j1,klow:khigh    ))
    allocate(vloc     (2:i1,2:j1,klow:khigh    ))
    allocate(utot     (2:i1,2:j1,klow:khigh    ))
    allocate(windfield(2:i1,2:j1,klow:khigh    ))
    allocate(uwhole   (2:i1,2:j1,klow:khigh    ))
    allocate(windfav  (          klow:khigh    ))
    allocate(uwholeav (          klow:khigh    ))
    allocate(uvar     (2:i1,2:j1,klow:khigh    ))
    allocate(vvar     (2:i1,2:j1,klow:khigh    ))
    allocate(wvar     (2:i1,2:j1,klow:khigh    ))
    allocate(utotvar  (2:i1,2:j1,klow:khigh    ))
    allocate(thloc    (2:i1,2:j1,klow:khigh    ))
    allocate(qtloc    (2:i1,2:j1,klow:khigh    ))
    allocate(thvar    (2:i1,2:j1,klow:khigh    ))
    allocate(qtvar    (2:i1,2:j1,klow:khigh    ))
    allocate(svloc    (2:i1,2:j1,klow:khigh,nsv))
    allocate(svvar    (2:i1,2:j1,klow:khigh,nsv))
    allocate(wures    (2:i1,2:j1,klow:khigh    ))
    allocate(wvres    (2:i1,2:j1,klow:khigh    ))
    allocate(wthlres  (2:i1,2:j1,klow:khigh    ))
    allocate(wqtres   (2:i1,2:j1,klow:khigh    ))
    allocate(wusub    (2:i1,2:j1,klow:khigh    ))
    allocate(wvsub    (2:i1,2:j1,klow:khigh    ))
    allocate(wthlsub  (2:i1,2:j1,klow:khigh    ))
    allocate(wqtsub   (2:i1,2:j1,klow:khigh    ))
    allocate(ugrad    (2:i1,2:j1,klow:khigh    ))
    allocate(vgrad    (2:i1,2:j1,klow:khigh    ))
    allocate(thgrad   (2:i1,2:j1,klow:khigh    ))
    allocate(qtgrad   (2:i1,2:j1,klow:khigh    ))
    allocate(wsvres   (2:i1,2:j1,klow:khigh,nsv))
    allocate(wsvsub   (2:i1,2:j1,klow:khigh,nsv))
    allocate(svgrad   (2:i1,2:j1,klow:khigh,nsv))
    allocate(thlqtcov (2:i1,2:j1,klow:khigh    ))
    allocate(ekhhalf  (2:i1,2:j1,klow:khigh    ))
    allocate(ekmhalf  (2:i1,2:j1,klow:khigh    ))

    mask      = .false.
    uloc      = 0.0
    vloc      = 0.0
    utot      = 0.0
    windfield = 0.0
    uwhole    = 0.0
    windfav   = 0.0
    uwholeav  = 0.0
    uvar      = 0.0
    vvar      = 0.0
    wvar      = 0.0
    utotvar   = 0.0
    thloc     = 0.0
    qtloc     = 0.0
    thvar     = 0.0
    qtvar     = 0.0
    svloc     = 0.0
    svvar     = 0.0
    wures     = 0.0
    wvres     = 0.0
    wthlres   = 0.0
    wqtres    = 0.0
    wusub     = 0.0
    wvsub     = 0.0
    wthlsub   = 0.0
    wqtsub    = 0.0
    ugrad     = 0.0
    vgrad     = 0.0
    thgrad    = 0.0
    qtgrad    = 0.0
    wsvres    = 0.0
    wsvsub    = 0.0
    svgrad    = 0.0
    thlqtcov  = 0.0
    ekhhalf   = 0.0
    ekmhalf   = 0.0

    do k=klow,khigh
      km              = k-1

      uloc(:,:,k)     = (dzf(k)*(u0(2:i1,2:j1,km)+u0(3:i2,2:j1,km))+dzf(km)*(u0(2:i1,2:j1,k)+u0(3:i2,2:j1,k)))/(4*dzh(k)) + cu
      vloc(:,:,k)     = (dzf(k)*(v0(2:i1,2:j1,km)+v0(2:i1,3:j2,km))+dzf(km)*(v0(2:i1,2:j1,k)+v0(2:i1,3:j2,k)))/(4*dzh(k)) + cv

      thloc(:,:,k)    = (dzf(k)*thl0(2:i1,2:j1,km  )+dzf(km)*thl0(2:i1,2:j1,k  ))/(2*dzh(k)) - 300.0   !To limit rounding errors
      qtloc(:,:,k)    = (dzf(k)* qt0(2:i1,2:j1,km  )+dzf(km)* qt0(2:i1,2:j1,k  ))/(2*dzh(k))
      svloc(:,:,k,:)  = (dzf(k)* sv0(2:i1,2:j1,km,:)+dzf(km)* sv0(2:i1,2:j1,k,:))/(2*dzh(k))

      ekhhalf(:,:,k)  = (dzf(k)* ekh(2:i1,2:j1,km  )+dzf(km)* ekh(2:i1,2:j1,k  ))/(2*dzh(k))
      ekmhalf(:,:,k)  = (dzf(k)* ekm(2:i1,2:j1,km  )+dzf(km)* ekm(2:i1,2:j1,k  ))/(2*dzh(k))

      ugrad(:,:,k)    = ((u0(2:i1,2:j1,k)+u0(3:i2,2:j1,k))-(u0(2:i1,2:j1,km)+u0(3:i2,2:j1,km)))/(2*dzh(k))
      vgrad(:,:,k)    = ((v0(2:i1,2:j1,k)+v0(2:i1,3:j2,k))-(v0(2:i1,2:j1,km)+v0(2:i1,3:j2,km)))/(2*dzh(k))

      thgrad(:,:,k)   = (thl0(2:i1,2:j1,k  )-thl0(2:i1,2:j1,km  ))/dzh(k)
      qtgrad(:,:,k)   = (qt0 (2:i1,2:j1,k  )-qt0 (2:i1,2:j1,km  ))/dzh(k)
      svgrad(:,:,k,:) = (sv0 (2:i1,2:j1,k,:)-sv0 (2:i1,2:j1,km,:))/dzh(k)
    end do ! klow - khigh

    utot      = sqrt(uloc**2+vloc**2)

    uvar      = uloc**2
    vvar      = vloc**2
    wvar      = w0(2:i1,2:j1,klow:khigh)**2
    utotvar   = uloc**2+vloc**2
    thvar     = thloc**2
    qtvar     = qtloc**2
    svvar     = svloc**2

    wures     = w0(2:i1,2:j1,klow:khigh)*uloc
    wvres     = w0(2:i1,2:j1,klow:khigh)*vloc
    wthlres   = w0(2:i1,2:j1,klow:khigh)*thloc
    wqtres    = w0(2:i1,2:j1,klow:khigh)*qtloc

    wusub     = - ekmhalf * (ugrad + (w0(3:i2,2:j1,klow:khigh)-w0(1:imax,2:j1  ,klow:khigh))/(2*dx))
    wvsub     = - ekmhalf * (vgrad + (w0(2:i1,3:j2,klow:khigh)-w0(2:i1  ,1:jmax,klow:khigh))/(2*dy))
    wthlsub   = - ekhhalf * thgrad
    wqtsub    = - ekhhalf * qtgrad

    do n=1,nsv
      wsvres(:,:,:,n) = w0(2:i1,2:j1,klow:khigh)*svloc(:,:,:,n)
      wsvsub(:,:,:,n) = - ekhhalf * svgrad(:,:,:,n)
    end do ! nsv

    thlqtcov  = thloc*qtloc

    select case (iwind)
      case (1)
        windfield = uloc
      case (2)
        windfield = vloc
      case (3)
        windfield = utot
    end select

    call slabsum(windfav,klow,khigh,windfield,2,i1,2,j1,klow,khigh,2,i1,2,j1,klow,khigh)
    windfav   = windfav / ijtot

!!  Quadrant 1: u' > 0, w' > 0
!!  Quadrant 2: u' < 0, w' > 0
!!  Quadrant 3: u' < 0, w' < 0
!!  Quadrant 4: u' > 0, w' < 0
    select case (isamp)
      case (1)
        do k=klow,khigh
          mask(:,:,k) = ( (w0(2:i1,2:j1,k) .gt. 0.0) .and. (windfield(2:i1,2:j1,k) .gt. windfav(k)) )
        end do
      case (2)
        do k=klow,khigh
          mask(:,:,k) = ( (w0(2:i1,2:j1,k) .gt. 0.0) .and. (windfield(2:i1,2:j1,k) .lt. windfav(k)) )
        end do
      case (3)
        do k=klow,khigh
          mask(:,:,k) = ( (w0(2:i1,2:j1,k) .lt. 0.0) .and. (windfield(2:i1,2:j1,k) .lt. windfav(k)) )
        end do
      case (4)
        do k=klow,khigh
          mask(:,:,k) = ( (w0(2:i1,2:j1,k) .lt. 0.0) .and. (windfield(2:i1,2:j1,k) .gt. windfav(k)) )
        end do
    end select

    if (hole .gt. 0.0) then
      do k=klow,khigh
        uwhole(:,:,k) = w0(2:i1,2:j1,k)*(windfield(2:i1,2:j1,k) - windfav(k))
      end do
      call slabsum(uwholeav,klow,khigh,uwhole,2,i1,2,j1,klow,khigh,2,i1,2,j1,klow,khigh)
      uwholeav        = uwholeav / ijtot
      do k=klow,khigh
        mask(:,:,k)   = mask(:,:,k) .and. ( abs(uwhole(:,:,k)).gt.(hole*abs(uwholeav(k))) )
      end do
    endif


    do k=klow,khigh
      nrsampl  (k,  isamp) = nrsampl  (k,  isamp) + count(mask    ( :  , :  ,k  ))
      uavl     (k,  isamp) = uavl     (k,  isamp) + sum  (uloc    ( :  , :  ,k  ),mask(:,:,k))
      vavl     (k,  isamp) = vavl     (k,  isamp) + sum  (vloc    ( :  , :  ,k  ),mask(:,:,k))
      wavl     (k,  isamp) = wavl     (k,  isamp) + sum  (w0      (2:i1,2:j1,k  ),mask(:,:,k))
      utotavl  (k,  isamp) = utotavl  (k,  isamp) + sum  (utot    ( :  , :  ,k  ),mask(:,:,k))
      thlavl   (k,  isamp) = thlavl   (k,  isamp) + sum  (thloc   ( :  , :  ,k  ),mask(:,:,k))
      qtavl    (k,  isamp) = qtavl    (k,  isamp) + sum  (qtloc   ( :  , :  ,k  ),mask(:,:,k))
      uvarl    (k,  isamp) = uvarl    (k,  isamp) + sum  (uvar    ( :  , :  ,k  ),mask(:,:,k))
      vvarl    (k,  isamp) = vvarl    (k,  isamp) + sum  (vvar    ( :  , :  ,k  ),mask(:,:,k))
      wvarl    (k,  isamp) = wvarl    (k,  isamp) + sum  (wvar    ( :  , :  ,k  ),mask(:,:,k))
      utotvarl (k,  isamp) = utotvarl (k,  isamp) + sum  (utotvar ( :  , :  ,k  ),mask(:,:,k))
      thlvarl  (k,  isamp) = thlvarl  (k,  isamp) + sum  (thvar   ( :  , :  ,k  ),mask(:,:,k))
      qtvarl   (k,  isamp) = qtvarl   (k,  isamp) + sum  (qtvar   ( :  , :  ,k  ),mask(:,:,k))
      do n=1,nsv
        svavl  (k,n,isamp) = svavl    (k,n,isamp) + sum  (svloc   ( :  , :  ,k,n),mask(:,:,k))
        svvarl (k,n,isamp) = svvarl   (k,n,isamp) + sum  (svvar   ( :  , :  ,k,n),mask(:,:,k))
      end do ! nsv
      wuresl   (k,  isamp) = wuresl   (k,  isamp) + sum  (wures   ( :  , :  ,k  ),mask(:,:,k))
      wvresl   (k,  isamp) = wvresl   (k,  isamp) + sum  (wvres   ( :  , :  ,k  ),mask(:,:,k))
      wthlresl (k,  isamp) = wthlresl (k,  isamp) + sum  (wthlres ( :  , :  ,k  ),mask(:,:,k))
      wqtresl  (k,  isamp) = wqtresl  (k,  isamp) + sum  (wqtres  ( :  , :  ,k  ),mask(:,:,k))
      wusubl   (k,  isamp) = wusubl   (k,  isamp) + sum  (wusub   ( :  , :  ,k  ),mask(:,:,k))
      wvsubl   (k,  isamp) = wvsubl   (k,  isamp) + sum  (wvsub   ( :  , :  ,k  ),mask(:,:,k))
      wthlsubl (k,  isamp) = wthlsubl (k,  isamp) + sum  (wthlsub ( :  , :  ,k  ),mask(:,:,k))
      wqtsubl  (k,  isamp) = wqtsubl  (k,  isamp) + sum  (wqtsub  ( :  , :  ,k  ),mask(:,:,k))
      do n=1,nsv
        wsvresl(k,n,isamp) = wsvresl  (k,n,isamp) + sum  (wsvres  ( :  , :  ,k,n),mask(:,:,k))
        wsvsubl(k,n,isamp) = wsvsubl  (k,n,isamp) + sum  (wsvsub  ( :  , :  ,k,n),mask(:,:,k))
      end do ! nsv
      thlqtcovl(k,  isamp) = thlqtcovl(k,  isamp) + sum  (thlqtcov( :  , :  ,k  ),mask(:,:,k))
    end do ! klow - khigh

    deallocate(mask                               )
    deallocate(uloc, vloc, utot, windfield,uwhole )
    deallocate(windfav, uwholeav                  )
    deallocate(uvar, vvar, wvar, utotvar          )
    deallocate(thloc, qtloc, thvar, qtvar         )
    deallocate(svloc, svvar                       )
    deallocate(wures, wvres, wthlres, wqtres      )
    deallocate(wusub, wvsub, wthlsub, wqtsub      )
    deallocate(ugrad, vgrad, thgrad, qtgrad       )
    deallocate(wsvres, wsvsub, svgrad             )
    deallocate(thlqtcov                           )
    deallocate(ekhhalf,ekmhalf                    )

  end subroutine doquadrant
!> Write the statistics to file
  subroutine writequadrant

    use modglobal, only : rtimee,zh,cexpnr,ifoutput,ijtot,nsv
    use modfields, only : presh
    use modmpi,    only : myid,comm3d,mpierr,mpi_sum,D_MPI_ALLREDUCE
    use modstat_nc, only: lnetcdf, writestat_nc,nc_fillvalue

    implicit none
    real, allocatable, dimension(:,:)          :: vars
    real, allocatable, dimension(:,:)          :: nrsamp
    real, allocatable, dimension(:,:)          :: uavg,vavg,wavg,utotavg,thlavg,qtavg
    real, allocatable, dimension(:,:)          :: uvar,vvar,wvar,utotvar,thlvar,qtvar
    real, allocatable, dimension(:,:,:)        :: svavg,svvar
    real, allocatable, dimension(:,:)          :: wures,wvres,wthlres,wqtres
    real, allocatable, dimension(:,:)          :: wusub,wvsub,wthlsub,wqtsub
    real, allocatable, dimension(:,:,:)        :: wsvres,wsvsub
    real, allocatable, dimension(:,:)          :: thlqtcov
    real, allocatable, dimension(:,:)          :: perc
    real, allocatable, dimension(:,:)          :: uavgnorm,vavgnorm,wavgnorm,utotavgnorm,thlavgnorm,qtavgnorm
    real, allocatable, dimension(:,:)          :: uvarnorm,vvarnorm,wvarnorm,utotvarnorm,thlvarnorm,qtvarnorm
    real, allocatable, dimension(:,:,:)        :: svavgnorm,svvarnorm
    real, allocatable, dimension(:,:)          :: wuresnorm,wvresnorm,wthlresnorm,wqtresnorm
    real, allocatable, dimension(:,:)          :: wusubnorm,wvsubnorm,wthlsubnorm,wqtsubnorm
    real, allocatable, dimension(:,:,:)        :: wsvresnorm,wsvsubnorm
    real, allocatable, dimension(:,:)          :: thlqtcovnorm

    integer                                    :: nsecs, nhrs, nminut, k, n
    integer                                    :: inorm

    character(3)                               :: csvname

    allocate(vars(klow:khigh,nvar))

    allocate(nrsamp    (klow:khigh,    4)                                                                                                                                     )
    allocate(uavg      (klow:khigh,    4),vavg      (klow:khigh,    4),wavg       (klow:khigh,4),utotavg    (klow:khigh,4),thlavg      (klow:khigh,4),qtavg    (klow:khigh,4) )
    allocate(uvar      (klow:khigh,    4),vvar      (klow:khigh,    4),wvar       (klow:khigh,4),utotvar    (klow:khigh,4),thlvar      (klow:khigh,4),qtvar    (klow:khigh,4) )
    allocate(svavg     (klow:khigh,nsv,4),svvar     (klow:khigh,nsv,4)                                                                                                        )
    allocate(wures     (klow:khigh,    4),wvres     (klow:khigh,    4),wthlres    (klow:khigh,4),wqtres     (klow:khigh,4),thlqtcov    (klow:khigh,4)                         )
    allocate(wusub     (klow:khigh,    4),wvsub     (klow:khigh,    4),wthlsub    (klow:khigh,4),wqtsub     (klow:khigh,4)                                                    )
    allocate(wsvres    (klow:khigh,nsv,4),wsvsub    (klow:khigh,nsv,4)                                                                                                        )

    allocate(perc      (klow:khigh    ,4)                                                                                                                                     )
    allocate(uavgnorm  (klow:khigh    ,4),vavgnorm  (klow:khigh    ,4),wavgnorm   (klow:khigh,4),utotavgnorm(klow:khigh,4),thlavgnorm  (klow:khigh,4),qtavgnorm(klow:khigh,4) )
    allocate(uvarnorm  (klow:khigh    ,4),vvarnorm  (klow:khigh    ,4),wvarnorm   (klow:khigh,4),utotvarnorm(klow:khigh,4),thlvarnorm  (klow:khigh,4),qtvarnorm(klow:khigh,4) )
    allocate(svavgnorm (klow:khigh,nsv,4),svvarnorm (klow:khigh,nsv,4)                                                                                                        )
    allocate(wuresnorm (klow:khigh    ,4),wvresnorm (klow:khigh    ,4),wthlresnorm(klow:khigh,4),wqtresnorm (klow:khigh,4),thlqtcovnorm(klow:khigh,4)                         )
    allocate(wusubnorm (klow:khigh    ,4),wvsubnorm (klow:khigh    ,4),wthlsubnorm(klow:khigh,4),wqtsubnorm (klow:khigh,4)                                                    )
    allocate(wsvresnorm(klow:khigh,nsv,4),wsvsubnorm(klow:khigh,nsv,4)                                                                                                        )

    nsecs   = nint(rtimee)
    nhrs    = int(nsecs/3600)
    nminut  = int(nsecs/60)-nhrs*60
    nsecs   = mod(nsecs,60)
    inorm   = nint(ijtot*timeav/dtav)

    call D_MPI_ALLREDUCE(nrsampl  ,nrsamp  ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(uavl     ,uavg    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(vavl     ,vavg    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wavl     ,wavg    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(utotavl  ,utotavg ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(thlavl   ,thlavg  ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(qtavl    ,qtavg   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(uvarl    ,uvar    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(vvarl    ,vvar    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wvarl    ,wvar    ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(utotvarl ,utotvar ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(thlvarl  ,thlvar  ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(qtvarl   ,qtvar   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(svavl    ,svavg   ,isamptot*knr*nsv,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(svvarl   ,svvar   ,isamptot*knr*nsv,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wuresl   ,wures   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wvresl   ,wvres   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlresl ,wthlres ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtresl  ,wqtres  ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wusubl   ,wusub   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wvsubl   ,wvsub   ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlsubl ,wthlsub ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtsubl  ,wqtsub  ,isamptot*knr    ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wsvresl  ,wsvres  ,isamptot*knr*nsv,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wsvsubl  ,wsvsub  ,isamptot*knr*nsv,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(thlqtcovl,thlqtcov,isamptot*knr    ,MPI_SUM,comm3d,mpierr)

!reset variables
    nrsampl   = 0.0
    uavl      = 0.0
    vavl      = 0.0
    wavl      = 0.0
    utotavl   = 0.0
    thlavl    = 0.0
    qtavl     = 0.0
    uvarl     = 0.0
    vvarl     = 0.0
    wvarl     = 0.0
    utotvarl  = 0.0
    thlvarl   = 0.0
    qtvarl    = 0.0
    svavl     = 0.0
    svvarl    = 0.0
    wuresl    = 0.0
    wvresl    = 0.0
    wthlresl  = 0.0
    wqtresl   = 0.0
    wusubl    = 0.0
    wvsubl    = 0.0
    wthlsubl  = 0.0
    wqtsubl   = 0.0
    wsvresl   = 0.0
    wsvsubl   = 0.0
    thlqtcovl = 0.0

    if (myid==0) then
      if (lnetcdf) then
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
      endif

      perc         = 0.0
      uavgnorm     = 0.0
      vavgnorm     = 0.0
      wavgnorm     = 0.0
      utotavgnorm  = 0.0
      thlavgnorm   = 0.0
      qtavgnorm    = 0.0
      uvarnorm     = 0.0
      vvarnorm     = 0.0
      wvarnorm     = 0.0
      utotvarnorm  = 0.0
      thlvarnorm   = 0.0
      qtvarnorm    = 0.0
      svavgnorm    = 0.0
      svvarnorm    = 0.0
      wuresnorm    = 0.0
      wvresnorm    = 0.0
      wthlresnorm  = 0.0
      wqtresnorm   = 0.0
      wusubnorm    = 0.0
      wvsubnorm    = 0.0
      wthlsubnorm  = 0.0
      wqtsubnorm   = 0.0
      wsvresnorm   = 0.0
      wsvsubnorm   = 0.0
      thlqtcovnorm = 0.0

!normalize variables
      perc         = nrsamp / inorm

      where (nrsamp .GT. 0)
        uavgnorm     = uavg     / nrsamp
        vavgnorm     = vavg     / nrsamp
        wavgnorm     = wavg     / nrsamp
        utotavgnorm  = utotavg  / nrsamp
        thlavgnorm   = thlavg   / nrsamp + 300.0
        qtavgnorm    = qtavg    / nrsamp
        uvarnorm     = uvar     / nrsamp - (uavg    / nrsamp)**2
        vvarnorm     = vvar     / nrsamp - (vavg    / nrsamp)**2
        wvarnorm     = wvar     / nrsamp - (wavg    / nrsamp)**2
        utotvarnorm  = utotvar  / nrsamp - (utotavg / nrsamp)**2
        thlvarnorm   = thlvar   / nrsamp - (thlavg  / nrsamp)**2
        qtvarnorm    = qtvar    / nrsamp - (qtavg   / nrsamp)**2
        wuresnorm    = wures    / nrsamp - (wavg    / nrsamp) * (uavg   / nrsamp)
        wvresnorm    = wvres    / nrsamp - (wavg    / nrsamp) * (vavg   / nrsamp)
        wthlresnorm  = wthlres  / nrsamp - (wavg    / nrsamp) * (thlavg / nrsamp)
        wqtresnorm   = wqtres   / nrsamp - (wavg    / nrsamp) * (qtavg  / nrsamp)
        wusubnorm    = wusub    / nrsamp
        wvsubnorm    = wvsub    / nrsamp
        wthlsubnorm  = wthlsub  / nrsamp
        wqtsubnorm   = wqtsub   / nrsamp
        thlqtcovnorm = thlqtcov / nrsamp - (thlavg  / nrsamp) * (qtavg  / nrsamp)
      end where

      do n=1,nsv
        where (nrsamp .GT. 0)
          svavgnorm (:,n,:) = svavg (:,n,:) / nrsamp
          svvarnorm (:,n,:) = svvar (:,n,:) / nrsamp - (svavg (:,n,:) / nrsamp)**2
          wsvresnorm(:,n,:) = wsvres(:,n,:) / nrsamp - (wavg          / nrsamp) * (svavg (:,n,:) / nrsamp)
          wsvsubnorm(:,n,:) = wsvsub(:,n,:) / nrsamp
        end where
      end do ! nsv


!write files
      do isamp = 1,isamptot
        open (ifoutput,file=trim(samplname(isamp))//'quadrant.'//cexpnr,position='append')
        write(ifoutput,'(//3A,/A,I5,A,I4,A,I2,A,I2,A)') &
          '#-------------------------- ',trim(longsamplname(isamp)),' ---------------------------'      &
          ,'#',nint(timeav),' --- AVERAGING TIMESTEP --- '      &
          ,nhrs,':',nminut,':',nsecs      &
          ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(A/4A)') &
           '#-------------------------------------------------------------------------------' &
           ,'# LEV   HEIGHT  PRESS.    SAMPLES    FRACTION ' &
           ,'       u_avg        v_avg        w_avg      |U|_avg  theta_l_avg       qt_avg        u_var ' &
           ,'       v_var        w_var      |U|_var  theta_l_var       qt_var  resolved wu   subgrid wu ' &
           ,' resolved wv   subgrid wv res. wthetal sub. wthetal  resolved wq   subgrid wq  cov(th_l,q)'
        do k=klow,khigh
          write(ifoutput,'(i5,1X,F8.2,1X,F7.2,1X,i10,1X,E11.4,21(1X,E12.5))') &
              k                      , &
              zh           (k)       , &
              presh        (k)/100.  , &
              nint(nrsamp  (k,isamp)), &
              perc         (k,isamp) , &
              uavgnorm     (k,isamp) , &
              vavgnorm     (k,isamp) , &
              wavgnorm     (k,isamp) , &
              utotavgnorm  (k,isamp) , &
              thlavgnorm   (k,isamp) , &
              qtavgnorm    (k,isamp) , &
              uvarnorm     (k,isamp) , &
              vvarnorm     (k,isamp) , &
              wvarnorm     (k,isamp) , &
              utotvarnorm  (k,isamp) , &
              thlvarnorm   (k,isamp) , &
              qtvarnorm    (k,isamp) , &
              wuresnorm    (k,isamp) , &
              wusubnorm    (k,isamp) , &
              wvresnorm    (k,isamp) , &
              wvsubnorm    (k,isamp) , &
              wthlresnorm  (k,isamp) , &
              wthlsubnorm  (k,isamp) , &
              wqtresnorm   (k,isamp) , &
              wqtsubnorm   (k,isamp) , &
              thlqtcovnorm (k,isamp)
        end do ! klow - khigh
        close(ifoutput)

        do n=1,nsv
          write (csvname(1:3),'(i3.3)') n
          open (ifoutput,file=trim(samplname(isamp))//'quadrantsv'//csvname//'.'//cexpnr,position='append')
          write(ifoutput,'(//3A,I4,A,/A,I5,A,I4,A,I2,A,I2,A)') &
            '#--------------------- ',trim(longsamplname(isamp)),' : scalar ',n,' ----------------------'      &
            ,'#',nint(timeav),' --- AVERAGING TIMESTEP --- '      &
            ,nhrs,':',nminut,':',nsecs      &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '

          write (ifoutput,'(A/2A)') &
             '#-------------------------------------------------------------------------------' &
             ,'# LEV   HEIGHT  PRESS.    SAMPLES    FRACTION ' &
             ,'     AVERAGE     VARIANCE RESOLV. FLUX   SUBG. FLUX'
          do k=klow,khigh
            write(ifoutput,'(i5,1X,F8.2,1X,F7.2,1X,i10,1X,E11.4,4(1X,E12.5))') &
                k                        , &
                zh           (k)         , &
                presh        (k)/100.    , &
                nint(nrsamp  (k,isamp))  , &
                perc         (k,isamp)   , &
                svavgnorm    (k,n,isamp) , &
                svvarnorm    (k,n,isamp) , &
                wsvresnorm   (k,n,isamp) , &
                wsvsubnorm   (k,n,isamp)
          end do ! klow - khigh
          close(ifoutput)
        end do ! nsv


        if (lnetcdf) then
          vars               = nc_fillvalue
          vars(:,1)          = nrsamp(:,isamp)

          where (nrsamp(:,isamp) .gt. 0)
            vars(:,       2) = uavgnorm     (:,isamp)
            vars(:,       3) = vavgnorm     (:,isamp)
            vars(:,       4) = wavgnorm     (:,isamp)
            vars(:,       5) = utotavgnorm  (:,isamp)
            vars(:,       6) = thlavgnorm   (:,isamp)
            vars(:,       7) = qtavgnorm    (:,isamp)
            vars(:,       8) = uvarnorm     (:,isamp)
            vars(:,       9) = vvarnorm     (:,isamp)
            vars(:,      10) = wvarnorm     (:,isamp)
            vars(:,      11) = utotvarnorm  (:,isamp)
            vars(:,      12) = thlvarnorm   (:,isamp)
            vars(:,      13) = qtvarnorm    (:,isamp)
            vars(:,2*nsv+14) = wuresnorm    (:,isamp)
            vars(:,2*nsv+15) = wvresnorm    (:,isamp)
            vars(:,2*nsv+16) = wthlresnorm  (:,isamp)
            vars(:,2*nsv+17) = wqtresnorm   (:,isamp)
            vars(:,2*nsv+18) = wusubnorm    (:,isamp)
            vars(:,2*nsv+19) = wvsubnorm    (:,isamp)
            vars(:,2*nsv+20) = wthlsubnorm  (:,isamp)
            vars(:,2*nsv+21) = wqtsubnorm   (:,isamp)
            vars(:,4*nsv+22) = thlqtcovnorm (:,isamp)
          end where

          do n=1,nsv
            where (nrsamp(:,isamp) .gt. 0)
              vars(:,      13+n) = svavgnorm (:,n,isamp)
              vars(:,  nsv+13+n) = svvarnorm (:,n,isamp)
              vars(:,2*nsv+21+n) = wsvresnorm(:,n,isamp)
              vars(:,3*nsv+21+n) = wsvsubnorm(:,n,isamp)
            end where
          end do ! nsv

          call writestat_nc(ncid,nvar,ncname(:,:,isamp),vars(klow:khigh,:),nrec,knr)
        end if ! lnetcdf

      end do ! isamp
    end if ! myid = 0

    deallocate(vars)

    deallocate(nrsamp                                                               )
    deallocate(uavg      ,vavg      ,wavg       ,utotavg    ,thlavg      ,qtavg     )
    deallocate(uvar      ,vvar      ,wvar       ,utotvar    ,thlvar      ,qtvar     )
    deallocate(svavg     ,svvar                                                     )
    deallocate(wures     ,wvres     ,wthlres    ,wqtres     ,thlqtcov               )
    deallocate(wusub     ,wvsub     ,wthlsub    ,wqtsub                             )
    deallocate(wsvres    ,wsvsub                                                    )

    deallocate(perc                                                                 )
    deallocate(uavgnorm  ,vavgnorm  ,wavgnorm   ,utotavgnorm,thlavgnorm  ,qtavgnorm )
    deallocate(uvarnorm  ,vvarnorm  ,wvarnorm   ,utotvarnorm,thlvarnorm  ,qtvarnorm )
    deallocate(svavgnorm ,svvarnorm                                                 )
    deallocate(wuresnorm ,wvresnorm ,wthlresnorm,wqtresnorm ,thlqtcovnorm           )
    deallocate(wusubnorm ,wvsubnorm ,wthlsubnorm,wqtsubnorm                         )
    deallocate(wsvresnorm,wsvsubnorm                                                )

  end subroutine writequadrant

end module modquadrant
