!> \file modsamptend.f90
!!  Calculates the tendencies of the main fields


!>
!!  Calculates the tendencies of the main fields
!>
!! Profiles of the individual terms of the prognostic equations.  Written to *tend.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Thijs Heus, MPI
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
module modsamptend
  use modglobal, only : longint
  use modsampdata
  implicit none
  private
  public :: initsamptend, samptend, exitsamptend, leibniztend, writesamptend
  save
!NetCDF variables
  character(80),allocatable,dimension(:,:,:) :: ncname
  character(80),dimension(1,4) :: tncname
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer,public,parameter :: tend_tot=1,tend_start=1,tend_hadv=2,tend_vadv=3,tend_subg=4,tend_force=5,tend_rad=6,&
                              tend_ls=7,tend_micro=8, tend_topbound=9,tend_pois=10,tend_addon=11, tend_coriolis=12, tend_totlb=13
  integer,parameter :: nrfields = 13
  character(20),dimension(10) :: samplname,longsamplname
  integer :: isamp,isamptot,ntsamp
  logical :: ldosamptendwrite = .false. !< write tendencies
  logical :: ldosamptendleib = .false. !< determine leibniz terms
  real :: lastrk3coef

  real, allocatable :: uptm(:,:,:),vptm(:,:,:),wptm(:,:,:),thlptm(:,:,:),qtptm(:,:,:),qrptm(:,:,:),nrptm(:,:,:)
  real, allocatable :: upav(:,:,:),vpav(:,:,:),wpav(:,:,:),thlpav(:,:,:),qtpav(:,:,:),qrpav(:,:,:),nrpav(:,:,:)
  real, allocatable :: upmn(:,:,:),vpmn(:,:,:),wpmn(:,:,:),thlpmn(:,:,:),qtpmn(:,:,:),qrpmn(:,:,:),nrpmn(:,:,:)
  real, allocatable :: ust(:,:),vst(:,:),wst(:,:),thlst(:,:),qtst(:,:),qrst(:,:),nrst(:,:)
  real, allocatable :: uwavr(:,:), vsavr(:,:), wavr(:,:), thlavr(:,:), qtavr(:,:), qravr(:,:), nravr(:,:)
  real, allocatable :: uwav(:,:), vsav(:,:), wav(:,:), thlav(:,:), qtav(:,:), qrav(:,:), nrav(:,:)
  real, allocatable :: thlwav(:,:), thlsav(:,:), qtwav(:,:), qtsav(:,:), qrwav(:,:), qrsav(:,:), nrwav(:,:), nrsav(:,:)
  real, allocatable :: uthlwav(:,:), vthlsav(:,:), uqtwav(:,:), vqtsav(:,:), uqrwav(:,:), vqrsav(:,:), unrwav(:,:), vnrsav(:,:)
  real, allocatable :: wthlav(:,:), wqtav(:,:), wqrav(:,:), wnrav(:,:)
  logical, allocatable :: tendmask(:,:,:,:)
  integer, allocatable :: nrsamptot(:,:),nrsamp(:,:),nrsamplast(:,:),nrsampnew(:,:)
  character(80) :: fname = 'samptend.xxx.nc'
  character(80) :: fname_block = 'samptend.xxxxyxxx.xxx.nc'
  integer :: ncid,nrec = 0

contains
!> Initialization routine, reads namelists and inits variables
subroutine initsamptend
    use modmpi,   only : myid
    use modglobal,only : dtmax,k1,ladaptive,&
                         btime,tres,j1,jh,i1,ih
    use modstat_nc, only : lnetcdf

    implicit none

    if (.not. lsamptend) return

    isamptot = 0
    if (lsampall) then
      isamptot = isamptot + 1
      samplname (isamptot) = 'all'
      longsamplname(isamptot) = 'All '
    endif
    if (lsampup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'upd'
      longsamplname(isamptot) = 'Updraft '
    end if
    if (lsampbuup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'buup'
      longsamplname(isamptot) = 'Buoyant Updraft '
    end if
    if (lsampcl) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cld'
      longsamplname(isamptot) = 'Cloud '
    end if
    if (lsampco) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cldcr'
      longsamplname(isamptot) = 'Cloud Core '
    end if
    if (lsampcldup) then
      isamptot = isamptot + 1
      samplname(isamptot) = 'cldup'
      longsamplname(isamptot) = 'Cloud Updraft '
    end if

    if(isamptot < 1) return
    if(.not.(lnetcdf)) return !only in netcdf at the moment

    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))
    tnext      = idtav   +btime
    tnextwrite = itimeav +btime


    if (abs(timeav/dtav-nint(timeav/dtav))>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    if (ltenddec .and. .not. lprocblock) then
      stop 'ltenddec is only intended to be used to complement processor-averaged budgets'
    end if

    if (lsamptendu) allocate (uptm(k1,nrfields,isamptot), upmn(k1,nrfields,isamptot), upav(k1,nrfields,isamptot), ust(k1,isamptot))
    if (lsamptendv) allocate (vptm(k1,nrfields,isamptot), vpmn(k1,nrfields,isamptot), vpav(k1,nrfields,isamptot), vst(k1,isamptot))
    if (lsamptendw) allocate (wptm(k1,nrfields,isamptot), wpmn(k1,nrfields,isamptot), wpav(k1,nrfields,isamptot), wst(k1,isamptot))
    if (lsamptendthl) allocate (thlptm(k1,nrfields,isamptot), thlpmn(k1,nrfields,isamptot), thlpav(k1,nrfields,isamptot), thlst(k1,isamptot))
    if (lsamptendqt) allocate (qtptm(k1,nrfields,isamptot), qtpmn(k1,nrfields,isamptot), qtpav(k1,nrfields,isamptot), qtst(k1,isamptot))
    if (lsamptendqr) allocate (qrptm(k1,nrfields,isamptot), qrpmn(k1,nrfields,isamptot), qrpav(k1,nrfields,isamptot), qrst(k1,isamptot))
    if (lsamptendnr) allocate (nrptm(k1,nrfields,isamptot), nrpmn(k1,nrfields,isamptot), nrpav(k1,nrfields,isamptot), nrst(k1,isamptot))

    ! Needed to decompose advective terms
    if (ltenddec) allocate (uwav(k1,isamptot), vsav(k1,isamptot), wav(k1,isamptot), uwavr(k1,isamptot), vsavr(k1,isamptot), & 
                            wavr(k1,isamptot), thlavr(k1,isamptot), qtavr(k1,isamptot), qravr(k1,isamptot), nravr(k1,isamptot))

    ! Only allocate these if you have the budget for a scalar and you want to decompose that scalar's budget
    if (ltenddec .and. lsamptendthl) allocate(thlav(k1,isamptot), thlwav(k1,isamptot), thlsav(k1,isamptot), wthlav(k1,isamptot), &
                                              uthlwav(k1,isamptot), vthlsav(k1,isamptot))
    if (ltenddec .and. lsamptendqt)  allocate(qtav(k1,isamptot), qtwav(k1,isamptot), qtsav(k1,isamptot), wqtav(k1,isamptot), &
                                              uqtwav(k1,isamptot), vqtsav(k1,isamptot))
    if (ltenddec .and. lsamptendqr)  allocate(qrav(k1,isamptot), qrwav(k1,isamptot), qrsav(k1,isamptot), wqrav(k1,isamptot), &
                                              uqrwav(k1,isamptot), vqrsav(k1,isamptot))
    if (ltenddec .and. lsamptendnr)  allocate(nrav(k1,isamptot), nrwav(k1,isamptot), nrsav(k1,isamptot), wnrav(k1,isamptot), &
                                              unrwav(k1,isamptot), vnrsav(k1,isamptot))

    allocate (tendmask(2-ih:i1+ih,2-jh:j1+jh,k1,isamptot))
    allocate (nrsamptot(k1,isamptot),nrsamp(k1,isamptot),nrsamplast(k1,isamptot),nrsampnew(k1,isamptot))
    
    ! Tendency terms
    if (lsamptendu) then
      uptm = 0.; upmn = 0.; upav = 0.; ust = 0.
    end if
    if (lsamptendv) then
     vptm = 0.; vpmn = 0.; vpav = 0.; vst = 0.
    end if
    if (lsamptendw) then
     wptm = 0.; wpmn = 0.; wpav = 0.; wst = 0.
    end if
    if (lsamptendthl) then
      thlptm = 0.; thlpmn = 0.; thlpav = 0.; thlst = 0.
    end if
    if (lsamptendqt) then
      qtptm = 0.; qtpmn = 0.; qtpav = 0.; qtst = 0.
    end if
    if (lsamptendqr) then 
      qrptm = 0.; qrpmn = 0.; qrpav = 0.; qrst = 0.
    end if
    if (lsamptendnr) then
     nrptm = 0.; nrpmn = 0.; nrpav = 0.; nrst = 0.
    end if

    ! Advection decomposition terms
    if (ltenddec) then
      uwav = 0.; vsav = 0; wav = 0.
    end if
    if (ltenddec .and. lsamptendthl) then
      thlav = 0.; thlwav = 0.; thlsav = 0.; wthlav = 0.; uthlwav = 0.; vthlsav = 0.
    end if
    if (ltenddec .and. lsamptendqt) then
      qtav = 0.; qtwav = 0.; qtsav = 0.; wqtav = 0.; uqtwav = 0.; vqtsav = 0.
    end if
    if (ltenddec .and. lsamptendqr) then
      qrav = 0.; qrwav = 0.; qrsav = 0.; wqrav = 0.; uqrwav = 0.; vqrsav = 0.
    end if
    if (ltenddec .and. lsamptendnr) then
      nrav = 0.; nrwav = 0.; nrsav = 0.; wnrav = 0.; unrwav = 0.; vnrsav = 0.
    end if

    tendmask=.false.
    nrsamp=0
    ntsamp=0
    nrsamptot=0
    nrsamplast=0
    nrsampnew=0

    if (lnetcdf) then
      if (lprocblock) then
        ! Each block writes its own budget to its own file
        call initnetcdf()
      else
        ! We average over all blocks and write to a single file
        if (myid==0) then
          call initnetcdf()
        end if
      end if      
    end if

  end subroutine initsamptend

  subroutine initnetcdf
  
    use modmpi,   only : cmyid
    use modglobal,only : cexpnr, kmax, output_prefix
    use modstat_nc, only : open_nc,define_nc,redefine_nc,ncinfo,nctiminfo,writestat_dims_nc
    implicit none
    logical :: proc = .true.
    character(80) :: dimst, dimsm, dimsu, dimsv
    integer :: nvar = 100 ! Current maximum number of fields
    integer :: ifield=0

    allocate(ncname(nvar,4,isamptot))

    call nctiminfo(tncname(1,:))
    if (lprocblock) then
      fname_block(10:17) = cmyid
      fname_block(19:21) = cexpnr
      dimst='tttt'
      dimsm='ttmt'
      dimsu='mttt'
      dimsv='tmtt'
      call open_nc(trim(output_prefix)//fname_block,ncid,nrec,n1=1,n2=1,n3=kmax)
      call define_nc( ncid,1,tncname)
      call writestat_dims_nc(ncid, 1, 1, proc)
    else
      fname(10:12) = cexpnr
      dimst='tt'
      dimsm='mt'
      call open_nc(fname,ncid,nrec,n3=kmax)
      call define_nc( ncid,1,tncname)
      call writestat_dims_nc(ncid)
    end if

    do isamp=1,isamptot
      ifield=0
      if (lsamptendu) then
        call ncinfo(ncname(ifield + 1,:,isamp),'utendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U horizontal advective tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 2,:,isamp),'utendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U vertical advective tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 3,:,isamp),'utenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U diffusive tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 4,:,isamp),'utendfor'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U tendency due to other forces','m/s^2',dimst)
        call ncinfo(ncname(ifield + 5,:,isamp),'utendcor','U coriolis tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 6,:,isamp),'utendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U large scale tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 7,:,isamp),'utendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U top boundary tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 8,:,isamp),'utendpois'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U pressure gradient tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 9,:,isamp),'utendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U in addons tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield + 10,:,isamp),'utendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U total tendency','m/s^2',dimst)
        ifield = ifield + 10
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'utendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U total tendency with leibniz terms','m/s^2',dimst)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendv) then
        call ncinfo(ncname(ifield +1,:,isamp),'vtendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V horizontal advective tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +2,:,isamp),'vtendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V vertical advective tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +3,:,isamp),'vtenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V diffusive tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +4,:,isamp),'vtendfor'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V tendency due to other forces','m/s^2',dimst)
        call ncinfo(ncname(ifield +5,:,isamp),'vtendcor'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V coriolis tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +6,:,isamp),'vtendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V large scale tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +7,:,isamp),'vtendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V top boundary tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +8,:,isamp),'vtendpois'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V pressure gradient tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +9,:,isamp),'vtendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V in addons tendency','m/s^2',dimst)
        call ncinfo(ncname(ifield +10,:,isamp),'vtendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V total tendency','m/s^2',dimst)
        ifield = ifield + 10
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'vtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V total tendency with leibniz terms','m/s^2',dimst)
          ifield = ifield +1
        end if
      end if
      if (lsamptendw) then
        call ncinfo(ncname(ifield +1,:,isamp),'wtendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W horizontal advective tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +2,:,isamp),'wtendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W vertical advective tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +3,:,isamp),'wtenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W diffusive tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +4,:,isamp),'wtendfor'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W tendency due to other forces','m/s^2',dimsm)
        call ncinfo(ncname(ifield +5,:,isamp),'wtendcor'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W coriolis tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +6,:,isamp),'wtendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W large scale tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +7,:,isamp),'wtendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W top boundary tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +8,:,isamp),'wtendpois'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W pressure gradient tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +9,:,isamp),'wtendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W in addons tendency','m/s^2',dimsm)
        call ncinfo(ncname(ifield +10,:,isamp),'wtendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W total tendency','m/s^2',dimsm)
        ifield = ifield + 10
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'wtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W total tendency with leibniz terms','m/s^2',dimst)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendthl) then
        call ncinfo(ncname(ifield +1,:,isamp),'thltendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l horizontal advective tendency','K/s',dimst)
        call ncinfo(ncname(ifield +2,:,isamp),'thltendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l vertical advective tendency','K/s',dimst)
        call ncinfo(ncname(ifield +3,:,isamp),'thltenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l diffusive tendency','K/s',dimst)
        call ncinfo(ncname(ifield +4,:,isamp),'thltendrad'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l radiative tendency','K/s',dimst)
        call ncinfo(ncname(ifield +5,:,isamp),'thltendmicro'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l microphysical tendency','K/s',dimst)
        call ncinfo(ncname(ifield +6,:,isamp),'thltendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l large scale tendency','K/s',dimst)
        call ncinfo(ncname(ifield +7,:,isamp),'thltendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l  top boundary tendency','K/s',dimst)
        call ncinfo(ncname(ifield +8,:,isamp),'thltendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l in addons tendency','K/s',dimst)
        call ncinfo(ncname(ifield +9,:,isamp),'thltendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'theta_l total tendency','K/s',dimst)
        ifield = ifield + 9
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'thltendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l total tendency with leibniz terms','K/s',dimst)
          ifield = ifield + 1
        end if
        if (ltenddec) then
          call ncinfo(ncname(ifield +1,:,isamp),'thlm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l block average','K',dimst)
          call ncinfo(ncname(ifield +2,:,isamp),'thlw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l west edge average','K',dimsu)
          call ncinfo(ncname(ifield +3,:,isamp),'thls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l south edge average','K',dimsv)
          call ncinfo(ncname(ifield +4,:,isamp),'uthlw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'u*theta_l west edge average','K m/s',dimsu)
          call ncinfo(ncname(ifield +5,:,isamp),'vthls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'v*theta_l south edge average','K m/s',dimsv)
          call ncinfo(ncname(ifield +6,:,isamp),'wthlm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'Vertical theta_l flux block average','K m/s',dimsm)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendqt) then
        call ncinfo(ncname(ifield +1,:,isamp),'qttendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content horizontal advective tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +2,:,isamp),'qttendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content vertical advective tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +3,:,isamp),'qttenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content diffusive tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +4,:,isamp),'qttendrad'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content radiative tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +5,:,isamp),'qttendmicro'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content microphysical tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +6,:,isamp),'qttendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content large scale tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +7,:,isamp),'qttendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content  top boundary tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +8,:,isamp),'qttendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content in addons tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +9,:,isamp),'qttendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'total water content total tendency','kg/kg/s',dimst)
        ifield = ifield + 9
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'qttendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content total tendency with leibniz terms','kg/kg/s',dimst)
          ifield = ifield + 1
        end if
        if (ltenddec) then
          call ncinfo(ncname(ifield +1,:,isamp),'qtm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content block average','kg/kg',dimst)
          call ncinfo(ncname(ifield +2,:,isamp),'qtw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content west edge average','kg/kg',dimsu)
          call ncinfo(ncname(ifield +3,:,isamp),'qts'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content south edge average','kg/kg',dimsv)
          call ncinfo(ncname(ifield +4,:,isamp),'uqtw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'u*qt west edge average','kg/kg m/s',dimsu)
          call ncinfo(ncname(ifield +5,:,isamp),'vqts'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'v*qt south edge average','kg/kg m/s',dimsv)
          call ncinfo(ncname(ifield +6,:,isamp),'wqtm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'Vertical total water content flux block average','kg/kg m/s',dimsm)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendqr) then
        call ncinfo(ncname(ifield +1,:,isamp),'qrtendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content horizontal advective tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +2,:,isamp),'qrtendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content vertical advective tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +3,:,isamp),'qrtenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content diffusive tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +4,:,isamp),'qrtendrad'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content radiative tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +5,:,isamp),'qrtendmicro'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'raom water content microphysical tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +6,:,isamp),'qrtendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content large scale tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +7,:,isamp),'qrtendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content  top boundary tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +8,:,isamp),'qrtendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content in addons tendency','kg/kg/s',dimst)
        call ncinfo(ncname(ifield +9,:,isamp),'qrtendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'rain water content total tendency','kg/kg/s',dimst)
        ifield = ifield + 9
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'qrtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'rain water content total tendency with leibniz terms','kg/kg/s',dimst)
          ifield = ifield + 1
        end if        
        if (ltenddec) then
          call ncinfo(ncname(ifield +1,:,isamp),'qrm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'rain water content block average','kg/kg',dimst)
          call ncinfo(ncname(ifield +2,:,isamp),'qrw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'rain water content west edge average','kg/kg',dimsu)
          call ncinfo(ncname(ifield +3,:,isamp),'qrs'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'rain water content south edge average','kg/kg',dimsv)
          call ncinfo(ncname(ifield +4,:,isamp),'uqrw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'u*qr west edge average','kg/kg m/s',dimsu)
          call ncinfo(ncname(ifield +5,:,isamp),'vqrs'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'v*qr south edge average','kg/kg m/s',dimsv)
          call ncinfo(ncname(ifield +6,:,isamp),'wqrm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'Vertical rain water content flux block average','kg/kg m/s',dimsm)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendnr) then
        call ncinfo(ncname(ifield +1,:,isamp),'nrtendhadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC horizontal advective tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +2,:,isamp),'nrtendvadv'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC vertical advective tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +3,:,isamp),'nrtenddif'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC diffusive tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +4,:,isamp),'nrtendrad'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC radiative tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +5,:,isamp),'nrtendmicro'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC microphysical tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +6,:,isamp),'nrtendls'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC large scale tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +7,:,isamp),'nrtendtop'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC top boundary tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +8,:,isamp),'nrtendaddon'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC addons tendency','/kg/s',dimst)
        call ncinfo(ncname(ifield +9,:,isamp),'nrtendtot'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'RDNC total tendency','/kg/s',dimst)
        ifield = ifield + 9
        if (ltendleib) then
          call ncinfo(ncname(ifield +1,:,isamp),'nrtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC total tendency with leibniz terms','/kg/s',dimst)
          ifield = ifield + 1
        end if
        if (ltenddec) then
          call ncinfo(ncname(ifield +1,:,isamp),'nrm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC block average','kg/kg',dimst)
          call ncinfo(ncname(ifield +2,:,isamp),'nrw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC west edge average','/kg',dimsu)
          call ncinfo(ncname(ifield +3,:,isamp),'nrs'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC south edge average','/kg',dimsv)
          call ncinfo(ncname(ifield +4,:,isamp),'unrw'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'u*nr west edge average','/kg m/s',dimsu)
          call ncinfo(ncname(ifield +5,:,isamp),'vnrs'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'v*nr south edge average','/kg m/s',dimsv)
          call ncinfo(ncname(ifield +6,:,isamp),'wnrm'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'Vertical RDNC flux block average','/kg m/s',dimsm)
          ifield = ifield + 6
        end if
      end if
      if (ltenddec) then
        call ncinfo(ncname(ifield +1,:,isamp),'wm'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'W block average','m/s',dimsm)
        call ncinfo(ncname(ifield +2,:,isamp),'uw'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'U western edge average','m/s',dimsu)
        call ncinfo(ncname(ifield +3,:,isamp),'vs'//samplname(isamp),&
        trim(longsamplname(isamp))//' '//'V southern edge average','m/s',dimsv)
        ifield = ifield + 3
      end if
      nvar = ifield ! total number of fields actually in use
      call define_nc(ncid,nvar,ncname(1:nvar,:,isamp)) ! Slice ncname up to nvar, which may have decreased
    enddo

  end subroutine initnetcdf

!> Performs the statistics, keeps track of what the tendencies were last time, and what they are this time.
  subroutine samptend(tendterm,firstterm,lastterm)
    use modmpi,    only : slabsum
    use modglobal, only : i1,j1,kmax,k1,ih,jh,&
                          cp,rv,rlv,rd,&
                          timee,rk3step,dt_lim,ijtot,nsv,rdt
    use modfields, only : up,vp,wp,thlp,qtp,svp,w0,thl0,ql0,exnf,qt0,u0,v0,sv0
    use modmicrodata, only : iqr,inr
    use modstat_nc, only : lnetcdf
    implicit none
    integer, intent(in)           :: tendterm !< name of the term to write down
    logical, intent(in), optional :: lastterm !< true if this is the last term of the equations; the write routine is entered.
    logical, intent(in), optional :: firstterm !< true if this is the first term of the equations
    real, allocatable, dimension(:,:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav
    integer :: i,j,k

    if (.not. lsamptend) return
    if(isamptot < 1) return
    if(.not.(lnetcdf)) return !only in netcdf at the moment
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if

    IF (present(firstterm)) THEN
    IF (firstterm) THEN
      ! Zero the process-contributed tendency fields (*ptm)
      ! Define sampling masks (tendmask) based on the current fields of buoyancy and w
      ! and sample the state variables over these isamptot masks
      nrsamplast=0
      tendmask=.false.
      if (lsamptendu) then 
        uptm = 0.; ust = 0.
      end if
      if (lsamptendv) then
        vptm = 0.; vst = 0.
      end if
      if (lsamptendw) then
        wptm = 0.; wst = 0.
      end if
      if (lsamptendthl) then
        thlptm = 0.; thlst = 0.
      end if
      if (lsamptendqt) then
        qtptm = 0.; qtst = 0.
      end if
      if (lsamptendqr) then
        qrptm = 0.; qrst = 0.
      end if
      if (lsamptendnr) then
        nrptm = 0.; nrst = 0.
      end if
      
      if (lsampco .or. lsampbuup) then
        allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
        allocate(thvav(k1))

        do k=1,k1
          thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                      *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
        enddo
        thvav = 0.0
        call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
        thvav = thvav/ijtot
      end if

      if (lsampup .or. lsampbuup .or. lsampcldup) then
        allocate(w0f(2-ih:i1+ih,2-jh:j1+jh,k1))
        do k=1,kmax
          w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
        end do
      end if

      do isamp=1,isamptot
        select case (samplname(isamp))
        case ('upd')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (w0f(i,j,k)>0.) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('buup')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if ((w0f(i,j,k)>0.0).and.(thv0(i,j,k) > thvav(k))) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cld')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cldcr')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0).and.thv0(i,j,k) > thvav(k)) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('cldup')
          do i=2,i1
          do j=2,j1
          do k=1,kmax
            if (ql0(i,j,k)>epsilon(1.0).and.w0f(i,j,k).gt.0.) then
                tendmask(i,j,k,isamp) = .true.
            endif
          enddo
          enddo
          enddo
        case ('all')
            tendmask(:,:,:,isamp)  = .true.
        end select
        do k=1,kmax
          nrsamp(k,isamp)= nrsamp(k,isamp)+count(tendmask(2:i1,2:j1,k,isamp))
          nrsamplast(k,isamp)= count(tendmask(2:i1,2:j1,k,isamp))
        end do
      enddo

      if (lsampco .or. lsampbuup) then
        deallocate(thv0)
        deallocate(thvav)
      end if
      if (lsampup .or. lsampbuup .or. lsampcldup) deallocate(w0f)

      do isamp=1,isamptot
      do k=1,kmax
        if (lsamptendu) ust(k,isamp) = sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if (lsamptendv) vst(k,isamp) = sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if (lsamptendw) wst(k,isamp) = sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if (lsamptendthl) thlst(k,isamp) = sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if (lsamptendqt) qtst(k,isamp) = sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if(nsv>1) then
          if (lsamptendqr) qrst(k,isamp) = sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
          if (lsamptendnr) nrst(k,isamp) = sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
        endif
      end do
      end do

      if (ltendleib) ldosamptendleib=.true.
      lastrk3coef = rdt / (4. - dble(rk3step))

    ENDIF
    ENDIF

    ! Strategy
    ! Call samptend after each new process has been added to the total tendency, which for each
    ! prognostic variable is stored in a single field which accumulates as processes are added.
    ! So each time a new process (tendterm) is added in program.f90, samptend is immediately called afterwards.
    ! Here, we keep track of both the total tendency before tendterm was added (uptm(k,tend_tot, isamp))
    ! and after tendterm was added (up), such that the difference is due to tendterm
    ! This is summed over cells in each level which satisfy tendmask, then stored in uptm(k,tendterm,isamp)

    if (lsamptendu) then
      do isamp=1,isamptot
      do k=1,kmax
        uptm(k,tendterm,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-uptm (k,tend_tot,isamp)
        uptm(k,tend_tot,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        upav(k,tendterm,isamp) = upav(k,tendterm,isamp)+uptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendv) then
      do isamp=1,isamptot
      do k=1,kmax
        vptm(k,tendterm,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-vptm (k,tend_tot,isamp)
        vptm(k,tend_tot,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        vpav(k,tendterm,isamp) = vpav(k,tendterm,isamp)+vptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendw) then
      do isamp=1,isamptot
      do k=1,kmax
        wptm(k,tendterm,isamp) = sum(wp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-wptm (k,tend_tot,isamp)
        wptm(k,tend_tot,isamp) = sum(wp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        wpav(k,tendterm,isamp) = wpav(k,tendterm,isamp)+wptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendthl) then
      do isamp=1,isamptot
      do k=1,kmax
        thlptm(k,tendterm,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-thlptm (k,tend_tot,isamp)
        thlptm(k,tend_tot,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        thlpav(k,tendterm,isamp) = thlpav(k,tendterm,isamp)+thlptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendqt) then
      do isamp=1,isamptot
      do k=1,kmax
        qtptm(k,tendterm,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-qtptm (k,tend_tot,isamp)
        qtptm(k,tend_tot,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        qtpav(k,tendterm,isamp) = qtpav(k,tendterm,isamp)+qtptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendqr) then
      do isamp=1,isamptot
      do k=1,kmax
        qrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))-qrptm (k,tend_tot,isamp)
        qrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
        qrpav(k,tendterm,isamp) = qrpav(k,tendterm,isamp)+qrptm(k,tendterm,isamp)
      end do
      end do
    end if
    if (lsamptendnr) then
      do isamp=1,isamptot
      do k=1,kmax
        nrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))-nrptm (k,tend_tot,isamp)
        nrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
        nrpav(k,tendterm,isamp) = nrpav(k,tendterm,isamp)+nrptm(k,tendterm,isamp)
      end do
      end do
    end if

    IF (present(lastterm)) THEN
    IF (lastterm) THEN
      ! Update the total tendency a final time
      if (lsamptendu) then
        do isamp=1,isamptot
        do k=1,kmax
          upav(k,tend_tot,isamp) = upav(k,tend_tot,isamp)+uptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendv) then
        do isamp=1,isamptot
        do k=1,kmax
          vpav(k,tend_tot,isamp) = vpav(k,tend_tot,isamp)+vptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendw) then
        do isamp=1,isamptot
        do k=1,kmax
          wpav(k,tend_tot,isamp) = wpav(k,tend_tot,isamp)+wptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendthl) then
        do isamp=1,isamptot
        do k=1,kmax
          thlpav(k,tend_tot,isamp) = thlpav(k,tend_tot,isamp)+thlptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendqt) then
        do isamp=1,isamptot
        do k=1,kmax
          qtpav(k,tend_tot,isamp) = qtpav(k,tend_tot,isamp)+qtptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendqr) then
        do isamp=1,isamptot
        do k=1,kmax
          qrpav(k,tend_tot,isamp) = qrpav(k,tend_tot,isamp)+qrptm(k,tend_tot,isamp)
        enddo
        enddo
      end if
      if (lsamptendnr) then
        do isamp=1,isamptot
        do k=1,kmax
          nrpav(k,tend_tot,isamp) = nrpav(k,tend_tot,isamp)+nrptm(k,tend_tot,isamp)
        enddo
        enddo
      end if

      ! Calculate the decomposed advection terms
      call decomposedadvtend

      tnext = tnext+idtav

      if (timee>=tnextwrite) then
        tnextwrite = tnextwrite+itimeav
        ldosamptendwrite = .true.
      end if
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
    END IF
    END IF

  end subroutine samptend

  subroutine decomposedadvtend
    ! Terms/variables needed to scale-decompose the advection terms, for each selected budget
    use modglobal, only : i1,imax,j1,jmax,k1,dzhi,dzf,dzh
    use modfields, only : w0,thl0,thl0h,qt0,qt0h,u0,v0,sv0
    use modmicrodata, only : iqr,inr
    use modsubgriddata, only : ekh
    use modsurfdata,only: thlflux,qtflux,svflux
    implicit none
    real :: ekhalf, thlhav, qthav, qrhav, qr0h, nrhav, nr0h
    real :: thlwavr, thlsavr, qtwavr, qtsavr
    real :: qrwavr, qrsavr, nrwavr, nrsavr
    real :: wthlavr, wqtavr, wqravr, wnravr
    integer :: i,j,k

    if(.not.(ltenddec)) return

      uwavr = 0.
      vsavr = 0.
      wavr = 0.
      thlavr = 0.
      thlwavr = 0.
      thlsavr = 0.
      qtavr = 0.
      qtwavr = 0.
      qtsavr = 0.
      qravr = 0.
      qrwavr = 0.
      qrsavr = 0.
      nravr = 0.
      nrwavr = 0.
      nrsavr = 0.

      do isamp=1,isamptot
      do k=1,k1
        ! uwav and vsav (averaged horizontal velocities over western/southern block edges)
        ! will not work for samp not equal all, because we do not sample object boundaries
        uwavr(k,isamp) = sum(u0(2   ,2:j1,k),tendmask(2   ,2:j1,k,isamp))/jmax
        vsavr(k,isamp) = sum(v0(2:i1,2   ,k),tendmask(2:i1,2   ,k,isamp))/imax
        wavr (k,isamp) = sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))/jmax/imax

        uwav(k,isamp) = uwav(k,isamp) + uwavr(k,isamp)
        vsav(k,isamp) = vsav(k,isamp) + vsavr(k,isamp)
        wav(k,isamp) = wav(k,isamp) + wavr(k,isamp)
      enddo
      enddo
      if (lsamptendthl) then
        do isamp=1,isamptot
        do k=1,k1
          ! Average over block, and over edges
          thlavr(k,isamp) = sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))/jmax/imax
          thlwavr = sum(thl0(2   ,2:j1,k),tendmask(2   ,2:j1,k,isamp))/jmax
          thlsavr = sum(thl0(2:i1,2   ,k),tendmask(2:i1,2   ,k,isamp))/imax

          thlav(k,isamp) = thlav(k,isamp) + thlavr(k,isamp)
          thlwav(k,isamp) = thlwav(k,isamp) + thlwavr
          thlsav(k,isamp) = thlsav(k,isamp) + thlsavr

          ! Horizontal flux on edges
          uthlwav(k,isamp) = uthlwav(k,isamp) + sum((u0(2   ,2:j1,k) - uwavr(k,isamp))*(thl0(2   ,2:j1,k) - thlwavr), tendmask(2,2:j1,k,isamp))/jmax
          vthlsav(k,isamp) = vthlsav(k,isamp) + sum((v0(2:i1,2   ,k) - vsavr(k,isamp))*(thl0(2:i1,2   ,k) - thlsavr), tendmask(2:i1,2,k,isamp))/imax

          ! Vertical flux (including the subgrid flux), excluding block-averaged contributions
          wthlavr = 0.
          if (k == 1) then
            do i=2,i1
            do j=2,j1
              wthlavr = wthlavr + thlflux(i,j)
            end do
            end do
          else
            thlhav = thlavr(k,isamp)*dzf(k-1)+thlavr(k-1,isamp)*dzf(k)/(2*dzh(k))
            do i=2,i1
            do j=2,j1
              ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
              wthlavr = wthlavr + (w0(i,j,k) - wavr(k,isamp))*(thl0h(i,j,k) - thlhav) &
                                - ekhalf*(thl0(i,j,k)-thl0(i,j,k-1))*dzhi(k)
            end do
            end do
          end if
          wthlav(k,isamp) = wthlav(k,isamp) + wthlavr/jmax/imax
        enddo
        enddo
      end if
      if (lsamptendqt) then
        do isamp=1,isamptot
        do k=1,k1
          qtavr(k,isamp)  = sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))/jmax/imax
          qtwavr = sum(qt0(2   ,2:j1,k),tendmask(2   ,2:j1,k,isamp))/jmax
          qtsavr = sum(qt0(2:i1,2   ,k),tendmask(2:i1,2   ,k,isamp))/imax

          qtav(k,isamp) = qtav(k,isamp) + qtavr(k,isamp)
          qtwav(k,isamp) = qtwav(k,isamp) + qtwavr
          qtsav(k,isamp) = qtsav(k,isamp) + qtsavr

          ! Horizontal flux on edges
          uqtwav(k,isamp) = uqtwav(k,isamp) + sum((u0(2   ,2:j1,k) - uwavr(k,isamp))*(qt0(2   ,2:j1,k) - qtwavr), tendmask(2,2:j1,k,isamp))/jmax
          vqtsav(k,isamp) = vqtsav(k,isamp) + sum((v0(2:i1,2   ,k) - vsavr(k,isamp))*(qt0(2:i1,2   ,k) - qtsavr), tendmask(2:i1,2,k,isamp))/imax

          ! Vertical flux (including the subgrid flux), excluding block-averaged contributions
          wqtavr = 0.
          if (k == 1) then
            do i=2,i1
            do j=2,j1
              wqtavr = wqtavr + qtflux(i,j)
            end do
            end do
          else
            qthav = qtavr(k,isamp)*dzf(k-1)+qtavr(k-1,isamp)*dzf(k)/(2*dzh(k))
            do i=2,i1
            do j=2,j1
              ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
              wqtavr = wqtavr + (w0(i,j,k) - wavr(k,isamp))*(qt0h(i,j,k) - qthav) &
                              - ekhalf*(qt0(i,j,k)-qt0(i,j,k-1))*dzhi(k)
            end do
            end do
          end if
          wqtav(k,isamp) = wqtav(k,isamp) + wqtavr/jmax/imax
        enddo
        enddo
      end if
      if (lsamptendqr) then
        do isamp=1,isamptot
        do k=1,k1
          qravr(k,isamp)  = sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))/jmax/imax
          qrwavr = sum(sv0(2   ,2:j1,k,iqr),tendmask(2   ,2:j1,k,isamp))/jmax
          qrsavr = sum(sv0(2:i1,2   ,k,iqr),tendmask(2:i1,2   ,k,isamp))/imax

          qrav(k,isamp) = qrav(k,isamp) + qravr(k,isamp)
          qrwav(k,isamp) = qrwav(k,isamp) + qrwavr
          qrsav(k,isamp) = qrsav(k,isamp) + qrsavr

          ! Horizontal flux on edges
          uqrwav(k,isamp) = uqrwav(k,isamp) + sum((u0(2   ,2:j1,k) - uwavr(k,isamp))*(sv0(2   ,2:j1,k,iqr) - qrwavr), tendmask(2,2:j1,k,isamp))/jmax
          vqrsav(k,isamp) = vqrsav(k,isamp) + sum((v0(2:i1,2   ,k) - vsavr(k,isamp))*(sv0(2:i1,2   ,k,iqr) - qrsavr), tendmask(2:i1,2,k,isamp))/imax

          ! Vertical flux (including the subgrid flux), excluding block-averaged contributions
          wqravr = 0.
          if (k == 1) then
            do i=2,i1
            do j=2,j1
              wqravr = wqravr + svflux(i,j,iqr)
            end do
            end do
          else
            qrhav = qravr(k,isamp)*dzf(k-1)+qravr(k-1,isamp)*dzf(k)/(2*dzh(k))
            do i=2,i1
            do j=2,j1
              ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
              qr0h = (sv0(i,j,k,iqr)*dzf(k-1)+sv0(i,j,k-1,iqr)*dzf(k))/(2*dzh(k))
              wqravr = wqravr + (w0(i,j,k) - wavr(k,isamp))*(qr0h - qrhav) &
                              - ekhalf*(sv0(i,j,k,iqr)-sv0(i,j,k-1,iqr))*dzhi(k)
            end do
            end do
          end if
          wqrav(k,isamp) = wqrav(k,isamp) + wqravr/jmax/imax
        enddo
        enddo
      end if
      if (lsamptendnr) then
        do isamp=1,isamptot
        do k=1,k1
          nravr(k,isamp)  = sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))/jmax/imax
          nrwavr = sum(sv0(2   ,2:j1,k,inr),tendmask(2   ,2:j1,k,isamp))/jmax
          nrsavr = sum(sv0(2:i1,2   ,k,inr),tendmask(2:i1,2   ,k,isamp))/imax

          nrav(k,isamp) = nrav(k,isamp) + nravr(k,isamp)
          nrwav(k,isamp) = nrwav(k,isamp) + nrwavr
          nrsav(k,isamp) = nrsav(k,isamp) + nrsavr

          ! Horizontal flux on edges
          unrwav(k,isamp) = unrwav(k,isamp) + sum((u0(2   ,2:j1,k) - uwavr(k,isamp))*(sv0(2   ,2:j1,k,inr) - nrwavr), tendmask(2,2:j1,k,isamp))/jmax
          vnrsav(k,isamp) = vnrsav(k,isamp) + sum((v0(2:i1,2   ,k) - vsavr(k,isamp))*(sv0(2:i1,2   ,k,inr) - nrsavr), tendmask(2:i1,2,k,isamp))/imax

          ! Vertical flux (including the subgrid flux), excluding block-averaged contributions
          wnravr = 0.
          if (k == 1) then
            do i=2,i1
            do j=2,j1
              wnravr = wnravr + svflux(i,j,inr)
            end do
            end do
          else
            nrhav = nravr(k,isamp)*dzf(k-1)+nravr(k-1,isamp)*dzf(k)/(2*dzh(k))
            do i=2,i1
            do j=2,j1
              ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
              nr0h = (sv0(i,j,k,inr)*dzf(k-1)+sv0(i,j,k-1,inr)*dzf(k))/(2*dzh(k))
              wnravr = wnravr + (w0(i,j,k) - wavr(k,isamp))*(nr0h - nrhav) &
                              - ekhalf*(sv0(i,j,k,inr)-sv0(i,j,k-1,inr))*dzhi(k)
            end do
            end do
          end if
          wnrav(k,isamp) = wnrav(k,isamp) + wnravr/jmax/imax
        enddo
        enddo
      end if
      ntsamp = ntsamp + 1

  end subroutine decomposedadvtend

  ! Make its own subroutine which, like leibniztend, does the decomposed variables
  ! Since (block) averaging is only done in writesamptend, but you need the actual averages to compute fluctuations,
  ! we actually need to compute averages for state variables here. Since leibniztend assumes division by nrsamp has not happened,
  ! this needs to occur after leibniztend would have been called, i.e. stick this subroutine after lebniztend in program.f90

  subroutine leibniztend
    ! If the sampling region changes in time, need to compute Leibniz terms to get total
    ! time derivative
    use modmpi,    only : slabsum
    use modglobal, only : i1,j1,kmax,k1,ih,jh,&
                          cp,rv,rlv,rd,&
                          ijtot
    use modfields, only : w0,thl0,ql0,exnf,qt0,u0,v0,sv0
    use modmicrodata, only : iqr,inr
    implicit none
    real, allocatable, dimension(:,:,:) :: w0f
    real, allocatable, dimension(:,:,:) :: thv0
    real, allocatable, dimension(:) :: thvav
    integer :: i,j,k


    if (.not. lsamptend) return
    if(.not.(ldosamptendleib)) return
    ldosamptendleib=.false.
    tendmask=.false.
    nrsampnew=0

    if (lsampco .or. lsampbuup) then
      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thvav(k1))

      do k=1,k1
        thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/ijtot
    end if

    if (lsampup .or. lsampbuup .or. lsampcldup) then
      allocate(w0f(2-ih:i1+ih,2-jh:j1+jh,k1))
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do
    end if

    do isamp=1,isamptot
      select case (samplname(isamp))
      case ('upd')
        do i=2,i1
        do j=2,j1
        do k=1,kmax
          if (w0f(i,j,k)>0.) then
              tendmask(i,j,k,isamp) = .true.
          endif
        enddo
        enddo
        enddo
      case ('buup')
        do i=2,i1
        do j=2,j1
        do k=1,kmax
          if ((w0f(i,j,k)>0.0).and.(thv0(i,j,k) > thvav(k))) then
              tendmask(i,j,k,isamp) = .true.
          endif
        enddo
        enddo
        enddo
      case ('cld')
        do i=2,i1
        do j=2,j1
        do k=1,kmax
          if (ql0(i,j,k)>epsilon(1.0)) then
              tendmask(i,j,k,isamp) = .true.
          endif
        enddo
        enddo
        enddo
      case ('cldcr')
        do i=2,i1
        do j=2,j1
        do k=1,kmax
          if (ql0(i,j,k)>epsilon(1.0).and.thv0(i,j,k) > thvav(k)) then
              tendmask(i,j,k,isamp) = .true.
          endif
        enddo
        enddo
        enddo
      case ('all')
          tendmask(:,:,:,isamp)  = .true.
      end select
      do k=1,kmax
        nrsampnew(k,isamp)= count(tendmask(2:i1,2:j1,k,isamp))
      end do
    enddo

    if (lsampco .or. lsampbuup) then
      deallocate(thv0)
      deallocate(thvav)
    end if

    do isamp=1,isamptot
    do k=1,kmax
      if((nrsampnew(k,isamp)>0).and.(nrsamplast(k,isamp)>0)) then! only do if sampling can be performed at both points in time
        if (lsamptendu) upav(k,tend_totlb,isamp) = upav(k,tend_totlb,isamp)+(ust(k,isamp)-&
        sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendv) vpav(k,tend_totlb,isamp) = vpav(k,tend_totlb,isamp)+(vst(k,isamp)-&
        sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendw) wpav(k,tend_totlb,isamp) = wpav(k,tend_totlb,isamp)+(wst(k,isamp)-&
        sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendthl) thlpav(k,tend_totlb,isamp) = thlpav(k,tend_totlb,isamp)+(thlst(k,isamp)-&
        sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendqt) qtpav(k,tend_totlb,isamp) = qtpav(k,tend_totlb,isamp)+(qtst(k,isamp)-&
        sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendqr) qrpav(k,tend_totlb,isamp) = qrpav(k,tend_totlb,isamp)+(qrst(k,isamp)-&
        sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
        if (lsamptendnr) nrpav(k,tend_totlb,isamp) = nrpav(k,tend_totlb,isamp)+(nrst(k,isamp)-&
        sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
      endif
    enddo
    enddo

    if (lsampup .or. lsampbuup .or. lsampcldup) deallocate(w0f)

  end subroutine leibniztend

!> Write the statistics to file
  subroutine writesamptend
    use modglobal, only : k1
    use modmpi,    only : myid,comm3d,mpierr,mpi_sum,D_MPI_ALLREDUCE
    use modstat_nc, only: lnetcdf
    implicit none
    integer :: field,k

    if (.not. lsamptend) return
    if (.not. ldosamptendwrite) return

    ldosamptendwrite=.false.

    if (lsamptendu) upmn = 0.
    if (lsamptendv) vpmn = 0.
    if (lsamptendw) wpmn = 0.
    if (lsamptendthl) thlpmn = 0.
    if (lsamptendqt) qtpmn = 0.
    if (lsamptendqr) qrpmn = 0.
    if (lsamptendnr) nrpmn = 0.
    nrsamptot=0

    if (lprocblock) then
      ! Keep sums on each processor
      nrsamptot = nrsamp
      if (lsamptendu) upmn = upav
      if (lsamptendv) vpmn = vpav
      if (lsamptendw) wpmn = wpav
      if (lsamptendthl) thlpmn = thlpav
      if (lsamptendqt) qtpmn = qtpav
      if (lsamptendqr) qrpmn = qrpav
      if (lsamptendnr) nrpmn = nrpav
    else
      ! Sum over all processors
      call D_MPI_ALLREDUCE(nrsamp   ,nrsamptot ,k1*isamptot         ,MPI_SUM,comm3d,mpierr)
      if (lsamptendu)   call D_MPI_ALLREDUCE(upav     ,upmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendv)   call D_MPI_ALLREDUCE(vpav     ,vpmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendw)   call D_MPI_ALLREDUCE(wpav     ,wpmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendthl) call D_MPI_ALLREDUCE(thlpav   ,thlpmn    ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendqt)  call D_MPI_ALLREDUCE(qtpav    ,qtpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendqr)  call D_MPI_ALLREDUCE(qrpav    ,qrpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
      if (lsamptendnr)  call D_MPI_ALLREDUCE(nrpav    ,nrpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    end if

    do field=1,nrfields
    do isamp=1,isamptot
    do k=1,k1
      if (nrsamptot(k,isamp)>0) then
        if (lsamptendu)   upmn  (k,field,isamp) = upmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendv)   vpmn  (k,field,isamp) = vpmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendw)   wpmn  (k,field,isamp) = wpmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendthl) thlpmn(k,field,isamp) = thlpmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendqt)  qtpmn (k,field,isamp) = qtpmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendqr)  qrpmn (k,field,isamp) = qrpmn (k,field,isamp)/nrsamptot(k,isamp)
        if (lsamptendnr)  nrpmn (k,field,isamp) = nrpmn (k,field,isamp)/nrsamptot(k,isamp)
      end if
    enddo
    enddo
    enddo

    ! Advection decomposition terms
    if (ltenddec) then
      do isamp=1,isamptot
      do k=1,k1
        uwav(k,isamp) = uwav(k,isamp)/ntsamp
        vsav(k,isamp) = vsav(k,isamp)/ntsamp
        wav(k,isamp) = wav(k,isamp)/ntsamp
        if (lsamptendthl) then
          thlav(k,isamp) = thlav(k,isamp)/ntsamp
          thlwav(k,isamp) = thlwav(k,isamp)/ntsamp
          thlsav(k,isamp) = thlsav(k,isamp)/ntsamp
          wthlav(k,isamp) = wthlav(k,isamp)/ntsamp
          uthlwav(k,isamp) = uthlwav(k,isamp)/ntsamp
          vthlsav(k,isamp) = vthlsav(k,isamp)/ntsamp
        end if
        if (lsamptendqt) then
          qtav(k,isamp) = qtav(k,isamp)/ntsamp
          qtwav(k,isamp) = qtwav(k,isamp)/ntsamp
          qtsav(k,isamp) = qtsav(k,isamp)/ntsamp
          wqtav(k,isamp) = wqtav(k,isamp)/ntsamp
          uqtwav(k,isamp) = uqtwav(k,isamp)/ntsamp
          vqtsav(k,isamp) = vqtsav(k,isamp)/ntsamp
        end if
        if (lsamptendqr) then
          qrav(k,isamp) = qrav(k,isamp)/ntsamp
          qrwav(k,isamp) = qrwav(k,isamp)/ntsamp
          qrsav(k,isamp) = qrsav(k,isamp)/ntsamp
          wqrav(k,isamp) = wqrav(k,isamp)/ntsamp
          uqrwav(k,isamp) = uqrwav(k,isamp)/ntsamp
          vqrsav(k,isamp) = vqrsav(k,isamp)/ntsamp
        end if
        if (lsamptendnr) then
          nrav(k,isamp) = nrav(k,isamp)/ntsamp
          nrwav(k,isamp) = nrwav(k,isamp)/ntsamp
          nrsav(k,isamp) = nrsav(k,isamp)/ntsamp
          wnrav(k,isamp) = wnrav(k,isamp)/ntsamp
          unrwav(k,isamp) = unrwav(k,isamp)/ntsamp
          vnrsav(k,isamp) = vnrsav(k,isamp)/ntsamp
        end if
      enddo
      enddo
    end if

    if (lnetcdf) then
      if (lprocblock) then
        ! Each block writes its own budget to its own file
        call writenetcdf_proc()
      else
        ! We have averaged over all blocks and write to a single file
        if (myid==0) then
          call writenetcdf()
        end if
      end if      
    end if

    if (lsamptendu) upav = 0.
    if (lsamptendv) vpav = 0.
    if (lsamptendw) wpav = 0.
    if (lsamptendthl) thlpav = 0.
    if (lsamptendqt) qtpav = 0.
    if (lsamptendqr) qrpav = 0.
    if (lsamptendnr) nrpav = 0.
    nrsamp = 0
    ntsamp = 0

    if (ltenddec) then
      uwav = 0.; vsav = 0; wav = 0.
    end if
    if (ltenddec .and. lsamptendthl) then
      thlav = 0.; thlwav = 0.; thlsav = 0.; wthlav = 0.; uthlwav = 0.; vthlsav = 0.
    end if
    if (ltenddec .and. lsamptendqt) then
      qtav = 0.; qtwav = 0.; qtsav = 0.; wqtav = 0.; uqtwav = 0.; vqtsav = 0.
    end if
    if (ltenddec .and. lsamptendqr) then
      qrav = 0.; qrwav = 0.; qrsav = 0.; wqrav = 0.; uqrwav = 0.; vqrsav = 0.
    end if
    if (ltenddec .and. lsamptendnr) then
      nrav = 0.; nrwav = 0.; nrsav = 0.; wnrav = 0.; unrwav = 0.; vnrsav = 0.
    end if

  end subroutine writesamptend


  subroutine writenetcdf

    use modstat_nc, only: writestat_nc
    use modglobal, only : kmax,k1,rtimee

    implicit none
    integer :: nvar = 73 ! Current maximum number of fields
    integer :: ifield=0
    real, allocatable :: vars(:,:)
    allocate(vars(1:k1,nvar))

    call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
    do isamp=1,isamptot
      ifield = 0
      vars=0.
      if (lsamptendu) then
        vars(:, ifield+1) = upmn(:,tend_hadv,isamp)
        vars(:, ifield+2) = upmn(:,tend_vadv,isamp)
        vars(:, ifield+3) = upmn(:,tend_subg,isamp)
        vars(:, ifield+4) = upmn(:,tend_force,isamp)
        vars(:, ifield+5) = upmn(:,tend_coriolis,isamp)
        vars(:, ifield+6) = upmn(:,tend_ls,isamp)
        vars(:, ifield+7) = upmn(:,tend_topbound,isamp)
        vars(:, ifield+8) = upmn(:,tend_pois,isamp)
        vars(:, ifield+9) = upmn(:,tend_addon,isamp)
        vars(:,ifield+10) = upmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(:,ifield+1) = upmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendv) then
        vars(:,ifield+1) = vpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = vpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = vpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = vpmn(:,tend_force,isamp)
        vars(:,ifield+5) = vpmn(:,tend_coriolis,isamp)
        vars(:,ifield+6) = vpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = vpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = vpmn(:,tend_pois,isamp)
        vars(:,ifield+9) = vpmn(:,tend_addon,isamp)
        vars(:,ifield+10) = vpmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(:,ifield+1) = vpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendw) then
        vars(:,ifield+1) = wpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = wpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = wpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = wpmn(:,tend_force,isamp)
        vars(:,ifield+5) = wpmn(:,tend_coriolis,isamp)
        vars(:,ifield+6) = wpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = wpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = wpmn(:,tend_pois,isamp)
        vars(:,ifield+9) = wpmn(:,tend_addon,isamp)
        vars(:,ifield+10) = wpmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(:,ifield+1) = wpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendthl) then
        vars(:,ifield+1) = thlpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = thlpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = thlpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = thlpmn(:,tend_rad,isamp)
        vars(:,ifield+5) = thlpmn(:,tend_micro,isamp)
        vars(:,ifield+6) = thlpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = thlpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = thlpmn(:,tend_addon,isamp)
        vars(:,ifield+9) = thlpmn(:,tend_tot,isamp)
        ifield = ifield+9
        if (ltendleib) then
          vars(:,ifield+1) = thlpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendqt) then
        vars(:,ifield+1) = qtpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = qtpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = qtpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = qtpmn(:,tend_rad,isamp)
        vars(:,ifield+5) = qtpmn(:,tend_micro,isamp)
        vars(:,ifield+6) = qtpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = qtpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = qtpmn(:,tend_addon,isamp)
        vars(:,ifield+9) = qtpmn(:,tend_tot,isamp)
        ifield = ifield + 9
        if (ltendleib) then
          vars(:,ifield+1) = qtpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendqr) then
        vars(:,ifield+1) = qrpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = qrpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = qrpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = qrpmn(:,tend_rad,isamp)
        vars(:,ifield+5) = qrpmn(:,tend_micro,isamp)
        vars(:,ifield+6) = qrpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = qrpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = qrpmn(:,tend_addon,isamp)
        vars(:,ifield+9) = qrpmn(:,tend_tot,isamp)
        ifield = ifield + 9
        if (ltendleib) then
          vars(:,ifield+1) = qrpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
      if (lsamptendnr) then
        vars(:,ifield+1) = nrpmn(:,tend_hadv,isamp)
        vars(:,ifield+2) = nrpmn(:,tend_vadv,isamp)
        vars(:,ifield+3) = nrpmn(:,tend_subg,isamp)
        vars(:,ifield+4) = nrpmn(:,tend_rad,isamp)
        vars(:,ifield+5) = nrpmn(:,tend_micro,isamp)
        vars(:,ifield+6) = nrpmn(:,tend_ls,isamp)
        vars(:,ifield+7) = nrpmn(:,tend_topbound,isamp)
        vars(:,ifield+8) = nrpmn(:,tend_addon,isamp)
        vars(:,ifield+9) = nrpmn(:,tend_tot,isamp)
        ifield = ifield + 9
        if (ltendleib) then
          vars(:,ifield+1) = nrpmn(:,tend_totlb,isamp)
          ifield = ifield + 1
        end if
      end if
    nvar = ifield
    call writestat_nc(ncid,nvar,ncname(:,:,isamp),vars(1:kmax,:),nrec,kmax)
    enddo

    deallocate(vars)

  end subroutine writenetcdf

  subroutine writenetcdf_proc

    use modstat_nc, only: writestat_nc
    use modglobal, only : kmax,k1,rtimee
    
    implicit none
    integer :: nvar = 73 ! Current maximum number of fields
    integer :: ifield=0
    real, allocatable :: vars(:,:,:,:)
    allocate(vars(1,1,1:k1,nvar))
    
    call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
    do isamp=1,isamptot
      ifield = 0
      vars=0.
      if (lsamptendu) then
        vars(1, 1, :, ifield+1) = upmn(:,tend_hadv,isamp)
        vars(1, 1, :, ifield+2) = upmn(:,tend_vadv,isamp)
        vars(1, 1, :, ifield+3) = upmn(:,tend_subg,isamp)
        vars(1, 1, :, ifield+4) = upmn(:,tend_force,isamp)
        vars(1, 1, :, ifield+5) = upmn(:,tend_coriolis,isamp)
        vars(1, 1, :, ifield+6) = upmn(:,tend_ls,isamp)
        vars(1, 1, :, ifield+7) = upmn(:,tend_topbound,isamp)
        vars(1, 1, :, ifield+8) = upmn(:,tend_pois,isamp)
        vars(1, 1, :, ifield+9) = upmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+10) = upmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = upmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
      end if
      if (lsamptendv) then
        vars(1, 1, :,ifield+1) = vpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = vpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = vpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = vpmn(:,tend_force,isamp)
        vars(1, 1, :,ifield+5) = vpmn(:,tend_coriolis,isamp)
        vars(1, 1, :,ifield+6) = vpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = vpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = vpmn(:,tend_pois,isamp)
        vars(1, 1, :,ifield+9) = vpmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+10) = vpmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = vpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
      end if
      if (lsamptendw) then
        vars(1, 1, :,ifield+1) = wpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = wpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = wpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = wpmn(:,tend_force,isamp)
        vars(1, 1, :,ifield+5) = wpmn(:,tend_coriolis,isamp)
        vars(1, 1, :,ifield+6) = wpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = wpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = wpmn(:,tend_pois,isamp)
        vars(1, 1, :,ifield+9) = wpmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+10) = wpmn(:,tend_tot,isamp)
        ifield = ifield+10
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = wpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
      end if
      if (lsamptendthl) then
        vars(1, 1, :,ifield+1) = thlpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = thlpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = thlpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = thlpmn(:,tend_rad,isamp)
        vars(1, 1, :,ifield+5) = thlpmn(:,tend_micro,isamp)
        vars(1, 1, :,ifield+6) = thlpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = thlpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = thlpmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+9) = thlpmn(:,tend_tot,isamp)
        ifield = ifield+9
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = thlpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
        if (ltenddec) then
          vars(1, 1, :,ifield+1) = thlav(:,isamp)
          vars(1, 1, :,ifield+2) = thlwav(:,isamp)
          vars(1, 1, :,ifield+3) = thlsav(:,isamp)
          vars(1, 1, :,ifield+4) = uthlwav(:,isamp)
          vars(1, 1, :,ifield+5) = vthlsav(:,isamp)
          vars(1, 1, :,ifield+6) = wthlav(:,isamp)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendqt) then
        vars(1, 1, :,ifield+1) = qtpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = qtpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = qtpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = qtpmn(:,tend_rad,isamp)
        vars(1, 1, :,ifield+5) = qtpmn(:,tend_micro,isamp)
        vars(1, 1, :,ifield+6) = qtpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = qtpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = qtpmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+9) = qtpmn(:,tend_tot,isamp)
        ifield = ifield+9
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = qtpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if        
        if (ltenddec) then
          vars(1, 1, :,ifield+1) = qtav(:,isamp)
          vars(1, 1, :,ifield+2) = qtwav(:,isamp)
          vars(1, 1, :,ifield+3) = qtsav(:,isamp)
          vars(1, 1, :,ifield+4) = uqtwav(:,isamp)
          vars(1, 1, :,ifield+5) = vqtsav(:,isamp)
          vars(1, 1, :,ifield+6) = wqtav(:,isamp)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendqr) then
        vars(1, 1, :,ifield+1) = qrpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = qrpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = qrpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = qrpmn(:,tend_rad,isamp)
        vars(1, 1, :,ifield+5) = qrpmn(:,tend_micro,isamp)
        vars(1, 1, :,ifield+6) = qrpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = qrpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = qrpmn(:,tend_addon,isamp)
        ifield = ifield+9
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = qrpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
        if (ltenddec) then
          vars(1, 1, :,ifield+1) = qrav(:,isamp)
          vars(1, 1, :,ifield+2) = qrwav(:,isamp)
          vars(1, 1, :,ifield+3) = qrsav(:,isamp)
          vars(1, 1, :,ifield+4) = uqrwav(:,isamp)
          vars(1, 1, :,ifield+5) = vqrsav(:,isamp)
          vars(1, 1, :,ifield+6) = wqrav(:,isamp)
          ifield = ifield + 6
        end if
      end if
      if (lsamptendnr) then
        vars(1, 1, :,ifield+1) = nrpmn(:,tend_hadv,isamp)
        vars(1, 1, :,ifield+2) = nrpmn(:,tend_vadv,isamp)
        vars(1, 1, :,ifield+3) = nrpmn(:,tend_subg,isamp)
        vars(1, 1, :,ifield+4) = nrpmn(:,tend_rad,isamp)
        vars(1, 1, :,ifield+5) = nrpmn(:,tend_micro,isamp)
        vars(1, 1, :,ifield+6) = nrpmn(:,tend_ls,isamp)
        vars(1, 1, :,ifield+7) = nrpmn(:,tend_topbound,isamp)
        vars(1, 1, :,ifield+8) = nrpmn(:,tend_addon,isamp)
        vars(1, 1, :,ifield+9) = nrpmn(:,tend_tot,isamp)
        ifield = ifield+9
        if (ltendleib) then
          vars(1, 1, :,ifield+1) = nrpmn(:,tend_totlb,isamp)
          ifield = ifield+1
        end if
        if (ltenddec) then
          vars(1, 1, :,ifield+1) = nrav(:,isamp)
          vars(1, 1, :,ifield+2) = nrwav(:,isamp)
          vars(1, 1, :,ifield+3) = nrsav(:,isamp)
          vars(1, 1, :,ifield+4) = unrwav(:,isamp)
          vars(1, 1, :,ifield+5) = vnrsav(:,isamp)
          vars(1, 1, :,ifield+6) = wnrav(:,isamp)
          ifield = ifield + 6
        end if
      end if
      if (ltenddec) then
        vars(1, 1, :,ifield+1) = wav(:,isamp)
        vars(1, 1, :,ifield+2) = uwav(:,isamp)
        vars(1, 1, :,ifield+3) = vsav(:,isamp)
        ifield = ifield+3
      end if
    nvar = ifield
    call writestat_nc(ncid,nvar,ncname(:,:,isamp),vars(:,:,1:kmax,:),nrec,1,1,kmax)
    enddo

    deallocate(vars)

  end subroutine writenetcdf_proc


!> Cleans up after the run
  subroutine exitsamptend
    use modstat_nc, only: exitstat_nc,lnetcdf
  implicit none

    if (.not. lsamptend) return
    if(isamptot == 0) return
    if(.not.(lnetcdf)) return
	call exitstat_nc(ncid)
    if (lsamptendu) deallocate (uptm, upmn, upav, ust)
    if (lsamptendv) deallocate (vptm, vpmn, vpav, vst)
    if (lsamptendw) deallocate (wptm, wpmn, wpav, wst)
    if (lsamptendthl) deallocate (thlptm, thlpmn, thlpav, thlst)
    if (lsamptendqt) deallocate (qtptm, qtpmn, qtpav, qtst)
    if (lsamptendqr) deallocate (qrptm, qrpmn, qrpav, qrst)
    if (lsamptendnr) deallocate (nrptm, nrpmn, nrpav, nrst)
    if (ltenddec) deallocate (wav, uwav, vsav)
    if (ltenddec .and. lsamptendthl) deallocate(thlav, thlwav, thlsav, uthlwav, vthlsav, wthlav)
    if (ltenddec .and. lsamptendqt) deallocate(qtav, qtwav, qtsav, uqtwav, vqtsav, wqtav)
    if (ltenddec .and. lsamptendqr) deallocate(qrav, qrwav, qrsav, uqrwav, vqrsav, wqrav)
    if (ltenddec .and. lsamptendnr) deallocate(nrav, nrwav, nrsav, unrwav, vnrsav, wnrav)

    deallocate (tendmask)
    deallocate (nrsamptot,nrsamp,nrsamplast,nrsampnew)

  end subroutine exitsamptend

end module
