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
  public :: initsamptend, samptend, exitsamptend, leibniztend
  save
!NetCDF variables
  integer, parameter :: nvar = 66
  character(80),allocatable,dimension(:,:,:) :: ncname
  character(80),dimension(1,4) :: tncname
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer,public,parameter :: tend_tot=1,tend_start=1,tend_adv=2,tend_subg=3,tend_force=4,tend_rad=5,&
                              tend_ls=6,tend_micro=7, tend_topbound=8,tend_pois=9,tend_addon=10, tend_coriolis=11, tend_totlb=12
  integer,parameter :: nrfields = 12
  character(20),dimension(10) :: samplname,longsamplname
  integer :: nsamples,isamp,isamptot
  logical :: ldosamptendwrite = .false. !< write tendencies after leibniz terms have been detemined
  logical :: ldosamptendleib = .false. !< determine leibniz terms
  real :: lastrk3coef

  real, allocatable :: uptm(:,:,:),vptm(:,:,:),wptm(:,:,:),thlptm(:,:,:),qtptm(:,:,:),qrptm(:,:,:),nrptm(:,:,:)
  real, allocatable :: upav(:,:,:),vpav(:,:,:),wpav(:,:,:),thlpav(:,:,:),qtpav(:,:,:),qrpav(:,:,:),nrpav(:,:,:)
  real, allocatable :: upmn(:,:,:),vpmn(:,:,:),wpmn(:,:,:),thlpmn(:,:,:),qtpmn(:,:,:),qrpmn(:,:,:),nrpmn(:,:,:)
  real, allocatable :: ust(:,:),vst(:,:),wst(:,:),thlst(:,:),qtst(:,:),qrst(:,:),nrst(:,:)
  logical, allocatable :: tendmask(:,:,:,:)
  integer, allocatable :: nrsamptot(:,:),nrsamp(:,:),nrsamplast(:,:),nrsampnew(:,:)
  character(80) :: fname = 'samptend.xxx.nc'
  integer :: ncid,nrec = 0

contains
!> Initialization routine, reads namelists and inits variables
subroutine initsamptend
    use modmpi,   only : myid
    use modglobal,only : cexpnr,dtmax,kmax,k1,ladaptive,&
                         btime,kmax,tres,cexpnr,j1,jh,i1,ih,kmax
    use modstat_nc, only : open_nc,define_nc,redefine_nc,ncinfo,nctiminfo,writestat_dims_nc,lnetcdf
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav

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

    idtav = dtav/tres
    itimeav = timeav/tres
    tnext      = idtav   +btime
    tnextwrite = itimeav +btime


    if (abs(timeav/dtav-nint(timeav/dtav))>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'dtav should be a integer multiple of dtmax'
    end if

    allocate (uptm(k1,nrfields,isamptot),vptm(k1,nrfields,isamptot),wptm(k1,nrfields,isamptot),thlptm(k1,nrfields,isamptot),&
    qtptm(k1,nrfields,isamptot),qrptm(k1,nrfields,isamptot),nrptm(k1,nrfields,isamptot))
    allocate (upmn(k1,nrfields,isamptot),vpmn(k1,nrfields,isamptot),wpmn(k1,nrfields,isamptot),thlpmn(k1,nrfields,isamptot),&
    qtpmn(k1,nrfields,isamptot),qrpmn(k1,nrfields,isamptot),nrpmn(k1,nrfields,isamptot))
    allocate (upav(k1,nrfields,isamptot),vpav(k1,nrfields,isamptot),wpav(k1,nrfields,isamptot),thlpav(k1,nrfields,isamptot),&
    qtpav(k1,nrfields,isamptot),qrpav(k1,nrfields,isamptot),nrpav(k1,nrfields,isamptot))
    allocate (tendmask(2-ih:i1+ih,2-jh:j1+jh,k1,isamptot))
    allocate (nrsamptot(k1,isamptot),nrsamp(k1,isamptot),nrsamplast(k1,isamptot),nrsampnew(k1,isamptot))
    allocate (ust(k1,isamptot),vst(k1,isamptot),wst(k1,isamptot),thlst(k1,isamptot),qtst(k1,isamptot),&
    qrst(k1,isamptot),nrst(k1,isamptot))
    uptm = 0.
    vptm = 0.
    wptm = 0.
    thlptm = 0.
    qtptm = 0.
    qrptm = 0.
    nrptm = 0.
    upmn = 0.
    vpmn = 0.
    wpmn = 0.
    thlpmn = 0.
    qtpmn = 0.
    qrpmn = 0.
    nrpmn = 0.
    upav = 0.
    vpav = 0.
    wpav = 0.
    thlpav = 0.
    qtpav = 0.
    qrpav = 0.
    nrpav = 0.
    ust = 0.
    vst = 0.
    wst = 0.
    thlst = 0.
    qtst = 0.
    qrst = 0.
    nrst = 0.

    tendmask=.false.
    nrsamp=0
    nrsamptot=0
    nrsamplast=0
    nrsampnew=0

    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = itimeav/idtav
      if (myid==0) then
        allocate(ncname(nvar,4,isamptot))
        call nctiminfo(tncname(1,:))
        fname(10:12) = cexpnr
        call open_nc(fname,ncid,nrec,n3=kmax)
        call define_nc( ncid,1,tncname)
        call writestat_dims_nc(ncid)
          do isamp=1,isamptot
          call ncinfo(ncname( 1,:,isamp),'utendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U advective tendency','m/s^2','tt')
          call ncinfo(ncname( 2,:,isamp),'utenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U diffusive tendency','m/s^2','tt')
          call ncinfo(ncname( 3,:,isamp),'utendfor'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U tendency due to other forces','m/s^2','tt')
          call ncinfo(ncname( 4,:,isamp),'utendcor','U coriolis tendency','m/s^2','tt')
          call ncinfo(ncname( 5,:,isamp),'utendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U large scale tendency','m/s^2','tt')
          call ncinfo(ncname( 6,:,isamp),'utendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U top boundary tendency','m/s^2','tt')
          call ncinfo(ncname( 7,:,isamp),'utendpois'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U pressure gradient tendency','m/s^2','tt')
          call ncinfo(ncname( 8,:,isamp),'utendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U in addons tendency','m/s^2','tt')
          call ncinfo(ncname( 9,:,isamp),'utendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U total tendency','m/s^2','tt')
          call ncinfo(ncname(10,:,isamp),'utendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'U total tendency with leibniz terms','m/s^2','tt')
          call ncinfo(ncname(11,:,isamp),'vtendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V advective tendency','m/s^2','tt')
          call ncinfo(ncname(12,:,isamp),'vtenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V diffusive tendency','m/s^2','tt')
          call ncinfo(ncname(13,:,isamp),'vtendfor'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V tendency due to other forces','m/s^2','tt')
          call ncinfo(ncname(14,:,isamp),'vtendcor'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V coriolis tendency','m/s^2','tt')
          call ncinfo(ncname(15,:,isamp),'vtendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V large scale tendency','m/s^2','tt')
          call ncinfo(ncname(16,:,isamp),'vtendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V top boundary tendency','m/s^2','tt')
          call ncinfo(ncname(17,:,isamp),'vtendpois'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V pressure gradient tendency','m/s^2','tt')
          call ncinfo(ncname(18,:,isamp),'vtendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V in addons tendency','m/s^2','tt')
          call ncinfo(ncname(19,:,isamp),'vtendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V total tendency','m/s^2','tt')
          call ncinfo(ncname(20,:,isamp),'vtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'V total tendency with leibniz terms','m/s^2','tt')
          call ncinfo(ncname(21,:,isamp),'wtendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W advective tendency','m/s^2','mt')
          call ncinfo(ncname(22,:,isamp),'wtenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W diffusive tendency','m/s^2','mt')
          call ncinfo(ncname(23,:,isamp),'wtendfor'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W tendency due to other forces','m/s^2','mt')
          call ncinfo(ncname(24,:,isamp),'wtendcor'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W coriolis tendency','m/s^2','mt')
          call ncinfo(ncname(25,:,isamp),'wtendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W large scale tendency','m/s^2','mt')
          call ncinfo(ncname(26,:,isamp),'wtendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W top boundary tendency','m/s^2','mt')
          call ncinfo(ncname(27,:,isamp),'wtendpois'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W pressure gradient tendency','m/s^2','mt')
          call ncinfo(ncname(28,:,isamp),'wtendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W in addons tendency','m/s^2','mt')
          call ncinfo(ncname(29,:,isamp),'wtendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W total tendency','m/s^2','mt')
          call ncinfo(ncname(30,:,isamp),'wtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'W total tendency with leibniz terms','m/s^2','tt')
          call ncinfo(ncname(31,:,isamp),'thltendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l advective tendency','K/s','tt')
          call ncinfo(ncname(32,:,isamp),'thltenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l diffusive tendency','K/s','tt')
          call ncinfo(ncname(33,:,isamp),'thltendrad'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l radiative tendency','K/s','tt')
          call ncinfo(ncname(34,:,isamp),'thltendmicro'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l microphysical tendency','K/s','tt')
          call ncinfo(ncname(35,:,isamp),'thltendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l large scale tendency','K/s','tt')
          call ncinfo(ncname(36,:,isamp),'thltendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l  top boundary tendency','K/s','tt')
          call ncinfo(ncname(37,:,isamp),'thltendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l in addons tendency','K/s','tt')
          call ncinfo(ncname(38,:,isamp),'thltendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l total tendency','K/s','tt')
          call ncinfo(ncname(39,:,isamp),'thltendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'theta_l total tendency with leibniz terms','K/s','tt')
          call ncinfo(ncname(40,:,isamp),'qttendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content advective tendency','kg/kg/s','tt')
          call ncinfo(ncname(41,:,isamp),'qttenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content diffusive tendency','kg/kg/s','tt')
          call ncinfo(ncname(42,:,isamp),'qttendrad'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content radiative tendency','kg/kg/s','tt')
          call ncinfo(ncname(43,:,isamp),'qttendmicro'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content microphysical tendency','kg/kg/s','tt')
          call ncinfo(ncname(44,:,isamp),'qttendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content large scale tendency','kg/kg/s','tt')
          call ncinfo(ncname(45,:,isamp),'qttendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content  top boundary tendency','kg/kg/s','tt')
          call ncinfo(ncname(46,:,isamp),'qttendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content in addons tendency','kg/kg/s','tt')
          call ncinfo(ncname(47,:,isamp),'qttendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content total tendency','kg/kg/s','tt')
          call ncinfo(ncname(48,:,isamp),'qttendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content total tendency with leibniz terms','kg/kg/s','tt')
          call ncinfo(ncname(49,:,isamp),'qrtendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content advective tendency','kg/kg/s','tt')
          call ncinfo(ncname(50,:,isamp),'qrtenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content diffusive tendency','kg/kg/s','tt')
          call ncinfo(ncname(51,:,isamp),'qrtendrad'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content radiative tendency','kg/kg/s','tt')
          call ncinfo(ncname(52,:,isamp),'qrtendmicro'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content microphysical tendency','kg/kg/s','tt')
          call ncinfo(ncname(53,:,isamp),'qrtendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content large scale tendency','kg/kg/s','tt')
          call ncinfo(ncname(54,:,isamp),'qrtendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content  top boundary tendency','kg/kg/s','tt')
          call ncinfo(ncname(55,:,isamp),'qrtendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content in addons tendency','kg/kg/s','tt')
          call ncinfo(ncname(56,:,isamp),'qrtendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content total tendency','kg/kg/s','tt')
          call ncinfo(ncname(57,:,isamp),'qrtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'total water content total tendency with leibniz terms','kg/kg/s','tt')
          call ncinfo(ncname(58,:,isamp),'nrtendadv'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC advective tendency','/kg/s','tt')
          call ncinfo(ncname(59,:,isamp),'nrtenddif'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC diffusive tendency','/kg/s','tt')
          call ncinfo(ncname(60,:,isamp),'nrtendrad'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC radiative tendency','/kg/s','tt')
          call ncinfo(ncname(61,:,isamp),'nrtendmicro'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC microphysical tendency','/kg/s','tt')
          call ncinfo(ncname(62,:,isamp),'nrtendls'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC large scale tendency','/kg/s','tt')
          call ncinfo(ncname(63,:,isamp),'nrtendtop'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC top boundary tendency','/kg/s','tt')
          call ncinfo(ncname(64,:,isamp),'nrtendaddon'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC addons tendency','/kg/s','tt')
          call ncinfo(ncname(65,:,isamp),'nrtendtot'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC total tendency','/kg/s','tt')
          call ncinfo(ncname(66,:,isamp),'nrtendleib'//samplname(isamp),&
          trim(longsamplname(isamp))//' '//'RDNC total tendency with leibniz terms','/kg/s','tt')
          call define_nc(ncid,nvar,ncname(:,:,isamp))
        enddo
      end if
    end if

  end subroutine initsamptend

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
    real, allocatable, dimension(:,:,:) :: w0f,wpf
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
      nrsamplast=0
      tendmask=.false.
      uptm = 0.
      vptm = 0.
      wptm = 0.
      thlptm = 0.
      qtptm = 0.
      qrptm = 0.
      nrptm = 0.
      ust = 0.
      vst = 0.
      wst = 0.
      thlst = 0.
      qtst = 0.
      qrst = 0.
      nrst = 0.

      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
                w0f(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thvav(k1))

      do k=1,k1
        thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do

      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/ijtot

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

      deallocate(thv0,w0f)
      deallocate(thvav)

      do isamp=1,isamptot
      do k=1,kmax
        ust(k,isamp) = sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        vst(k,isamp) = sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        wst(k,isamp) = sum(w0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        thlst(k,isamp) = sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        qtst(k,isamp) = sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
        if(nsv>1) then
        qrst(k,isamp) = sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
        nrst(k,isamp) = sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
        endif
      end do
      end do

      ldosamptendleib=.true.
      lastrk3coef = rdt / (4. - dble(rk3step))

    ENDIF
    ENDIF

    allocate(wpf(2-ih:i1+ih,2-jh:j1+jh,k1))

    do k=1,kmax
      wpf (2:i1,2:j1,k) = 0.5*(wp (2:i1,2:j1,k) + wp  (2:i1,2:j1,k+1))
    end do

    do isamp=1,isamptot
    do k=1,kmax
      uptm(k,tendterm,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-uptm (k,tend_tot,isamp)
      vptm(k,tendterm,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-vptm (k,tend_tot,isamp)
      wptm(k,tendterm,isamp) = sum(wpf (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-wptm (k,tend_tot,isamp)
      thlptm(k,tendterm,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-thlptm (k,tend_tot,isamp)
      qtptm(k,tendterm,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))-qtptm (k,tend_tot,isamp)
      if(nsv>1) then
      qrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))-qrptm (k,tend_tot,isamp)
      nrptm(k,tendterm,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))-nrptm (k,tend_tot,isamp)
      endif
      uptm(k,tend_tot,isamp) = sum(up (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      vptm(k,tend_tot,isamp) = sum(vp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      wptm(k,tend_tot,isamp) = sum(wpf (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      thlptm(k,tend_tot,isamp) = sum(thlp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      qtptm(k,tend_tot,isamp) = sum(qtp (2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))
      if(nsv>1) then
      qrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))
      nrptm(k,tend_tot,isamp) = sum(svp (2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))
      endif
      upav(k,tendterm,isamp) = upav(k,tendterm,isamp)+uptm(k,tendterm,isamp)
      vpav(k,tendterm,isamp) = vpav(k,tendterm,isamp)+vptm(k,tendterm,isamp)
      wpav(k,tendterm,isamp) = wpav(k,tendterm,isamp)+wptm(k,tendterm,isamp)
      thlpav(k,tendterm,isamp) = thlpav(k,tendterm,isamp)+thlptm(k,tendterm,isamp)
      qtpav(k,tendterm,isamp) = qtpav(k,tendterm,isamp)+qtptm(k,tendterm,isamp)
      qrpav(k,tendterm,isamp) = qrpav(k,tendterm,isamp)+qrptm(k,tendterm,isamp)
      nrpav(k,tendterm,isamp) = nrpav(k,tendterm,isamp)+nrptm(k,tendterm,isamp)
    end do
    end do

    deallocate(wpf)

    IF (present(lastterm)) THEN
    IF (lastterm) THEN
      do isamp=1,isamptot
      do k=1,kmax
        upav(k,tend_tot,isamp) = upav(k,tend_tot,isamp)+uptm(k,tend_tot,isamp)
        vpav(k,tend_tot,isamp) = vpav(k,tend_tot,isamp)+vptm(k,tend_tot,isamp)
        wpav(k,tend_tot,isamp) = wpav(k,tend_tot,isamp)+wptm(k,tend_tot,isamp)
        thlpav(k,tend_tot,isamp) = thlpav(k,tend_tot,isamp)+thlptm(k,tend_tot,isamp)
        qtpav(k,tend_tot,isamp) = qtpav(k,tend_tot,isamp)+qtptm(k,tend_tot,isamp)
        qrpav(k,tend_tot,isamp) = qrpav(k,tend_tot,isamp)+qrptm(k,tend_tot,isamp)
        nrpav(k,tend_tot,isamp) = nrpav(k,tend_tot,isamp)+nrptm(k,tend_tot,isamp)
      enddo
      enddo
      tnext = tnext+idtav

      if (timee>=tnextwrite) then
        tnextwrite = tnextwrite+itimeav
        ldosamptendwrite = .true.
      end if
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
    END IF
    END IF

  end subroutine samptend

  subroutine leibniztend
    use modmpi,    only : slabsum
    use modglobal, only : i1,j1,kmax,k1,ih,jh,&
                          cp,rv,rlv,rd,&
                          ijtot,nsv
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

      allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1),&
                w0f(2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(thvav(k1))

      do k=1,k1
        thv0(2:i1,2:j1,k) = (thl0(2:i1,2:j1,k)+rlv*ql0(2:i1,2:j1,k)/(cp*exnf(k))) &
                    *(1+(rv/rd-1)*qt0(2:i1,2:j1,k)-rv/rd*ql0(2:i1,2:j1,k))
      enddo
      do k=1,kmax
        w0f (2:i1,2:j1,k) = 0.5*(w0 (2:i1,2:j1,k) + w0  (2:i1,2:j1,k+1))
      end do

      thvav = 0.0
      call slabsum(thvav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      thvav = thvav/ijtot

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

      deallocate(thv0)
      deallocate(thvav)

      do isamp=1,isamptot
      do k=1,kmax
        if((nrsampnew(k,isamp)>0).and.(nrsamplast(k,isamp)>0)) then! only do if sampling can be performed at both points in time
          upav(k,tend_totlb,isamp) = upav(k,tend_totlb,isamp)+(ust(k,isamp)-&
          sum(u0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          vpav(k,tend_totlb,isamp) = vpav(k,tend_totlb,isamp)+(vst(k,isamp)-&
          sum(v0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          wpav(k,tend_totlb,isamp) = wpav(k,tend_totlb,isamp)+(wst(k,isamp)-&
          sum(w0f(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          thlpav(k,tend_totlb,isamp) = thlpav(k,tend_totlb,isamp)+(thlst(k,isamp)-&
          sum(thl0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          qtpav(k,tend_totlb,isamp) = qtpav(k,tend_totlb,isamp)+(qtst(k,isamp)-&
          sum(qt0(2:i1,2:j1,k),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          if(nsv>1) then
            qrpav(k,tend_totlb,isamp) = qrpav(k,tend_totlb,isamp)+(qrst(k,isamp)-&
            sum(sv0(2:i1,2:j1,k,iqr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
            nrpav(k,tend_totlb,isamp) = nrpav(k,tend_totlb,isamp)+(nrst(k,isamp)-&
            sum(sv0(2:i1,2:j1,k,inr),tendmask(2:i1,2:j1,k,isamp))*nrsamplast(k,isamp)/nrsampnew(k,isamp))/lastrk3coef
          endif
        endif
      enddo
      enddo

      deallocate(w0f)

      if(ldosamptendwrite) then
        ldosamptendwrite=.false.
        call writesamptend
        upav = 0.
        vpav = 0.
        wpav = 0.
        thlpav = 0.
        qtpav = 0.
        qrpav = 0.
        nrpav = 0.
        nrsamp = 0
      endif

  end subroutine leibniztend

!> Write the statistics to file
  subroutine writesamptend
    use modglobal, only : kmax,k1,rtimee
    use modmpi,    only : mpi_integer,myid,comm3d,mpierr,mpi_sum &
                        , D_MPI_ALLREDUCE
    use modstat_nc, only: lnetcdf,writestat_nc
    implicit none
    integer :: field,k
    real, allocatable :: vars(:,:)

    if (.not. lsamptend) return

    allocate(vars(1:k1,nvar))
    upmn = 0.
    vpmn = 0.
    wpmn = 0.
    thlpmn = 0.
    qtpmn = 0.
    qrpmn = 0.
    nrpmn = 0.
    nrsamptot=0

    call D_MPI_ALLREDUCE(nrsamp   ,nrsamptot ,k1*isamptot         ,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(upav     ,upmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(vpav     ,vpmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(wpav     ,wpmn      ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(thlpav   ,thlpmn    ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(qtpav    ,qtpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(qrpav    ,qrpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)
    call D_MPI_ALLREDUCE(nrpav    ,nrpmn     ,k1*nrfields*isamptot,MPI_SUM,comm3d,mpierr)

    do field=1,nrfields
    do isamp=1,isamptot
    do k=1,k1
      if (nrsamptot(k,isamp)>0) then
        upmn  (k,field,isamp) = upmn (k,field,isamp)/nrsamptot(k,isamp)
        vpmn  (k,field,isamp) = vpmn (k,field,isamp)/nrsamptot(k,isamp)
        wpmn  (k,field,isamp) = wpmn (k,field,isamp)/nrsamptot(k,isamp)
        thlpmn(k,field,isamp) = thlpmn (k,field,isamp)/nrsamptot(k,isamp)
        qtpmn (k,field,isamp) = qtpmn (k,field,isamp)/nrsamptot(k,isamp)
        qrpmn (k,field,isamp) = qrpmn (k,field,isamp)/nrsamptot(k,isamp)
        nrpmn (k,field,isamp) = nrpmn (k,field,isamp)/nrsamptot(k,isamp)
      endif
    enddo
    enddo
    enddo

    if(myid == 0) then
      if (lnetcdf) then
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
        do isamp=1,isamptot
          vars=0.
          vars(:, 1) = upmn(:,tend_adv,isamp)
          vars(:, 2) = upmn(:,tend_subg,isamp)
          vars(:, 3) = upmn(:,tend_force,isamp)
          vars(:, 4) = upmn(:,tend_coriolis,isamp)
          vars(:, 5) = upmn(:,tend_ls,isamp)
          vars(:, 6) = upmn(:,tend_topbound,isamp)
          vars(:, 7) = upmn(:,tend_pois,isamp)
          vars(:, 8) = upmn(:,tend_addon,isamp)
          vars(:, 9) = upmn(:,tend_tot,isamp)
          vars(:,10) = upmn(:,tend_totlb,isamp)
          vars(:,11) = vpmn(:,tend_adv,isamp)
          vars(:,12) = vpmn(:,tend_subg,isamp)
          vars(:,13) = vpmn(:,tend_force,isamp)
          vars(:,14) = vpmn(:,tend_coriolis,isamp)
          vars(:,15) = vpmn(:,tend_ls,isamp)
          vars(:,16) = vpmn(:,tend_topbound,isamp)
          vars(:,17) = vpmn(:,tend_pois,isamp)
          vars(:,18) = vpmn(:,tend_addon,isamp)
          vars(:,19) = vpmn(:,tend_tot,isamp)
          vars(:,20) = vpmn(:,tend_totlb,isamp)
          vars(:,21) = wpmn(:,tend_adv,isamp)
          vars(:,22) = wpmn(:,tend_subg,isamp)
          vars(:,23) = wpmn(:,tend_force,isamp)
          vars(:,24) = wpmn(:,tend_coriolis,isamp)
          vars(:,25) = wpmn(:,tend_ls,isamp)
          vars(:,26) = wpmn(:,tend_topbound,isamp)
          vars(:,27) = wpmn(:,tend_pois,isamp)
          vars(:,28) = wpmn(:,tend_addon,isamp)
          vars(:,29) = wpmn(:,tend_tot,isamp)
          vars(:,30) = wpmn(:,tend_totlb,isamp)
          vars(:,31) = thlpmn(:,tend_adv,isamp)
          vars(:,32) = thlpmn(:,tend_subg,isamp)
          vars(:,33) = thlpmn(:,tend_rad,isamp)
          vars(:,34) = thlpmn(:,tend_micro,isamp)
          vars(:,35) = thlpmn(:,tend_ls,isamp)
          vars(:,36) = thlpmn(:,tend_topbound,isamp)
          vars(:,37) = thlpmn(:,tend_addon,isamp)
          vars(:,38) = thlpmn(:,tend_tot,isamp)
          vars(:,39) = thlpmn(:,tend_totlb,isamp)
          vars(:,40) = qtpmn(:,tend_adv,isamp)
          vars(:,41) = qtpmn(:,tend_subg,isamp)
          vars(:,42) = qtpmn(:,tend_rad,isamp)
          vars(:,43) = qtpmn(:,tend_micro,isamp)
          vars(:,44) = qtpmn(:,tend_ls,isamp)
          vars(:,45) = qtpmn(:,tend_topbound,isamp)
          vars(:,46) = qtpmn(:,tend_addon,isamp)
          vars(:,47) = qrpmn(:,tend_tot,isamp)
          vars(:,48) = qrpmn(:,tend_totlb,isamp)
          vars(:,49) = qrpmn(:,tend_adv,isamp)
          vars(:,50) = qrpmn(:,tend_subg,isamp)
          vars(:,51) = qrpmn(:,tend_rad,isamp)
          vars(:,52) = qrpmn(:,tend_micro,isamp)
          vars(:,53) = qrpmn(:,tend_ls,isamp)
          vars(:,54) = qrpmn(:,tend_topbound,isamp)
          vars(:,55) = qrpmn(:,tend_addon,isamp)
          vars(:,56) = qrpmn(:,tend_tot,isamp)
          vars(:,57) = qrpmn(:,tend_totlb,isamp)
          vars(:,58) = nrpmn(:,tend_adv,isamp)
          vars(:,59) = nrpmn(:,tend_subg,isamp)
          vars(:,60) = nrpmn(:,tend_rad,isamp)
          vars(:,61) = nrpmn(:,tend_micro,isamp)
          vars(:,62) = nrpmn(:,tend_ls,isamp)
          vars(:,63) = nrpmn(:,tend_topbound,isamp)
          vars(:,64) = nrpmn(:,tend_addon,isamp)
          vars(:,65) = nrpmn(:,tend_tot,isamp)
          vars(:,66) = nrpmn(:,tend_totlb,isamp)
        call writestat_nc(ncid,nvar,ncname(:,:,isamp),vars(1:kmax,:),nrec,kmax)
        enddo
      end if
    end if
    deallocate(vars)

  end subroutine writesamptend

!> Cleans up after the run
  subroutine exitsamptend
    use modstat_nc, only: exitstat_nc,lnetcdf
    use modmpi, only: myid
  implicit none

    if (.not. lsamptend) return
    if(isamptot == 0) return
    if(.not.(lnetcdf)) return

    if (myid==0) then
       call exitstat_nc(ncid)
    end if
    
    deallocate (uptm,vptm,wptm,thlptm,qtptm,qrptm,nrptm)
    deallocate (upmn,vpmn,wpmn,thlpmn,qtpmn,qrpmn,nrpmn)
    deallocate (upav,vpav,wpav,thlpav,qtpav,qrpav,nrpav)
    deallocate (ust,vst,wst,thlst,qtst,qrst,nrst)
    deallocate (tendmask)
    deallocate (nrsamptot,nrsamp,nrsamplast,nrsampnew)

  end subroutine exitsamptend

end module
