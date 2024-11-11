!> \file modgenstat.f90
!!  Genstat calculates slab averages of several variables

!>
!!  Genstat calculates slab averages of several variables
!>
!!  Written to fields.expnr, moments.expnr and flux1.expnr and flux2.expnr
!! If netcdf is true, this module leads the profiles.expnr.nc output
!!  \author Hans Cuijpers, K.N.M.I.
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
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



module modgenstat

    !-----------------------------------------------------------------|
    !                                                                 |
    !*** *stattend*  calculates generic slabaveraged statistics       |
    !                                                                 |
    !      Pier Siebesma   K.N.M.I.     12/05/1995                    |
    !      Hans Cuijpers   I.M.A.U.     21/06/1995                    |
    !                                                                 |
    !     purpose.                                                    |
    !     --------                                                    |
    !                                                                 |
    !     stattend.f calculates:                                      |
    !                                                                 |
    !     * time averaged fieldsof theta_l, theta, theta_v, qt, qv   |
    !        ql, u and v                                              |
    !     * time averaged tendencies  of theta_l, theta, qt, qv, ql,  |
    !        u and v                                                  |
    !     * time averaged turbulent fluxes of theta_l, theta_v,       |
    !       theta, qt, qv, ql, u and v                                |
    !     * variances of qt, w, u, theta_l and theta_v                |
    !     * skewness of qt and w                                      |
    !*** *genstat*  calculates timeseries of several variables       |
    !                                                                 |
    !____________________SETTINGS_AND_SWITCHES________________________|
    !                     IN &NAMTIMESTAT                             |
    !                                                                 |
    !    dtav           SAMPLING INTERVAL                             |
    !                                                                 |
    !    timeav         INTERVAL OF WRITING                           |
    !                                                                 |
    !    lstat      SWITCH TO ENABLE TIMESERIES                       |
    !-----------------------------------------------------------------|
  use modprecision
  use modtimer

  implicit none
  ! private
  PUBLIC :: initgenstat, genstat, exitgenstat
  save

!NetCDF variables
  integer :: nvar = 48
  integer :: ncid,nrec = 0
  character(80) :: fname = 'profiles.xxx.nc'
  character(80),allocatable, dimension(:,:) :: ncname
  character(80),dimension(1,4) :: tncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  logical :: lstat= .false. ! switch for conditional sampling cloud (on/off)
  integer :: nsamples
!     ----  total fields  ---

  real, allocatable  :: umn   (:)       ,vmn   (:),  wmn  (:)
  real, allocatable  :: thlmn (:)       ,thvmn (:)
  real, allocatable  :: qtmn  (:)       ,qlmn  (:),  qlhmn(:),cfracmn(:),hurmn(:),tamn(:)
  real, allocatable  :: clwmn(:), climn(:), plwmn(:), plimn(:)

! real, allocatable  ::     --- fluxes (resolved, subgrid and total) ---
  real, allocatable  :: wthlsmn (:),wthlrmn (:),wthltmn(:)
  real, allocatable  :: wthvsmn (:),wthvrmn (:),wthvtmn(:)
  real, allocatable  :: wqlsmn (:),wqlrmn (:),wqltmn(:)
  real, allocatable  :: wqtsmn (:),wqtrmn (:),wqttmn(:)
  real, allocatable  :: wsvsmn (:,:),wsvrmn(:,:),wsvtmn(:,:)
  real, allocatable  :: uwtmn  (:),vwtmn  (:) !total    uw, vw
  real, allocatable  :: uwsmn  (:),vwsmn  (:) !resolved uw, vw
  real, allocatable  :: uwrmn  (:),vwrmn  (:) !subgrid  uw, vw
! real, allocatable  ::     --- various moments ---

  real, allocatable  :: w2mn       (:), skewmn (:)
  real, allocatable  :: w2submn    (:)
  real, allocatable  :: u2mn       (:), v2mn  (:),     qt2mn(:)
  real, allocatable  :: thl2mn     (:), thv2mn(:),     th2mn(:),     ql2mn(:)
  real, allocatable  :: svmmn(:,:),svptmn(:,:),svplsmn(:,:),svpmn(:,:)
  real, allocatable  :: sv2mn(:,:)

  real(field_r), allocatable :: umav (:)     ! slab averaged u_0    at full level
  real(field_r), allocatable :: vmav (:)     ! slab averaged v_0    at full level
  real(field_r), allocatable :: wmav (:)     ! slab averaged w_0    at full level
  real(field_r), allocatable :: thlmav (:)     ! slab averaged thl_0    at full level
  real(field_r), allocatable :: thmav (:)     ! slab averaged th_0    at full level
  real(field_r), allocatable :: qtmav (:)     ! slab averaged qt_0    at full level
  real(field_r), allocatable :: qlmav (:)     ! slab averaged ql_0    at full level
  real, allocatable :: cfracav (:)     ! slab averaged ql_0    at full level
  real, allocatable :: hurav (:)
  real, allocatable :: clwav(:), cliav(:), plwav(:), pliav(:)
  real(field_r), allocatable :: svmav (:,:)     ! slab averaged ql_0    at full level
  real(field_r), allocatable :: taav (:)
  real, allocatable :: svpav(:,:)                  !  slab average total tendency of sv(n)
  real, allocatable :: svptav(:,:)                 !  slab average tendency of sv(n) due to turb.

  real, allocatable :: uptav(:)                      !  slab averaged tendency for u
  real, allocatable :: vptav(:)                      !  slab averaged tendency for v

  real, allocatable :: thptav(:)     ! slab averaged turbulence tendency of theta
  real, allocatable :: qlptav(:)     ! slab averaged turbulence tendency of q_liq
  real, allocatable :: uwtot (:)     ! slab averaged tot w-u flux at half levels
  real, allocatable :: vwtot (:)     ! slab averaged tot w-v flux at half levels
  real, allocatable :: uwsub (:)     ! slab averaged sgs w-u flux at half levels
  real, allocatable :: vwsub (:)     ! slab averaged sgs w-v flux at half levels
  real, allocatable :: uwres (:)     ! slab averaged res w-u flux at half levels
  real, allocatable :: vwres (:)     ! slab averaged res w-v flux at half levels


  real, allocatable :: thvhav(:)
  real, allocatable :: th0av(:)
  real, allocatable :: wthlsub(:)     ! slab averaged sub w-theta_l flux at half levels
  real, allocatable :: wthlres(:)     ! slab averaged res w-theta_l flux at half levels
  real, allocatable :: wthltot(:)     ! slab averaged tot w-theta_l flux at half levels

  real, allocatable :: wqtsub(:)    ! slab averaged sub w-qtot    flux at half levels
  real, allocatable :: wqtres(:)    ! slab averaged res w-qtot    flux at half levels
  real, allocatable :: wqttot(:)    ! slab averaged tot w-qtot    flux at half levels

  real, allocatable :: wqlsub(:)    ! slab averaged sub w-ql      flux at half levels
  real, allocatable :: wqlres(:)    ! slab averaged tot w-ql      flux at half levels
  real, allocatable :: wqltot(:)    ! slab averaged tot w-ql      flux at half levels

  real, allocatable :: wthvtot(:)    !  slab averaged total wthv-flux
  real, allocatable :: wthvres(:)    !  slab averaged res   wthv-flux
  real, allocatable :: wthvsub(:)    !  slab averaged sub   wthv-flux

  real, allocatable :: wsvsub(:,:)! slab averaged sub w-sv(n)  flux
  real, allocatable :: wsvres(:,:)! slab averaged res w-sv(n)  flux
  real, allocatable :: wsvtot(:,:)! slab averaged tot w-sv(n)  flux

  real, allocatable :: cszmn(:)    ! Smagorinsky constant
  real, allocatable :: cszav(:)    ! Smagorinsky constant

  real, allocatable :: qlmnlast(:)
  real, allocatable :: wthvtmnlast(:)

    ! Local fields

  real,allocatable, dimension(:):: qlptavl   ! slab averaged turbulence tendency of q_liq
  real,allocatable, dimension(:):: uwsubl
  real,allocatable, dimension(:):: vwsubl
  real,allocatable, dimension(:):: uwresl
  real,allocatable, dimension(:):: vwresl
  real,allocatable, dimension(:):: qlhav
  real,allocatable, dimension(:):: u2av    , &
            v2av    , &
            w2av    , &
            w3av    , &
            w2subav , &
            qt2av   , &
            thl2av  , &
            thv2av  , &
            th2av   , &
            ql2av
  real, allocatable, dimension(:,:) :: sv2av
  real(field_r),allocatable, dimension(:,:,:)::  thv0
  real(field_r),allocatable, dimension(:)::   thvmav
  real(field_r),allocatable, dimension(:,:,:):: sv0h

contains

  subroutine initgenstat
    use modmpi, only : myid,mpierr, comm3d, D_MPI_BCAST
    use modglobal, only : i1, ih, j1, jh, kmax, k1, nsv, ifnamopt, fname_options, ifoutput, &
                          cexpnr, dtav_glob, timeav_glob, dt_lim, btime, tres, &
                          lwarmstart, checknamelisterror
    use modstat_nc, only : lnetcdf, open_nc, ncinfo, define_nc, nctiminfo, writestat_dims_nc
    use modsurfdata, only : isurf, ksoilmax
    use modlsm, only : kmax_soil
    use modtracers, only : tracer_prop

    implicit none

    integer n, ierr
    character(40) :: name

    namelist/NAMGENSTAT/ &
    dtav,timeav,lstat

    call timer_tic('modgenstat/initgenstat', 0)

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMGENSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMGENSTAT')
      write(6 ,NAMGENSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(timeav     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lstat      ,1,0,comm3d,mpierr)
    idtav = int(dtav/tres,kind=longint)
    itimeav = int(timeav/tres,kind=longint)

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = int(itimeav/idtav)
    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
       write (*,*) timeav, dtav, nsamples, timeav/dtav-nsamples
       write (*,*) itimeav, idtav, tres
      stop 'NAMGENSTAT: timeav must be a integer multiple of dtav'
    end if

    allocate(umn(k1),vmn(k1),wmn(k1))
    allocate(thlmn (k1)       ,thvmn (k1))
    allocate(qtmn  (k1)       ,qlmn  (k1),  qlhmn(k1),cfracmn(k1), hurmn(k1), tamn(k1))
    allocate(clwmn(k1), climn(k1), plwmn(k1), plimn(k1))
    allocate(wthlsmn (k1),wthlrmn (k1),wthltmn(k1))
    allocate(wthvsmn (k1),wthvrmn (k1),wthvtmn(k1))
    allocate(wqlsmn (k1),wqlrmn (k1),wqltmn(k1))
    allocate(wqtsmn (k1),wqtrmn (k1),wqttmn(k1))
    allocate(wsvsmn (k1,nsv),wsvrmn(k1,nsv),wsvtmn(k1,nsv))
    allocate(uwtmn  (k1),vwtmn  (k1),uwrmn  (k1),  vwrmn(k1),uwsmn(k1),vwsmn(k1))
    allocate(w2mn       (k1), skewmn (k1))
    allocate(w2submn    (k1))
    allocate(u2mn       (k1), v2mn  (k1),     qt2mn(k1))
    allocate(thl2mn     (k1), thv2mn(k1),     th2mn(k1),     ql2mn(k1))
    allocate(svmmn(k1,nsv),svptmn(k1,nsv),svplsmn(k1,nsv),svpmn(k1,nsv))
    allocate(sv2mn(k1,nsv))

    allocate(umav (k1))
    allocate(vmav (k1))
    allocate(wmav (k1))
    allocate(thlmav (k1))
    allocate(thmav (k1))
    allocate(qtmav (k1))
    allocate(qlmav (k1))
    allocate(cfracav(k1))
    allocate(hurav(k1))
    allocate(clwav(k1), cliav(k1), plwav(k1), pliav(k1))
    allocate(taav(k1))
    allocate(svmav (k1,nsv))
    allocate(uptav(k1))
    allocate(vptav(k1))

    allocate(thptav(k1))
    allocate(qlptav(k1))
    allocate(uwtot (k1))
    allocate(vwtot (k1))
    allocate(uwsub (k1))
    allocate(vwsub (k1))
    allocate(uwres (k1))
    allocate(vwres (k1))
    allocate(wthlsub(k1))
    allocate(wthlres(k1))
    allocate(wthltot(k1))
    allocate(wqtsub(k1))
    allocate(wqtres(k1))
    allocate(wqttot(k1))
    allocate(wqlsub(k1))
    allocate(wqlres(k1))
    allocate(wqltot(k1))
    allocate(wthvtot(k1))
    allocate(wthvres(k1))
    allocate(wthvsub(k1))
    allocate(wsvsub(k1,nsv))
    allocate(wsvres(k1,nsv))
    allocate(wsvtot(k1,nsv))
    allocate(thvhav(k1))
    allocate(th0av(k1))
    allocate(svptav(k1,nsv))
    allocate(svpav(k1,nsv))

    allocate(cszmn(k1), cszav(k1))

    allocate(qlmnlast(k1))
    allocate(wthvtmnlast(k1))

    allocate( qlhav(k1))
    allocate( u2av    (k1), &
              v2av    (k1), &
              w2av    (k1), &
              w3av    (k1), &
              w2subav (k1), &
              qt2av   (k1), &
              thl2av  (k1), &
              thv2av  (k1), &
              th2av   (k1), &
              ql2av   (k1), &
              sv2av   (k1,nsv))
    allocate(thv0(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(thvmav(k1))
    allocate(sv0h(2-ih:i1+ih,2-jh:j1+jh,k1))

    umn      = 0.
    vmn      = 0.
    wmn      = 0.
    thlmn    = 0.
    thvmn    = 0.
    qtmn     = 0.
    qlmn     = 0.
    qlhmn    = 0.
    cfracmn  = 0.
    hurmn    = 0.
    clwmn    = 0.
    climn    = 0.
    plwmn    = 0.
    plimn    = 0.
    tamn     = 0.

    wthlsmn =  0.
    wthlrmn =  0.
    wthltmn =  0.

    wthvsmn =  0.
    wthvrmn =  0.
    wthvtmn =  0.

    wqtsmn =  0.
    wqtrmn =  0.
    wqttmn =  0.

    wqlsmn =  0.
    wqlrmn =  0.
    wqltmn =  0.

    uwtmn  = 0.
    vwtmn  = 0.
    uwrmn  = 0.
    vwrmn  = 0.
    uwsmn  = 0.
    vwsmn  = 0.

    u2mn     = 0.
    v2mn     = 0.
    w2mn     = 0.
    w2submn  = 0.
    skewmn   = 0.
    qt2mn    = 0.
    thl2mn   = 0.
    thv2mn   = 0.
    th2mn    = 0.
    ql2mn    = 0.

    svmmn   = 0.
    svpmn   = 0.
    svpav   = 0.
    svptmn  = 0.
    svptav  = 0.

    sv2mn = 0.

    wsvsmn = 0.
    wsvrmn = 0.
    wsvtmn = 0.

    cszav = 0.
    cszmn = 0.

    qlmnlast = 0.
    wthvtmnlast = 0.

    if(myid==0)then
       if (.not. lwarmstart) then
          open (ifoutput,file='field.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='flux1.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='flux2.'//cexpnr,status='replace')
          close (ifoutput)
          open (ifoutput,file='moments.'//cexpnr,status='replace')
          close (ifoutput)
          do n=1,nsv
             name = 'svnnnfld.'//cexpnr
             write (name(3:5),'(i3.3)') n ! TODO apply tracer_props
             open (ifoutput,file=name,status='replace')
             close (ifoutput)
          end do
          do n=1,nsv
             name = 'svnnnflx.'//cexpnr
             write (name(3:5),'(i3.3)') n ! TODO apply tracer_props
             open (ifoutput,file=name,status='replace')
             close (ifoutput)
          end do
       end if

      if (lnetcdf) then
        fname(10:12) = cexpnr
        nvar = nvar + 7*nsv
        allocate(ncname(nvar,4))
        call nctiminfo(tncname(1,:))
        call ncinfo(ncname( 1,:),'rhof','Full level slab averaged density','kg/m^3','tt')
        call ncinfo(ncname( 2,:),'rhobf','Full level base-state density','kg/m^3','tt')
        call ncinfo(ncname( 3,:),'rhobh','Half level base-state density','kg/m^3','mt')
        call ncinfo(ncname( 4,:),'presh','Pressure at cell center','Pa','tt')
        call ncinfo(ncname( 5,:),'u','West-East velocity','m/s','tt')
        call ncinfo(ncname( 6,:),'v','South-North velocity','m/s','tt')
        call ncinfo(ncname( 7,:),'w','Vertical velocity','m/s','mt')
        call ncinfo(ncname( 8,:),'thl','Liquid water potential temperature','K','tt')
        call ncinfo(ncname( 9,:),'thv','Virtual potential temperature','K','tt')
        call ncinfo(ncname(10,:),'qt','Total water specific humidity','kg/kg','tt')
        call ncinfo(ncname(11,:),'ql','Liquid water specific humidity','kg/kg','tt')
        call ncinfo(ncname(12,:),'wthls','SFS-Theta_l flux','Km/s','mt')
        call ncinfo(ncname(13,:),'wthlr','Resolved Theta_l flux','Km/s','mt')
        call ncinfo(ncname(14,:),'wthlt','Total Theta_l flux','Km/s','mt')
        call ncinfo(ncname(15,:),'wthvs','SFS-buoyancy flux','Km/s','mt')
        call ncinfo(ncname(16,:),'wthvr','Resolved buoyancy flux','Km/s','mt')
        call ncinfo(ncname(17,:),'wthvt','Total buoyancy flux','Km/s','mt')
        call ncinfo(ncname(18,:),'wqts','SFS-moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(19,:),'wqtr','Resolved moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(20,:),'wqtt','Total moisture flux','kg/kg m/s','mt')
        call ncinfo(ncname(21,:),'wqls','SFS-liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(22,:),'wqlr','Resolved liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(23,:),'wqlt','Total liquid water flux','kg/kg m/s','mt')
        call ncinfo(ncname(24,:),'uws','SFS-momentum flux (uw)','m^2/s^2','mt')
        call ncinfo(ncname(25,:),'uwr','Resolved momentum flux (uw)','m^2/s^2','mt')
        call ncinfo(ncname(26,:),'uwt','Total momentum flux (uw)','m^2/s^2','mt')
        call ncinfo(ncname(27,:),'vws','SFS-momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(28,:),'vwr','Resolved momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(29,:),'vwt','Total momentum flux (vw)','m^2/s^2','mt')
        call ncinfo(ncname(30,:),'w2s','SFS-TKE','m^2/s^2','mt')
        call ncinfo(ncname(31,:),'w2r','Resolved vertical velocity variance','m^2/s^2','mt')
        !call ncinfo(ncname(31,:),'w2t','Total vertical velocity variance','m^2/s^2','mt')
        call ncinfo(ncname(32,:),'skew','vertical velocity skewness','-','mt')
        call ncinfo(ncname(33,:),'u2r','Resolved horizontal velocity variance (u)','m^2/s^2','tt')
        call ncinfo(ncname(34,:),'v2r','Resolved horizontal velocity variance (v)','m^2/s^2','tt')
        call ncinfo(ncname(35,:),'thl2r','Resolved theta_l variance','K^2','tt')
        call ncinfo(ncname(36,:),'thv2r','Resolved buoyancy variance','K^2','tt')
        call ncinfo(ncname(37,:),'th2r','Resolved theta variance','K^2','tt')
        call ncinfo(ncname(38,:),'qt2r','Resolved total water variance','(kg/kg)^2','tt')
        call ncinfo(ncname(39,:),'ql2r','Resolved liquid water variance','(kg/kg)^2','tt')
        call ncinfo(ncname(40,:),'cs','Smagorinsky constant','-','tt')
        call ncinfo(ncname(41,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname(42,:),'hur','Relative humidity','%','tt')
        call ncinfo(ncname(43,:),'hus','Specific humidity','kg/kg','tt')
        call ncinfo(ncname(44,:),'ta', 'Temperature','K','tt')
        call ncinfo(ncname(45,:),'clw', 'Specific cloud liquid water content','kg/kg','tt')
        call ncinfo(ncname(46,:),'cli', 'Specific cloud ice content','kg/kg','tt')
        call ncinfo(ncname(47,:),'plw', 'Specific precipitation liquid water content','kg/kg','tt')
        call ncinfo(ncname(48,:),'pli', 'Specific precipitation ice content','kg/kg','tt')
        do n = 1, nsv
          call ncinfo(ncname(48+7*(n-1)+1,:),trim(tracer_prop(n)%tracname), trim(tracer_prop(n)%traclong)//' specific mixing ratio', trim(tracer_prop(n)%unit),'tt')
          call ncinfo(ncname(48+7*(n-1)+2,:),trim(tracer_prop(n)%tracname)//'p', trim(tracer_prop(n)%traclong)//' tendency',trim(tracer_prop(n)%unit)//'/s)','tt')
          call ncinfo(ncname(48+7*(n-1)+3,:),trim(tracer_prop(n)%tracname)//'pt', trim(tracer_prop(n)%traclong)//' turbulence tendency',trim(tracer_prop(n)%unit)//'/s','tt')
          call ncinfo(ncname(48+7*(n-1)+4,:),trim(tracer_prop(n)%tracname)//'2r','Resolved '//trim(tracer_prop(n)%traclong)//' variance','('//trim(tracer_prop(n)%unit)//')^2','tt')
          call ncinfo(ncname(48+7*(n-1)+5,:),'w'//trim(tracer_prop(n)%tracname)//'s','SFS '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','mt')
          call ncinfo(ncname(48+7*(n-1)+6,:),'w'//trim(tracer_prop(n)%tracname)//'r','Resolved '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','mt')
          call ncinfo(ncname(48+7*(n-1)+7,:),'w'//trim(tracer_prop(n)%tracname)//'t','Total '//trim(tracer_prop(n)%traclong)//' flux',trim(tracer_prop(n)%unit)//' m/s','mt')
        end do

        if (isurf==1) then
          call open_nc(fname,  ncid,nrec,n3=kmax,ns=ksoilmax)
        else if (isurf==11) then
          call open_nc(fname,  ncid,nrec,n3=kmax,ns=kmax_soil)
        else
          call open_nc(fname,  ncid,nrec,n3=kmax)
        endif
        if (nrec == 0) then
          call define_nc( ncid, 1, tncname)
          call writestat_dims_nc(ncid)
        end if
        call define_nc( ncid, NVar, ncname)
      end if

    end if

    !$acc enter data copyin(umn, vmn, wmn, thlmn, thvmn, qtmn, qlmn, qlhmn, cfracmn, wthlsmn, wthlrmn, wthltmn, &
    !$acc&                  wthvsmn, wthvrmn, wthvtmn, wqtsmn, wqtrmn, wqttmn, wqlsmn, wqlrmn, wqltmn, &
    !$acc&                  uwtmn, vwtmn, uwrmn, vwrmn, uwsmn, vwsmn, u2mn, v2mn, w2mn, w2submn, skewmn, &
    !$acc&                  qt2mn, thl2mn, thv2mn, th2mn, ql2mn, svmmn, svpmn, svpav, svptmn, svptav, &
    !$acc&                  sv2mn, wsvsmn, wsvrmn, wsvtmn, cszav, cszmn, qlmnlast, wthvtmnlast, qlhav, &
    !$acc&                  wqlsub, wqlres, wthlsub, wthlres, wqtsub, wqtres, wthvsub, wthvres, wqttot, &
    !$acc&                  wqltot, wthltot, wthvtot, wsvsub, wsvres, wsvtot, uwres, vwres, uwsub, vwsub, &
    !$acc&                  uwtot, vwtot, umav, vmav, wmav, thvmav, thlmav, qtmav, qlmav, cfracav, u2av, v2av, &
    !$acc&                  w2av, w2subav, qt2av, thl2av, thv2av, th2av, svmav, svpav, svptav, sv2av, w3av, &
    !$acc&                  ql2av, thvmav, thmav, qlptav, thv0, sv0h, hurav, clwav, cliav, plwav, pliav, taav, &
    !$acc&                  hurmn, clwmn, climn, plwmn, plimn, tamn)

    call timer_toc('modgenstat/initgenstat')

  end subroutine initgenstat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine genstat

    use modglobal, only : rk3step, timee, dt_lim
    implicit none
    if (.not. lstat) return
    if (rk3step/=3) return

    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if

    call timer_tic('modgenstat/genstat', 0)

    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_genstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writestat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

    call timer_toc('modgenstat/genstat')

  end subroutine genstat

  subroutine do_genstat

    use modfields, only : u0,v0,w0,thl0,qt0,qt0h, &
                          ql0,ql0h,thl0h,thv0h,sv0,exnf,exnh,tmp0,presf, &
                          um, vm, wm, svm, qtm, thlm, e12m  
    use modsurfdata,only: thls,qts,ustar,thlflux,qtflux,svflux
    use modsubgriddata,only : ekm, ekh, csz
    use modglobal, only : i1,ih,j1,jh,k1,kmax,nsv,dzf,dzh,rlv,rv,rd,cp,dzhi, &
                          ijtot,cu,cv,iadv_sv,iadv_kappa,eps1,dxi,dyi,tup,tdn,lopenbc
    use modmpi,    only : comm3d,mpi_sum,mpierr,slabsum,D_MPI_ALLREDUCE
    use advec_kappa, only : halflev_kappa
    use modmicrodata, only: tuprsg, tdnrsg, iqr
    use modthermodynamics, only: qsat_tab
    implicit none

    real :: cthl,cqt,den

    integer :: i, j, k, n
    real :: tsurf, c1, c2
    real :: qs0h, t0h, ekhalf, euhalf, evhalf
    real :: wthls, wthlr, wqts, wqtr, wqls, wqlr, wthvs, wthvr
    real :: uws,vws,uwr,vwr
    real(field_r) :: upcu, vpcv
    real :: qsat, qls, ilratio
    real :: qlhav_s, wthlsub_s, wqtsub_s, wqlsub_s, wthvsub_s
    real :: uwsub_s, vwsub_s, uwres_s, vwres_s
    real :: wqlres_s, wthlres_s, wthvres_s, wqtres_s
    real :: hurav_s, clwav_s, cliav_s, plwav_s, pliav_s
    real :: wsvsub_s, wsvres_s
    real :: a_dry, b_dry, a_moist, b_moist

    call timer_tic('modgenstat/do_genstat', 1)

    !--------------------------------------------------------
    !1.0    RESET ARRAYS FOR SLAB AVERAGES
    !---    ------------------------------
    !--------------------------------------------------------

    !$acc kernels default(present)
    qlhav = 0.0
    u2av = 0.0
    v2av = 0.0
    w2av = 0.0
    w3av = 0.0
    w2subav = 0.0
    qt2av = 0.0
    thl2av = 0.0
    thv2av = 0.0
    th2av = 0.0
    ql2av = 0.0
    thvmav = 0.0
    wqtsub = 0.0
    wqtres = 0.0
    wqlsub = 0.0
    wqlres = 0.0
    umav = 0.0
    vmav = 0.0
    wmav = 0.0
    thlmav = 0.0
    thmav = 0.0
    qtmav = 0.0
    qlmav = 0.0
    wthlsub = 0.0
    wthlres = 0.0
    wthvsub = 0.0
    wthvres = 0.0
    sv2av = 0.0
    uwres = 0.0
    vwres = 0.0
    uwsub = 0.0
    vwsub = 0.0
    svmav = 0.0
    qlptav = 0.0
    wqltot = 0.0
    wqttot = 0.0
    wthvtot = 0.0
    wthltot = 0.0
    uwtot = 0.0
    vwtot = 0.0
    wsvsub = 0.0
    wsvres = 0.0
    sv2av = 0.0
    cfracav= 0.0
    cszav = 0.0
    taav = 0.0
    hurav = 0.0
    clwav = 0.0
    cliav = 0.0
    plwav = 0.0
    pliav = 0.0

    !$acc end kernels

    !-------------------------------------------------------------
    ! 2     CALCULATE SLAB AVERAGED OF FLUXES AND SEVERAL MOMENTS
    !-------------------------------------------------------------
    !------------------------------------------
    ! 2.1 SLAB AVERAGES OF PROGNOSTIC VARIABLES
    !------------------------------------------

    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          thv0(i,j,k) = (thl0(i,j,k)+rlv*ql0(i,j,k)/(cp*exnf(k))) &
                        *(1+(rv/rd-1)*qt0(i,j,k)-rv/rd*ql0(i,j,k))
        enddo
      enddo
    enddo

    !$acc parallel loop default(present) async(1)
    do k = 1, k1
      cfracav(k) = cfracav(k)+count(ql0(2:i1,2:j1,k)>0)
    end do

    !$acc wait(1)

    !$acc host_data use_device(umav, um, vmav, vm, wmav, wm, thlmav, thlm, qtmav, qtm, &
    !$acc&                     qlmav, ql0, thvmav, thv0, taav, tmp0)
    call slabsum(umav  ,1,k1,um  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(vmav  ,1,k1,vm  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(wmav  ,1,k1,wm  ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thlmav,1,k1,thlm,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qtmav ,1,k1,qtm ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(qlmav ,1,k1,ql0 ,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(thvmav,1,k1,thv0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    call slabsum(taav  ,1,k1,tmp0,2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
    !$acc end host_data
    if (nsv > 0) then
      !$acc host_data use_device(svmav, svm)
      do n = 1, nsv
        call slabsum(svmav(1:1,n),1,k1,svm(:,:,:,n),2-ih,i1+ih,2-jh,j1+jh,1,k1,2,i1,2,j1,1,k1)
      enddo
      !$acc end host_data
    end if

    !$acc kernels default(present) async(1)
    umav    = umav    / ijtot + cu
    vmav    = vmav    / ijtot + cv
    wmav    = wmav    / ijtot
    thlmav  = thlmav  / ijtot
    qtmav   = qtmav   / ijtot
    qlmav   = qlmav   / ijtot
    cfracav = cfracav / ijtot
    taav    = taav    / ijtot
    thvmav  = thvmav  / ijtot
    svmav   = svmav   / ijtot
    thmav   = thlmav + (rlv/cp)*qlmav/exnf
    cszav   = csz
    !$acc end kernels

    !-------------------------------------------
    ! 2.2 SLAB AVERAGES OF DIAGNOSTICS VARIABLES
    !-------------------------------------------
    !     TREAT LOWEST LAYER FIRST
    !-------------------------------------------

    !$acc update self(exnh(1))

    qls   = 0.0 ! hj: no liquid water at the surface
    tsurf = thls*exnh(1)+(rlv/cp)*qls
    qsat  = qts - qls

    a_dry = 1.+(rv/rd-1)*qts
    b_dry = (rv/rd-1)
    a_moist = (1.-qts+rv/rd*qsat*(1.+rlv/(rv*tsurf))) / (1.+rlv/(rv*tsurf)*rlv/(cp*tsurf)*qsat)
    b_moist = a_moist*rlv/(tsurf*cp)-1.

    if (qls < eps1) then
      c1 = a_dry
      c2 = b_dry
    else
      c1 = a_moist
      c2 = b_moist
    end if

    den = 1. + (rlv**2)*qsat/(rv*cp*(tsurf**2))
    cthl = (exnh(1)*cp/rlv)*((1-den)/den)
    cqt = 1./den

    !$acc parallel loop collapse(2) default(present) private(upcu, vpcv) &
    !$acc& reduction(+: qlhav(1), wthlsub(1), wqtsub(1), wthvsub(1), uwsub(1), vwsub(1), hurav(1), clwav(1), cliav(1), plwav(1), pliav(1)) async(1)
    do j = 2, j1
      do i = 2, i1
        qlhav(1) = qlhav(1) + ql0h(i,j,1)
        wthlsub(1) = wthlsub(1) + thlflux(i,j)
        wqtsub(1) = wqtsub(1) + qtflux (i,j)
        wthvsub(1) = wthvsub(1) + ( c1*thlflux(i,j)+c2*thls*qtflux(i,j) ) !hj: thv0 replaced by thls
        wqlsub(1) = 0.0

        !Momentum flux
        upcu = um(i, j, 1) + cu
        upcu = sign(1._field_r, upcu) * max(abs(upcu), eps1)

        uwsub(1) = uwsub(1) - (0.5 * (ustar(i,j) + ustar(i-1,j)))**2 &
                    * upcu / sqrt(upcu**2 + ((vm(i,j,1) + vm(i-1,j,1) + vm(i,j+1,1) + vm(i-1,j+1,1)) / 4. + cv)**2)

        vpcv = vm(i, j, 1) + cv
        vpcv = sign(1._field_r, vpcv) * max(abs(vpcv), eps1)

        vwsub(1) = vwsub(1) - (0.5 * (ustar(i,j) + ustar(i,j-1)))**2 &
                    * vpcv / sqrt(vpcv**2 + ((um(i,j,1) + um(i+1,j,1) + um(i,j-1,1) + um(i+1,j-1,1)) / 4. + cu)**2)

        hurav(1) = hurav(1) + 100 * (qt0(i,j,1) - ql0(i,j,1)) / qsat_tab(tmp0(i,j,1), presf(1))

        ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,1)-tdn) / (tup-tdn)))
        clwav(1) = clwav(1) + ql0(i,j,1) * ilratio
        cliav(1) = cliav(1) + ql0(i,j,1) * (1-ilratio)

        if (iqr > 0) then
           ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,1)-tdnrsg)/(tuprsg-tdnrsg)))
           plwav(1) = plwav(1) + sv0(i,j,1,iqr) * ilratio
           pliav(1) = pliav(1) + sv0(i,j,1,iqr) * (1-ilratio)
        end if
      end do
    end do

    !-------------------------------------------
    !     HIGHER LAYERS
    !-------------------------------------------

    !$acc parallel loop gang default(present) &
    !$acc& private(qlhav_s, wqlsub_s, wqlres_s, wthlsub_s, wthlres_s, wthvsub_s, wthvres_s, &
    !$acc&         wqtsub_s, wqtres_s, uwres_s, vwres_s, uwsub_s, vwsub_s, &
    !$acc&         hurav_s, clwav_s, cliav_s, plwav_s, pliav_s) async(1)
    do k = 2, kmax
      qlhav_s = 0.0
      wthlsub_s = 0.0
      wqtsub_s = 0.0
      wqlsub_s = 0.0
      wthvsub_s = 0.0
      uwsub_s = 0.0
      vwsub_s = 0.0
      uwres_s = 0.0
      vwres_s = 0.0
      wqlres_s = 0.0
      wthlres_s = 0.0
      wthvres_s = 0.0
      wqtres_s = 0.0
      clwav_s = 0.0
      cliav_s = 0.0
      plwav_s = 0.0
      pliav_s = 0.0
      hurav_s = 0.0

      !$acc loop collapse(2) &
      !$acc& private(qs0h, t0h, den, cthl, cqt, a_dry, b_dry, a_moist, b_moist, ekhalf, euhalf, evhalf, wthls, wthlr, &
      !$acc&         wqts, wqtr, wqls, wqlr, wthvs, wthvr, uwr, vwr, uws, vws) &
      !$acc& reduction(+:qlhav_s, wqlsub_s, wqlres_s, wthlsub_s, wthlres_s, wthvsub_s, wthvres_s, &
      !$acc&             wqtsub_s, wqtres_s, uwres_s, vwres_s, uwsub_s, vwsub_s, &
      !$acc&             hurav_s, clwav_s, cliav_s, plwav_s, pliav_s)
      do j = 2, j1
        do i = 2, i1
          !------------------------------------------------------
          ! Calculate ql and thv at time t0 at full and half level
          !----------------------------------------------------
          qlhav_s = qlhav_s  + ql0h(i,j,k)

          !-----------------------------------------------------------
          ! Calculate prefactors for subgrid wthv and wql fluxes
          ! at half levels
          !-----------------------------------------------------------
          qs0h  =  (qt0h(i,j,k) - ql0h(i,j,k))
          t0h   =  exnh(k)*thl0h(i,j,k) + (rlv/cp)*ql0h(i,j,k)

          den   = 1. + (rlv**2)*qs0h/(rv*cp*(t0h**2))
          cthl  = (exnh(k)*cp/rlv)*((1-den)/den)
          cqt   =  (1./den)

          a_dry = 1. + (rv/rd-1)*qt0h(i,j,k)
          b_dry = (rv/rd-1)
          a_moist = (1.-qt0h(i,j,k)+rv/rd*qs0h * (1.+rd/rv*rlv/(rd*t0h)))/den
          b_moist =  a_moist*rlv/(t0h*cp)-1.

          c1 = merge(a_moist, a_dry, ql0h(i,j,k) > 0.)
          c2 = merge(b_moist, b_dry, ql0h(i,j,k) > 0.)

          !-----------------------------------------------------------
          ! Calculate resolved and subgrid fluxes at half levels
          !-----------------------------------------------------------
          ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
          euhalf = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i-1,j,k  ) )  + &
                     dzf(k  ) * ( ekm(i,j,k-1) + ekm(i-1,j,k-1) ) ) * &
                      ( 0.25 * dzhi(k) )
          evhalf = ( dzf(k-1) * ( ekm(i,j,k  ) + ekm(i,j-1,k  ) )  + &
                     dzf(k  ) * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * &
                      ( 0.25 * dzhi(k) )

          wthls  = -ekhalf*(thl0(i,j,k)-thl0(i,j,k-1))*dzhi(k)
          wthlr  = (w0(i,j,k)-wmav(k))*thl0h(i,j,k)

          wqts   = -ekhalf*(qt0(i,j,k)-qt0(i,j,k-1))*dzhi(k)
          wqtr   = (w0(i,j,k)-wmav(k))*qt0h(i,j,k)

          wqls   = cthl*wthls+ cqt*wqts
          wqlr   = (w0(i,j,k)-wmav(k))*ql0h(i,j,k)

          wthvs  = c1*wthls + c2*thl0h(i,j,k)*wqts
          wthvr  = (w0(i,j,k)-wmav(k))*thv0h(i,j,k)

          uwr    = (w0(i,j,k)+w0(i-1,j,k)-2*wmav(k)) &
                   *((u0(i,j,k-1)+cu)*dzf(k)+(u0(i,j,k)+cu)*dzf(k-1))*(0.25*dzhi(k))
          vwr    = (w0(i,j,k)+w0(i,j-1,k)-2*wmav(k)) &
                   *((v0(i,j,k-1)+cv)*dzf(k)+(v0(i,j,k)+cv)*dzf(k-1))*(0.25*dzhi(k))
          uws    = -euhalf &
                   *((u0(i,j,k)-u0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i-1,j,k))*dxi)
          vws    = -evhalf &
                   *((v0(i,j,k)-v0(i,j,k-1))/dzh(k)+(w0(i,j,k)-w0(i,j-1,k))*dyi)

          hurav_s = hurav_s + 100 * (qt0(i,j,k) - ql0(i,j,k)) / qsat_tab(tmp0(i,j,k), presf(k))

          ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdn) / (tup-tdn)))
          clwav_s = clwav_s + ql0(i,j,k) * ilratio
          cliav_s = cliav_s + ql0(i,j,k) * (1-ilratio)

          if (iqr > 0) then
            ilratio = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))
            plwav_s = plwav_s + sv0(i,j,k,iqr) * ilratio
            pliav_s = pliav_s + sv0(i,j,k,iqr) * (1-ilratio)
          end if

          if (ql0h(i,j,k)>0) then
            wqlsub_s = wqlsub_s + wqls
          end if

          wqlres_s = wqlres_s + wqlr

          wthlsub_s = wthlsub_s + wthls
          wthlres_s = wthlres_s + wthlr

          wthvsub_s = wthvsub_s + wthvs
          wthvres_s = wthvres_s + wthvr

          wqtsub_s = wqtsub_s + wqts
          wqtres_s = wqtres_s + wqtr

          uwres_s = uwres_s + uwr
          vwres_s = vwres_s + vwr
          uwsub_s = uwsub_s + uws
          vwsub_s = vwsub_s + vws
        end do
      end do

      wqlsub(k) = wqlsub_s
      wqlres(k) = wqlres_s
      wthlsub(k) = wthlsub_s
      wthlres(k) = wthlres_s
      wthvsub(k) = wthvsub_s
      wthvres(k) = wthvres_s
      wqtsub(k) = wqtsub_s
      wqtres(k) = wqtres_s
      uwres(k) = uwres_s
      vwres(k) = vwres_s
      uwsub(k) = uwsub_s
      vwsub(k) = vwsub_s
      hurav(k) = hurav_s
      clwav(k) = clwav_s
      cliav(k) = cliav_s
      plwav(k) = plwav_s
      pliav(k) = pliav_s
    end do

    !------------
    ! 2.3 MOMENTS
    !------------
    call calc_moment(u2av, um, 2, umav, cu)
    call calc_moment(v2av, vm, 2, vmav, cv)
    call calc_moment(w2av, wm, 2, wmav)
    call calc_moment(w3av, wm, 3, wmav)
    call calc_moment(w2subav, e12m, 2, wmav)
    call calc_moment(qt2av, qtm, 2, qtmav)
    call calc_moment(thl2av, thlm, 2, thlmav)
    call calc_moment(thv2av, thv0, 2, thvmav)
    call calc_moment(th2av, thlm, 2, thmav)
    call calc_moment(ql2av, ql0, 2, qlmav)
    if (nsv > 0) then
      do n = 1, nsv
        call calc_moment(sv2av(:, n), svm(:, :, :, n), 2, svmav(:, n))
      end do
    do n = 1, nsv
        if (iadv_sv==iadv_kappa .and. .not. lopenbc) then
           call halflev_kappa(sv0(:,:,:,n),sv0h)
        else
          !$acc parallel loop collapse(3) default(present) async(1)
          do k = 2, k1
            do j = 2, j1
              do i = 2, i1    ! note: sv0h only defined and only used for k=2...
                sv0h(i,j,k) = (sv0(i,j,k,n)*dzf(k-1)+sv0(i,j,k-1,n)*dzf(k))/(2*dzh(k))
              enddo
            enddo
          enddo
        end if

        !$acc wait(1)

        !$acc parallel loop default(present) private(wsvres_s)
        do k = 2, kmax
          wsvres_s = 0.0
          !$acc loop collapse(2) reduction(+: wsvres_s)
          do j = 2, j1
            do i = 2, i1
              wsvres_s = wsvres_s + (w0(i,j,k)-wmav(k))*sv0h(i,j,k)
            end do
          end do
          wsvres(k,n) = wsvres_s
        end do


        !$acc parallel loop gang private(wsvsub_s) async(1)
        do k = 1, kmax
          wsvsub_s = 0.0
          if (k == 1) then
            !$acc loop collapse(2) reduction(+:wsvsub_s)
            do j = 2, j1
              do i = 2, i1
                wsvsub_s = wsvsub_s + svflux(i,j,n)
              end do
            end do
          else
            !$acc loop collapse(2) private(ekhalf) reduction(+: wsvsub_s)
            do j = 2, j1
              do i= 2, i1
                ekhalf = (ekh(i,j,k)*dzf(k-1)+ekh(i,j,k-1)*dzf(k))/(2*dzh(k))
                wsvsub_s= wsvsub_s-ekhalf*(sv0(i,j,k,n)-sv0(i,j,k-1,n)) / dzh(k)
              end do
            end do
          endif
          wsvsub(k,n) = wsvsub_s
        end do
      end do
    end if

    !---------------------------
    ! 3 GATHER ACCROSS MPI RANKS
    !---------------------------

    ! MPI communication
    !$acc wait(1)
    !$acc host_data use_device(qlhav, wqlsub, wqlres, wthlsub, wthlres, wthvsub, &
    !$acc&                     wthvres, uwsub, vwsub, uwres, vwres, u2av, v2av, &
    !$acc&                     w2av, w3av, w2subav, qt2av, thl2av, thv2av, th2av, &
    !$acc&                     ql2av, qlptav, sv2av, wsvsub, wsvres, cfracav, &
    !$acc&                     hurav, clwav, cliav, plwav, pliav)
    call D_MPI_ALLREDUCE(qlhav, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqlsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqlres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthlres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wqtres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthvsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(wthvres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(uwsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(vwsub, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(uwres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(vwres, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(u2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(v2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w3av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(w2subav, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qt2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(thl2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(thv2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(th2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(ql2av, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(qlptav, k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(cfracav,k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(hurav,k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(clwav,k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(cliav,k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(plwav,k1, MPI_SUM, comm3d,mpierr)
    call D_MPI_ALLREDUCE(pliav,k1, MPI_SUM, comm3d,mpierr)

    if (nsv > 0) then
      do n = 1, nsv
        call D_MPI_ALLREDUCE(sv2av(:,n),k1, MPI_SUM, comm3d,mpierr)
        call D_MPI_ALLREDUCE(wsvsub(:,n), k1, MPI_SUM, comm3d,mpierr)
        call D_MPI_ALLREDUCE(wsvres(:,n), k1, MPI_SUM, comm3d,mpierr)
      end do
    end if
    !$acc end host_data

    !------------
    ! 4 NORMALIZE
    !------------
    !$acc kernels default(present) async(1)
    cfracav = cfracav / ijtot
    qlhav   = qlhav  /ijtot

    wqlsub  = wqlsub /ijtot
    wqlres  = wqlres /ijtot

    wthlsub  = wthlsub /ijtot
    wthlres  = wthlres /ijtot

    wqtsub  = wqtsub /ijtot
    wqtres  = wqtres /ijtot

    wthvsub  = wthvsub /ijtot
    wthvres  = wthvres /ijtot

    wqttot  = wqtres + wqtsub
    wqltot  = wqlres + wqlsub
    wthltot  = wthlres + wthlsub
    wthvtot  = wthvres + wthvsub

    if (nsv > 0) then
      wsvsub = wsvsub / ijtot
      wsvres = wsvres / ijtot
      wsvtot = wsvsub + wsvres
    end if

    uwres    = uwres / ijtot
    vwres    = vwres / ijtot
    uwsub    = uwsub / ijtot
    vwsub    = vwsub / ijtot
    uwtot    = uwres + uwsub
    vwtot    = vwres + vwsub

    hurav = hurav / ijtot
    clwav = clwav / ijtot
    cliav = cliav / ijtot
    plwav = plwav / ijtot
    pliav = pliav / ijtot
    !$acc end kernels

    !---------------------------------
    ! 5 ADD SLAB AVERAGES TO TIME MEAN
    !---------------------------------
    !$acc kernels default(present) async(1)
    umn     = umn     + umav
    vmn     = vmn     + vmav
    wmn     = wmn     + wmav
    thvmn   = thvmn   + thvmav
    thlmn   = thlmn   + thlmav
    qtmn    = qtmn    + qtmav
    qlmn    = qlmn    + qlmav
    cfracmn = cfracmn + cfracav
    hurmn   = hurmn   + hurav
    clwmn   = clwmn   + clwav
    climn   = climn   + cliav
    plwmn   = plwmn   + plwav
    plimn   = plimn   + pliav
    tamn    = tamn    + taav
    qlhmn   = qlhmn   + qlhav
    wthlsmn = wthlsmn + wthlsub
    wthlrmn = wthlrmn + wthlres
    wthltmn = wthltmn + wthltot
    wthvsmn = wthvsmn + wthvsub
    wthvrmn = wthvrmn + wthvres
    wthvtmn = wthvtmn + wthvtot
    wqtsmn  = wqtsmn  + wqtsub
    wqtrmn  = wqtrmn  + wqtres
    wqttmn  = wqttmn  + wqttot
    wqlsmn  = wqlsmn  + wqlsub
    wqlrmn  = wqlrmn  + wqlres
    wqltmn  = wqltmn  + wqltot
    uwtmn   = uwtmn   + uwtot
    vwtmn   = vwtmn   + vwtot
    uwrmn   = uwrmn   + uwres
    vwrmn   = vwrmn   + vwres
    uwsmn   = uwsmn   + uwsub
    vwsmn   = vwsmn   + vwsub
    u2mn    = u2mn    + u2av
    v2mn    = v2mn    + v2av
    w2mn    = w2mn    + w2av
    w2submn = w2submn + w2subav
    qt2mn   = qt2mn   + qt2av
    thl2mn  = thl2mn  + thl2av
    thv2mn  = thv2mn  + thv2av
    th2mn   = th2mn   + th2av
    ql2mn   = ql2mn   + ql2av
    skewmn  = skewmn  + w3av/max(w2av**1.5,epsilon(w2av(1))) !
    cszmn   = cszmn   + cszav
    if (nsv > 0) then
      svmmn  = svmmn  + svmav
      svpmn  = svpmn  + svpav
      svptmn = svptmn + svptav
      sv2mn  = sv2mn  + sv2av
      wsvsmn = wsvsmn + wsvsub
      wsvrmn = wsvrmn + wsvres
      wsvtmn = wsvtmn + wsvtot
    end if
    !$acc end kernels

    call timer_toc('modgenstat/do_genstat')
  end subroutine do_genstat

  subroutine calc_moment(prof, var, n, mean, c_in)
    use modglobal, only: ijtot, i1,j1,k1,ih,jh

    implicit none

    integer, intent(in) :: n
    real, intent(out) :: prof(1:k1)
    real(field_r), intent(in) :: var(2-ih:i1+ih, 2-jh:j1+jh, 1:k1)
    real(field_r), optional, intent(in) :: mean(1:k1)
    real(field_r), optional, intent(in) :: c_in !< Translational velocity
    real(field_r) :: c
    real(field_r) :: prof_s
    integer :: i, j, k

    if (.not.present(c_in)) then
      c = 0.
    else
      c = c_in
    end if

    !$acc parallel loop default(present) private(prof_s) async
    do k = 1, k1
      prof_s = 0.0
      !$acc loop collapse(2) reduction(+: prof_s)
      do j = 2, j1
        do i = 2, i1
          prof_s = prof_s + (var(i, j, k) + c - mean(k))**n
        end do
      end do
      prof(k) = prof_s / ijtot
    end do
    !$acc wait

  end subroutine calc_moment

  subroutine writestat
      use modglobal, only : kmax,k1,nsv, zh,zf,rtimee,rlv,cp,cexpnr,ifoutput
      use modfields, only : presf,presh,exnf,exnh,rhof,rhobf,rhobh
      use modsubgriddata, only : csz
      use modmpi,    only : myid
      use modstat_nc, only: lnetcdf, writestat_nc
      implicit none


      real,dimension(k1,nvar) :: vars
      real,allocatable, dimension(:) :: tmn, thmn
      integer nsecs, nhrs, nminut,k,n
      real convt, convq
      character(40) :: name

      call timer_tic('modgenstat/writestat', 1)

      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)
      convt   = 86400.
      convq   = 86400*1000.
      allocate(tmn   (k1), thmn  (k1))

      !$acc kernels default(present)
      umn    = umn    /nsamples
      vmn    = vmn    /nsamples
      wmn    = wmn    /nsamples
      thvmn  = thvmn  /nsamples
      thlmn  = thlmn  /nsamples
      qtmn   = qtmn   /nsamples
      qlmn   = qlmn   /nsamples
      cfracmn= cfracmn/nsamples
      hurmn  = hurmn /nsamples
      clwmn  = clwmn /nsamples
      climn  = climn /nsamples
      plwmn  = plwmn /nsamples
      plimn  = plimn /nsamples
      tamn  =  tamn   /nsamples
      qlhmn  = qlhmn  /nsamples


      wthlsmn = wthlsmn/nsamples
      wthlrmn = wthlrmn/nsamples
      wthltmn = wthltmn/nsamples

      wqtsmn = wqtsmn/nsamples
      wqtrmn = wqtrmn/nsamples
      wqttmn = wqttmn/nsamples

      wqlsmn = wqlsmn/nsamples
      wqlrmn = wqlrmn/nsamples
      wqltmn = wqltmn/nsamples

      wthvsmn = wthvsmn/nsamples
      wthvrmn = wthvrmn/nsamples
      wthvtmn = wthvtmn/nsamples

      uwtmn  = uwtmn /nsamples
      vwtmn  = vwtmn /nsamples
      uwrmn  = uwrmn /nsamples
      vwrmn  = vwrmn /nsamples
      uwsmn  = uwsmn /nsamples
      vwsmn  = vwsmn /nsamples

      w2mn     = w2mn   /nsamples
      skewmn   = skewmn /nsamples
      w2submn  = w2submn/nsamples
      qt2mn    = qt2mn  /nsamples
      v2mn     = v2mn   /nsamples
      u2mn     = u2mn   /nsamples
      thl2mn   = thl2mn /nsamples
      thv2mn   = thv2mn /nsamples
      th2mn    = th2mn  /nsamples
      ql2mn    = ql2mn  /nsamples

      if (nsv > 0) then
        svmmn  = svmmn  / nsamples
        svpmn  = svpmn  / nsamples
        svptmn = svptmn / nsamples
        sv2mn  = sv2mn  / nsamples
        wsvsmn = wsvsmn / nsamples
        wsvrmn = wsvrmn / nsamples
        wsvtmn = wsvtmn / nsamples
      end if

      cszmn = cszmn / nsamples
      !$acc end kernels


  !     ------------------------------------------
  !     2.0  Construct other time averaged fields
  !     ------------------------------------------

      !$acc kernels default(present) copy(thmn, tmn)
      thmn = thlmn + (rlv/cp)*qlmn/exnf
      tmn  = thmn*exnf
      !$acc end kernels

      !$acc update self(umn, vmn, wmn, thvmn, thlmn, qtmn, qlmn, cfracmn, qlhmn, &
      !$acc&            wthlsmn, wthlrmn, wthltmn, wqtsmn, wqtrmn, wqttmn, &
      !$acc&            wqlsmn, wqlrmn, wqltmn, wthvsmn, wthvrmn, wthvtmn, &
      !$acc&            uwtmn, vwtmn, uwrmn, vwrmn, uwsmn, vwsmn, w2mn, skewmn, &
      !$acc&            w2submn, qt2mn, v2mn, u2mn, thl2mn, thv2mn, th2mn, ql2mn, &
      !$acc&            cszmn, cfracmn, hurmn, clwmn, climn, plwmn, plimn, tamn)

      !$acc update self(svmmn, svpmn, svptmn, sv2mn, wsvsmn, wsvrmn, wsvtmn) if(nsv > 0)

  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='field.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/2A/2A)') &
          '#--------------------------------------------------------' &
          ,'#LEV  HGHT    PRES    TEMP       TH_L     THETA      TH_V     ' &
          ,'  QT_AV      QL_AV      U       V   CLOUD FRACTION  CS' &
          ,'#      (M)    (MB)   (----------- (KELVIN) ---------------)    ' &
          ,'(----(G/KG)------)  (--- (M/S ---)   (-----------)  (---)'
      do k=1,kmax
        write(ifoutput,'(I3,F10.2,F7.1,5F10.4,F12.5,3F11.4,F11.5)') &
            k, &
            zf    (k), &
            presf (k)/100., &
            tmn   (k), &
            thlmn (k), &
            thmn  (k), &
            thvmn (k), &
            qtmn  (k)*1000., &
            qlmn  (k)*1000., &
            umn   (k), &
            vmn   (k), &
            cfracmn(k), &
            cszmn(k)
      end do
      close (ifoutput)

  !     -------------------------
  !     6.5   write the fluxes
  !     -----------------------------------------------------------------------

      open (ifoutput,file='flux1.'//cexpnr,position='append')

      write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
          ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(2A/A/2A/2A/2A/2A)') &
            '#---------------------------------------------------------------' &
          ,'---------------------------------)' &
          ,'#                                                               ' &
          ,'#                  |                                  TURBULENT ' &
          ,'FLUXES                           |                             ' &
          ,'#LEV HEIGHT  PRES  |   WTHL_SUB    WTHL_RES    WTHL_TOT        WQ' &
          ,'T_SUB      WQT_RES     WQT_TOT   ' &
          ,'#     (M)    (MB) | (----------   (K M/S)   --------------)     ' &
          ,'(---------- (M/S)   -------------)                             ' &
          ,'#---------------------------------------------------------------' &
          ,'---------------------------------)'

      write(ifoutput,'(I3,F10.2,F7.1,6E13.5)') &
            (k, &
            zh       (k), &
            presh     (k)/100., &
            wthlsmn   (k)                              , &
            wthlrmn   (k)                              , &
            wthltmn   (k)                              , &
            wqtsmn   (k)                              , &
            wqtrmn   (k)                              , &
            wqttmn   (k),&
             k=1,kmax)



      close(ifoutput)

      open (ifoutput,file='flux2.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#--------------------------------------------------------' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A/A/3A/4A/2A)') &
            '#(------------------------------------------------------------)' &
          ,'#                                                             ' &
          ,'#                                     TURBULENT FLUXES        ' &
          ,'#LEV HEIGHT  PRES  |    UW_TOT      VW_TOT       UW_SGS   ' &
          ,'   VW_SGS       UW_RES       VW_RES ' &
          ,'      WTH_TOT       WQ_L       WTHV_SUB     WTHV_RES     WTHV_TOT' &
          ,'#     (M)    (MB) |   ' &
          ,'(---------------------------- (M/S)^2  ---------------------------------)' &
          ,'   (-(K M/S)-)   (-(' &
          ,'M/S)-)     (----------  (K M/S)  -----------)' &
          ,'#(------------------------------------------------------------' &
          ,'-----------------------------------------)'

      write(ifoutput,'(I3,F10.2,F7.1,11E13.5)') &
          (k, &
            zh       (k), &
            presh     (k)/100., &
            uwtmn    (k)                              , &
            vwtmn    (k)                              , &
            uwsmn    (k)                              , &
            vwsmn    (k)                              , &
            uwrmn    (k)                              , &
            vwrmn    (k)                              , &
            wthltmn   (k) + wqltmn(k)*(rlv/cp)/exnh(k) , &
            wqltmn   (k)                              , &
            wthvsmn   (k)                              , &
            wthvrmn   (k)                              , &
            wthvtmn   (k)                              , &
            k=1,kmax)


      close(ifoutput)

      open (ifoutput,file='moments.'//cexpnr,position='append')

      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
        '#--------------------------------------------------------'      &
        ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
        ,nhrs,':',nminut,':',nsecs      &
        ,'   HRS:MIN:SEC AFTER INITIALIZATION '

      write (ifoutput,'(3A/3A/3A)') &
            '#----------------------------------------------------' &
            ,'---------------------------------------------------' &
            ,'------------------------------' &
            ,'#  LEV   HGHT   PRES     THL**2       THV**2         ' &
            ,'TH**2       QT**2       U*U    ' &
            ,'  V*V     HGHT     W*W     SKEWW     SFS-TKE' &
            ,'#        (M)   (MB)     (--------------(K*K)--------' &
            ,'------)     (G/KG)^2     (-(M/S)' &
            ,'^2)-)      (M)     (M/S)^2  ()        (M/S)^2'

      write(ifoutput,'(I5,F10.2,F7.1,4E13.5,2F9.4,F10.2,3F9.4)') &
            (k, &
            zf    (k), &
            presf (k)/100., &
            thl2mn(k), &
            thv2mn(k), &
            th2mn (k), &
            qt2mn (k)*1e6, &
!             qs2mn (k)*1e6, &
            u2mn  (k), &
            v2mn  (k), &
            zh    (k), &
            w2mn  (k), &
            skewmn(k), &
            w2submn(k), &
            k=1,kmax)

      close(ifoutput)


  !----   Write information about scalar field and its tendencies ------

      do n=1,nsv
        name = 'svnnnfld.'//cexpnr
        write (name(3:5),'(i3.3)') n ! TODO apply tracer_props
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '

        write (ifoutput,'(2A/A/A/A,I2.2,A,/A/A/2A)') &
            '#--------------------------------------------------------- ' &
            ,'--------------' &
            ,'#               --FIELD &  T E N D E N C I E S  --------    ' &
            ,'#                                                           ' &
            ,'#                  |        SCALAR(' &
            ,n & ! TODO apply tracer_props
            ,')                  |' &
            ,'# LEV HEIGHT  PRES      SV(1)      TURB    TOTAL  |' &
            ,'#      (M)   (MB)  |       -----  (KG/KG/DAY) ----- |' &
            ,'#----------------------------------------------------------' &
            ,'-------------'
        write(ifoutput,'(I4,2F10.2,E14.5e3,2F10.2,E14.5E3)') &
            (k, &
              zf       (k), &
              presf    (k)/100., &
              svmmn    (k,n), &
              svptmn   (k,n)  *convt, &
              svpmn    (k,n)  *convt, &
              sv2mn    (k,n), &
              k=1,kmax)

        close(ifoutput)



  !        -----------------------
  !        Write the scalar fluxes
  !        -----------------------
        name = 'svnnnflx.'//cexpnr
        write (name(3:5),'(i3.3)') n ! TODO apply tracer_props
        open (ifoutput,file=name,position='append')

        write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
            '#-------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#',(timeav),'--- AVERAGING TIMESTEP --- ' &
            ,nhrs,':',nminut,':',nsecs &
            ,'   HRS:MIN:SEC AFTER INITIALIZATION '
        write (ifoutput,'(2A/A/A,I2.2,A,/A/A/2A)') &
            '#---------------------------------------------------------------' &
            ,'---------------------------------)' &
            ,'#                                                               ' &
            ,'#                 |         TURBULENT FLUXES (SV=',n,')       |' & ! TODO apply tracer_props
            ,'#LEV HEIGHT  PRES       WSV_SUB     WSV_RES     WSV_TOT' &
            ,'#     (M)    (MB) | (----------   (KG/KG M/S)   --------------) ' &
            ,'#---------------------------------------------------------------' &
            ,'---------------------------------)'

        write(ifoutput,'(I3,2F10.2,3E14.5E3)') &
            (k, &
              zh       (k), &
              presh    (k)/100., &
              wsvsmn   (k,n)                            , &
              wsvrmn   (k,n)                            , &
              wsvtmn   (k,n)                            , &
              k=1,kmax)

         close(ifoutput)

      end do
      if (lnetcdf) then
        vars(:, 1)=rhof
        vars(:, 2)=rhobf
        vars(:, 3)=rhobh
        vars(:, 4)=presh
        vars(:, 5)=umn
        vars(:, 6)=vmn
        vars(:, 7)=wmn
        vars(:, 8)=thlmn
        vars(:, 9)=thvmn
        vars(:,10)=qtmn
        vars(:,11)=qlmn
        vars(:,12)=wthlsmn
        vars(:,13)=wthlrmn
        vars(:,14)=wthltmn
        vars(:,15)=wthvsmn
        vars(:,16)=wthvrmn
        vars(:,17)=wthvtmn
        vars(:,18)=wqtsmn
        vars(:,19)=wqtrmn
        vars(:,20)=wqttmn
        vars(:,21)=wqlsmn
        vars(:,22)=wqlrmn
        vars(:,23)=wqltmn
        vars(:,24)=uwsmn
        vars(:,25)=uwrmn
        vars(:,26)=uwtmn
        vars(:,27)=vwsmn
        vars(:,28)=vwrmn
        vars(:,29)=vwtmn
        vars(:,30)=w2submn
        vars(:,31)=w2mn
        !vars(:,31)=w2submn+w2mn
        vars(:,32)=skewmn
        vars(:,33)=u2mn
        vars(:,34)=v2mn
        vars(:,35)=thl2mn
        vars(:,36)=thv2mn
        vars(:,37)=th2mn
        vars(:,38)=qt2mn
        vars(:,39)=ql2mn
        vars(:,40)=csz
        vars(:,41)=cfracmn
        vars(:,42)=hurmn
        vars(:,43)=qtmn-qlmn
        vars(:,44)=tamn
        vars(:,45)=clwmn
        vars(:,46)=climn
        vars(:,47)=plwmn
        vars(:,48)=plimn

        do n=1,nsv
          vars(:,48+7*(n-1)+1)=svmmn(:,n)
          vars(:,48+7*(n-1)+2)=svpmn(:,n)
          vars(:,48+7*(n-1)+3)=svptmn(:,n)
          vars(:,48+7*(n-1)+4)=sv2mn(:,n)
          vars(:,48+7*(n-1)+5)=wsvsmn(:,n)
          vars(:,48+7*(n-1)+6)=wsvrmn(:,n)
          vars(:,48+7*(n-1)+7)=wsvtmn(:,n)
        end do
        call writestat_nc(ncid,1,tncname,(/rtimee/),nrec,.true.)
        call writestat_nc(ncid,nvar,ncname,vars(1:kmax,:),nrec,kmax)
      end if

    end if ! end if(myid==0)

      !$acc kernels default(present)
      qlmnlast=qlmn
      wthvtmnlast=wthvtmn

      umn      = 0.
      vmn      = 0.
      wmn      = 0.
      thlmn    = 0.
      thvmn    = 0.
      qtmn     = 0.
      qlmn     = 0.
      qlhmn    = 0.
      cfracmn  = 0.
      hurmn    = 0.
      clwmn    = 0.
      climn    = 0.
      plwmn    = 0.
      plimn    = 0.
      tamn     = 0.

      wthlsmn =  0.
      wthlrmn =  0.
      wthltmn =  0.

      wthvsmn =  0.
      wthvrmn =  0.
      wthvtmn =  0.

      wqtsmn =  0.
      wqtrmn =  0.
      wqttmn =  0.

      wqlsmn =  0.
      wqlrmn =  0.
      wqltmn =  0.

      uwtmn  =  0.
      vwtmn  =  0.
      uwrmn  =  0.
      vwrmn  =  0.
      uwsmn  =  0.
      vwsmn  =  0.


      u2mn     = 0.
      v2mn     = 0.
      w2mn     = 0.
      w2submn  = 0.
      skewmn   = 0.
      qt2mn    = 0.
      thl2mn   = 0.
      thv2mn   = 0.
      th2mn    = 0.
      ql2mn    = 0.

      svmmn   = 0.
      svpmn   = 0.
      svptmn  = 0.

      sv2mn = 0.

      wsvsmn = 0.
      wsvrmn = 0.
      wsvtmn = 0.

      cszmn  = 0.
      !$acc end kernels

      deallocate(tmn, thmn)

      call timer_toc('modgenstat/writestat')
  end subroutine writestat

  subroutine exitgenstat
    use modmpi, only : myid
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none

    if(.not.(lstat)) return

    if(lnetcdf .and. myid==0) call exitstat_nc(ncid)

    deallocate(umn       ,vmn, wmn)
    deallocate(thlmn        ,thvmn )
    deallocate(qtmn         ,qlmn  ,  qlhmn, cfracmn, hurmn)
    deallocate(clwmn,climn,plwmn,plimn)
    deallocate(wthlsmn ,wthlrmn ,wthltmn)
    deallocate(wthvsmn ,wthvrmn ,wthvtmn)
    deallocate(wqlsmn ,wqlrmn ,wqltmn)
    deallocate(wqtsmn ,wqtrmn ,wqttmn)
    deallocate(wsvsmn ,wsvrmn,wsvtmn)
    deallocate(uwtmn,vwtmn,uwrmn,vwrmn,uwsmn,vwsmn )
    deallocate(w2mn       , skewmn )
    deallocate(w2submn    )
    deallocate(u2mn       , v2mn  ,     qt2mn)
    deallocate(thl2mn     , thv2mn,     th2mn,     ql2mn)
    deallocate(svmmn,svptmn,svplsmn,svpmn)
    deallocate(sv2mn)

    deallocate(umav )
    deallocate(vmav )
    deallocate(wmav )
    deallocate(thlmav )
    deallocate(thmav )
    deallocate(qtmav )
    deallocate(qlmav )
    deallocate(cfracav )
    deallocate(hurav )
    deallocate(clwav,cliav,plwav,pliav)
    deallocate(svmav )
    deallocate(taav )
    deallocate(uptav)
    deallocate(vptav)

    deallocate(thptav)
    deallocate(qlptav)
    deallocate(uwtot )
    deallocate(vwtot )
    deallocate(uwres )
    deallocate(vwres )
    deallocate(uwsub )
    deallocate(vwsub )
    deallocate(wthlsub)
    deallocate(wthlres)
    deallocate(wthltot)
    deallocate(wqtsub)
    deallocate(wqtres)
    deallocate(wqttot)
    deallocate(wqlsub)
    deallocate(wqlres)
    deallocate(wqltot)
    deallocate(wthvtot)
    deallocate(wthvres)
    deallocate(wthvsub)
    deallocate(wsvsub)
    deallocate(wsvres)
    deallocate(wsvtot)
    deallocate(thvhav)
    deallocate(th0av)
    deallocate(svpav)
    deallocate(svptav)

    deallocate(cszmn)
    deallocate(cszav)

    deallocate(qlmnlast)
    deallocate(wthvtmnlast)

    deallocate( qlhav)
    deallocate( u2av    , &
              v2av    , &
              w2av    , &
              w3av    , &
              w2subav , &
              qt2av   , &
              thl2av  , &
              thv2av  , &
              th2av   , &
              ql2av   , &
              sv2av   )
    deallocate(thv0)
    deallocate(thvmav)
    deallocate(sv0h)

  end subroutine exitgenstat

end module modgenstat
