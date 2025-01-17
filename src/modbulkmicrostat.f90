!> \file modbulkmicrostat.f90
!!  Calculates profiles coming from the bulkmicrophysics


!>
!!  Calculates profiles coming from the bulkmicrophysics
!>
!! Profiles coming from the bulkmicrophysics. Written to precep.expnr for the
!! rain rates etc., and to qlptend.expnr and nptend.expnr for the
!! tendencies is rain water content, droplet number, respectively.
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Olivier Geoffroy, KNMI
!!  \author Johan van de Dussen, TU Delft
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
module modbulkmicrostat
  use modprecision, only : longint, field_r
  use modtimer

implicit none
private
PUBLIC  :: initbulkmicrostat, bulkmicrostat, exitbulkmicrostat, bulkmicrotend
save
!NetCDF variables
  integer,parameter :: nvar = 20
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real          :: dtav, timeav
  integer(kind=longint):: idtav, itimeav, tnext, tnextwrite
  integer          :: nsamples
  logical          :: lmicrostat = .false.
  integer, parameter      :: nrfields = 6  , &
                 iauto    = 2 , &
                  iaccr    = 3 , &
               ievap    = 4 , &
               ised      = 5, &
               iclip = 6
  real, allocatable, dimension(:,:)  :: Npav    , &
               Npmn    , &
               qrpav  , &
               qrpmn
  real, allocatable, dimension(:)    :: &
               precav  , &
               precmn  , &
               preccountav  , &
               preccountmn  , &
               prec_prcav  , &
               prec_prcmn  , &
               cloudcountav  , &
               cloudcountmn  , &
               raincountav  , &
               raincountmn  , &
               Nrrainav  , &
               Nrrainmn  , &
               qrav    , &
               qrmn    , &
               Dvrav  , &
               Dvrmn

  real(field_r), allocatable, dimension(:) :: tend_np, &
                                              tend_qrp

contains
!> Initialization routine, reads namelists and inits variables
subroutine initbulkmicrostat
    use modmpi,    only  : myid, comm3d, mpierr, D_MPI_BCAST
    use modglobal, only  : ifnamopt, fname_options, cexpnr, ifoutput, &
         dtav_glob, timeav_glob, ladaptive, k1, dtmax,btime,tres,lwarmstart,checknamelisterror
    use modstat_nc, only : lnetcdf,define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modmicrodata,only: imicro, imicro_bulk, imicro_sice, imicro_bulk3 !#sb3
    use modbulkmicrostat3,only:initbulkmicrostat3  ! #sb3
    implicit none
    integer      :: ierr

    namelist/NAMBULKMICROSTAT/ &
    lmicrostat, dtav, timeav
    
    ! #sb3 START
    if (imicro.eq.imicro_bulk3) then
       call initbulkmicrostat3
       return
    endif
    ! #sb3 END

    if ((imicro /=imicro_bulk) .and. (imicro /= imicro_sice)) return

    
    
    dtav  = dtav_glob
    timeav  = timeav_glob
    if(myid==0)then
      open (ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMBULKMICROSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMBULKMICROSTAT')
      write(6,NAMBULKMICROSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(lmicrostat,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav      ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(timeav    ,1,0,comm3d,mpierr)
    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = int(itimeav / idtav)

    if (.not. lmicrostat) return
    if (abs(timeav/dtav - nsamples) > 1e-4) then
      stop 'timeav must be an integer multiple of dtav (NAMBULKMICROSTAT)'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax - nint(dtav/dtmax)) > 1e-4) then
      stop 'dtav must be an integer multiple of dtmax (NAMBULKMICROSTAT)'
    end if

    allocate(Npav    (k1, nrfields), &
             Npmn    (k1, nrfields), &
             qrpav   (k1, nrfields), &
             qrpmn   (k1, nrfields))
    allocate(&
             precav    (k1)    , &
             precmn    (k1)    , &
             preccountav  (k1)    , &
             preccountmn  (k1)    , &
             prec_prcav  (k1)    , &
             prec_prcmn  (k1)    , &
             cloudcountav  (k1)    , &
             cloudcountmn  (k1)    , &
             raincountav  (k1)    , &
             raincountmn  (k1)    , &
             Nrrainav  (k1)    , &
             Nrrainmn  (k1)    , &
             qrav    (k1)    , &
             qrmn    (k1)    , &
             Dvrav    (k1)    , &
             Dvrmn    (k1))
    Npmn    = 0.0
    qrpmn    = 0.0
    precmn    = 0.0
    preccountmn  = 0.0
    prec_prcmn  = 0.0
    cloudcountmn  = 0.0
    raincountmn  = 0.0
    Nrrainmn  = 0.0
    qrmn    = 0.0
    Dvrmn    = 0.0

    allocate(tend_np(k1))
    allocate(tend_qrp(k1))
    tend_np(:) = 0.0
    tend_qrp(:) = 0.0

    !$acc enter data copyin(tend_np, tend_qrp, Npmn, qrpmn, &
    !$acc&                  Npav, qrpav, precav, preccountav, prec_prcav, &
    !$acc&                  cloudcountav, raincountav, Nrrainav, qrav, Dvrav, &
    !$acc&                  preccountmn, prec_prcmn, &
    !$acc&                  precmn, cloudcountmn, raincountmn, Nrrainmn, qrmn, Dvrmn)

    if (myid == 0 .and. .not. lwarmstart) then
      open (ifoutput,file = 'precep.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'nptend.'//cexpnr ,status = 'replace')
      close(ifoutput)
      open (ifoutput,file = 'qlptend.'//cexpnr,status = 'replace')
      close(ifoutput)
    end if

    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = int(itimeav / idtav)
      if (myid==0) then
        call nctiminfo(tncname(1,:))
        call ncinfo(ncname( 1,:),'cfrac','Cloud fraction','-','tt')
        call ncinfo(ncname( 2,:),'rainrate','Echo rain rate','W/m^2','tt')
        call ncinfo(ncname( 3,:),'preccount','Precipitation flux area fraction','-','tt')
        call ncinfo(ncname( 4,:),'nrrain','Rain droplet number concentration','#/m3','tt')
        call ncinfo(ncname( 5,:),'raincount','Rain water content area fraction','-','tt')
        call ncinfo(ncname( 6,:),'precmn','Rain rate','W/m^2','tt')
        call ncinfo(ncname( 7,:),'dvrmn','Precipitation mean diameter','m','tt')
        call ncinfo(ncname( 8,:),'qrmn','Precipitation specific humidity','kg/kg','tt')
        call ncinfo(ncname( 9,:),'npauto','Autoconversion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(10,:),'npaccr','Accretion rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(11,:),'npsed','Sedimentation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(12,:),'npevap','Evaporation rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(13,:),'npclip','Rain drop tendency due to clipping','#/m3/s','tt')
        call ncinfo(ncname(14,:),'nptot','Total rain drop tendency','#/m3/s','tt')
        call ncinfo(ncname(15,:),'qrpauto','Autoconversion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(16,:),'qrpaccr','Accretion rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(17,:),'qrpsed','Sedimentation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(18,:),'qrpevap','Evaporation rain water content tendency','kg/kg/s','tt')
        call ncinfo(ncname(19,:),'qrpclip','Rain water content tendency due to clipping','kg/kg/s','tt')
        call ncinfo(ncname(20,:),'qrptot','Total rain water content tendency','kg/kg/s','tt')
        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initbulkmicrostat

!------------------------------------------------------------------------------!
!> General routine, does the timekeeping
  subroutine bulkmicrostat
    use modglobal,    only  : rk3step, timee, dt_lim
    use modmicrodata,only: imicro, imicro_bulk3             ! #sb3
    implicit none

    
    if (imicro.eq.imicro_bulk3) return  ! #sb3 treated separately
    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return
    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if

    call timer_tic('modbulkmicrostat/bulkmicrostat', 1)

    if (timee >= tnext) then
      tnext = tnext + idtav
      call dobulkmicrostat
    end if
    if (timee >= tnextwrite) then
      tnextwrite = tnextwrite + itimeav
      call writebulkmicrostat
    end if

    call timer_toc('modbulkmicrostat/bulkmicrostat')

  end subroutine bulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for rainrate etc.
  subroutine dobulkmicrostat
    use modglobal,    only  : i1, j1, k1, ijtot
    use modmicrodata,  only  : qr,precep,Dvr,Nr,epscloud,epsqr,epsprec,imicro,imicro_bulk
    use modfields,  only  : ql0
    use modmpiinterface
    use modgpumpiinterface
    use modmpi
#if defined(_OPENACC)
    use openacc
    use modgpumpiinterface
#endif
    implicit none

    integer :: i, j, k
    real :: c_count, r_count, p_count, p_sum_cl
    real :: Nr_sum, p_sum, qr_sum, Dvr_sum_cl

    !$acc parallel loop gang default(present) private(c_count, r_count, p_count, p_sum_cl,&
    !$acc&                                            Nr_sum, p_sum, qr_sum, Dvr_sum_cl)
    do k = 1, k1
      c_count = 0.0
      r_count = 0.0
      p_count = 0.0
      p_sum_cl = 0.0
      Nr_sum = 0.0
      p_sum = 0.0
      qr_sum = 0.0
      Dvr_sum_cl = 0.0
      !$acc loop collapse(2) reduction(+:c_count, r_count, p_count, p_sum_cl,&
      !$acc&                             Nr_sum, p_sum, qr_sum, Dvr_sum_cl)
      do j = 2, j1
        do i = 2, i1
          if (ql0(i,j,k) > epscloud) then
            c_count = c_count + 1.0
          endif
          if (qr(i,j,k) > epsqr) then
            r_count = r_count + 1.0
          endif
          if (precep(i,j,k) > epsprec) then
            p_count = p_count + 1.0
            p_sum_cl = p_sum_cl + precep(i,j,k)
          endif
          Nr_sum = Nr_sum + Nr(i,j,k)
          p_sum = p_sum + precep(i,j,k)
          qr_sum = qr_sum + qr(i,j,k)
          if (imicro==imicro_bulk .and. qr(i,j,k) > epsqr) then
            Dvr_sum_cl = Dvr_sum_cl + Dvr(i,j,k)
          end if
        end do
      end do
      cloudcountav(k) = c_count
      raincountav (k) = r_count
      preccountav (k) = p_count
      prec_prcav  (k) = p_sum_cl
      Nrrainav    (k) = Nr_sum
      precav      (k) = p_sum
      qrav        (k) = qr_sum
      if (imicro==imicro_bulk) then
        Dvrav     (k) = Dvr_sum_cl
      end if
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, cloudcountav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, raincountav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, preccountav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, prec_prcav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, Dvrav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, Nrrainav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, precav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, qrav, k1, MPI_REAL8, MPI_SUM, comm3d, mpierr)

    !$acc kernels default(present)
    cloudcountmn(:) = cloudcountmn(:) +  cloudcountav(:) / ijtot
    raincountmn(:)  = raincountmn(:)  +  raincountav(:)  / ijtot
    preccountmn(:)  = preccountmn(:)  +  preccountav(:)  / ijtot
    prec_prcmn(:)   = prec_prcmn(:)   +  prec_prcav(:)   / ijtot
    Dvrmn(:)        = Dvrmn(:)        +  Dvrav(:)        / ijtot
    Nrrainmn(:)     = Nrrainmn(:)     +  Nrrainav(:)     / ijtot
    precmn(:)       = precmn(:)       +  precav(:)       / ijtot
    qrmn(:)         = qrmn(:)         +  qrav(:)         / ijtot
    !$acc end kernels

  end subroutine dobulkmicrostat

!------------------------------------------------------------------------------!
!> Performs the calculations for the tendencies etc.
  subroutine bulkmicrotend
    use modmpi,    only  : slabsum, slabsum_multi
    use modglobal,    only  : rk3step, timee, dt_lim, k1, ih, i1, jh, j1, ijtot
    use modmicrodata,  only  : qrp, Nrp
    implicit none

    integer        :: ifield = 0

    if (.not. lmicrostat)  return
    if (rk3step /= 3)  return
    if (timee == 0)    return

    if (timee < tnext .and. timee < tnextwrite) then
      dt_lim  = minval((/dt_lim, tnext - timee, tnextwrite - timee/))
      return
    end if

    call timer_tic('modbulkmicrostat/bulkmicrotend', 1)

    ifield    = mod(ifield, nrfields) + 1

    !$acc kernels default(present)
    tend_np(:) = 0.0
    tend_qrp(:) = 0.0
    !$acc end kernels

    !$acc host_data use_device(tend_np, Nrp, tend_qrp, qrp)
    call slabsum_multi(tend_np , 1,k1,Nrp  ,2,i1,2,j1,1,k1,2,i1,2,j1,1,k1, &
                       tend_qrp      ,qrp)
    !$acc end host_data

    !$acc kernels default(present)
    Npav(:,ifield)  = tend_np(:)  - sum(Npav (:,1:ifield-1),2)
    qrpav(:,ifield) = tend_qrp(:) - sum(qrpav(:,1:ifield-1),2)
    !$acc end kernels

    if (ifield == nrfields) then
      !$acc kernels default(present)
      Npmn(:,:)  = Npmn(:,:)  + Npav(:,:)  / nsamples / ijtot
      qrpmn(:,:) = qrpmn(:,:) + qrpav(:,:) / nsamples / ijtot
      Npav(:,:)  = 0.0
      qrpav(:,:) = 0.0
      !$acc end kernels
    end if

    call timer_toc('modbulkmicrostat/bulkmicrotend')

  end subroutine bulkmicrotend

!------------------------------------------------------------------------------!
!> Write the stats to file
  subroutine writebulkmicrostat
    use modmpi,    only  : myid
    use modglobal,    only  : rtimee, ifoutput, cexpnr, k1,kmax, &
              rlv, zf
    use modfields,    only  : presf,rhof
    use modstat_nc, only: lnetcdf, writestat_nc
    use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec

    implicit none
    real,dimension(k1,nvar) :: vars

    integer    :: nsecs, nhrs, nminut
    integer    :: k

    nsecs    = nint(rtimee)
    nhrs    = int (nsecs/3600)
    nminut    = int (nsecs/60)-nhrs*60
    nsecs    = mod (nsecs,60)

    !$acc update self(Npmn, qrpmn, cloudcountmn, raincountmn, preccountmn,&
    !$acc&            prec_prcmn, Dvrmn, Nrrainmn, precmn, qrmn)

    cloudcountmn(:) = cloudcountmn(:) / nsamples
    raincountmn(:)  = raincountmn(:)  / nsamples
    preccountmn(:)  = preccountmn(:)  / nsamples
    prec_prcmn(:)   = prec_prcmn(:)   / nsamples
    Dvrmn(:)        = Dvrmn(:)        / nsamples
    Nrrainmn(:)     = Nrrainmn(:)     / nsamples
    precmn(:)       = precmn(:)       / nsamples
    qrmn(:)         = qrmn(:)         / nsamples

    where (raincountmn > 0.)
      Dvrmn = Dvrmn / raincountmn
    elsewhere
      Dvrmn = 0.0
    end where
    where (preccountmn > 0.)
      prec_prcmn = prec_prcmn/preccountmn
    elsewhere
      prec_prcmn = 0.0
    end where

    if (myid == 0) then
    open (ifoutput,file='precep.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/2A/2A)')             &
      '#------------------------------------------------------------'     &
      ,'------------'                 &
      ,'#               --------   PRECIPITATION ------    '       &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   RHO(k)  PRES  |CLOUDCOVER  ECHORAINRATE  PRECCOUNT '   &
      ,'    NRRAIN      RAINCOUNT     PREC(k)     <Dvr(k)>     <qr(k)>'   &
      ,'#      (M)             (MB)  |----------  ---W/M2----   --------- '   &
      ,'    ------      ---------     -------     --------    ---------'   &
      ,'#-----------------------------------------------------------------'   &
      ,'---------------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F8.3,F7.1,8E13.5)') &
      (k          , &
      zf    (k)      , &
      rhof    (k)      , &
      presf    (k)/100.    , &
      cloudcountmn  (k)      , &
      prec_prcmn  (k)*rhof(k)*rlv  , &
      preccountmn  (k)      , &
      Nrrainmn  (k)      , &
      raincountmn  (k)      , &
      precmn    (k)*rhof(k)*rlv  , &
      Dvrmn    (k)      , &
      qrmn    (k)      , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='nptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S NRAIN ------    '     &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (#/M3/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      Npmn    (k,iauto)    , &
      Npmn    (k,iaccr)    , &
      Npmn    (k,ised)    , &
      Npmn    (k,ievap)    , &
      sum(Npmn  (k,2:nrfields))    , &
      k=1,kmax)
    close(ifoutput)

    open (ifoutput,file='qlptend.'//cexpnr,position='append')
    write(ifoutput,'(//2A,/A,F5.0,A,I4,A,I2,A,I2,A)')         &
      '#-------------------------------------------------------------'   &
      ,'---------------------------------)'           &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '         &
      ,nhrs,':',nminut,':',nsecs             &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
    write (ifoutput,'(2A/A/A/2A/A/A)')             &
      '#------------------------------------------------------------'     &
      , '------------'               &
      ,'#               --------   T E N D E N C I E S QRAIN ------    '   &
      ,'#                                                           '     &
      ,'# LEV HEIGHT   PRES  |  AUTO         ACCR          SEDIM    '     &
      ,'     EVAP         TOT '             &
      ,'#      (M)   (MB)  |  ---------   (KG/KG/S)      ----------'     &
      ,'#-----------------------------------------------------------'
    write(ifoutput,'(I4,F10.2,F7.1,5E13.5)') &
      (k          , &
      zf    (k)      , &
      presf    (k)/100.    , &
      qrpmn    (k,iauto)    , &
      qrpmn    (k,iaccr)    , &
      qrpmn    (k,ised)    , &
      qrpmn    (k,ievap)    , &
      sum(qrpmn  (k,2:nrfields))    , &
                        k=1,kmax)
    close(ifoutput)
      if (lnetcdf) then
        vars(:, 1) = cloudcountmn
        vars(:, 2) = prec_prcmn  (:)*rhof(:)*rlv
        vars(:, 3) = preccountmn  (:)
        vars(:, 4) = Nrrainmn  (:)
        vars(:, 5) = raincountmn  (:)
        vars(:, 6) = precmn    (:)*rhof(:)*rlv
        vars(:, 7) = Dvrmn    (:)
        vars(:, 8) = qrmn    (:)
        vars(:, 9) =Npmn    (:,iauto)
        vars(:,10) =Npmn    (:,iaccr)
        vars(:,11) =Npmn    (:,ised)
        vars(:,12) =Npmn    (:,ievap)
        vars(:,13) =Npmn    (:,iclip)
        do k=1,k1
        vars(k,14) =sum(Npmn  (k,2:nrfields))
        enddo
        vars(:,15) =qrpmn    (:,iauto)
        vars(:,16) =qrpmn    (:,iaccr)
        vars(:,17) =qrpmn    (:,ised)
        vars(:,18) =qrpmn    (:,ievap)
        vars(:,19) =qrpmn    (:,iclip)
        do k=1,k1
        vars(k,20) =sum(qrpmn  (k,2:nrfields))
        enddo
        call writestat_nc(ncid_prof,nvar,ncname,vars(1:kmax,:),nrec_prof,kmax)
      end if

    end if

    !$acc kernels default(present)
    cloudcountmn(:) = 0.0
    raincountmn(:)  = 0.0
    preccountmn(:)  = 0.0
    prec_prcmn(:)   = 0.0
    Dvrmn(:)        = 0.0
    Nrrainmn(:)     = 0.0
    precmn(:)       = 0.0
    qrmn(:)         = 0.0
    Npmn(:,:)         = 0.0
    qrpmn(:,:)        = 0.0
    !$acc end kernels

  end subroutine writebulkmicrostat

!------------------------------------------------------------------------------!

  subroutine exitbulkmicrostat
    use modmicrodata,only: imicro, imicro_bulk3     ! #sb3
    use modbulkmicrostat3, only:exitbulkmicrostat3  ! #sb3
    implicit none
    ! #sb3 START
    if (imicro.eq.imicro_bulk3) then
       call exitbulkmicrostat3
       return
    endif
    ! #sb3 END
    
    if (.not. lmicrostat)  return

    !$acc exit data delete(tend_np, tend_qrp, Npmn, qrpmn, &
    !$acc&                 Npav, qrpav, precav, preccountav, prec_prcav, &
    !$acc&                 cloudcountav, raincountav, Nrrainav, qrav, Dvrav, &
    !$acc&                 preccountmn, prec_prcmn, &
    !$acc&                 precmn, cloudcountmn, raincountmn, Nrrainmn, qrmn, Dvrmn)

    deallocate(Npav     , &
               Npmn     , &
               qrpav    , &
               qrpmn)
    deallocate(&
         precav    , &
         precmn    , &
         preccountav    , &
         preccountmn    , &
         prec_prcav    , &
         prec_prcmn    , &
         cloudcountav    , &
         cloudcountmn    , &
         raincountav    , &
         raincountmn    , &
         Nrrainav    , &
         Nrrainmn    , &
         qrav      , &
         qrmn      , &
         Dvrav    , &
         Dvrmn)

    deallocate(tend_np, tend_qrp)

  end subroutine exitbulkmicrostat

!------------------------------------------------------------------------------!

end module
