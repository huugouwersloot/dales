!> \file modlsmstat.f90
!!  Calculates the land surface model statistics


!>
!!  Calculates the land surface model statistics
!>
!! Profiles of the LSM statistics. Written to lsmstat.expnr
!! If netcdf is true, this module also writes in the profiles.expnr.nc output
!!  \author Stephan de Roode, TU Delft
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
module modlsmstat

  use modglobal, only : longint

implicit none
!private
PUBLIC :: initlsmstat, lsmstat, exitlsmstat
save
!NetCDF variables
  integer,parameter :: nvar = 5
  character(80),dimension(nvar,4) :: ncname

  real    :: dtav, timeav
  integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
  integer :: nsamples
  logical :: lstat= .false. !< switch to enable the lsmiative statistics (on/off)

!     ------

!   --------------
  real, allocatable :: gammasav(:)
  real, allocatable :: phiwav(:)
  real, allocatable :: tsoilav(:)
  real, allocatable :: lambdaav(:)
  real, allocatable :: lambdasav(:)

!
  real, allocatable :: gammasmn(:)
  real, allocatable :: phiwmn(:)
  real, allocatable :: tsoilmn(:)
  real, allocatable :: lambdamn(:)
  real, allocatable :: lambdasmn(:)

contains
!> Initialization routine, reads namelists and inits variables
  subroutine initlsmstat
    use modmpi,    only : myid,mpierr, comm3d, D_MPI_BCAST
    use modglobal, only : dtmax, ifnamopt,fname_options, ifoutput, cexpnr,dtav_glob,timeav_glob,ladaptive,dt_lim,btime,tres,lwarmstart,checknamelisterror
    use modstat_nc, only : lnetcdf,define_nc,ncinfo
    use modgenstat, only : idtav_prof=>idtav, itimeav_prof=>itimeav,ncid_prof=>ncid
    use modsurfdata,only : ksoilmax,isurf
    use modlsm, only : kmax_soil
    implicit none

    integer ierr, kdim_soil
    namelist/NAMLSMSTAT/ &
    dtav,timeav,lstat

    dtav=dtav_glob;timeav=timeav_glob
    lstat = .false.
    if (.not. (isurf ==1 .or. isurf==11)) return
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMSTAT,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMLSMSTAT')
      write(6 ,NAMLSMSTAT)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(timeav     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lstat      ,1,0,comm3d,mpierr)
    idtav = int(dtav / tres, kind=kind(idtav))
    itimeav = int(timeav / tres, kind=kind(itimeav))

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = int(itimeav / idtav)


    if(.not.(lstat)) return
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'NAMLSMSTAT: timeav must be a integer multiple of dtav'
    end if
    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'NAMLSMSTAT: dtav should be a integer multiple of dtmax'
    end if

    if (isurf == 1) then
        kdim_soil = ksoilmax
    else if (isurf == 11) then
        kdim_soil = kmax_soil
    end if

    allocate(phiwav(kdim_soil))
    allocate(tsoilav(kdim_soil))
    allocate(lambdaav(kdim_soil))
    allocate(lambdasav(kdim_soil))
    allocate(gammasav(kdim_soil))

    allocate(phiwmn(kdim_soil))
    allocate(tsoilmn(kdim_soil))
    allocate(lambdamn(kdim_soil))
    allocate(lambdasmn(kdim_soil))
    allocate(gammasmn(kdim_soil))

    phiwmn = 0.0
    tsoilmn = 0.0
    lambdamn = 0.0
    lambdasmn = 0.0
    gammasmn = 0.0

    if(myid==0 .and. .not. lwarmstart)then
      open (ifoutput,file='lsmstat.'//cexpnr,status='replace')
      close (ifoutput)
    end if

    if (lnetcdf) then
      idtav = idtav_prof
      itimeav = itimeav_prof
      tnext      = idtav+btime
      tnextwrite = itimeav+btime
      nsamples = int(itimeav / idtav)

      if (myid==0) then
        call ncinfo(ncname( 1,:),'tsoil','Soil temperature','W/m^2','tts')
        call ncinfo(ncname( 2,:),'phiw','Soil moisture content','W/m^2','tts')
        call ncinfo(ncname( 3,:),'lambda','Heat conductivity soil layer','W/m/K','tts')
        call ncinfo(ncname( 4,:),'lambdas','Soil moisture diffusivity soil layer','m^2/s','tts')
        call ncinfo(ncname( 5,:),'gammas','Soil moisture conductivity soil layer','M/s','tts')

        call define_nc( ncid_prof, NVar, ncname)
      end if

   end if

  end subroutine initlsmstat
!> General routine, does the timekeeping
  subroutine lsmstat
    use modglobal, only : rk3step,timee,dt_lim
    implicit none
    if (.not. lstat) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then
      tnext = tnext+idtav
      call do_lsmstat
    end if
    if (timee>=tnextwrite) then
      tnextwrite = tnextwrite+itimeav
      call writelsmstat
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))

  end subroutine lsmstat

!> Calculates the statistics
  subroutine do_lsmstat

    use modmpi,    only :  slabsum
    use modglobal, only : ijtot,i1,j1,i2,j2
    use modsurfdata, only : isurf,ksoilmax,tsoil,phiw,lambda,lambdas,gammas
    use modlsm, only : kmax_soil

    implicit none

    integer kdim_soil
    if (isurf == 1) then
        kdim_soil = ksoilmax
    else if (isurf == 11) then
        kdim_soil = kmax_soil
    end if

    tsoilav  = 0.
    phiwav  = 0.
    lambdaav  = 0.
    lambdasav  = 0.
    gammasav = 0.

    call slabsum(tsoilav ,1,kdim_soil,tsoil ,1,i2,1,j2,1,kdim_soil,2,i1,2,j1,1,kdim_soil)
    call slabsum(phiwav ,1,kdim_soil,phiw ,1,i2,1,j2,1,kdim_soil,2,i1,2,j1,1,kdim_soil)
    call slabsum(lambdaav ,1,kdim_soil,lambda ,1,i2,1,j2,1,kdim_soil,2,i1,2,j1,1,kdim_soil)
    call slabsum(lambdasav ,1,kdim_soil,lambdas ,1,i2,1,j2,1,kdim_soil,2,i1,2,j1,1,kdim_soil)
    call slabsum(gammasav ,1,kdim_soil,gammas ,1,i2,1,j2,1,kdim_soil,2,i1,2,j1,1,kdim_soil)
 !    ADD SLAB AVERAGES TO TIME MEAN

    phiwmn = phiwmn + phiwav/ijtot
    tsoilmn = tsoilmn + tsoilav/ijtot
    lambdamn = lambdamn + lambdaav/ijtot
    lambdasmn = lambdasmn + lambdasav/ijtot
    gammasmn = gammasmn + gammasav/ijtot

  end subroutine do_lsmstat


!> Write the statistics to file
  subroutine writelsmstat
      use modmpi,    only : myid
      use modglobal, only : cexpnr,ifoutput,rtimee
      use modstat_nc, only: lnetcdf, writestat_nc
      use modgenstat, only: ncid_prof=>ncid,nrec_prof=>nrec
      use modsurfdata,only : isurf, ksoilmax, zsoilc
      use modlsm, only : kmax_soil, z_soil

      implicit none
      real, allocatable :: vars(:,:)
      integer nsecs, nhrs, nminut,k, kdim_soil

      if (isurf == 1) then
          kdim_soil = ksoilmax
      else if (isurf == 11) then
          kdim_soil = kmax_soil
      end if

      allocate(vars(kdim_soil, nvar))

      nsecs   = nint(rtimee)
      nhrs    = int(nsecs/3600)
      nminut  = int(nsecs/60)-nhrs*60
      nsecs   = mod(nsecs,60)

      phiwmn   = phiwmn    /nsamples
      tsoilmn   = tsoilmn    /nsamples
      lambdamn   = lambdamn    /nsamples
      lambdasmn   = lambdasmn    /nsamples
      gammasmn   = gammasmn   /nsamples
  !     ----------------------
  !     2.0  write the fields
  !           ----------------

    if(myid==0)then
      open (ifoutput,file='lsmstat.'//cexpnr,position='append')
      write(ifoutput,'(//A,/A,F5.0,A,I4,A,I2,A,I2,A)') &
      '#--------------------------------------------------------'      &
      ,'#',(timeav),'--- AVERAGING TIMESTEP --- '      &
      ,nhrs,':',nminut,':',nsecs      &
      ,'   HRS:MIN:SEC AFTER INITIALIZATION '
      write (ifoutput,'(A/A/A)') &
          '#--------------------------------------------------------------------------' &
          ,'#LEV HEIGHT  T_SOIL      SOIL MOIST   HEAT COND.   MOIST DIFF.  MOIST COND.' &
          ,'#    (M)    (K)          (M^3/M^3)    (W/M/K)      (M^2/S)       (M/S)      '
      if (isurf == 1) then
        do k=1,kdim_soil
          write(ifoutput,'(I3,F8.4,F10.4,4E13.4)') &
              k, zsoilc(k),&
              tsoilmn(k),&
              phiwmn(k),&
              lambdamn(k),&
              lambdasmn(k),&
              gammasmn(k)
        end do
      else if (isurf == 11) then
        do k=1,kdim_soil
          write(ifoutput,'(I3,F8.4,F10.4,4E13.4)') &
              k, z_soil(k),&
              tsoilmn(k),&
              phiwmn(k),&
              lambdamn(k),&
              lambdasmn(k),&
              gammasmn(k)
        end do
      end if
      close (ifoutput)
      if (lnetcdf) then
        vars(:, 1) = tsoilmn
        vars(:, 2) = phiwmn
        vars(:, 3) = lambdamn
        vars(:, 4) = lambdasmn
        vars(:, 5) = gammasmn
       call writestat_nc(ncid_prof,nvar,ncname,vars(1:kdim_soil,:),nrec_prof,kdim_soil)
      end if
    end if ! end if(myid==0)

    phiwmn = 0.0
    tsoilmn = 0.0
    lambdamn = 0.0
    lambdasmn = 0.0
    gammasmn  = 0.0

    deallocate(vars)


  end subroutine writelsmstat

!> Cleans up after the run
  subroutine exitlsmstat
    implicit none

    !deallocate variables that are needed in modlsmiation

    if(.not.(lstat)) return

    deallocate(phiwav,tsoilav,lambdaav,lambdasav,gammasav)
    deallocate(phiwmn,tsoilmn,lambdamn,lambdasmn,gammasmn)



  end subroutine exitlsmstat


end module modlsmstat
