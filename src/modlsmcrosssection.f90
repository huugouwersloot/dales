!> \file modlsmcrosssection.f90
!!   Dumps an instantenous lsmcrosssection of the field

!>
!! Dumps an instantenous lsmcrosssection of the field.
!>
!! lsmcrosssections in the yz-plane and in the xy-plane            |
    !        of u,v,w,thl,thv,qt,ql. Written to movv_*.expnr and movh_*.expnr
!! If netcdf is true, this module leads the cross.myid.expnr.nc output

!!  \par Revision list
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
module modlsmcrosssection


  use modglobal, only : longint
  use modsurfdata,only : ksoilmax

implicit none
private
PUBLIC :: initlsmcrosssection, lsmcrosssection,exitlsmcrosssection
save
!NetCDF variables
  integer,parameter :: nvar=2
  integer :: nvar3
  integer :: ncid1 = 0
  integer :: ncid2 = 0
  integer :: ncid3 = 0
  integer :: nrec1 = 0
  integer :: nrec2 = 0
  integer :: nrec3 = 0
  integer :: crossheight
  character(4) :: cheight
  character(80) :: fname1 = 'lsmcrossxz.xxxxyxxx.xxx.nc'
  character(80) :: fname2 = 'lsmcrossxy.xxxx.xxxxyxxx.xxx.nc'
  character(80) :: fname3 = 'surfcross.xxxxyxxx.xxx.nc'
  character(80), dimension(nvar,4) :: ncname1
  character(80), dimension(1,4) :: tncname1
  character(80), dimension(nvar,4) :: ncname2
  character(80), dimension(1,4) :: tncname2
  character(80), allocatable, dimension(:,:) :: ncname3
  character(80), dimension(1,4) :: tncname3

  real    :: dtav
  integer(kind=longint) :: idtav, tnext
  logical :: lcross = .false.     !< switch for doing the lsmcrosssection (on/off)
  logical :: lcrosssoil = .false. !< switch for doing vertical soil crosssection (on/off)
  integer :: crossplane = 2       !< Location of the xz lsmcrosssection

contains
!> Initializing lsmcrosssection. Read out the namelist, initializing the variables
  subroutine initlsmcrosssection
    use modmpi,      only : myid,myidy,mpierr,comm3d,cmyid,D_MPI_BCAST
    use modglobal,   only : imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,j1,dt_lim,cexpnr,tres,btime,checknamelisterror,&
                            output_prefix
    use modstat_nc,  only : lnetcdf,open_nc, define_nc,ncinfo,nctiminfo,writestat_dims_nc
    use modsurfdata, only : isurf
    use modlsm,      only : lags
    implicit none

    integer :: ierr, ii

    namelist/NAMLSMCROSSSECTION/ &
    lcross, lcrosssoil, dtav, crossheight, crossplane

    crossheight=2
    ncid2=2
    nrec2=0

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMLSMCROSSSECTION,iostat=ierr)
       call checknamelisterror(ierr, ifnamopt, 'NAMLSMCROSSSECTION')
      write(6 ,NAMLSMCROSSSECTION)
      close(ifnamopt)
    end if

    if (lcross .and. .not. (isurf == 1 .or. isurf == 2 .or. isurf == 11)) then
       lcross = .FALSE.
       write (6,*) "Ignoring lcross, lsmcrossection currently implemented only for isurf==1, 2, or 11."
    endif

    if (lcrosssoil .and. .not. (isurf == 1 .or. isurf == 11)) then
       lcrosssoil = .FALSE.
       write (6,*) "Ignoring lcrosssoil, lsm soil crossection currently implemented only for isurf==1 or 11."
    endif

    call D_MPI_BCAST(dtav       ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcross     ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lcrosssoil ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossheight,1,0,comm3d,mpierr)
    call D_MPI_BCAST(crossplane ,1,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lcross .or. lcrosssoil)) return
    dt_lim = min(dt_lim,tnext)

    if (lcrosssoil) then
       if((crossheight>ksoilmax) .or. crossplane>j1) then
          stop 'lsmcrosssection: lsmcrosssection out of range'
       end if
       if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
          stop 'lsmcrosssection: dtav should be a integer multiple of dtmax'
       end if
       if (lnetcdf) then
          if (myidy==0) then
             fname1(12:19) = cmyid
             fname1(21:23) = cexpnr
             call nctiminfo(tncname1(1,:))
             call ncinfo(ncname1( 1,:),'tsoil', 'xz crosssection of the Soil temperature','K','t0tts')
             call ncinfo(ncname1( 2,:),'phiw', 'xz crosssection of the Soil moisture','m3/m3','t0tts')
             call open_nc(trim(output_prefix)//fname1,  ncid1,nrec1,n1=imax,ns=ksoilmax)
             if (nrec1==0) then
                call define_nc( ncid1, 1, tncname1)
                call writestat_dims_nc(ncid1)
                call define_nc( ncid1, NVar, ncname1)
             end if
          end if
          write(cheight,'(i4.4)') crossheight
          fname2(12:15) = cheight
          fname2(17:24) = cmyid
          fname2(26:28) = cexpnr
          call nctiminfo(tncname2(1,:))
          call ncinfo(ncname2( 1,:),'tsoil', 'xy crosssection of the Soil temperature','K','tt0t')
          call ncinfo(ncname2( 2,:),'phiw', 'xy crosssection of the Soil moisture','m3/m3','tt0t')
          call open_nc(trim(output_prefix)//fname2,  ncid2,nrec2,n1=imax,n2=jmax)
          if (nrec2==0) then
             call define_nc( ncid2, 1, tncname2)
             call writestat_dims_nc(ncid2)
             call define_nc( ncid2, NVar, ncname2)
          end if
       end if
    end if

    if (lnetcdf .and. lcross) then
        ! Surface values
        fname3(11:18) = cmyid
        fname3(20:22) = cexpnr
        if (isurf == 1) then
            nvar3 = 12
         else if (isurf == 2) then
            nvar3 = 8
        else if (isurf == 11) then
            nvar3 = 13
            if (lags) nvar3 = nvar3 + 2
        end if

        allocate(ncname3(nvar3,4))

        if (isurf == 1) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'Qnet','Net radiation','W/m^2','tt0t')
            call ncinfo(ncname3( 2,:),'H','Sensible heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 3,:),'LE','Latent heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 4,:),'G0','Ground heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 5,:),'tskin','Skin temperature','K','tt0t')
            call ncinfo(ncname3( 6,:),'tendskin','Skin tendency','W/m^2','tt0t')
            call ncinfo(ncname3( 7,:),'rs','Surface resistance','s/m','tt0t')
            call ncinfo(ncname3( 8,:),'ra','Aerodynamic resistance','s/m','tt0t')
            call ncinfo(ncname3( 9,:),'cliq','Fraction of vegetated surface covered with liquid water','-','tt0t')
            call ncinfo(ncname3(10,:),'Wl','Liquid water reservoir','m','tt0t')
            call ncinfo(ncname3(11,:),'rssoil','Soil evaporation resistance','s/m','tt0t')
            call ncinfo(ncname3(12,:),'rsveg','Vegetation resistance','s/m','tt0t')
            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
              call define_nc(ncid3, 1, tncname3)
              call writestat_dims_nc(ncid3)
            end if
            call define_nc(ncid3, nvar3, ncname3)
        else if (isurf == 2) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'hfss','Surface upward sensible heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 2,:),'hfls','Surface upward latent heat flux','W/m^2','tt0t')
            call ncinfo(ncname3( 3,:),'obuk', 'Obukhov length', 'm', 'tt0t')
            call ncinfo(ncname3( 4,:),'ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
            call ncinfo(ncname3( 5,:),'Cs', 'Drag coefficient for scalars', '-', 'tt0t')
            call ncinfo(ncname3( 6,:),'Cm', 'Drag coefficient for momentum','-', 'tt0t')
            call ncinfo(ncname3( 7,:),'z0h','Surface roughness length for heat','m', 'tt0t')
            call ncinfo(ncname3( 8,:),'z0m','Surface roughness length for momentum','m', 'tt0t')
            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
               call define_nc(ncid3, 1, tncname3)
               call writestat_dims_nc(ncid3)
            end if
            call define_nc(ncid3, nvar3, ncname3)
        else if (isurf == 11) then
            call nctiminfo(tncname3(1,:))
            call ncinfo(ncname3( 1,:),'H', 'Sensible heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 2,:),'LE', 'Latent heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 3,:),'G0', 'Ground heat flux', 'W/m^2', 'tt0t')
            call ncinfo(ncname3( 4,:),'tskin', 'Skin temperature', 'K', 'tt0t')
            call ncinfo(ncname3( 5,:),'obuk', 'Obukhov length', 'm', 'tt0t')
            call ncinfo(ncname3( 6,:),'ustar', 'Friction velocity', 'm/s^-1', 'tt0t')
            call ncinfo(ncname3( 7,:),'cliq', 'Fraction of vegetated surface covered with liquid water', '-', 'tt0t')
            call ncinfo(ncname3( 8,:),'wl', 'Liquid water reservoir', 'm', 'tt0t')
            call ncinfo(ncname3( 9,:),'ra', 'Aerodynamic resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(10,:),'rssoil', 'Soil evaporation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(11,:),'rsveg', 'Vegetation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(12,:),'f1', 'f1(SWD) function vegetation resistance', 's/m', 'tt0t')
            call ncinfo(ncname3(13,:),'f2_b', 'f2(theta) function soil resistance', 's/m', 'tt0t')
            if (lags) then
              call ncinfo(ncname3(14,:),'an_co2', 'Net CO2 assimilation', 'ppb m s-1', 'tt0t')
              call ncinfo(ncname3(15,:),'resp_co2', 'CO2 respiration soil', 'ppb m s-1', 'tt0t')
            end if

            call open_nc(trim(output_prefix)//fname3,  ncid3,nrec3,n1=imax,n2=jmax)
            if (nrec3==0) then
              call define_nc(ncid3, 1, tncname3)
              call writestat_dims_nc(ncid3)
            end if
            call define_nc( ncid3, nvar3, ncname3)
        end if
    end if
  end subroutine initlsmcrosssection
!>Run lsmcrosssection. Mainly timekeeping
  subroutine lsmcrosssection
    use modglobal, only : rk3step,timee,dt_lim
    use modstat_nc, only : writestat_nc
    implicit none


    if (.not. (lcross .or. lcrosssoil)) return
    if (rk3step/=3) return
    if(timee<tnext) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    tnext = tnext+idtav
    dt_lim = minval((/dt_lim,tnext-timee/))

    if (lcrosssoil) then
       call wrtvert
       call wrthorz
    end if
    if (lcross) then
       call wrtsurf
    end if
  end subroutine lsmcrosssection


!> Do the xz lsmcrosssections and dump them to file
  subroutine wrtvert
  use modglobal, only : imax,i1,cexpnr,ifoutput,rtimee
  use modsurfdata, only : tsoil, phiw
  use modmpi,    only : myidy
  use modstat_nc, only : lnetcdf, writestat_nc
  implicit none

  integer i,k

  real, allocatable :: vars(:,:,:)

  if( myidy /= 0 ) return

    open(ifoutput,file='movv_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)

    open(ifoutput,file='movv_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,crossplane,k),i=2,i1),k=1,ksoilmax)
    close(ifoutput)
    if (lnetcdf) then
      allocate(vars(1:imax,1:ksoilmax,2))
      vars(:,:,1) = tsoil(2:i1,crossplane,1:ksoilmax)
      vars(:,:,2) = phiw(2:i1,crossplane,1:ksoilmax)
      call writestat_nc(ncid1,1,tncname1,(/rtimee/),nrec1,.true.)
      call writestat_nc(ncid1,2,ncname1(1:2,:),vars,nrec1,imax,ksoilmax)
      deallocate(vars)
    end if

  end subroutine wrtvert

!> Do the xy lsmcrosssections and dump them to file
  subroutine wrthorz
    use modglobal, only : imax,jmax,i1,j1,cexpnr,ifoutput,rtimee
    use modsurfdata, only : tsoil,phiw
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    integer i,j
    real, allocatable :: vars(:,:,:)

    write(cheight,'(i4.4)') crossheight
    open(ifoutput,file='movh_tsoil.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((tsoil(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    open(ifoutput,file='movh_phiw.'//cexpnr,position='append',action='write')
    write(ifoutput,'(es12.5)') ((phiw(i,j,crossheight),i=2,i1),j=2,j1)
    close(ifoutput)

    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,2))
        vars(:,:,1) = tsoil(2:i1,2:j1,crossheight)
        vars(:,:,2) = phiw(2:i1,2:j1,crossheight)
        call writestat_nc(ncid2,1,tncname2,(/rtimee/),nrec2,.true.)
        call writestat_nc(ncid2,2,ncname2(1:2,:),vars,nrec2,imax,jmax)
        deallocate(vars)
    end if

  end subroutine wrthorz

  !> Do the xy lsmcrosssections and dump them to file
  subroutine wrtsurf
    use modglobal, only : imax,jmax,i1,j1,rtimee,cp,rlv
    use modfields, only : rhof
    use modsurfdata, only : Qnet, H, LE, G0, rs, ra, tskin, tendskin, &
                            cliq, rsveg, rssoil, Wl, isurf, obl, ustar, &
                            Cs, Cm, z0h, z0m, qtflux, thlflux
    use modlsm, only : f1, f2b, lags, an_co2, resp_co2
    use modstat_nc, only : lnetcdf, writestat_nc
    implicit none

    ! LOCAL
    real, allocatable :: vars(:,:,:)
    if (lnetcdf) then
        allocate(vars(1:imax,1:jmax,nvar3))

        if (isurf == 1) then
            vars(:,:,1) = qnet(2:i1,2:j1)
            vars(:,:,2) = h(2:i1,2:j1)
            vars(:,:,3) = le(2:i1,2:j1)
            vars(:,:,4) = g0(2:i1,2:j1)
            vars(:,:,5) = tskin(2:i1,2:j1)
            vars(:,:,6) = tendskin(2:i1,2:j1)
            vars(:,:,7) = rs(2:i1,2:j1)
            vars(:,:,8) = ra(2:i1,2:j1)
            vars(:,:,9) = cliq(2:i1,2:j1)
            vars(:,:,10) = Wl(2:i1,2:j1)
            vars(:,:,11) = rssoil(2:i1,2:j1)
            vars(:,:,12) = rsveg(2:i1,2:j1)
        else if (isurf == 2) then
            vars(:,:, 1) = rhof(1) * cp * thlflux(2:i1,2:j1)  ! H not allocated for isurf=2
            vars(:,:, 2) = rhof(1) * rlv * qtflux (2:i1,2:j1) ! LE not allocated for isurf=2
            vars(:,:, 3) = obl(2:i1,2:j1)
            vars(:,:, 4) = ustar(2:i1,2:j1)
            vars(:,:, 5) = Cs(2:i1,2:j1)
            vars(:,:, 6) = Cm(2:i1,2:j1)
            vars(:,:, 7) = z0h(2:i1,2:j1)
            vars(:,:, 8) = z0m(2:i1,2:j1)
        else if (isurf == 11) then
            vars(:,:, 1) = H(2:i1,2:j1)
            vars(:,:, 2) = LE(2:i1,2:j1)
            vars(:,:, 3) = G0(2:i1,2:j1)
            vars(:,:, 4) = tskin(2:i1,2:j1)
            vars(:,:, 5) = obl(2:i1,2:j1)
            vars(:,:, 6) = ustar(2:i1,2:j1)
            vars(:,:, 7) = cliq(2:i1,2:j1)
            vars(:,:, 8) = Wl(2:i1,2:j1)
            vars(:,:, 9) = ra(2:i1,2:j1)
            vars(:,:,10) = rssoil(2:i1,2:j1)
            vars(:,:,11) = rsveg(2:i1,2:j1)
            vars(:,:,12) = f1(2:i1,2:j1)
            vars(:,:,13) = f2b(2:i1,2:j1)
            if (lags) then
              vars(:,:,14) = an_co2(2:i1,2:j1)
              vars(:,:,15) = resp_co2(2:i1,2:j1)
            endif
        end if

        call writestat_nc(ncid3, 1, tncname3, (/rtimee/), nrec3, .true.)
        call writestat_nc(ncid3, nvar3, ncname3(1:nvar3,:), vars, nrec3, imax, jmax)

        deallocate(vars)
    end if


  end subroutine wrtsurf
!> Clean up when leaving the run
  subroutine exitlsmcrosssection
    use modstat_nc, only : exitstat_nc,lnetcdf
    use modmpi, only : myidy
    implicit none

    if(lcrosssoil .and. lnetcdf) then
      if (myidy==0) then
        call exitstat_nc(ncid1) ! xz soil
      end if
      call exitstat_nc(ncid2)   ! xy soil
    end if
    if(lcross .and. lnetcdf) then
       call exitstat_nc(ncid3) ! surface
       deallocate(ncname3)
    end if
  end subroutine exitlsmcrosssection

end module modlsmcrosssection
