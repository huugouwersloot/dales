!> \file modbulkmicro.f90

!>
!!  Bulk microphysics.
!>
!! Calculates bulk microphysics using a two moment scheme.
!! \see  Seifert and Beheng (Atm. Res., 2001)
!! \see  Seifert and Beheng (Met Atm Phys, 2006)
!! \see  Stevens and Seifert (J. Meteorol. Soc. Japan, 2008)  (rain sedim, mur param)
!! \see  Seifert (J. Atm Sc., 2008) (rain evap)
!! \see  Khairoutdinov and Kogan (2000) (drizzle param : auto, accr, sedim, evap)
!!  \author Olivier Geoffroy, K.N.M.I.
!!  \author Margreet van Zanten, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \par Revision list
!! \todo documentation
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


module modbulkmicro
!   Amount of liquid water is splitted into cloud water and precipitable
!   water (as such it is a two moment scheme). Cloud droplet number conc. is
!   fixed in place and time.
!
!   same rhof value used for some diagnostics calculation (in modbulkmicrostat, modtimestat)
!
!   Cond. sampled timeav averaged profiles are weighted with fraction of condition,
!   similarly as is done in sampling.f90
!
!   bulkmicro is called from *modmicrophysics*
!*********************************************************************
  use modprecision, only : field_r
  use modtimer
  use modmicrodata, only: qrbase, qrroof, qcbase, qcroof
  implicit none
  private
  public initbulkmicro, exitbulkmicro, bulkmicro

  real :: gamma25
  real :: gamma3
  real :: gamma35
  contains

!> Initializes and allocates the arrays
  subroutine initbulkmicro
    use modglobal, only : i1,j1,k1,ih,jh
    use modmicrodata, only : lacz_gamma, Nr, Nrp, qr, qrp, thlpmcr, &
                             qtpmcr, Dvr, xr, mur, &
                             lbdr, iqr, inr, &
                             precep, qrmask, qcmask
    use modtracers,   only: add_tracer
    implicit none

    ! Setup two tracers for precipitation
    call add_tracer("qr", long_name="rain water mixing ratio", &
                    unit="kg/kg", lmicro=.true., isv=iqr) 
    
    call add_tracer("Nr", long_name="rain droplet number concentration", &
                    unit="1/m^3", lmicro=.true., isv=inr)

                                        ! Fields accessed by:
    allocate(Nr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,qr       (2:i1,2:j1,k1)  & ! dobulkmicrostat, dosimpleicestat
            ,Nrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,qrp      (2:i1,2:j1,k1)  & ! bulkmicrotend, simpleicetend
            ,Dvr      (2:i1,2:j1,k1)  & ! dobulkmicrostat
            ,precep   (2:i1,2:j1,k1)  ) ! dobulkmicrostat, dosimpleicestat, docape

    allocate(thlpmcr  (2:i1,2:j1,k1)  & !
            ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1) & ! ghost cells added here for modvarbudget
            ,xr       (2:i1,2:j1,k1)  & !
            ,mur      (2:i1,2:j1,k1)  & !
            ,lbdr     (2:i1,2:j1,k1)  & !
            ,qrmask   (2:i1,2:j1,k1)  & !
            ,qcmask   (2:i1,2:j1,k1)  )

    gamma25=lacz_gamma(2.5)
    gamma3=2.
    gamma35=lacz_gamma(3.5)

    !$acc enter data copyin(Nr, qr, Nrp, qrp, Dvr, precep, &
    !$acc&                  thlpmcr, qtpmcr, xr, mur, lbdr, qrmask, qcmask)

  end subroutine initbulkmicro

!> Cleaning up after the run
  subroutine exitbulkmicro
  !*********************************************************************
  ! subroutine exitbulkmicro
  !*********************************************************************
    use modmicrodata, only : Nr,Nrp,qr,qrp,thlpmcr,qtpmcr, &
                             Dvr,xr,mur,lbdr, &
                             precep,qrmask,qcmask
    implicit none

    !$acc exit data delete(Nr, qr, Nrp, qrp, Dvr, precep, &
    !$acc&                 thlpmcr, qtpmcr, xr, mur, lbdr, qrmask, qcmask)

    deallocate(Nr,Nrp,qr,qrp,thlpmcr,qtpmcr)
    deallocate(Dvr,xr,mur,lbdr)
    deallocate(precep,qrmask,qcmask)

  end subroutine exitbulkmicro

!> Calculates the microphysical source term.
  subroutine bulkmicro
    use modglobal, only : i1,j1,kmax,k1,rdt,rk3step,timee,rlv,cp
    use modfields, only : sv0,svm,svp,qtp,thlp,ql0,exnf,rhof
    use modbulkmicrostat, only : bulkmicrotend
    use modmpi,    only : myid
    use modmicrodata, only : Nr, qr, Nrp, qrp, thlpmcr, qtpmcr, delt, &
                             l_sedc, l_mur_cst, l_lognormal, l_rain, &
                             qrmask, qrmin, qcmask, qcmin, &
                             mur_cst, inr, iqr, l_sb
    use bulkmicro_sb, only: do_bulkmicro_sb
    use bulkmicro_kk, only: do_bulkmicro_kk
    implicit none
    integer :: i, j, k
    real :: qrtest,nr_cor,qr_cor
    real :: qrsum_neg, qrsum, Nrsum_neg, Nrsum

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          Nr(i,j,k) = sv0(i,j,k,inr)
          qr(i,j,k) = sv0(i,j,k,iqr)
          Nrp(i,j,k)     = 0.0
          qrp(i,j,k)     = 0.0
          thlpmcr(i,j,k) = 0.0
          qtpmcr(i,j,k)  = 0.0
        enddo
      enddo
    enddo

    delt = rdt/ (4. - dble(rk3step))

    if (timee.eq.0 .and. rk3step.eq.1 .and. myid.eq.0) then
      write(*,*) 'l_lognormal',l_lognormal
      write(*,*) 'rhof(1)', rhof(1),' rhof(10)', rhof(10)
      write(*,*) 'l_mur_cst',l_mur_cst,' mur_cst',mur_cst
      write(*,*) 'nuc = param'
    endif

    !*********************************************************************
    ! remove neg. values of Nr and qr
    !*********************************************************************
    if (l_rain) then
      qrsum_neg = 0.0
      qrsum = 0.0
      Nrsum_neg = 0.0
      Nrsum = 0.00
      !$acc parallel loop collapse(3) default(present) reduction(+: qrsum_neg, qrsum, Nrsum_neg, Nrsum)
      do k = 1, k1
        do j = 2, j1
          do i = 2, i1
            qrsum = qrsum + qr(i,j,k)
            Nrsum = Nrsum + Nr(i,j,k)
            if (qr(i,j,k) < 0.0) then
              qrsum_neg = qrsum_neg + qr(i,j,k)
              qr(i,j,k) = 0.0
            end if
            if (Nr(i,j,k) < 0.0) then
              Nrsum_neg = Nrsum_neg + Nr(i,j,k)
              Nr(i,j,k) = 0.0
            end if
          enddo
        enddo
      enddo

      ! LE: Commenting those out for now, popping up too often.
      !if ( -qrsum_neg > 0.000001*qrsum) then
      !  write(*,*)'amount of neg. qr thrown away is too high  ',timee,' sec'
      !end if
      !if ( -Nrsum_neg > 0.000001*Nrsum) then
      !   write(*,*)'amount of neg. Nr thrown away is too high  ',timee,' sec'
      !end if
    end if   ! l_rain

    !*********************************************************************
    ! Find gridpoints where the microphysics scheme should run
    !*********************************************************************

#if defined(DALES_GPU)
    ! Faster with OpenACC acceleration as it enables collapse(3)
    qrbase = k1 + 1
    qrroof = 1 - 1
    qcbase = k1 + 1
    qcroof = 1 - 1
    !$acc parallel loop collapse(3) default(present) reduction(min:qrbase,qcbase)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          ! Update mask prior to using it
          qrmask(i,j,k) = (qr(i,j,k) > qrmin .and. Nr(i,j,k) > 0.0)
          qcmask(i,j,k) = ql0(i,j,k) > qcmin
          if (qrmask(i,j,k)) then
            qrbase = min(k, qrbase)
          endif
          if (qcmask(i,j,k)) then
            qcbase = min(k, qcbase)
          endif
        enddo
      enddo
    enddo
    qrbase = max(1, qrbase)
    qcbase = max(1, qcbase)

    if (qrbase.le.k1 .or. qcbase.le.k1) then
      !$acc parallel loop collapse(3) default(present) reduction(max:qrroof,qcroof)
      do k = min(qrbase,qcbase), k1
        do j = 2, j1
          do i = 2, i1
            if (qrmask(i,j,k)) then
              qrroof = max(k, qrroof)
            endif
            if (qcmask(i,j,k)) then
              qcroof = max(k, qcroof)
            endif
          enddo
        enddo
      enddo
      qrroof = min(k1, qrroof)
      qcroof = min(k1, qcroof)
    endif
#else
    qrmask = qr.gt.qrmin.and.Nr.gt.0
    qrbase = k1 + 1
    qrroof = 1 - 1
    do k=1,kmax
      if (any(qrmask(:,:,k))) then
        qrbase = max(1, k)
        exit
      endif
    enddo
    if (qrbase.le.k1) then
      do k=kmax,qrbase,-1
        if (any(qrmask(:,:,k))) then
          qrroof = min(kmax, k)
          exit
        endif
      enddo
    endif

    qcmask = ql0(2:i1,2:j1,1:k1).gt.qcmin
    qcbase = k1 + 1
    qcroof = 1 - 1
    do k=1,kmax
      if (any(qcmask(:,:,k))) then
        qcbase = max(1, k)
        exit
      endif
    enddo
    if (qcbase.le.k1) then
      do k=kmax,qcbase,-1
        if (any(qcmask(:,:,k))) then
          qcroof = min(kmax, k)
          exit
        endif
      enddo
    endif
#endif

    ! if there is nothing to do, we can return at this point
    ! if (min(qrbase,qcbase).gt.max(qrroof,qcroof)) return

    if (l_sedc) then
      call sedimentation_cloud
    endif

    !*********************************************************************
    ! call microphysical processes subroutines
    !*********************************************************************
    if (l_rain) then
      if (l_sb) then
        call do_bulkmicro_sb
      else
        call do_bulkmicro_kk
      end if
    end if

    !*********************************************************************
    ! remove negative values and non physical low values
    !*********************************************************************
    !$acc parallel loop collapse(3) default(present) private(qr_cor, Nr_cor)
    do k = min(qrbase, qcbase), max(qrroof, qcroof)
      do j = 2, j1
        do i = 2, i1
          qr_cor = min(svp(i,j,k,iqr) + qrp(i,j,k) + (svm(i,j,k,iqr) / delt), &
                       0.0_field_r)
          Nr_cor = min(svp(i,j,k,iNr) + Nrp(i,j,k) + (svm(i,j,k,iNr) / delt), &
                       0.0_field_r)

          qrp(i,j,k) = qrp(i,j,k) - qr_cor
          Nrp(i,j,k) = Nrp(i,j,k) - Nr_cor
        end do
      end do
    end do

    call bulkmicrotend

    ! qcbase/qcroof are based on ql0.gt.qcmin and
    ! qrbase/qrroof are based on qr.gt.qrmin
    ! but we need boundaries to update qtp/thlp.
    !
    ! The difference between them comes from:
    !  * sedimentation_cloud updated qtpmcr/thlpmcr at qcbase-1
    !  * sedimentation_rain updated qrbase/qrroof,
    !    but at those levels qtmpcr/thlpmcr are either zero or
    !    already in the qc boundaries.
    if (qcbase.le.k1) qcbase = max(1, qcbase - 1)

    if (min(qrbase,qcbase) .gt. max(qrroof, qcroof)) return

    !$acc parallel loop collapse(3) default(present)
    do k = min(qrbase,qcbase), max(qrroof, qcroof)
      do j = 2, j1
        do i = 2, i1
          qtp (i,j,k) = qtp (i,j,k) + qtpmcr (i,j,k)
          thlp(i,j,k) = thlp(i,j,k) + thlpmcr(i,j,k)

          svp(i,j,k,iqr) = svp(i,j,k,iqr) + qrp(i,j,k)
          svp(i,j,k,inr) = svp(i,j,k,inr) + Nrp(i,j,k)
        enddo
      enddo
    enddo
  end subroutine bulkmicro

  !> Sedimentation of cloud water ((Bretherton et al,GRL 2007))
  !!
  !!   The sedimentation of cloud droplets assumes a lognormal DSD in which the
  !!   geometric std dev. is assumed to be fixed at 1.3.
  !! sedimentation of cloud droplets
  !! lognormal CDSD is assumed (1 free parameter : sig_g)
  !! terminal velocity : Stokes velocity is assumed (v(D) ~ D^2)
  !! flux is calc. anal.
  subroutine sedimentation_cloud
    use modglobal, only : i1,j1,rlv,cp,dzf,pi
    use modfields, only : rhof,exnf,ql0
    use modmicrodata, only : csed,c_St,rhow,sig_g,Nc_0, &
                             qtpmcr,thlpmcr,qcmask
    implicit none
    integer :: i, j, k
    real :: sedc

    call timer_tic('modbulkmicro/sedimentation_cloud', 1)

    if (qcbase .gt. qcroof) return

    csed = c_St*(3./(4.*pi*rhow))**(2./3.)*exp(5.*log(sig_g)**2.)

    !$acc parallel loop collapse(3) default(present)
    do k = qcbase, qcroof
      do j = 2, j1
        do i = 2, i1
          if (qcmask(i,j,k)) then
            sedc = csed*Nc_0**(-2./3.)*(ql0(i,j,k)*rhof(k))**(5./3.)

            !$acc atomic update
            qtpmcr(i,j,k)  = qtpmcr (i,j,k) - sedc /(dzf(k)*rhof(k))
            !$acc atomic update
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + sedc * (rlv/(cp*exnf(k)))/(dzf(k)*rhof(k))

            if (k > 1) then
              !$acc atomic update
              qtpmcr(i,j,k-1)  = qtpmcr(i,j,k-1) + sedc / (dzf(k-1)*rhof(k-1))
              !$acc atomic update
              thlpmcr(i,j,k-1) = thlpmcr(i,j,k-1) - sedc * (rlv/(cp*exnf(k-1)))/(dzf(k-1)*rhof(k-1))
            end if
          endif
        enddo
      enddo
    enddo

    call timer_toc('modbulkmicro/sedimentation_cloud')

  end subroutine sedimentation_cloud

end module modbulkmicro
