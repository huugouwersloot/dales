!> \file modsimpleice.f90

!>
!!  Ice microphysics.
!>
!! Calculates ice microphysics in a cheap scheme without prognostic nr
!!  simpleice is called from *modmicrophysics*
!! \see  Grabowski, 1998, JAS
!! and Khairoutdinov and Randall, 2006, JAS
!!  \author Steef B\"oing, TU Delft
!!  \par Revision list
!
! FJ: seems the ref should be Khairoutdinov and Randall, 2003, JAS
!
! http://dx.doi.org/10.1175/1520-0469(1998)055%3C3283:TCRMOL%3E2.0.CO;2
! http://dx.doi.org/10.1175/JAS3810.1
! http://dx.doi.org/10.1175/1520-0469(2003)060<0607:CRMOTA>2.0.CO;2   - 2003
!
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


module modsimpleice2
  use modfields, only : rhobf
  use modprecision, only : field_r
  implicit none
  private
  public initsimpleice2, exitsimpleice2, simpleice2
  real :: gamb1r
  real :: gambd1r
  real :: gamb1s
  real :: gambd1s
  real :: gamb1g
  real :: gambd1g
  real :: gam2dr
  real :: gam2ds
  real :: gam2dg
  real :: gammaddr3
  real :: gammadds3
  real :: gammaddg3
  contains

!> Initializes and allocates the arrays
  subroutine initsimpleice2
    use modglobal, only : ih,i1,jh,j1,k1,lacz_gamma
    use modmicrodata, only : qr, qrp, nr, nrp, thlpmcr, qtpmcr, sed_qr, qr_spl, &
                             ilratio, rsgratio, sgratio, &
                             lambdar, lambdas, lambdag, &
                             qrmask, qcmask, precep, &
                             ccrz,ccsz,ccgz,ccrz2,ccsz2,ccgz2,&
                             bbg,bbr,bbs,ddg,ddr,dds, iqr, inr
    use modtracers, only: add_tracer
    implicit none
    integer:: i, j, k

    ! Nr is not used with simpleice. Added here temporarily for compatibility with old restart files
    call add_tracer("Nr", long_name="rain droplet number concentration", &
                    unit="1/m^3", lmicro=.true., isv=inr)

    call add_tracer("qr", long_name="Total precipitation mixing ratio", &
                    unit="kg/kg", lmicro=.true., isv=iqr) 

    allocate (qr(2:i1,2:j1,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,qrp(2:i1,2:j1,k1)       & ! qr tendency due to microphysics only, for statistics
             ,nr(2:i1,2:j1,k1)        & ! qr (total precipitation!) converted from a scalar variable
             ,nrp(2:i1,2:j1,k1)       & ! qr tendency due to microphysics only, for statistics
             ,thlpmcr(2:i1,2:j1,k1)   & ! thl tendency due to microphysics only, for statistics
             ,qtpmcr(2-ih:i1+ih,2-jh:j1+jh,k1) & ! qt tendency due to microphysics only, for statistics. Ghost cells for modvarbudget.
             ,sed_qr(2:i1,2:j1,k1)    & ! sedimentation rain droplets mixing ratio
             ,qr_spl(2:i1,2:j1,k1)    & ! time-splitting substep qr
             ,ilratio(2:i1,2:j1,k1)   & ! partition ratio cloud water vs cloud ice
             ,rsgratio(2:i1,2:j1,k1)  & ! partition ratio rain vs. snow/graupel
             ,sgratio(2:i1,2:j1,k1)   & ! partition ratio snow vs graupel
             ,lambdar(2:i1,2:j1,k1)   & ! slope parameter for rain
             ,lambdas(2:i1,2:j1,k1)   & ! slope parameter for snow
             ,lambdag(2:i1,2:j1,k1))    ! slope parameter for graupel

    allocate (qrmask(2:i1,2:j1,k1)    & ! mask for rain water
             ,qcmask(2:i1,2:j1,k1))     ! mask for cloud water

    allocate(precep(2:i1,2:j1,k1))      ! precipitation for statistics

    allocate(ccrz(k1),ccsz(k1),ccgz(k1))
    allocate(ccrz2(k1),ccsz2(k1),ccgz2(k1))

    gamb1r=lacz_gamma(bbr+1)
    gambd1r=lacz_gamma(bbr+ddr+1)
    gamb1s=lacz_gamma(bbs+1)
    gambd1s=lacz_gamma(bbs+dds+1)
    gamb1g=lacz_gamma(bbg+1)
    gambd1g=lacz_gamma(bbg+ddg+1)
    gam2dr=lacz_gamma(2.5+0.5*ddr)
    gam2ds=lacz_gamma(2.5+0.5*dds)
    gam2dg=lacz_gamma(2.5+0.5*ddg)
    gammaddr3=lacz_gamma(3.+ddr)
    gammadds3=lacz_gamma(3.+dds)
    gammaddg3=lacz_gamma(3.+ddg)

    !$acc enter data create(nrp, nr, qrp, qr, thlpmcr, qtpmcr, sed_qr, qr_spl, &
    !$acc&                  ilratio, rsgratio, sgratio, lambdar, lambdas, &
    !$acc&                  lambdag, qrmask, qcmask, precep, &
    !$acc&                  ccrz, ccsz, ccgz, ccrz2, ccsz2, ccgz2)

    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          nrp(i,j,k)=0. ! not used in this scheme
          nr (i,j,k)=0. ! set to 0 here in case the statistics use them
        enddo
      enddo
    enddo

  end subroutine initsimpleice2

!> Cleaning up after the run
  subroutine exitsimpleice2
    use modmicrodata, only : nr,nrp,qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl, &
                             ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag, &
                             qrmask,qcmask,precep, &
                             ccrz,ccsz,ccgz,ccrz2,ccsz2,ccgz2
    implicit none

    !$acc exit data delete(nrp, nr, &
    !$acc&                 qrp, qr, thlpmcr, qtpmcr, sed_qr, qr_spl, &
    !$acc&                 ilratio, rsgratio, sgratio, lambdar, lambdas, &
    !$acc&                 lambdag, qrmask, qcmask, precep, &
    !$acc&                 ccrz, ccsz, ccgz, ccrz2, ccsz2, ccgz2)

    deallocate(nr,nrp,qr,qrp,thlpmcr,qtpmcr,sed_qr,qr_spl,ilratio,rsgratio,sgratio,lambdar,lambdas,lambdag)
    deallocate(qrmask,qcmask)
    deallocate(precep)
    deallocate(ccrz,ccsz,ccgz)
    deallocate(ccrz2,ccsz2,ccgz2)
  end subroutine exitsimpleice2


!> Calculates the microphysical source term.
  subroutine simpleice2
    use modglobal, only : i1,ih,j1,jh,k1,rdt,rk3step,timee,rlv,cp,tup,tdn,pi,tmelt,kmax,dzf,dzh
    use modfields, only : sv0,svm,svp,qtp,thlp,qt0,ql0,exnf,rhof,tmp0,rhobf,qvsl,qvsi,esl,surf_rain
    use modmicrodata, only : sed_qr,qrp,&
                             aag,aar,aas,bbg,bbr,bbs,betag,betar,betas,ccg,ccr,ccs,&
                             ccgz2,ccrz2,ccsz2,ddg,ddr,dds,&
                             ccgz,ccrz,ccsz,&
                             n0rg,n0rr,n0rs,&
                             betakessi,ceffgl,ceffri,ceffrl,ceffsi,ceffsl,ceffgi,courantp,delt,&
                             qr,qtpmcr,thlpmcr,evapfactor,iqr,n0rg,n0rs,Nc_0,&
                             qr_spl,precep,&
                             qcmin,qrmin,qli0,qll0,tdnrsg,tdnsg,tuprsg,tupsg,&
                             l_berry,l_graupel,l_rain,l_warm,timekessl

    use modsimpleicestat, only : simpleicetend

    implicit none

    integer:: i,j,k
    real(field_r) :: qrsmall, qrsum
    real(field_r) :: qll,qli,ddisp,lwc,autl,tc,times,auti,aut ! autoconvert
    real(field_r) :: qrr,qrs,qrg, gaccrl,gaccsl,gaccgl,gaccri,gaccsi,gaccgi,accr,accs,accg,acc  !accrete
    real(field_r) :: ssl,ssi,ventr,vents,ventg,thfun,evapdepr,evapdeps,evapdepg,devap  !evapdep
    real(field_r) :: dt_spl,wfallmax,vtr,vts,vtg,vtf ! precipitation
    real(field_r) :: tmp_lambdar, tmp_lambdas, tmp_lambdag
    integer :: jn
    integer :: n_spl      !<  sedimentation time splitting loop

    real(field_r) :: ilratio_,lambdar_,lambdas_,lambdag_, rsgratio_, sgratio_  ! local values instead of global arrays
    logical :: qrmask_, qcmask_
    logical :: rain_present, snow_present, graupel_present   ! logicals for presence of different forms of water in the current cell
    logical :: any_qr, any_snow_graupel                      ! logicals for precense of any precipitation, and for presense of snow/graupel in the whole system

    delt = rdt/ (4. - dble(rk3step))

    wfallmax = 9.9 ! cap for fall velocity
    n_spl = ceiling(wfallmax*delt/(minval(dzf)*courantp)) ! number of sub-timesteps for precipitation
    dt_spl = delt/real(n_spl)                              ! fixed time step for precipitation sub-stepping!

    sed_qr = 0. ! reset sedimentation fluxes

    ! Density corrected fall speed parameters, see Tomita 2008
    ! rhobf is constant in time. could initialize tables in initsimpleice2, but that would requires that
    ! microphysics is initialized AFTER readinitfiles (for rhobf def), which breaks chemistry or surface
    do k = 1, k1
      ccrz(k)=ccr*(1.29/rhobf(k))**0.5
      ccsz(k)=ccs*(1.29/rhobf(k))**0.5
      ccgz(k)=ccg*(1.29/rhobf(k))**0.5

      ! these coefficients are used in evapdep - tabulated because of sqrt
      ccrz2(k) = gam2dr*.27*n0rr*sqrt(ccrz(k)/2.e-5)
      ccsz2(k) = gam2ds*.39*n0rs*sqrt(ccsz(k)/2.e-5) ! NOTE: .27 here is suspect -> .39
      ccgz2(k) = gam2dg*.27*n0rg*sqrt(ccgz(k)/2.e-5)
    end do

    ! used to check on negative qr and nr
    qrsum = 0.
    qrsmall = 0.

    ! reset microphysics tendencies (different bounds on the two loops)
    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2, j1
        do i = 2, i1
          qrp(i,j,k) = 0.
          thlpmcr(i,j,k) = 0.
        enddo
      enddo
    enddo
    !$acc parallel loop collapse(3) default(present)
    do k = 1, k1
      do j = 2-jh, j1+jh
        do i = 2-ih, i1+ih
          qtpmcr(i,j,k) = 0.
        enddo
      enddo
    enddo

    any_qr = .false.
    any_snow_graupel = .false.
    qrmask_ = .false. ! needed if l_rain is false

    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          rain_present = .false.
          snow_present = .false.
          graupel_present = .false.

          ! initialise qr
          qr(i,j,k) = sv0(i,j,k,iqr)

          ! initialise qc mask
          if (ql0(i,j,k) > qcmin) then
            qcmask_ = .true.
          else
            qcmask_ = .false.
          end if

          ! initialise qr mask and check if we are not throwing away too much rain
          if (l_rain) then
            qrsum = qrsum + qr(i,j,k)
            if (qr(i,j,k) <= qrmin) then
              qrmask_ = .false.
              if (qr(i,j,k)<0.) then
                qrsmall = qrsmall - qr(i,j,k)
                qr(i,j,k) = 0.
              end if
            else
              qrmask_= .true. ! this cell
              any_qr = .true. ! whole system
            endif
          endif


          ! logic
          !
          ! qrmask: true if cell contains rain. qr > threshold
          ! rsgratio, sgratio, lambda* calculated
          !
          ! qcmask: true if cell contains cloud - condensed water - ql > threshold
          ! ilratio, qll, qli calculated
          !
          ! qr, qrp  - rain, tendency
          ! qtpmcr   - qt tendency from microphysics
          ! thlpmcr  - thl tendency from microphysics



          !partitioning and determination of intercept parameter

          if (qrmask_ .eqv. .true.) then
            if (l_warm) then !partitioning and determination of intercept parameter
              rsgratio_=1.   ! rain vs snow/graupel partitioning
              rain_present = .true.

              sgratio_=0.   ! snow versus graupel partitioning
              lambdar_=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k))))**(1./(1.+bbr)) ! lambda rain
              !lambdas_=lambdar_ ! lambda snow    ! probably not right but they will not be used
              !lambdag_=lambdar_ ! lambda graupel
            elseif(l_graupel) then
              rsgratio_=max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg))) ! rain vs snow/graupel partitioning   rsg = 1 if t > tuprsg
              sgratio_=max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnsg)/(tupsg-tdnsg))) ! snow versus graupel partitioning    sg = 1 -> only graupel
              if (rsgratio_ > 0) then                                     ! sg = 0 -> only snow
                rain_present = .true.
                lambdar_=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k)*rsgratio_)))**(1./(1.+bbr)) ! lambda rain
              endif
              if (rsgratio_ < 1) then
                any_snow_graupel = .true.     ! whole system  ! note: this may miss snow that appears during the precipitation substepping
                if (sgratio_ > 0) then
                  graupel_present = .true.
                  lambdag_=(aag*n0rg*gamb1g/(rhof(k)*(qr(i,j,k)*(1.-rsgratio_)*sgratio_)))**(1./(1.+bbg)) ! graupel
                endif
                if (sgratio_ < 1) then
                  snow_present = .true.
                  lambdas_=(aas*n0rs*gamb1s/(rhof(k)*(qr(i,j,k)*(1.-rsgratio_)*(1.-sgratio_))))**(1./(1.+bbs)) ! snow
                endif
              endif

              ! no snow/graupel -> large lambda
            else ! rain, snow but no graupel
              rsgratio_=max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg)))   ! rain vs snow/graupel partitioning
              sgratio_=0.
              if (rsgratio_ > 0) then
                rain_present = .true.
                lambdar_=(aar*n0rr*gamb1r/(rhof(k)*(qr(i,j,k)*rsgratio_)))**(1./(1.+bbr)) ! lambda rain
              endif
              if (rsgratio_ < 1) then
                snow_present = .true.
                any_snow_graupel = .true.     ! whole system
                lambdas_=(aas*n0rs*gamb1s/(rhof(k)*(qr(i,j,k)*(1.-rsgratio_))))**(1./(1.+bbs)) ! lambda snow
              endif
              ! lambdag_=lambdas_ ! FJ: probably wrong - routines below don't always check sgratio
            end if
          endif

          ! Autoconvert
          if (qcmask_ .eqv. .true.) then
            if(l_warm) then
              ilratio_=1.
            else
              ilratio_=max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdn)/(tup-tdn)))! cloud water vs cloud ice partitioning
            endif

            ! ql partitioning - used here and in Accrete
            qll = ql0(i,j,k)*ilratio_
            qli = ql0(i,j,k) - qll

            if (l_berry .eqv. .true.) then ! Berry/Hsie autoconversion
              ! ql partitioning
              ! qll=ql0(i,j,k)*ilratio(i,j,k)
              ! qli=ql0(i,j,k)-qll

              ddisp=0.146-5.964e-2*log(Nc_0/2.e9) ! Relative dispersion coefficient for Berry autoconversion
              lwc=1.e3*rhof(k)*qll ! Liquid water content in g/kg
              autl=1./rhof(k)*1.67e-5*lwc*lwc/(5. + .0366*Nc_0/(1.e6*ddisp*(lwc+1.e-6)))
              tc=tmp0(i,j,k)-tmelt ! Temperature wrt melting point
              times=min(1.e3,(3.56*tc+106.7)*tc+1.e3) ! Time scale for ice autoconversion
              auti=qli/times
              aut = min(autl + auti,ql0(i,j,k)/delt)
              qrp(i,j,k) = qrp(i,j,k)+aut
              qtpmcr(i,j,k) = qtpmcr(i,j,k)-aut
              thlpmcr(i,j,k) = thlpmcr(i,j,k)+(rlv/(cp*exnf(k)))*aut
            else  ! Lin/Kessler autoconversion as in Khairoutdinov and Randall, 2006
              ! ql partitioning
              ! qll=ql0(i,j,k)*ilratio(i,j,k)
              ! qli=ql0(i,j,k)-qll

              autl=max(0._field_r,timekessl*(qll-qll0))
              tc=tmp0(i,j,k)-tmelt
              auti=max(0._field_r,betakessi*exp(0.025_field_r*tc)*(qli-qli0))
              aut = min(autl + auti,ql0(i,j,k)/delt)
              qrp(i,j,k) = qrp(i,j,k) + aut
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - aut
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k))) * aut
            endif
          endif

          ! Accrete
          if (qrmask_ .eqv. .true.) then
            if (qcmask_ .eqv. .true.) then ! apply mask
              ! ql partitioning - calculated in Autoconvert
              !qll=ql0(i,j,k)*ilratio(i,j,k)
              !qli=ql0(i,j,k)-qll

              ! qr partitioning
              qrr = qr(i,j,k)*rsgratio_
              qrs = qr(i,j,k)*(1.-rsgratio_)*(1.-sgratio_)
              qrg = qr(i,j,k)*(1.-rsgratio_)*sgratio_
              ! collection of cloud water by rain etc.

              accr = 0
              accs = 0
              accg = 0
              if (rain_present) then
                gaccrl=pi/4.*ccrz(k)*ceffrl*rhof(k)*qll*qrr*lambdar_**(bbr-2.-ddr)*gammaddr3/(aar*gamb1r)
                gaccri=pi/4.*ccrz(k)*ceffri*rhof(k)*qli*qrr*lambdar_**(bbr-2.-ddr)*gammaddr3/(aar*gamb1r)
                accr=(gaccrl+gaccri) !*qrr/(qrr+1.e-9)
              endif

              if (snow_present) then
                gaccsl=pi/4.*ccsz(k)*ceffsl*rhof(k)*qll*qrs*lambdas_**(bbs-2.-dds)*gammadds3/(aas*gamb1s)
                gaccsi=pi/4.*ccsz(k)*ceffsi*rhof(k)*qli*qrs*lambdas_**(bbs-2.-dds)*gammadds3/(aas*gamb1s)
                accs=(gaccsl+gaccsi) !*qrs/(qrs+1.e-9)  ! why this division? makes accr small if qr* << 1e-9
              endif    ! disable accretion if no snow present - now done with if.

              if (graupel_present) then
                gaccgl=pi/4.*ccgz(k)*ceffgl*rhof(k)*qll*qrg*lambdag_**(bbg-2.-ddg)*gammaddg3/(aag*gamb1g)
                gaccgi=pi/4.*ccgz(k)*ceffgi*rhof(k)*qli*qrg*lambdag_**(bbg-2.-ddg)*gammaddg3/(aag*gamb1g)
                accg=(gaccgl+gaccgi) !*qrg/(qrg+1.e-9)
              endif

              acc= min(accr + accs + accg, ql0(i,j,k)/delt)  ! total growth by accretion
              qrp(i,j,k) = qrp(i,j,k) + acc
              qtpmcr(i,j,k) = qtpmcr(i,j,k) - acc
              thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*acc
            end if
          end if

          ! evapdep
          if (qrmask_ .eqv. .true.) then
            ! saturation ratios
            ssl=(qt0(i,j,k)-ql0(i,j,k))/qvsl(i,j,k)
            ssi=(qt0(i,j,k)-ql0(i,j,k))/qvsi(i,j,k)
            !integration over ventilation factors and diameters, see e.g. seifert 2008
            evapdepr = 0
            evapdeps = 0
            evapdepg = 0
            thfun=1.e-7/(2.2*tmp0(i,j,k)/esl(i,j,k)+2.2e2/tmp0(i,j,k))  ! thermodynamic function


            if (rain_present) then
              !ventr=.78*n0rr/lambdar_**2 + gam2dr*.27*n0rr*sqrt(ccrz(k)/2.e-5)*lambdar_**(-2.5-0.5*ddr)
              ventr=.78*n0rr/lambdar_**2 + ccrz2(k) * lambdar_**(-2.5-0.5*ddr)
              evapdepr=(4.*pi/(betar*rhof(k)))*(ssl-1.)*ventr*thfun
            endif
            if (snow_present) then
              !vents=.78*n0rs/lambdas_**2 + gam2ds*.27*n0rs*sqrt(ccsz(k)/2.e-5)*lambdas_**(-2.5-0.5*dds)
              !vents=.78*n0rs/lambdas_**2 + ccsz2(k)*lambdas_**(-2.5-0.5*dds)
              vents=.65*n0rs/lambdas_**2 + ccsz2(k)*lambdas_**(-2.5-0.5*dds)  ! FJ corrected coefficient
              evapdeps=(4.*pi/(betas*rhof(k)))*(ssi-1.)*vents*thfun
            endif
            if (graupel_present) then
              !ventg=.78*n0rg/lambdag_**2 + gam2dg*.27*n0rg*sqrt(ccgz(k)/2.e-5)*lambdag_**(-2.5-0.5*ddg)
              ventg=.78*n0rg/lambdag_**2 + ccgz2(k)*lambdag_**(-2.5-0.5*ddg)
              evapdepg=(4.*pi/(betag*rhof(k)))*(ssi-1.)*ventg*thfun
            endif

            ! total growth by deposition and evaporation
            ! limit with qr and ql after accretion and autoconversion
            devap= max(min(evapfactor*(evapdepr+evapdeps+evapdepg),ql0(i,j,k)/delt-qrp(i,j,k)),-qr(i,j,k)/delt-qrp(i,j,k))
            qrp(i,j,k) = qrp(i,j,k) + devap
            qtpmcr(i,j,k) = qtpmcr(i,j,k) - devap
            thlpmcr(i,j,k) = thlpmcr(i,j,k) + (rlv/(cp*exnf(k)))*devap

            ! Grabowski 1998 has different coefficients here for snow
            ! also ccsz2 table needs to be updated
          end if

          ! precipitate - part 1
          qr_spl(i,j,k) = qr(i,j,k) ! prepare for sub-timestepping precipitation
                                    ! this is the first substep, using lambdas already calculated
          if (qrmask_ .eqv. .true.) then
            vtf = 0
            if (rain_present) then
              vtr=ccrz(k)*(gambd1r/gamb1r)/(lambdar_**ddr)  ! terminal velocity rain
              vtf = vtf + rsgratio_*vtr
            endif
            if (snow_present) then
              vts=ccsz(k)*(gambd1s/gamb1s)/(lambdas_**dds)  ! terminal velocity snow
              vtf = vtf + (1.-rsgratio_)*(1.-sgratio_)*vts
            endif
            if (graupel_present) then
              vtg=ccgz(k)*(gambd1g/gamb1g)/(lambdag_**ddg)  ! terminal velocity graupel
              vtf = vtf + (1.-rsgratio_)*sgratio_*vtg
            endif
            ! vtf=rsgratio_*vtr+(1.-rsgratio_)*(1.-sgratio_)*vts+(1.-rsgratio_)*sgratio(i,j,k)*vtg ! weighted
            vtf = min(wfallmax,vtf)
            !write(*,*) 'vtf', vtf

            precep(i,j,k) = vtf*qr_spl(i,j,k)
            sed_qr(i,j,k) = precep(i,j,k)*rhobf(k) ! convert to flux
          else
            precep(i,j,k) = 0.
            sed_qr(i,j,k) = 0.
          end if
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax
      do j = 2, j1
        do i = 2, i1
          qr_spl(i,j,k) = qr_spl(i,j,k) + &
                          (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl / (dzh(k+1)*rhobf(k))
        enddo
      enddo
    enddo

    surf_rain = surf_rain + sed_qr(:,:,1)*dt_spl

    ! precipitate part 2

!    write (*,*) 'any_qr:', any_qr
!    write (*,*) 'any_snow_graupel:', any_snow_graupel
!    write(*,*) 'n_spl', n_spl

    ! begin time splitting loop
    if (n_spl > 1 .and. any_qr) then
      do jn = 2, n_spl
        ! reset fluxes at each step of loop
        ! sed_qr = 0. ! not needed !?

        if (any_snow_graupel) then
          !$acc parallel loop collapse(3) default(present)
          do k = 1, kmax
            do j= 2 , j1
              do i = 2, i1
                if (qr_spl(i,j,k) > qrmin) then
                  ! re-evaluate lambda and rsgratios
                  rsgratio_ = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnrsg)/(tuprsg-tdnrsg))) ! rain vs snow/graupel partitioning   rsg = 1 if t > tuprsg
                  sgratio_ = max(0._field_r,min(1._field_r,(tmp0(i,j,k)-tdnsg) /(tupsg -tdnsg))) ! snow versus graupel partitioning    sg = 1 -> only graupel

                  !these ifs are here to avoid performing the power calculations unless they are going to be used
                  if (rsgratio_ > 0) then
                    tmp_lambdar = (aar*n0rr*gamb1r/(rhof(k)*(qr_spl(i,j,k)*rsgratio_)))**(1./(1.+bbr)) ! lambda rain
                    vtr = ccrz(k)*(gambd1r/gamb1r)/(tmp_lambdar**ddr)  ! terminal velocity rain
                  else
                    vtr = 0
                  end if

                  if ( (1.-rsgratio_)*(1.-sgratio_) > 0 ) then
                    tmp_lambdas = (aas*n0rs*gamb1s/(rhof(k)*(qr_spl(i,j,k)*(1.-rsgratio_)*(1.-sgratio_))))**(1./(1.+bbs)) ! lambda snow
                    vts = ccsz(k)*(gambd1s/gamb1s)/(tmp_lambdas**dds)  ! terminal velocity snow
                  else
                    vts = 0
                  end if

                  if ( (1.-rsgratio_)*sgratio_ > 0 ) then
                    tmp_lambdag = (aag*n0rg*gamb1g/(rhof(k)*(qr_spl(i,j,k)*(1.-rsgratio_)*sgratio_)))**(1./(1.+bbg)) ! lambda graupel
                    vtg = ccgz(k)*(gambd1g/gamb1g)/(tmp_lambdag**ddg)  ! terminal velocity graupel
                  else
                    vtg = 0
                  end if

                  vtf = rsgratio_*vtr+(1.-rsgratio_)*(1.-sgratio_)*vts+(1.-rsgratio_)*sgratio_*vtg  ! mass-weighted terminal velocity
                  vtf = min(wfallmax,vtf)
                  sed_qr(i,j,k) = vtf*qr_spl(i,j,k)*rhobf(k)
                else
                  sed_qr(i,j,k) = 0.
                endif
              enddo
            enddo
          enddo

        else ! alternative loops when there is only rain
          !$acc parallel loop collapse(3) default(present)
          do k = 1, kmax
            do j = 2, j1
              do i = 2, i1
                if (qr_spl(i,j,k) > qrmin) then
                  !*rsgratio(i,j,k) removed from here, since it is 1
                  tmp_lambdar = (aar*n0rr*gamb1r/(rhof(k)*(qr_spl(i,j,k))))**(1./(1.+bbr)) ! lambda rain
                  vtf = ccrz(k)*(gambd1r/gamb1r)/(tmp_lambdar**ddr)  ! terminal velocity rain

                  vtf = min(wfallmax,vtf)
                  sed_qr(i,j,k) = vtf*qr_spl(i,j,k)*rhobf(k)
                else
                  sed_qr(i,j,k) = 0.
                endif
              enddo
            enddo
          enddo
        endif

        !$acc parallel loop collapse(3) default(present)
        do k = 1, kmax
          do j = 2, j1
            do i = 2, i1
              qr_spl(i,j,k) = qr_spl(i,j,k) + &
                              (sed_qr(i,j,k+1) - sed_qr(i,j,k))*dt_spl/(dzh(k+1)*rhobf(k))
            enddo
          enddo
        enddo
        surf_rain = surf_rain + sed_qr(:,:,1)*dt_spl

      enddo ! end time splitting loop
    endif ! if n>1

    ! no thl and qt tendencies build in, implying no heat transfer between precipitation and air

    if (qrsmall > 0.000001*qrsum) then
      write(*,*)'amount of neg. qr thrown away is too high  ',timee, ' sec', qrsmall, qrsum
    end if

    !$acc parallel loop collapse(3) default(present)
    do k = 1, kmax !was k1
      do j = 2, j1
        do i = 2, i1
          qrp(i,j,k) = qrp(i,j,k) + (qr_spl(i,j,k) - qr(i,j,k))/delt
          thlp(i,j,k) = thlp(i,j,k) + thlpmcr(i,j,k)
          qtp(i,j,k) = qtp(i,j,k) + qtpmcr(i,j,k)
          svp(i,j,k,iqr) = svp(i,j,k,iqr) + qrp(i,j,k)
          ! clip the qr tendency so that qr >= 0 next step
          svp(i,j,k,iqr) = max(svp(i,j,k,iqr), -svm(i,j,k,iqr)/delt)
        enddo
      enddo
    enddo

!    if (l_rain) then
!      call simpleicetend !after corrections
!    endif
  end subroutine simpleice2

end module modsimpleice2
