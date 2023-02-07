module integration
   use sysparam
   use globals
   use psidef
   use derivatives
   use dvr
   use timingmod
   implicit none
   private
   public :: rkqs, Adams, SetupAdams, bsstep, SetupBS
   public :: h_prev, dA_prev, dpsi_prev
   ! For Adams' method
   double precision :: pow_shrink, pow_grow
   double complex, dimension(:,:), allocatable :: dA_prev, dpsi_prev  ! derivatives at previous times (re-shifted at the end)
   integer, dimension(:), allocatable :: ipiv
   double precision, dimension(:), allocatable :: h_prev   ! previous time steps (re-shifted at the end)
   double precision, dimension(:), allocatable :: BB, XX
   double precision, dimension(:,:), allocatable :: MM, MMcopy
   double complex, dimension(:), allocatable :: XXc
   ! For Bulirsch-Stoer method
   integer, parameter :: kMaxx = 8, iMax = kMaxx + 1
   integer, dimension(iMax), parameter :: nSeq = (/ 2, 4, 6, 8, 10, 12, 14, 16, 18 /)
   integer :: kMax, kOpt
   double precision :: tNew
   double precision, parameter :: safe1 = 0.25d0, safe2 = 0.7d0
   double precision, dimension(:) :: Ak(iMax), tStore(iMax)
   double precision, dimension(:,:) :: alf(kMaxx,kMaxx)
   double complex, dimension(:,:), allocatable :: qA, qpsi
   logical :: first
   !
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! 4th-5th-order Runge-Kutta with variable step size

   subroutine rkqs(A0,psi0,dA,dpsi,tt,htry,hdid,hnext,int_method)
! Adapted from numerical recipes
   implicit none
   double precision, intent(in) :: htry
   double precision, intent(out) :: hdid, hnext
   integer, intent(in) :: int_method
   double precision :: tt, dt, errmax
   double precision :: ovl_11, ovl_22
   double complex :: ovl_12
   double complex, dimension(nAConftot * npackets) :: A0, dA, Atemp, Aerr, dANew
   double complex, dimension(dimpsi) :: psi0, dpsi, psitemp, psierr, dpsiNew
   double precision, parameter :: safety = 0.9d0, pgrow = -0.2d0, pshrnk = -0.25d0, &
                                  minscal = 0.1d0, maxscal = 5.d0
   double precision, parameter :: errcon = (5.d0 / safety) ** (1.d0 / pgrow)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   ! Set the stepsize to the initial trial value
   dt = htry
   do 
      if (int_method == 1) then          ! Cash-Karp parameters
         call rkck(A0,psi0,dA,dpsi,tt,dt,Atemp, psitemp, Aerr, psierr)
      else if (int_method == 2) then     ! Dormand-Prince parameters
         call rkdp(A0,psi0,dA,dpsi,tt,dt,Atemp, psitemp, Aerr, psierr, dANew, dpsiNew)
      end if
      ! Error 
      call DeltaPsi(Atemp,psitemp,Aerr,psierr,errmax)
      errmax = errmax / eps_int
!
      if (errmax .lt. dOne .or. dt .le. min_stepsize) exit
      dt = max(minscal * dt, safety * dt * errmax ** pshrnk)
      dt = max(dt,min_stepsize)
      cycle
      if (dt .lt. min_stepsize) then 
         write(1988,*)
         write(1988,*) ovl_11, ovl_22, ovl_12, errmax
         write(0,'(a,f12.10,a)') 'The stepsize has become too small: ', dt * au2fs, ' fs'
         call StepSizeDrop(A0,psi0,tt)
         stop
      end if
   end do
   if (errmax .gt. errcon) then
      hnext = safety * dt * (errmax ** pgrow)
   else
      hnext = dt * maxscal
   end if
   hnext = max(hnext,min_stepsize)
   hnext = min(hnext,writetime)
   hdid = dt
   call zlacpy('All',nAConftot,npackets,Atemp(1),nAConftot,A0(1),nAConftot)
   call zcopy(dimpsi,psitemp,1,psi0,1)
   if (int_method == 2) then    ! The derivatives for the next step have been already calculated 
      call zlacpy('All',nAConftot,npackets,dANew(1),nAConftot,dA(1),nAConftot)
      call zcopy(dimpsi,dpsiNew,1,dpsi,1)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine rkqs



   subroutine rkck(A0,psi0,dA,dpsi,tt,dt, &
                   Atemp,psitemp,Aerr,psierr)
   ! Runge-Kutta step with Cash-Karp parameters
   implicit none
   double precision :: tt, dt
   double complex, dimension(:), intent(in) :: A0(nAconftot * npackets), &
                                               dA(nAconftot * npackets), &
                                               psi0(dimpsi), dpsi(dimpsi)
   double complex, dimension(:), intent(out) :: Atemp(nAconftot * npackets), &
                                                Aerr(nAconftot * npackets), &
                                                psitemp(dimpsi), psierr(dimpsi)
   double complex, dimension(nAconftot * npackets) :: Ak2, Ak3, Ak4, Ak5, Ak6
   double complex, dimension(dimpsi) :: psik2, psik3, psik4, psik5, psik6
   double precision, parameter ::  &
      A2 = 0.2d0, A3 = 0.3d0, A4 = 0.6d0, A5 = 1.0d0, A6 = 0.875d0, &
      B21 = 0.2d0, &
      B31 = 0.075d0, B32 = 0.225d0, &
      B41 = 0.3d0  , B42 = -0.9d0 , B43 = 1.2d0, &
      B51 = -11.0d0 / 54.0d0, B52 = 2.5d0, B53 = -70.0d0 / 27.0d0, B54 = 35.0d0 / 27.0d0, &
      B61 = 1631.0d0 / 55296.0d0, B62 = 175.d0 / 512.d0, B63 = 575.0d0 / 13824.0d0, &
      B64 = 44275.0d0 / 110592.0d0, B65 = 253.0d0 / 4096.0d0, &
      C1 = 37.0d0 / 378.0d0, C3 = 250.0d0 / 621.0d0, C4 = 125.0d0 / 594.0d0, C6 = 512.0d0 / 1771.0d0, &
      D1 = 2825.0d0 / 27648.0d0, D3 = 18575.0d0 / 48384.0d0, D4 = 13525.0d0 / 55296.0d0, &
      D5 = 277.0d0 / 14336.0d0 , D6 = 0.25d0
      !DC1 = C1 - 2825.0d0 / 27648.0d0, DC3 = C3 - 18575.0d0 / 48384.0d0, &
      !DC4 = C4 - 13525.0d0 / 55296.0d0, DC5 = -277.0d0 / 14336.0d0, DC6 = C6 - 0.25d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   ! First step
   Atemp   = A0   + B21 * dt * dA
   psitemp = psi0 + B21 * dt * dpsi
   ! Second step
   call Calc_Derivs(Atemp,psitemp, Ak2,psik2, tt + A2 * dt)
   Atemp   = A0   + dt * (B31 * dA   + B32 * Ak2) 
   psitemp = psi0 + dt * (B31 * dpsi + B32 * psik2) 
   ! Third step
   call Calc_Derivs(Atemp,psitemp, Ak3,psik3, tt + A3 * dt)
   Atemp   = A0   + dt * (B41 * dA   + B42 * Ak2   + B43 * Ak3)
   psitemp = psi0 + dt * (B41 * dpsi + B42 * psik2 + B43 * psik3)
   ! Fourth step
   call Calc_Derivs(Atemp,psitemp, Ak4,psik4, tt + A4 * dt)
   Atemp   = A0   + dt * (B51 * dA   + B52 * Ak2   + B53 * Ak3   + B54 * Ak4)
   psitemp = psi0 + dt * (B51 * dpsi + B52 * psik2 + B53 * psik3 + B54 * psik4)
   ! Fifth step
   call Calc_Derivs(Atemp,psitemp, Ak5,psik5, tt + A5 * dt)
   Atemp   = A0   + dt * (B61 * dA   + B62 * Ak2   + B63 * Ak3   + B64 * Ak4   + B65 * Ak5)
   psitemp = psi0 + dt * (B61 * dpsi + B62 * psik2 + B63 * psik3 + B64 * psik4 + B65 * psik5)
   ! Sixth step
   call Calc_Derivs(Atemp,psitemp, Ak6,psik6, tt + A6 * dt)
   ! Accumulate increments with proper weights
   Atemp   = A0   + dt * (C1 * dA   + C3 * Ak3   + C4 * Ak4   + C6 * Ak6)      ! 5th order estimate
   psitemp = psi0 + dt * (C1 * dpsi + C3 * psik3 + C4 * psik4 + C6 * psik6)
   Aerr    = A0   + dt * (D1 * dA   + D3 * Ak3   + D4 * Ak4   + D5 * Ak5   + D6 * Ak6)
   psierr  = psi0 + dt * (D1 * dpsi + D3 * psik3 + D4 * psik4 + D5 * psik5 + D6 * psik6)
   ! Estimate error as difference between fourth and fifth order methods
   !Aerr = Atemp - dt * (DC1 * dA + DC3 * Ak3 &
   !                   + DC4 * Ak4 + DC5 * Ak5 + DC6 * Ak6)
   !psierr = psitemp - dt * (DC1 * dpsi + DC3 * psik3 &
   !                       + DC4 * psik4 + DC5 * psik5 + DC6 * psik6)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine rkck

   

   subroutine rkdp(A0,psi0,dA,dpsi,tt,dt, &
                   Atemp,psitemp,Aerr,psierr, Ak7, psik7)
   ! Runge-Kutta step with Dormand-Prince parameters
   implicit none
   double precision :: tt, dt
   double complex, dimension(:), intent(in) :: A0(nAconftot * npackets), &
                                               dA(nAconftot * npackets), &
                                               psi0(dimpsi), dpsi(dimpsi)
   double complex, dimension(:), intent(out) :: Atemp(nAconftot * npackets), &
                                                Aerr(nAconftot * npackets), &
                                                psitemp(dimpsi), psierr(dimpsi)
   double complex, dimension(nAconftot * npackets) :: Ak2, Ak3, Ak4, Ak5, Ak6, Ak7
   double complex, dimension(dimpsi) :: psik2, psik3, psik4, psik5, psik6, psik7
   double precision, parameter ::  &
      A2 = 0.2d0, A3 = 0.3d0, A4 = 0.8d0, A5 = 8.0d0 / 9.0d0, A6 = dOne, A7 = dOne, &
      B21 = 0.2d0, &
      B31 = 0.075d0, B32 = 0.225d0, &
      B41 = 44.0d0 / 45.0d0  , B42 = -56.0d0 / 15.0d0 , B43 = 32.0d0 / 9.0d0, &
      B51 = 19372.0d0 / 6561.0d0, B52 = -25360.0d0 / 2187.0d0, &
      B53 = 64448.0d0 / 6561.0d0, B54 = -212.0d0 / 729.0d0, &
      B61 = 9017.0d0 / 3168.0d0, B62 = -355.0d0 / 33.0d0, B63 = 46732.0d0 / 5247.0d0, &
      B64 = 49.0d0 / 176.0d0, B65 = -5103.0d0 / 18656.0d0, &
      C1 = 35.0d0 / 384.0d0, C3 = 500.0d0 / 1113.0d0, C4 = 125.0d0 / 192.0d0, &
      C5 = -2187.0d0 / 6784.0d0, C6 = 11.0d0 / 84.0d0, &
      D1 = 5179.0d0 / 57600.0d0, D3 = 7571.0d0 / 16695.0d0, D4 = 393.0d0 / 640.0d0, &
      D5 = -92097.0d0 / 339200.0d0, D6 = 187.0d0 / 2100.0d0, D7 = 0.025d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   ! First step
   Atemp   = A0   + B21 * dt * dA
   psitemp = psi0 + B21 * dt * dpsi
   ! Second step
   call Calc_Derivs(Atemp,psitemp, Ak2,psik2, tt + A2 * dt)
   Atemp   = A0   + dt * (B31 * dA   + B32 * Ak2) 
   psitemp = psi0 + dt * (B31 * dpsi + B32 * psik2) 
   ! Third step
   call Calc_Derivs(Atemp,psitemp, Ak3,psik3, tt + A3 * dt)
   Atemp   = A0   + dt * (B41 * dA   + B42 * Ak2   + B43 * Ak3)
   psitemp = psi0 + dt * (B41 * dpsi + B42 * psik2 + B43 * psik3)
   ! Fourth step
   call Calc_Derivs(Atemp,psitemp, Ak4,psik4, tt + A4 * dt)
   Atemp   = A0   + dt * (B51 * dA   + B52 * Ak2   + B53 * Ak3   + B54 * Ak4)
   psitemp = psi0 + dt * (B51 * dpsi + B52 * psik2 + B53 * psik3 + B54 * psik4)
   ! Fifth step
   call Calc_Derivs(Atemp,psitemp, Ak5,psik5, tt + A5 * dt)
   Atemp   = A0   + dt * (B61 * dA   + B62 * Ak2   + B63 * Ak3   + B64 * Ak4   + B65 * Ak5)
   psitemp = psi0 + dt * (B61 * dpsi + B62 * psik2 + B63 * psik3 + B64 * psik4 + B65 * psik5)
   ! Sixth step
   call Calc_Derivs(Atemp,psitemp, Ak6,psik6, tt + A6 * dt)
   ! Accumulate increments with proper weights
   Atemp   = A0   + dt * (C1 * dA   + C3 * Ak3   + C4 * Ak4   + C5 * Ak5   + C6 * Ak6)      ! 5th order estimate
   psitemp = psi0 + dt * (C1 * dpsi + C3 * psik3 + C4 * psik4 + C5 * psik5 + C6 * psik6)
   ! Seventh step
   call Calc_Derivs(Atemp,psitemp, Ak7,psik7, tt + A7 * dt)
   Aerr    = A0   + dt * (D1 * dA   + D3 * Ak3   + D4 * Ak4   + D5 * Ak5   + D6 * Ak6   + D7 * Ak7)
   psierr  = psi0 + dt * (D1 * dpsi + D3 * psik3 + D4 * psik4 + D5 * psik5 + D6 * psik6 + D7 * psik7)
   ! Estimate error as difference between fourth and fifth order methods
   !Aerr = Atemp - dt * (DC1 * dA + DC3 * Ak3 &
   !                   + DC4 * Ak4 + DC5 * Ak5 + DC6 * Ak6)
   !psierr = psitemp - dt * (DC1 * dpsi + DC3 * psik3 &
   !                       + DC4 * psik4 + DC5 * psik5 + DC6 * psik6)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine rkdp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Fixed-order Adams PECE formulas with variable step size

   subroutine SetupAdams()
   implicit none
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   allocate(h_prev(int_order), dA_prev(nAconftot * npackets,int_order), dpsi_prev(dimpsi,int_order))
   allocate(MM(0:int_order + 1,-1:int_order), MMcopy(0:int_order + 1,-1:int_order))
   allocate(BB(0:int_order + 1),XX(-1:int_order),XXc(-1:int_order))
   MM(0,:) = dOne
   MM(1:int_order + 1,0) = dZero
   allocate(ipiv(0:int_order + 1))
   !pow_shrink = - dOne / dble(int_order + 2)
   pow_shrink = - dOne / dble(int_order + 1)
   pow_grow = - dOne / dble(int_order + 3)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   end subroutine SetupAdams



   subroutine Adams(A0,psi0,dA,dpsi,tt,htry,hdid,hnext)
! Integration using Adams' PECE formulas with variable step size
! (Theory from the book of Hairer, Norsett, Wanner)
! - A0 and psi0 contain the wavefunction at time t_n  (A_n and psi_n)
! - dA and dpsi contain the derivative at time t_n    (dA_n and dpsi_n)
! - h_prev(j) contains the time intervals t_(n - j + 1) - t_(n - j)
! - dA_prev(:,j) contains dA_(n - j)
! - dpsi_prev(:,j) contains dpsi_(n - j)
   implicit none
   double complex, dimension(nAconftot * npackets) :: A0, dA, Atemp, &
                                                      ANew, dANew, dANew1, ANew1
   double complex, dimension(dimpsi) :: psi0, dpsi, psitemp, psiNew, dpsiNew, dpsiNew1, psiNew1
   double precision :: tt, htry, hdid, hnext, errmax, dt
   double precision :: hsum, dummy
   integer :: i, j, info, k
   double precision, parameter :: minscal = 0.3d0, maxscal = 4.d0, safety = 0.95d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Columns >= 1 of the matrix M(j,i) = [x_(n-i) - x_n] ** j
   hsum = dZero
   do i = 1, int_order
      hsum = hsum - h_prev(i)
      dummy = hsum
      do j = 1, int_order + 1
         MM(j,i) = dummy
         dummy = dummy * hsum
      end do
   end do
   dt = htry
   do   ! This cycle is repeated as long as the step size dt needs to be decreased
      ! '-1' column of the matrix M(j,i) = [x_(n-i) - x_n] ** j
      BB(0) = dt
      dummy = dt
      do j = 1, int_order + 1
           MM(j,-1) = dummy
           dummy = dummy * dt
           BB(j) = dummy / dble(j + 1)
      enddo 
      ! Calculate B(j) = [ x_(n+1) - x_n ] ** (j + 1) / (j + 1)
      ! ! ! ! ! ! !
      ! (1) Predictor (order k)
      ! ! ! ! ! ! !
      call dcopy(int_order,BB(0),1,XX(0),1)
      call dlacpy('',int_order,int_order,MM(0,0),int_order + 2,MMcopy(0,0),int_order + 2)
      call dgesv(int_order, 1, MMcopy(0,0), int_order + 2, &
                 ipiv, XX(0), int_order + 2,info)
      call zlacp2('All',int_order,1,XX(0),int_order + 2,XXc(0),int_order + 2)
         ! Atemp = A0 + dA * XXc(0) + matmul(dA_prev(:,1:int_order - 1),XXc(1:int_order - 1))
      call zcopy(nAConfTot * npackets,A0,1,Atemp,1)
      call zaxpy(nAConfTot * npackets,XXc(0),dA,1,Atemp,1)
      call zgemv('N',nAConfTot * npackets,int_order - 1,cOne,dA_prev(1,1),nAConfTot * npackets,XXc(1),1,cOne,Atemp(1),1)
      call zcopy(dimpsi,psi0,1,psitemp,1)
      call zaxpy(dimpsi,XXc(0),dpsi,1,psitemp,1)
      call zgemv('N',dimpsi,int_order - 1,cOne,dpsi_prev(1,1),dimpsi,XXc(1),1,cOne,psitemp(1),1)
      ! Estimate the derivative at t + dt
      call Calc_Derivs(Atemp,psitemp,dANew,dpsiNew, tt + dt)
      ! ! ! ! ! ! ! 
      ! (2) Corrector (order k + 1)
      ! ! ! ! ! ! !
      call dcopy(int_order + 1,BB(0),1,XX(-1),1)
      call dlacpy('',int_order + 1,int_order + 1,MM(0,-1),int_order + 2,MMcopy(0,-1),int_order + 2)
      call dgesv(int_order + 1, 1, MMcopy(0,-1), int_order + 2, &
                 ipiv, XX(-1), int_order + 2, info)
      call zlacp2('All',int_order + 1,1,XX(-1),int_order + 2,XXc(-1),int_order + 2)
      call zcopy(nAConfTot * npackets,A0,1,Atemp,1)
      call zaxpy(nAConfTot * npackets,XXc(0),dA,1,Atemp,1)
      call zgemv('N',nAConfTot * npackets,int_order - 1,cOne,dA_prev(1,1),nAConfTot * npackets,XXc(1),1,cOne,Atemp(1),1)
      call zcopy(dimpsi,psi0,1,psitemp,1)
      call zaxpy(dimpsi,XXc(0),dpsi,1,psitemp,1)
      call zgemv('N',dimpsi,int_order - 1,cOne,dpsi_prev(1,1),dimpsi,XXc(1),1,cOne,psitemp(1),1)
      ANew = dANew * XXc(-1) + Atemp
      psiNew = dpsiNew * XXc(-1) + psitemp
         ! Other EC iterations
      do k = 1, int_par1
         call Calc_Derivs(ANew,psiNew,dANew1,dpsiNew1, tt + dt)
         ANew = dANew1 * XXc(-1) + Atemp
         psiNew = dpsiNew1 * XXc(-1) + psitemp
      end do
      ! ! ! ! ! ! !
      ! (3) Corrector (order k + 2)
      ! ! ! ! ! ! !
      call dcopy(int_order + 2,BB(0),1,XX(-1),1)
      call dlacpy('',int_order + 2,int_order + 2,MM(0,-1),int_order + 2,MMcopy(0,-1),int_order + 2)
      call dgesv(int_order + 2, 1, MMcopy(0,-1), int_order + 2, &
                 ipiv, XX(-1), int_order + 2, info)
      call zlacp2('All',int_order + 2,1,XX(-1),int_order + 2,XXc(-1),int_order + 2)
      call zcopy(nAConfTot * npackets,A0,1,Atemp,1)
      call zaxpy(nAConfTot * npackets,XXc(0),dA,1,Atemp,1)
      call zgemv('N',nAConfTot * npackets,int_order,cOne,dA_prev(1,1),nAConfTot * npackets,XXc(1),1,cOne,Atemp(1),1)
      call zcopy(dimpsi,psi0,1,psitemp,1)
      call zaxpy(dimpsi,XXc(0),dpsi,1,psitemp,1)
      call zgemv('N',dimpsi,int_order,cOne,dpsi_prev(1,1),dimpsi,XXc(1),1,cOne,psitemp(1),1)
      ANew1 = dANew * XXc(-1) + Atemp
      psiNew1 = dpsiNew * XXc(-1) + psitemp
         ! Other EC iterations
      do k = 1, int_par1
         call Calc_Derivs(ANew1,psiNew1,dANew1,dpsiNew1, tt + dt)
         ANew1 = dANew1 * XXc(-1) + Atemp
         psiNew1 = dpsiNew1 * XXc(-1) + psitemp
      end do
      !!!!
      ! Error control
      !!!!
      call DeltaPsi(ANew,psiNew,ANew1,psiNew1,errmax)
      errmax = errmax / eps_int
      if (errmax .lt. dOne) exit
      dt = max(minscal * dt, safety * dt * errmax ** pow_shrink)
      if (dt .lt. min_stepsize) then
         write(0,'(a,f12.10,a)') 'The stepsize has become too small: ', dt * au2fs, ' fs'
         call StepSizeDrop(ANew,psiNew,tt)
         stop
      end if
   end do
   ! Increase step size
   hnext = min(maxscal * dt, dt * errmax ** pow_shrink)
   hnext = min(hnext,writetime)
   hdid = dt
   ! Shift the data of previous points
   dA_prev   = cshift(dA_prev, shift = -1, dim = 2)
   dpsi_prev = cshift(dpsi_prev, shift = -1, dim = 2)
   h_prev    = cshift(h_prev, shift = -1, dim = 1)
   call zlacpy('All',nAConftot,npackets,dA(1),nAConftot,dA_prev(1,1),nAConftot)
   !call zcopy(nAConftot * npackets,dA,1,dA_prev(1,1),1)
   call zcopy(dimpsi,dpsi,1,dpsi_prev(1,1),1)
   h_prev(1) = hdid
   call zlacpy('All',nAConftot,npackets,ANew1(1),nAConftot,A0(1),nAConftot)
   !call zcopy(nAConftot * npackets,ANew1,1,A0,1)
   call zcopy(dimpsi,psiNew1,1,psi0,1)
   ! Update the derivatives 
   call Calc_Derivs(A0,psi0,dA,dpsi,tt + hdid)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine Adams



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SetupBS
   integer :: k, iq, kOpt
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   tNew = -1.d16  ! This should be set to an "impossible" value
   ! Compute the work coefficients Ak  (Shouldn't they be defined as integers???)
   Ak(1) = nSeq(1) + dOne
   do k = 1, kMaxx
      Ak(k + 1) = Ak(k) + nSeq(k + 1)
   end do
   ! Compute the tableau alpha(k,q)
   do iq = 2, kMaxx
      do k = 1, iq - 1 
         alf(k,iq) = (safe1 * eps_int) ** ((Ak(k + 1) - Ak(iq + 1)) &
                   / ((Ak(iq + 1) - Ak(1) + dOne) * dble(2 * k + 1)))
      end do
   end do
   ! Determine the optimal row number for convergence
   do kOpt = 2, kMaxx - 1
      if (Ak(kOpt + 1) .gt. Ak(kOpt) * alf(kOpt - 1,kOpt)) exit
   end do
   kMax = kOpt
   ! Allocate qA and qpsi (used for polynomial extrapolation)
   allocate(qA(nAconftot * npackets,iMax), qpsi(dimpsi,iMax))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine SetupBS


   
   subroutine bsstep(A0,psi0,dA,dpsi,tt,htry,hdid,hnext)
   double precision, intent(in) :: htry
   double precision, intent(out) :: hdid, hnext
   double precision :: tt, dt, tEst, errmax
   double complex, dimension(nAconftot * npackets) :: A0, dA, ASave, ASeq, A1
   double complex, dimension(dimpsi) :: psi0, dpsi, psiSave, psiSeq, psi1
   integer :: k, kk
   double precision :: red, work, workMin, fact, scal
   double precision, parameter :: RedMax = 1.d-5, RedMin = 0.7d0, scalMax = 0.1d0
   double precision, dimension(:) :: error(kMaxx)
   logical :: reduct, converged
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   dt = htry
   ! Save the starting values
   ASave = A0
   psiSave = psi0
   ! Check if the step size changed. In this case the order window is re-established
   if (dt .ne. hnext .or. tt .ne. tNew) then
      first = .true.
      kOpt = kMax
   end if
   reduct = .false.
   converged = .false.
   ! 
   do    ! Repeat until convergence is reached for a certain step size
      do k = 1, kMax
         tNew = tt + dt
         ! Midpoint integration
         call mmid(ASave,psiSave,dA,dpsi,tt,dt,nSeq(k),ASeq,psiSeq)
         tEst = (dt / nSeq(k)) ** 2
         ! Extrapolation
         call pzextr(k, tEst, ASeq, psiSeq, A0, psi0, A1, psi1)
         if (k .ne. 1) then
            !
            ! Error evaluation
            !
            call DeltaPsi(A0,psi0,A1,psi1,errmax)
            errmax = errmax / eps_int
            error(k - 1) = (errmax / safe1) ** (dOne / dble(2 * k - 1))
         end if
         if (k .ne. 1 .and. (k .ge. kOpt - 1 .or. first)) then
            if (errmax .lt. dOne) then  ! Converged
               converged = .true.
               exit  
            end if
            ! Stepsize reduction
            if (k .eq. kmax .or. k .eq. kOpt + 1) then
               red = safe2 / error(k - 1)
               exit
            else if (k .eq. kOpt) then
               if (alf(kOpt - 1,kOpt) .lt. error(k - 1)) then
                  red = dOne / error(k - 1)
                  exit
               endif
            else if (kOpt .eq. kMax) then
               if (alf(k - 1,kMax - 1) .lt. error(k - 1)) then
                  red = alf(k - 1,kMax - 1) * safe2 / error(k - 1)
                  exit
               end if
            else if (alf(k - 1,kOpt) .lt. error(k - 1)) then
               red = alf(k - 1,kOpt - 1) / error(k - 1)
               exit
            end if
         end if
      end do
      if (converged) exit   ! Step accepted
      red = max(min(red, RedMin), RedMax)
      dt = dt * red
      reduct = .true.
   end do
   ! Successful step taken
   !tt = tNew
   hdid = dt
   first = .false.
   ! Compute optimal row for convergence and corresponding step size
   workMin = 1.d16
   do kk = 1, k - 1
      fact = max(error(kk), scalMax)
      work = fact * Ak(kk + 1)
      if (work .lt. workMin) then
         scal = fact
         workMin = work
         kOpt = kk + 1
      endif
   end do
   hnext = dt / scal
   ! Check for possible order increase but not if the step size was just reduced
   if (kOpt .ge. k .and. kOpt .ne. kMax .and. .not. reduct) then
      fact = max(scal / alf(kOpt - 1,kOpt), scalMax)
      if (Ak(kOpt + 1) * fact .le. workMin) then
         hnext = dt / fact
         kOpt = kOpt + 1
      end if
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine bsstep


   subroutine mmid(A0,psi0,dA,dpsi,tt,htot,nStep,AOut,psiOut)
   integer, intent(in) :: nStep
   double precision, intent(in) :: tt, htot
   double complex, dimension(:), intent(in) :: A0(nAconftot * npackets), dA(nAconftot * npackets), &
                                               psi0(dimpsi), dpsi(dimpsi)
   double complex, dimension(:), intent(out) :: AOut(nAconftot * npackets), psiOut(dimpsi)
   double complex, dimension(:) :: Am(nAconftot * npackets), An(nAconftot * npackets), &
                                   psim(dimpsi) , psin(dimpsi) , &
                                   Aswap(nAconftot * npackets), psiswap(dimpsi)
   integer :: n
   double precision :: dt, tx, dt2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   dt = htot / nStep
   Am   = A0
   psim = psi0
   An   = A0 + dt * dA
   psin = psi0 + dt * dpsi
   tx = tt + dt
   call Calc_Derivs(An,psin,AOut,psiOut,tx)   ! AOut and psiOut are used for temporary storage of derivatives
   dt2 = dt * 2
   do n = 2, nStep
      Aswap   = Am   + dt2 * AOut
      psiswap = psim + dt2 * psiOut
      Am      = An
      psim    = psin
      An      = Aswap
      psin    = psiswap
      tx = tx + dt
      call Calc_Derivs(An,psin,AOut,psiOut,tx)
   end do
   AOut   = cOneHalf * (Am + An + dt * AOut)
   psiOut = cOneHalf * (psim + psin + dt * psiOut)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine mmid


   subroutine pzextr(iEst, tEst, AEst, psiEst, A0, psi0, A1, psi1)
   integer, intent(in) :: iEst
   double precision, intent(in) :: tEst
   double complex, dimension(nAconftot * npackets) :: AEst, A0, A1, Atemp
   double complex, dimension(dimpsi) :: psiEst, psi0, psi1, psitemp
   integer :: k1, i
   double precision :: delta, f1, f2
   double complex :: cdelta, cdummy
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   tStore(iEst) = tEst
   A0 = AEst
   psi0 = PsiEst
   A1 = AEst
   psi1 = psiEst
   if (iEst .eq. 1) then
      qA(:,1) = AEst
      qpsi(:,1) = psiEst
      return
   end if
   Atemp = AEst
   psitemp = psiEst
   do k1 = 1, iEst - 1
      delta = dOne / (tStore(iEst - k1) - tEst)
      f1 = tEst * delta
      f2 = tStore(iEst - k1) * delta
      do i = 1, nAconftot * npackets
         cdummy   = qA(i,k1)
         qA(i,k1) = A1(i)
         cdelta   = Atemp(i) - cdummy
         A1(i)    = f1 * cdelta
         Atemp(i) = f2 * cdelta
         A0(i)    = A0(i) + A1(i)
      end do
      do i = 1, dimpsi
         cdummy     = qpsi(i,k1)
         qpsi(i,k1) = psi1(i)
         cdelta     = psitemp(i) - cdummy
         psi1(i)    = f1 * cdelta
         psitemp(i) = f2 * cdelta
         psi0(i)    = psi0(i) + psi1(i)
      end do
   end do
   qA(:,iEst) = A1
   qpsi(:,iEst) = psi1
   ! At this stage, A1 and psi1 are the errors
   A1 = A0 + A1
   psi1 = psi0 + psi1
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine pzextr

end module integration



