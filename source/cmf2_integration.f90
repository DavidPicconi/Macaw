module cmf2_integration
   use sysparam
   use globals
   use psidef
   use storage, only : S00, nS00max
   use dvr, only : ngpm
   use timingmod
   implicit none
   ! For Lanczos integration
   integer :: nKrylov                                      
   integer, parameter :: nKrylovMax = 20                      ! maximum number of Krylov vectors
   double precision :: dnorm0
   double precision, parameter :: eps_B = 1e-8                ! integration accuracy (might go to the input)
   double precision, dimension(:) :: alphaLan(nKrylovMax)     ! Lanczos energies (after diagonalisation)
   double precision, dimension(:,:) :: ULan(nKrylovMax,nKrylovMax)
   double complex, dimension(:,:), allocatable :: BK, BU   ! Krylov space, product between BK and the eigenvectors of the reduced Hamiltonian

contains

   subroutine SetupLanczos
   implicit none
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
   allocate(BK(nAconftot,nKrylovMax))
   allocate(BU(nAconftot,nKrylovMax))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
   return
   end subroutine SetupLanczos


   subroutine Lanczos(B0, BT, dt, BuildKrylov)
   ! Lanczos integration starting from B0 for a time interval dt
   ! BT is the propagated vector at the end of the step
   ! BuildKrylov is a flag to tell whether the Krylov space has to be constructed or not
   implicit none
   double precision, intent(in) :: dt
   double complex, dimension(:), intent(in)  :: B0(nAconftot)
   double complex, dimension(:), intent(out) :: BT(nAconftot)
   logical, intent(in) :: BuildKrylov
   integer :: j, info
   double precision :: error
   double precision, dimension(:) :: betaLan(nKrylovMax), &   ! upper diagonal of the Lanczos matrix
                                     WORK(2 * nKrylovMax - 2), &
                                     RWORK(2 * nAconftot * nKrylovMax)
   double complex, dimension(:) :: aux(nKrylovMax)
   double precision, external :: dznrm2
   double complex, external :: zdotc
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
   call continue_timer(8)
   !
   if (BuildKrylov) then    ! Krylov space vectors are not already available 
      call zcopy(nAconftot, B0, 1, BK(1,1), 1)    ! BISOGNA NORMALIZZARE B0?
      ! 
      dnorm0 = dznrm2(nAconftot, BK(1,1), 1)
      call zdrscl(nAconftot, dnorm0, BK(1,1), 1)
      !
      call Calc_DHD_B(BK(1,1), BK(1,2))
      alphaLan(1) = dble(zdotc(nAConftot, BK(1,1), 1, BK(1,2), 1))
      call zaxpy(nAconftot, complex(-alphaLan(1),dZero), BK(1,1), 1, BK(1,2), 1)
      betaLan(1) = dznrm2(nAconftot, BK(1,2), 1)
      call zdrscl(nAconftot, betaLan(1), BK(1,2), 1)
      error = betaLan(1) * abs(dt)
      j = 2
      do 
         if (j == nKrylovMax) then
            write(0,*) alphaLan * au2eV
            write(0,*) betaLan * au2eV
            write(0,*) eps_B, error
            write(0,*) 'ERROR: Maximum number of Krylov vectors reached'
            stop
         end if
         call Calc_DHD_B(BK(1,j), BK(1,j + 1))
         alphaLan(j) = real(zdotc(nAConftot, BK(1,j), 1, BK(1,j + 1), 1))
         call zaxpy(nAconftot, complex(-betaLan(j - 1),dZero), BK(1,j - 1), 1, BK(1,j + 1), 1)
         call zaxpy(nAconftot, complex(-alphaLan(j),dZero) , BK(1,j)    , 1, BK(1,j + 1), 1)
         betaLan(j) = dznrm2(nAconftot, BK(1,j + 1), 1)
         call zdrscl(nAconftot, betaLan(j), BK(1,j + 1), 1)
         error = error * betaLan(j) / dble(j) * abs(dt)
         if (error .lt. eps_B) exit
         j = j + 1   ! At the moment we create a vector without using it: it might be worth to estimate alphaLan(j + 1)
      end do
      nKrylov = j
      !
      ! ENERGIES
      !
      call dstev('V',nKrylov,alphaLan,betaLan,ULan,nKrylovMax,WORK,info)   ! alphaLan is rewritten here
      ! BU
      call zlacrm(nAconftot, nKrylov, BK(1,1), nAconftot, ULan, nKrylovMax, BU(1,1), nAconftot, RWORK)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
! The Krylov space has been built and the Lanczos energies have been obtained
! Propagate  B_j(t + dt) = sum_r (BU)_jr * exp(-ii * alpha_r * dt) * ULan_1r*
   !call zgemm('C', 'N', nKrylov, nKrylov, nAconftot, cOne, BK(1,1), nAconftot, BK(1,1), nAconftot, cZero, tmp(1,1), nKrylovMax)
   do j = 1, nKrylov
      aux(j) = ULan(1,j) * exp(der_factor * alphaLan(j) * dt)    ! modify for imaginary time propagation
   end do
   call zgemv('N', nAconftot, nKrylov, cOne, BU(1,1), nAconftot, aux, 1, cZero, BT, 1)
   ! Normalization (IN PRINCIPLE IS NECESSARY ONLY FOR RELAXATION CALCULATIONS)
   dnorm0 = dznrm2(nAconftot, BT, 1)
   call zdrscl(nAconftot, dnorm0, BT, 1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
   call suspend_timer(8)
   return
   end subroutine Lanczos



   subroutine rk_SPF_ck(psi0, DD0, i, typeinfo, htry, hdid, hnext)
   implicit none
   double complex, dimension(:) :: psi0(dimpsi)
   double complex, dimension(:,:,:,:) :: DD0(nS00max,nS00max,nmodGWP,nstates)
   integer, intent(in) :: i
   character(len = 3), intent(in) :: typeinfo  ! 'DVR' or 'GWP'
   double precision, intent(in) :: htry
   double precision, intent(out) :: hdid, hnext
   double precision :: dt, error
   integer :: im, nPar, iel, i1, i2, nS
   double precision, parameter :: safety = 0.9d0, pgrow = -0.2d0, pshrnk = -0.25d0, &
                                  minscal = 0.1d0, maxscal = 5.d0
   double precision, parameter :: errcon = (5.d0 / safety) ** (1.d0 / pgrow)
   ! Dormand-Prince parameters   
   double complex, dimension(dimpsi) :: psitemp, psierr
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates) :: DDtemp
   double complex, dimension(dimpsi) :: dpsi, psik2, psik3, psik4, psik5, psik6
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates) :: dDD, DDk2, DDk3, DDk4, DDk5, DDk6
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
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
   ! Get the mode and number of SPF parameters
   if (typeinfo == 'DVR') then
      im = imodDVR(i)
      nPar = ngpm(i)
   else
      im = imodGWP(i)
      nPar = nmodeDOF(im)
   end if
   !
   dt = htry
   do 
      !
      ! Start Cash-Karp steps
      !
      ! Step 0
      call Calc_Der_SPF(psi0, DD0, dpsi, dDD, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + B21 * dt * dpsi(i1:i2)
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) + B21 * dt * dDD(1:nS,1:nS,i,iel)
      end do
      ! Step 1
      call Calc_Der_SPF(psitemp, DDtemp, psik2, DDk2, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B31 * dpsi(i1:i2) + B32 * psik2(i1:i2))
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B31 * dDD(1:nS,1:nS,i,iel) + B32 * DDk2(1:nS,1:nS,i,iel))
      end do
      ! Step 2
      call Calc_Der_SPF(psitemp, DDtemp, psik3, DDk3, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B41 * dpsi(i1:i2) + B42 * psik2(i1:i2) + B43 * psik3(i1:i2))
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B41 *  dDD(1:nS,1:nS,i,iel) + B42 * DDk2(1:nS,1:nS,i,iel) &
                                          + B43 * DDk3(1:nS,1:nS,i,iel) )
      end do
      ! Step 3
      call Calc_Der_SPF(psitemp, DDtemp, psik4, DDk4, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B51 *  dpsi(i1:i2) + B52 * psik2(i1:i2) + B53 * psik3(i1:i2) &
                                            + B54 * psik4(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B51 *  dDD(1:nS,1:nS,i,iel) + B52 * DDk2(1:nS,1:nS,i,iel) &
                                          + B53 * DDk3(1:nS,1:nS,i,iel) + B54 * DDk4(1:nS,1:nS,i,iel) )
      end do
      ! Step 4
      call Calc_Der_SPF(psitemp, DDtemp, psik5, DDk5, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B61 *  dpsi(i1:i2) + B62 * psik2(i1:i2) + B63 * psik3(i1:i2) &
                                            + B64 * psik4(i1:i2) + B65 * psik5(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B61 *  dDD(1:nS,1:nS,i,iel) + B62 * DDk2(1:nS,1:nS,i,iel) &
                                          + B63 * DDk3(1:nS,1:nS,i,iel) + B64 * DDk4(1:nS,1:nS,i,iel) &
                                          + B65 * DDk5(1:nS,1:nS,i,iel) )
      end do
      ! Step 5
      call Calc_Der_SPF(psitemp, DDtemp, psik6, DDk6, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (C1 *  dpsi(i1:i2) + C3 * psik3(i1:i2) + C4 * psik4(i1:i2) &
                                            + C6 * psik6(i1:i2) )
         psierr(i1:i2)  = psi0(i1:i2) + dt * (D1 *  dpsi(i1:i2) + D3 * psik3(i1:i2) + D4 * psik4(i1:i2) &
                                            + D5 * psik5(i1:i2) + D6 * psik6(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (C1 *  dDD(1:nS,1:nS,i,iel) + C3 * DDk3(1:nS,1:nS,i,iel) &
                                          + C4 * DDk4(1:nS,1:nS,i,iel) + C6 * DDk6(1:nS,1:nS,i,iel))
      end do
      !
      ! Evaluate the error
      call Calc_Err_SPF(psitemp, psierr, i, typeinfo, error)
      error = error / eps_psi
      if (error .lt. dOne) exit
      dt = max(minscal * dt, safety * dt * error ** pshrnk)
   end do
   !
   ! Successful step
   !
   if (error .gt. errcon) then
      hnext = safety * dt * (error ** pgrow)
   else
      hnext = dt * maxscal
   end if
   !hnext = max(hnext,min_stepsize)
   !hnext = min(hnext,writetime)
   hdid = dt
   do iel = 1, merge(nstates, 1, MultiSet(im))       
      nS = nSPF(im,iel)
      i1 = ipsiSPF(im,iel)
      i2 = i1 + nPar * nS - 1
      psi0(i1:i2) = psitemp(i1:i2) 
      if (typeinfo == 'GWP') DD0(1:nS,1:nS,i,iel) = DDtemp(1:nS,1:nS,i,iel)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
   return
   end subroutine rk_SPF_ck



   subroutine rk_SPF_dp(psi0, DD0, i, typeinfo, htry, hdid, hnext)
   implicit none
   double complex, dimension(:) :: psi0(dimpsi)
   double complex, dimension(:,:,:,:) :: DD0(nS00max,nS00max,nmodGWP,nstates)
   integer, intent(in) :: i
   character(len = 3), intent(in) :: typeinfo  ! 'DVR' or 'GWP'
   double precision, intent(in) :: htry
   double precision, intent(out) :: hdid, hnext
   double precision :: dt, error
   integer :: im, nPar, iel, i1, i2, nS
   double precision, parameter :: safety = 0.9d0, pgrow = -0.2d0, pshrnk = -0.25d0, &
                                  minscal = 0.1d0, maxscal = 5.d0
   double precision, parameter :: errcon = (5.d0 / safety) ** (1.d0 / pgrow)
   ! Dormand-Prince parameters   
   double complex, dimension(dimpsi) :: psitemp, psierr
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates) :: DDtemp
   double complex, dimension(dimpsi) :: dpsi, psik2, psik3, psik4, psik5, psik6, psik7
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates) :: dDD, DDk2, DDk3, DDk4, DDk5, DDk6, DDk7
   double precision, parameter ::  &
           A2 = 0.2d0, A3 = 0.3d0, A4 = 0.8d0, A5 = 8.0d0 / 9.0d0, A6 = dOne, A7 = dOne, &   ! The A coefficients are actually not used
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
   ! Get the mode and number of SPF parameters
   if (typeinfo == 'DVR') then
      im = imodDVR(i)
      nPar = ngpm(i)
   else
      im = imodGWP(i)
      nPar = nmodeDOF(im)
   end if
   !
   dt = htry
   do
      !
      ! Start Dormand-Prince steps
      !
      ! Step 0
      call Calc_Der_SPF(psi0, DD0, dpsi, dDD, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + B21 * dt * dpsi(i1:i2)
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) + B21 * dt * dDD(1:nS,1:nS,i,iel)
      end do
      ! Step 1
      call Calc_Der_SPF(psitemp, DDtemp, psik2, DDk2, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B31 * dpsi(i1:i2) + B32 * psik2(i1:i2))
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B31 * dDD(1:nS,1:nS,i,iel) + B32 * DDk2(1:nS,1:nS,i,iel))
      end do
      ! Step 2
      call Calc_Der_SPF(psitemp, DDtemp, psik3, DDk3, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B41 * dpsi(i1:i2) + B42 * psik2(i1:i2) + B43 * psik3(i1:i2))
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B41 *  dDD(1:nS,1:nS,i,iel) + B42 * DDk2(1:nS,1:nS,i,iel) &
                                          + B43 * DDk3(1:nS,1:nS,i,iel) )
      end do
      ! Step 3
      call Calc_Der_SPF(psitemp, DDtemp, psik4, DDk4, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B51 *  dpsi(i1:i2) + B52 * psik2(i1:i2) + B53 * psik3(i1:i2) &
                                            + B54 * psik4(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B51 *  dDD(1:nS,1:nS,i,iel) + B52 * DDk2(1:nS,1:nS,i,iel) &
                                          + B53 * DDk3(1:nS,1:nS,i,iel) + B54 * DDk4(1:nS,1:nS,i,iel) )
      end do
      ! Step 4
      call Calc_Der_SPF(psitemp, DDtemp, psik5, DDk5, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (B61 *  dpsi(i1:i2) + B62 * psik2(i1:i2) + B63 * psik3(i1:i2) &
                                            + B64 * psik4(i1:i2) + B65 * psik5(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (B61 *  dDD(1:nS,1:nS,i,iel) + B62 * DDk2(1:nS,1:nS,i,iel) &
                                          + B63 * DDk3(1:nS,1:nS,i,iel) + B64 * DDk4(1:nS,1:nS,i,iel) &
                                          + B65 * DDk5(1:nS,1:nS,i,iel) )
      end do
      ! Step 5
      call Calc_Der_SPF(psitemp, DDtemp, psik6, DDk6, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psitemp(i1:i2) = psi0(i1:i2) + dt * (C1 *  dpsi(i1:i2) + C3 * psik3(i1:i2) + C4 * psik4(i1:i2) &
                                            + C5 * psik5(i1:i2) + C6 * psik6(i1:i2) )
         if (typeinfo == 'GWP') &
            DDtemp(1:nS,1:nS,i,iel) = DD0(1:nS,1:nS,i,iel) &
                                    + dt * (C1 *  dDD(1:nS,1:nS,i,iel) + C3 * DDk3(1:nS,1:nS,i,iel) &
                                          + C4 * DDk4(1:nS,1:nS,i,iel) + C5 * DDk5(1:nS,1:nS,i,iel) + C6 * DDk6(1:nS,1:nS,i,iel))
      end do
      ! Step 6
      call Calc_Der_SPF(psitemp, DDtemp, psik7, DDk7, i, typeinfo)
      do iel = 1, merge(nstates, 1, MultiSet(im))       
         nS = nSPF(im,iel)
         i1 = ipsiSPF(im,iel)
         i2 = i1 + nPar * nS - 1
         psierr(i1:i2) = psi0(i1:i2) + dt * (D1 *  dpsi(i1:i2) + D3 * psik3(i1:i2) + D4 * psik4(i1:i2) &
                                           + D5 * psik5(i1:i2) + D6 * psik6(i1:i2) + D7 * psik7(i1:i2) )
      end do
      !
      ! Evaluate the error
      call Calc_Err_SPF(psitemp, psierr, i, typeinfo, error)
      error = error / eps_psi
      if (error .lt. dOne) exit
      dt = max(minscal * dt, safety * dt * error ** pshrnk)
   end do
   !
   ! Successful step
   !
   if (error .gt. errcon) then
      hnext = safety * dt * (error ** pgrow)
   else
      hnext = dt * maxscal
   end if
   !hnext = max(hnext,min_stepsize)
   !hnext = min(hnext,writetime)
   hdid = dt
   do iel = 1, merge(nstates, 1, MultiSet(im))       
      nS = nSPF(im,iel)
      i1 = ipsiSPF(im,iel)
      i2 = i1 + nPar * nS - 1
      psi0(i1:i2) = psitemp(i1:i2) 
      if (typeinfo == 'GWP') DD0(1:nS,1:nS,i,iel) = DDtemp(1:nS,1:nS,i,iel)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
   return
   end subroutine rk_SPF_dp


end module cmf2_integration
