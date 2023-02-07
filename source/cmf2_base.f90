module cmf2_base
   use sysparam
   use globals
   use psidef
   use storage, only : S00, nS00max, energyMF
   use hamilton, only : Dh1modeD, nOpMax, Calc_h1mode, Calc_DHD, nHamTD
   use cmf2_integration
   implicit none
   private
   public :: SetupCMF2, cmf2

contains


   subroutine SetupCMF2
!  Allocate the matrices D and D^H * H * D
!  needed with the constant mean field approach
!
!  Initialise the D matrix (it is assumed that the overlap matrix is already calculated)
!
!  Initialise the integrators
   implicit none
   integer :: i, im, iel, nS
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
! D matrix and D^H * H * D
   allocate(DD(nS00max,nS00max,nmodGWP,nstates))
   allocate(Dh1modeD(nS00max,nS00max,nOpMax,nmodGWP))
! Initialise D to S^(-1/2)
   do i = 1, nmodGWP
      im = imodGWP(i)
      do iel = 1, merge(nstates, 1, MultiSet(im))
         nS = nSPF(im,iel)
         !call InverseSquareRootU(nS, S00(1,1,i,iel,iel), nS00max, DD(1,1,i,iel), nS00max) 
         call InverseSquareRoot(nS, S00(1,1,i,iel,iel), nS00max, DD(1,1,i,iel), nS00max, dZero) 
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
! INITIALISE THE INTEGRATORS
   ! Lanczos
   if (int_type_A == 5) call SetupLanczos
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   return
   end subroutine SetupCMF2



   subroutine cmf2(A0,psi0,DD0,tt,htry,hdid,hnext)
   implicit none
   double precision, intent(in) :: htry
   double precision, intent(out) :: hdid, hnext
   integer :: i, im, iel, nS
   double precision, intent(in) :: tt
   double complex, dimension(:) :: A0(nAconftot), AHalf(nAconftot), &
                                   B0(nAconftot), BHalf(nAconftot), B02(nAconftot), &
                                   psi0(dimpsi), psiHalf(dimpsi), psiHalf2(dimpsi)
   double complex, dimension(:,:,:,:) :: DD0(nS00max,nS00max,nmodGWP,nstates), &
                                      DDHalf(nS00max,nS00max,nmodGWP,nstates), &
                                     DDHalf2(nS00max,nS00max,nmodGWP,nstates)
   double precision :: err_psi, err_A, error
   double precision, external :: dznrm2
   double precision, parameter :: safe = 0.95d0, maxscal = 1.5d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!  
   ! Recalculate D00 at the beginning of the CMF2 step
   ! Initialise D to S^(-1/2)
   do i = 1, nmodGWP
      call Calc_S00(psi0,i)
      im = imodGWP(i)
      do iel = 1, merge(nstates, 1, MultiSet(im))
         nS = nSPF(im,iel)
         call InverseSquareRootU(nS, S00(1,1,i,iel,iel), nS00max, DD0(1,1,i,iel), nS00max) 
         !call InverseSquareRoot(nS, S00(1,1,i,iel,iel), nS00max, DD0(1,1,i,iel), nS00max) 
      end do
   end do
   !
   hdid = htry
   do
      ! Calculate the operators and the mean fields
      do i = 1, nmodGWP
         call Calc_S00(psi0,i)
      end do
      do im = 1, nmodes
         call Calc_rho1mode(im, A0)
         call Calc_h1mode(psi0,im)
      end do
      call Calc_MeanField(A0,psi0,tt)
      call Calc_rhom1MF
      call Calc_DHD(DD0)
   !
   ! STEP 1: Propagate the coefficients from 0 to tau/2
   !
      ! Calculate the initial vector B = D^H * S * A
      call A_to_B(A0,B0,S00,DD0)
      call UpdateTDcoeff(tt)
      ! Propagate
      select case(int_type_A)
         case(5)  ! Lanczos
            call Lanczos(B0, BHalf, dOneHalf * hdid, .true.)
         case default
            write(0,*) 'Error: integrator not yet implemented for the A coefficients'
            stop
      end select
   !
   ! STEP 2: Propagate the SPFs from 0 to tau/2 (different propagations for different modes)
   !
      ! Copy the initial values
      psiHalf = psi0
      DDHalf = DD0
      ! DVR modes
      do i = 1, nmodDVR
         call propagate_SPF(dOneHalf * hdid, i, 'DVR', psiHalf, DDHalf)    ! integrate from t0 to t0 + 0.5*htry
      end do
      ! GWP modes
      do i = 1, nmodGWP
         call propagate_SPF(dOneHalf * hdid, i, 'GWP', psiHalf, DDHalf) 
      end do
   !
   ! STEP 3: Calculate the mean fields at the half step and re-propagate the SPFs from 0 to tau/2
   !
      ! Transform the coefficient vector
      call B_to_A(AHalf,BHalf,DDHalf)   ! more correct than call B_to_A(AHalf,BHalf,DD0) , right?
      ! Calculate the operators and the mean fields
      do i = 1, nmodGWP
         call Calc_S00(psiHalf,i)
      end do
      do im = 1, nmodes
         call Calc_rho1mode(im, AHalf)
         call Calc_h1mode(psiHalf,im)
      end do
      call Calc_MeanField(AHalf,psiHalf,tt + dOneHalf * hdid)
      call Calc_rhom1MF
      ! Copy the initial values
      psiHalf2 = psi0
      DDHalf2 = DD0
      ! DVR modes
      do i = 1, nmodDVR
         call propagate_SPF(dOneHalf * hdid, i, 'DVR', psiHalf2, DDHalf2)    ! integrate from t0 to t0 + 0.5*htry
      end do
      ! GWP modes
      do i = 1, nmodGWP
         call propagate_SPF(dOneHalf * hdid, i, 'GWP', psiHalf2, DDHalf2) 
      end do
      ! Evaluate the CMF error for the SPFs
      call DeltaPsi(AHalf,psiHalf,AHalf,psiHalf2,err_psi)   
   !
   ! STEP 4: Propagate the SPFs till the end of the step using the mean fields at the half step
   !
      ! DVR modes
      do i = 1, nmodDVR
         call propagate_SPF(dOneHalf * hdid, i, 'DVR', psiHalf2, DDHalf2)    ! integrate from t0 to t0 + 0.5*htry
      end do
      ! GWP modes
      do i = 1, nmodGWP
         call propagate_SPF(dOneHalf * hdid, i, 'GWP', psiHalf2, DDHalf2) 
      end do
   !
   ! STEP 5: Calculate the Hamiltonian matrix at the final step and back-propagate the coefficients from h/2 to 0
   !
      do i = 1, nmodGWP
         call Calc_S00(psiHalf2,i)
      end do
      do im = 1, nmodes
         call Calc_h1mode(psiHalf2,im)
      end do
      call Calc_DHD(DDHalf2)
      call UpdateTDcoeff(tt + hdid)
      ! Propagate
      select case(int_type_A)
         case(5)  ! Lanczos
            call Lanczos(BHalf, B02, - dOneHalf * hdid, .true.)
         case default
            write(0,*) 'Error: integrator not yet implemented for the A coefficients'
            stop
      end select
      ! Evaluate the error in the coefficients
      call zaxpy(nAConftot, -cOne, B0,1, B02,1)
      err_A = dznrm2(nAconftot, B02, 1) * dOneHalf
      !
      ! Check whether the step has to be repeated
      !
      error = max(err_A, err_psi) / eps_int
      if (error .gt. dOne) then
      !if (error .gt. dOne * 2) then
         hdid = hdid / sqrt(error) * safe
         cycle
      else
         ! Terminate the propagation of the coefficients
         select case(int_type_A)
            case(5)  ! Lanczos
               call Lanczos(BHalf, B0, dOneHalf * hdid, .false.)
            case default
               write(0,*) 'Error: integrator not yet implemented for the A coefficients'
               stop
         end select 
         ! Update
         call B_to_A(A0,B0,DDHalf2)
         psi0 = psiHalf2
         DD0  = DDHalf2
         ! Set the next step size and exit
         hnext = min(hdid / sqrt(error), hdid * maxscal)
         hnext = min(hnext, writetime)
         exit
      end if
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   return
   end subroutine cmf2



   subroutine propagate_SPF(dt, i, typeinfo, psiHalf, DDHalf)
   ! Propagate the SPFs of a single mode overl the interval dt
   ! psiHalf and DDHalf contain the initial values of the SPFs and the D matrices
   ! and are updated at the end of the step
   implicit none
   double precision, intent(in) :: dt
   integer, intent(in) :: i
   character(len = 3), intent(in) :: typeinfo   ! 'DVR' or 'GWP'
   double complex, dimension(:) :: psiHalf(dimpsi)
   double complex, dimension(:,:,:,:) :: DDHalf(nS00max,nS00max,nmodGWP,nstates)
   !
   double precision :: dt_done, stepsize, step_done, step_next
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
   select case(int_type_psi)
      case(1)   ! Runge-Kutta / Cash-Karp
         dt_done = dZero
         stepsize = dt
         do
            call rk_SPF_ck(psiHalf, DDHalf, i, typeinfo, stepsize, step_done, step_next)
            if (stepsize < min_stepsize) then
               write(2989,'(a,a,a,i0)') 'The step size for the ', typeinfo, ' mode n. ', i
               write(2989,'(a,f15.12)') 'has reached the value ', stepsize * au2fs, ' fs'
            end if
            dt_done = dt_done + step_done
            if (dt_done .ge. dt) exit   ! the propagation is done
            step_next = max(step_next, min_stepsize * 2)   !!! This is not really an adaptive step size !!
            stepsize = (dt - dt_done) / (int((dt - dt_done) / step_next) + 1)
            stepsize = min(stepsize, dt - dt_done)
         end do
      case(2)   ! Runge-Kutta / Dormand - Prince
         dt_done = dZero
         stepsize = dt
         do
            call rk_SPF_dp(psiHalf, DDHalf, i, typeinfo, stepsize, step_done, step_next)
            if (stepsize < min_stepsize) then
               write(2989,'(a,a,a,i0)') 'The step size for the ', typeinfo, ' mode n. ', i
               write(2989,'(a,f15.12)') 'has reached the value ', stepsize * au2fs, ' fs'
            end if
            dt_done = dt_done + step_done
            if (dt_done .ge. dt) exit   ! the propagation is done
            step_next = max(step_next, min_stepsize * 2)   !!! This is not really an adaptive step size !!
            stepsize = (dt - dt_done) / (int((dt - dt_done) / step_next) + 1)
            stepsize = min(stepsize, dt - dt_done)
         end do
      case default
         write(0,*) 'Error: integrator not yet implemented for the SPFs'
         stop
   end select
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   return
   end subroutine propagate_SPF


   
   subroutine A_to_B(A0,B0,ovl,DD0)
   ! Calculate B = D^H * S * A
   double complex, dimension(nAConftot), intent(in) :: A0
   double complex, dimension(nAConftot), intent(out) :: B0
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates), intent(in) :: DD0
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates,nstates), intent(in) :: ovl
   integer :: iel, iA, nAC, im, jel, i, nS
   double complex, dimension(:,:) :: SD(nS00max,nS00max)
   double complex, dimension(:) :: Baux1(nAconftot), Baux2(nAconftot)
   double complex, dimension(:), pointer :: pA1, pA2, pswap
   target :: Baux1, Baux2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   if (nmodGWP .gt. 0) then
      do iel = 1, nstates
         iA = AStateOffset(iel)
         nAC = nAConf(iel)
         call zcopy(nAC, A0(iA + 1), 1, Baux2, 1)
         call zlacgv(nAC, Baux2, 1)
         pA1 => Baux1
         pA2 => Baux2
         do i = 1, nmodGWP
            im = imodGWP(i)
            jel = merge(iel, 1, MultiSet(im))
            nS = nSPF(im,iel)
            ! Form the matrix S * D
            call zhemm('L','U',nS,nS,cOne,ovl(1,1,i,jel,jel),nS00max, &
                       DD0(1,1,i,jel),nS00max,cZero,SD,nS00max)
            ! Multiply A^H * S * D
            call Mult_GM_MatV(nSPF(1,iel),im,'N',nS,nS,SD,nS00max,pA2,pA1)
            pswap => pA2
            pA2   => pA1
            pA1   => pswap
         end do
         call zcopy(nAC, pA2, 1, B0(iA + 1), 1)    
         call zlacgv(nAC, B0(iA + 1), 1)
      end do
   else
      call zcopy(nAconftot, A0, 1, B0, 1)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   return
   end subroutine A_to_B

   
   
   subroutine B_to_A(A0,B0,DD0)
   ! Calculate A = D * B
   double complex, dimension(nAConftot), intent(out) :: A0
   double complex, dimension(nAConftot), intent(in)  :: B0
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates), intent(in) :: DD0
   integer :: iel, iA, nAC, im, jel, i, nS
   double complex, dimension(:) :: Baux1(nAconftot), Baux2(nAconftot)
   double complex, dimension(:), pointer :: pA1, pA2, pswap
   target :: Baux1, Baux2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   if (nmodGWP .gt. 0) then
      do iel = 1, nstates
         iA = AStateOffset(iel)
         nAC = nAConf(iel)
         call zcopy(nAC, B0(iA + 1), 1, Baux2, 1)
         call zlacgv(nAC, Baux2, 1)
         pA1 => Baux1
         pA2 => Baux2
         do i = 1, nmodGWP
            im = imodGWP(i)
            jel = merge(iel, 1, MultiSet(im))
            nS = nSPF(im,iel)
            ! Multiply B^H * D^H
            call Mult_GM_MatV(nSPF(1,iel),im,'C',nS,nS,DD0(1,1,i,jel),nS00max,pA2,pA1)
            pswap => pA2
            pA2   => pA1
            pA1   => pswap
         end do
         call zcopy(nAC, pA2, 1, A0(iA + 1), 1)    
         call zlacgv(nAC, A0(iA + 1), 1)
      end do
   else
      call zcopy(nAconftot, B0, 1, A0, 1)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!                           
   return
   end subroutine B_to_A


end module cmf2_base
