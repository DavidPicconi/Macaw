module derivatives
   use sysparam
   use globals
   use psidef
   use storage
   use hamilton
   use dvr
   use timingmod
   implicit none 
   private
   public :: Calc_Derivs, StepSizeDrop
!
   interface
      subroutine Calc_S00(psi0, i)
         use sysparam
         use psidef
         use storage
         implicit none
         integer, intent(in) :: i
         integer :: iel1, iel2, im, nmD, nS1, nS2, j1, j2, iP, iP1, iP2
         double complex, dimension(dimpsi), intent(in) :: psi0
         target :: psi0
      end subroutine Calc_S00
!   end interface
!
!
!   interface
      subroutine Calc_CY(psi0, i, DD0)     ! This interface is necessary because the argument DD0 is optional
         use psidef
         use storage, only : nS00max
         double complex, dimension(:), intent(in) :: psi0(dimpsi)
         integer, intent(in) :: i
         double complex, dimension(nS00max,nS00max,nmodGWP,nstates), optional :: DD0
         end subroutine Calc_CY
   end interface
contains


   subroutine Calc_Derivs(A0,psi0,dA,dpsi,tt)
   implicit none 
   integer :: info = 0
   integer :: i, k, im, iel, jel, iel1, iel2, iHam, nMD, kdof, iW, &
              nG, nS, nS1, nS2, nn1, kg, iOD, j1, j2, iP1, iP2, iA, nAC, iD
   integer, dimension(nmodes) :: sizeTens
   double precision, intent(in) :: tt
   double precision :: dummy
   double complex, dimension(:), intent(in) :: A0(nAConftot * npackets), psi0(dimpsi)
   double complex, dimension(:), intent(out) :: dA(nAConftot * npackets), dpsi(dimpsi)
   double complex, dimension(:,:) :: tau(nS00max,nS00max)
   ! For the SPF projector
   double complex, dimension(:,:) :: ovl(nSPFmax,nSPFmax), ovldot(nSPFmax,nSPFmax)
   ! Lapack
   double complex, external :: zdotc
   !
   double complex, dimension(:) :: A0copy(nAConftot * npackets)
   double complex, dimension(:), pointer :: pA1, pA2, pswap
   ! For the harmonic oscillator bath
   integer :: jQ, jP, jQQ, jQP, jPQ, jPP, lQ, lP, j, l
   double precision, dimension(npackets) :: WPpopm1
   double complex, dimension(npackets, npackets) :: Cdot, Caux, FC, Caux1, Caux2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   do i = 1, nmodGWP
      call Calc_S00(psi0,i)
      call Calc_S00m1(i)
   end do
   do im = 1, nmodes
      call Calc_rho1mode(im, A0)
      call Calc_h1mode(psi0, im)
   end do
   ! 
   if (npackets > 1) then
      iA = 1
      do iW = 1, npackets
         call norm_wf(A0(iA), psi0, WPpop(iW))
         iA = iA + nAConftot
      end do
      ! Regularize the wave packet populations to avoid division by zero
      ! ??? (This is done also into Calc_MeanField .... 
      do iW = 1, npackets
         dummy = eps_pop * exp(- (WPpop(iW) / eps_pop) ** 2)
         WPpopm1(iW) = WPpop(iW) / (WPpop(iW) ** 2 + dummy ** 2)
      enddo
   end if
   !
   call Calc_MeanField(A0,psi0,tt)
   call Calc_rhom1MF
   do i = 1, nmodGWP
      call Calc_CY(psi0,i)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Start the timer
   call continue_timer(1)
   !
! 1) Derivatives of primary SPFs 
   dpsi = cZero    ! important initialization
   do i = 1, nmodDVR
      im = imodDVR(i)
      nG = ngpm(i)
      nMD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         if (IsIdentity(im,iHam)) cycle    
         iel1 = Oper(-1,iHam,im)
         iel2 = Oper( 0,iHam,im)
         nS2 = nSPF(im,iel2)
         iP2 = ipsiSPF(im,iel2)
         call zcopy(nS2 * nG, psi0(iP2),1, psi_sect,1)      ! copy the section of psi0
                                                          ! corresponding to the SPF of the mode im for the state iel2
         ! Operate with the operator iHam on the SPFs of mode im
         do k = 1, nMD
            iOD = Oper(k,iHam,im)
            if (iOD .eq. 0) cycle
            kdof = iDOF(k,im)
            nn1 = ngp(kdof)
            iD = indDVR(kdof)
            select case(iOD)
               case(-1)
                  call MultTens(d1DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)  ! MultTens should work also for non-symmetric matrices
                  call zscal(nS2 * nG,ciOne,psi_sect(1),1)   ! +i because MultTens multiplies from the right
               case(-2)
                  call MultTens(d2DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)
               case(-3)
                  call MultTens(xd1DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)  ! MultTens should work also for non-symmetric matrices
                  call zscal(nS2 * nG,ciOne,psi_sect(1),1)   ! +i because MultTens multiplies from the right
               case(-4)
                  call MultTens(operDVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)
               case(1:20)
                  call MultTensDiag(gridPow(1,iOD,iD),iel2,im,k,psi_sect(1))
               case(201:)
                  call MultTensDiag(gridOper(1, iOD - 200,iD),iel2,im,k,psi_sect(1))
               case(:-200)
                  call MultTensProj(gridOper(1,-iOD - 200,iD),iel2,im,k,psi_sect(1),nS2 * nG)
            end select            
         end do
         nS1 = nSPF(im,iel1)
         iP1 = ipsiSPF(im,iel1)
         call zgemm('N','T',nG,nS1,nS2, der_factor, psi_sect(1),nG, rhom1MF(1,1,iHam,i), nSPFmax, &
                    cOne, dpsi(iP1), nG)
      end do
      ! Project out the previous SPFs
      do iel = 1, merge(nstates,1,MultiSet(im))   ! 1 if the mode is treated with a single-set of SPFs, otherwise nstates
         nS = nSPF(im,iel)
         iP1 = ipsiSPF(im,iel)
         ! Calculate the overlap matrix between SPFs ovl
         call zherk('U','C',nS,nG,dOne,psi0(iP1),nG,dZero,ovl,nSPFmax)
         ! Calculate the overlap between phi and d(phi)/dt
         call zgemm('C','N',nS,nS,nG,cOne, psi0(iP1),nG, dpsi(iP1),nG, &
                    cZero,ovldot,nSPFmax)
         ! Calculate ovl^-1 * ovldot and store the result on ovldot
         call zposv('U',nS,nS,ovl,nSPFmax,ovldot,nSPFmax,info)  
         ! Multiply psi0 by ovldot
         call zgemm('N','N',nG,nS,nS, -cOne, psi0(iP1),nG, ovldot,nSPFmax, cOne, dpsi(iP1),nG)
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (nmodGWP .eq. 0) then 
      call zcopy(nAConfTot * npackets,Adot_tmp(1),1,dA(1),1)
   else
! 2) Derivatives of the Gaussian parameters
      do i = 1, nmodGWP
         im = imodGWP(i)
         do iel = 1, nstates
            if (.not. MultiSet(im) .and. iel > 1) cycle
            nG = nSPF(im,iel) * nmodeDOF(im)
! 'Standard' regularisation 
            iP1 = ipsiSPF(im,iel)
            !call DirectSol(nG,CMat(1,1,i,iel),nCYmax,YVec(1,i,iel),dpsi(iP1))
            !call DirectTikhSol(nG,CMat(1,1,i,iel),nCYmax,YVec(1,i,iel),dpsi(iP1),eps_C)
            call StdRegSol(nG,CMat(1,1,iel,i),nCYmax,YVec(1,iel,i),dpsi(iP1),eps_C,'tikh')
            call zscal(nG,der_factor,dpsi(iP1),1)
            !call zscal(nG,-ciOne,dpsi(iP1),1)
         end do
      end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! 3) Derivatives of the configuration coefficients
   ! At this point: Adot_tmp = (-i * H * A)^H  (complex conjugate, calculated in Calc_MeanField) 
      ! Operate with the inverse of the overlap matrix (dA --> S^(-1) * dA)
      call zcopy(nAConftot * npackets, A0,1, A0copy,1)    ! Is this copy necessary?
      call zlacgv(nAConftot * npackets, A0copy,1)
      do iel = 1, nstates
         !iA = AStateOffset(iel)
         nAC = nAConf(iel)
         do iW = 1, npackets
            iA = APacketOffset(iW) + AStateOffset(iel)
            pA1 => Adot_tmp(iA + 1:iA + nAC)
            pA2 => Aux2
            sizeTens = nSPF(:,iel)
            do i = 1, nmodGWP
               im = imodGWP(i)
               jel = merge(iel,1,MultiSet(im))
               call Mult_HeM_MatV(sizeTens,im,S00m1(1,1,i,jel),nS00max,pA1,pA2)
               pswap => pA1
               pA1 => pA2
               pA2 => pswap
            end do   
            call zcopy(nAC,pA1(1),1,dA(iA + 1),1)
         end do
      ! At this point: dA = pA1 = (-i * S^(-1) * H * A)^H   (complex conjugate)
      ! Operate with the tau matrix (dA --> dA - S^(-1) * tau * A
         !call zcopy(nAC,A0(iA + 1),1,Aux1(1),1)
         !call zlacgv(nAC,Aux1(1),1)
         do i = 1, nmodGWP
            im = imodGWP(i)
            jel = merge(iel,1,MultiSet(im))
            nS = nSPF(im,iel)
            nMD = nmodeDOF(im)
            ! Multiply by the tau matrix
            kg = 1
            iP2 = ipsiSPF(im,iel) 
            do j2 = 1, nS
               do j1 = 1, nS
                  tau(j1,j2) = zdotc(nMD,Sa0Sm1(kg,j1,i,jel),1,dpsi(iP2),1)
               end do
               kg = kg + nMD
               iP2 = iP2 + nMD
            end do
            ! tau contains the derivatives of the phases ...
            iP2 = ipsiSPF(im,iel) - 1
            do j2 = 1, nS
               dummy = dZero
               do k = 1, nMD
                  dummy = dummy &
                     + real(psi0(iP2 + k)) * real(dpsi(iP2 + k)) / GWPa(k,j2,i)
               end do
               tau(j2,j2) = tau(j2,j2) + dummy * dOneHalf
               iP2 = iP2 + nMD
            end do
            !
            do iW = 1, npackets
               iA = APacketOffset(iW) + AStateOffset(iel)
               call Mult_GM_MatV(sizeTens,im,'C',nS,nS,tau,nS00max,A0copy(iA + 1),Aux2)
               ! Update dA
               call zaxpy(nAC,-cOne,Aux2(1),1,dA(iA + 1),1)
            end do
         end do
      end do
   end if
   call zlacgv(nAConftot * npackets,dA,1)
   ! Multiply by the Lagrange multipliers
   if (npackets > 1) then
      call zaxpy(nAConftot * npackets, complex(multLagDiag,dZero), A0, 1, dA, 1)
      call zhemm('R','L',nAConfTot,npackets,ciOne,multLag,npackets,A0,nAConfTot,cOne,dA,nAConfTot)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Derivatives of the cross-correlation functions (maybe move this to another subroutine)
   if (lHObath) then
      ! i * sum_a F_a * C_a
      FC = cZero
      lQ = iCQ
      lP = iCP
      do l = 1, nHObath
         call zhemm('L','U',npackets,npackets,complex(dZero, cos(wBath(l) * tt)), &
                    FQbath(1,1,l),npackets,psi0(lQ),npackets,cOne,FC,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero, sin(wBath(l) * tt)), &
                    FQbath(1,1,l),npackets,psi0(lP),npackets,cOne,FC,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero,-sin(wBath(l) * tt)), &
                    FPbath(1,1,l),npackets,psi0(lQ),npackets,cOne,FC,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero, cos(wBath(l) * tt)), &
                    FPbath(1,1,l),npackets,psi0(lP),npackets,cOne,FC,npackets)
         lQ = lQ + npackets**2
         lP = lP + npackets**2
      end do
      !
      ! CQ matrix
      !
      jQ  = iCQ
      jP  = iCP
      jQQ = iCQQ
      jQP = iCQP
      jPQ = iCPQ
      jPP = iCPP
      do j = 1, nHObath
         call zhemm('L','U',npackets,npackets,complex(dZero,-cos(wBath(j) * tt)), &
                    FQbath(1,1,j),npackets,psi0(jQQ),npackets,cZero,Cdot,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero,-sin(wBath(j) * tt)), &
                    FQbath(1,1,j),npackets,psi0(jQP),npackets,cOne,Cdot,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero, sin(wBath(j) * tt)), &
                    FPbath(1,1,j),npackets,psi0(jQQ),npackets,cOne,Cdot,npackets)    !!!! PERCHE' C'E' cZero QUI???
         call zhemm('L','U',npackets,npackets,complex(dZero,-cos(wBath(j) * tt)), &
                    FPbath(1,1,j),npackets,psi0(jQP),npackets,cOne,Cdot,npackets)
         !
         call zhemm('R','U',npackets,npackets,cOne, &
                    psi0(jQ),npackets, FC,npackets, cOne,Cdot,npackets)
         !call zgemm('N','N',npackets,npackets,npackets,cOne, &
         !           FC,npackets,psi0(jQ),npackets,cOne,Cdot,npackets)
         do concurrent (j1 = 1:npackets, j2 = 1:npackets)
            Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
         end do
         !do j2 = 1, npackets
         !   do j1 = 1, npackets
         !      Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
         !   end do
         !end do
         call zhemm('L','L',npackets,npackets,-ciOne, &
                    multLag,npackets,psi0(jQ),npackets,cOne,Cdot,npackets)
         do concurrent (j1 = 1:npackets, j2 = 1:npackets)
            Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))
         end do
         !do j2 = 1, npackets
         !   do j1 = 1, npackets
         !      Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))
         !   end do
         !end do
         call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jQ),npackets)
         !
         jQ  = jQ  + npackets**2
         jP  = jP  + npackets**2
         jQQ = jQQ + npackets**2
         jQP = jQP + npackets**2
         jPQ = jPQ + npackets**2
         jPP = jPP + npackets**2
      end do
      !
      ! CP matrix
      !
      jQ  = iCQ
      jP  = iCP
      jQQ = iCQQ
      jQP = iCQP
      jPQ = iCPQ
      jPP = iCPP
      do j = 1, nHObath
         call zhemm('L','U',npackets,npackets,complex(dZero,-cos(wBath(j) * tt)), &
                    FQbath(1,1,j),npackets,psi0(jPQ),npackets,cZero,Cdot,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero,-sin(wBath(j) * tt)), &
                    FQbath(1,1,j),npackets,psi0(jPP),npackets,cOne,Cdot,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero, sin(wBath(j) * tt)), &
                    FPbath(1,1,j),npackets,psi0(jPQ),npackets,cOne,Cdot,npackets)
         call zhemm('L','U',npackets,npackets,complex(dZero,-cos(wBath(j) * tt)), &
                    FPbath(1,1,j),npackets,psi0(jPP),npackets,cOne,Cdot,npackets)
         !
         call zhemm('R','U',npackets,npackets,cOne, &
                    psi0(jP),npackets, FC,npackets, cOne,Cdot,npackets)
         !call zgemm('N','N',npackets,npackets,npackets,cOne, &
         !           FC,npackets,psi0(jP),npackets,cOne,Cdot,npackets)
         do concurrent (j1 = 1:npackets, j2 = 1:npackets)
            Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
         end do
         !do j2 = 1, npackets
         !   do j1 = 1, npackets
         !      Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
         !   end do
         !end do
         call zhemm('L','L',npackets,npackets,-ciOne, &
                    multLag,npackets,psi0(jP),npackets,cOne,Cdot,npackets)
         do concurrent (j1 = 1:npackets, j2 = 1:npackets)
            Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))
         end do
         !do j2 = 1, npackets
         !   do j1 = 1, npackets
         !      Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))
         !   end do
         !end do
         call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jP),npackets)
         !
         jQ  = jQ  + npackets**2
         jP  = jP  + npackets**2
         jQQ = jQQ + npackets**2
         jQP = jQP + npackets**2
         jPQ = jPQ + npackets**2
         jPP = jPP + npackets**2
      end do
      if (HObathPTorder > 2) then
         ! QUALCOSA ANCORA NON TORNA QUI
         !
         ! CQQ matrix
         !
         jQQ = iCQQ
         do j = 1, nHObath
            call zgemm('N','N',npackets,npackets,npackets,cOne, &
                       FC,npackets,psi0(jQQ),npackets,cZero,Cdot,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
               end do
            end do
            call zhemm('L','L',npackets,npackets,-ciOne, &
                       multLag,npackets,psi0(jQQ),npackets,cOne,Cdot,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))
               end do
            end do
            call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jQQ),npackets)
            !
            jQQ = jQQ + npackets**2
         end do
         !!!!! QUESTA PARTE (CQP E CPQ) E' SBAGLIATA 
         !!!!  NON SONO HERMITIANI
         !
         ! CQP matrix
         !
         jQP = iCQP
         do j = 1, nHObath
            call zgemm('N','N',npackets,npackets,npackets,cOne, &
                       FC,npackets,psi0(jQP),npackets,cZero,Caux1,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Caux1(j1,j2) = WPpopm1(j1) * Caux1(j1,j2)
               end do
            end do
            call zhemm('L','L',npackets,npackets,-ciOne, &
                       multLag,npackets,psi0(jQP),npackets,cOne,Caux1,npackets)
            !
            call zgemm('N','C',npackets,npackets,npackets, cOne, &   ! +cOne here??
                       psi0(jQP),npackets,FC,npackets,cZero,Caux2,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Caux2(j1,j2) = WPpopm1(j2) * Caux2(j1,j2)
               end do
            end do
            call zhemm('R','L',npackets,npackets,ciOne, &
                       multLag,npackets,psi0(jQP),npackets,cOne,Caux2,npackets)
            !
            Caux = Caux1 + Caux2
            call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jQP),npackets)
            !
            jQP = jQP + npackets**2
         end do
         !
         ! CPQ matrix
         !
         jPQ = iCPQ
         do j = 1, nHObath
            call zgemm('N','N',npackets,npackets,npackets,cOne, &
                       FC,npackets,psi0(jPQ),npackets,cZero,Caux1,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Caux1(j1,j2) = WPpopm1(j1) * Caux1(j1,j2)
               end do
            end do
            call zhemm('L','L',npackets,npackets,-ciOne, &
                       multLag,npackets,psi0(jPQ),npackets,cOne,Caux1,npackets)
            !
            call zgemm('N','C',npackets,npackets,npackets, cOne, &    !!! + cOne here???
                       psi0(jPQ),npackets,FC,npackets,cZero,Caux2,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Caux2(j1,j2) = WPpopm1(j2) * Caux2(j1,j2)
               end do
            end do
            call zhemm('R','L',npackets,npackets,ciOne, &
                       multLag,npackets,psi0(jPQ),npackets,cOne,Caux2,npackets)
            !
            Caux = Caux1 + Caux2
            call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jPQ),npackets)
            !
            jPQ = jPQ + npackets**2
         end do
         !
         ! CPP matrix
         !
         jPP = iCPP
         do j = 1, nHObath
            call zgemm('N','N',npackets,npackets,npackets,cOne, &
                       FC,npackets,psi0(jPP),npackets,cZero,Cdot,npackets)
            do j2 = 1, npackets
               do j1 = 1, npackets
                  Cdot(j1,j2) = WPpopm1(j1) * Cdot(j1,j2)
               end do
            end do
            call zhemm('L','L',npackets,npackets,-ciOne, &
                       multLag,npackets,psi0(jPP),npackets,cOne,Cdot,npackets)
            do j2 = 1, npackets 
               do j1 = 1, npackets
                  Caux(j1,j2) = Cdot(j1,j2) + conjg(Cdot(j2,j1))    
               end do
            end do
            call zlacpy('All',npackets,npackets,Caux,npackets,dpsi(jPP),npackets)
            !
            jPP = jPP + npackets**2
         end do
      end if
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Suspend the timer
   call suspend_timer(1)
   return
   end subroutine Calc_Derivs



   subroutine StepSizeDrop(A0,psi0,tt)
! It is called when the step size goes below the minimal allowed value
! The density matrices, Gaussian overlaps and the C-matrices are written in the log files
   double precision, intent(in) :: tt
   double complex, dimension(:), intent(in) :: A0(nAConftot * npackets), psi0(dimpsi)
   integer :: i, iel, im, nS, j, LWORK, info
   double precision, dimension(:), allocatable :: Wv, RWORK
   double complex, dimension(:,:), allocatable :: UMat
   double complex, dimension(:), allocatable :: WORK
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Calculate the S, rho and C matrices
   do i = 1, nmodGWP
      call Calc_S00(psi0,i)
   end do
   do i = 1, nmodes
      call Calc_rho1mode(i,A0)
   end do
   do i = 1, nmodGWP
      call Calc_CY(psi0,i)
   end do
! Write information about the interruption of the propagation
   write(ilog_unit,*) 
   write(ilog_unit,'(a,f7.2,a)') 'Propagation interrupted at ', tt * au2fs, ' fs'
   write(ilog_unit,'(a)') '   (the integration step size became too small)'
! Write the rho matrices
   write(ilog_unit,*) 
   write(ilog_unit,*) 'rho matrices'
   LWORK = 3 * nSPFmax - 2
   allocate(UMat(nSPFmax,nSPFmax), Wv(nSPFmax), WORK(LWORK), RWORK(LWORK))
   do iel = 1, nstates
      do im = 1, nmodes
         nS = nSPF(im,iel)
         UMat = cZero
         forall(j = 1:nS) UMat(j,j:nS) = rho(j,j:nS,im,iel)
         ! write information
         write(ilog_unit,'(a,i0,a,i0)') '   State ', iel, ';   Mode: ', im
         do j = 1, nS
            write(ilog_unit,'(20(f8.4,",",f8.4,3x))') UMat(j,1:nS)
         end do
         ! write eigenvalues
         call zheev('V','U',nS,UMat,nS00max,Wv,WORK,LWORK,RWORK,info)
         write(ilog_unit,'(a,20(2x,f8.4))') '   Eigenvalues:', Wv(1:nS)
      end do
   end do
   deallocate(UMat,Wv,WORK,RWORK)
! Write the S matrices
   write(ilog_unit,*) 
   write(ilog_unit,*) 'S matrices'
   LWORK = 3 * nS00max - 2
   allocate(UMat(nS00max,nS00max), Wv(nS00max), WORK(LWORK), RWORK(LWORK))
   do iel = 1, nstates
      do i = 1, nmodGWP
         im = imodGWP(i)
         nS = nSPF(im,iel)
         UMat = cZero
         forall(j = 1:nS) UMat(j,j:nS) = S00(j,j:nS,i,iel,iel)
         ! write information
         write(ilog_unit,'(a,i0,a,i0)') '   State ', iel, ';   Mode: ', im
!         write(0,*) '(', nS, '("(",f10.7,",",f10.7,")",2x))'
!         write(formatstr,'(a,i0,a)') '(', nS, '("(",f10.7,",",f10.7,")",2x))'
!         write(0,*) formatstr
         do j = 1, nS
            write(ilog_unit,'(20(f7.4,",",f7.4,3x))') UMat(j,1:nS)
         end do
         ! write eigenvalues
         call zheev('V','U',nS,UMat,nS00max,Wv,WORK,LWORK,RWORK,info)
         write(ilog_unit,'(a,20(2x,f7.4))') '   Eigenvalues:', Wv(1:nS)
      end do
   end do
   deallocate(UMat,Wv,WORK,RWORK)
! Write the C matrices
   write(ilog_unit,*) 
   write(ilog_unit,*) 'C matrices'
   LWORK = 3 * nCYmax - 2
   allocate(UMat(nCYmax,nCYmax), Wv(nCYmax), WORK(LWORK), RWORK(LWORK))
   do iel = 1, nstates
      do i = 1, nmodGWP
         im = imodGWP(i)
         nS = nSPF(im,iel) * nmodeDOF(im)
         UMat = cZero
         forall(j = 1:nS) UMat(j,j:nS) = CMat(j,j:nS,i,iel)
         ! write information
         write(ilog_unit,'(a,i0,a,i0)') '   State ', iel, ';   Mode: ', im
         do j = 1, nS
            write(ilog_unit,'(20(f7.4,",",f7.4,3x))') UMat(j,1:nS)
         end do
         ! write eigenvalues
         call zheev('V','U',nS,UMat,nCYmax,Wv,WORK,LWORK,RWORK,info)
         write(ilog_unit,'(a,20(2x,f7.4))') '   Eigenvalues:', Wv(1:nS)
      end do
   end do
   deallocate(UMat,Wv,WORK,RWORK)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine StepSizeDrop

end module derivatives
