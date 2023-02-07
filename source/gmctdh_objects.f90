subroutine alloc_AVec
use sysparam
use psidef
implicit none
integer :: a_status, i
integer :: nTens
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Allocate Aux1 and Aux2
nTens = 1
do i = 1, nmodes
   nTens = nTens * maxval(nSPF(i,:))
end do
allocate(Aux1(nTens), Aux2(nTens), stat = a_status)
if (a_status .ne. 0) call err_alloc('Aux1/Aux2','alloc_AVec',a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine alloc_Avec



subroutine norm_wf(A0,psi0,wfnorm)
! It calculates the squared norm of a wave packet
use sysparam
use psidef
use dvr
use storage
use timingmod
implicit none
integer :: i, j, j1, j2, iel, im, iP, nmD, nS, nAC, iA, nG
double precision, intent(out) :: wfnorm
integer, dimension(:), pointer :: sizeTens
double complex, dimension(:), intent(in) :: A0(nAConftot), psi0(dimpsi)
double complex, dimension(:,:) :: ovl(nSPFmax,nSPFmax)
double complex, external :: zdotu   ! Lapack 
double complex, dimension(:), pointer :: pA1, pA2, pswap
!
interface
   pure double complex function GauOvlMultiD(A1,A2,xi1,xi2,nmD)
      integer, intent(in) :: nmD
      double precision, dimension(nmD), intent(in) :: A1, A2
      double complex, dimension(nmD), intent(in) :: xi1, xi2
   end function GauOvlMultiD
end interface
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Continue the timer
call continue_timer(2)
!
wfnorm = dZero
do iel = 1, nstates
   iA = AStateOffset(iel) 
   nAC = nAConf(iel)
   call zcopy(nAC,A0(iA + 1),1,Aux1(1),1)
   call zlacgv(nAC,Aux1(1),1)
   sizeTens => nSPF(:,iel)
   pA1 => Aux1
   pA2 => Aux2
   ! DVR modes
   do i = 1, nmodDVR
      im = imodDVR(i)
      nS = nSPF(im,iel)
      nG = ngpm(i)
      ! Calculate the overlap matrix
      iP = ipsiSPF(im,iel)    ! If the mode is single-set then ipsiSPF(im,iel) = ipsiSPF(im,1)
      call zherk('U','C',nS,nG,dOne,psi0(iP),nG,dZero,ovl,nSPFmax)
      ! Multiply the A vector by the overlap matrix
      call Mult_HeM_MatV(sizeTens,im,ovl,nSPFmax,pA1,pA2)
      pswap => pA1
      pA1 => pA2
      pA2 => pswap
   end do
   ! GWP modes
   do concurrent (j = 1:nS00max) 
      ovl(j,j) = cOne
   end do
   do i = 1, nmodGWP
      im = imodGWP(i)
      nmD = nmodeDOF(im)
      nS = nSPF(im,iel)
      ! Calculate the overlap matrix
      iP = ipsiSPF(im,iel) - nmD  
      do concurrent (j1 = 1:nS, j2 = 1:nS, j1 .lt. j2) 
         ovl(j1,j2) = GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi0(iP + j1 * nmD),psi0(iP + j2 * nmD),nmD)
      end do
      ! Multiply the A vector by the overlap matrix
      call Mult_HeM_MatV(sizeTens,im,ovl,nSPFmax,pA1,pA2)
      pswap => pA1
      pA1 => pA2
      pA2 => pswap
   end do
   wfnorm = wfnorm + real(zdotu(nAC,pA1(1),1,A0(iA + 1),1))
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Suspend the timer
call suspend_timer(2)
return
end subroutine norm_wf



subroutine WF_Ovl(A1,psi1,A2,psi2,overlap)
! It calculates the overlap between two wavefunctions
use sysparam
use psidef
use dvr
use timingmod
implicit none
integer :: iel, iA, nAC, i, im, nS, j1, j2, iP, nG, nmD
integer, dimension(:), pointer :: sizeTens
double complex, intent(out) :: overlap
double complex, dimension(:), intent(in) :: A1(nAConftot), A2(nAConftot), &
                                            psi1(dimpsi), psi2(dimpsi)
double complex, dimension(:,:) :: ovl(nSPFmax,nSPFmax)
double complex, external :: zdotu   ! Lapack 
double complex, dimension(:), pointer :: pA1, pA2, pswap
!
interface
   pure double complex function GauOvlMultiD(A1,A2,xi1,xi2,nmD)
      integer, intent(in) :: nmD
      double precision, dimension(nmD), intent(in) :: A1, A2
      double complex, dimension(nmD), intent(in) :: xi1, xi2
   end function GauOvlMultiD
end interface
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Continue the timer
call continue_timer(2)
!
overlap = cZero
do iel = 1, nstates
   iA = AStateOffset(iel) 
   nAC = nAConf(iel)
   call zcopy(nAC,A2(iA + 1),1,Aux1(1),1)
   call zlacgv(nAC,Aux1(1),1)
   sizeTens => nSPF(:,iel)
   pA1 => Aux1
   pA2 => Aux2
   !
   ! DVR modes
   !
   do i = 1, nmodDVR
      im = imodDVR(i)
      nS = sizeTens(im)
      nG = ngpm(i)
      ! Calculate the overlap matrix
      iP = ipsiSPF(im,iel)
      call zgemm('C','N',nS,nS,nG, cOne,psi1(iP),nG, psi2(iP),nG, cZero,ovl,nSPFmax)
      ! Multiply by the overlap matrix
      call Mult_GM_MatV(sizeTens,im,'C',nS,nS,ovl,nSPFmax,pA1,pA2)
      pswap => pA1
      pA1 => pA2
      pA2 => pswap
   end do
   !
   ! GWP modes
   !
   do i = 1, nmodGWP
      im = imodGWP(i)
      nS = nSPF(im,iel)
      nmD = nmodeDOF(im)
      ! Calculate the overlap matrix
      iP = ipsiSPF(im,iel) - nmD
      do concurrent (j1 = 1:nS, j2 = 1:nS)
         ovl(j1,j2) = GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi1(iP + j1 * nmD),psi2(iP + j2 * nmD),nmD)
      end do
      ! Multiply by the overlap matrix
      call Mult_GM_MatV(sizeTens,im,'C',nS,nS,ovl,nSPFmax,pA1,pA2)
      pswap => pA1
      pA1 => pA2
      pA2 => pswap
   end do
   ! 
   overlap = overlap + zdotu(nAC,pA1(1),1,A1(iA + 1),1)
end do
overlap = conjg(overlap)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Suspend the timer
call suspend_timer(2)
return
end subroutine WF_Ovl



subroutine Trace_wf(A0,psi0,trace)
! Calculate the sum of the norms of the wave packets        
use sysparam
use psidef
double complex, dimension(:), intent(in) :: A0(nAConftot * npackets), &
                                            psi0(dimpsi)
double precision, intent(out) :: trace
integer :: iW, iA
double precision :: dummy
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
trace = dZero
do iW = 1, npackets
   iA = APacketOffset(iW) + 1
   call norm_WF(A0(iA), psi0, dummy)
   trace = trace + dummy
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Trace_wf



subroutine DeltaPsi(A1,psi1,A2,psi2,delta)
!! OLD:        
!! Calculate delta = Tr(drho,drho)^(1/4)
!! where drho = rho1 / Tr(rho1) - rho2 / Tr(rho2)
! 
! NEW:
! Calculate delta = sum_J (2 - 2 * Re <psi_1J|psi_2J> / sqrt(<psi_1J|psi_1J> * <psi_2J|psi_2J>))^(1/2)
use sysparam
use psidef
double complex, dimension(nAConftot * npackets), intent(in) :: A1, A2
double complex, dimension(dimpsi) :: psi1, psi2
double precision, intent(out) :: delta 
integer :: iW, iA
double precision :: pop, ovl_11, ovl_22, dummy
double complex :: ovl_12
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
if (npackets == 1) then
   call norm_WF(A1, psi1, ovl_11)
   call norm_WF(A2, psi2, ovl_22)
   call WF_ovl(A1, psi1, A2, psi2, ovl_12)
   delta = sqrt(abs(2 * (dOne - dble(ovl_12) / sqrt(ovl_11 * ovl_22)) ))
else        
   delta = dZero
   pop = dZero
   iA = 1
   do iW = 1, npackets
      call norm_WF(A1(iA), psi1, ovl_11)
      call norm_WF(A2(iA), psi2, ovl_22)
      call WF_ovl(A1(iA), psi1, A2(iA), psi2, ovl_12)
      !delta = delta + sqrt(abs(2 * (dOne - dble(ovl_12) / sqrt(ovl_11 * ovl_22)) )) * ovl_11
      !pop1 = pop1 + ovl_11
      dummy = sqrt(ovl_11 * ovl_22)
      pop = pop + dummy
      delta = delta + abs(dummy - real(ovl_12))
      !
      iA = iA + nAConftot
   end do
   delta = sqrt(2.d0 * delta / pop)
   !delta = delta / pop1
   return
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ! OLD VERSION (it does not work well)
!!! ! Calculate the wave packet norms
!!! do iW = 1, npackets
!!!    iA = APacketOffset(iW) + 1
!!!    call norm_WF(A1(iA), psi1, ovl_11(iW))
!!!    call norm_WF(A2(iA), psi2, ovl_22(iW))
!!! end do
!!! ! Population and purities
!!! pop1 = sum(ovl_11)
!!! pop2 = sum(ovl_22)
!!! pur1 = sum(ovl_11 ** 2)
!!! pur2 = sum(ovl_22 ** 2)
!!! ! Tr(rho1 * rho2)
!!! delta = dZero
!!! do iW1 = 1, npackets
!!!    iA1 = APacketOffset(iW1) + 1
!!!    do iW2 = 1, npackets
!!!       iA2 = APacketOffset(iW2) + 1
!!!       call WF_ovl(A1(iA1),psi1, A2(iA2),psi2,ovl_12)
!!!       delta = delta - 2 * dble(conjg(ovl_12) * ovl_12)
!!!    end do
!!! end do
!!! !
!!! delta = delta / (pop1 * pop2) + pur1  / pop1**2 + pur2 / pop2**2 
!!! delta = sqrt(sqrt(delta))   ! Maybe one sqrt???
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
!return
end subroutine DeltaPsi  



subroutine DeltaPsi_WF(A1,psi1,A2,psi2,delta)
! Calculate delta = ||Psi1 - Psi2|| / (||Psi1|| * ||Psi2||)
use sysparam
use psidef, only : nAConftot, dimpsi
double complex, dimension(:), intent(in) :: A1(nAConftot), A2(nAConftot), &
                                            psi1(dimpsi), psi2(dimpsi)
double precision, intent(out) :: delta
double precision :: ovl_11, ovl_22
double complex :: ovl_12
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
call norm_WF(A1,psi1,ovl_11)
call norm_WF(A2,psi2,ovl_22)
call WF_ovl(A1,psi1,A2,psi2,ovl_12)
delta = 2.d0 * (dOne - real(ovl_12) / sqrt(ovl_11 * ovl_22))
delta = sqrt(abs(delta))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DeltaPsi_WF



subroutine Calc_SWP(A0,psi0,SWP)
! Calculate the overlap matrix between wave packets      
use sysparam
use psidef
use storage
double complex, dimension(nAConftot * npackets), intent(in) :: A0
double complex, dimension(dimpsi), intent(in) :: psi0
double complex, dimension(npackets,npackets), intent(out) :: SWP
integer :: iW1, iA1, iW2, iA2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do iW1 = 1, npackets
   iA1 = ApacketOffset(iW1) + 1
   do iW2 = 1, iW1
      iA2 = ApacketOffset(iW2) + 1
      call WF_Ovl(A0(iA1),psi0, A0(iA2),psi0, SWP(iW1,iW2))
      SWP(iW2,iW1) = conjg(SWP(iW1,iW2))
   end do
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_SWP


subroutine Calc_S00(psi0, i)
! It calculates the overlap matrices for GWP mode i
use sysparam
use psidef
use storage
implicit none
double complex, dimension(dimpsi), intent(in) :: psi0
integer, intent(in) :: i
integer :: iel1, iel2, im, nmD, nS1, nS2, j1, j2, iP, iP1, iP2
!
interface
   pure double complex function GauOvlMultiD(A1,A2,xi1,xi2,nmD)
      integer, intent(in) :: nmD
      double precision, dimension(nmD), intent(in) :: A1, A2
      double complex, dimension(nmD), intent(in) :: xi1, xi2
   end function GauOvlMultiD
end interface
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
im = imodGWP(i)
nmD = nmodeDOF(im)
! Electronic diagonal
do iel1 = 1, merge(nstates, 1, MultiSet(im))
   nS1 = nSPF(im,iel1)
   iP = ipsiSPF(im,iel1) - nmD
   do concurrent (j1 = 1:nS1, j2 = 1:nS1, j1 < j2) 
      S00(j1,j2,i,iel1,iel1) = GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                            psi0(iP + j1 * nmD),psi0(iP + j2 * nmD),nmD)
   end do
   do concurrent (j2 = 1:nS1, j1 = 1:nS1, j1 > j2) 
      S00(j1,j2,i,iel1,iel1) = conjg(S00(j2,j1,i,iel1,iel1))
   enddo
end do
if (.not. MultiSet(im)) return
! Different electronic states
do iel2 = 1, nstates
   do iel1 = 1, iel2 -1
      nS1 = nSPF(im,iel1)
      nS2 = nSPF(im,iel2)
      iP1 = ipsiSPF(im,iel1) - nmD
      iP2 = ipsiSPF(im,iel2) - nmD
      do concurrent (j1 = 1:nS1, j2 = 1:nS2) 
         S00(j1,j2,i,iel1,iel2) = GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                               psi0(iP1 + j1 * nmD),psi0(iP2 + j2 * nmD),nmD)
      end do
   end do
end do
! Hermitian conjugate
do iel2 = 1, nstates
   do iel1 = iel2 + 1, nstates
      nS1 = nSPF(im,iel1)
      nS2 = nSPF(im,iel2)
      do concurrent (j1 = 1:nS1, j2 = 1:nS2)
         S00(j1,j2,i,iel1,iel2) = conjg(S00(j2,j1,i,iel2,iel1))
      enddo
   end do
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_S00



subroutine Calc_S00m1(i)
! It calculates the inverse overlap matrices for the GWP mode i
! It follows the regularisation used for the density matrix
use sysparam
use psidef
use storage
use globals
implicit none
integer, intent(in) :: i
integer :: iel, nS
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do iel = 1, merge(nstates, 1, MultiSet(imodGWP(i)))
   nS = nSPF(imodGWP(i),iel)
   ! Standard regularised inverse
   !call DirectTikhInv(nS,S00(1,1,i,iel,iel),nS00max,S00m1(1,1,i,iel),nS00max,eps_ovl)
   call StdRegInv(nS,S00(1,1,i,iel,iel),nS00max,S00m1(1,1,i,iel),nS00max,eps_ovl,'S00 ','tikh')
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_S00m1



subroutine Calc_rho1mode(im,A0)
! It calculates the density matrix for the mode im
! - SPFs for non-GWP modes are assumed orthonormal
! - overlap matrices for GWP modes are read from the array S00 (which must be pre-calculated)
use sysparam
use psidef
use storage
implicit none
integer, intent(in) :: im
double complex, dimension(nAConftot * npackets), intent(in) :: A0
integer :: iel, jel, i, im2, nS, iA, nAC, j1, j2, iW, &
           loop_length, nloops, iA1, iA2, iAsize   ! (for the final building of the density matrix)
integer, dimension(:), pointer :: sizeTens
double complex, dimension(:), pointer :: pA1, pA2, pswap
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do iel = 1, nstates
   nAC = nAConf(iel)
   sizeTens => nSPF(:,iel)
   nS = sizeTens(im)
   ! Initialise the density matrix to zero
   call zlaset('U',nS,nS,cZero,cZero,rho(1,1,im,iel),nSPFmax)
   loop_length = product(sizeTens(1:im - 1))
   nloops = product(sizeTens(im + 1:nmodes))
   iAsize = nS * loop_length
   ! Loop over wave packets
   do iW = 1, npackets
      iA = APacketOffset(iW) + AStateOffset(iel)
      call zcopy(nAC, A0(iA + 1),1, Aux1(1),1)
      call zlacgv(nAC,Aux1(1),1)
      pA1 => Aux1
      pA2 => Aux2
      do i = 1, nmodGWP
         im2 = imodGWP(i)
         if (im2 .eq. im) cycle
         jel = merge(iel, 1, MultiSet(im2))
         call Mult_HeM_MatV(sizeTens,im2,S00(1,1,i,jel,jel),nS00max,pA1,pA2)
         pswap => pA1
         pA1 => pA2
         pA2 => pswap
      enddo
      ! Update the density matrix
      iA1 = iA + 1
      iA2 = 1
      call zlacgv(nAC,pA1(1),1)
      do i = 1, nloops
         call zher2k('U','C',nS,loop_length,cOneHalf,pA1(iA2),loop_length,A0(iA1),loop_length, &
                     dOne,rho(1,1,im,iel),nSPFmax)
         iA1 = iA1 + iAsize
         iA2 = iA2 + iAsize
      enddo
   enddo
   !
enddo
!
! For the single-set modes, calculate the trace of rho over the electronic states
!
if (.not. MultiSet(im)) then
   nS = nSPF(im,1)
   call zlacpy('U', nS, nS, rho(1,1,im,1), nSPFmax, rho_ave(1,1,im), nSPFmax)
   do iel = 2, nstates
      do j2 = 1, nS
         do j1 = 1, j2
            rho_ave(j1,j2,im) = rho_ave(j1,j2,im) + rho(j1,j2,im,iel)
         end do
      end do
   end do
end if  
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_rho1mode



subroutine Calc_MeanField(A0,psi0,tt)
! It calculates all mean fields
use sysparam
use globals, only : der_factor, eps_pop
use psidef
use storage
use hamilton
use timingmod
use dvr, only : dvrType
implicit none
integer :: iHam, iel1, iel2, iA, iA1, iA2, nAC1, nAC2, im, im2, jm, iW, iW2, &
           nS1, nS2, iOp, iB1, iB2, i, j, iB1size, iB2size, nloops, loop_length, &
           iDissOp, iExtraOp, iT, jQ, jP, iExtraOpQ, iExtraOpP
integer, dimension(nmodes) :: sizeTens
double precision, intent(in) :: tt
double precision :: dummy, dw, cwt, swt
double complex :: aet, cdummy
double complex, dimension(nAConftot * npackets), intent(in) :: A0
double complex, dimension(dimpsi), intent(in) :: psi0    ! needed for the harmonic perturbative bath
double complex, dimension(nAConftot,0:nModes - 1) :: AL
double complex, dimension(nAConftot,0:nModes) :: AR
double complex, external :: zdotu, zdotc
! For the dissipative terms
double precision, dimension(npackets) :: WPpopm1
double complex, dimension(nAConfTot * npackets) :: VA, AuxA, VQA, VPA
double complex, dimension(npackets,npackets) :: AVA, AVAAVpA, cFQsFP, sFQcFP
!
double complex, dimension(:), pointer :: pA1, pA2, pswap
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Start timer
call continue_timer(3)
!
if (nHamTD .gt. 0) call UpdateTDcoeff(tt)   ! Update time-dependent coefficients 
!
MeanField = cZero  ! This initialisation is important!
Adot_tmp = cZero   ! This initialisation is important!
energyMF = dZero 
do iHam = 1, nHam
   iel1 = lhs_el(iHam)
   iel2 = rhs_el(iHam)
   nAC1 = nAConf(iel1)
   nAC2 = nAConf(iel2)
   !
   do iW = 1, npackets
      iA1 = APacketOffset(iW) + AStateOffset(iel1)
      iA2 = APacketOffset(iW) + AStateOffset(iel2)
      !
      ! Get the right-multiplied vectors: H_f x A, H_(f-1) x H_f x A, ...
      !
      sizeTens = nSPF(:,iel2)
      call zcopy(nAC2, A0(iA2 + 1), 1, AR(1,0), 1)
      call zlacgv(nAC2, AR(1,0), 1)
      do im = 1, nmodes
         jm = nmodes - im + 1
         if (SkipTerm(jm,iHam)) then  ! Do not multiply by the identity matrix
             call zcopy(nAC2, AR(1,im - 1), 1, AR(1,im), 1)
             cycle
         end if
         ! Multiply by the 1-mode Hamiltonian matrix
         nS1 = nSPF(jm,iel1)
         nS2 = nSPF(jm,iel2)
         iOp = iOper(jm,iHam)
         if (iel1 == iel2 .or. .not. MultiSet(jm)) then
            call Mult_HeM_MatV(sizeTens,jm,h1mode(1,1,iOp,jm),nSPFmax,AR(1,im - 1),AR(1,im))
         else
            call Mult_GM_MatV(sizeTens,jm,'C',nS2,nS1,h1mode(1,1,iOp,jm),nSPFmax,AR(1,im - 1),AR(1,im))
            sizeTens(jm) = nS1
            nAC2 = nAC2 * nS1 / nS2
         end if
      end do
      !
      ! Get the left-multiplied vectors: A x H_1, A x H_1 x H_2, ...
      !
      sizeTens = nSPF(:,iel1)
      call zcopy(nAC1,A0(iA1 + 1), 1, AL(1,0), 1)
      call zlacgv(nAC1, AL(1,0), 1)
      do im = 1, nmodes - 1
         if (SkipTerm(im,iHam)) then
            call zcopy(nAC1, AL(1,im - 1), 1, AL(1,im), 1)
            cycle
         end if
         ! Multiply by the 1-mode Hamiltonian matrix
         nS1 = nSPF(im,iel1)
         nS2 = nSPF(im,iel2)
         iOp = iOper(im,iHam)
         if (iel1 == iel2 .or. .not. MultiSet(im)) then
            call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,AL(1,im - 1),AL(1,im))
         else
            call Mult_GM_MatV(sizeTens,im,'N',nS1,nS2,h1mode(1,1,iOp,im),nSPFmax,AL(1,im - 1),AL(1,im))
            sizeTens(im) = nS2
            nAC1 = nAC1 * nS2 / nS1
         end if
      end do
      !
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      ! Calculate Adot_tmp (the complex conjugate of -i * H * A)
      call zaxpy(nAC1, conjg(der_factor * ae(iHam)), AR(1,nModes), 1, Adot_tmp(iA1 + 1),1)
      ! Calculate the energy
      energyMF = energyMF + real(conjg(ae(iHam)) * zdotu(nAC1, A0(iA1 + 1), 1, AR(1,nModes), 1))
      ! Calculate the mean field
      do im = 1, nmodes
         if (IsIdentity(im,iOper(im,iHam))) cycle  ! terms associated with the identity are projected out
         jm = nmodes - im
         iOp = iOper(im,iHam)
         nS1 = nSPF(im,iel1)
         nS2 = nSPF(im,iel2)
         loop_length = product(nSPF(1:im - 1,iel2))
         nloops = product(nSPF(im + 1:nmodes,iel1))
         iB1size = nS1 * loop_length 
         iB2size = nS2 * loop_length
         call zlacgv(loop_length * nS2 * nloops, AR(1,jm), 1)
         iB1 = 1
         iB2 = 1
         do j = 1, nloops
            call zgemm('T','N',nS1,nS2,loop_length,ae(iHam),AL(iB1,im - 1),loop_length, &
                       AR(iB2,jm),loop_length,cOne,MeanField(1,1,iOp,im),nSPFmax)
            iB1 = iB1 + iB1size 
            iB2 = iB2 + iB2size
         end do
      end do
   !
   end do   ! close the loop over the wave packets (iW)
end do
!
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
!
! DISSIPATIVE TERMS
!
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
multLagDiag = dZero
call zlaset('L',npackets,npackets,cZero,cZero,multLag,npackets)
!multLag = cZero      ! Lagrange multipliers to impose population conservation and orthogonality
                     ! They are calculated here for being used in the derivative module
! Regularize the wave packet populations to avoid division by zero
do j = 1, npackets
   dummy = eps_pop * exp(- (WPpop(j) / eps_pop) ** 2)
   WPpopm1(j) = WPpop(j) / (WPpop(j) ** 2 + dummy ** 2)
enddo
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Type 1 dissipation terms:  alpha * (V rho V^+)
!    (the global prefactor of the operator has already been divided by two)
do iDissOp = 1, nDiss1
   iExtraOp = iDiss1(iDissOp)
   ! Calculate the matrix [A^+ * V * A]
   VA = cZero
   do iT = 1, nOpTerms(iExtraOp)   
      aet = conjg(aeExtra(iT,iExtraOp))
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      do iW = 1, npackets
         iA2 = APacketOffset(iW) + AStateOffset(iel2)
         call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1), 1)
         call zlacgv(nAC2,Aux2(1),1)
         pA1 => Aux1
         pA2 => Aux2
         sizeTens = nSPF(:,iel2)
         do im = 1, nmodes
            if (dvrType(iDOF(1,im)) .ne. 'gwp' .and. IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle    ! SkipTerm
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            ! Multiply by the 1-mode-Hamiltonian matrix
            iOp = iOperExtra(im,iT,iExtraOp)
            if (iel1 .eq. iel2 .or. .not. MultiSet(im)) then  ! h1mode is assumed to be Hermitian
               call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            else
               call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            end if
            pswap => pA2
            pA2   => pA1
            pA1   => pswap
            sizeTens(im) = nS1
         end do
         iA1 = APacketOffset(iW) + AStateOffset(iel1)
         call zaxpy(nAC1,aet,pA2(1),1,VA(iA1 + 1),1)
      end do
   end do
   ! Build the matrix
   call zgemm('T','N',npackets,npackets,nAConfTot,cOne,A0,nAConftot,VA,nAConfTot,cZero,AVA,npackets)
   call zlacgv(npackets * npackets,AVA(1,1),1)
   !
   ! Update the Lagrange multiplier, which ensures the conservation of population
   dummy = dZero
   do j = 1, npackets
      do i = 1, npackets
         dummy = dummy - (real(AVA(i,j))**2 + aimag(AVA(i,j))**2) * WPpopm1(i)
         !dummy = dummy - real(conjg(AVA(i,j)) * AVA(i,j)) * WPpopm1(i)
      enddo
   enddo
   multLagDiag = multLagDiag + dummy * aeD1(iDissOp)
   ! Calculate the product [A^+ V A] * [A^+ V^+ A]
   call zherk('L','N',npackets,npackets,aeD1(iDissOp),AVA,npackets,dZero,AVAAVpA,npackets)
   forall (j = 1:npackets,i = 1:npackets, i .gt. j) &
      multLag(i,j) = multLag(i,j) + AVAAVpA(i,j) * (WPpopm1(i) + WPpopm1(j))
   ! Update Adot_tmp  (the complex conjugate of S * dA/dt)
   call zgemm('N','T',nAConfTot,npackets,npackets, &
              complex(aeD1(iDissOp),dZero),VA,nAConfTot,AVA,npackets,cZero,AuxA,nAConfTot)
   iA = 1
   do iW = 1, npackets
      call zaxpy(nAConfTot, complex(WPpopm1(iW),dZero), AuxA(iA),1, Adot_tmp(iA),1)
      iA = iA + nAConfTot
   end do
   !
   call zgemm('N','C',nAConfTot,npackets,npackets,complex(dZero,aeD1(iDissOp)), &
              A0,nAConftot, AVA,npackets, cZero,AuxA,nAConfTot)
   iA = 1
   do iW = 1, npackets
      call zdscal(nAConfTot, WPpopm1(iW), AuxA(iA), 1)
      iA = iA + nAConfTot
   end do
   !
   ! Evaluate the mean field
   !
   do iT = 1, nOpTerms(iExtraOp)
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)      
      do im = 1, nmodes
         if (IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle
         loop_length = product(nSPF(1:im - 1,iel1))
         nloops = product(nSPF(im + 1:nmodes,iel1))
         iB1size = nSPF(im,iel1) * loop_length
         iB2size = nSPF(im,iel2) * loop_length      
         do iW = 1, npackets
            iA2 = APacketOffset(iW) + AStateOffset(iel2)
            call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1),1)
            call zlacgv(nAC2,Aux2(1),1)
            pA1 => Aux1
            pA2 => Aux2
            sizeTens = nSPF(:,iel2)
            do im2 = 1, nmodes
               if (im2 == im) cycle    ! do not operate on the MeanField mode or with the identity matrix
               if (dvrType(iDOF(1,im2)) .ne. 'gwp' .and. IsIdentity(im2,iOperExtra(im2,iT,iExtraOp))) cycle   ! Skip term
               nS1 = nSPF(im2,iel1)
               nS2 = nSPF(im2,iel2)
               ! Multiply by the 1-mode-Hamiltonian matrix
               iOp = iOperExtra(im2,iT,iExtraOp)
               if (iel1 == iel2 .or. .not. MultiSet(im2)) then   ! Assume Hermitian h1mode
                  call Mult_HeM_MatV(sizeTens,im2,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               else
                  call Mult_GM_MatV(sizeTens,im2,'C',nS2,nS1,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               end if
               pswap => pA2
               pA2   => pA1
               pA1   => pswap
               sizeTens(im2) = nS1
            end do
            ! At this stage <pA2| = < psi_iW^(im) | V^+ , where psi_iwp^(im) is a single-hole function
            ! Update the mean-field matrix
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            iOp = iOperExtra(im,iT,iExtraOp)
            !
            iB1 = APacketOffset(iW) + AStateOffset(iel1) + 1
            iB2 = 1
            do j = 1, nloops
               call zgemm('T','N',nS2,nS1,loop_length, aeExtra(iT,iExtraOp), pA2(iB2),loop_length, &
                          AuxA(iB1), loop_length,cOne, MeanField(1,1,iOp,im), nSPFmax)
               iB1 = iB1 + iB1size
               iB2 = iB2 + iB2size
            end do
         end do
      end do
   end do
end do
!
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Type 2 dissipation terms: beta * (W * rho + rho * W^+)
!    the same code for the standard Hamiltonian can be used
!
do iDissOp = 1, nDiss2
   iExtraOp = iDiss2(iDissOp)
   do iT = 1, nOpTerms(iExtraOp)
      aet = aeD2(iDissOp) * aeExtra(iT,iExtraOp)
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      !
      do iW = 1, npackets
         iA1 = APacketOffset(iW) + AStateOffset(iel1)
         iA2 = APacketOffset(iW) + AStateOffset(iel2)
         !
         ! Get the right-multiplied vectors: H_f x A, H_(f-1) x H_f x A, ...
         !
         sizeTens = nSPF(:,iel2)
         call zcopy(nAC2, A0(iA2 + 1), 1, AR(1,0), 1)
         call zlacgv(nAC2, AR(1,0), 1)
         do im = 1, nmodes
            jm = nmodes - im + 1
            if (dvrType(iDOF(1,jm)) .ne. 'gwp' .and. IsIdentity(jm,iOperExtra(jm,iT,iExtraOp))) then  ! Do not multiply by the identity matrix
                call zcopy(nAC2, AR(1,im - 1), 1, AR(1,im), 1)
                cycle
            end if
            ! Multiply by the 1-mode Hamiltonian matrix
            nS1 = nSPF(jm,iel1)
            nS2 = nSPF(jm,iel2)
            iOp = iOperExtra(jm,iT,iExtraOp)
            if (iel1 == iel2 .or. .not. MultiSet(jm)) then
               call Mult_HeM_MatV(sizeTens,jm,h1mode(1,1,iOp,jm),nSPFmax,AR(1,im - 1),AR(1,im))
            else
               call Mult_GM_MatV(sizeTens,jm,'C',nS2,nS1,h1mode(1,1,iOp,jm),nSPFmax,AR(1,im - 1),AR(1,im))
               sizeTens(jm) = nS1
               nAC2 = nAC2 * nS1 / nS2
            end if
         end do
         !
         ! Get the left-multiplied vectors: A x H_1, A x H_1 x H_2, ...
         !
         sizeTens = nSPF(:,iel1)
         call zcopy(nAC1,A0(iA1 + 1), 1, AL(1,0), 1)
         call zlacgv(nAC1, AL(1,0), 1)
         do im = 1, nmodes - 1
            if (dvrType(iDOF(1,im)) .ne. 'gwp' .and. IsIdentity(im,iOperExtra(im,iT,iExtraOp))) then
               call zcopy(nAC1, AL(1,im - 1), 1, AL(1,im), 1)
               cycle
            end if
            ! Multiply by the 1-mode Hamiltonian matrix
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            iOp = iOperExtra(im,iT,iExtraOp)
            if (iel1 == iel2 .or. .not. MultiSet(im)) then
               call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,AL(1,im - 1),AL(1,im))
            else
               call Mult_GM_MatV(sizeTens,im,'N',nS1,nS2,h1mode(1,1,iOp,im),nSPFmax,AL(1,im - 1),AL(1,im))
               sizeTens(im) = nS2
               nAC1 = nAC1 * nS2 / nS1
            end if
         end do
         !
         nAC1 = nAConf(iel1)
         nAC2 = nAConf(iel2)
         ! Update Adot_tmp (the complex conjugate of -i * H * A) and the Lagrange multipliers
         cdummy = der_factor * ciOne * aet
         call zaxpy(nAC1, conjg(cdummy), AR(1,nModes), 1, Adot_tmp(iA1 + 1),1)
         dummy = - real(conjg(cdummy) * zdotu(nAC1, AR(1,nModes),1, A0(iA1 + 1),1) )
         multLagDiag = multLagDiag + dummy
         do iW2 = 1, iW - 1
            iB2 = APacketOffset(iW2) + AStateOffset(iel1) + 1
            multLag(iW,iW2) = multLag(iW,iW2) &
                            + conjg(cdummy) * zdotu(nAC1, AR(1,nModes),1, A0(iB2),1)
         end do
         do iW2 = iW + 1, npackets
            iB2 = APacketOffset(iW2) + AStateOffset(iel1) + 1
            multLag(iW2,iW) = multLag(iW2,iW) &
                            + cdummy * conjg(zdotu(nAC1, AR(1,nModes),1, A0(iB2),1))
         end do
         ! Update the mean field
         do im = 1, nmodes
            if (IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle  ! terms associated with the identity are projected out
            jm = nmodes - im
            iOp = iOperExtra(im,iT,iExtraOp)
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            loop_length = product(nSPF(1:im - 1,iel2))
            nloops = product(nSPF(im + 1:nmodes,iel1))
            iB1size = nS1 * loop_length
            iB2size = nS2 * loop_length
            call zlacgv(loop_length * nS2 * nloops, AR(1,jm), 1)
            iB1 = 1
            iB2 = 1
            do j = 1, nloops
               call zgemm('T','N',nS1,nS2,loop_length,ciOne * aet,AL(iB1,im - 1),loop_length, &
                          AR(iB2,jm),loop_length,cOne,MeanField(1,1,iOp,im),nSPFmax)
               iB1 = iB1 + iB1size
               iB2 = iB2 + iB2size
            end do
         end do
         !
      end do   ! Close the loop on the wave packet index (iW)
   end do
end do
!
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Terms related to the perturbative harmonic bath
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
!     The matrices FQbath and FPbath are calculated here
!
jQ = iCQ
jP = iCP
do iDissOp = 1, nHObath
   ! Calculate the matrix FQ * A 
   iExtraOp = iHObathQ(iDissOp)
   VQA = cZero
   do iT = 1, nOpTerms(iExtraOp)
      aet = conjg(aeExtra(iT,iExtraOp))
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      do iW = 1, npackets
         iA2 = APacketOffset(iW) + AStateOffset(iel2)
         call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1), 1)
         call zlacgv(nAC2,Aux2(1),1)
         pA1 => Aux1
         pA2 => Aux2
         sizeTens = nSPF(:,iel2)
         do im = 1, nmodes
            if (dvrType(iDOF(1,im)) .ne. 'gwp' .and. IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle    ! SkipTerm
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            ! Multiply by the 1-mode-Hamiltonian matrix
            iOp = iOperExtra(im,iT,iExtraOp)
            if (iel1 .eq. iel2 .or. .not. MultiSet(im)) then  ! h1mode is assumed to be Hermitian
               call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            else
               call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            end if
            pswap => pA2
            pA2   => pA1
            pA1   => pswap
            sizeTens(im) = nS1
         end do
         iA1 = APacketOffset(iW) + AStateOffset(iel1)
         call zaxpy(nAC1,aet,pA2(1),1,VQA(iA1 + 1),1)
      end do
   end do
   ! Calculate the matrix FP * A 
   iExtraOp = iHObathP(iDissOp)
   VPA = cZero
   do iT = 1, nOpTerms(iExtraOp)
      aet = conjg(aeExtra(iT,iExtraOp))
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      do iW = 1, npackets
         iA2 = APacketOffset(iW) + AStateOffset(iel2)
         call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1), 1)
         call zlacgv(nAC2,Aux2(1),1)
         pA1 => Aux1
         pA2 => Aux2
         sizeTens = nSPF(:,iel2)
         do im = 1, nmodes
            if (dvrType(iDOF(1,im)) .ne. 'gwp' .and. IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle    ! SkipTerm
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            ! Multiply by the 1-mode-Hamiltonian matrix
            iOp = iOperExtra(im,iT,iExtraOp)
            if (iel1 .eq. iel2 .or. .not. MultiSet(im)) then  ! h1mode is assumed to be Hermitian
               call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            else
               call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
            end if
            pswap => pA2
            pA2   => pA1
            pA1   => pswap
            sizeTens(im) = nS1
         end do
         iA1 = APacketOffset(iW) + AStateOffset(iel1)
         call zaxpy(nAC1,aet,pA2(1),1,VPA(iA1 + 1),1)
      end do
   end do
   ! Build the representation matrices of the operators FQ and FP
   ! in the basis of the wave packets
   call zgemm('T','N',npackets,npackets,nAConfTot,cOne, &
              VQA,nAConftot,A0,nAConfTot,cZero,FQbath(1,1,iDissOp),npackets)
   !call zgemm('T','N',npackets,npackets,nAConfTot,cOne, &
   !           A0,nAConftot,VQA,nAConfTot,cZero,FQbath(1,1,iDissOp),npackets)
   !call zlacgv(npackets * npackets,FQbath(1,1,iDissOp),1)
   !
   call zgemm('T','N',npackets,npackets,nAConfTot,cOne, &
              VPA,nAConftot,A0,nAConfTot,cZero,FPbath(1,1,iDissOp),npackets)
   !call zgemm('T','N',npackets,npackets,nAConfTot,cOne, &
   !           A0,nAConftot,VPA,nAConfTot,cZero,FPbath(1,1,iDissOp),npackets)
   !call zlacgv(npackets * npackets,FPbath(1,1,iDissOp),1)
   ! Symmetrize
   do iW = 1, npackets
      do iW2 = 1, iW - 1
         cdummy = cOneHalf * (FQbath(iW2,iW,iDissOp) + conjg(FQbath(iW,iW2,iDissOp)))
         FQbath(iW2,iW,iDissOp) = cdummy
         FQbath(iW,iW2,iDissOp) = conjg(cdummy)
         cdummy = cOneHalf * (FPbath(iW2,iW,iDissOp) + conjg(FPbath(iW,iW2,iDissOp)))
         FPbath(iW2,iW,iDissOp) = cdummy
         FPbath(iW,iW2,iDissOp) = conjg(cdummy)
      end do
      FQbath(iW,iW,iDissOp) = complex(real(FQbath(iW,iW,iDissOp)),dZero)
      FPbath(iW,iW,iDissOp) = complex(real(FPbath(iW,iW,iDissOp)),dZero)
   end do
   ! cos(w * t) and sin(w * t)
   cwt = cos(wBath(iDissOp) * tt)
   swt = sin(wBath(iDissOp) * tt)
   ! Update the Lagrange multiplier matrix
   cFQsFP = complex(dZero,-cwt) * FQbath(:,:,iDissOp) &
          + complex(dZero, swt) * FPbath(:,:,iDissOp)
   sFQcFP = complex(dZero,-swt) * FQbath(:,:,iDissOp) &
          + complex(dZero,-cwt) * FPbath(:,:,iDissOp)
   call zher2k('L','N',npackets,npackets,cOne, &
               cFQsFP,npackets,psi0(jQ),npackets,dOne,multLag,npackets)
   call zher2k('L','N',npackets,npackets,cOne, &
               sFQcFP,npackets,psi0(jP),npackets,dOne,multLag,npackets)
   !call zher2k('L','N',npackets,npackets,complex(dZero,-cwt), &
   !        FQbath(1,1,iDissOp),npackets,psi0(jQ),npackets,dOne,multLag,npackets)
   !call zher2k('L','N',npackets,npackets,complex(dZero,-swt), &
   !        FQbath(1,1,iDissOp),npackets,psi0(jP),npackets,dOne,multLag,npackets)
   !call zher2k('L','N',npackets,npackets,complex(dZero, swt), &
   !        FPbath(1,1,iDissOp),npackets,psi0(jQ),npackets,dOne,multLag,npackets)
   !call zher2k('L','N',npackets,npackets,complex(dZero,-cwt), &
   !        FPbath(1,1,iDissOp),npackets,psi0(jP),npackets,dOne,multLag,npackets)
   !
   ! Update Adot_tmp  (the complex conjugate of S * dA/dt)
   call zgemm('N','T',nAConfTot,npackets,npackets, complex(dZero, cwt), &
              VQA,nAConfTot,psi0(jQ),npackets,cOne,Adot_tmp,nAConfTot)
   call zgemm('N','T',nAConfTot,npackets,npackets, complex(dZero, swt), &
              VQA,nAConfTot,psi0(jP),npackets,cOne,Adot_tmp,nAConfTot)
   call zgemm('N','T',nAConfTot,npackets,npackets, complex(dZero,-swt), &
              VPA,nAConfTot,psi0(jQ),npackets,cOne,Adot_tmp,nAConfTot)
   call zgemm('N','T',nAConfTot,npackets,npackets, complex(dZero, cwt), &
              VPA,nAConfTot,psi0(jP),npackets,cOne,Adot_tmp,nAConfTot)
   !
   ! Evaluate the mean field
   !
   ! coordinate coupling
   iExtraOp = iHObathQ(iDissOp)
   call zhemm('R','U',nAConfTot,npackets,complex(cwt,dZero), &
              psi(jQ),npackets, A0,nAConfTot, cZero,AuxA,nAConfTot)
   call zhemm('R','U',nAConfTot,npackets,complex(swt,dZero), &
              psi(jP),npackets, A0,nAConfTot,  cOne,AuxA,nAConfTot)
   !call zgemm('N','N',nAConfTot,npackets,npackets,complex(cos(wBath(iDissOp) * tt),dZero), &
   !           A0,nAConftot, psi(jQ),npackets, cZero,AuxA,nAConfTot)
   !call zgemm('N','N',nAConfTot,npackets,npackets,complex(sin(wBath(iDissOp) * tt),dZero), &
   !           A0,nAConftot, psi(jP),npackets, cOne ,AuxA,nAConfTot)
   !
   do iT = 1, nOpTerms(iExtraOp)
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      do im = 1, nmodes
         if (IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle
         loop_length = product(nSPF(1:im - 1,iel1))
         nloops = product(nSPF(im + 1:nmodes,iel1))
         iB1size = nSPF(im,iel1) * loop_length
         iB2size = nSPF(im,iel2) * loop_length
         do iW = 1, npackets
            iA2 = APacketOffset(iW) + AStateOffset(iel2)
            call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1),1)
            call zlacgv(nAC2,Aux2(1),1)
            pA1 => Aux1
            pA2 => Aux2
            sizeTens = nSPF(:,iel2)
            do im2 = 1, nmodes
               if (im2 == im) cycle    ! do not operate on the MeanField mode or with the identity matrix
               if (dvrType(iDOF(1,im2)) .ne. 'gwp' .and. IsIdentity(im2,iOperExtra(im2,iT,iExtraOp))) cycle   ! Skip term
               nS1 = nSPF(im2,iel1)
               nS2 = nSPF(im2,iel2)
               ! Multiply by the 1-mode-Hamiltonian matrix
               iOp = iOperExtra(im2,iT,iExtraOp)
               if (iel1 == iel2 .or. .not. MultiSet(im2)) then   ! Assume Hermitian h1mode
                  call Mult_HeM_MatV(sizeTens,im2,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               else
                  call Mult_GM_MatV(sizeTens,im2,'C',nS2,nS1,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               end if
               pswap => pA2
               pA2   => pA1
               pA1   => pswap
               sizeTens(im2) = nS1
            end do
            ! At this stage <pA2| = < psi_iW^(im) | V^+ , where psi_iW^(im) is a single-hole function
            ! Update the mean-field matrix
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            iOp = iOperExtra(im,iT,iExtraOp)
            !
            iB1 = APacketOffset(iW) + AStateOffset(iel1) + 1
            iB2 = 1
            do j = 1, nloops
               call zgemm('T','N',nS2,nS1,loop_length, aeExtra(iT,iExtraOp), pA2(iB2),loop_length, &
                          AuxA(iB1), loop_length,cOne, MeanField(1,1,iOp,im), nSPFmax)
               iB1 = iB1 + iB1size
               iB2 = iB2 + iB2size
            end do
         end do
      end do
   end do
   ! momentum coupling
   iExtraOp = iHObathP(iDissOp)
   call zhemm('R','U',nAConfTot,npackets,-complex(swt,dZero), &
              psi(jQ),npackets, A0,nAConfTot, cZero,AuxA,nAConfTot)
   call zhemm('R','U',nAConfTot,npackets, complex(cwt,dZero), &
              psi(jP),npackets, A0,nAConfTot,  cOne,AuxA,nAConfTot)
   !call zgemm('N','N',nAConfTot,npackets,npackets,complex(-sin(wBath(iDissOp) * tt),dZero), &
   !           A0,nAConftot, psi(jQ),npackets, cZero,AuxA,nAConfTot)
   !call zgemm('N','N',nAConfTot,npackets,npackets,complex( cos(wBath(iDissOp) * tt),dZero), &
   !           A0,nAConftot, psi(jP),npackets, cOne ,AuxA,nAConfTot)
   !
   do iT = 1, nOpTerms(iExtraOp)
      iel1 = lhs_elExtra(iT,iExtraOp)
      iel2 = rhs_elExtra(iT,iExtraOp)
      nAC1 = nAConf(iel1)
      nAC2 = nAConf(iel2)
      do im = 1, nmodes
         if (IsIdentity(im,iOperExtra(im,iT,iExtraOp))) cycle
         loop_length = product(nSPF(1:im - 1,iel1))
         nloops = product(nSPF(im + 1:nmodes,iel1))
         iB1size = nSPF(im,iel1) * loop_length
         iB2size = nSPF(im,iel2) * loop_length
         do iW = 1, npackets
            iA2 = APacketOffset(iW) + AStateOffset(iel2)
            call zcopy(nAC2, A0(iA2 + 1),1, Aux2(1),1)
            call zlacgv(nAC2,Aux2(1),1)
            pA1 => Aux1
            pA2 => Aux2
            sizeTens = nSPF(:,iel2)
            do im2 = 1, nmodes
               if (im2 == im) cycle    ! do not operate on the MeanField mode or with the identity matrix
               if (dvrType(iDOF(1,im2)) .ne. 'gwp' .and. IsIdentity(im2,iOperExtra(im2,iT,iExtraOp))) cycle   ! Skip term
               nS1 = nSPF(im2,iel1)
               nS2 = nSPF(im2,iel2)
               ! Multiply by the 1-mode-Hamiltonian matrix
               iOp = iOperExtra(im2,iT,iExtraOp)
               if (iel1 == iel2 .or. .not. MultiSet(im2)) then   ! Assume Hermitian h1mode
                  call Mult_HeM_MatV(sizeTens,im2,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               else
                  call Mult_GM_MatV(sizeTens,im2,'C',nS2,nS1,h1mode(1,1,iOp,im2),nSPFmax,pA2,pA1)
               end if
               pswap => pA2
               pA2   => pA1
               pA1   => pswap
               sizeTens(im2) = nS1
            end do
            ! At this stage <pA2| = < psi_iW^(im) | V^+ , where psi_iW^(im) is a single-hole function
            ! Update the mean-field matrix
            nS1 = nSPF(im,iel1)
            nS2 = nSPF(im,iel2)
            iOp = iOperExtra(im,iT,iExtraOp)
            !
            iB1 = APacketOffset(iW) + AStateOffset(iel1) + 1
            iB2 = 1
            do j = 1, nloops
               call zgemm('T','N',nS2,nS1,loop_length, aeExtra(iT,iExtraOp), pA2(iB2),loop_length, &
                          AuxA(iB1), loop_length,cOne, MeanField(1,1,iOp,im), nSPFmax)
               iB1 = iB1 + iB1size
               iB2 = iB2 + iB2size
            end do
         end do
      end do
   end do
   !
   !
   jQ = jQ + npackets**2
   jP = jP + npackets**2
end do
!!!!
!!!!
!!!!
do j = 1, npackets
   do i = j + 1, npackets
      dW = WPpop(i) - WPpop(j)
      dummy = eps_pop * exp(- (dW / eps_pop) ** 2)
      multLag(i,j) = multLag(i,j) * complex(dZero, dw / (dw ** 2 + dummy ** 2))
      !multLag(i,j) = ciOne * multLag(i,j) * dw / (dw ** 2 + eps_pop ** 2)
   enddo
enddo
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Suspend timer
call suspend_timer(3)
!
return
end subroutine Calc_MeanField



subroutine Calc_rhom1MF
! It calculates the product matrix rho^-1 * MeanField
! for the non-GWP modes
! - Before matrix inversion, the density matrix is regularised
use sysparam
use globals
use psidef
use hamilton
use storage
implicit none
integer :: i, iel, iHam
integer :: im, nS, nS2
double complex, dimension(:,:) :: rhoinv(nSPFmax,nSPFmax)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do i = 1, nmodDVR
   im = imodDVR(i)
   if (MultiSet(im)) then
      do iel = 1, nstates
         nS = nSPF(im,iel)
         ! Regularisation
         !call DirectRegInv(nS,rho(1,1,im,iel),nSPFmax,rhoinv,nSPFmax,eps_rho)
         !call DirectTikhInv(nS,rho(1,1,im,iel),nSPFmax,rhoinv,nSPFmax,eps_rho)
         call StdRegInv(nS,rho(1,1,im,iel),nSPFmax,rhoinv,nSPFmax,eps_rho,'rho ','poly')
            ! Multiply rhoinv by the MeanField
         do iHam = 1, nOper(im)
            if (Oper(-1,iHam,im) .ne. iel .or. IsIdentity(im,iHam)) cycle
            nS2 = nSPF(im,Oper(0,iHam,im))
            call zhemm('L','U',nS,nS2, cOne, rhoinv,nSPFmax, &
                       MeanField(1,1,iHam,im),nSPFmax, cZero, &
                       rhom1MF(1,1,iHam,i), nSPFmax)
         end do
         !
      end do
   else  
      ! Single set formalism
      nS = nSPF(im,1)
      call StdRegInv(nS,rho_ave(1,1,im),nSPFmax,rhoinv,nSPFmax,eps_rho,'rho ','poly')
      do iHam = 1, nOper(im)
         if (IsIdentity(im,iHam)) cycle
         call zhemm('L','U',nS,nS, cOne, rhoinv,nSPFmax, &
                    MeanField(1,1,iHam,im),nSPFmax, cZero, &
                    rhom1MF(1,1,iHam,i), nSPFmax)
      end do
   end if
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_rhom1MF



subroutine Calc_CY(psi0, i, DD0)
! It calculates the matrices C and Y for i-th GWP mode,
! The C matrix is not regularised here
use sysparam
use globals
use psidef
use hamilton
use storage
use timingmod
implicit none
!
integer, intent(in) :: i
double complex, dimension(:), intent(in) :: psi0(dimpsi)
!
integer :: iel, im, nS, nMD, j1, j2, k, iP, iP1, iP2, &
           iC1, k1, k2, iel2, iHam, nS2
double complex :: cdummy
double complex, dimension(:) :: gInt(nDOF)
double complex, dimension(:,:) :: Sa0(nCYmax,nS00max), Sab(nCYmax,nCYmax), Ha0(nCYmax,nS00max)
double complex, dimension(:,:,:) :: gIntx(nDOF,nS00max,nS00max)
double complex, dimension(:,:), pointer :: Sa0S00m1, Sa0D0
double complex, dimension(nS00max,nS00max,nmodGWP,nstates), optional :: DD0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Start timer
call continue_timer(4)
!
im = imodGWP(i)
nMD = nmodeDOF(im)
do iel = 1, merge(nstates, 1, MultiSet(im))
   nS = nSPF(im,iel)
   YVec(1:nS * nMD,iel,i) = cZero   ! Important initialization
   ! Calculate the matrix elements of the x operator between Gaussians
   iP = ipsiSPF(im,iel) - 1
   do concurrent (j2 = 1:nS, j1 = 1:nS, k = 1:nMD, j1 .le. j2) 
      gIntx(k,j1,j2) = GauInt1D(GWPa(k,j1,i),GWPa(k,j2,i), &
                                psi0(iP + (j1 - 1) * nMD + k),psi0(iP + (j2 - 1) * nMD + k),1)
   end do
   ! Calculate the matrix Sa0 (for the GWP-mode i and the state iel)
   do concurrent (j2 = 1:nS, j1 = 1:nS, k = 1:nMD, j1 .le. j2) 
      Sa0((j1 - 1) * nMD + k,j2) = S00(j1,j2,i,iel,iel) * gIntx(k,j1,j2)
   end do
   do concurrent (j2 = 1:nS, j1 = 1:nS, k = 1:nMD, j1 .gt. j2) 
      Sa0((j1 - 1) * nMD + k,j2) = conjg(Sa0((j2 - 1) * nMD + k,j1))
   end do
   !
   Sa0S00m1 => Sa0Sm1(:,:,i,iel)   ! Sa0Sm1 is used in the calculation of derivatives
   if (.not. present(DD0)) then
      call zhemm('R','U',nS * nMD,nS,cOne,S00m1(1,1,i,iel),nS00max,Sa0,nCYmax,cZero,Sa0S00m1,nCYmax)
   else 
      Sa0D0 => Sa0D(:,:,i,iel)        ! Sa0D is used to calculate the derivatives of the D matrix
      call zgemm('N','N',nS * nMD,nS,nS,cOne,Sa0,nCYmax,DD0(1,1,i,iel),nS00max,cZero,Sa0D0,nCYmax)
      call zgemm('N','C',nS * nMD,nS,nS,cOne,Sa0D0,nCYmax,DD0(1,1,i,iel),nS00max,cZero,Sa0S00m1,nCYmax)
   end if
   !
   ! Calculate the Sab matrix
      ! k1 .ne. k2
   do concurrent (j2 = 1:nS, k2 = 1:nMD, j1 = 1:nS, k1 = 1:nMD, j1 .le. j2) 
      Sab((j1 - 1) * nMD + k1,(j2 - 1) * nMD + k2) = S00(j1,j2,i,iel,iel) &
         * merge(GauInt1D(GWPa(k1,j1,i),GWPa(k2,j2,i),psi0(iP + (j1 - 1) * nMD + k1),psi0(iP + (j2 - 1) * nMD + k2),2), &
                 gIntx(k1,j1,j2) * gIntx(k2,j1,j2), k1 == k2)
   end do
      ! k1 .eq. k2
   !do concurrent (j2 = 1:nS, j1 = 1:nS, k = 1:nMD, j1 .le. j2) 
   !   Sab((j1 - 1) * nMD + k,(j2 - 1) * nMD + k) = S00(j1,j2,i,iel,iel) &
   !      * GauInt1D(GWPa(k,j1,i),GWPa(k,j2,i),psi0(iP + (j1 - 1) * nMD + k),psi0(iP + (j2 - 1) * nMD + k),2)
   !end do
   !
   if (.not. present(DD0)) then
      call zher2k('U','N',nS * nMD,nS,-cOneHalf,Sa0S00m1,nCYmax,Sa0,nCYmax,dOne,Sab,nCYmax)
   else
      call zherk('U','N',nS * nMD,nS,-dOne,Sa0D0,nCYmax,dOne,Sab,nCYmax)
   end if
   !
   ! Now the C-matrix can be calculated
   !    (only the upper diagonal)
   if (MultiSet(im)) then
      do concurrent (j2 = 1:nS, k2 = 1:nMD, j1 = 1:nS, k1 = 1:nMD, j1 .le. j2) 
         CMat((j1 - 1) * nMD + k1,(j2 - 1) * nMD + k2,iel,i) = &
            rho(j1,j2,im,iel) * Sab((j1 - 1) * nMD + k1,(j2 - 1) * nMD + k2)
      end do
   else
      do concurrent (j2 = 1:nS, k2 = 1:nMD, j1 = 1:nS, k1 = 1:nMD, j1 .le. j2)
         CMat((j1 - 1) * nMD + k1,(j2 - 1) * nMD + k2,iel,i) = &
            rho_ave(j1,j2,im) * Sab((j1 - 1) * nMD + k1,(j2 - 1) * nMD + k2)
      end do
   end if
   ! The calculation of the Y-vector starts here
   do iHam = 1, nOper(im)
      if (IsIdentity(im,iHam)) cycle
      if (Oper(-1,iHam,im) .ne. iel .and. MultiSet(im)) cycle   
      iel2 = Oper(0,iHam,im)
      nS2 = nSPF(im,iel2)
      ! Calculate the matrix Ha0
      iP2 = ipsiSPF(im,iel2) - 1
      do j2 = 1, nS2
         iC1 = 0
         iP1 = ipsiSPF(im,iel) - 1
         do j1 = 1, nS
            cdummy = S00(j1,j2,i,iel,iel2)
            do concurrent (k = 1:nMD)
               gInt(k) = GauInt1D(GWPa(k,j1,i),GWPa(k,j2,i), &
                                 psi0(iP1 + k),psi0(iP2 + k),Oper(k,iHam,im))
            end do
            do concurrent (k = 1:nMD) 
               Ha0(iC1 + k,j2) = cdummy * product(gInt(1:k - 1)) * product(gInt(k + 1:nMD)) &
                                    * Gau10Int1D(GWPa(k,j1,i),GWPa(k,j2,i), &
                                                 psi0(iP1 + k),psi0(iP2 + k),Oper(k,iHam,im))
            end do
            iC1 = iC1 + nMD
            iP1 = iP1 + nMD
         end do
         iP2 = iP2 + nMD
      end do
      if (iel2 .eq. iel) then  ! h1mode is assumed Hermitian
         call zhemm('R','U',nS * nMD,nS, -cOne,h1mode(1,1,iHam,im),nSPFmax, &
                    Sa0S00m1(1,1),nCYmax, cOne,Ha0(1,1),nCYMax)
      else
         call zgemm('N','N',nS * nMD,nS2,nS, -cOne,Sa0S00m1(1,1),nCYmax,h1mode(1,1,iHam,im),nSPFmax, &
                 cOne, Ha0(1,1),nCYMax)
      end if
      ! Update Yvec
      do j2 = 1, nS2
         iC1 = 1
         do j1 = 1, nS
            call zaxpy(nMD, MeanField(j1,j2,iHam,im), Ha0(iC1,j2), 1, Yvec(iC1,iel,i), 1)
            iC1 = iC1 + nMD
         end do
      end do
   end do
!
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Suspend timer
call suspend_timer(4)
return
end subroutine Calc_CY



subroutine Calc_PrecC(psi0)
! Calculate the preconditioner useful to invert the C matrix
use globals, only : eps_C
use psidef
use storage
implicit none
!
double complex, dimension(:), intent(in) :: psi0(dimpsi)
integer :: i, im, nMD, nS, iel, nG
!
interface
   subroutine Calc_CY(psi0, i, DD0)     ! This interface is necessary because the argument DD0 is optional
   use psidef
   use storage, only : nS00max
   double complex, dimension(:), intent(in) :: psi0(dimpsi)
   integer, intent(in) :: i
   double complex, dimension(nS00max,nS00max,nmodGWP,nstates), optional :: DD0
   end subroutine Calc_CY
end interface
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do i = 1, nmodGWP
   call Calc_CY(psi0,i)
   im = imodGWP(i)
   nMD = nmodeDOF(im)
   do iel = 1, merge(nstates, 1, MultiSet(im))
      nS = nSPF(im,iel)
      nG = nS * nMD
      call InverseSquareRoot(nG, CMat(1,1,iel,i), nCYmax, PrecC(1,1,iel,i), nCYmax, eps_C)
   end do
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Calc_PrecC
 



subroutine UpdateTDcoeff(tt)
! It calculates the time-dependent coefficients at the time tt
use sysparam
use globals
use hamilton
double precision :: tt
integer :: i, iHam
double precision :: c0, T0, deltaT, omega, phase
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do i = 1, nHamTD
   iHam = iHamTD(i)
   select case(typeHamTD(i))
      case(1)   ! Gaussian
         c0 = parHamTD(0,i)
         T0 = parHamTD(1,i)
         deltaT = parHamTD(2,i)
         omega = parHamTD(3,i)
         phase = parHamTD(4,i)
         ae(iHam) = c0 * exp(- dOneHalf * ((tt - T0) / deltaT) ** 2 &
                             + ciOne * omega * (tt - T0) + ciOne * phase)
      case default
         stop
   end select
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine UpdateTDCoeff


!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!



subroutine Mult_HeM_MatV(sizeTens,im,HM,LDHM,A1,A2)
! Multiplication between an Hermitian matrix HM
! and the matricisation of the vector A1 along the mode im.
! The vectorized result is written into A2
use sysparam
use psidef
implicit none
integer, intent(in) :: im, LDHM
integer, dimension(*), intent(in) :: sizeTens  ! nTens (size of each tensor dimension)
integer :: iloop, loop_length, nloops, nS, iAsize, iA
double complex, dimension(*), intent(in) :: A1
double complex, dimension(LDHM,*), intent(in) :: HM  ! only the upper diagonal needs to be referenced
double complex, dimension(*), intent(out) :: A2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
loop_length = product(sizeTens(1:im - 1))
nS = sizeTens(im)
nloops = product(sizeTens(im + 1:nmodes))
! Multiply by the matrix HM
iA = 1
iAsize = nS * loop_length
do iloop = 1, nloops
   call zhemm('R','U',loop_length,nS,cOne,HM,LDHM, &
              A1(iA),loop_length, cZero,A2(iA),loop_length)
   iA = iA + iAsize
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Mult_HeM_MatV



subroutine Mult_GM_MatV(sizeTens,im,transG,n1,n2,GM,LDGM,A1,A2)
! Multiplication between the matricisation of the vector A1 along the mode im.
! and the Hermitian conjugate of the n2xn1 matrix GM
! The vectorized result is written into A2
use sysparam
use psidef
implicit none
integer, intent(in) :: im, n1, n2, LDGM
integer, dimension(*), intent(in) :: sizeTens  ! nTens (size of each tensor dimension)
character(len = 1), intent(in) :: transG   ! 'C' or 'N' according to Lapack convention
integer :: iloop, loop_length, nloops, iA1, iA2, iA1size, iA2size
double complex, dimension(*), intent(in) :: A1
double complex, dimension(LDGM,*), intent(in) :: GM  
double complex, dimension(*), intent(out) :: A2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
loop_length = product(sizeTens(1:im - 1))
nloops = product(sizeTens(im + 1:nmodes))
! Multiply by the matrix GM
iA1 = 1
iA2 = 1
iA1size = n1 * loop_length
iA2size = n2 * loop_length
do iloop = 1, nloops
   call zgemm('N',transG,loop_length,n2,n1,cOne,A1(iA1),loop_length,GM,LDGM,cZero,A2(iA2),loop_length)
   !call zgemm('N','C',loop_length,n2,n1,cOne,A1(iA1),loop_length,GM,LDGM,cZero,A2(iA2),loop_length)
   iA1 = iA1 + iA1size
   iA2 = iA2 + iA2size
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Mult_GM_MatV



!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Functions and subroutines related to Gaussian integrals
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!

pure double complex function GauOvl1D(A1,A2,xi1,xi2)
! It yields the overlap <g_1|g_2>, where g_i is a normalised Gaussian function:
! g_i = N(A_i,xi_i) * exp(A_i * x^2 + xi_i * x)
! WARNING: In this function, it is assumed that A_i < 0
implicit none
double precision, intent(in) :: A1, A2
double complex, intent(in) :: xi1, xi2
double precision :: dummy
double complex :: cdummy
double precision, parameter :: OneFourth = 0.25d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
cdummy = real(xi1) ** 2 / A1 + real(xi2) ** 2 / A2 &
         - (conjg(xi1) + xi2) ** 2 / (A1 + A2)
dummy = sqrt(4 * A1 * A2 / (A1 + A2) ** 2)
GauOvl1D = exp(OneFourth * cdummy) * sqrt(dummy)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end function GauOvl1D



pure double complex function GauOvlMultiD(A1,A2,xi1,xi2,nmD)
! It yields  the overlap <g_1|g_2>, where g_i is a normalised nmD-dimensional Gaussian function
! g_i = N(A_i,xi_i) * exp(A_i * x^2 + xi_i * x)
use sysparam
implicit none
integer, intent(in) :: nmD
double precision, dimension(nmD), intent(in) :: A1, A2
double complex, dimension(nmD), intent(in) :: xi1, xi2
double precision :: dummy
double complex :: cdummy
double precision, parameter :: OneFourth = 0.25d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
if (all(abs(A1 - A2) < epsilon(dOne))) then
   cdummy = sum((real(xi1) ** 2 + real(xi2) ** 2 - cOneHalf * (conjg(xi1) + xi2) ** 2) / A1)
   GauOvlMultiD = exp(OneFourth * cdummy)
   return
end if        
cdummy = sum(real(xi1) ** 2 / A1 + real(xi2) ** 2 / A2 - (conjg(xi1) + xi2) ** 2 / (A1 + A2))
dummy = sqrt(product(4 * A1 * A2 / (A1 + A2) ** 2))
GauOvlMultiD = exp(OneFourth * cdummy) * sqrt(dummy)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end function GauOvlMultiD
