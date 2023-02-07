subroutine Calc_DHD_B(B0, HB0)
! Calculate the product between the (transformed) Hamiltonian D^H * H * D 
! and the coefficient vector
use sysparam
use globals
use hamilton
use psidef
use storage, only : nS00max
implicit none
double complex, dimension(:), intent(in)  :: B0(nAconftot)
double complex, dimension(:), intent(out) :: HB0(nAconftot)
integer :: iel1, iel2, iA1, iA2, nAC1, nAC2, nS1, nS2, iOp, im, i, iHam
integer, dimension(nmodes) :: sizeTens
double complex, dimension(:) :: Baux1(nAconftot), Baux2(nAconftot), B0cgv(nAconftot)
double complex, dimension(:), pointer :: pA1, pA2, pswap
target :: Baux1, Baux2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
!HB0 = cZero   ! this initialisation is important
call zlaset('All',nAconftot,1,cZero,cZero,HB0(1),nAconftot)
call zcopy(nAconftot, B0(1), 1, B0cgv(1), 1)
call zlacgv(nAconftot, B0cgv(1), 1)
!
do iHam = 1, nHam
   iel1 = lhs_el(iHam)
   iel2 = rhs_el(iHam)
   iA1 = AStateOffset(iel1)
   iA2 = AStateOffset(iel2)
   nAC1 = nAConf(iel1)
   nAC2 = nAConf(iel2)
   call zcopy(nAC2, B0cgv(iA2 + 1), 1, Baux2(1), 1)
   !call zcopy(nAC2, B0(iA2 + 1), 1, Baux2(1), 1)
   !call zlacgv(nAC2, Baux2(1), 1)
   pA1 => Baux1
   pA2 => Baux2
   ! Loop separately over the DVR and GWP modes
   ! because the DVR modes use the matrix h1mode
   ! whereas the GWP modes use Dh1modeD
   sizeTens = nSPF(:,iel2)
   ! DVR modes
   do i = 1, nmodDVR
      im = imodDVR(i)
      if (IsIdentity(im,iOper(im,iHam))) cycle
      nS1 = nSPF(im,iel1)
      nS2 = nSPF(im,iel2)
      iOp = iOper(im,iHam)
      if (iel1 == iel2 .or. .not. MultiSet(im)) then
         call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
      else
         call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,h1mode(1,1,iOp,im),nSPFmax,pA2,pA1)
      end if
      pswap => pA2
      pA2 => pA1
      pA1 => pswap
      sizeTens(im) = nS1
   end do
   ! GWP modes
   do i = 1, nmodGWP
      im = imodGWP(i)
      if (IsIdentity(im,iOper(im,iHam))) cycle    ! because D^H * S * D = 1
      nS1 = nSPF(im,iel1)
      nS2 = nSPF(im,iel2)
      iOp = iOper(im,iHam)
      if (iel1 == iel2 .or. .not. MultiSet(im)) then
         call Mult_HeM_MatV(sizeTens,im,Dh1modeD(1,1,iOp,i),nS00max,pA2,pA1)
      else
         call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,Dh1modeD(1,1,iOp,i),nS00max,pA2,pA1)
      end if
      pswap => pA2
      pA2 => pA1
      pA1 => pswap
      sizeTens(im) = nS1
   end do
   ! Update HB0
   call zaxpy(nAC1,conjg(ae(iHam)),pA2,1,HB0(iA1 + 1),1)
end do
! Up to now the complex conjugate of H * B has been calculated
call zlacgv(nAconftot, HB0, 1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
return      
end subroutine Calc_DHD_B



subroutine Calc_Der_SPF(psi0,DD0,dpsi,dDD,i,typeinfo)
use sysparam
use globals
use hamilton
use psidef
use storage
use dvr
use timingmod
double complex, dimension(dimpsi), intent(in)  :: psi0
double complex, dimension(dimpsi), intent(out) :: dpsi
double complex, dimension(nS00max,nS00max,nmodGWP,nstates), intent(in)  :: DD0
double complex, dimension(nS00max,nS00max,nmodGWP,nstates), intent(out) :: dDD
integer, intent(in) :: i
character(len = 3), intent(in) :: typeinfo
integer :: k, im, nG, nMD, iHam, iel, iel1, iel2, nS, nS1, nS2, &
           iP1, iP2, iOD, kdof, iD, info, kg, j2
! For the SPF projector and the derivatives of the D matrix
double complex, dimension(:,:) :: ovl(nSPFmax,nSPFmax), ovldot(nSPFmax,nSPFmax), &
                                  tau(nS00max,nS00max)
!
double precision :: dummy
double complex, dimension(nstates) :: cFac
double complex, external :: zdotc
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
if (typeinfo == 'DVR') then   ! SPFs on DVR grids
   call continue_timer(6)
   im = imodDVR(i)
   nG = ngpm(i)
   nMD = nmodeDOF(im)
   !
   cFac = cZero   ! This basically initialises the derivatives to zero
   do iHam = 1, nOper(im)
      if (IsIdentity(im,iHam)) cycle
      iel1 = Oper(-1,iHam,im)
      iel2 = Oper(0,iHam,im)
      nS2 = nSPF(im,iel2)
      iP2 = ipsiSPF(im,iel2)
      call zcopy(nS2 * nG,psi0(iP2),1,psi_sect,1)      ! copy the section of psi0
                                                       ! corresponding to the SPF of the mode im for the state iel2
      ! Operate with the operator iHam on the SPFs of mode im
      do k = 1, nMD
         iOD = Oper(k,iHam,im)
         if (iOD .eq. 0) cycle
         kdof = iDOF(k,im)
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
            case(-4)  ! sin^-2
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
      ! Multiply by the mean field matrix
      call zgemm('N','T',nG,nS1,nS2, der_factor, psi_sect(1),nG, rhom1MF(1,1,iHam,i), nSPFmax, cFac(iel1), dpsi(iP1), nG)
      cFac(iel1) = cOne   ! next iteration, add
   end do
   ! Project out the current SPFs
   do iel = 1, merge(nstates,1,MultiSet(im))   ! 1 if the mode is treated with a single-set of SPFs, otherwise nstates
      nS = nSPF(im,iel)
      iP1 = ipsiSPF(im,iel)
      call zherk('U','C',nS,nG,dOne,psi0(iP1),nG,dZero,ovl,nSPFmax)  ! overlap matrix
      call zgemm('C','N',nS,nS,nG,cOne, psi0(iP1),nG, dpsi(iP1),nG, cZero,ovldot,nSPFmax)
      call zposv('U',nS,nS,ovl,nSPFmax,ovldot,nSPFmax,info)
      call zgemm('N','N',nG,nS,nS, -cOne, psi0(iP1),nG, ovldot,nSPFmax, cOne, dpsi(iP1),nG)
   end do
   call suspend_timer(6)
!   
else  ! Gaussian wave packets
!  
   call continue_timer(7)
   !
   im = imodGWP(i)
   ! Calculate the matrix C and the vector Y
   call Calc_S00(psi0,i)
   call Calc_S00m1(i)
   !
   call suspend_timer(7)
   call Calc_h1mode(psi0,im)
   !call Check_DSD(nSPF(im,1), S00(1,1,i,1,1), DD0(1,1,i,1), nS00max)
   call Calc_CY(psi0, i)
   !call Calc_CY(psi0, i, DD0)
   call continue_timer(7)
   !
   do iel = 1, merge(nstates, 1, MultiSet(im))
      !
      ! Calculate the derivatives of the Gaussian parameters
      !
      nS = nSPF(im,iel)
      nMD = nmodeDOF(im)
      nG = nS * nMD
      iP1 = ipsiSPF(im,iel)
      !call DirectTikhSol(nG,CMat(1,1,iel,i),nCYmax,YVec(1,iel,i),dpsi(iP1),eps_C)
      !call CGPrecSol(nG,CMat(1,1,iel,i),nCYmax,YVec(1,iel,i),dpsi(iP1),PrecC(1,1,iel,i),nCYMax)
      !call CGSol(nG,CMat(1,1,iel,i),nCYmax,YVec(1,iel,i),dpsi(iP1),eps_C)
      call StdRegSol(nG,CMat(1,1,iel,i),nCYmax,YVec(1,iel,i),dpsi(iP1),eps_C,'tikh')
      call zscal(nG,der_factor,dpsi(iP1),1)
      !
      ! Calculate the derivative of the D matrix
      !
      kg = 1
      iP2 = ipsiSPF(im,iel)
      do j2 = 1, nS
         call zgemv('C',nMD,nS,cOne,Sa0Sm1(kg,1,i,iel),nCYmax,dpsi(iP2),1,cZero,tau(1,j2),1)
         !do j1 = 1, nS
         !   tau(j1,j2) = zdotc(nMD,Sa0Sm1(kg,j1,i,iel),1,dpsi(iP2),1)
         !end do
         kg = kg + nMD
         iP2 = iP2 + nMD
      end do
      ! tau contains the derivatives of the phases ...
      iP2 = ipsiSPF(im,iel) - 1
      do j2 = 1, nS
         dummy = dZero
         do k = 1, nMD
            dummy = dummy + dble(psi0(iP2 + k)) * dble(dpsi(iP2 + k)) / GWPa(k,j2,i)
         end do
         tau(j2,j2) = tau(j2,j2) + dummy * cOneHalf
         iP2 = iP2 + nMD
      end do
      !
      call zgemm('N','N',nS,nS,nS,-cOne,tau,nS00max,DD0(1,1,i,iel),nS00max,cZero,dDD(1,1,i,iel),nS00max)
   end do
   !
   call suspend_timer(7)
end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
return      
end subroutine Calc_Der_SPF



subroutine Calc_Err_SPF(psi1,psi2,i,typeinfo,error)
! Calculate the difference norm between two wavefunctions
! where only the SPFs of the mode i of type typeinfo ('DVR' or 'GWP') change
! error = Tr(rho * (S11 + S22 - S12 - S21)), where rho is the reduced density matrix
use sysparam
use globals
use psidef
use storage, only : rho, rho_ave
use dvr, only : ngpm
double complex, dimension(dimpsi), intent(in)  :: psi1, psi2
integer, intent(in) :: i
character(len = 3), intent(in) :: typeinfo
double precision, intent(out) :: error
integer :: im, nS, nG, iel, iP, j1, j2, nMD
double complex, dimension(dimpsi) :: deltapsi
double complex, dimension(nSPFmax,nSPFmax) :: ovl
!
interface
   pure double complex function GauOvlMultiD(A1,A2,xi1,xi2,nmD)
   integer, intent(in) :: nmD
   double precision, dimension(nmD), intent(in) :: A1, A2
   double complex, dimension(nmD), intent(in) :: xi1, xi2
   end function GauOvlMultiD
end interface
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
if (typeinfo == 'DVR') then   ! SPFs on DVR grids
   im = imodDVR(i)
   nG = ngpm(i)
   error = dZero
   do iel = 1, merge(nstates, 1, MultiSet(im))
      nS = nSPF(im,iel)
      iP = ipsiSPF(im,iel)  ! If the mode is single-set then ipsiSPF(im,iel) = ipsiSPF(im,1)
      call zcopy(nS * nG, psi1(iP), 1, deltapsi(iP), 1)
      call zaxpy(nS * nG, -cOne, psi2(iP), 1, deltapsi(iP), 1)
      call zherk('U', 'C', nS,nG,dOne, deltapsi(iP), nG, dZero, ovl,nSPFmax)
      !call zherk('U','C',nS,nG,dOne,psi1(iP),nG,dZero,ovl,nSPFmax)
      !call zherk('U','C',nS,nG,dOne,psi2(iP),nG,dOne, ovl,nSPFmax)
      !call zher2k('U','C',nS,nG,-cOne,psi1(iP),nG,psi2(iP),nG,dOne,ovl,nSPFmax)
      ! Trace
      if (MultiSet(im)) then
         do j2 = 1, nS
            do j1 = 1, j2 - 1
               error = error + 2 * dble(rho(j1,j2,im,iel) * conjg(ovl(j1,j2)))
            end do
            error = error + dble(rho(j2,j2,im,iel) * ovl(j2,j2))    !
         end do
      else ! single-set: use the density matrix traced over the electronic states
         do j2 = 1, nS
            do j1 = 1, j2 - 1
               error = error + 2 * dble(rho_ave(j1,j2,im) * conjg(ovl(j1,j2)))
               end do
            error = error + dble(rho_ave(j2,j2,im) * ovl(j2,j2))    
         end do
      end if
   end do
else     ! GWP
   im = imodGWP(i)
   nmD = nmodeDof(im)
   error = dZero
   do iel = 1, merge(nstates, 1, MultiSet(im))
      nS = nSPF(im,iel)
      iP = ipsiSPF(im,iel) - nmD
      do j1 = 1, nS
         ovl(j1,j1) = 2 * cOne &
                    - 2 * dble(GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j1,i), &
                                            psi1(iP + j1 * nmD),psi2(iP + j1 * nmD),nmD))
      end do
      do concurrent (j1 = 1:nS, j2 = 1:nS, j1 .lt. j2)
         ovl(j1,j2) = GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi1(iP + j1 * nmD),psi1(iP + j2 * nmD),nmD) &
                    + GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi2(iP + j1 * nmD),psi2(iP + j2 * nmD),nmD) &
                    - GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi1(iP + j1 * nmD),psi2(iP + j2 * nmD),nmD) &
                    - GauOvlMultiD(GWPa(1,j1,i),GWPa(1,j2,i), &
                                   psi2(iP + j1 * nmD),psi1(iP + j2 * nmD),nmD) 
      end do
      ! Trace
      if (MultiSet(im)) then
         do j2 = 1, nS
            do j1 = 1, j2 - 1
               error = error + 2 * dble(rho(j1,j2,im,iel) * conjg(ovl(j1,j2)))
            end do
            error = error + dble(rho(j2,j2,im,iel) * ovl(j2,j2))    ! For single-set it is better to sum the rho_s and then muliply
         end do
      else ! single-set: use the density matrix traced over the electronic states
         do j2 = 1, nS
            do j1 = 1, j2 - 1
               error = error + 2 * dble(rho_ave(j1,j2,im) * conjg(ovl(j1,j2)))
            end do
            error = error + dble(rho_ave(j2,j2,im) * ovl(j2,j2))    ! For single-set it is better to sum the rho_s and then muliply
         end do
      end if
   end do
end if
!
error = sqrt(abs(error))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
return      
end subroutine Calc_Err_SPF



subroutine Check_DSD(nn, S, D, LD)
! Check whether D^H * S * D = 1
use sysparam
integer, intent(in) :: nn, LD
double complex, dimension(LD,*), intent(in) :: S, D
double complex, dimension(:,:) :: SD(nn,nn), DSD(nn,nn)
integer :: i, j
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!  
write(2908,*)
write(2908,'(a,i2)') 'Check D and S matrices. nn = ', nn
call zhemm('L', 'U', nn, nn, cOne, S, LD, D, LD, cZero, SD, nn)  
call zher2k('U', 'C', nn, nn, cOneHalf, D, LD, SD, nn, cZero, DSD, nn)
! Set lower to have it Hermitian
do j = 1, nn
   do i = j + 1, nn
      DSD(i,j) = conjg(DSD(j,i))
   end do
end do
write(2908,*) 'D^H * S * D'
do i = 1, min(4,nn)
   write(2908,'(4(f8.5,1x,f8.5,3x))') DSD(i,1:min(4,nn))
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!           
return      
end subroutine Check_DSD
