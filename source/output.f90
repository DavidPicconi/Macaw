module output
   use sysparam
   use globals
   use psidef
   use storage
   use hamilton
   use dvr
   implicit none
   private
   public :: energyStart, Density, DensityPop, DensityMatrix
   public :: WriteOutput
   double precision :: energyStart
   double precision, dimension(:,:), allocatable :: Density, DensityPop
   double complex, dimension(:,:), allocatable :: DensityMatrix
contains

   subroutine WriteOutput(tt)
! Write information on the output files
   implicit none
   integer :: i, j, iel, jel, im, nS, iW, iA, &
             j1, j2, iP1, iP2, nG, ires_unit, iC
   double precision :: energy, tt, entropy, poptot
   double precision, dimension(:) :: pop(nstates)
   double complex :: auto, cdummy
   ! Variables for the calculation of natural populations
   integer :: info, LWORK
   double precision, dimension(:) :: Wv(nSPFmax), RWORK(3 * nSPFmax - 2)
   double precision, dimension(:,:,:) :: natpop(nSPFmax,nmodes,nstates), &
                                         GGP(nS00max,nmodGWP,nstates), PopErr(nS00max,nmodGWP,nstates)
   double complex, dimension(:) :: WORK(nSPFmax * nSPFmax)
   double complex, dimension(:,:) :: rhot(nSPFmax,nSPFmax), S12(nS00max,nS00max), &
                                     Srho(nS00max,nS00max), SWP(npackets,npackets)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Calculate all objects
   do i = 1, nmodGWP
      call Calc_S00(psi,i)
   end do
   do im = 1, nmodes
      call Calc_rho1mode(im,AVector)
      call Calc_h1mode(psi,im)
   end do
   call Calc_MeanField(AVector,psi,tt)   ! Used to evaluate the energy
! Calculate here the preconditioner for the C matrix (used only for CG iterations -- experimental)
!   call Calc_PrecC(psi)
! Calculate the entropy
   poptot = dZero
   iA = 1
   do iW = 1, npackets
      call norm_wf(AVector(iA), psi, WPpop(iW))
      poptot = poptot + WPpop(iW)
      iA = iA + nAConftot
   end do
   WPpop = WPpop / poptot
   entropy = dZero
   do iW = 1, npackets
      if (WPpop(iW) < 1d-8) cycle
      entropy = entropy - WPpop(iW) * log(WPpop(iW))
   end do
   entropy = entropy / log(dble(npackets))   ! respace to the maximum entropy
   call dlasrt('D', npackets, WPpop, info)
! Calculate electronic state populations
   ! The density matrix of mode 1 is used
   select case(dvrtype(iDOF(1,1)))
      case('gwp')   
         do iel = 1, nstates
            jel = merge(iel, 1, MultiSet(1))
            nS = nSPF(1,jel)
            pop(iel) = dZero
            do j2 = 1, nS
               do j1 = 1, j2 - 1
                  pop(iel) = pop(iel) + 2.d0 * real(rho(j1,j2,1,iel) * S00(j1,j2,1,jel,jel))
               end do
               pop(iel) = pop(iel) + real(rho(j2,j2,1,iel) * S00(j2,j2,1,jel,jel))
            end do
         end do
      case default
         do iel = 1, nstates
            jel = merge(iel, 1, MultiSet(1))
            nS = nSPF(1,jel)
            pop(iel) = dZero
            do i = 1, nS
               pop(iel) = pop(iel) + real(rho(i,i,1,iel))
            end do
         end do
   end select
   ! Energy
   energy = energyMF / sum(pop)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NATURAL POPULATIONS
   LWORK = nSPFmax * nSPFmax
! Calculate natural populations and natural orbitals for DVR modes
   do i = 1, nmodDVR
      im = imodDVR(i)
      do iel = 1, nstates
         nS = nSPF(im,iel)
         call zlacpy('U',nS,nS,rho(1,1,im,iel),nSPFmax,rhot(1,1),nSPFmax)
         call zheev('V','U',nS,rhot,nSPFmax,Wv, &
                    WORK,LWORK,RWORK,info)
         ! Natural population are set in decreasing order
         do j = 1, nS
            natpop(j,im,iel) = Wv(nS - j + 1)
         end do
      end do
   end do
! Calculate natural populations, gross Gaussian populations 
! and population errors for the GWP modes
   GGP    = dZero
   PopErr = dZero
   do i = 1, nmodGWP
      im = imodGWP(i)
      do iel = 1, nstates
         jel = merge(iel, 1, MultiSet(im))
         nS = nSPF(im,jel)
         call SquareRoot(nS, S00(1,1,i,jel,jel), nS00max, S12, nS00max)
         call zlacgv(nS00max**2, S12(1,1), 1)
         call zhemm('R','U',nS,nS,cOne, rho(1,1,im,iel),nSPFmax, &
                     S12,nS00max, cZero, Srho,nS00max)
         call zhemm('R','U',nS,nS,cOne, S12,nS00max, &
                     Srho,nS00max, cZero, rhot,nSPFmax)
         call zheev('V','U',nS,rhot,nSPFmax,Wv, &
                    WORK,LWORK,RWORK,info)
         ! Natural population are set in decreasing order
         do j = 1, nS
            natpop(j,im,iel) = Wv(nS - j + 1)
         end do
         ! Gross Gaussian population and population errors
         do j1 = 1, nS
            do j2 = 1, j1 - 1
               GGP(j1,i,iel)    = GGP(j1,i,iel)    + real(rho(j2,j1,im,iel) * S00(j2,j1,i,jel,jel))
               PopErr(j1,i,iel) = PopErr(j1,i,iel) + real(rho(j2,j1,im,iel) * S00(j2,j1,i,jel,jel)) * 2
            end do
            PopErr(j1,i,iel) = PopErr(j1,i,iel) - real(rho(j1,j1,im,iel))
            do j2 = j1, nS
               GGP(j1,i,iel)    = GGP(j1,i,iel)    + real(rho(j1,j2,im,iel) * S00(j1,j2,i,jel,jel))
               PopErr(j1,i,iel) = PopErr(j1,i,iel) + real(rho(j1,j2,im,iel) * S00(j1,j2,i,jel,jel)) * 2
            end do
            PopErr(j1,i,iel) = abs(PopErr(j1,i,iel))
         end do
         call dlasrt('D', nS, PopErr(1,i,iel), info)
      end do
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the (statistical average of the) autocorrelation function
   auto = cZero
   iA = 1
   do iW = 1, npackets
      call WF_Ovl(AVector_Start(iA),psi_Start,AVector(iA),psi,cdummy)
      auto = auto + cdummy
      iA = iA + nAConftot
   end do
!
! Write output 
!
   write(iout_unit,*)
   write(iout_unit,'(a,f9.3,a)') 'Time = ', time * au2fs, ' fs'
   write(iout_unit,'(a,f12.8,a,4x,a,f9.4,a,4x,a,f9.4,a)') &
      ' Norm = ', sum(pop), ' ,', &
      ' Energy = ', energy * au2eV, ' eV,', &
      ' Delta-E = ', (energy - energyStart) * 1000.0d0 * au2eV, ' meV'
   !
   if (npackets > 1) then
      write(iout_unit,*) 'Natural populations of the statistical mixture:'
      j = 1
      do i = 1, npackets / 6 + 1
         write(iout_unit,'(1x,6(2x,f9.7))') WPpop(j:min(j + 5,npackets))
         j = j + 6
      end do
      write(iout_unit,'(a,f9.6)') ' Relative entropy S / S_max: ', entropy
   end if
   !
   if (nstates > 1) write(iout_unit,'(a,30(2x,f9.6))') ' Electronic populations:', pop
   ! 
   write(iout_unit,*) ' Natural populations:'
   do iel = 1, nstates
      write(iout_unit,'(a,i0)') '   State ', iel
      do im = 1, nmodes
         nS = nSPF(im,iel)
         write(iout_unit,'(a,i3,a,6(2x,f8.5))') &
            '    Mode ', im, ':', (natpop(j1,im,iel), j1 = 1, min(nS,6))
         j2 = 7 
         do 
            if (j2 .gt. nS) exit
            write(iout_unit,'(a,6(2x,f8.5))') '             ', &
               (natpop(j1,im,iel), j1 = j2, min(nS,j2 + 5))
            j2 = j2 + 6
         end do
      end do
   end do
   !
   if (nmodGWP > 0) then
      write(iout_unit,*) ' Gross Gaussian populations:'
      do iel = 1, nstates
         write(iout_unit,'(a,i0)') '   State ', iel
         do i = 1, nmodGWP
            im = imodGWP(i)
            nS = nSPF(im,iel)
            write(iout_unit,'(a,i3,a,6(2x,f8.5))') &
               '    Mode ', im, ':', (GGP(j1,i,iel), j1 = 1, min(nS,6))
            j2 = 7
            do
               if (j2 .gt. nS) exit
               write(iout_unit,'(a,6(2x,f8.5))') '             ', &
                  (GGP(j1,i,iel), j1 = j2, min(nS,j2 + 5))
               j2 = j2 + 6
            end do
         end do
      end do
      write(iout_unit,*) ' Errors on the population by neglecting one Gaussian:'
      do iel = 1, nstates
         write(iout_unit,'(a,i0)') '   State ', iel
         do i = 1, nmodGWP
            im = imodGWP(i)
            nS = nSPF(im,iel)
            write(iout_unit,'(a,i3,a,6(2x,f8.5))') &
               '    Mode ', im, ':', (PopErr(j1,i,iel), j1 = 1, min(nS,6))
            j2 = 7
            do
               if (j2 .gt. nS) exit
               write(iout_unit,'(a,6(2x,f8.5))') '             ', &
                  (PopErr(j1,i,iel), j1 = j2, min(nS,j2 + 5))
               j2 = j2 + 6
            end do
         end do
      end do
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!! 
! Write auto
   write(iauto_unit,'(2x,f10.4,5x,3(1x,f18.14))') &
             time * au2fs, real(auto), aimag(auto), abs(auto)
! Write restart
   ires_unit = freeunit()
   open(ires_unit, file = 'restart', status = 'replace', form = 'unformatted', action = 'write')
   write(ires_unit) nstates, nmodes
   do i = 1, nstates
      write(ires_unit) nSPF(1:nmodes,i)
   end do
   write(ires_unit) time
   write(ires_unit) AVector
   write(ires_unit) psi
   write(ires_unit) GWPa
   close(ires_unit)
! Write spop
   if (nstates .gt. 1) then
      write(ispop_unit,'(f10.4,2x,30(2x,f10.7))') time * au2fs, pop
      flush(ispop_unit)
   end if
! Write psi
   if (wr_psi) then
      write(ipsi_unit) time, AVector, psi
      wr_psi = .false.
      flush(ipsi_unit)
   end if
! Write densmat
   if (wr_densmat) then
      write(idm_unit) time
      do i = 1, nmodes_densmat
         im = imode_densmat(i)
         select case(dvrtype(iDOF(1,im)))
            case ('gwp')
               nG = nmodeDOF(im)
            case default
               do j = 1, nmodDVR
                  if (imodDVR(j) .eq. im) exit
               end do
               nG = ngpm(j)
         end select
         do iel = 1, nstates
            nS = nSPF(im,iel)
            iP1 = ipsiSPF(im,iel)
            iP2 = iP1 + nS * nG - 1
            write(idm_unit) im, iel, rho(1:nS,1:nS,im,iel), psi(iP1:iP2)   ! Use only the upper diagonal!! The lower diagonal contains rubbish
         end do
      end do
      flush(idm_unit)
   end if
! Write density
   if (wr_density) then
      write(idens_unit) time
      if (nmodDVR > 0) then
         write(iout_unit,*) 'Population of the basis set underlying the DVR'
      end if
      do i = 1, nDOF
         nG = ngp(i)
         call Calc_Density(i)
         select case(dvrType(i))
            case('gwp')
               do iel = 1, nstates
                  write(idens_unit) i, Density(1:nG,iel)
               end do
            case default
               write(iout_unit,'(a,i2)') '   DOF:', i
               do iel = 1, nstates
                  write(idens_unit) i, Density(1:nG,iel), DensityPop(1:nG,iel)
                  if (nG == 2) then
                     write(iout_unit,'(a,i2,a,3x,f7.4,2x,f7.4)') &
                        '     State ', iel, ':', DensityPop(1,iel), DensityPop(nG,iel)
                  else
                     write(iout_unit,'(a,i2,a,3x,f7.4,2x,f7.4,a,f7.4,2x,f7.4)') &
                        '     State ', iel, ':', DensityPop(1:2,iel), ' ... ', DensityPop(nG - 1:nG,iel)
                  end if
               end do
         end select
      end do
      flush(idens_unit)
   end if 
! Write wave packet overlaps
   if (wr_wpOvl) then
      call Calc_SWP(AVector,psi,SWP)   ! It is assumed that the Gaussian overlap matrix has been already calculated
      write(iwpovl_unit,'(a,f7.2,a)') 'Time: ', time * au2fs, ' fs'
      iC = 1
      do j = 1, (npackets + 5) / 6
         write(iwpovl_unit,'(7x,6(8x,i3,8x))') (i, i = iC, min(npackets,iC + 5))
         do i = 1, npackets
            write(iwpovl_unit,'(i3,2x,6(2x,"(",f7.4,",",f7.4,")"))') i, (SWP(i,iW), iW = iC, min(npackets,iC + 5))
         enddo
         iC = iC + 6
      enddo
      flush(iwpovl_unit)
   end if
! Calculate expectation values
   if (nExpect .gt. 0) then
      call Expectation(tt)
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   flush(iout_unit)
   flush(iauto_unit)
   return
   end subroutine WriteOutput
   
   
   
   subroutine Calc_Density(id)
! It calculates the vibrational density along a specific degree of freedom (id) for all electronic states
   implicit none
   integer, intent(in) :: id
   integer :: im, iel, k, nS, kg, kg1, kg2, nG, nmD, nn, &
              j, j1, j2, iS1, iS2, igwp
   double precision :: delta, ASum, A1, A2, xi1, xi2
   double complex :: cdummy, GauOvl1D, xiSum
   double precision, dimension(:) :: gridtmp(maxval(ngp))
   double precision, parameter :: pi24 = dOne / asin(dOne) ** 2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   im = modeOfDOF(id)
   do k = 1, nmodeDOF(im)
      if (iDOF(k,im) .eq. id) exit
   end do
   nG = ngp(id)
   select case(dvrType(id)) 
      case ('gwp')
         do igwp = 1, nmodGWP 
            if (imodGWP(igwp) .eq. im) exit
         end do
         delta = (dvrpar(2,id) - dvrpar(1,id)) / dble(nG - 1)
         do kg = 1, nG
            gridtmp(kg) = dvrpar(1,id) + delta * dble(kg - 1)
         end do
         nmD = nmodeDOF(im)
         Density = dZero
         do iel = 1, nstates
            nS = nSPF(im,iel)
            iS1 = ipsiSPF(im,iel) - 1
            do j1 = 1, nS
               A1 = GWPa(k,j1,igwp)
               xi1 = real(psi(iS1 + k)) ** 2
               iS2 = ipsiSPF(im,iel) - 1
               do j2 = 1, j1 - 1
                  A2 = GWPa(k,j2,igwp)
                  xi2 = real(psi(iS2 + k)) ** 2
                  cdummy = conjg(rho(j2,j1,im,iel)) * (pi24 * A1 * A2 * exp(xi1 / A1 + xi2 / A2)) ** 0.25d0
                  do kg = 1, nmD
                     if (kg .eq. k) cycle
                     cdummy = cdummy * GauOvl1D(GWPa(kg,j1,igwp),GWPa(kg,j2,igwp), &
                                                psi(iS1 + kg),psi(iS2 + kg))
                  end do
                  ASum = A1 + A2
                  xiSum = conjg(psi(iS1 + k)) + psi(iS2 + k)
                  do kg = 1, nG
                     Density(kg,iel) = Density(kg,iel) &
                        + real(cdummy * exp(gridtmp(kg) ** 2 * ASum + gridtmp(kg) * xiSum))
                  end do
                  iS2 = iS2 + nMD
               end do
               do j2 = j1, nS
                  A2 = GWPa(k,j2,igwp)
                  xi2 = real(psi(iS2 + k)) ** 2
                  cdummy = rho(j1,j2,im,iel) * (pi24 * A1 * A2 * exp(xi1 / A1 + xi2 / A2)) ** 0.25d0
                  do kg = 1, nmD
                     if (kg .eq. k) cycle
                     cdummy = cdummy * GauOvl1D(GWPa(kg,j1,igwp),GWPa(kg,j2,igwp), &
                                                psi(iS1 + kg),psi(iS2 + kg))
                  end do
                  ASum = A1 + A2
                  xiSum = conjg(psi(iS1 + k)) + psi(iS2 + k)
                  forall (kg = 1:nG) &
                     Density(kg,iel) = Density(kg,iel) &
                        + real(cdummy * exp(gridtmp(kg) ** 2 * ASum + gridtmp(kg) * xiSum))
                  iS2 = iS2 + nMD
               end do
               iS1 = iS1 + nMD
            end do
         end do
         do concurrent (iel = 1:nstates, kg = 1:nG) 
            Density(kg,iel) = Density(kg,iel) * delta
         end do
         return
      case default
         do j = 1, nmodDVR
            if (imodDVR(j) .eq. im) exit
         end do
         nn = ngpm(j) / nG
         do iel = 1, nstates
            DensityMatrix = cZero
            nS = nSPF(im,iel)
            iS1 = ipsiSPF(im,iel)
            call MatricisePsi(iel,im,k,psi(iS1:iS1 + ngpm(j) * nS - 1),psiMat)   ! CONTROLLARE DI PASSARE LA SEZIONE DI psi GIUSTA
            ! Calculate the density matrix in coordinate representation
            do j1 = 1, nS
               iS1 = (j1 - 1) * nn         
               ! j2 .eq. j1
               cdummy = rho(j1,j1,im,iel)
               forall (kg1 = 1:nG, kg2 = 1:nG) &
                  DensityMatrix(kg1,kg2) = DensityMatrix(kg1,kg2) &
                     + cdummy * dot_product(psiMat(iS1 + 1:iS1 + nn,kg2),psiMat(iS1 + 1:iS1 + nn,kg1))
               ! j2 .ne .j1
               do j2 = j1 + 1, nS
                  iS2 = (j2 - 1) * nn
                  cdummy = rho(j1,j2,im,iel)
                  forall (kg1 = 1:nG, kg2 = 1:nG) &
                     DensityMatrix(kg1,kg2) = DensityMatrix(kg1,kg2) &   ! Forse kg1 e kg2 vanno invertiti
                        + cdummy * dot_product(psiMat(iS1 + 1:iS1 + nn,kg2),psiMat(iS2 + 1:iS2 + nn,kg1)) &
                        + conjg(cdummy) * dot_product(psiMat(iS2 + 1:iS2 + nn,kg2),psiMat(iS1 + 1:iS1 + nn,kg1))
               end do
            end do
            ! Density
            !call zlacrm(nG,nG,DensityMatrix(1,1),maxval(ngp),trafoDVR(1,1,indDVR(id)),maxval(ngp),DT(1,1),maxval(ngp))
            forall (kg = 1:nG) 
               Density(kg,iel) = real(DensityMatrix(kg,kg))
               !DensityPop(kg,iel) = dot_product(trafoDVR(1:nG,kg,indDVR(id)), DT(1:nG,kg))
               DensityPop(kg,iel) = real(dot_product(trafoDVR(kg,1:nG,indDVR(id)), &
                                                matmul( DensityMatrix(1:nG,1:nG),trafoDVR(kg,1:nG,indDVR(id)) ) ))
            end forall
         end do
   end select
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   end subroutine Calc_Density



   subroutine Expectation(tt)
   ! It calculates the requested expectation values for the current wavefunction
   ! It assumes that the integrals on the array h1mode have been already calculated
   implicit none
   double precision, intent(in) :: tt
   integer :: i, iOp, iT, iel1, iel2, iA1, iA2, nAC1, nAC2, im, nS1, nS2, jOp, iW
   integer, dimension(nmodes) :: sizeTens
   double complex, dimension(nExpect) :: expect
   double complex, dimension(:), pointer :: pA1, pA2, pswap
   double complex, external :: zdotu
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   expect = cZero
   do iW = 1, npackets
      do i = 1, nExpect
         iOp = iExpect(i)
         do iT = 1, nOpTerms(iOp)
            iel1 = lhs_elExtra(iT,iOp)
            iel2 = rhs_elExtra(iT,iOp)
            iA1 = APacketOffset(iW) + AStateOffset(iel1)
            iA2 = APacketOffset(iW) + AStateOffset(iel2)
            nAC1 = nAConf(iel1)
            nAC2 = nAConf(iel2)
            sizeTens = nSPF(:,iel2)
            call zcopy(nAC2,AVector(iA2 + 1),1,Aux2(1),1)
            call zlacgv(nAC2,Aux2(1),1)
            pA1 => Aux1
            pA2 => Aux2
            do im = 1, nmodes
               nS1 = nSPF(im,iel1)
               nS2 = nSPF(im,iel2)
               ! Multiply by the 1-mode-operator matrix
               jOp = iOperExtra(im,iT,iOp)
               if (iel1 .eq. iel2 .or. .not. MultiSet(im)) then
               ! h1mode is assumed to be Hermitian
                  call Mult_HeM_MatV(sizeTens,im,h1mode(1,1,jOp,im),nSPFmax,pA2,pA1)
               else
                  call Mult_GM_MatV(sizeTens,im,'C',nS2,nS1,h1mode(1,1,jOp,im),nSPFmax,pA2,pA1)
               end if
               pswap => pA2
               pA2 => pA1
               pA1 => pswap
               sizeTens(im) = nS1
            end do
            expect(i) = expect(i) + zdotu(nAC1,pA2(1),1,AVector(iA1 + 1),1) * conjg(aeExtra(iT,iOp))
         end do
      end do
   end do
   ! Warning: complex conjugation
   write(iexpect_unit,'(1x,f8.2,2x,800(3x, f14.9,1x,f14.9))') tt * au2fs, conjg(expect)
   flush(iexpect_unit)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   end subroutine Expectation


end module output
