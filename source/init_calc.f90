subroutine init_calc
! It initialises all arrays which are used in the propagation
use sysparam
use globals
use psidef
use dvr
use hamilton
use storage
implicit none
integer :: a_status
integer :: i, j, k, iC, nmD
double precision :: cnorm, dummy
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Write general information
write(iout_unit,'(a,i0,a)') 'Propagation runs on ', nstates, ' electronic states'
write(iout_unit,'(a,i0,a)') '  and includes ', ndof, ' degrees of freedom,'
write(iout_unit,'(a,i0,a)') '  combined in ', nmodes, ' logical modes'
write(iout_unit,*)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Calculate the number of configurations for each state
! and allocate the coefficient vector
allocate(nAconf(nstates), APacketOffset(npackets), AStateOffset(nstates), stat = a_status)
if (a_status .ne. 0) call err_alloc('nAconf/AStateOffset', 'init_calc', a_status)
!
nAConf = product(nSPF, dim = 1)
write(iout_unit,'(a)') ' el. state |  n. of configurations '
write(iout_unit,'(a)') '-----------------------------------'
iC = 0
do i = 1, nstates
   AStateOffset(i) = iC
   write(iout_unit,'(5x,i3,3x,a,4x,i7)') i, '|', nAconf(i)
   iC = iC + nAConf(i)
enddo
nAconfTot = iC
! Check that the number of configurations is larger than the number of wave packets
if (nAConfTot .lt. npackets) then
   write(0,*)
   write(0,*) 'ERROR: With the given G-MCTDH specifications'
   write(0,*) '       the number of wave packets is larger'
   write(0,*) '       than the available number of configurations:'
   write(0,'(8x,i0,a,i0,a)') npackets, ' wave packets, ', nAConfTot, ' configurations'
   write(0,*)
   write(0,*) 'Increase the number of SPFs for one or mode modes'
   stop
endif
!
allocate(WPpop(npackets), multLag(npackets,npackets))
allocate(Avector(nAconfTot * npackets), stat = a_status)
if (a_status .ne. 0) call err_alloc('Avector', 'init_calc', a_status)
do i = 1, npackets
   APacketOffset(i) = nAconfTot * (i - 1)
end do
Avector = cZero
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Set-up DVR grids and derivative matrices
   ! Get which modes use DVR and which use GWP representation
nmodDVR = 0
nmodGWP = 0
do i = 1, nmodes
   if (trim(dvrtype(iDOF(1,i))) .eq. 'gwp') then
      nmodGWP = nmodGWP + 1
   else
      nmodDVR = nmodDVR + 1
   endif
enddo
allocate(imodDVR(nmodDVR), imodGWP(nmodGWP))
i = 0
j = 0
nS00max = 0
do k = 1, nmodes
   if (trim(dvrtype(iDOF(1,k))) .eq. 'gwp') then
      i = i + 1
      imodGWP(i) = k
      nS00max = max(nS00max, maxval(nSPF(k,:)))
   else
      j = j + 1
      imodDVR(j) = k
   endif
enddo
   ! Get which dofs use DVR and which use GWP representation
ndofDVR = count(dvrtype .ne. 'gwp')
ndofGWP = count(dvrtype .eq. 'gwp')
allocate(idofDVR(ndofDVR), idofGWP(ndofGWP), GWPa(ndofGWP,nS00max,ndofGWP), indDVR(ndof))
i = 0
j = 0
do k = 1, ndof
   if (trim(dvrtype(k)) .eq. 'gwp') then
      i = i + 1
      idofGWP(i) = k
   else
      j = j + 1
      idofDVR(j) = k
      indDVR(k) = j
   endif
enddo
   ! calculate DVR derivative matrices
allocate(operDVR(ngpMax,ngpMax,ndofDVR), d1DVR(ngpMax,ngpMax,ndofDVR), &
         d2DVR(ngpMax,ngpMax,ndofDVR), xd1DVR(ngpMax,ngpMax,ndofDVR), &
         trafoDVR(ngpMax,ngpMax,ndofDVR), grid(ngpMax,ndofDVR), wgrid(ngpMax,ndofDVR))
do i = 1, ndofDVR
   k = idofDVR(i)
   select case(trim(dvrtype(k)))
      case('sin')
         call sin_DVR(ngp(k), dvrpar(1,k), dvrpar(2,k), &
                      grid(1,i), wgrid(1,i), d1DVR(1,1,i), d2DVR(1,1,i), xd1DVR(1,1,i), &
                      trafoDVR(1,1,i), ngpMax)
      case('HO')
         call harm_DVR(ngp(k), dvrpar(1,k), dvrpar(2,k), &
                       grid(1,i), wgrid(1,i), d1DVR(1,1,i), d2DVR(1,1,i), xd1DVR(1,1,i), &
                       trafoDVR(1,1,i), ngpMax)
      case('Leg')
         call Leg_DVR(ngp(k), idvrpar(k), grid(1,i), wgrid(1,i), operDVR(1,1,i), d1DVR(1,1,i), d2DVR(1,1,i), &
                      trafoDVR(1,1,i), ngpMax)
      case('2pi')
         call twopi_DVR(ngp(k), grid(1,i), wgrid(1,i), d1DVR(1,1,i), d2DVR(1,1,i), &
                      trafoDVR(1,1,i), ngpMax)
      case('par')
         dummy = (dvrpar(2,k) - dvrpar(1,k)) / dble(ngp(k) - 1)
         do j = 1, ngp(k)
            grid(j,i) = dvrpar(1,k) + dummy * dble(j - 1)
         enddo
   end select
enddo
! Set up multi-dimensional direct tensor product DVR grids
   ! Calculate ngpm
allocate(ngpm(nmodDVR))
do i = 1, nmodDVR
   j = imodDVR(i)
   nmD = nmodeDOF(j)
   ngpm(i) = product(ngp(iDOF(1:nmD,j)))
enddo
ngpmMax = maxval(ngpm)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Initialise the logical-mode-structure of the Hamiltonian
! Initialise Hamiltonian arrays needed in the propagation
call alloc_h1mode
call alloc_IsIdentity
!call scale_operators    ! For the parametrized operators this does not work!
   ! Write information about the operator basis in the log file
write(ilog_unit,*)
write(ilog_unit,*) 'Operators for the individual modes'
do j = 1, nmodes
   write(ilog_unit,'(a,i0)') '   Mode ', j
   write(ilog_unit,*) '  #  el. states    elementary operators'
   do i = 1, nOper(j)
      write(ilog_unit,'(i4,4x,2(1x,i2),6x,30(i0,1x))') i, Oper(-1:nmodeDOF(j),i,j)
   enddo
   write(ilog_unit,*)
enddo
flush(ilog_unit)
! Allocate Mean Fields arrays
call alloc_MeanField
call alloc_rhom1MF
! Allocate the arrays which are needed for propagation
call alloc_S_matrices
call alloc_CY
call alloc_rho
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Initialise psi vector
allocate(ipsiSPF(nmodes,nstates))
dimPsi = 0
do i = 1, nstates
   do j = 1, nmodes
      if (.not. MultiSet(j) .and. i > 1) then
         ipsiSPF(j,i) = ipsiSPF(j,1)
         cycle
      end if
      select case(trim(dvrtype(idof(1,j))))
         case('gwp')
            ipsiSPF(j,i) = dimPsi + 1
            dimPsi = dimPsi + nmodeDOF(j) * nSPF(j,i)
         case default
            iC = 1
            do k = 1, nmodeDOF(j)
               iC = iC * ngp(iDOF(k,j))
            enddo
            ipsiSPF(j,i) = dimPsi + 1
            dimPsi = dimPsi + iC * nSPF(j,i)
      end select
   enddo
enddo
!
write(ilog_unit,*)
write(ilog_unit,'(a,i0,a)') 'The SPF-vector contains ', dimPsi, ' variational parameters'
!!!!
! In case of a harmonic oscillator bath the cross-correlation functions
! are written in the psi vector (better for the integrator)
if (lHObath) then
   allocate(FQbath(npackets,npackets,nHObath))
   allocate(FPbath(npackets,npackets,nHObath))
   iCQ = dimpsi + 1                  ! pointer to the beginning of the CQ vector
   dimpsi = dimpsi + nHObath * npackets**2
   iCP = dimpsi + 1                  ! pointer to the beginning of the CP vector
   dimpsi = dimpsi + nHObath * npackets**2
   iCQQ = dimpsi + 1
   dimpsi = dimpsi + nHObath * npackets**2
   iCQP = dimpsi + 1
   dimpsi = dimpsi + nHObath * npackets**2
   iCPQ = dimpsi + 1
   dimpsi = dimpsi + nHObath * npackets**2
   iCPP = dimpsi + 1
   dimpsi = dimpsi + nHObath * npackets**2
end if
!!!!
! Read initial wavefunction from the input file
allocate(psi(dimpsi), psi_Start(dimpsi), stat = a_status)
if (a_status .ne. 0) call err_alloc('psi', 'init_calc', a_status)
psi = cZero
call alloc_AVec   ! init_wf needs Aux1 and Aux2 to be allocated
call init_wf
! Write information about the initial Gaussian distribution
if (nmodGWP .gt. 0) call write_gau_pos
write(ilog_unit,*) ''
! Normalise the initial wavefunction mixture
call trace_wf(AVector,psi,cnorm)
!call norm_wf(AVector,psi,cnorm)
write(ilog_unit,*)
write(ilog_unit,'(a,f11.6)') 'Wavefunction norm before normalisation: ', cnorm
AVector = AVector / sqrt(cnorm)
AVector_Start = AVector
psi_Start = psi
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
flush(iout_unit)
flush(ilog_unit)
return
end subroutine init_calc



subroutine init_wf
! It sets the initial wavefunctions
! Three possibilities are available:
! 1) 'restart'
!    The initial wavefunction is read from a restart file (not yet implemented)
! 2) 'readpsi'
!    The initial wavefunction is read from a psi file (not yet implemented)
! 3) 'single-conf'
!    The first configuration of a given state is populated.
!    The unpopulated SPF functions are automatically orthogonalised
use sysparam
use globals
use psidef
use dvr
use Hamilton, only : lHObath, nHObath
implicit none
integer :: io_inp, io_psi, i, j, k, im, i2, iC, iP, nS, iA, iASize, iW
integer :: iel, kdof, nG, id, LL, MM
character(len = 200) :: line, line2, readPsiFile
double precision :: dummy,  xx
double precision :: x0, p0, dx, cR, cI, TKin
double precision :: cR1, cI1, x1, alpha1, cR2, cI2, x2, alpha2   ! for the 'cat' superposition
double complex :: cdummy
integer, dimension(nmodes)  :: jSPF
integer, dimension(ndofGWP) :: iDistGWP
double precision, dimension(3,ndofGWP) :: ParGWP
double complex, dimension(:,:), allocatable :: phi
character(20) :: muFile   ! dipole moment file
integer :: io_mu
double precision, dimension(:) :: tdm(1000)
logical :: exists
double precision, external :: dznrm2
! For reading from restart file...
integer :: nstates_t, dimpsi_t, iP_t, loop_length, nloops, iloop, iel_Start
integer, dimension(:), allocatable :: associatedstate, nAConf_t
integer, dimension(:,:), allocatable :: nSPF_t, Jindex_t, Jindex
integer, dimension(:,:,:), allocatable :: ipsiSPF_t
double complex, dimension(:), allocatable :: AVector_t, psi_t
! Auxiliary
integer, parameter :: LWORK = 200
double precision, external :: AssociatedLegendreP
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
io_inp = freeunit()
open(unit = io_inp, file = trim(adjustl(input_file_psi)), form = 'formatted', action = 'read')
do
   line = ''
   read(io_inp, '(a200)') line
   line2 = line
   i = index(line,'#')
   if (i .ne. 0) line2 = line(1:i - 1)
   if(trim(adjustl(line2)) .eq. 'INIT-WF') exit
enddo
! Read whether the wavefunction should be read from an external file
! or a single configuration is given
read(io_inp, '(a200)') line
line2 = line
i = index(line,'#')
if (i .ne. 0) line2 = line(1:i - 1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! SPFs AND COEFFICIENTS FROM RESTART FILE
if (trim(adjustl(line2)) .eq. 'restart') then
   if (npackets > 1) then
      write(0,*) 'ERROR'
      write(0,*) '   restart not yet implemented for statistical mixtures'
      write(0,*)
      stop
   end if
   read(io_inp, '(a200)') line
   readPsiFile = trim(adjustl(line))
   io_psi = freeunit()
   open(unit = io_psi, file = readPsiFile, form = 'unformatted', status = 'old', action = 'read')
   ! Read the specifications of the wavefunction on the psi file
   read(io_psi) nstates_t, i   ! i = nmodes
   if (i .ne. nmodes) then
      write(0,*) 'Error: '
      write(0,'(a,a,a,i0,a,i0)') &
         '   The wavefunction on the file ', trim(adjustl(readPsiFile)), &
         ' should consist of ', nmodes, ' modes instead of ', i
      stop
   endif
   allocate(nSPF_t(nmodes,nstates_t))
   do i = 1, nstates_t
      read(io_psi) nSPF_t(1:nmodes,i)
   enddo
   ! Read the electronic state association
   allocate(associatedstate(nstates))
   read(io_inp, *) iel_Start   ! initial state (all coefficients for the other states are made zero)
                            ! if iel_Start = 0 just copy the coefficients
   do i = 1, nstates
      read(io_inp, '(a200)') line
      j = index(line,'<--')
      if (j .eq. 0) then
         write(0,*) 'Error in INIT-WF:'
         write(0,*) '   Electronic state associations must be specified'
      endif
      read(line(1:j - 1),*) iel
      read(line(j + 3:),*) associatedstate(iel)
   enddo
   ! check that for the GWP modes enough Gaussians are available
   do iel = 1, nstates
      j = associatedstate(iel)
      do i = 1, nmodGWP
         im = imodGWP(i)
         if (nSPF_t(im,j) .ne. nSPF(im,iel)) then
            write(0,*) 'Error in INIT-WF:'
            write(0,'(a,i0)') '  Mode ', im
            write(0,*) '  For GWPs the number of SPFs should be the same'
            write(0,*) nSPF_t(im,j), nSPF(im,iel)
         endif
      enddo
   enddo
   ! Set up and read the temporary A and psi vectors
   allocate(nAConf_t(nstates_t))
   forall(i = 1:nstates_t) nAconf_t(i) = product(nSPF_t(:,i))
   allocate(Avector_t(sum(nAconf_t)))
   allocate(ipsiSPF_t(maxval(nSPF_t),nmodes,nstates_t))
   dimpsi_t = 0
   do i = 1, nstates_t
      do j = 1, nmodes
         if (.not. MultiSet(j) .and. i > 1) then
            forall (k = 1:nSPF_t(j,1)) ipsiSPF_t(k,j,i) = ipsiSPF_t(k,j,1)
            cycle
         end if
         select case(trim(dvrtype(idof(1,j))))
            case('gwp')
               do k = 1, nSPF_t(j,i)
                  ipsiSPF_t(k,j,i) = dimPsi_t + 1
                  dimPsi_t = dimPsi_t + nmodeDOF(j)
               enddo
            case default
            iC = 1
            do k = 1, nmodeDOF(j)
               iC = iC * ngp(iDOF(k,j))
            enddo
            do k = 1, nSPF_t(j,i)
               ipsiSPF_t(k,j,i) = dimPsi_t + 1
               dimPsi_t = dimPsi_t + iC
            enddo
         end select
      enddo
   enddo
   allocate(psi_t(dimpsi_t))
   read(io_psi) dummy
   read(io_psi) AVector_t
   read(io_psi) psi_t
   read(io_psi) GWPa
   ! Set up the single particle functions
   do iel = 1, nstates
      do i = 1, nmodDVR
         im = imodDVR(i)
         nG = ngpm(i)
         k  = min(nSPF(im,iel),nSPF_t(im,associatedstate(iel)))
         iP = ipsiSPF(im,iel)
         do j = 1, k
            iP_t = ipsiSPF_t(j,im,associatedstate(iel))
            psi(iP:iP + nG - 1) = psi_t(iP_t:iP_t + nG - 1)
            iP = iP + nG
         enddo
         ! Fill the remaining SPFs and orthogonalize
         allocate(phi(nG,1))
         do j = k + 1, nSPF(im,iel)
            iP = ipsiSPF(im,iel) + (j - 2) * nG
            forall(i2 = 1:nG) phi(i2,1) = 0.001d0 * cOne * (i2 - nG / 2) * psi(iP + i2 - 1)
            do i2 = 1, j - 1
               iP = ipsiSPF(im,iel) + (i2 - 1) * nG
               cdummy = dot_product(psi(iP:iP + nG - 1),phi(1:nG,1))
               phi(1:nG,1) = phi(1:nG,1) - psi(iP:iP + nG - 1) * cdummy
            enddo
            dummy = dznrm2(nG, phi(1,1), 1)
            call zdscal(nG, dOne / dummy, phi(1,1), 1)
            ! orthogonalize a second time (might improve the orthogonalization if 'dummy' is very small...
            do i2 = 1, j - 1
               iP = ipsiSPF(im,iel) + (i2 - 1) * nG
               cdummy = dot_product(psi(iP:iP + nG - 1),phi(1:nG,1))
               phi(1:nG,1) = phi(1:nG,1) - psi(iP:iP + nG - 1) * cdummy
            enddo
            dummy = real(dot_product(phi(1:nG,1),phi(1:nG,1)))
            phi(1:nG,1) = phi(1:nG,1) / sqrt(dummy)
            !
            iP = ipsiSPF(im,iel) + (j - 1) * nG
            psi(iP:iP + nG - 1) = phi(1:nG,1)
         enddo
         deallocate(phi)
         !
      enddo
      do i = 1, nmodGWP
         im = imodGWP(i)
         nG = nmodeDOF(im)
         iP = ipsiSPF(im,iel)
         do j = 1, nSPF(im,iel)
            iP_t = ipsiSPF_t(j,im,associatedstate(iel))
            psi(iP:iP + nG - 1) = psi_t(iP_t:iP_t + nG - 1)
            !GWPa(1:nG,j,i) = - dOneHalf   ! temporary
            iP = iP + nG
         enddo
      enddo
   enddo
   ! Set up the coefficients
      ! temporary Jindex
   allocate(Jindex_t(nmodes,sum(nAConf_t)))
   allocate(Jindex(nmodes,nAconfTot))
      ! Jindex
   do i = 1, nstates
      do j = 1, nmodes
         loop_length = product(nSPF(:j - 1,i))
         nloops = product(nSPF(j + 1:,i))
         iC = sum(nAconf(:i - 1))
         do iloop = 1, nloops
            do k = 1, nSPF(j,i)
               Jindex(j,iC + 1:iC + loop_length) = k
               iC = iC + loop_length
            enddo
         enddo
      enddo
   enddo
      ! Jindex_t
   do i = 1, nstates_t
      do j = 1, nmodes
         loop_length = product(nSPF_t(:j - 1,i))
         nloops = product(nSPF_t(j + 1:,i))
         iC = sum(nAconf_t(:i - 1))
         do iloop = 1, nloops
            do k = 1, nSPF_t(j,i)
               Jindex_t(j,iC + 1:iC + loop_length) = k
               iC = iC + loop_length
            enddo
         enddo
      enddo
   enddo
   AVector = cZero
   do iel = 1, nstates
      iC = AStateOffset(iel)
      i2 = associatedstate(iel)
      do i = iC + 1, iC + nAConf(iel)
         k = sum(nAconf_t(:i2 - 1))
         iP = 0
         do j = k + 1, k + nAconf_t(i2)
            if (all(Jindex_t(:,j) .eq. Jindex(:,i))) then
               iP = j
               exit
            endif
         enddo
         if (iP .ne. 0 .and. (iel_Start == iel .or. iel_Start == 0)) AVector(i) = AVector_t(iP)
      enddo
   enddo
   !
   close(io_psi)
   deallocate(nSPF_t,associatedstate)
   deallocate(AVector_t,psi_t)
   deallocate(ipsiSPF_t,Jindex_t,Jindex)
   return
endif
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! SINGLE CONFIGURATION SPECIFICATION
if (trim(adjustl(line2)) .eq. 'single-conf') then  ! Initially populated configuration is read
   allocate(phi(max(1,maxval(ngp)),ndof))
   read(io_inp, '(a200)') line
   i = index(line, '#')
   if (i .ne. 0) line = line(1:i - 1)
   !!!!!
   if (index(line, 'el') .ne. 0) then
      ! Read the electronic state and initialize the 'virtual' wave packets automatically
      ! with populations given by eps_pop * n, with n = 2,3,4,...
      i = index(line, 'el')
      read(line(i + 2:), *) iel
      if (all(dvrtype(1:nDOF) .ne. 'par')) then
         Avector = cZero
         do iW = 1, npackets
            iA = APacketOffset(iW) + AStateOffset(iel)
            AVector(iA + iW) = merge(cOne, cOne * sqrt(eps_pop * iW), iW == 1)
         end do
      else  ! parametric coordinates
         if (npackets > 1) then
            write(0,*) 'ERROR:'
            write(0,*) 'Multi-packet wavefunctions with parametric coordinates'
            write(0,*) ' not yet implemented'
            stop
         end if
         iA = AStateOffset(iel)
         ! First all coefficients are set to one
         AVector(iA + 1:iA + nAConf(iel)) = cOne
         ! Then, the coefficients not associated to the SPF 1 are set to zero
         do im = 1, nmodes
            if (dvrtype(iDOF(1,im)) .eq. 'par') cycle
            loop_length = product(nSPF(1:im - 1,iel))
            nS = nSPF(im,iel)
            nloops = product(nSPF(im + 1:nmodes,iel))
            iA = AStateOffset(iel)
            iASize = nS * loop_length
            do iloop = 1, nloops
               AVector(iA + 1 + loop_length:iA + iASize) = cZero
               iA = iA + iASize
            enddo
         enddo
         !
      end if
      !!!!!
   elseif (index(line, 'packets') .ne. 0) then
      ! The states of the mixture are declared explicitly
      AVector = cZero
      do iW = 1, npackets
         read(io_inp,*) iel, jSPF(1:nmodes), cR, cI
         iA = APacketOffset(iW) + AStateOffset(iel) + 1
         do k = nmodes, 1, -1
            iA = iA + (jSPF(k) - 1) * product(nSPF(:k - 1,iel))
         end do
         AVector(iA) = complex(cR,cI)
         write(ilog_unit,'(i7,a,"(",f7.4,",",f7.4,")",i2,2x,20(2x,i2))') &
            iW, '   |   ', cR, cI, iel, jSPF
      end do
   endif
   !
   if (npackets == 1) then
      write(ilog_unit,*)
      write(ilog_unit,'(a,i0)') 'The initial wavefunction is set on the electronic state ', iel
   end if
   !
   ! Read the initial SPFs
   !
   do k = 1, ndof
      read(io_inp, '(a200)') line
      if (dvrtype(k) .eq. 'par') cycle
      i = index(line, '#')
      if (i .ne. 0) line = line(1:i - 1)
      !
      i = index(line, 'gauss')
      if (i .ne. 0)  then ! Gaussian function
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .eq. 'gwp') stop 'Input Error: use "gwp" instead of "gauss"'
         read(line(i + 5:), *) x0, p0, dx
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. kdof) exit
         enddo
         if (abs(x0) < epsilon(x0) .and. abs(p0) < epsilon(p0)) then  ! symmetric version
            do j = 1, ngp(kdof) / 2 + 1
               phi(j,kdof) = complex(exp(- (dOneHalf * grid(j,i2) / dx) ** 2) * sqrt(wgrid(j,i2)), dZero)
               phi(ngp(kdof) - j + 1,kdof) = phi(j,kdof)
            end do
         else
            do j = 1, ngp(kdof)
               xx = grid(j,i2) - x0
               phi(j,kdof) = exp(complex(- (dOneHalf * xx / dx) ** 2, p0 * xx)) &
                           * sqrt(wgrid(j,i2))
            enddo
         end if
         cycle
      endif
      !
      i = index(line, 'Leg')
      if (i .ne. 0) then  ! associated Legendre function
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .ne. 'Leg') stop 'Input Error: Legendre functions available only for Leg-DVR'
         do id = 1, nmodDVR
            if (idofDVR(id) == kdof) exit
         end do
         read(line(i + 3:), *) LL, MM   ! associated legendre function
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. kdof) exit
         enddo
         do j = 1, ngp(kdof)
            phi(j,kdof) = complex(AssociatedLegendreP(LL,MM,grid(j,i2)) * sqrt(wgrid(j,i2)), dZero)
            !phi(j,kdof) = complex(trafoDVR(i2 + 1,j,id), dZero)
         enddo
         cycle
      endif
      !
      i = index(line, '2pi')
      if (i .ne. 0) then  ! ring-periodic function
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .ne. '2pi') stop 'Input Error: 2pi functions available only for 2pi-DVR'
         do id = 1, nmodDVR
            if (idofDVR(id) == kdof) exit
         end do
         read(line(i + 3:), *) i2   ! eigenstate number (0 is the ground state)
         do j = 1, ngp(kdof)
            phi(j,kdof) = complex(trafoDVR(i2 + 1,j,id), dZero)
         enddo
         cycle
      endif
      !
      i = index(line, 'cat')
      if (i .ne. 0) then ! 'Cat' superposition of two Gaussians
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .eq. 'gwp') stop 'Input Error: cat initial SPFs are implemented only for grid dvr'
         read(line(i + 3:),*) cR1, cI1, x1, alpha1, cR2, cI2, x2, alpha2
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. kdof) exit
         enddo
         do j = 1, ngp(kdof)
            xx = grid(j,i2)
            phi(j,kdof) = exp(- alpha1 * (xx - x1) ** 2) * complex(cR1,cI1) &
                        + exp(- alpha2 * (xx - x2) ** 2) * complex(cR2,cI2)
            phi(j,kdof) = phi(j,kdof) * sqrt(wgrid(j,i2))
         enddo
         cycle
      endif
      !
      i = index(line, 'mu-gau')
      if (i .ne. 0) then ! TDM function, on the file 'mu_dof' applied to a Gaussian
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .eq. 'gwp') stop 'Input Error: mu*Gaussian initial SPFs are implemented only for grid dvr'
         read(line(i + 6:), *) x0, p0, dx
         ! read the dipole function
         write(muFile,'(a,i0)') 'TDM_', kdof
         inquire(file = muFile, exist = exists)
         if (.not. exists) then
            write(0,*) 'INPUT ERROR: dipole function file ', muFile, ' not found'
            stop
         endif
         io_mu = freeunit()
         open(unit = io_mu, file = muFile, status = 'old', form = 'formatted', action = 'read')
         do j = 1, ngp(kdof)
            read(io_mu,*) xx, tdm(j)
         enddo
         close(io_mu)
         !
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. kdof) exit
         enddo
         do j = 1, ngp(kdof)
            xx = grid(j,i2) - x0
            phi(j,kdof) = exp(complex(- (dOneHalf * xx / dx) ** 2, p0 * xx)) &
                          * sqrt(wgrid(j,i2)) * tdm(j)
         enddo
         cycle
      endif
      !
      i = index(line, 'eigenf1D')
      if (i .ne. 0) then  ! eigenfunction of a 1D potential
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .eq. 'gwp') stop 'Input Error: eigenf1D initial SPFs are implemented only for grid dvr'
         read(line(i + 8:), *) Tkin, j   ! Kinetic energy prefactor, selected eigenstates
         write(muFile,'(a,i0)') 'Pot1D_', kdof
         inquire(file = muFile, exist = exists)
         if (.not. exists) then
            write(0,*) 'INPUT ERROR: 1D potential file ', muFile, ' not found'
            stop
         endif
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. kdof) exit
         enddo
         ! Call a function to evaluate the eigenstate
         call Eigenf1DPot(ngp(kdof),grid(1,i2), Tkin,d2DVR(1,1,i2),ngpmax, muFile, j,phi(1,kdof))
         cycle
      end if
      !
      i = index(line, 'gwp')
      if (i .ne. 0) then  ! gwp distribution
         read(line(:i - 1), *) kdof
         if (dvrtype(kdof) .ne. 'gwp') &
            stop 'Input Error: gwp initial SPFs are allowed only for "gwp representation"'
         do j = 1, ndofGWP
            if (idofGWP(j) .eq. kdof) exit
         enddo
         read(line(i + 3:), *) i2, x0, p0, dx !, ovl
         ! check for overlap specification
         iC = index(line,'ovl')
         if (iC .eq. 0) then
            ParGWP(2,j) = 0.7d0
         else
            line2 = line(iC:)
            iC = index(line2, '=')
            if (iC .eq. 0) stop 'Input Error: keyword ovl must be followed by ='
            read(line2(iC + 1:),*) ParGWP(2,j)
         endif
         ! check for width scaling specification
         iC = index(line,'scw')
         if (iC .eq. 0) then
            ParGWP(3,j) = dOne
         else
            line2 = line(iC:)
            iC = index(line2, '=')
            if (iC .eq. 0) stop 'Input Error: keyword scw must be followed by ='
            read(line2(iC + 1:),*) ParGWP(3,j)
         endif
         !
         ParGWP(1,j) = dx
         iDistGWP(j) = i2
         if (max(abs(x0), abs(p0)) < epsilon(x0)) then
            phi(1,kdof) = cZero
         else
            phi(1,kdof) = complex(dOneHalf * x0 / dx ** 2, p0)
         end if
         cycle
      endif
      !
      write(0,*) 'Input error:'
      write(0,*) '   The chosen initial SPF is unknown'
      write(0,*) line
      stop
   enddo
   !
   call SetSPF(phi,iDistGWP,ParGWP)
   ! Orthogonalize the wave packets so to obtained population-weighted natural states
   call OrthoWP(Avector,psi)
   close(io_inp)
   !
   ! If a harmonic oscillator bath is included, initialize the values of the cross-correlation functions
   if (lHObath) then
      psi(iCQ:iCQQ - 1) = cZero
      i = iCQQ
      do i2 = 1, nHObath
         call zlaset('All',npackets,npackets,cZero,complex(dOneHalf,dZero),psi(i),npackets)
         i = i + npackets**2
      end do
      i = iCQP
      do i2 = 1, nHObath
      call zlaset('All',npackets,npackets,cZero,complex(dZero,dOneHalf),psi(i),npackets)
         i = i + npackets**2
      end do
      i = iCPQ
      do i2 = 1, nHObath
      call zlaset('All',npackets,npackets,cZero,complex(dZero,-dOneHalf),psi(i),npackets)
         i = i + npackets**2
      end do
      i = iCPP
      do i2 = 1, nHObath
      call zlaset('All',npackets,npackets,cZero,complex(dOneHalf,dZero),psi(i),npackets)
         i = i + npackets**2
      end do
   end if
   !
   return
endif
!
write(0,*) 'Input error: '
write(0,*) '   The INIT-WF block must start with the keyword read or single-conf'
stop
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
end subroutine init_wf



subroutine SetSPF(phi,iDistGWP, ParGWP)
! It initialises the SPF functions for the initial state
! in which only the first configuration of the state iel is populated.
! The SPFs of the populated modes are defined in the array phi.
! The array iDistGWP contains flags which are associated with different distributions for the unpopulated Gaussians.
! The array ParGWP contains parameters which defines the distributions of unpopulated Gaussians (1 = width, 2 = overlap)
use sysparam
use globals
use psidef
use dvr
use hamilton, only : ll_psi, nl_psi
implicit none
integer :: i, j, k, l, kg, im, iP, ngpt, nmD, jel, nG,  &
           i1, i2, i3, i4, j1, j2, nS, iDG, io_grid, kg_save, jj, idummy
!integer, intent(in) :: ngpmax, &  ! maximum number of grid points
!                             iel        ! initially populated electronic state
integer, dimension(ndofGWP), intent(in) :: iDistGWP
integer, dimension(:,:), allocatable :: index_lat
double precision :: dummy, dx, scw
double precision, dimension(:), intent(in) :: ParGWP(3,ndofGWP)
double precision, dimension(:) :: grid_op1D(ngpmax)
double precision, dimension(:,:), allocatable :: w_values, grid_op
double complex :: cdummy
double complex, dimension(:,:), intent(in) :: phi(ngpmax,ndof)
double complex, dimension(:), allocatable :: psi2_t
double complex, dimension(:,:), allocatable :: psi_t, xi_values
double precision, external :: dznrm2  ! Lapack
double complex, external :: zdotc
! for grid_index
integer :: iC, loop_length, nloops, iloop
integer, dimension(:,:,:), allocatable :: grid_index  ! grid_index(i,j,k) = grid point for the dof i corresponding to the element j
                                                             ! of the direct product grid for mode k
character(len = 20) :: GWPgridfile ! file containing the initial gwp distribution
logical :: exists
double precision, parameter :: sqrt2 = sqrt(2.d0)
! Auxiliary
integer :: info
integer, parameter :: LWORK = 200
double complex, dimension(:) :: tau(100), WORK(LWORK)
double complex, dimension(:,:) :: ovl(100,100)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! grid_index and ovl are constructed and deallocated after building the initial SPFs
allocate(grid_index(maxval(nmodeDOF), maxval(ngpm), nmodDVR))
do i = 1, nmodDVR
   im = imodDVR(i)
   nmD = nmodeDOF(im)
   do j = 1, nmD
      loop_length = product(ngp(iDOF(1:j - 1,im)))
      nloops = product(ngp(iDOF(j + 1:nmD,im)))
      iC = 0
      do iloop = 1, nloops
         do k = 1, ngp(iDOF(j,im))
            grid_index(j,iC + 1:iC + loop_length,i) = k
            iC = iC + loop_length
         enddo
      enddo
   enddo
enddo
! 1) Modes described by DVR grids
do i = 1, nmodDVR
   im = imodDVR(i)
   ! run over the n-dimensional grid
   nmD = nmodeDOF(im)
   ngpt = ngpm(i)
   nS = maxval(nSPF(im,:))
   !!!
   ! Procedure for the 'parameter coordinates': the primitive grid points are the SPFs
   if (all(dvrtype(iDOF(1:nmD,im)) .eq. 'par')) then
      ! Check that the parametric coordinates are not combined with the non-parametric ones
      if (any(dvrtype(iDOF(1:nmD,im)) .ne. 'par')) then
         write(0,'(a,1x,i0)') ' Error with mode number:', im
         write(0,*) ' A combined mode cannot contain parametric and non-parametric coordinates'
      endif
      do jel = 1, nstates
         iP = ipsiSPF(im,jel)
         do j = 1, nSPF(im,jel)
            psi(iP + j - 1) = cOne
            iP = iP + ngpt
         enddo
      enddo
      cycle
   endif
   !!!
   ! Procedure for the other DVR modes
   !!!
   allocate(grid_op(ngpt,nmD), psi_t(ngpt,maxval(nSPF(im,:))), psi2_t(ngpt))
   do j = 1, ngpt
      cdummy = cOne
      do k = 1, nmD
         kg = grid_index(k,j,i)
         cdummy = cdummy * phi(kg,iDOF(k,im))
      enddo
      psi_t(j,1) = cdummy
   enddo
   ! Normalize
   dummy = dznrm2(ngpt, psi_t(1,1), 1)
   call zdrscl(ngpt, dummy, psi_t(1,1), 1)
   ! Add further SPFs using Gram-Schmidt orthogonalisation:
   !    Subsequent SPFs are obtained by multiplying previous SPFs
   !    by the 'q' operators of the different dofs
   !
   ! Define the operator
   grid_op = dOne
   do k = 1, nmD
      do j = 1, ngpt
         kg = grid_index(k,j,i)
         do i2 = 1, ndofDVR
            if (idofDVR(i2) .eq. iDOF(k,im)) exit
         enddo
         grid_op(j,k) = merge(cos(grid(kg,i2)), grid(kg,i2), &
                              dvrtype(idofDVR(i2)) == 'Leg')
         !grid_op(j,k) = grid(kg,i2)
      enddo
   enddo
   ! Initialize the remaining SPFs
   kg = 0
   l  = 1  ! reference SPF to apply the operator
   idummy = 0
   jj = 1
   do j = 2, maxval(nSPF(im,:))
      jj = jj + 1
      kg_save = kg
      kg = mod(kg, nMD) + 1
      l  = (jj - 2) / nMD + 1
      ! Copy the SPF phi_l into the phi_j
      call zcopy(ngpt, psi_t(1,l), 1, psi_t(1,j), 1)
      ! Define the operator which multiplies phi_j
      do i2 = 1, ndofDVR
         if (idofDVR(i2) .eq. iDOF(kg,im)) exit
      enddo
      select case(dvrtype(iDOF(kg,im)))
         case('Leg')
            forall (k = 1:ngp(iDOF(kg,im))) grid_op1D(k) = cos(grid(k,i2))
         case('2pi')
            forall (k = 1:ngp(iDOF(kg,im))) grid_op1D(k) = merge(cos(grid(k,i2)), &
                                                                 sin(grid(k,i2)), idummy == 0)
         case default
            forall (k = 1:ngp(iDOF(kg,im))) grid_op1D(k) = grid(k,i2)
      end select
      ! Multiply
      iP = 1
      do i2 = 1, nl_psi(kg,im)
         do k = 1, ngp(iDOF(kg,im))
            call zdscal(ll_psi(kg,im),grid_op1D(k),psi_t(iP,j),1)
            iP = iP + ll_psi(kg,im)
         end do
      end do
      ! For the 2pi-DVR this step has to be repeated twice
      if (dvrtype(iDOF(kg,im)) == '2pi' .and. idummy == 0) then
         kg = kg_save
         jj = jj - 1
         idummy = 1
      end if
   end do
   !
   ! Gram-Schmidt orthogonalization
   !
   ! (1) Project out the first SPF from the other ones
   do j = 2, nS
      cdummy = - zdotc(ngpt, psi_t(1,1), 1, psi_t(1,j), 1)
      call zaxpy(ngpt, cdummy, psi_t(1,1), 1, psi_t(1,j), 1)
   end do
   ! (2) QR decomposition to orthogonalize the virtual SPFs
   call zgeqrf(ngpt,nS - 1, psi_t(1,2),ngpt,tau(1),WORK,LWORK,info)
   call zungqr(ngpt,nS - 1, nS - 1, psi_t(1,2),ngpt,tau(1),WORK,LWORK,info)
   !
   do jel = 1, nstates
      iP = ipsiSPF(im,jel)
      call zcopy(ngpt * nSPF(im,jel), psi_t(1,1), 1, psi(iP), 1)
   enddo
   !
!!!   iP = ipsiSPF(im,1)
!!!   call zgemm('C','N',nSPF(im,1),nSPF(im,1),ngpt,cOne,psi(iP),ngpt,psi(iP),ngpt,cZero,ovl,100)
!!!   write(0,*) nSPFmax
!!!   write(0,*)
!!!   write(0,*) ovl(1,1), ovl(1,2), ovl(1,3)
!!!   write(0,*) ovl(2,1), ovl(2,2), ovl(2,3)
!!!   write(0,*) ovl(3,1), ovl(3,2), ovl(3,3)
!!!   iP = ipsiSPF(im,2)
!!!   do j = 1, ngpt
!!!      write(2222,*) j, real(psi(iP + j - 1)), aimag(psi(iP + j - 1)), &
!!!              real(psi(iP + j - 1 + ngpt)), aimag(psi(iP + j - 1 + ngpt)), &
!!!              real(psi(iP + j - 1 + ngpt*2)), aimag(psi(iP + j - 1 + ngpt*2)), &
!!!              real(psi(iP + j - 1 + ngpt*3)), aimag(psi(iP + j - 1 + ngpt*3)), &
!!!              real(psi(iP + j - 1 + ngpt*4)), aimag(psi(iP + j - 1 + ngpt*4)), &
!!!              real(psi(iP + j - 1 + ngpt*5)), aimag(psi(iP + j - 1 + ngpt*5)), &
!!!              real(psi(iP + j - 1 + ngpt*6)), aimag(psi(iP + j - 1 + ngpt*6))
!!!   end do
   !
   deallocate(grid_op, psi_t, psi2_t)
   write(iout_unit,'(a,i0,a,20(1x,i3))') '   Mode ', im, '. Degrees of freedom: ', iDOF(1:nmD,im)
   write(iout_unit,'(a)') '      Discrete Variable Representation'
enddo
deallocate(grid_index)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! 2) Modes described by GWP distributions
do i = 1, nmodGWP
   im = imodGWP(i)
   nmD = nmodeDOF(im)
   allocate(xi_values(-nSPFmax:nSPFmax,nmD))
   allocate(w_values(-nSPFmax:nSPFmax,nmD))
   do j = 1, nmD
      k = iDOF(j,im)
      do kg = 1, ndofGWP
         if (idofGWP(kg) .eq. k) exit
      enddo
      ! define the widths
      ! ESEGUIRE UN CHECK: 2 * scw / (1 + scw ** 2) > overlap  (CREDO)
      dx = ParGWP(1,kg)
      scw = ParGWP(3,kg)
      forall (l = -nSPFmax:nSPFmax) w_values(l,j) = dx * scw ** abs(l)
      ! define the xi's
      !xi_values(0,j) = phi(1,k)
      xi_values(0,j) = dZero
      do l = 1, nSPFmax
         dummy = - log(ParGWP(2,kg)) + dOneHalf * log(2 * scw / (dOne + scw ** 2))
         dummy = dummy * (dOne + scw ** 2) / w_values(l - 1,j) ** 2
         dummy = (real(xi_values(l - 1,j)) + sqrt(dummy)) / scw ** 2
         xi_values( l,j) =   dummy
         xi_values(-l,j) = - dummy
      enddo
      if (abs(phi(1,k)) > epsilon(dOne)) &
         forall (l = -nSPFmax:nSPFmax) xi_values(l,j) = xi_values(l,j) + phi(1,k)
      cycle
      !
      ! (temporary) equidistant grid in the coordinate space
      dx = ParGWP(1,kg)
      w_values(:,j) = dx
      ! difference between xi values
      dummy = sqrt2 / dx * sqrt(- log(ParGWP(2,kg)))
      ! set the grid
      forall (l = -nSPFmax:nSPFmax) xi_values(l,j) = phi(1,k) + dble(l) * dummy
   enddo
   ! Set the initial parameters of the Gaussians
      ! Check that all DOFs are associated with the same distribution
   k = iDOF(1,im)
   do kg = 1, ndofGWP
      if (idofGWP(kg) .eq. k) exit
   enddo
   iDG = iDistGWP(kg)
   do j = 2, nmD
      k = iDOF(j,im)
      do kg = 1, ndofGWP
         if (idofGWP(kg) .eq. k) exit
      enddo
      if (iDistGWP(kg) .ne. iDG) then
         write(0,*) 'Input error:'
         write(0,*) '   the same GWP distribution must be used'
         write(0,'(a,i2)') '   for all degrees of freedom of mode ', im
         stop
      endif
   enddo
      !  Calculate lattice indices
   nS = maxval(nSPF(im,:))
   allocate(index_lat(nmD,nS))
   index_lat = 0
   if (nS .gt. 1) then   ! start the procedures only if there are more than one Gaussian
      select case(iDG)
         case(0)
            i1 = 1
            i2 = 1
            l = 2
            do
               do i3 = i1, i2
                  do j = 1, nmD
                     do i4 = -1, 1, 2
                        index_lat(:,l) = index_lat(:,i3)
                        index_lat(j,l) = index_lat(j,l) + i4
                        ! check that the index sequence is not existing already
                        j2 = 0
                        do j1 = 1, l - 1
                           if(all(index_lat(:,l) .eq. index_lat(:,j1))) then
                              j2 = 1
                              exit
                           endif
                        enddo
                        if (j2 .eq. 1) cycle   ! cycle because the index sequence was already existing
                        l = l + 1
                        if (l .gt. nS) exit
                     enddo
                     if (l .gt. nS) exit
                  enddo
                  if (l .gt. nS) exit
               enddo
               i1 = i2 + 1
               i2 = l - 1
               if (l .gt. nS) exit
            enddo
         case(1)   ! distribution read from external file
            write(GWPgridfile,'(a,i0)') 'gwpgrid_', im
            inquire(file = GWPgridfile, exist = exists)
            if (.not. exists) then
               write(0,*) 'INPUT ERROR: Gaussian grid file ', GWPgridfile, ' not found'
               stop
            endif
            io_grid = freeunit()
            open(unit = io_grid, file = GWPgridfile, status = 'old', form = 'formatted', action = 'read')
            do l = 1, nS
               read(io_grid,*) index_lat(1:nMD,l)
            enddo
            close(io_grid)
         case default
            write(0,*) 'Input error:'
            write(0,'(a,i3)') '   Unknown distribution type:', iDG
            stop
      end select
   endif
      !
   allocate(psi_t(nmD,nS))
   forall(l = 1:nS, j = 1:nmD)
      psi_t(j,l) = xi_values(index_lat(j,l),j)
      GWPa(j,l,i) = - (dOneHalf / w_values(index_lat(j,l),j)) ** 2
   end forall
   do jel = 1, nstates
      iP = ipsiSPF(im,jel)
      call zcopy(nMD * nSPF(im,jel), psi_t(1,1), 1, psi(iP), 1)
   enddo
   !
   deallocate(psi_t,index_lat)
   deallocate(xi_values, w_values)
   write(iout_unit,'(a,i0,a,20(1x,i3))') '   Mode ', im, '. Degrees of freedom: ', iDOF(1:nmD,im)
   write(iout_unit,'(a)') '      GWP representation'
enddo
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine SetSPF



subroutine OrthoWP(A0,psi0)
! Replace the wave packets with the population-weighted natural states
use sysparam
use globals
use psidef
implicit none
double complex, dimension(nAConftot * npackets), intent(inout) :: A0
double complex, dimension(dimpsi), intent(inout) :: psi0
integer :: iA, jA, i, j, info
double precision, dimension(:) :: RWORK(3 * npackets - 2)
double complex, dimension(:) :: WORK(3 * npackets)
double complex, dimension(nAConftot * npackets) :: A0tmp
double complex, dimension(:,:) :: ovl(npackets,npackets)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Calculate the overlap matrix between the wave packets
!write(0,*)
!write(0,*) 'Before'
do i = 1, npackets
   iA = APacketOffset(i) + 1
   do j = i, npackets
      jA = APacketOffset(j) + 1
      call WF_Ovl(A0(iA),psi0, A0(jA),psi0, ovl(i,j))
   end do
!   write(0,'(i3,3x,f11.8,1x,f11.8)') i, ovl(i,i)
end do
! Calculate the natural populations
!write(0,*) 'Eigenvalues'
call zheev('V','U',npackets,ovl,npackets,WPpop,WORK,3 * npackets, RWORK, info)
!do i = 1, npackets
!   write(0,'(i3,3x,f11.8,1x,f11.8)') i, WPpop(i)
!end do
! Transform the coefficients
call zgemm('N','N',nAConftot,npackets,npackets,cOne, &
           A0,nAconftot, ovl,npackets, cZero, A0tmp,nAConftot)
! Replace the original coefficients with the transformed ones
call zcopy(nAConftot * npackets, A0tmp,1, A0,1)
return
!
! Check
!write(0,*) 'After'
do i = 1, npackets
   iA = APacketOffset(i) + 1
   do j = i, npackets
      jA = APacketOffset(j) + 1
      call WF_Ovl(A0(iA),psi0, A0(jA),psi0, ovl(i,j))
   end do
!   write(0,'(i3,3x,f11.8,1x,f11.8)') i, ovl(i,i)
end do
!stop
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine OrthoWP



subroutine write_Gau_pos
! It writes the initial positions of the Gaussians in the log file
use sysparam
use globals
use psidef
implicit none
integer :: i, j, k, im, iP, nmD, iel, ielMax
double precision, dimension(:) :: qq(nDOF), pp(nDOF)
double complex, dimension(:) :: xi(nDOF)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
write(ilog_unit,*)
do i = 1, nmodGWP
   im = imodGWP(i)
   nmD = nmodeDOF(im)
   write(ilog_unit,'(a,i0)') 'Initial positions and exponents of the Gaussians of mode ', im
   write(ilog_unit,'(a,6x,20(8x,i3,11x))') '   DOF', iDOF(1:nmD,im)
   write(ilog_unit,'(a,20(6x,a,11x))') '   SPF no.    ', ('(q,A)', j = 1, nmD)
   ! Find the electronic state having more Gaussians
   ielMax = 1
   do iel = 2, nstates
      if (nSPF(im,iel) .gt. nSPF(im,ielMax)) ielMax = iel
   enddo
   iP = ipsiSPF(im,ielMax)
   do j = 1, nSPF(im,ielMax)
      xi(1:nmD) = psi(iP:iP + nmD - 1)
      pp(1:nmD) = imag(xi(1:nmD))
      qq(1:nmD) = dOneHalf * real(xi(1:nmD)) / GWPa(1:nmD,j,i)
      !write(ilog_unit,'(3x,i3,4x,20(3x,a,f8.4,a,f8.4,a))') &
      !   j, ('(',qq(k),',', pp(k),')', k = 1, nmD)
      write(ilog_unit,'(3x,i3,4x,20(3x,a,f8.4,a,f8.4,a))') &
         j, ('(',qq(k),',', GWPa(k,j,i),')', k = 1, nmD)
      iP = iP + nMD
   enddo
enddo
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine write_Gau_pos
