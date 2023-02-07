module Hamilton
   use sysparam
   use globals
   use psidef
   use dvr
   use storage
   use timingmod
   implicit none
   private
   public :: nHam, ae, iOpDof, lhs_el, rhs_el, h1mode, Dh1modeD, GauInt1D, Gau10Int1D, &
             nOper, iOper, Oper, nOpMax, psiMat, IsIdentity, SkipTerm, psi_sect, nPsiMat, &
             ll_psi, nl_psi
   public :: nHamTD, iHamTD, typeHamTD, parHamTD, nOpParam, OpParam, OpParamFunc
   public :: alloc_h1mode, alloc_IsIdentity, alloc_MeanField, alloc_rhom1MF
   public :: Calc_h1mode, Calc_DHD
   public :: MatricisePsi, BackMatricisePsi, MultTens, MultTensDiag, MultTensProj
   public :: nOpTerms, nExtraOpMax, nExtraOpTermsMax, nExtraOp, nExpect, iExpect, &
             lhs_elExtra, rhs_elExtra, aeExtra, iOpDofExtra, iOperExtra
   ! Dissipative terms
   public :: nDiss1, nDiss2, iDiss1, iDiss2, aeD1, aeD2
   ! Harmonic oscillator bath
   public :: lHObath, HObathPTorder, nHObath, wBath, iHObathQ, iHObathP
   !
   integer :: nHam   !  number of Hamiltonian terms
   integer :: nOpMax
   integer, dimension(:), allocatable :: lhs_el, rhs_el    ! xhs_el(i) = left/right el. state label of the i-th Hamiltonian term
   integer, dimension(:,:), allocatable :: iOpDof          ! iOpDof(i,j) = primitive operator for the i-th dof and the j-th Ham. term
                                                                  !  0 : identity
                                                                  !  i : q^i
                                                                  ! -2 : dq^2
   double complex, dimension(:), allocatable :: ae   ! ae(i) = coefficient of the i-th Hamiltonian term 
   double complex, dimension(:,:,:,:), allocatable :: h1mode ! h1mode(i,j,k,l) = matrix element between the SPFs i,j of the operator k of the mode l
                                                             !    The Hamiltonian term l defines the electronic states associated with i,j
   double complex, dimension(:,:,:,:), allocatable :: Dh1modeD    ! D^h * H * D (used only with the CMF2 scheme)
   logical, dimension(:,:), allocatable :: IsIdentity   ! IsIdentity(i,j) is .true. when the operator for the mode i and the Ham. term j is the identity
   integer, dimension(:), allocatable :: nOper          ! nOper(i) = number of operators for the mode i
   integer, dimension(:,:), allocatable :: iOper        ! iOper(i,j) = operator for the mode i in the j-th Hamiltonian term
   integer, dimension(:,:,:), allocatable :: Oper       ! Oper(i,j,k) = primitive operators for the j-th operator of the mode k.
                                                        !               The positions i = -1 and i = 0 contain the electronic labels
                                                        !               The positions i = 1:nmodeDOF(k) contain the primitive operator label (same as in iOpDof)
   target :: Oper
   ! For parameter-dependent operators
   integer :: nOpParam        ! number of parameterized elementary operators
   integer, dimension(:), allocatable :: OpParamFunc  ! OpParam(i) = parametrized function for the i-th parametrized operator
                                                      !    -  0: Gaussian, exp(p1*x + p2*x**2)
                                                      !    -  1: interpolation 
   double precision, dimension(:,:), allocatable :: OpParam      ! OpParam(i,j) = i-th parameter for the j-th parametrized operator      
   ! For time-dependent Hamiltonian terms
   integer :: nHamTD
   integer, dimension(1000) :: iHamTD                        ! iHamTD(i) = i-th time-dependent Hamiltonian term 
   integer, dimension(1000) :: typeHamTD                     ! typeHamTD(i) = time-dependent function for the i-th time-dependent Hamiltonian term
   double precision, dimension(0:4,1000) :: parHamTD                ! parHamTD(j,i) = j-th parameter defining the i-th time-dependent Hamiltonian term
   ! Auxiliary matrices, allocated in alloc_h1mode
   integer, dimension(:,:), allocatable :: ll_psi, nl_psi
   double complex, dimension(:), allocatable :: psi_sect
   integer :: nPsiMat
   double complex, dimension(:,:), allocatable :: psiMat
   ! 
   logical, dimension(:,:), allocatable :: SkipTerm            ! dictates whether the term should be skipped in mean field calculation
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Additional operators (for time-dependent expectation values, for operating on the initial state, etc.)
   integer, parameter :: nExtraOpMax = 800, nExtraOpTermsMax = 8000
   integer :: nExtraOp                          ! number of extra operators 
   integer, dimension(nExtraOpMax * 3) :: nOpTerms  ! nOpTerms(i) = number of terms in the i-th extra operator
   integer, dimension(nExtraOpTermsMax,nExtraOpMax) :: lhs_elExtra, rhs_elExtra 
   integer, dimension(:,:,:), allocatable :: iOpDofExtra
   double complex, dimension(nExtraOpTermsMax,nExtraOpMax) :: aeExtra  ! aeExtra(i,j) = coefficient of the operator term i for the operator j
   integer, dimension(:,:,:), allocatable :: iOperExtra  ! iOperExtra(i,j,k) = operator for the mode i, of the operator term j for the operator k
   integer :: nExpect    ! number of operators for which the expectation value should be calculated
   integer, dimension(nExtraOpMax) :: iExpect   ! iExpect(i) = index of the i-th operator for which the expectation value should be calculated 
   !
   ! Dissipation
   !
   integer :: nDiss1, nDiss2    ! Number of dissipation operators of type 1 and 2
   integer, dimension(nExtraOpMax) :: iDiss1, iDiss2   ! iDiss1(i) = index of the i-th dissipation operator of type 1 
   double precision, dimension(:), allocatable :: aeD1, aeD2       ! global coefficients of the dissipative superoperators (real)
   ! Harmonic perturbative bath
   logical :: lHObath         ! is there a harmonic bath?
   integer :: HObathPTorder, nHObath   ! perturbative order
   integer, dimension(nExtraOpMax) :: iHObathQ, iHObathP
   double precision, dimension(:), allocatable :: wBath    ! Bath frequency
! 
contains



!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Functions and subroutines related to allocation of arrays
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!

   subroutine alloc_h1mode
   ! It allocates the array h1mode
   ! The matrices h1mode(:,:,i,j) correspondent to the identity are set here,
   !    in order to avoid to set these matrices every time Calc_h1mode is called
   implicit none
   integer :: a_status
   integer :: iHam, i, j, im, nMD, iel, nS, nG, nn1, nn2, &
              kdof, k, nOp, nn, nmDmax, iOp, iT, iFile, nLines, idummy
   integer, dimension(:), allocatable :: OpVec
   integer, dimension(:,:,:), allocatable :: Oper_t
   double precision :: p1, p2
   logical :: inserted
   ! For 1D-interpolants
   double precision, dimension(:), allocatable :: xdummy, ydummy
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
!  Define the 'operator basis' for each mode
   nmDmax = maxval(nmodeDOF)
   allocate(nOper(nmodes),iOper(nmodes,nHam))
   allocate(Oper_t(-1:nmDmax,nHam + sum(nOpTerms(1:nExtraOp)),nmodes))
   allocate(OpVec(-1:nmDmax))
   allocate(iOperExtra(nmodes,nExtraOpTermsMax,nExtraOpMax))
   do im = 1, nmodes
      ! Starting index in the arrays Oper_t and OpVec to check whether an operator has already inserted
      nmD = nmodeDOF(im)
      nOp = 0   ! This counter corresponds to the number of operators for the mode im
      do iHam = 1, nHam
         if (MultiSet(im)) then
            OpVec(-1) = lhs_el(iHam)    ! different operators for different electronic states
            OpVec(0) = rhs_el(iHam)
         else
            OpVec(-1:0) = 1             ! the same operator for all electronic states
         end if
         OpVec(1:nMD) = iOpDOF(iDOF(1:nMD,im),iHam)
         ! check whether the operator is already present
         inserted = .false.
         do i = 1, nOp
            if ( all(Oper_t(-1:nmD,i,im) .eq. OpVec(-1:nmD)) ) then
               inserted = .true.
               iOper(im,iHam) = i
               exit
            end if
         end do
         if (.not. inserted) then
            nOp = nOp + 1
            Oper_t(-1:nmD,nOp,im) = OpVec(-1:nmD)
            iOper(im,iHam) = nOp
         end if
      end do
      ! Extra operators
      do iOp = 1, nExtraOp
         do iT = 1, nOpTerms(iOp)
            if (MultiSet(im)) then
               OpVec(-1) = lhs_elExtra(iT,iOp)    ! different operators for different electronic states
               OpVec(0) = rhs_elExtra(iT,iOp) 
            else 
               OpVec(-1:0) = 1                    ! the same operator for all the electronic states
            end if
            OpVec(1:nMD) = iOpDOFExtra(iDOF(1:nMD,im),iT,iOp)
            ! check whether the operator is already present
            inserted = .false.
            do i = 1, nOp
               if ( all(Oper_t(-1:nmD,i,im) .eq. OpVec(-1:nmD)) ) then
                  inserted = .true.
                  iOperExtra(im,iT,iOp) = i
                  exit
               end if
            end do
            if (.not. inserted) then
               nOp = nOp + 1
               Oper_t(-1:nmD,nOp,im) = OpVec(-1:nmD)
               iOperExtra(im,iT,iOp) = nOp
            end if
         end do
      end do
      !
      nOper(im) = nOp
   end do
   nOpMax = maxval(nOper)
   allocate(Oper(-1:nmDmax,nOpMax,nmodes))
   ! copy the values of Oper_t into Oper, and deallocate Oper_t
   do im = 1, nmodes
      nMD = nmodeDOF(im)
      nOp = nOper(im)
      Oper(-1:nMD,1:nOp,im) = Oper_t(-1:nMD,1:nOp,im) 
   end do
   deallocate(Oper_t)
!
! Set up the array gridPow, whcih contains the values of q^n on the grid points
!
   ! Get the maximum power needed
   nn = 0
   do i = 1, nmodDVR
      im = imodDVR(i)
      nMD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         do k = 1, nMD
            if (Oper(k,iHam,im) .le. 20 .and. Oper(k,iHam,im) .gt. nn) nn = Oper(k,iHam,im)
         end do
      end do
   end do
   allocate(gridPow(ngpMax,0:nn,ndofDVR))
   gridPow = dOne
   do i = 1, nmodDVR
      im = imodDVR(i)
      nMD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         do k = 1, nMD
            nOp = Oper(k,iHam,im)
            if (nOp .le. 0 .or. nOp .gt. 200) cycle     ! (negative powers not yet implemented)
            kdof = iDOF(k,im)
            if (nOp .eq. 1) then
               call dcopy(ngp(kdof),grid(1,indDVR(kdof)),1,gridPow(1,nOp,indDVR(kdof)),1)
            else
               gridPow(1:ngp(kdof),nOp,indDVR(kdof)) = grid(1:ngp(kdof),indDVR(kdof)) ** nOp
            end if
         end do
      end do
   end do
!
! Set up the array gridOper, which contains the value of the parametrized operator on the DVR grid
!
   allocate(gridOper(ngpMax,nOpParam,ndofDVR))
   gridOper = dZero
   do i = 1, nmodDVR
      im = imodDVR(i)
      nMD = nmodeDOF(im)
      do iHam = 1, nOper(im)     
         OpVec(-1:nMD) = Oper(-1:nMD,iHam,im)
         do k = 1, nMD
            iOp = OpVec(k)
            if (iOp .le. 200 .and. iOp .gt. -200) cycle
            if (iOp .gt. 200)  iOp =   iOp - 200
            if (iOp .le. -200) iOp = - iOp - 200
            kdof = iDOF(k,im)
            select case(OpParamFunc(iOp))
               case(0)   ! Gaussian exp(p1*x + p2*x**2)
                  p1 = OpParam(1,iOp)
                  p2 = OpParam(2,iOp)
                  gridOper(1:ngp(kdof),iOp,indDVR(kdof)) &
                     = exp(p1 * grid(1:ngp(kdof),indDVR(kdof)) + p2 * grid(1:ngp(kdof),indDVR(kdof)) ** 2)
               case(1)   ! One-dimensional interpolation
                  iFile = int(OpParam(1,iOp))
                  iT = freeunit()
                  open(unit = iT, file = trim(adjustl(strings(iFile))), form = 'formatted', status = 'old', action = 'read')
                  ! Get the number of lines
                  nLines = 0
                  idummy = 0
                  do 
                     read(iT,*, iostat = idummy) 
                     if (idummy .ne. 0) exit
                     nlines = nlines + 1
                  end do
                  ! Read the data
                  allocate(xdummy(nLines), ydummy(nLines))
                  rewind(iT)
                  do j = 1, nLines
                     read(iT,*) xdummy(j), ydummy(j)
                  end do
                  ! Interpolation
                  call Spline(nLines, xdummy, ydummy, ngp(kdof), grid(1,indDVR(kdof)), gridOper(1,iOp,indDVR(kdof)))
                  deallocate(xdummy, ydummy)
                  close(iT)
               case(2)   ! indicator function for the interval [x0-dx/2,x0+dx/2]
                  p1 = OpParam(1,iOp)
                  p2 = OpParam(2,iOp)
                  do j = 1, ngp(kdof)
                     if (abs(grid(j,indDVR(kdof)) - p1) .lt. dOneHalf * p2) then
                        gridOper(j,iOp,indDVR(kdof)) = dOne
                     else
                        gridOper(j,iOp,indDVR(kdof)) = dZero
                     endif
                  enddo
               case(3)   ! Hermite function x**r * exp(p1 * x + p2 * x**2)
                  p1 = OpParam(1,iOp)
                  p2 = OpParam(2,iOp)
                  idummy = int(OpParam(3,iOp))
                  do concurrent (j = 1:ngp(kdof))
                     gridOper(j,iOp,indDVR(kdof)) = grid(j,indDVR(kdof)) ** idummy &
                                                  * exp(p1 * grid(j,indDVR(kdof)) + p2 * grid(j,indDVR(kdof))**2)
                  end do
               case(4)   ! Powers of sine ( sin(p1 * x)^p2 )
                  p1 = OpParam(1,iOp)
                  p2 = OpParam(2,iOp)
                  do concurrent (j = 1:ngp(kdof))
                     gridOper(j,iOp,indDVR(kdof)) = sin(p1 * grid(j,indDVR(kdof))) ** p2
                  end do
               case(5)   ! Powers of cosine ( cos(p1 * x)^p2 )
                  p1 = OpParam(1,iOp)
                  p2 = OpParam(2,iOp)
                  do concurrent (j = 1:ngp(kdof))
                     gridOper(j,iOp,indDVR(kdof)) = cos(p1 * grid(j,indDVR(kdof))) ** p2
                  end do
               case(6)   ! Projector onto a 1D eigenfunction
                  iFile = int(OpParam(1,iOp))
                  ! Evaluation of the eigenfunction
                  call Eigenf1DPot_real(ngp(kdof), grid(1,indDVR(kdof)), &
                                        OpParam(2,iOp), d2DVR(1,1,indDVR(kdof)),ngpMax, &
                                        strings(iFile), int(OpParam(3,iOp)), gridOper(1,iOp,indDVR(kdof)))
            end select
         end do    
      end do
   end do
!
   allocate(h1mode(nSPFmax,nSPFmax,nOpMax,nmodes), stat = a_status)
   if (a_status .ne. 0) call err_alloc('h1mode', 'alloc_h1mode', a_status)
   h1mode = cZero
! Initialise the identity matrices for the non-GWP modes
   do i = 1, nmodDVR
      im = imodDVR(i)
      nmD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         if (Oper(-1,iHam,im) .ne. Oper(0,iHam,im)) cycle  ! different electronic states
         if (any(Oper(1:nmD,iHam,im) .ne. 0)) cycle
         iel = Oper(-1,iHam,im)
         do j = 1, nSPF(im,iel)
            h1mode(j,j,iHam,im) = cOne
         end do
      end do
   end do
! Allocate the auxiliary matrix psiMat
   nn1 = 0
   nn2 = 0
   do i = 1, nmodDVR 
      im = imodDVR(i)
      nG = ngpm(i)
      nS = maxval(nSPF(im,:))
      do k = 1, nmodeDOF(im)
         kdof = iDOF(k,im)
         nn1 = max(nn1,ngp(kdof))
         nn2 = max(nn2,nG * nS / ngp(kdof))
      end do
   end do   
   allocate(psiMat(nn2,nn1))
   nPsiMat = nn2
! Initialise the auxiliary arrays ll_psi and nl_psi
   nMD = 0
   do i = 1, nmodDVR
      nMD = max(nMD,nmodeDOF(imodDVR(i)))
   end do
   allocate(ll_psi(nMD,nmodes), nl_psi(nMD,nmodes))
   do i = 1, nmodDVR
      im = imodDVR(i)
      nMD = nmodeDOF(im)
      do k = 1, nMD
         ll_psi(k,im) = product(ngp(iDOF(1:k - 1,im)))
         nl_psi(k,im) = product(ngp(iDOF(k + 1:nMD,im)))
      end do
   end do 
! Allocate psi_sect
   nn = 0
   do i = 1, nmodDVR
      im = imodDVR(i)
      nn = max(nn, ngpm(i) * maxval(nSPF(im,:)))
   end do
   allocate(psi_sect(nn))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   deallocate(OpVec) 
   return
   end subroutine alloc_h1mode



   subroutine alloc_IsIdentity
   ! Sets up the array IsIdentity
   implicit none
   integer :: a_status
   integer :: iHam, im, nMD
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   allocate(IsIdentity(nmodes,nOpMax), stat = a_status)
   if (a_status .ne. 0) call err_alloc('IsIdentity', 'alloc_IsIdentity', a_status)
   !
   IsIdentity = .false.
   do im = 1, nmodes
      nmD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         if (Oper(-1,iHam,im) .ne. Oper(0,iHam,im)) cycle   ! terms for different electronic states are not identity
         if (any(Oper(1:nmD,iHam,im) .ne. 0)) cycle
         IsIdentity(im,iHam) = .true.
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine alloc_IsIdentity



   subroutine alloc_MeanField
   implicit none
   integer :: a_status, i, j
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   allocate(MeanField(nSPFmax,nSPFmax,nOpMax,nmodes), stat = a_status)
   if (a_status .ne. 0) call err_alloc('MeanField', 'alloc_MeanField', a_status)
   allocate(Adot_tmp(nAconftot * npackets), stat = a_status)
   if (a_status .ne. 0) call err_alloc('Adot_tmp', 'alloc_MeanField', a_status)
   ! Allocate SkipTerm
   allocate(SkipTerm(nmodes,nHam))
   SkipTerm = .false.
   do j = 1, nHam
      do i = 1, nmodes
         if (dvrType(iDOF(1,i)) .ne. 'gwp' .and. IsIdentity(i,iOper(i,j))) &
            SkipTerm(i,j) = .true.
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine alloc_MeanField


   subroutine alloc_rhom1MF
   implicit none
   integer :: a_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   allocate(rhom1MF(nSPFmax,nSPFmax,nOpMax,nmodDVR), stat = a_status)
   if (a_status .ne. 0) call err_alloc('rhom1MF', 'alloc_rhom1MF', a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine alloc_rhom1MF



!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Functions and subroutines related to calculation of Hamiltonian matrix elements and mean fields
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!


   subroutine Calc_h1mode(psi0, im)
   ! It calculates the array h1mode which contains the matrix elements 
   ! of the operators for the mode im between SPFs or GWPs
   ! - It is assumed that the (intra-state) overlap matrices S00 have been already calculated
   implicit none
   integer, intent(in) :: im   ! the mode
   integer :: iHam, iel1, iel2, i, nMD, iD, &
              nS1, nS2, nG, k, iOD, j1, j2, iP1, iP2
   integer, dimension(:), pointer :: OpVec
   double complex, dimension(dimpsi), intent(in) :: psi0
   double complex, dimension(:) :: gInt(nDOF)
   double complex, external :: zdotc
   logical :: EqualStates
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Start timer
   call continue_timer(5)
   !
   if (dvrtype(iDOF(1,im)) .ne. 'gwp') then
   ! DVR modes
      do i = 1, nmodDVR
         if (imodDVR(i) == im) exit
      end do
      nMD = nmodeDOF(im)
      nG = ngpm(i)
      do iHam = 1, nOper(im)   
         if (IsIdentity(im,iHam)) cycle   ! In this case h1mode is set up in the initialization (subroutine alloc_h1mode)
         OpVec(-1:nMD) => Oper(-1:nMD,iHam,im) 
         iel1 = OpVec(-1)
         iel2 = OpVec(0)
         nS1 = nSPF(im,iel1)
         nS2 = nSPF(im,iel2)
         iP2 = ipsiSPF(im,iel2)
         call zcopy(nS2 * nG,psi0(iP2),1,psi_sect,1)     ! copy the section of psi0
                                                         ! corresponding to the SPF of the mode im for the state iel2
         ! Operate on the SPFs with the individual DOF operators
         do k = 1, nMD
            iOD = OpVec(k)
            if (iOD .eq. 0) cycle   ! Do not operate with the identity
            iD = indDVR(iDOF(k,im))
            select case(iOD)
               case(-1)
                  call MultTens(d1DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)   ! MultTens works also for non-symmetric matrices
                  call zscal(nS2 * nG,ciOne,psi_sect(1),1)  ! + i because MultTens multiplies from the right
               case(-2)  ! dq^2
                  call MultTens(d2DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)
               case(-3)  ! - i d/dx * x - i x * d/dx
                  call MultTens(xd1DVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)   ! MultTens works also for non-symmetric matrices
                  call zscal(nS2 * nG,ciOne,psi_sect(1),1)  ! + i because MultTens multiplies from the right
               case(-4)  ! sin^-2
                  call MultTens(operDVR(1,1,iD),ngpmax,iel2,im,k,psi_sect(1),nS2 * nG)
               case(1:20)
                  call MultTensDiag(gridPow(1,iOD,iD),iel2,im,k,psi_sect(1))
               case(201:)
                  call MultTensDiag(gridOper(1, iOD - 200,iD),iel2,im,k,psi_sect(1))
               case(:-200)   ! 1D-projector
                  call MultTensProj(gridOper(1,-iOD - 200,iD),iel2,im,k,psi_sect(1),nS2 * nG)
            end select
         end do
         !
         if (iel1 .eq. iel2) then   ! h1mode is Hermitian
            call zher2k('U','C',nS2,nG,cOneHalf, psi0(iP2),nG, psi_sect(1),nG, &
                        dZero,h1mode(1,1,iHam,im),nSPFmax)
         else
            iP1 = ipsiSPF(im,iel1)
            call zgemm('C','N',nS1,nS2,nG,cOne, psi0(iP1),nG, psi_sect(1),nG, &
                       cZero,h1mode(1,1,iHam,im),nSPFmax)
         end if
         ! Make the very small values equal to zero (perhaps it enforces symmetry ?) 
         !!!do j2 = 1, nS2
         !!!   do j1 = 1, merge(j2, nS1, iel1 == iel2)
         !!!      h1mode(j1,j2,iHam,im) = merge(h1mode(j1,j2,iHam,im), &
         !!!                                    cZero, &
         !!!                                    abs(h1mode(j1,j2,iHam,im)) > epsilon(dOne))
         !!!   end do
         !!!end do
         !
      end do
   !
   else
   ! GWP modes
   do i = 1, nmodGWP
      if (imodGWP(i) == im) exit
   end do
      nMD = nmodeDOF(im)
      do iHam = 1, nOper(im)
         OpVec(-1:nMD) => Oper(-1:nMD,iHam,im) 
         iel1 = OpVec(-1)
         iel2 = OpVec(0)
         nS1 = nSPF(im,iel1)
         nS2 = nSPF(im,iel2)
         if (IsIdentity(im,iHam)) then
            call zlacpy('U',nS1,nS2,S00(1,1,i,iel1,iel2),nS00max,h1mode(1,1,iHam,im),nSPFmax)   ! copy the overlap matrix into h1mode
            cycle
         end if
         !
         EqualStates = iel1 .eq. iel2
         iP2 = ipsiSPF(im,iel2) - 1
         do j2 = 1, nS2
            iP1 = ipsiSPF(im,iel1) - 1
            do j1 = 1, merge(j2,nS1,EqualStates)  ! j2 if (iel1 .eq. iel2) else nS1   [h1mode is assumed Hermitian for iel1 .eq. iel2]
               do concurrent (k = 1:nMD)
                  gInt(k) = GauInt1D(GWPa(k,j1,i), GWPa(k,j2,i),  psi0(iP1 + k), psi0(iP2 + k), OpVec(k))
               end do
               h1mode(j1,j2,iHam,im) = S00(j1,j2,i,iel1,iel2) * product(gInt(1:nMD))
               iP1 = iP1 + nMD
            end do
            iP2 = iP2 + nmD
         end do
      end do
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Suspend timer
call suspend_timer(5)
!
   return
   end subroutine Calc_h1mode



   subroutine Calc_DHD(DD0)
   implicit none
   double complex, dimension(:,:,:,:), intent(in) :: DD0(nS00max,nS00max,nmodGWP,nstates)
   integer :: i, im, iel1, iel2, nS1, nS2, iHam, j1, j2
   double complex, dimension(:,:) :: aux(nS00max,nS00max)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   do i = 1, nmodGWP
      im = imodGWP(i)
      do iHam = 1, nOper(im)
         iel1 = Oper(-1,iHam,im)
         iel2 = Oper( 0,iHam,im)
         nS1  = nSPF(im,iel1)
         nS2  = nSPF(im,iel2)
         if (IsIdentity(im,iHam)) then  ! D^H * S * D = 1
            call zlaset('All', nS1, nS2, cZero, cOne, Dh1modeD(1,1,iHam,i), nS00max)
            cycle
         end if
         if (iel1 == iel2) then  ! Hermitian algebra
            call zhemm('L','U',nS1,nS2,cOne,h1mode(1,1,iHam,im),nSPFmax, &
                       DD0(1,1,i,iel2),nS00max,cZero,aux,nS00max)                     ! aux = H * D
            call zher2k('U','C',nS1,nS1,cOneHalf,DD0(1,1,i,iel1),nS00max, &
                        aux,nS00max,cZero,Dh1modeD(1,1,iHam,i),nS00max)             ! DHD = 1/2 * (D^H * aux + aux^H * D)
         else
            call zgemm('C','N',nS1,nS2,nS1,cOne,DD0(1,1,i,iel1),nS00max, &
                       h1mode(1,1,iHam,im),nSPFmax,cZero,aux,nS00max)                ! aux = D^H * H
            call zgemm('N','N',nS1,nS2,nS2,cOne,aux,nS00max, &
                       DD0(1,1,i,iel2),nS00max,cZero,Dh1modeD(1,1,iHam,i),nS00max)   ! DHD = aux * D
         end if
         ! Make the very small values equal to zero (perhaps it enforces symmetry ?) 
         do j2 = 1, nS2
            do j1 = 1, merge(j2, nS1, iel1 == iel2)
               Dh1modeD(j1,j2,iHam,i) = merge(Dh1modeD(j1,j2,iHam,i), cZero, abs(Dh1modeD(j1,j2,iHam,i)) > epsilon(dOne))
            end do
         end do
         !
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine Calc_DHD


   subroutine MultTens(AM,LDA,iel,im,k,psiV,npsi)
!  Multiplication along the DOF k of mode im between the 
!  tensor psi and the matrix AM
!  psiV = psiV *_k AM   (k-DOF tensor multiplication)
!  - psiV contains the SPFs for mode im (which includes - in principle - many combined DOFs)
!  - AM is a matrix which operates only on the DOF k
   implicit none
   integer, intent(in) :: LDA, iel, im, k, npsi
   integer :: iP, loop_length, ngpt, j, nSize
   double complex, dimension(*), intent(inout) :: psiV
   double complex, dimension(npsi) :: psiV2
   double precision, dimension(2 * npsi) :: RWORK 
   double precision, dimension(LDA,*), intent(in) :: AM
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   loop_length = ll_psi(k,im)   ! product(ngp(iDOF(1:k - 1,im)))
   ngpt = ngp(iDOF(k,im))
   nSize = loop_length * ngpt
   iP = 1
   do j = 1,nSPF(im,iel) * nl_psi(k,im)
      call zlacrm(loop_length,ngpt,psiV(iP),loop_length,AM,LDA,psiV2(iP),loop_length,RWORK)
      iP = iP + nSize
   end do
   call zcopy(npsi,psiV2,1,psiV,1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine MultTens


  
   subroutine MultTensDiag(AV,iel,im,k,psiV)
!  Multiplication along the DOF k of mode im between the 
!  tensor psi and the diagonal matrix AV
!  psiV2 = psiV1 *_k AM   (k-DOF tensor multiplication)
!  - psiV1 contains the SPFs for mode im (which includes - in principle - many combined DOFs)
!  - AV is a diagonal matrix which operates only on the DOF k
   implicit none
   integer, intent(in) :: iel, im, k
   integer :: iP, loop_length, ngpt, j, kg
   double complex, dimension(*), intent(inout) :: psiV
   double precision, dimension(*), intent(in) :: AV
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   loop_length = ll_psi(k,im)   ! product(ngp(iDOF(1:k - 1,im)))
   ngpt = ngp(iDOF(k,im))
   iP = 1
   do j = 1, nSPF(im,iel) * nl_psi(k,im)
      do kg = 1, ngpt
         call zdscal(loop_length,AV(kg),psiV(iP),1)
         iP = iP + loop_length
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine MultTensDiag


   subroutine MultTensProj(AV,iel,im,k,psiV,npsi)
!  Multiplication along the DOF k of mode im between the 
!  tensor psi and the projector matrix AV * AV^T
!  - psiV1 contains the SPFs for mode im (which includes - in principle - many combined DOFs)
!  - AV is the projector wavefunction
   implicit none
   integer, intent(in) :: iel, im, k, npsi
   integer :: iP, loop_length, ngpt, j, kg, nSize
   double complex, dimension(*), intent(inout) :: psiV
   double complex, dimension(npsi) :: psiV2
   double precision, dimension(2 * npsi) :: RWORK 
   double precision, dimension(*), intent(in) :: AV
   double precision, dimension(ngpmax,ngpmax) :: AM
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   loop_length = ll_psi(k,im)   ! product(ngp(iDOF(1:k - 1,im)))
   ngpt = ngp(iDOF(k,im))
   nSize = loop_length * ngpt
! Construct AV * AV^T
   call dlaset('All',ngpt,ngpt,dZero,dZero,AM,ngpmax)
   call dger(ngpt,ngpt,dOne,Av(1),1,AV(1),1, AM,ngpmax)
! Multiply
   iP = 1
   do j = 1,nSPF(im,iel) * nl_psi(k,im)
      call zlacrm(loop_length,ngpt,psiV(iP),loop_length,AM,ngpmax,psiV2(iP),loop_length,RWORK)
      iP = iP + nSize
   end do
   call zcopy(npsi,psiV2,1,psiV,1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine MultTensProj


   subroutine MatricisePsi(iel,im,k,psiV,psiM)
!  Matricisation of the SPF vector for the electronic state iel and the mode im
!  - psiV (input) is the complete psi vector
!  - psiM (output) is the psi vector for the mode im and the state iel,
!    in order that each column corresponds to a grid point for the k-th DOF of mode im
   implicit none
   integer, intent(in) :: iel, im, k
   integer :: iP, loop_length, nloops, ngpt, iC, j, iloop, kg, kk
   double complex, dimension(:), intent(in) :: psiV
   double complex, dimension(:,:), intent(out) :: psiM
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!        ! CONTROLLARE BENE TUTTA QUESTA SUBROUTINE
   loop_length = ll_psi(k,im)   ! product(ngp(iDOF(1:k - 1,im)))
   nloops = nl_psi(k,im)        ! product(ngp(iDOF(k + 1:nmD,im)))
   ngpt = ngp(iDOF(k,im))
   iC = 0
   iP = 0
   do j = 1, nSPF(im,iel)
      do iloop = 1, nloops
         ! compact form
         do concurrent (kg = 1:ngpt, kk = 1:loop_length) 
            psiM(iC + kk,kg) = psiV(iP + loop_length * (kg - 1) + kk)
         end do
         iP = iP + loop_length * ngpt
         ! extendend form
!         do kg = 1, ngpt
!            call zcopy(loop_length,psiV(iP + 1),1,psiM(iC + 1,kg),1)
!            !psiM(iC + 1:iC + loop_length,kg) = psiV(iP + 1:iP + loop_length)
!            iP = iP + loop_length
!         end do
         !
         iC = iC + loop_length
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine MatricisePsi



   subroutine BackMatricisePsi(iel,im,k,psiV,psiM)
   !subroutine BackMatricisePsi(iel,im,k,psiV,psiM,nn1,nn2)
!  It construct the psi vector psiV from its matricised form
!  (Inverse of MatricisePsi)
   implicit none
   integer, intent(in) :: iel, im, k
   integer :: iP, loop_length, nloops, ngpt, iC, j, iloop, kg, kk
   double complex, dimension(:), intent(out) :: psiV
   double complex, dimension(:,:), intent(in) :: psiM
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!        ! CONTROLLARE BENE TUTTA QUESTA SUBROUTINE
   loop_length = ll_psi(k,im)   ! product(ngp(iDOF(1:k - 1,im)))
   nloops = nl_psi(k,im)        ! product(ngp(iDOF(k + 1:nmD,im)))
   ngpt = ngp(iDOF(k,im))
   iC = 0
   iP = 0
!   iC = 1
!   iP = 1
   do j = 1, nSPF(im,iel)
      do iloop = 1, nloops
         ! compact form
         do concurrent (kg = 1:ngpt, kk = 1:loop_length) 
            psiV(iP + loop_length * (kg - 1) + kk) = psiM(iC + kk,kg)
         end do
         iP = iP + loop_length * ngpt
         ! extended form
!         do kg = 1, ngpt
!            call zcopy(loop_length,psiM(iC,kg),1,psiV(iP),1)
!            iP = iP + loop_length
!         end do
         iC = iC + loop_length
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine BackMatricisePsi



!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Functions and subroutines related to Gaussian integrals
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!

   pure recursive double complex function GauInt1D(A1,A2,xi1,xi2,iOp) result(gInt)
! It yields the integral <g_1|Operator|g_2> / <g_1|g_2>, where g_i is a Gaussian function:
! g_i = exp(A_i * x^2 + xi_i * x)
!    The operator is defined by the flag iOp
! WARNING: In this function, it is assumed that A_i < 0
   integer :: k, idummy
   integer, intent(in) :: iOp
   double precision, intent(in) :: A1, A2
   double precision ::  A1A2
   double complex, intent(in) :: xi1, xi2
   double precision :: p1, p2
   double complex :: xi1c, cdummy1, cdummy2
   double precision :: dummy
   double precision, parameter :: dOneFourth = 0.25d0
   ! The following table contains coefficients
   !    T(k,n) = Binomial(n,2k) * |(2k-1)!!|
   double precision, dimension(0:10,0:20), parameter :: &
      Coeff = reshape(   &    
            (/ &
             1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. , &
             1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. , &
             1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0. , &
             1.,3.,0.,0.,0.,0.,0.,0.,0.,0.,0. , &
             1.,6.,3.,0.,0.,0.,0.,0.,0.,0.,0. , &
             1.,10.,15.,0.,0.,0.,0.,0.,0.,0.,0. ,&
             1.,15.,45.,15.,0.,0.,0.,0.,0.,0.,0. ,&
             1.,21.,105.,105.,0.,0.,0.,0.,0.,0.,0. ,&
             1.,28.,210.,420.,105.,0.,0.,0.,0.,0.,0. ,&
             1.,36.,378.,1260.,945.,0.,0.,0.,0.,0.,0. ,&
             1.,45.,630.,3150.,4725.,945.,0.,0.,0.,0.,0. ,&
             1.,55.,990.,6930.,17325.,10395.,0.,0.,0.,0.,0. ,&
             1.,66.,1485.,13860.,51975.,62370.,10395.,0.,0.,0.,0. ,&
             1.,78.,2145.,25740.,135135.,270270.,135135.,0.,0.,0.,0. ,&
             1.,91.,3003.,45045.,315315.,945945.,945945.,135135.,0.,0.,0. ,&
             1.,105.,4095.,75075.,675675.,2837835.,4729725.,2027025.,0.,0.,0. ,&
             1.,120.,5460.,120120.,1351350.,7567560.,18918900.,16216200.,2027025.,0.,0. ,&
             1.,136.,7140.,185640.,2552550.,18378360.,64324260.,91891800.,34459425.,0.,0. ,&
             1.,153.,9180.,278460.,4594590.,41351310.,192972780.,413513100.,310134825.,34459425.,0. ,&
             1.,171.,11628.,406980.,7936110.,87297210.,523783260.,1571349780.,1964187225.,654729075.,0. ,&
             1.,190.,14535.,581400.,13226850.,174594420.,1309458150.,5237832600.,9820936125.,6547290750.,654729075. &
             /) , (/ 11, 21  /))
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   xi1c = conjg(xi1)
   select case(iOp)
      case(0)     ! identity
         gInt = cOne
      case(-1)    ! p = - i * dq
         gInt = - ciOne * (A1 * xi2 - A2 * xi1c) / (A1 + A2)
      case(-2)    ! dq^2
         A1A2 = A1 + A2
         gInt = ((A2 * xi1c - A1 * xi2) / A1A2) ** 2 + 2.d0 * A1 * A2 / A1A2
         !gInt = (A2 * xi1c - A1 * xi2) ** 2 + 2 * A1 * A2 * A1A2
         !gInt = gInt / A1A2 ** 2
      case(-3)    ! <g_i| qp + pq |g_j> = - i * <g_i| 4*A_j*q^2 + 2*xi_j*q + 1 |g_j>
         gInt = A1 ** 2 - A2 ** 2 &
              + (A2 * xi1c - A1 * xi2) * (xi1c + xi2)
         gInt = - ciOne * gInt / (A1 + A2) ** 2
      case(1)
         gInt = (xi1c + xi2) / (-2.d0 * (A1 + A2))
      case(2)
         A1A2 = 2.d0 * (A1 + A2)
         gInt = ((xi1c + xi2) / A1A2) ** 2 - dOne / A1A2
         !gInt = ((xi1c + xi2) / A1A2) ** 2 - cOne / A1A2
         !gInt = (xi1c + xi2) ** 2 - A1A2
         !gInt = gInt / A1A2 ** 2
      case(3,5,7,9,11,13,15,17,19)  ! q^iOp, with iOp odd
         gInt = cZero
         cdummy1 = xi1c + xi2
         A1A2 = - 2.d0 * (A1 + A2)
         do k = 0, iOp / 2
            gInt = gInt + Coeff(k,iOp) &
                                  / A1A2 ** (iOp - k) &
                                  * cdummy1 ** (iOp - k - k)
         end do
      case(4,6,8,10,12,14,16,18,20)   ! q^iOp, with iOp even: (xi1c + xi2) ** 0 is set to 1 even for xi1c + xi2 = 0
         gInt = cZero
         A1A2 = - 2 * (A1 + A2)
         cdummy2 = (xi1c + xi2) ** 2 / A1A2
         cdummy1 = cOne
         do k = iOp / 2, 0, -1
            gInt = gInt + Coeff(k,iOp) * cdummy1
            cdummy1 = cdummy1 * cdummy2
         end do
         gInt = gInt / A1A2 ** (iOp / 2)
      case(201:)  ! Parameter-dependent operators
         select case(OpParamFunc(iOp - 200))
            case(0)    ! Gaussian exp(p1*x + p2*x**2)
               p1 = OpParam(1,iOp - 200)
               p2 = OpParam(2,iOp - 200)
               dummy = (A1 + A2) / (A1 + A2 + p2)
               cdummy1 = (xi1c + xi2) ** 2 / (A1 + A2) - (xi1c + xi2 + p1) ** 2 / (A1 + A2 + p2)
               gInt = sqrt(dummy) * exp(dOneFourth * cdummy1)
            case(3)    ! Hermite function x**r * exp(p1*x + p2*x**2)
               p1 = OpParam(1,iOp - 200)
               p2 = OpParam(2,iOp - 200)
               idummy = int(OpParam(3,iOp - 200))   ! power r
               cdummy1 = GauInt1D(A1,A2 + p2,xi1,xi2 + p1,idummy)
               cdummy1 = cdummy1 * exp(- dOneFourth * (xi1c + p1 + xi2)**2 / (A1 + A2 + p2) &
                                   + dOneFourth * (xi1c + xi2)**2 / (A1 + A2))
               gInt = cdummy1 * sqrt((A1 + A2) / (A1 + A2 + p2))
         end select
   end select
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end function GauInt1D



   pure double complex function Gau10Int1D(A1,A2,xi1,xi2,iOp)
! It yields the integral <dg_1/dxi_1 |Operator|g_2> / <g_1|g_2>, where g_i is a Gaussian function:
! g_i = exp(A_i * x^2 + xi_i * x)
!    The operator is defined by the flag iOp
   integer, intent(in) :: iOp
   double precision, intent(in) :: A1, A2
   double complex, intent(in) :: xi1, xi2
   integer :: idummy
   double precision :: p1, p2
   double complex :: xi1c
   double precision, parameter :: dOneFourth = 0.25d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   ! CONTROLLARE BENE LE FORMULE
   select case(iOp)
      case(-1) 
         Gau10Int1D = - ciOne * xi2 * GauInt1D(A1,A2,xi1,xi2,1) &
                      - complex(dZero,2.d0 * A2) * GauInt1D(A1,a2,xi1,xi2,2)
      case(-2)
         Gau10Int1D = (A2 + A2 + xi2 ** 2) * GauInt1D(A1,A2,xi1,xi2,1) &
                    + 4.d0 * A2 * xi2 * GauInt1D(A1,A2,xi1,xi2,2) &
                    + 4.d0 * A2 ** 2 * GauInt1D(A1,A2,xi1,xi2,3) 
      case(-3)
         Gau10Int1D = - ciOne &
                    * (GauInt1D(A1,A2,xi1,xi2,1) &
                    + 2.d0 * xi2 * GauInt1D(A1,A2,xi1,xi2,2) &
                    + 4.d0 * A2 * GauInt1D(A1,A2,xi1,xi2,3))
      case(0:19)
         Gau10Int1D = GauInt1D(A1,A2,xi1,xi2,iOp + 1) 
      case(201:)    ! Parameter-dependent operator
         select case(OpParamFunc(iOp - 200))
            case(0)    ! Gaussian exp(p1*x + p2*x**2)
               p1 = OpParam(1,iOp - 200)
               p2 = OpParam(2,iOp - 200)
               Gau10Int1D = - dOneHalf * (conjg(xi1) + xi2 + p1) / (A1 + A2 + p2) &
                          * GauInt1D(A1,A2,xi1,xi2,iOp)
            case(1)    ! Interpolating function
               stop
            case(3)    ! Hermite function x**r * exp(p1*x + p2*x**2)
               xi1c = conjg(xi1)
               p1 = OpParam(1,iOp - 200)
               p2 = OpParam(2,iOp - 200)
               idummy = int(OpParam(3,iOp - 200)) + 1
               Gau10Int1D = GauInt1D(A1,A2 + p2,xi1,xi2 + p1,idummy)
               Gau10Int1D = Gau10Int1D * exp(- dOneFourth * (xi1c + p1 + xi2)**2 / (A1 + A2 + p2) &
                                             + dOneFourth * (xi1c + xi2)**2 / (A1 + A2))
               Gau10Int1D = Gau10Int1D * sqrt((A1 + A2) / (A1 + A2 + p2))
         end select
   end select 
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end function Gau10Int1D


   double complex function GauIntND(n,A1,A2,xi1,xi2,p1,p2)
! It calculates the integral < g1 | exp(p1*x + p2*x**2 | g2 >
! where g1 and g2 are normalized multi-dimensional Gaussians
   integer, intent(in) :: n
   double precision, dimension(n), intent(in) :: A1, A2, p1, p2
   double complex, dimension(n), intent(in) :: xi1, xi2
   double precision :: dummy
   double complex :: cdummy
!   double precision, parameter :: dFour = 4.d0, dOneFourth = 0.25d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Prefactor
   dummy = sqrt(product(4 * A1 * A2 / (A1 + A2 + p2) ** 2))
   cdummy = sum(real(xi1) ** 2 / (16 * A1) + real(xi2) ** 2 / (16 * A2) &
          - (conjg(xi1) + xi2 + p1) ** 2 / (4 * (A1 + A2 + p2)))
   GauIntND = sqrt(dummy) * exp(cdummy)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end function GauIntND


end module Hamilton
