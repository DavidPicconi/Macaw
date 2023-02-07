module storage
   use sysparam
   use psidef
   implicit none
   private
   public :: S00, S00m1, Sa0Sm1, Sa0D, rho, rho_ave, MeanField, rhom1MF, energyMF, &
             Cmat, PrecC, nCYmax, Yvec, nS00max, Adot_tmp, multLag, multLagDiag, FQBath,FPBath
   public :: alloc_S_matrices, alloc_rho, alloc_CY
   !
   integer :: nS00max  ! size of dimensions 1 and 2 of array S00
   integer :: nCYmax   ! size of dimensions 1 (1 and 2) of the arrays Yvec (Cmat)
   double precision :: energyMF    ! energy calculated during the mean field evaluation
   double complex, dimension(:,:,:,:,:), allocatable :: S00        ! S00(i,j,k,l,l) = overlap between the SPFs i,j of the GWP-mode k for the electronic states l,m
   double complex, dimension(:,:,:,:), allocatable :: S00m1      ! S00(i,j,k,l) = element of S^-1 matrix (see definition of S00, with only one el. state)
   double complex, dimension(:,:,:,:), allocatable :: Sa0Sm1     ! Sa0Sm1(i,j,k,l) = element of matrix Sa0 * S^-1 (see definition of S00, with only one el. state)
   double complex, dimension(:,:,:,:), allocatable :: Sa0D     ! Sa0Sm1(i,j,k,l) = element of matrix Sa0 * D (see definition of S00, with only one el. state)
                                                               !                 The index j refers to both the SPF and the DOF
   double complex, dimension(:,:,:,:), allocatable :: Cmat       ! Cmat(i,j,k,l) = element of C matrix for the GWP-mode k and the el. state l
                                                               !                 i,j are indices for the SPF and the dof
   double complex, dimension(:,:,:,:), allocatable :: PrecC      ! Preconditioner for CG iterations on the C matrix
   double complex, dimension(:,:,:), allocatable :: Yvec         ! Yvec(i,k,l) = element of Y vector for the GWP-mode k and the el. state l
                                                               !                 i is an index which includes SPFs and the dof
   double complex, dimension(:,:,:,:), allocatable :: rho          ! rho(i,j,k,l) = (i,j) element of the density matrix of the mode k for the electronic state l
   double complex, dimension(:,:,:), allocatable :: rho_ave        ! Trace of rho over the electronic states
   double complex, dimension(:,:,:,:), allocatable :: MeanField    ! MeanField(i,j,k,l) = (i,j) element of the mean field of the mode l for the Hamiltonian term k
   double complex, dimension(:,:,:,:), allocatable :: rhom1MF      ! rhom1MF(i,j,k,l) = (i,j) element of the matrix (rho^-1 * MeanField) 
                                                                 ! for the DVR-mode k for the Hamiltonian term l
   double complex, dimension(:), allocatable :: Adot_tmp         ! Auxiliary array used to calculate the derivative of coefficients. Calculated in Calc_MeanField
   double complex, dimension(:,:), allocatable :: multLag          ! Lagrange multiplier to enforce the orthogonality between wave packets
   double precision :: multLagDiag                                 ! Lagrange multiplier to impose the conservation of the total population
!!!   
   !
   target :: Adot_tmp, Sa0Sm1, Sa0D
!!!
! For the Harmonic operator bath
   double complex, dimension(:,:,:), allocatable :: FQbath       ! FQbath(i,j,k) = transition matrix element
                                                                 !                 of the system-bath operator F_k  ( * Q_k)
                                                                 !                 between the wave packets i and j
   double complex, dimension(:,:,:), allocatable :: FPbath       ! FPbath(i,j,k) = transition matrix element
                                                                 !                 of the system-bath operator F_k  ( * P_k)
                                                                 !                 between the wave packets i and j
!
contains



   subroutine alloc_S_matrices
! It allocates the matrices S00 and S00m1
! The matrix S0am is allocated into alloc_CY
   implicit none
   integer :: a_status, iel, i, im, nmD, nS, j
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! S00
   allocate(S00(nS00max,nS00max,nmodGWP,nstates,nstates), stat = a_status)
   if (a_status .ne. 0) call err_alloc('S00', 'alloc_S_matrices', a_status)
      ! Initialise diagonal elements as 1 (they are kept equal to 1 for the whole simulation)
   do iel = 1, nstates
      do i = 1, nmodGWP
         im = imodGWP(i)
         nmD = nmodeDOF(im)
         nS = nSPF(im,iel)
         forall (j = 1:nS) S00(j,j,i,iel,iel) = cOne
      enddo
   enddo
   ! S00m1
   allocate(S00m1(nS00max,nS00max,nmodGWP,nstates), stat = a_status)
   if (a_status .ne. 0) call err_alloc('S00m1', 'alloc_S_matrices', a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!! 
   return
   end subroutine alloc_S_matrices


   
   subroutine alloc_rho
   implicit none
   integer :: a_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   allocate(rho(nSPFmax,nSPFmax,nmodes,nstates), stat = a_status)
   allocate(rho_ave(nSPFmax,nSPFmax,nmodes), stat = a_status)
   if (a_status .ne. 0) call err_alloc('rho', 'alloc_rho', a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine alloc_rho



   subroutine alloc_CY
   implicit none
   integer :: a_status
   integer :: i, im
!!!i!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   nCYmax = 0
   do i = 1, nmodGWP
      im = imodGWP(i)
      nCYmax = max(nCYmax, maxval(nSPF(im,:)) * nmodeDOF(im))
   enddo
   allocate(Yvec(nCYmax,nstates,nmodGWP), &
            Cmat(nCYmax,nCYmax,nstates,nmodGWP), &
            PrecC(nCYmax,nCYmax,nstates,nmodGWP), stat = a_status)
   if (a_status .ne. 0) call err_alloc('C, Y', 'alloc_CY', a_status)
   Cmat = cZero
   Yvec = cZero
! Allocate S0am
   allocate(Sa0Sm1(nCYmax,nS00max,nmodGWP,nstates), stat = a_status)
   if (a_status .ne. 0) call err_alloc('S0am', 'alloc_CY', a_status)
   allocate(Sa0D(nCYmax,nS00max,nmodGWP,nstates), stat = a_status)
   if (a_status .ne. 0) call err_alloc('Sa0D', 'alloc_CY', a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine alloc_CY



end module storage
