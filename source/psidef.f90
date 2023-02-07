module psidef
   use sysparam
   implicit none
   private
   public :: npackets, nstates, nmodes, ndof, nmodeDOF, modeOfDOF, nSPF, iDOF, ndofDVR, ndofGWP, &
             nAconf, nAconftot, WPpop, Avector, APacketOffset, AStateOffset, idofDVR, idofGWP, indDVR, &
             dimPsi, psi, DD, ipsiSPF, input_file_psi, GWPa, AVector_Start, psi_Start, &
             nmodDVR, nmodGWP, imodDVR, imodGWP, nSPFmax, Aux1, Aux2, MultiSet, &
             iCQ, iCP, iCQQ, iCQP, iCPQ, iCPP
   ! 
   integer :: npackets   ! number of wave packets expanded in the same set of configurations
   integer :: nstates    ! number of electronic states
   integer :: nmodes     ! number of logical modes 
   integer :: nmodDVR    ! number of logical modes for which a grid representation is used
   integer :: nmodGWP    ! number of logical modes for which a Gaussian representation is used
   integer :: ndof       ! number of degrees of freedom
   integer :: ndofDVR    ! number of dofs for which a grid representation is used
   integer :: ndofGWP    ! number of dofs for which a Gaussian representation is used
   integer :: nAconftot  ! total number of configurations
   integer :: dimpsi     ! number of SPF parameters to be propagated
   integer :: nSPFmax    ! maximum number of SPFs for any mode and any electronic state
   !
   character(len = 30) :: input_file_psi   ! alias for the input file
   ! 
   integer, dimension(:),  allocatable :: nmodeDOF, &     ! nmodeDOF(i) = no. of dof for the mode i 
                                          modeOfDOF, &    ! modeOfDOF(i) = mode to which the DOF i belongs
                                          nAconf, &       ! nAconf(i) = no. of configurations for the el. state i
                                          APacketOffset,& ! APacketOffset(i) = last index in the A-Vector first config. for the packet i
                                          AStateOffset, & ! AStateOffset(i)  = last index in the A-Vector before the first config. on the state i
                                          imodDVR, &      ! imodDVR(i) = i-th mode for which a grid representation is used
                                          imodGWP, &      ! imodGWP(i) = i-th mode for which a Gaussian representation is used
                                          idofDVR, &      ! idofDVR(i) = i-th dof for which a grid representation is used
                                          idofGWP, &      ! idofGWP(i) = i-th dof for which a Gaussian representation is used
                                          indDVR          ! indDVR(i) = j, such that idofDVR(j) = i
   integer, dimension(:,:), allocatable :: nSPF, &   ! nSPF(i,j) = no. of SPFs for the mode i and the el. state j
                                           iDOF      ! iDOF(i,j) = i-th dof of mode j
   integer, dimension(:,:), allocatable :: ipsiSPF   ! ipsiSPF(j,i) = index of the psi vector corresponding to the beginning
                                                     ! of the block of SPFs of the mode j of the state i
   ! Multi-set or single-set character for multi-state dynamics
   logical, dimension(:), allocatable :: MultiSet    ! MultiSet(i) = .false. if the single-set formalism is used for mode i
   !
   double precision, dimension(:,:,:), allocatable :: GWPa  ! GWPa(i,j,k) = GWP A-parameter for the i-th dof of the GWP-mode k for the SPF j
   !
   double complex, dimension(:), allocatable :: AVector,       & ! vector of the A-Coefficients
                                                psi,           & ! vector of the SPFs
                                                AVector_Start, & ! initial AVector
                                                psi_Start        ! initial SPFs vector
   double precision, dimension(:), allocatable :: WPpop     ! WPpop(i) = population (norm) of the wave packet i
!!!
! For the perturbative harmonic bath
   integer :: iCQ, iCP, iCQQ, iCQP, iCPQ, iCPP      ! Initial indices for the cross-correlation matrices
                                                    ! of the bath operators Q, P, QQ, QP, PQ, and PP
!!!
! For the constant mean field scheme
   double complex, dimension(:,:,:,:), allocatable :: DD
   ! Auxiliary arrays
   double complex, dimension(:), allocatable :: Aux1, Aux2
   ! Targets for pointers
   target :: psi, nSPF, Aux1, Aux2

end module psidef
