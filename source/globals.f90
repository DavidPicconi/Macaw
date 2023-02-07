module globals
   use sysparam
   implicit none
   private
   public :: pi, au2eV, au2fs, au2cmm1, & ! Physical constants
             ilog_unit, iout_unit, iauto_unit, ispop_unit, isteps_unit, &
             ipsi_unit, idm_unit, idens_unit, iexpect_unit, iwpovl_unit, &
             read_integer, read_float, & !I/O-Routines
             time, stepsize, min_stepsize, mintime, maxtime, &
             int_type, int_order, int_type_A, int_type_psi, int_par1, writetime, tpsi, & !propagation parameters
             eps_rho, eps_ovl, eps_c, eps_int, eps_A, eps_psi, eps_pop, &  ! regularization parameters
             wr_psi, relaxation, inorm_psi, der_factor, &
             wr_densmat, nmodes_densmat, imode_densmat, wr_density, wr_wpOvl, wr_steps, & 
             nStrings, nStringsMax, strings
   ! Physical and mathematical constants
   double precision, parameter :: pi = 4.0d0 * datan(dOne), &
                           au2eV = 27.2113961d0,      au2fs = 0.02418884326505d0, &
                           au2cmm1 = 219474.63d0
   ! Propagation Settings
   double precision :: time, stepsize, min_stepsize, mintime, maxtime, writetime
   double precision :: tpsi
   logical :: relaxation    ! whether performing a relaxation run
   integer :: inorm_psi     ! = 0 : never normalize the wavefunction
                            ! = 1 : normalize the wavefunction only up to the first output time
                            ! = 2 : normalize the wavefunction every output time
                            ! = 3 : normalize the wavefunction every integration step
   double complex :: der_factor  ! factor in front of the r.h.s.:  dPsi/dt = f * H * Psi
                                 ! f = -i for propagation run
                                 ! f = -1 for relaxation run
   ! Integration Properties
   integer :: int_type, int_order, int_par1, int_type_A, int_type_psi
   ! Various parameters
   double precision :: eps_rho    ! 'epsilon' parameter used in the regularisation of the density matrix
   double precision :: eps_ovl    ! 'epsilon' parameter used in the regularisation of the overlap matrix
   double precision :: eps_C      ! 'epsilon' parameter used in the regularisation of the C matrix
   double precision :: eps_int    ! 'epsilon' parameter used in the accuracy of the integrator
   double precision :: eps_A      ! accuracy parameter for the propagation of the coefficients in the CMF scheme
   double precision :: eps_psi    ! accuracy parameter for the propagation of the SPFs in the CMF scheme
   double precision :: eps_pop    ! regularisation parameters for open system dynamics
   ! Output
   integer :: ilog_unit, iout_unit, iauto_unit, ispop_unit, isteps_unit, &
              ipsi_unit, idm_unit , idens_unit, iexpect_unit, iwpovl_unit      ! write units
   logical :: wr_psi, wr_densmat, wr_density, wr_wpOvl, wr_steps  ! write flags 
   integer :: nmodes_densmat
   integer, dimension(:), allocatable :: imode_densmat
   ! Array holding strings
   integer :: nStrings                       ! number of strings: updated every time a string needs to be added
   integer, parameter :: nStringsMax = 2000  ! maximum number of strings allowed (the array strings is preallocated)
   character(len = 100), dimension(nStringsMax) :: strings
       
!
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines used for reading from input
! Maybe they should be moved into input.f90

   integer function read_integer(varname,value)
   character(len = *), intent(in) :: varname, value
   integer :: r_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   read(value, *,iostat = r_status) read_integer
   if (r_status .ne. 0) then
      write(6,'(5a)') 'Error reading "', trim(varname), '". ', trim(value), ' is not an integer!'
      stop 
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end function


   double precision function read_float(varname,value)
   character(len = *), intent(in) :: varname, value
   integer :: r_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   read(value, *, iostat = r_status) read_float
   if (r_status .ne. 0) then
      write(6,'(5a)') 'Error reading "', trim(varname), '". ', trim(value), ' is not a float!'
      stop
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   end function


   character(len = 1) function read_char(varname,value)
   character(len = *), intent(in) :: varname, value
   integer :: r_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   read (value, '(a1)', iostat = r_status) read_char
   if (r_status .ne. 0) then
     write(6,'(5a)') 'Error reading "', trim(varname), '". ', trim(value), ' is not a character!'
     stop 
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   end function

end module globals
