module sysparam
! It contains:
! - system parameters such as KINDS and predefined variables for 1 and 0.
! - a function to get the next free io-unit
! - subroutine to report allocation and i/o errors
   implicit none
   private
   public :: dOne, dZero, dOneHalf, cOne, ciOne, cZero, cOneHalf
   public :: freeunit, err_alloc, err_read
   !Predefined 0 and 1 as double and complex
   double precision, parameter :: dOne = 1.d0, dZero = 0.d0, dOneHalf = 0.5d0
   double complex, parameter :: cOne = (1.d0,0.d0), ciOne = (0.d0,1.d0), &
                               cZero = (0.d0,0.d0), cOneHalf = (0.5d0,0.d0)
!
contains


   integer function freeunit()
   ! Searches for the next available I/O unit
   logical :: used
   integer, parameter :: iinp = 5, iout = 6, &
                         if90inp = 100, if90out = 101, if90err = 102, &
                         imax = 299
   integer, dimension(5), parameter :: busy = (/ iinp, iout, if90inp, if90out, if90err /)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   freeunit = 0
   used = .true.
   do while (used)
     freeunit = freeunit + 1
     if (freeunit .gt. imax) stop 'No further I/O units are available'
     if (any(busy .eq. freeunit)) cycle
     inquire(unit = freeunit, opened = used)    
   enddo
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end function freeunit

 
   subroutine err_alloc(array, sub, a_status)
   implicit none
   character(len = *) :: array, sub
   integer :: a_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   write(0,'(a,1x,i0)') '--- ALLOCATION ERROR', a_status
   write(0,'(7x,2a,5x,2a)') 'array: ', array, 'subroutine: ', sub
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   stop
   end subroutine err_alloc

 
   subroutine err_read(rfile, sub, f_status)
   implicit none
   character(len = *) :: rfile, sub
   integer :: f_status
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   write(0,'(a,1x,i0)') '--- READ ERROR', f_status
   write(0,'(7x,2a,5x,2a)') 'file: ', trim(adjustl(rfile)), 'subroutine: ', sub
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   stop
   end subroutine err_read

end module sysparam

