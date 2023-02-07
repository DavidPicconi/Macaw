module timingmod
   use sysparam
   implicit none
   private
   public :: init_timing, suspend_timer, continue_timer, write_timings
   !
   integer, save :: CountStart, CountRate
   integer, parameter :: n_timers = 8    ! number of subroutines for which the timing is done
   character(len = 20), dimension(1:n_timers) :: routine_names
   logical, dimension(n_timers), save :: timer_running = .false.
   integer, dimension(n_timers) :: routine_calls
   integer, dimension(n_timers), save :: starttime, routine_net_total
contains


   subroutine init_timing()
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Get initial time and clock count rate
   call system_clock(count = CountStart, count_rate = CountRate)
   routine_names(1) = 'CalcDerivs'
   routine_names(2) = 'normWF/wf_ovl'
   routine_names(3) = 'Calc_MeanField'
   routine_names(4) = 'Calc_CY'
   routine_names(5) = 'Calc_h1mode'
   routine_names(6) = 'Calc_Der_SPF_DVR'
   routine_names(7) = 'Calc_Der_SPF_GWP'
   routine_names(8) = 'Lanczos'
   starttime = 0
   routine_net_total = 0
   routine_calls = 0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine init_timing



   subroutine suspend_timer(ri)
   integer, intent(in) :: ri
   integer :: endtime
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (.not. timer_running(ri)) then
      write(0,*) 'Timer', ri
      stop '   you tried to suspend a non-running timer'
   else
      timer_running(ri) = .false.
      call system_clock(count = endtime)
      routine_net_total(ri) = routine_net_total(ri) + endtime - starttime(ri)
      routine_calls(ri) = routine_calls(ri) + 1
   endif
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine suspend_timer

 

   subroutine continue_timer(ri)
   integer, intent(in)::ri
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (.not. timer_running(ri)) then
      timer_running(ri) = .true.
      call system_clock(count = starttime(ri))
   else
      write(0,*) 'Timer', ri
      stop '   You tried to start a running timer'
   endif
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine continue_timer



   subroutine write_timings()
   integer :: endtime, walltime_h, walltime_m, walltime_s
   double precision :: walltime, subtime
   integer :: itime_unit, ri, ncalls
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   itime_unit = freeunit()
   open(unit = itime_unit, file = 'timing')
   call system_clock(count = endtime)
   ! Write the total walltime
   walltime = (endtime - CountStart) / dble(CountRate)
   walltime_h = int(walltime / 3600.0d0)
   walltime_m = int((walltime - (walltime_h * 3600.0d0)) / 60.0d0)
   walltime_s = int(walltime - walltime_h * 3600.0d0 - walltime_m * 60.0d0)
   write(itime_unit,'(a,i0,a,i0,a,i0,a)') &
     'Total time:   ', walltime_h, ' h   ', walltime_m, ' m   ', walltime_s, ' s'
   write(itime_unit,*)
   write(itime_unit,*)
   ! Write time information for individual subroutines
   write(itime_unit,*) '         Subroutine    Calls N    time/N [ms]', &
                       '      time [min]    %time'
   write(itime_unit,*) 
   do ri = 1, n_timers
      ncalls = routine_calls(ri)
      subtime = routine_net_total(ri) / dble(CountRate)
      write(itime_unit,'(a20,4x,i7,3x,f9.2,3x,f11.2,6x,f6.2)') &
         trim(adjustl(routine_names(ri))), ncalls, subtime / dble(ncalls) * 1000.d0, &
         subtime / 60.d0, subtime / walltime * 100.d0
   enddo   
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   close(itime_unit)
   return
   end subroutine write_timings


end module timingmod
