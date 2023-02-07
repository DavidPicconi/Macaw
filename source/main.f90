program gMCTDH
use sysparam
use globals
use timingmod
use input
implicit none
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Initialise timing
call init_timing()
! Open the log file
ilog_unit = freeunit()
open(unit = ilog_unit, file = 'log', form = 'formatted', status = 'replace', action = 'write')
! Read the input file
call read_input()
!!!!
!
! Initialise the calculation
!
! Open the main output file
iout_unit = freeunit()
open(iout_unit, file = 'output', status = 'replace', form = 'formatted', action = 'write')
! Initialise wavefunction and variables needed in the propagation 
call init_calc
! Start the propagation
call propagate_wf
! Write final message and the timing file
write(iout_unit,*)
write(iout_unit,'(a)') 'Propagation was succesful'
call write_timings()
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Close the output files
close(ilog_unit)
close(iout_unit)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
end program gMCTDH
