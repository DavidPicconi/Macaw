subroutine propagate_wf
use sysparam
use globals
use integration
use cmf2_base
use storage, only : energyMF
use dvr, only : ngp
use psidef
use derivatives
use hamilton
use output
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
implicit none
integer :: a_status
integer :: i, im, iProp
integer :: ipsi, npsi
double precision :: tnext, step_done, step_next
double precision :: wfnorm
double complex, dimension(:), allocatable :: dAVector, dpsi
logical :: wr_psi_in = .false.
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Open additional output files
   ! Autocorrelation function file
iauto_unit = freeunit()
open(iauto_unit, file = 'auto', status = 'replace', form = 'formatted', action = 'write')
write(iauto_unit,'(a)') '#    time[fs]         Re(autocorrel)     Im(autocorrel)     Abs(autocorrel)'
   ! If multiple electronic states are present, open the spop file
if (nstates .gt. 1) then
   ispop_unit = freeunit()
   open(ispop_unit, file = 'spop', status = 'replace', form = 'formatted', action = 'write')
   write(ispop_unit,'(a)') '#  time[fs]      State populations'
end if
   ! Steps file
if (wr_steps) then
   isteps_unit = freeunit()
   open(isteps_unit, file = 'steps', status = 'replace', form = 'formatted', action = 'write')
   write(isteps_unit,'(a)') '#                Integration steps  [fs]'
   write(isteps_unit,'(a)') '#  time [fs]      attempted       done       suggested'
endif
   ! Wavefunction file
if (wr_psi) then
   wr_psi_in = .true.
   ipsi_unit = freeunit()
   open(ipsi_unit, file = 'psi', status = 'replace', form = 'unformatted', action = 'write')
end if
   ! Density matrix file
if (wr_densmat) then
   idm_unit = freeunit()
   open(idm_unit, file = 'densmat', status = 'replace', form = 'unformatted', action = 'write')
end if
   ! 1D density file
if (wr_density) then
   idens_unit = freeunit()
   open(idens_unit, file = 'density', status = 'replace', form = 'unformatted', action = 'write')
   allocate(Density(maxval(ngp),nstates), DensityPop(maxval(ngp),nstates))
   allocate(DensityMatrix(maxval(ngp),maxval(ngp)))
end if
   ! Wave packet overlap file
if (wr_wpOvl) then
   iwpovl_unit = freeunit()
   open(iwpovl_unit, file = 'wpOvl', status = 'replace', form = 'formatted', action = 'write')
   write(iwpovl_unit,'(a)') 'Wave packet overlap matrix'
   write(iwpovl_unit,*)
endif
   ! Expectation values file
if (nExpect .gt. 0) then
   iexpect_unit = freeunit()
   open(iexpect_unit, file = 'expectation', status = 'replace', form = 'formatted', action = 'write')
end if
   !
write(iout_unit,*)
write(iout_unit,'(a,f9.3,a,f9.3,a)') 'Propagation from ', &
         mintime * au2fs, ' fs to ', maxtime * au2fs, ' fs'
npsi = -1
if (wr_psi) then
   wr_psi_in = .true.
   npsi = int(tpsi / writetime)
   write(iout_unit,'(a,f7.3,a)') 'The wavefunction is stored every ', npsi * writetime * au2fs, ' fs'
end if
if (wr_densmat) then
   write(iout_unit,'(a,20(i0,1x))') &
      'The file densmat contains the density matrices and SPFs for the modes: ', &
      imode_densmat(1:nmodes_densmat)
end if
write(iout_unit,*)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Calculate initial energy
do i = 1, nmodGWP
   call Calc_S00(psi,i)
end do
do im = 1, nmodes
   call Calc_rho1mode(im,AVector)
   call Calc_h1mode(psi,im)
end do
call Calc_MeanField(AVector,psi,mintime)
energyStart = energyMF
! Allocate dAvector and dpsi
allocate(dAVector(nAconftot * npackets), dpsi(dimpsi), stat = a_status)
if (a_status .ne. 0) call err_alloc('dAVector/dpsi', 'propagate_wf',a_status)
!
maxtime = maxtime + writetime * 0.1d0
time = mintime
stepsize = writetime * eps_int    ! initial step size
!stepsize = writetime * 0.1d0
! Set counter for writing the wavefunction
ipsi = 0
!
select case(int_type)  ! Different procedures for different integrators
   case(1:2) ! Runge-Kutta
      do while(time .lt. maxtime)
         ! Write output files (populations, autocorrelation function, etc.)
         call WriteOutput(time)
         !
         tnext = time + writetime
         do
            if (int_type == 1) call Calc_Derivs(AVector,psi,dAVector,dpsi,time)   ! only for Cash-Karp parameters
            call rkqs(AVector,psi,dAVector,dpsi,time,stepsize,step_done,step_next,int_type)
            if (wr_steps) then
               write(isteps_unit,'(f12.6,3(2x,f12.6))') time * au2fs, stepsize * au2fs, step_done * au2fs, step_next * au2fs
               flush(isteps_unit)
            endif
            time = time + step_done
            if (time .ge. tnext - min_stepsize) exit
            stepsize = merge(tnext - time, &
                            (tnext - time) / (int((tnext - time) / step_next) + 1), &
                            int((tnext - time) / step_next) == 0)
            !stepsize = (tnext - time) / (int((tnext - time) / step_next) + 1)
            !stepsize = min(stepsize, tnext - time)
            if (inorm_psi == 3) then
               call trace_wf(AVector,psi,wfnorm)
               !call norm_wf(AVector,psi,wfnorm)
               AVector = AVector / sqrt(wfnorm)
            endif
         end do
         ipsi = ipsi + 1
         if (ipsi .eq. npsi .and. wr_psi_in) then
            ipsi = 0
            wr_psi = .true.
         end if
         ! Check if the wavefunction needs to be normalized
         if (inorm_psi .ne. 0) then
            call trace_wf(AVector,psi,wfnorm)
            AVector = AVector / sqrt(wfnorm)
            inorm_psi = merge(0, inorm_psi, inorm_psi == 1)
         endif
      end do
   case(3) ! Adams PECE formula
      ! Allocate variables needed in Adams method
      call SetupAdams
      ! Calculate derivatives at the initial time
      call Calc_Derivs(AVector,psi,dAVector,dpsi,time)
      !
      iProp = 0
      do while(time .lt. maxtime)
         ! Write output files (populations, autocorrelation function, etc.)
         call WriteOutput(time)
         !
         tnext = time + writetime
         do
            if (iProp .lt. int_order) then ! the first int_order steps are done with Runge-Kutta (Cash-Karp)
               dA_prev(:,int_order - iProp) = dAVector
               dpsi_prev(:,int_order - iProp) = dpsi
               if (inorm_psi == 3) then
                  call trace_wf(AVector,psi,wfnorm)
                  AVector = AVector / sqrt(wfnorm)
               endif
               call rkqs(AVector,psi,dAVector,dpsi,time,stepsize,step_done,step_next,1)  ! Avector and psi are updated here
               call Calc_Derivs(AVector,psi,dAVector,dpsi,time)
               if (wr_steps) then
                  write(isteps_unit,'(f12.6,3(2x,f12.6))') time * au2fs, stepsize * au2fs, step_done * au2fs, step_next * au2fs
                  flush(isteps_unit)
               endif
               ! Store the step sizes and the derivatives
               h_prev(int_order - iProp) = step_done
               iProp = iProp + 1
               ! Update time
               time = time + step_done
               if (time .ge. tnext - min_stepsize) exit
               stepsize = merge(tnext - time, &
                               (tnext - time) / (int((tnext - time) / step_next) + 1), &
                                int((tnext - time) / step_next) == 0)
               !if (time .ge. tnext) exit
               !stepsize = (tnext - time) / (int((tnext - time) / step_next) + 1)
               !stepsize = min(stepsize, tnext - time)
               cycle
            end if
            ! Adams' integration
            call Adams(AVector,psi,dAVector,dpsi,time,stepsize,step_done,step_next)
            if (wr_steps) then
               write(isteps_unit,'(f12.6,3(2x,f12.6))') time * au2fs, stepsize * au2fs, step_done * au2fs, step_next * au2fs
               flush(isteps_unit)
            endif
            time = time + step_done
            if (time .ge. tnext) exit
            stepsize = (tnext - time) / (int((tnext - time) / step_next) + 1)
            stepsize = min(stepsize, tnext - time)
         end do
         ipsi = ipsi + 1
         if (ipsi .eq. npsi .and. wr_psi_in) then
            ipsi = 0
            wr_psi = .true.
         end if
         ! Check if the wavefunction needs to be normalized
         if (inorm_psi .ne. 0) then
            call trace_wf(AVector,psi,wfnorm)
            !call norm_wf(AVector,psi,wfnorm)
            AVector = AVector / sqrt(wfnorm)
            inorm_psi = merge(0, inorm_psi, inorm_psi == 1)
         endif
      end do
   case (4)   ! Bulirsch-Stoer
      ! Allocate variables needed in BS method
      call SetupBS
      step_next = -1.d15   ! this needs to be initialized to some "impossible" value
      do while(time .lt. maxtime)
         ! Write output files (populations, autocorrelation function, etc.)
         call WriteOutput(time)
         !
         tnext = time + writetime
         do
            call Calc_Derivs(AVector,psi,dAVector,dpsi,time)
            call bsstep(AVector,psi,dAVector,dpsi,time,stepsize,step_done,step_next)  ! Avector and psi are updated here
            if (wr_steps) then
               write(isteps_unit,'(f12.6,3(2x,f12.6))') time * au2fs, stepsize * au2fs, step_done * au2fs, step_next * au2fs
               flush(isteps_unit)
            endif
            time = time + step_done
            if (time .ge. tnext - min_stepsize) exit
               stepsize = merge(tnext - time, &
                               (tnext - time) / (int((tnext - time) / step_next) + 1), &
                                int((tnext - time) / step_next) == 0)
            !if (time .ge. tnext) then                 ! next output point has been reached
            !   stepsize = min(step_next, writetime)   ! don't ask for a step size larger than the output time
            !   exit
            !else
            !   stepsize = (tnext - time) / (int((tnext - time) / step_next) + 1)   ! try to make the step sizes uniform
            !   stepsize = min(stepsize, tnext - time)
            !endif
            if (inorm_psi == 3) then
               call trace_wf(AVector,psi,wfnorm)
               !call norm_wf(AVector,psi,wfnorm)
               AVector = AVector / sqrt(wfnorm)
            endif
         end do
         ipsi = ipsi + 1
         if (ipsi .eq. npsi .and. wr_psi_in) then
            ipsi = 0
            wr_psi = .true.
         end if
         ! Check if the wavefunction needs to be normalized
         if (inorm_psi .ne. 0) then
            call trace_wf(AVector,psi,wfnorm)
            AVector = AVector / sqrt(wfnorm)
            inorm_psi = merge(0, inorm_psi, inorm_psi == 1)
         endif
      end do
   case (10)   ! Constant mean field scheme
      if (npackets > 1) then
         write(0,*) ' ERROR:'
         write(0,*) '    Constant-mean field integrator not available for multi-packet run'
         write(0,*) '    (it might become available for non-interacting wave packets)'
         stop
      end if
      call SetupCMF2
      do while (time .lt. maxtime)
         ! Write output files (populations, autocorrelation function, etc.)
         call WriteOutput(time)
         !
         tnext = time + writetime
         do
            call cmf2(AVector,psi,DD,time,stepsize,step_done,step_next)  ! Avector and psi are updated here
            if (wr_steps) then
               write(isteps_unit,'(f12.6,3(2x,f12.6))') time * au2fs, stepsize * au2fs, step_done * au2fs, step_next * au2fs
               flush(isteps_unit)
            endif
            time = time + step_done
            if (time .ge. tnext - min_stepsize) exit
               stepsize = merge(tnext - time, &
                               (tnext - time) / (int((tnext - time) / step_next) + 1), &
                                int((tnext - time) / step_next) == 0)
            !if (time .ge. tnext) then                 ! next output point has been reached
            !   stepsize = min(step_next, writetime)   ! don't ask for a step size larger than the output time
            !   exit
            !else
            !   stepsize = (tnext - time) / (int((tnext - time) / step_next) + 1)   ! try to make the step sizes uniform
            !   stepsize = min(stepsize, tnext - time)
            !endif
            if (inorm_psi == 3) then
               call trace_wf(AVector,psi,wfnorm)
               !call norm_wf(AVector,psi,wfnorm)
               AVector = AVector / sqrt(wfnorm)
            endif
         end do
         ipsi = ipsi + 1
         if (ipsi .eq. npsi .and. wr_psi_in) then
            ipsi = 0
            wr_psi = .true.
         end if
         ! Check if the wavefunction needs to be normalized
         if (inorm_psi .ne. 0) then
            call trace_wf(AVector,psi,wfnorm)
            !call norm_wf(AVector,psi,wfnorm)
            AVector = AVector / sqrt(wfnorm)
            inorm_psi = merge(0, inorm_psi, inorm_psi == 1)
         endif
      end do
end select
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Close output files
close(iauto_unit)
if (nstates .gt. 1) close(ispop_unit)
if (wr_psi_in) close(ipsi_unit)
if (wr_densmat) close(idm_unit)
if (wr_steps) close(isteps_unit)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine propagate_wf
