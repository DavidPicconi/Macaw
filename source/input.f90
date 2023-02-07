module input
! This module contains routines to read the input file and the Hamiltonian file
   use sysparam 
   use globals
   use psidef
   use dvr
   use Hamilton
   implicit none
   private
   public :: read_input, nTermsEF, lEF, rEF, coeffEF
   !
   integer :: a_status, f_status
   ! auxiliary matrices
   integer, dimension(:,:), allocatable :: iDOF_t
   ! for the parameter-dependent operators
   integer, dimension(:), allocatable :: OpParamFunc_t
   double precision, dimension(:,:), allocatable :: OpParam_t
   ! Electronic functions (maximum 50 with 50 expansion terms)
   integer, parameter :: maxTermsEF = 50, maxExpEF = 50
   integer, dimension(:) :: nTermsEF(maxTermsEF)            ! nTerms(i) = number of terms of the i-th electronic function
   integer, dimension(:,:) :: lEF(maxExpEF,maxTermsEF), &
                              rEF(maxExpEF,maxTermsEF)      ! lEF(i,j) = left electronic state of the i-th term of the j-th electronic function
   double complex, dimension(:,:) :: coeffEF(maxExpEF,maxTermsEF)  ! coeffEF(i,j) = expansion coefficient
   character(len = 100) :: Hamiltonfile, Dissipationfile
   character(len = 100), dimension(:) :: ExtraOpFile(nExtraOpMax)
!
contains

   
   subroutine read_input()
   ! It processes the input file line by line.
   ! After the input file has been completely processed, the wavefunction arrays are allocated.
   implicit none
   logical :: exists, lDiss
   character(len = 30), save :: inp_fn
   character(len = 60) :: wrfmt, IntName
   character(len = 200) :: line, line_sub, varName, varValue
   integer :: io_inp, i_unit
   integer :: i, im, sindex, nmD, iType
   integer :: nclargs, iargc
   character(len = 30), dimension(:), allocatable :: clargs
   character(len = 30) :: chardummy
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Default values
   npackets = 1
   nstates = 1
   ndof = -1
   nmodes = 0
   Hamiltonfile    = ''
   Dissipationfile = ''
   lDiss   = .false.
   nDiss1  = 0
   nDiss2  = 0
   lHOBath = .false.
   nHOBath = 0
   mintime = dZero
   eps_C   = 1.d-8
   eps_rho = 1.d-8
   eps_ovl = 1.d-8
   eps_int = 1.d-8
   eps_A   = 1.d-8
   eps_psi = 1.d-8
   eps_pop = 1.d-6
   writetime    = dOne / au2fs
   min_stepsize = 1.d-9 / au2fs
   wr_psi     = .false.
   wr_densmat = .false.
   wr_wpOvl   = .false.
   wr_steps   = .false.
   relaxation = .false.
   der_factor = -ciOne   ! The default is a propagation, not a relaxation
   inorm_psi = 0
   nmodes_densmat = 0
   int_type  = 1    ! Default integrator: Runge-Kutta
   int_order = 4    ! It is not referenced if int_type = 1
   !
   nExtraOp = 0
   nExpect  = 0
   iExpect  = 0
   !
   nStrings = 0
   !
! Read the options 
   nclargs = iargc()
   if (nclargs .eq. 0) stop 'Input file must be given!'
   allocate(clargs(nclargs), stat = a_status)
   if (a_status .ne. 0) call err_alloc('clargs', 'read_input', a_status)
   do i = 1, nclargs
      call getarg(i, clargs(i))   ! get the i-th option
   end do 
   inp_fn = clargs(nclargs)       ! the last argument is the input file
   do i = 1, nclargs - 1  ! the last argument is the input file
      select case(trim(adjustl(clargs(i))))
         case default
            write(0,*) '--- INPUT ERROR'
            write(0,'(6x,a,1x,a)') 'Unknown option', trim(adjustl(clargs(i)))
      end select
   end do
   deallocate(clargs, stat = a_status)
   if (a_status .ne. 0) call err_alloc('dealloc. clargs', 'read_input', a_status)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Read the input file
   inquire(file = trim(adjustl(inp_fn)), exist = exists)   ! check that the file exists
   if (.not. exists) then
      write(0,'(3a)') 'Input file ', inp_fn, ' not found!'
      stop
   end if   
   input_file_psi = inp_fn    ! the input file name is stored because it is used also in psimod.f90
   io_inp = freeunit()
   open(unit = io_inp, file = trim(adjustl(inp_fn)), status = 'old', action = 'read', iostat = f_status)
   if (f_status .ne. 0) call err_read(inp_fn, 'read_input', f_status)
   ! Process the input line by line
   do 
      line = ''
      read(io_inp,'(a200)',iostat = f_status) line
      if (f_status .ne. 0) exit
      ! Read variable name and (when applicable) value
      varName = ''
      varValue = ''
      sindex = index(line, '#')
      if (sindex .ne. 0) line = line(1:sindex - 1)
      sindex = 0
      ! Don't read the initial wavefunction now
      if (trim(adjustl(line)) .eq. 'INIT-WF') then
         do 
            read(io_inp,'(a200)',iostat = f_status) line
            sindex = index(line, '#')
            if (sindex .ne. 0) line = line(1:sindex - 1)
            if(trim(adjustl(line)) .eq. 'END-INIT-WF') exit
            if (f_status .ne. 0) then
               write(0,*) 'Input error:'
               write(0,*) '   Closing instruction END-INIT-WF missing'
               stop
            end if
         end do
         read(io_inp,'(a200)',iostat = f_status) line
         if (f_status .ne. 0) exit
         sindex = index(line, '#')
         if (sindex .ne. 0) line = line(1:sindex - 1)
      end if
      !
      sindex = index(line, '=')
      if (sindex .ne. 0) then
         varName = trim(adjustl(line(1:sindex - 1)))
         call lc(varName)
         varValue = trim(adjustl(line(sindex + 1:)))
         call lc(varValue)
         sindex = 0
      else
         varName = trim(adjustl(line))
         call lc(varName)
         varValue = 'blank'
      end if
      ! varName and value have been read. Now read the content
      select case(trim(adjustl(varName)))
         case('npackets')
            npackets = read_integer(varName,varValue)
         case('nstates')
            nstates = read_integer(varName,varValue)
         case('ndof') 
            ndof = read_integer(varName,varValue)
            ! If ndof > 99 there will be problems in reading the Hamiltonian
            if (ndof > 99) stop 'Error: the maximum number of degrees of freedom is 99'
            if (ndof < 0 ) stop 'Error: the number of degrees of freedom should be positive'
         case('nmodes')
            nmodes = read_integer(varName,varValue)
         case('spf')
            if (nmodes .eq. 0) stop 'Input error: with the keyword SPF, nmodes must be > 0'
            if (ndof .eq. 0) stop 'Input error: with the keyword SPF, ndof must be > 0'
            if (nmodes .gt. ndof) stop 'Input error: nmodes > ndof'
            allocate(nSPF(nmodes,nstates), nmodeDOF(nmodes), MultiSet(nmodes), &
                     iDOF_t(ndof,nmodes), modeOfDOF(ndof), stat = a_status)
            if (a_status .ne. 0) call err_alloc('spf arrays', 'read_input', a_status)
            MultiSet = .true.  ! Multi-set is the default
            do im = 1, nmodes
               read(io_inp,'(a200)',iostat = f_status) line_sub
               if (f_status .ne. 0) stop 'Input error: expecting SPF declaration'
               call read_declare_spf(line_sub,im)      
            end do
         case('grid')
            ndofDVR = 0
            ndofGWP = 0
            if (ndof .eq. 0) stop 'Input error: with the keyword GRID, ndof must be > 0'
            allocate(ngp(ndof), dvrtype(ndof), idvrpar(ndof), dvrpar(2,ndof), stat = a_status)
            if (a_status .ne. 0) call err_alloc('dvr arrays', 'read_input', a_status)
            idvrpar = 0
            ngp = 0
            do i = 1, ndof
               read(io_inp,'(a200)',iostat = f_status) line_sub
               if (f_status .ne. 0) stop 'Input error: expecting DVR declaration'
               call read_declare_dvr(line_sub)
            end do 
            ngpMax = maxval(ngp)   
         case('hamiltonian')
            Hamiltonfile = trim(adjustl(varValue))
         case('expect', 'dissipation', 'lcharmbathpt')
            if (.not. allocated(iOpDOFExtra))  then
               if (ndof .eq. -1) then
                  write(0,*) 'Error: ndof should be defined before the keyword expect'
                  stop
               else
                  allocate(iOpDOFExtra(ndof,nExtraOpTermsMax,nExtraOpMax))    
               end if
            end if
            nExtraOp = nExtraOp + 1
            if (nExtraOp .gt. nExtraOpMax) then
               write(0,*) 'Too many extra operators found in the input'
               write(0,*) 'The maximum is ', nExtraOpMax
               stop
            end if
            ExtraOpFile(nExtraOp) = trim(adjustl(varValue))
            select case(trim(adjustl(varName)))
               case('expect')
                  nExpect = nExpect + 1
                  iExpect(nExpect) = nExtraOp
               case('dissipation')
                  lDiss = .true.
                  ! count the number of dissipative operators of type 1 and 2
                  i_unit = freeunit()
                  open(unit = i_unit, file = ExtraOpFile(nExtraOp), status = 'old', action = 'read', form = 'formatted')
                  read(i_unit,*) iType
                  close(i_unit)
                  if (iType == 1) then
                     nDiss1 = nDiss1 + 1
                     iDiss1(nDiss1) = nExtraOp
                  else if (iType == 2) then
                     nDiss2 = nDiss2 + 1
                     iDiss2(nDiss2) = nExtraOp
                  end if
               case('lcharmbathpt')
                  lHObath = .true.
                  nHObath = nHObath + 1
                  iHObathQ(nHObath) = nExtraOp    
                  nExtraOp = nExtraOp + 1
                  ExtraOpFile(nExtraOp) = ExtraOpFile(nExtraOp - 1)   
                  iHObathP(nHObath) = nExtraOp
            end select
         case('harmbathptorder')
            HObathPTorder = read_integer(varName,varValue)
            if (HObathPTorder < 2 .or. HObathPTorder > 4) then
               write(0,*) 'ERROR: at the moment the order in perturbation theory'
               write(0,*) '       for the harmonic bath needs to be 2 or 3'
               stop
            end if
         case('tinit')
            mintime = read_float(varName,varValue)
            mintime = mintime / au2fs
         case('tfinal')
            maxtime = read_float(varName,varValue)
            maxtime = maxtime / au2fs
         case('tout')
            writetime = read_float(varName,varValue)
            writetime = writetime / au2fs
         case('eps_rho')
            eps_rho = read_float(varName,varValue)
         case('eps_c')
            eps_C = read_float(varName,varValue)
         case('eps_int')
            eps_int = read_float(varName,varValue)
         case('eps_s')
            eps_ovl = read_float(varName,varValue)
         case('eps_pop')
            eps_pop = read_float(varName,varValue)
         case('minstep')
            min_stepsize = read_float(varName,varValue)
            min_stepsize = min_stepsize / au2fs
         case('psi')
            wr_psi = .true.
            tpsi = read_float(varName,varValue)
            tpsi = tpsi / au2fs
         case ('relaxation')
            relaxation = .true. 
            inorm_psi = 2
            der_factor = -cOne
         case ('normalize')
            inorm_psi = read_integer(varName,varValue)   
         case('densmat')
            if (ndof .eq. 0) stop 'Input error: with the keyword DENSMAT, ndof must be > 0'
            if (nmodes_densmat .eq. 0) allocate(imode_densmat(nDOF))
            wr_densmat = .true.
            nmodes_densmat = nmodes_densmat + 1
            imode_densmat(nmodes_densmat) = read_integer(varName,varValue)
         case('density')
            wr_density = .true.
         case('wp_ovl')
            wr_wpOvl = .true.
         case('steps')
            wr_steps = .true.
         case('integration')
            IntName = trim(adjustl(varValue))
            if (IntName .eq. 'rk45') then
               int_type = 1
               int_order = 4
            elseif (IntName .eq. 'dopri') then
               int_type = 2
               int_order = 4
            elseif (IntName(1:3) .eq. 'abm') then
               int_type = 3
               read(IntName(4:5),*) int_order
               if (trim(adjustl(IntName(6:))) .eq. '') then
                  int_par1 = 0
               else
                  read(IntName(6:),*) int_par1
               end if
            elseif (IntName(1:2) .eq. 'bs') then
               int_type = 4
            elseif (IntName(1:4) .eq. 'cmf2') then
               int_type = 10
               read(io_inp,*,iostat = f_status) chardummy, eps_A
               if (f_status .ne. 0) stop 'Input error: expecting integrator declaration'
               call lc(chardummy)
               select case(trim(adjustl(chardummy)))
                  case('lanczos')
                     int_type_A = 5
                  case('rk45')
                     int_type_A = 1
                  case('dopri')
                     int_type_A = 2
                  case('bs')
                     int_type_A = 3
                  case default
                     write(0,*) 'Integrator ', trim(adjustl(chardummy)), &
                                ' not possible for the A coefficients'
                     stop
               end select
               read(io_inp,*,iostat = f_status) chardummy, eps_psi
               if (f_status .ne. 0) stop 'Input error: expecting integrator declaration'
               call lc(chardummy)
               select case(trim(adjustl(chardummy)))
                  case('rk45')
                     int_type_psi = 1
                  case('dopri')
                     int_type_psi = 2
                  case('bs')
                     int_type_psi = 3
                  case default
                      write(0,*) 'Integrator ', trim(adjustl(chardummy)), &
                                 ' not possible for the SPFs'
                      stop
               end select
            else
               write(0,*) 'Input Error:'
               write(0,'(2a)') '   unknown integrator: ', IntName
               stop
            end if            
         case('')
            continue
         case default
            write(0,*) 'Input instruction:'
            write(0,*) trim(line)
            write(0,*) 'cannot be understood'
            stop
      end select
   end do
   close(io_inp)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Check the consistency of the parameters and initialise arrays
   if (nstates .lt. 1) stop 'Input error: nstates < 1'
   do im = 1, nmodes
      do i = 1, nmodeDOF(im)
         if (iDOF_t(i,im) .gt. ndof) then
            write(0,'(a,i0)') 'Input error: Attempt to assign the dof. ', iDOF_t(i,im)
            write(0,'(a,i0)') '   The maximum dof number is', ndof
            stop
         end if
      end do
   end do
   if (Hamiltonfile .eq. '') stop 'Input error: Hamiltonian file not given'
   nSPFmax = maxval(nSPF)
   allocate(iDOF(maxval(nmodeDOF),nmodes))
   iDOF = iDOF_t(1:maxval(nmodeDOF),:nmodes)
   deallocate(iDOF_t)
   ! a mode cannot have simultaneously dvr and gwp representation
   do im = 1, nmodes
      nmD = nmodeDOF(im)
      if (all(dvrtype(iDOF(1:nmD,im)) .eq. 'gwp') .or. all(dvrtype(iDOF(1:nmD,im)) .ne. 'gwp')) cycle
      write(0,'(a,i0)') 'Input error for mode ', im
      write(0,*) '   A mode cannot have simultaneously dvr and gwp representation'
      stop
   end do
      ! ... other checks might be required
!
   ! Define the array modeOfDOF
   do im = 1, nmodes
      do i = 1, nmodeDOF(im)
         modeOfDOF(iDOF(i,im)) = im
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   ! Write data in the log file
   write(ilog_unit,'(a,i0)') 'Number of electronic states: ', nstates
   write(ilog_unit,'(a,i0)') 'Number of degrees of freedom: ', ndof
   write(ilog_unit,'(a,i0)') 'Number of logical modes: ', nmodes
   !
   write(ilog_unit,*)
   write(ilog_unit,*) 'SPF structure'
   do im = 1, nmodes
      nmD = nmodeDOF(im)
      write(wrfmt,'(a,i0,a,i0,a)') '(3x,a,i0,a,', nmD, '(i0,2x),a,', nstates, '(i0,2x))'
      write(ilog_unit,trim(wrfmt)) &
         'Mode ', im, ' -     DOFs: ', (iDOF(i,im), i = 1, nmD), &
         ' ;    SPFs: ', (nSPF(im,i), i = 1, nstates)
   end do
   !
   write(ilog_unit,*)
   write(ilog_unit,*) 'Primitive basis representations'
   do i = 1, ndof
      write(ilog_unit,'(3x,a,i0,a,a,3x,i0,2(4x,f10.4))') &
         'DOF ', i, ' - ', dvrtype(i), ngp(i), dvrpar(:,i)
   end do
   !
   write(ilog_unit,*)
   write(ilog_unit,*) 'The Hamiltonian is read from file: ', Hamiltonfile
   ! Read the operators
      ! Allocate the dissipative operators
   if (lDiss) allocate(aeD1(nDiss1), aeD2(nDiss2))
      ! Allocate the Harmonic bath frequencies
   if (lHObath) allocate(wBath(nHObath))
      ! Allocate iOpDOFExtra 'just in case'  (????)
   if (.not. allocated(iOpDOFExtra)) allocate(iOpDOFExtra(ndof,nExtraOpTermsMax,nExtraOpMax))
   call Read_Operators()
   !
   flush(ilog_unit)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine read_input



   subroutine lc(string)
   ! It re-writes a string in lowercase
   character(len = *), intent(inout) :: string
   integer :: i
   integer :: capdiff, lclim, uclim, lslim, sav
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   lclim = iachar('A')
   uclim = iachar('Z')
   lslim = iachar('a')
   capdiff = lclim - lslim
   do i = 1, len(string)
      sav = iachar(string(i:i))
      if ((sav .ge. lclim) .and. (sav.le.uclim)) &
         string(i:i) = achar(sav - capdiff)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine lc



   subroutine read_declare_spf(line,imode)
   ! It read a SPF declaration
   integer :: imode, j, j1, j2, iC
   character(len = 200) :: line, sDOF, sSPF
   character(len = 1) :: set_type
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   j = index(line, '#')
   if (j .ne. 0) line = line(1:j - 1)
   j = index(line, '=')
   if (j .eq. 0) then 
      write(0,*) 'Input error. Wrong SPF declaration:'
      write(0,*) trim(adjustl(line))
      stop
   end if
   sDOF = line(:j - 1)
   sSPF = line(j + 1:)
! Read whether the single-set or the multi-set formalism should be used
   j2 = index(sDOF, ',')
   read(sDOF(:j2 - 1), *, iostat = f_status) set_type
   if (f_status .ne. 0) stop 'Input error: Wrong SPF declaration (MS/SS)'
   if (set_type .eq. 'S') MultiSet(imode) = .false.
! Read the dofs of the i-th mode
   iC = 1   ! iC is the DOF number for the mode imode
   j1 = j2 + 1
   do
      j2 = index(sDOF(j1:), ',')
      if (j2 .eq. 0) then
         read(sDOF(j1:), *, iostat = f_status) iDOF_t(iC,imode)
         if (f_status .ne. 0) stop 'Input error: Wrong SPF declaration (DOF)'
         nmodeDOF(imode) = iC
         exit
      end if
      read(sDOF(j1:j1 + j2 - 1), *, iostat = f_status) iDOF_t(iC,imode)
      if (f_status .ne. 0) stop 'Input error: Wrong SPF declaration (DOF)'
      j1 = j1 + j2  
      iC = iC + 1
   end do
! Read the number of spfs for each state
   iC = 1  ! iC denotes the electronic state
   j1 = 1
   do 
      j2 = index(sSPF(j1:), ',')
      if (j2 .eq. 0) then
         read(sSPF(j1:), *, iostat = f_status) nSPF(imode,iC) 
         if (f_status .ne. 0) stop 'Input error: Wrong SPF declaration (SPF)'
         exit
      end if
      read(sSPF(j1:j1 + j2 - 1), *, iostat = f_status) nSPF(imode,iC)
      if (f_status .ne. 0) stop 'Input error: Wrong SPF declaration (SPF)'
      j1 = j1 + j2
      iC = iC + 1
   end do
! If the modes is described with the single-set formalism, set an equal number of SPFs
! for all the electronic states
   if (.not. MultiSet(imode)) then
      do iC = 2, nstates
         nSPF(imode,iC) = nSPF(imode,1)
      end do
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine read_declare_spf



   subroutine Read_Operators()
   ! Read the Hamiltonian and the additional operators (e.g. for expectation values)
   integer :: i, j, i_unit, f_status, iEF, iExp, iOp, nOpT, iHam, incre, idummy
   double precision :: cR, cI, dummy
   logical :: exists
   character(len = 500) :: line
   character(len = 100), dimension(:) :: OperatorFile(nExtraOpMax)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Array with all operator files
   OperatorFile(1) = Hamiltonfile
   OperatorFile(2:nExtraOp + 1) = ExtraOpFile(1:nExtraOp)
! Preliminary: Check that the operator files exist
   ! Additional operators
   do iOp = 1, nExtraOp + 1
      inquire(file = OperatorFile(iOp), exist = exists)
      if (.not. exists) then
         write(0,*) 'Input error:'
         write(0,*) 'Operator file, ', OperatorFile(iOp), ' not found'
         stop
      end if
   end do 
!
! First loop on the operator files. Look for electronic functions (EF)
!
   do iOp = 1, nExtraOp + 1
      i_unit = freeunit()
      open(unit = i_unit, file = OperatorFile(iOp), status = 'old', action = 'read', form = 'formatted')
      ! Go after the end of the Hamiltonian
      do
         read(i_unit, '(a500)', iostat = f_status) line
         if (f_status .ne. 0) exit
         if (line(1:4) .eq. '----') exit
      end do
      ! Look for electronic functions
      do
         read(i_unit, '(a500)', iostat = f_status) line
         if (f_status .ne. 0) exit
         if (trim(adjustl(line)) .eq. '') cycle
         i = index(line, '#')
         if (i .ne. 0 .and. trim(adjustl(line(1:i - 1))) .eq. '') cycle
         if (i .ne. 0) line = line(1:i - 1)
         i = index(line, 'EF')
         if (i .eq. 0) cycle
         read(line(i + 2:),*) iEF, iExp
         nTermsEF(iEF) = iExp
         do j = 1, iExp
            read(i_unit,*) cR, cI, lEF(j,iEF), rEF(j,iEF)
            if (lEF(j,iEF) .lt. 1 .or. rEF(j,iEF) .lt. 1) stop 'Input error: non-positive electronic label'
            if (lEF(j,iEF) .gt. nstates .or. rEF(j,iEF) .gt. nstates) &
               stop 'Input error: too large electronic label'
            coeffEF(j,iEF) = cR * cOne + cI * ciOne
         end do
      end do
      close(i_unit)
   end do
!
! Second loop on the operator files. Count the number of Hamiltonian terms
!
   ! Hamiltonian
   i_unit = freeunit()
   open(unit = i_unit, file = Hamiltonfile, status = 'old', action = 'read', form = 'formatted')
   nHam = 0
   nHamTD = 0  ! number of time-dependent terms in the Hamiltonian
   do
      line = ''
      read(i_unit, '(a500)', iostat = f_status) line
      if (f_status .ne. 0) exit
      if (line(1:4) .eq. '----') exit
      if (trim(adjustl(line)) .eq. '') cycle
      i = index(line, '#')
      if (i .ne. 0 .and. trim(adjustl(line(1:i - 1))) .eq. '') cycle
      if (i .ne. 0) line = line(1:i - 1)
      call count_ham_terms(line,nHam)   ! update the number of Hamiltonian terms
   end do 
   close(i_unit)
   ! Allocate Hamiltonian variables
   allocate(ae(nHam), lhs_el(nHam), rhs_el(nHam))
   allocate(iOpDof(ndof,nHam))
   iOpDof = 0     ! Set the identity as default for all dofs
   ! Additional operators
   nOpTerms = 0
   do iOp = 1, nExtraOp
      i_unit = freeunit()
      open(unit = i_unit, file = ExtraOpFile(iOp), status = 'old', action = 'read', form = 'formatted')
      ! If this is a dissipative operator, the first line contains a global coefficient
      if (any(iDiss1 == iOp) .or. any(iDiss2 == iOp)) then
         read(i_unit,*) idummy, dummy
         if (idummy == 1) then
            do j = 1, nDiss1
               if (iDiss1(j) == iOp) exit
            end do
            aeD1(j) = dummy * dOneHalf   ! Divide by two here. Check Eq. (34) and (37) of JPC 150, 224106 (2019)
         else if (idummy == 2) then
            do j = 1, nDiss2
               if (iDiss2(j) == iOp) exit
            end do
            aeD2(j) = dummy
         end if
      end if
      ! If this is a Harmonic bath operator coupled to a coordinate Q, the first line contains the frequency
      if (any(iHObathQ == iOp)) then
         read(i_unit,*) dummy
         do j = 1, nHObath
            if (iHObathQ(j) == iOp) exit
         end do
         wBath(j) = dummy
      end if
      ! If this is a Harmonic bath operator coupled to a momentum P, the Q operator should be skipped
      if (any(iHObathP == iOp)) then
         do
            line = ''
            read(i_unit, '(a500)', iostat = f_status) line
            if (line(1:4) .eq. '----') exit
         end do
      end if
      !
      nOpT = 0
      do
         line = ''
         read(i_unit, '(a500)', iostat = f_status) line
         if (f_status .ne. 0) exit
         if (line(1:4) .eq. '----') exit
         if (trim(adjustl(line)) .eq. '') cycle
         i = index(line, '#')
         if (i .ne. 0 .and. trim(adjustl(line(1:i - 1))) .eq. '') cycle
         if (i .ne. 0) line = line(1:i - 1)
         call count_ham_terms(line, nOpT)   ! update the number of operator terms
      end do
      nOpTerms(iOp) = nOpT
      ! If this is a Harmonic bath operator, read a second operator
      close(i_unit)
   end do
   iOpDofExtra = 0   ! Set the identity as default for all dofs
!
! Third loop on the operator files. Read the operators
!
   ! Allocate temporary variables for the parametrized operators
   allocate(OpParam_t(3,(nHam + sum(nOpTerms)) * ndof), OpParamFunc_t((nHam + sum(nOpTerms)) * ndof))
   nOpParam = 0    ! number of elementary parametrized operators
   OpParamFunc_t = 0
   OpParam_t = dZero
   ! Read the Hamiltonian
   i_unit = freeunit()
   open(unit = i_unit, file = Hamiltonfile, status = 'old', action = 'read', form = 'formatted')
   iHam = 1
   do
      line = ''
      read(i_unit, '(a500)', iostat = f_status) line
      if (f_status .ne. 0) exit
      if (line(1:4) .eq. '----') exit
      if (trim(adjustl(line)) .eq. '') cycle      ! cycle because this line is empty
      i = index(line, '#')
      if (i .ne. 0 .and. trim(adjustl(line(1:i - 1))) .eq. '') cycle   ! cycle because this line contains only a comment
      if (i .ne. 0) line = line(1:i - 1)
      call read_declare_ham(line,iHam,incre)
      iHam = iHam + incre
   end do
   close(i_unit)
   ! Read the additional operators
   do iOp = 1, nExtraOp
      i_unit = freeunit()
      open(unit = i_unit, file = ExtraOpFile(iOp), status = 'old', action = 'read', form = 'formatted')
      if (any(iDiss1 == iOp) .or. any(iDiss2 == iOp)) read(i_unit,*) idummy, dummy
      if (any(iHObathQ == iOp)) read(i_unit,*) dummy    ! the first line is the frequency
      !
      if (any(iHObathP == iOp)) then   
         ! skip the operator associated with the coordinate,
         ! so to read the one associated with the momentum
         do
            line = ''
            read(i_unit, '(a500)', iostat = f_status) line
            if (line(1:4) .eq. '----') exit
         end do
      end if
      !
      iHam = 1
      do
         line = ''
         read(i_unit, '(a500)', iostat = f_status) line
         if (f_status .ne. 0) exit
         if (line(1:4) .eq. '----') exit
         if (trim(adjustl(line)) .eq. '') cycle      ! cycle because this line is empty
         i = index(line, '#')
         if (i .ne. 0 .and. trim(adjustl(line(1:i - 1))) .eq. '') cycle   ! cycle because this line contains only a comment
         if (i .ne. 0) line = line(1:i - 1)
         call read_declare_op(line,iHam,incre,iOp)
         iHam = iHam + incre
      end do 
      close(i_unit)
   end do
!
! Copy the parameters into OpParam
!
   allocate(OpParam(3,nOpParam), OpParamFunc(nOpParam))
   call dlacpy('All', 3, nOpParam, OpParam_t, 3, OpParam, 3)
   OpParamFunc = OpParamFunc_t(1:nOpParam) 
   deallocate(OpParam_t, OpParamFunc_t)
!
!  Write the Hamiltonian log
!
   write(ilog_unit,*) '   Hamiltonian terms'
   write(ilog_unit,*)
   write(ilog_unit,'(a)') '  Term  |    Coeff.    |  el. operator  |  dof operators '
   write(ilog_unit,'(a)') '---------------------------------------------------------'
   do iHam = 1, nHam
      write(ilog_unit,'(2x,i4,2x,a1,1x,f12.8,1x,a1,3x,a1,i2,1x,a2,i2,1x,a1,3x,a1,3x,20(i3,1x))') &
         iHam, '|', dble(ae(iHam)), '|', '|', lhs_el(iHam), '><', rhs_el(iHam), '|', &
         '|', (iOpDof(i,iHam), i = 1, ndof)
   end do
   if (nHamTD .gt. 0) then
      write(ilog_unit,*)
      write(ilog_unit,'(a,1x,30(1x,i0))') &
        '  The following Hamiltonian terms are time-dependent:', iHamTD(1:nHamTD)
   end if
!
! Write the additional operator log
!
   do iOp = 1, nExtraOp
      write(ilog_unit,*)
      write(ilog_unit,*) '   Terms for the operator of file ', trim(ExtraOpFile(iOp))
      write(ilog_unit,*)
      write(ilog_unit,'(a)') '  Term  |   Re[Coeff.]    Im[Coeff.] |  el. operator  |  dof operators '
      write(ilog_unit,'(a)') '---------------------------------------------------------'
      do iHam = 1, nOpTerms(iOp)
         write(ilog_unit,'(2x,i4,2x,a1,1x,f12.8,2x,f12.8,1x,a1,3x,a1,i2,1x,a2,i2,1x,a1,3x,a1,3x,20(i4,1x))') &
            iHam, '|', real(aeExtra(iHam,iOp)), imag(aeExtra(iHam,iOp)), &
            '|', '|', lhs_elExtra(iHam,iOp), '><', rhs_elExtra(iHam,iOp), '|', &
            '|', (iOpDofExtra(i,iHam,iOp), i = 1, ndof)
      end do
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine Read_Operators



   subroutine count_ham_terms(line,nterms)
   integer, intent(inout) :: nterms
   integer :: i, j, k, iel1, iel2, iEF
   character(len = 500), intent(in) :: line
   character(len = 500) :: line2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   line2 = line
   i = index(line2, '|')
   if (i .eq. 0) then  ! only the coefficient must be read --> identity operator for the electronic dof
      nterms = nterms + nstates
      return
   end if
   line2 = line2(i + 1:)
   i = index(line2, '|')
   if (i .ne. 0) line2 = line2(1:i - 1)
   if (trim(adjustl(line2)) .eq. '1') then   ! identity operator
      nterms = nterms + nstates
      return
   end if
   j = index(line2, 'S')
   if (j .ne. 0) then
      k = index(line2, '&')
      read(line2(j + 1:k - 1), *) iel1
      read(line2(k + 1:), *) iel2
      if (iel1 .eq. iel2) then
         nterms = nterms + 1
      else
         nterms = nterms + 2
      end if
      return
   end if
   j = index(line2, 'Z')
   if (j .ne. 0) then
      k = index(line2, '&')
      read(line2(j + 1:k - 1), *) iel1
      read(line2(k + 1:), *) iel2
      if (iel1 .eq. iel2) then
         write(0,*) 'Input Error:'
         write(0,*) ' Electronic operators "Za&b" are allowed only between different states'
      else
         nterms = nterms + 1
      end if
      return
   end if
   j = index(line2, 'EF')
   if (j .ne. 0) then
      read(line2(j + 2:),*) iEF
      nterms = nterms + nTermsEF(iEF)
   end if
   if (j .eq. 0) then
      write(0,*) 'Input Error:'
      write(0,*) 'Only operator labels "1", "Sa&b" or "Za&b" are allowed'
      stop
   end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine count_ham_terms



   subroutine read_declare_ham(line,iHam,incre)
   ! It reads a line from the Hamiltonian input
   ! It generates a number of Hamiltonian terms, and gives the increment for iHam
   integer, intent(in) :: iHam
   integer, intent(out) :: incre
   integer :: i, j, i1, i2, j1, j2, kdof, iC, ind, iEF, idummy, iC1, iC2
   character(len = 500), intent(in) :: line
   character(len = 500) :: line2, line3
   double precision :: coeff, p1, p2, p3, p4
   character(len = 100) :: InterpDataFile
   double complex :: cdummy
!   integer, dimension(:), allocatable :: OpParamFunc_t
!   double precision, dimension(:,:), allocatable :: OpParam_t
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   line2 = line
! Read the coefficient
   i = index(line2, '|')
   if (i .eq. 0) then
      line3 = line2 
   else
      line3 = line2(1:i - 1)
   end if
   read(line3, *, iostat = f_status) coeff
   ! Check if it is imaginary
   i1 = index(line3, 'i')
   if (i1 .eq. 0) then
      cdummy = cOne
   else
      cdummy = ciOne
   endif
   !
   if (f_status .ne. 0) call err_read('Hamiltonian coeff.', 'read_declare_ham', f_status)
   if (i .eq. 0) then   !  only the coefficient must be read
      incre = nstates
      forall(j = 1:nstates)
         lhs_el(iHam + j - 1) = j
         rhs_el(iHam + j - 1) = j
         ae(iHam + j - 1) = coeff * cdummy
      end forall
      return
   end if
   line2 = line2(i + 1:)
! Read the electronic state operator
   i = index(line2, '|')
   if (i .eq. 0) then
      line3 = line2
   else
      line3 = line2(1:i - 1)
   end if
   if (trim(adjustl(line3)) .eq. '1') then
      incre = nstates
      forall(j = 1:nstates)
         lhs_el(iHam + j - 1) = j
         rhs_el(iHam + j - 1) = j
         ae(iHam + j - 1) = coeff * cdummy
      end forall
   elseif(index(line3,'S') .ne. 0) then
      i1 = index(line3,'S')
      i2 = index(line3,'&')
      read(line3(i1 + 1:i2 - 1),*) j1
      read(line3(i2 + 1:),*) j2
      if (j1 .eq. j2) then
         incre = 1
         lhs_el(iHam) = j1
         rhs_el(iHam) = j2
         ae(iHam) = coeff * cdummy
      else
         incre = 2
         lhs_el(iHam) = j1
         rhs_el(iHam) = j2
         lhs_el(iHam + 1) = j2
         rhs_el(iHam + 1) = j1
         ae(iHam) = coeff * cdummy
         ae(iHam + 1) = coeff * cdummy
      end if 
   elseif(index(line3,'Z') .ne. 0) then
      i1 = index(line3,'Z')
      i2 = index(line3,'&')
      read(line3(i1 + 1:i2 - 1),*) j1
      read(line3(i2 + 1:),*) j2
      incre = 1
      lhs_el(iHam) = j1
      rhs_el(iHam) = j2
      ae(iHam) = coeff * cdummy
   elseif(index(line3,'EF') .ne. 0) then
      i1 = index(line3,'EF')
      read(line3(i1 + 2:),*) iEF
      incre = nTermsEF(iEF)
      do j = 0, incre - 1
         lhs_el(iHam + j) = lEF(j + 1,iEF)
         rhs_el(iHam + j) = rEF(j + 1,iEF)
         ae(iHam + j) = coeff * cdummy * coeffEF(j + 1,iEF)
      end do   
   else 
      write(0,*) 'Input Error:'
      write(0,*) 'Only operator labels "1", "Sa&b" or "Za&b" are allowed'
      stop
   end if  
   if (i .eq. 0) return
   line2 = line2(i:)
! Read the dof operators
   do 
      i = index(line2,'|')
      j = index(line2(i + 1:), '|')
      if (j .eq. 0) then
         line3 = line2(i + 1:)
      else
         line3 = line2(i + 1:i + j - 1)
      end if
      read(line3(1:2),*, iostat = f_status) kdof
      if (f_status .ne. 0) call err_read('Hamiltonian dof', 'read_declare_ham', f_status)
      if (kdof .eq. 0) then ! Time dependent coefficient
         ind = index(line3(3:),'gau')
         if (ind .eq. 0) then
            write(0,*) 'Input Error: only time-dependent function "gau" allowed'
            stop
         end if
         read(line3(ind + 5:),*) p1, p2, p3, p4
         do iC = 1, incre
            nHamTD = nHamTD + 1 
            iHamTD(nHamTD) = iHam + iC - 1
            typeHamTD(nHamTD) = 1
            parHamTD(0,nHamTD) = coeff
            parHamTD(1,nHamTD) = p1 / au2fs 
            parHamTD(2,nHamTD) = p2 / au2fs
            parHamTD(3,nHamTD) = p3 / au2cmm1
            parHamTD(4,nHamTD) = p4 / 180 * pi
         end do
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      if (trim(adjustl(line3(3:))) .eq. '1') then   ! identity is the default
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      !
      ! Check if parametrized functions are used
      !
         ! Gaussian
      ind = index(line3(3:), 'gau')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
         read(line3(ind + 5:),*) p1, p2   ! exp(p1*x + p2*x**2)
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 0
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
         ! 1D-interpolating function
      ind = index(line3(3:), 'interp')
      if (ind .ne. 0) then
         ! Check whether this function has been already read
         read(line3(ind + 8:),*) InterpDataFile
         if (any(strings == InterpDataFile)) then
            do iC1 = 1, nStrings
               if (strings(iC1) == InterpDataFile) exit
            end do
            do iC2 = 1, nOpParam
               if (OpParamFunc_t(iC2) == 1 .and. OpParam_t(1,iC2) == dble(iC1)) exit
            end do
            iOpDof(kdof,iHam:iHam + incre - 1) = 200 + iC2 
         else   ! It is a new interpolation file
            nOpParam = nOpParam + 1
            iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
            nStrings = nStrings + 1
            if (nStrings .gt. nStringsMax) then
                write(0,*) 'ERROR: too many strings to be saved'
                stop
            end if
            strings(nStrings) = InterpDataFile
            OpParam_t(1,nOpParam) = dble(nStrings)
            OpParamFunc_t(nOpParam) = 1
         end if
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
         ! indicator function Ind(x0 - dx/2, x0 + dx/2)
      ind = index(line3(3:), 'square')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
         read(line3(ind + 8:),*) p1, p2   ! x0, dx
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 2
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
         ! "Hermite" function  x**r * exp(p1*x + p2*x**2)
      ind = index(line3(3:), 'hermite')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
         read(line3(ind + 9:),*) idummy, p1, p2   ! idummy is the power r
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParam_t(3,nOpParam) = dble(idummy)
         OpParamFunc_t(nOpParam) = 3
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! sin^-2
      ind = index(line3(3:), 'sin^-2')
      if (ind .ne. 0) then
         if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: sin^-2 operator allowed only for Leg-DVR'
         iOpDof(kdof,iHam:iHam + incre - 1) = -4
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! Powers of sine ( sin(p1 * x)^p2 )
      ind = index(line3(3:), 'sin')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
         if (len(trim(line3(ind + 5:))) == 0) then
            p1 = 1.0
            p2 = 1.0
         else
            read(line3(ind + 5:),*) p1, p2
         end if
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 4
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! Powers of sine ( cos(p1 * x)^p2 )
      ind = index(line3(3:), 'cos')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDof(kdof,iHam:iHam + incre - 1) = 200 + nOpParam
         if (len(trim(line3(ind + 5:))) == 0) then
            p1 = 1.0
            p2 = 1.0
         else
            read(line3(ind + 5:),*) p1, p2
         end if
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 5
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
         ! Projector onto an eigenstate
      ind = index(line3(3:), 'eig_proj')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         read(line3(ind + 10:),*) InterpDataFile, OpParam_t(2,nOpParam), idummy    ! Potential file, kinetic energy prefactor, energy level
         OpParam_t(3,nOpParam) = dble(idummy)
         iOpDof(kdof,iHam:iHam + incre - 1) = -200 - nOpParam
         nStrings = nStrings + 1
         if (nStrings .gt. nStringsMax) then
             write(0,*) 'ERROR: too many strings to be saved'
             stop
         end if
         strings(nStrings) = InterpDataFile
         OpParam_t(1,nOpParam) = dble(nStrings)
         OpParamFunc_t(nOpParam) = 6
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
!         ! Copy the parameters into OpParam
!      allocate(OpParam(3,nOpParam), OpParamFunc(nOpParam))
!      call dlacpy('All', 3, nOpParam, OpParam_t(1,1), 3, OpParam(1,1), 3)
!      OpParamFunc = OpParamFunc_t(1:nOpParam)
!         ! Deallocate the large array OpParam_t
!      deallocate(OpParam_t, OpParamFunc_t)
      !
      select case(trim(adjustl(line3(3:))))
         case('p')
            iOpDof(kdof,iHam:iHam + incre - 1) = -1
         case('dq^2')
            iOpDof(kdof,iHam:iHam + incre - 1) = -2
         case('L0^2')  ! Second derivative for the Legendre DVR
            if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: L0^2 operator allowed only for Leg-DVR'
            iOpDof(kdof,iHam:iHam + incre - 1) = -2
         case('sin^-2')  !
            if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: sin^-2 operator allowed only for Leg-DVR'
            iOpDof(kdof,iHam:iHam + incre - 1) = -4
         case('qp')  ! qp + pq
            iOpDof(kdof,iHam:iHam + incre - 1) = -3
         case('q')
            iOpDof(kdof,iHam:iHam + incre - 1) = 1
         case('q^2')
            iOpDof(kdof,iHam:iHam + incre - 1) = 2
         case('q^3')
            iOpDof(kdof,iHam:iHam + incre - 1) = 3
         case('q^4')
            iOpDof(kdof,iHam:iHam + incre - 1) = 4
         case('q^5')
            iOpDof(kdof,iHam:iHam + incre - 1) = 5
         case('q^6')
            iOpDof(kdof,iHam:iHam + incre - 1) = 6
         case('q^7')
            iOpDof(kdof,iHam:iHam + incre - 1) = 7
         case('q^8')
            iOpDof(kdof,iHam:iHam + incre - 1) = 8
         case('q^9')
            iOpDof(kdof,iHam:iHam + incre - 1) = 9
         case('q^10')
            iOpDof(kdof,iHam:iHam + incre - 1) = 10
         case default
            write(0,*) 'Input Error:'
            write(0,*) line
            write(0,*) '  dof ', kdof
            write(0,'(3a)') 'Operator ', trim(adjustl(line3(3:))), ' not implemented'
            stop
      end select
      !
      if (j .eq. 0) return
      line2 = line2(i + j:)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine read_declare_ham

   


   subroutine read_declare_op(line,iOp,incre,iOpExtra)
   ! It reads a line from the Hamiltonian input
   ! It generates a number of Hamiltonian terms, and gives the increment for iHam
   integer, intent(in) :: iOp, iOpExtra
   integer, intent(out) :: incre
   character(len = 500), intent(in) :: line
   integer :: i, j, i1, i2, j1, j2, kdof, ind, idummy, iC1, iC2
   character(len = 500) :: line2, line3
   double precision :: dcoeff, p1, p2
   double complex :: coeff, cdummy
   character(len = 100) :: InterpDataFile
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   line2 = line
! Read the coefficient
   i = index(line2, '|')
   if (i .eq. 0) then
      line3 = line2 
   else
      line3 = line2(1:i - 1)
   end if
   read(line3, *, iostat = f_status) dcoeff
   ! Check if it is imaginary
   i1 = index(line3, 'i')
   if (i1 .eq. 0) then
      cdummy = cOne
   else
      cdummy = ciOne
   endif
   coeff = dcoeff * cdummy
   !
   if (i .eq. 0) then   !  only the coefficient must be read
      incre = nstates
      forall(j = 1:nstates)
         lhs_elExtra(iOp + j - 1,iOpExtra) = j
         rhs_elExtra(iOp + j - 1,iOpExtra) = j
         aeExtra(iOp + j - 1,iOpExtra) = coeff
      end forall
      return
   end if
   line2 = line2(i + 1:)
! Read the electronic state operator
   i = index(line2, '|')
   if (i .eq. 0) then
      line3 = line2
   else
      line3 = line2(1:i - 1)
   end if
   if (trim(adjustl(line3)) .eq. '1') then
      incre = nstates
      forall(j = 1:nstates)
         lhs_elExtra(iOp + j - 1,iOpExtra) = j
         rhs_elExtra(iOp + j - 1,iOpExtra) = j
         aeExtra(iOp + j - 1,iOpExtra) = coeff
      end forall
   elseif(index(line3,'S') .ne. 0) then
      i1 = index(line3,'S')
      i2 = index(line3,'&')
      read(line3(i1 + 1:i2 - 1),*) j1
      read(line3(i2 + 1:),*) j2
      if (j1 .eq. j2) then
         incre = 1
         lhs_elExtra(iOp,iOpExtra) = j1
         rhs_elExtra(iOp,iOpExtra) = j2
         aeExtra(iOp,iOpExtra) = coeff
      else
         incre = 2
         lhs_elExtra(iOp,iOpExtra) = j1
         rhs_elExtra(iOp,iOpExtra) = j2
         lhs_elExtra(iOp + 1,iOpExtra) = j2
         rhs_elExtra(iOp + 1,iOpExtra) = j1
         aeExtra(iOp,iOpExtra) = coeff 
         aeExtra(iOp + 1,IOpExtra) = coeff
      end if 
   elseif(index(line3,'Z') .ne. 0) then
      i1 = index(line3,'Z')
      i2 = index(line3,'&')
      read(line3(i1 + 1:i2 - 1),*) j1
      read(line3(i2 + 1:),*) j2
      incre = 1
      lhs_elExtra(iOp,iOpExtra) = j1
      rhs_elExtra(iOp,iOpExtra) = j2
      aeExtra(iOp,iOpExtra) = coeff
   else 
      write(0,*) 'Input Error:'
      write(0,*) 'Only operator labels "1", "Sa&b" or "Za&b" are allowed'
      stop
   end if  
   if (i .eq. 0) return
   line2 = line2(i:)
! Read the dof operators
   do 
      i = index(line2,'|')
      j = index(line2(i + 1:), '|')
      if (j .eq. 0) then
         line3 = line2(i + 1:)
      else
         line3 = line2(i + 1:i + j - 1)
      end if
      read(line3(1:2),*, iostat = f_status) kdof
      if (f_status .ne. 0) call err_read('Hamiltonian dof', 'read_declare_ham', f_status)
      if (trim(adjustl(line3(3:))) .eq. '1') then   ! identity is the default
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      !
      ! Check if parametrized functions are used
      !
      ! Gaussian
      ind = index(line3(3:), 'gau')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
         read(line3(ind + 5:),*) p1, p2   ! exp(p1*x + p2*x**2)
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 0
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! 1D-interpolating function
      ind = index(line3(3:), 'interp')
      if (ind .ne. 0) then
         read(line3(ind + 8:),*) InterpDataFile
         ! check whether the function has already beed found
         if (any(strings == InterpDataFile)) then
            do iC1 = 1, nStrings
               if (strings(iC1) == InterpDataFile) exit
            end do
            do iC2 = 1, nOpParam
               if (OpParamFunc_t(iC2) == 1 .and. OpParam_t(1,iC2) == dble(iC1)) exit
            end do
            iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + iC2 
         else   ! It is a new interpolation file
            nOpParam = nOpParam + 1
            iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
            nStrings = nStrings + 1
            if (nStrings .gt. nStringsMax) then
                write(0,*) 'ERROR: too many strings to be saved'
                stop
            end if
            strings(nStrings) = InterpDataFile
            OpParam_t(1,nOpParam) = dble(nStrings)
            OpParamFunc_t(nOpParam) = 1
         end if
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! indicator function for the interval [x0-dx/2,x0+dx/2]
      ind = index(line3(3:), 'square')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
         read(line3(ind + 8:),*) p1, p2   ! x0, dx
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 2
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! "Hermite" function  x**r * exp(p1*x + p2*x**2)
      ind = index(line3(3:), 'hermite')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
         read(line3(ind + 9:),*) idummy, p1, p2   ! idummy is the power r
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParam_t(3,nOpParam) = dble(idummy)
         OpParamFunc_t(nOpParam) = 3
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! sin^-2
      ind = index(line3(3:), 'sin^-2')
      if (ind .ne. 0) then
         if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: sin^-2 operator allowed only for Leg-DVR'
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -4
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! Powers of sine ( sin(p1 * x)^p2 )
      ind = index(line3(3:), 'sin')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
         if (len(trim(line3(ind + 5:))) == 0) then
            p1 = 1.0
            p2 = 1.0
         else
            read(line3(ind + 5:),*) p1, p2
         end if
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 4
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! Powers of sine ( cos(p1 * x)^p2 )
      ind = index(line3(3:), 'cos')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 200 + nOpParam
         if (len(trim(line3(ind + 5:))) == 0) then
            p1 = 1.0
            p2 = 1.0
         else
            read(line3(ind + 5:),*) p1, p2
         end if
         OpParam_t(1,nOpParam) = p1
         OpParam_t(2,nOpParam) = p2
         OpParamFunc_t(nOpParam) = 5
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      ! Projector onto an eigenstate
      ind = index(line3(3:), 'eig_proj')
      if (ind .ne. 0) then
         nOpParam = nOpParam + 1
         read(line3(ind + 10:),*) InterpDataFile, OpParam_t(2,nOpParam), idummy    ! Potential file, kinetic energy prefactor, energy level
         OpParam_t(3,nOpParam) = dble(idummy)
         iOpDofExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -200 - nOpParam
         nStrings = nStrings + 1
         if (nStrings .gt. nStringsMax) then
             write(0,*) 'ERROR: too many strings to be saved'
             stop
         end if
         strings(nStrings) = InterpDataFile
         OpParam_t(1,nOpParam) = dble(nStrings)
         OpParamFunc_t(nOpParam) = 6
         !
         if (j .eq. 0) return
         line2 = line2(i + j:)
         cycle
      end if
      !
      select case(trim(adjustl(line3(3:))))
         case('p')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -1
         case('dq^2')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -2
         case('L0^2')  ! Second derivative for the Legendre DVR
            if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: L0^2 operator allowed only for Leg-DVR'
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -2
         case('sin^-2') 
            if (dvrType(kdof) .ne. 'Leg') stop 'ERROR: sin^-2 operator allowed only for Leg-DVR'
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -4
         case('qp')  ! qp + pq
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = -3
         case('q')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 1
         case('q^2')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 2
         case('q^3')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 3
         case('q^4')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 4
         case('q^5')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 5
         case('q^6')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 6
         case('q^7')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 7
         case('q^8')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 8
         case('q^9')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 9
         case('q^10')
            iOpDOFExtra(kdof,iOp:iOp + incre - 1,iOpExtra) = 10
         case default
            write(0,*) 'Input Error:'
            write(0,*) line
            write(0,*) '  dof ', kdof
            write(0,'(3a)') 'Operator ', trim(adjustl(line3(3:))), ' not implemented'
            stop
      end select
      !
      if (j .eq. 0) return
      line2 = line2(i + j:)
   end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine read_declare_op


end module input
