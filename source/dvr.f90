module dvr
   use sysparam
   use globals
   use psidef
   implicit none
   private
   public :: dvrType, ngp, grid, wgrid, idvrpar, dvrpar, operDVR, d1DVR, xd1DVR, d2DVR, trafoDVR, gridPow, gridOper
   public :: ngpm, ngpMax, ngpmMax
   public :: read_declare_dvr, sin_DVR, harm_DVR, Leg_DVR, twopi_DVR, Eigenf1DPot, Eigenf1DPot_real
   !
   character(6), dimension(:), allocatable :: dvrType
   integer :: ngpMax, ngpmMax   ! maximum values of ngp and ngpm
   integer, dimension(:), allocatable :: ngp, &  ! ngp(i) = no. of grid points for the dof i
                                         ngpm    ! ngpm(i) = no. of direct product grid points for the DVR mode i
   integer, dimension(:), allocatable :: idvrpar    ! dvrpar(i) = value of the additional integer parameter for the DVR of the DOF i 
   double precision, dimension(:,:), allocatable :: grid  , &  ! grid(i,j) = value i-th grid point for the DVR-dof j
                                                    wgrid , &  ! wgrid(i,j) = DVR weight at the i-th grid point for the DVr-dof j
                                                    dvrpar     ! dvrpar(i,j) = value of the i-th parameter for the dvr of dof j
   double precision, dimension(:,:,:), allocatable :: d1DVR    ! d1DVR(:,:,i)  = dvr first derivative matrix for the DVR-dof i
   double precision, dimension(:,:,:), allocatable :: xd1DVR   ! xd1DVR(:,:,i) = dvr matrix for the operator x*dx + dx*x for the DVR-dof i
   double precision, dimension(:,:,:), allocatable :: d2DVR    ! d2DVR(:,:,i)  = dvr second derivative matrix for the DVR-dof i
   double precision, dimension(:,:,:), allocatable :: operDVR  ! extraDVR(:,:,i)  = dvr matrix for additional special operators (e.g. sin^-2)
   double precision, dimension(:,:,:), allocatable :: trafoDVR ! trafoDVR(:,:,i) = DVR-to-FBR transformation matrix for the DVR-dof i
   !
   double precision, dimension(:,:,:), allocatable :: gridPow  !  gridPow(i,j,k) = q(i)^j for the DVR-dof k
   double precision, dimension(:,:,:), allocatable :: gridOper ! gridOper(i,j,k) = Operator(j) evaluated on the grid point i for the DVR-dof k
   !   
contains



   subroutine read_declare_dvr(line)
   ! It reads a line of DVR declaration in the input file
   use sysparam
   use psidef
   implicit none
   character(len = 200) :: line
   integer :: j, kdof
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   j = index(line, '#')
   if (j .ne. 0) line = line(1:j - 1)
! Look for the available dvr types
   ! sin
   j = index(line, 'sin')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = 'sin'
      read(line(j + 3:),*) ngp(kdof), dvrpar(1,kdof), dvrpar(2,kdof)
      return
   endif
   ! harmonic oscillator
   j = index(line, 'HO')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = 'HO'
      read(line(j + 3:),*) ngp(kdof), dvrpar(1,kdof), dvrpar(2,kdof)
      return
   endif
   ! Legendre
   j = index(line, 'Leg')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = 'Leg'
      read(line(j + 3:),*) idvrpar(kdof), ngp(kdof)
      idvrpar(kdof) = abs(idvrpar(kdof))
      return
   end if
   ! 2-pi periodic
   j = index(line, '2pi')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = '2pi'
      read(line(j + 3:),*) ngp(kdof)
      return
   end if
   ! parameters on a grid
   j = index(line, 'par')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = 'par'
      read(line(j + 3:),*) ngp(kdof), dvrpar(1,kdof), dvrpar(2,kdof)
      return
   endif
   ! gwp
   j = index(line, 'gwp')
   if (j .ne. 0) then
      read(line(1:j - 1),*) kdof
      dvrType(kdof) = 'gwp'
      if (trim(adjustl(line(j + 3:))) .eq. '') then
         dvrpar(1,kdof) = -10.d0
         dvrpar(2,kdof) = 10.d0
         ngp(kdof) = 101
      else
         read(line(j + 3:),*) ngp(kdof), dvrpar(1,kdof), dvrpar(2,kdof)
      endif
      return
   endif
   ! no dvr has been found
   write(0,*) 'Input error. Wrong DVR declaration:'
   write(0,*) trim(adjustl(line))
   stop
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   end subroutine



   subroutine sin_DVR(nP,xi,xf,gridP,wP,d1Mat,d2Mat,xd1Mat,trafo,LDD)
   implicit none
   integer :: i, j
   integer, intent(in) :: nP, LDD
   double precision, intent(in) :: xi, xf
   double precision, dimension(*), intent(out) :: gridP, wP
   double precision, dimension(LDD,*), intent(out) :: d1Mat, d2Mat, xd1Mat, trafo
   double precision :: dx, dummy, dNP1
   double precision, parameter :: piHalf = dOneHalf * pi
   double precision, dimension(nP,nP) :: aux, aux1
   double precision, dimension(nP) :: energ
   double precision, dimension(3 * nP - 1) :: WORK
   integer :: info
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (nP .lt. 2) stop 'too few grid points selected for a DVR'
   ! Generate the grid and the weights
   dx = (xf - xi) / dble(nP - 1)
   dummy = xi
   do i = 1, nP
      gridP(i) = dummy
      wP(i) = dx
      dummy = dummy + dx
   enddo
   ! Symmetrize in case the grid is centered at zero
   if (abs(abs(xf) - abs(xi)) < epsilon(xf)) then
      do i = 1, nP / 2
         dummy = dOneHalf * (- gridP(i) + gridP(nP - i + 1))
         gridP(i) = - dummy
         gridP(nP - i + 1) =   dummy
      end do
      if (mod(nP,2) == 1) gridP(i) = dZero
   end if
   !
   dNP1 = dble(nP + 1)
   ! Generate the matrix of first derivatives
   forall(i = 1:nP, j = 1:nP, i .ne. j)
      trafo(i,j) = sin(i * j * pi / dNP1) * sqrt(2.0d0 / dNP1)
      aux(i,j) = mod(abs(i - j),2) * i * j / dble(i ** 2 - j ** 2)
   end forall
   forall (i = 1:nP) 
      trafo(i,i) = sin(i ** 2 * pi / dNP1) * sqrt(2.0d0 / dNP1)
      aux(i,i) = dZero
   end forall
   call dsymm('R','U',nP,nP,4.d0 / (xf - xi),trafo,LDD,aux,nP,dZero,aux1,nP)
   call dsymm('L','U',nP,nP,dOne,trafo,LDD,aux1,nP,dZero,d1Mat,LDD)
   ! Calculate the matrix of the operator x * dx + dx * x
   forall (i = 1:nP, j = 1:nP) xd1Mat(i,j) = d1Mat(i,j) * (gridP(i) + gridP(j))
   ! Generate the matrix of second derivatives
   forall(i = 1:nP) d2Mat(i,i) = &
      dble(2 * nP ** 2 + 4 * nP + 3) / 3.d0 - dOne / sin(pi * i / dNP1) ** 2
   do i = 1, nP - 1
      dummy = - dOne
      do j = i + 1, nP
         d2Mat(i,j) = ( dOne / sin(piHalf * (i - j) / dNP1) ** 2  &
                      - dOne / sin(piHalf * (i + j) / dNP1) ** 2) &
                      * dummy
         !d2Mat(j,i) = d2Mat(i,j)
         dummy = - dummy
      enddo
   enddo
   dummy = dOneHalf * (pi / (xf - xi + 2 * dx)) ** 2
   forall(i = 1:nP, j = 1:nP, j >= i) d2Mat(i,j) = - d2Mat(i,j) * dummy
   forall(i = 1:nP, j = 1:nP, j <  i) d2Mat(i,j) = d2Mat(j,i)
   return
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Test
   dummy = 0.13
   forall (i = 1:nP, j = 1:nP) aux(i,j) = -dOneHalf * dummy * d2Mat(i,j)
   forall (i = 1:nP) aux(i,i) = aux(i,i) + 1.7 + dOneHalf * dummy * (gridP(i) - 2.d0) ** 2 !- 0.002 * gridP(i) + 0.0001*gridP(i)**4
   call dsyev('V', 'U', nP, aux, nP, energ, WORK, 3 * nP - 1, info)
   write(0,*)i
   write(0,*)
   do i = 1, 10
      write(0,'(i3,3x,f9.6)') i, energ(i)
   enddo
   stop
   return
   end subroutine sin_DVR


   subroutine twopi_DVR(nP,gridP,wP,d1Mat,d2Mat,trafo,LDD)
   ! 2pi-periodic DVR with nP points
   implicit none
   integer :: i, j
   integer, intent(in) :: nP, LDD
   double precision, dimension(*), intent(out) :: gridP, wP
   double precision, dimension(LDD,*), intent(out) :: d1Mat, d2Mat, trafo
   double precision :: dummy, dx, alpha, dNP
   double precision, dimension(nP,nP) :: aux
   double precision, external :: dnrm2
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (nP .lt. 2) stop 'too few grid points selected for a DVR'
   if (mod(nP,2) == 0) stop 'An odd number of points must be requested for the 2pi-DVR'
   ! Generate the grid and the weights
   dNP = dble(nP)
   dx = (pi + pi) / dNP
   dummy = dx * dOneHalf - pi
   do i = 1, (nP - 1) / 2
      gridP(i) = dummy 
      wP(i) = dx
      gridP(nP - i + 1) = - dummy
      wP(nP - i + 1) = dx
      dummy = dummy + dx
   end do
   gridP((nP + 1) / 2) = dZero
   wP((nP + 1) / 2) = dx 
   ! Matrix of second and first derivatives
   forall(i = 1:nP) d1Mat(i,i) = dZero
   forall(i = 1:nP) d2Mat(i,i) = - dble(nP ** 2  - 1) / 12.d0
   do i = 1, nP - 1
      dummy = dOne
      do j = i + 1, nP
         alpha = pi * dble(i - j) / dNP
         d1Mat(i,j) = - dOneHalf * dummy / sin(alpha)
         d2Mat(i,j) = dOneHalf * cos(alpha) / sin(alpha)**2 * dummy
         dummy = - dummy
      enddo
   enddo
   forall(i = 1:nP, j = 1:nP, j < i) 
      d1Mat(i,j) = - d1Mat(j,i)
      d2Mat(i,j) =   d2Mat(j,i)
   end forall
   ! Transformation matrix
   forall (j = 1:nP) aux(j,1) = dOne / sqrt(dNP)
   do i = 1, (nP - 1) / 2
      forall(j = 1:nP) 
         aux(j,2 * i    ) = sin(gridP(j) * dble(i))
         aux(j,2 * i + 1) = cos(gridP(j) * dble(i))
      end forall
      dummy = dnrm2(nP, aux(1,2 * i), 1)
      call drscl(nP,dummy,aux(1,2 * i),1)
      dummy = dnrm2(nP, aux(1,2 * i + 1), 1)
      call drscl(nP,dummy,aux(1,2 * i + 1),1)
   end do
   !forall(i = 1:nP, j = 1:nP, j >= i) aux(i,j) = - d2Mat(i,j)
   !call dsyev('V','U',nP,aux,LDD,W,WORK,3 * nP - 1,info)
   !! Transpose
   do i = 1, nP
      do j = 1, nP
         trafo(i,j) = aux(j,i)
      end do
   end do
   !write(0,*)
   !write(0,*) W(1:nP)
   !write(0,*)
  !write(0,*) gridP(1:nP)
   !do i = 1, nP
   !   write(1111,*) gridP(i), trafo(1,i), trafo(2,i), trafo(3,i), trafo(4,i), trafo(5,i)
   !end do
   !stop
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   return
   end subroutine twopi_DVR



   subroutine harm_DVR(nP,x0,mw,gridP,wP,d1Mat,d2Mat,xd1Mat,trafo,LDD)
   ! Harmonic DVR
   implicit none
   integer :: i, info, j
   integer, intent(in) :: nP, LDD
   double precision, intent(in) :: x0, mw   ! center of the grid, product mass * frequency
   double precision, dimension(*), intent(out) :: gridP, wP
   double precision, dimension(LDD,*), intent(out) :: d1Mat, d2Mat, xd1Mat, trafo
   double precision :: dummy
   double precision, parameter :: piHalf = dOneHalf * pi
   double precision, dimension(3 * nP - 1) :: WORK
   double precision, dimension(nP,nP) :: aux, aux1
   double precision, dimension(nP) :: energ
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
   if (nP .lt. 2) stop 'too few grid points selected for a DVR'
   ! Generate tihe grid and the weights
      ! aux is the matrix representation of the x operator
   call dlaset('U', nP, nP, dZero, dZero, trafo(1,1), LDD)
   forall (i = 1:nP - 1) trafo(i,i + 1) = sqrt(dOneHalf * i / mw)
   call dsyev('V', 'U', nP, trafo(1,1), LDD, gridP(1), WORK, 3 * nP - 1, info)
       ! Symmetrize the grid
   do i = 1, nP / 2
      dummy = dOneHalf * (- gridP(i) + gridP(nP - i + 1))
      gridP(i) = - dummy
      gridP(nP - i + 1) =   dummy
   end do
   if (mod(nP,2) == 1) gridP(i) = dZero
   forall (i = 1:nP) gridP(i) = gridP(i) + x0
      ! Weights
      ! Make the weights positive
   forall (i = 1:nP, trafo(1,i) .lt. dZero) trafo(:,i) = - trafo(:,i)
   forall (i = 1:nP)  &
      wP(i) = trafo(1,i) ** 2 * exp(mw * (gridP(i) - x0) ** 2) * sqrt(pi / mw)
   ! Calculate the matrix of the first derivatives
   call dlaset('All', nP, nP, dZero, dZero, aux(1,1), nP)
   forall (i = 1:nP - 1) 
      aux(i,i + 1) =  sqrt(dble(i))
      aux(i + 1,i) = -sqrt(dble(i))
   end forall
   aux = aux * sqrt(dOneHalf * mw)
   call dgemm('N', 'N', nP, nP, nP, dOne, aux(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1), nP)
   call dgemm('T', 'N', nP, nP, nP, dOne, trafo(1,1), LDD, aux1(1,1), nP, dZero, d1Mat(1,1), LDD)
   ! Calculate the matrix representation of the operator x * dx + dx * x
   call dlaset('U', nP, nP, dZero, dZero, aux(1,1), nP)
   forall (i = 1:nP - 2) 
      aux(i, i + 2) =   sqrt(dble(i * (i + 1)))
      aux(i + 2, i) = - sqrt(dble(i * (i + 1)))
   end forall
   call dgemm('N','N', nP,nP,nP, dOne, aux,nP, trafo(1,1),LDD, dZero, aux1,nP)
   call dgemm('T','N', nP,nP,nP, dOne, trafo(1,1),LDD, aux1,nP, dZero, xd1Mat(1,1),LDD)
   ! Anti-Symmetrize
   do i = 1, nP
      do j = 1, i - 1
         dummy = dOneHalf * (d1Mat(i,j) - d1Mat(j,i))
         d1Mat(i,j) =   dummy
         d1Mat(j,i) = - dummy
         dummy = dOneHalf * (xd1Mat(i,j) - xd1Mat(j,i))
         xd1Mat(i,j) =   dummy
         xd1Mat(j,i) = - dummy
      end do
      d1Mat(i,i)  = dZero
      xd1Mat(i,i) = dZero
   end do
   !
   ! Calculate the matrix of second derivatives
   !
   call dlaset('U', nP, nP, dZero, dZero, aux(1,1), nP)
   forall (i = 1:nP) aux(i,i) = - dble(2 * i - 1)
   forall (i = 1:nP - 2) aux(i, i + 2) = sqrt(dble(i * (i + 1)))
   call dsymm('L', 'U', nP, nP, dOneHalf, aux(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1),nP)
   call dsyr2k('U', 'T', nP, nP, dOneHalf * mw, trafo(1,1), LDD, aux1(1,1), nP, dZero, d2Mat(1,1), LDD)
   !call dgemm('T', 'N', nP, nP, nP, mw, trafo(1,1), LDD, aux1(1,1), nP, dZero, d2Mat(1,1), LDD)
   ! Symmetrize
   do i = 1, nP
      do j = i + 1, nP
         d2Mat(j,i) = d2Mat(i,j)
         !dummy = dOneHalf * (d2Mat(i,j) + d2Mat(j,i))
         !d2Mat(i,j) = dummy
         !d2Mat(j,i) = dummy
      end do
   end do
   return
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! TEST
   forall (i = 1:nP, j = 1:nP) aux(i,j) = -dOneHalf * 0.001 * d2Mat(i,j)
   forall (i = 1:nP) aux(i,i) = aux(i,i) + dOneHalf * 0.001 * gridP(i) ** 2 - 0.002 * gridP(i) + 0.0001 * gridP(i)**4
   call dsyev('V', 'U', nP, aux, nP, energ, WORK, 3 * nP - 1, info)
   write(0,*)
   do i = 1, 10
      write(0,'(i3,3x,f9.6)') i, energ(i)
   enddo
   stop
   !
   return
   end subroutine harm_DVR



   subroutine Leg_DVR(nP,MM,gridP,wP,sinm2Mat,d1Mat,d2Mat,trafo,LDD)
   ! Legendre DVR with nP points
   implicit none
   integer :: i, info, j, k, idummy, jdummy, li, lj
   integer, intent(in) :: nP, LDD, MM
   double precision, dimension(*), intent(out) :: gridP, wP
   double precision, dimension(LDD,*), intent(out) :: sinm2Mat, d1Mat, d2Mat, trafo
   double precision :: dummy
   double precision, dimension(3 * nP - 1) :: WORK
   !
   integer, parameter :: n_XL = 3001
   double precision, dimension(n_XL) :: gridP_XL
   double precision, dimension(n_XL,n_XL) :: trafo_XL
   double precision, dimension(3 * n_XL - 1) :: WORK_XL
   double precision, dimension(nP,nP) :: sinm1Mat
   !
   double precision, dimension(nP,nP) :: aux, aux1
   double precision, parameter, dimension(11) :: &
      weight_factor = (/ 2.d0, 4.d0/3.d0, 16.d0/15.d0, 32.d0/35.d0, &
                        256.d0/315.d0, 512.d0/693.d0, 2048.d0/3003.d0, &
                        4096.d0/6435.d0, 65536.d0/109395.d0, &
                        131072.d0/230945.d0, 524288.d0/969969.d0 /)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   if (MM > 10) stop 'Leg DVR implemented only for M <= 10'
   if (nP .lt. 2) stop 'too few grid points selected for a DVR'
   ! Generate the grid and the weights
   ! aux is the matrix representation of the x = cos(theta) operator
   call dlaset('U', nP, nP, dZero, dZero, trafo(1,1), LDD)
   forall (i = 1:nP - 1) trafo(i,i + 1) = sqrt(dble((MM + i)**2 - MM**2) / dble(4 * (MM + i)**2 - 1))
   !forall (i = 1:nP - 1) trafo(i,i + 1) = sqrt(dble(i**2) / dble(4 * i**2 - 1))
   call dsyev('V', 'U', nP, trafo(1,1), LDD, gridP(1), WORK, 3 * nP - 1, info)
   forall (i = 1:nP) gridP(i) = acos(gridP(i))
   ! Symmetrize the grid about pi/2
   do i = 1, nP / 2
      dummy = dOneHalf * (gridP(i) + pi - gridP(nP - i + 1))
      gridP(i) = dummy
      gridP(nP - i + 1) =  pi - dummy
   end do
   if (mod(nP,2) == 1) gridP(i) = dOneHalf * pi
   ! Weights
      ! Make the square root of the weights positive
   forall (i = 1:nP, trafo(1,i) < dZero) trafo(:,i) = - trafo(:,i)
   forall (i = 1:nP) wP(i) = (trafo(1,i) / sin(gridP(i))**MM) ** 2 * weight_factor(MM + 1)
   !forall (i = 1:nP) wP(i) = trafo(1,i) ** 2 * 2
      ! note that trafo(i,1:nP) is the representation of the i-th eigenstate of L^2 in the DVR basis
   !
   ! Calculate the matrix of the operator d/dtheta sin(theta)   (antisymmetric)
   !
   call dlaset('U', nP, nP, dZero, dZero, aux(1,1), nP)
   forall (i = 1:nP - 1) 
      aux(i,i + 1) =   sqrt(dble((MM+i)**4 - (MM+i)**2*MM**2) / dble(4 * (MM + i)**2 - 1))
      aux(i + 1,i) = - sqrt(dble((MM+i)**4 - (MM+i)**2*MM**2) / dble(4 * (MM + i)**2 - 1))
   end forall
   call dgemm('N', 'N', nP, nP, nP, dOne, aux(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1), nP)
   call dgemm('T', 'N', nP, nP, nP, dOne, trafo(1,1), LDD, aux1(1,1), nP, dZero, d1Mat(1,1), LDD)
   ! Enforce anti-symmetry
   do i = 1, nP
      do j = 1, i - 1
         dummy = dOneHalf * (d1Mat(i,j) - d1Mat(j,i))
         d1Mat(i,j) =   dummy
         d1Mat(j,i) = - dummy
      end do
      d1Mat(i,i)  = dZero
   end do
   !
   ! Calculate the matrix of the operator 1/sin^2
   !
   if (MM == 0) then   ! approximate 1/sin^2 as (1/sin)^2  [see JCP 154, 104115 (2021)]
      ! Set up a large grid for the evaluation of 1/sin(theta)
      call dlaset('U', n_XL, n_XL, dZero, dZero, trafo_XL(1,1), n_XL)
      forall (i = 1:n_XL - 1) trafo_XL(i,i + 1) = dble(i) / sqrt(dble(4 * i**2 - 1))
      call dsyev('V', 'U', n_XL, trafo_XL(1,1), n_XL, gridP_XL(1), WORK_XL, 3 * n_XL - 1, info)
      forall (i = 1:n_XL) gridP_XL(i) = acos(gridP_XL(i))
      ! Symmetrize the grid
      do i = 1, n_XL / 2
         dummy = dOneHalf * (gridP_XL(i) + pi - gridP_XL(n_XL - i + 1))
         gridP_XL(i) = dummy
         gridP_XL(n_XL - i + 1) =  pi - dummy
      end do
      gridP_XL(i) = dOneHalf * pi
      ! Set up the sin^-1 matrix in the basis of the nP lowest Legendre polynomials
      call dlaset('U', nP, nP, dZero, dZero, sinm1Mat(1,1), nP)   !! the generic aux matrix could play the role of sinm1Mat...
      do i = 1, n_XL
         call dsyr('U',nP, dOne / sin(gridP_XL(i)), trafo_XL(1,i),1, sinm1Mat, nP)
      end do
      ! Transform back to the DVR basis
      call dsymm('L', 'U', nP, nP, dOne, sinm1Mat(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1),nP)
      call dsyr2k('U', 'T', nP, nP, dOneHalf, trafo(1,1), LDD, aux1(1,1), nP, dZero, sinm1Mat(1,1), nP)
      forall (i = 1:nP, j = 1:nP, i > j) sinm1Mat(i,j) = sinm1Mat(j,i)
      ! Square the sin^-1 matrix
      write(0,*)
      do i = 1, 8
         write(0,'(8(2x,f8.2))') sinm1Mat(i,1:8)
      end do
      call dsyrk('U','N',nP,nP,dOne,sinm1Mat,nP,dZero,sinm2Mat,LDD)
      forall (i = 1:nP, j = 1:nP, i > j) sinm2Mat(i,j) = sinm2Mat(j,i)
      write(0,*)
      do i = 1, 8
         write(0,'(8(2x,f8.2))') sinm2Mat(i,1:8)
      end do
      ! Transform back to the DVR basis
!!!      call dsymm('L', 'U', nP, nP, dOne, sinm2Mat(1,1), LDD, trafo(1,1), LDD, dZero, aux1(1,1),nP)
!!!      call dsyr2k('U', 'T', nP, nP, dOneHalf, trafo(1,1), LDD, aux1(1,1), nP, dZero, sinm2Mat(1,1), LDD)
!!!      forall (i = 1:nP, j = 1:nP, i > j) sinm2Mat(i,j) = sinm2Mat(j,i)
!!!      write(0,*)
!!!      do i = 1, 8
!!!         write(0,'(8(2x,f8.2))') sinm2Mat(i,1:8)
!!!      end do
      !! Check (diagonalize)
      !call dsyev('V','U',nP,sinm2Mat,LDD,gridP,WORK,3*nP-1,info)
      !write(0,*)
      !do i = 1, nP
      !   write(0,*) i, gridP(i)
      !end do
      !stop
   !
   else  ! MM > 0
      call dlaset('U', nP, nP, dZero, dZero, aux(1,1), nP)
      do i = 1, nP
         li = i - 1 + MM
         idummy = 1
         do k = 0, 2 * MM - 1
            idummy = idummy * (li + MM - k)
         end do
         do j = i, nP, 2
            lj = j - 1 + MM
            jdummy = 1
            do k = 0, 2 * MM - 1
               jdummy = jdummy * (lj + MM - k)
            end do
            aux(i,j) = sqrt(dble((2 * li + 1) * (2 * lj + 1))) / dble(2 * MM) &
                     * sqrt(dble(idummy) / dble(jdummy))
         end do
      end do
      call dsymm('L', 'U', nP, nP, dOne, aux(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1),nP)
      call dsyr2k('U', 'T', nP, nP, dOneHalf, trafo(1,1), LDD, aux1(1,1), nP, dZero, sinm2Mat(1,1), LDD)
      do j = 1, nP
         do i = j + 1, nP
            sinm2Mat(i,j) = sinm2Mat(j,i)
         end do
      end do
   end if
!   write(0,*)
!   do i = 1, 8
!      write(0,'(8(2x,f8.2))') aux(i,1:8)
!   end do
!   write(0,*)
!   do i = 1, 8
!      write(0,'(8(2x,f8.2))') sinM2Mat(i,1:8)
!   end do
!   stop
   !
   ! Calculate the matrix of the operator L0^2 = -1/sin(t) * d/dt sin(t) d/dt
   !
   call dlaset('U', nP, nP, dZero, dZero, aux(1,1), nP)
   forall (i = 1:nP) aux(i,i) = dble((MM+i)**2 - (MM+i))
   call dsymm('L', 'U', nP, nP, dOne, aux(1,1), nP, trafo(1,1), LDD, dZero, aux1(1,1),nP)
   call dlacpy('U',nP,nP,sinm2Mat,LDD,d2Mat,LDD)
   call dsyr2k('U', 'T', nP, nP, dOneHalf, trafo(1,1), LDD, aux1(1,1), nP, -dble(MM)**2, d2Mat(1,1), LDD)
   ! Symmetrize
   do i = 1, nP
      do j = i + 1, nP
         d2Mat(j,i) = d2Mat(i,j)
      end do
   end do
   ! Remove M^2/sin^2
  !! do i = 1, nP
  !!    do j = 1, nP
  !!       d2Mat(i,j) = d2Mat(i,j) - dble(MM)**2 * sinm2Mat(i,j)
  !!    end do
  !! end do
   !do i = 1, nP
   !   d2Mat(i,i) = d2Mat(i,i) - (dble(MM) / sin(gridP(i)))**2
   !end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine Leg_DVR


 
   subroutine Eigenf1DPot(nn, grid, facT, D2, LDD, PotFile, iEig, phi)
   ! Calculate one of the eigenstates of the operator
   !   H = - facT * d^2/dx^2  + V(x)
   ! 
   ! nn      : number of grid points
   ! grid    : DVR grid points
   ! facT    : kinetic energy prefactor
   ! D2      : DVR matrix of the second derivative operator
   ! LDD     : leading dimension of D2
   ! PotFile : file containing the potential, to be interpolated to the DVR grid
   ! iEig    : desired eigenstate (iEig = 0) is the ground state
   ! phi     : eigenstate, as output
   implicit none
   integer, intent(in) :: nn, LDD, iEig
   double precision, intent(in) :: facT
   double precision, dimension(nn), intent(in) :: grid
   double precision, dimension(LDD,*), intent(in) :: D2
   character(20) :: Potfile
   double complex, dimension(nn), intent(out) :: phi
   !
   integer :: nLines, idummy, iT, j, k, info
   double precision, dimension(nn) :: gridPot, W
   double precision, dimension(:), allocatable :: xdummy, ydummy
   double precision, dimension(nn,nn) :: UM
   double precision, dimension(3 * nn) :: WORK
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
! Read the potential
   iT = freeunit()
   open(unit = iT, file = trim(adjustl(PotFile)), form = 'formatted', status = 'old', action = 'read')
   ! Get the number of lines
   nLines = 0
   idummy = 0
   do
      read(iT,*, iostat = idummy)
      if (idummy .ne. 0) exit
      nLines = nLines + 1
   end do
   ! Read the data
   allocate(xdummy(nLines), ydummy(nLines))
   rewind(iT)
   do j = 1, nLines
      read(iT,*) xdummy(j), ydummy(j)
   end do
   ! Interpolation
   call Spline(nLines, xdummy, ydummy, nn, grid, gridPot)
   deallocate(xdummy,ydummy)
   close(iT)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
! Diagonalize
   forall (j = 1:nn, k = 1:nn) UM(j,k) = - facT * D2(j,k)
   forall (j = 1:nn) UM(j,j) = UM(j,j) + gridPot(j)
   call dsyev('V', 'U', nn, UM, nn, W, WORK, 3 * nn, info)
! Set up the wavefunction
   forall (j = 1:nn) phi(j) = complex(UM(j,iEig + 1), dZero)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine Eigenf1DPot


   subroutine Eigenf1DPot_real(nn, grid, facT, D2, LDD, PotFile, iEig, phi)
   ! Calculate one of the eigenstates of the operator
   !   H = - facT * d^2/dx^2  + V(x)
   ! 
   ! nn      : number of grid points
   ! grid    : DVR grid points
   ! facT    : kinetic energy prefactor
   ! D2      : DVR matrix of the second derivative operator
   ! LDD     : leading dimension of D2
   ! PotFile : file containing the potential, to be interpolated to the DVR grid
   ! iEig    : desired eigenstate (iEig = 0) is the ground state
   ! phi     : eigenstate, as output
   implicit none
   integer, intent(in) :: nn, LDD, iEig
   double precision, intent(in) :: facT
   double precision, dimension(nn), intent(in) :: grid
   double precision, dimension(LDD,*), intent(in) :: D2
   character(20) :: Potfile
   double precision, dimension(nn), intent(out) :: phi
   !
   integer :: nLines, idummy, iT, j, k, info
   double precision, dimension(nn) :: gridPot, W
   double precision, dimension(:), allocatable :: xdummy, ydummy
   double precision, dimension(nn,nn) :: UM
   double precision, dimension(3 * nn) :: WORK
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
! Read the potential
   iT = freeunit()
   open(unit = iT, file = trim(adjustl(PotFile)), form = 'formatted', status = 'old', action = 'read')
   ! Get the number of lines
   nLines = 0
   idummy = 0
   do
      read(iT,*, iostat = idummy)
      if (idummy .ne. 0) exit
      nLines = nLines + 1
   end do
   ! Read the data
   allocate(xdummy(nLines), ydummy(nLines))
   rewind(iT)
   do j = 1, nLines
      read(iT,*) xdummy(j), ydummy(j)
   end do
   ! Interpolation
   call Spline(nLines, xdummy, ydummy, nn, grid, gridPot)
   deallocate(xdummy,ydummy)
   close(iT)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
! Diagonalize
   forall (j = 1:nn, k = 1:nn) UM(j,k) = - facT * D2(j,k)
   forall (j = 1:nn) UM(j,j) = UM(j,j) + gridPot(j)
   call dsyev('V', 'U', nn, UM, nn, W, WORK, 3 * nn, info)
! Set up the wavefunction
   forall (j = 1:nn) phi(j) = UM(j,iEig + 1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   
   return
   end subroutine Eigenf1DPot_real


end module dvr
