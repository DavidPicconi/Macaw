!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mathematics/Algebra subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Standard regularized matrix inversion

subroutine StdRegInv(nn,AM,LDA,BM,LDB,eps,MatName,RegType)
! Regularized inversion of a matrix A
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix to be inverted
!    LDA : leading dimension of AM, as declared in the calling program
!    BM  : (output) inverse of the AM matrix
!    LDB : leading dimension of BM, as declared in the calling program
!    eps : (input) regularization paramenter
! MatName: (input) name of the matrix (used in the error message)
! RegType: (input) Type of regularization
integer, intent(in) :: nn, LDA, LDB
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:,:), intent(out) :: BM(LDB,*)
double precision, intent(in) :: eps
character(4), intent(in) :: MatName, RegType
integer :: info
integer :: j
double precision :: Wv2e2, expW
double precision, dimension(:) :: Wv(nn), RWORK(3 * nn - 2)
double complex, dimension(:) :: WORK(2 * nn - 1)
double complex, dimension(:,:) :: UM(nn,nn) 
double precision, parameter :: dOne = 1.d0, dOneHalf = 0.5d0, StdRegThresh = 64.d0
double complex, parameter :: cZero = (0.d0,0.d0)
double precision, external :: dasum
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
call zlacpy('U',nn,nn,AM,LDA,UM,nn)
! Diagonalisation (necessary to identify zero eigenvalues)
info = 0
call zheev('V','U',nn,UM,nn,Wv,WORK,2 * nn - 1,RWORK,info)
if (info .ne. 0) then
   write(0,*) 'Error in subroutine StdRegInv'
   write(0,*) 'for the regularization of the matrix ', MatName
   write(0,*) '   Lapack subroutine zheev returns a value of info =', info
   write(0,*)
   do j = 1, nn
      write(0,*) j, AM(j,j:nn)
   end do
   stop
end if
! Regularisation
select case (RegType)
   case('tikh')
      Wv = Wv / (Wv ** 2 + eps ** 2)  ! Tikhonov
   case('tik2') ! 2-iterated Tikhonov
      do j = 1, nn
         Wv2e2 = (Wv(j) ** 2 + eps ** 2)
         Wv(j) = Wv(j) * (Wv2e2 + eps ** 2) / Wv2e2 ** 2
      end do
   case('poly')
      where (Wv .lt. eps)
         Wv = dOneHalf * (eps + Wv ** 2 / eps)
      end where
      Wv = dOne / Wv
   case('expo')
      do j = 1, nn
         expW = eps * exp(- (Wv(j) / eps) ** 2)
         Wv(j) = Wv(j) / (Wv(j) ** 2 + expW ** 2)
      end do
   case('strd')  ! standard MCTDH 
      where (Wv .lt. eps * StdRegThresh)
         Wv = Wv + eps * exp(-Wv / eps)
      end where
      Wv = dOne / Wv
   case default
      write(0,*) 'Error in subroutine StdRegInv:'
      write(0,*) ' no regularization type specified for the matrix ', MatName
      write(0,*) RegType
      stop
end select
call zlaset('U',nn,nn,cZero,cZero,BM,LDB)
do j = 1, nn
   call zher('U',nn,Wv(j),UM(1,j),1,BM,LDB)
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return  
end subroutine StdRegInv



subroutine StdRegSol(nn,AM,LDA,YV,XV,eps,RegType)
! Solution of the linear system A * X = Y via direct regularized inverse of the matrix AM
! Regularized inversion of a matrix A
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
!    eps : (input) regularization paramenter
! RegType: (input) Type of regularization
use sysparam
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
double precision, intent(in) :: eps
character(4), intent(in) :: RegType
integer :: info
integer :: j
double precision :: Wv2e2, expW
double precision, dimension(:) :: Wv(nn), RWORK(3 * nn - 2)
double complex, dimension(:) :: WORK(2 * nn - 1)
double complex, dimension(:,:) :: UM(nn,nn)
double complex, dimension(:) :: UYV(nn)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
call zlacpy('U',nn,nn,AM,LDA,UM,nn)
! Diagonalisation (necessary to identify zero eigenvalues)
info = 0
call zheev('V','U',nn,UM,nn,Wv,WORK,2 * nn - 1,RWORK,info)
! temp: print conditioning number
!eps2 = max(eps**2, Wv(1) * Wv(nn) * (Wv(nn) - 1000.0 * Wv(1)) / (1000.0 * Wv(nn) - Wv(1)))
!write(1989,*) Wv(1), Wv(nn), Wv(nn) / Wv(1), &
!        (Wv(nn)**2 + eps2) / (Wv(1)**2 + eps2) * Wv(1) / Wv(nn), sqrt(eps2)
!
if (info .ne. 0) then
   write(0,*) nn
   write(0,*) 'Error in subroutine StdRegSol'
   write(0,*) '   Lapack subroutine zheev returns a value of info =', info
   write(0,*) AM(1,1), AM(1,2)
   write(0,*) AM(2,1), AM(2,2)
   stop
end if
! Transform the vector YV -> UM^H YV   
call zgemv('C',nn,nn,cOne,UM,nn,YV,1,cZero,UYV,1)
! Solution with Regularisation
select case (RegType)
   case('tikh')
      UYV = UYV * Wv / (Wv ** 2 + eps ** 2)
   case('tik2')
      do j = 1, nn
         Wv2e2 = Wv(j) ** 2 + eps ** 2
         UYV(j) = UYV(j) * Wv(j) * (Wv2e2 + eps ** 2) / Wv2e2 ** 2
      end do
   case('poly')
      where (Wv .lt. eps)
         UYV = UYV * 2 * eps / (eps ** 2 + Wv ** 2)
      elsewhere
         UYV = UYV / Wv
      end where
   case('expo')
      do j = 1, nn
         expW = eps * exp(- (Wv(j) / eps) ** 2)
         UYV(j) = UYV(j) * Wv(j) / (Wv(j) ** 2 + expW ** 2)
      end do
   case default
      write(0,*) 'Error in subroutine StdRegSol:'
      write(0,*) ' no regularization type specified '
      write(0,*) RegType
      stop
end select
!
! Back-transform and write the result on the solution vector
call zgemv('N',nn,nn,cOne,UM,nn,UYV,1,cZero,XV,1)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return  
end subroutine StdRegSol



subroutine DirectSol(nn,AM,LDA,YV,XV)
! Solution of the linear system A * X = Y via direct inversion the matrix AM
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(inout) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
integer :: info
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
call zcopy(nn,YV(1),1,XV(1),1)
call zposv('U',nn,1,AM(1,1),LDA,XV(1),nn,info)
if (info .gt. 0) then
   write(0,*) 'info = ', info
   write(0,*) 'Factorization not completed'
   !stop
end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DirectSol



subroutine DirectLSSol(nn,AM,LDA,YV,XV)
! Least-square solution of the linear system A * X = Y
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!           (only upper diagonal is set up)
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(inout) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
integer :: info, i, j
double complex, dimension(:) :: WORK(2 * nn)
double precision, parameter :: eps = 1.d-14
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
call zcopy(nn,YV(1),1,XV(1),1)
! Set up the lower diagonal of AM
forall(i = 1:nn, j = 1:nn, i .gt. j) AM(i,j) = conjg(AM(j,i))
! Add a tiny positive terms on the diagonal
forall(i = 1:nn) AM(i,i) = AM(i,i) + eps
call zgels('N',nn,nn,1,AM(1,1),LDA,XV(1),nn,WORK,2 * nn,info)
if (info .gt. 0) then
   write(0,*) 'info = ', info
   write(0,*) 'The matrix does not have the full rank,'
   write(0,*) 'therefore it could not be factorized'
   stop
end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DirectLSSol



subroutine DirectExtrapSol(nn,AM,LDA,YV,XV,eps)
! Solution of the linear system A * X = Y via direct inversion the matrix AM
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
integer, intent(in) :: nn, LDA
double precision, intent(in) :: eps
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
double complex, dimension(:,:) :: AMeps(nn,nn), Xeps(nn,4)
integer :: info, n, i
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
do n = 0, 3
   call zlacpy('U',nn,nn,AM,LDA,AMeps,nn)
   forall (i = 1:nn) AMeps(i,i) = AMeps(i,i) + eps * 2 ** n
   call zcopy(nn,YV(1),1,Xeps(1,n + 1),1)
   call zposv('U',nn,1,AMeps(1,1),LDA,Xeps(1,n + 1),nn,info)
   if (info .gt. 0) then
      write(0,*) 'DirectExtrapSol:'
      write(0,*) 'n = ', n
      write(0,*) 'info = ', info
      write(0,*) 'Factorization not completed'
      stop
   end if
end do
! Extrapolation
forall (i = 1:nn) &
   XV(i) = (64 * Xeps(i,1) - 56 * Xeps(i,2) + 14 * Xeps(i,3) - Xeps(i,4)) / 21.d0
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DirectExtrapSol



subroutine DirectTikhInv(nn,AM,LDA,BM,LDB,eps)
! Matrix inversion of a Hermitian matrix AM with Tikhonov regularization
!    nn  : dimension of AM
!    AM  : (input) nn x nn (full!!) matrix
!    LDA : (input) leading dimension of AM, as declared in the calling program
!    BM  : (output) regularized inverse of AM (upper diagonal)
!    LDB : (input) leading dimension of BM, as declared in the calling program
!    eps : (input) regularization paramenter
use sysparam
integer, intent(in) :: nn, LDA, LDB
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:,:), intent(out) :: BM(LDB,*)
double precision, intent(in) :: eps
integer :: info
integer :: i, j
double complex, dimension(:,:) :: AM2(nn,nn)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
!call zlacpy('all',nn,nn,AM,LDA,BM,LDB)
! Set up the lower part of BM
!do concurrent (i = 1:nn, j = 1:nn, i .gt. j) 
!   BM(i,j) = conjg(AM(j,i))
!end do
! Square of AM
call zlaset('U',nn,nn,cZero,complex(eps**2,dZero),AM2(1,1),nn)
call zherk('U','N',nn,nn,cOne,AM(1,1),LDA,cOne,AM2(1,1),nn)
!call zherk('U','N',nn,nn,cOne,BM(1,1),LDB,cZero,AM2(1,1),nn)
   ! Tikhonov regularization
!do concurrent (i = 1:nn) 
!   AM2(i,i) = AM2(i,i) + eps ** 2
!end do
!call zposv('U',nn,nn,AM2(1,1),nn,BM(1,1),LDB,info)
call zpotrf('U',nn,AM2(1,1),nn,info)   ! Cholesky factorization
call zpotri('U',nn,AM2(1,1),nn,info)   ! inverse (AM**2 + eps**2)^-1
! Symmetric rank-2 construction of AM^-1
do concurrent (i = 1:nn, j = 1:nn, i .gt. j) 
   AM2(i,j) = conjg(AM2(j,i))
end do
call zher2k('U','N',nn,nn,cOneHalf,AM2(1,1),nn,AM(1,1),LDA,dZero,BM(1,1),LDB)
! Symmetrize the inverse matrix
!BM(1:nn,1:nn) = 0.5d0 * (BM(1:nn,1:nn) + conjg(transpose(BM(1:nn,1:nn))))
!if (info .gt. 0) then
!   write(0,*) 'info = ', info
!   write(0,*) 'Factorization not completed'
!   stop
!end if
!write(0,*)
!do i = 1, 4
!   write(0,'(4(3x, "(", f8.3, ",", f8.3, ")"))') AM(i,1:4) 
!enddo
!write(0,*)
!do i = 1, 4
!   write(0,'(4(3x, "(", f8.3, ",", f8.3, ")"))') BM(i,1:4) 
!enddo
!stop
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DirectTikhInv



subroutine DirectTikhSol(nn,AM,LDA,YV,XV,eps)
! Solution of the linear system A * X = Y via direct inversion the matrix AM
! with Tikhonov regularization, i.e. solve (A*A + eps**2) * X = A * Y 
! Regularized inversion of a matrix A
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
!    eps : (input) regularization paramenter
use sysparam
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(inout) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
double precision, intent(in) :: eps
integer :: info
integer :: i, j
double complex, dimension(:,:) :: AM2(nn,nn)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Set up the lower part of AM
do j = 1, nn
   do i = j + 1, nn
      AM(i,j) = conjg(AM(j,i))
   end do
end do
! Square of AM
call zlaset('U',nn,nn,cZero,cOne,AM2(1,1),nn)
call zherk('U','N',nn,nn,cOne,AM(1,1),LDA,eps**2,AM2(1,1),nn)
   ! Tikhonov regularization
!forall (i = 1:nn) AM2(i,i) = AM2(i,i) + eps ** 2
call zhemv('U',nn,cOne,AM(1,1),LDA,YV(1),1,cZero,XV(1),1)
call zposv('U',nn,1,AM2(1,1),nn,XV(1),nn,info)
if (info .gt. 0) then
   write(0,*) 'info = ', info
   write(0,*) 'Factorization not completed'
   stop
end if
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine DirectTikhSol


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conjugate gradient

subroutine CGSol(nn,AM,LDA,YV,XV)
! Solution of the linear system A * X = Y via conjugate gradient iterations
! Regularized inversion of a matrix A
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
use sysparam
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
integer, parameter :: maxit = 20000
integer :: i, iit
integer, save :: counter = 0
integer, save, dimension(maxit + 1) :: itHist = 0
double precision, parameter :: tol = 1.d-8
double precision :: delta, alpha, beta, delta_old
double complex, dimension(nn) :: r0, p0, Ap0
double precision, external :: dznrm2, dzasum
double complex, external :: zdotc
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
counter = counter + 1
!!! Initialize the X vector (divide Y by the diagonal terms of the AM matrix - with regularization)
!forall (i = 1:nn) XV(i) = cOne*0.1d0 ! YV(i) * AM(i,i) / (AM(i,i) ** 2 + eps ** 2)
!call zhemv('U',nn,-cOne,AM(1,1),LDA,XV(1),1,cZero,r0(1),1)
!call zaxpy(nn,cOne,YV(1),1,r0(1),1)
!delta = dznrm2(nn,r0(1),1)
!!! Initialize the X vector to zero
forall (i = 1:nn) 
   XV(i) = cZero
   r0(i) = YV(i)
end forall
delta = real(zdotc(nn,r0,1,r0,1))
!!!
if (maxval(abs(r0)) < tol) then 
!if (sqrt(delta) .lt. tol) then 
   !write(3002,'(i6,3x,a,i3,3x,a,f9.6)') counter, 'Iterations: ', 0, 'Maximum rel. error: ', maxval(abs(r0)) / tol
   !write(3002,'(a,i3,3x,a,e10.4)') 'Iterations: ', 0, 'Residual: ', sqrt(delta)
   itHist(1) = itHist(1) + 1
   call WriteHist(counter,itHist)
   return
end if
p0 = r0
iit = 0
do
   iit = iit + 1
   if (iit .gt. maxit) then
      write(0,*)  'Conjugate gradient: Maximum number of iterations reached'
      write(0,*) '  delta = ', delta
      write(0,*) '  Maximum rel. error = ', maxval(abs(r0)) / tol
      stop
   end if
   call zhemv('U',nn,cOne,AM(1,1),LDA,p0(1),1,cZero,Ap0(1),1)
   alpha = delta / real(zdotc(nn,p0(1),1,Ap0(1),1))
   call zaxpy(nn, complex(alpha,dZero), p0(1),1,XV(1),1)
   call zaxpy(nn,-complex(alpha,dZero),Ap0(1),1,r0(1),1)
   !r0 = r0 - alpha * Ap0
   delta_old = delta
   if (maxval(abs(r0)) < tol) then 
      !write(3002,'(i6,3x,a,i3,3x,a,f9.6)') counter, 'Iterations: ', iit, 'Maximum rel. error: ', maxval(abs(r0)) / tol
      itHist(iit + 1) = itHist(iit + 1) + 1
      call WriteHist(counter,itHist)
      return
   end if
   !delta = real(zdotc(nn,r0,1,r0,1))
   !if (sqrt(delta) .lt. tol) return
   delta = real(zdotc(nn,r0,1,r0,1))
   beta = delta / delta_old
   p0 = r0 + beta * p0
end do
!write(3002,'(a,i3,3x,a,e10.4)') 'Iterations: ', iit, 'Residual: ', sqrt(delta)
!flush(3002)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine CGSol


subroutine CGPrecSol(nn,AM,LDA,YV,XV)
! Solution of the linear system A * X = Y via conjugate gradient iterations
! Regularized inversion of a matrix A
!    nn  : dimension of AM
!    AM  : (input) Hermitian nn x nn positive definite matrix at the l.h.s
!    LDA : leading dimension of AM, as declared in the calling program
!    YV  : (input) vector at the r.h.s. of the equation
!    XV  : (output) solution vector
!    eps : (input) regularization parameter
use sysparam
integer, intent(in) :: nn, LDA
double complex, dimension(:,:), intent(in) :: AM(LDA,*)
double complex, dimension(:), intent(in) :: YV(*)
double complex, dimension(:), intent(out) :: XV(*)
integer, parameter :: maxit = 20000
integer :: i, iit
integer, save :: counter = 0
integer, save, dimension(maxit + 1) :: itHist = 0
double precision, parameter :: tol = 1.d-4
double precision :: delta, alpha, beta, delta_old
double complex, dimension(nn) :: r0, p0, Ap0
double precision, external :: dznrm2, dzasum
double complex, external :: zdotc
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
counter = counter + 1
!!! Initialize the X vector 
!call zhemv('U',nn,cOne,Prec(1,1),LDP,YV(1),1,cZero,r0(1),1)
!call zhemv('U',nn,cOne,Prec(1,1),LDP,r0(1),1,cZero,XV(1),1)
!
!call zlacpy('U',nn,nn,AM(1,1),LDA,aux,nn)
!forall (i = 1:nn) 
!   aux(i,i) = aux(i,i) + 1d-6
!   XV(i) = YV(i)
!end forall
!call zposv('U',nn,1,aux,nn,XV(1),nn,info)
!call zhemv('U',nn,-cOne,AM(1,1),LDA,XV(1),1,cZero,r0(1),1)
!call zaxpy(nn,cOne,YV(1),1,r0(1),1)
!
forall (i = 1:nn) 
   XV(i) = cZero
   r0(i) = YV(i)
end forall
delta = real(zdotc(nn,r0,1,r0,1))
!!!
!if (maxval(abs(r0)) < tol) then 
if (sqrt(delta) < tol) then 
   !write(3002,'(i6,3x,a,i3,3x,a,f9.6)') counter, 'Iterations: ', 0, 'Maximum rel. error: ', maxval(abs(r0)) / tol
   !write(3002,'(a,i3,3x,a,e10.4)') 'Iterations: ', 0, 'Residual: ', sqrt(delta)
   itHist(1) = itHist(1) + 1
   call WriteHist(counter,itHist)
   return
end if
p0 = r0
iit = 0
do
   iit = iit + 1
   if (iit .gt. maxit) then
      write(0,*)  'Conjugate gradient: Maximum number of iterations reached'
      write(0,*) '  delta = ', delta
      write(0,*) '  Maximum rel. error = ', maxval(abs(r0)) / tol
      stop
   end if
   call zhemv('U',nn,cOne,AM(1,1),LDA,p0(1),1,cZero,Ap0(1),1)
   alpha = delta / real(zdotc(nn,p0(1),1,Ap0(1),1))
   XV(1:nn) = XV(1:nn) + alpha * p0
   !call zaxpy(nn,cZero + alpha,p0(1) ,1,XV(1),1)
   !call zaxpy(nn,cZero - alpha,Ap0(1),1,r0(1),1)
   r0 = r0 - alpha * Ap0
   delta_old = delta
   delta = real(zdotc(nn,r0,1,r0,1))
   if (sqrt(delta) < tol) then
   !if (maxval(abs(r0)) < tol) then 
      !write(3002,'(i6,3x,a,i3,3x,a,f9.6)') counter, 'Iterations: ', iit, 'Maximum rel. error: ', maxval(abs(r0)) / tol
      itHist(iit + 1) = itHist(iit + 1) + 1
      call WriteHist(counter,itHist)
      return
   end if
   !if (sqrt(delta) .lt. tol) return
   beta = delta / delta_old
   p0 = r0 + beta * p0
end do
!write(3002,'(a,i3,3x,a,e10.4)') 'Iterations: ', iit, 'Residual: ', sqrt(delta)
!flush(3002)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine CGPrecSol


subroutine WriteHist(counter,itHist)
integer, dimension(80), intent(in) :: itHist
integer, intent(in) :: counter
integer :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit = 5151, file = 'CG.hist')
write(5151,*) '# Distribution after', counter, ' CG applications'
do i = 1, 80
   write(5151,*) i - 1, itHist(i), dble(itHist(i)) / real(counter)
end do
close(5151)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
return
end subroutine WriteHist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! (Inverse) square root of a matrix

subroutine SquareRoot(n, AA, LDAA, Am1, LDA)
! Calculate the square root of the Hermitian matrix AA
! i.e.  Am1*Am1 = AA
use sysparam
integer, intent(in) :: n, LDAA, LDA   ! size and leading dimensions of AA and Am1
double complex, dimension(LDAA,*), intent(in) :: AA
double complex, dimension(LDA ,*), intent(out) :: Am1
integer :: info, j, j1, j2
double precision, dimension(n) :: Wv
double precision, dimension(3 * n - 2) :: RWORK
double complex, dimension(2 * n - 1) :: WORK
double complex, dimension(n,n) :: aux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call zlacpy('U', n, n, AA, LDAA, aux, n)
info = 0
call zheev('V', 'U', n, aux, n, Wv, WORK, 2 * n - 1, RWORK, info)
call zlaset('U', n, n, cZero, cZero, Am1, LDA)
do j = 1, n
   call zher('U', n, sqrt(Wv(j)), aux(1,j), 1, Am1, LDA)
end do
! Fill the lower part
do j2 = 1, n
   do j1 = j2 + 1, n
      Am1(j1,j2) = conjg(Am1(j2,j1))
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
return
end subroutine SquareRoot



subroutine InverseSquareRoot(n, AA, LDAA, Am1, LDA, eps)
! Calculate the inverse square root of the Hermitian matrix AA
! i.e.  Am1*Am1 = AA^-1
use sysparam
integer, intent(in) :: n, LDAA, LDA   ! size and leading dimensions of AA and Am1
double complex, dimension(LDAA,*), intent(in) :: AA
double complex, dimension(LDA ,*), intent(out) :: Am1
double precision, intent(in) :: eps
integer :: info, j, j1, j2
double precision, dimension(n) :: Wv
double precision, dimension(3 * n - 2) :: RWORK
double complex, dimension(2 * n - 1) :: WORK
double complex, dimension(n,n) :: aux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call zlacpy('U', n, n, AA, LDAA, aux, n)
info = 0
call zheev('V', 'U', n, aux, n, Wv, WORK, 2 * n - 1, RWORK, info)
call zlaset('U', n, n, cZero, cZero, Am1, LDA)
do j = 1, n
   call zher('U', n, sqrt(Wv(j) / (Wv(j)**2 + eps**2)), aux(1,j), 1, Am1, LDA)
   !call zher('U', n, dOne / sqrt(abs(Wv(j))), aux(1,j), 1, Am1, LDA)
end do
! Fill the lower part
do j2 = 1, n
   do j1 = j2 + 1, n
      Am1(j1,j2) = conjg(Am1(j2,j1))
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
return
end subroutine InverseSquareRoot


subroutine InverseSquareRootU(n, AA, LDAA, Am1, LDA)
! Calculate a (improper) inverse square root of the Hermitian matrix AA
! i.e.  Am1^H * AA * Am1 = 1 and Am1 * Am1^H = AA^-1
! using Cholesky decomposition
use sysparam
integer, intent(in) :: n, LDAA, LDA   ! size and leading dimensions of AA and Am1
double complex, dimension(LDAA,*), intent(in) :: AA
double complex, dimension(LDA ,*), intent(out) :: Am1
integer :: info, j1, j2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call zlacpy('U', n, n, AA, LDAA, Am1, LDA)
! Cholesky decomposition
call zpotrf('U', n, Am1, LDA, info)
! Inversion
call ztrtri('U', 'N', n, Am1, LDA, info)
! Fill the lower diagonal
do j2 = 1, n
   do j1 = j2 + 1, n
      Am1(j1,j2) = cZero
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
return
end subroutine InverseSquareRootU




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cubic spline
subroutine Spline(n, x, y, nn, xx, yy)
! Construct a cubic spline interpolation of a function defined on a grid
! onto another grid. 
! Quadratic extrapolation is used for the points outside the input grid.
!
!    n  : number of input grid points
!    x  : input grid
!    y  : values of the function on the input grid
!    nn : number of output grid points
!    xx : output grid
!    yy : values of the function on the output grid
use sysparam
integer, intent(in) :: n, nn
double precision, dimension(:), intent(in) :: x(n), y(n), xx(nn)
double precision, dimension(:), intent(out) :: yy(nn)
integer :: i, info, iEq, iABC, j, iC
integer, dimension(:) :: ipiv(3 * n - 3)
double precision :: dx, xt
double precision, dimension(:) :: Coeff(3 * n - 3)
double precision, dimension(:,:) :: SplineMat(3 * n - 3, 3 * n - 3)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! In each interval [x_(i-1), x_i], the function is approximated as 
! f_i(x) = y_i + a_i * (x - x_i) + b_i * (x - x_i)**2 + c_i * (x - x_i)**3
! The parameters are arranged as Coeff = (a1,b1,c1,a2,b2,c2,...,aN-1,bN-1,cN-1)
! and c1 = cN-1 = 0
SplineMat = dZero
Coeff = dZero
! Boundary: zero third derivative in the first interval (c3 = 0)
SplineMat(1,3) = dOne
iEq = 1
iABC = 0
do i = 1, n - 2
   dx = x(i + 1) - x(i)
   ! Eqs. (1): continuity of the function
   SplineMat(iEq + 1,iABC + 1) = dx
   SplineMat(iEq + 1,iABC + 2) = dx ** 2
   SplineMat(iEq + 1,iABC + 3) = dx ** 3
   Coeff(iEq + 1) = y(i + 1) - y(i)
   ! Eqs. (2): continuity of the first derivative
   SplineMat(iEq + 2,iABC + 1) = dOne
   SplineMat(iEq + 2,iABC + 2) = 2 * dx
   SplineMat(iEq + 2,iABC + 3) = 3 * dx ** 2
   SplineMat(iEq + 2,iABC + 4) = - dOne
   ! Eqs. (3): continuity of the second derivative
   SplineMat(iEq + 3,iABC + 2) = dOne
   SplineMat(iEq + 3,iABC + 3) = 3 * dx
   SplineMat(iEq + 3,iABC + 5) = - dOne
   iEq = iEq + 3
   iABC = iABC + 3
end do
! Boundary: passage for the last point
iEq = iEq + 1
dx = x(n) - x(n - 1)
SplineMat(iEq,iABC + 1) = dx
SplineMat(iEq,iABC + 2) = dx ** 2
SplineMat(iEq,iABC + 3) = dx ** 3
Coeff(iEq) = y(n) - y(n - 1)
! Boundary: zero third derivative in the last interval
SplineMat(iEq + 1,iABC + 3) = dOne
! Get the spline coefficients
call dgesv(3 * n - 3,1,SplineMat,3 * n - 3,ipiv,Coeff,3 * n - 3,info)
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
! Interpolate the function on the output grid
do i = 1, nn
   xt = xx(i)
   ! Get the interval
   do j = 1, n - 1
      if (xt .lt. x(j + 1)) exit
   end do
   if (xt .gt. x(n - 1)) j = n - 1
   ! Evaluate the function using the spline coefficients of the j-th interval
   dx = xt - x(j)
   iC = 3 * (j - 1) + 1
   yy(i) = y(j) + Coeff(iC) * dx + Coeff(iC + 1) * dx**2 + Coeff(iC + 2) * dx**3
end do
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end subroutine Spline




!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
double precision function AssociatedLegendreP(L,M,theta)
use sysparam        
integer, intent(in) :: L, M
double precision, intent(in) :: theta
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
if (abs(M) > L) then
   write(0,*) 'ERROR: requested a Legendre function with M > L'
   stop
end if
!
select case(L)
   case(0)  ! M = 0
      AssociatedLegendreP = sqrt(0.5d0)
   case(1)
      select case(M)
         case(0) 
            AssociatedLegendreP = sqrt(1.5d0) * cos(theta)
         case(1)
            AssociatedLegendreP = - sqrt(0.75) * sin(theta)
         case(-1)
            AssociatedLegendreP =   sqrt(0.75) * sin(theta)
      end select
   case(2)
      select case(M)
         case(0)
            AssociatedLegendreP = sqrt(0.15625) * (dOne + 3.d0 * cos(2.d0 * theta))
         case(1)
            AssociatedLegendreP = - sqrt(0.9375) * sin(2.d0 * theta)
         case(-1)
            AssociatedLegendreP =   sqrt(0.9375) * sin(2.d0 * theta)
         case(2)
            AssociatedLegendreP =   sqrt(0.9375) * sin(theta)**2
         case(-2)
            AssociatedLegendreP =   sqrt(0.9375) * sin(theta)**2
      end select
   case default
      write(0,*) 'Associated Legendre function not implemented for L,M = ', L, M
      stop
end select
!!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!   !!!!
return
end function AssociatedLegendreP
