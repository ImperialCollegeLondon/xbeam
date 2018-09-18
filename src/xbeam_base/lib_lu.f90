!->Copyright by The University of Michigan, Aerospace Department. 2003
!
!->Module LIB_LU. Rafa Palacios. 22Jun2003
!
!->Description.-
!
!  Matrix operations using LU decomposition.
!
!->Subroutines.-
!
!   lu_bksubs:    LU back-substitution.
!   lu_decomp:    LU decomposition of a general square matrix.
!   lu_determ:    Determinant of a square matrix.
!   lu_invers:    Inversion of a square matrix.
!   lu_sparse:    Inverse of Sparse matrix/vector product
!
!-> Modifications.-
! 20120317 A.Da Ronch Subroutine lu_sparse added
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_lu
 implicit none
!
!   (There are no public variables in this module).
!
 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_BKSUBS.
!
!->Description.-
!
!   Solves the set of N linear equations A�X=B. A is given by its LU decomposition.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_bksubs (A, INDX, B)
!
!-> I/O Variables.
!
  real(8), intent(in)   :: A(:,:)
  integer, intent(in)   :: INDX(:)
  real(8), intent(inout):: B(:)
!
!-> Local Variables.
!
  integer :: N               ! Size of the matrices.
  integer :: I, II, LL, J
  real(8) :: Sum
!
      n= size(A,1)
!
      II = 0
      DO I = 1, N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0)THEN
          DO J = II, I-1
            SUM = SUM - A(I, J)*B(J)
          END DO
      ELSEIF(abs(SUM)>1e-20)THEN
          II = I
        ENDIF
        B(I) = SUM
      END DO
!
      DO I = N, 1, -1
        SUM = B(I)
        IF(I.LT.N)THEN
          DO J = I+1, N
            SUM = SUM - A(I, J)*B(J)
          END DO
        ENDIF
        B(I) = SUM/A(I, I)
      END DO
!
  return
 end subroutine lu_bksubs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_DECOMP.
!
!->Description.-
!
!   LU Decomposition of a NxN Matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_decomp (A, INDX, D)
!
!-> I/O Variables.
!
  real(8), intent(inout):: A(:,:)
  integer, intent(out)  :: INDX(:)
  real(8), intent(out)  :: D
!
!-> Local Variables.
!
  integer :: i,j,k                   ! Counters.
  integer :: imax                    ! Max Counter.
  integer :: N                       ! Size of the matrices.
  real(8) :: Aamax, Sum, Dum
  real(8), allocatable :: VV (:)     ! Stores the implicit scaling of each row.
!
  real(8), parameter :: Tiny=1.0d-20 ! Small number.
!
!
      n= size(A,1)
      allocate (VV(N))
      imax = 0
!
      D = 1.0
      DO I = 1, N
        AAMAX = 0.0
        DO J = 1, N
          IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        END DO
        IF(abs(AAMAX)<Tiny) then
            ! call backtrace()
            STOP 'Singular matrix'
        end if
        VV(I:I) = 1.d0/AAMAX
      END DO
!
      DO J = 1, N
        DO I = 1, J-1
          SUM = A(I, J)
          DO K = 1, I-1
            SUM = SUM - A(I, K)*A(K, J)
          END DO
          A(I, J) = SUM
        END DO
!
        AAMAX = 0.0
        DO I = J, N
          SUM = A(I, J)
          DO K = 1, J-1
            SUM = SUM - A(I, K)*A(K, J)
          END DO
          A(I, J) = SUM
          DUM = VV(I)*ABS(SUM)
          IF(DUM.GE.AAMAX)THEN
            IMAX = I
            AAMAX = DUM
          ENDIF
        END DO
        IF(J.NE.IMAX)THEN
          DO K = 1, N
            DUM = A(IMAX, K)
            A(IMAX, K) = A(J, K)
            A(J, K) = DUM
          END DO
          D = -D
          VV(IMAX) = VV(J)
        ENDIF
        INDX(J) = IMAX
        IF(abs(A(J, J))<Tiny) A(J, J) = TINY
        IF(J.NE.N)THEN
          DUM = 1.0/A(J, J)
          DO I = J+1, N
            A(I, J) = A(I, J)*DUM
          END DO
        ENDIF
      END DO
!
  deallocate (VV)
  return
 end subroutine lu_decomp

 function inv(A) result(Ainv)
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  integer :: i, j, m

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A, 1)
  m = size(A, 2)
 !  if (maxval(abs(A)) > 1e12) then
 !      open(unit=2, file='matrix.txt', ACTION="write", STATUS="replace")
 !      print*, 'UPPPSSSSS'
 !      do i=1, n
 !         write(2, '(*(2X e18.7))')((A(i,j)) ,j=1,M)
 !      end do
 !     close(2)
 !     read(*,*)
 ! end if
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv



!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_INVERS.
!
!->Description.-
!
!   Invert a matrix using LU Decomposition. It can also return the determinant.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_invers (Matrix, InvMatrix, Det)
!
!-> I/O Variable.
!
  real(8),intent(in) :: Matrix(:,:)
  real(8),intent(out):: InvMatrix(:,:)
  real(8),optional,intent(out):: Det

!-> Local variables.
!
  integer::i,N
  real(8) :: D
  real(8), allocatable :: MatrixLU(:,:)
  integer, allocatable :: Indexes(:)
!
! Set dimension and the auxiliary vector.
!
  N= size(Matrix,1)
  allocate (Indexes(N)); Indexes = 0
  allocate (MatrixLU(N,N))
!
! Initialize InvMatrix as the identity matrix and MatrixLU as Matrix.
!
  MatrixLU=Matrix
!
  InvMatrix=0.d0
  do i=1,N
    InvMatrix(i,i)=1.d0
  end do
!
! LU Decomposition.
!
  call lu_decomp (MatrixLU,Indexes,D)
!
! The inverse is the solution to units vectors.
!
  do i=1,N
    call lu_bksubs (MatrixLU,Indexes,InvMatrix(:,i))
  end do
!
! Return the determinant of the matrix.
!
  if (present(Det)) then
    Det=D
    do i=1,N
      Det=Det*MatrixLU(i,i)
    end do
  end if
!
  return
 end subroutine lu_invers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_DETERM.
!
!->Description.-
!
!   Invert a matrix using LU Decomposition.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lu_determ (Matrix, Det)
!
!-> I/O Variable.
!
  real(8), intent(in) :: Matrix(:,:)
  real(8), intent(out):: Det
!
!-> Local variables.
!
  integer::i,N
  real(8), allocatable :: MatrixLU(:,:)
  integer, allocatable :: Indexes(:)
!
! Set dimension and the auxiliary vector.
!
  N= size(Matrix,1)
  allocate (Indexes(N))
  allocate (MatrixLU(N,N)); MatrixLU=Matrix
!
! LU Decomposition.
!
  call lu_decomp (MatrixLU,Indexes,Det)
!
! The determinant of a LU decomposed matrix is just the product of the
! diagonal elements.
!
  do i=1,N
    Det=Det*MatrixLU(i,i)
  end do

  deallocate(MatrixLU)
  deallocate(Indexes)
!
  return
 end subroutine lu_determ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LU_SPARSE
!
!->Description.-
!
!   Inverse of sparse matrix/vector product.
!
!-> Modifications.-
! 20120317 A.Da Ronch Subroutine added for NOLAPACK conditional compilation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lu_sparse(dimSprMat, SprMat, b, X)
        use lib_sparse
        integer,     intent(in) :: dimSprMat     ! Current storage dimension
        type(sparse),intent(in) :: SprMat(:)     ! Sparse matrix
        real(8),     intent(in) :: b(:)          ! Forcing vector
        real(8),     intent(out):: X(:)          ! Solution vector
        real(8),allocatable     :: FulMat(:,:)   ! Full matrix
        real(8),allocatable     :: invFulMat(:,:)! Full matrix
        integer:: dimb

        dimb=size(b)

        ! allocate(invFulMat(dimb,dimb))
        !
        ! ! Calculate the inverse
        ! ! call lu_invers(FulMat, invFulMat)
        ! invFulMat = inv(SprMat(1)%a)
        !
        ! ! Calculate matrix-vector product
        ! ! X=0.0d0
        ! ! do i1=1,dimb
        ! ! !   do i2=1,dimb
        ! !         ! original
        ! !         X(i1) = X(i1) + dot_product(invFulMat(i1,:), b(:))
        ! !     ! end do
        ! ! end do
        ! X = MATMUL(invFulMat, b)
        !
        ! deallocate(invFulMat)
        call lu_solve(size(SprMat(1)%a, dim=1), SprMat(1)%a, b, X)
        return
    end subroutine lu_sparse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lu_solve(numdof,A, b, X, balancing)
        use, intrinsic :: IEEE_ARITHMETIC
        use, intrinsic :: ISO_C_BINDING
        integer, intent(IN)     :: numdof
        real(8), intent(IN)     :: A(numdof, numdof)
        real(8), intent(in)     :: b(numdof)          ! Forcing vector
        real(8), intent(out)    :: X(numdof)          ! Solution vector
        ! real(8),allocatable     :: invFulMat(:,:)! Full matrix
        logical(c_bool), optional       :: balancing
        integer:: dimb
        real(8)                 :: r(numdof)
        real(8)                 :: c(numdof)
        real(8)                 :: rcond
        real(8)                 :: ferr(1)
        real(8)                 :: berr(1)
        real(8)                 :: work(4*numdof)
        integer                 :: iwork(numdof)
        integer                 :: info
        integer                 :: ipivot(numdof)
        character               :: equed
        character               :: fact
        character               :: trans
        integer                 :: unit
        integer                 :: i

        ! real(8),allocatable     :: A_copy(:,:)
        real(8)                 :: A_copy(numdof, numdof)
        real(8)                 :: b_copy(numdof)
        real(8)                 :: AF(numdof, numdof)
        ! real(8)                 :: work(size(A, dim=1)*4)
        ! integer                 :: iwork(size(A, dim=1))
        integer                 :: nrows
        logical                 :: with_balancing

        if (any(ieee_is_nan(A))) then
             open(newunit=unit, file='debug_failed_Asys.txt')
             do i=1, size(A(:,1))
                 write(unit,*)A(i, :)
             end do
             close(unit)
            STOP 'NaN in lu_solve'
        end if

        A_copy = A
        X = b
        b_copy = b
        nrows = size(b)

        ! solve system
        with_balancing = .FALSE.
        if (present(balancing)) then
            if (balancing) then
                with_balancing = .TRUE.
            end if
        end if
        if (.NOT. with_balancing) then
            call dgesv(numdof,&
                       1,&
                       A_copy,&
                       numdof,&
                       ipivot,&
                       x,&
                       numdof,&
                       info)
        else
            fact = 'E'
            trans = 'N'
            call dgesvx(fact,&
                        trans,&
                        numdof,&
                        1,&
                        A_copy,&
                        numdof,&
                        AF,&
                        numdof,&
                        ipivot,&
                        equed,&
                        r,&
                        c,&
                        b_copy,&
                        numdof,&
                        x,&
                        numdof,&
                        rcond,&
                        ferr,&
                        berr,&
                        work,&
                        iwork,&
                        info)
        end if

        if (info /= 0) then
            print*, 'Info in DGESV = ', info
             open(newunit=unit, file='debug_failed_Asys.txt')
             do i=1, size(A(:,1))
                 write(unit,*)A(i, :)
             end do
             close(unit)
            stop
        end if
        ! deallocate(A_copy)
        return
    end subroutine lu_solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module lib_lu
