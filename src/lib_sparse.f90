!->Module LIB_SPARSE. Rafa Palacios. 30Aug2006
!                     Salvatore Maraniello 29 Sep 2014
!
!->Description.-
!
!  This module defines tools for handling sparse matrix.
!
!->Subroutines:
!
!   sparse_zero:      Initialize to zero a sparse matrix.
!   sparse_addmat:    Add submatrix to sparse matrix.
!   sparse_addval:    Add element to sparse matrix.
!   sparse_addsparse: Add sparse submatrix to a larger sparse matrix.
!   sparse_bandwidth: Compute bandwidth of sparse matrix.
!   sparse_copy2band: Write sparse matrix in banded form.
!   sparse_matvmul  : Multiply a sparse matrix by a vector.
!   sparse_precond:   Precondition linear system with sparse matrix.
!   sparse_precond2:  Precondition linear dynamic system with sparse matrix.
!
!-> New Subroutines
!   sparse_set_row_zero: given a sparse matrix, sets the specified rows to zero
!   sparse_set_rows_zero: delete a list of rows (optimised)
!   sparse_set_colrows_zero: delete a list of rows and columns (optimised)
!   sparse print_nonzero: print all elements in a sparse matrix in list form
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_sparse
 implicit none

! Define derived types.
  type sparse
    integer i
    integer j
    real(8) a
  end type

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ZERO.
!
!->Description.-
!
!   Initialize to zero a sparse matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_zero (DimArray,Matrix)
  integer,     intent(out):: DimArray  ! Current storage dimension.
  type(sparse),intent(out):: Matrix(:) ! Sparse matrix.

  DimArray=0
  Matrix(:)%i=0
  Matrix(:)%j=0
  Matrix(:)%a=0.d0

 return 
 end subroutine sparse_zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ADDVAL.
!
!->Description.-
!
!   Add element to sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_addval (i,j,Val,DimArray,Matrix,Factor)

  integer,     intent(in)   :: i, j      ! Coefficient in the matrix.
  real(8),     intent(in)   :: Val       ! Value to be added.
  integer,     intent(inout):: DimArray  ! Current storage dimension.
  type(sparse),intent(inout):: Matrix(:) ! Sparse matrix.
  real(8),optional,intent(in)::Factor    ! Factor that multiplies the new terms.

  ! Local variables.
  logical:: Flag     ! Logical flag.
  integer:: k        ! Counter of element
  real(8):: Val2     ! Modified val

  if (present(Factor)) then
    Val2=Val*Factor
  else
    Val2=Val
  end if

! If there is already a (i,j) term in the matrix, add the new value.
  Flag=.false.
  do k=1,DimArray
    if ((Matrix(k)%i.eq.i).and.(Matrix(k)%j.eq.j)) then
      Matrix(k)%a=Matrix(k)%a+Val2
      Flag=.true.
      exit
    end if
  end do

! If (i,j) was empty, create a new entry.

  if (.not.Flag) then
    if (DimArray.eq.size(Matrix)) STOP 'ERROR: Not enough memory for sparse matrix allocation defined in solver routine'

    DimArray=DimArray+1
    Matrix(DimArray)%i= i
    Matrix(DimArray)%j= j
    Matrix(DimArray)%a= Val2
  end if

  return
 end subroutine sparse_addval


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ADDMAT.
!
!->Description.-
!
!   Add submatrix to sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_addmat (i1,j1,Mat,DimArray,SprMat)

  integer,     intent(in)   :: i1,j1     ! First coefficient minus one on the SprMat.
  real(8),     intent(in)   :: Mat(:,:)  ! Value to be added.
  integer,     intent(inout):: DimArray  ! Current storage dimension.
  type(sparse),intent(inout):: SprMat(:) ! Sparse matrix.

  integer:: i,j    ! Counter on the element of the matrix.

  do i=1,size(Mat,DIM=1)
    do j=1,size(Mat,DIM=2)
      if (Mat(i,j).ne.0.d0) then
        call sparse_addval (i1+i,j1+j,Mat(i,j),DimArray,SprMat)
      end if
    end do
  end do

  return
 end subroutine sparse_addmat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_SPARSE2FULL.
!
!->Description.-
!
!   Convert sparse matrix to full matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_sparse2full (DimArray,SprMat,FullMat)

  integer,     intent(in)   :: DimArray      ! Current storage dimension.
  type(sparse),intent(in)   :: SprMat(:)     ! Sparse matrix.
  real(8),     intent(out)   :: FullMat(:,:)  ! Full matrix.

  integer:: k    ! Counter on the element of the matrix.

  do k=1,DimArray
    FullMat(SprMat(k)%i,SprMat(k)%j)=SprMat(k)%a
  end do

  return
 end subroutine sparse_sparse2full


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ADDSPARSE.
!
!->Description.-
!
!   Add sparse submatrix to sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_addsparse (i1,j1,DimSubmat,Submat,DimMat,Mat,Factor,FilterI,FilterJ)

  integer,     intent(in)   ::  i1,j1       ! First coefficient on the Mat.
  integer,     intent(in)   ::  DimSubmat   ! Current dimension of the submatrix.
  type(sparse),intent(in)   ::  Submat(:)   ! Matrix to be added.
  integer,     intent(inout)::  DimMat      ! Current storage dimension.
  type(sparse),intent(inout)::  Mat(:)      ! Sparse matrix.
  real(8),optional,intent(in):: Factor      ! Factor that multiplies the new terms.
  integer,optional,intent(in):: FilterI(:)  ! Filter on the rows of Submat.
  integer,optional,intent(in):: FilterJ(:)  ! Filter on the columns of Submat.

  integer:: i,j,k       ! Counters on the elements of the matrix.
  real(8)::ActualFactor ! Added RPN 25.04.2011 to avoid Factor being used without definition.

  if (present(Factor)) then
    ActualFactor=Factor
  else
    ActualFactor=1.d0
  end if

  do k=1,DimSubmat
    if (Submat(k)%a.ne.0.d0) then

      if (present(FilterI)) then
        i=FilterI(Submat(k)%i)
      else
        i=Submat(k)%i
      end if

      if (present(FilterJ)) then
        j=FilterJ(Submat(k)%j)
      else
        j=submat(k)%j
      end if

      if ((i.ne.0).and.(j.ne.0)) &
&       call sparse_addval (i1+i,j1+j,Submat(k)%a,DimMat,Mat,Factor=ActualFactor)
    end if
  end do

  return
 end subroutine sparse_addsparse



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_GETBANDWIDTH.
!
!->Description.-
!
!   Get lower and upper bandwidths of sparse matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_getbandwidth (DimArray,Matrix,KL,KU)

  integer,     intent(in) :: DimArray  ! Current storage dimension.
  type(sparse),intent(in) :: Matrix(:) ! Sparse matrix.
  integer,     intent(out):: KL        ! Number of subdiagonals of the banded form.
  integer,     intent(out):: KU        ! Number of superdiagonals of the banded form.

  integer :: k  ! Counter on the elements of the array.

  KL=0
  KU=0

  do k=1,DimArray
    if ((Matrix(k)%i-Matrix(k)%j).gt.KL) KL=Matrix(k)%i-Matrix(k)%j
    if ((Matrix(k)%j-Matrix(k)%i).gt.KU) KU=Matrix(k)%j-Matrix(k)%i
  end do

  return
 end subroutine sparse_getbandwidth


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_COPY2BAND.
!
!->Description.-
!
!   Copy matrix to bandwidth form used in LAPACK.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_copy2band (DimArray,Matrix,KL,KU,BandedMat)

  integer,     intent(in) :: DimArray       ! Current storage dimension.
  type(sparse),intent(in) :: Matrix(:)      ! Sparse matrix.
  integer,     intent(in) :: KL             ! Number of superdiagonals of the banded form.
  integer,     intent(in) :: KU             ! Number of subdiagonals of the banded form.
  real(8),     intent(out):: BandedMat(:,:) ! Number of superdiagonals of the banded form.

  integer :: k  ! Counter on the elements of the array.

  do k=1,DimArray
    BandedMat(KL+KU+1+Matrix(k)%i-Matrix(k)%j,Matrix(k)%j) = Matrix(k)%a
  end do

  return
 end subroutine sparse_copy2band


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_MATVMUL.
!
!->Description.-
!
!   Multiply a sparse matrix by a vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function sparse_matvmul (DimArray,SprMat,DimVector,Vector)

! I/O Variables.
  integer,     intent(in):: DimArray    ! Current storage dimension.
  type(sparse),intent(in):: SprMat(:)   ! Sparse matrix.
  integer,     intent(in):: DimVector   ! Dimension of the output vector.
  real(8),     intent(in):: Vector(:)   ! Column vector to be postmultiplied.
  real(8)                :: sparse_matvmul(DimVector)

! Local variables.
  integer:: i,j,k    ! Counters on the element of the matrix.

  sparse_matvmul=0.d0

  do k=1,DimArray
    i=SprMat(k)%i
    j=SprMat(k)%j
    sparse_matvmul(i)= sparse_matvmul(i) + SprMat(k)%a*Vector(j)
  end do

  return
 end function sparse_matvmul



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine SPARSE_PRECOND
!
!-> Description:
!
!    Equation preconditioner.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_precond (ks,Kglobal,Qglobal,MinDelta)

! I/O Variables.
  integer,intent(inout)::      ks             ! Current storage size of stiffness matrix.
  type(sparse),intent(inout):: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  real(8),intent(inout)::      Qglobal(:)     ! Global vector of discrete generalize forces.
  real(8),intent(in)   ::      MinDelta       ! Threshold for coefficients set to zero.

! Local variables.
  integer             :: k,k1                 ! Counters.
  real(8),allocatable :: Kmax(:)              ! Maximum values of each row in K.
  real(8)             :: Val                  ! Auxiliar number.

! Obtain maximum value of each row of the stiffness.
  allocate(Kmax(size(Qglobal))); Kmax=0.d0
  do k=1,ks
    if (abs(Kglobal(k)%a).gt.abs(Kmax(Kglobal(k)%i))) Kmax(Kglobal(k)%i)=Kglobal(k)%a
  end do

! Normalize each row of K by its maximum and remove elements below a threshold.
  k=1
  do while (k.le.ks)
    Val=Kglobal(k)%a/Kmax(Kglobal(k)%i)

    if (abs(Val).gt.MinDelta) then
      Kglobal(k)%a=Val
      k=k+1
    else
      do k1=k+1,ks
        Kglobal(k1-1)%i=Kglobal(k1)%i
        Kglobal(k1-1)%j=Kglobal(k1)%j
        Kglobal(k1-1)%a=Kglobal(k1)%a
      end do
      ks=ks-1
    end if
  end do

! Normalize force vector.
  do k=1,size(Qglobal)
    Qglobal(k)=Qglobal(k)/Kmax(k)
  end do

  return
 end subroutine sparse_precond



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine SPARSE_PRECOND2
!
!-> Description:
!
!    Dynamic equation preconditioner.
!
!-> Remarks.-
!
!    1) Equations are normalised by the max value in the stiffness matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_precond2 (ks,Kglobal,ms,Mglobal,Qglobal,MinDelta)

! I/O Variables.
  integer,intent(inout)::      ks             ! Current storage size of stiffness matrix.
  type(sparse),intent(inout):: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  integer,intent(inout)::      ms             ! Current storage size of mass matrix.
  type(sparse),intent(inout):: Mglobal(:)     ! Global mass matrix in sparse storage.
  real(8),intent(inout)::      Qglobal(:)     ! Global vector of discrete generalize forces.
  real(8),intent(in)   ::      MinDelta       ! Threshold for coefficients set to zero.

! Local variables.
  integer             :: j,j1                 ! Counters.
  real(8),allocatable :: Kmax(:)              ! Maximum values of each row in K.
  real(8)             :: Val                  ! Auxiliar number.

! Obtain maximum value of each row of the stiffness.
  allocate(Kmax(size(Qglobal))); Kmax=0.d0
  do j=1,ks
    if (abs(Kglobal(j)%a).gt.abs(Kmax(Kglobal(j)%i))) Kmax(Kglobal(j)%i)=Kglobal(j)%a
  end do

! Normalize each row of K by its maximum and remove elements below a threshold.
  j=1
  do while (j.le.ks)
    Val=Kglobal(j)%a/Kmax(Kglobal(j)%i)

    if (abs(Val).gt.MinDelta) then
      Kglobal(j)%a=Val
      j=j+1
    else
      do j1=j+1,ks
        Kglobal(j1-1)%i=Kglobal(j1)%i
        Kglobal(j1-1)%j=Kglobal(j1)%j
        Kglobal(j1-1)%a=Kglobal(j1)%a
      end do
      ks=ks-1
    end if
  end do

! Normalize each row of M by its maximum and remove elements below a threshold.
  j=1
  do while (j.le.ms)
    Val=Mglobal(j)%a/Kmax(Mglobal(j)%i)
    j=j+1
  end do

! Normalize force vector.
  do j=1,size(Qglobal)
    Qglobal(j)=Qglobal(j)/Kmax(j)
  end do

  return
 end subroutine sparse_precond2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine sparse_set_rows_unit
!
!-> Description:
!
!    Given a sparse matrix, the routine delete allocates unit values along
!    the diagonal corresponding to the rows specified in rows_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_set_rows_unit (rows_list,DimArray,SprMat)

  integer,     intent(in)   :: rows_list(:)! rows for which all entries of sparse matrix will be set to zero
  integer,     intent(inout):: DimArray      ! Storage dimension of sparse matrix
  type(sparse),intent(inout):: SprMat(:)     ! Sparse matrix.

  integer:: nn, rr

  do nn=1,size(rows_list)
    rr = rows_list(nn)
    call sparse_addval(rr,rr,1.d0,DimArray,SprMat,1.d0)
  end do

  return
 end subroutine sparse_set_rows_unit





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine sparse_set_rows_zero
!
!-> Description:
!
!    Given a sparse matrix, the routine delete all the elements from a certain
!    list of rows
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_set_rows_zero (rows_to_del,DimArray,SprMat)

  integer,     intent(in)   :: rows_to_del(:)! rows for which all entries of sparse matrix will be set to zero
  integer,     intent(inout):: DimArray      ! Storage dimension of sparse matrix
  type(sparse),intent(inout):: SprMat(:)     ! Sparse matrix.

  integer:: nn, jj
  logical:: DelRow


  nn=1
  do while (nn<=DimArray)
    DelRow=.false.
    ! check if row is in the list
    jj=1
    do while (jj<=size(rows_to_del))
      if (SprMat(nn)%i==rows_to_del(jj)) then
        DelRow=.true.
      end if
      jj=jj+1
    end do

    !if (DelRow .eqv. .true.) then
    if (DelRow) then
      jj=nn+1
      do while (jj<=DimArray)
        SprMat(jj-1)=SprMat(jj)
        jj=jj+1
      end do
      SprMat(DimArray)%i=0; SprMat(DimArray)%j=0; SprMat(DimArray)%a=0.d0
      DimArray=DimArray-1
    else
      nn=nn+1
    end if
  end do


  return
 end subroutine sparse_set_rows_zero



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine sparse_set_colrows_zero
!
!-> Description:
!
!    Given a sparse matrix, the routine delete all the elements from the rows
!    and columns listed in colrows_to_del
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_set_colrows_zero (colrows_to_del,DimArray,SprMat)

  integer,     intent(in)   :: colrows_to_del(:) ! col & rows for which all entries of sparse matrix will be set to zero
  integer,     intent(inout):: DimArray          ! Storage dimension of sparse matrix
  type(sparse),intent(inout):: SprMat(:)         ! Sparse matrix.

  integer:: nn, jj
  logical:: Del


  nn=1
  do while (nn<=DimArray)
    Del=.false.
    ! check if row/col is in the list
    jj=1
    do while (jj<=size(colrows_to_del))
      if ( (SprMat(nn)%i==colrows_to_del(jj)) .or.  (SprMat(nn)%j==colrows_to_del(jj)) ) then
        Del=.true.
      end if
      jj=jj+1
    end do

    if (Del) then
      jj=nn+1
      do while (jj<=DimArray)
        SprMat(jj-1)=SprMat(jj)
        jj=jj+1
      end do
      SprMat(DimArray)%i=0; SprMat(DimArray)%j=0; SprMat(DimArray)%a=0.d0
      DimArray=DimArray-1
    else
      nn=nn+1
    end if
  end do


  return
 end subroutine sparse_set_colrows_zero






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine sparse_set_row_zero
!
!-> Description:
!
!    Given a sparse matrix, the routine delete all the elements from a certain
!    rows
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_set_row_zero (row_to_del,DimArray,SprMat)

  integer,     intent(in)   :: row_to_del! rows for which all entries of sparse matrix will be set to zero
  integer,     intent(inout):: DimArray  ! Storage dimension of sparse matrix
  type(sparse),intent(inout):: SprMat(:) ! Sparse matrix.

  integer:: nn, jj
  logical:: DelEntry(DimArray)


  !DelEntry(DimArray)=.false.

  nn=1
  do while (nn<=DimArray)
    if (SprMat(nn)%i == row_to_del) then
      !SprMat(nn)%a = 0.d0
      !DelEntry(nn)=.true.
      jj=nn+1
      do while (jj<=DimArray)
        SprMat(jj-1)=SprMat(jj)
        jj=jj+1
      end do
      SprMat(DimArray)%i=0; SprMat(DimArray)%j=0; SprMat(DimArray)%a=0.d0
      DimArray=DimArray-1
    else
      nn=nn+1
    end if
  end do


  return
 end subroutine sparse_set_row_zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine print_sparse
!
!-> Description:
!
!    Prints non-zero elements inside a sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_print_nonzero (DimArray,SprMat)

  integer,     intent(in):: DimArray  ! Storage dimension of sparse matrix
  type(sparse),intent(in):: SprMat(:) ! Sparse matrix.

  integer:: nn

  write (*,'(A4,A4,A12)'), 'i', 'j', 'val'

  nn=1
  do nn=1,DimArray
    write (*,'(1I4,1I4,1PE12.3)'), SprMat(nn)%i, SprMat(nn)%j, SprMat(nn)%a
  end do

  return
 end subroutine sparse_print_nonzero



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine print_mat
!
!-> Description:
!
!    Prints all elements (also zero) of a sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_print_mat (DimArray,SprMat)

  integer,     intent(in):: DimArray  ! Storage dimension of sparse matrix
  type(sparse),intent(in):: SprMat(:) ! Sparse matrix.

  integer:: nn, mm, Imax, Jmax, ii, jj, Nfound
  logical :: Found

  !write (*,'(1I4,1I4,1PE12.3)'), SprMat(nn)%i, SprMat(nn)%j, SprMat(nn)%a
  !write (*,'(A4,A4,A12)'), 'i', 'j', 'val'

  ! get matrix size:
  Imax = maxval(SprMat(1:DimArray)%i)
  Jmax = maxval(SprMat(1:DimArray)%j)


  Nfound=0

  do ii=1,Imax
    do jj=1,Jmax

      ! find element
      Found=.false.
      nn=1
      do while ( (nn<=DimArray) .and. (Found .eqv. .false.) )
        if ( (SprMat(nn)%i == ii) .and. (SprMat(nn)%j == jj) ) then
          Found=.true.
          Nfound=nn
        end if
        nn=nn+1
      end do

      ! print
      if ( Found .eqv. .true. ) then
        write(*,'(1PE12.3,2X,$)'), SprMat(Nfound)%a
      else
        write(*,'(1PE12.3,2X,$)'), 0.d0
      end if

    end do

    ! go next row
    write(*,'(1X)')

  end do

  ! go next row
  write(*,'(1X)')

  return
 end subroutine sparse_print_mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_sparse
