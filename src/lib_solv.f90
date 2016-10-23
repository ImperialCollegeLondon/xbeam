!->Module LIB_SOLV. Salvatore Maraniello. 04Sep2014
!
!->Description.-
!
!  This module contains tools for the solution of non-linear systems
!
!->Subroutines:
! - residual_check: checks Newtonm iteration convergence via residual
! - error_check: checks convergence making an estimation of the error
! - delta_check: original convergence check, based on the maximum residual value
!   and the maximum incremental value
! - solv_set_vec_rows_zero: (analogous to lib_sparse-sparse_set_rows_to_zero)
!      Given a vector/matrix, sets rows specified in rows_list to zero.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module lib_solv

  implicit none

 contains


! subroutine residual_check
! ------------------------------------------------------------------------------
    ! Checks the convergence of a Newton iteration using the criteria
    !
    ! Res < Tr Ro + Ta
    !
    ! where Res and Ro are the current and initial residual. The check is such
    ! that:
    ! 1. if the relative residual Res/Res0 < Tr, the test is always passed
    ! 2. if the relative residual is Res/Res0 > Tr, the tolerance for the
    !    absolute residual is increased of a factor as much higher as closer
    !    the relative residual gets to the relative tolerance Tr.
    !    The test is in fact passed if:
    !    Res < Ta / ( 1 - Tr/(Res/Res0) )
    ! 3. In the event that Res/Res0 ~ 1 (e.g. initial guess very close to
    !    solution) the test is passed if Res < Ta
    ! 4. If v is dimensional, Ta should be a dimensional quantity as well
    !
    !
    ! Remarks:
    ! 1. Res0 is calculated only if Iter = 1. Otherwise the input value is used.
    ! --------------------------------------------------------------------------
  subroutine residual_check(Iter,v,Res,Res0,passed,Tr,Ta,print_info)

    ! input:
    integer, intent(in)           :: Iter        ! iteration number
    real(8), intent(in)           :: v(:)
    ! in/out:
    real(8), intent(inout)        :: Res0        ! Residual at iteration 1.
    ! optional
    real(8),             optional :: Tr, Ta      ! relative and absolute tolerance
    logical,             optional :: print_info  ! print update or not
    ! output:
    logical, intent(out)          :: passed      ! outcome of convergence test
    real(8), intent(out)          :: Res ! absolute residual


    ! default values for optional input
    if ( present(Tr).eqv. .false.) then
        Tr = 1e-6
    end if
    if ( present(Ta).eqv. .false.) then
        Ta = 1e-6
    end if
    if ( present(print_info).eqv. .false.) then
        print_info = .true.
    end if

    ! current residual
    Res = vector_norm(v)      ! sqrt( dot_product(v,v)/dble(NumDof) )

    ! initial residual: shared variable
    ! The .le. condition allows to recompute Res0 for quantities for which
    ! the initial Ro may be too close to zero.
    ! if (Iter .eq. 1) then
    if (Iter .eq. 1) then
        Res0 = Res ! sqrt( dot_product(v,v)/dble(NumDof) )
    end if

    ! check convergence
    passed = .false.
    if ( Res < Tr*Res0 + Ta ) then
        passed = .true.
    end if

    !if (print_info .eqv. .true.) write (*,'(2X,A,1PE10.3,A,1PE10.3,$)') &
    !&                           'Res=',Res, ' ResRel=',ResRel
    if (print_info .eqv. .true.) write (*,'(2X,1PE10.3,2X,1PE10.3,$)') ,Res, Res/Res0

  end subroutine residual_check


! subroutine error_check
! ------------------------------------------------------------------------------
    ! Check the convergence of a Newton iteration using an approximation of the
    ! error. The check is
    ! passed if:
    !
    ! E_i < Ta
    !
    ! where E is an approximation of the error at the i-th iteration and Ta is
    ! the tolerance level (note that the eq. above is dimensional, therefore Ta
    ! is an absolute tolerance). The error is approximated as:
    !
    ! E_i ~ |dx_{i-1}|**2 / |dx_i|
    !
    ! where:
    ! dx_i: is the delta step in the solution x at iteriation i
    ! |.| : norm computed as per vector_norm function
    !
    ! Remarks:
    ! 1. The approximation holds for superlinear methods only.
    ! 2. The methods passes in output the norm dx_i and required in input dx_{i-1}
    ! 3. The first check is possible only at the 2nd iteration
    !
    ! --------------------------------------------------------------------------
    subroutine error_check(Iter,dx,DX_old,DX_now,Er,passed,Ta,print_info)

    ! input:
    integer, intent(in)   :: Iter        ! iteration number
    real(8), intent(in)   :: dx(:)       ! solution delta
    real(8), intent(in)   :: DX_old      ! norm of dx at previous iteration
    ! optional
    real(8), optional :: Ta               ! tolerance level
    logical, optional :: print_info      ! print update or not
    ! output:
    real(8), intent(out)   :: DX_now     ! norm of dx at current iteration
    real(8), intent(out)   :: Er         ! Error at current iteration
    logical, intent(out)   :: passed     ! outcome of convergence test


    ! default values for optional input
    if ( present(Ta).eqv. .false.) then
        Ta = 1e-9
    end if
    if ( present(print_info).eqv. .false.) then
        print_info = .true.
    end if

    ! vector norm
    DX_now = vector_norm(dx)

    ! check convergence
    passed = .false.
    if (Iter > 1) then
        Er = DX_old**2 / DX_now
        if ( Er < Ta ) then
            passed = .true.
        end if
        if (print_info .eqv. .true.) write (*,'(2X,1PE10.3,$)') Er
    else
        if (print_info .eqv. .true.) write (*,'(2X,9X,$)')
    end if

    end subroutine error_check


! subroutine delta_check
! ------------------------------------------------------------------------------
    ! Original convergence check used for non-linear solvers of cbem3_solv.
    ! The chekc is passed if:
    !
    ! delta = max(abs(r)) + max(abs(dx)) < T
    !
    ! where ideally:
    ! 1. r and dx are residual and delta x vectors (where x is the  solution).
    !    Ideally, r_vec and dx_vec would be non-dimensional, but the original
    !    check does not take care of that.
    ! 2. T tolerance level
    ! --------------------------------------------------------------------------
  subroutine delta_check(r,dx,delta,passed,T,print_info)

    ! input
    real(8), intent(in)  :: r(:), dx(:) ! residual and delta solution vectors
    real(8), optional    :: T           ! tolerance
    logical, optional    :: print_info  ! print on screenupdate or not
    ! output
    logical, intent(out) :: passed
    real(8), intent(out) :: delta

    ! default values for optional input
    if ( present(T).eqv. .false.) then
        T = 1e-5
    end if
    if ( present(print_info).eqv. .false.) then
        print_info = .true.
    end if

    ! check convergence
    passed = .false.
    delta = maxval(abs(r))+maxval(abs(dx))
    if (  delta < T ) then
        passed = .true.
    end if

    !if (print_info .eqv. .true.) write (*,'(2X,A,1PE10.3,A,1PE10.3,$)') &
    ! &                           'DeltaF=',maxval(abs(Qglobal)), ' DeltaX=',maxval(abs(DeltaX))
    if (print_info .eqv. .true.) write (*,'(2X,1PE10.3,2X,1PE10.3,$)') maxval(abs(r)), maxval(abs(dx))

  end subroutine delta_check



! function vector_norm
! ------------------------------------------------------------------------------
  ! Computes the vector norm n such that:
  ! n**2 = int_0^1{v**2}dx
  ! ----------------------------------------------------------------------------
  function vector_norm(v)

    ! input:
    real(8), intent(in)  :: v(:)  ! vector
    ! output:
    real(8) :: vector_norm        ! norm of the vector
    ! internal
    integer  :: N                 ! size of v

    N=size(v)
    vector_norm = sqrt( dot_product(v,v)/dble(N) )

  end function vector_norm



! subroutine separate_dofs(v,avec,bvec,va,vb)
! ------------------------------------------------------------------------------
  ! Given a vector v, the routine separates the elements of v into 2 different
  ! arrays, va and vb. avec and bvec are arrays of arbitrary length and such
  ! that:
  ! va( nn*(len(avec)-1) + avec(ii) ) = v( nn*(len(avec)+len(bvec))-1 + avec(ii))
  ! with nn = 1 .... len(v)/( len(avec)+len(bvec) )
  !
  ! e.g.
  ! v: displacement solution
  ! va: solution related to the translational dof
  ! vb: solution related to rotational dof
  ! avec: (/1,2,3/) - given the dof associated to the node nn, the translational dof are the first 3
  ! bvec: (/4,5,6/) - given the dof associated to the node nn, the rotational dof are the last 3
  !
  ! remark:
  ! - output vectors va and vb need to be allocated outside the subroutine
  ! ----------------------------------------------------------------------------
  subroutine separate_dofs(v,avec,bvec,va,vb)

    ! input
    real(8), intent(in) :: v(:)             ! input vector
    integer, intent(in) :: avec(:), bvec(:) ! positions vectors
    ! output
    real(8), intent(inout) :: va(:), vb(:)
    ! internal
    integer :: Na, Nb, Ntot ! length of avec, bvec and repetition pattern in v
    integer :: kk, ii       ! counters

    Na = size(avec)
    Nb = size(bvec)
    Ntot = Na + Nb

     kk=0
     do while (kk < size(v)/Ntot) ! Divide Displ from Rotations

        do ii=1,Na
            va(Na*kk+ii)=v(Ntot*kk+avec(ii))
        end do

        do ii=1,Nb
            vb(Nb*kk+ii)=v(Ntot*kk+bvec(ii))
        end do

        kk = kk+1
     end do

  end subroutine separate_dofs




! subroutine solv_set_vec_rows_zero
! ------------------------------------------------------------------------------
  ! Given a vector v, the routine sets to zero the values on the rows specified
  ! in rows_list
  !-----------------------------------------------------------------------------

 subroutine solv_set_vec_rows_zero (rows_list,v)

  integer,     intent(in)   :: rows_list(:)! rows for which all entries of sparse matrix will be set to zero
  real(8),     intent(inout) :: v(:)

  integer:: nn, rr

  do nn=1,size(rows_list)
    rr = rows_list(nn)
    v(rr)=0.d0
  end do

  return
 end subroutine solv_set_vec_rows_zero



! subroutine print_mat
! ------------------------------------------------------------------------------
  ! Print entried of rank 2 array (for debugging).
  !-----------------------------------------------------------------------------

 subroutine print_mat(mat,Imax,Jmax)
  real(8), intent(in) :: mat(:,:)
  integer :: ii, jj, Imax, Jmax

  do ii=1,Imax
    do jj=1,Jmax
      write (*,'(2X,F5.2,$)') mat(ii,jj)
    end do
    print *, ' '
  end do

 end subroutine print_mat


end module lib_solv
