module lib_mat


contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !->Subroutine mat_ADDMAT.
    !
    !->Description.-
    !
    !   Add submatrix to sparse matrix
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine mat_addmat (i1,j1,in_Mat,Mat)
    integer,     intent(in)   :: i1,j1     ! First coefficient minus one on the mat
    real(8),     intent(in)   :: in_Mat(:,:)  ! Value to be added.
    real(8),target,     intent(inout):: mat(:, :) ! mat matrix.

    real(8), pointer          :: original(:,:)

    original => mat(i1+1:i1 + size(in_Mat,1), j1+1:j1 + size(in_Mat,2))
    original = original + in_Mat
    nullify(original)
end subroutine mat_addmat

end module lib_mat
