module debug_utils
    use, intrinsic          :: iso_c_binding
    use                     :: xbeam_shared
    implicit none

    interface print_matrix
        module procedure print_1d_matrix_double, print_1d_matrix_int, print_3d_matrix_double, print_2d_matrix_double, &
                            print_2d_matrix_int,&
                            print_node, print_elem
    end interface
contains

    subroutine print_1d_matrix_double(name, matrix)
        character(len=*), intent(IN)        :: name
        real(c_double), intent(IN)           :: matrix(:)

        integer(c_int)                          :: unit
        integer(c_int)                          :: i


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix)
            write(unit,*) matrix(i)
        end do
        close(unit)
    end subroutine
    subroutine print_1d_matrix_int(name, matrix)
        character(len=*), intent(IN)        :: name
        integer(c_int), intent(IN)           :: matrix(:)

        integer(c_int)                          :: unit
        integer(c_int)                          :: i


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix)
            write(unit,*) matrix(i)
        end do
        close(unit)
    end subroutine
    subroutine print_2d_matrix_int(name, matrix)
        character(len=*), intent(IN)        :: name
        integer(c_int), intent(IN)           :: matrix(:,:)

        integer(c_int)                          :: unit
        integer(c_int)                          :: i


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix(:,1))
            write(unit,*) matrix(i, :)
        end do
        close(unit)
    end subroutine
    subroutine print_2d_matrix_double(name, matrix)
        character(len=*), intent(IN)        :: name
        real(c_double), intent(IN)           :: matrix(:,:)

        integer(c_int)                          :: unit
        integer(c_int)                          :: i


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix(:,1))
            write(unit,*) matrix(i, :)
        end do
        close(unit)
    end subroutine
    subroutine print_3d_matrix_double(name, matrix)
        character(len=*), intent(IN)        :: name
        real(c_double), intent(IN)           :: matrix(:,:,:)

        integer(c_int)                          :: unit
        integer(c_int)                          :: i
        integer(c_int)                          :: j


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix(:,1,1))
            write(unit,*) 'i = ', i
            do j=1, size(matrix(i,:,1))
				write(unit,*) matrix(i, j, :)
			end do
        end do
        close(unit)
    end subroutine
    subroutine print_node(name, matrix)
        character(len=*), intent(IN)        :: name
        type(xbnode), intent(IN)           :: matrix(:)
        integer(c_int)                          :: unit
        integer(c_int)                          :: i


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix)
            write(unit,*) '---'
            write(unit,*) 'node = ', i
            write(unit,*) 'Master'
            write(unit,*) matrix(i)%Master
            write(unit,*) 'Vdof'
            write(unit,*) matrix(i)%Vdof
            write(unit,*) 'Fdof'
            write(unit,*) matrix(i)%Fdof
        end do
        close(unit)
    end subroutine
    subroutine print_elem(name, matrix)
        character(len=*), intent(IN)        :: name
        type(xbelem), intent(IN)           :: matrix(:)
        integer(c_int)                          :: unit
        integer(c_int)                          :: i, j


        open(newunit=unit, file='debug_'//name//'.txt')
        do i=1, size(matrix)
            write(unit,*) '---'
            write(unit,*) 'elem = ', i
            write(unit,*) 'NumNodes'
            write(unit,*) matrix(i)%NumNodes
            write(unit,*) 'MemNo'
            write(unit,*) matrix(i)%MemNo
            write(unit,*) 'Conn'
            write(unit,*) matrix(i)%Conn
            write(unit,*) 'Master'
            do j=1,matrix(i)%NumNodes
				write(unit,*) matrix(i)%Master(j, :)
            end do
            write(unit,*) 'Length'
            write(unit,*) matrix(i)%Length
            write(unit,*) 'PreCurv'
            write(unit,*) matrix(i)%PreCurv
            write(unit,*) 'Psi'
            write(unit,*) matrix(i)%Psi
            write(unit,*) 'Vector'
            write(unit,*) matrix(i)%Vector
            write(unit,*) 'Mass'
            do j=1,6
				write(unit,*) matrix(i)%Mass(j,:)
			end do
			write(unit,*) 'InvStiff'
			do j=1,6
				write(unit,*) matrix(i)%InvStiff(j,:)
		    end do

        end do
        close(unit)
    end subroutine

end module debug_utils
