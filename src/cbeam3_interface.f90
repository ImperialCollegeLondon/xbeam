module cbeam3_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: cbeam3_solv
    use                                 :: debug_utils
    use                                 :: lib_sparse

    implicit none

    integer(c_int), parameter, private  :: max_elem_node = MaxElNod

contains
    subroutine cbeam3_solv_nlnstatic_python(n_elem,&
                                            n_node,&
                                            num_nodes,&
                                            mem_number,&
                                            conn,&
                                            master,&
                                            n_mass,&
                                            mass_db,&
                                            mass_indices,&
                                            n_stiffness,&
                                            stiffness_db,&
                                            inv_stiffness_db,&
                                            stiffness_indices,&
                                            for_delta,&
                                            rbmass,&
                                            master_node,&
                                            vdof,&
                                            fdof,&
                                            options,&
                                            pos_ini,&
                                            psi_ini,&
                                            pos_def,&
                                            psi_def,&
                                            app_forces&
                                            ) bind(C)
        integer(c_int), intent(IN)      :: n_elem
        integer(c_int), intent(IN)      :: n_node

        ! elem data
        integer(c_int), intent(IN)      :: num_nodes(n_elem)
        integer(c_int), intent(IN)      :: mem_number(n_elem)
        integer(c_int), intent(IN)      :: conn(n_elem, max_elem_node)
        integer(c_int), intent(IN)      :: master(n_elem, max_elem_node, 2)
        integer(c_int), intent(IN)      :: n_mass
        real(c_double), intent(IN)      :: mass_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: mass_indices(n_elem)
        integer(c_int), intent(IN)      :: n_stiffness
        real(c_double), intent(IN)      :: stiffness_db(n_stiffness, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_stiffness, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(n_elem, max_elem_node, 3)
        real(c_double), intent(IN)      :: rbmass(n_elem, max_elem_node, 6, 6)

        ! node data
        integer(c_int), intent(IN)      :: master_node(n_node, 2)
        integer(c_int), intent(IN)      :: vdof(n_node)
        integer(c_int), intent(IN)      :: fdof(n_node)

        type(xbelem)                    :: elements(n_elem)
        type(xbnode)                    :: nodes(n_node)
        type(xbopts), intent(INOUT)     :: options

        real(c_double), intent(IN)      :: pos_ini(n_node, 3)
        real(c_double), intent(IN)      :: psi_ini(n_elem, max_elem_node, 3)
        real(c_double), intent(INOUT)   :: pos_def(n_node, 3)
        real(c_double), intent(INOUT)   :: psi_def(n_elem, max_elem_node, 3)

        real(c_double), intent(IN)      :: app_forces(n_node, 6)
        ! ADC XXX: careful with forces in master FoR

        integer(c_int)                  :: num_dof
        real(c_double)                  :: applied_forces(n_node, 6)!legacy var

        integer(c_int)                  :: nodes_per_elem
        integer(c_int)                  :: i

        ! call print_matrix('conn',conn)
        ! call print_matrix('pos_ini', pos_ini)
        ! call print_matrix('psi_ini1', psi_ini(:, 1, :))
        ! call print_matrix('psi_ini2', psi_ini(:, 2, :))
        ! call print_matrix('psi_ini3', psi_ini(:, 3, :))
        ! call print_matrix('app_forces', app_forces)
        ! call print_matrix('fdof',fdof)
        ! call print_matrix('vdof',vdof)
        ! call print_matrix('master1',master(:,:,1))
        ! call print_matrix('master2',master(:,:,2))
        ! call print_matrix('masternode',master_node)

        num_dof = count(vdof > 0)*6
        applied_forces = app_forces
        ! do i=1, n_node
        !     print*, applied_forces(i,:)
        ! end do
        ! do i=1, num_app_forces
        !     ! applied_forces(node_app_forces(i), :) = app_forces(i, :)
        !     applied_forces(node_app_forces(i), :) = &
        !             applied_forces(node_app_forces(i), :) + app_forces(i, :)
        ! end do

        ! gaussian nodes
        nodes_per_elem = count(conn(1,:) /= 0)
        options%NumGauss = nodes_per_elem - 1

        ! do i=1, n_node
        !     print*, pos_def(i,:)
        ! end do
        ! print*, '---'
        ! do i=1, n_node
        !     print*, pos_ini(i,:)
        ! end do

        elements = generate_xbelem(n_elem,&
                                   num_nodes,&
                                   mem_number,&
                                   conn,&
                                   master,&
                                   n_mass,&
                                   mass_db,&
                                   mass_indices,&
                                   n_stiffness,&
                                   stiffness_db,&
                                   inv_stiffness_db,&
                                   stiffness_indices,&
                                   for_delta,&
                                   psi_ini,&
                                   rbmass)

        nodes = generate_xbnode(n_node,&
                                master_node,&
                                vdof,&
                                fdof)

        ! call print_xbelem(elements)
        ! call print_xbnode(nodes)
       ! overloaded function
       call cbeam3_solv_nlnstatic(num_dof,&
                                  elements,&
                                  nodes,&
                                  applied_forces,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  options &
                                  )
        ! call output_elems (elements,pos_def,psi_def)
        !
        ! call print_matrix('pos_def', pos_def)
        ! call print_matrix('psi_def1', psi_def(:, 1, :))
        ! call print_matrix('psi_def2', psi_def(:, 2, :))
        ! call print_matrix('psi_def3', psi_def(:, 3, :))
    end subroutine cbeam3_solv_nlnstatic_python


    subroutine cbeam3_solv_nlndyn_python   (n_elem,&
                                            n_node,&
                                            n_tsteps,&
                                            time,&
                                            num_nodes,&
                                            mem_number,&
                                            conn,&
                                            master,&
                                            n_mass,&
                                            mass_db,&
                                            mass_indices,&
                                            n_stiffness,&
                                            stiffness_db,&
                                            inv_stiffness_db,&
                                            stiffness_indices,&
                                            for_delta,&
                                            rbmass,&
                                            master_node,&
                                            vdof,&
                                            fdof,&
                                            options,&
                                            pos_ini,&
                                            psi_ini,&
                                            pos_def,&
                                            psi_def,&
                                            app_forces,&
                                            dynamic_forces_amplitude,&
                                            dynamic_forces_time,&
                                            forced_vel,&
                                            forced_acc,&
                                            pos_def_history,&
                                            psi_def_history,&
                                            pos_dot_def_history,&
                                            psi_dot_def_history)bind(C)

        integer(c_int), intent(IN)      :: n_elem
        integer(c_int), intent(IN)      :: n_node
        integer(c_int), intent(IN)      :: n_tsteps
        real(c_double), intent(IN)      :: time(n_tsteps)

        ! elem data
        integer(c_int), intent(IN)      :: num_nodes(n_elem)
        integer(c_int), intent(IN)      :: mem_number(n_elem)
        integer(c_int), intent(IN)      :: conn(n_elem, max_elem_node)
        integer(c_int), intent(IN)      :: master(n_elem, max_elem_node, 2)
        integer(c_int), intent(IN)      :: n_mass
        real(c_double), intent(IN)      :: mass_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: mass_indices(n_elem)
        integer(c_int), intent(IN)      :: n_stiffness
        real(c_double), intent(IN)      :: stiffness_db(n_mass, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(n_elem, max_elem_node, 3)
        real(c_double), intent(IN)      :: rbmass(n_elem, max_elem_node, 6, 6)

        ! node data
        integer(c_int), intent(IN)      :: master_node(n_node, 2)
        integer(c_int), intent(IN)      :: vdof(n_node)
        integer(c_int), intent(IN)      :: fdof(n_node)

        ! data structures to be reconstructed
        type(xbelem)                    :: elements(n_elem)
        type(xbnode)                    :: nodes(n_node)
        type(xbopts), intent(INOUT)     :: options

        real(c_double), intent(IN)      :: pos_ini(n_node, 3)
        real(c_double), intent(IN)      :: psi_ini(n_elem, max_elem_node, 3)
        real(c_double), intent(OUT)     :: pos_def_history(n_tsteps, n_node, 3)
        real(c_double), intent(OUT)     :: psi_def_history(n_tsteps, n_elem, max_elem_node, 3)
        real(c_double), intent(OUT)     :: pos_dot_def_history(n_tsteps, n_node, 3)
        real(c_double), intent(OUT)     :: psi_dot_def_history(n_tsteps, n_elem, max_elem_node, 3)

        real(c_double), intent(IN)      :: app_forces(n_node, 6)
        ! ADC: careful, forces in master FoR

        ! dynamic forces
        real(c_double), intent(IN)      :: dynamic_forces_amplitude(n_node, 6)
        real(c_double), intent(IN)      :: dynamic_forces_time(n_tsteps)
        real(c_double), intent(IN)      :: forced_vel(n_tsteps, 6)
        real(c_double), intent(IN)      :: forced_acc(n_tsteps, 6)

        integer(c_int)                  :: num_dof
        real(c_double)                  :: applied_forces(n_node, 6)! static
                                                                    ! loads

        real(c_double)                  :: pos_def(n_node, 3)
        real(c_double)                  :: pos_dot_def(n_node, 3)
        real(c_double)                  :: psi_def(n_elem, max_elem_node, 3)
        real(c_double)                  :: psi_dot_def(n_elem, max_elem_node, 3)
        integer(c_int)                  :: i_out
        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem


        num_dof = count(vdof > 0)*6
        ! applied_forces = 0.0_c_double
        ! do i=1, num_app_forces
        !     applied_forces(node_app_forces(i), :) = app_forces(i, :)
        ! end do
        applied_forces = app_forces

        ! gaussian nodes
        nodes_per_elem = count(conn(1,:) /= 0)
        options%NumGauss = nodes_per_elem - 1

        elements = generate_xbelem(n_elem,&
                                   num_nodes,&
                                   mem_number,&
                                   conn,&
                                   master,&
                                   n_mass,&
                                   mass_db,&
                                   mass_indices,&
                                   n_stiffness,&
                                   stiffness_db,&
                                   inv_stiffness_db,&
                                   stiffness_indices,&
                                   for_delta,&
                                   psi_ini,&
                                   rbmass)

        nodes = generate_xbnode(n_node,&
                                master_node,&
                                vdof,&
                                fdof)


        i_out = 10  ! random unit for output
        ! updating position vectors
        ! pos_def = pos_ini
        ! psi_def = psi_ini
        pos_def_history = 0.0_c_double
        pos_def_history(1, :, :) = pos_def
        psi_def_history = 0.0_c_double
        psi_def_history(1, :, :, :) = psi_def
        pos_dot_def_history = 0.0_c_double
        psi_dot_def_history = 0.0_c_double
        pos_dot_def = 0.0_c_double
        psi_dot_def = 0.0_c_double

        call cbeam3_solv_nlndyn (i_out,&
                                   num_dof,&
                                   Time,&
                                   elements,&
                                   nodes,&
                                   applied_forces,&
                                   dynamic_forces_amplitude,&
                                   dynamic_forces_time,&
                                   forced_vel,&
                                   forced_acc,&
                                   pos_ini,&
                                   psi_ini,&
                                   pos_def,&
                                   psi_def,&
                                   pos_dot_def,&
                                   psi_dot_def,&
                                   Options,&
                                   pos_def_history,&
                                   psi_def_history,&
                                   pos_dot_def_history,&
                                   psi_dot_def_history)

    end subroutine cbeam3_solv_nlndyn_python



    subroutine cbeam3_solv_modal_python(n_elem,&
                                            n_node,&
                                            num_dof,&
                                            num_nodes,&
                                            mem_number,&
                                            conn,&
                                            master,&
                                            n_mass,&
                                            mass_db,&
                                            mass_indices,&
                                            n_stiffness,&
                                            stiffness_db,&
                                            inv_stiffness_db,&
                                            stiffness_indices,&
                                            for_delta,&
                                            rbmass,&
                                            master_node,&
                                            vdof,&
                                            fdof,&
                                            options,&
                                            pos_ini,&
                                            psi_ini,&
                                            pos_def,&
                                            psi_def,&
                                            fullM,&
                                            fullK) bind(C)
        integer(c_int), intent(IN)      :: n_elem
        integer(c_int), intent(IN)      :: n_node
        integer(c_int), intent(IN)      :: num_dof

        ! elem data
        integer(c_int), intent(IN)      :: num_nodes(n_elem)
        integer(c_int), intent(IN)      :: mem_number(n_elem)
        integer(c_int), intent(IN)      :: conn(n_elem, max_elem_node)
        integer(c_int), intent(IN)      :: master(n_elem, max_elem_node, 2)
        integer(c_int), intent(IN)      :: n_mass
        real(c_double), intent(IN)      :: mass_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: mass_indices(n_elem)
        integer(c_int), intent(IN)      :: n_stiffness
        real(c_double), intent(IN)      :: stiffness_db(n_mass, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(n_elem, max_elem_node, 3)
        real(c_double), intent(IN)      :: rbmass(n_elem, max_elem_node, 6, 6)

        ! node data
        integer(c_int), intent(IN)      :: master_node(n_node, 2)
        integer(c_int), intent(IN)      :: vdof(n_node)
        integer(c_int), intent(IN)      :: fdof(n_node)

        type(xbelem)                    :: elements(n_elem)
        type(xbnode)                    :: nodes(n_node)
        type(xbopts), intent(INOUT)     :: options

        real(c_double), intent(IN)      :: pos_ini(n_node, 3)
        real(c_double), intent(IN)      :: psi_ini(n_elem, max_elem_node, 3)
        real(c_double), intent(OUT)     :: pos_def(n_node, 3)
        real(c_double), intent(OUT)     :: psi_def(n_elem, max_elem_node, 3)


        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem
        real(c_double)                  :: Vrel(1, 6)

        type(sparse), allocatable       :: Mglobal(:)
        type(sparse), allocatable       :: Cglobal(:)
        type(sparse), allocatable       :: Kglobal(:)
        integer(c_int)                  :: ms
        integer(c_int)                  :: cs
        integer(c_int)                  :: ks

        real(c_double), intent(OUT)     :: fullM(num_dof, num_dof)
        real(c_double), intent(OUT)     :: fullK(num_dof, num_dof)

        Vrel = 0.0_c_double

        ! gaussian nodes
        nodes_per_elem = count(conn(1,:) /= 0)
        options%NumGauss = nodes_per_elem - 1



        elements = generate_xbelem(n_elem,&
                                   num_nodes,&
                                   mem_number,&
                                   conn,&
                                   master,&
                                   n_mass,&
                                   mass_db,&
                                   mass_indices,&
                                   n_stiffness,&
                                   stiffness_db,&
                                   inv_stiffness_db,&
                                   stiffness_indices,&
                                   for_delta,&
                                   psi_ini,&
                                   rbmass)

        nodes = generate_xbnode(n_node,&
                                master_node,&
                                vdof,&
                                fdof)
       ! overloaded function
       call cbeam3_solv_modal    (num_dof,&
                                  elements,&
                                  nodes,&
                                  Vrel,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  options,&
                                  Mglobal,&
                                  Cglobal,&
                                  Kglobal,&
                                  ms,&
                                  cs,&
                                  ks&
                                  )


      ! inflate Mglobal and Kglobal
      fullM = Mglobal(1)%a
      fullK = Kglobal(1)%a

    end subroutine cbeam3_solv_modal_python



    function generate_xbelem(n_elem,&
                               num_nodes,&
                               mem_number,&
                               conn,&
                               master,&
                               n_mass,&
                               mass_db,&
                               mass_indices,&
                               n_stiffness,&
                               stiffness_db,&
                               inv_stiffness_db,&
                               stiffness_indices,&
                               for_delta,&
                               psi_ini,&
                               RBMass) result(elements)
        ! elem data
        integer(c_int), intent(IN)      :: n_elem
        integer(c_int), intent(IN)      :: num_nodes(n_elem)
        integer(c_int), intent(IN)      :: mem_number(n_elem)
        integer(c_int), intent(IN)      :: conn(n_elem, max_elem_node)
        integer(c_int), intent(IN)      :: master(n_elem, max_elem_node, 2)
        integer(c_int), intent(IN)      :: n_mass
        real(c_double), intent(IN)      :: mass_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: mass_indices(n_elem)
        integer(c_int), intent(IN)      :: n_stiffness
        real(c_double), intent(IN)      :: stiffness_db(n_stiffness, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_stiffness, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(:,:,:)
        real(c_double), intent(IN)      :: psi_ini(:,:,:)
        real(c_double), intent(IN)      :: RBMass(n_elem, max_elem_node, 6, 6)

        type(xbelem)                    :: elements(n_elem)

        integer(c_int)                  :: i
        integer(c_int)                  :: inode_global
        integer(c_int)                  :: inode_local


        do i=1, n_elem
            elements(i)%NumNodes    = num_nodes(i)
            elements(i)%MemNo       = 0*mem_number(i)
            elements(i)%Conn        = conn(i, :)
            elements(i)%Master      = master(i, :, :)
            elements(i)%Length      = 0.0d0  ! dont think it is necessary
            select case (elements(i)%NumNodes)
            case(3)
                inode_local = 3
            case(2)
                inode_local = 1
            case default
                stop 'Not supported numNodes'
            end select
            inode_global = elements(i)%Conn(inode_local)
            elements(i)%Psi         = psi_ini(i, inode_local, :)
            elements(i)%Vector      = for_delta(i, inode_local, :)
            elements(i)%Mass        = mass_db(mass_indices(i), :, :)
            elements(i)%Stiff       = stiffness_db(stiffness_indices(i), :, :)
            elements(i)%InvStiff    = inv_stiffness_db(stiffness_indices(i),:,:)
            elements(i)%RBMass      = RBMass(i, :, :, :)
        end do

    end function generate_xbelem



    subroutine print_xbelem(input)
        type(xbelem), intent(IN)            :: input(:)
        integer                             :: unit, i, ii

        open(newunit=unit, file='debug_elem.txt')

        do i=1, size(input)
            write(unit,*) "-----------------------------------"
            write(unit,*) "NumNodes = ", input(i)%NumNodes
            write(unit,*) "MemNo = ", input(i)%MemNo
            write(unit,*) "Conn = ", input(i)%Conn
            write(unit,*) "Master = ", input(i)%Master
            write(unit,*) "Psi = ", input(i)%Psi
            write(unit,*) "Vector = ", input(i)%Vector
            write(unit,*) "Mass = "
            do ii=1, 6
                write(unit,*) input(i)%Mass(ii,:)
            end do
            write(unit,*) "Stiff = "
            do ii=1, 6
                write(unit,*) input(i)%Stiff(ii,:)
            end do
            write(unit,*) ""
        end do
        close (unit)
    end subroutine print_xbelem


    subroutine print_xbnode(input)
        type(xbnode), intent(IN)            :: input(:)
        integer                             :: unit, i

        open(newunit=unit, file='debug_node.txt')

        do i=1, size(input)
            write(unit, *) "-----------------------------------"
            write(unit, *) "master= ", input(i)%master
            write(unit, *) "vdof= ", input(i)%vdof
            write(unit, *) "fdof= ", input(i)%fdof
            write(unit, *) ""
        end do
        close(unit)
    end subroutine print_xbnode

    function generate_xbnode(n_node,&
                             master,&
                             vdof,&
                             fdof) result(nodes)
        ! node data
        integer(c_int), intent(IN)      :: n_node
        integer(c_int), intent(IN)      :: master(n_node, 2)
        integer(c_int), intent(IN)      :: vdof(n_node)
        integer(c_int), intent(IN)      :: fdof(n_node)

        type(xbnode)                    :: nodes(n_node)

        integer(c_int)                  :: i

        do i=1, n_node
            nodes(i)%Master      = master(i, :)
            nodes(i)%Vdof        = vdof(i)
            nodes(i)%fdof        = fdof(i)
        end do
    end function generate_xbnode


subroutine output_elems (Elem,Coords,Psi)
  use lib_fem

! I/O Variables.
  type(xbelem), intent(in)   :: Elem   (:)        ! Element information.
  real(8),      intent(in)   :: Coords (:,:)      ! Coordinates of the grid points.
  real(8),      intent(in)   :: Psi    (:,:,:)    ! CRV of the nodes in the elements.

! Local variables.
  integer:: i                    ! Counter.
  integer:: iElem                ! Counter on the finite elements.
  integer:: NumE                 ! Number of elements in the model.
  integer:: NumNE                ! Number of nodes in an element.
  real(8):: PosElem (MaxElNod,3) ! Coordinates/CRV of nodes in the element.
  integer:: unit

  open(newunit=unit, file='output_def.txt', status='replace')
  NumE=size(Elem)

! Loop in the elements in the model.
  do iElem=1,NumE
    ! associates the coordinates of the nn-th node (global) to the ii-th node
    ! (local) of the element.
    call fem_glob2loc_extract (Elem(iElem)%Conn,Coords,PosElem,NumNE)

    do i=1,NumNE
      write (unit,'(2I4,1P6E15.6)') iElem,i,PosElem(i,:),Psi(iElem,i,:)
    end do
  end do
  close(unit)

end subroutine output_elems


end module cbeam3_interface
