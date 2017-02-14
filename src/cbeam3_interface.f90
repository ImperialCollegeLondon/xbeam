module cbeam3_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: cbeam3_solv
    use                                 :: debug_utils
    use                                 :: interface_lapack
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
                                            num_app_forces,&
                                            app_forces,&
                                            node_app_forces) bind(C)
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
        real(c_double), intent(IN)      :: stiffness_db(n_mass, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(n_node, 3)
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

        integer(c_int), intent(IN)      :: num_app_forces
        real(c_double), intent(IN)      :: app_forces(num_app_forces, 6)
        ! ADC XXX: careful with forces in master FoR
        integer(c_int), intent(IN)      :: node_app_forces(num_app_forces)

        integer(c_int)                  :: num_dof
        real(c_double)                  :: applied_forces(n_node, 6)!legacy var

        integer(c_int)                  :: i
        integer(c_int)                  :: j
        integer(c_int)                  :: nodes_per_elem

        num_dof = count(vdof > 0)*6
        applied_forces = 0.0_c_double
        do i=1, num_app_forces
            applied_forces(node_app_forces(i), :) = app_forces(i, :)
        end do

        print*, options%gravity_dir_x
        print*, options%gravity_dir_y
        print*, options%gravity_dir_z

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

    do i=1, n_elem
        !print*, 'i=', i
        do j=1,3
            !print*, 'j=', j
            !print*, psi_ini(i, j, :)
        end do
    end do

        nodes = generate_xbnode(n_node,&
                                master_node,&
                                vdof,&
                                fdof)

    do i=1, n_node
        !do j=1, 3
            !print*, nodes(i)%Vdof
            !print*, nodes(i)%Fdof
        !end do
        !print*, '--'
    end do

        call print_matrix('Node', nodes)
        call print_matrix('Elem', elements)
        call print_matrix('PsiIni', psi_ini)
        call print_matrix('PosDef', pos_def)
        call print_matrix('PsiDef', psi_def)
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

    
        call print_matrix('PosDef2', pos_def)
        call print_matrix('PsiDef2', psi_def)
        print*, 'Leaving fortran'

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
                                            num_app_forces,&
                                            app_forces,&
                                            node_app_forces,&
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
        real(c_double), intent(IN)      :: for_delta(n_node, 3)
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

        integer(c_int), intent(IN)      :: num_app_forces
        real(c_double), intent(IN)      :: app_forces(num_app_forces, 6)
        ! ADC: careful, forces in master FoR
        integer(c_int), intent(IN)      :: node_app_forces(num_app_forces)

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
        integer(c_int)                  :: j
        integer(c_int)                  :: nodes_per_elem


        num_dof = count(vdof > 0)*6
        applied_forces = 0.0_c_double
        do i=1, num_app_forces
            applied_forces(node_app_forces(i), :) = app_forces(i, :)
        end do

        ! gaussian nodes
        nodes_per_elem = count(conn(1,:) /= 0)
        options%NumGauss = nodes_per_elem - 1

        call print_matrix('forced_vel', forced_vel)
        call print_matrix('dynamic_forces_amplitude',dynamic_forces_amplitude)
        call print_matrix('time', time)

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


        call print_matrix('Node', nodes)
        call print_matrix('Elem', elements)
        call print_matrix('PsiIni', psi_ini)
        call print_matrix('PosDef', pos_def)
        call print_matrix('PsiDef', psi_def)
        i_out = 10  ! random unit for output
        ! updating position vectors
        pos_def = pos_ini
        psi_def = psi_ini
        pos_def_history = 0.0_c_double
        pos_def_history(1, :, :) = pos_def
        psi_def_history = 0.0_c_double
        psi_def_history(1, :, :, :) = psi_def
        pos_dot_def_history = 0.0_c_double
        psi_dot_def_history = 0.0_c_double
        pos_dot_def = 0.0_c_double
        psi_dot_def = 0.0_c_double

        ! Signature for cbeam3_solv_nlndyn
          !integer,      intent(in)   :: iOut              ! Output file.
          !integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
          !real(8),      intent(in)   :: Time      (:)     ! Time steps.
          !type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
          !type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
          !real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
          !real(8),      intent(in)   :: Fa        (:,:)  ! Amplitude of the dynamic nodal forces.
          !real(8),      intent(in)   :: Ftime     (:)    ! Time history of the applied forces.
          !real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
          !real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
          !real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
          !real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
          !real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
          !real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
          !real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
          !real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
          !type(xbopts),intent(in)    :: Options           ! Solver parameters.

          !real(8), intent(OUT)       :: pos_def_history(size(Time) + 1, size(Node), 3)
          !real(8), intent(OUT)       :: pos_dot_def_history(size(Time) + 1, size(Node), 3)
          !real(8), intent(OUT)       :: psi_def_history(size(Time) + 1, size(Elem), MaxElNod, 3)
          !real(8), intent(OUT)       :: psi_dot_def_history(size(Time) + 1, size(Elem), MaxElNod, 3)

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
        real(c_double), intent(IN)      :: for_delta(n_node, 3)
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
        integer(c_int)                  :: j
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
      fullM = 0.0_c_double
      do i=1, ms
          fullM(Mglobal(i)%i, Mglobal(i)%j) = Mglobal(i)%a
      end do
      fullK = 0.0_c_double
      do i=1, ks
          fullK(Kglobal(i)%i, Kglobal(i)%j) = Kglobal(i)%a
      end do


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
        real(c_double), intent(IN)      :: stiffness_db(n_mass, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_mass, 6, 6)
        integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
        real(c_double), intent(IN)      :: for_delta(:,:)
        real(c_double), intent(IN)      :: psi_ini(:,:,:)
        real(c_double), intent(IN)      :: RBMass(n_elem, max_elem_node, 6, 6)

        type(xbelem)                    :: elements(n_elem)

        integer(c_int)                  :: i
        integer(c_int)                  :: inode_global
        integer(c_int)                  :: inode_local


        do i=1, n_elem
            elements(i)%NumNodes    = num_nodes(i)
            elements(i)%MemNo       = mem_number(i)
            elements(i)%Conn        = conn(i, :)
            elements(i)%Master      = master(i, :, :)
            !elements(i)%Length      = 0.5d0 ! TODO quick and ugly solyution
            elements(i)%Length      = 0.0d0  ! dont think it is necessary
            select case (elements(i)%NumNodes)
            case(3)
                inode_local = 3
            case(2)
                inode_local = 1
            end select
            inode_global = elements(i)%Conn(inode_local)
            elements(i)%Psi         = psi_ini(i, inode_local, :)
            elements(i)%Vector      = for_delta(inode_global,:)
            elements(i)%Mass        = mass_db(mass_indices(i), :, :)
            elements(i)%Stiff       = stiffness_db(stiffness_indices(i), :, :)
            elements(i)%InvStiff    = inv_stiffness_db(stiffness_indices(i),:,:)
            elements(i)%RBMass      = RBMass(i, :, :, :)
        end do

    end function generate_xbelem



    subroutine print_xbelem(input)
        type(xbelem), intent(IN)            :: input

        print*, "-----------------------------------"
        print*, "NumNodes = ", input%NumNodes
        print*, "MemNo = ", input%MemNo
        print*, "Conn = ", input%Conn
        print*, "Master = ", input%Master
        print*, "Psi = ", input%Psi
        print*, "Vector = ", input%Vector
        print*, "Mass = ", input%Mass
        print*, "Stiff = ", input%Stiff
        print*, ""
    end subroutine print_xbelem



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


end module cbeam3_interface
