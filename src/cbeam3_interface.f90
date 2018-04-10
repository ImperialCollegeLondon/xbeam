module cbeam3_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: cbeam3_solv
    use                                 :: debug_utils
    use                                 :: lib_sparse
    use                                 :: xbeam_asbly
    use                                 :: lib_xbeam

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
                                            applied_forces,&
                                            gravity_forces&
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

        real(c_double), intent(IN)      :: applied_forces(n_node, 6)
        real(c_double), intent(INOUT)   :: gravity_forces(n_node, 6)
        ! ADC XXX: careful with forces in master FoR

        integer(c_int)                  :: num_dof

        integer(c_int)                  :: nodes_per_elem
        integer(c_int)                  :: i

        integer                         :: unit
        Logical                         :: halt

        num_dof = count(vdof > 0)*6

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
        ! call print_matrix('conn',conn)
        ! call print_matrix('pos_ini', pos_ini)
        ! call print_matrix('pos_def', pos_def)
        ! call print_matrix('psi_ini1', psi_ini(:, 1, :))
        ! call print_matrix('psi_ini2', psi_ini(:, 2, :))
        ! call print_matrix('psi_ini3', psi_ini(:, 3, :))
        ! call print_matrix('psi_def1', psi_def(:, 1, :))
        ! call print_matrix('psi_def2', psi_def(:, 2, :))
        ! call print_matrix('psi_def3', psi_def(:, 3, :))
        ! call print_matrix('fdof',fdof)
        ! call print_matrix('app_forces',applied_forces)
        ! ! call print_matrix('gravity_forces',gravity_forces)
        ! call print_matrix('vdof',vdof)
        ! call print_xbelem(elements)
        ! call print_xbnode(nodes)
        ! print*, 'RBMass:'
        ! print*, elements(1)%RBMASS(1, 1, :)
        ! print*, elements(1)%RBMASS(1, 2, :)
        ! print*, elements(1)%RBMASS(1, 3, :)
        ! print*, elements(1)%RBMASS(1, 4, :)
        ! print*, elements(1)%RBMASS(1, 5, :)
        ! print*, elements(1)%RBMASS(1, 6, :)

       ! overloaded function
       gravity_forces = 0.0d0
       call cbeam3_solv_nlnstatic(num_dof,&
                                  n_elem,&
                                  n_node, &
                                  elements,&
                                  nodes,&
                                  applied_forces,&
                                  gravity_forces,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  options&
                                  )
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


    subroutine cbeam3_solv_nlndyn_step_python   (n_elem,&
                                                 n_node,&
                                                 dt,&
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
                                                 pos_def_dot,&
                                                 psi_def,&
                                                 psi_def_dot,&
                                                 steady_app_forces,&
                                                 dynamic_app_forces,&
                                                 quat,&
                                                 forced_vel,&
                                                 forced_acc)bind(C)

        integer(c_int), intent(IN)      :: n_elem
        integer(c_int), intent(IN)      :: n_node
        real(c_double), intent(IN)      :: dt

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
        real(c_double), intent(INOUT)   :: pos_def(n_node, 3)
        real(c_double), intent(INOUT)   :: psi_def(n_elem, max_elem_node, 3)
        real(c_double), intent(INOUT)   :: pos_def_dot(n_node, 3)
        real(c_double), intent(INOUT)   :: psi_def_dot(n_elem, max_elem_node, 3)

        real(c_double), intent(IN)      :: steady_app_forces (n_node, 6)
        ! ADC: careful, forces in master FoR

        ! dynamic
        real(c_double), intent(IN)      :: dynamic_app_forces(n_node, 6)
        real(c_double), intent(IN)      :: quat(4)
        real(c_double), intent(IN)      :: forced_vel(6)
        real(c_double), intent(IN)      :: forced_acc(6)

        integer(c_int)                  :: num_dof

        integer(c_int)                  :: i_out
        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem

        integer(c_int), parameter        :: n_tsteps = 1

        ! aux variables
        real(c_double)                  :: time(n_tsteps + 1)
        real(c_double)                  :: forced_vel_mat(2, 6)
        real(c_double)                  :: forced_acc_mat(2, 6)


        forced_vel_mat(1,:) = forced_vel
        forced_acc_mat(1,:) = forced_acc
        forced_vel_mat(2,:) = forced_vel
        forced_acc_mat(2,:) = forced_acc
        ! print*, 'cbeam3_interface, line 417'
        ! print*, quat
        ! fa_fake = 0
        ! qquat = quat(:,1)
        ! print*, 'cbeam3_interface, line 417'
        ! ftime_fake = 0
        ! print*, 'cbeam3_interface, line 417'
        ! print*, size(forced_vel)
        ! print*, (forced_vel)
        ! print*, 'cbeam3_interface, line 417'
        !
        !
        ! print*, 'cbeam3_interface, line 423'
        ! print*,  n_node
        ! call print_matrix('conn',conn)
        ! call print_matrix('pos_ini', pos_ini)
        ! call print_matrix('pos_def', pos_def)
        ! call print_matrix('pos_def_dot', pos_def_dot)
        ! call print_matrix('psi_ini1', psi_ini(:, 1, :))
        ! call print_matrix('psi_ini2', psi_ini(:, 2, :))
        ! call print_matrix('psi_ini3', psi_ini(:, 3, :))
        ! call print_matrix('psi_def1', psi_def(:, 1, :))
        ! call print_matrix('psi_def2', psi_def(:, 2, :))
        ! call print_matrix('psi_def3', psi_def(:, 3, :))
        ! call print_matrix('fdof',fdof)
        ! call print_matrix('steady_app_forces',steady_app_forces)
        ! call print_matrix('dynamic_app_forces',dynamic_app_forces)
        ! call print_matrix('vdof',vdof)
        ! call print_matrix('master1',master(:,:,1))
        ! call print_matrix('master2',master(:,:,2))
        ! call print_matrix('masternode',master_node)
        ! call print_matrix('psi_def_dot1',psi_def_dot(:, 1, :))
        ! call print_matrix('psi_def_dot2',psi_def_dot(:, 2, :))
        ! call print_matrix('forced_vel', forced_vel)
        ! call print_matrix('forced_acc', forced_acc)
        num_dof = count(vdof > 0)*6

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


        ! call print_xbelem(elements)
        ! call print_xbnode(nodes)
        i_out = 10  ! random unit for output
        time(1) = 0
        time(2) = dt

        ! updating position vectors
        call cbeam3_solv_nlndyn_step (i_out,&
                                      num_dof,&
                                      time,&
                                      elements,&
                                      nodes,&
                                      steady_app_forces+dynamic_app_forces,&
                                      steady_app_forces*0.0_c_double,&
                                      time*0.0_c_double,&
                                      forced_vel_mat,&
                                      forced_acc_mat,&
                                      pos_ini,&
                                      psi_ini,&
                                      pos_def,&
                                      psi_def,&
                                      pos_def_dot,&
                                      psi_def_dot,&
                                      Options,&
                                      quat)

    end subroutine cbeam3_solv_nlndyn_step_python

    ! subroutine cbeam3_solv_modal_python(n_elem,&
    !                                         n_node,&
    !                                         num_dof,&
    !                                         num_nodes,&
    !                                         mem_number,&
    !                                         conn,&
    !                                         master,&
    !                                         n_mass,&
    !                                         mass_db,&
    !                                         mass_indices,&
    !                                         n_stiffness,&
    !                                         stiffness_db,&
    !                                         inv_stiffness_db,&
    !                                         stiffness_indices,&
    !                                         for_delta,&
    !                                         rbmass,&
    !                                         master_node,&
    !                                         vdof,&
    !                                         fdof,&
    !                                         options,&
    !                                         pos_ini,&
    !                                         psi_ini,&
    !                                         pos_def,&
    !                                         psi_def,&
    !                                         fullM,&
    !                                         fullK) bind(C)
    !     integer(c_int), intent(IN)      :: n_elem
    !     integer(c_int), intent(IN)      :: n_node
    !     integer(c_int), intent(IN)      :: num_dof
    !
    !     ! elem data
    !     integer(c_int), intent(IN)      :: num_nodes(n_elem)
    !     integer(c_int), intent(IN)      :: mem_number(n_elem)
    !     integer(c_int), intent(IN)      :: conn(n_elem, max_elem_node)
    !     integer(c_int), intent(IN)      :: master(n_elem, max_elem_node, 2)
    !     integer(c_int), intent(IN)      :: n_mass
    !     real(c_double), intent(IN)      :: mass_db(n_mass, 6, 6)
    !     integer(c_int), intent(IN)      :: mass_indices(n_elem)
    !     integer(c_int), intent(IN)      :: n_stiffness
    !     real(c_double), intent(IN)      :: stiffness_db(n_mass, 6, 6)
    !     real(c_double), intent(IN)      :: inv_stiffness_db(n_mass, 6, 6)
    !     integer(c_int), intent(IN)      :: stiffness_indices(n_elem)
    !     real(c_double), intent(IN)      :: for_delta(n_elem, max_elem_node, 3)
    !     real(c_double), intent(IN)      :: rbmass(n_elem, max_elem_node, 6, 6)
    !
    !     ! node data
    !     integer(c_int), intent(IN)      :: master_node(n_node, 2)
    !     integer(c_int), intent(IN)      :: vdof(n_node)
    !     integer(c_int), intent(IN)      :: fdof(n_node)
    !
    !     type(xbelem)                    :: elements(n_elem)
    !     type(xbnode)                    :: nodes(n_node)
    !     type(xbopts), intent(INOUT)     :: options
    !
    !     real(c_double), intent(IN)      :: pos_ini(n_node, 3)
    !     real(c_double), intent(IN)      :: psi_ini(n_elem, max_elem_node, 3)
    !     real(c_double), intent(OUT)     :: pos_def(n_node, 3)
    !     real(c_double), intent(OUT)     :: psi_def(n_elem, max_elem_node, 3)
    !
    !
    !     integer(c_int)                  :: i
    !     integer(c_int)                  :: nodes_per_elem
    !     real(c_double)                  :: Vrel(1, 6)
    !
    !     type(sparse), allocatable       :: Mglobal(:)
    !     type(sparse), allocatable       :: Cglobal(:)
    !     type(sparse), allocatable       :: Kglobal(:)
    !     integer(c_int)                  :: ms
    !     integer(c_int)                  :: cs
    !     integer(c_int)                  :: ks
    !
    !     real(c_double), intent(OUT)     :: fullM(num_dof, num_dof)
    !     real(c_double), intent(OUT)     :: fullK(num_dof, num_dof)
    !
    !     Vrel = 0.0_c_double
    !
    !     ! gaussian nodes
    !     nodes_per_elem = count(conn(1,:) /= 0)
    !     options%NumGauss = nodes_per_elem - 1
    !
    !
    !
    !     elements = generate_xbelem(n_elem,&
    !                                num_nodes,&
    !                                mem_number,&
    !                                conn,&
    !                                master,&
    !                                n_mass,&
    !                                mass_db,&
    !                                mass_indices,&
    !                                n_stiffness,&
    !                                stiffness_db,&
    !                                inv_stiffness_db,&
    !                                stiffness_indices,&
    !                                for_delta,&
    !                                psi_ini,&
    !                                rbmass)
    !
    !     nodes = generate_xbnode(n_node,&
    !                             master_node,&
    !                             vdof,&
    !                             fdof)
    !    ! overloaded function
    !    call cbeam3_solv_modal    (num_dof,&
    !                               elements,&
    !                               nodes,&
    !                               Vrel,&
    !                               pos_ini,&
    !                               psi_ini,&
    !                               pos_def,&
    !                               psi_def,&
    !                               options,&
    !                               Mglobal,&
    !                               Cglobal,&
    !                               Kglobal,&
    !                               ms,&
    !                               cs,&
    !                               ks&
    !                               )
    !
    !
    !   ! inflate Mglobal and Kglobal
    !   fullM = Mglobal(1)%a
    !   fullK = Kglobal(1)%a
    !
    ! end subroutine cbeam3_solv_modal_python



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
        real(c_double), intent(IN)      :: for_delta(n_elem, max_elem_node, 3)
        real(c_double), intent(IN)      :: psi_ini(n_elem,3,3)
        real(c_double), intent(IN)      :: RBMass(n_elem, max_elem_node, 6, 6)

        type(xbelem)                    :: elements(n_elem)

        integer(c_int)                  :: i
        integer(c_int)                  :: j
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
            ! do j=1, 3
            !     if (any(elements(i)%RBMass(j, :, :) /= 0.0)) then
            !         print*, 'RBMASS'
            !         print*, i, j
            !         print*, elements(i)%RBMass(j, 1, :)
            !         print*, elements(i)%RBMass(j, 2, :)
            !         print*, elements(i)%RBMass(j, 3, :)
            !         print*, elements(i)%RBMass(j, 4, :)
            !         print*, elements(i)%RBMass(j, 5, :)
            !         print*, elements(i)%RBMass(j, 6, :)
            !         print*, '-----'
            !     end if
            ! end do
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
            nodes(i)%Sflag       = 0
        end do
    end function generate_xbnode

    subroutine cbeam3_solv_state2disp_python (&
                                              n_elem,&
                                              n_node,&
                                              numdof,&
                                              pos_ini,&
                                              psi_ini,&
                                              pos_def,&
                                              psi_def,&
                                              pos_dot,&
                                              psi_dot,&
                                              master_node,&
                                              vdof,&
                                              num_nodes,&
                                              master_elem,&
                                              q,&
                                              dqdt) bind(C)
    use lib_fem
    use lib_rotvect
    use lib_bgeom
    use lib_cbeam3
    use iso_c_binding

    ! I/O Variables.
    integer(c_int), intent(IN)              :: n_elem
    integer(c_int), intent(IN)              :: n_node
    integer(c_int), intent(IN)              :: numdof
    real(c_double), intent(IN)              :: pos_ini(n_node, 3)
    real(c_double), intent(IN)              :: psi_ini(n_elem, 3, 3)
    real(c_double), intent(OUT)             :: pos_def(n_node, 3)
    real(c_double), intent(OUT)             :: psi_def(n_elem, 3, 3)
    real(c_double), intent(OUT)             :: pos_dot(n_node, 3)
    real(c_double), intent(OUT)             :: psi_dot(n_elem, 3, 3)
    integer(c_int), intent(IN)              :: master_node(n_node, 2)
    integer(c_int), intent(IN)              :: vdof(n_node)
    integer(c_int), intent(IN)              :: num_nodes(n_elem)
    integer(c_int), intent(IN)              :: master_elem(n_elem, 3, 2)
    real(c_double), intent(IN)              :: q(numdof + 10)
    real(c_double), intent(IN)              :: dqdt(numdof + 10)

    ! Local variables.
    integer :: i,j,k               ! Counters.
    integer :: ix                  ! Counter on the degrees of freedom.
    integer :: iElem               ! Counter on the elements.
    integer :: iNode               ! Counter on the nodes.

    ! Store current displacement and its time derivative at all nodes and the
    ! rotations and its first derivatives at the master nodes of each element.
    ix = 0
    do iNode = 1, n_node
        iElem = master_node(iNode, 1)
        k = master_node(iNode, 2)
            ! Constrained nodes.
            if (Vdof(iNode) == 0) then
              pos_def(iNode,:) = pos_ini(iNode,:)
              pos_dot(iNode,:) = 0.d0

              psi_def(iElem,k,:) = psi_ini(iElem, k, :)
              psi_dot(iElem,k,:) = 0.d0

            ! Unconstrained nodes.
            else
              pos_def(iNode, :) = q   (ix+1:ix+3)
              pos_dot(iNode, :) = dqdt(ix+1:ix+3)

              psi_def(iElem,k,:) = q   (ix+4:ix+6)
              psi_dot(iElem,k,:) = dqdt(ix+4:ix+6)

              ix=ix+6
            end if
        end do

        ! Compute rotation vector and time derivative at slave nodes within elements.
        do i=1, n_elem
            do j=1, num_nodes(i)
            ! Copy rotation and derivative from master node for each slave node.
              if (master_elem(i, j ,1) /= 0) then
                psi_def(i,j,:) = psi_def(master_elem(i, j, 1), master_elem(i, j, 2), :)
                psi_dot(i,j,:) = psi_dot(master_elem(i, j, 1), master_elem(i, j, 2), :)
              end if
            end do
        ! Include master-to-slave initial rotations from the undeformed configuration.
        call cbeam3_projm2s (num_nodes(i), master_elem(i, :, :), psi_ini(i,:,:), psi_ini, psi_def(i,:,:))
    end do

end subroutine cbeam3_solv_state2disp_python

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_DISP2STATE_PYTHON
!
!-> Description:
!
!    Write current state vector from current displacements and rotations.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_disp2state_python (&
                                           n_elem,&
                                           n_node,&
                                           numdof,&
                                           pos_def,&
                                           psi_def,&
                                           pos_dot,&
                                           psi_dot,&
                                           vdof,&
                                           master_node,&
                                           q,&
                                           dqdt) bind(C)
  use lib_fem

! I/O Variables.
    integer(c_int), intent(IN)              :: n_elem
    integer(c_int), intent(IN)              :: n_node
    integer(c_int), intent(IN)              :: numdof
    real(c_double), intent(IN)              :: pos_def(n_node, 3)
    real(c_double), intent(IN)              :: psi_def(n_elem, 3, 3)
    real(c_double), intent(IN)              :: pos_dot(n_node, 3)
    real(c_double), intent(IN)              :: psi_dot(n_elem, 3, 3)
    integer(c_int), intent(IN)              :: vdof(n_node)
    integer(c_int), intent(IN)              :: master_node(n_node, 2)
    real(c_double), intent(INOUT)             :: q(numdof + 10)
    real(c_double), intent(INOUT)             :: dqdt(numdof + 10)


! Local variables.
  integer :: k                   ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Loop in all nodes in the model.
  ix=0
  do iNode=1, n_node
    iElem=master_node(iNode, 1)
    if (Vdof(iNode) /= 0) then
      k=master_node(iNode, 2)

! Current nodal displacements and derivatives.
      q   (ix+1:ix+3)= pos_def(iNode,:)
      dqdt(ix+1:ix+3)= pos_dot(iNode,:)

! Cartesian rotation vector at master nodes.
      q   (ix+4:ix+6)=psi_def(iElem,k,:)
      dqdt(ix+4:ix+6)=psi_dot(iElem,k,:)
      ix=ix+6
    end if
  end do

  return
 end subroutine cbeam3_solv_disp2state_python



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
