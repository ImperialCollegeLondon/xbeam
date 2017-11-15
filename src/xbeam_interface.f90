module xbeam_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: xbeam_solv
    use                                 :: debug_utils
    use                                 :: cbeam3_interface

implicit none

    integer(c_int), parameter, private  :: max_elem_node = MaxElNod

contains
subroutine xbeam_solv_couplednlndyn_python(n_elem,&
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
                                        !    pos_def,&
                                        !    psi_def,&
                                           app_forces,&
                                           dynamic_forces,&
                                           for_vel,&
                                           for_acc,&
                                           pos_def_history,&
                                           psi_def_history,&
                                           pos_dot_def_history,&
                                           psi_dot_def_history,&
                                           quat_history,&
                                           success) bind(C)

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
   real(c_double), intent(OUT)     :: quat_history(n_tsteps, 4)

   real(c_double), intent(IN)      :: app_forces(n_node, 6)

   ! dynamic forces
   real(c_double), intent(IN)      :: dynamic_forces(n_node, 6, n_tsteps)
   real(c_double), intent(OUT)     :: for_vel(n_tsteps, 6)
   real(c_double), intent(OUT)     :: for_acc(n_tsteps, 6)

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

   logical(c_bool), intent(INOUT)  :: success

   success = .TRUE.
   quat_history = 0.0_c_double
   quat_history(1, 1) = 1.0_c_double

   ! number of nodes per element
   nodes_per_elem = 0
   nodes_per_elem = count(conn(1, :) /= 0)
   options%NumGauss = nodes_per_elem - 1

   num_dof = count(vdof > 0)*6
   applied_forces = app_forces
   ! do i=1, num_app_forces
   !     applied_forces(node_app_forces(i), :) = app_forces(i, :)
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


   ! call print_matrix('Node', nodes)
   ! call print_matrix('Elem', elements)
   ! call print_matrix('PsiIni', psi_ini)
   ! call print_matrix('PosDef', pos_def)
   ! call print_matrix('PsiDef', psi_def)
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

   call  xbeam_solv_couplednlndyn (i_out,&
                                   num_dof,&
                                   Time,&
                                   elements,&
                                   nodes,&
                                   applied_forces,&
                                   dynamic_forces,&
                                   for_vel,&
                                   for_acc,&
                                   pos_ini,&
                                   psi_ini,&
                                   pos_def,&
                                   psi_def,&
                                   pos_dot_def,&
                                   psi_dot_def,&
                                   pos_def_history,&
                                   psi_def_history,&
                                   pos_dot_def_history,&
                                   psi_dot_def_history,&
                                   quat_history,&
                                   Options,&
                                   success)



end subroutine xbeam_solv_couplednlndyn_python

    subroutine xbeam_solv_nlndyn_step_python   (n_elem,&
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
                                                for_vel,&
                                                for_acc&
                                                )bind(C)

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
        real(c_double), intent(INOUT)   :: quat(4)
        real(c_double), intent(INOUT)   :: for_vel(6)
        real(c_double), intent(INOUT)   :: for_acc(6)

        integer(c_int)                  :: num_dof
        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem

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
        ! print*, 'IN'
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

        ! updating position vectors
        call xbeam_solv_couplednlndyn_step_updated(num_dof,&
                                                   dt,&
                                                   n_node,&
                                                   n_elem,&
                                                   elements,&
                                                   nodes,&
                                                   pos_ini,&
                                                   pos_def,&
                                                   psi_ini,&
                                                   psi_def,&
                                                   pos_def_dot,&
                                                   psi_def_dot,&
                                                   steady_app_forces,&
                                                   dynamic_app_forces,&
                                                   for_vel,&
                                                   for_acc,&
                                                   quat,&
                                                   options)

    end subroutine xbeam_solv_nlndyn_step_python
end module xbeam_interface
