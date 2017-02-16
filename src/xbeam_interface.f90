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
                                           pos_def,&
                                           psi_def,&
                                           num_app_forces,&
                                           app_forces,&
                                           node_app_forces,&
                                           dynamic_forces_amplitude,&
                                           dynamic_forces_time,&
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
   real(c_double), intent(OUT)     :: quat_history(n_tsteps, 4)

   integer(c_int), intent(IN)      :: num_app_forces
   real(c_double), intent(IN)      :: app_forces(num_app_forces, 6)
   ! ADC: careful, forces in master FoR
   integer(c_int), intent(IN)      :: node_app_forces(num_app_forces)

   ! dynamic forces
   real(c_double), intent(IN)      :: dynamic_forces_amplitude(n_node, 6)
   real(c_double), intent(IN)      :: dynamic_forces_time(n_tsteps)
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
   integer(c_int)                  :: j
   integer(c_int)                  :: nodes_per_elem

   logical(c_bool), intent(INOUT)  :: success

   success = .TRUE.
   quat_history = 0.0_c_double
   quat_history(1, 1) = 1.0_c_double

   ! number of nodes per element
   nodes_per_elem = 0
   nodes_per_elem = count(conn(1, :) /= 0)
   print*, nodes_per_elem
   options%NumGauss = nodes_per_elem - 1

   num_dof = count(vdof > 0)*6
   applied_forces = 0.0_c_double
   do i=1, num_app_forces
       applied_forces(node_app_forces(i), :) = app_forces(i, :)
   end do

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


   !call print_matrix('Node', nodes)
   !call print_matrix('Elem', elements)
   !call print_matrix('PsiIni', psi_ini)
   !call print_matrix('PosDef', pos_def)
   !call print_matrix('PsiDef', psi_def)
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
                                   dynamic_forces_amplitude,&
                                   dynamic_forces_time,&
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

end module xbeam_interface
