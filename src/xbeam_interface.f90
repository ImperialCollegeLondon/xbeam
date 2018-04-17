module xbeam_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: xbeam_solv
    use                                 :: debug_utils
    use                                 :: cbeam3_interface

implicit none

    integer(c_int), parameter, private  :: max_elem_node = MaxElNod
    real(8),private,parameter,dimension(4,4):: Unit4= &       ! 4x4 Unit matrix.
& reshape((/1.d0,0.d0,0.d0,0.d0,&
            0.d0,1.d0,0.d0,0.d0,&
            0.d0,0.d0,1.d0,0.d0,&
            0.d0,0.d0,0.d0,1.d0/),(/4,4/))


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
   ! quat_history = 0.0_c_double
   ! quat_history(1, 1) = 1.0_c_double

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

    subroutine xbeam_solv_nlndyn_step_python   (numdof,&
                                                iter,&
                                                n_elem,&
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
                                                gravity_forces,&
                                                quat,&
                                                for_vel,&
                                                for_acc,&
                                                q,&
                                                dqdt,&
                                                dqddt&
                                                )bind(C)

        integer(c_int), intent(IN)      :: numdof
        integer(c_int), intent(IN)      :: iter
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
        real(c_double), intent(IN)      :: stiffness_db(n_stiffness, 6, 6)
        real(c_double), intent(IN)      :: inv_stiffness_db(n_stiffness, 6, 6)
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
        real(c_double), intent(INOUT)   :: gravity_forces(n_node, 6)
        real(c_double), intent(INOUT)   :: quat(4)
        real(c_double), intent(INOUT)   :: for_vel(6)
        real(c_double), intent(INOUT)   :: for_acc(6)
        real(c_double), intent(INOUT)   :: q(numdof + 10)
        real(c_double), intent(INOUT)   :: dqdt(numdof + 10)
        real(c_double), intent(INOUT)   :: dqddt(numdof + 10)

        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem


    ! if (any(isnan(steady_app_forces))) then
    !     print*, 'xbeamint, 291'
    !     stop
    ! end if
    !
    ! if (any(isnan(dynamic_app_forces))) then
    !     print*, 'xbeamint, 296'
    !     stop
    ! end if
        ! print*, quat
        ! fa_fake = 0
        ! qquat = quat(:,1)
        ! print*, 'cbeam3_interface, line 417'
        ! ftime_fake = 0
        ! print*, 'cbeam3_interface, line 417'
        ! print*, size(forced_vel)l
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
        ! call print_matrix('steady_app_forces',steady_app_forces + dynamic_app_forces)
        ! call print_matrix('vdof',vdof)
        ! call print_matrix('master1',master(:,:,1))
        ! call print_matrix('master2',master(:,:,2))
        ! call print_matrix('masternode',master_node)
        ! call print_matrix('psi_def_dot1',psi_def_dot(:, 1, :))
        ! call print_matrix('psi_def_dot2',psi_def_dot(:, 2, :))
        ! call print_matrix('q',q)
        ! call print_matrix('dqdt',dqdt)
        ! call print_matrix('dqddt',dqddt)
        ! call print_matrix('forced_vel', forced_vel)
        ! call print_matrix('forced_acc', forced_acc)
        ! print*, 'IN'

        ! print*, 'for_vel, init: ', for_vel(1:3)

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
        ! call print_elem('elem_before_step', elements)
        ! call print_node('node_before_step', nodes)
        ! options%NumGauss = nodes_per_elem - 1

        gravity_forces = 0.0d0
        call xbeam_solv_couplednlndyn_step_updated(numdof,&
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
                                                   gravity_forces,&
                                                   for_vel,&
                                                   for_acc,&
                                                   quat,&
                                                   q,&
                                                   dqdt,&
                                                   dqddt,&
                                                   options)

        ! print*, 'grav'
        ! print*, 'XBEAM'
        ! print*, gravity_forces(1, :)
        call correct_gravity_forces(n_node, n_elem, gravity_forces, psi_def, elements, nodes)
        ! print*, gravity_forces(1, :)
        ! print*, '--'
        ! print*, 'for_vel, finish: ', for_vel(1:3)
        ! call print_matrix('q', q)
        ! call print_matrix('dqdt', dqdt)
        ! call print_matrix('dqddt', dqddt)
    end subroutine xbeam_solv_nlndyn_step_python



    subroutine xbeam_solv_nlndyn_init_python   (numdof,&
                                                iter,&
                                                n_elem,&
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
                                                for_acc,&
                                                q,&
                                                dqdt,&
                                                dqddt&
                                                )bind(C)

        use lib_fem
        use lib_rot
        use lib_rotvect
        use lib_mat
        use lib_xbeam
        use lib_lu
        use cbeam3_asbly
        use cbeam3_solv
        use xbeam_asbly
        use xbeam_solv
        use iso_c_binding
        integer(c_int), intent(IN)      :: numdof
        integer(c_int), intent(IN)      :: iter
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
        real(c_double), intent(IN)      :: stiffness_db(n_stiffness, 6, 6)
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

        ! dynamic
        real(c_double), intent(IN)      :: dynamic_app_forces(n_node, 6)
        real(c_double), intent(INOUT)   :: quat(4)
        real(c_double), intent(INOUT)   :: for_vel(6)
        real(c_double), intent(INOUT)   :: for_acc(6)
        real(c_double), intent(OUT)     :: q(numdof + 10)
        real(c_double), intent(OUT)     :: dqdt(numdof + 10)
        real(c_double), intent(OUT)     :: dqddt(numdof + 10)

        integer(c_int)                  :: i
        integer(c_int)                  :: nodes_per_elem
        real(8)                         :: CSS(numdof, numdof)     ! Sparse damping matrix.
        real(8)                         :: KSS(numdof, numdof)     ! Elast stiffness matrix in sparse storage.
        real(8)                         :: MSS(numdof, numdof)     ! Elast mass matrix in sparse storage.
        real(8)                         :: Felast(numdof, numdof)  ! Applied external forces on structure
        real(8)                         :: Qelast(numdof)  ! Elast vector of discrete generalize forces.
        real(8)                         :: MSR(numdof, 6)   ! Mass and damping from the motions of reference system.
        real(8)                         :: CSR(numdof, 6)   ! Mass and damping from the motions of reference system.

        ! Define variables for rigid system matrices.
        real(8)                         :: CRS(6, numdof)      ! rigid Sparse damping matrix.
        real(8)                         :: KRS(6, numdof)      ! rigid stiffness matrix in sparse storage.
        real(8)                         :: MRS(6, numdof)   ! Mass and damping from the motions of reference system.
        real(8)                         :: MRS_gravity(6, numdof + 6)   ! Mass and damping from the motions of reference system.
        real(8)                         :: MSS_gravity(numdof + 6, numdof + 6)   ! Mass and damping from the motions of reference system.
        real(8)                         :: MRR_gravity(6, 6)   ! Mass and damping from the motions of reference system.
        real(8)                         :: Frigid(6, numdof + 6)   ! rigid matrix of applied forces in sparse format
        real(8)                         :: Qrigid(6)   ! rigid vector of discrete generalize forces.
        real(8)                         :: MRR(6, 6)    ! rigid Mass and damping from the motions of reference system.
        real(8)                         :: CRR(6, 6)    ! rigid Mass and damping from the motions of reference system.
        real(8)                         :: CQR(4, 6)
        real(8)                         :: CQQ(4, 4)  ! Tangent matrices from linearisation of quaternion equation.

        ! Define variables for rigid-body motion
        real(8)                         :: Cao (3, 3)
                    ! Rotation operator from reference to inertial frame

        ! Define variables for complete system matrices.
        real(8)                         :: Ctotal(numdof + 10, numdof + 10)    ! Total Sparse damping matrix.
        real(8)                         :: Ktotal(numdof + 10, numdof + 10)    ! Total stiffness matrix in sparse storage.
        real(8)                         :: Mtotal(numdof + 10, numdof + 10)    ! Total mass matrix in sparse storage.
        real(8)                         :: Asys(numdof + 10, numdof + 10)    ! System matrix in sparse storage.
        real(8)                         :: Qtotal(numdof + 10)    ! Total vector of discrete generalize forces.
        real(8)                         :: DQ(numdof + 10)
        real(8)                         :: beta
        real(8)                         :: gamma
        integer(c_int)                  :: listin(n_node)
        real(8)                         :: mindelta
        integer(c_int)                  :: i_iter
        logical                         :: converged

        ! print*, 'Init, ', 510
        ! print*, 'numdof', numdof
        ! print*, 'size of q: ', size(q)
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
        ! call print_matrix('steady_app_forces',steady_app_forces + dynamic_app_forces)
        ! call print_matrix('vdof',vdof)
        ! call print_matrix('master1',master(:,:,1))
        ! call print_matrix('master2',master(:,:,2))
        ! call print_matrix('masternode',master_node)
        ! call print_matrix('psi_def_dot1',psi_def_dot(:, 1, :))
        ! call print_matrix('psi_def_dot2',psi_def_dot(:, 2, :))

        Qelast = 0.0d0
        Qrigid = 0.0d0
        Qtotal = 0.0d0

        Frigid = 0.0d0
        Felast = 0.0d0

        MSS = 0.0d0
        MSR = 0.0d0
        MRS = 0.0d0
        MRR = 0.0d0

        CSS = 0.0d0
        CSR = 0.0d0
        CRS = 0.0d0
        CRR = 0.0d0
        CQR = 0.0d0
        CQQ = 0.0d0

        KSS = 0.0d0
        KRS = 0.0d0

        Mtotal = 0.0d0
        Ctotal = 0.0d0
        Ktotal = 0.0d0
        Asys = 0.0d0

        Q = 0.0d0
        dQdt = 0.0d0
        dQddt = 0.0d0

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

        ! ListIN = vdof
        ListIN = 0
        do i=1,n_node
            ListIN(i)=Nodes(i)%Vdof
        end do

        gamma = 0.5d0 + Options%NewmarkDamp
        beta = 0.25d0*(gamma + 0.5d0)*(gamma + 0.5d0)

        call cbeam3_solv_disp2state(nodes,&
                                    pos_def,&
                                    psi_def,&
                                    0.0d0*pos_def,&
                                    0.0d0*psi_def,&
                                    Q(1:numdof),&
                                    dQdt(1:numdof))
        dQdt(numdof+7:numdof+10) = quat
        Cao = xbeam_rot(dqdt(numdof+7:numdof+10))

        call cbeam3_asbly_dynamic(numdof,&
                                  n_node,&
                                  n_elem,&
                                  elements,&
                                  nodes,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  0.0d0*pos_def,&
                                  0.0d0*psi_def,&
                                  0.0d0*pos_def,&
                                  0.0d0*psi_def,&
                                  steady_app_forces + dynamic_app_forces,&
                                  dQdt(numdof+1:numdof+6),&
                                  0.0d0*dQddt(numdof+1:numdof+6),&
                                  MSS,&
                                  MSR,&
                                  CSS,&
                                  CSR,&
                                  KSS,&
                                  Felast,&
                                  Qelast,&
                                  options,&
                                  Cao)

        Qelast = Qelast - matmul(Felast, fem_m2v(steady_app_forces + &
                                                 dynamic_app_forces, &
                                                 numdof,&
                                                 Filter=ListIN))

        call xbeam_asbly_dynamic(numdof,&
                                 n_node,&
                                 n_elem,&
                                 elements,&
                                 nodes,&
                                 pos_ini,&
                                 psi_ini,&
                                 pos_def,&
                                 psi_def,&
                                 0.0d0*pos_def,&
                                 0.0d0*psi_def,&
                                 0.0d0*pos_def,&
                                 0.0d0*psi_def,&
                                 dQdt(numdof+1:numdof+6),&
                                 0.0d0*dQddt(numdof+1:numdof+6),&
                                 dQdt(numdof+7:numdof+10),&
                                 MRS,&
                                 MRR,&
                                 CRS,&
                                 CRR,&
                                 CQR,&
                                 CQQ,&
                                 KRS,&
                                 Frigid,&
                                 Qrigid,&
                                 options,&
                                 Cao)
        Qrigid = Qrigid - matmul(Frigid, fem_m2v(steady_app_forces + &
                                                 dynamic_app_forces, numdof+6))

        if (options%gravity_on) then
            call xbeam_asbly_M_gravity(&
                                         numdof,&
                                         n_node,&
                                         n_elem,&
                                         Elements,&
                                         Nodes,&
                                         pos_ini,&
                                         psi_ini,&
                                         pos_def,&
                                         psi_def,&
                                         MRS_gravity,&
                                         MSS_gravity,&
                                         MRR_gravity,&
                                         options)
            Qrigid = Qrigid + matmul(MRS_gravity, &
                                     cbeam3_asbly_gravity_dynamic(NumDof + 6,&
                                                                  options,&
                                                                  Cao))

            Qelast = Qelast + matmul(MSS, &
                                     cbeam3_asbly_gravity_dynamic(NumDof,&
                                                                  options,&
                                                                  Cao))
        end if

        Qtotal(1:numdof) = Qelast
        Qtotal(numdof+1:numdof+6) = Qrigid
        Qtotal(numdof+7:numdof+10) = matmul(CQQ, dQdt(numdof+7:numdof+10))

        Mtotal(1:numdof, 1:numdof) = MSS
        Mtotal(1:numdof, numdof+1:numdof+6) = MSR
        Mtotal(numdof+1:numdof+6, 1:numdof) = MRS
        Mtotal(numdof+1:numdof+6, numdof+1:numdof+6) = MRR
        Mtotal(numdof+7:numdof+10, numdof+7:numdof+10) = unit4

        call lu_solve(numdof, Mtotal, -Qtotal, dQddt)
        ! ! Predictor step.
        ! Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
        ! ! Q    = initialQ + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
        ! dQdt = dQdt + (1.d0-gamma)*dt*dQddt
        ! ! dQdt = initialdQdt + (1.d0-gamma)*dt*dQddt
        ! dQddt= 0.d0
        !
        ! converged = .FALSE.
        ! do i_iter=1, 1000
        !     print*, 'ITER ', i_iter
        !     Qelast = 0.0d0
        !     Qrigid = 0.0d0
        !     Qtotal = 0.0d0
        !
        !     Frigid = 0.0d0
        !     Felast = 0.0d0
        !
        !     MSS = 0.0d0
        !     MSR = 0.0d0
        !     MRS = 0.0d0
        !     MRR = 0.0d0
        !
        !     CSS = 0.0d0
        !     CSR = 0.0d0
        !     CRS = 0.0d0
        !     CRR = 0.0d0
        !     CQR = 0.0d0
        !     CQQ = 0.0d0
        !
        !     KSS = 0.0d0
        !     KRS = 0.0d0
        !
        !     Mtotal = 0.0d0
        !     Ctotal = 0.0d0
        !     Ktotal = 0.0d0
        !     Asys = 0.0d0
        !
        !
        !     ! call print_matrix('init_q',q)
        !     ! call print_matrix('init_dq',dqdt)
        !
        !     call cbeam3_solv_state2disp(&
        !                                 elements,&
        !                                 nodes,&
        !                                 pos_ini,&
        !                                 psi_ini,&
        !                                 Q(1:numdof),&
        !                                 dQdt(1:numdof),&
        !                                 pos_def,&
        !                                 psi_def,&
        !                                 pos_def_dot,&
        !                                 psi_def_dot)
        !     Cao = xbeam_Rot(dQdt(numdof+7:numdof+10))
        !
        !     !
        !     ! system matrix generation
        !     call cbeam3_asbly_dynamic(numdof,&
        !                               n_node,&
        !                               n_elem,&
        !                               elements,&
        !                               nodes,&
        !                               pos_ini,&
        !                               psi_ini,&
        !                               pos_def,&
        !                               psi_def,&
        !                               0.0d0*pos_def,&
        !                               0.0d0*psi_def,&
        !                               0.0d0*pos_def,&
        !                               0.0d0*psi_def,&
        !                               steady_app_forces + dynamic_app_forces,&
        !                               dQdt(numdof+1:numdof+6),&
        !                               0.0d0*dQddt(numdof+1:numdof+6),&
        !                               MSS,&
        !                               MSR,&
        !                               CSS,&
        !                               CSR,&
        !                               KSS,&
        !                               Felast,&
        !                               Qelast,&
        !                               options,&
        !                               Cao)
        !     mindelta = options%mindelta*max(1.0d0, maxval(abs(Qelast)))
        !     Qelast = Qelast - matmul(Felast, fem_m2v(steady_app_forces + &
        !                                              dynamic_app_forces, &
        !                                              numdof,&
        !                                              Filter=ListIN))
        !
        !     call xbeam_asbly_dynamic(numdof,&
        !                              n_node,&
        !                              n_elem,&
        !                              elements,&
        !                              nodes,&
        !                              pos_ini,&
        !                              psi_ini,&
        !                              pos_def,&
        !                              psi_def,&
        !                              0.0d0*pos_def,&
        !                              0.0d0*psi_def,&
        !                              0.0d0*pos_def,&
        !                              0.0d0*psi_def,&
        !                              dQdt(numdof+1:numdof+6),&
        !                              0.0d0*dQddt(numdof+1:numdof+6),&
        !                              dQdt(numdof+7:numdof+10),&
        !                              MRS,&
        !                              MRR,&
        !                              CRS,&
        !                              CRR,&
        !                              CQR,&
        !                              CQQ,&
        !                              KRS,&
        !                              Frigid,&
        !                              Qrigid,&
        !                              options,&
        !                              Cao)
        !     Qrigid = Qrigid - matmul(Frigid, fem_m2v(steady_app_forces + &
        !                                              dynamic_app_forces, numdof+6))
        !
        !     if (options%gravity_on) then
        !         call xbeam_asbly_MRS_gravity(&
        !                                      numdof,&
        !                                      n_node,&
        !                                      n_elem,&
        !                                      Elements,&
        !                                      Nodes,&
        !                                      pos_ini,&
        !                                      psi_ini,&
        !                                      pos_def,&
        !                                      psi_def,&
        !                                      MRS_gravity,&
        !                                      options)
        !         Qrigid = Qrigid + matmul(MRS_gravity, &
        !                                  cbeam3_asbly_gravity_dynamic(NumDof + 6,&
        !                                                               options,&
        !                                                               Cao))
        !
        !         Qelast = Qelast + matmul(MSS, &
        !                                  cbeam3_asbly_gravity_dynamic(NumDof,&
        !                                                               options,&
        !                                                               Cao))
        !     end if
        !
        !     Mtotal(1:numdof, 1:numdof) = MSS
        !     Mtotal(1:numdof, numdof+1:numdof+6) = MSR
        !     Mtotal(numdof+1:numdof+6, 1:numdof) = MRS
        !     Mtotal(numdof+1:numdof+6, numdof+1:numdof+6) = MRR
        !     Mtotal(numdof+7:numdof+10, numdof+7:numdof+10) = unit4
        !
        !     Qtotal(1:numdof) = Qelast
        !     Qtotal(numdof+1:numdof+6) = Qrigid
        !     Qtotal(numdof+7:numdof+10) = matmul(CQQ, dQdt(numdof+7:numdof+10))
        !
        !     Qtotal = Qtotal + matmul(Mtotal, dQddt)
        !
        !     ! convergence check
        !     print*, maxval(abs(Qtotal)), '  ', mindelta
        !     if (maxval(abs(Qtotal)) < mindelta) then
        !         print*, 'converged in ', iter
        !         converged = .TRUE.
        !     end if
        !
        !     if (converged) then
        !         exit
        !     endif
        !
        !     ! call lu_solve(Mtotal, -Qtotal, dQddt)
        !
        !     ! damping and stiffness matrices
        !     call mat_addmat(0, 0, CSS, Ctotal)
        !     call mat_addmat(0, numdof, CSR, Ctotal)
        !     call mat_addmat(numdof, 0, CRS, Ctotal)
        !     call mat_addmat(numdof, numdof, CRR, Ctotal)
        !     call mat_addmat(numdof+6, numdof, CQR, Ctotal)
        !     call mat_addmat(numdof+6, numdof+6, CQQ, Ctotal)
        !
        !     call mat_addmat(0, 0, KSS, Ktotal)
        !     call mat_addmat(numdof, 0, KRS, Ktotal)
        !
        !     ! assembly of A matrix
        !     Asys = Ktotal
        !     Asys = Asys + Ctotal*gamma/(beta*dt)
        !     Asys = Asys + Mtotal/(beta*dt*dt)
        !
        !     ! calculation of the correction
        !     DQ = 0.0d0
        !     ! call print_matrix('Qtotal',Qtotal)
        !     ! call print_matrix('Frigid',Frigid)
        !     ! call print_matrix('Asys',Asys)
        !     ! call print_matrix('Mtotal',Mtotal)
        !     call lu_solve(Asys, -Qtotal, DQ)
        !     ! call print_matrix('DQ',DQ)
        !
        !
        !     ! reconstruction of state vectors
        !     Q = Q + DQ
        !     dQdt = dQdt + DQ*gamma/(beta*dt)
        !     dQddt = dQddt + DQ/(beta*dt*dt)
        !     ! dQddt = 0.0d0
        ! end do
        ! Predictor step.
        ! Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
        ! ! Q    = initialQ + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
        ! dQdt = dQdt + (1.d0-gamma)*dt*dQddt
        ! ! dQdt = initialdQdt + (1.d0-gamma)*dt*dQddt
        ! dQddt= 0.d0
        call cbeam3_solv_state2disp(elements,&
                                    nodes,&
                                    pos_ini,&
                                    psi_ini,&
                                    Q(1:numdof),&
                                    dQdt(1:numdof),&
                                    pos_def,&
                                    psi_def,&
                                    pos_def_dot,&
                                    psi_def_dot)
        ! dQdt(1:numdof+6) = 0.0d0
        ! dQddt(numdof+1:) = 0.0d0
        for_vel = dQdt(numdof+1:numdof+6)
        for_acc = dQddt(numdof+1:numdof+6)
        dQdt(numdof+7:numdof+10) = rot_unit(dQdt(numdof+7:numdof+10))
        quat = dQdt(numdof+7:numdof+10)
        ! dQddt = 0.0d0
    end subroutine xbeam_solv_nlndyn_init_python




end module xbeam_interface
