module cbeam3_interface
    use, intrinsic                      :: iso_c_binding
    use                                 :: xbeam_shared
    use                                 :: cbeam3_solv
    use                                 :: debug_utils

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
                                            RBMass,&
                                            master_node,&
                                            vdof,&
                                            fdof,&
                                            options,&
                                            applied_forces,&
                                            pos_ini,&
                                            psi_ini,&
                                            pos_def,&
                                            psi_def) bind(C)
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
        real(c_double), intent(IN)      :: RBMass(n_elem, max_elem_node, 6, 6)

        ! node data
        integer(c_int), intent(IN)      :: master_node(n_node, 2)
        integer(c_int), intent(IN)      :: vdof(n_node)
        integer(c_int), intent(IN)      :: fdof(n_node)

        type(xbelem)                    :: elements(n_elem)
        type(xbnode)                    :: nodes(n_node)
        type(xbopts), intent(IN)        :: options

        real(c_double), intent(IN)      :: applied_forces(n_node, 6)
        real(c_double), intent(IN)      :: pos_ini(n_node, 3)
        real(c_double), intent(IN)      :: psi_ini(n_elem, max_elem_node, 3)
        real(c_double), intent(INOUT)   :: pos_def(n_node, 3)
        real(c_double), intent(INOUT)   :: psi_def(n_elem, max_elem_node, 3)

        integer(c_int)                  :: num_dof


        integer(c_int)                  :: i
        integer(c_int)                  :: j

        num_dof = count(vdof > 0)*6

        print*, 'Num_dof = ', num_dof


        !subroutine cbeam3_solv_nlnstatic (NumDof,Elem,Node,AppForces,Coords,Psi0, &
        !&                                  PosDefor,PsiDefor,Options)
        !use lib_fem
        !use lib_sparse
        !!<<<<<<< HEAD
        !use lib_solv
        !!#ifdef NOLAPACK
        !use lib_lu
        !use cbeam3_asbly

        !! I/O Variables.
        !integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
        !type(xbelem),intent(in)    :: Elem(:)           ! Element information.
        !type(xbnode),intent(in)    :: Node(:)           ! Nodal information.
        !real(8),      intent(in)   :: AppForces (:,:)   ! Applied nodal forces.
        !real(8),      intent(in)   :: Coords   (:,:)    ! Initial coordinates of the grid points.
        !real(8),      intent(in)   :: Psi0     (:,:,:)  ! Initial CRV of the nodes in the elements.
        !real(8),      intent(inout):: PosDefor (:,:)    ! Current coordinates of the grid points
        !real(8),      intent(inout):: PsiDefor (:,:,:)  ! Current CRV of the nodes in the elements.
        !type(xbopts),intent(in)    :: Options           ! Solver parameters.



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
                                   RBMass)

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
       call cbeam3_solv_nlnstatic(num_dof,&
                                  elements,&
                                  nodes,&
                                  applied_forces,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  options)



    end subroutine cbeam3_solv_nlnstatic_python


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
        real(c_double), intent(IN)      :: RBMass(n_elem, max_elem_node, 6, 6)

        type(xbelem)                    :: elements(n_elem)

        integer(c_int)                  :: i


        do i=1, n_elem
            elements(i)%NumNodes    = num_nodes(i)
            elements(i)%MemNo       = mem_number(i)
            elements(i)%Conn        = conn(i, :)
            elements(i)%Master      = master(i, :, :)
            elements(i)%Length      = 0.0d0
            elements(i)%Psi         = [0.0d0, 0.0d0, 0.0d0]
            elements(i)%Vector      = [0.0d0, 1.0d0, 0.0d0]
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


    !function generate_xbopts(follower_force,&
                             !follower_force_rig,&
                             !print_info,&
                             !out_in_a_frame,&
                             !out_in_b_frame,&
                             !elem_proj,&
                             !max_iterations,&
                             !num_load_steps,&
                             !num_gauss,&
                             !solution,&
                             !delta_curved,&
                             !min_delta,&
                             !newmark_damp) result(options)
        !! options
        !logical(c_bool), intent(IN)     :: follower_force
        !logical(c_bool), intent(IN)     :: follower_force_rig
        !logical(c_bool), intent(IN)     :: print_info
        !logical(c_bool), intent(IN)     :: out_in_a_frame
        !logical(c_bool), intent(IN)     :: out_in_b_frame
        !integer(c_int), intent(IN)      :: elem_proj
        !integer(c_int), intent(IN)      :: max_iterations
        !integer(c_int), intent(IN)      :: num_load_steps
        !integer(c_int), intent(IN)      :: num_gauss
        !integer(c_int), intent(IN)      :: solution
        !real(c_double), intent(IN)      :: delta_curved
        !real(c_double), intent(IN)      :: min_delta
        !real(c_double), intent(IN)      :: newmark_damp

        !type(xbopts)                    :: options

        !options%FollowerForce       = follower_force
        !options%FollowerForceRig    = follower_force_rig
        !options%PrintInfo           = print_info
        !options%OutInAFrame         = out_in_a_frame
        !options%OutInBFrame         = out_in_b_frame
        !options%ElemProj            = elem_proj
        !options%MaxIterations       = max_iterations
        !options%NumLoadSteps        = num_load_steps
        !options%NumGauss            = num_gauss
        !options%Solution            = solutions
        !options%DeltaCurved         = delta_curved
        !options%MinDelta            = min_delta
        !options%NewmarkDamp         = newmark_damp
    !end function generate_xbopts






end module cbeam3_interface
