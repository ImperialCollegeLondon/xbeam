!-> Copyright by Imperial College London, 2008
!
!-> Module.- XBEAM_SHARED Rafa Palacios. 15Jul2008 - Last Update 07Jan2011 Henrik Hesse
!
!-> Description.-
!
!  This module defines shared parameters and type definitions shared by all
!  XBeam routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_shared
    use, intrinsic          :: iso_c_binding
 implicit none


! Problem constants.
 integer,parameter:: MaxElNod=3                ! Max number of nodes per element.
 real(8),parameter:: Pi=3.14159265358979
 real(8),parameter,private,dimension(3,3):: Unit= &    ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))


! Define derived types with default values:

! Element information.
 type xbelem
  integer:: NumNodes = 0               ! Number of nodes in the element.
  integer:: MemNo = 0                 ! Member to which the element belongs.
  integer:: Conn   (MaxElNod) = 0
    ! Connectivities (i.e., for each node, contains the local number to global number association)
  integer:: Master (MaxElNod,2) = 0   ! Master node for each node j in the element.
                                   ! (j,1): # master elem  (or 0 if current is master).
                                   ! (j,2): Node within master element (or 0 if current is master).
  real(8):: Length = 0.0d0                ! Length in the undeformed configuration.
  real(8):: PreCurv (3) = 0.0d0            ! Initial curvature of the element.
  real(8):: Psi     (3) = 0.0d0            ! Rotation vector of (undeformed) element frame.
  real(8):: Vector  (3) = 0.0d0            ! Element orientatIon vector. It goes along the local Y axis.
                                   ! sm: this defines the element y axis and is in the a frame components
  real(8):: Mass    (6,6) = 0.0d0          ! Mass matrix (constant along the element).
  real(8):: Stiff   (6,6) = 0.0d0          ! Element stiffness (constant along the element).
  real(8):: InvStiff(6,6) = 0.0d0          ! Inverse of the element stiffness (constant along the element).
  real(8):: RBMass  (MaxElNod,6,6) = 0.0d0 ! Non-Structural (lumped) mass at element nodes.
 end type xbelem

! Nodal information.
 type xbnode
  integer:: Master (2) = 0         ! Master node for current node.
                                   ! (1): # master elem.
                                   ! (2): # node within master element.
  integer:: Vdof  = 0                  ! Number of node for which displacements and rotations are unknown
                                   ! (free -internal and not - and clamped) in the velocity/displacements vector
                                   ! For clamped node, this signifies weak enforcement of BCs.
  integer:: Fdof = 0                  ! Number of node in the force vector for which all the forces/moments need to be computed
                                   ! (clamped BC, internal nodes)
  integer:: Sflag=0                ! Flag for spherical joint at the node (1: hinged, 0: no hinge)
                                   ! This will imply weak enforcement of BCs.
                                   ! Only working for static solver and solution 912

 end type xbnode

! Simulation options (with default values).
 type, bind(C) :: xbopts
  logical(c_bool):: FollowerForce   =.true.   ! =T: Follower force.
  logical(c_bool):: FollowerForceRig=.true.   ! =T: Follower force in the body-fixed frame.
  logical(c_bool):: PrintInfo    =.true.      ! =T: Print information on screen.
  logical(c_bool):: OutInBframe  =.false.      ! =T print velocities in B-frame (if not, use a-frame)
  logical(c_bool):: OutInaframe  =.true.     ! =T print velocities in a-frame (if not, Inertial frame)
  integer(c_int):: ElemProj     = 0          ! =0: Element info computed in the global frame.
                                      ! =1: Element info computed in a fixed element frame.
                                      ! =2: Element info computed in a moving element frame.
  integer(c_int):: MaxIterations=99          ! Maximum number of iterations.
  integer(c_int):: NumLoadSteps=1            ! Number of load increments.
  integer(c_int):: NumGauss=2                ! Number of Gauss points in the integration.
  integer(c_int):: Solution=111              ! Solution process:
                                      ! =102/112: cbeam3 linear/nonlinear static.
                                      ! =202/212: cbeam3 linear/nonlinear structural dynamic
                                      ! =302/312: cbeam3 linear/nonlinear static + structural dynamic
                                      ! =900/910:        linear/nonlinear rigid-body dynamic
                                      ! =902/912: cbeam3 linear/nonlinear flexible-body dynamic
                                      ! =    922: cbeam3 nonlinear static + flexible-body dynamic
                                      ! =    952: cbeam3 linear flexible with nonlinear rigid-body dynamic
  real(c_double):: DeltaCurved=1d-2          ! Minimum angle for two unit vectors to be parallel.
  real(c_double):: MinDelta=1d-8             ! Relative convergence parameter for Newton-Raphson iterations.
  real(c_double):: abs_threshold=1d-13       ! Absolute convergence parameter for Newton-Raphson iterations.
  real(c_double):: NewmarkDamp=1.d-4         ! Numerical damping in the Newmark integration scheme.
  logical(c_bool):: gravity_on = .FALSE.
  real(c_double):: gravity = 0.0d0
  real(c_double):: gravity_dir_x = 0
  real(c_double):: gravity_dir_y = 0
  real(c_double):: gravity_dir_z = 1
  logical(c_bool):: balancing = .false.
  real(c_double):: relaxation_factor = 0.3
  logical(c_bool):: load_ramping_conv = .false.
 end type

end module xbeam_shared
