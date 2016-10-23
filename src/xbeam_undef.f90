!-> Module.- XBEAM_UNDEF Rafa Palacios. 15Jul2008
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Compute properties of undeformed structure.
!
!-> Subroutines.-
!
!    -xbeam_undef_geom      : Compute geometric properties of elements.
!    |-xbeam_undef_relatree : Compute tree of nodes shared between elements.
!                              Identify master and slave nodes.
!    -xbeam_undef_dofs      : Identify degrees of freedom in the problem.
!    |-xbeam_undef_nodeindep: Identify independent nodes.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module xbeam_undef
  use xbeam_shared
  implicit none

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_GEOM
!
!-> Description:
!
!    Compute basic geometric parameters on the undeformed configuration.
!
!-> Remarks.-
!
!  1) The undeformed frame (b) is defined within each element. For nodes
!     belonging to two or more elements, there will be as many undeformed
!     frames.
!
!  2) The displaced coordinate system (a) is given as the cartesian rotation
!     vector that defines its orientation with respect to the global frame.
!
! -> sm Example: call from main_xxx.f90:
!   allocate(PsiIni(NumElems,MaxElNod,3)); PsiIni=0.d0
!   call xbeam_undef_geom (Elem,PosIni,PhiNodes,PsiIni,Options)
!        xbeam_undef_geom(inout, in  ,  in    , out  ,  in   )
!   - PosIni   (in): Coords in unput_xxx.f90
!   - PhiNodes (in): pretwist
!   - PsiIni  (out): CRV at the node
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_geom (Elem,Coords,PhiNodes,Psi0,Options)
  use lib_fem
  use lib_bgeom

! I/O Variables.
  type(xbelem), intent(inout):: Elem     (:)     ! Element information.
  real(8),      intent(in)   :: Coords   (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: PhiNodes (:)     ! Pretwist at the nodes.
  real(8),      intent(out)  :: Psi0     (:,:,:) ! Initial CRV at element nodes Psi0(NumElems,MaxElNod,3)
  type(xbopts), intent(in)   :: Options          ! Solver parameters.

! Local variables.
  real(8)            :: ElemPhi(MaxElNod)      ! Pretwist at the nodes of the element.
  integer            :: j                      ! Counter on the elements in the model.
  real(8)            :: LocCoords(MaxElNod,3)  ! Global coordinates of nodes in the element.
 
! Loop in all elements in the model.
  do j=1,size(Elem)

  ! Extract element information.
  ! sm: given the node nn (global numbering), associated to the node ii (local) of
  ! the element j, the function allocates in LocCoords(ii,:) the row Coords(nn,:)
  ! Elem%NumNodes is here recomputed and allocated
    call fem_glob2loc_extract (Elem(j)%Conn,Coords,LocCoords,Elem(j)%NumNodes)
  ! similar call for the pretwist. No need to insert in a function as
  ! Elme%NumNodes is available
    ElemPhi(1:Elem(j)%NumNodes)= PhiNodes(Elem(j)%Conn(1:Elem(j)%NumNodes))

  ! Compute initial displaced coordinate system (a).
  ! sm: input/purpose recap
  ! the 1st and last node of the element [r <-> LocCoords(1:2,:)] define the x axis
  ! the local y axis (in global coordinates) is set in input under Elem%Vector
  ! the function computes the last direction to define the local element frame,
  ! the transformation matrix and the CRV associated.
  !
  ! sm: from lib_bgeom:
  !      bgeom_elemframe (          r (in),        V (in),  Psi (out), Delta(in))
    call bgeom_elemframe (LocCoords(1:2,:),Elem(j)%Vector,Elem(j)%Psi, &
&                         Options%DeltaCurved)

! Compute undeformed frame (B) at all nodes of each element. Given by CRV in Psi0.
    call bgeom_nodeframe (Elem(j)%NumNodes,LocCoords,ElemPhi,Elem(j)%Vector,&
&                         Psi0(j,:,:),Options%DeltaCurved)

! Compute element length.
    call bgeom_elemlength (Options%NumGauss,Elem(j)%NumNodes,LocCoords,Elem(j)%Length)

! Compute curvature in the initial configuration.
    if (Elem(j)%NumNodes.eq.2) then
      call bgeom_elemcurv2 (LocCoords(1:2,:),Psi0(j,:,:),Elem(j)%PreCurv)
    end if

  end do

! Determine tree of connectivities. For each element determine if their nodes are master or
! slaves. In the last case, determine its master node.
  call xbeam_undef_relatree (Elem)

  return
 end subroutine xbeam_undef_geom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_DOFS
!
!-> Description:
!
!    Initialize the solution process.
!
!-> Remarks.-
!
!    1. Nodes(nn) contains information of the nn-th node (global numbering).
!    2. Nodes(nn)%Master=(jj,kk) means that the nn-th node (global numbering) is
!    associated to the kk-th node (local numbering) of the jj-th element.
!    3. Nodes%Master only points to master nodes: this ensure uniqueness in the
!    association
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_dofs (Elem,BoundConds,Node,NumDof, Solution)
  use lib_fem
  use lib_bgeom

! I/O Variables.
  type(xbelem),intent(in)    :: Elem      (:)      ! Element information.
  integer,     intent(in)    :: BoundConds(:)      ! Boundary conditions.
  type(xbnode),intent(inout) :: Node      (:)      ! Nodal information.
  integer,     intent(out)   :: NumDof             ! Number of independent degrees of freedom.
  integer,     intent(inout), optional :: Solution ! Solution number. Used for Sperical Joint BCs only

! Local variables.
  integer:: NumE                           ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.
  integer :: j,k                           ! Counter on the elements.
  integer,allocatable:: ListFr    (:)      ! List of non-free nodes.
  integer,allocatable:: ListIN    (:)      ! List of independent nodes.
  integer,allocatable:: ListSflag (:)      ! Flag for spherical joint BC (1: hinge / 0: no hinge)

! Initialize
  NumE=size(Elem)
  NumN=size(Node)

! Determine ID of master for each node in the model.
  do j=1,NumE
    do k=1,Elem(j)%NumNodes
        if (Elem(j)%Master(k,1).eq.0) then
          Node(Elem(j)%Conn(k))%Master(1)= j
          Node(Elem(j)%Conn(k))%Master(2)= k
        end if
    end do
  end do

! Get list of independent nodes.
  allocate (ListIN(NumN)); ListIN= 0
  allocate (ListFr(NumN)); ListFr= 0
  allocate (ListSflag(NumN)); ListSflag= 0
  ! NumDof is here the number of independent nodes

  if (present(Solution)) then
    call xbeam_undef_nodeindep (NumN,BoundConds,NumDof,ListIN,ListFr,ListSflag,Solution)
  else
    call xbeam_undef_nodeindep (NumN,BoundConds,NumDof,ListIN,ListFr,ListSflag)
  end if
  NumDof=NumDof*6

  do j=1,NumN
    Node(j)%Vdof=ListIN(j)
    Node(j)%Fdof=ListFr(j)
    Node(j)%Sflag=ListSflag(j)
  end do

  deallocate (ListIn,ListFr,ListSflag)

  !print *, 'NumDof: ', NumDof
  !print *, 'Node%Vdof: ', Node%Vdof
  !print *, 'Node%Fdof: ', Node%Fdof
  !print *, 'Node%Sflag:', Node%Sflag

  return
 end subroutine xbeam_undef_dofs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_RELATREE
!
!-> Description:
!
!    Compute tree of nodes shared between elements.
!
!-> Remarks.-
!
!  1) Elem(i)%Master(j,1)=0 -> Node j of element i is a master one.
!
!     Elem(i)%Master(j,1)=i0 & Elem(i)%Master(j,2)=j0 -> Node j of element i
!     is slave node of node j0 of element i0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_relatree (Elem)

! I/O Variables.
  type(xbelem),intent(inout):: Elem (:)     ! Element information

! Local variables.
  integer :: iElem,jElem        ! Counters on the elements in the model.
  integer :: iNode,jNode        ! Counters on the nodes in the elements.

! Loop in the nodes of all elements.
  do iElem=1,size(Elem)

! Initalize to zero (i.e., by default all nodes are master nodes).
    Elem(iElem)%Master=0

! Loop in the nodes within the element
    do iNode=1,Elem(iElem)%NumNodes

! For each node, see if it has already appeared. For that purpose, scan
! all elements from 1 to iElem, and compare their connectivites to the current
! node.
      jElem=1
      do while ((Elem(iElem)%Master(iNode,1).eq.0).and.(jElem.lt.iElem))
        do jNode=1,Elem(iElem)%NumNodes
          if (Elem(jElem)%Conn(jNode) .eq. Elem(iElem)%Conn(iNode)) then
            Elem(iElem)%Master(iNode,1)=jElem
            Elem(iElem)%Master(iNode,2)=jNode
          end if
        end do
        jElem=jElem+1
      end do
    end do
  end do

  return
 end subroutine xbeam_undef_relatree


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_UNDEF_NODEINDEP
!
!-> Description:
!
!    Identify independent nodes in the problem.
!
!-> Remarks.-
!   1. NumFr and NumIN are counters; they count the number of nodes and forces
!      for which to solve.
!      a. If no BCs are applied to a node, both forces and displacements have to
!         be found (both have a +1)
!      b. If a clamped BC is applied, the forces have to be solved, but no the
!         displacements (NumFr = NumFr + 1). Vice-versa for free end, for which
!         only the displacements need to be solved for (NumIN = NumIN + 1).
!   2. ListIN and ListFr, are arrays of length equal to the total number of
!      nodes.
!      a. If ListIN(ii)=0 : the nn-th node (global numbering) is not independent
!      b. If ListIN(ii)=jj: the nn-th node (global numbering) is associated to the
!         jj-th independent nodal dispacement.
!   3. ListSflag BEHAVES DIFFERENTLY FROM ListIN and ListFr: ListHg will return
!      1 if the node is hinged, 0 if not. This will add only 6 dof to the global
!      velocity/displacement/force vector instead of 6
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_undef_nodeindep (NumN,BoundCond,NumIN,ListIN,ListFr,ListSflag,Solution)

! I/O Variables.
  integer,intent(in)::  NumN               ! Number of nodes in the model.
  integer,intent(in)::  BoundCond(:)       ! Boundary conditions: 1: Clamped; -1: free; 0: nothing
  integer,intent(out):: NumIN              ! Number of independent nodes.
  integer,intent(out):: ListIN (:)         ! List of independent nodes.
  integer,intent(out):: ListFr (:)         ! List of nodes with independent force vector (not free nodes).

  integer,intent(out):: ListSflag (:)      ! 1 if node is hinged, 0 if not
  integer,intent(inout), optional :: Solution ! Solution number. Used for Sperical Joint BCs only
  integer :: SolutionCode ! copy of solution number

! Local variables.
  integer:: iNode                          ! Counter on the nodes.
  integer:: NumFr                          ! Counter on the force vector.

! Solution not present: if not present (wrapper) rigid-flex body dynamic solution
! is assumed
  if (present(Solution)) then
      SolutionCode=Solution
  else
      ! for wrapper
      print *, 'xbeam_undef: setting Solution=912!!!'
      SolutionCode=912
  end if

! Loop on the nodes and remove then from the final list if they are constrained.
  NumIN=0
  NumFr=0
  ListIN=0
  ListFr=0
  ListSflag=0

  do iNode=1,NumN
    select case (BoundCond(iNode))
      case (0) ! no BCs (independent nodes)
        NumIN=NumIN+1
        NumFr=NumFr+1
        ListIN(iNode)=NumIN
        ListFr(iNode)=NumFr
      case (-1) ! free end
        NumIN=NumIN+1
        ListIN(iNode)=NumIN
      case (1) ! clamp
        NumFr=NumFr+1
        ListFr(iNode)=NumFr
      case (2)
        select case (SolutionCode)
          case default
            print *, 'Spherical Joint not implemented for Solution: ', Solution
            stop 'Program terminated - xbeam_undef_nodeindep'
          case (102) ! spherical joint implemented
            NumIN=NumIN+1
            NumFr=NumFr+1
            ListIN(iNode)=NumIN
            ListFr(iNode)=NumFr
            ListSflag(iNode)=1
          case (912,932) ! treat the spherical joint as clamped
            !print *, 'Treating spherical joint as a clamp'
            NumFr=NumFr+1
            ListFr(iNode)=NumFr
            ListSflag(iNode)=1
        end select
    end select
  end do

  return
 end subroutine xbeam_undef_nodeindep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module xbeam_undef
