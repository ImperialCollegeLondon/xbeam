!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Module.- XBEAM_SOLV  Henrik Hesse. 07/01/2011 - Last Update 07/01/2011
!                        Salvatore Maraniello. 20/10/2014 - Last Update 23/10/2014
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!   Solve beam equations.
!
!-> Subroutines.-
!
!   -xbeam_solv_coupledlindyn:      Coupled linear dynamic solution
!   -xbeam_solv_couplednlindyn:     Coupled nonlinear dynamic solution
!   -xbeam_solv_rigidlindyn:        Linear rigid-body dynamic solution
!   -xbeam_solv_rigidnlindyn:       Nonlinear rigid-body dynamic solution
!
!-> Remark:
! - SM: xbeam_solv_couplednlindyn modified to account for spherical BCs.
! - SM: convergence check improved for sol932        ?
!   sm: uneffective for sol932. To be investigated   ?
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module xbeam_solv
 use xbeam_shared
 use lib_solv
 implicit none

! Shared variables within the module.
 integer,private,parameter:: MaxNodCB3=3               ! Max number of nodes per element is 3.
 integer,private,parameter:: DimMat=18                 ! Memory index for sparse matrices.

 real(8),private,parameter,dimension(4,4):: Unit4= &       ! 4x4 Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0/),(/4,4/))

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine xbeam_SOLV_COUPLEDLINDYN
!
!-> Description:
!
!    Linear dynamic solution of coupled multibeam problem for given applied forces
!    accounting for interaction between structural and rigid-body dynamics.
!
!-> Remarks:
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_solv_coupledlindyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,           &
&                                     Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,  &
&                                     PosDotDefor,PsiDotDefor,DynOut,Options)
  use lib_rot
  use lib_fem
  use lib_sparse
  use lib_out
  use lib_xbeam
  use lib_lu
  use cbeam3_asbly
  use cbeam3_solv
  use xbeam_asbly

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Quat      (4)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Undeformed coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Undeformed CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Initial/final position vector of grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep                          ! Current time step.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: DX(:),DXDt(:),DXDDt(:)     ! Generalized coordinates and derivatives for structure.
  real(8),allocatable:: DQ(:),DQDt(:),DQDDt(:)     ! Generalized coordinates and derivatives for coupled system.

  integer,allocatable:: ListIN     (:)     ! List of independent nodes.
  real(8),allocatable:: DeltaPos   (:,:)   ! Initial/final position vector of grid points
  real(8),allocatable:: DeltaPsi   (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),allocatable:: DeltaPosDot(:,:)   ! Current time derivatives of the coordinates of the grid points
  real(8),allocatable:: DeltaPsiDot(:,:,:) ! Current time derivatives of the CRV of the nodes in the elements.

  ! Define vaiables for structural system matrices
  integer:: as,cs,ks,ms,fs,mr,cr,kr,fr,ctot,ktot,mtot
  type(sparse),allocatable:: Asys   (:)     ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: CSS(:)         ! Sparse damping matrix.
  type(sparse),allocatable:: KSS(:)         ! Elast stiffness matrix in sparse storage.
  type(sparse),allocatable:: MSS(:)         ! Elast mass matrix in sparse storage.
  type(sparse),allocatable:: Felast(:)      ! Applied external forces on structure
  real(8),allocatable::      Qelast(:)      ! Elast vector of discrete generalize forces.
  real(8),allocatable::      MSR(:,:)       ! Mass and damping from the motions of reference system.
  real(8),allocatable::      CSR(:,:)       ! Mass and damping from the motions of reference system.

  ! Define variables for rigid system matrices.
  type(sparse),allocatable:: MRS(:)             ! Sparse damping matrix.
  type(sparse),allocatable:: CRS(:)             ! Sparse damping matrix.
  type(sparse),allocatable:: KRS(:)             ! elast stiffness matrix in sparse storage.
  type(sparse),allocatable:: Frigid(:)          ! Applied external forces on rigid-body motion
  real(8),allocatable::      Qrigid(:)          ! Elast vector of discrete generalize forces.
  real(8),allocatable::      MRR(:,:)           ! Mass and damping from the motions of reference system.
  real(8),allocatable::      CRR(:,:)           ! Mass and damping from the motions of reference system.
  real(8),allocatable::      CQR(:,:),CQQ(:,:)  ! Tangent matrices from linearisation of quaternion equation.

  ! Define variables for rigid-body motion
  real(8),allocatable:: Cao (:,:)                   ! Rotation operator from reference to inertial frame
  real(8),allocatable:: ACoa(:,:)                   ! Rotation operator from reference to inertial frame

  ! Define variables for complete system matrices.
  type(sparse),allocatable:: Ctotal(:)    ! Total Sparse damping matrix.
  type(sparse),allocatable:: Ktotal(:)    ! Total stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mtotal(:)    ! Total mass matrix in sparse storage.
  real(8),allocatable::      Qtotal(:)    ! Total vector of discrete generalize forces.

  character(len=80)  :: Text          ! Text with current time information to print out.
  real(8),allocatable:: Displ(:,:),Veloc(:,:)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (4*DimMat*NumDof)); call sparse_zero (as,Asys)

  allocate (MSS(DimMat*NumDof)); call sparse_zero (ms,MSS)
  allocate (CSS(DimMat*NumDof)); call sparse_zero (cs,CSS)
  allocate (KSS(DimMat*NumDof)); call sparse_zero (ks,KSS)
  allocate (Felast(DimMat*NumDof)); call sparse_zero (fs,Felast)
  allocate (Qelast(NumDof));   Qelast = 0.d0
  allocate (MSR   (NumDof,6)); MSR    = 0.d0
  allocate (CSR   (NumDof,6)); CSR    = 0.d0

  allocate (MRS(DimMat*NumDof)); call sparse_zero (mr,MRS)
  allocate (CRS(DimMat*NumDof)); call sparse_zero (cr,CRS)
  allocate (KRS(DimMat*NumDof)); call sparse_zero (kr,KRS)
  allocate (Frigid(DimMat*NumDof)); call sparse_zero (fr,Frigid)
  allocate (Qrigid   (6));      Qrigid      = 0.d0
  allocate (MRR(6,6));    MRR   = 0.d0
  allocate (CRR(6,6));    CRR   = 0.d0

  allocate (CQR(4,6));    CQR   = 0.d0
  allocate (CQQ(4,4));    CQQ   = 0.d0

  allocate (Mtotal(DimMat*NumDof)); call sparse_zero (mtot,Mtotal)
  allocate (Ctotal(DimMat*NumDof)); call sparse_zero (ctot,Ctotal)
  allocate (Ktotal(DimMat*NumDof)); call sparse_zero (ktot,Ktotal)
  allocate (Qtotal(NumDof+10));   Qtotal= 0.d0

! Updated state vector with rigid body states
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (DXDt  (NumDof)); DXDt   = 0.d0
  allocate (DXDDt (NumDof)); DXDDt  = 0.d0

  allocate (DQ    (NumDof+10)); DQ     = 0.d0
  allocate (DQDt  (NumDof+10)); DQDt   = 0.d0
  allocate (DQDDt (NumDof+10)); DQDDt  = 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  allocate(Cao      (3,3));   Cao       = Unit
  allocate(ACoa     (6,6));   ACoa      = 0.d0

! Compute system information at initial condition.
  allocate (Displ      (NumN,            6)); Displ       = 0.d0
  allocate (Veloc      (NumN,            6)); Veloc       = 0.d0
  allocate (DeltaPos   (NumN,            3)); DeltaPos    = 0.d0
  allocate (DeltaPsi   (NumE(1),MaxElNod,3)); DeltaPsi    = 0.d0
  allocate (DeltaPosDot(NumN,            3)); DeltaPosDot = 0.d0
  allocate (DeltaPsiDot(NumE(1),MaxElNod,3)); DeltaPsiDot = 0.d0

! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,DX,DXDt)

! Check here how to enforce BC for DQDt
  DQDt(NumDof+1:NumDof+6) = Vrel(1,:)
  DQDt(NumDof+7:NumDof+10)= Quat

  Cao = xbeam_Rot(DQDt(NumDof+7:NumDof+10))

! Compute tangent matrices at initial time.
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                   &
&                            0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),0.d0*VrelDot(1,:),  &
&                            ms,MSS,MSR,cs,CSS,CSR,ks,KSS,fs,Felast,Qelast,Options,Cao)

  call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                    &
&                           0.d0*PosDefor,0.d0*PsiDefor,Vrel(1,:),0.d0*VrelDot(1,:),DQDt(NumDof+7:NumDof+10),   &
&                           mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

! Assembly of total matrices combining structural as well as rigid-body matrices.
  ! Mass matrix
  call sparse_addsparse (0,0,ms,MSS,mtot,Mtotal)
  call sparse_addmat    (0,NumDof,MSR,mtot,Mtotal)
  call sparse_addmat    (NumDof,0,transpose(MSR),mtot,Mtotal)
  call sparse_addmat    (NumDof,NumDof,MRR,mtot,Mtotal)
  ! Damping matrix
  call sparse_addsparse (0,0,cs,CSS,ctot,Ctotal)
  call sparse_addmat    (0,NumDof,CSR,ctot,Ctotal)
  call sparse_addsparse (NumDof,0,cr,CRS,ctot,Ctotal)
  call sparse_addmat    (NumDof,NumDof,CRR,ctot,Ctotal)
  ! Stiffness matrix
  call sparse_addsparse (0,0,ks,KSS,ktot,Ktotal)
  call sparse_addsparse (NumDof,0,kr,KRS,ktot,Ktotal)
  ! Contribution of quaternions to mass and damping matrices
  call sparse_addmat    (NumDof+6,NumDof+6,Unit4,mtot,Mtotal)
  call sparse_addmat    (NumDof+6,NumDof  ,CQR,ctot,Ctotal)
  call sparse_addmat    (NumDof+6,NumDof+6,CQQ,ctot,Ctotal)

! Find initial acceleration. Account for initial conditions with velocities and acceleration of the a frame
  Qelast  = sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN)) &
&           - matmul(MSR,VrelDot(1,:)) - matmul(CSR,Vrel(1,:))
  Qrigid  = sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Ftime(1)*Fa,NumDof+6)) &
&           - matmul(MRR,VrelDot(1,:)) - matmul(CRR,Vrel(1,:))

  Qtotal(1:NumDof)          = Qelast
  Qtotal(NumDof+1:NumDof+6) = Qrigid
  Qtotal(NumDof+7:NumDof+10)=-matmul(CQQ,DQDt(NumDof+7:NumDof+10))

  call lu_sparse(mtot,Mtotal,Qtotal,DQDDt)

! Loop in the time steps. 
  do iStep=1,size(Time)-1
    dt     =Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A)') trim(Text)

! Predictor step.
    DQ   = DQ   + dt*DQDt + (0.5d0-beta)*dt*dt*DQDDt
    DQDt = DQDt + (1.d0-gamma)*dt*DQDDt

! Update quaternions with new states
    Cao = xbeam_Rot(DQDt(NumDof+7:NumDof+10))

! Compute system functionals and matrices. Only updating Felast, Frigid, CQR and CQQ to update the orientation of the body-fixed FoR
! Need to make sure that CSS, CSR, CRS, CRR are not overwritten.
    CQR = 0.d0; CQQ = 0.d0
    call sparse_zero (fs,Felast); call sparse_zero (fr,Frigid)
    
    call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor+DeltaPos,PsiDefor+DeltaPsi,                               &
&                              PosDotDefor+DeltaPosDot,PsiDotDefor+DeltaPsiDot,0.d0*PosDefor,0.d0*PsiDefor,             &
&                              F0+Ftime(iStep+1)*Fa,dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),              &
&                              ms,MSS,MSR,ms,MSS,MSR,ms,MSS,fs,Felast,Qelast,Options,Cao)

    call xbeam_asbly_orient (Elem,Node,PosDefor+DeltaPos,PsiDefor+DeltaPsi,dQdt(NumDof+1:NumDof+6),dQdt(NumDof+7:NumDof+10),   &
&                            CQR,CQQ,fr,Frigid,Options,Cao)

! Damping matrix with updated CQR, CQQ
    call sparse_zero      (ctot,Ctotal)
    call sparse_addsparse (0,0,cs,CSS,ctot,Ctotal)
    call sparse_addmat    (0,NumDof,CSR,ctot,Ctotal)
    call sparse_addsparse (NumDof,0,cr,CRS,ctot,Ctotal)
    call sparse_addmat    (NumDof,NumDof,CRR,ctot,Ctotal)
    call sparse_addmat    (NumDof+6,NumDof  ,CQR,ctot,Ctotal)
    call sparse_addmat    (NumDof+6,NumDof+6,CQQ,ctot,Ctotal)

! Compute right-hand side of the equation.
    Qelast  = sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
    Qrigid  = sparse_matvmul(fr,Frigid,6     ,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof+6))

    Qtotal=0.d0
    Qtotal(1:NumDof)          = Qelast
    Qtotal(NumDof+1:NumDof+6) = Qrigid
    Qtotal(NumDof+7:NumDof+10)= matmul(CQQ,dQdt(NumDof+7:NumDof+10))

    Qtotal= Qtotal - sparse_matvmul(ctot,Ctotal,NumDof+10,DQDt) &
&                  - sparse_matvmul(ktot,Ktotal,NumDof+10,DQ)

! Compute left-hand side of the equation (if dt=constant then Asys=constant too, but this was not
! assumed here).
    call sparse_zero (as,Asys)
    call sparse_addsparse(0,0,mtot,Mtotal,as,Asys,Factor=1.d0)
    call sparse_addsparse(0,0,ctot,Ctotal,as,Asys,Factor=gamma*dt)
    call sparse_addsparse(0,0,ktot,Ktotal,as,Asys,Factor=beta*dt*dt)

! Solve equation.
    call lu_sparse(as,Asys,Qtotal,DQDDt)

    DQ    = DQ   + beta *dt*dt*DQDDt
    DQDt  = DQDt + gamma*dt   *DQDDt

! Update positions and velocities on the current converged time step.
    DX  = DQ  (1:NumDof)
    DXDt= DQDt(1:NumDof)

! Update solution and store relevant information at current time step.
    call cbeam3_solv_update_lindyn (Elem,Node,PsiDefor,DX,DXDt,DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)

! Postprocessing
!    DQDt(NumDof+7:NumDof+10)=DQDt(NumDof+7:NumDof+10)/xbeam_2norm(DQDt(NumDof+7:NumDof+10))
    Cao = xbeam_Rot(DQDt(NumDof+7:NumDof+10))
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)
    
! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = DQDt (NumDof+1:NumDof+6)
        VrelDot(iStep+1,:) = DQDDt(NumDof+1:NumDof+6)
! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,DQDt (NumDof+1:NumDof+6))
        VrelDot(iStep+1,:) = matmul(ACoa,DQDDt(NumDof+1:NumDof+6))
    end if

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:) + DeltaPos(k,:)
    end do

  end do

! Write information at last time step.
  PosDefor= PosDefor + DeltaPos
  PsiDefor= PsiDefor + DeltaPsi

  deallocate (MSS,MSR,MRS,MRR,CSS,CSR,CRS,CRR,KSS,KRS)
  deallocate (Mtotal,Ctotal,Ktotal,Qtotal)
  deallocate (Asys,Felast,Qelast,Frigid,Qrigid)
  deallocate (DX,DXDt,DXDDt,DQ,DQDt,DQDDt)
  deallocate (ListIn,Displ,Veloc)
  deallocate (DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)
  return
 end subroutine xbeam_solv_coupledlindyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_SOLV_COUPLEDNLNDYN
!
!-> Description:
!
!    Nonlinear dynamic solution of multibeam problem under applied forces,
!    coupled with rigid-body dynamics
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_solv_couplednlndyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,           &
&                                     Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,  &
&                                     PosDotDefor,PsiDotDefor,DynOut,Options,SUCCESS)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_out
  use lib_sparse
  use lib_xbeam
  use lib_lu
  use cbeam3_asbly
  use cbeam3_solv
  use xbeam_asbly

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  !real(8),      intent(in)   :: Fdyn      (:,:,:) ! applied dynamic force of size (NumNodes, 6, NumSteps+1)
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(inout):: Quat      (4)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

 logical, intent(inout), optional :: SUCCESS  ! Variable to allow python wrapper to handle ecceptions.
                                              ! If the solution does not converge, the variable is set to .false.


! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)   ! Generalized coordinates and derivatives for structure.
  real(8),allocatable:: X(:)
  real(8),allocatable:: dQdt(:),dQddt(:)   ! Generalized coordinates and derivatives for coupled system.
  real(8),allocatable:: Q(:), DQ(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for structure system matrices.

  integer:: as,cs,ks,ms,fs,cr,kr,mr,fr,ctot,ktot,mtot ! size of sparse matrices
  type(sparse),allocatable:: Asys(:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: CSS(:)     ! Sparse damping matrix.
  type(sparse),allocatable:: KSS(:)     ! Elast stiffness matrix in sparse storage.
  type(sparse),allocatable:: MSS(:)     ! Elast mass matrix in sparse storage.
  type(sparse),allocatable:: Felast(:)  ! Applied external forces on structure
  real(8),allocatable::      Qelast(:)  ! Elast vector of discrete generalize forces.
  real(8),allocatable::      MSR(:,:)   ! Mass and damping from the motions of reference system.
  real(8),allocatable::      CSR(:,:)   ! Mass and damping from the motions of reference system.
  
  ! Define variables for rigid system matrices.
  type(sparse),allocatable:: CRS(:)      ! rigid Sparse damping matrix.
  type(sparse),allocatable:: KRS(:)      ! rigid stiffness matrix in sparse storage.
  type(sparse),allocatable:: MRS(:)      ! rigid mass matrix in sparse storage.
  type(sparse),allocatable:: Frigid(:)   ! rigid matrix of applied forces in sparse format
  real(8),allocatable::      Qrigid(:)   ! rigid vector of discrete generalize forces.
  real(8),allocatable::      MRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CQR(:,:),CQQ(:,:)  ! Tangent matrices from linearisation of quaternion equation.

  ! Define variables for rigid-body motion
  real(8),allocatable:: Cao (:,:)                   ! Rotation operator from reference to inertial frame
  real(8),allocatable:: ACoa(:,:)                   ! Rotation operator from reference to inertial frame

  ! Define variables for complete system matrices.
  type(sparse),allocatable:: Ctotal(:)    ! Total Sparse damping matrix.
  type(sparse),allocatable:: Ktotal(:)    ! Total stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mtotal(:)    ! Total mass matrix in sparse storage.
  real(8),allocatable::      Qtotal(:)    ! Total vector of discrete generalize forces.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! =Define variables for sperical joint BC
  logical :: SphFlag
  integer :: sph_rows(3) ! rows to be modified to include a spherical joint BC

  ! Parameters to Check Convergence
  logical :: converged       = .false.
  logical :: passed_delta    = .false.! true if the subiteration (newton) converged according to delta check
  logical :: passed_res      = .false.! true if the subiteration (newton) converged according to residual check
  logical :: passed_resfrc   = .false.! true if the subiteration (newton) converged according to residual check (forces only)
  logical :: passed_resmmt   = .false.! true if the subiteration (newton) converged according to residual check (Moments only)
  logical :: passed_err      = .false.! true if the subiteration (newton) converged according to error check
  logical :: passed_errpos   = .false.! true if the subiteration (newton) converged according to error check (position only)
  logical :: passed_errpsi   = .false.! true if the subiteration (newton) converged according to error check (rotations only)

  ! variables for residual check
  real(8) :: Fsc, Msc                   ! scaling factor for forces and moments
  real(8) :: QglFrc(NumDof/2), QglMmt(NumDof/2) ! residual related to forces and moments
  real(8) :: Res, Res0                  ! current and  initial residual
  real(8) :: ResFrc, ResFrc0, ResFrcRel ! absolute, initial and relative residual
  real(8) :: ResMmt, ResMmt0, ResMmtRel ! absolute, initial and relative residual
  real(8) :: TaRes, TaResFrc, TaResMmt  ! absolute tolerance for Res, ResMmt, ResFrc

 ! variables for error based check
  real(8)   :: Possc, Psisc         ! scaling factor for position and rotations (Psi and quaternion)
  real(8)   :: DeltaPos(NumDof/2)   ! delta displacements at current and previous iteration
  real(8)   :: DeltaPsi(NumDof/2)   ! delta rotations at current and previous iteration
  real(8)   :: ErrX, ErrPos, ErrPsi ! Error extimation for all DoF, displacements and rotations
  real(8)   :: DX_now, DX_old       ! Norm of DeltaX at current and old iteration
  real(8)   :: DPos_now, DPos_old   ! Norm of translational dofs of DeltaX at current and old iteration
  real(8)   :: DPsi_now, DPsi_old   ! Norm of rotational dofs of DeltaX at current and old iteration
  real(8)   :: TaX, TaPos, TaPsi    ! Absolute tolerance for DeltaX, DeltaPos and DeltaPsi

 ! compatibility
  real(8), allocatable  :: Fdyn      (:,:,:) ! applied dynamic force of size (NumNodes, 6, NumSteps+1)
  allocate (Fdyn( size(Fa(:,1)),6,size(Time) )); Fdyn=0.0_8;
  do iStep=1,size(Time)
    Fdyn(:,:,iStep)=Ftime(iStep)*Fa;
  end do


  ! Determine scaling factors for convergence test (absolute tolerances)
  Psisc = 1.0_8 ! used for quaternion scaling as well
  Possc = maxval(abs(Coords))
  Fsc = 1.0_8; Msc = 1.0_8;
  do iStep=1,size(Time)-1
    Fsc = max( maxval(abs( F0(:,1:3)+Fdyn(:,1:3,iStep+1) )), Fsc);
    Msc = max( maxval(abs( F0(:,4:6)+Fdyn(:,4:6,iStep+1) )), Msc);
  end do
  !print*, 'xbeam_solv_couplednlndyn:'
  !print*, 'scaling factor for simulation:'
  !print*, 'Fsc=', Fsc, 'Msc=', Msc

  ! Initialise
  NumN=size(Node)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do

  ! Initialize (sperical joint) sm
  SphFlag=.false.
  do k=1,NumN
      if (Node(k)%Sflag == 1) then
        SphFlag=.true.
      end if
  end do
  sph_rows = (/1,2,3/)+NumDof


  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys(4*DimMat*NumDof)); call sparse_zero (as,Asys)

  allocate (MSS(DimMat*NumDof));    call sparse_zero (ms,MSS)
  allocate (CSS(DimMat*NumDof));    call sparse_zero (cs,CSS)
  allocate (KSS(DimMat*NumDof));    call sparse_zero (ks,KSS)
  allocate (Felast(DimMat*NumDof)); call sparse_zero (fs,Felast)
  allocate (Qelast(NumDof));        Qelast = 0.d0
  allocate (MSR   (NumDof,6));      MSR    = 0.d0
  allocate (CSR   (NumDof,6));      CSR    = 0.d0

  allocate (MRS(DimMat*NumDof));    call sparse_zero (mr,MRS)
  allocate (CRS(DimMat*NumDof));    call sparse_zero (cr,CRS)
  allocate (KRS(DimMat*NumDof));    call sparse_zero (kr,KRS)
  allocate (Frigid(DimMat*NumDof)); call sparse_zero (fr,Frigid)
  allocate (Qrigid(6));             Qrigid = 0.d0
  allocate (MRR   (6,6));           MRR    = 0.d0
  allocate (CRR   (6,6));           CRR    = 0.d0

  allocate (CQR(4,6));    CQR   = 0.d0
  allocate (CQQ(4,4));    CQQ   = 0.d0

  allocate (Mtotal(2*DimMat*NumDof));   call sparse_zero (mtot,Mtotal)
  allocate (Ctotal(2*DimMat*NumDof));   call sparse_zero (ctot,Ctotal)
  allocate (Ktotal(2*DimMat*NumDof));   call sparse_zero (ktot,Ktotal)
  allocate (Qtotal(NumDof+6+4));        Qtotal= 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Updated state vector with rigid body states and quaternions
  allocate (Q     (NumDof+6+4)); Q      = 0.d0
  allocate (dQdt  (NumDof+6+4)); dQdt   = 0.d0
  allocate (dQddt (NumDof+6+4)); dQddt  = 0.d0
  allocate (DQ    (NumDof+6+4)); DQ     = 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  allocate(Cao      (3,3));   Cao       = Unit
  allocate(ACoa     (6,6));   ACoa      = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

  ! sm exceptions handling for python wrapper
  if (present(SUCCESS)) then
      SUCCESS=.true.
  end if

  ! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)

  Q(1:NumDof)             = X(:)
  Q(NumDof+1:NumDof+6)    = 0
  dQdt(1:NumDof)          = dXdt(:)
  dQdt(NumDof+1:NumDof+6) = Vrel(1,:)
  dQdt(NumDof+7:NumDof+10)= Quat

  ! (spherical joint - sm)
  if (SphFlag) then
    dQdt(NumDof+1:NumDof+3) = 0.d0
  end if

  Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

  ACoa(1:3,1:3) = transpose(Cao)
  ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
  if (Options%OutInaframe) then
      Vrel   (1,:) = dQdt (NumDof+1:NumDof+6)
      VrelDot(1,:) = dQddt(NumDof+1:NumDof+6)

! Export velocities and accelerations in inertial frame
  else
      Vrel   (1,:) = matmul(ACoa,dQdt (NumDof+1:NumDof+6))
      VrelDot(1,:) = matmul(ACoa,dQddt(NumDof+1:NumDof+6))
  end if

  do k=1,NumN
      DynOut(k,:) = PosDefor(k,:)
  end do

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
! sm: why we are neglecting it? Let's set it to be on-zero
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                            F0+Fdyn(:,:,1),dQdt(NumDof+1:NumDof+6),            &
&                            0.d0*dQddt(NumDof+1:NumDof+6),                     &
&                            ms,MSS,MSR,cs,CSS,CSR,ks,KSS,fs,Felast,Qelast,Options,Cao)

  Qelast= Qelast - sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Fdyn(:,:,1),NumDof,Filter=ListIN))

  call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,    &
&                           dQdt(NumDof+1:NumDof+6),                            &
&                           0.d0*dQddt(NumDof+1:NumDof+6),                      &
&                           dQdt(NumDof+7:NumDof+10),                           &
&                           mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

  ! Sperical Joint - sm:
  !  Qrigid: Total generalized forces on the element: gyroscopic, lumped masses.
  !          This is non-linear and acts like residual
  !  Frigid: influence coeff. matrix for applied forces.
  ! All terms related to the translational dof and that will influence the RHS
  ! of Newmark-beta time stepping are set to zero.
  if (SphFlag) then
    call sparse_set_rows_zero((/1,2,3/),fr,Frigid)
    Qrigid(1:3)=0.d0
  end if
  Qrigid= Qrigid - sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Fdyn(:,:,1),NumDof+6))

  !!!print *, 'MRR:'
  !!!write (*,'(5X,1P6E12.4)') MRR

! Assemble coupled system
  Qtotal(1:NumDof)          = Qelast
  Qtotal(NumDof+1:NumDof+6) = Qrigid
  Qtotal(NumDof+7:NumDof+10)= matmul(CQQ,dQdt(NumDof+7:NumDof+10))

  call sparse_addsparse (0,0,ms,MSS,mtot,Mtotal)
  call sparse_addmat    (0,NumDof,MSR,mtot,Mtotal)
  call sparse_addmat    (NumDof,0,transpose(MSR),mtot,Mtotal)
  call sparse_addmat    (NumDof,NumDof,MRR,mtot,Mtotal)
  call sparse_addmat    (NumDof+6,NumDof+6,Unit4,mtot,Mtotal)

! Spherical Joint - sm:
! unit diag term in Mtotal to the translational dof of the a FoR. The aim is to
! implement M*dBetadt = 0
! Note: also the columns deleted (though unecessary if the beam root is not
! moving) to avoid numerical error. If a moving hinge/spherical joint is implemented,
! the columns related to the translational dof can't be set to zero!
  if (SphFlag) then
    call sparse_set_colrows_zero(sph_rows,mtot,Mtotal)
    call sparse_set_rows_unit(sph_rows,mtot,Mtotal)
  end if

! Solve matrix system
  call lu_sparse(mtot,Mtotal,-Qtotal,dQddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    if (Options%PrintInfo) then
      write (*,'(X)')
      write (*,'(5X,A,$)') trim(Text)
    end if

! Predictor step.
    Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
    dQdt = dQdt + (1.d0-gamma)*dt*dQddt
    dQddt= 0.d0

    ! Spherical Joint - sm
    ! Predictor step is forced to be zero for translational dofs.
    if (SphFlag) then
      Q(NumDof+1:NumDof+3)=0.d0
      dQdt(NumDof+1:NumDof+3)=0.d0
      dQddt(NumDof+1:NumDof+3)=0.d0
    end if

! Iteration until convergence.
  converged=.false.
! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      !Spherical Joint - sm
      ! reminder of possible issues
      if (Iter.gt.Options%MaxIterations) then
        print *, '\N success=', SUCCESS
        print *, 'Reminders:'
        print *, '1. Old convergence criteria still used. Consider upgrading as per static solver.'
        print *, '2. Spherical Joint, not hinge!'
        !!! sm: stop commented to allow python wrapper to handle exceptions
        if (present(SUCCESS)) then
            SUCCESS=.false.
            ! always print last iteration delta if crash occurrs
            if (.not.(Options%PrintInfo)) then
              write (*,'(5X,A,$)') trim(Text)
              write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qtotal))
            end if
            print *, 'Solution did not converge (18235)'
            exit
            ! stop the iterations
        else
            STOP 'Solution did not converge (18235)'
        end if
      end if

! Update nodal positions and velocities .
      X   (:) = Q   (1:NumDof)
      dXdt(:) = dQdt(1:NumDof)
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Update quaternions with new states
      Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qelast = 0.d0
      MSR    = 0.d0
      CSR    = 0.d0
      call sparse_zero (ms,MSS)
      call sparse_zero (cs,CSS)
      call sparse_zero (ks,KSS)
      call sparse_zero (fs,Felast)
      Qrigid = 0.d0
      MRR    = 0.d0
      CRR    = 0.d0
      CQR    = 0.d0
      CQQ    = 0.d0
      call sparse_zero (mr,MRS)
      call sparse_zero (cr,CRS)
      call sparse_zero (kr,KRS)
      call sparse_zero (fr,Frigid)
      Qtotal = 0.d0
      call sparse_zero (mtot,Mtotal)
      call sparse_zero (ctot,Ctotal)
      call sparse_zero (ktot,Ktotal)

      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                                F0+Fdyn(:,:,iStep+1),dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),                    &
&                                ms,MSS,MSR,cs,CSS,CSR,ks,KSS,fs,Felast,Qelast,Options,Cao)

      call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                               dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),dQdt(NumDof+7:NumDof+10),                &
&                               mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

      ! Spherical Joint - sm:
      ! for spherical joint doesn't matter the orientation of a, we want all the
      ! global external force to be zero.
      ! Note: for hinge BCs no rotations should be required to understand which
      ! moment has to be set to zero.
      if (SphFlag) then
        call sparse_set_rows_zero((/1,2,3/),fr,Frigid)
        Qrigid(1:3)=0.d0
      end if

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qelast)))

! Compute the residual.
      Qelast = Qelast - sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Fdyn(:,:,iStep+1),NumDof,Filter=ListIN))
      Qrigid = Qrigid - sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Fdyn(:,:,iStep+1),NumDof+6))

      Qtotal(1:NumDof)          = Qelast
      Qtotal(NumDof+1:NumDof+6) = Qrigid
      Qtotal(NumDof+7:NumDof+10)= matmul(CQQ,dQdt(NumDof+7:NumDof+10))

      call sparse_addsparse (0,0,ms,MSS,mtot,Mtotal)
      call sparse_addmat    (0,NumDof,MSR,mtot,Mtotal)
!      call sparse_addsparse (NumDof,0,mr,MRS,mtot,Mtotal)
      call sparse_addmat    (NumDof,0,transpose(MSR),mtot,Mtotal)
      call sparse_addmat    (NumDof,NumDof,MRR,mtot,Mtotal)
      call sparse_addmat    (NumDof+6,NumDof+6,Unit4,mtot,Mtotal)

      ! Spherical Joint - sm:
      ! set Mass matrix to unit. If the hinge/spherical joint translates, the
      ! coulmn part cannot be set to zero (see also above).
      if (SphFlag) then
        call sparse_set_colrows_zero(sph_rows,mtot,Mtotal)
        call sparse_set_rows_unit(sph_rows,mtot,Mtotal)
      end if

      Qtotal= Qtotal + sparse_matvmul(mtot,Mtotal,NumDof+6+4,dQddt)

      ! ------------------------------------------------------ Check convergence.

      if (Options%PrintInfo) then
        ! delta check (original code)
        !write (*,'(2X,1PE10.3,$)') maxval(abs(Qtotal))
        write (*,'(5X,A,I4,A,1PE12.3,$)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qtotal))
      end if

      if (maxval(abs(Qtotal)).lt.MinDelta) then
          converged=.true. ! sm
      end if

      ! exit do loop: kept outside in case of other criteria will be introduced
      if (converged.eqv..true.) then
        exit
      end if

! ----------------------------------------------------- end check of convergence


! Compute total damping and stiffness matrices linear system
      call sparse_addsparse (0,0,cs,CSS,ctot,Ctotal)
      call sparse_addmat    (0,NumDof,CSR,ctot,Ctotal)
      call sparse_addsparse (NumDof,0,cr,CRS,ctot,Ctotal)
      call sparse_addmat    (NumDof,NumDof,CRR,ctot,Ctotal)

      call sparse_addsparse (0,0,ks,KSS,ktot,Ktotal)
      call sparse_addsparse (NumDof,0,kr,KRS,ktot,Ktotal)

! Contribution of quaternions to damping matrix
      call sparse_addmat    (NumDof+6,NumDof  ,CQR,ctot,Ctotal)
      call sparse_addmat    (NumDof+6,NumDof+6,CQQ,ctot,Ctotal)

      ! Spherical Joint - sm:
      ! set damping and stiffness terms columns and rows to zero. this can also
      ! be set to unit (won't change the nature of the equations at the joint).
      if (SphFlag .eqv. .true.) then
        call sparse_set_colrows_zero(sph_rows,ctot,Ctotal)
        call sparse_set_colrows_zero(sph_rows,ktot,Ktotal)
      end if

! Compute Jacobian
      ! Spherical Joint - note:
      ! here for the translational dofs of the hinge/spherical joint we will
      ! solve 1.d0/(beta*dt*dt*mass*acceleration = 0.0
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ktot,Ktotal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,ctot,Ctotal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,mtot,Mtotal,as,Asys,Factor=1.d0/(beta*dt*dt))
      !!!! sm modify
      ! call sparse_print_nonzero(as,Asys)

! Calculation of the correction.
      call lu_sparse(as,Asys,-Qtotal,DQ)

      Q     = Q     + DQ
      dQdt  = dQdt  + gamma/(beta*dt)*DQ
      dQddt = dQddt + 1.d0/(beta*dt*dt)*DQ

      ! Spherical Joint - sm
      ! unsure what's the beast approach. The code works with this part being
      ! commented out.
      !if (SphFlag) then
        !Q(NumDof+1:NumDof+3)=0.d0
        !dQdt(NumDof+1:NumDof+3)=0.d0
        !dQddt(NumDof+1:NumDof+3)=0.d0
      !end if

    end do

! Update nodal positions and velocities on the current converged time step.
    X(:)    = Q(1:NumDof)
    dXdt(:) = dQdt(1:NumDof)
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Postprocessing
    Quat= dQdt(NumDof+7:NumDof+10)
    Cao = xbeam_Rot(Quat)
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = dQdt (NumDof+1:NumDof+6)
        VrelDot(iStep+1,:) = dQddt(NumDof+1:NumDof+6)

! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,dQdt (NumDof+1:NumDof+6))
        VrelDot(iStep+1,:) = matmul(ACoa,dQddt(NumDof+1:NumDof+6))
    end if

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do

 ! check crash - error handlying in python
    if (present(SUCCESS)) then
       if (SUCCESS .eqv. .false.) then
         print *, 'Returning to wrapper with SUCCESS=.false.'
         exit
       end if
     end if

  end do

  deallocate (MSS,MSR,CSS,CSR,CRS,CRR,CQR,CQQ,KSS,KRS)
  deallocate (Asys,Felast,Frigid,Qelast,Qrigid)
  deallocate (Mtotal,Ktotal,Ctotal,Qtotal)
  deallocate (Q,DQ,dQdt,dQddt)
  deallocate (X,dXdt,dXddt)
  deallocate (ListIN,Displ,Veloc)

  return
 end subroutine xbeam_solv_couplednlndyn



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_SOLV_COUPLEDNLNDYN
!
!-> Description:
!
!    Nonlinear dynamic solution of multibeam problem under applied forces,
!    coupled with rigid-body dynamics suitable for optimal control
!
!-> Remarks.-
!    Main changes in respect to xbeam_solv_couplednlndyn are:
!       - Fa, Ftime substituted by Fdyn
!       - SUCCESS flag variable added
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_solv_couplednlndyn_opt_control (iOut,NumDof,Time,Elem,Node,F0,Fdyn,               &
&                                     Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,  &
&                                     PosDotDefor,PsiDotDefor,DynOut,Options,SUCCESS)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_out
  use lib_sparse
  use lib_xbeam
  use lib_lu
  !use interface_lapack
  use cbeam3_asbly
  use cbeam3_solv
  use xbeam_asbly

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fdyn      (:,:,:) ! applied dynamic force of size (NumNodes, 6, NumSteps+1)
  !real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  !real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(inout):: Quat      (4)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

 logical, intent(inout), optional :: SUCCESS  ! Variable to allow python wrapper to handle ecceptions.
                                              ! If the solution does not converge, the variable is set to .false.


! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)   ! Generalized coordinates and derivatives for structure.
  real(8),allocatable:: X(:)
  real(8),allocatable:: dQdt(:),dQddt(:)   ! Generalized coordinates and derivatives for coupled system.
  real(8),allocatable:: Q(:), DQ(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for structure system matrices.

  integer:: as,cs,ks,ms,fs,cr,kr,mr,fr,ctot,ktot,mtot ! size of sparse matrices
  type(sparse),allocatable:: Asys(:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: CSS(:)     ! Sparse damping matrix.
  type(sparse),allocatable:: KSS(:)     ! Elast stiffness matrix in sparse storage.
  type(sparse),allocatable:: MSS(:)     ! Elast mass matrix in sparse storage.
  type(sparse),allocatable:: Felast(:)  ! Applied external forces on structure
  real(8),allocatable::      Qelast(:)  ! Elast vector of discrete generalize forces.
  real(8),allocatable::      MSR(:,:)   ! Mass and damping from the motions of reference system.
  real(8),allocatable::      CSR(:,:)   ! Mass and damping from the motions of reference system.

  ! Define variables for rigid system matrices.
  type(sparse),allocatable:: CRS(:)      ! rigid Sparse damping matrix.
  type(sparse),allocatable:: KRS(:)      ! rigid stiffness matrix in sparse storage.
  type(sparse),allocatable:: MRS(:)      ! rigid mass matrix in sparse storage.
  type(sparse),allocatable:: Frigid(:)   ! rigid matrix of applied forces in sparse format
  real(8),allocatable::      Qrigid(:)   ! rigid vector of discrete generalize forces.
  real(8),allocatable::      MRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CQR(:,:),CQQ(:,:)  ! Tangent matrices from linearisation of quaternion equation.

  ! Define variables for rigid-body motion
  real(8),allocatable:: Cao (:,:)                   ! Rotation operator from reference to inertial frame
  real(8),allocatable:: ACoa(:,:)                   ! Rotation operator from reference to inertial frame

  ! Define variables for complete system matrices.
  type(sparse),allocatable:: Ctotal(:)    ! Total Sparse damping matrix.
  type(sparse),allocatable:: Ktotal(:)    ! Total stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mtotal(:)    ! Total mass matrix in sparse storage.
  real(8),allocatable::      Qtotal(:)    ! Total vector of discrete generalize forces.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! =Define variables for sperical joint BC
  logical :: SphFlag
  integer :: sph_rows(3) ! rows to be modified to include a spherical joint BC

  !!!! hinge in arbitrary position:
  !real(8) :: possph0(3)  ! Initial Position of sperical joint in global FoR
  !real(8) :: Rsph0(3)   ! Position of sperical joint in FoR A. Constant in time
  !real(8) :: posA(3),velA(3),accA(3) ! Position, velocity, accelerations of origin of FoR A when spherical joint is applied

  ! Parameters to Check Convergence
  logical :: converged       = .false.
  logical :: passed_delta    = .false.! true if the subiteration (newton) converged according to delta check
  logical :: passed_res      = .false.! true if the subiteration (newton) converged according to residual check
  logical :: passed_resfrc   = .false.! true if the subiteration (newton) converged according to residual check (forces only)
  logical :: passed_resmmt   = .false.! true if the subiteration (newton) converged according to residual check (Moments only)
  logical :: passed_err      = .false.! true if the subiteration (newton) converged according to error check
  logical :: passed_errpos   = .false.! true if the subiteration (newton) converged according to error check (position only)
  logical :: passed_errpsi   = .false.! true if the subiteration (newton) converged according to error check (rotations only)

  ! variables for residual check
  real(8) :: Fsc, Msc                   ! scaling factor for forces and moments
  real(8) :: QglFrc(NumDof/2), QglMmt(NumDof/2) ! residual related to forces and moments
  real(8) :: Res, Res0                  ! current and  initial residual
  real(8) :: ResFrc, ResFrc0, ResFrcRel ! absolute, initial and relative residual
  real(8) :: ResMmt, ResMmt0, ResMmtRel ! absolute, initial and relative residual
  real(8) :: TaRes, TaResFrc, TaResMmt  ! absolute tolerance for Res, ResMmt, ResFrc

 ! variables for error based check
  real(8)   :: Possc, Psisc         ! scaling factor for position and rotations (Psi and quaternion)
  real(8)   :: DeltaPos(NumDof/2)   ! delta displacements at current and previous iteration
  real(8)   :: DeltaPsi(NumDof/2)   ! delta rotations at current and previous iteration
  real(8)   :: ErrX, ErrPos, ErrPsi ! Error extimation for all DoF, displacements and rotations
  real(8)   :: DX_now, DX_old       ! Norm of DeltaX at current and old iteration
  real(8)   :: DPos_now, DPos_old   ! Norm of translational dofs of DeltaX at current and old iteration
  real(8)   :: DPsi_now, DPsi_old   ! Norm of rotational dofs of DeltaX at current and old iteration
  real(8)   :: TaX, TaPos, TaPsi    ! Absolute tolerance for DeltaX, DeltaPos and DeltaPsi

  ! Determine scaling factors for convergence test (absolute tolerances)
  Psisc = 1.0_8 ! used for quaternion scaling as well
  Possc = maxval(abs(Coords))
  Fsc = 1.0_8; Msc = 1.0_8;
  do iStep=1,size(Time)-1
    Fsc = max( maxval(abs( F0(:,1:3)+Fdyn(:,1:3,iStep+1) )), Fsc);
    Msc = max( maxval(abs( F0(:,4:6)+Fdyn(:,4:6,iStep+1) )), Msc);
  end do

  ! Initialize (sperical joint) sm
  NumN=size(Node)
  SphFlag=.false.
  do k=1,NumN
      if (Node(k)%Sflag == 1) then
        SphFlag=.true.
        !Rsph0=Coords(k,:) ! Local position of spherical joint
        !possph0=Rsph0     ! assume FoRs A and G are coincident at t=0
      end if
  end do
  sph_rows = (/1,2,3/)+NumDof

  ! Initialise
  NumN=size(Node)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do

  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys(4*DimMat*NumDof)); call sparse_zero (as,Asys)

  allocate (MSS(DimMat*NumDof));    call sparse_zero (ms,MSS)
  allocate (CSS(DimMat*NumDof));    call sparse_zero (cs,CSS)
  allocate (KSS(DimMat*NumDof));    call sparse_zero (ks,KSS)
  allocate (Felast(DimMat*NumDof)); call sparse_zero (fs,Felast)
  allocate (Qelast(NumDof));        Qelast = 0.d0
  allocate (MSR   (NumDof,6));      MSR    = 0.d0
  allocate (CSR   (NumDof,6));      CSR    = 0.d0

  allocate (MRS(DimMat*NumDof));    call sparse_zero (mr,MRS)
  allocate (CRS(DimMat*NumDof));    call sparse_zero (cr,CRS)
  allocate (KRS(DimMat*NumDof));    call sparse_zero (kr,KRS)
  allocate (Frigid(DimMat*NumDof)); call sparse_zero (fr,Frigid)
  allocate (Qrigid(6));             Qrigid = 0.d0
  allocate (MRR   (6,6));           MRR    = 0.d0
  allocate (CRR   (6,6));           CRR    = 0.d0

  allocate (CQR(4,6));    CQR   = 0.d0
  allocate (CQQ(4,4));    CQQ   = 0.d0

  allocate (Mtotal(2*DimMat*NumDof));   call sparse_zero (mtot,Mtotal)
  allocate (Ctotal(2*DimMat*NumDof));   call sparse_zero (ctot,Ctotal)
  allocate (Ktotal(2*DimMat*NumDof));   call sparse_zero (ktot,Ktotal)
  allocate (Qtotal(NumDof+6+4));        Qtotal= 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Updated state vector with rigid body states and quaternions
  allocate (Q     (NumDof+6+4)); Q      = 0.d0
  allocate (dQdt  (NumDof+6+4)); dQdt   = 0.d0
  allocate (dQddt (NumDof+6+4)); dQddt  = 0.d0
  allocate (DQ    (NumDof+6+4)); DQ     = 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  allocate(Cao      (3,3));   Cao       = Unit
  allocate(ACoa     (6,6));   ACoa      = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

  ! sm exceptions handling for python wrapper
  if (present(SUCCESS)) then
      SUCCESS=.true.
  end if

  ! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)

  Q(1:NumDof)             = X(:)
  Q(NumDof+1:NumDof+6)    = 0
  dQdt(1:NumDof)          = dXdt(:)
  dQdt(NumDof+1:NumDof+6) = Vrel(1,:)
  dQdt(NumDof+7:NumDof+10)= Quat

  ! (spherical joint - sm)
  if (SphFlag) then
    dQdt(NumDof+1:NumDof+3) = 0.d0
  end if

  Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

  ACoa(1:3,1:3) = transpose(Cao)
  ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
  if (Options%OutInaframe) then
      Vrel   (1,:) = dQdt (NumDof+1:NumDof+6)
      VrelDot(1,:) = dQddt(NumDof+1:NumDof+6)

! Export velocities and accelerations in inertial frame
  else
      Vrel   (1,:) = matmul(ACoa,dQdt (NumDof+1:NumDof+6))
      VrelDot(1,:) = matmul(ACoa,dQddt(NumDof+1:NumDof+6))
  end if

  do k=1,NumN
      DynOut(k,:) = PosDefor(k,:)
  end do

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                            F0+Fdyn(:,:,1),dQdt(NumDof+1:NumDof+6),            &
&                            0.d0*dQddt(NumDof+1:NumDof+6),                     &
&                            ms,MSS,MSR,cs,CSS,CSR,ks,KSS,fs,Felast,Qelast,Options,Cao)

  Qelast= Qelast - sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Fdyn(:,:,1),NumDof,Filter=ListIN))

  call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,    &
&                           dQdt(NumDof+1:NumDof+6),                            &
&                           0.d0*dQddt(NumDof+1:NumDof+6),                      &
&                           dQdt(NumDof+7:NumDof+10),                                                                       &
&                           mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

  ! Sperical Joint - sm:
  !  Qrigid: Total generalized forces on the element: gyroscopic, lumped masses.
  !          This is non-linear and acts like residual
  !  Frigid: influence coeff. matrix for applied forces.
  ! All terms related to the translational dof and that will influence the RHS
  ! of Newmark-beta time stepping are set to zero.
  if (SphFlag) then
    call sparse_set_rows_zero((/1,2,3/),fr,Frigid)
    Qrigid(1:3)=0.d0
  end if
  Qrigid= Qrigid - sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Fdyn(:,:,1),NumDof+6))

  !!!print *, 'MRR:'
  !!!write (*,'(5X,1P6E12.4)') MRR

! Assemble coupled system
  Qtotal(1:NumDof)          = Qelast
  Qtotal(NumDof+1:NumDof+6) = Qrigid
  Qtotal(NumDof+7:NumDof+10)= matmul(CQQ,dQdt(NumDof+7:NumDof+10))

  call sparse_addsparse (0,0,ms,MSS,mtot,Mtotal)
  call sparse_addmat    (0,NumDof,MSR,mtot,Mtotal)
  call sparse_addmat    (NumDof,0,transpose(MSR),mtot,Mtotal)
  call sparse_addmat    (NumDof,NumDof,MRR,mtot,Mtotal)
  call sparse_addmat    (NumDof+6,NumDof+6,Unit4,mtot,Mtotal)

! Spherical Joint - sm:
! unit diag term in Mtotal to the translational dof of the a FoR. The aim is to
! implement M*dBetadt = 0
! Note: also the columns deleted (though unecessary if the beam root is not
! moving) to avoid numerical error. If a moving hinge/spherical joint is implemented,
! the columns related to the translational dof can't be set to zero!
  if (SphFlag) then
    call sparse_set_colrows_zero(sph_rows,mtot,Mtotal)
    call sparse_set_rows_unit(sph_rows,mtot,Mtotal)
  end if


! Solve matrix system
!  call lapack_sparse (mtot,Mtotal,-Qtotal,dQddt)
  call lu_sparse (mtot,Mtotal,-Qtotal,dQddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    if (Options%PrintInfo) then
      write (*,'(X)')
      write (*,'(5X,A,$)') trim(Text)
    end if

! Predictor step.
    Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
    dQdt = dQdt + (1.d0-gamma)*dt*dQddt
    dQddt= 0.d0

    ! Spherical Joint - sm
    ! Predictor step is forced to be zero for translational dofs.
    if (SphFlag) then
      Q(NumDof+1:NumDof+3)=0.d0
      dQdt(NumDof+1:NumDof+3)=0.d0
      dQddt(NumDof+1:NumDof+3)=0.d0
    end if

! Iteration until convergence.
  converged=.false.
! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      !Spherical Joint - sm
      ! reminder of possible issues
      if (Iter.gt.Options%MaxIterations) then
        print *, '\N success=', SUCCESS
        print *, 'Reminders:'
        print *, '1. Old convergence criteria still used. Consider upgrading as per static solver.'
        print *, '2. Spherical Joint, not hinge!'
        !!! sm: stop commented to allow python wrapper to handle exceptions
        if (present(SUCCESS)) then
            SUCCESS=.false.
            ! always print last iteration delta if crash occurrs
            if (.not.(Options%PrintInfo)) then
              write (*,'(5X,A,$)') trim(Text)
              write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qtotal))
            end if
            print *, 'Solution did not converge (18235)'
            exit
            ! stop the iterations
        else
            STOP 'Solution did not converge (18235)'
        end if
      end if

! Update nodal positions and velocities .
      X   (:) = Q   (1:NumDof)
      dXdt(:) = dQdt(1:NumDof)
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Update quaternions with new states
      Cao = xbeam_Rot(dQdt(NumDof+7:NumDof+10))

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qelast = 0.d0
      MSR    = 0.d0
      CSR    = 0.d0
      call sparse_zero (ms,MSS)
      call sparse_zero (cs,CSS)
      call sparse_zero (ks,KSS)
      call sparse_zero (fs,Felast)
      Qrigid = 0.d0
      MRR    = 0.d0
      CRR    = 0.d0
      CQR    = 0.d0
      CQQ    = 0.d0
      call sparse_zero (mr,MRS)
      call sparse_zero (cr,CRS)
      call sparse_zero (kr,KRS)
      call sparse_zero (fr,Frigid)
      Qtotal = 0.d0
      call sparse_zero (mtot,Mtotal)
      call sparse_zero (ctot,Ctotal)
      call sparse_zero (ktot,Ktotal)

      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                                F0+Fdyn(:,:,iStep+1),dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),                    &
&                                ms,MSS,MSR,cs,CSS,CSR,ks,KSS,fs,Felast,Qelast,Options,Cao)

      call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,   &
&                               dQdt(NumDof+1:NumDof+6),0.d0*dQddt(NumDof+1:NumDof+6),dQdt(NumDof+7:NumDof+10),                &
&                               mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

      ! Spherical Joint - sm:
      ! for spherical joint doesn't matter the orientation of a, we want all the
      ! global external force to be zero.
      ! Note: for hinge BCs no rotations should be required to understand which
      ! moment has to be set to zero.
      if (SphFlag) then
        call sparse_set_rows_zero((/1,2,3/),fr,Frigid)
        Qrigid(1:3)=0.d0
      end if

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qelast)))

! Compute the residual.
      Qelast = Qelast - sparse_matvmul(fs,Felast,NumDof,fem_m2v(F0+Fdyn(:,:,iStep+1),NumDof,Filter=ListIN))
      Qrigid = Qrigid - sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Fdyn(:,:,iStep+1),NumDof+6))

      Qtotal(1:NumDof)          = Qelast
      Qtotal(NumDof+1:NumDof+6) = Qrigid
      Qtotal(NumDof+7:NumDof+10)= matmul(CQQ,dQdt(NumDof+7:NumDof+10))

      call sparse_addsparse (0,0,ms,MSS,mtot,Mtotal)
      call sparse_addmat    (0,NumDof,MSR,mtot,Mtotal)
!      call sparse_addsparse (NumDof,0,mr,MRS,mtot,Mtotal)
      call sparse_addmat    (NumDof,0,transpose(MSR),mtot,Mtotal)
      call sparse_addmat    (NumDof,NumDof,MRR,mtot,Mtotal)
      call sparse_addmat    (NumDof+6,NumDof+6,Unit4,mtot,Mtotal)

      ! Spherical Joint - sm:
      ! set Mass matrix to unit. If the hinge/spherical joint translates, the
      ! coulmn part cannot be set to zero (see also above).
      if (SphFlag) then
        call sparse_set_colrows_zero(sph_rows,mtot,Mtotal)
        call sparse_set_rows_unit(sph_rows,mtot,Mtotal)
      end if

      Qtotal= Qtotal + sparse_matvmul(mtot,Mtotal,NumDof+6+4,dQddt)

      ! ------------------------------------------------------ Check convergence.

      if (Options%PrintInfo) then
        ! delta check (original code)
        !write (*,'(2X,1PE10.3,$)') maxval(abs(Qtotal))
        write (*,'(5X,A,I4,A,1PE12.3,$)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qtotal))
      end if

      if (maxval(abs(Qtotal)).lt.MinDelta) then
          converged=.true. ! sm
      end if

      ! exit do loop: kept outside in case of other criteria will be introduced
      if (converged.eqv..true.) then
        exit
      end if

! ----------------------------------------------------- end check of convergence

! Compute total damping and stiffness matrices linear system
      call sparse_addsparse (0,0,cs,CSS,ctot,Ctotal)
      call sparse_addmat    (0,NumDof,CSR,ctot,Ctotal)
      call sparse_addsparse (NumDof,0,cr,CRS,ctot,Ctotal)
      call sparse_addmat    (NumDof,NumDof,CRR,ctot,Ctotal)

      call sparse_addsparse (0,0,ks,KSS,ktot,Ktotal)
      call sparse_addsparse (NumDof,0,kr,KRS,ktot,Ktotal)

! Contribution of quaternions to damping matrix
      call sparse_addmat    (NumDof+6,NumDof  ,CQR,ctot,Ctotal)
      call sparse_addmat    (NumDof+6,NumDof+6,CQQ,ctot,Ctotal)


      ! Spherical Joint - sm:
      ! set damping and stiffness terms columns and rows to zero. this can also
      ! be set to unit (won't change the nature of the equations at the joint).
      if (SphFlag .eqv. .true.) then
        call sparse_set_colrows_zero(sph_rows,ctot,Ctotal)
        call sparse_set_colrows_zero(sph_rows,ktot,Ktotal)
      end if

! Compute Jacobian
      ! Spherical Joint - note:
      ! here for the translational dofs of the hinge/spherical joint we will
      ! solve 1.d0/(beta*dt**2*mass*acceleration) = 0.0
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ktot,Ktotal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,ctot,Ctotal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,mtot,Mtotal,as,Asys,Factor=1.d0/(beta*dt*dt))
      !!!! sm modify
      ! call sparse_print_nonzero(as,Asys)

! Calculation of the correction.
      !call lapack_sparse (as,Asys,-Qtotal,DQ)
      call lu_sparse (as,Asys,-Qtotal,DQ)

      Q     = Q     + DQ
      dQdt  = dQdt  + gamma/(beta*dt)*DQ
      dQddt = dQddt + 1.d0/(beta*dt*dt)*DQ

      ! Spherical Joint - sm
      ! unsure what's the beast approach. The code works with this part being
      ! commented out.
      ! Without hardcoding the terms to zero, we make sure the system, durint the
      ! iteration, will converge to the 0.0 solution. Also, non zero values will
      ! not affect other dofs (columns zero in all matrices).
      !if (SphFlag) then
        !Q(NumDof+1:NumDof+3)=0.d0
        !dQdt(NumDof+1:NumDof+3)=0.d0
        !dQddt(NumDof+1:NumDof+3)=0.d0
      !end if

    end do

! Update nodal positions and velocities on the current converged time step.
    X(:)    = Q(1:NumDof)
    dXdt(:) = dQdt(1:NumDof)
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Postprocessing
    Quat= dQdt(NumDof+7:NumDof+10)
    Cao = xbeam_Rot(Quat)
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = dQdt (NumDof+1:NumDof+6)
        VrelDot(iStep+1,:) = dQddt(NumDof+1:NumDof+6)

! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,dQdt (NumDof+1:NumDof+6))
        VrelDot(iStep+1,:) = matmul(ACoa,dQddt(NumDof+1:NumDof+6))
    end if

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do



 ! check crash - error handlying in python
    if (present(SUCCESS)) then
       if (SUCCESS .eqv. .false.) then
         print *, 'Returning to wrapper with SUCCESS=.false.'
         exit
       end if
     end if


  end do

  deallocate (MSS,MSR,CSS,CSR,CRS,CRR,CQR,CQQ,KSS,KRS)
  deallocate (Asys,Felast,Frigid,Qelast,Qrigid)
  deallocate (Mtotal,Ktotal,Ctotal,Qtotal)
  deallocate (Q,DQ,dQdt,dQddt)
  deallocate (X,dXdt,dXddt)
  deallocate (ListIN,Displ,Veloc)

  return
 end subroutine xbeam_solv_couplednlndyn_opt_control




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_SOLV_RIGIDlndyn
!
!-> Description:
!
!    Linear rigid-body dynamic solution of multibeam problem for given applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_solv_rigidlndyn  (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,             &
&                                   Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,    &
&                                   PosDotDefor,PsiDotDefor,Options)
  use lib_rot
  use lib_fem
  use lib_sparse
  use lib_out
  use lib_xbeam
  use lib_lu
  use cbeam3_asbly
  use cbeam3_solv
  use xbeam_asbly

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Quat      (4)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Undeformed coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Undeformed CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Initial/final position vector of grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep                          ! Current time step.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8)               :: DQ(10),DQDt(10),DQDDt(10)
  integer,allocatable   :: ListIN     (:)     ! List of independent nodes.
  real(8),allocatable   :: DeltaPos   (:,:)   ! Initial/final position vector of grid points
  real(8),allocatable   :: DeltaPsi   (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),allocatable   :: DeltaPosDot(:,:)   ! Current time derivatives of the coordinates of the grid points
  real(8),allocatable   :: DeltaPsiDot(:,:,:) ! Current time derivatives of the CRV of the nodes in the elements.

  ! Define vaiables for structural system matrices
  integer                   :: as,cr,kr,mr,fr,ctot,ktot,mtot
  type(sparse),allocatable  :: Asys(:)            ! System matrix for implicit Newmark method.
  real(8),allocatable       :: MSR(:,:)           ! Mass and damping from the motions of reference system.
  real(8),allocatable       :: CSR(:,:)           ! Mass and damping from the motions of reference system.
  real(8)                   :: CQR(4,6),CQQ(4,4)  ! Tangent matrices from linearisation of quaternion equation.

  ! Define variables for rigid system matrices.
  type(sparse),allocatable  :: CRS(:)      ! rigid Sparse damping matrix.
  type(sparse),allocatable  :: KRS(:)      ! rigid stiffness matrix in sparse storage.
  type(sparse),allocatable  :: MRS(:)      ! rigid mass matrix in sparse storage.
  type(sparse),allocatable  :: Frigid(:)   ! rigid matrix of applied forces in sparse format
  real(8),allocatable       :: Qrigid(:)   ! rigid vector of discrete generalize forces.
  real(8)                   :: MRR(6,6)    ! rigid Mass and damping from the motions of reference system.
  real(8)                   :: CRR(6,6)    ! rigid Mass and damping from the motions of reference system.

  ! Define variables for rigid-body motion
  real(8)               :: Cao (3,3)                   ! Rotation operator from reference to inertial frame
  real(8)               :: ACoa(6,6)                   ! Rotation operator from reference to inertial frame

  ! Define variables for complete system matrices.
  type(sparse),allocatable  :: Ctotal(:)    ! Total Sparse damping matrix.
  type(sparse),allocatable  :: Ktotal(:)    ! Total stiffness matrix in sparse storage.
  type(sparse),allocatable  :: Mtotal(:)    ! Total mass matrix in sparse storage.
  real(8)                   :: Qtotal(10)   ! Total vector of discrete generalize forces.

  character(len=80)     :: Text          ! Text with current time information to print out.
  real(8),allocatable   :: Displ(:,:),Veloc(:,:)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (4*DimMat*NumDof)); call sparse_zero (as,Asys)

  allocate (MSR   (NumDof,6));      MSR = 0.d0
  allocate (CSR   (NumDof,6));      CSR = 0.d0

  allocate (MRS(DimMat*NumDof));    call sparse_zero (mr,MRS)
  allocate (CRS(DimMat*NumDof));    call sparse_zero (cr,CRS)
  allocate (KRS(DimMat*NumDof));    call sparse_zero (kr,KRS)
  allocate (Frigid(DimMat*NumDof)); call sparse_zero (fr,Frigid)
  allocate (Qrigid(6));             Qrigid = 0.d0

  MRR = 0.d0; CRR = 0.d0

  allocate (Mtotal(DimMat*NumDof)); call sparse_zero (mtot,Mtotal)
  allocate (Ctotal(DimMat*NumDof)); call sparse_zero (ctot,Ctotal)
  allocate (Ktotal(DimMat*NumDof)); call sparse_zero (ktot,Ktotal)
  Qtotal = 0.d0

! Updated state vector with rigid body states
  DQ = 0.d0; DQDt = 0.d0; DQDDt = 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  DQDt(1:6) = Vrel(1,:)
  DQDt(7:10)= Quat

  Cao = (xbeam_Rot(DQDt(7:10)))

! Compute system information at initial condition.
  allocate (Displ      (NumN,            6)); Displ       = 0.d0
  allocate (Veloc      (NumN,            6)); Veloc       = 0.d0
  allocate (DeltaPos   (NumN,            3)); DeltaPos    = 0.d0
  allocate (DeltaPsi   (NumE(1),MaxElNod,3)); DeltaPsi    = 0.d0
  allocate (DeltaPosDot(NumN,            3)); DeltaPosDot = 0.d0
  allocate (DeltaPsiDot(NumE(1),MaxElNod,3)); DeltaPsiDot = 0.d0

  call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,    &
&                           0.d0*PosDefor,0.d0*PsiDefor,DQDt(1:6),0.d0*DQDDt(1:6),DQDt(7:10),   &
&                           mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

! Find initial acceleration.
  Qtotal(1:6) = sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Ftime(1)*Fa,NumDof+6))
  Qtotal(7:10)= matmul(CQQ,DQDt(7:10))

  call sparse_addmat (0,0,MRR,mtot,Mtotal)
  call sparse_addmat (6,6,Unit4,mtot,Mtotal)

  call lu_sparse(mtot,Mtotal,Qtotal,DQDDt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt     =Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A)') trim(Text)

! Predictor step.
    DQ   = DQ   + dt*DQDt + (0.5d0-beta)*dt*dt*DQDDt
    DQDt = DQDt + (1-gamma)*dt*DQDDt

! Update quaternions with new states
    Cao = xbeam_Rot(DQDt(7:10))

! Compute system functionals and matrices. Only updating Frigid, CQR and CQQ to update the orientation of the body-fixed FoR
! Need to make sure that CRR is not overwritten.
    call sparse_zero (ctot,Ctotal); CSR = 0.d0; call sparse_zero (cr,CRS); CQR = 0.d0; CQQ = 0.d0
    call sparse_zero (fr,Frigid)

    call xbeam_asbly_orient (Elem,Node,PosDefor,PsiDefor,DQDt(1:6),DQDt(7:10),CQR,CQQ,fr,Frigid,Options,Cao)

! Compute right-hand side of the equation.
    Qtotal(1:6) = sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof+6))
    Qtotal(7:10)= matmul(CQQ,DQDt(7:10))
 
    ! Compute total damping matrix of linear system
    call sparse_addmat (0,0,CRR,ctot,Ctotal)

    ! Contribution of quaternions to damping matrix
    call sparse_addmat (6,0,CQR,ctot,Ctotal)
    call sparse_addmat (6,6,CQQ,ctot,Ctotal)

    Qtotal= Qtotal - sparse_matvmul(ctot,Ctotal,10,DQDt)

! Compute left-hand side of the equation (if dt=constant then Asys=constant too, but this was not assumed here).
    call sparse_zero (as,Asys)
    call sparse_addsparse(0,0,mtot,Mtotal,as,Asys,Factor=1.d0)
    call sparse_addsparse(0,0,ctot,Ctotal,as,Asys,Factor=gamma*dt)

! Solve equation.
    call lu_sparse(as,Asys,Qtotal,DQDDt)

    DQ    = DQ   + beta *dt*dt*DQDDt
    DQDt  = DQDt + gamma*dt   *DQDDt

! Postprocessing
    Cao = xbeam_Rot(DQDt(7:10))
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = DQDt (1:6)
        VrelDot(iStep+1,:) = DQDDt(1:6)
! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,DQDt (1:6))
        VrelDot(iStep+1,:) = matmul(ACoa,DQDDt(1:6))
    end if

  end do

! Write information at last time step.
  PosDefor= 0.d0
  PsiDefor= 0.d0

  deallocate (ListIn,Asys,MSR,CSR)
  deallocate (Displ,Veloc,DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)
  return
 end subroutine xbeam_solv_rigidlndyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine XBEAM_SOLV_RIGIDNLNDYN
!
!-> Description:
!
!    Nonlinear rigid-body dynamic solution of multibeam problem under applied forces
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine xbeam_solv_rigidnlndyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,                 &
&                                   Vrel,VrelDot,Quat,Coords,Psi0,PosDefor,PsiDefor,        &
&                                   PosDotDefor,PsiDotDefor,DynOut,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_out
  use lib_sparse
  use lib_xbeam
  use lib_lu
  use cbeam3_asbly
  use cbeam3_solv
  use xbeam_asbly

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(out)  :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(out)  :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(inout):: Quat      (4)     ! Quaternions to describes motion of reference system.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  type(xbopts), intent(in)   :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dQdt(:),dQddt(:)   ! Generalized coordinates and derivatives for coupled system.
  real(8),allocatable:: Q(:), DQ(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define vaiables for system matrices
  integer:: as,cr,kr,mr,fr,ctot,ktot,mtot
  type(sparse),allocatable:: Asys(:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: CRS(:)      ! rigid Sparse damping matrix.
  type(sparse),allocatable:: KRS(:)      ! rigid stiffness matrix in sparse storage.
  type(sparse),allocatable:: MRS(:)      ! rigid mass matrix in sparse storage.
  type(sparse),allocatable:: Frigid(:)   ! rigid matrix of applied forces in sparse format
  real(8),allocatable::      Qrigid(:)   ! rigid vector of discrete generalize forces.
  real(8),allocatable::      MRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CRR(:,:)    ! rigid Mass and damping from the motions of reference system.
  real(8),allocatable::      CQR(:,:),CQQ(:,:)  ! Tangent matrices from linearisation of quaternion equation.

  ! Define variables for rigid-body motion
  real(8),allocatable:: Cao (:,:)                   ! Rotation operator from reference to inertial frame
  real(8),allocatable:: ACoa(:,:)                   ! Rotation operator from reference to inertial frame

  ! Define variables for complete system matrices.
  type(sparse),allocatable:: Ctotal(:)    ! Total Sparse damping matrix.
  type(sparse),allocatable:: Ktotal(:)    ! Total stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mtotal(:)    ! Total mass matrix in sparse storage.
  real(8),allocatable::      Qtotal(:)    ! Total vector of discrete generalize forces.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

! Initialize.
  NumN=size(Node)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=1.d0/2.d0+Options%NewmarkDamp
  beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Asys   (4*DimMat*NumDof)); call sparse_zero (as,Asys)

  allocate (MRS(DimMat*NumDof));    call sparse_zero (mr,MRS)
  allocate (CRS(DimMat*NumDof));    call sparse_zero (cr,CRS)
  allocate (KRS(DimMat*NumDof));    call sparse_zero (kr,KRS)
  allocate (Frigid(DimMat*NumDof)); call sparse_zero (fr,Frigid)
  allocate (Qrigid(6));             Qrigid = 0.d0
  allocate (MRR   (6,6));           MRR    = 0.d0
  allocate (CRR   (6,6));           CRR    = 0.d0

  allocate (CQR(4,6));    CQR   = 0.d0
  allocate (CQQ(4,4));    CQQ   = 0.d0

  allocate (Mtotal(DimMat*NumDof)); call sparse_zero (mtot,Mtotal)
  allocate (Ctotal(DimMat*NumDof)); call sparse_zero (ctot,Ctotal)
  allocate (Ktotal(DimMat*NumDof)); call sparse_zero (ktot,Ktotal)
  allocate (Qtotal(6+4));           Qtotal= 0.d0

! Updated state vector with rigid body states
  allocate (Q     (6+4)); Q      = 0.d0
  allocate (dQdt  (6+4)); dQdt   = 0.d0
  allocate (dQddt (6+4)); dQddt  = 0.d0
  allocate (DQ    (6+4)); DQ     = 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  allocate(Cao      (3,3));   Cao       = Unit
  allocate(ACoa     (6,6));   ACoa      = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

  dQdt(1:6) = Vrel(1,:)
  dQdt(7:10)= Quat

  Cao = xbeam_Rot(dQdt(7:10))

  do k=1,NumN
      DynOut(k,:) = PosDefor(k,:)
  end do

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,    &
&                           0.d0*PosDefor,0.d0*PsiDefor,dQdt(1:6),0.d0*dQddt(1:6),dQdt(7:10),   &
&                           mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

  write (*,'(1P6E12.4)') MRR

  Qrigid= Qrigid - sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Ftime(1)*Fa,NumDof+6))

! Assemble coupled system
  Qtotal(1:6) = Qrigid
  Qtotal(7:10)= matmul(CQQ,dQdt(7:10))

  call sparse_addmat (0,0,MRR,mtot,Mtotal)
  call sparse_addmat (6,6,Unit4,mtot,Mtotal)

  call lu_sparse(mtot,Mtotal,-Qtotal,dQddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A,$)') trim(Text)

! Predictor step.
    Q    = Q + dt*dQdt + (0.5d0-beta)*dt*dt*dQddt
    dQdt = dQdt + (1.d0-gamma)*dt*dQddt
    dQddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) STOP 'Solution did not converge (18235)'

! Update quaternions with new states
      Cao = xbeam_Rot(dQdt(7:10))

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qrigid      = 0.d0
      MRR   = 0.d0
      CRR   = 0.d0
      CQR   = 0.d0
      CQQ   = 0.d0
      call sparse_zero (cr,CRS)
      call sparse_zero (kr,KRS)
      call sparse_zero (fr,Frigid)
      Qtotal = 0.d0
      call sparse_zero (mtot,Mtotal)
      call sparse_zero (ctot,Ctotal)
      call sparse_zero (ktot,Ktotal)

      call xbeam_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                &
&                               0.d0*PosDefor,0.d0*PsiDefor,dQdt(1:6),0.d0*dQddt(1:6),(/1.d0,0.d0,0.d0,0.d0/),  &
&                               mr,MRS,MRR,cr,CRS,CRR,CQR,CQQ,kr,KRS,fr,Frigid,Qrigid,Options,Cao)

! Assemble total matrices
      Qtotal(1:6) = Qrigid-sparse_matvmul(fr,Frigid,6,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof+6))
      Qtotal(7:10)= matmul(CQQ,dQdt(7:10))

      call sparse_addmat (0,0,MRR,mtot,Mtotal)
      call sparse_addmat (6,6,Unit4,mtot,Mtotal)

! Compute the residual.
      Qtotal= Qtotal + sparse_matvmul(mtot,Mtotal,6+4,dQddt)

      ! Check convergence.
      if (maxval(abs(Qtotal)).lt.Options%MinDelta) then
        write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qtotal))
        exit
      end if

      ! Compute total damping matrix of linear system
      call sparse_addmat    (0,0,CRR,ctot,Ctotal)

      ! Contribution of quaternions to damping matrix
      call sparse_addmat    (6,0,CQR,ctot,Ctotal)
      call sparse_addmat    (6,6,CQQ,ctot,Ctotal)

      ! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ctot,Ctotal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,mtot,Mtotal,as,Asys,Factor=1.d0/(beta*dt*dt))

! Calculation of the correction.
      call lu_sparse(as,Asys,-Qtotal,DQ)

      Q     = Q     + DQ
      dQdt  = dQdt  + gamma/(beta*dt)*DQ
      dQddt = dQddt + 1.d0/(beta*dt*dt)*DQ

    end do

! Update nodal positions and velocities on the current converged time step.
    Cao = xbeam_Rot(dQdt(7:10))
    ACoa(1:3,1:3) = transpose(Cao)
    ACoa(4:6,4:6) = transpose(Cao)

! Export velocities and accelerations in body frame
    if (Options%OutInaframe) then
        Vrel   (iStep+1,:) = dQdt (1:6)
        VrelDot(iStep+1,:) = dQddt(1:6)
! Export velocities and accelerations in inertial frame
    else
        Vrel   (iStep+1,:) = matmul(ACoa,dQdt (1:6))
        VrelDot(iStep+1,:) = matmul(ACoa,dQddt(1:6))
    end if

    do k=1,NumN
        DynOut((iStep)*NumN+k,:) = PosDefor(k,:)
    end do

  end do

  deallocate (MRS,CRS,CRR,CQR,CQQ,KRS)
  deallocate (Asys,Frigid,Qrigid)
  deallocate (Mtotal,Ktotal,Ctotal,Qtotal)
  deallocate (Q,DQ,dQdt,dQddt)
  deallocate (ListIN,Displ,Veloc)
  
  return
 end subroutine xbeam_solv_rigidnlndyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module xbeam_solv
