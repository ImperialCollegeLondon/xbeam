!-> Copyright by Imperial College London, 2009
!
!-> Module.- CBEAM3_SOLV Rafa Palacios. 27Aug2009.
!
!-> Language: FORTRAN90, Free format.
!
!-> Description.-
!
!  Solve beam equations.
!
!-> Subroutines.-
!
!    -cbeam3_solv_nlnstatic     : Nonlinear static solution.
!    -cbeam3_solv_linstatic     : Linear static solution.
!    -cbeam3_solv_modal         : Find natural modes.
!    -cbeam3_solv_nlndyn        : Solve nonlinear dynamic problem.
!    -cbeam3_solv_lindyn        : Solve linear dynamic problem.
!    -cbeam3_solv_update_static : Update solution variables in static problem.
!    -cbeam3_solv_update_lindyn : Update solution variables in linear dynamic problem.
!    -cbeam3_solv_state2disp    : Write state vector for displacements/rotations.
!    -cbeam3_solv_disp2state    : Write state vector for displacements/rotations.
!
!-> Remarks.-
!
!-> Modifications.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module cbeam3_solv
  use, intrinsic    :: iso_c_binding
  use debug_utils
  use xbeam_shared
  use                                 :: debug_utils
  use xbeam_asbly
  use lib_rot
  implicit none

! Private variables.
  integer,private,parameter:: DimMat=24    ! Memory index for sparse matrices.

  real(8),private,parameter,dimension(3,3):: Unit= &       ! 4x4 Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
  real(8),private,parameter,dimension(4,4):: Unit4= &       ! 4x4 Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,0.d0,1.d0/),(/4,4/))

 interface cbeam3_solv_nlnstatic
    module procedure :: cbeam3_solv_nlnstatic_old
 end interface cbeam3_solv_nlnstatic

 interface cbeam3_solv_nlndyn
    module procedure :: cbeam3_solv_nlndyn_updated, cbeam3_solv_nlndyn_old
 end interface cbeam3_solv_nlndyn

 !interface cbeam3_solv_modal
!    module procedure :: cbeam3_solv_modal_old, cbeam3_solv_modal_updated
 !end interface cbeam3_solv_modal

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNSTATIC
!
!-> Description:
!
!    Steady-state solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlnstatic_old (NumDof,n_elem,n_node,Elem,Node,AppForces,nodal_gravity_forces,Coords,Psi0, &
&                                  PosDefor,PsiDefor,Options)
  use lib_fem
  use lib_sparse
  use lib_solv
!#ifdef NOLAPACK
  use lib_lu
  use cbeam3_asbly

! I/O Variables.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  integer,      intent(in)   :: n_elem ! Number of independent DoFs.
  integer,      intent(in)   :: n_node ! Number of independent DoFs.
  type(xbelem),intent(in)    :: Elem(n_elem)           ! Element information.
  type(xbnode),intent(in)    :: Node(n_node)           ! Nodal information.
  real(8),      intent(in)   :: AppForces (n_node,6)   ! Applied nodal forces.
  real(8),      intent(OUT)  :: nodal_gravity_forces(n_node, 6)   ! Applied nodal forces.
  real(8),      intent(in)   :: Coords   (n_node,3)    ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0     (n_elem,3,3)  ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor (n_node, 3)    ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor (n_elem, 3, 3)  ! Current CRV of the nodes in the elements.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: Delta                          ! Flag for convergence in iterations.
  integer:: fs                             ! Current storage size of force matrix.
  integer:: ms                             ! Current storage size of mass matrix.
  integer:: iLoadStep                      ! Counter in the load steps.
  integer:: Iter                           ! Counter on iterations.
  integer:: i,j,k                          ! Auxiliary integer.
  integer:: ks                             ! Current storage size of stiffness matrix.

  integer::      ListIN (n_node)    ! List of independent nodes.
  real(8)::      Qglobal(numdof)    ! Global vector of discrete generalize forces.
  real(8)::      DeltaX (numdof)    ! Unknown in the linearized system.
  real(8):: Kglobal(numdof,numdof)    ! Global stiffness matrix in sparse storage.
  real(8):: Fglobal(numdof,numdof)    ! Influence coefficients matrix for applied forces.
  real(8):: Mglobal(numdof + 6,numdof + 6)    ! Global stiffness matrix in sparse storage.

  real(8):: MRR(6, 6)

  ! Parameters to Check Convergence
  logical :: converged = .false.
  logical :: passed_delta = .false.! true if the subiteration (newton) converged according to delta check
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
  real(8)   :: Possc, Psisc   ! scaling factor for position and rotations
  real(8)   :: DeltaPos(NumDof/2)   ! delta displacements at current and previous iteration
  real(8)   :: DeltaPsi(NumDof/2)   ! delta rotations at current and previous iteration
  real(8)   :: ErrX, ErrPos, ErrPsi ! Error extimation for all DoF, displacements and rotations
  real(8)   :: DX_now, DX_old       ! Norm of DeltaX at current and old iteration
  real(8)   :: DPos_now, DPos_old   ! Norm of translational dofs of DeltaX at current and old iteration
  real(8)   :: DPsi_now, DPsi_old   ! Norm of rotational dofs of DeltaX at current and old iteration
  real(8)   :: TaX, TaPos, TaPsi    ! Absolute tolerance for DeltaX, DeltaPos and DeltaPsi

  integer   :: ii
  real(8)   :: gravity_forces(6)

 ! Determine scaling factors for convergence test (absolute tolerances)
!   Psisc = 1.0_8
!   Possc = maxval(abs(Coords))
!   Fsc = maxval(abs( AppForces(:,1:3) ));
!   Msc = maxval(abs( AppForces(:,4:6) ));
!
!   ! Correct Msc accounting for forces contribution. It is assumed that the
!   ! beam is constrained at Node 1 and that here the coordinates are (0,0,0)
!   do i =1,3
!       !print *, 'i=', 1, ' Msc=', Msc
!       do j = 1,3
!           if (i /= j) then
!               MSc = max(Msc,maxval(abs( AppForces(:,i)*Coords(:,j) )))
!           end if
!       end do
!   end do
!
!   ! avoid zero scaling factors
!   if (abs(Fsc) < 1e-8) then
!       Fsc=1.0_8
!   end if
!   if (abs(Msc) < 1e-8) then
!       Msc=1.0_8
!   end if
!
! ! Determine Absolute tolerances
! TaRes    = max(Fsc,Msc) *Options%MinDelta ! this could be made more severe picking the minimum.
! TaResFrc =     Fsc      *Options%MinDelta
! TaResMmt =         Msc  *Options%MinDelta
!
! TaX   = max(Possc,Psisc)*Options%MinDelta
! TaPos =     Possc       *Options%MinDelta
! TaPsi =           Psisc *Options%MinDelta

DX_old = 1.0d0*options%mindelta

! Initialize geometric constants and system state.
  ListIN = 0
  do k=1,size(Node)
    ListIN(k)=Node(k)%Vdof
  end do

  Kglobal = 0.0d0
  Mglobal = 0.0d0
  gravity_forces = 0.0d0
  nodal_gravity_forces = 0.d00
  Qglobal= 0.d0
  DeltaX = 0.d0
  Fglobal = 0.0d0

! Apply loads in substeps.
  do iLoadStep=1,Options%NumLoadSteps
    Iter  = 0
    Delta = Options%MinDelta+1.d0

! Iteration until convergence.
  converged=.false.
    do while (converged .eqv. .false.)!(Delta.gt.Options%MinDelta)
      Iter= Iter+1
      if (Iter.gt.Options%MaxIterations) then
          print*, 'Residual is: ', maxval(abs(DeltaX))
          STOP 'Static equations did not converge (17235)'
      end if

      if ((Options%PrintInfo) .AND. (Iter.eq.1)) then
          write (*,'(17X,A12,A12,A11,A12,A12,A12,A12,A12,A12,A13,A11)')      &
              & 'DeltaF ','DeltaX ',                                         &
              & 'Res','ResRel', 'ResFrc','ResRelFrc','ResMmt','ResRelMmt',   &
              & 'ErX','ErPos ','ErPsi'

          ! write (*,'(A16,$)') 'Tolerance'
          ! write (*,'(2X,1PE10.3,2X,1PE10.3,2X,$)') Options%MinDelta,Options%MinDelta
          ! write (*,'(2X,1PE10.3,2X,1PE10.3,2X,$)') TaRes,   Options%MinDelta
          ! write (*,'(2X,1PE10.3,2X,1PE10.3,2X,$)') TaResFrc,Options%MinDelta
          ! write (*,'(2X,1PE10.3,2X,1PE10.3,2X,$)') TaResMmt,Options%MinDelta
          ! write (*,'(2X,1PE10.3,2X,1PE10.3,2X,1PE10.3,2X)') TaX,TaPos,TaPsi

          write (*,'(A8,A8,A13,A11,A12,A12,A12,A12,A12,A12,A12,A13,A11)')      &
              & 'LoadStep','Subiter',                                          &
              & 'DeltaF ','DeltaX ',                                           &
              & 'Res','ResRel', 'ResFrc','ResRelFrc','ResMmt','ResRelMmt', &
              & 'ErX','ErPos ','ErPsi'

      end if
      ! if (Options%PrintInfo) write(*,'(I8,I8,$)')  iLoadStep, Iter

! Assembly matrices and functional.
      Qglobal=0.d0
      Kglobal=0.0d0
      Fglobal=0.0d0
      Mglobal = 0.0d0
      call cbeam3_asbly_static (numdof, n_elem, n_node,Elem,Node,Coords,Psi0,&
                                PosDefor,PsiDefor,&
                                AppForces*dble(iLoadStep)/dble(Options%NumLoadSteps), &
                                Kglobal,Fglobal,Qglobal,Options,Mglobal,MRR)

      Qglobal= Qglobal - dble(iLoadStep)/dble(Options%NumLoadSteps) * &
      &              MATMUL(Fglobal,fem_m2v(AppForces,NumDof,Filter=ListIN))

      if (options%gravity_on) then
          ! print*, 'Mglobal cbeam3'
          ! print*, Mglobal(1, 1:6)
          ! print*, Mglobal(2, 1:6)
          ! print*, Mglobal(3, 1:6)
          ! print*, Mglobal(4, 1:6)
          ! print*, Mglobal(5, 1:6)
          ! print*, Mglobal(6, 1:6)
          ! print*, '--'
          ! print*, 'gravity cbeam3'
          ! print*, cbeam3_asbly_gravity_static(6,options)
        nodal_gravity_forces = -fem_v2m(MATMUL(Mglobal,&
                                         cbeam3_asbly_gravity_static(NumDof + 6,&
                                                                     options)),&
                                  n_node, 6)
        gravity_forces = -MATMUL(MRR, cbeam3_asbly_gravity_static(6, options))
        ! print*, 'MRR'
        ! print*, MRR(1, :)
        ! print*, MRR(2, :)
        ! print*, MRR(3, :)
        ! print*, MRR(4, :)
        ! print*, MRR(5, :)
        ! print*, MRR(6, :)
        ! print*, '--'
        ! print*, 'product'
        ! print*, total_gravity_acceleration
        Qglobal = Qglobal - dble(iLoadStep)/dble(Options%NumLoadSteps)*&
                  fem_m2v(nodal_gravity_forces, numdof, filter=ListIN)
      end if

! Solve equation and update the global vectors.
      ! call lu_sparse(ks,Kglobal,-Qglobal,DeltaX)
      DeltaX = 0.0d0
      call lu_solve(size(Kglobal, dim=1), Kglobal,-Qglobal,DeltaX)
      DeltaX = 0.7d0*DeltaX
      call cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)

      if (iter > 1) then
          if (maxval(abs(DeltaX)) < DX_old) then
              converged = .TRUE.

            !   print*, 'APPFORCES = ', sum(appforces, dim=1)
            !   print*, 'gravity_forces = ', sum(gravity_forces, dim=1)
            !   print*, sum(gravity_forces + AppForces, dim=1)
          end if
      end if

    if (iter == 1) then
      DX_old = max(1.0d0, maxval(abs(DeltaX)))*options%mindelta
    end if
! Convergence parameter delta (original):
 !      call delta_check(Qglobal,DeltaX,Delta,passed_delta,Options%MinDelta,Options%PrintInfo)
 ! ! Check convergence using the residual:
 !      call separate_dofs(Qglobal,(/1,2,3/),(/4,5,6/),QglFrc,QglMmt)
 !
 !      call residual_check(Iter,Qglobal,Res,Res0,passed_res,Options%MinDelta,&
 !                         &TaRes,Options%PrintInfo  )
 !
 !      if ( (iLoadStep .eq. 1) .and. (Iter.eq.2) ) then
 !          ! update forcesd and moments residual at 2nd iteration to avoid zero
 !          ! due to trivial solution
 !          if (maxval(abs( AppForces(:,1:3))) < 1e-8) then
 !              ! jump first iteration
 !              call residual_check(1,QglFrc,ResFrc,ResFrc0,passed_resfrc,Options%MinDelta,&
 !                                 &TaResFrc, Options%PrintInfo)
 !          else
 !              call residual_check(Iter  ,QglFrc,ResFrc,ResFrc0,passed_resfrc,Options%MinDelta,&
 !                                 &TaResFrc, Options%PrintInfo)
 !          end if
 !
 !          if (maxval(abs( AppForces(:,4:6))) < 1e-8) then
 !              call residual_check(1,QglMmt,ResMmt,ResMmt0,passed_resmmt,Options%MinDelta,&
 !                                 &TaResMmt, Options%PrintInfo)
 !          else
 !              call residual_check(Iter,QglMmt,ResMmt,ResMmt0,passed_resmmt,Options%MinDelta,&
 !                                 &TaResMmt, Options%PrintInfo)
 !          end if
 !
 !      else
 !          call residual_check(Iter  ,QglFrc,ResFrc,ResFrc0,passed_resfrc,Options%MinDelta,&
 !                             &TaResFrc, Options%PrintInfo)
 !          call residual_check(Iter,QglMmt,ResMmt,ResMmt0,passed_resmmt,Options%MinDelta,&
 !                             &TaResMmt, Options%PrintInfo)
 !      end if
 !
 ! ! SuperLinear Convergence Test
 !      call separate_dofs(DeltaX,(/1,2,3/),(/4,5,6/),DeltaPos,DeltaPsi)
 !
 !      call error_check(Iter,DeltaX  ,DX_old  ,DX_now  ,ErrX  ,passed_err   ,TaX  , Options%PrintInfo)
 !      call error_check(Iter,DeltaPos,DPos_old,DPos_now,ErrPos,passed_errpos,TaPos, Options%PrintInfo)
 !      call error_check(Iter,DeltaPsi,DPsi_old,DPsi_now,ErrPsi,passed_errpsi,TaPsi, Options%PrintInfo)
 !
 !      DX_old   = DX_now
 !      DPos_old = DPos_now
 !      DPsi_old = DPsi_now
 !
 !      if (Options%PrintInfo) write (*,'(1X)')
 !
 ! ! Global Convergence
 !    if (passed_res .eqv. .true.) then
 !        converged=.true.
 !        if (Options%PrintInfo) write (*,'(A)') 'Global residual converged!'
 !    end if
 !
 !    if  (passed_err .eqv. .true. ) then
 !        converged=.true.
 !        if (Options%PrintInfo) write (*,'(A)') 'Global error converged!'
 !    end if
 !
 !
 !    if ( (passed_resfrc .eqv. .true.) .and. (passed_resmmt .eqv. .true.) ) then
 !        converged=.true.
 !        if (Options%PrintInfo) write (*,'(A)') 'Forces and Moments residual converged!'
 !    end if
 !
 !    if ( (passed_errpos .eqv. .true.) .and. (passed_errpsi .eqv. .true.) ) then
 !        converged=.true.
 !        if (Options%PrintInfo) write (*,'(A)') 'Displacements and Rotations error converged!'
 !    end if

    end do
  end do
 end subroutine cbeam3_solv_nlnstatic_old


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_LINSTATIC
!
!-> Description:
!
!    Linear steady-state solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!    a. NumDof is 6*Number of Independent Nodes (see xbeam_undef_dofs &
!    xbeam_undef_nodeindep)
!    b. ForceStatic is a matrix (row: global numbering; columns: forces and Moments)
!    c. Psi0 (PsiIni is main): CRV at the nodes [Psi0(NumElems,MaxElNod,3)].
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine cbeam3_solv_linstatic (NumDof,Elem,Node,AppForces,Coords,Psi0, &
! &                                  PosDefor,PsiDefor,Options)
!
!   use lib_fem
!   use lib_sparse
!   use lib_lu
!   use cbeam3_asbly
!
! ! I/O Variables.
!   integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
!   type(xbelem),intent(in)    :: Elem(:)           ! Element information.
!   type(xbnode),intent(in)    :: Node(:)           ! Nodal information.
!   real(8),      intent(in)   :: AppForces (:,:)   ! Applied nodal forces.
!   real(8),      intent(in)   :: Coords   (:,:)    ! Initial coordinates of the grid points.
!   real(8),      intent(in)   :: Psi0     (:,:,:)  ! Initial CRV of the nodes in the elements.
!   real(8),      intent(inout):: PosDefor (:,:)    ! Current coordinates of the grid points
!   real(8),      intent(inout):: PsiDefor (:,:,:)  ! Current CRV of the nodes in the elements.
!   type(xbopts),intent(in)    :: Options           ! Solver parameters.
!
! ! Local variables.
!   integer:: fs                             ! Current storage size of force matrix.
!   integer:: k                              ! Auxiliary integer.
!   integer:: ks                             ! Current storage size of stiffness matrix.
!   integer:: NumN                           ! Number of nodes in the model.
!
!   integer,allocatable::      ListIN (:)    ! List of independent nodes.
!   real(8),allocatable::      DeltaX (:)    ! Unknown in the linearized system.
!   real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
!   type(sparse),allocatable:: Fglobal(:)    ! Influence coefficients matrix for applied forces.
!   type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
!
!   print*, "In cbeam3_solv_linstatic"
!
! ! Initialize.
!   NumN=size(Node)
!   allocate (ListIN (NumN));
!   do k=1,NumN
!     ListIN(k)=Node(k)%Vdof
!   end do
!
! ! Allocate memory for solver (Use a conservative estimate of the size of the matrix Kglobal).
! ! - DimMat=24 is a parameter
! ! - ks and fs are defined by the sparse_zero call and initialised to zero.
! ! - Kglobal and Fglobal are empty (0 at (0,0) element)
! ! - They are both output
!   allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
!   allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
!   allocate (Qglobal(NumDof));        Qglobal= 0.d0
!   allocate (DeltaX (NumDof));        DeltaX = 0.d0
!
!  do k=1, size(AppForces(:,1))
!      print*, AppForces(k, :)
!  end do
!
! ! Assembly matrices and functional.
!   Qglobal=0.d0
!   call sparse_zero (ks,Kglobal) ! sm BUG: isn't this repeated?
!   call sparse_zero (fs,Fglobal)
!
!   call cbeam3_asbly_static (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,AppForces, & ! input
! &                           Kglobal,Fglobal,Qglobal,Options)               ! output (except for Options)
!
!   ! Check Kglobal is filled correctly
!   ! do k=1,size(Kglobal) ! ADC: CHANGED (was size(Kglobal) + 1)
!   !   if ((Kglobal(k)%i>NumDof) .or. (Kglobal(k)%j>NumDof)) then
!   !     print *, 'Out of Bounds!!! Allocated: (', Kglobal(k)%i,',',Kglobal(k)%j,')'
!   !     stop 'Execution terminated!'
!   !   end if
!   ! end do
!
! ! Forces on the unconstrained nodes.
! ! sm: AppForces has shape (Nodes,6), where the columns contain forces and moments.
! ! fem_m2v reorders them into a vector (i.e. for node ii, (ii-1)+1 will be the x
! ! force and (ii-1)+6 the z moment. Nodes for which the solution has not to be
! ! found will not be counted.
!   Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(AppForces,NumDof,Filter=ListIN))
!
! !stop 'sparse_matvmul ok!'
!
! ! Solve equation and update the global vectors.
! ! Kglobal * deltaX = Qglobal
! !#ifdef NOLAPACK
! !  call lu_sparse(ks,Kglobal,Qglobal,DeltaX)
! !#else
! !  call lapack_sparse (ks,Kglobal,Qglobal,DeltaX)
! !#endif
!   call lu_sparse(ks,Kglobal,Qglobal,DeltaX)
!
!   ! sm: PosDefor and PsiDefor are the only one being updated
!   call cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)
!   print*, "RUN"
!
!   deallocate (Kglobal,Qglobal,DeltaX)
!   print*, Coords(size(Coords(:,1)),:)
!   print*, PosDefor(size(Coords(:,1)),:)
!   return
!  end subroutine cbeam3_solv_linstatic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_MODAL_UPDATED
!   Find linear vibration modes.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_modal_updated  (num_dof,&
                                        n_elem,&
                                        n_node,&
                                        Elem,&
                                        Node,&
                                        Vrel,&
                                        Coords,&
                                        Psi0,&
                                        PosDefor,&
                                        PsiDefor,&
                                        Options,&
                                        FullMglobal,&
                                        FullCglobal,&
                                        FullKglobal&
                                        )

  ! Libraries
  use cbeam3_asbly

  ! I/O Variables.
  integer,      intent(in) :: num_dof               ! Degrees of freedom
  integer,      intent(in) :: n_elem                 ! Number of elements
  integer,      intent(in) :: n_node                  ! Number of nodes
  type(xbelem), intent(in) :: Elem      (n_elem)     ! Element information.
  type(xbnode), intent(in) :: Node      (n_node)     ! Nodal information.
  real(8),      intent(in) :: Vrel      (6)       ! Velocity of the reference frame.
  real(8),      intent(in) :: Coords    (n_node,3)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (n_elem,3,3) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (n_node,3)    ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (n_elem,3,3) ! Current CRV of the nodes in the elements.
  type(xbopts), intent(in) :: Options           ! Solver parameters.

  ! System matrices
  real(8), intent(inout) :: FullCglobal(num_dof,num_dof)
  real(8), intent(inout) :: FullMglobal(num_dof,num_dof)
  real(8), intent(inout) :: FullKglobal(num_dof,num_dof)

  ! Local variables.
  integer:: k                            ! Counters.

  ! Initialize matrices
  FullCglobal = 0.d0
  FullMglobal = 0.d0
  FullKglobal = 0.d0

  ! Assembly the martices
  call cbeam3_asbly_modal_updated (num_dof,n_elem,n_node,Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Vrel, &
        &                          FullMglobal,FullCglobal,FullKglobal,Options)

  return
 end subroutine cbeam3_solv_modal_updated


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_MODAL_OLD
!
!-> Description:
!
!   Find linear vibration modes.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_modal_old  (iOut,NumDof,Elem,Node,Vrel,Coords,Psi0, &
&                               PosDefor,PsiDefor,Options)
  use lib_fem
  use lib_sparse
  use cbeam3_asbly

! I/O Variables.
  integer,      intent(in) :: iOut              ! Output file.
  integer,      intent(in) :: NumDof            ! Number of independent DoFs.
  type(xbelem), intent(in) :: Elem      (:)     ! Element information.
  type(xbnode), intent(in) :: Node      (:)     ! Nodal information.
  real(8),      intent(in) :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in) :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in) :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(in) :: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(in) :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  type(xbopts), intent(in) :: Options           ! Solver parameters.

! Local variables.
  integer:: k                            ! Counters.
  integer:: cs,ks,ms
  type(sparse),allocatable:: Cglobal(:)     ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)     ! Global mass matrix in sparse storage.


! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)

! Compute tangent matrices at initial time.
  call cbeam3_asbly_modal_old (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,Vrel(1,:), &
&                          ms,Mglobal,cs,Cglobal,ks,Kglobal,Options)

! Write matrices in files.
  ! open (unit=71,file='Msparse',status='replace')
  ! do k=1,ms
  !   write (71,'(2I4,1PE18.10)') Mglobal(k)%i,Mglobal(k)%j,Mglobal(k)%a
  ! end do
  ! close (71)
  !
  ! open (unit=72,file='Csparse',status='replace')
  ! do k=1,cs
  !   write (72,'(2I4,1PE18.10)') Cglobal(k)%i,Cglobal(k)%j,Cglobal(k)%a
  ! end do
  ! close (72)
  !
  ! open (unit=73,file='Ksparse',status='replace')
  ! do k=1,ks
  !   write (73,'(2I4,1PE18.10)') Kglobal(k)%i,Kglobal(k)%j,Kglobal(k)%a
  ! end do
  ! close (73)

! End of routine.
  deallocate (Mglobal,Kglobal,Cglobal)
  return
 end subroutine cbeam3_solv_modal_old


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN_UPDATED
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem under applied forces.
!
!-> Remarks.-
!    ADC: this is a modified version, which outputs the whole history
!         of displacements and rotations for every timestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn_updated (iOut,&
                                        NumDof,&
                                        Time,&
                                        Elem,&
                                        Node,&
                                        F0,&
                                        Fa,&
                                        Ftime,&
                                        Vrel,&
                                        VrelDot,&
                                        Coords,&
                                        Psi0,&
                                        PosDefor,&
                                        PsiDefor,&
                                        PosDotDefor,&
                                        PsiDotDefor,&
                                        Options,&
                                        pos_def_history,&
                                        psi_def_history,&
                                        pos_dot_def_history,&
                                        psi_dot_def_history,&
                                        in_quat)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
  use lib_lu
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)  ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)    ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),intent(IN), optional   :: in_quat(:)
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

  real(8), intent(OUT)       :: pos_def_history(:,:,:)
  real(8), intent(OUT)       :: pos_dot_def_history(:,:,:)
  real(8), intent(OUT)       :: psi_def_history(:,:,:,:)
  real(8), intent(OUT)       :: psi_dot_def_history(:,:,:,:)

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
  real(8),allocatable:: X(:), DX(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for system matrices.
  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
  real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  type(outopts)      ::  OutOptions    ! Output options.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Cao(3,3)
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

! Initialize.
  NumN=size(Node)
  NumE(1)=size(Elem)
  allocate (ListIN (NumN));
  do k=1,NumN
    ListIN(k)=Node(k)%Vdof
  end do
  gamma=0.5d0+Options%NewmarkDamp
  beta =0.25d0*(gamma+0.5d0)*(gamma+0.5d0)

! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
  call sparse_allocate (Asys   , NumDof, NumDof)
  call sparse_allocate (Mglobal, NumDof, NumDof)
  call sparse_allocate (Cglobal, NumDof, NumDof)
  call sparse_allocate (Kglobal, NumDof, NumDof)
  call sparse_allocate (Fglobal, NumDof, NumDof)
  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  if (.NOT. present(in_quat)) then
      Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4
  else
      Quat = in_quat
  end if
  Cao = xbeam_Rot(Quat)

! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  if (options%gravity_on .eqv. .TRUE.) then
      call cbeam3_asbly_dynamic (Elem,&
                                 Node,&
                                 Coords,&
                                 Psi0,&
                                 PosDefor,&
                                 PsiDefor,&
                                 PosDotDefor,&
                                 PsiDotDefor,&
                                 0.d0*PosDefor,&
                                 0.d0*PsiDefor,&
                                 F0+Ftime(1)*Fa,Vrel(1,:),&
                                 VrelDot(1,:),&
                                 ms,Mglobal,&
                                 Mvel,&
                                 cs,Cglobal,&
                                 Cvel,&
                                 ks,Kglobal,&
                                 fs,Fglobal,&
                                 Qglobal,&
                                 Options,&
                                 Cao)

      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
                                        NumDof,&
                                        cbeam3_asbly_gravity_dynamic(NumDof, options, Cao), '749')
  else
      call cbeam3_asbly_dynamic (Elem,&
                                 Node,&
                                 Coords,&
                                 Psi0,&
                                 PosDefor,&
                                 PsiDefor,&
                                 PosDotDefor,&
                                 PsiDotDefor,&
                                 0.d0*PosDefor,&
                                 0.d0*PsiDefor,&
                                 F0+Ftime(1)*Fa,&
                                 Vrel(1,:),&
                                 VrelDot(1,:),&
                                 ms,Mglobal,&
                                 Mvel,&
                                 cs,Cglobal,&
                                 Cvel,&
                                 ks,Kglobal,&
                                 fs,Fglobal,&
                                 Qglobal,&
                                 Options,&
                                 Cao)
  end if

  Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))

  call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
  call lu_sparse(as,Asys,-Qglobal,dXddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    if (Options%PrintInfo) then
      call out_time(iStep,Time(iStep+1),Text)
      ! write (*,'(5X,A,$)') trim(Text)
    end if
! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
    dXdt = dXdt + (1.d0-gamma)*dt*dXddt
    dXddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) then
        write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        STOP 'Solution did not converge (18235)'
      end if

! Update nodal positions and velocities .
      call cbeam3_solv_state2disp (Elem,&
                                   Node,&
                                   Coords,&
                                   Psi0,&
                                   X,&
                                   dXdt,&
                                   PosDefor,&
                                   PsiDefor,&
                                   PosDotDefor,&
                                   PsiDotDefor)

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qglobal= 0.d0
      Mvel   = 0.d0
      Cvel   = 0.d0
      call sparse_zero (ms,Mglobal)
      call sparse_zero (cs,Cglobal)
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)

  if (options%gravity_on .eqv. .TRUE.) then
      call cbeam3_asbly_dynamic (Elem,&
                                 Node,&
                                 Coords,&
                                 Psi0,&
                                 PosDefor,&
                                 PsiDefor,&
                                 PosDotDefor,&
                                 PsiDotDefor,&
                                 0.d0*PosDefor,&
                                 0.d0*PsiDefor,&
                                 F0+Ftime(iStep+1)*Fa,&
                                 Vrel(iStep+1,:),&
                                 VrelDot(iStep+1,:),&
                                 ms,Mglobal,&
                                 Mvel,&
                                 cs,Cglobal,&
                                 Cvel,&
                                 ks,Kglobal,&
                                 fs,Fglobal,&
                                 Qglobal,&
                                 Options,&
                                 Cao)
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
                                        NumDof,&
                                        cbeam3_asbly_gravity_dynamic(NumDof,&
                                                                     options,&
                                                                     Cao), &
                                        '807')

  else
      call cbeam3_asbly_dynamic (Elem,&
                                 Node,&
                                 Coords,&
                                 Psi0,&
                                 PosDefor,&
                                 PsiDefor,&
                                 PosDotDefor,&
                                 PsiDotDefor,&
                                 0.d0*PosDefor,&
                                 0.d0*PsiDefor,&
                                 F0+Ftime(iStep+1)*Fa,&
                                 Vrel(iStep+1,:),&
                                 VrelDot(iStep+1,:),&
                                 ms,Mglobal,&
                                 Mvel,&
                                 cs,Cglobal,&
                                 Cvel,&
                                 ks,Kglobal,&
                                 fs,Fglobal,&
                                 Qglobal,&
                                 Options,&
                                 Cao)
  end if

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))

! Compute the residual.
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
      Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))


! Check convergence.
      if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
        if (Options%PrintInfo) then
          write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        end if
        exit
      end if

! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))


! Calculation of the correction.
      call lu_sparse(as,Asys,-Qglobal,DX)

      X    = X     + DX
      dXdt = dXdt  + gamma/(beta*dt)*DX
      dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
    end do

! Update nodal positions and velocities on the current converged time step.
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

    pos_def_history(iStep + 1, :, :) = PosDefor
    pos_dot_def_history(iStep + 1, :, :) = PosDotDefor
    psi_def_history(iStep + 1, :, :, :) = PsiDefor
    psi_dot_def_history(iStep + 1, :, :, :) = PsiDotDefor

  end do
  deallocate (ListIN,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (X,DX,dXdt,dXddt)
  deallocate (Displ,Veloc)
 end subroutine cbeam3_solv_nlndyn_updated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN_STEP
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem under applied forces.
!
!-> Remarks.-
!    ADC: this is a modified version, which outputs the whole history
!         of displacements and rotations for every timestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn_step(num_dof,&
                                    num_elem,&
                                    num_node,&
                                    dt,&
                                    elem,&
                                    node,&
                                    static_forces,&
                                    dynamic_forces,&
                                    gravity_forces,&
                                    quat,&
                                    for_vel,&
                                    for_acc,&
                                    pos_ini,&
                                    psi_ini,&
                                    pos_def,&
                                    psi_def,&
                                    pos_dot_def,&
                                    psi_dot_def,&
                                    q,&
                                    dqdt,&
                                    options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
  use lib_lu
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)          :: num_dof! Number of independent DoFs.
  integer,      intent(IN)          :: num_elem
  integer,      intent(IN)          :: num_node
  real(8),      intent(IN)          :: dt
  type(xbelem), intent(in)          :: elem      (num_elem)     ! Element information.
  type(xbnode), intent(in)          :: node      (num_node)     ! Nodal information.
  real(8),      intent(IN)          :: static_forces(num_node, 6)
  real(8),      intent(IN)          :: dynamic_forces(num_node, 6)
  real(8),      intent(OUT)         :: gravity_forces(num_node, 6)
  real(8),      intent(INOUT)       :: quat(4)
  real(8),      intent(in)          :: for_vel(6)
  real(8),      intent(in)          :: for_acc(6)
  real(8),      intent(in)          :: pos_ini(num_node, 3)   ! Initial coordinates of the grid points.
  real(8),      intent(in)          :: psi_ini(num_elem, 3, 3) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout)       :: pos_def(num_node, 3)   ! Current coordinates of the grid points
  real(8),      intent(inout)       :: psi_def(num_elem, 3, 3) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout)       :: pos_dot_def(num_node, 3)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout)       :: psi_dot_def(num_elem, 3, 3)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)         :: q(num_dof)
  real(8),      intent(out)         :: dqdt(num_dof)
  type(xbopts) ,intent(in)          :: options           ! solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  integer:: k                              ! Counters.
  integer:: Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8)                           :: dXdt(num_dof),dXddt(num_dof)  ! Generalized coordinates and derivatives.
  real(8)                           :: X(num_dof), DX(num_dof)
  real(8)                           :: residual
  real(8)                           :: previous_residual

  integer                           ::  ListIN     (num_node)    ! List of independent nodes.

  ! Define variables for system matrices.
  real(8)                           :: Asys(num_dof, num_dof)
  real(8)                           :: Mglobal(num_dof, num_dof)
  real(8)                           :: Cglobal(num_dof, num_dof)
  real(8)                           :: Kglobal(num_dof, num_dof)
  real(8)                           :: Fglobal(num_dof, num_dof)
  real(8)                           :: Qglobal(num_dof)
  real(8)                           :: Mvel(num_dof, 6)
  real(8)                           :: Cvel(num_dof, 6)
  real(8)                           :: MSS_gravity(num_dof+6, num_dof+6)
  real(8)                           :: MRS_gravity(6, num_dof+6)
  real(8)                           :: MRR_gravity(6, 6)

  ! Rotation operator for body-fixed frame using quaternions
  real(8)                           :: Cao(3,3)
  real(8)                           :: Temp(4, 4)
  logical                           :: converged
  real(8)                           :: old_x
  real(8)                           :: old_dx

    ListIN = 0
    do k=1,num_node
        ListIN(k)=Node(k)%Vdof
    end do

    gamma = 0.5d0 + options%NewmarkDamp
    beta = 0.25d0*(gamma + 0.5d0)*(gamma + 0.5d0)

    ! update state vector
    X = 0.0d0
    dXdt = 0.0d0
    dXddt = 0.0d0
    call cbeam3_solv_disp2state(node,pos_def,psi_def,pos_dot_def,psi_dot_def,X,dXdt)

    X = X + dt*dXdt + (0.5d0 - beta)*dt*dt*dXddt
    dXdt = dXdt + (1.0d0 - gamma)*dt*dXddt
    dXddt = 0.0d0

    ! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(for_vel(4:6))*dt),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(for_vel(4:6))*dt),Quat))
    Cao = xbeam_Rot(Quat)

    mindelta = 0
    old_x = 1.0d0
    converged = .FALSE.
    old_DX = 1.0d0
    ! Iteration loop -----------------------------------------
    do iter = 1, options%maxiterations + 1
        if (iter == options%maxiterations + 1) then
            print*, 'Solver did not converge in ', iter, ' iterations.'
            exit
        end if

        ! update positions and velocities
        call cbeam3_solv_state2disp(Elem,&
                                    Node,&
                                    pos_ini,&
                                    psi_ini,&
                                    X,&
                                    dXdt,&
                                    pos_def,&
                                    psi_def,&
                                    pos_dot_def,&
                                    psi_dot_def)
        ! call cbeam3_solv_state2accel(Elem, Node, dqddt(1:numdof), pos_ddot_def, psi_ddot_def)

        ! system functionals and matrices initialisation
        Mvel = 0.0d0
        Cvel = 0.0d0
        Asys = 0.0d0
        Mglobal = 0.0d0
        Cglobal = 0.0d0
        Kglobal = 0.0d0
        Qglobal = 0.0d0
        Fglobal = 0.0d0

        ! system matrix generation
        call cbeam3_asbly_dynamic(num_dof,&
                                  num_node,&
                                  num_elem,&
                                  elem,&
                                  node,&
                                  pos_ini,&
                                  psi_ini,&
                                  pos_def,&
                                  psi_def,&
                                  pos_dot_def,&
                                  psi_dot_def,&
                                  ! pos_ddot_def,&
                                  ! psi_ddot_def,&
                                  0.0d0*pos_dot_def,&
                                  0.0d0*psi_dot_def,&
                                  static_forces + dynamic_forces,&
                                  for_vel,&
                                  for_acc,&
                                  Mglobal,&
                                  Mvel,&
                                  Cglobal,&
                                  Cvel,&
                                  Kglobal,&
                                  Fglobal,&
                                  Qglobal,&
                                  options,&
                                  Cao)

        ! compute residual
        Qglobal = Qglobal + MATMUL(Mglobal, dXddt)
        Qglobal = Qglobal + MATMUL(Mvel, for_acc)
        ! Qglobal = Qglobal - MATMUL(Fglobal, fem_m2v(static_forces + dynamic_forces, num_dof, filter=ListIN))
        Qglobal = Qglobal - MATMUL(Fglobal, fem_m2v(static_forces + dynamic_forces, num_dof, filter=ListIN))

        if (options%gravity_on) then
            call xbeam_asbly_M_gravity(num_dof,&
                                        num_node,&
                                        num_elem,&
                                        Elem,&
                                        Node,&
                                        pos_ini,&
                                        psi_ini,&
                                        pos_def,&
                                        psi_def,&
                                        MRS_gravity,&
                                        MSS_gravity,&
                                        MRR_gravity,&
                                        options)
            gravity_forces = -fem_v2m(MATMUL(MSS_gravity,&
                                      cbeam3_asbly_gravity_dynamic(Num_Dof + 6,options, Cao)),&
                                      num_node, 6)
            Qglobal = Qglobal - fem_m2v(gravity_forces, num_dof, filter=ListIN)
        end if

        ! Qtotal(1:numdof) = Qelast
        ! Qtotal(numdof+1:numdof+6) = Qrigid
        ! Qtotal(numdof+7:numdof+10) = matmul(CQQ, dQdt(numdof+7:numdof+10))

        ! call mat_addmat(0, 0, MSS, Mtotal)
        ! call mat_addmat(0, numdof, MSR, Mtotal)
        ! call mat_addmat(numdof, 0, MRS, Mtotal)
        ! call mat_addmat(numdof, numdof, MRR, Mtotal)
        ! call mat_addmat(numdof + 6, numdof + 6, unit4, Mtotal)

        ! Qtotal = Qtotal + matmul(Mtotal, dqddt)

        ! convergence check
        ! print*, 'maxval(abs(Qelast)) = ', maxval(abs(Qelast))
        ! if (maxval(abs(Qelast)) < mindelta .AND.&
        !     maxval(abs(Qrigid)) < MinDeltarigid) then
        !     ! print*, 'converged in ', iter
        !     converged = .TRUE.
        ! end if

        ! ! damping and stiffness matrices
        ! call mat_addmat(0, 0, CSS, Ctotal)
        ! call mat_addmat(0, numdof, CSR, Ctotal)
        ! call mat_addmat(numdof, 0, CRS, Ctotal)
        ! call mat_addmat(numdof, numdof, CRR, Ctotal)
        ! call mat_addmat(numdof+6, numdof, CQR, Ctotal)
        ! call mat_addmat(numdof+6, numdof+6, CQQ, Ctotal)
        !
        ! call mat_addmat(0, 0, KSS, Ktotal)
        ! call mat_addmat(numdof, 0, KRS, Ktotal)

        ! assembly of A matrix
        ! if (any(abs(Ktotal) > 1e11)) then
        !     print*, 'KTOTAL---------------------------------'
            ! call print_matrix('Mtotal', Mtotal)
        !     call print_matrix('pos_def', pos_def)
        !     call print_elem('failed_elements', elem)
        !     stop
        ! end if

        ! print*, 'Test'
        ! if (any(isnan(Ktotal))) then
        !     print*, 'Ktotal'
        !     stop
        ! end if
        ! if (any(isnan(Ctotal))) then
        !     print*, 'Ctotal'
        !     stop
        ! end if
        ! if (any(isnan(Mtotal))) then
        !     print*, 'Mtotal'
        !     stop
        ! end if
!
        ! call print_matrix('Kglobal', Kglobal)
        ! call print_matrix('Cglobal', Cglobal)
        ! call print_matrix('Mglobal', Mglobal)
        ! stop
        Asys = Kglobal
        Asys = Asys + Cglobal*gamma/(beta*dt)
        Asys = Asys + Mglobal/(beta*dt*dt)

        ! print*, size(Asys, dim=1), size(Asys, dim=2)
        ! if (any(isnan(Asys))) then
        !     print*, 'Asys'
        !     stop
        ! end if


        ! calculation of the correction
        DX = 0.0d0
        call lu_solve(num_dof, Asys, -Qglobal, DX)

        if (Iter > 1) then
            ! print*, (maxval(abs(DX))/old_DX)
            ! if ((sqrt(dot_product(q, q))-old_q)/old_q < options%MinDelta) then
            ! if (abs(maxval(abs(DQ)) - old_DQ)/old_DQ < options%MinDelta) then
            if (maxval(abs(DX))/old_DX < options%MinDelta) then
                converged = .TRUE.
            end if
        end if


        if (converged) then
            exit
        endif

        ! reconstruction of state vectors
        X = X + DX
        dXdt  = dXdt  + gamma/(beta*dt)*DX
        dXddt = dXddt + 1.d0/(beta*dt*dt)*DX

        if (iter == 1) then
            old_DX = max(maxval(abs(DX)), 1.0d0)
        end if
        ! old_q = sqrt(dot_product(q, q))
        ! end of convergence check

        ! call cbeam3_solv_state2disp(elem,&
        !                             node,&
        !                             pos_ini,&
        !                             psi_ini,&
        !                             Q(1:numdof),&
        !                             dQdt(1:numdof),&
        !                             pos_def,&
        !                             psi_def,&
        !                             pos_dot_def,&
        !                             psi_dot_def)
    end do

    ! dQdt(numdof+7:numdof+10) = rot_unit(dQdt(numdof+7:numdof+10))
    ! update of the nodal positions
    call cbeam3_solv_state2disp(elem,&
                                node,&
                                pos_ini,&
                                psi_ini,&
                                X,&
                                dXdt,&
                                pos_def,&
                                psi_def,&
                                pos_dot_def,&
                                psi_dot_def)


  q(1:size(q) - 10) = X
  dqdt(1:size(q) - 10) = dXdt

    ! quat = dQdt(numdof+7:numdof+10)
    ! if (options%OutInaframe) then
    !     for_vel = dQdt(numdof+1:numdof+6)
    !     for_acc = dQddt(numdof+1:numdof+6)
    ! else
    !     print*, 'outinaframe is not TRUE, please check!'
    !     stop
    ! end if

! ! Initialize.
!   do k=1, num_node
!     ListIN(k)=Node(k)%Vdof
!   end do
!   gamma=0.5d0+Options%NewmarkDamp
!   beta =0.25d0*(gamma+0.5d0)*(gamma+0.5d0)
!
!   Cao = xbeam_Rot(Quat)
!
! ! Extract initial displacements and velocities.
!   call cbeam3_solv_disp2state (node,pos_def,psi_def,pos_dot_def,psi_dot_def,X,dXdt)
!
! ! Compute initial acceleration (we are neglecting qdotdot in Kmass).
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(1)*Fa,&
!                                  Vrel(1,:),&
!                                  VrelDot(1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof, options, Cao), '749')
!   else
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(1)*Fa,&
!                                  Vrel(1,:),&
!                                  VrelDot(1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!   end if
!
!   Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))
!
!   call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
!   call lu_sparse(as,Asys,-Qglobal,dXddt)
!
! ! Loop in the time steps.
!   do iStep=1,size(Time)-1
!     dt= Time(iStep+1)-Time(iStep)
!     if (Options%PrintInfo) then
!       call out_time(iStep,Time(iStep+1),Text)
!       ! write (*,'(5X,A,$)') trim(Text)
!     end if
! ! Update transformation matrix for given angular velocity
!     call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
!     Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
!     Cao = xbeam_Rot(Quat)
!
! ! Predictor step.
!     X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
!     dXdt = dXdt + (1.d0-gamma)*dt*dXddt
!     dXddt= 0.d0
!
! ! Iteration until convergence.
!
!     previous_residual = 1e10
!     do Iter=1,Options%MaxIterations+1
!       if (Iter.gt.Options%MaxIterations) then
!         ! write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         STOP 'Solution did not converge (18235)'
!       end if
!
! ! Update nodal positions and velocities .
!       call cbeam3_solv_state2disp (Elem,&
!                                    Node,&
!                                    Coords,&
!                                    Psi0,&
!                                    X,&
!                                    dXdt,&
!                                    PosDefor,&
!                                    PsiDefor,&
!                                    PosDotDefor,&
!                                    PsiDotDefor)
!
! ! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
!       Qglobal= 0.d0
!       Mvel   = 0.d0
!       Cvel   = 0.d0
!       call sparse_zero (ms,Mglobal)
!       call sparse_zero (cs,Cglobal)
!       call sparse_zero (ks,Kglobal)
!       call sparse_zero (fs,Fglobal)
!
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(iStep+1)*Fa,&
!                                  Vrel(iStep+1,:),&
!                                  VrelDot(iStep+1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof,&
!                                                                      options,&
!                                                                      Cao), &
!                                         '807')
!
!   else
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(iStep+1)*Fa,&
!                                  Vrel(iStep+1,:),&
!                                  VrelDot(iStep+1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!   end if
!
! ! Compute admissible error.
!       MinDelta=Options%MinDelta
!
! ! Compute the residual.
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
!       Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
!
!
! ! Check convergence.
!     !   if (maxval(abs(DX)).lt.MinDelta) then
!     residual = maxval(abs(DX))
!     ! if (Iter > 1) then
!       if (abs((residual - previous_residual)/previous_residual).lt.MinDelta) then
!         if (Options%PrintInfo) then
!           write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         end if
!         exit
!       end if
!     end if
!     previous_residual = residual
!
! ! Compute Jacobian
!       call sparse_zero (as,Asys)
!       call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
!       call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
!       call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))
!
! ! Calculation of the correction.
!       call lu_sparse(as,Asys,-Qglobal,DX)
!
!       X    = X     + DX
!       dXdt = dXdt  + gamma/(beta*dt)*DX
!       dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
!     end do
!
! ! Update nodal positions and velocities on the current converged time step.
!     call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!   end do
!
!   deallocate (ListIN,Mvel,Cvel)
!   deallocate (Asys,Fglobal,Mglobal)
!   deallocate (Kglobal,Cglobal,Qglobal)
!   deallocate (X,DX,dXdt,dXddt)
!   deallocate (Displ,Veloc)
 end subroutine cbeam3_solv_nlndyn_step
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !-> Subroutine CBEAM3_SOLV_NLNDYN_STEP
! !
! !-> Description:
! !
! !    Linear dynamic solution of multibeam problem under applied forces.
! !
! !-> Remarks.-
! !    ADC: this is a modified version, which outputs the whole history
! !         of displacements and rotations for every timestep
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine cbeam3_solv_nlndyn_step(iOut,&
!                                     NumDof,&
!                                     Time,&
!                                     Elem,&
!                                     Node,&
!                                     F0,&
!                                     Fa,&
!                                     Ftime,&
!                                     Vrel,&
!                                     VrelDot,&
!                                     Coords,&
!                                     Psi0,&
!                                     PosDefor,&
!                                     PsiDefor,&
!                                     PosDotDefor,&
!                                     PsiDotDefor,&
!                                     Options,&
!                                     in_quat)
!   use lib_fem
!   use lib_rot
!   use lib_rotvect
!   use lib_lu
!   use lib_out
!   use lib_sparse
!   use lib_lu
!   use cbeam3_asbly
!   use lib_xbeam
!
! ! I/O Variables.
!   integer,      intent(in)       :: iOut              ! Output file.
!   integer,      intent(in)       :: NumDof            ! Number of independent DoFs.
!   real(8),      intent(in)       :: Time      (:)     ! Time steps.
!   type(xbelem), intent(in)       :: Elem      (:)     ! Element information.
!   type(xbnode), intent(in)       :: Node      (:)     ! Nodal information.
!   real(8),      intent(in)       :: F0        (:,:)   ! Applied static nodal forces.
!   real(8),      intent(in)       :: Fa        (:,:)  ! Amplitude of the dynamic nodal forces.
!   real(8),      intent(in)       :: Ftime     (:)    ! Time history of the applied forces.
!   real(8),      intent(in)       :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
!   real(8),      intent(in)       :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
!   real(8),      intent(in)       :: Coords    (:,:)   ! Initial coordinates of the grid points.
!   real(8),      intent(in)       :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
!   real(8),      intent(inout)    :: PosDefor  (:,:)   ! Current coordinates of the grid points
!   real(8),      intent(inout)    :: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
!   real(8),      intent(inout)    :: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
!   real(8),      intent(inout)    :: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
!   type(xbopts) ,intent(in)       :: Options           ! Solver parameters.
!   real(8),intent(IN), optional   :: in_quat(:)
!
! ! Local variables.
!   real(8):: beta,gamma                     ! Newmark coefficients.
!   real(8):: dt                             ! Time step
!   integer:: k                              ! Counters.
!   integer:: iStep,Iter                     ! Counters on time steps and subiterations.
!   real(8):: MinDelta                       ! Value of Delta for convergence.
!   integer:: NumE(1)                        ! Number of elements in the model.
!   integer:: NumN                           ! Number of nodes in the model.
!
!   real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
!   real(8),allocatable:: X(:), DX(:)
!   real(8)               :: residual
!   real(8)               :: previous_residual
!
!   integer,allocatable::  ListIN     (:)    ! List of independent nodes.
!
!   ! Define variables for system matrices.
!   integer:: as,cs,ks,ms,fs
!   type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
!   type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
!   type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
!   type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
!   type(sparse),allocatable:: Fglobal(:)
!   real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
!   real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
!   real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.
!
!   ! Define vaiables for output information.
!   character(len=80)  ::  Text          ! Text with current time information to print out.
!   type(outopts)      ::  OutOptions    ! Output options.
!   real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
!   real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.
!
!   ! Rotation operator for body-fixed frame using quaternions
!   real(8) :: Cao(3,3)
!   real(8) :: Quat(4)
!   real(8) :: Temp(4,4)
!
! ! Initialize.
!   NumN=size(Node)
!   NumE(1)=size(Elem)
!   allocate (ListIN (NumN));
!   do k=1,NumN
!     ListIN(k)=Node(k)%Vdof
!   end do
!   gamma=0.5d0+Options%NewmarkDamp
!   beta =0.25d0*(gamma+0.5d0)*(gamma+0.5d0)
!
! ! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
!   call sparse_allocate (Asys   , NumDof, NumDof)
!   call sparse_allocate (Mglobal, NumDof, NumDof)
!   call sparse_allocate (Cglobal, NumDof, NumDof)
!   call sparse_allocate (Kglobal, NumDof, NumDof)
!   call sparse_allocate (Fglobal, NumDof, NumDof)
!   allocate (Qglobal(NumDof));   Qglobal= 0.d0
!   allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
!   allocate (Cvel   (NumDof,6)); Cvel   = 0.d0
!
!   allocate (X     (NumDof)); X      = 0.d0
!   allocate (DX    (NumDof)); DX     = 0.d0
!   allocate (dXdt  (NumDof)); dXdt   = 0.d0
!   allocate (dXddt (NumDof)); dXddt  = 0.d0
!
! ! Compute system information at initial condition.
!   allocate (Veloc(NumN,6)); Veloc=0.d0
!   allocate (Displ(NumN,6)); Displ=0.d0
!
! ! Allocate quaternions and rotation operator for initially undeformed system
!   if (.NOT. present(in_quat)) then
!       Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4
!   else
!       Quat = in_quat
!   end if
!   Cao = xbeam_Rot(Quat)
!
! ! Extract initial displacements and velocities.
!   call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)
!
! ! Compute initial acceleration (we are neglecting qdotdot in Kmass).
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(1)*Fa,&
!                                  Vrel(1,:),&
!                                  VrelDot(1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof, options, Cao), '749')
!   else
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(1)*Fa,&
!                                  Vrel(1,:),&
!                                  VrelDot(1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!   end if
!
!   Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))
!
!   call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
!   call lu_sparse(as,Asys,-Qglobal,dXddt)
!
! ! Loop in the time steps.
!   do iStep=1,size(Time)-1
!     dt= Time(iStep+1)-Time(iStep)
!     if (Options%PrintInfo) then
!       call out_time(iStep,Time(iStep+1),Text)
!       ! write (*,'(5X,A,$)') trim(Text)
!     end if
! ! Update transformation matrix for given angular velocity
!     call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
!     Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
!     Cao = xbeam_Rot(Quat)
!
! ! Predictor step.
!     X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
!     dXdt = dXdt + (1.d0-gamma)*dt*dXddt
!     dXddt= 0.d0
!
! ! Iteration until convergence.
!
!     previous_residual = 1e10
!     do Iter=1,Options%MaxIterations+1
!       if (Iter.gt.Options%MaxIterations) then
!         ! write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         STOP 'Solution did not converge (18235)'
!       end if
!
! ! Update nodal positions and velocities .
!       call cbeam3_solv_state2disp (Elem,&
!                                    Node,&
!                                    Coords,&
!                                    Psi0,&
!                                    X,&
!                                    dXdt,&
!                                    PosDefor,&
!                                    PsiDefor,&
!                                    PosDotDefor,&
!                                    PsiDotDefor)
!
! ! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
!       Qglobal= 0.d0
!       Mvel   = 0.d0
!       Cvel   = 0.d0
!       call sparse_zero (ms,Mglobal)
!       call sparse_zero (cs,Cglobal)
!       call sparse_zero (ks,Kglobal)
!       call sparse_zero (fs,Fglobal)
!
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(iStep+1)*Fa,&
!                                  Vrel(iStep+1,:),&
!                                  VrelDot(iStep+1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof,&
!                                                                      options,&
!                                                                      Cao), &
!                                         '807')
!
!   else
!       call cbeam3_asbly_dynamic (Elem,&
!                                  Node,&
!                                  Coords,&
!                                  Psi0,&
!                                  PosDefor,&
!                                  PsiDefor,&
!                                  PosDotDefor,&
!                                  PsiDotDefor,&
!                                  0.d0*PosDefor,&
!                                  0.d0*PsiDefor,&
!                                  F0+Ftime(iStep+1)*Fa,&
!                                  Vrel(iStep+1,:),&
!                                  VrelDot(iStep+1,:),&
!                                  ms,Mglobal,&
!                                  Mvel,&
!                                  cs,Cglobal,&
!                                  Cvel,&
!                                  ks,Kglobal,&
!                                  fs,Fglobal,&
!                                  Qglobal,&
!                                  Options,&
!                                  Cao)
!   end if
!
! ! Compute admissible error.
!       MinDelta=Options%MinDelta
!
! ! Compute the residual.
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
!       Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
!
!
! ! Check convergence.
!     !   if (maxval(abs(DX)).lt.MinDelta) then
!     residual = maxval(abs(DX))
!     ! if (Iter > 1) then
!       if (abs((residual - previous_residual)/previous_residual).lt.MinDelta) then
!         if (Options%PrintInfo) then
!           write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         end if
!         exit
!       end if
!     end if
!     previous_residual = residual
!
! ! Compute Jacobian
!       call sparse_zero (as,Asys)
!       call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
!       call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
!       call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))
!
! ! Calculation of the correction.
!       call lu_sparse(as,Asys,-Qglobal,DX)
!
!       X    = X     + DX
!       dXdt = dXdt  + gamma/(beta*dt)*DX
!       dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
!     end do
!
! ! Update nodal positions and velocities on the current converged time step.
!     call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!   end do
!
!   deallocate (ListIN,Mvel,Cvel)
!   deallocate (Asys,Fglobal,Mglobal)
!   deallocate (Kglobal,Cglobal,Qglobal)
!   deallocate (X,DX,dXdt,dXddt)
!   deallocate (Displ,Veloc)
!  end subroutine cbeam3_solv_nlndyn_step
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !-> Subroutine CBEAM3_SOLV_NLNDYN_STEP
! !
! !-> Description:
! !
! !    Linear dynamic solution of multibeam problem under applied forces.
! !
! !-> Remarks.-
! !    ADC: this is a modified version, only runs one step
! !         it is useful for coupled simulations, where forces
! !         change every tstep
! !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine cbeam3_solv_nlndyn_step (iOut,&
!                                      NumDof,&
!                                      dt,&
!                                      Elem,&
!                                      Node,&
!                                      F0,&
!                                      Fdyn,&
!                                      Vrel,&
!                                      VrelDot,&
!                                      Coords,&
!                                      Psi0,&
!                                      PosDefor,&
!                                      PsiDefor,&
!                                      PosDotDefor,&
!                                      PsiDotDefor,&
!                                      Options)
!   use lib_fem
!   use lib_rot
!   use lib_rotvect
!   use lib_lu
!   use lib_out
!   use lib_sparse
!   use lib_lu
!   use cbeam3_asbly
!   use lib_xbeam
!
! ! I/O Variables.
!   integer,      intent(in)   :: iOut              ! Output file.
!   integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
!   real(8),      intent(in)   :: dt                ! time step
!   type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
!   type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
!   real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
!   real(8),      intent(in)   :: Fdyn      (:,:)  ! Amplitude of the dynamic nodal forces.
!   real(8),      intent(in)   :: Vrel      (:)   ! Time history of the velocities of the reference frame.
!   real(8),      intent(in)   :: VrelDot   (:)   ! Time history of the accelerations of the reference frame.
!   real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
!   real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
!   real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
!   real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
!   real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
!   real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
!   type(xbopts),intent(in)    :: Options           ! Solver parameters.
!
! ! Local variables.
!   real(8):: beta,gamma                     ! Newmark coefficients.
!   integer:: k                              ! Counters.
!   integer:: iStep,Iter                     ! Counters on time steps and subiterations.
!   real(8):: MinDelta                       ! Value of Delta for convergence.
!   integer:: NumE(1)                        ! Number of elements in the model.
!   integer:: NumN                           ! Number of nodes in the model.
!
!   real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
!   real(8),allocatable:: X(:), DX(:)
!
!   integer,allocatable::  ListIN     (:)    ! List of independent nodes.
!
!   ! Define variables for system matrices.
!   integer:: as,cs,ks,ms,fs
!   type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
!   type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
!   type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
!   type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
!   type(sparse),allocatable:: Fglobal(:)
!   real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
!   real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
!   real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.
!
!   ! Define vaiables for output information.
!   character(len=80)  ::  Text          ! Text with current time information to print out.
!   type(outopts)      ::  OutOptions    ! Output options.
!   real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
!   real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.
!
!   ! Rotation operator for body-fixed frame using quaternions
!   real(8) :: Cao(3,3)
!   real(8) :: Quat(4)
!   real(8) :: Temp(4,4)
!
! ! Initialize.
!   NumN=size(Node)
!   NumE(1)=size(Elem)
!   allocate (ListIN (NumN));
!   do k=1,NumN
!     ListIN(k)=Node(k)%Vdof
!   end do
!   gamma=0.5d0+Options%NewmarkDamp
!   beta =0.25d0*(gamma+0.5d0)*(gamma+0.5d0)
!
! ! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
!   call sparse_allocate (Asys   , NumDof, NumDof)
!   call sparse_allocate (Mglobal, NumDof, NumDof)
!   call sparse_allocate (Cglobal, NumDof, NumDof)
!   call sparse_allocate (Kglobal, NumDof, NumDof)
!   call sparse_allocate (Fglobal, NumDof, NumDof)
!   allocate (Qglobal(NumDof));   Qglobal= 0.d0
!   allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
!   allocate (Cvel   (NumDof,6)); Cvel   = 0.d0
!
!   allocate (X     (NumDof)); X      = 0.d0
!   allocate (DX    (NumDof)); DX     = 0.d0
!   allocate (dXdt  (NumDof)); dXdt   = 0.d0
!   allocate (dXddt (NumDof)); dXddt  = 0.d0
!
! ! Compute system information at initial condition.
!   allocate (Veloc(NumN,6)); Veloc=0.d0
!   allocate (Displ(NumN,6)); Displ=0.d0
!
! ! Allocate quaternions and rotation operator for initially undeformed system
!   Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4
!
! ! Extract initial displacements and velocities.
!   call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)
!
! ! Compute initial acceleration (we are neglecting qdotdot in Kmass).
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
!     &                            0.d0*PosDefor,0.d0*PsiDefor,F0+Fdyn,Vrel(:),VrelDot(:),         &
!     &                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)
!
!       !print*, 'size of Mglobal:', sparse_max_index(ms, Mglobal)
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof, options, Cao), '749')
!   else
!       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
!     &                            0.d0*PosDefor,0.d0*PsiDefor,F0+Fdyn,Vrel,VrelDot,         &
!     &                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao&
!                                  )
!   end if
!
!   Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Fdyn,NumDof,Filter=ListIN))
!
!   call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
!   call lu_sparse(as,Asys,-Qglobal,dXddt)
!
! ! Loop in the time steps.
!   ! do iStep=1,size(Time)-1
!     ! dt= Time(iStep+1)-Time(iStep)
!     ! if (Options%PrintInfo) then
!     !   call out_time(iStep,Time(iStep+1),Text)
!     !   write (*,'(5X,A,$)') trim(Text)
!     ! end if
! ! Update transformation matrix for given angular velocity
!     call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(4:6))*dt),Temp)
!     Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(4:6))*dt),Quat))
!     Cao = xbeam_Rot(Quat)
!
! ! Predictor step.
!     X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
!     dXdt = dXdt + (1.d0-gamma)*dt*dXddt
!     dXddt= 0.d0
!
! ! Iteration until convergence.
!     do Iter=1,Options%MaxIterations+1
!       if (Iter.gt.Options%MaxIterations) then
!         write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         STOP 'Solution did not converge (18235)'
!       end if
!
! ! Update nodal positions and velocities .
!       call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!
! ! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
!       Qglobal= 0.d0
!       Mvel   = 0.d0
!       Cvel   = 0.d0
!       call sparse_zero (ms,Mglobal)
!       call sparse_zero (cs,Cglobal)
!       call sparse_zero (ks,Kglobal)
!       call sparse_zero (fs,Fglobal)
!
!   if (options%gravity_on .eqv. .TRUE.) then
!       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
!                        0.d0*PosDefor,0.d0*PsiDefor,Fdyn+F0,Vrel,VrelDot, &
!     &                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao&
!                                  )
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,&
!                                         NumDof,&
!                                         cbeam3_asbly_gravity_dynamic(NumDof, options, Cao), '807')
!   else
!       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
!                        0.d0*PosDefor,0.d0*PsiDefor,F0+Fdyn,Vrel,VrelDot, &
!     &                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao&
!                                  )
!   end if
!
! ! Compute admissible error.
!       MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))
!
! ! Compute the residual.
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot)
!       Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Fdyn,NumDof,Filter=ListIN))
!
!
! ! Check convergence.
!       if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
!         if (Options%PrintInfo) then
!           write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         end if
!         exit
!       end if
!
! ! Compute Jacobian
!       call sparse_zero (as,Asys)
!       call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
!       call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
!       call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))
!
! ! Calculation of the correction.
!       call lu_sparse(as,Asys,-Qglobal,DX)
!
!       X    = X     + DX
!       dXdt = dXdt  + gamma/(beta*dt)*DX
!       dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
!     end do
!
! ! Update nodal positions and velocities on the current converged time step.
!     call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!
!   deallocate (ListIN,Mvel,Cvel)
!   deallocate (Asys,Fglobal,Mglobal)
!   deallocate (Kglobal,Cglobal,Qglobal)
!   deallocate (X,DX,dXdt,dXddt)
!   deallocate (Displ,Veloc)
!   return
!  end subroutine cbeam3_solv_nlndyn_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem under applied forces.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn_old (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,              &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
  use lib_lu
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)  ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)    ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical(c_bool),      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
  real(8),allocatable:: X(:), DX(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for system matrices.
  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
  real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  type(outopts)      ::  OutOptions    ! Output options.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Cao(3,3)
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

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
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

! sm:
  !allocate( Fdyn(NumN,6,size(Time)) ) ! temporary
  !do iStep=1,size(Time)
  !  Fdyn(:,:,iStep)=Ftime(iStep)*Fa
  !end do

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

! Extract initial displacements and velocities.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)

! Compute initial acceleration (we are neglecting qdotdot in Kmass).
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
&                            0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),VrelDot(1,:),         &
&                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

  Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))



  call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
  call lu_sparse(as,Asys,-Qglobal,dXddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    if (Options%PrintInfo) then
      call out_time(iStep,Time(iStep+1),Text)
      ! write (*,'(5X,A,$)') trim(Text)
    end if
! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
    dXdt = dXdt + (1.d0-gamma)*dt*dXddt
    dXddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) then
        write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        STOP 'Solution did not converge (18235)'
      end if

! Update nodal positions and velocities .
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qglobal= 0.d0
      Mvel   = 0.d0
      Cvel   = 0.d0
      call sparse_zero (ms,Mglobal)
      call sparse_zero (cs,Cglobal)
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)


       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                     &
 &                                0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:), &
 &                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)


! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))

! Compute the residual.
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
      Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))


! Check convergence.
      if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
        if (Options%PrintInfo) then
          write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        end if
        exit
      end if

! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))

! Calculation of the correction.
      call lu_sparse(as,Asys,-Qglobal,DX)

      X    = X     + DX
      dXdt = dXdt  + gamma/(beta*dt)*DX
      dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
    end do

! Update nodal positions and velocities on the current converged time step.
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Write output to export to main program (obsolete!).
    PosPsiTime(iStep+1,1:3)= PosDefor(NumN,:)
    PosPsiTime(iStep+1,4:6)= PsiDefor(NumE(1),Elem(NumE(1))%NumNodes,:)

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do

  end do
  VelocTime = 0.0d0

  deallocate (ListIN,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (X,DX,dXdt,dXddt)
  deallocate (Displ,Veloc)
  return
 end subroutine cbeam3_solv_nlndyn_old



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem under applied forces for
!    optimal control
!
!-> Remarks.-
!
!    The solver is identical to cbeam3_solv_nlndyn but the input Fa, Ftime
!    have been replaced by a single 3 rank array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine cbeam3_solv_nlndyn_opt_control (iOut,NumDof,Time,Elem,Node,F0,Fdyn,      &
! &                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
! &                               PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,DynOut, &
! &                               OutGrids,Options)
!   use lib_fem
!   use lib_rot
!   use lib_rotvect
!   use lib_lu
!   use lib_out
!   use lib_sparse
! !#ifdef NOLAPACK
!   use lib_lu
! !#else
! !  use interface_lapack
! !#endif
!   use cbeam3_asbly
!   use lib_xbeam
!
!
! ! I/O Variables.
!   integer,      intent(in)   :: iOut              ! Output file.
!   integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
!   real(8),      intent(in)   :: Time      (:)     ! Time steps.
!   type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
!   type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
!   real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
!   !real(8),      intent(in)   :: Fa        (:,:)  ! Amplitude of the dynamic nodal forces.
!   !real(8),      intent(in)   :: Ftime     (:)    ! Time history of the applied forces.
!   real(8),      intent(in)   :: Fdyn      (:,:,:) ! applied dynamic force of size (NumNodes, 6, NumSteps+1)
!   real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
!   real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
!   real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
!   real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
!   real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
!   real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
!   real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
!   real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
!   real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
!   real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
!   real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
!   logical,      intent(inout):: OutGrids  (:)     ! Output grids.
!   type(xbopts),intent(in)    :: Options           ! Solver parameters.
!
! ! Local variables.
!   real(8):: beta,gamma                     ! Newmark coefficients.
!   real(8):: dt                             ! Time step
!   integer:: k                              ! Counters.
!   integer:: iStep,Iter                     ! Counters on time steps and subiterations.
!   real(8):: MinDelta                       ! Value of Delta for convergence.
!   integer:: NumE(1)                        ! Number of elements in the model.
!   integer:: NumN                           ! Number of nodes in the model.
!
!   real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
!   real(8),allocatable:: X(:), DX(:)
!
!   integer,allocatable::  ListIN     (:)    ! List of independent nodes.
!
!   ! Define variables for system matrices.
!   integer:: as,cs,ks,ms,fs
!   type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
!   type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
!   type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
!   type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
!   type(sparse),allocatable:: Fglobal(:)
!   real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
!   real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
!   real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.
!
!   ! Define vaiables for output information.
!   character(len=80)  ::  Text          ! Text with current time information to print out.
!   type(outopts)      ::  OutOptions    ! Output options.
!   real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
!   real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.
!
!   ! Rotation operator for body-fixed frame using quaternions
!   real(8) :: Cao(3,3)
!   real(8) :: Quat(4)
!   real(8) :: Temp(4,4)
!
! ! Initialize.
!   NumN=size(Node)
!   NumE(1)=size(Elem)
!   allocate (ListIN (NumN));
!   do k=1,NumN
!     ListIN(k)=Node(k)%Vdof
!   end do
!   gamma=1.d0/2.d0+Options%NewmarkDamp
!   beta =1.d0/4.d0*(gamma+0.5d0)*(gamma+0.5d0)
!
! ! Allocate memory for solver (Use a conservative estimate of the size of the matrices).
!   allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
!   allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
!   allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
!   allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
!   allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
!   allocate (Qglobal(NumDof));   Qglobal= 0.d0
!   allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
!   allocate (Cvel   (NumDof,6)); Cvel   = 0.d0
!
!   allocate (X     (NumDof)); X      = 0.d0
!   allocate (DX    (NumDof)); DX     = 0.d0
!   allocate (dXdt  (NumDof)); dXdt   = 0.d0
!   allocate (dXddt (NumDof)); dXddt  = 0.d0
!
! ! Compute system information at initial condition.
!   allocate (Veloc(NumN,6)); Veloc=0.d0
!   allocate (Displ(NumN,6)); Displ=0.d0
!
! ! sm:
!   !allocate( Fdyn(NumN,6,size(Time)) ) ! temporary
!   !do iStep=1,size(Time)
!   !  Fdyn(:,:,iStep)=Ftime(iStep)*Fa
!   !end do
!
! ! Allocate quaternions and rotation operator for initially undeformed system
!   Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4
!
! ! Extract initial displacements and velocities.
!   call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)
!
! ! sm: substitute applied force term with Fdyn
! ! Compute initial acceleration (we are neglecting qdotdot in Kmass).
! !  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
! !&                            0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),VrelDot(1,:),         &
! !&                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)
! !
! !  Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))
!   call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,           &
! &                            0.d0*PosDefor,0.d0*PsiDefor,F0+Fdyn(:,:,1),Vrel(1,:),VrelDot(1,:),         &
! &                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)
!
!   Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Fdyn(:,:,1),NumDof,Filter=ListIN))
!
!
!   call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
! !#ifdef NOLAPACK
!   call lu_sparse(as,Asys,-Qglobal,dXddt)
! !#else
! !  call lapack_sparse (as,Asys,-Qglobal,dXddt)
! !#endif
!
! ! Loop in the time steps.
!   do iStep=1,size(Time)-1
!     dt= Time(iStep+1)-Time(iStep)
!     if (Options%PrintInfo) then
!       call out_time(iStep,Time(iStep+1),Text)
!       write (*,'(5X,A,$)') trim(Text)
!     end if
! ! Update transformation matrix for given angular velocity
!     call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
!     Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
!     Cao = xbeam_Rot(Quat)
!
! ! Predictor step.
!     X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
!     dXdt = dXdt + (1.d0-gamma)*dt*dXddt
!     dXddt= 0.d0
!
! ! Iteration until convergence.
!     do Iter=1,Options%MaxIterations+1
!       if (Iter.gt.Options%MaxIterations) then
!         write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         STOP 'Solution did not converge (18235)'
!       end if
!
! ! Update nodal positions and velocities .
!       call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!
! ! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
!       Qglobal= 0.d0
!       Mvel   = 0.d0
!       Cvel   = 0.d0
!       call sparse_zero (ms,Mglobal)
!       call sparse_zero (cs,Cglobal)
!       call sparse_zero (ks,Kglobal)
!       call sparse_zero (fs,Fglobal)
!
!  ! sm: substitute applied force term with Fdyn
!  !      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                  &
!  !&                                0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),    &
!  !&                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)
!       call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                  &
! &                                0.d0*PosDefor,0.d0*PsiDefor,F0+Fdyn(:,:,iStep+1),Vrel(iStep+1,:),VrelDot(iStep+1,:),    &
! &                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)
!
! ! Compute admissible error.
!       MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))
!
! ! Compute the residual.
!       Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
!       ! sm: substitute applied force term with Fdyn
!       ! Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
!       Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Fdyn(:,:,iStep+1),NumDof,Filter=ListIN))
!
! ! Check convergence.
!       if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
!         if (Options%PrintInfo) then
!           write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
!         end if
!         exit
!       end if
!
! ! Compute Jacobian
!       call sparse_zero (as,Asys)
!       call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
!       call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
!       call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))
!
! ! Calculation of the correction.
! !#ifdef NOLAPACK
!       call lu_sparse(as,Asys,-Qglobal,DX)
! !#else
! !      call lapack_sparse (as,Asys,-Qglobal,DX)
! !#endif
!       X    = X     + DX
!       dXdt = dXdt  + gamma/(beta*dt)*DX
!       dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
!     end do
!
! ! Update nodal positions and velocities on the current converged time step.
!     call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
!
! ! Write output to export to main program (obsolete!).
!     PosPsiTime(iStep+1,1:3)= PosDefor(NumN,:)
!     PosPsiTime(iStep+1,4:6)= PsiDefor(NumE(1),Elem(NumE(1))%NumNodes,:)
!
!     do k=1,NumN
!         DynOut(iStep*NumN+k,:) = PosDefor(k,:)
!     end do
!
!   end do
!
!   deallocate (ListIN,Mvel,Cvel)
!   deallocate (Asys,Fglobal,Mglobal)
!   deallocate (Kglobal,Cglobal,Qglobal)
!   deallocate (X,DX,dXdt,dXddt)
!   deallocate (Displ,Veloc)
!   return
!  end subroutine cbeam3_solv_nlndyn_opt_control



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_NLNDYN_ACCEL
!
!-> Description:
!
!    Nonlinear dynamic solution of multibeam problem under applied forces
!    with accelerations added as input/output variable.
!
!-> Remarks: Initial accelerations should be calculated outside this routine.
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_nlndyn_accel (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,        &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosDDot,PsiDDot,             &
&                               PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_fem
  use lib_rot
  use lib_rotvect
  use lib_lu
  use lib_out
  use lib_sparse
  use lib_lu
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem), intent(in)   :: Elem      (:)     ! Element information.
  type(xbnode), intent(in)   :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Initial coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Initial CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Current coordinates of the grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Current CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDDot(:,:)  ! Current accelerations of the coordinates of the grid points
  real(8),      intent(inout):: PsiDDot(:,:,:)! Current accelerations of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical(c_bool),      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma                     ! Newmark coefficients.
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep,Iter                     ! Counters on time steps and subiterations.
  real(8):: MinDelta                       ! Value of Delta for convergence.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: dXdt(:),dXddt(:)  ! Generalized coordinates and derivatives.
  real(8),allocatable:: X(:), DX(:)

  integer,allocatable::  ListIN     (:)    ! List of independent nodes.

  ! Define variables for system matrices.
  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable:: Asys   (:)    ! System matrix for implicit Newmark method.
  type(sparse),allocatable:: Cglobal(:)    ! Sparse damping matrix.
  type(sparse),allocatable:: Kglobal(:)    ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable:: Mglobal(:)    ! Global mass matrix in sparse storage.
  type(sparse),allocatable:: Fglobal(:)
  real(8),allocatable::      Qglobal(:)    ! Global vector of discrete generalize forces.
  real(8),allocatable::      Mvel(:,:)     ! Mass and damping from the motions of reference system.
  real(8),allocatable::      Cvel(:,:)     ! Mass and damping from the motions of reference system.

  ! Define vaiables for output information.
  character(len=80)  ::  Text          ! Text with current time information to print out.
  type(outopts)      ::  OutOptions    ! Output options.
  real(8),allocatable::  Displ(:,:)    ! Current nodal displacements/rotations.
  real(8),allocatable::  Veloc(:,:)    ! Current nodal velocities.

  ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Cao(3,3)
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

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
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)
  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X     (NumDof)); X      = 0.d0
  allocate (DX    (NumDof)); DX     = 0.d0
  allocate (dXdt  (NumDof)); dXdt   = 0.d0
  allocate (dXddt (NumDof)); dXddt  = 0.d0

! Compute system information at initial condition.
  allocate (Veloc(NumN,6)); Veloc=0.d0
  allocate (Displ(NumN,6)); Displ=0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

! Extract initial displacements, velocities AND accelerations.
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,X,dXdt)
  call cbeam3_solv_accel2state (Node,PosDDot,PsiDDot,dXddt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt= Time(iStep+1)-Time(iStep)
    if (Options%PrintInfo) then
      call out_time(iStep,Time(iStep+1),Text)
      ! write (*,'(5X,A,$)') trim(Text)
    end if
! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    X    = X + dt*dXdt + (0.5d0-beta)*dt*dt*dXddt
    dXdt = dXdt + (1.d0-gamma)*dt*dXddt
    dXddt= 0.d0

! Iteration until convergence.
    do Iter=1,Options%MaxIterations+1
      if (Iter.gt.Options%MaxIterations) STOP 'Solution did not converge (18235)'

! Update nodal positions and velocities .
      call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)

! Compute system functionals and matrices. (Use initial accelerations for Kgyr).
      Qglobal= 0.d0
      Mvel   = 0.d0
      Cvel   = 0.d0
      call sparse_zero (ms,Mglobal)
      call sparse_zero (cs,Cglobal)
      call sparse_zero (ks,Kglobal)
      call sparse_zero (fs,Fglobal)

      call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                  &
&                                0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),    &
&                                ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

! Compute admissible error.
      MinDelta=Options%MinDelta*max(1.d0,maxval(abs(Qglobal)))

! Compute the residual.
      Qglobal= Qglobal + sparse_matvmul(ms,Mglobal,NumDof,dXddt) + matmul(Mvel,Vreldot(iStep+1,:))
      Qglobal= Qglobal - sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))

! Check convergence.
      if (maxval(abs(DX)+abs(Qglobal)).lt.MinDelta) then
        if (Options%PrintInfo) then
          write (*,'(5X,A,I4,A,1PE12.3)') 'Subiteration',Iter, '  Delta=', maxval(abs(Qglobal))
        end if
        exit
      end if

! Compute Jacobian
      call sparse_zero (as,Asys)
      call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=1.d0)
      call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma/(beta*dt))
      call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0/(beta*dt*dt))

! Calculation of the correction.
      call lu_sparse(as,Asys,-Qglobal,DX)

      X    = X     + DX
      dXdt = dXdt  + gamma/(beta*dt)*DX
      dXddt= dXddt + 1.d0/(beta*dt*dt)*DX
    end do

! Update nodal positions and velocities on the current converged time step.
    call cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor)
    call cbeam3_solv_state2accel(Elem,Node,dXddt,PosDDot,PsiDDot)

!!! Postprocesing (for single cantilever beams) !!!
! Store data in output variables (V_B, Omega_B).
!    if (any(OutGrids)) then
!      Veloc=0.d0
!      Displ=0.d0
!
!      do k=1,NumN
!        Displ(k,1:3)= PosDefor(k,1:3)-Coords(k,1:3)
!      end do
!
!      Veloc(1,1:3)= PosDotDefor(1,:)+rot_cross(Vrel(iStep,4:6), PosDefor(1,:))+Vrel(iStep,1:3)
!      do k=2,NumN
!        Displ(k,4:6)=PsiDefor(k-1,2,:)
!        CBa=rotvect_psi2mat(PsiDefor(k-1,2,:))
!        Veloc(k,1:3)= matmul(CBa, PosDotDefor(k,:) + Vrel(iStep,1:3)      &
!&                               + rot_cross(Vrel(iStep,4:6),PosDefor(k,:)))
!        Veloc(k,4:6)= matmul(rotvect_psi2rot(PsiDefor(k-1,2,:)),PsiDotDefor(k-1,2,:)) &
!&                   + matmul(CBa,Vrel(iStep,4:6))
!      end do
!
!!  Write output information in output file.
!      OutOptions%PrintDispl=.true.
!      OutOptions%PrintVeloc=.true.
!      call out_title   (iOut,trim(Text))
!      call out_outgrid (iOut,'NODE',OutOptions,1,NumE,6,OutGrids,DISPL=Displ,VELOC=Veloc)
!
!! Write output to export to main program (obsolete!).
!      VelocTime (iStep+1,:)= Veloc (1:NumN,3)
!    end if

! Write output to export to main program (obsolete!).
    PosPsiTime(iStep+1,1:3)= PosDefor(NumN,:)
    PosPsiTime(iStep+1,4:6)= PsiDefor(NumE(1),Elem(NumE(1))%NumNodes,:)

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:)
    end do

  end do

  VelocTime = 0.0d0
  deallocate (ListIN,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (X,DX,dXdt,dXddt)
  deallocate (Displ,Veloc)
  return
 end subroutine cbeam3_solv_nlndyn_accel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_LINDYN
!
!-> Description:
!
!    Linear dynamic solution of multibeam problem for given applied forces.
!
!-> Remarks.-
!
!   1) This routine assumes that the motion of the reference system is
!      prescribed.
!
!   2) Coords/Psi0 defines the unloaded geometry, while PosDefor/PsiDefor
!      brings the initial geometry (after static equilibrium). F0 should
!      include the static loading for that equilibrium.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_lindyn (iOut,NumDof,Time,Elem,Node,F0,Fa,Ftime,              &
&                               Vrel,VrelDot,Coords,Psi0,PosDefor,PsiDefor,          &
&                               PosDotDefor,PsiDotDefor,PosPsiTime,VelocTime,DynOut, &
&                               OutGrids,Options)
  use lib_rot
  use lib_fem
  use lib_sparse
  use lib_out
  use lib_lu
  use cbeam3_asbly
  use lib_xbeam

! I/O Variables.
  integer,      intent(in)   :: iOut              ! Output file.
  integer,      intent(in)   :: NumDof            ! Number of independent DoFs.
  real(8),      intent(in)   :: Time      (:)     ! Time steps.
  type(xbelem),intent(in)    :: Elem      (:)     ! Element information.
  type(xbnode),intent(in)    :: Node      (:)     ! Nodal information.
  real(8),      intent(in)   :: F0        (:,:)   ! Applied static nodal forces.
  real(8),      intent(in)   :: Fa        (:,:)   ! Amplitude of the dynamic nodal forces.
  real(8),      intent(in)   :: Ftime     (:)     ! Time history of the applied forces.
  real(8),      intent(in)   :: Vrel      (:,:)   ! Time history of the velocities of the reference frame.
  real(8),      intent(in)   :: VrelDot   (:,:)   ! Time history of the accelerations of the reference frame.
  real(8),      intent(in)   :: Coords    (:,:)   ! Undeformed coordinates of the grid points.
  real(8),      intent(in)   :: Psi0      (:,:,:) ! Undeformed CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDefor  (:,:)   ! Initial/final position vector of grid points
  real(8),      intent(inout):: PsiDefor  (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),      intent(inout):: PosDotDefor(:,:)  ! Current time derivatives of the coordinates of the grid points
  real(8),      intent(inout):: PsiDotDefor(:,:,:)! Current time derivatives of the CRV of the nodes in the elements.
  real(8),      intent(out)  :: PosPsiTime(:,:)   ! Time-history of the position/rotation at selected nodes.
  real(8),      intent(out)  :: VelocTime (:,:)   ! Time-history of the time derivatives at selected nodes.
  real(8),      intent(out)  :: DynOut    (:,:)   ! Time-history of displacement of all nodes.
  logical,      intent(inout):: OutGrids  (:)     ! Output grids.
  type(xbopts),intent(in)    :: Options           ! Solver parameters.

! Local variables.
  real(8):: beta,gamma
  real(8):: dt                             ! Time step
  integer:: k                              ! Counters.
  integer:: iStep                          ! Current time step.
  integer:: NumE(1)                        ! Number of elements in the model.
  integer:: NumN                           ! Number of nodes in the model.

  real(8),allocatable:: X0(:),DX(:),DXDt(:),DXDDt(:)
  integer,allocatable:: ListIN     (:)     ! List of independent nodes.
  real(8),allocatable:: DeltaPos   (:,:)   ! Initial/final position vector of grid points
  real(8),allocatable:: DeltaPsi   (:,:,:) ! Initial/final CRV of the nodes in the elements.
  real(8),allocatable:: DeltaPosDot(:,:)   ! Current time derivatives of the coordinates of the grid points
  real(8),allocatable:: DeltaPsiDot(:,:,:) ! Current time derivatives of the CRV of the nodes in the elements.

  integer:: as,cs,ks,ms,fs
  type(sparse),allocatable  :: Asys   (:)     ! System matrix for implicit Newmark method.
  type(sparse),allocatable  :: Cglobal(:)     ! Sparse damping matrix.
  type(sparse),allocatable  :: Kglobal(:)     ! Global stiffness matrix in sparse storage.
  type(sparse),allocatable  :: Mglobal(:)     ! Global mass matrix in sparse storage.
  type(sparse),allocatable  :: Fglobal(:)
  real(8),allocatable       :: Qglobal(:)     ! Global vector of discrete generalize forces.
  real(8),allocatable       :: Mvel(:,:)      ! Mass and damping from the motions of reference system.
  real(8),allocatable       :: Cvel(:,:)      ! Mass and damping from the motions of reference system.

  real(8) :: Cao(3,3)           ! Rotation operator for body-fixed frame using quaternions
  real(8) :: Quat(4)
  real(8) :: Temp(4,4)

  character(len=80)  :: Text          ! Text with current time information to print out.
  type(outopts)      :: OutOptions    ! Output options.
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
  allocate (Asys   (DimMat*NumDof)); call sparse_zero (as,Asys)
  allocate (Mglobal(DimMat*NumDof)); call sparse_zero (ms,Mglobal)
  allocate (Cglobal(DimMat*NumDof)); call sparse_zero (cs,Cglobal)
  allocate (Kglobal(DimMat*NumDof)); call sparse_zero (ks,Kglobal)
  allocate (Fglobal(DimMat*NumDof)); call sparse_zero (fs,Fglobal)

  allocate (Qglobal(NumDof));   Qglobal= 0.d0
  allocate (Mvel   (NumDof,6)); Mvel   = 0.d0
  allocate (Cvel   (NumDof,6)); Cvel   = 0.d0

  allocate (X0     (NumDof));   X0     = 0.d0
  allocate (DX     (NumDof));   DX     = 0.d0
  allocate (DXDt   (NumDof));   DXDt   = 0.d0
  allocate (DXDDt  (NumDof));   DXDDt  = 0.d0

  allocate (Displ      (NumN,         6));      Displ   =    0.d0
  allocate (Veloc      (NumN,         6));      Veloc   =    0.d0
  allocate (DeltaPos   (NumN,         3));      DeltaPos=    0.d0
  allocate (DeltaPsi   (NumE(1),MaxElNod,3));   DeltaPsi=    0.d0
  allocate (DeltaPosDot(NumN,         3));      DeltaPosDot= 0.d0
  allocate (DeltaPsiDot(NumE(1),MaxElNod,3));   DeltaPsiDot= 0.d0

! Allocate quaternions and rotation operator for initially undeformed system
  Quat = (/1.d0,0.d0,0.d0,0.d0/); Cao = Unit; Temp = Unit4

!Find initial conditions
  call cbeam3_solv_disp2state (Node,Coords,Psi0,PosDotDefor,PsiDotDefor,X0,DXDt)
  call cbeam3_solv_disp2state (Node,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,DX,DXDt)

  DX=DX-X0

! Compute tangent matrices at initial time.
  call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,Coords,Psi0,PosDotDefor,                                   &
&                            PsiDotDefor,0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(1)*Fa,Vrel(1,:),VrelDot(1,:),   &
&                            ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Unit)

! Find initial acceleration.
  Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(1)*Fa,NumDof,Filter=ListIN))   &
&          - matmul(Mvel,VrelDot(1,:)) - matmul(Cvel,Vrel(1,:))                             &
&          - sparse_matvmul(cs,Cglobal,NumDof,DXDt) - sparse_matvmul(ks,Kglobal,NumDof,DX)

  call sparse_addsparse(0,0,ms,Mglobal,as,Asys)
  call lu_sparse(as,Asys,Qglobal,DXDDt)

! Loop in the time steps.
  do iStep=1,size(Time)-1
    dt     =Time(iStep+1)-Time(iStep)
    call out_time(iStep,Time(iStep+1),Text)
    write (*,'(5X,A)') trim(Text)

! Update transformation matrix for given angular velocity
    call lu_invers ((Unit4+0.25d0*xbeam_QuadSkew(Vrel(iStep+1,4:6))*(Time(iStep+1)-Time(iStep))),Temp)
    Quat=matmul(Temp,matmul((Unit4-0.25d0*xbeam_QuadSkew(Vrel(iStep,4:6))*(Time(iStep+1)-Time(iStep))),Quat))
    Cao = xbeam_Rot(Quat)

! Predictor step.
    DX= DX + dt*DXDt + (0.5d0-beta)*dt*dt*DXDDt
    DXDt=DXDt + (1-gamma)*dt*DXDDt

! Compute system functionals and matrices. Only updating Felast.
    call sparse_zero (ms,Mglobal); Mvel=0.d0
    call sparse_zero (cs,Cglobal); Cvel=0.d0
    call sparse_zero (ks,Kglobal);
    call sparse_zero (fs,Fglobal);

    call cbeam3_asbly_dynamic (Elem,Node,Coords,Psi0,PosDefor,PsiDefor,PosDotDefor,PsiDotDefor,                &
&                              0.d0*PosDefor,0.d0*PsiDefor,F0+Ftime(iStep+1)*Fa,Vrel(iStep+1,:),VrelDot(iStep+1,:),  &
&                              ms,Mglobal,Mvel,cs,Cglobal,Cvel,ks,Kglobal,fs,Fglobal,Qglobal,Options,Cao)

! Compute right-hand side of the equation.
    Qglobal= sparse_matvmul(fs,Fglobal,NumDof,fem_m2v(F0+Ftime(iStep+1)*Fa,NumDof,Filter=ListIN))
    Qglobal= Qglobal - matmul(Mvel,Vreldot(iStep+1,:)) &
&                    - matmul(Cvel,Vrel   (iStep+1,:))
    Qglobal= Qglobal - sparse_matvmul(cs,Cglobal,NumDof,DXDt) &
&                    - sparse_matvmul(ks,Kglobal,NumDof,DX)

! Compute left-hand side of the equation (if dt=constant then Asys=constant too, but this was not
! assumed here).
    call sparse_zero (as,Asys)
    call sparse_addsparse(0,0,ms,Mglobal,as,Asys,Factor=1.d0)
    call sparse_addsparse(0,0,cs,Cglobal,as,Asys,Factor=gamma*dt)
    call sparse_addsparse(0,0,ks,Kglobal,as,Asys,Factor=beta*dt*dt)

! Solve equation.
    call lu_sparse(as,Asys,Qglobal,DXDDt)

    DX    = DX   + beta *dt*dt*DXDDt
    DXDt  = DXDt + gamma*dt   *DXDDt

! Update solution and store relevant information at current time step.
    call cbeam3_solv_update_lindyn (Elem,Node,PsiDefor,DX,DXDt,DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)

! Store data in output variables.
    do k=1,NumN
      Veloc(k,1:3)= DeltaPosDot(k,:) + Vrel(iStep,1:3) &
&                 + rot_cross(Vrel(iStep,4:6),PosDefor(k,:)+DeltaPos(k,:))
    end do
    VelocTime (iStep+1, : )= Veloc   (1:NumN,3)

    PosPsiTime(iStep+1,1:3)= Coords(NumN,:) + DeltaPos(NumN,:)
    PosPsiTime(iStep+1,4:6)= Psi0(NumE(1),Elem(NumE(1))%NumNodes,:) &
&                          + DeltaPsi(NumE(1),Elem(NumE(1))%NumNodes,:)

!  Write output information in output file.
    if (any(OutGrids)) then
      do k=1,NumN
        if (OutGrids(k)) then
          Displ(k,1:3)= DeltaPos(k,:)
          Displ(k,4:6)= DeltaPsi(Node(k)%Master(1),Node(k)%Master(2),:)
        end if
      end do
      call out_title   (iOut,trim(Text))
      call out_outgrid (iOut,'NODE',OutOptions,1,NumE,6,OutGrids,DISPL=Displ,VELOC=Veloc)
    end if

    do k=1,NumN
        DynOut(iStep*NumN+k,:) = PosDefor(k,:) + DeltaPos(k,:)
    end do

  end do

! Write information at last time step.
  PosDefor= Coords + DeltaPos
  PsiDefor= Psi0 + DeltaPsi

  deallocate (ListIn,Mvel,Cvel)
  deallocate (Asys,Fglobal,Mglobal)
  deallocate (Kglobal,Cglobal,Qglobal)
  deallocate (DX,DXdt,DXDDt)
  deallocate (Displ,Veloc)
  deallocate (DeltaPos,DeltaPsi,DeltaPosDot,DeltaPsiDot)
  return
 end subroutine cbeam3_solv_lindyn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_UPDATE_STATIC
!
!-> Description:
!
!    Update results from static solution.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_update_static (Elem,Node,Psi0,DeltaX,PosDefor,PsiDefor)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem      (:)       ! Element information.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Psi0      (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: DeltaX    (:)       ! Incremental State Vector.
  real(8),intent(inout)   :: PosDefor  (:,:)     ! Current nodal position.
  real(8),intent(inout)   :: PsiDefor  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.
! Store current state in CBEAM3 elements.
  ix=0
  do iNode=1,size(PosDefor,DIM=1)
      !print*, '**'
      !print*,  ' iNode = ', iNode
    iElem=Node(iNode)%Master(1)
      !print*, ' iElem = ', iElem
    ! if ((Node(iNode)%Vdof.ne.0).and.(Elem(iElem)%MemNo.eq.0)) then
    if (Node(iNode)%Vdof.ne.0) then
      k=Node(iNode)%Master(2)

! Nodal displacements.
      PosDefor(iNode,:)   = PosDefor(iNode,:)   + DeltaX(ix+1:ix+3)

! Cartesian rotation vector at master nodes.
      PsiDefor(iElem,k,:) = PsiDefor(iElem,k,:) + DeltaX(ix+4:ix+6)
      ix=ix+6
    end if
  end do

!!! Post-processing.
! Compute rotation vector at slave nodes of CBEAM3 elements.
  do i=1,size(Elem)
    ! if (Elem(i)%MemNo.eq.0) then

! Copy rotation from master node for each slave node.
      do j=1,Elem(i)%NumNodes
        if (Elem(i)%Master(j,1).ne.0) then
          PsiDefor(i,j,:)=PsiDefor(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        end if
      end do

! Include master-to-slave initial rotations from the undeformed configuration.
      call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,PsiDefor(i,:,:))
    ! end if
  end do

  return
 end subroutine cbeam3_solv_update_static

 subroutine cbeam3_solv_update_static_python(num_dof,&
                                             n_node,&
                                             node_master,&
                                             vdof,&
                                             n_elem,&
                                             elem_master,&
                                             elem_num_nodes,&
                                             psi0,&
                                             pos_def,&
                                             psi_def,&
                                             deltax) bind(C)
    use, intrinsic                          :: iso_c_binding
    use lib_fem
    use lib_rotvect
    use lib_bgeom
    use lib_cbeam3

    integer(c_int), intent(IN)              :: num_dof
    integer(c_int), intent(IN)              :: n_node
    integer(c_int), intent(IN)              :: node_master(n_node, 2)
    integer(c_int), intent(IN)              :: vdof(n_node)
    integer(c_int), intent(IN)              :: n_elem
    integer(c_int), intent(IN)              :: elem_master(n_elem, MaxElNod, 2)
    integer(c_int), intent(IN)              :: elem_num_nodes(n_elem)
    real(c_double), intent(IN)              :: psi0(n_elem, MaxElNod, 3)
    real(c_double), intent(INOUT)           :: pos_def(n_node, 3)
    real(c_double), intent(INOUT)           :: psi_def(n_elem, MaxElNod, 3)
    real(c_double), intent(IN)              :: deltax(num_dof+6)

    ! Local variables.
    integer :: i,j,k               ! Counters.
    integer :: ix                  ! Counter on the degrees of freedom.
    integer :: iElem               ! Counter on the elements.
    integer :: iNode               ! Counter on the nodes.

  ! Store current state in CBEAM3 elements.
    ix=0
    do iNode=1, n_node
        print*, 'node'
        print*, iNode
      iElem=node_master(iNode, 1)
      if (vdof(iNode) /= 0) then
        k=node_master(iNode, 2)
        print*, vdof(inode)
        ix = (vdof(inode) - 1)*6
        print*, ix

  ! Nodal displacements.
        pos_def(iNode,:)   = pos_def(iNode,:)   + deltax(ix+1:ix+3)

  ! Cartesian rotation vector at master nodes.
        psi_def(iElem,k,:) = psi_def(iElem,k,:) + deltax(ix+4:ix+6)
        !ix=ix+6
      end if
    end do

  !!! Post-processing.
  ! Compute rotation vector at slave nodes of CBEAM3 elements.
    do i=1, n_elem
      ! Copy rotation from master node for each slave node.
        do j=1, elem_num_nodes(i)
          if (elem_master(i, j, 1) /= 0) then
            psi_def(i,j,:) = psi_def(elem_master(i, j, 1),elem_master(i, j, 2),:)
          end if
        end do

  ! Include master-to-slave initial rotations from the undeformed configuration.
        !call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,PsiDefor(i,:,:))
        call cbeam3_projm2s (elem_num_nodes(i),elem_master(i, :, :),psi0(i,:,:),psi0,psi_def(i,:,:))
    end do

 end subroutine cbeam3_solv_update_static_python


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_UPDATE_lindyn
!
!-> Description:
!
!    Update results from linear dynamic solution.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_update_lindyn (Elem,Node,Psi0,DX,DXDt,Pos,Psi,PosDot,PsiDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem      (:)       ! Element information.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Psi0      (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: DX        (:)       ! Incremental State Vector.
  real(8),intent(in)      :: DXDt      (:)       ! Time derivative of DX.
  real(8),intent(inout)   :: Pos       (:,:)     ! Delta nodal position.
  real(8),intent(inout)   :: Psi       (:,:,:)   ! Delta rotation at element nodes.
  real(8),intent(inout)   :: PosDot    (:,:)     ! Delta nodal velocities.
  real(8),intent(inout)   :: PsiDot    (:,:,:)   ! Delta nodal angular velocities.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current displacement and its time derivative at all nodes and the
! rotations and its first derivatives at the master nodes of each element.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then
      Pos   (iNode,:)  = 0.d0
      Psi   (iElem,k,:)= 0.d0
      PosDot(iNode,:)  = 0.d0
      PsiDot(iElem,k,:)= 0.d0

    ! Unconstrained nodes.
    else
      Pos   (iNode,:)= DX(ix+1:ix+3)
      PosDot(iNode,:)= DXDt(ix+1:ix+3)

      Psi   (iElem,k,:)= DX(ix+4:ix+6)
      PsiDot(iElem,k,:)= DXDt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

!!! Post-processing.
! Compute rotation vector at slave nodes of CBEAM3 elements.
  do i=1,size(Elem)

! Copy rotation and derivative from master node for each slave node.
    do j=1,Elem(i)%NumNodes
      if (Elem(i)%Master(j,1).ne.0) then
        Psi   (i,j,:)= Psi0 (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:) &
&                    + Psi  (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        PsiDot(i,j,:)=PsiDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do

    ! Include master-to-slave initial rotations from the undeformed configuration.
    call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,Psi(i,:,:))

    ! Compute the delta value.
    do j=1,Elem(i)%NumNodes
      if (Elem(i)%Master(j,1).ne.0) then
        Psi   (i,j,:)= Psi (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:) &
&                    - Psi0(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do
  end do


  return
 end subroutine cbeam3_solv_update_lindyn



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_STATE2DISP
!
!-> Description:
!
!    Extract current positiion, orientation and velocities from the
!    state vector.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_state2disp (Elem,Node,Coords,Psi0,X,dXdt,Pos,Psi,PosDot,PsiDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem    (:)       ! Element information.
  type(xbnode),intent(in) :: Node    (:)       ! Nodal information.
  real(8),intent(in)      :: Coords  (:,:)     ! Initial coordinates of the grid points.
  real(8),intent(in)      :: Psi0    (:,:,:)   ! Initial rotation vector at element grids.
  real(8),intent(in)      :: X       (:)       ! Current generalized coordinates.
  real(8),intent(in)      :: dXdt    (:)       ! Time derivatives of X.
  real(8),intent(out)     :: Pos     (:,:)     ! Current nodal position.
  real(8),intent(out)     :: Psi     (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(out)     :: PosDot  (:,:)     ! Current nodal position.
  real(8),intent(out)     :: PsiDot  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current displacement and its time derivative at all nodes and the
! rotations and its first derivatives at the master nodes of each element.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then
      Pos   (iNode,:)= Coords(iNode,:)
      PosDot(iNode,:)= 0.d0

      Psi   (iElem,k,:)= Psi0(iElem,k,:)
      PsiDot(iElem,k,:)= 0.d0

    ! Unconstrained nodes.
    else
      Pos   (iNode,:)= X   (ix+1:ix+3)
      PosDot(iNode,:)= dXdt(ix+1:ix+3)

      Psi   (iElem,k,:)= X   (ix+4:ix+6)
      PsiDot(iElem,k,:)= dXdt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

  !do i=1, size(Elem)
    !do j=1,Elem(i)%NumNodes
      !print*, Elem(i)%Master(j, :)
    !end do
  !end do

! Compute rotation vector and time derivative at slave nodes within elements.
  do i=1,size(Elem)
    do j=1,Elem(i)%NumNodes

! Copy rotation and derivative from master node for each slave node.
      if (Elem(i)%Master(j,1).ne.0) then
        Psi(i,j,:)   =Psi   (Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
        PsiDot(i,j,:)=PsiDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do

! Include master-to-slave initial rotations from the undeformed configuration.
    call cbeam3_projm2s (Elem(i)%NumNodes,Elem(i)%Master,Psi0(i,:,:),Psi0,Psi(i,:,:))
  end do

  return
 end subroutine cbeam3_solv_state2disp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_STATE2ACCEL
!
!-> Description:
!
!    Extract current acceleration from the state vector.
!
!-> Remarks:
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_state2accel (Elem,Node,dXddt,PosDDot,PsiDDot)
  use lib_fem
  use lib_rotvect
  use lib_bgeom
  use lib_cbeam3

! I/O Variables.
  type(xbelem),intent(in) :: Elem    (:)       ! Element information.
  type(xbnode),intent(in) :: Node    (:)       ! Nodal information.
  real(8),intent(in)      :: dXddt    (:)       ! Time derivatives of X.
  real(8),intent(out)     :: PosDDot  (:,:)     ! Current nodal position.
  real(8),intent(out)     :: PsiDDot  (:,:,:)   ! Current rotation at element nodes.

! Local variables.
  integer :: i,j,k               ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Store current accelerations at all nodes and the
! rotational accel. at the master nodes of each element.
  ix=0
  do iNode=1,size(PosDDot,DIM=1)
    iElem=Node(iNode)%Master(1)
    k=Node(iNode)%Master(2)

    ! Constrained nodes.
    if (Node(iNode)%Vdof.eq.0) then

      PosDDot(iNode,:) = 0.d0
      PsiDDot(iElem,k,:) = 0.d0

    ! Unconstrained nodes.
    else
      PosDDot(iNode,:)=dXddt(ix+1:ix+3)
      PsiDDot(iElem,k,:) = dXddt(ix+4:ix+6)

      ix=ix+6
    end if
  end do

! Compute rotation vector and time derivative at slave nodes within elements.
  do i=1,size(Elem)
    do j=1,Elem(i)%NumNodes

! Copy rotational accel. from master node for each slave node.
      if (Elem(i)%Master(j,1).ne.0) then
        PsiDDot(i,j,:)=PsiDDot(Elem(i)%Master(j,1),Elem(i)%Master(j,2),:)
      end if
    end do
  end do

  return
 end subroutine cbeam3_solv_state2accel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_DISP2STATE
!
!-> Description:
!
!    Write current state vector from current displacements and rotations.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_disp2state (Node,Pos,Psi,PosDot,PsiDot,X,dXdt)
  use lib_fem

! I/O Variables.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: Pos       (:,:)     ! Current nodal position.
  real(8),intent(in)      :: Psi       (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(in)      :: PosDot    (:,:)     ! Current nodal position.
  real(8),intent(in)      :: PsiDot    (:,:,:)   ! Current rotation at element nodes.
  real(8),intent(out)     :: X         (:)       ! Current displacements/rotations.
  real(8),intent(out)     :: dXdt      (:)       ! Current time derivative of X.


! Local variables.
  integer :: k                   ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Loop in all nodes in the model.
  ix=0
  do iNode=1,size(Pos,DIM=1)
    iElem=Node(iNode)%Master(1)
    if (Node(iNode)%Vdof.ne.0) then
      k=Node(iNode)%Master(2)

! Current nodal displacements and derivatives.
      X   (ix+1:ix+3)= Pos(iNode,:)
      dXdt(ix+1:ix+3)= PosDot(iNode,:)

! Cartesian rotation vector at master nodes.
      X   (ix+4:ix+6)= Psi(iElem,k,:)
      dXdt(ix+4:ix+6)= PsiDot(iElem,k,:)
      ix=ix+6
    end if
  end do

  return
 end subroutine cbeam3_solv_disp2state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine CBEAM3_SOLV_ACCEL2STATE
!
!-> Description:
!
!    Extract 2nd time-derivative of state vector from accelerations.
!
!-> Remarks:
!-> Author: Rob Simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine cbeam3_solv_accel2state (Node,PosDDot,PsiDDot,dXddt)
  use lib_fem

! I/O Variables.
  type(xbnode),intent(in) :: Node      (:)       ! Nodal information.
  real(8),intent(in)      :: PosDDot    (:,:)     ! Current nodal acceleration.
  real(8),intent(in)      :: PsiDDot    (:,:,:)   ! Current rotational accel. at element nodes.
  real(8),intent(out)     :: dXddt      (:)       ! Current second time-derivative of X.

! Local variables.
  integer :: k                   ! Counters.
  integer :: ix                  ! Counter on the degrees of freedom.
  integer :: iElem               ! Counter on the elements.
  integer :: iNode               ! Counter on the nodes.

! Loop in all nodes in the model.
  ix=0
  do iNode=1,size(PosDDot,DIM=1)
    iElem=Node(iNode)%Master(1)
    if (Node(iNode)%Vdof.ne.0) then
      k=Node(iNode)%Master(2)

! Current nodal displacements and derivatives.
      dXddt(ix+1:ix+3)= PosDDot(iNode,:)

! Cartesian rotation vector at master nodes.
      dXddt(ix+4:ix+6)= PsiDDot(iElem,k,:)
      ix=ix+6
    end if
  end do

  return
 end subroutine cbeam3_solv_accel2state

 subroutine dprint(line, ifile)
    character(len=*), intent(IN)        :: line
    character(len=*), intent(IN)        :: ifile
    print*, 'Line: ', line
    print*, 'File: ', ifile
 end subroutine dprint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module cbeam3_solv
