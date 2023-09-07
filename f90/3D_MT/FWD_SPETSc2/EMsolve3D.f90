!**********************************************************************
! driver modules for solving the forward EM problem, including setup and
! solver
! this module is modified to use sparse matrix

module EMsolve3D
  use sg_boundary! work between different data types
  ! (between boundary conditions and
  ! complex vectors)
  use sg_sparse_vector, only: add_scvector
  use modelOperator3D  ! Maxwell operator module for sp
  use vectranslate     ! translate back and forth between Cvec and vec
  use solver           ! generic solvers rewrite for sp
  use solnspace

  implicit none
  public        :: FWDSolve3D
  public        :: deallSolverDiag, deallEMsolveControl
  public        :: createSolverDiag, getEMsolveDiag, setEMsolveControl

  interface FWDSolve3D
     MODULE PROCEDURE FWDSolve3Ds
     MODULE PROCEDURE FWDSolve3Dp
  end interface
  private       :: SdivCorr

  type :: emsolve_control
    ! Values of solver control parameters, e.g., read in from file
    ! plus other information on how the solver is to be initialized, 
    ! called, etc.
    ! idea is that this is the public access version of this info, which is
    ! copied into private version for actual solver control
    integer                   ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
    real(kind = 8)            ::      tolEMfwd, tolEMadj, tolDivCor
    logical                   ::      E0fromFile
    logical                   ::      UseDefaults
    logical                   ::      read_E0_from_File=.false.
    character (len=80)        ::      E0fileName
    integer                   ::      ioE0
    character (len=80)        ::      AirLayersMethod
    integer                   ::      AirLayersNz
    real(kind = 8)            ::      AirLayersMaxHeight, AirLayersAlpha, AirLayersMinTopDz
    real(kind = 8), pointer, dimension(:)   :: AirLayersDz
    logical                   ::      AirLayersPresent=.false.
    character (len=10)        ::      solver_name="QMR"		
    character (len=50) , public      ::   get_1D_from="Geometric_mean"	
  end type emsolve_control

  type :: emsolve_diag
    ! Solver diagnostic arrays, computed during run of forward solver.
    ! idea is that this is the public access version of this info, which is
    ! copied from the private version in module em_solve where this info is
    ! initially stored
    logical                   :: diagOut
    character (len=80)        :: fn_diagn
    integer                   :: ioDiag
    integer                   :: nIterTotal, nDivCor
    real(kind = 8), pointer, dimension(:)      ::      EMrelErr
    real(kind = 8), pointer, dimension(:,:)    ::      divJ
    real(kind = 8), pointer, dimension(:,:)    ::      DivCorRelErr
  end type emsolve_diag

  ! Default solver control parameters
  ! number of iterations for each call to divergence correction:
  integer, parameter    ::              IterPerDivCorDef = 120
  ! maximum number of divergence correction calls allowed
  integer, parameter    ::              MaxDivCorDef = 7
  ! maximum number of PCG iterations for divergence correction
  integer, parameter    ::              MaxIterDivCorDef = 100
  ! misfit tolerance for convergence of EMsolve algorithm
  real(kind=prec), parameter       ::      tolEMDef = 1E-7
  ! misfit tolerance for convergence of divergence correction solver
  ! note this is deprecated in SP2/SPETSc2
  real(kind=prec), parameter       ::      tolDivCorDef = 1E-7
  !Solver name, by default we use QMR
  character (len=10)  		   ::   solver_name="QMR"
  character (len=50) , public      ::   get_1D_from="Geometric_mean"

  save

  type(timer_t), private :: timer

  ! Actual values of control parameters must be set before first use,
  !     by call to setEMsolveControl
  !  of em_solve; are saved between calls, private to this module
  integer,  private        ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
  integer,  private        ::      MaxIterTotal ! = MaxDivCor*IterPerDivCor
  real(kind=prec), private ::      tolEMfwd, tolEMadj, tolDivCor

  ! EMsolve diagnostics: these are computed during execution of em_solve
  !   can be retrieved by call to getEmsolveDiag
  integer, private   :: nIterTotal, nDivCor
  logical, private   :: failed
  ! nIterTotal keeps tally on number of iterations so far
  ! nDivCor keeps tally on number of divergence correction so far
  real(kind=prec), pointer, dimension(:), private	::	EMrelErr
  real(kind=prec), pointer, dimension(:,:), private	::	divJ
  real(kind=prec), pointer, dimension(:,:), private	::	DivCorRelErr

Contains

!**********************************************************************
! main subroutine to Solve the forward EM problem;
! modified to use the sparse matrix data structure defined in
! sp modelOperator3D module
! renamed to "3Ds(olo)" to distinguish from the petsc version "3Dp(ara)"
! 
! If bRHS%adj = 'TRN' solves transposed problem  A^T x = b
!
! below is Anna's comment copied from the MF equivalent subroutine
!
! Note [AK 2018-05-10]:
! Any physical source has already been pre-multiplied by
!       [- ISIGN i\omega\mu_0]
! to yield
!       [- ISIGN i\omega\mu_0 j]
! on input to this routine. Note that this also holds for the secondary field
! formulation, where
!       j = dsigma * e,
! as well as for the tidal forcing, where
!       j = sigma (v x B).
! However, we still want to pre-compute the source RHS outside of this
! routine, for "generality".
! Specifically, Jmult supplies an interior source on the RHS that is
! not physical and is not pre-multiplied by that factor (except in Pmult).
! So it's cleaner to pass on the complete interior forcing in bRHS.
! For divergence correction, we divide by [+ ISIGN i\omega\mu_0] to get
!       i[- Div(j)].
! The plus sign is needed because we're taking the divergence of
!       curl(curl(E)) + ISIGN i\omega\mu_0 sigma E = f - curl(curl(b))
! Terms 1 and 4 vanishes, leaving:
!       Div(sigma E) - Div(f)/(+ ISIGN i\omega\mu_0) =  0.
! For a physical source j, this is equivalent to Div(sigma E) + Div(j) = 0;
! but the divergence correction may be applied also for non-physical sources,
! such  as in Jmult ('FWD') and JmultT ('TRN').
  subroutine FWDSolve3Ds(bRHS,omega,eSol,device_id)

    ! redefine some of the interfaces (locally) for our convenience
    use sg_vector !, only: copy => copy_cvector, &
    use vectranslate     ! translate back and forth between Cvec and vec
    use solver
    use spoptools
    ! generic routines for vector operations on the edge/face nodes
    ! in a staggered grid
    ! in cvec copy, remember the order is copy(new, old) i.e new = old

    implicit none
    !  INPUTS:
    type (RHS_t), intent(in)      :: bRHS
    real(kind=prec), intent(in)   :: omega
    !dummy parameter for compatibility
    integer, intent(in),optional  :: device_id
    !  OUTPUTS:
    !  eSol must be allocated before calling this routine
    type (cvector), intent(inout) :: eSol

    ! LOCAL VARIABLES
    logical                     :: converged,trans
    integer                     :: iter, fid
    integer                     :: Ne,Nei,Nni,Nn,i
    complex(kind=prec)          :: iOmegaMuInv
    ! e(lectric field) s(ource) b(rhs) phi0(div(s))
    complex(kind=prec), pointer, dimension (:) :: e,s,b
    complex(kind=prec), allocatable, dimension (:) :: ei,si,phi0
    complex(kind=prec), allocatable, dimension (:) :: temp, stemp
    character(80)                                  :: cfile
    !  band-aid cvector ...
    type (cvector)              :: tvec
    !  diagnostic structure for Krylov Subspace Solvers(KSS)
    type (solverControl_t)      :: KSSiter
    ! some warnings 
    if (present(device_id)) then
        if (device_id.ge.0) then
            write(6,*) "WARNING: we do not have GPU support without PETSc"
            write(6,*) "(in this branch, at least)"
            write(6,*) "fall back to CPU calculations"
        end if
    end if
    !  initialize solver diagnostic variables
    nIterTotal = 0
    nDivCor = 0
    EMrelErr = R_ZERO
    divJ = R_ZERO
    DivCorRelErr = R_ZERO
    failed = .false.
    trans = (bRHS%adj .eq. TRN)
    iOmegaMuInv = C_ONE/cmplx(0.0,1.0d0*ISIGN*omega*MU_0,kind=prec)
    if (.not.eSol%allocated) then
       write(0,*) 'eSol in EMsolve not allocated yet'
       stop
    else
    !   determine the edge numbers of the mesh
    !   need to write a interface for these
    !   since these will be private after debugging
        Nei = size(EDGEi,1)
        Ne = size(EDGEb,1)+Nei
        if (output_level > 3) then
            write(*,'(a36,i8,a4,i8)') 'FWDsolve3D model grid #edges: Nei=', &
                Nei,' Ne=',Ne
        end if
    end if
    ! allocate/initialize local data structures
    ! cboundary is a quite complex type...
    ! *essentially it should be e(EDGEb)
    ! for now we don't have an interface to deal with cboundary
    ! so just use cvectors to deliver the value...
    call create_cvector(bRHS%grid, tvec, eSol%gridtype)
    allocate(e(Ne))
    allocate(ei(Nei))
    allocate(s(Ne))
    allocate(si(Nei))
    allocate(b(Nei))
    allocate(temp(Ne))
    allocate(stemp(Nei))
    ! at this point e should be all zeros if there's no initial guess
    call getCVector(eSol,e)
    if(bRHS%nonZero_Source) then ! source (TRN)
        !   this is for *all* nodes
        Nni = size(NODEi,1)
        Nn  = size(NODEb,1) + Nni
        if (output_level > 3) then
            write(*,'(a36,i8,a4,i8)') 'FWDsolve3D source grid #nodes:           Nni=',    Nni,' Nn=',Nn
        end if
    ! uncomment the following line to try divergence correction in CCGD
    !    allocate(phi0(Nn)) ! make sure you *WANT* to do this, first!
    endif
    ! Using boundary condition and sources from rHS data structure
    ! construct vector b (defined only on interior nodes) for rHS of
    ! reduced (interior nodes only) linear system of equations

    if(trans) then ! TRN, trans=.true.
       !  In this case boundary conditions do not enter into forcing
       !    for reduced (interior node) linear system; solution on
       !    boundary is determined after solving for interior nodes
       if (bRHS%nonZero_Source) then
          if (bRHS%sparse_Source) then
             ! sparse source
             ! not sure how to do it efficiently with normal array
             ! for now it is just a walkaround, probably not going to
             ! be used by most
              call add_scvector(C_ONE,bRHS%sSparse,tvec)
              call getVector(tvec,s) !s is of size nEdge (all edges)
          else
             ! normal source
              call getVector(bRHS%s,s)
          endif
       else
          ! doesn't need to tamper with this part for sparse matrix
          ! let it go with cvectors...
          call zero(eSol)
          write(0,*) 'Warning: no sources for adjoint problem'
          write(0,*) 'Solution is identically zero'
          if(bRHS%nonzero_BC) then
          ! just copy input BC into boundary nodes of solution and return
             Call setBC(bRHS%bc, eSol)
          endif ! otherwise the eSol should be all zeros
          return
       endif
       ! NOTE that here we DO NOT divide the source by volume weights before
       ! the Div as the divcorr operation in SP is using VDiv instead of Div
       si = s(EDGEi) ! taking only the interior edges
       ! note that Div is formed from inner edges to all nodes
       ! divide by iOmegaMu and Volume weight to get the source term j
       ! uncomment the following line, to do divergence correction
       ! call Div(si,phi0)
       si = si * iOmegaMuInv / Vedge(EDGEi)
       ! calculate the modification term V_E GD_II j
       call RMATxCVEC(GDii,si,stemp)
       ! now i\omega\mu_0 V_E j + V_E GD_II j
       b = s(EDGEi) + stemp * Vedge(EDGEi)
    else ! trans = .false.
       ! In the usual forward model case BC does enter into forcing
       ! First compute contribution of BC term to RHS of reduced interior
       ! node system of equations : - A_IB*b
       if (bRHS%nonzero_BC) then
          !   copy from rHS structure into zeroed complex edge vector
          !   note that bRHS%bc is a cboundary type
          Call setBC(bRHS%bc, tvec) ! setBC -> copy_bcvector
          !   get info form BC
          call getVector(tvec,s)
          !   but only the boundary parts
          e(EDGEb) = s(EDGEb)
          !   Then multiply by curl_curl operator (use Mult_Aib ...
          !   Note that Mult_Aib is already multiplies by volume weights
          !   required to symetrize the problem, so the result is V*A_IB*b)
          !   essentially b = A(i,b)*e(b)
          Call Mult_Aib(e(EDGEb), trans, b)
       endif
       ! Add internal sources if appropriate: 
       ! Note that these must be multiplied explictly by volume weights
       !     [V_E^{-1} C_II^T V_F] C_II e + i\omega\mu_0\sigma e
       !   = - i\omega\mu_0  j - V_E^{-1} G_IA \Lambda D_AI j
       !     - [V_E^{-1} C_II^T V_F] C_IB b
       ! here we multiply by V_E throughout to obtain a symmetric system:
       !      V_E A_II e = - i\omega\mu_0 V_E j - V_E GD_II j - V_E A_IB b
       ! where
       !      A_II = V_E^{-1} C_II^T V_F C_II + i\omega\mu_0 \sigma,
       ! while
       !      A_IB = V_E^{-1} C_II^T V_F C_IB,
       ! and
       !      GD_II = V_E{-1} G_IA \Lambda D_AI.
       if (bRHS%nonzero_Source) then
          if (bRHS%sparse_Source) then
             ! sparse source
             call zero(tvec)
             call add_scvector(C_ONE,bRHS%sSparse,tvec)
             call getVector(tvec,s)
          else
             ! normal source
             call getVector(bRHS%s, s)
          endif
          ! At this point, s = - ISIGN * i\omega\mu_0 j
          ! Now Div(s) - will later divide by i_omega_mu to get the general
          ! divergence correction (j)
          temp = s*Vedge
          ! uncomment the following line to do divergence correction 
          ! call Div(temp(EDGEi), phi0)
          ! now temp = - ISIGN * i\omega\mu_0 V_E j
          ! divide by iOmegaMu to get the source term (j)
          si = s(EDGEi) * iOmegaMuInv 
          ! calculate the modification term GD_II j
          call RMATxCVEC(GDii,si,stemp)
          ! i\omega\mu_0 V_E j + V_E GD_II j
          stemp = temp(EDGEi) + stemp * Vedge(EDGEi)
          ! now add the V_E A_IB b term
          if(bRHS%nonzero_BC) then
             b = stemp - b
          else
             b = stemp
          endif
       else ! there is no source
           b = -b
       endif
    endif
    ! uncomment the following 3 lines to do divergence correction 
    ! if (bRHS%nonzero_Source) then
    !     phi0 = phi0 * iOmegaMuInv ! 1/i_omega_mu
    ! endif
    ! Outer part of KSS loop ... alternates between Calls to KSS solver
    ! and Calls to divcor  ... this will be part of EMsolve
    !
    ! e = current best solution (only on interior edges)
    ! b = rHS
    ! resetting
    nIterTotal = 0
    nDivCor = 0
    call reset_time(timer)
    ! Initialize iteration control/diagnostic structure for KSS
    if (trans) then
       KSSiter%tol = tolEMadj
    else
      if (bRHS%nonzero_BC) then
        KSSiter%tol = tolEMfwd
      else
        KSSiter%tol = tolEMadj
      end if
    end if

    KSSiter%niter = 0
    KSSiter%maxIt = IterPerDivCor
    allocate(KSSiter%rerr(IterPerDivCor))
    KSSiter%rerr = 0.0
    converged = .false.
    failed = .false.
    !   just take the interior elements
    !   Note: e here can be used for some initial guess
    ei = e(EDGEi)
    loop: do while ((.not.converged).and.(.not.failed))
       if (trim(solver_name) .eq. 'QMR') then
           write(*,*) 'I am using QMR with initial relative error ',KSSiter%rerr(1)
           Call QMR(b, ei, KSSiter)
       elseif (trim(solver_name) .eq. 'TFQMR') then
           write(*,*) 'I am using TFQMR with initial relative error ',KSSiter%rerr(1)
           Call TFQMR(b, ei, KSSiter)
       elseif (trim(solver_name) .eq. 'BICG') then
           write(*,*) 'I am using BICG with initial relative error ',KSSiter%rerr(1)
           Call BICG(b, ei, KSSiter)
       else
           write(*,*) 'Unsupported Forward Solver Method'
       end if
       ! algorithm is converged when the relative error is less than tolerance
       ! (in which case KSSiter%niter will be less than KSSiter%maxIt)
       converged = KSSiter%niter .lt. KSSiter%maxIt
       ! there are two ways of failing:
       !    1) the specific KSS did not work or
       !    2) total number of divergence corrections exceeded
       failed = failed .or. KSSiter%failed
       ! update diagnostics output from KSS
       do iter = 1,KSSiter%niter
           EMrelErr(nIterTotal+iter) = KSSiter%rerr(iter)
       end do
       if (KSSiter%niter.eq.0) then ! in case a initial guess is good enough
           KSSiter%niter = 1
           EMrelErr(KSSiter%niter) = KSSiter%rerr(1)
       endif
       nIterTotal = nIterTotal + KSSiter%niter
       nDivCor = nDivCor+1
       if( nDivCor < MaxDivCor) then
          ! do divergence correction
          ! NOTE: this is obsolete in SP2/SPETSc2
       else
          ! max number of divergence corrections exceeded; convergence failed
          failed = .true.
       endif
    end do loop
    if (output_level > 2) then
       write (*,'(a12,a20,i8,g15.7)') node_info, 'finished solving:',     &
   &            nIterTotal, EMrelErr(nIterTotal)
       write (*,'(a12,a22,f12.6)') node_info, 'solving time (sec): ',     &
   &          elapsed_time(timer)
    end if
    e(EDGEi) = ei
    !  After solving symetrized system, need to do different things for
    !   transposed, standard cases
    if(trans) then ! trans = .true.
       !   compute solution on boundary nodes: first  A_IB^T eSol
       call Mult_Aib(ei ,trans, s)
       !   Multiply solution on interior nodes by volume weights
       !   but after filling the solution on boundary
       temp = Vedge*e
       e = temp
       ! then b - A_IB^T eSol, where b is input boundary values (if any)
       if(bRHS%nonzero_BC) then
           e(EDGEb) = e(EDGEb) - s(EDGEb)
       else
           e(EDGEb) = -s(EDGEb)
       endif
    else ! trans = .false.
       ! just copy input BC into boundary nodes of solution
       if(.not.bRHS%nonzero_BC) then
           e(EDGEb) = 0
       endif
    endif
    call setVector(e,eSol)
    ! deallocate local temporary arrays
    deallocate(b)
    deallocate(s)
    deallocate(e)
    deallocate(ei)
    deallocate(temp)
    deallocate(stemp)
    deallocate(si)
    Call deall(tvec)
    deallocate(KSSiter%rerr)

  end subroutine FWDSolve3Ds

!**********************************************************************
! main subroutine to Solve the forward EM problem with petsc
! modified to use the petsc VEC and MAT distributed structure
! the ultimate idea is to parallel-ly prepare everything locally
! before assembly and solve the system once for all
!
! for now this is just (another) layer of translation
! to feed the system matrix to KSP solver.
!
! needs work to convert the spModelOperator module into a new module
! using VEC and MAT
! 
! I have merged the normal SP2 version and the SPETSc2 version
! the petsc version should work just fine (I think)
! if you have encountered any problems, please let me know:
! - send the feedback to donghao@cugb.edu.cn
! 
! If bRHS%adj = 'TRN' solves transposed problem  A^T x = b
! below is Anna's comment copied from the MF equivalent subroutine
!
! Note [AK 2018-05-10]:
! Any physical source has already been pre-multiplied by
!       [- ISIGN i\omega\mu_0]
! to yield
!       [- ISIGN i\omega\mu_0 j]
! on input to this routine. Note that this also holds for the secondary field
! formulation, where
!       j = dsigma * e,
! as well as for the tidal forcing, where
!       j = sigma (v x B).
! However, we still want to pre-compute the source RHS outside of this
! routine, for "generality".
! Specifically, Jmult supplies an interior source on the RHS that is
! not physical and is not pre-multiplied by that factor (except in Pmult).
! So it's cleaner to pass on the complete interior forcing in bRHS.
! For divergence correction, we divide by [+ ISIGN i\omega\mu_0] to get
!       i[- Div(j)].
! The plus sign is needed because we're taking the divergence of
!       curl(curl(E)) + ISIGN i\omega\mu_0 sigma E = f - curl(curl(b))
! Terms 1 and 4 vanishes, leaving:
!       Div(sigma E) - Div(f)/(+ ISIGN i\omega\mu_0) =  0.
! For a physical source j, this is equivalent to Div(sigma E) + Div(j) = 0;
! but the divergence correction may be applied also for non-physical sources,
! such  as in Jmult ('FWD') and JmultT ('TRN').

  subroutine FWDSolve3Dp(bRHS,omega,eSol,device_id, comm_local)
!----------------------------------------------------------------------
     ! redefine some of the interfaces (locally) for our convenience
     use sg_vector !, only: copy => copy_cvector, &
     use vectranslate     ! translate back and forth between Cvec and vec
     use solver
     use spoptools
     ! generic routines for vector operations on the edge/face nodes
     ! in a staggered grid
     ! in cvec copy, remember the order is copy(new, old) i.e new = old
!----------------------------------------------------------------------
!                   dependencies related to Petsc
!----------------------------------------------------------------------
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscviewer.h>
     use petscksp
     implicit none
!----------------------------------------------------------------------
!                   Variables related to Petsc
!----------------------------------------------------------------------
     integer, dimension(:), allocatable :: idx
     integer, dimension(:), pointer     :: isizes,isubs
     integer, dimension(3)              :: iedges
     integer                            :: ni, nj, na
     integer                            :: block_size
     integer                            :: rank_local,size_local
     integer                            :: ierr
     integer                            :: Nsub,Ntotal,istart,iend,csize
     integer, allocatable, dimension(:) :: ilocal,jlocal
     complex(kind=prec), allocatable, dimension (:) :: vlocal
     real(kind=prec)                    :: ptol=1e-2
     real                               :: normu
     type(spMatCSR_Cmplx)               :: Aii, Alocal
     Vec                                :: xvec,bvec,uvec,vvec
     Vec                                :: sxvec
     VecScatter                         :: xscat
     Mat                                :: Amat
     KSP                                :: ksp_local
     KSP,allocatable,dimension(:)       :: ksp_sub
     KSPType                            :: stype
     PC                                 :: pc_local,pc_sub
     PCType                             :: ptype,psubtype
     PetscViewer                        :: pviewer
     !  INPUTS:
     type (RHS_t), intent(in)           :: bRHS
     real(kind=prec), intent(in)        :: omega
     integer, intent(in)                :: device_id
     integer, intent(in)                :: comm_local
     !  OUTPUTS:
     !  eSol must be allocated before calling this routine
     type (cvector), intent(inout)      :: eSol
!----------------------------------------------------------------------
     ! LOCAL VARIABLES
     logical                                        :: converged,trans
     integer                                        :: iter, fid
     integer                                        :: Ne,Nei,Neb
     integer                                        :: Nn,Nni,Nnb,i
     complex(kind=prec)                             :: iOmegaMuInv
     ! e(lectric field) s(ource) b(rhs) phi0(div(s))
     complex(kind=prec), pointer, dimension (:)     :: e,s,b
     complex(kind=prec), allocatable, dimension (:) :: ei,si,phi0,phii
     complex(kind=prec), allocatable, dimension (:) :: temp,stemp
     
     character(80)                                  :: cfile
     ! band-aid cvector ...
     type (cvector)                                 :: tmpvec
     ! diagnostic structure for Krylov Subspace Solvers(KSS)
     type (solverControl_t)                         :: KSSiter
     ! initialize solver diagnostic variables
#ifndef PETSC_USE_COMPLEX
     write(6,*) 'PETSC IS NOT CONFIGURED WITH COMPLEX NUMBERS'
     write(6,*) 'NONE OF THE COMPLEX OPERATIONS WILL PROCEED, ABORTING'
     STOP
#endif
     ! now initialize solver diagnostic variables
     nIterTotal = 0
     nDivCor = 0
     EMrelErr = R_ZERO
     divJ = R_ZERO
     DivCorRelErr = R_ZERO
     failed = .false.
     iOmegaMuInv = C_ONE/cmplx(0.0,1.0d0*ISIGN*omega*MU_0,kind=prec)
     call MPI_COMM_RANK(comm_local,rank_local,ierr)
     call MPI_COMM_SIZE(comm_local,size_local,ierr)

     if (output_level > 3) then
         if ((device_id.ge.0).and.(size_local.eq.1)) then 
             write(6,*) 'GPU ACC STARTS...'
         else
             write(6,*) 'FALL BACK TO CPU...'
         end if
     end if
     ! PETSc initialization 
     PETSC_COMM_WORLD = comm_local
     call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
     !------------------------------------------------------------------------
     !     check the input variables
     !------------------------------------------------------------------------
     if (rank_local.eq.0) then ! leader
         trans = (bRHS%adj .eq. TRN)
         if (.not.eSol%allocated) then
             write(0,*) 'eSol in EMsolve not allocated yet'
             stop
         end if 
         ! determine the edge numbers of the mesh
         ! need to write a interface for these
         ! since these will be private after debugging
         Nei = size(EDGEi,1)
         Neb = size(EDGEb,1)
         Ne = Nei+Neb
         if(bRHS%nonZero_Source) then ! source (TRN)
             !   this is for *all* nodes
             Nni = size(NODEi,1)
             Nnb = size(NODEb,1)
             Nn = Nni+Nnb
             if (output_level > 3) then
                 write(*,'(a36,i8,a4,i8)') 'FWDsolve3D source grid #nodes:           Nni=',    Nni,' Nn=',Nn
             end if
         ! uncomment the following line to try divergence correction in CCGD
         !    allocate(phi0(Nn)) ! make sure you *WANT* to do this, first!
         endif
         iedges(1) = eSol%Nx*(eSol%Ny-1)*(eSol%Nz-1)
         iedges(2) = eSol%Ny*(eSol%Nx-1)*(eSol%Nz-1)
         iedges(3) = eSol%Nz*(eSol%Nx-1)*(eSol%Ny-1)
     end if
     ! common part for leader and workers
     ! broadcast those parameters to all workers
     call MPI_BCAST(Nei,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Neb,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Ne ,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nni,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nnb,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nn ,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(iedges ,3, MPI_INTEGER,0, comm_local,ierr)
     ! calculate the sizes of each local row matrix
     ! as well as the size of block preconditioners
     call calc_dist(size_local,iedges,isizes,isubs)
     if (rank_local .eq. 0) then !leader
         ! leader does all the works to allocate/initialize 
         ! local data structures
         if (output_level > 3) then
             ! for debug
             write(6,*) 'system iedges =', iedges
             write(6,*) 'system isizes =', isizes
             write(6,*) 'system isubs =', isubs
         endif 
         ! cboundary is a quite complex type...
         ! *essentially it should be e(EDGEb)
         ! for now we don't have an interface to deal with cboundary
         ! so just use cvectors to deliver the value...
         call create_cvector(bRHS%grid, tmpvec, eSol%gridtype)
         allocate(e(Ne))
         allocate(ei(Nei))
         allocate(s(Ne))
         allocate(si(Nei))
         allocate(b(Nei))
         allocate(temp(Ne))
         allocate(stemp(Nei))
         ! at this point e should be all zeros if there's no initial guess
         call getCVector(eSol,e)
         ! Using boundary condition and sources from rHS data structure
         ! construct vector b (defined only on interior nodes) for rHS of
         ! reduced (interior nodes only) linear system of equations
         if(trans) then ! TRN, trans=.true.
             !  In this case boundary conditions do not enter into forcing
             !  for reduced (interior node) linear system; solution on
             !  boundary is determined after solving for interior nodes
             if (bRHS%nonZero_Source) then
                 if (bRHS%sparse_Source) then
                     ! sparse source
                     ! not sure how to do it efficiently with normal array
                     ! for now it is just a walkaround, probably not going to
                     ! be used by most
                     call add_scvector(C_ONE,bRHS%sSparse,tmpvec)
                     call getVector(tmpvec,s)
                 else
                     ! normal source
                     call getVector(bRHS%s,s)
                 endif
             else
                 ! doesn't need to tamper with this part for sparse matrix
                 ! let it go with cvectors...
                 call zero(eSol)
                 write(0,*) 'Warning: no sources for adjoint problem'
                 write(0,*) 'Solution is identically zero'
                 if (bRHS%nonzero_BC) then
                     ! just copy input BC into boundary nodes of solution 
                     ! and return
                     Call setBC(bRHS%bc, eSol)
                 endif ! otherwise the eSol should be all zeros
                 return
             endif
             ! NOTE that here we DO NOT divide the source by volume weights
             ! before the Div as the divcorr operation in SP is using VDiv
             ! instead of Div
             si = s(EDGEi) ! taking only the interior edges
             ! note that Div is formed from inner edges to all nodes
             ! uncomment the following line, to do divergence correction
             ! call Div(si,phi0)
             ! divide by iOmegaMu and Volume weight to get the source term j
             si = si * iOmegaMuInv / Vedge(EDGEi)
             ! calculate the modification term V_E GD_II j
             call RMATxCVEC(GDii,si,stemp)
             ! now i\omega\mu_0 V_E j + V_E GD_II j
             b = s(EDGEi) + stemp * Vedge(EDGEi)
         else ! trans = .false.
             ! In the usual forward model case BC does enter into forcing
             ! First compute contribution of BC term to RHS of reduced 
             ! interior node system of equations : - A_IB*b
             if (bRHS%nonzero_BC) then
                 !   copy from rHS structure into zeroed complex edge vector
                 !   note that bRHS%bc is a cboundary type
                 Call setBC(bRHS%bc, tmpvec) ! setBC -> copy_bcvector
                 !   get info form BC
                 call getVector(tmpvec,s)
                 !   but only the boundary parts
                 e(EDGEb) = s(EDGEb)
                 !   Then multiply by curl_curl operator (use Mult_Aib ...
                 !   Note that Mult_Aib already multiplies by volume weights
                 !   required to symetrize problem, so the result is V*A_IB*b)
                 !   essentially b = A(i,b)*e(b)
                 Call Mult_Aib(e(EDGEb), trans, b)
             endif
        ! Add internal sources if appropriate: 
        ! Note that these must be multiplied explictly by volume weights
        !     [V_E^{-1} C_II^T V_F] C_II e + i\omega\mu_0\sigma e
        !   = - i\omega\mu_0  j - V_E^{-1} G_IA \Lambda D_AI j
        !     - [V_E^{-1} C_II^T V_F] C_IB b
        ! here we multiply by V_E throughout to obtain a symmetric system:
        !      V_E A_II e = - i\omega\mu_0 V_E j - V_E GD_II j - V_E A_IB b
        ! where
        !      A_II = V_E^{-1} C_II^T V_F C_II + i\omega\mu_0 \sigma,
        ! while
        !      A_IB = V_E^{-1} C_II^T V_F C_IB,
        ! and
        !      GD_II = V_E{-1} G_IA \Lambda D_AI.
             if (bRHS%nonzero_Source) then
                 if (bRHS%sparse_Source) then
                     ! sparse source
                     call zero(tmpvec)
                     call add_scvector(C_ONE,bRHS%sSparse,tmpvec)
                     call getVector(tmpvec,s)
                 else
                     ! normal source
                     call getVector(bRHS%s, s)
                 endif
                 ! At this point, s = - ISIGN * i\omega\mu_0 j
                 ! Now Div(s) - will later divide by i_omega_mu to get the 
                 ! general divergence correction (j)
                 temp = s*Vedge
                 ! now temp = - ISIGN * i\omega\mu_0 V_E j
                 ! divide by iOmegaMu to get the source (j)
                 si = s(EDGEi) * iOmegaMuInv
                 ! calculate the modification term GD_II j
                 call RMATxCVEC(GDii,si,stemp)
                 ! i\omega\mu_0 V_E j + V_E GD_II j
                 stemp = temp(EDGEi) + stemp * Vedge(EDGEi)
                 ! now add the V_E A_IB b term
                 if(bRHS%nonzero_BC) then
                     b = stemp - b
                 else
                     b = stemp
                 endif
             else! there is no source
                 b = -b
             endif
         endif ! trans 
     endif ! if rank == 0
     !-------------------------------------------------------------------------
     !    now beginning PETSc translation session
     !-------------------------------------------------------------------------
     if (rank_local.eq.0) then ! leader
         ! firstly leader constructs the system matrix Aii
         ! this will not be necesary if we build a module to contruct 
         ! the petsc version of operators at some stage
         call CSR_R2Cdiag(AAii,VOmegaMuSig,Aii)
         call deall_spMatCSR(AAii) ! release the CCGD matrix
         ! now calculate the rows that should be stored locally...
         ! Aii (Nei x Nei)
         do i = 2,size_local !now send those to your fellow workers
             call splitCMAT(Aii,i-1,size_local,Alocal,isizes)
             ni = size(Alocal%row)
             nj = size(Alocal%col)
             na = nj
             ! send the info to workers
             call MPI_SEND(ni,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(nj,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(na,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             ! now send the local matrix (in CSR format) that will be 
             ! stored in that worker process
             call MPI_SEND(Alocal%row,ni, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Alocal%col,nj, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Alocal%val,nj, MPI_DOUBLE_COMPLEX,i-1,2,         &
    &             comm_local, ierr)
             call deall_spMatCSR(Alocal) ! and release the temp sp matrix
         end do
         ! now deal with the local rows in leader
         call splitCMAT(Aii,rank_local,size_local,Alocal,isizes)
         block_size=Alocal%nRow
         call MatCreateMPIAIJWithArrays(comm_local,block_size,block_size, &
    &         Nei, Nei, Alocal%row,Alocal%col,Alocal%val,Amat,ierr)
         call deall_spMatCSR(Aii) ! release the temp sparse matrix
     else !workers
         ! get the info from leader
         call MPI_RECV(ni,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(nj,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(na,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         ! now get the local matrix (in CSR format) 
         allocate(ilocal(ni))
         allocate(jlocal(nj))
         allocate(vlocal(na))
         call MPI_RECV(ilocal,ni, MPI_INTEGER,0,2, comm_local,       &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(jlocal,nj, MPI_INTEGER,0,2, comm_local,       & 
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(vlocal,nj, MPI_DOUBLE_COMPLEX,0,2, comm_local,&
    &        MPI_STATUS_IGNORE,ierr)
         block_size=ni-1
         call MatCreateMPIAIJWithArrays(comm_local,block_size,block_size, &
    &         Nei, Nei, ilocal ,jlocal ,vlocal ,Amat,ierr)
     end if
     ! common part for leader and workers
     if (device_id.ge.0) then
         ! set as cuSPARSE format
         call MatSetType(Amat,MATMPIAIJCUSPARSE,ierr)
     else
         ! fall back to CPUs
         call MatSetType(Amat,MATMPIAIJ,ierr)
     end if
     call MatSetFromOptions(Amat,ierr)
     call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
     call deall_spMatCSR(Alocal) ! release the temp sp matrix
     call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
     !-------------------------------------------------------------------------
     !    PETSc vectors
     !-------------------------------------------------------------------------
     ! these part should be thread safe
     ! call VecCreateMPI(comm_local,block_size,Nei,bvec,ierr)
     ! create from matrix size -
     call MatCreateVecs(Amat,bvec,PETSC_NULL_VEC,ierr)
     call VecSetFromOptions(bvec,ierr)
     call VecSet(bvec,C_ZERO,ierr)
     call VecDuplicate(bvec,xvec,ierr)
     call VecSet(xvec,C_ZERO,ierr)
     ! call VecDuplicate(bvec,uvec,ierr)
     ! call VecSet(uvec,C_ONE,ierr)
     ! call VecDuplicate(bvec,vvec,ierr)
     ! create a vector to harvest the final result...
     call VecScatterCreateToZero(xvec,xscat,sxvec,ierr)
     !------------------------------------------------------------------------!
     ! Outer part of KSS loop ... alternates between Calls to KSS solver
     ! and Calls to divcor  ... this will be part of EMsolve
     ! these part should be thread safe
     !
     ! e = current best solution (only on interior edges)
     ! b = rHS
     !
     ! at present we don't really have the option to skip
     ! the divergence correction.  Not sure how/if this should
     ! be done.
     ! resetting
     nIterTotal = 0
     nDivCor = 0
     call reset_time(timer)
     ! Initialize iteration control/diagnostic structure for KSS
     if (rank_local.eq.0) then
         if (trans) then
             KSSiter%tol = tolEMadj
         else
             if (bRHS%nonzero_BC) then
                 KSSiter%tol = tolEMfwd
             else
                 KSSiter%tol = tolEMadj
             end if
         end if
         call MPI_BCAST(KSSiter%tol,1, MPI_DOUBLE,0, comm_local,ierr)
     else
         call MPI_BCAST(KSSiter%tol,1, MPI_DOUBLE,0, comm_local,ierr)
     end if

     KSSiter%niter = 0
     KSSiter%maxIt = IterPerDivCor
     allocate(KSSiter%rerr(IterPerDivCor))
     KSSiter%rerr = 0.0
     converged = .false.
     failed = .false.
     ! just take the interior elements
     ! construct the x and b as petsc VEC
     if (rank_local.eq.0) then ! leader prepare the vectors
         ei = e(EDGEi)
         call VecSetValues(bvec, Nei, (/(i,i=0,Nei-1)/), b, INSERT_VALUES, &
     &           ierr)
         call VecSetValues(xvec, Nei, (/(i,i=0,Nei-1)/), ei,               &
     &           INSERT_VALUES,ierr)
     end if 
     call VecAssemblyBegin(xvec,ierr)
     call VecAssemblyBegin(bvec,ierr)
     call VecAssemblyEnd(xvec,ierr)        
     call VecAssemblyEnd(bvec,ierr)
     !------------------------------------------------------------------------
     ! idea to test: for non-zero source START with divergence
     !               correction
     ! if (bRHS%nonzero_Source) then
     !     nDivCor = 1
     !     Call SdivCorr(ei,phi0)
     ! endif
     !-------------------------------------------------------------------------
     !     initialize KSP solver
     !-------------------------------------------------------------------------
     call KSPCreate(comm_local,ksp_local,ierr)
     call KSPSetOperators(ksp_local,Amat,Amat,ierr)
     !-------------------------------------------------------------------------
     !     configure PC preconditioner
     !-------------------------------------------------------------------------
     stype = KSPBCGS
     ! stype = KSPPREONLY
     ! stype = KSPTFQMR
     ptype = PCBJACOBI
     psubtype = PCILU
     ! psubtype = PCSOR
     call KSPGetPC(ksp_local,pc_local,ierr)
     call PCSetType(pc_local,ptype,ierr)
     call KSPSetType(ksp_local,stype,ierr)
     Ntotal = size(isubs)
     csize = 0
     Nsub = 0
     istart = 1
     do i=1,Ntotal
         csize = csize + isubs(i)
         if (csize.le.sum(isizes(1:rank_local))) then ! go on
             istart = istart + 1
         elseif (csize.le.sum(isizes(1:rank_local+1))) then 
             ! count sub-blocks in the local block
             Nsub = Nsub + 1
         else 
             exit
         end if
     end do
     iend = istart+Nsub-1
     if (output_level > 3) then
     ! for debug
         write(6,*) 'number of sub blocks =', Nsub, rank_local
         write(6,*) 'from #', istart, 'to #', iend, rank_local
     end if
     call PCBJacobiSetTotalBlocks(pc_local,Ntotal,            &
    &     isubs,ierr)
     call PCBJacobiSetLocalBlocks(pc_local, Nsub,             &
    &    isubs(istart:iend),ierr)
     ! set up the block jacobi structure
     call KSPSetup(ksp_local,ierr)
     ! allocate sub ksps
     allocate(ksp_sub(Nsub))
     call PCBJacobiGetSubKSP(pc_local,Nsub,istart,        &
    &     ksp_sub,ierr)
     do i=1,Nsub
         call KSPGetPC(ksp_sub(i),pc_sub,ierr)
         !ILU preconditioner
         call PCSetType(pc_sub,psubtype,ierr)
         call PCFactorSetLevels(pc_sub,0,ierr) ! use ILU(0) here
         call KSPSetType(ksp_sub(i),KSPPREONLY,ierr)
         ! this following lines are for MUMPS (use with caution)
         ! call PCFactorSetMatSolverType(pc_sub,MATSOLVERMUMPS,ierr)
         ! call PCFactorSetUpMatSolverType(pc_sub,ierr)
         ! call KSPSetTolerances(ksp_sub(i),ptol,PETSC_DEFAULT_REAL, &
    !&     PETSC_DEFAULT_REAL,KSSiter%maxit,ierr)
     end do
     call KSPSetFromOptions(ksp_local,ierr)
     call KSPSetup(ksp_local,ierr)
     call KSPSetTolerances(ksp_local,KSSiter%tol,PETSC_DEFAULT_REAL, &
    &     PETSC_DEFAULT_REAL,KSSiter%maxit,ierr)
     ! this is not yet supported by BCGS
     ! call KSPSetNormType(ksp_local, KSP_NORM_UNPRECONDITIONED, ierr)
     call KSPSetResidualHistory(ksp_local,KSSiter%rerr,KSSiter%maxit,    &
    &                PETSC_TRUE,ierr)
     call KSPSetInitialGuessNonzero(ksp_local, PETSC_TRUE,ierr)
!------------------------------------------------------------------------------
!     iteration starts
!------------------------------------------------------------------------------
     loop: do while ((.not.converged).and.(.not.failed))
         call KSPsolve(ksp_local,bvec,xvec,ierr)
         call KSPGetIterationNumber(ksp_local,KSSiter%niter,ierr)
         ! algorithm is converged when the relative error is less than 
         ! tolerance
         ! (in which case KSSiter%niter will be less than KSSiter%maxIt)
         converged = KSSiter%niter .lt. KSSiter%maxIt
         ! there are two ways of failing: 
         !    1) the specific KSS did not work or
         !    2) total number of divergence corrections exceeded
         failed = failed .or. KSSiter%failed
         ! update diagnostics output from KSS
         do iter = 1,KSSiter%niter
             EMrelErr(nIterTotal+iter) = KSSiter%rerr(iter)
         end do
         if (KSSiter%niter.eq.0) then ! in case a initial guess is good enough
             KSSiter%niter = 1
             EMrelErr(KSSiter%niter) = KSSiter%rerr(1)
         endif
         nIterTotal = nIterTotal + KSSiter%niter
         nDivCor = nDivCor+1
         if (nDivCor < MaxDivCor) then
             if (rank_local.eq.0) then ! leader
                 if (output_level > 3) then
                     ! note that the relative residual is "preconditioned"
                     ! (i.e. after applying the ILU preconditioner)
                     ! which is different from what you are experiencing 
                     ! in the SP/SP2/MF verstion
                     write(6,*) 'iter: ',nIterTotal,' residual: ',       &
    &                     EMrelErr(nIterTotal)
                 endif
             endif
         else
             ! max number of divergence corrections exceeded; 
             ! convergence failed
             failed = .true.
         endif
     end do loop
     call VecScatterBegin(xscat, xvec, sxvec, INSERT_VALUES,             &
    &     SCATTER_FORWARD,ierr)
     call VecScatterEnd(xscat, xvec, sxvec, INSERT_VALUES,               &
    &     SCATTER_FORWARD,ierr)
     if (rank_local.eq.0) then ! leader
         if (output_level > 2) then
             write (*,'(a12,a20,i8,g15.7)') node_info,'finished solving:',& 
    &             nIterTotal, EMrelErr(nIterTotal)
             write (*,'(a12,a22,f12.6)') node_info,'solving time (sec): ',&
    &             elapsed_time(timer)
         end if
     end if
     
     if (rank_local.eq.0) then ! leader prepare the results...
         call VecGetValues(sxvec,Nei,(/(i,i=0,Nei-1)/),ei,ierr)
         e(EDGEi) = ei
         ! After solving symetrized system, need to do different things for
         ! transposed, standard cases
         if (trans) then ! trans = .true.
             ! compute solution on boundary nodes: first  A_IB^T eSol
             call Mult_Aib(ei ,trans, s)
             ! Multiply solution on interior nodes by volume weights
             ! but after filling the solution on boundary
             temp = Vedge*e
             e = temp
             ! then b - A_IB^T eSol, where b is input boundary values (if any)
             if(bRHS%nonzero_BC) then
                 e(EDGEb) = e(EDGEb) - s(EDGEb)
             else
                 e(EDGEb) = -s(EDGEb)
             endif
         else ! trans = .false.
             ! just copy input BC into boundary nodes of solution
             if(.not.bRHS%nonzero_BC) then
                 e(EDGEb) = 0
             endif
         endif
         call setVector(e,eSol)
     else !while others sit back and watch 
         ! write(6,*) 'hanging around...'
     end if
     call MPI_BARRIER(comm_local,ierr)
     if (rank_local .eq. 0) then ! Leader 
         ! deallocate local temporary arrays
         deallocate(b)
         deallocate(s)
         deallocate(si)
         deallocate(e)
         deallocate(ei)
         deallocate(temp)
         deallocate(stemp)
         call deall(tmpvec)
     else
         deallocate(ilocal)
         deallocate(jlocal)
         deallocate(vlocal)
     end if
     if (allocated(ksp_sub)) then
         deallocate(ksp_sub)
     endif
     deallocate(KSSiter%rerr)
     call VecScatterDestroy(xscat,ierr);
     call VecDestroy(sxvec,ierr)
     call VecDestroy(xvec,ierr)
     call VecDestroy(bvec,ierr)
     call KSPDestroy(ksp_local,ierr)
     call MatDestroy(Amat,ierr)
     call PetscFinalize(ierr)

  end subroutine FWDSolve3Dp

!**********************************************************************
! solver_divcorrr contains the subroutine that would solve the divergence
! correction. Solves the divergene correction using pre-conditioned
! conjuagte gradient
! NOTE: this is obsolete in SP2/SPETSc2
subroutine SdivCorr(inE,phi0)
  ! Purpose: driver routine to compute divergence correction for input/output
  ! electric field vector inE
  ! Optional argument phi0 is scaled divergence of source term
  ! to be subtracted from current divergence

  use modelOperator3D
  ! rename routines for linear algebra operations; change to apply
  ! PCG to different problem

  implicit none

  complex(kind=prec), intent(inout),        dimension (:) :: inE
  complex(kind=prec), intent(in), optional, dimension (:) :: phi0

  !  local variables
  type (solverControl_t)                         :: PCGiter
  complex(kind=prec), allocatable, dimension (:) :: phiSol, phiRHS, phiAll
  complex(kind=prec), allocatable, dimension (:) :: cE
  integer                                        :: Nni,Nn,Nei,fid
  logical                                        :: SourceTerm

  SourceTerm = present(phi0)

  ! initialize PCGiter (maximum iterations allowed per set of diveregence
  ! correction, error tolerence, and relative error book keeping)
  PCGiter%maxIt = MaxIterDivCor
  PCGiter%tol = tolDivCor
  allocate(PCGiter%rerr(PCGiter%maxIt))
  PCGiter%rerr = 0.0

  ! alocating phiSol, phiRHS, phiAll
  Nni = size(NODEi)
  Nn = size(NODEb)+Nni
  Nei = size(inE)
  allocate(phiSol(Nni))
  allocate(phiRHS(Nni))
  allocate(phiAll(Nn))
  allocate(cE(Nei))
  phiSol=C_ZERO ! some initializations
  phiAll=C_ZERO

  ! compute divergence of currents for input electric field
  ! phi = V*Div*sigma*e
  call DivC(inE,phiRHS)
  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     phiRHS = phiRHS - phi0(NODEi)
  endif
  ! compute current Divergence before correction (using dot product)
  ! note this is the phi with volumn weights - not really just divJ
  divJ(1,nDivCor) = sqrt(dot_product(phiRHS,phiRHS))
  if (divJ(1,nDivCor).eq.0.0)  then
      ! initial current divergence is already zero
      ! no point to do the correction
      divJ(2,nDivCor) = divJ(1,nDivCor)
      if (output_level > 2) then
          write (*,'(a12,a25,g15.7,a25)')  node_info,                   &
    &           'current divergence is: ', divJ(2,nDivCor),             &
    &           'no need to do correction'
      endif
      return
  endif
  ! PCG is a (sort of) generic pre-conditioned CG algorithm
  Call PCG(phiRHS,phiSol,PCGiter)
  DivCorRelErr(:,nDivCor) = PCGiter%rerr
  if (output_level > 2) then
     write (*,'(a12,a32,i5,g15.7)') node_info,                            &
        'finished divergence correction:', PCGiter%niter,                 &
         PCGiter%rerr(PCGiter%niter)
  end if
  ! compute gradient of phiSol (Divergence correction for inE)
  phiAll(NODEi) = phiSol
  call Grad(phiAll,cE)

  ! subtract Divergence correction from inE
  inE = inE - cE

  ! divergence of the corrected output electrical field
  Call DivC(inE,phiRHS)

  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     phiRHS = phiRHS - phi0(NODEi)
  endif

  ! as in WS code, compute the current Divergence after correction
  ! (using the dot product)
  divJ(2,nDivCor) = sqrt(dot_product(phiRHS,phiRHS))

  ! output level defined in basic file_units module
  ! write(6,*) divJ(1,nDivCor), divJ(2,nDivCor)
  if (output_level > 3) then
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents before correction: ', divJ(1, nDivCor)
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents  after correction: ', divJ(2, nDivCor)
  end if

  ! deallocate the temporary work arrays
  deallocate(phiRHS)
  deallocate(phiSol)
  deallocate(phiAll)
  deallocate(cE)
  deallocate(PCGiter%rerr)

end subroutine SdivCorr ! SdivCorr

  !**********************************************************************
  ! Bundles the Inits that are used for an EM problem. These Inits can be
  ! used separately as well.
  subroutine ModelOperatorSetUp()

    ! COPIED OVER FROM MATRIX FREE VERSION FOR CONSISTENCY. CURRENTLY NOT USED
    ! BUT NECESSARY FOR REUSE OF THE SAME FORWARD SOLVER MODULE

  end subroutine ModelOperatorSetUp

  !**********************************************************************
  ! Deallocate the model operators after an EM problem is finished
  subroutine ModelOperatorCleanUp()

    ! COPIED OVER FROM MATRIX FREE VERSION FOR CONSISTENCY. CURRENTLY NOT USED
    ! BUT NECESSARY FOR REUSE OF THE SAME FORWARD SOLVER MODULE


  end subroutine ModelOperatorCleanUp

  !**********************************************************************
  ! setEMsolveControl sets actual solver control parameters, using info
  !  in structure solverControl, and allocates diagnostic arrays
  subroutine setEMsolveControl(solverControl,tolEM)

     type (emsolve_control), intent(in) ::  solverControl
     real (8), intent(in), optional     ::  tolEM

     if(solverControl%UseDefaults) then
        IterPerDivCor = IterPerDivCorDef
        MaxDivCor = MaxDivCorDef
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = MaxIterDivCorDef
        tolEMfwd = tolEMDef
        tolEMadj = tolEMDef
        tolDivCor = tolDivCorDef
        solver_name="QMR"
        get_1D_from="Geometric_mean"
     else
        IterPerDivCor = solverControl%IterPerDivCor
        MaxDivCor = solverControl%MaxDivCor
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = solverControl%MaxIterDivCor
        tolEMfwd = solverControl%tolEMfwd
        tolEMadj = solverControl%tolEMadj
        tolDivCor = solverControl%tolDivCor
        solver_name=solverControl%solver_name
        get_1D_from=solverControl%get_1D_from
     endif

     if (present(tolEM)) then
        tolEMfwd = tolEM
        tolEMadj = tolEM
     endif

     !  first check to see if diagnostic arrays are allocated
     !     ... if so deallocate first
     if(associated(EMrelErr)) then
        deallocate(EMrelErr)
     endif
     if(associated(divJ)) then
        deallocate(divJ)
     endif
     if(associated(DivCorRelErr)) then
        deallocate(DivCorRelErr)
     endif
     !   then allocate all arrays
     allocate(EMrelErr(MaxIterTotal))
     allocate(divJ(2,MaxDivCor))
     allocate(DivCorRelErr(MaxIterDivCor,MaxDivCor))

  end subroutine setEMsolveControl

   ! ***************************************************************************
   ! * readEMsolveControl reads the EM solver configuration from file
   subroutine readEMsolveControl(solverControl,rFile,fileExists,tolEM)

    type(emsolve_control), intent(inout)    :: solverControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    real(8), intent(in), optional           :: tolEM
    integer                                 :: ios
	logical                             :: exists
    character(80)                           :: string
    integer                                 :: istat

    ! Initialize inverse solver configuration

    inquire(FILE=rFile,EXIST=exists)
    if (present(fileExists)) then
       fileExists = exists
    end if

    if (.not. exists) then
       solverControl%UseDefaults = .true.
       if (present(tolEM)) then
          call setEMsolveControl(solverControl,tolEM)
       else
          call setEMsolveControl(solverControl)
       end if
       return
    else
       solverControl%UseDefaults = .false.
       write(*,*) node_info,'Reading EM solver configuration from file ',trim(rFile)
    end if

    !open (unit=ioFwdCtrl,file=rFile,status='old',iostat=ios)
     open (unit=ioFwdCtrl,file=rFile,form='formatted',status='old',iostat=ios)
    if(ios/=0) then
       write(0,*) node_info,'Error opening file: ', rFile
    end if

    ! This is the list of options specified in the startup file

    read (ioFwdCtrl,'(a48,i5)') string,solverControl%IterPerDivCor
    if (output_level > 2) then
       write (*,*)
       write (*,'(a12,a48,i5)') node_info,string,solverControl%IterPerDivCor
    end if
    read (ioFwdCtrl,'(a48,i5)') string,solverControl%MaxDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,i5)') node_info,string,solverControl%MaxDivCor
    end if
    read (ioFwdCtrl,'(a48,i5)') string,solverControl%MaxIterDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,i5)') node_info,string,solverControl%MaxIterDivCor
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolEMfwd
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolEMfwd
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolEMadj
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolEMadj
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolDivCor
    end if



    ! Check if there is an additional line.
    ! if yes, it corresponds to the larger E field solution.

    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a80)',iostat=istat) solverControl%E0fileName
    if (istat .eq. 0 ) then
        inquire(FILE=solverControl%E0fileName,EXIST=exists)
        if (index(string,'#')>0) then
            ! This is a comment line
            solverControl%read_E0_from_File=.false.
        else if (.not.exists) then
            write(*,*) node_info,'Nested E-field solution file not found and will not be used. '
            solverControl%read_E0_from_File=.false.
        else
            if (output_level > 2) then
                write (*,'(a12,a48,a80)') node_info,string,adjustl(solverControl%E0fileName)
            end if
            solverControl%read_E0_from_File=.true.
            solverControl%ioE0=ioE
        end if

    else
        solverControl%read_E0_from_File=.false.
    end if

    ! Now keep on reading for the air layers info
    ! If the number of air layers conflicts with that from the model file, we
    ! update the grid to use the controls
    ! Method options are: mirror; fixed height; read from file
    ! For backwards compatibility, defaults to what was previously hardcoded

    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a80)',iostat=istat) solverControl%AirLayersMethod
    !if (output_level > 2) then
    !    write (*,'(a12,a48,a80)') node_info,string,adjustl(solverControl%AirLayersMethod)
    !end if
    if (istat .eq. 0 ) then
        solverControl%AirLayersPresent = .true.
        if (index(solverControl%AirLayersMethod,'mirror')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersNz,solverControl%AirLayersAlpha,solverControl%AirLayersMinTopDz
        else if (index(solverControl%AirLayersMethod,'fixed height')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersNz,solverControl%AirLayersMaxHeight
        else if (index(solverControl%AirLayersMethod,'read from file')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,'(i3)',advance='no',iostat=istat) solverControl%AirLayersNz
            allocate(solverControl%AirLayersDz(solverControl%AirLayersNz), STAT=istat)
            if (solverControl%AirLayersNz > 0) then
                read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersDz
            end if
        else
            solverControl%AirLayersPresent = .false.
            call warning('Unknown air layers method option in readEMsolveControl')
        end if
    else
        solverControl%AirLayersPresent = .false.
    end if

    if (solverControl%AirLayersNz <= 0) then
        write(*,*) node_info,'Problem reading the air layers. Resort to defaults '
        solverControl%AirLayersPresent = .false.
    end if


    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a10)',iostat=istat) solverControl%solver_name
    if (istat .ne. 0) then
       solverControl%solver_name = 'QMR' ! default
    elseif (output_level > 2) then
       write (*,'(a12,a48,a)') node_info,string,adjustl(solverControl%solver_name)
    end if

    ! For any secondary field calculation approach...
    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a50)',iostat=istat) solverControl%get_1D_from
    if (istat .ne. 0) then
       solverControl%get_1D_from = 'Geometric_mean' ! default
    elseif (output_level > 2) then
       write (*,'(a12,a48,a)') node_info,string,adjustl(solverControl%get_1D_from)
    end if

    close(ioFwdCtrl)

    call setEMsolveControl(solverControl)


   end subroutine readEMsolveControl

  !**********************************************************************
  !   deallEMsolveControl deallocate
  subroutine  deallEMsolveControl(solverControl)
     type(emsolve_control), intent(inout),optional    :: solverControl

     integer istat

     if (present(solverControl)) then
        deallocate(solverControl%AirLayersDz, STAT=istat)
     end if

     deallocate(EMrelErr, STAT=istat)
     deallocate(divJ, STAT=istat)
     deallocate(DivCorRelErr, STAT=istat)

  end subroutine deallEMsolveControl

  !**********************************************************************
  ! getEMsolveDiag retrieves solver diagnositics
  subroutine getEMsolveDiag(solverDiagnostics)

     type (emsolve_diag), intent(inout)    ::      solverDiagnostics

     solverDiagnostics%nIterTotal = nIterTotal
     solverDiagnostics%nDivCor = nDivCor
     solverDiagnostics%EMrelErr = EMrelErr
     solverDiagnostics%divJ = divJ
     solverDiagnostics%DivCorRelErr = DivCorRelErr

  end subroutine getEMsolveDiag

  !***************************************************************************
  ! * createSolverDiag initializes emsolve_diag structure
  subroutine createSolverDiag(solverParams,solverDiag)

    implicit none
    type (emsolve_control), intent(in)  :: solverParams
    type (emsolve_diag), intent(inout)  :: solverDiag
    integer                             :: maxIterTotal

    maxIterTotal = solverParams%MaxDivCor*solverParams%IterPerDivCor
    allocate(solverDiag%EMrelErr(maxIterTotal))
    allocate(solverDiag%divJ(2,solverParams%MaxDivCor))
    allocate(solverDiag%DivCorRelErr(solverParams%MaxIterDivCor,  &
        solverParams%MaxDivCor))

  end subroutine createSolverDiag

  !***************************************************************************
  ! * deallSolverDiag deallocates emsolve_diag structure
  subroutine deallSolverDiag(solverDiag)
    type (emsolve_diag), intent(inout)  :: solverDiag

    deallocate(solverDiag%EMrelErr)
    deallocate(solverDiag%divJ)
    deallocate(solverDiag%DivCorRelErr)

  end subroutine deallSolverDiag

     !*************************************************************************
     ! * calculate the parallel data structure distribution among the processes
     ! * basically this is trying to figure out how you are going to split 
     ! * a matrix into a number of row-matrices to distribute onto different 
     ! * processes
     ! *
     ! * the preconditioner (Diagnol Block ILU) will always divide the 
     ! * matrix into 3 blocks, according to the iedges(x,y,z) length
     ! * 
     ! * the basic idea is to assign these 3 blocks to different processors
     ! * 
     ! * a new feature in petsc is more than one processor can be assigned 
     ! * for one single block
     ! *
     ! * the actual partition of row-matrices should follow this rule
     ! * for example: 
     ! * 1 processor -> 3 blocks
     ! * block x p1 
     ! * block y p1
     ! * block z p1
     ! * 2 processors -> 3 blocks
     ! * block x p1 
     ! * block y p1,p2
     ! * block z p2
     ! * 3 processors -> 3 blocks
     ! * block x p1 
     ! * block y p2
     ! * block z p3
     ! * 4 processors -> 3 blocks
     ! * block x p1
     ! * block y p2,p3
     ! * block z p4
     ! * 5 processors -> 3 blocks
     ! * block x p1 
     ! * block y p2,p3
     ! * block z p4,p5
     ! * 6 processors -> 3 blocks
     ! * block x p1,p2 
     ! * block y p3,p4
     ! * block z p5,p6
     ! * 7 processors -> 3 blocks
     ! * block x p1,p2 
     ! * block y p3,p4,p5
     ! * block z p6,p7
     !......
     subroutine calc_dist(Nproc,iedges,isizes,isubs)
     implicit none
     integer,intent(in)                            :: Nproc
     integer,intent(in), dimension(3)              :: iedges
     integer,intent(out), pointer,dimension(:)     :: isizes
     integer,intent(out), pointer,dimension(:)     :: isubs
     ! local variables
     integer                                       :: Nsub=1
     integer                                       :: i,j,k
     allocate(isizes(Nproc))

     if (Nproc.eq. 1) then ! 1 proc, 3 sub blocks
         k = 0
         allocate(isubs(3*Nproc))
         Nsub = 1
         do i = 1,3  ! loop through x/y/z
             do j = 1,Nsub !loop through these sub blocks (if any)
                 isubs(k+j) = iedges(i)/Nsub
             end do
             isubs((k+Nsub)) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
             k = k + Nsub
         end do
         isizes = sum(iedges)
    ! else if (Nproc.eq.2) then ! 2 proc, 3+3 sub blocks
    !     Nsub = 2
    !     allocate(isubs(3*Nsub))
    !     do i= 1,3 ! loop for x,y,z
    !         do j = 1,Nsub !loop through these sub blocks (if any)
    !             isubs((i-1)*Nsub+j) = iedges(i)/Nsub
    !         end do
    !         isubs((i*Nsub)) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
    !     end do
    !     isizes(1) = sum(isubs(1:3))
    !     isizes(2) = sum(isubs(4:6))
     else if (Nproc.eq.2) then ! 2 proc, 1+2 sub blocks
         Nsub = 1
         allocate(isubs(3*Nsub))
         do i= 1,3
            do j = 1,Nsub !loop through these sub blocks (if any)
                  isubs((i-1)*Nsub+j) = iedges(i)/Nsub
              end do
            isubs((i*Nsub)) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
         end do
         isizes(1) = isubs(1)
         isizes(2) = sum(isubs(2:3))
     else if (MOD(Nproc,3).eq.0) then ! 3 6 9...
         ! each side of edges (x,y,z) has Nsub blocks
         ! each process has exactly 1 block
         allocate(isubs(Nproc))
         Nsub = Nproc/3
         do i = 1,3  ! loop through x/y/z
             do j = 1,Nsub-1 !loop through these sub blocks (if any)
                 isubs((i-1)*Nsub+j) = iedges(i)/Nsub
             end do
             isubs(i*Nsub) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
         end do
         isizes = isubs

     elseif (MOD(Nproc,3).eq.1) then ! 4 7 10...
         ! x edges have Nsub+1 processes, while y and z have Nsub blocks
         k = 0
         allocate(isubs(Nproc))
         do i = 1,3  ! loop through x/y/z
             if (i .eq. 1) then ! x
                 Nsub = Nproc/3+1 
             else
                 Nsub = Nproc/3 
             endif
             do j = 1,Nsub !loop through these sub blocks (if any)
                 isubs(k+j) = iedges(i)/Nsub
             end do
             isubs((k+Nsub)) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
             k = k + Nsub
         end do
         isizes = isubs

     elseif (MOD(Nproc,3).eq.2) then ! 5 8 11...
         ! each process has exactly one sub block
         ! x and y have Nsub+1 processes, while z has Nsub processes
         k = 0
         allocate(isubs(Nproc))
         do i = 1,3  ! loop through x/y/z
             if ((i .eq. 1).or.(i.eq.3)) then ! x
                 Nsub = Nproc/3+1 
             else
                 Nsub = Nproc/3 
             endif
             do j = 1,Nsub !loop through these sub blocks (if any)
                 isubs(k+j) = iedges(i)/Nsub
             end do
             isubs((k+Nsub)) = iedges(i)-(iedges(i)/Nsub)*(Nsub-1)
             k = k + Nsub
         end do
         isizes = isubs
     endif
     end subroutine calc_dist
end module EMsolve3D
