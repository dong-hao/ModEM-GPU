
Module Declaration_MPI
#ifdef MPI
! include 'mpif.h'
     use mpi
     implicit none

! Declaration of general MPI stuff
!********************************************************************
Integer        :: taskid,total_number_of_Proc,number_of_workers
Integer        :: MASTER, FROM_MASTER, FROM_WORKER,TAG,ierr,dest
Integer        :: STATUS(MPI_STATUS_SIZE)
parameter         (MASTER=0,FROM_MASTER=1,FROM_WORKER=2,Tag=1)
!********************************************************************
! additional parameters needed by two-layered parallelization
!********************************************************************
integer        :: comm_world, comm_leader, comm_local
integer        :: rank_world, rank_leader, rank_local
integer        :: size_world, size_leader, size_local
integer        :: info_world, info_leader, info_local
integer        :: group_world, group_leader
! this is used to store the name of proc/cpu for current rank and identify 
! ranks from different hosts, useful when grouping cpus according to 
! topology - note hostname_len cannot exceed 40 (hard coded here)
character*(40) :: hostname_MPI
integer        :: hostname_len, ngroup
! one master switch for parallel paradigm - currently not supported for 
! on-the-fly modification, 
! TODO: need to think about it - if we really need an on-the-fly config
! 0 : simple grouping, 1 proc per each fwd/trn task (default)
! 1 : topology grouping, 1 node per each fwd/trn task
! 2 : equal grouping, n procs per each fwd/trn task
! 3 : dynamic grouping, variable number of procs per each fwd/trn task, 
!     for load-balancing
integer        :: para_method = 0
!********************************************************************
! additional parameters needed by CUDA acceleration
!********************************************************************
integer        :: size_gpu = 0 
! change the cpus_per_gpu if you want to use more than one cpus to 
! feed one gpu, modify at your own risk! 
integer        :: cpus_per_gpu = 3 ! hard coded here 
integer        :: device_id = -1
! this is used to store the timer of each mpi sub-process
DOUBLE PRECISION    :: previous_time
integer, allocatable, dimension(:) :: prev_group_sizes

!********************************************************************
! Parameters required to create an MPI derived data types.
!********************************************************************
Integer        :: cvector_mpi_3D,gridDef3D_mpi,eAll_mpi,dvecMTX_mpi
Integer        :: dvec_mpi,grid_t_mpi,worker_job_task_mpi
Integer        :: modelParam_t_mpi_sing,userdef_control_MPI
Integer        :: modelParam_t_mpi,extent
integer        :: oldtypes(0:20), blockcounts(0:20),offsets(0:20)
integer        :: block_lengths(0:20)
integer        :: displacements(0:20)
integer        :: address(0:21)
integer        :: typelist(0:21)
!********************************************************************



! Parameters used in communication
!********************************************************************
Integer        :: answers_to_receive,received_answers,recv_loop
Integer        :: who, which_stn,which_per,which_dt,which_pol,orginal_nPol
Integer , pointer, dimension(:)  :: eAll_location
logical                          :: eAll_exist=.false.
real*8,   pointer, dimension(:)  :: model_para_vec
character, pointer, dimension(:) :: eAll_para_vec       !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: e_para_vec          !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: sigma_para_vec      !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: data_para_vec       !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: worker_job_package  !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: userdef_control_package !! needed for MPI_pack/MPI_unpack; counted in bytes

Integer                          :: Nbytes              !! used in all MPI_pack/MPI_unpack
!********************************************************************


Integer                          :: nPol_MPI

! Time measuring
!********************************************************************
DOUBLE PRECISION    :: starttime,endtime,time_used
DOUBLE PRECISION    :: starttime_total,endtime_total
!********************************************************************


! A drived data type to distribute Info. between Processors.
!********************************************************************
type :: define_worker_job
     SEQUENCE
     character*80  :: what_to_do='NOTHING'
     Integer       :: per_index,Stn_index,pol_index,data_type_index,data_type
     Integer       :: taskid
     logical       :: keep_E_soln=.false.
     logical       :: several_Tx=.false.
     logical       :: create_your_own_e0=.false.
     ! 2022.10.06, Liu Zhongyin, add iSite storing the site index in rx of dataBlock_t
     Integer       :: iSite
 end type define_worker_job
type(define_worker_job), save :: worker_job_task
!********************************************************************

#ifdef CUDA

interface

   ! kernelc_getDevNum --> call cuda c kernel code with iso_c_binding
   integer(c_int) function kernelc_getDevNum() & 
    &              bind (C, name="kernelc_getDevNum" )
     ! interface to get the number of GPU devices - 
     ! only useful if you have any!   
     use iso_c_binding
     implicit none
   end function kernelc_getDevNum
    
end interface

#endif

Contains

!##########################################################################
subroutine create_worker_job_task_place_holder

     implicit none
     integer index,Nbytes1,Nbytes2,Nbytes3

       CALL MPI_PACK_SIZE(80, MPI_CHARACTER, MPI_COMM_WORLD, Nbytes1,  ierr)
       ! 2019.05.10, Liu Zhongyin, replace 6 with 7, for the iSite
       ! CALL MPI_PACK_SIZE(6, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(7, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(3, MPI_LOGICAL, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=(Nbytes1+Nbytes2+Nbytes3)+1

         if(.not. associated(worker_job_package)) then
            allocate(worker_job_package(Nbytes))
         end if

end subroutine create_worker_job_task_place_holder
!*******************************************************************************

subroutine Pack_worker_job_task
implicit none
integer index

index=1

        call MPI_Pack(worker_job_task%what_to_do,80, MPI_CHARACTER, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%per_index ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%Stn_index ,1 ,	 	MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%pol_index ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%data_type_index ,1 , 	MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%data_type ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%taskid ,1 , 			MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

        call MPI_Pack(worker_job_task%keep_E_soln,1, 		MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%several_Tx,1, 		MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%create_your_own_e0,1, MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

        ! 2019.05.08, Liu Zhongyin, add iSite for rx in dataBlock_t
        call MPI_Pack(worker_job_task%iSite, 1,             MPI_INTEGER, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine Pack_worker_job_task

subroutine Unpack_worker_job_task
implicit none
integer index
index=1
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%what_to_do,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)

        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%per_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%Stn_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%pol_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%data_type_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%data_type ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%taskid ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%keep_E_soln,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%several_Tx,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%create_your_own_e0,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)

        ! 2019.05.08, Liu Zhongyin, add iSite for rx in dataBlock_t
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%iSite, 1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

end subroutine Unpack_worker_job_task

subroutine gather_runtime(comm_current,time_passed,time_buff)
! simple subroutine to get the runtime of each sub tasks to access the
! parallel efficiency 
! collective on comm_current
      implicit none
      double precision,intent(in)                :: time_passed
      integer,intent(in)                         :: comm_current
      double precision,intent(out),pointer,dimension(:)   :: time_buff
      integer                                    :: current_rank
      integer                                    :: current_size,root=0
      call MPI_COMM_RANK(comm_current,current_rank,ierr)
      call MPI_COMM_SIZE(comm_current,current_size,ierr)
      allocate(time_buff(current_size))
      call MPI_Gather(time_passed, 1, MPI_DOUBLE_PRECISION, time_buff, 1,     &
     &     MPI_DOUBLE_PRECISION, root,comm_current,ierr) 
      return
end subroutine gather_runtime

#endif

end module Declaration_MPI
