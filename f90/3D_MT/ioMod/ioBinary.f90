!   Derived from file_io.f90, with some additions to make
!   this more like Modular2D IO, but at the same time to allow
!   previously used 3D files formats to be handled.
!   Modified by Egbert, Nov. 2007, following previous modifications
!    to codes originally written by Kush Tandon.
!
!   Comments following are from file_io:
!   CHANGES FROM KUSH VERION:
!    1) Input is input, does not ever compare to
!   expected array sizes etc.  (as it was the input routine can only be
!    used when you already know sizes; it's useful to be able to just read the
!    file;  any consistency checking is kept separate.
!    2)  comments cleaned up
!    3) nMode is fixed at 2 ... if we get more general, this will all
!      have to be changed anyway.
!    4)  order of arguments changed most places
!    5)  got rid of routine for closing files: one executable stmt, make no sense
!    6)  BC IO only supported for BC data type: these routines don't
!        need to know about type excite
!    7) eliminate skipping and rewinding for output routines: in general
!         this makes no sense for writing a sequential binary file!
!      NOTE: we might want to at least allow the option to just read the next
!       records, assuming that the file is positioned correctly (as it will
!        be in all applications we now have)
!    8) format of impedance file changed: all site locations are in a second header
!          record, all impedances for one frequency are in one sequential record

! *****************************************************************************
module ioBinary
  ! This module contains io routines for boundary nodes,
  ! electric field solutions, impedances, and solver diagnostics

  use griddef
  use math_constants
  use emsolve3d
  use dataspace
  use datafunc
  use ForwardSolver, only: solnVectorMTX
  use transmitters
  use receivers
  use datatypes

  implicit none

  public                   :: FileReadInit, FileWriteInit
  public                   :: ZfileReadInit, ZfileWriteInit
  public                   :: BCfileRead, BCfileWrite
  public                   :: EfileRead, EfileWrite
  public                   :: DfileWrite
  public                   :: ZfileRead, ZfileWrite
  public                   :: read_Z,write_Z

Contains

  ! ***************************************************************************
  ! open file, read header record; used for reading header of
  !     all output files except impedance ... this just reads, returns
  !     what is in header.  Check consistency elsewhere if desired.
  !    If you don't want to check consistency, or get the header info,
  !    just open the file directly, don't use this routine!
   subroutine FileReadInit(inFile,ioNum,fileGrid,filePer,fileMode,fileVersion,ios)

    implicit none
    !  file name and iouint
    character (len=80), intent(in)              :: inFile
    integer, intent(in)                 	:: ioNum

    !  file header contents
    character (len=20), intent(out)             :: fileVersion
    type (grid_t), intent(out)                :: fileGrid
    integer, intent(out)                        :: filePer
    integer, intent(out)                        :: fileMode
    integer, intent(out)                        :: ios

    open (unit=ioNum,file=inFile,status='old', form ='unformatted',iostat=ios)

    IF( ios/=0) THEN
       WRITE(0,*) 'Error opening file: ', inFile
    else
       ! read the header
       read(ioNum) fileVersion, filePer, fileMode, &
	fileGrid%nx, fileGrid%ny, fileGrid%nz, fileGrid%nzAir, &
	fileGrid%ox, fileGrid%oy, fileGrid%oz, fileGrid%rotdeg

       ! the basic grid information is also stored in the header
       ! (don't we need to allocate first???)
       read(ioNum) fileGrid%dx
       read(ioNum) fileGrid%dy
       read(ioNum) fileGrid%dz

    endif
    return
    end subroutine FileReadInit

  !**********************************************************************
  ! check consistency between grids, number of periods
  subroutine CheckgridConsist(runGrid, fileGrid, nPer, filePer, ier)
    type (grid_t), intent(in)         :: fileGrid, runGrid
    integer, intent(in)			:: nPer, filePer
    integer, intent(out)		:: ier
    logical				:: ok = .true.

    ier = 0
    ok = ok .and. fileGrid%nx .eq. runGrid%nx
    ok = ok .and. fileGrid%ny .eq. runGrid%ny
    ok = ok .and. fileGrid%nz .eq. runGrid%nz
    ok = ok .and. fileGrid%nzAir .eq. runGrid%nzAir
    if(.not. ok) then
       write(0,*) 'Mismatch in grid size '
       ier = -1
    endif
    ok = ok .and. fileGrid%ox .eq. runGrid%ox
    ok = ok .and. fileGrid%oy .eq. runGrid%oy
    ok = ok .and. fileGrid%oz .eq. runGrid%oz
    if(.not. ok) then
       write(0,*) 'Mismatch in grid origin'
       ier = -1
    endif

    if(fileGrid%rotdeg /= runGrid%rotdeg) then
       write(0,*) 'Mismatch in grid rotation'
       ier = -1
    endif

    if(nPer /= filePer) then
       write(0,*) 'Mismatch in number of periods'
       ier = -1
    endif

  end subroutine CheckGridConsist

  ! ***************************************************************************
  !  open MT impedance file and read header (no reason for consistency check!)
  subroutine ZfileReadInit(inFile, ioNum, nPer, nSites, sites, ifBzTF,ios)

    implicit none
    !  input file name and io unit number
    character (len=80), intent(in)		:: inFile
    integer, intent(in)				:: ioNum

    ! output info from header
    integer, intent(out)			:: nPer
    integer, intent(out)			:: nSites
    integer, intent(out)			:: ios
    logical, intent(out)			:: ifBzTF
    real(kind=prec), pointer, dimension(:,:) :: sites

    integer		:: k

    open (unit=ioNum,file=inFile,status='old', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in ZfileReadInit: ', inFile
    endif

    ! read the header  ... why no version here?????
    read(ioNum) nPer, nSites, ifBzTF

    !  read in site locations
    allocate(sites(3,nSites))
    read(ioNum) (sites(:,k),k=1,nSites)

  end subroutine ZfileReadInit

  ! ***************************************************************************
  ! write header for MT_fwd output files (except impedance file)
  subroutine FileWriteInit(version, outFile, ioNum, inGrid, nPer, nMode,ios)

    implicit none
    character (len=20), intent(in)		:: version
    ! the version of the software in use
    character (len=80), intent(in)		:: outFile
    ! the filename
    type(grid_t), intent(in)			:: inGrid
    integer, intent(in)				:: nPer
    integer, intent(in)				:: nMode
    integer, intent(in)				:: ioNum
    integer, intent(out)			:: ios

    open (unit=ioNum,file=outFile,STATUS='unknown', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in FileWriteInit: ', outFile
    else
       ! write the header (contains the basic information for the forward
       ! modeling). the header is 4 lines
       write(ioNum) version,nPer,nMode,inGrid%nx,inGrid%ny,inGrid%nz,inGrid%nzAir, &
       inGrid%ox, inGrid%oy, inGrid%oz, inGrid%rotdeg
       write(ioNum) inGrid%dx
       write(ioNum) inGrid%dy
       write(ioNum) inGrid%dz
    endif
  end subroutine FileWriteInit

  ! ***************************************************************************
  ! write header for MT impedance file ... why no version?????
  subroutine ZfileWriteInit(inFile, ioNum, nPer, nSites, sites, ifBzTF,ios)

    implicit none
    ! file name and io unit
    character (len=80), intent(in)		:: inFile
    integer, intent(in)				:: ioNum

    ! header contents
    integer, intent(in)				:: nPer,nSites
    logical, intent(in)				:: ifBzTF
    real(kind=prec), intent(inout)	:: sites(3,nSites)

    !  iostatus
    integer, intent(out)			:: ios

    ! local variables
    integer	:: k

    open (unit=ioNum,file=inFile,status='unknown', form ='unformatted',&
         iostat=ios)

    if( ios/=0) then
       write(0,*) 'Error opening file in ZfileWriteInit: ', inFile
    else
       ! write the header
       write(ioNum) nPer, nSites, ifBzTF
    endif

    write(ioNum) (sites(:,k),k=1,nSites)

  end subroutine ZfileWriteInit

  ! ***************************************************************************
  ! read BCs for one mode, frequency; open unit ioNum first ...
  subroutine BCfileRead(ioNum, ifreq, imode, fileOmega, inBC)

    implicit none
    integer, intent(in)				:: ioNum,ifreq,imode

    real (kind=prec), intent(out)			:: fileOmega
    type (cboundary), intent(inout)		:: inBC

    !  local variables
    integer					:: nRecSkip, iskip
!    integer					:: iskip

    !  following could be optional oututs
    integer				:: fileIfreq, fileMode
    character (len = 20)            	:: ModeName

    !  hard code number of modes for now
    integer		:: nMode = 2

    ! calculate number of records before the header for this frequency/mode
    nRecSkip = ((ifreq-1)*nMode+(imode-1))*13+4

    ! rewind file, and skip to header record
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read the frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq, fileMode, ModeName

    ! read BC data for one frequency/mode - 12 records
    read(ioNum) inBC%xYMax
    read(ioNum) inBC%zYMax
    read(ioNum) inBC%xYMin
    read(ioNum) inBC%zYMin
    read(ioNum) inBC%yXMax
    read(ioNum) inBC%zXMax
    read(ioNum) inBC%yXMin
    read(ioNum) inBC%zXMin
    read(ioNum) inBC%xZMin
    read(ioNum) inBC%yZMin
    read(ioNum) inBC%xZMax
    read(ioNum) inBC%yZMax

  end subroutine BCfileRead

  ! ***************************************************************************
  ! write BCs for one mode, frequency; open unit ioNum first ...
  !  Note: for a sequential binary file skipping around makes no sense in general
  !   Just write next record (only thing that works, as far as I know!)
  subroutine BCfileWrite(ioNum, Omega, iFreq, iMode, ModeName, outBC)

    implicit none
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ioNum
    integer, intent(in)				:: ifreq,imode
    character (len=20), intent(in)		:: ModeName
    type (cboundary), intent(inout)		:: outBC

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode, ModeName

    ! write the BC data - 12 records
    write(ioNum) outBC%xYMax
    write(ioNum) outBC%zYMax
    write(ioNum) outBC%xYMin
    write(ioNum) outBC%zYMin
    write(ioNum) outBC%yXMax
    write(ioNum) outBC%zXMax
    write(ioNum) outBC%yXMin
    write(ioNum) outBC%zXMin
    write(ioNum) outBC%xZMin
    write(ioNum) outBC%yZMin
    write(ioNum) outBC%xZMax
    write(ioNum) outBC%yZMax

  end subroutine BCfileWrite

  ! ***************************************************************************
  ! read electric field solution for one mode, frequency; open unit ioNum first ...
  subroutine EfileRead(ioNum, ifreq, imode, fileOmega, inE)

    implicit none
    integer, intent(in)				:: ioNum,ifreq,imode
    type (cvector), intent(inout)		:: inE
    real (kind=prec), intent(out)	:: fileOmega

    !  local variables
    integer					:: nRecSkip, iskip
!    integer					:: iskip

    !  following could be optional oututs
    integer				:: fileIfreq, fileMode
    character (len = 20)            	:: ModeName

    !  hard code number of modes for now
    integer		:: nMode = 2

    ! calculate number of records before the header for this frequency/mode
    nRecSkip = ((ifreq-1)*nMode+(imode-1))*4+4

    ! rewind the file, skip to header record
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq, fileMode, ModeName

    ! read electrical field data - 3 records
    read(ioNum) inE%x
    read(ioNum) inE%y
    read(ioNum) inE%z

  end subroutine EfileRead

  ! ***************************************************************************
  ! write electrical field solution for one frequency/mode
  !    ... again no sense to skipping, rewinding for sequential binary write
  subroutine EfileWrite(ioNum,Omega, iFreq, iMode, ModeName, outE)

    implicit none
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ioNum,iFreq,iMode
    character (len=20), intent(in)		:: ModeName
    type (cvector), intent(in)			:: outE

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode,ModeName

    ! write the electrical field data - 3 records
    write(ioNum) outE%x
    write(ioNum) outE%y
    write(ioNum) outE%z

  end subroutine EfileWrite

  ! ***************************************************************************
  ! read impedances for nSites for one frequency
  subroutine ZfileRead(ioNum,ifreq,ifBzTF,nSites,Z,fileOmega)

    implicit none
    integer, intent(in)				:: ioNum,iFreq
    logical, intent(in)				:: ifBzTF
    integer, intent(in)				:: nSites
    complex (kind=prec), intent(out)	:: Z(3,2,nSites)
    real (kind=prec), intent(out)	:: fileOmega

    integer					:: nRecSkip,fileIfreq
    integer					:: iskip,i,j,k

    ! rewind file, skip past records for other frequencies
    nRecSkip = ifreq*2
    rewind(ioNum)
    do iskip = 1,nRecSkip
       read(ioNum)
    enddo

    ! read the frequency header - 1 record
    read(ioNum) fileOmega, fileIfreq

    if (ifBzTF) then
        read(ioNum) Z
    else
        read(ioNum)(((Z(i,j,k),i=1,2),j=1,2),k=1,nSites)
    end if

  end subroutine ZfileRead

  ! ***************************************************************************
  ! writes impedances for nSites for one frequency
  subroutine ZfileWrite(ioNum,Omega,ifreq,nSites,Z,ifBzTF)
  !   NOTE: impedances are stored in real dataVec structure:
  !    first dimension is mode, second component (real/imag for
  !     two impedance components + BzTF  ... rows and
  !     columns are switched compared to usual impedance
  !     storage.  Output in the usual way


    implicit none
    integer, intent(in)				:: ioNum
    real (kind=prec), intent(in)	:: Omega
    integer, intent(in)				:: ifreq
    integer, intent(in)				:: nSites
    real (kind=prec), intent(in)	:: Z(2,6,nSites)
    logical, intent(in)				:: ifBzTF
    integer					:: i,j,k,nComp

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq

    if (ifBzTF) then
       nComp = 6
    else
       nComp = 4
    endif

    write(ioNum) (((Z(i,j,k),j=1,nComp),i=1,2),k=1,nSites)

  end subroutine ZfileWrite

  ! ***************************************************************************
  ! writes solution diagnostics from unit opened by FileWriteInit. The prefix D
  ! is for diagnostics.
  subroutine DfileWrite(ioNum,Omega, iFreq, iMode, ModeName,solverDiag)

    implicit none
    real (kind=prec), intent(in)		:: Omega
    integer, intent(in)				:: ifreq,imode
    character (len=20), intent(in)		:: ModeName
    integer, intent(in)				:: ioNum
    type(emsolve_diag), intent(in)		:: solverDiag
    integer					:: nRecSkip, nMode
    integer					:: iskip
    integer					:: i, j

    ! write the frequency header - 1 record
    write(ioNum) Omega, iFreq, iMode, ModeName

    ! write the solution diagnostics - 3 records
    write(ioNum) solverDiag%nIterTotal, solverDiag%nDivCor
    write(ioNum) (solverDiag%EMrelErr(i),i = 1,solverDiag%nIterTotal)
    write(ioNum) ((solverDiag%divJ(i,j), i = 1,2), j = 1,solverDiag%nDivCor)

    !  why no output of divcor convergence?

  end subroutine DfileWrite ! DfileWrite

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!   Addtional routines, to allow greater compatability with Modular2D
!    code structure
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!   Routines to read and write dataVectorMTX objects; these are generalized
!   versions of write_Z and read_Z from Modular2D  ...
!    BUT NOTE: we still have to work on a more generic way to deal with
!    IO of dataVec objects (issue is meta-data, which is stored in transmitter
!    and receiver dictionaries for use by the program, but which here is
!    kept with the data in the file
!
   !**********************************************************************
   subroutine write_Z(fid,cfile,nTx,periods,nSites,sites,allData)
   ! writes impedance file, including list of periods, siteLocations
   !   NOTE: this assumes that the arrays "sites" and "periods" are
     !    essentially identical to the receiver and transmitter dictionaries
     !    (in which case, why have both these arrays and the dicts?)
     !   Also, we just get ncomp from each dataVec, and infer dataType
     !    from this.  Not at all general with regard to dataVec objects.

      integer, intent(in)			:: fid
      character(*), intent(in)			:: cfile
      integer, intent(in)			:: nTx,nSites
      real(kind=prec),intent(in)	:: periods(nTx)
      real(kind=prec),intent(in)	:: sites(3,nSites)
      type(dataVectorMTX_t),intent(in)			:: allData
      real(kind=prec), dimension(:,:), pointer :: siteTemp

     ! local variables
      integer   :: ns,iTx,k,j,nComp

      open(unit=fid,file=cfile,form='unformatted',status='unknown')
      write(fid) allData%nTx
      ! loop over periods
      do iTx = 1,allData%nTx
         ns = allData%d(iTx)%nSite
         ncomp = allData%d(iTx)%nComp
         ! write period, number of sites for this period
         write(fid) periods(iTx),nComp,ns
         ! write site locations for this period
         allocate(siteTemp(3,ns))
         do k = 1,ns
            siteTemp(:,k) = sites(:,allData%d(iTx)%rx(k))
         enddo
         write(fid) siteTemp
         !  note that each data field in a dataVec (i.e., allData%d(iTx)%data)
         !   is a real array of size (nComp,ns)
         write(fid) allData%d(iTx)%data
         write(fid) allData%d(iTx)%err
         deallocate(siteTemp)
      enddo
      close(fid)
      end subroutine write_Z
     !**********************************************************************
     subroutine read_Z(fid,cfile,nTx,periods,nSites,sites,allData)
     ! reads in data file, returns list of periods, , siteLocations, and
     !   sets up data vector structure, including data and error bars
     !   Also returns a list of periods, and sites ... not very general!
      integer, intent(in)       			:: fid
      character(*), intent(in)  			:: cfile
      integer, intent(out)      			:: nTx,nSites
      real(kind=prec),dimension(:), pointer     :: periods
      real(kind=prec),dimension(:,:), pointer   :: sites
      type(dataVectorMTX_t), intent(inout)   			:: allData

     ! local variables
      integer   	:: nComp,ns,iTx,k,l,j,Ndata
      real(kind=prec), pointer, dimension(:,:) :: siteTemp,siteTempAll
      logical		:: newSite

      open(unit=fid,file=cfile,form='unformatted',status='old')
      read(fid) nTx
      allocate(periods(nTx))
      allocate(allData%d(nTx))
      allData%allocated = .true.
      allData%errorBar = .true.
      allData%nTx = nTx

     ! loop over dataVec instances
      Ndata = 0
      do iTx = 1,nTx
         ! read in number of sites for this dataVec
         !  nTx is number of dataVecs ... might not all be for different periods!
         !   really should clean up list of periods (effectively the
	 !    transmitter dictionary
         read(fid) periods(iTx),nComp,ns
         ! read in site locations
         allocate(siteTemp(3,ns))
         read(fid) siteTemp
         ! create dataVec object, read in data
         allData%d(iTx)%errorBar = .true.
         call create_dataBlock(nComp,ns,allData%d(iTx))
         Ndata  = Ndata + nComp*ns
         allData%d(iTx)%tx = iTx

         selectcase(nComp)
            case(8)
               allData%d(iTx)%datatype =  Full_Impedance
            case(4)
               allData%d(iTx)%datatype =  Off_Diagonal_Impedance
         endselect

         read(fid) allData%d(iTx)%data
         read(fid) allData%d(iTx)%err

         if(iTx .eq. 1) then
           ! allocate temporary storage for full sites list
           ! (this might not always work ... I am assuming that the
           !        max number of sites is ns from period one * nTx)
            allocate(siteTempAll(3,ns*nTx))
            nSites = ns
            do k = 1,ns
               siteTempAll(:,k) = siteTemp(:,k)
               allData%d(iTx)%rx(k) = k
            enddo
         else
            ! check to see if site locations are already in list
            !  if not add to list; in any event set reciever "pointer" rx
            do k = 1,ns
               newSite = .true.
               do l = 1,nSites
                  if((siteTemp(1,k).eq.siteTempAll(1,l)).and.  &
                        (siteTemp(2,k).eq.siteTempAll(2,l)).and. &
                        (siteTemp(3,k).eq.siteTempAll(3,l))) then
                     newSite = .false.
                     allData%d(iTx)%rx(k) = l
                     exit
                  endif
               enddo
               if(newSite) then
                  nSites = nSites+1
                  siteTempAll(:,nSites) = siteTemp(:,k)
                  allData%d(iTx)%rx(k) = nSites
               endif
            enddo
         endif
         deallocate(siteTemp)
      enddo

      ! copy list of unique sites into "sites" array
      allocate(sites(3,nSites))
      do k = 1,nSites
         sites(:,k) = siteTempAll(:,k)
      enddo
      close(fid)
      end subroutine read_Z

!******************************************************************
      subroutine write_solnVectorMTX(fid,cfile,eAll)

      !  open cfile on unit fid, writes out object of
      !   type cvector in standard format (readable by matlab
      !   routine readcvector.m), closes file
      !  NOT coded at present to specifically write out TE/TM
      !    solutions, periods, etc. (can get this infor from
      !    eAll%solns(j)%tx, but only with access to TXdict.

      integer, intent(in)               :: fid
      character(*), intent(in)          :: cfile
      type(solnVectorMTX_t), intent(in)               :: eAll

      !   local variables
      integer           :: j,k,nMode = 2, ios
      character (len=20) 		:: version = '',ModeNames(2)
      real (kind=prec)    :: omega

      ModeNames(1) = 'Ey'
      ModeNames(2) = 'Ex'

      call FileWriteInit(version,cfile,fid,eAll%solns(1)%grid &
               ,eAll%nTX, nMode,ios)

      do j = 1,eAll%nTx
         do k = 1,2
           omega = txDict(eAll%solns(j)%tx)%omega
           call EfileWrite(fid, omega, j,  &
             k, ModeNames(k), eAll%solns(j)%pol(k))
         enddo
      enddo
      close(fid)
      end subroutine write_solnVectorMTX

end module ioBinary
