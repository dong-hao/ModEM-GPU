! *****************************************************************************
module DataIO
  ! This module contains io routines for reading and writing the data vectors
  ! Version: 3D MT-CSEM-SFF-TIDE-GLOBAL

  use math_constants
  use file_units
  use utilities
  use dataspace
  use gridcalc
  use transmitters
  use receivers
  use datatypes

  implicit none

  private

  ! switch between data formats by leaving uncommented one of the options below
  interface read_dataVectorMTX
	MODULE PROCEDURE read_Z_list
  end interface

  interface write_dataVectorMTX
	MODULE PROCEDURE write_Z_list
  end interface

  interface deall_dataFileInfo
    MODULE PROCEDURE deall_fileInfo
  end interface

  public     :: read_dataVectorMTX, write_dataVectorMTX, deall_dataFileInfo

  type :: data_file_block

      ! this block of information constitutes user preferences about the data format;
      ! there is one entry per each transmitter type and data type... (iTxt,iDt)
      ! if there are multiple data blocks of the same transmitter & data types,
      ! the last value is used.
      character(200) :: info_in_file
      character(20)  :: sign_info_in_file
      integer        :: sign_in_file
      character(20)  :: units_in_file
      real           :: origin_in_file(2)
      real           :: geographic_orientation

     ! these lists contain the indices into the data vector for each data type;
     ! they make it possible to sort the data by receiver for output.
     ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
     ! these indices are typically allocated as we read the data file
     integer, pointer, dimension(:)   :: tx_index
     integer, pointer, dimension(:)   :: dt_index
     integer, pointer, dimension(:,:) :: rx_index

     ! some transmitter types and data types don't go together
     logical         :: defined

  end type data_file_block

  ! private dictionary of data block info dimension (nTxt,nDt)
  ! where nTxt = number of all possible transmitter types
  !       nDt  = number of all possible data types
  ! number of transmitter types comes from the DICT/txTypes module
  ! and defines the number of conceptually different types of sources
  type (data_file_block), pointer, save, private, dimension(:,:) :: fileInfo

  ! we are converting from an "old format" to a "new format"
  ! the only difference being that in the new format, there is
  ! an additional line in the head that indicates transmitter type.
  ! on output, use the same format as on input. AK 25 May 2018
  logical, save, private  :: old_data_file_format = .true.

Contains

!**********************************************************************
! Sorts out the data block header

  function DataBlockHeader(txType,dataType) result (header)

    integer, intent(in)         :: txType
    integer, intent(in)         :: dataType
    character(200)              :: header

    select case (dataType)

       case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error HxAzi HyAzi ExAzi EyAzi'

       case(Full_Interstation_TF)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Ref_Code Ref_Lat '// &
                   'Ref_Lon Ref_X(m) Ref_Y(m) Ref_Z(m) Component Real Imag Error '// &
                   'HxAzi HyAzi HxRefAzi HyRefAzi'

       case(Off_Diagonal_Rho_Phase,Phase_Tensor)
          header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Value Error HxAzi HyAzi ExAzi EyAzi'

       case(Ex_Field, Ey_Field, Bx_Field, By_Field, Bz_Field)
              if (txType == CSEM) then
                header = 'Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) '// &
                         'Code X(m) Y(m) Z(m) Component Real Imag Error'
              else if (txType == TIDE) then
                header = 'Tx_Name Tx_Period(s) Tx_Amplitude Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'
              else
                header = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error'
              end if
	      
       case(Exy_Ampli_Phase)
          header = 'Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(x) Z(m) Component Real Error Rx_Azi'
	  
       case(Exy_Field)
          header = 'Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(x) Z(m) Component Real Imag Error Rx_Azi'

    end select

  end function DataBlockHeader

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine init_fileInfo(nTxt,nDt,nTx,nRx)

    integer, intent(in) :: nTxt,nDt
    integer, intent(in), optional :: nTx,nRx
    integer     :: istat,iTxt,iDt

    allocate(fileInfo(nTxt,nDt),STAT=istat)


     do iTxt = 1,nTxt
       do iDt = 1,nDt
         fileInfo(iTxt,iDt)%defined = .false.
         if (present(nTx) .and. present(nRx)) then
           allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
           allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
         end if
       end do
     end do

  end subroutine init_fileInfo

  ! **************************************************************************
  ! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_fileInfo()

    integer     :: i,j, istat

    if (associated(fileInfo)) then

     do i = 1,size(fileInfo,1)
       do j = 1,size(fileInfo,2)
          if (associated(fileInfo(i,j)%tx_index)) then
             deallocate(fileInfo(i,j)%tx_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%dt_index)) then
             deallocate(fileInfo(i,j)%dt_index,STAT=istat)
          end if
          if (associated(fileInfo(i,j)%rx_index)) then
             deallocate(fileInfo(i,j)%rx_index,STAT=istat)
          end if
       end do
     end do

     deallocate(fileInfo,STAT=istat)

    end if

  end subroutine deall_fileInfo

!**********************************************************************
! writes data in the ASCII list data file; it is convenient to work
! with data ordered by site (NOT by frequency), so we are going into
! some pains here to write them out in this order ...

   subroutine write_Z_list(allData,cfile)

    character(*), intent(in)                  :: cfile
    type(dataVectorMTX_t), intent(in)         :: allData
    ! local variables
    integer                         :: nTx,nRx,nDt,ncomp
    integer                         :: countData
    real(8), allocatable            :: value(:) ! (ncomp)
    real(8), allocatable            :: error(:) ! (ncomp)
    logical, allocatable            :: exist(:) ! (ncomp)
    character(2)                    :: temp = '> '
    character(8)                    :: today
    character(50)                   :: siteid,ref_siteid,compid
    character(20)                   :: sitename
    character(1000)                 :: strtemp
    integer                         :: iTxt,iTx,iRx,iDt,icomp,i,j,k,istat,ios,nBlocks
    real(8)                         :: x(3),ref_x(3),Xx(3),Period,SI_factor,large
    real(8)                         :: lat,lon,ref_lat,ref_lon,rx_azimuth
    logical                         :: conjugate, isComplex
    type(orient_t)                  :: azimu ! 2022.09.28, Liu Zhongyin, add azimuth for MT
    
    !===========================================================================
    !======================================================== New Local Variable
    !===========================================================================
    character(8)                    :: Dipole
    character(40)              	    :: Txid=''
    real(8) 			    :: Moment, Azi, Dip, LatTx, LongTx, Tx(3)
    real(8)                         :: Omega, Amplitude

    open(unit=ioDat,recl=4096,file=cfile,form='formatted',status='unknown')

    ! For each data type in dictionary, if data of this type exists, write it out.
    WRITE_TX_TYPE: do iTxt = 1,5

      ! subset to extract the data of one transmitter type only & work with that;
      ! call subset_dataVectorMTX(allData,txType_allData,iTxt)

      ! ideally we only want to loop over relevant data types but this hasn't been
      ! set up: would need a separate dictionary of data types per transmitter type...
      WRITE_DATA_TYPE: do iDt = 1,size(typeDict)

        ! skip those data types that don't belong to our transmitter type at all
        ! or even those that have not been defined on input
        if (.not. fileInfo(iTxt,iDt)%defined) then
            cycle WRITE_DATA_TYPE
        end if
    
      nBlocks = countDataBlock(allData,iDt,iTxt)
      if (nBlocks == 0) then
	! no data for this data type; skip it - this shouldn't happen anymore
	! since the "defined" logical deals with this on input
        cycle WRITE_DATA_TYPE
      else
        ! count the number of transmitters and receivers
        nTx = 0
        nRx = 0
        do iTx = 1,size(txDict)
            if (fileInfo(iTxt,iDt)%tx_index(iTx) > 0) then
                nTx = nTx + 1
            end if
        end do
        do iRx = 1,size(rxDict)
            do iTx = 1,size(txDict)
                if (fileInfo(iTxt,iDt)%rx_index(iTx,iRx) > 0) then
                    nRx = nRx + 1
                    exit
                end if
            end do
        end do
      end if

      ! write the data type header
      call date_and_time(today)
      call compact(fileInfo(iTxt,iDt)%info_in_file)
      write(strtemp,*) adjustl(trim(fileInfo(iTxt,iDt)%info_in_file)),' [',trim(today(1:4)),'-',trim(today(5:6)),'-',trim(today(7:8)),']'
      write(ioDat,'(a2)',advance='no') '# '
      write(ioDat,*,iostat=ios) adjustl(trim(strtemp))
      write(ioDat,'(a2)',advance='no') '# '
      write(ioDat,*,iostat=ios) adjustl(trim(DataBlockHeader(iTxt,iDt)))

      ! the new format is critical for JOINT modeling and inversion; otherwise, can stick
      ! to the old format for backwards compatibility. Will always write in the same format
      ! as the input data file
      if (.not. old_data_file_format) then
            write(ioDat,'(a2)',advance='no') '+ '
            write(ioDat,*,iostat=ios) trim(tx_type_name(iTxt))
      end if

      ! write the remainder of data type header
      call compact(typeDict(iDt)%name)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(typeDict(iDt)%name)
      call compact(fileInfo(iTxt,iDt)%sign_info_in_file)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(fileInfo(iTxt,iDt)%sign_info_in_file)
      call compact(fileInfo(iTxt,iDt)%units_in_file)
      write(ioDat,'(a2)',advance='no') temp
      write(ioDat,*,iostat=ios) trim(fileInfo(iTxt,iDt)%units_in_file)
      write(ioDat,'(a2,f8.2)',iostat=ios) temp,fileInfo(iTxt,iDt)%geographic_orientation
      write(ioDat,'(a2,2f9.3)',iostat=ios) temp,fileInfo(iTxt,iDt)%origin_in_file(1),fileInfo(iTxt,iDt)%origin_in_file(2)
      write(ioDat,'(a2,2i6)',iostat=ios) temp,nTx,nRx

      if (fileInfo(iTxt,iDt)%sign_in_file == ISIGN) then
          conjugate = .false.
      else if (abs(fileInfo(iTxt,iDt)%sign_in_file) == 1) then
          conjugate = .true.
      end if
      SI_factor = ImpUnits(typeDict(iDt)%units,fileInfo(iTxt,iDt)%units_in_file)

      ncomp = typeDict(iDt)%nComp
      allocate(value(ncomp),error(ncomp),exist(ncomp),STAT=istat)
      isComplex = typeDict(iDt)%isComplex
      countData = 0

      ! write data in order that is consistent with all previous versions of ModEM
      do iRx = 1,size(rxDict)
        do iTx = 1,size(txDict)

            k = fileInfo(iTxt,iDt)%rx_index(iTx,iRx)
            i = fileInfo(iTxt,iDt)%dt_index(iTx)
            j = fileInfo(iTxt,iDt)%tx_index(iTx)
            if (k == 0) then
                cycle
            end if
            value = SI_factor * allData%d(j)%data(i)%value(:,k)
            if (allData%d(j)%data(i)%errorBar) then
                error = SI_factor * allData%d(j)%data(i)%error(:,k)
            else
                error = LARGE_REAL
            end if
            exist = allData%d(j)%data(i)%exist(:,k)

            if (iTxt == CSEM) then
                Dipole = txDict(iTx)%Dipole
                Moment = txDict(iTx)%moment
                Azi = txDict(iTx)%AzimuthTx
                Dip = txDict(iTx)%dipTx
                Tx = txDict(iTx)%xyzTx
                !Txid = txDict(iTx)%id
            end if

            if (iTxt == TIDE) then
                Omega = txDict(iTx)%omega
                Amplitude = txDict(iTx)%amplitude
                Txid = txDict(iTx)%id
            end if

            Period = txDict(iTx)%period
            siteid = rxDict(iRx)%id
            x = rxDict(iRx)%x

            ! 2022.10.04, Liu Zhongyin, assign azimuth to azimu
            azimu = allData%d(j)%data(i)%orient(k)

            select case (iDt)

                case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)

                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(es13.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        if (FindStr(gridCoords, CARTESIAN)>0) then
                            write(ioDat,'(a50,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                   
                        elseif (FindStr(gridCoords, SPHERICAL)>0) then
                            read(siteid,'(a20,2f15.3)',iostat=ios) sitename,Xx(1),Xx(2)
                            write(ioDat,'(a20,5f15.3)',iostat=ios,advance='no') trim(sitename),x(1),x(2),Xx(1),Xx(2),x(3)
                            
                        end if
                        

                        ! 2022.09.28, Liu Zhongyin, add azimu while writing
                        if (conjugate) then
                                write(ioDat,'(a8,3es15.6,4f9.3)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp), &
                                azimu%azimuth%Hx-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Ex-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Ey-fileInfo(iTxt,iDt)%geographic_orientation
                        else
                                write(ioDat,'(a8,3es15.6,4f9.3)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp), &
                                azimu%azimuth%Hx-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Ex-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Ey-fileInfo(iTxt,iDt)%geographic_orientation
                        end if
                        countData = countData + 1
                    end do

                case(Full_Interstation_TF)

                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        ref_siteid = rxDict(iRx)%id_ref
                        ref_x = rxDict(iRx)%r
                        write(ioDat,'(es13.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(ref_siteid),ref_x(:)
                        ! 2022.09.28, Liu Zhongyin, add azimu while writing
                        if (conjugate) then
                                write(ioDat,'(a8,3es15.6,4f9.3)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp), &
                                azimu%azimuth%Hx-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hx_ref-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy_ref-fileInfo(iTxt,iDt)%geographic_orientation
                        else
                                write(ioDat,'(a8,3es15.6,4f9.3)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp), &
                                azimu%azimuth%Hx-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hx_ref-fileInfo(iTxt,iDt)%geographic_orientation, &
                                azimu%azimuth%Hy_ref-fileInfo(iTxt,iDt)%geographic_orientation
                        end if
                        countData = countData + 1
                    end do

                case(Off_Diagonal_Rho_Phase,Phase_Tensor)

                    do icomp = 1,ncomp
                        if (.not. exist(icomp)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        ! For apparent resistivities only, log10 of the values was used internally in the program;
                        ! writing out the linear apparent resistivity
                        if (index(compid,'RHO')>0) then
                            value(icomp) = 10**value(icomp)
                            ! Avoid Inf for FWD calculation
                            if (error(icomp) .ge. LARGE_REAL) then
                                error(icomp) = LARGE_REAL
                            else
                                error(icomp) = 10**error(icomp)
                            endif
                        end if

                        ! For Phase only, now using radians but output degrees [LiuZhongyin 2017.05.27]
                        if (index(compid,'PHS')>0) then
                            if (conjugate) then
                                value(icomp) = value(icomp)*R2D
                            else
                                value(icomp) = -value(icomp)*R2D
                            endif
                        end if
			
                        write(ioDat,'(es12.6)',    iostat=ios,advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        ! 2022.09.07, Liu Zhongyin, add azimu while writing
                            write(ioDat,'(a8,3es15.6,4f9.3)',iostat=ios) trim(compid),value(icomp),error(icomp), &
                            azimu%azimuth%Hx-fileInfo(iTxt,iDt)%geographic_orientation, &
                            azimu%azimuth%Ex-fileInfo(iTxt,iDt)%geographic_orientation, &
                            azimu%azimuth%Hy-fileInfo(iTxt,iDt)%geographic_orientation, &
                            azimu%azimuth%Ey-fileInfo(iTxt,iDt)%geographic_orientation
                        countData = countData + 1
                    end do

                case(Ex_Field,Ey_Field,Bx_Field,By_Field,Bz_Field)
					
                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        if (iTxt == CSEM) then
                            call compact(Dipole)
                            write(ioDat,'(a8)',iostat=ios,advance='no') trim(Dipole)
                            write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                            write(ioDat,'(2es12.6)', iostat=ios, advance='no') Period
                            write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                            write(ioDat,'(2es12.6)', iostat=ios, advance='no') Moment
                            write(ioDat,'(2f9.3, 3f12.3)',iostat=ios,advance='no') Azi, Dip, Tx
                        else if (iTxt == TIDE) then
                            call compact(Txid)
                            write(ioDat,'(a15)',  iostat=ios,advance='no') trim(Txid)
                            write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                            write(ioDat,'(es13.6)',  iostat=ios,advance='no') Period
                            write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                            write(ioDat,'(es13.6)',  iostat=ios,advance='no') Amplitude
                            write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        else
                            write(ioDat,'(es13.6)',  iostat=ios,advance='no') Period
                        end if
                        write(ioDat,'(a40,3f15.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
!                        write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp), error(2*icomp)
!                        error(2*icomp) = sqrt(value(2*icomp-1)*value(2*icomp-1) + value(2*icomp)*value(2*icomp))/100.0;
                        if (conjugate) then
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp)
                        else
                            write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp)
                        end if
                        countData = countData + 1
                    end do

                case(Exy_Ampli_Phase)
		
                    rx_azimuth = rxDict(iRx)%Rx_Azi
                    do icomp = 1,ncomp
                        if (.not. exist(icomp)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(a8)',iostat=ios,advance='no') trim(Dipole)
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(2es12.6)', iostat=ios, advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(2es12.6)', iostat=ios, advance='no') Moment
                        write(ioDat,'(2f9.3, 3f12.3)',iostat=ios,advance='no') Azi, Dip, Tx
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a15,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        if (icomp .eq. 1) then
                            write(ioDat,'(a10,3es15.6)',iostat=ios) trim(compid),value(icomp),0.1,rx_azimuth
                        else
                            write(ioDat,'(a10,3es15.6)',iostat=ios) trim(compid),value(icomp),5.0,rx_azimuth
                        end if
							
                        countData = countData + 1
                    end do	
						
                case(Exy_Field)
		
                    rx_azimuth = rxDict(iRx)%Rx_Azi
                    do icomp = 1,ncomp/2
                        if (.not. exist(2*icomp-1)) then
                            cycle
                        end if
                        compid = typeDict(iDt)%id(icomp)
                        write(ioDat,'(a8)',iostat=ios,advance='no') trim(Dipole)
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(2es12.6)', iostat=ios, advance='no') Period
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(2es12.6)', iostat=ios, advance='no') Moment
                        write(ioDat,'(2f9.3, 3f12.3)',iostat=ios,advance='no') Azi, Dip, Tx
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
                        write(ioDat,'(a15,3f12.3)',iostat=ios,advance='no') trim(siteid),x(:)
                        write(ioDat, '(a1)', iostat=ios,advance='no') ' '
!                        write(ioDat,'(a8,3es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp), error(2*icomp)
!                        error(2*icomp) = sqrt(value(2*icomp-1)*value(2*icomp-1) + value(2*icomp)*value(2*icomp))/100.0;
                        if (conjugate) then
                            write(ioDat,'(a10,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),-value(2*icomp),error(2*icomp),rx_azimuth
                        else
                            write(ioDat,'(a10,4es15.6)',iostat=ios) trim(compid),value(2*icomp-1),value(2*icomp),error(2*icomp),rx_azimuth
                        end if
                        countData = countData + 1
                    end do	
			
            end select

        end do  ! transmitters
      end do  ! receivers

      if (output_level > 4) then
	write(0,*) 'Written ',countData,' data values of type ',tx_type_name(iTxt),': ',trim(typeDict(iDt)%name),' to file'
      end if
      deallocate(value, error, exist, STAT=istat)

     end do WRITE_DATA_TYPE ! data types

    end do WRITE_TX_TYPE

    close(ioDat)

   end subroutine write_Z_list


!**********************************************************************
! reads in the ASCII list data file, sets up all dictionaries
! and the allData structure, including data and error bars.
! logic here is quite complicated, but once written can be used
! to read any kind of data, by adding a new case statement.

   subroutine read_Z_list(allData,cfile)

    character(*), intent(in)               :: cfile
    type(dataVectorMTX_t), intent(inout)   :: allData
    ! local variables
    type(dataVectorMTX_t)           :: newData
    integer                         :: nTx,nRx,nDt,ncomp,iRx,iTx,icomp
    integer                         :: countData,countRx
    complex(8), allocatable         :: value(:,:,:) ! (nTx,nRx,ncomp)
    real(8), allocatable            :: error(:,:,:) ! (nTx,nRx,ncomp)
    logical, allocatable            :: exist(:,:,:) ! (nTx,nRx,ncomp)
    integer, allocatable            :: new_TxType(:) ! contains txType indices (nTx)
    integer, allocatable            :: new_Tx(:) ! contains txDict indices (nTx)
    integer, allocatable            :: new_Rx(:) ! contains rxDict indices (nRx)
    character(2)                    :: temp
    character(200)                  :: txTypeName,typeName,typeInfo,typeHeader
    character(50)                   :: siteid,ref_siteid,compid
    integer                         :: nTxt,iTxt,iDt,i,j,k,istat,ios
    character(40)                   :: code,ref_code
    real(8)                         :: x(3),ref_x(3), Period,SI_factor
    real(8)                         :: lat,lon,ref_lat,ref_lon,rx_azimuth
    real(8)                         :: Zreal, Zimag, Zerr
    logical                         :: conjugate, errorBar, isComplex

    ! 2022.09.28, Liu Zhongyin, add azimu variable
    real(kind=prec), allocatable    :: HxAzimuth(:,:), ExAzimuth(:,:), HxAzimuth_ref(:,:) !(nTx,nRx)
    real(kind=prec), allocatable    :: HyAzimuth(:,:), EyAzimuth(:,:), HyAzimuth_ref(:,:) !(nTx,nRx)
    real(kind=prec)                 :: Hxangle, Exangle, Hxangle_ref, Hyangle, Eyangle, Hyangle_ref

    integer                         :: ncount
    character(1000)                 :: tmpline

    !===========================================================================
    !======================================================== New Local Variable
    !===========================================================================
    character(8)                    :: Dipole
    character(40)              	    :: Txid=''
    real(8) 			    :: Moment, Azi, Dip, LatTx, LongTx, Tx(3)
    real(8)                         :: Omega, Amplitude
    type(transmitter_t)             :: aTx

    ! First, set up the data type dictionary, if it's not in existence yet
    call setup_typeDict()

    ! Save the user preferences
    nTxt = 5
    nDt = size(typeDict)
    call init_fileInfo(nTxt,nDt)

    ! Now, read the data file
    open(unit=ioDat,file=cfile,form='formatted',status='old')
      
    ! Read the data blocks for each data type
    READ_DATA_TYPE: do
      
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeInfo
    	read(ioDat,'(a2,a200)',iostat=ios) temp,typeHeader
    	read(ioDat,'(a2,a100)',iostat=ios) temp,typeName

        ! If transmitter name exists, it precedes the typeName
        if (temp(1:1) == '+') then
            txTypeName = typeName
            read(ioDat,'(a2,a100)',iostat=ios) temp,typeName
            old_data_file_format = .false.
        else
            txTypeName = 'MT'
        end if
        iTxt = tx_type_index(txTypeName)
    	if (ios /= 0) exit
    
    	! Read new data type
    	call compact(typeName)
    	iDt = ImpType(typeName)
    	ncomp = typeDict(iDt)%nComp
    	if (typeDict(iDt)%isComplex) then
        	ncomp = ncomp/2
    	end if

        call compact(typeInfo)
    	fileInfo(iTxt,iDt)%defined = .true.
    	fileInfo(iTxt,iDt)%info_in_file = typeInfo
    	
    	! Sort out the sign convention
    	read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%sign_info_in_file
    	if(index(fileInfo(iTxt,iDt)%sign_info_in_file,'-')>0) then
      		fileInfo(iTxt,iDt)%sign_in_file = - 1
    	else
      		fileInfo(iTxt,iDt)%sign_in_file = 1
    	end if
    	if (fileInfo(iTxt,iDt)%sign_in_file == ISIGN) then
      		conjugate = .false.
    	else
      		conjugate = .true.
    	end if

        read(ioDat,'(a2,a20)',iostat=ios) temp,fileInfo(iTxt,iDt)%units_in_file
        SI_factor = ImpUnits(fileInfo(iTxt,iDt)%units_in_file,typeDict(iDt)%units)

        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%geographic_orientation
        read(ioDat,*,iostat=ios) temp,fileInfo(iTxt,iDt)%origin_in_file(1),fileInfo(iTxt,iDt)%origin_in_file(2)
        read(ioDat,*,iostat=ios) temp,nTx,nRx
        !write(0,'(a6,i5,a18,i8,a24)') 'Found ',nTx,' transmitters and ',nRx,' receivers in data block'


        if (output_level > 3) then
            write(0,*) node_info,'Reading data type: ',trim(typeName)
            write(0,*) node_info,'Sign convention in file: ',trim(fileInfo(iTxt,iDt)%sign_info_in_file)
            write(0,*) node_info,'Units in file: ',trim(fileInfo(iTxt,iDt)%units_in_file)
            write(0,*) node_info,'Number of transmitters: ',nTx
            write(0,*) node_info,'Number of receivers: ',nRx
        end if


        ! Allocate temporary data arrays
        allocate(new_TxType(nTx),new_Tx(nTx),new_Rx(nRx),STAT=istat)
        allocate(value(nTx,nRx,ncomp),error(nTx,nRx,ncomp),exist(nTx,nRx,ncomp),STAT=istat)

        ! 2022.09.28, Liu Zhongyin, add azimu allocation
        allocate(HxAzimuth(nTx,nRx),stat=istat)
        allocate(HyAzimuth(nTx,nRx),stat=istat)
        allocate(ExAzimuth(nTx,nRx),stat=istat)
        allocate(EyAzimuth(nTx,nRx),stat=istat)
        allocate(HxAzimuth_ref(nTx,nRx),stat=istat)
        allocate(HyAzimuth_ref(nTx,nRx),stat=istat)

        new_TxType(:) = 0
        new_Tx(:) = 0
        new_Rx(:) = 0
        value(:,:,:) = dcmplx(0.0d0,0.0d0)
        error(:,:,:) = LARGE_REAL
        exist(:,:,:) = .FALSE.
        countData = 0

        ! 2022.09.28, Liu Zhongyin, add azimu initial
        HxAzimuth(:,:) = R_ZERO
        HyAzimuth(:,:) = R_ZERO
        ExAzimuth(:,:) = R_ZERO
        EyAzimuth(:,:) = R_ZERO
        HxAzimuth_ref(:,:) = R_ZERO
        HyAzimuth_ref(:,:) = R_ZERO

        READ_DATA_LINE: Do

            select case (iDt)
            case(Ex_Field,Ey_Field,Bx_Field,By_Field,Bz_Field)
                if (iTxt == CSEM) then
                    read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0 ) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='CSEM'
                    aTx%nPol=1
                    aTx%Dipole = Dipole
                    aTx%period = Period
                    aTx%omega = 2.0d0*PI/Period
                    aTx%xyzTx = Tx
                    aTx%azimuthTx = Azi
                    aTx%dipTx = Dip
                    atx%moment = Moment
                    !aTx%id = Txid

                else if (iTxt == TIDE) then
                    read(ioDat,*,iostat=ios) Txid, Period, Amplitude, code, lat, lon, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0 ) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='TIDE'
                    aTx%nPol=1
                    aTx%omega = 2.0d0*PI/Period
                    aTx%period = Period
                    aTx%amplitude = Amplitude
                    aTx%id = Txid
                    call compact(aTx%id)

                else
                    read(ioDat,*,iostat=ios) Period, code, lat, lon, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr
                    if (ios /= 0) then
                        backspace(ioDat)
                        exit
                    end if
                    ! Find component id for this value
                    icomp = ImpComp(compid,iDt)
                    aTx%Tx_type='MT'
                    aTx%nPol=2
                    aTx%Dipole =''
                    aTx%period = Period
                    aTx%omega = 2.0d0*PI/Period

                end if

                ! Now overwrite aTx%Tx_type with txTypeName... allows for general SFF computation
                aTx%Tx_type = tx_type_name(iTxt)

                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code
                iRx = update_rxDict(x,siteid)

            case(Exy_Ampli_Phase)
                read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zerr, rx_azimuth
                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if
            			
                aTx%Tx_type='CSEM'
                aTx%nPol=1
                aTx%Dipole = Dipole
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period            		
                aTx%xyzTx = Tx
                aTx%azimuthTx = Azi
                aTx%dipTx = Dip
                atx%moment = Moment
						
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a20,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid,rx_azimuth)
		
            case(Exy_Field)
                read(ioDat,*,iostat=ios) Dipole, Period, Moment, Azi, Dip, Tx(1), Tx(2), Tx(3), code, x(1), x(2), x(3), compid, Zreal, Zimag, Zerr, rx_azimuth
                if (ios /= 0) then
                    backspace(ioDat)
                    exit
                end if
            			
                aTx%Tx_type='CSEM'
                aTx%nPol=1
                aTx%Dipole = Dipole
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period            		
                aTx%xyzTx = Tx
                aTx%azimuthTx = Azi
                aTx%dipTx = Dip
                atx%moment = Moment
						
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a20,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid,rx_azimuth)
		
            case(Full_Impedance,Off_Diagonal_Impedance,Full_Vertical_Components)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(12)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(13)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(14)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Exangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(15)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Exangle,Eyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Eyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle_ref = Hxangle_ref + 90.0
                case default
                    exit
                end select
		
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                ! Now overwrite aTx%Tx_type with txTypeName... allows for general SFF computation
                aTx%Tx_type = tx_type_name(iTxt)
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                
                if (FindStr(gridCoords, CARTESIAN)>0) then
                    write(siteid,'(a20,2f9.3)') code,lat,lon
                   
                elseif (FindStr(gridCoords, SPHERICAL)>0) then
                write(siteid,'(a20,2f15.3)') code,x(1),x(2)
                    x(1) = lat
                    x(2) = lon
                end if
                iRx = update_rxDict(x,siteid)

            case(Full_Interstation_TF)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(17)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(18)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(19)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hxangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(20)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Hxangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(21)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3), &
                        ref_code,ref_lat,ref_lon,ref_x(1),ref_x(2),ref_x(3),compid,Zreal,Zimag,Zerr,Hxangle,Hyangle,Hxangle_ref,Hyangle_ref
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = 0.0
                    Hxangle_ref = Hxangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hyangle_ref + fileInfo(iTxt,iDt)%geographic_orientation
                case default
                    exit
                end select
		
                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
		! Note that rx_azimuth is NOT used for MT: instead, we're supporting the
		! possibility that all fields components have different azimuths (stored in allData)
                write(siteid,'(a22,2f9.3)') code,lat,lon
                write(ref_siteid,'(a22,2f9.3)') ref_code,ref_lat,ref_lon
		rx_azimuth = R_ZERO
                iRx = update_rxDict(x,siteid,rx_azimuth,ref_x,ref_siteid)


            case(Off_Diagonal_Rho_Phase,Phase_Tensor)
                read(ioDat,'(a)',iostat=ios) tmpline

                if ((ios /= 0) .or. (tmpline(1:1)=='#')) then
                    backspace(ioDat)
                    exit
                end if

                ! Liu Zhongyin, 2019.08.27, add new codes for reading data
                backspace(ioDat)
                call strcount(tmpline, ' ', ncount)
                select case (ncount)
                case(10)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr
                    Hxangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(11)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hxangle + 90.0
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(12)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Hxangle
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(13)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle,Exangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Exangle + 90.0
                    Hyangle_ref = Hxangle_ref + 90.0
                case(14)
                    read(ioDat,*,iostat=ios) Period,code,lat,lon,x(1),x(2),x(3),compid,Zreal,Zerr,Hxangle,Hyangle,Exangle,Eyangle
                    Hxangle = Hxangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Exangle = Exangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hxangle_ref = 0.0
                    Hyangle = Hyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Eyangle = Eyangle + fileInfo(iTxt,iDt)%geographic_orientation
                    Hyangle_ref = Hxangle_ref + 90.0
                case default
                    exit
                end select

                ! Find component id for this value
                icomp = ImpComp(compid,iDt)

                ! For apparent resistivities only, use log10 of the values
                if (index(compid,'RHO')>0) then
                    Zerr  = Zerr/Zreal/dlog(10.0d0) ! Propagation of error
                    Zreal = log10(Zreal)
                end if

            	! For Phase only, using radians but reading degrees [LiuZhongyin 2017.05.27]
            	if (index(compid,'PHS')>0) then
                    if (conjugate) then
                	Zreal = Zreal*D2R
                    else
                        Zreal = -Zreal*D2R
                    endif
                    Zerr  = Zerr*D2R
                end if

                ! Update the transmitter dictionary and the index (sets up if necessary)
                aTx%Tx_type='MT'
                aTx%nPol=2
                aTx%period = Period
                aTx%omega = 2.0d0*PI/Period
                iTx = update_txDict(aTx)

                ! Update the receiver dictionary and index (sets up if necessary)
                ! For now, make lat & lon part of site ID; could use directly in the future
                write(siteid,'(a22,2f9.3)') code,lat,lon
                iRx = update_rxDict(x,siteid)

            end select

            ! complete transmitter dictionary update
            do i = 1,nTx
                if ((new_Tx(i) == iTx) .or. (new_Tx(i) == 0)) then
                    exit
                end if
            end do
            new_Tx(i) = iTx
            new_TxType(i) = iTxt

            ! complete receiver dictionary update
            do j = 1,nRx
                if ((new_Rx(j) == iRx) .or. (new_Rx(j) == 0)) then
                    exit
                end if
            end do
            new_Rx(j) = iRx

            ! record the value for storage in the data vector
            if (typeDict(iDt)%isComplex) then
                if (conjugate) then
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,-Zimag)
                else
                    value(i,j,icomp) = SI_factor * dcmplx(Zreal,Zimag)
                end if
            else
                value(i,j,icomp) = SI_factor * Zreal
            end if
            error(i,j,icomp) = SI_factor * Zerr
            exist(i,j,icomp) = .TRUE.

            ! 2022.09.28, Liu Zhongyin, assign angle to azimu
            HxAzimuth(i,j) = Hxangle
            HyAzimuth(i,j) = Hyangle
            ExAzimuth(i,j) = Exangle
            EyAzimuth(i,j) = Eyangle
            HxAzimuth_ref(i,j) = Hxangle_ref
            HyAzimuth_ref(i,j) = Hyangle_ref

            countData = countData + 1

        end do READ_DATA_LINE

	write(0,*) 'Read ',countData,' data values of ',trim(tx_type_name(iTxt)),' type ',trim(typeDict(iDt)%name),' from file'
	call create_dataVectorMTX(nTx,newData)
	newData%allocated = .TRUE.
	errorBar = .TRUE.
        SAVE_DATA: do i = 1,nTx

	       ! Count how many receivers we really have for this transmitter
	       countRx = 0
	       do j = 1,nRx
	        if(count(exist(i,j,:))>0) then
	            countRx = countRx + 1
	        end if
	       end do

	       ! Create a data vector for this transmitter and data type
	       call create_dataVector(1,newData%d(i))
	       newData%d(i)%tx = new_Tx(i)
	       newData%d(i)%txType = new_TxType(i)
	       newData%d(i)%allocated = .TRUE.
	       call create_dataBlock(typeDict(iDt)%nComp,countRx,newData%d(i)%data(1),typeDict(iDt)%isComplex,errorBar)
	       k = 1
	       do j = 1,nRx
	           ! If no data for this receiver, skip it
	           if(count(exist(i,j,:))==0) then
	            cycle
	           end if
	           ! Otherwise, write all components to data vector
	           do icomp = 1,ncomp
	            if(typeDict(iDt)%isComplex) then
	               newData%d(i)%data(1)%value(2*icomp-1,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%value(2*icomp  ,k) = imag(value(i,j,icomp))
	               newData%d(i)%data(1)%error(2*icomp-1,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%error(2*icomp  ,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%exist(2*icomp-1,k) = exist(i,j,icomp)
	               newData%d(i)%data(1)%exist(2*icomp  ,k) = exist(i,j,icomp)

	               ! 2022.09.28, Liu Zhongyin, add azimuth
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx = HxAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy = HyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ex = ExAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ey = EyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx_ref = HxAzimuth_ref(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy_ref = HyAzimuth_ref(i,j)
	            else
	               newData%d(i)%data(1)%value(icomp,k) = real(value(i,j,icomp))
	               newData%d(i)%data(1)%error(icomp,k) = error(i,j,icomp)
	               newData%d(i)%data(1)%exist(icomp,k) = exist(i,j,icomp)

	               ! 2022.09.28, Liu Zhongyin, add azimuth
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx = HxAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy = HyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ex = ExAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Ey = EyAzimuth(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hx_ref = HxAzimuth_ref(i,j)
	               newData%d(i)%data(1)%orient(k)%azimuth%Hy_ref = HyAzimuth_ref(i,j)
	            end if
	           end do
	           newData%d(i)%data(1)%rx(k) = new_Rx(j)
	           k = k+1
	       end do
	       newData%d(i)%data(1)%dataType = iDt
	       newData%d(i)%data(1)%tx = new_Tx(i)
	       newData%d(i)%data(1)%txType = new_TxType(i)
	       newData%d(i)%data(1)%allocated = .TRUE.

        end do SAVE_DATA

	! Merge the new data into the main data vector
	call merge_dataVectorMTX(allData,newData,allData)
    
	! 2022.09.28, Liu Zhongyin, deallocate azimu
	deallocate(HxAzimuth,stat=istat)
	deallocate(HyAzimuth,stat=istat)
	deallocate(ExAzimuth,stat=istat)
	deallocate(EyAzimuth,stat=istat)
	deallocate(HxAzimuth_ref,stat=istat)
	deallocate(HyAzimuth_ref,stat=istat)

	deallocate(value,error,exist,STAT=istat)
	deallocate(new_TxType,new_Tx,new_Rx,STAT=istat)
	call deall_dataVectorMTX(newData)

    end do READ_DATA_TYPE

    close(ioDat)

    ! Finished reading the data: write an empty line to screen
    write(0,*)

    ! Finally, set up the index vectors in the data type dictionary - used for output
    nTxt = 5
    nTx = size(txDict)
    nRx = size(rxDict)
    do iTxt = 1,nTxt
    	do iDt = 1,nDt
		allocate(fileInfo(iTxt,iDt)%tx_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%dt_index(nTx),STAT=istat)
	        allocate(fileInfo(iTxt,iDt)%rx_index(nTx,nRx),STAT=istat)
	        call index_dataVectorMTX(allData,iTxt,iDt,fileInfo(iTxt,iDt)%tx_index,fileInfo(iTxt,iDt)%dt_index,fileInfo(iTxt,iDt)%rx_index)
	end do
    end do

   end subroutine read_Z_list

end module DataIO
