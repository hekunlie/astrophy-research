	program main
        use mpi
        include 'para.inc'

!> The purpose of this program is to measure cosmic shear

        integer ierr, my_id, num_procs
        common /MPIpar/my_id, num_procs

        character*(strl) SOURCE_LIST,DIR_OUT

        !----- MPI initialization ---------
        call MPI_Init( ierr )
        call MPI_comm_rank(MPI_COMM_WORLD, my_id, ierr )
        call MPI_comm_size(MPI_COMM_WORLD, num_procs, ierr )
	! my_id = 0,1,...,num_procs-1
	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

	call getarg(1,SOURCE_LIST)

	call getarg(2,DIR_OUT)

c	call pipeline(SOURCE_LIST,DIR_OUT)

c	call astrometry_multi_expo_MPI(SOURCE_LIST,DIR_OUT)
	
c	call shear_field_distortion_test_MPI(SOURCE_LIST,DIR_OUT)
	
	call anomaly_detection_expo_MPI(SOURCE_LIST,DIR_OUT)

c        !----- MPI finalization ------------------
 50	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 	call MPI_Finalize (ierr)
  
	stop
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pipeline(SOURCE_LIST,DIR_OUT)
	use mpi
	implicit none
	include 'para.inc'
	
        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs

        character*(strl) SOURCE_LIST,DIR_OUT,filename
	character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

        real anomaly(NMAX_EXPO)
        common /anomaly_pass/ anomaly
	
	integer ierr,status(mpi_status_size)

	external process
	
	if (my_id.eq.0) call initialize(SOURCE_LIST)

	call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(N_CHIP,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
c	call MPI_Bcast(eye_check,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(IMAGE_FILE,NMAX_EXPO*NMAX_CHIP*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_OUTPUT,NMAX_EXPO*strl,mpi_character
     .,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(ASTROMETRY_CAT,NMAX_EXPO*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(SOURCE_CAT,NMAX_EXPO*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_PSF,NMAX_EXPO*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)

	if (mod(PROCESS_stage,17).eq.0) then
          if (my_id.eq.0) then
  	    filename=trim(DIR_OUT)//'/anomaly_expo.dat'
            call read_anomaly(filename)
          endif
	  call MPI_Bcast(anomaly,NMAX_EXPO,mpi_float,0
     .,MPI_COMM_WORLD,ierr)
	endif
	
	call mpi_distribute(N_EXPO,process,'main process')
  	      
	write(*,*) my_id,'All finished !!'

	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initialize(SOURCE_LIST)
	implicit none
	include 'para.inc'
	
        character*(strl) SOURCE_LIST,IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)

	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	integer eye_check(NMAX_EXPO)
	common /eye_pass/ eye_check

	character*(strl) imagename,dirname,astroname,sourcename
        integer ierror,iexpo,i,nchip,eye_c

c Initialize the parameters and data array:	
	N_EXPO=0
	do i=1,NMAX_EXPO
	  N_CHIP(i)=0
        enddo

c read in the contents of file "source_list": 	
        open(unit=10,file=SOURCE_LIST,status='old',iostat=ierror)
        rewind 10

        if (ierror.ne.0) then
	  pause 'SOURCE_LIST reading error!!'
        endif
	
        do while (ierror.eq.0)	   
          read(10,*,iostat=ierror) iexpo,nchip,eye_c
	  if (ierror.eq.0) then
            N_EXPO=N_EXPO+1
	    N_CHIP(N_EXPO)=nchip
	    eye_check(N_EXPO)=eye_c
	    read(10,'(A)',iostat=ierror) DIR_OUTPUT(N_EXPO)
	    read(10,'(A)',iostat=ierror) ASTROMETRY_CAT(N_EXPO)
            read(10,'(A)',iostat=ierror) SOURCE_CAT(N_EXPO)	
            read(10,'(A)',iostat=ierror) DIR_PSF(N_EXPO)	
	    do i=1,N_CHIP(N_EXPO)
              read(10,'(A)',iostat=ierror) IMAGE_FILE(N_EXPO,i)	
	    enddo
	  endif
	enddo
        close(10)
        write(*,*) 'Total number of EXPOSURE: ',N_EXPO


        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine process(iexpo)
	implicit none
        include 'para.inc'

	integer iexpo,nchip
	character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	character*(strl) ASTROMETRY_CAT(NMAX_EXPO),SOURCE_CAT(NMAX_EXPO)	
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	character*(strl) IMAGE_FILE_pass(NMAX_CHIP)
	character*(strl) DIR_OUTPUT_pass,DIR_PSF_pass
	character*(strl) ASTROMETRY_CAT_pass,SOURCE_CAT_pass	
	
	integer i,ichip

c copy the source file information:
	
	DIR_OUTPUT_pass=DIR_OUTPUT(iexpo)
	DIR_PSF_pass=DIR_PSF(iexpo)
	ASTROMETRY_CAT_pass=ASTROMETRY_CAT(iexpo)
	SOURCE_CAT_pass=SOURCE_CAT(iexpo)

	nchip=N_CHIP(iexpo)
	
	do i=1,nchip
          IMAGE_FILE_pass(i)=IMAGE_FILE(iexpo,i)
	enddo

	if (mod(PROCESS_stage,2).eq.0) then
	  do ichip=1,nchip	 
	    call chip_pre_process(IMAGE_FILE_pass(ichip)
     .,ASTROMETRY_CAT_pass,SOURCE_CAT_pass,DIR_OUTPUT_pass)
	  enddo
	endif
cc
c generate the solution for the astrometric calibration	
	if (mod(PROCESS_stage,3).eq.0)
     .call chip_process_astrometry(IMAGE_FILE_pass
     .,nchip,DIR_OUTPUT_pass)
cc
        if (mod(PROCESS_stage,5).eq.0) then	
	  do ichip=1,nchip	 
            call chip_process_source(
     .IMAGE_FILE_pass,ichip,SOURCE_CAT_pass,DIR_OUTPUT_pass)
	  enddo
	endif
cc
	if (mod(PROCESS_stage,7).eq.0) then 	
	  do ichip=1,nchip	 
	    call chip_process_Fourier_T(
     .IMAGE_FILE_pass(ichip),DIR_OUTPUT_pass)
	  enddo
	endif
c!
	if (mod(PROCESS_stage,11).eq.0.and.ext_PSF.ne.1)
     . call PSF_reconstruction(IMAGE_FILE_pass,nchip,DIR_OUTPUT_pass)	 
	
cc
	if (mod(PROCESS_stage,13).eq.0)
     . call gen_shear_cat(IMAGE_FILE_pass,nchip,DIR_OUTPUT_pass
     .,DIR_PSF_pass)	  

	if (mod(PROCESS_stage,17).eq.0.and.ext_cat.eq.1) 	  
     . call combine_expo_catalog(nchip,IMAGE_FILE_pass
     .,DIR_OUTPUT_pass,iexpo)
	
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	


