        subroutine anomaly_detection_expo_MPI(SOURCE_LIST,DIR_OUT)
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

	integer i,j,complete,iexpo,jexpo,kexpo,npp,u
	integer source,tag,ierr,status(mpi_status_size)

	real para(np_anmly)
        

	if (my_id.eq.0) call initialize(SOURCE_LIST)

	call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(N_CHIP,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(IMAGE_FILE,NMAX_EXPO*NMAX_CHIP*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_OUTPUT,NMAX_EXPO*strl,mpi_character
     .,0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

        do i=1,np_anmly
          para(i)=0.
        enddo
       
	iexpo=0
	if (my_id.eq.0) then
	  iexpo=1
	  jexpo=N_EXPO
	  filename=trim(DIR_OUT)//'/anomaly_expo.dat'
          open(unit=10,file=filename,status='replace')
          rewind 10           
c	  filename=trim(DIR_OUT)//'/signature_expo.dat'
c          open(unit=20,file=filename,status='replace')
c          rewind 20           
        endif
       
	complete=0
	do while (complete.eq.0)
	  if (my_id.ne.0) then
	    call MPI_SEND(para,np_anmly,mpi_float,0,0,MPI_COMM_WORLD,ierr)
	    call MPI_RECV(iexpo,1,mpi_int,0,0,MPI_COMM_WORLD,status,ierr)
	    if (iexpo.eq.0) then
	      complete=1
	    else
	      call collect_expo_info_expo(iexpo,para)
	    endif
	  else
	    call MPI_RECV(para,np_anmly,mpi_float,MPI_ANY_SOURCE
     .,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
	    tag=status(MPI_TAG)
	    source=status(MPI_SOURCE)
	    call MPI_SEND(iexpo,1,mpi_int,source,tag,MPI_COMM_WORLD,ierr)
	    if (para(1).gt.0.5) then
              jexpo=jexpo-1
              kexpo=int(para(1)+0.5)
c              do i=1,N_CHIP(kexpo)
              write(10,*) para(1),para(2)
     .,trim(IMAGE_FILE(int(para(1)+0.5),1))
c                if (para(i+1).gt.0.02) write(20,*) kexpo,i 
c              enddo   
	    endif
	    write(*,*) 'processing expo: ',source,iexpo,jexpo
	    if (iexpo.ne.0) iexpo=iexpo+1
	    if (iexpo.gt.N_EXPO) iexpo=0
	    if (jexpo.eq.0) complete=1
	  endif	  
	enddo            

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

        if (my_id.eq.0) then
          close(10)
c          close(20)
        endif
        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        subroutine collect_expo_info_expo(iexpo,para)
        implicit none
        include 'para.inc'

        integer iexpo
        real para(np_anmly)
        
        character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	character*(strl) PREFIX,filename
	real aa(npara),star(NMAX_CHIP*nstar_max,ns,ns)
        real star_para(NMAX_CHIP*nstar_max,npara)
        real stamp(ns,ns),ave_p,p,stamp2(ns,ns),norm1,norm2
	integer ichip,ierror,i,j,k,u,v,w,nn1,nn2,nstar

        integer len_sam
        parameter (len_sam=100)
       
        do i=1,np_anmly
          para(i)=0.
        enddo
        para(1)=iexpo
       
        call get_PREFIX(IMAGE_FILE(iexpo,1),PREFIX)
	PREFIX=trim(DIR_OUTPUT(iexpo))//'/stamps/'//trim(PREFIX)
  	filename=trim(PREFIX)//'_star_info_expo.dat'
	nstar=0
        
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
	  write(*,*) filename
	  pause 'Catalog file error!!'
        endif   
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,7)
 	  if (ierror.lt.0) cycle
	  nstar=nstar+1
	  do i=1,7
	    star_para(nstar,i)=aa(i)
	  enddo
        enddo
        close(10)

        nn1=ns*len_sam
   	nn2=ns*(int(nstar/len_sam)+1)
        filename=trim(PREFIX)//'_star_power_expo.fits'
	call read_stamps(nstar_max*NMAX_CHIP,1,nstar,ns,ns
     .,star,nn1,nn2,filename)

        if (nstar.ge.2) then
          ave_p=0.
          do i=1,nstar-1
             do w=i+1,nstar
                norm1=0.
                norm2=0. 
                do j=1,ns
                   do k=1,ns
                      stamp(j,k)=star(i,j,k)  
                      stamp2(j,k)=star(w,j,k)
                      norm1=norm1+stamp(j,k)
                      norm2=norm2+stamp2(j,k)
                   enddo   
                enddo
                norm1=1./norm1
                norm2=1./norm2
                do j=1,ns
                   do k=1,ns
                      stamp(j,k)=stamp(j,k)*norm1  
                      stamp2(j,k)=stamp2(j,k)*norm2
                   enddo   
                enddo
                call ana_chi2(ns,stamp,stamp2,p)
                ave_p=ave_p+p
             enddo
          enddo
          para(2)=ave_p*2./((nstar-1.)*nstar)
        else
          para(2)=99.
        endif

         
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        subroutine anomaly_detection_MPI(SOURCE_LIST,DIR_OUT)
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

	integer i,j,complete,iexpo,jexpo,kexpo,npp,u
	integer source,tag,ierr,status(mpi_status_size)

	real para(np_anmly)
        

	if (my_id.eq.0) call initialize(SOURCE_LIST)

	call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(N_CHIP,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(IMAGE_FILE,NMAX_EXPO*NMAX_CHIP*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_OUTPUT,NMAX_EXPO*strl,mpi_character
     .,0,MPI_COMM_WORLD,ierr)

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

        do i=1,np_anmly
          para(i)=0.
        enddo
       
	iexpo=0
	if (my_id.eq.0) then
	  iexpo=1
	  jexpo=N_EXPO
	  filename=trim(DIR_OUT)//'/anomaly.dat'
          open(unit=10,file=filename,status='replace')
          rewind 10           
c	  filename=trim(DIR_OUT)//'/signature.dat'
c          open(unit=20,file=filename,status='replace')
c          rewind 20           
        endif
       
	complete=0
	do while (complete.eq.0)
	  if (my_id.ne.0) then
	    call MPI_SEND(para,np_anmly,mpi_float,0,0,MPI_COMM_WORLD,ierr)
	    call MPI_RECV(iexpo,1,mpi_int,0,0,MPI_COMM_WORLD,status,ierr)
	    if (iexpo.eq.0) then
	      complete=1
	    else
	      call collect_expo_info(iexpo,para)
	    endif
	  else
	    call MPI_RECV(para,np_anmly,mpi_float,MPI_ANY_SOURCE
     .,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
	    tag=status(MPI_TAG)
	    source=status(MPI_SOURCE)
	    call MPI_SEND(iexpo,1,mpi_int,source,tag,MPI_COMM_WORLD,ierr)
	    if (para(1).gt.0.5) then
              jexpo=jexpo-1
              kexpo=int(para(1)+0.5)
              do i=1,N_CHIP(kexpo)
                write(10,*) para(1),para(i+1)
c                if (para(i+1).gt.0.02) write(20,*) kexpo,i 
              enddo   
	    endif
	    write(*,*) 'processing expo: ',source,iexpo,jexpo
	    if (iexpo.ne.0) iexpo=iexpo+1
	    if (iexpo.gt.N_EXPO) iexpo=0
	    if (jexpo.eq.0) complete=1
	  endif	  
	enddo            

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

        if (my_id.eq.0) then
          close(10)
c          close(20)
        endif
        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        subroutine collect_expo_info(iexpo,para)
        implicit none
        include 'para.inc'

        integer iexpo
        real para(np_anmly)
        
        character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	character*(strl) PREFIX,filename
	real aa(npara),star(nstar_max,ns,ns),star_para(nstar_max,npara)
        real stamp(ns,ns),ave_p,p,stamp2(ns,ns)
	integer ichip,ierror,i,j,k,u,v,w,nn1,nn2,nstar

       
        do i=1,np_anmly
          para(i)=0.
        enddo
        para(1)=iexpo
       
        do ichip=1,N_CHIP(iexpo)
          call get_PREFIX(IMAGE_FILE(iexpo,ichip),PREFIX)
	  PREFIX=trim(DIR_OUTPUT(iexpo))//'/stamps/'//trim(PREFIX)
  	  filename=trim(PREFIX)//'_star_info.dat'
	  nstar=0          
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
	    write(*,*) filename
	    pause 'Catalog file error!!'
          endif   
c   	  read(10,*) 'ig xp yp SNR' 
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (aa(i),i=1,4)
 	    if (ierror.lt.0) cycle
	    nstar=nstar+1
	    do i=1,4
	      star_para(nstar,i)=aa(i)
	    enddo
          enddo
          close(10)

	  nn1=ns*len_s
   	  nn2=ns*(int(nstar/len_s)+1)
          filename=trim(PREFIX)//'_star_power.fits'
	  call read_stamps(nstar_max,1,nstar,ns,ns
     .,star,nn1,nn2,filename)

          if (nstar.ge.2) then
            ave_p=0.
            do i=1,nstar-1
              do j=1,ns
                do k=1,ns
                  stamp(j,k)=star(i,j,k)  
                enddo   
              enddo
c              call get_star_property(ns,stamp,p)
c              ave_p=ave_p+p
              do w=i+1,nstar
                do u=1,ns
                  do v=1,ns
                     stamp2(u,v)=star(w,u,v)
                  enddo
                enddo
                call ana_chi2(ns,stamp,stamp2,p)
                ave_p=ave_p+p
              enddo
            enddo
            para(ichip+1)=ave_p*2./((nstar-1.)*nstar)
          else
            para(ichip+1)=-1.
          endif
        enddo

         
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
