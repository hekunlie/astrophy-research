        subroutine shear_field_distortion_test_MPI(SOURCE_LIST,DIR_OUT)
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

c	integer eye_check(NMAX_EXPO)
c	common /eye_pass/ eye_check

        real anomaly(NMAX_EXPO)
        common /anomaly_pass/ anomaly

	integer i,j,complete,iexpo,jexpo,kexpo,npp,u
	integer source,tag,ierr,status(mpi_status_size)

	real x1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real x2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real y1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real y2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de(ngfieldmax*NMAX_EXPO/num_procs_min)
        real yy1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real yy2(ngfieldmax*NMAX_EXPO/num_procs_min)
        integer ng
        common /shear_data_pass/ x1,x2,y1,y2,de,de1,de2,yy1,yy2,ng        
        
	real xmin,xmax
	common /bound_pass/ xmin,xmax
        real yeqx
	external yeqx

        external read_shear_cat
 
	if (my_id.eq.0) call initialize(SOURCE_LIST)

        
	call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(N_CHIP,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
c	call MPI_Bcast(eye_check,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)

	call MPI_Bcast(IMAGE_FILE,NMAX_EXPO*NMAX_CHIP*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_OUTPUT,NMAX_EXPO*strl,mpi_character
     .,0,MPI_COMM_WORLD,ierr)

        if (my_id.eq.0) then
	  filename=trim(DIR_OUT)//'/anomaly_expo.dat'
          call read_anomaly(filename)
        endif
	call MPI_Bcast(anomaly,NMAX_EXPO,mpi_float,0,MPI_COMM_WORLD,ierr)
        
        call mpi_distribute(N_EXPO,read_shear_cat,'reading cat')
        
        npp=ngfieldmax*NMAX_EXPO/num_procs_min
	xmax=0.005
	xmin=-0.005

        
	filename=trim(DIR_OUT)//'/full_shear_field_test_g1_opt'
	call plot_comparison_MPI(npp,ng,x1,y1,de1,40,yeqx,filename)

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) 

	xmax=0.005
	xmin=-0.005

        filename=trim(DIR_OUT)//'/full_shear_field_test_g2_opt'
	call plot_comparison_MPI(npp,ng,x2,y2,de2,40,yeqx,filename)

        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine read_anomaly(filename)
        implicit none
        include 'para.inc'

        character*(strl) filename
        real anomaly(NMAX_EXPO)
        common /anomaly_pass/ anomaly
        integer i,j,ierror,iexpo
        real v_chi2,tmp

        do i=1,NMAX_EXPO
          anomaly(i)=100.
        enddo
        
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) trim(filename),' is missing!'
          pause
        endif
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) tmp,v_chi2
          if (ierror.lt.0) cycle
          iexpo=int(tmp+0.5)
          anomaly(iexpo)=v_chi2    
	enddo
        close(10)
        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        subroutine read_shear_cat(iexpo)
        implicit none
        include 'para.inc'

        integer iexpo

        character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

c	integer eye_check(NMAX_EXPO)
c	common /eye_pass/ eye_check
        
	character*(strl) PREFIX,filename

	integer ichip,ierror,i,j,k,u
	real cat(npara),temp,gf,SNR
        integer n,m,mis
        save n
        data n /0/

	real x1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real x2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real y1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real y2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de2(ngfieldmax*NMAX_EXPO/num_procs_min)
        real de(ngfieldmax*NMAX_EXPO/num_procs_min)
        real yy1(ngfieldmax*NMAX_EXPO/num_procs_min)
        real yy2(ngfieldmax*NMAX_EXPO/num_procs_min)
        integer ng
        common /shear_data_pass/ x1,x2,y1,y2,de,de1,de2,yy1,yy2,ng        
        
        real anomaly(NMAX_EXPO)
        common /anomaly_pass/ anomaly
       
        m=0
        mis=0

c        if (eye_check(iexpo).ne.1) goto 50
 
        if (anomaly(iexpo).gt.chi2_thresh) goto 50
 
        do ichip=1,N_CHIP(iexpo)

          call get_PREFIX(IMAGE_FILE(iexpo,ichip),PREFIX)
	  filename=trim(DIR_OUTPUT(iexpo))
     .//'/result/'//trim(PREFIX)//'_shear.dat'        
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            cycle
          endif
          read(10,*)

          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (cat(u),u=1,ih2)
            if (ierror.lt.0) cycle	
            if (cat(istar).lt.20) cycle
            if (cat(iflux_alt).lt.2.) cycle

c            SNR=cat(ih_flux)/sqrt(cat(ih_area))
c            if (SNR.lt.30) cycle

c            if (cat(ih_area).lt.30) cycle

            if (abs(cat(igf1)).gt.0.005) cycle
            if (abs(cat(igf2)).gt.0.005) cycle
c            gf=sqrt(cat(igf1)**2+cat(igf2)**2)            
c            if (gf.gt.0.005) cycle
            
            if (cat(i_imax).ge.ns.or.cat(i_jmax).ge.ns) then
              mis=mis+1
              cycle
            endif
           
            m=m+1              
	    n=n+1
	    temp=cat(iflux_alt)**2
            temp=1./temp
              
  	    x1(n)=cat(igf1)
	    y1(n)=cat(ig1)
            yy1(n)=cat(ig1)*temp
	    de1(n)=cat(ide)-cat(ih1)              
  	    x2(n)=cat(igf2)
	    y2(n)=cat(ig2)
            yy2(n)=cat(ig2)*temp
	    de2(n)=cat(ide)+cat(ih1)              
	    de(n)=cat(ide)*temp
              
	    if (isnan(y1(n)).or.isnan(y2(n))
     ..or.isnan(de1(n)).or.isnan(de2(n))) then
              write(*,*) 'Something is wrong!'
              n=n-1
              m=m-1
	    endif
	    if (isnan(x1(n)).or.isnan(x2(n))) then
              write(*,*) 'Something is wrong!'
              n=n-1
              m=m-1
	    endif
	  enddo
          close(10)
        enddo

 50     write(*,*) iexpo,trim(IMAGE_FILE(iexpo,1)),m,mis,n
        ng=n
        
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
	subroutine plot_comparison_MPI(np,n,x,y,de,nb
     .,func,filename)
        use mpi        
	implicit none

        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs
	
	integer n,np,nbin,i,j,is,nb,ntot,ierr
        parameter (nbin=5)
	real func
	external func
	real x(np),y(np),de(np),yy(np),dd(np)
	character filename*(*),fname*200
	real xmin,xmax,dx
	common /bound_pass/ xmin,xmax

	real arr(nb,4),c,dc

        if (my_id.eq.0) then
          fname=trim(filename)//'.dat'
          open(unit=50,file=fname)
          rewind 50
        endif
       
	dx=(xmax-xmin)/nb

	do i=1,nb
	  arr(i,1)=xmin+dx*(i-0.5)
	  is=0	
	  do j=1,n
	    if (x(j).ge.xmin+dx*(i-1).and.x(j).le.xmin+dx*i) then
	      is=is+1
	      yy(is)=y(j)
	      dd(is)=de(j)
	    endif
	  enddo
	  call MPI_BARRIER(MPI_COMM_WORLD, ierr ) 

	  call statis_MPI(is,np,yy,dd,nbin,c,dc)
          arr(i,2)=c
          arr(i,3)=dc
	  call MPI_Reduce(is,ntot,1,mpi_int,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)

          if (my_id.eq.0) then
	    arr(i,4)=ntot
	    write(*,*) i,(arr(i,j),j=1,4)
	    write(50,*) i,(arr(i,j),j=1,4)
          endif
	enddo

        if (my_id.eq.0) then
	  close(50)
          fname=trim(filename)//'.fits'
	  call map_function(func,nb,nb,arr,500,500,fname)
        endif
       
        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine shear_field_distortion_test()
        implicit none
        include 'para.inc'

	character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)	

	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT

	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	character*(strl) PREFIX,filename

	integer iexpo,ichip

	integer ierror,i,j,k
	real cat(npara),yeqx,temp
	external yeqx
	integer u,v

	integer npp
	real x1(ngfieldmax*NMAX_EXPO),x2(ngfieldmax*NMAX_EXPO)
        real y1(ngfieldmax*NMAX_EXPO),y2(ngfieldmax*NMAX_EXPO)
        real de1(ngfieldmax*NMAX_EXPO),de2(ngfieldmax*NMAX_EXPO)
        real de(ngfieldmax*NMAX_EXPO)
        real yy1(ngfieldmax*NMAX_EXPO),yy2(ngfieldmax*NMAX_EXPO)
	integer n,nn,ig,m,mis

	real aa,bb

	real xmin,xmax
	common /bound_pass/ xmin,xmax

        npp=ngfieldmax*NMAX_EXPO

        n=0

        do iexpo=1,N_EXPO
           m=0
           mis=0
           
          do ichip=1,N_CHIP(iexpo)

            call get_PREFIX(IMAGE_FILE(iexpo,ichip),PREFIX)
	    filename=trim(DIR_OUTPUT(iexpo))
     .//'/result/'//trim(PREFIX)//'_shear.dat'
        
            open(unit=10,file=filename,status='old',iostat=ierror)
            rewind 10

            if (ierror.ne.0) then
              write(*,*) trim(filename),' is missing!'
              cycle
            endif
            read(10,*)

            do while (ierror.ge.0)
              read(10,*,iostat=ierror) (cat(u),u=1,ih2)
              if (ierror.lt.0) cycle	
c             if (cat(iflag).lt.5.) cycle

              if (cat(istar).lt.20.) cycle
              if (cat(iflux_alt).lt.2.) cycle
              if (cat(i_imax).ge.ns) then
                 mis=mis+1
                 cycle
              endif
              if (cat(i_jmax).ge.ns) then
                 mis=mis+1
                 cycle
              endif
              m=m+1
              
c-------correct for offset------------------------
c	      aa=-102.329*(cat(ipeak)+0.1021)**(-5.39557)
c     .+151.28*cat(ipeak)**(-5.957)
c	      bb=0.023115*cat(ipeak)**(-1.43919)
c             cat(ig1)=cat(ig1)-aa*cat(ide)
c     .+(aa*cat(ih1)+bb*cat(ih2))
c             cat(ig2)=cat(ig2)-bb*cat(ide)
c     .+(aa*cat(ih2)-bb*cat(ih1))
c--------------------------------------------------------------------------

c-------correct for field distortion------------------------
c             cat(ig1)=cat(ig1)-cat(igf1)*cat(ide)
c     .+(cat(igf1)*cat(ih1)+cat(igf2)*cat(ih2))
c             cat(ig2)=cat(ig2)-cat(igf2)*cat(ide)
c     .+(cat(igf1)*cat(ih2)-cat(igf2)*cat(ih1))
c--------------------------------------------------------------------------
c	      temp=cat(ipeak)**(-2)
c	      temp=cat(iflux_alt)**(-1)
	      temp=cat(iflux_alt)**2
              temp=1./temp
              
	      n=n+1

  	      x1(n)=cat(igf1)
	      y1(n)=cat(ig1)
              yy1(n)=cat(ig1)*temp
	      de1(n)=cat(ide)-cat(ih1)
              
  	      x2(n)=cat(igf2)
	      y2(n)=cat(ig2)
              yy2(n)=cat(ig2)*temp
	      de2(n)=cat(ide)+cat(ih1)
              
	      de(n)=cat(ide)*temp
              
	      if (isnan(y1(n)).or.isnan(y2(n))
     ..or.isnan(de1(n)).or.isnan(de2(n))) then
                write(*,*) 'Something is wrong!'
                n=n-1
	      endif
	      if (isnan(x1(n)).or.isnan(x2(n))) then
                write(*,*) 'Something is wrong!'
                n=n-1
	      endif
	    enddo
            close(10)
          enddo
          write(*,*) iexpo,n,filename,m,mis
        enddo


	xmax=0.005
	xmin=-0.005

	filename='full_shear_field_test_g1_ave'
	call plot_comparison(npp,n,x1,yy1,de,40,0,yeqx,filename)
	filename='full_shear_field_test_g1_opt'
	call plot_comparison(npp,n,x1,y1,de1,40,1,yeqx,filename)

	xmax=0.005
	xmin=-0.005
	filename='full_shear_field_test_g2_ave'
        call plot_comparison(npp,n,x2,yy2,de,40,0,yeqx,filename)
	filename='full_shear_field_test_g2_opt'
	call plot_comparison(npp,n,x2,y2,de2,40,1,yeqx,filename)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






