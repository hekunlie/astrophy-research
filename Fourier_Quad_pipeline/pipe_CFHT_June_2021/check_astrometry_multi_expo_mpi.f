        subroutine astrometry_multi_expo_MPI(SOURCE_LIST,DIR_OUT)
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

	integer ierr,i,j,k

        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs
     .,PUs,validity

        external fit_CDs,fit_CDs_only
        
        do i=1,NMAX_EXPO
          validity(i)=0 
          do j=1,2
            do k=1,2  
              CDs(i,j,k)=0d0
            enddo
            CRVALs(i,j)=0d0
          enddo
        enddo

        do i=1,NMAX_CHIP
          do j=1,2 
            CRPIXs(i,j)=0d0
          enddo
        enddo

        do i=1,2
          do j=1,npd
            PUs(i,j)=0d0
          enddo
        enddo
           
	if (my_id.eq.0) call initialize(SOURCE_LIST)
        
	call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(N_CHIP,NMAX_EXPO,mpi_int,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(IMAGE_FILE,NMAX_EXPO*NMAX_CHIP*strl
     .,mpi_character,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(DIR_OUTPUT,NMAX_EXPO*strl,mpi_character
     .,0,MPI_COMM_WORLD,ierr)

        call mpi_distribute(N_EXPO,fit_CDs,'fitting CDs')
        call fit_PUs()
        call mpi_distribute(N_EXPO,fit_CDs_only,'fitting CDs_only')
        
       
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine readin_astro_cat(iexpo)
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
        
	integer ichip        
	character*(strl) PREFIX,filename
	integer nss(NMAX_CHIP),n_user,n_ref,i
	integer nss_max
	parameter (nss_max=10000)
	double precision ra2(NMAX_CHIP,nss_max),dec2(NMAX_CHIP,nss_max)
	double precision x2(NMAX_CHIP,nss_max),y2(NMAX_CHIP,nss_max)
	double precision CRVAL2(2),tmp1,tmp2
        common /astro_cat_data_pass/ ra2,dec2,x2,y2,nss,CRVAL2
        
        do ichip=1,N_CHIP(iexpo)

	  call get_PREFIX(IMAGE_FILE(iexpo,ichip),PREFIX)
	  filename=trim(DIR_OUTPUT(iexpo))//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
	   
          open(unit=10,file=filename,status='old')
          rewind 10
	  read(10,*) tmp1,tmp2,CRVAL2(1),CRVAL2(2)
	  read(10,*) 
          read(10,*) nss(ichip),n_user,n_ref
	  do i=1,nss(ichip)
	    read(10,*) ra2(ichip,i),dec2(ichip,i)
     .,x2(ichip,i),y2(ichip,i)
	  enddo
	  close(10)

        enddo

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_PUs(iexpo)
        use mpi
        implicit none
        include 'para.inc'

        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs

        integer iexpo
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO
        
        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs
     .,PUs,validity
        integer tmp(NMAX_EXPO),i,j,k,ierr,nC,nn,order,px,py
        double precision tmp2(NMAX_EXPO,2),tmp3(NMAX_EXPO,2,2)
        external fit_PUs_accre

	double precision bvec(ndim),matx(ndim,ndim)
        common /astro_matx_pass/ matx,bvec

        double precision matx_tot(ndim,ndim),bvec_tot(ndim),vec(ndim)
        double precision matx_1(ndim,ndim)
        
	call MPI_Allreduce(validity,tmp,NMAX_EXPO,mpi_int,MPI_SUM
     .,MPI_COMM_WORLD,ierr)
	call MPI_Allreduce(CRVALs,tmp2,NMAX_EXPO*2,mpi_double
     .,MPI_SUM,MPI_COMM_WORLD,ierr)
	call MPI_Allreduce(CDs,tmp3,NMAX_EXPO*4,mpi_double,MPI_SUM
     .,MPI_COMM_WORLD,ierr)

        do i=1,NMAX_EXPO
          do j=1,2
            do k=1,2
              CDs(i,j,k)=tmp3(i,j,k)
            enddo
            CRVALs(i,j)=tmp2(i,j)
          enddo   
          validity(i)=tmp(i)
        enddo

        do i=1,ndim
          do j=1,ndim
            matx(i,j)=0d0
            matx_tot(i,j)=0d0 
          enddo
          bvec(i)=0d0
          bvec_tot(i)=0d0
          vec(i)=0d0
        enddo
       
        call mpi_distribute(N_EXPO,fit_PUs_accre,'adding chi2:PUs')
	call MPI_Reduce(matx,matx_tot,ndim*ndim,mpi_double,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)
	call MPI_Reduce(bvec,bvec_tot,ndim,mpi_double,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)

        nC=nchip_max
        
        if (my_id.eq.0) then

          j=npd*2+nC*2 

          do i=npd*2+1,npd*2+nC
            matx_tot(i,j+1)=0.5d0
            matx_tot(j+1,i)=0.5d0
          enddo
          do i=npd*2+nC+1,npd*2+nC+nC
            matx_tot(i,j+2)=0.5d0
            matx_tot(j+2,i)=0.5d0
          enddo

          k=j+2
	  call matrix_inverse_doub(matx_tot,k,ndim,matx_1)
	  do i=1,k
	    vec(i)=0d0
	    do j=1,k
	      vec(i)=vec(i)+matx_1(i,j)*bvec_tot(j)
	    enddo
	  enddo

          px=0
	  py=1
  	  order=1
	  nn=0
	  do while (nn.lt.npd)
            if (py.eq.order) then
	      order=order+1
	      px=order
	      py=0
	    else
	      px=px-1
	      py=py+1
	    endif
	    nn=nn+1
	    PUs(1,nn)=vec(nn)
	    PUs(2,nn)=vec(npd+nn+order-py*2)
          enddo
          do i=1,nC
            CRPIXs(i,1)=vec(npd*2+i)
            CRPIXs(i,2)=vec(npd*2+nC+i)
          enddo
          
        endif

	call MPI_Bcast(PUs,npd*2,mpi_double,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(CRPIXs,NMAX_CHIP*2,mpi_double,0
     .,MPI_COMM_WORLD,ierr)
        
                
c chi^2=\sum_i{CD11*x(i)+CD12*y(i)-CD11*CRPIX(k,1)-CD12*CRPIX(k,2)-xi(i)+PU(1,1)*r
c +PU(1,2)*xi(i)^2+PU(1,3)*xi(i)*eta(i)+PU(1,4)*eta(i)^2
c +PU(1,5)*xi(i)^3+PU(1,6)*xi(i)^2*eta(i)+PU(1,7)*xi(i)*eta(i)^2
c +PU(1,8)*eta(i)^3+PU(1,9)*r^3+... ...}^2

c +\lambda_1\sum_k{CRPIX(k,1)}+\lambda_2\sum_k{CRPIX(k,2)}      

c +\sum_i{CD21*x(i)+CD22*y(i)-CD21*CRPIX(k,1)-CD22*CRPIX(k,2)-eta(i)+PU(2,1)*r
c +PU(2,2)*eta(i)^2+PU(2,3)*xi(i)*eta(i)+PU(2,4)*xi(i)^2
c +PU(2,5)*eta(i)^3+PU(2,6)*eta(i)^2*xi(i)+PU(2,7)*eta(i)*xi(i)^2
c +PU(2,8)*xi(i)^3+PU(2,9)*r^3+... ...}^2

c vector order:
c PU(1,1),PU(1,2),...,PU(2,1),PU(2,2),...,CRPIX(1,1),CRPIX(2,1),...,CRPIX(1,2),CRPIX(2,2),...,\lambda_1,\lambda_2

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_PUs_accre(iexpo)
        implicit none
        include 'para.inc'

        integer iexpo
        
        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs
     .,PUs,validity
        
        if (validity(iexpo).ne.1) return
       
        call readin_astro_cat(iexpo)
        
        call update_PUs(iexpo)


        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_CDs(iexpo)
        implicit none
        include 'para.inc'

        integer iexpo
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

      
        if (N_CHIP(iexpo).ne.nchip_max) return
       
        call readin_astro_cat(iexpo)
	        
        call update_CDs(iexpo)


        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_CDs_only(iexpo)
        implicit none
        include 'para.inc'

        integer iexpo
    
        call readin_astro_cat(iexpo)
	        
        call update_CDs_only(iexpo)


        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_PUs(iexpo)
	implicit none
        include 'para.inc'
        
	integer iexpo
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	integer nss(NMAX_CHIP),nss_max
	parameter (nss_max=10000)
	double precision ra(NMAX_CHIP,nss_max),dec(NMAX_CHIP,nss_max)
	double precision x(NMAX_CHIP,nss_max),y(NMAX_CHIP,nss_max)
	double precision CRVAL(2)
        common /astro_cat_data_pass/ ra,dec,x,y,nss,CRVAL
        
        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs,PUs,validity

        double precision vec(ndim),vec1(ndim),vec2(ndim)       
	double precision matx_1(ndim,ndim),r

	double precision bvec(ndim),matx(ndim,ndim)
        common /astro_matx_pass/ matx,bvec

        integer tot_source,ichip,nC,i,j,k,l,u,v,nn,ntot,order,px,py
	double precision xi,eta,PU(2,npd),s_expo(2),CD(2,2)

        nC=N_CHIP(iexpo)

        CRVAL(1)=CRVALs(iexpo,1)
        CRVAL(2)=CRVALs(iexpo,2)
        CD(1,1)=CDs(iexpo,1,1)
        CD(1,2)=CDs(iexpo,1,2)
        CD(2,1)=CDs(iexpo,2,1)
        CD(2,2)=CDs(iexpo,2,2)
        
	do k=1,nC
	  do i=1,nss(k)
	    call ra_dec_to_xi_eta(ra(k,i),dec(k,i),xi,eta
     .,CRVAL(1),CRVAL(2))

            do j=1,ndim
              vec(j)=0d0
            enddo
            
	    px=0
	    py=1
	    order=1
	    nn=0

	    do while (nn.lt.npd)
	      if (py.eq.order) then
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      nn=nn+1
	      vec(nn)=xi**px*eta**py
	    enddo
            
            vec(npd*2+k)=-CD(1,1)
            vec(npd*2+nC+k)=-CD(1,2)

	    do u=1,ndim
	      do v=1,ndim
	        matx(u,v)=matx(u,v)+vec(u)*vec(v)
	      enddo
	      bvec(u)=bvec(u)+(xi-CD(1,1)*x(k,i)-CD(1,2)*y(k,i))*vec(u)
	    enddo

            do j=1,npd
              vec(j+npd)=vec(j)
              vec(j)=0d0 
            enddo

            vec(npd*2+k)=-CD(2,1)
            vec(npd*2+nC+k)=-CD(2,2)

	    do u=1,ndim
	      do v=1,ndim
	        matx(u,v)=matx(u,v)+vec(u)*vec(v)
	      enddo
	      bvec(u)=bvec(u)+(eta-CD(2,1)*x(k,i)-CD(2,2)*y(k,i))*vec(u)
	    enddo
            
	  enddo
	enddo	    	


        
c r=sqrt(xi(i)^2+eta(i)^2)

c chi^2=\sum_i{CD11*x(i)+CD12*y(i)-CD11*CRPIX(k,1)-CD12*CRPIX(k,2)-xi(i)+PU(1,1)*r
c +PU(1,2)*xi(i)^2+PU(1,3)*xi(i)*eta(i)+PU(1,4)*eta(i)^2
c +PU(1,5)*xi(i)^3+PU(1,6)*xi(i)^2*eta(i)+PU(1,7)*xi(i)*eta(i)^2
c +PU(1,8)*eta(i)^3+PU(1,9)*r^3+... ...}^2

c     +\lambda_1\sum_k{CRPIX(k,1)}+\lambda_2\sum_k{CRPIX(k,2)}      

c      +\sum_i{CD21*x(i)+CD22*y(i)-CD21*CRPIX(k,1)-CD22*CRPIX(k,2)-eta(i)+PU(2,1)*r
c +PU(2,2)*eta(i)^2+PU(2,3)*xi(i)*eta(i)+PU(2,4)*xi(i)^2
c +PU(2,5)*eta(i)^3+PU(2,6)*eta(i)^2*xi(i)+PU(2,7)*eta(i)*xi(i)^2
c +PU(2,8)*xi(i)^3+PU(2,9)*r^3+... ...}^2

c     vector order:
c     PU(1,1),PU(1,2),...,PU(2,1),PU(2,2),...,CRPIX(1,1),CRPIX(2,1),...,CRPIX(1,2),CRPIX(2,2),...,\lambda_1,\lambda_2
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_CDs_only(iexpo)
	implicit none
        include 'para.inc'
        
	integer iexpo
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO
        character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT        
	character*(strl) PREFIX,filename

	integer nss(NMAX_CHIP),nss_max
	parameter (nss_max=10000)
	double precision ra(NMAX_CHIP,nss_max),dec(NMAX_CHIP,nss_max)
	double precision x(NMAX_CHIP,nss_max),y(NMAX_CHIP,nss_max)
	double precision CRVAL(2)
        common /astro_cat_data_pass/ ra,dec,x,y,nss,CRVAL
        
        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs,PUs,validity

	double precision vec(nchip_max+2),vec1(nchip_max+2)
        double precision vec2(nchip_max+2),tmp1,tmp2
	double precision bvec1(nchip_max+2),bvec2(nchip_max+2)
	double precision matx(nchip_max+2,nchip_max+2)
	double precision matx_1(nchip_max+2,nchip_max+2),r

        integer tot_source,ichip,nC,i,j,k,l,u,v,nn,ntot,order,px,py
	double precision xi,eta,s_expo(2),CD(2,2)
        integer valid(nchip_max),iC
        
        tot_source=0
	do ichip=1,N_CHIP(iexpo)
  	  if (nss(ichip).ge.10) then
            tot_source=tot_source+nss(ichip)
            valid(ichip)=1 
 	  else
 	    valid(ichip)=0
	  endif
	enddo
        if (tot_source.lt.nchip_max*3) then
          do ichip=1,N_CHIP(iexpo)
            valid(ichip)=0
          enddo
          goto 80
        endif
        
        nC=nchip_max        
	ntot=nC+2

        do u=1,ntot
          do v=1,ntot
	    matx(u,v)=0d0
	  enddo
	  bvec1(u)=0d0
	  bvec2(u)=0d0
        enddo

        iC=0
	do k=1,nC
          if (valid(k).eq.0) cycle
          iC=iC+1 
	  do i=1,nss(k)
	    call ra_dec_to_xi_eta(ra(k,i),dec(k,i),xi,eta
     .,CRVAL(1),CRVAL(2))

            tmp1=0d0
            tmp2=0d0           
	    px=0
	    py=1
	    order=1
	    nn=0
	    do while (nn.lt.npd)
	      if (py.eq.order) then
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      nn=nn+1
	      tmp1=tmp1+PUs(1,nn)*xi**px*eta**py
              tmp2=tmp2+PUs(2,nn)*eta**px*xi**py
	    enddo

	    do j=1,ntot
	      vec(j)=0d0
	    enddo
            
	    vec(1)=x(k,i)
	    vec(2)=y(k,i)
	    vec(2+iC)=1d0

	    do u=1,2+iC
	      do v=1,2+iC
	        matx(u,v)=matx(u,v)+vec(u)*vec(v)
	      enddo
	      bvec1(u)=bvec1(u)+(xi-tmp1)*vec(u)
	      bvec2(u)=bvec2(u)+(eta-tmp2)*vec(u)
	    enddo

	  enddo
	enddo	    	

        k=2+iC
        
	call matrix_inverse_doub(matx,k,ntot,matx_1)

	do i=1,k
	  vec1(i)=0d0
	  vec2(i)=0d0
	  do j=1,k
	    vec1(i)=vec1(i)+matx_1(i,j)*bvec1(j)
	    vec2(i)=vec2(i)+matx_1(i,j)*bvec2(j)
	  enddo
	enddo

	CD(1,1)=vec1(1)
	CD(1,2)=vec1(2)
	CD(2,1)=vec2(1)
	CD(2,2)=vec2(2)

        s_expo(1)=0d0
        s_expo(2)=0d0
        iC=0
	do k=1,nC
          if (valid(k).eq.0) cycle
          iC=iC+1 
          s_expo(1)=s_expo(1)+vec1(2+iC)
     .+CD(1,1)*CRPIXs(k,1)+CD(1,2)*CRPIXs(k,2)
          s_expo(2)=s_expo(2)+vec2(2+iC)
     .+CD(2,1)*CRPIXs(k,1)+CD(2,2)*CRPIXs(k,2)
        enddo

        s_expo(1)=s_expo(1)/iC
        s_expo(2)=s_expo(2)/iC

        call update_CRVAL(s_expo,CRVAL,PUs,npd)

 80     call get_PREFIX(IMAGE_FILE(iexpo,1),PREFIX)
        filename=trim(DIR_OUTPUT(iexpo))//'/astrometry/'
     .//trim(PREFIX)//'.head2'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) CRVAL(1),CRVAL(2)
        do i=1,npd
	  write(10,*) PUs(1,i),PUs(2,i)
	enddo
	do k=1,N_CHIP(iexpo)
	  write(10,*) k,valid(k),CRPIXs(k,1),CRPIXs(k,2)
     .,CD(1,1),CD(1,2),CD(2,1),CD(2,2)
        enddo
        close(10)        
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_CDs(iexpo)
	implicit none
        include 'para.inc'
        
	integer iexpo
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO

	integer nss(NMAX_CHIP),nss_max
	parameter (nss_max=10000)
	double precision ra(NMAX_CHIP,nss_max),dec(NMAX_CHIP,nss_max)
	double precision x(NMAX_CHIP,nss_max),y(NMAX_CHIP,nss_max)
	double precision CRVAL(2)
        common /astro_cat_data_pass/ ra,dec,x,y,nss,CRVAL
        
        double precision CDs(NMAX_EXPO,2,2),PUs(2,npd)
        double precision CRVALs(NMAX_EXPO,2),CRPIXs(NMAX_CHIP,2)
        integer validity(NMAX_EXPO)
        common /astrometry_multi_pass/ CDs,CRVALs,CRPIXs,PUs,validity

	double precision vec(npd+nchip_max+2),vec1(npd+nchip_max+2)
        double precision vec2(npd+nchip_max+2)
	double precision bvec1(npd+nchip_max+2),bvec2(npd+nchip_max+2)
	double precision matx(npd+nchip_max+2,npd+nchip_max+2)
	double precision matx_1(npd+nchip_max+2,npd+nchip_max+2),r

        integer tot_source,ichip,nC,i,j,k,l,u,v,nn,ntot,order,px,py
	double precision xi,eta,PU(2,npd),s_expo(2),CD(2,2)

        tot_source=0
	do ichip=1,N_CHIP(iexpo)
  	  if (nss(ichip).ge.10) then
	    tot_source=tot_source+nss(ichip)
 	  else
 	    return
	  endif
	enddo
        if (tot_source.lt.(npd+nchip_max)*3) return

        nC=nchip_max
        
	ntot=npd+nC+2

        do u=1,ntot
          do v=1,ntot
	    matx(u,v)=0d0
	  enddo
	  bvec1(u)=0d0
	  bvec2(u)=0d0
        enddo

	do k=1,nC
	  do i=1,nss(k)
	    call ra_dec_to_xi_eta(ra(k,i),dec(k,i),xi,eta
     .,CRVAL(1),CRVAL(2))

	    do j=1,ntot
	      vec(j)=0d0
	    enddo

c	    r=sqrt(xi(k,i)**2+eta(k,i)**2)
	    px=0
	    py=1
	    order=1
	    nn=0

	    do while (nn.lt.npd)
	      if (py.eq.order) then
c	        if (mod(order,2).eq.1) then
c		  nn=nn+1
c	          vec(nn)=r**order
c	          if (nn.eq.npd) cycle
c	        endif
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      nn=nn+1
	      vec(nn)=xi**px*eta**py
	    enddo

	    vec(npd+1)=x(k,i)
	    vec(npd+2)=y(k,i)
	    vec(npd+2+k)=1d0

	    j=npd+2+k

	    do u=1,j
	      do v=1,j
	        matx(u,v)=matx(u,v)+vec(u)*vec(v)
	      enddo
	      bvec1(u)=bvec1(u)+xi*vec(u)
	      bvec2(u)=bvec2(u)+eta*vec(u)
	    enddo

	  enddo
	enddo	    	

	call matrix_inverse_doub(matx,ntot,ntot,matx_1)

	do i=1,ntot
	  vec1(i)=0d0
	  vec2(i)=0d0
	  do j=1,ntot
	    vec1(i)=vec1(i)+matx_1(i,j)*bvec1(j)
	    vec2(i)=vec2(i)+matx_1(i,j)*bvec2(j)
	  enddo
	enddo


        px=0
	py=1
	order=1
	nn=0
	do while (nn.lt.npd)
          if (py.eq.order) then
	    order=order+1
	    px=order
	    py=0
	  else
	    px=px-1
	    py=py+1
	  endif
	  nn=nn+1
	  PU(1,nn)=vec1(nn)
	  PU(2,nn)=vec2(nn+order-py*2)
	enddo
	CD(1,1)=vec1(npd+1)
	CD(1,2)=vec1(npd+2)
	CD(2,1)=vec2(npd+1)
	CD(2,2)=vec2(npd+2)

        s_expo(1)=0d0
        s_expo(2)=0d0
	do k=1,nC
          s_expo(1)=s_expo(1)+vec1(npd+2+k)
          s_expo(2)=s_expo(2)+vec2(npd+2+k)
        enddo

        s_expo(1)=s_expo(1)/nC
        s_expo(2)=s_expo(2)/nC

        call update_CRVAL(s_expo,CRVAL,PU,npd)

        validity(iexpo)=1
        do i=1,2
          do j=1,2
            CDs(iexpo,i,j)=CD(i,j)
          enddo
          CRVALs(iexpo,i)=CRVAL(i)
        enddo
        
c r=sqrt(xi(i)^2+eta(i)^2)

c chi^2=\sum_i{CD11(k)*x(i)+CD12(k)*y(i)+S1(k)-xi(i)+PU(1,1)*r
c +PU(1,2)*xi(i)^2+PU(1,3)*xi(i)*eta(i)+PU(1,4)*eta(i)^2
c +PU(1,5)*xi(i)^3+PU(1,6)*xi(i)^2*eta(i)+PU(1,7)*xi(i)*eta(i)^2
c +PU(1,8)*eta(i)^3+PU(1,9)*r^3+... ...}^2

c      +\sum_i{CD21(k)*x(i)+CD22(k)*y(i)+S2(k)-eta(i)+PU(2,1)*r
c +PU(2,2)*eta(i)^2+PU(2,3)*xi(i)*eta(i)+PU(2,4)*xi(i)^2
c +PU(2,5)*eta(i)^3+PU(2,6)*eta(i)^2*xi(i)+PU(2,7)*eta(i)*xi(i)^2
c +PU(2,8)*xi(i)^3+PU(2,9)*r^3+... ...}^2


	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine update_CRVAL(s_expo,CRVAL,PU,npd)
	implicit none
	
	integer npd
	double precision s_expo(2),CRVAL(2),PU(2,npd)

	double precision xxx,yyy,a,d,da,dd,const1,xi,eta
	double precision cosda,tandc,sumra,diffra
	double precision pi
	parameter (pi=3.1415926d0)

	const1=pi/180d0
	tandc=tan(CRVAL(2)*const1)

	call mapping_PU(s_expo(1),s_expo(2),xi,eta,npd,PU,1)

        xxx=xi*const1
	yyy=eta*const1

	da=xxx/(cos(CRVAL(2)*const1)*(1.-yyy*tandc))
	da=da-da*da*da*0.33333333333
	a=sumra(da/const1,CRVAL(1))

	cosda=1.-da*da*0.5+da**4/24.
	dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)	  
	dd=dd-dd*dd*dd*0.3333333333	  	 
	d=dd/const1+CRVAL(2)

        CRVAL(1)=a
        CRVAL(2)=d

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
