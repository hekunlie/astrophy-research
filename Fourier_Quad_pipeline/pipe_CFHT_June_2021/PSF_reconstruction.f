	subroutine PSF_reconstruction(IMAGE_FILE,nchip,DIR_OUTPUT)
        implicit none   
	include 'para.inc'

	integer ichip,nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT
	
	integer nstar(NMAX_CHIP)
	double precision star_para(NMAX_CHIP,nstar_max,npara)
	real star_test(NMAX_CHIP,nstar_max,ns,ns)
	common /star_info_pass/ star_test,star_para,nstar

	real star_orig(NMAX_CHIP,nstar_max,ns,ns)
	real star(NMAX_CHIP*nstar_max,ns,ns)
	
	real ee(2),source_p(ns,ns),aa(npara),model(ns,ns),size
	integer ierror,i,j,k,u,v,nn1,nn2,area,nums
	character*(strl) PREFIX,filename,headname

	double precision CRPIX(2),CD(2,2),PU(2,npd),CRVAL(2)
	double precision x,y,xx,yy,step

	integer nc
	real p_chip(nchip,4)
	
	integer ntot,status
	double precision PSF_coe(ns,ns,npo),PSF_coe_l(ns,ns,npl)
	double precision posi(NMAX_CHIP*nstar_max,2)

	integer nm
	parameter (nm=1000)
	real PSFmap(nm,nm),sk(NMAX_CHIP*nstar_max,5)

	integer nmax_stamp,len_sam,nbad
	parameter (nmax_stamp=5000)
	parameter (len_sam=100)

	
        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        headname=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'

	nc=0
	do ichip=1,nchip

	  nstar(ichip)=0
	  ierror=0

          call read_astrometry_para(headname,ichip
     .,CRPIX,CD,CRVAL,PU,npd,ierror)	       
	  if (ierror.eq.1) cycle

	  nc=nc+1
	  x=1d0
	  y=1d0
	  call xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	  p_chip(nc,1)=xx
	  p_chip(nc,2)=yy
	  x=2000d0
	  y=4000d0
	  call xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	  p_chip(nc,3)=xx
	  p_chip(nc,4)=yy
	  

	  
          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
          PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

	  filename=trim(PREFIX)//'_star_can_info.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
	    write(*,*) filename
	    pause 'Catalog file error!!'
          endif   
          read(10,*) 
c 	  read(10,*) 'ig xp yp SNR' 
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (aa(i),i=1,4)
 	    if (ierror.lt.0) cycle
	    nstar(ichip)=nstar(ichip)+1
	    do i=1,4
	      star_para(ichip,nstar(ichip),i)=aa(i)
	    enddo

          enddo
          close(10)

	  nn1=ns*len_s
 	  nn2=ns*(int(nstar(ichip)/len_s)+1)
          filename=trim(PREFIX)//'_star_can_power.fits'
	  call read_stamps(NMAX_CHIP*nstar_max,1,nstar(ichip),ns,ns
     .,star,nn1,nn2,filename)

	  nums=nstar(ichip)
	  nstar(ichip)=0

	  do i=1,nums
	    if (star_para(ichip,i,4).ge.SNR_PSF) then
	      nstar(ichip)=nstar(ichip)+1
	      do j=1,4
	        star_para(ichip,nstar(ichip),j)=star_para(ichip,i,j)
	      enddo
	      star_para(ichip,nstar(ichip),5)=1.
	      do u=1,ns
	        do v=1,ns
	  	  star_orig(ichip,nstar(ichip),u,v)=star(i,u,v)
		  source_p(u,v)=star(i,u,v)
	        enddo
	      enddo

              call smooth_image55_hole(ns,ns,source_p)
	      do u=1,ns
	        do v=1,ns
  	  	  star_test(ichip,nstar(ichip),u,v)=source_p(u,v)
	        enddo
	      enddo
	      
	      x=star_para(ichip,i,2)
	      y=star_para(ichip,i,3)
	      call xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	      star_para(ichip,nstar(ichip),6)=xx
	      star_para(ichip,nstar(ichip),7)=yy

c	      call get_power_all(ns,ns,source_p,ee,size,0.2)
c	      star_para(ichip,nstar(ichip),8)=size
	      call get_power_area(ns,ns,source_p,area,0.2)
	      call get_power_e(ns,ns,source_p,ee,0.2)
	      star_para(ichip,nstar(ichip),8)=area
	      star_para(ichip,nstar(ichip),9)=ee(1)
	      star_para(ichip,nstar(ichip),10)=ee(2)
	      

	    endif
	  enddo
c	  write(*,*) 'reading star cans:',nstar(ichip)
	enddo
c----------------------------------------------------------------------------
c!
	call star_selection_2018(nchip,status)

	filename=trim(PREFIX)//'_star_can_info_expo.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
	write(10,*) '#  star_id  area  e1  e2'
        do k=1,nchip
	  do i=1,nstar(k)
	    write(10,*) star_para(k,i,5),(star_para(k,i,u),u=8,10)    
          enddo
	enddo
	close(10)

	
	ntot=0
	do k=1,nchip
	  nums=0
	  do i=1,nstar(k)
	    if (star_para(k,i,5).lt.0) cycle
	    ntot=ntot+1
	    nums=nums+1
	    posi(nums,1)=star_para(k,i,2)
	    posi(nums,2)=star_para(k,i,3)
	    do u=1,ns
	      do v=1,ns
	        star_test(k,i,u,v)=star_orig(k,i,u,v)
	        star(nums,u,v)=star_orig(k,i,u,v)
	      enddo
	    enddo	    
	  enddo
          call get_PREFIX(IMAGE_FILE(k),PREFIX)
          PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
  	  filename=trim(PREFIX)//'_PSF_coe_local.dat'
          open(unit=10,file=filename,status='replace')
          rewind 10
	  if (nums.ge.nstar_min_local) then
            call interpolate_PSF(nums,NMAX_CHIP*nstar_max,star
     .,posi,ns,npl,PSF_coe_l)
	    write(10,*) nums,1
	    do i=1,ns
	      do j=1,ns
	        write(10,*) (PSF_coe_l(i,j,u),u=1,npl)
	      enddo
	    enddo		 
	  else
	    write(10,*) nums,-1
	  endif
	  close(10)
	enddo

	ntot=0
	do k=1,nchip
	  do i=1,nstar(k)
	    if (star_para(k,i,5).lt.0) cycle
	    ntot=ntot+1
	    posi(ntot,1)=star_para(k,i,6)
	    posi(ntot,2)=star_para(k,i,7)
	    do u=1,ns
	      do v=1,ns
	        star(ntot,u,v)=star_orig(k,i,u,v)
	      enddo
	    enddo	    
	  enddo
	enddo 
        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        filename=trim(PREFIX)//'_PSF_coe.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
        if (ntot.ge.nstar_min) then
          call interpolate_PSF(ntot,NMAX_CHIP*nstar_max,star
     .,posi,ns,npo,PSF_coe)
	  write(10,*) ntot,1
	  do i=1,ns
	    do j=1,ns
	      write(10,*) (PSF_coe(i,j,u),u=1,npo)
	    enddo
	  enddo		 
        else
	  write(10,*) ntot,-1
        endif
	close(10)
	
c--------------------------------------------------------------------	
	filename=trim(PREFIX)//'_star_info_expo.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10

	ntot=0
        do k=1,nchip
          do i=1,nstar(k)
	    if (star_para(k,i,5).lt.0) cycle
	    ntot=ntot+1
	    sk(ntot,1)=star_para(k,i,6)
	    sk(ntot,2)=star_para(k,i,7)
	    sk(ntot,3)=star_para(k,i,8)
	    sk(ntot,4)=star_para(k,i,9)
	    sk(ntot,5)=star_para(k,i,10)
            if (ntot.le.nmax_stamp)
     . write(10,*) k,(star_para(k,i,j),j=1,7)
          enddo
	enddo
        write(*,*) trim(PREFIX),' n,status :',ntot,'/',status

	close(10)

	nums=min(ntot,nmax_stamp)
	
        nn1=ns*len_sam
        nn2=ns*(int(nums/len_sam)+1)
        filename=trim(PREFIX)//'_star_power_expo.fits'
        call write_stamps(nstar_max*NMAX_CHIP,1,nums,ns,ns
     .,star,nn1,nn2,filename)

	nbad=0
	do k=1,nchip
	   do i=1,nstar(k)
	      if (star_para(k,i,5).gt.0) cycle
	      nbad=nbad+1
	      do u=1,ns
		 do v=1,ns
		    star(nbad,u,v)=star_orig(k,i,u,v)
		 enddo
	      enddo	    
	   enddo
	enddo 
	nums=min(nbad,nmax_stamp)	
	nn1=ns*len_sam
	nn2=ns*(int(nums/len_sam)+1)
	filename=trim(PREFIX)//'_star_bad_expo.fits'
        call write_stamps(nstar_max*NMAX_CHIP,1,nums,ns,ns
     .,star,nn1,nn2,filename)

	
c----------make plots--------------------------------------------	

	call draw_shear_expo(nm,PSFmap,nchip,nc,p_chip
     .,NMAX_CHIP*nstar_max,ntot,sk,200.,1.)
	
        filename=trim(PREFIX)//'_PSF_source.fits'
	call writeimage(filename,nm,nm,nm,nm,PSFmap)

	step=500d0
	ntot=0
        do k=1,nchip
	  if (nstar(k).eq.0) cycle
	  ierror=0
          call read_astrometry_para(headname,k
     .,CRPIX,CD,CRVAL,PU,npd,ierror)	       
	  if (ierror.eq.1) cycle
          do i=1,3
	    do j=1,7	
	      x=i*step
	      y=j*step
	      call xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	      call get_PSF_model(ns,npo,PSF_coe,xx,yy,model)
	      call get_power_area(ns,ns,model,area,0.2)
	      call get_power_e(ns,ns,model,ee,0.1)
	      ntot=ntot+1
	      sk(ntot,1)=xx
	      sk(ntot,2)=yy
	      sk(ntot,3)=area
	      sk(ntot,4)=ee(1)
	      sk(ntot,5)=ee(2)
c	      write(*,*) i,j,(sk(ntot,u),u=1,5)
	    enddo
	  enddo
	enddo

	call draw_shear_expo(nm,PSFmap,nchip,nc,p_chip
     .,NMAX_CHIP*nstar_max,ntot,sk,200.,1.)	
        filename=trim(PREFIX)//'_PSF_model.fits'
	call writeimage(filename,nm,nm,nm,nm,PSFmap)

	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine star_selection_2018(nchip,status)
	implicit none
	include 'para.inc'

	integer nchip
	integer nstar(NMAX_CHIP)
        double precision star_para(NMAX_CHIP,nstar_max,npara)
        real star_test(NMAX_CHIP,nstar_max,ns,ns)
        common /star_info_pass/ star_test,star_para,nstar

	real star(nstar_max,ns,ns),posi(nstar_max,2)

	integer i,j,k,u,v,status,w,miss,nums,n1,n2
	real tmp(nstar_max),temp,source_p(ns,ns)
	real model(ns,ns),coe(npl,ns,ns),xy(2)
	real mean,chi2(nstar_max)

	
	status=1

	n1=ns_2+1-ns/4
	n2=ns_2+1+ns/4

	do k=1,nchip
	   do u=n1,n2
	      do v=n1,n2
		 do i=1,nstar(k)
		    tmp(i)=star_test(k,i,u,v)
		 enddo
		 call sort(nstar(k),nstar_max,tmp)
		 model(u,v)=tmp(nstar(k)*3/4)
	      enddo
	   enddo

	   do i=1,nstar(k)
	      do u=1,ns
		 do v=1,ns
		    source_p(u,v)=star_test(k,i,u,v)
		 enddo
	      enddo
	      call get_chi2(ns,n1,n2,source_p,model,chi2(i))
	      tmp(i)=chi2(i)
	   enddo
	   call sort(nstar(k),nstar_max,tmp)
	   do i=1,nstar(k)
	      if (chi2(i).gt.9.*tmp(nstar(k)/4)) star_para(k,i,5)=-1.
	   enddo

	   miss=1
	   do while (miss.gt.0)
	      nums=0
	      mean=0.
	      do i=1,nstar(k)
		 if (star_para(k,i,5).gt.0.) then
		    nums=nums+1
		    posi(nums,1)=star_para(k,i,2)
		    posi(nums,2)=star_para(k,i,3)
		    do u=1,ns
		       do v=1,ns
			  star(nums,u,v)=star_test(k,i,u,v)
		       enddo
		    enddo
		 endif
	      enddo
	      if (nums.lt.nstar_min_local) then
		 nstar(k)=0
		 goto 66
	      endif
	      call fit_map(nums,nstar_max,ns,n1,n2,star,posi,npl,coe)
	      do i=1,nstar(k)
		 if (star_para(k,i,5).gt.0.) then
		    do u=1,ns
		       do v=1,ns
			  source_p(u,v)=star_test(k,i,u,v)
		       enddo
		    enddo
		    xy(1)=star_para(k,i,2)
		    xy(2)=star_para(k,i,3)
		    call get_model(ns,n1,n2,xy,npl,coe,model)
		    call get_chi2(ns,n1,n2,source_p,model,chi2(i))
		    mean=mean+chi2(i)
		 endif
	      enddo
	      mean=mean/nums
	      miss=0
	      do i=1,nstar(k)
		 if (star_para(k,i,5).gt.0.) then
		    if (chi2(i).gt.4.*mean) then
		       star_para(k,i,5)=-1.
		       miss=miss+1
		    endif
		 endif
	      enddo
	   enddo

 66	enddo

 
	return
	END 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine star_selection(nchip,status)
	implicit none
	include 'para.inc'

	integer nchip
	integer nstar(NMAX_CHIP)
        double precision star_para(NMAX_CHIP,nstar_max,npara)
        real star_test(NMAX_CHIP,nstar_max,ns,ns)
        common /star_info_pass/ star_test,star_para,nstar

	integer npsam
	parameter (npsam=NMAX_CHIP*nstar_max*nstar_max)

	integer i,j,k,u,v,ntot,status,w
	real tmp(npsam),temp,map1(ns,ns),map2(ns,ns)
	real thresh,peak,sig,chimin(NMAX_CHIP,nstar_max)
	real chi_d(NMAX_CHIP,nstar_max,nstar_max)

	integer id(nstar_max),group_size(nstar_max)
	integer group_id(nstar_max),current_group_id
	
	status=1
	ntot=0
        do k=1,nchip
          ntot=ntot+nstar(k)
        enddo
	if (ntot.lt.nstar_min*2) then
	   nstar=0
	   status=-1
	   return
	endif

        do k=1,nchip	   
	  do i=1,nstar(k)
	    chimin(k,i)=1000.  

            temp=0.  
	    do u=1,ns
	       do v=1,ns
		  temp=temp+star_test(k,i,u,v)
	       enddo
	    enddo
	    temp=1./temp
	    do u=1,ns
	       do v=1,ns
		  star_test(k,i,u,v)=star_test(k,i,u,v)*temp
	       enddo
	    enddo
	  enddo
	enddo

	ntot=0
        do k=1,nchip 
	   do i=1,nstar(k)
	      ntot=ntot+1
	      tmp(ntot)=star_para(k,i,8)
	   enddo
	enddo
        call sort(ntot,npsam,tmp)
	thresh=tmp((ntot*2)/3)
	
c----Determine the threshold for chi^2 ------------------------------
	
	ntot=0
        do k=1,nchip 
	   do i=1,nstar(k)-1
	      do j=i+1,nstar(k)
		 do u=1,ns
		    do v=1,ns
		       map1(u,v)=star_test(k,i,u,v)
		       map2(u,v)=star_test(k,j,u,v)
		    enddo
		 enddo
		 call ana_chi2(ns,map1,map2,temp)
		 temp=sqrt(temp)
		 chimin(k,i)=min(temp,chimin(k,i))
		 chimin(k,j)=min(temp,chimin(k,j))
		 chi_d(k,i,j)=temp
		 chi_d(k,j,i)=temp
		 if (star_para(k,i,8).lt.thresh) cycle   
		 if (star_para(k,j,8).lt.thresh) cycle 
		 ntot=ntot+1	      
		 tmp(ntot)=temp
	      enddo 
	   enddo
	enddo

        call get_peak_width_low_side(npsam,ntot,tmp,peak,sig)

	thresh=peak+3.*sig

        do k=1,nchip
	   ntot=0
	   do i=1,nstar(k)
	      if (chimin(k,i).gt.thresh) then
		 star_para(k,i,5)=-1
		 cycle
	      endif
	      ntot=ntot+1
	   enddo
	   if (ntot.lt.nstar_min_local) then
	      nstar(k)=0
	      cycle
	   endif

	   ntot=0
	   do i=1,nstar(k)
	      if (star_para(k,i,5).lt.0) cycle
	      ntot=ntot+1
	      id(ntot)=i
	   enddo
	   do i=1,ntot
	      group_id(i)=0
	      group_size(i)=0
	   enddo
	   current_group_id=0
	   do i=1,ntot
	      if (group_id(i).eq.0) then
		 current_group_id=current_group_id+1
		 group_id(i)=current_group_id
	      endif	 
	      do j=i+1,ntot
		 if (group_id(j).eq.group_id(i)) cycle	     
		 if (chi_d(k,id(i),id(j)).gt.thresh) cycle
		 if (group_id(j).eq.0) then
		    group_id(j)=group_id(i)
		 else
		    u=group_id(j)
		    v=group_id(i)
		    do w=1,ntot
		       if (group_id(w).eq.u) group_id(w)=v
		    enddo		    
		 endif
	      enddo
	   enddo	
c-----------------Find the largest group ---------------------------
	   do i=1,ntot
	      u=group_id(i)
	      group_size(u)=group_size(u)+1
	   enddo
	   u=group_size(1)
	   v=1
	   do i=2,current_group_id
	      if (group_size(i).gt.u) then
		 u=group_size(i)
		 v=i
	      endif
	   enddo
	   do i=1,ntot
	      j=id(i)	   
	      if (group_id(i).ne.v) star_para(k,j,5)=-1
	   enddo	   
	enddo
  
	return
	END 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine interpolate_PSF(nsam,npsam,image,posi,ns,npp,PSF_coe)
        implicit none   

	integer ns,npp,nsam,npsam
	double precision posi(npsam,2),PSF_coe(ns,ns,npp)
	real image(npsam,ns,ns)
	
	integer i,j,k,u,v,px,py,order,nn
	double precision matx(ns,ns,npp,npp),bvec(ns,ns,npp)
	double precision vec(npp),matx_tmp(npp,npp),matx_1(npp,npp)
	double precision xx,yy

	
        do i=1,ns
	  do j=1,ns
	    do u=1,npp
              PSF_coe(i,j,u)=0d0		
	      do v=1,npp
	        matx(i,j,u,v)=0d0
	      enddo
	      bvec(i,j,u)=0d0
	    enddo
	  enddo
	enddo	    	

        do i=1,ns
	  do j=1,ns
            do k=1,nsam
	      xx=posi(k,1)
	      yy=posi(k,2)

	      px=0
	      py=-1
	      order=-1
	      do nn=1,npp
	        if (py.eq.order) then
	          order=order+1
	          px=order
	          py=0
	        else
	          px=px-1
	          py=py+1
	        endif
	        vec(nn)=xx**px*yy**py
	      enddo

	      do u=1,npp
	        do v=1,npp
	          matx(i,j,u,v)=matx(i,j,u,v)+vec(u)*vec(v)
	        enddo
	        bvec(i,j,u)=bvec(i,j,u)+image(k,i,j)*vec(u)
	      enddo
	    enddo
	  enddo	    	
	enddo


	do i=1,ns
	  do j=1,ns	  
	    do u=1,npp
	      do v=1,npp
	        matx_tmp(u,v)=matx(i,j,u,v)
	      enddo
	    enddo
	    call matrix_inverse_doub(matx_tmp,npp,npp,matx_1)
	    do u=1,npp
	      do v=1,npp
	        PSF_coe(i,j,u)=PSF_coe(i,j,u)+matx_1(u,v)*bvec(i,j,v)
	      enddo
	    enddo
	  enddo
	enddo


 
        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_PSF_model(ns,np,PSF_coe,xx,yy,model)
	implicit none

	integer ns,np,i,j,nn,px,py,order
	double precision xx,yy,PSF_coe(ns,ns,np)
	real model(ns,ns)

	do i=1,ns
	  do j=1,ns
	    model(i,j)=0.	    
	    px=0
	    py=-1
	    order=-1
	    do nn=1,np
	      if (py.eq.order) then
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      model(i,j)=model(i,j)+xx**px*yy**py*PSF_coe(i,j,nn)
	    enddo
	  enddo
	enddo	    	

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_population_boundary(np,n,dat,bound,ierror)
	implicit none

	integer np,n,ierror
	real dat(np),bound

	real dat_sort(n),sum,sum2,mean,sig,arr(n-1)
	integer i,nmax,indx(n-1),m

	do i=1,n
	  dat_sort(i)=dat(i)
	enddo

	call sort(n,n,dat_sort)

	do i=1,n-1
	   arr(i)=dat_sort(i+1)-dat_sort(i)
c	   write(*,*) i,arr(i)
	enddo

	arr(1)=0.
	call indexx(n-1,n-1,arr,indx)

	m=n-1
	nmax=indx(m)-1
	do while (nmax.le.10)
	  m=m-1
	  nmax=indx(m)-1
	enddo  

c	write(*,*) nmax
	
c	if (nmax.ge.2) then
	  bound=dat_sort(nmax)
	  ierror=0
c	else
c	  ierror=1
c	endif


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_chi2_in(ns,n1,n2,image,model,chi2)
	implicit none

	integer ns,n1,n2
	real image(ns,ns)
	integer i,j
	real model(ns,ns),cc,chi2

	cc=ns/2.+1.

	chi2=0.
	do i=n1,n2
	  do j=n1,n2
	    chi2=chi2+(image(i,j)-model(i,j))**2
	  enddo
	enddo	    		
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_chi2_out(ns,n1,n2,image,model,chi2)
	implicit none

	integer ns,n1,n2
	real image(ns,ns)
	integer i,j
	real model(ns,ns),cc,chi2

	cc=ns/2.+1.

	chi2=0.
	do i=n1,n2
	  do j=n1,n2
	    chi2=chi2
     .+((image(i,j)-model(i,j))*((i-cc)**2+(j-cc)**2))**2
	  enddo
	enddo	    		
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
