	subroutine chip_process_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
        implicit none   
	include 'para.inc'

	integer nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT
	integer proc_error,ichip
	character PREFIX*(strl),filename*(strl)
	double precision CRPIX(2),CD(2,2),PU(2,npd),CRVAL(2)
	double precision CRPIX_b(2),CD_b(2,2),PU_b(2,npd),CRVAL_b(2)
	integer n,i,j,k
	double precision ra,dec,x,y,ra2,dec2
	
  	if (ASTROMETRY_trivial.eq.1) then	
	  call get_astrometry_trivial(IMAGE_FILE,nchip,DIR_OUTPUT)
        else
	  call get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
c	  if (ext_cat.eq.1) call get_astrometry_gal(IMAGE_FILE,nchip
c     .,DIR_OUTPUT)
        endif

	if (ASTROMETRY_trivial.eq.1) return
	
        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_check.dat' 

 	open(unit=50,file=trim(filename),status='replace')
        rewind 50
	
	do ichip=1,nchip	  
	  proc_error=0
	  
          call get_PREFIX(IMAGE_FILE(1),PREFIX)
          filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head' 

     	  call read_astrometry_para(filename,ichip
     .,CRPIX,CD,CRVAL,PU,npd,proc_error)
	  
	  if (proc_error.eq.0) then

            call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
            filename=trim(DIR_OUTPUT)//'/stamps/'
     .//trim(PREFIX)//'_norm.fits'	
            call update_para(filename,CRPIX,CD)

            filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
    	    open(unit=40,file=trim(filename),status='old')
            rewind 40
	    read(40,*) 
	    read(40,*) 
	    read(40,*) n,j,k
	    do i=1,n
	      read(40,*) ra,dec,x,y
	      call coordinate_transfer_PU(ra2,dec2,x,y,1
     .,CRPIX,CD,CRVAL,PU,npd)
	      write(50,*) ra,dec,ra2,dec2
	    enddo
	    close(40)	    
	  endif
	enddo

	close(50)

c	if (ext_cat.eq.0) return
	
c        call get_PREFIX(IMAGE_FILE(1),PREFIX)
c        filename=trim(DIR_OUTPUT)//'/astrometry/'
c     .//trim(PREFIX)//'_check_gal.dat' 
c 	open(unit=50,file=trim(filename))
c        rewind 50

c        filename=trim(DIR_OUTPUT)//'/astrometry/'
c     .//trim(PREFIX)//'_compare.dat' 
c 	open(unit=60,file=trim(filename))
c        rewind 60
	
c	do ichip=1,nchip	  
c	  proc_error=0
	  
c         call get_PREFIX(IMAGE_FILE(1),PREFIX)
c          filename=trim(DIR_OUTPUT)//'/astrometry/'
c     .//trim(PREFIX)//'_gal.head' 
c     	  call read_astrometry_para(filename,ichip
c     .,CRPIX,CD,CRVAL,PU,npd,proc_error)

c          filename=trim(DIR_OUTPUT)//'/astrometry/'
c     .//trim(PREFIX)//'.head' 
c     	  call read_astrometry_para(filename,ichip
c     .,CRPIX_b,CD_b,CRVAL_b,PU_b,npd,proc_error)
	  
c	  if (proc_error.eq.0) then

c            call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
c            filename=trim(DIR_OUTPUT)//'/astrometry/'
c     .//trim(PREFIX)//'_astro_gal.dat'
c    	    open(unit=40,file=trim(filename))
c            rewind 40
c	    read(40,*) 
c	    read(40,*) 
c	    read(40,*) n,j,k
c	    do i=1,n
c	      read(40,*) ra,dec,x,y
c	      call coordinate_transfer_PU(ra2,dec2,x,y,1
c     .,CRPIX,CD,CRVAL,PU,npd)
c	      write(50,*) ra,dec,ra2,dec2
c	      call coordinate_transfer_PU(ra,dec,x,y,1
c     .,CRPIX_b,CD_b,CRVAL_b,PU_b,npd)
c	      write(60,*) ra,dec,ra2,dec2
c	   enddo
c	    close(40)	    
c	  endif
c	enddo
c
c	close(50)
c	close(60)
	
        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_astrometry_trivial(IMAGE_FILE,nchip,DIR_OUTPUT)
	implicit none
	include 'para.inc'
	
	integer nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT

	character PREFIX*(strl),filename*(strl)
	integer i,valid,ichip,k

	double precision CRPIX2(nchip,2),CD2(nchip,2,2)
	double precision PU(2,npd),CRVAL2(2)

	do ichip=1,nchip

	  call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
	  filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
	   
          open(unit=10,file=filename,status='old')
          rewind 10
	  read(10,*) CRPIX2(ichip,1),CRPIX2(ichip,2)
     .,CRVAL2(1),CRVAL2(2)
	  read(10,*) CD2(ichip,1,1),CD2(ichip,1,2)
     .,CD2(ichip,2,1),CD2(ichip,2,2)
	  close(10)

	enddo

	do i=1,npd
	  PU(1,i)=0d0
	  PU(2,i)=0d0
	enddo

	valid=1

        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) CRVAL2(1),CRVAL2(2)
        do i=1,npd
	  write(10,*) PU(1,i),PU(2,i)
	enddo
	do k=1,nchip
	  write(10,*) k,valid,CRPIX2(k,1),CRPIX2(k,2)
     .,CD2(k,1,1),CD2(k,1,2),CD2(k,2,1),CD2(k,2,2)
        enddo
        close(10)

	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
	implicit none
	include 'para.inc'
	
	integer nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT
	integer ichip

	character*(strl) PREFIX,filename
	integer nss(nchip),n_user,n_ref,i,k
	integer nss_max
	parameter (nss_max=1000000)
	double precision ra2(nchip,nss_max),dec2(nchip,nss_max)
	double precision x2(nchip,nss_max),y2(nchip,nss_max)

	double precision CRPIX2(nchip,2),CD2(nchip,2,2)
	double precision PU(2,npd),CRVAL2(2)
	integer valid(nchip),tot_valid,tot_source
   
	tot_valid=0
	tot_source=0
	
	do ichip=1,nchip

	  call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
	  filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
	   
          open(unit=10,file=filename,status='old')
          rewind 10
	  read(10,*) CRPIX2(ichip,1),CRPIX2(ichip,2)
     .,CRVAL2(1),CRVAL2(2)
	  read(10,*) CD2(ichip,1,1),CD2(ichip,1,2)
     .,CD2(ichip,2,1),CD2(ichip,2,2)
          read(10,*) nss(ichip),n_user,n_ref
	  do i=1,nss(ichip)
	    read(10,*) ra2(ichip,i),dec2(ichip,i)
     .,x2(ichip,i),y2(ichip,i)
	  enddo
	  close(10)
  	  if (nss(ichip).ge.10) then
	    valid(ichip)=1
	    tot_valid=tot_valid+1
	    tot_source=tot_source+nss(ichip)
 	  else
 	    valid(ichip)=0
	  endif
	enddo

	 
	if (tot_source.ge.(npd+tot_valid*3)*3) then
	  call measure_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,CRPIX2,CD2,CRVAL2,PU,npd,valid)
          call check_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,CRPIX2,CD2,CRVAL2,PU,npd,valid)	  
        else
	  do i=1,nchip
	    valid(i)=0
	  enddo
	endif


        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) CRVAL2(1),CRVAL2(2)
        do i=1,npd
	  write(10,*) PU(1,i),PU(2,i)
	enddo
	do k=1,nchip
	  write(10,*) k,valid(k),CRPIX2(k,1),CRPIX2(k,2)
     .,CD2(k,1,1),CD2(k,1,2),CD2(k,2,1),CD2(k,2,2)
        enddo
        close(10)

	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine read_astrometry_para(filename,ichip
     .,CRPIX,CD,CRVAL,PU,npd,proc_error)
	implicit none

	integer ichip,npd,proc_error
	character filename*(*)
	double precision CRPIX(2),CD(2,2),PU(2,npd),CRVAL(2)
	integer valid,j,k,i

	if (proc_error.eq.1) return
	
	open(unit=10,file=filename,status='old')
        rewind 10
        read(10,*) CRVAL(1),CRVAL(2)
  	do i=1,npd
	  read(10,*) PU(1,i),PU(2,i)
	enddo
	do k=1,ichip-1
	   read(10,*)
	enddo
        read(10,*) j,valid,CRPIX(1),CRPIX(2)
     .,CD(1,1),CD(1,2),CD(2,1),CD(2,2)
	
        close(10)

	if (valid.eq.0) proc_error=1
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_astrometry_data_gal(cat_gal,nx,ny,npx,npy
     .,map,sigmap,CRPIX,CD,CRVAL,filename,proc_error,saturation_thresh)
        implicit none   
	
	real astrometry_shift_ratio
	parameter (astrometry_shift_ratio=0.3)

	character*(*) cat_gal,filename
	integer nx,ny,npx,npy,ierror,n_user,n_ref,proc_error
	real map(npx,npy),sigmap(npx,npy),flux_min,saturation_thresh

	integer i,j,u,v,k
	double precision ra,dra,diffra,dec(2),a,d,flux,mag,tmp
	integer nss_max,nss
	parameter (nss_max=1000000)
	double precision xs(nss_max),ys(nss_max)
	double precision xr(nss_max),yr(nss_max)
	double precision ra_r(nss_max),dec_r(nss_max)
	integer box(nss_max),astrometry_shift_range
	double precision ra2(nss_max),dec2(nss_max)
	double precision x2(nss_max),y2(nss_max)

	double precision CRPIX(2),CD(2,2),CRVAL(2)


	if (proc_error.eq.1) then
	   nss=0
	   n_user=0
	   n_ref=0
	  goto 40
	endif

	
c Get the ranges of ra & dec:
	call get_ra_dec_range(nx,ny,ra,dec,dra,CRPIX,CD,CRVAL
     .,astrometry_shift_ratio)

c------------------------------------------------------------------
c Get the positions of the stars in the reference catalog:

	n_ref=0
        open(unit=20,file=trim(cat_gal),status='old',iostat=ierror)
        rewind 20
        if (ierror.ne.0) then
	  write(*,'(A)') cat_gal
	  write(*,*) ierror
	  pause 'Catalog file error!!'
        endif   
	read(20,*)
        do while (ierror.ge.0)
          read(20,*,iostat=ierror) a,d
 	  if (ierror.lt.0) cycle
	  if (abs(diffra(a,ra)).gt.dra*0.5) cycle
	  if (d.lt.dec(1).or.d.gt.dec(2)) cycle
	  n_ref=n_ref+1
	  ra_r(n_ref)=a
	  dec_r(n_ref)=d
	  call coordinate_transfer_simple(a,d,xr(n_ref),yr(n_ref),-1
     .,CRPIX,CD,CRVAL)
        enddo
        close(20)


	call get_astrometry_catalog(nx,ny,npx,npy,map
     .,sigmap,n_ref,nss_max,nss,xs,ys,saturation_thresh*1.5)


	n_user=nss

c Match the point sources in the "user" and "ref" catalogs:
	astrometry_shift_range=int(max(nx,ny)*astrometry_shift_ratio)

	call pattern_matching(nss_max,n_ref,xr,yr,nss_max,n_user,xs,ys
     .,astrometry_shift_range,box)

c Get the Astrometric Calibration parameters:


	nss=0	
        do i=1,n_ref
	  if (box(i).eq.0) cycle
	  nss=nss+1
	  ra2(nss)=ra_r(i)
	  dec2(nss)=dec_r(i)
	  j=box(i)
	  x2(nss)=xs(j)
	  y2(nss)=ys(j) 
        enddo	

	write(*,*) nss,n_ref,n_user,trim(filename)

 40	open(unit=10,file=trim(filename),status='replace')
        rewind 10
	write(10,*) CRPIX(1),CRPIX(2),CRVAL(1),CRVAL(2)
	write(10,*) CD(1,1),CD(1,2),CD(2,1),CD(2,2)
	write(10,*) nss,n_user,n_ref
	do i=1,nss
	  write(10,*) ra2(i),dec2(i),x2(i),y2(i)
	enddo
	close(10)


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine get_astrometry_gal(IMAGE_FILE,nchip,DIR_OUTPUT)
	implicit none
	include 'para.inc'
	
	integer nchip,npchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT
	integer ichip

	character*(strl) PREFIX,filename
	integer nss(nchip),n_user,n_ref,i,k
	integer nss_max
	parameter (nss_max=1000000)
	double precision ra2(nchip,nss_max),dec2(nchip,nss_max)
	double precision x2(nchip,nss_max),y2(nchip,nss_max)

	double precision CRPIX2(nchip,2),CD2(nchip,2,2)
	double precision PU(2,npd),CRVAL2(2)
	integer valid(nchip),tot_valid,tot_source
   
	tot_valid=0
	tot_source=0
	
	do ichip=1,nchip

	  call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
	  filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro_gal.dat'
	   
          open(unit=10,file=filename,status='old')
          rewind 10
	  read(10,*) CRPIX2(ichip,1),CRPIX2(ichip,2)
     .,CRVAL2(1),CRVAL2(2)
	  read(10,*) CD2(ichip,1,1),CD2(ichip,1,2)
     .,CD2(ichip,2,1),CD2(ichip,2,2)
          read(10,*) nss(ichip),n_user,n_ref
	  do i=1,nss(ichip)
	    read(10,*) ra2(ichip,i),dec2(ichip,i)
     .,x2(ichip,i),y2(ichip,i)
	  enddo
	  close(10)
  	  if (nss(ichip).ge.10) then
	    valid(ichip)=1
	    tot_valid=tot_valid+1
	    tot_source=tot_source+nss(ichip)
 	  else
 	    valid(ichip)=0
	  endif
	enddo

	 
	if (tot_source.ge.(npd+tot_valid*3)*3) then
	  call measure_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,CRPIX2,CD2,CRVAL2,PU,npd,valid)
        else
	  do i=1,nchip
	    valid(i)=0
	  enddo
	endif


        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_gal.head'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) CRVAL2(1),CRVAL2(2)
        do i=1,npd
	  write(10,*) PU(1,i),PU(2,i)
	enddo
	do k=1,nchip
	  write(10,*) k,valid(k),CRPIX2(k,1),CRPIX2(k,2)
     .,CD2(k,1,1),CD2(k,1,2),CD2(k,2,1),CD2(k,2,2)
        enddo
        close(10)

	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_astrometry_data_trivial(CRPIX,CD,CRVAL,filename)
        implicit none   
	
	character filename*(*)
	double precision CRPIX(2),CD(2,2),CRVAL(2)
	
	open(unit=10,file=trim(filename),status='replace')
        rewind 10
	write(10,*) CRPIX(1),CRPIX(2),CRVAL(1),CRVAL(2)
	write(10,*) CD(1,1),CD(1,2),CD(2,1),CD(2,2)
	close(10)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gen_astrometry_data(cat_standard,nx,ny,npx,npy
     .,map,sigmap,CRPIX,CD,CRVAL,filename,proc_error,saturation_thresh)
        implicit none   
	
	real astrometry_shift_ratio
	parameter (astrometry_shift_ratio=0.2)

	character cat_standard*(*),filename*(*)
	integer nx,ny,npx,npy,ierror,n_user,n_ref,proc_error
	real map(npx,npy),sigmap(npx,npy),flux_min,saturation_thresh

	integer i,j,u,v,k
	double precision ra,dra,diffra,dec(2),a,d,flux,mag
	integer nss_max,nss,n_user_max
	parameter (nss_max=10000)
	parameter (n_user_max=200)
	
	double precision xs(nss_max),ys(nss_max)
	double precision xr(nss_max),yr(nss_max)
	double precision ra_r(nss_max),dec_r(nss_max)
	integer box(nss_max),astrometry_shift_range
	double precision ra2(nss_max),dec2(nss_max)
	double precision x2(nss_max),y2(nss_max)

	double precision CRPIX(2),CD(2,2),CRVAL(2)


	if (proc_error.eq.1) then
	   nss=0
	   n_user=0
	   n_ref=0
	  goto 40
	endif

	
c Get the ranges of ra & dec:
	call get_ra_dec_range(nx,ny,ra,dec,dra,CRPIX,CD,CRVAL
     .,astrometry_shift_ratio)

c------------------------------------------------------------------
c Get the positions of the stars in the reference catalog:

	n_ref=0
        open(unit=20,file=trim(cat_standard)
     .,status='old',iostat=ierror)
        rewind 20
        if (ierror.ne.0) then
	  write(*,'(A)') cat_standard
	  write(*,*) ierror
	  pause 'Catalog file error!!'
        endif   
	read(20,*)
        do while (ierror.ge.0)
          read(20,*,iostat=ierror) a,d
 	  if (ierror.lt.0) cycle
	  if (abs(diffra(a,ra)).gt.dra*0.5) cycle
	  if (d.lt.dec(1).or.d.gt.dec(2)) cycle
	  n_ref=n_ref+1
	  ra_r(n_ref)=a
	  dec_r(n_ref)=d
	  call coordinate_transfer_simple(a,d,xr(n_ref),yr(n_ref),-1
     .,CRPIX,CD,CRVAL)
        enddo
        close(20)


	call get_astrometry_catalog(nx,ny,npx,npy,map
     .,sigmap,n_user_max,nss_max,nss,xs,ys,saturation_thresh*1.5)


	n_user=nss

c Match the point sources in the "user" and "ref" catalogs:
	astrometry_shift_range=int(max(nx,ny)*astrometry_shift_ratio)

	call pattern_matching(nss_max,n_ref,xr,yr,nss_max,n_user,xs,ys
     .,astrometry_shift_range,box)

c Get the Astrometric Calibration parameters:


	nss=0	
        do i=1,n_ref
	  if (box(i).eq.0) cycle
	  nss=nss+1
	  ra2(nss)=ra_r(i)
	  dec2(nss)=dec_r(i)
	  j=box(i)
	  x2(nss)=xs(j)
	  y2(nss)=ys(j) 
        enddo	

	write(*,*) nss,n_ref,n_user,trim(filename)

 40	open(unit=10,file=trim(filename),status='replace')
        rewind 10
	write(10,*) CRPIX(1),CRPIX(2),CRVAL(1),CRVAL(2)
	write(10,*) CD(1,1),CD(1,2),CD(2,1),CD(2,2)
	write(10,*) nss,n_user,n_ref
	do i=1,nss
	  write(10,*) ra2(i),dec2(i),x2(i),y2(i)
	enddo
	close(10)


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine get_astrometry_catalog(nx,ny,npx,npy,image
     .,sigmap,n_user_max,ns_max,ns,xs,ys,saturation)
	implicit none

	integer nx,ny,npx,npy,n_user_max
	real image(npx,npy),sigmap(npx,npy),flux,saturation
	integer mark(nx,ny)

	integer sizemax
	parameter (sizemax=10000)
	integer nb,nbb,buffer(sizemax,2)
	integer i,j,ix,iy,jx,jy,u,v,k1,k2,k,xp,yp
	real temp,peak,xc,yc

	integer ns_max,ns,toobig,nsb
	double precision xs(ns_max),ys(ns_max)
	double precision xsb(ns_max),ysb(ns_max)

	real flux_array(ns_max)
	integer order(ns_max)

	integer area,area_limit
	parameter (area_limit=400)


	do i=1,nx
	  do j=1,ny
	    if (image(i,j).ge.sigmap(i,j)*5.) then
              mark(i,j)=1
	    else
	      mark(i,j)=0
	    endif
	  enddo
	enddo

	nsb=0

	do i=1,nx
	  do j=1,ny
	    if (mark(i,j).eq.1) then
	      nbb=0
	      nb=1
	      buffer(nb,1)=i
	      buffer(nb,2)=j
	      mark(i,j)=-1
	      xp=i
	      yp=j
	      peak=image(i,j)
	      flux=0.
	      toobig=0
	      do while (nb.gt.nbb)
	        k1=nbb+1
	        k2=nb
	        nbb=nb
	        do k=k1,k2
	          ix=buffer(k,1)
	          iy=buffer(k,2)
	          do u=max(ix-3,1),min(ix+3,nx)
	            do v=max(iy-3,1),min(iy+3,ny)
	              if (mark(u,v).eq.1) then
	                nb=nb+1
	                buffer(nb,1)=u
	                buffer(nb,2)=v
	                mark(u,v)=-1
                        flux=flux+image(u,v)
	                if (image(u,v).gt.peak) then	
	                  peak=image(u,v)
	                  xp=u
	                  yp=v
	                endif
	                if (nb.eq.sizemax) then
                          toobig=1
	                  goto 20
	                endif
		      elseif (mark(u,v).eq.2) then
                        toobig=1
			goto 20
	              endif
	            enddo
	          enddo	                 
	        enddo  	        	   
	      enddo
20            if (toobig.eq.1.or.peak.gt.saturation) then
	        do k=1,nb
	          ix=buffer(k,1)
	          iy=buffer(k,2)
		  mark(ix,iy)=2
		enddo
	      else
	        area=0
	        flux=0.
	        xc=0.
	        yc=0.
	        temp=peak*0.5
                do k=1,nb
	          ix=buffer(k,1)
	          iy=buffer(k,2)
	          if (image(ix,iy).ge.temp) then
	            area=area+1
	            flux=flux+image(ix,iy)
	            xc=xc+image(ix,iy)*ix
	            yc=yc+image(ix,iy)*iy
	          endif
	        enddo
	        xc=xc/flux
	        yc=yc/flux

	        if (area.le.area_limit) then
	          nsb=nsb+1
	          xsb(nsb)=xc
	          ysb(nsb)=yc
		  flux_array(nsb)=flux
		  order(nsb)=nsb
                  if (nsb.eq.ns_max) then
	            ns=0
	            return
	          endif
	        endif
	      endif
 	    endif
	  enddo
	enddo	

	call indexx(nsb,ns_max,flux_array,order)

	ns=0
	i=nsb
	do while (ns.lt.n_user_max.and.i.ge.1)
	  ns=ns+1
	  xs(ns)=xsb(order(i))
	  ys(ns)=ysb(order(i))
c	  write(*,*) ns,flux_array(order(i)),xs(ns),ys(ns)
	  i=i-1
	enddo
c	pause

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine check_astrometry_global(np,n,nC,ra,dec,x,y
     .,CRPIX,CD,CRVAL,PU,npd,valid)
	implicit none

	integer np,nC,npd,iC
	integer valid(nC),n(nC),i
	double precision x(nC,np),y(nC,np),ra(nC,np),dec(nC,np)
	double precision xx(np),yy(np),a(np),d(np),ra_c,dra,dec_c,ddec
	double precision CRPIX(nC,2),CD(nC,2,2),CRVAL(2),PU(2,npd)
	double precision tolerate_shift,aa,dd,crp(2),cdd(2,2),diffra
	parameter (tolerate_shift=0.1d0)

	do iC=1,nC

	  if (valid(iC).eq.0) cycle 

	  do i=1,n(iC)
	    a(i)=ra(iC,i)
	    d(i)=dec(iC,i)
	  enddo
	 
	  call get_ra_dec_bound(np,n(iC),a,d,ra_c,dra,dec_c,ddec)

	  crp(1)=CRPIX(iC,1)
	  crp(2)=CRPIX(iC,2)
	  
	  cdd(1,1)=CD(iC,1,1)
	  cdd(1,2)=CD(iC,1,2)
	  cdd(2,1)=CD(iC,2,1)
	  cdd(2,2)=CD(iC,2,2)

	  do i=1,n(iC)
	     
	    call coordinate_transfer_PU(aa,dd,x(iC,i),y(iC,i)
     .,1,crp,cdd,CRVAL,PU,npd)
	    if (abs(diffra(aa,ra_c)).gt.dra*(0.5+tolerate_shift) 
     ..or.abs(dd-dec_c).gt.ddec*(0.5+tolerate_shift)
     ..or.isnan(aa).or.isnan(dd)) then
	      valid(iC)=0
  	      exit
	    endif
	  enddo

	enddo
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine measure_astrometry_global(np,n,nC,ra,dec,x,y
     .,CRPIX,CD,CRVAL,PU,npd,valid)
	implicit none

	integer np,nC,npd,iC
	integer valid(nC),n(nC)
	double precision x(nC,np),y(nC,np),ra(nC,np),dec(nC,np)
	double precision xi(nC,np),eta(nC,np)
	double precision CRPIX(nC,2),CD(nC,2,2),CRVAL(2),PU(2,npd)

	integer i,j,k,l,u,v,nn,ntot,order,px,py


	double precision vec(npd+nC*3),vec1(npd+nC*3),vec2(npd+nC*3)
	double precision bvec1(npd+nC*3),bvec2(npd+nC*3)
	double precision matx(npd+nC*3,npd+nC*3)
	double precision matx_1(npd+nC*3,npd+nC*3),r
	

	ntot=npd+nC*3

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
	  do i=1,n(k)
	    call ra_dec_to_xi_eta(ra(k,i),dec(k,i),xi(k,i),eta(k,i)
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
	      vec(nn)=xi(k,i)**px*eta(k,i)**py
	    enddo



	    j=npd+(iC-1)*3
	    vec(j+1)=x(k,i)
	    vec(j+2)=y(k,i)
	    vec(j+3)=1d0

	    j=npd+iC*3

	    do u=1,j
	      do v=1,j
	        matx(u,v)=matx(u,v)+vec(u)*vec(v)
	      enddo
	      bvec1(u)=bvec1(u)+xi(k,i)*vec(u)
	      bvec2(u)=bvec2(u)+eta(k,i)*vec(u)
	    enddo

	  enddo
	enddo	    	

        k=npd+iC*3

	call matrix_inverse_doub(matx,k,ntot,matx_1)
	do i=1,k
	  vec1(i)=0d0
	  vec2(i)=0d0
	  do j=1,k
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
c	    if (mod(order,2).eq.1) then
c              nn=nn+1
c	      PU(1,nn)=vec1(nn)
c	      PU(2,nn)=vec2(nn)
c	      if (nn.eq.npd) cycle
c	    endif
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


	iC=0
	do k=1,nC
	  if (valid(k).eq.0) cycle
	  iC=iC+1
	  j=npd+(iC-1)*3
	  CD(k,1,1)=vec1(j+1)
	  CD(k,1,2)=vec1(j+2)

	  CD(k,2,1)=vec2(j+1)
	  CD(k,2,2)=vec2(j+2)

c	  CD(k,1,1)*CRPIX(k,1)+CD(k,1,2)*CRPIX(k,2)=-vec1(j+3)
c	  CD(k,2,1)*CRPIX(k,1)+CD(k,2,2)*CRPIX(k,2)=-vec2(j+3)
	
	  CRPIX(k,1)=-(vec1(j+3)*CD(k,2,2)-vec2(j+3)*CD(k,1,2))
     ./(CD(k,1,1)*CD(k,2,2)-CD(k,2,1)*CD(k,1,2))
	  CRPIX(k,2)=-(vec2(j+3)*CD(k,1,1)-vec1(j+3)*CD(k,2,1))
     ./(CD(k,1,1)*CD(k,2,2)-CD(k,2,1)*CD(k,1,2))

c	  write(*,*) k,CD(k,1,1),CD(k,1,2),CD(k,2,1),CD(k,2,2)
c     .,CRPIX(k,1),CRPIX(k,2)
	enddo
c	pause

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
