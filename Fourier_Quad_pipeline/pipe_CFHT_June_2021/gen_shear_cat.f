	subroutine gen_shear_cat(IMAGE_FILE,nchip
     .,DIR_OUTPUT,DIR_PSF)
 	implicit none
	include 'para.inc'

	integer ichip,nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT,DIR_PSF
	
	if (expo_global.eq.0.or.ext_PSF.eq.1) then
	  do ichip=1,nchip	 
	    call chip_shear(IMAGE_FILE,ichip
     .,DIR_OUTPUT,DIR_PSF)
	  enddo
        else
	  call expo_shear(IMAGE_FILE,nchip,DIR_OUTPUT)	  
	endif

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine chip_shear(IMAGE_FILE,ichip,DIR_OUTPUT,DIR_PSF)
        implicit none   
	include 'para.inc'

	integer ichip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT,DIR_PSF
	integer proc_error,nstar,status
	character*(strl) PREFIX1,PREFIX2,PREFIX_head
	character*(strl) PREFIX,headname,filename

	integer ngal,ierror,nn1,nn2,i,j,u,v,nx,ny,k
	real g1,g2,de,h1,h2
	real gal_p_coll(ngal_max,ns,ns),gal_para(ngal_max,npara)

	double precision CRPIX(2),CD(2,2),CRVAL(2)
	double precision PU(2,npd)

	real gal_p(ns,ns),psf_model(ns,ns),aa(npara)
	double precision ra,dec,gf1,gf2,temp,cos2,sin2,cos4,sin4
	double precision x,y
	integer parity

	double precision local_coe(ns,ns,npl)
	real ePSF(ns,ns),ePSF_p(ns,ns)

	proc_error=0
c----------------------------------------------------------------
	if (ext_PSF.eq.1) then

	  filename=trim(DIR_PSF)//'/PSF.fits'
          call readimage(filename,nx,ny,ns,ns,ePSF)
          call get_power(ns,ns,ePSF,ePSF_p,0)
	  do i=1,ns
	    do j=1,ns
	      local_coe(i,j,1)=ePSF_p(i,j)
	    enddo
	  enddo	

	else
          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
          PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
	  filename=trim(PREFIX)//'_PSF_coe_local.dat'
          open(unit=10,file=filename)
          rewind 10
	  read(10,*) nstar,status
	  if (status.eq.-1) then
	    proc_error=1
	  else
	    do i=1,ns
	      do j=1,ns
	        read(10,*) (local_coe(i,j,k),k=1,npl)
	      enddo
	    enddo
	  endif
	  close(10)

	endif
c----------------------------------------------------------------


	
	call get_PREFIX(IMAGE_FILE(1),PREFIX_head)
	headname=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX_head)//'.head'

        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        PREFIX1=trim(DIR_OUTPUT)//'/stamps/'
     .//trim(PREFIX)
	PREFIX2=trim(DIR_OUTPUT)//'/result/'
     .//trim(PREFIX)

        if (proc_error.eq.1) then
	  ngal=0
	  goto 50
	endif

        call read_astrometry_para(headname,ichip
     .,CRPIX,CD,CRVAL,PU,npd,proc_error)	     

        if (proc_error.eq.1) then
	  ngal=0
	  goto 50
	endif

	ngal=0
	filename=trim(PREFIX1)//'_source_info.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
	  write(*,*) filename
	  pause 'Catalog file error in shear_proc!!'
        endif   
        read(10,*) 
c 	read(10,*) 'ig xc yc sigma peak imax jmax half_light_flux'
c     .,' half_light_area flag flux2 flux_alt' 
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,iclass)
 	  if (ierror.lt.0) cycle
	  ngal=ngal+1
	  do i=1,iclass
	    gal_para(ngal,i)=aa(i)
	  enddo
        enddo
        close(10)

	nn1=ns*len_g
 	nn2=ns*(int(ngal/len_g)+1)
        filename=trim(PREFIX1)//'_source_p.fits'
	call read_stamps(ngal_max,1,ngal,ns,ns
     .,gal_p_coll,nn1,nn2,filename)

50	filename=trim(PREFIX2)//'_shear.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10

 	write(10,*) 'ig xc yc sigma nstar imax jmax '
     .,'half_light_flux half_light_area flag flux2 flux_alt ' 
     .,'ra dec gf1 gf2 g1 g2 de h1 h2' 

	do i=1,ngal

          do u=1,ns
            do v=1,ns
              gal_p(u,v)=gal_p_coll(i,u,v)
            enddo
          enddo	

      	  x=gal_para(i,2)
	  y=gal_para(i,3)

	  if (ext_PSF.eq.1) then
	    call get_PSF_model(ns,1,local_coe,x,y,psf_model)
	  else
	    call get_PSF_model(ns,npl,local_coe,x,y,psf_model)
	  endif
	 
	  gal_para(i,istar)=nstar

	  call coordinate_transfer_PU(ra,dec,x,y,1
     .,CRPIX,CD,CRVAL,PU,npd)

	  gal_para(i,ira)=ra
	  gal_para(i,idec)=dec

	  call field_distortion_PU(x,y,npd,PU,CD,CRPIX
     .,gf1,gf2,cos2,sin2,parity)

	  gal_para(i,igf1)=gf1
	  gal_para(i,igf2)=gf2

	  call get_shear(ns,gal_p,psf_model,g1,g2,de,h1,h2)

	  gal_para(i,igf2+1)=g1*cos2+g2*sin2
	  gal_para(i,igf2+2)=g2*cos2-g1*sin2
	  gal_para(i,igf2+3)=de

	  cos4=cos2*cos2-sin2*sin2
	  sin4=2d0*sin2*cos2

	  gal_para(i,igf2+4)=h1*cos4+h2*sin4
	  gal_para(i,igf2+5)=h2*cos4-h1*sin4
	
	  if (parity.eq.-1) then
	    gal_para(i,igf2+2)=-gal_para(i,igf2+2)
	    gal_para(i,igf2+5)=-gal_para(i,igf2+5)
	  endif

          write(10,*) (gal_para(i,j),j=1,ih2)

        enddo

	close(10)

        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE get_shear(n,gal,psf,g1,g2,de,h1,h2)
        implicit none

	integer n
        real g1,g2,de,h1,h2,gal(n,n),psf(n,n)
        integer i,j,n_2,cc
	real ks,peak,thresh,area,ks_2,kx,kx2,ky,ky2,k2,k,temp,temp1	
	real filter_deriv,filter,ff,norm
	
	real pi,PSFr_ratio
	parameter (pi=3.1415926) 
	parameter (PSFr_ratio=0.75)

	peak=psf(1,1)
	do i=1,n
	  do j=1,n
	    if (psf(i,j).gt.peak) peak=psf(i,j)
	  enddo
	enddo

        thresh=exp(-1.)*peak

        area=0.
        do i=1,n
          do j=1,n
	    if (psf(i,j).ge.thresh) area=area+1.
	  enddo
	enddo
	
	ks=sqrt(area/pi)	 
	ks_2=(ks*PSFr_ratio)**(-2)
	
	thresh=peak*1e-5
	n_2=n/2
	cc=1+n_2

        g1=0.
        g2=0.
        de=0.
	h1=0.
	h2=0.

	do i=1,n
	  kx=i-cc
	  kx2=kx*kx
	  do j=1,n
	    ky=j-cc
	    ky2=ky*ky
	    k2=kx2+ky2
	    k=sqrt(k2)
	    if (psf(i,j).gt.thresh) then
	      ff=k2*ks_2
	      temp=exp(-ff)/psf(i,j)
	      temp1=temp*gal(i,j)

	      g1=g1-temp1*(kx2-ky2)
	      g2=g2-temp1*2.*kx*ky
	      de=de+temp1*k2*(2.-ff)
	      h1=h1+temp1*ks_2*(k2*k2-8.*kx2*ky2)
	      h2=h2+temp1*ks_2*4.*kx*ky*(kx2-ky2)

	    endif
	  enddo
	enddo


        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine expo_shear(IMAGE_FILE,nchip,DIR_OUTPUT)
        implicit none   
	include 'para.inc'

	integer ichip,nchip
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT
	integer proc_error,nstar_tot,status

	character*(strl) PREFIX1,PREFIX2,PREFIX,headname,filename

	integer ngal,ierror,nn1,nn2,i,j,u,v,k
	real g1,g2,de,h1,h2
	real gal_p_coll(ngal_max,ns,ns),gal_para(ngal_max,npara)

	double precision CRPIX(2),CD(2,2),CRVAL(2)
	double precision PU(2,npd)

	real gal_p(ns,ns),psf_model(ns,ns),aa(npara)
	double precision ra,dec,gf1,gf2,temp,cos2,sin2,cos4,sin4
	double precision x,y,xx,yy
	integer parity

	double precision PSF_coe(ns,ns,npo)

	proc_error=0

        call get_PREFIX(IMAGE_FILE(1),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

	filename=trim(PREFIX)//'_PSF_coe.dat'
        open(unit=10,file=filename)
        rewind 10
	read(10,*) nstar_tot,status
	if (status.eq.-1) then
	  proc_error=1
	else
	  do i=1,ns
	    do j=1,ns
	      read(10,*) (PSF_coe(i,j,k),k=1,npo)
	    enddo
	  enddo
	endif
	close(10)
	
	call get_PREFIX(IMAGE_FILE(1),PREFIX)
	headname=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'

	do ichip=1,nchip

          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)

          if (proc_error.eq.1) then
  	    ngal=0
	    goto 50
	  endif

	  ierror=0
          call read_astrometry_para(headname,ichip
     .,CRPIX,CD,CRVAL,PU,npd,ierror)	     

          if (ierror.eq.1) then
	    ngal=0
	    goto 50
	  endif

	  ngal=0
	  filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_source_info.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
	    write(*,*) filename
	    pause 'Catalog file error!!'
          endif   
          read(10,*)
 
c 	  read(10,*) 'ig xc yc sigma peak imax jmax'
c     .,' half_light_flux half_light_area flag flux2 flux_alt' 

          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (aa(i),i=1,iclass)
 	    if (ierror.lt.0) cycle
	    ngal=ngal+1
	    do i=1,iclass
	      gal_para(ngal,i)=aa(i)
	    enddo
          enddo
          close(10)

	  nn1=ns*len_g
 	  nn2=ns*(int(ngal/len_g)+1)
          filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_source_p.fits'
	  call read_stamps(ngal_max,1,ngal,ns,ns
     .,gal_p_coll,nn1,nn2,filename)

50	  filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)
     .//'_shear.dat'
          open(unit=10,file=filename,status='replace')
          rewind 10

 	  write(10,*) 'ig xc yc sigma nstar imax jmax half_light_flux'
     .,' half_light_area flag flux2 flux_alt ' 
     .,'ra dec gf1 gf2 g1 g2 de h1 h2' 

	  do i=1,ngal

            do u=1,ns
              do v=1,ns
                gal_p(u,v)=gal_p_coll(i,u,v)
              enddo
            enddo	

      	    x=gal_para(i,2)
	    y=gal_para(i,3)
            call xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	    call get_PSF_model(ns,npo,PSF_coe,xx,yy,psf_model)
	
	    gal_para(i,istar)=nstar_tot

	    call coordinate_transfer_PU(ra,dec,x,y,1
     .,CRPIX,CD,CRVAL,PU,npd)

	    gal_para(i,ira)=ra
	    gal_para(i,idec)=dec

	    call field_distortion_PU(x,y,npd,PU,CD,CRPIX
     .,gf1,gf2,cos2,sin2,parity)

	    gal_para(i,igf1)=gf1
	    gal_para(i,igf2)=gf2

	    call get_shear(ns,gal_p,psf_model,g1,g2,de,h1,h2)

	    gal_para(i,igf2+1)=g1*cos2+g2*sin2
	    gal_para(i,igf2+2)=g2*cos2-g1*sin2
	    gal_para(i,igf2+3)=de

	    cos4=cos2*cos2-sin2*sin2
	    sin4=2d0*sin2*cos2

	    gal_para(i,igf2+4)=h1*cos4+h2*sin4
	    gal_para(i,igf2+5)=h2*cos4-h1*sin4
	
	    if (parity.eq.-1) then
	      gal_para(i,igf2+2)=-gal_para(i,igf2+2)
	      gal_para(i,igf2+5)=-gal_para(i,igf2+5)
	    endif

            write(10,*) (gal_para(i,j),j=1,ih2)
	
          enddo	
	  close(10)

	enddo	

        return
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


