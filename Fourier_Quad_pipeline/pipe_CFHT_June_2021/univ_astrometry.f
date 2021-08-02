	subroutine pattern_matching(np0,n0,x0,y0,np1,n1,x1,y1
     .,shift_range,box_final)
	implicit none

	integer np0,np1,n0,n1,shift_range,n_match
	double precision x0(np0),y0(np0),x1(np1),y1(np1)
	
	integer radius,r_fof 
	parameter (radius=40)
	parameter (r_fof=1)

	integer box(n0,0:n1),i,j,k
	integer box_final(np0),u,v,dx,dy,xp,yp

	real mark1(-shift_range:shift_range,-shift_range:shift_range)
	real mark2(-shift_range:shift_range,-shift_range:shift_range)
	integer mark3(-shift_range:shift_range,-shift_range:shift_range)

	integer nn,dpl,changed,xmin,xmax,ymin,ymax
	real thresh,peak,temp,norm

	real map(2*shift_range+1,2*shift_range+1)
	real arr((2*shift_range+1)*(2*shift_range+1))

c	character filename*100
	double precision const
	double precision xx(np0,2),xxt(np0,2),coe(2,3)

	integer crowdy
	parameter (crowdy=0)

	nn=shift_range
	do i=-nn,nn
	  do j=-nn,nn
	    mark1(i,j)=0.
	    mark2(i,j)=0.
	  enddo
	enddo	

        do i=1,n0
	  box(i,0)=0
	  do j=1,n1
	    dx=int(x1(j)-x0(i)+0.5)
	    dy=int(y1(j)-y0(i)+0.5)
            if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
              mark1(dx,dy)=mark1(dx,dy)+1.
	      box(i,0)=box(i,0)+1
	      box(i,box(i,0))=j
	    endif
	  enddo
        enddo

	do i=-nn,nn
	  do j=-nn,nn
	    mark2(i,j)=mark1(i,j)
	  enddo
	enddo	

c	do i=1,2*nn+1
c	  do j=1,2*nn+1
c	    map(i,j)=mark1(i-nn-1,j-nn-1)
c	  enddo
c	enddo
c	filename='./pattern.fits'
c	call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)	
c	pause


	if (radius.gt.0) then
	  const=(radius*0.25d0)**(-2)
  	  do i=-nn,nn
  	    do j=-nn,nn
	      if (mark1(i,j).gt.0) then
	    do u=max(i-radius,-shift_range),min(i+radius,shift_range)
	      do v=max(j-radius,-shift_range),min(j+radius,shift_range)
	        mark2(u,v)=mark2(u,v)
     .+mark1(i,j)/(((i-u)**2+(j-v)**2)*const+1.)
	      enddo
	    enddo
 	      endif
	    enddo
	  enddo	

c	  do i=1,2*nn+1
c	    do j=1,2*nn+1
c	      map(i,j)=mark2(i-nn-1,j-nn-1)
c	    enddo
c	  enddo
c	  filename='./pattern_smoothed.fits'
c	  call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)	
	endif

	peak=0.
	xp=0
	yp=0
	u=0
	do i=-nn+radius,nn-radius
	  do j=-nn+radius,nn-radius
	    u=u+1
	    arr(u)=mark2(i,j)
	    if (mark2(i,j).gt.peak) then
	      peak=mark2(i,j)
	      xp=i
	      yp=j
	    endif  
	  enddo
	enddo	
	
	call sort(u,(2*nn+1)**2,arr)

c	thresh=int(arr(u/2))*3+1

	thresh=(arr(u/2)+peak)*0.5

	changed=1
	xmin=xp
	xmax=xp
	ymin=yp
	ymax=yp
	mark3(xp,yp)=1
	do while (changed.eq.1)
	  changed=0
	  do i=xmin,xmax
	    do j=ymin,ymax
	      if (mark3(i,j).eq.1) then
	        do u=i-r_fof,i+r_fof
	          do v=j-r_fof,j+r_fof
	            if (mark3(u,v).eq.0.and.mark2(u,v).gt.thresh) then
	              mark3(u,v)=1
	              xmin=min(xmin,u)
	              xmax=max(xmax,u)
	              ymin=min(ymin,v)
	              ymax=max(ymax,v)
	              changed=1
	            endif
	          enddo
	        enddo
	      endif
	    enddo
	  enddo
	enddo


        do i=1,n0
          box_final(i)=0	    
	  do j=1,box(i,0)
	    k=box(i,j)

	    dx=int(x1(k)-x0(i)+0.5)
	    dy=int(y1(k)-y0(i)+0.5)
            if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
	      if (mark3(dx,dy).eq.1) then
	        if (box_final(i).eq.0) then
	          box_final(i)=k
	        else
	          box_final(i)=0
	          goto 30
	        endif
	      endif
	    endif
	  enddo
30      enddo

	n_match=0
        do i=1,n0
	  if (box_final(i).eq.0) cycle
	  dpl=0
	  do j=i+1,n0
	    if (box_final(j).eq.0) cycle
	    if (box_final(i).eq.box_final(j)) then
	      box_final(j)=0
	      dpl=1
	    endif
	  enddo	
	  if (dpl.eq.1) then
	    box_final(i)=0
	    cycle
	  endif
	  n_match=n_match+1
	  xx(n_match,1)=x0(i)
	  xx(n_match,2)=y0(i)
	  xxt(n_match,1)=x1(box_final(i))
	  xxt(n_match,2)=y1(box_final(i))
	  
	enddo

	if (n_match.lt.6) then
          do i=1,n0
	    box_final(i)=0	  
	  enddo
	  n_match=0
	  return
	endif

	call fit_linear_2D(n_match,np0,xx,xxt,coe)

	do i=-nn,nn
	  do j=-nn,nn
	    mark1(i,j)=0.
	    mark3(i,j)=0
	  enddo
	enddo	

        do i=1,n0
	  do j=1,n1
	    dx=int(x1(j)-(x0(i)*coe(1,1)+y0(i)*coe(1,2)+coe(1,3))+0.5)
	    dy=int(y1(j)-(x0(i)*coe(2,1)+y0(i)*coe(2,2)+coe(2,3))+0.5)
            if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
              mark1(dx,dy)=mark1(dx,dy)+1.
	    endif
	  enddo
        enddo

c	do i=1,2*nn+1
c	  do j=1,2*nn+1
c	    map(i,j)=mark1(i-nn-1,j-nn-1)
c	  enddo
c	enddo
c	filename='./pattern_corrected.fits'
c	call writeimage(filename,2*nn+1,2*nn+1,2*nn+1,2*nn+1,map)	
c	pause

	if (crowdy.eq.1) then
c--------------------------------------------------------------
c    FOR HIGH SOURCE DENSITY:
	  u=0
	  do i=-nn,nn
	    do j=-nn,nn
	      u=u+1
	      arr(u)=mark1(i,j)
	    enddo
	  enddo	
	  call sort(u,u,arr)

	  peak=mark1(0,0)
	  thresh=(arr(u/2)+peak)*0.5
c---------------------------------------------------------------
	else
c--------------------------------------------------------------
c    FOR LOW SOURCE DENSITY:
	  thresh=0.5
c--------------------------------------------------------------
	endif

	changed=1
	xmin=0
	xmax=0
	ymin=0
	ymax=0
	mark3(0,0)=1
	do while (changed.eq.1)
	  changed=0
	  do i=xmin,xmax
	    do j=ymin,ymax
	      if (mark3(i,j).eq.1) then
	        do u=i-r_fof,i+r_fof
	          do v=j-r_fof,j+r_fof
	            if (mark3(u,v).eq.0.and.mark1(u,v).gt.thresh) then
	              mark3(u,v)=1
	              xmin=min(xmin,u)
	              xmax=max(xmax,u)
	              ymin=min(ymin,v)
	              ymax=max(ymax,v)
	              changed=1
	            endif
	          enddo
	        enddo
	      endif
	    enddo
	  enddo
	enddo
cc--------------------------------------------------------------------------
        do i=1,n0
          box_final(i)=0	    
	enddo
	n_match=0
        do i=1,n0
	  do j=1,n1
	    dx=int(x1(j)-(x0(i)*coe(1,1)+y0(i)*coe(1,2)+coe(1,3))+0.5)
	    dy=int(y1(j)-(x0(i)*coe(2,1)+y0(i)*coe(2,2)+coe(2,3))+0.5)
            if (dx.ge.-nn.and.dx.le.nn.and.dy.ge.-nn.and.dy.le.nn) then
              if (mark3(dx,dy).eq.1) then
	        n_match=n_match+1
	        box_final(i)=j
	        goto 50
	      endif
	    endif
	  enddo
50      enddo


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine mapping_PU(xx,yy,xi,eta,npd,PU,direc)
        implicit none

	integer npd,direc
	double precision xi,eta,xx,yy,PU(2,npd),dxi,deta,r

	integer n,px,py,order,i



	if (direc.eq.1) then
	  xi=xx
	  eta=yy
c         r=sqrt(xi**2+eta**2)

	  do i=1,3
	    dxi=0d0
	    deta=0d0

	    px=0
	    py=1
	    order=1
	    n=0
	    do while (n.lt.npd)
	      if (py.eq.order) then
c	        if (mod(order,2).eq.1) then
c	    	  n=n+1
c	          dxi=dxi+PU(1,n)*r**order
c	          deta=deta+PU(2,n)*r**order
c	          if (n.eq.npd) cycle
c	        endif
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      n=n+1
	      dxi=dxi+PU(1,n)*xi**px*eta**py
	      deta=deta+PU(2,n)*eta**px*xi**py
	    enddo

  	    xi=xx+dxi
	    eta=yy+deta
c            r=sqrt(xi**2+eta**2)
	    
	  enddo

	else

	  xx=xi
	  yy=eta
c         r=sqrt(xi**2+eta**2)

	  px=0
	  py=1
	  order=1
	  n=0
	  do while (n.lt.npd)
	    if (py.eq.order) then
c	      if (mod(order,2).eq.1) then
c		n=n+1
c	        xx=xx-PU(1,n)*r**order
c               yy=yy-PU(2,n)*r**order
c	        if (n.eq.npd) cycle
c	      endif
	      order=order+1
	      px=order
	      py=0
	    else
	      px=px-1
	      py=py+1
	    endif
	    n=n+1
	    xx=xx-PU(1,n)*xi**px*eta**py
            yy=yy-PU(2,n)*eta**px*xi**py
	  enddo

	endif

        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine ra_dec_to_xi_eta(ra,dec,xi,eta,CRVAL1,CRVAL2)
	implicit none
	
	double precision ra,dec,xi,eta,CRVAL1,CRVAL2,diffra

	double precision da,dd,const1,cosda,tandc,tandd,tanda,x,y
	double precision pi
	parameter (pi=3.1415926d0)

	const1=pi/180d0
	tandc=tan(CRVAL2*const1)

        da=diffra(ra,CRVAL1)*const1
        dd=(dec-CRVAL2)*const1
	tandd=tan(dd)
	cosda=cos(da)
	tanda=tan(da)

        y=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
	y=y/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
	x=tanda*(cos(CRVAL2*const1)*(1.-y*tandc))

	xi=x/const1
	eta=y/const1

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function diffra(ra1,ra2)
	implicit none

	double precision ra1,ra2,diffra

	diffra=ra1-ra2
	if (diffra.lt.-180d0) then
	  diffra=diffra+360d0
	elseif (diffra.gt.180d0) then
	  diffra=diffra-360d0
	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function sumra(dra,ra)
	implicit none

	double precision dra,ra,sumra

	sumra=ra+dra
	if (sumra.ge.360d0) then
	  sumra=sumra-360d0
	elseif (sumra.lt.0d0) then
	  sumra=sumra+360d0
	endif

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine xy_to_xxyy(x,y,xx,yy,CRPIX,CD)
	implicit none
	
	double precision x,y,xx,yy
	double precision CRPIX(2),CD(2,2)

	xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
	yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine coordinate_transfer_PU(a,d,x,y,direc
     .,CRPIX,CD,CRVAL,PU,npd)
	implicit none
	
	integer direc,npd
	double precision a,d,x,y,xx,yy,xxx,yyy,rr,dxx,dyy
	double precision CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)
	double precision PU(2,npd)

	double precision da,dd,const1,const2,ds,xi,eta,diffra
	double precision cosda,tandc,tandd,tanda,temp,sumra
	double precision pi
	parameter (pi=3.1415926d0)

c	integer tmp_sig
c	common /temp_pass/ tmp_sig

	const1=pi/180d0
	tandc=tan(CRVAL(2)*const1)

	if (direc.eq.1) then
	  xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
	  yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	  call mapping_PU(xx,yy,xi,eta,npd,PU,1)

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
 
	else
	
	  da=diffra(a,CRVAL(1))*const1
	  dd=(d-CRVAL(2))*const1
	  tandd=tan(dd)
	  cosda=cos(da)
	  tanda=tan(da)

	  yy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
	  yy=yy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
	  xx=tanda*(cos(CRVAL(2)*const1)*(1.-yy*tandc))

	  xi=xx/const1
	  eta=yy/const1

c	  if (tmp_sig.eq.1.and.direc.eq.-1) then
c	    write(*,*) da/const1,dd/const1,xi,eta
c	  endif
  
	  call mapping_PU(xx,yy,xi,eta,npd,PU,2)

	  temp=CD(1,1)*CD(2,2)-CD(1,2)*CD(2,1)
	  temp=1./temp

 	  CD_1(1,1)=CD(2,2)*temp
 	  CD_1(2,2)=CD(1,1)*temp
	  CD_1(1,2)=-CD(1,2)*temp
	  CD_1(2,1)=-CD(2,1)*temp

	  x=xx*CD_1(1,1)+yy*CD_1(1,2)+CRPIX(1)
	  y=xx*CD_1(2,1)+yy*CD_1(2,2)+CRPIX(2)

	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine field_distortion_PU(x,y,npd
     .,PU,CD,CRPIX,g1,g2,cos2,sin2,parity)
	implicit none
	
	integer npd
	double precision xi,eta,x,y,xx,yy,g1,g2,cos2,sin2,aa,bb
	double precision CRPIX(2),CD(2,2),PU(2,npd),cos1,sin1
	integer parity

	double precision mat(2,2),dm(2,2),temp,det,sqrt_det
	double precision dx_dxx,dx_dyy,dy_dxx,dy_dyy
	double precision dxx_dxi,dxx_deta,dyy_dxi,dyy_deta
	double precision r,dr_dxi,dr_deta
	integer px,py,n,order


	xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
        yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	call mapping_PU(xx,yy,xi,eta,npd,PU,1)

	temp=1d0/(CD(1,1)*CD(2,2)-CD(1,2)*CD(2,1))

	dx_dxx=CD(2,2)*temp
	dx_dyy=-CD(1,2)*temp
	dy_dxx=-CD(2,1)*temp
	dy_dyy=CD(1,1)*temp

	dxx_dxi=1d0
	dyy_deta=1d0
	dxx_deta=0d0
	dyy_dxi=0d0

c        r=sqrt(xi**2+eta**2)
c	if (r.gt.0d0) then 
c	  dr_dxi=xi/r
c	  dr_deta=eta/r
c	else  
c	  dr_dxi=0d0
c	  dr_deta=0d0
c	endif

	px=0
	py=1
	order=1
	n=0
	do while (n.lt.npd)
	  if (py.eq.order) then
c	    if (mod(order,2).eq.1) then
c              n=n+1
c              dxx_dxi=dxx_dxi-PU(1,n)*order*r**(order-1)*dr_dxi
c              dxx_deta=dxx_deta-PU(1,n)*order*r**(order-1)*dr_deta
c              dyy_deta=dyy_deta-PU(2,n)*order*r**(order-1)*dr_deta
c	      dyy_dxi=dyy_dxi-PU(2,n)*order*r**(order-1)*dr_dxi
c	      if (n.eq.npd) cycle
c	    endif
	    order=order+1
	    px=order
	    py=0
	  else
	    px=px-1
	    py=py+1
	  endif
	  n=n+1
          dxx_dxi=dxx_dxi-PU(1,n)*px*xi**(px-1)*eta**py
          dxx_deta=dxx_deta-PU(1,n)*xi**px*py*eta**(py-1)
          dyy_deta=dyy_deta-PU(2,n)*px*eta**(px-1)*xi**py
	  dyy_dxi=dyy_dxi-PU(2,n)*eta**px*py*xi**(py-1)

	enddo


	mat(1,1)=dx_dxx*dxx_dxi+dx_dyy*dyy_dxi
	mat(1,2)=dx_dxx*dxx_deta+dx_dyy*dyy_deta
	mat(2,1)=dy_dxx*dxx_dxi+dy_dyy*dyy_dxi
	mat(2,2)=dy_dxx*dxx_deta+dy_dyy*dyy_deta

	det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
	temp=1d0/det

	dm(1,1)=mat(2,2)*temp
	dm(1,2)=-mat(1,2)*temp
	dm(2,1)=-mat(2,1)*temp
	dm(2,2)=mat(1,1)*temp

	parity=1
	if (det.lt.0) then
	  dm(1,1)=-dm(1,1)
	  dm(1,2)=-dm(1,2)
	  parity=-1
	endif

	sqrt_det=sqrt(dm(1,1)*dm(2,2)-dm(1,2)*dm(2,1))

	dm(1,1)=dm(1,1)/sqrt_det
	dm(1,2)=dm(1,2)/sqrt_det
	dm(2,1)=dm(2,1)/sqrt_det
	dm(2,2)=dm(2,2)/sqrt_det

	cos1=0.5d0*(dm(1,1)+dm(2,2))
	sin1=0.5d0*(dm(1,2)-dm(2,1))

	aa=-0.5d0*(dm(1,2)+dm(2,1))
	bb=0.5d0*(dm(2,2)-dm(1,1))

	g1=aa*sin1+bb*cos1
	g2=aa*cos1-bb*sin1

	cos2=cos1*cos1-sin1*sin1
	sin2=2d0*sin1*cos1	

	if (parity.eq.-1) g2=-g2
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_ra_dec_bound(np,n,a,d,ra,dra,dec,ddec)
	implicit none

	integer np,n	
	double precision ra,dec,dra,ddec
	double precision a(np),d(np),ra1,ra2,dec1,dec2
	double precision tmp(n)
	integer i

	dec1=d(1)
	dec2=d(1)
	
	do i=1,n	   
	  dec1=min(dec1,d(i))
	  dec2=max(dec2,d(i))
	enddo
	dec=(dec1+dec2)*0.5d0	
	ddec=dec2-dec1

	ra1=a(1)
	ra2=a(1)
	do i=1,n
	  ra1=min(ra1,a(i))
	  ra2=max(ra2,a(i))
	enddo

	ra=0.5*(ra1+ra2)
	dra=ra2-ra1

	if (dra.gt.180d0) then

	  do i=1,n	      
	    if (a(i).lt.180d0) then
	      tmp(i)=a(i)+360d0
	    else
	      tmp(i)=a(i)
	    endif
	  enddo

	  ra1=tmp(1)
	  ra2=tmp(1)
	  do i=1,n
	    ra1=min(ra1,tmp(i))
	    ra2=max(ra2,tmp(i))
	  enddo

	  dra=ra2-ra1
	  ra=0.5*(ra1+ra2)
	  if (ra.ge.360d0) ra=ra-360d0

	endif			

	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_ra_dec_range_fine(nx,ny,ra,dec,dra
     .,CRPIX,CD,CRVAL,PU,npd,astrometry_shift_ratio)
	implicit none

	integer nx,ny,npd
	double precision ra,dec(2),dra,dec_c,ddec
	double precision x,y,a(4),d(4)

	double precision CRPIX(2),CD(2,2),CRVAL(2)
	double precision PU(2,npd)
	
	real astrometry_shift_ratio

	x=1d0
	y=1d0
	call coordinate_transfer_PU(a(1),d(1),x,y,1,CRPIX,CD,CRVAL,PU,npd)
	x=nx
	y=1d0	
	call coordinate_transfer_PU(a(2),d(2),x,y,1,CRPIX,CD,CRVAL,PU,npd)
	x=nx
	y=ny	
	call coordinate_transfer_PU(a(3),d(3),x,y,1,CRPIX,CD,CRVAL,PU,npd)
	x=1d0
	y=ny	
	call coordinate_transfer_PU(a(4),d(4),x,y,1,CRPIX,CD,CRVAL,PU,npd)

	call get_ra_dec_bound(4,4,a,d,ra,dra,dec_c,ddec)	

	dec(1)=dec_c-0.5*ddec
	dec(2)=dec_c+0.5*ddec
	
	ddec=ddec*astrometry_shift_ratio

	dec(1)=dec(1)-ddec
	dec(2)=dec(2)+ddec

	dra=dra*(1d0+2d0*astrometry_shift_ratio)

	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_ra_dec_range(nx,ny,ra,dec,dra
     .,CRPIX,CD,CRVAL,astrometry_shift_ratio)
	implicit none

	integer nx,ny
	double precision ra,dec(2),dra,dec_c,ddec
	double precision x,y,a(4),d(4)

	double precision CRPIX(2),CD(2,2),CRVAL(2)
	real astrometry_shift_ratio

	x=1d0
	y=1d0	
	call coordinate_transfer_simple(a(1),d(1),x,y,1,CRPIX,CD,CRVAL)
	x=nx
	y=1d0	
	call coordinate_transfer_simple(a(2),d(2),x,y,1,CRPIX,CD,CRVAL)
	x=nx
	y=ny	
	call coordinate_transfer_simple(a(3),d(3),x,y,1,CRPIX,CD,CRVAL)
	x=1d0
	y=ny	
	call coordinate_transfer_simple(a(4),d(4),x,y,1,CRPIX,CD,CRVAL)

	call get_ra_dec_bound(4,4,a,d,ra,dra,dec_c,ddec)	

	dec(1)=dec_c-0.5*ddec
	dec(2)=dec_c+0.5*ddec
	
	ddec=ddec*astrometry_shift_ratio

	dec(1)=dec(1)-ddec
	dec(2)=dec(2)+ddec

	dra=dra*(1d0+2d0*astrometry_shift_ratio)
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine coordinate_transfer_simple(a,d,x,y,direc
     .,CRPIX,CD,CRVAL)
	implicit none
	
	integer direc
	double precision a,d,x,y,xx,yy,xxx,yyy
	double precision CRPIX(2),CD(2,2),CD_1(2,2),CRVAL(2)

	double precision da,dd,const1,const2,ds,diffra,sumra
	double precision cosda,tandc,tandd,tanda,temp
	double precision pi
	parameter (pi=3.1415926d0)

	const1=pi/180d0
	tandc=tan(CRVAL(2)*const1)

	if (direc.eq.1) then
	  xx=CD(1,1)*(x-CRPIX(1))+CD(1,2)*(y-CRPIX(2))
	  yy=CD(2,1)*(x-CRPIX(1))+CD(2,2)*(y-CRPIX(2))	

	  xxx=xx*const1
	  yyy=yy*const1

	  da=xxx/(cos(CRVAL(2)*const1)*(1.-yyy*tandc))
	  da=da-da*da*da*0.33333333333
	  a=sumra(da/const1,CRVAL(1))
	  cosda=1.-da*da*0.5+da**4/24.

	  dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)	  
	  dd=dd-dd*dd*dd*0.3333333333	  	 
	  d=dd/const1+CRVAL(2)
 
	else
	
	  da=diffra(a,CRVAL(1))*const1
	  dd=(d-CRVAL(2))*const1
	  tandd=tan(dd)
	  cosda=cos(da)
	  tanda=tan(da)

	  yy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
	  yy=yy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
	  xx=tanda*(cos(CRVAL(2)*const1)*(1.-yy*tandc))

	  xx=xx/const1
	  yy=yy/const1

	  temp=CD(1,1)*CD(2,2)-CD(1,2)*CD(2,1)
	  temp=1./temp

 	  CD_1(1,1)=CD(2,2)*temp
 	  CD_1(2,2)=CD(1,1)*temp
	  CD_1(1,2)=-CD(1,2)*temp
	  CD_1(2,1)=-CD(2,1)*temp

	  x=xx*CD_1(1,1)+yy*CD_1(1,2)+CRPIX(1)
	  y=xx*CD_1(2,1)+yy*CD_1(2,2)+CRPIX(2)

	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

