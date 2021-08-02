	subroutine get_peak_width_low_side(np,n,arr,p,sig)
	implicit none

	integer np,n
	real arr(np),p,sig,den(n),den2(n),thresh
	integer i,j,ip
	
        call sort(n,np,arr)

        do i=2,n-1
	   den(i)=((arr(i+1)-arr(i-1))/2.)**2
	enddo

	den(1)=(arr(2)-arr(1))**2
	den(n)=(arr(n)-arr(n-1))**2

	do i=1,n
	   den2(i)=0.
	enddo
	
	do i=3,n-2
	   do j=i-2,i+2
	      den2(i)=den2(i)+den(j)
	   enddo
	   den2(i)=1./sqrt(den2(i)/5.)
	enddo

	den2(1)=1./sqrt((den(1)+den(2)+den(3))/3.)
	den2(2)=1./sqrt((den(1)+den(2)+den(3)+den(4))/4.)
	den2(n)=1./sqrt((den(n)+den(n-1)+den(n-2))/3.)
	den2(n-1)=1./sqrt((den(n)+den(n-1)+den(n-2)+den(n-3))/4.)

	ip=1
	p=den2(1)

	do i=2,n
	   if (den2(i).gt.p) then
	      p=den2(i)
	      ip=i
	   endif
	enddo

	thresh=p*0.5
	do i=ip+1,n
	   if (den2(i).lt.thresh) exit
	enddo

	p=arr(ip)
	sig=arr(i)-arr(ip)
	
	return
	end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_peak_width(np,n,arr,p,sig,status,direc)
	implicit none

	integer np,n,status,direc
	real arr(np),p,sig,d
	integer i,j,ip,u,v,b1,b2
	integer nb,smooth
	parameter (nb=200)
	parameter (smooth=3)
	real den(0:nb+1),posi(nb),vol1,vol2,tmp,thresh,vol
	integer mark(0:nb+1),bound1,bound2,num,bb(2,2)
	
	call sort(n,np,arr)

	if (direc.gt.0) n=(n*99)/100
	
	d=(arr(n)-arr(1))/(nb-1.)

	do i=1,nb
	  posi(i)=arr(1)+(i-1.)*d
	  den(i)=0.
	  mark(i)=0
	enddo
        den(0)=0.
	den(nb+1)=0.
	mark(0)=0
	mark(nb+1)=0
	
	do i=1,n
	  tmp=(arr(i)-arr(1))/d+1 
	  ip=int(tmp+0.5)
	  do j=max(ip-4*smooth,1),min(ip+4*smooth,nb)
	    den(j)=den(j)+exp(-0.5*((tmp-j)/smooth)**2)
	  enddo
        enddo

	vol1=0
	vol2=0
        bb(2,1)=0
        bb(2,2)=0
	bb(1,1)=0
	bb(1,2)=0
	
	do i=1,nb
	  if (mark(i).gt.0) cycle 
	  mark(i)=1 
	  if (den(i).gt.den(i-1).and.den(i).gt.den(i+1)) then
   	    bound1=i
	    bound2=i
	    ip=i
	    thresh=den(i)*0.5
	    do while ((den(bound2+1).gt.thresh
     ..or.den(bound2+1).lt.den(bound2)).and.bound2.lt.nb)
	      bound2=bound2+1
	      mark(bound2)=2 
 	      if (den(bound2)*0.5.gt.thresh) then
		thresh=den(bound2)*0.5
		ip=bound2
	      endif
	    enddo
	    do while ((den(bound1-1).gt.thresh
     ..or.den(bound1-1).lt.den(bound1)).and.bound1.gt.1
     ..and.mark(bound1-1).ne.2)
	      bound1=bound1-1
	      if (den(bound1)*0.5.gt.thresh) then
		thresh=den(bound1)*0.5
		ip=bound1
	      endif
	    enddo

	    thresh=thresh*0.5
	    b2=ip
	    do while (b2.lt.bound2.and.den(b2).gt.thresh)
	      b2=b2+1
	    enddo
	    b1=ip
	    do while (b1.gt.bound1.and.den(b1).gt.thresh)
	      b1=b1-1
	    enddo
	    bound1=b1
	    bound2=b2
	    
	    vol=0.
	    do j=bound1,bound2
	       vol=vol+den(j)
	    enddo
	       
	    if (vol.gt.vol1) then
	      vol2=vol1
	      bb(2,1)=bb(1,1)
	      bb(2,2)=bb(1,2)
              vol1=vol
	      bb(1,1)=bound1
	      bb(1,2)=bound2
            elseif (vol.gt.vol2) then
	      vol2=vol
              bb(2,1)=bound1
	      bb(2,2)=bound2
	    endif
	  endif
	enddo

	if (direc.gt.0) then
	  if (bb(1,1).gt.bb(2,1)) then
	     bound1=bb(1,1)
	     bound2=bb(1,2)
	     status=1
	  else
	     bound1=bb(2,1)
	     bound2=bb(2,2)
	     status=2
	  endif
	else
	  if (bb(1,1).gt.bb(2,1)) then
	     bound1=bb(2,1)
	     bound2=bb(2,2)
	     status=2
	  else
	     bound1=bb(1,1)
	     bound2=bb(1,2)
	     status=1
	  endif
	endif
	
	p=0
	sig=0
	num=0
	do i=1,n
	  if (arr(i).lt.posi(bound1).or.arr(i).gt.posi(bound2)) cycle 
	  p=p+arr(i)
	  sig=sig+arr(i)**2
	  num=num+1
	enddo
	p=p/num
	sig=max(sqrt(sig/num-p*p),0.02*(posi(bound2)-posi(bound1)))

	
c	write(*,*) bound1,bound2,num,posi(bound1),posi(bound2)	
c        do i=1,nb
c	  write(*,*) 'den(posi):',i,posi(i),den(i)
c	enddo 
c	write(*,*) bb(1,1),bb(1,2),vol1
c	write(*,*) bb(2,1),bb(2,2),vol2
c	write(*,*) 'final:',p,sig,p-3*sig,p+3*sig,status
	
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine get_med_sig(np,n,dat,med,sig)
	implicit none

	integer n,np
	real dat(np),med,sig,arr(n)
	integer i
	
	do i=1,n
	  arr(i)=dat(i)
	enddo

	call sort(n,n,arr)
	
	med=arr(n/2)
	sig=0.5*(arr(n*5/6)-arr(n/6))

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_mean_sig(np,n,dat,mean,sig)
	implicit none

	integer n,np
	real dat(np),mean,sig
	integer i

	mean=0.
	sig=0.
	do i=1,n
	  mean=mean+dat(i)
	  sig=sig+dat(i)**2
	enddo
	mean=mean/n
	sig=sqrt(sig/n-mean**2)	


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_PREFIX(imagefile,PREFIX)
	implicit none

	character*(*) imagefile,PREFIX
	integer i,p_dot,p_slash,n

	p_dot=0
	p_slash=0

	n=len(trim(imagefile))
	do i=n,1,-1
	  if (imagefile(i:i).eq.'.'.and.p_dot.eq.0) p_dot=i
	  if (imagefile(i:i).eq.'/'.and.p_slash.eq.0) p_slash=i
	enddo
	if (p_dot.eq.0.or.p_slash.eq.0..or.p_dot.le.p_slash+1) 
     .pause 'Image_file name is NOT normal !'		
	
	PREFIX=imagefile(p_slash+1:p_dot-1)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_field_name(DIR_OUTPUT,field)
	implicit none

	character DIR_OUTPUT*(*),field*(*)
	integer i,p_slash,n

	p_slash=0

        n=len(trim(DIR_OUTPUT))
	do i=n,1,-1
	  if (DIR_OUTPUT(i:i).eq.'/'.and.p_slash.eq.0) then
            p_slash=i
            exit
          endif
	enddo
	
	field=DIR_OUTPUT(p_slash+1:LEN_TRIM(DIR_OUTPUT))

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_expo_name(image_file,expo)
	implicit none

	character*(*) image_file,expo
	integer i,p_slash,p_us,n

	p_slash=0
        p_us=0
        n=len(trim(image_file))
	do i=n,1,-1
	  if (image_file(i:i).eq.'_'.and.p_us.eq.0) p_us=i
	  if (image_file(i:i).eq.'/'.and.p_slash.eq.0) p_slash=i
	enddo

	if (p_us.eq.0.or.p_slash.eq.0..or.p_us.le.p_slash+1) 
     .pause 'Image_file name is NOT normal !'		
	
	expo=image_file(p_slash+1:p_us-1)
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine fit_poly_2D_1(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y
	integer np,n
	real arr(np,4),c(3),s_2
	real coe(3,3),coe_1(3,3),vec(3),temp(3)
	real x,y,f
	integer i,j,u,v

	do i=1,3
	  do j=1,3
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  x=arr(i,1)
	  y=arr(i,2)
	  f=arr(i,3)
	  s_2=arr(i,4)

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y

	  do u=1,3
	    do v=1,3
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)*s_2
	    enddo
	    vec(u)=vec(u)+f*temp(u)*s_2		  
	  enddo
    
	enddo

        call matrix_inverse(coe,3,coe_1)

	do j=1,3
	  c(j)=0.
	  do i=1,3
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function yeqx(x)
	implicit none
	real yeqx,x	
	yeqx=x
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function stepf(a,b)
	implicit none

	integer stepf
	real a,b

	if (a.gt.b) then
	  stepf=1
	elseif (a.lt.b) then
	  stepf=-1
	else
	  stepf=0
	endif

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine simple_quadratic_fitting(np,n,x,a)
	implicit none
	
	integer np,n,i,j
	real x(2,np),a(3),u(n),cc,b
	real c(3,4),denomi,temp(3,3),determinant3

c to find the best fitting for x(2,i)=a(1)*x(1,i)**2+a(2)*x(1,i)+a(3)

c to define x(1,i)=b*u(i)+cc

	do i=1,n
          u(i)=x(1,i)
	enddo

	call sort(n,n,u)
	cc=u(n/2)
	b=u(n)-u(1)

	do i=1,n
          u(i)=(x(1,i)-cc)/b
	enddo
			

	do i=1,3
	  do j=1,4
	    c(i,j)=0.
	  enddo
	enddo
	
	do i=1,n
	  c(1,1)=c(1,1)+u(i)**4
	  c(2,1)=c(2,1)+u(i)**3
	  c(3,1)=c(3,1)+u(i)**2

	  c(1,2)=c(1,2)+u(i)**3
	  c(2,2)=c(2,2)+u(i)**2
	  c(3,2)=c(3,2)+u(i)

	  c(1,3)=c(1,3)+u(i)**2
	  c(2,3)=c(2,3)+u(i)
	  c(3,3)=c(3,3)+1.

	  c(1,4)=c(1,4)+u(i)*u(i)*x(2,i)
	  c(2,4)=c(2,4)+u(i)*x(2,i)
	  c(3,4)=c(3,4)+x(2,i)
	enddo

	do i=1,3
	  do j=1,3
	    temp(i,j)=c(i,j)
	  enddo
	enddo

	denomi=determinant3(temp)

	do i=1,3
	  do j=2,3
	    temp(i,j)=c(i,j)
	  enddo
	enddo
	temp(1,1)=c(1,4)
	temp(2,1)=c(2,4)
	temp(3,1)=c(3,4)
	a(1)=determinant3(temp)/denomi

	do i=1,3
	  do j=1,3,2
	    temp(i,j)=c(i,j)
	  enddo
	enddo
	temp(1,2)=c(1,4)
	temp(2,2)=c(2,4)
	temp(3,2)=c(3,4)
	a(2)=determinant3(temp)/denomi

	do i=1,3
	  do j=1,2
	    temp(i,j)=c(i,j)
	  enddo
	enddo
	temp(1,3)=c(1,4)
	temp(2,3)=c(2,4)
	temp(3,3)=c(3,4)
	a(3)=determinant3(temp)/denomi

	a(1)=a(1)/b/b
	a(2)=a(2)/b-2.*a(1)*cc
	a(3)=a(3)-a(1)*cc*cc-a(2)*cc

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function determinant3(c)
	implicit none

	real c(3,3),determinant3
	
	determinant3=c(1,1)*c(2,2)*c(3,3)+c(2,1)*c(3,2)*c(1,3)
     .+c(3,1)*c(1,2)*c(2,3)
	determinant3=determinant3-c(1,3)*c(2,2)*c(3,1)
     .-c(1,2)*c(2,1)*c(3,3)-c(1,1)*c(2,3)*c(3,2)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function quad_func(a,x)
	implicit none

	real quad_func,x,a(3)
	
	quad_func=a(1)*x*x+a(2)*x+a(3)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine matrix_inverse2(ma,n,np,vec)
	implicit none

	integer n,np
	real ma(np,np),vec(np),ma_1(n,n),d,b(n),tmp(n)
	integer indx(n),i,j

        call ludcmp(ma,n,np,indx,d)

	do i=1,n
	  do j=1,n
	    if (j.eq.i) then
	      b(j)=1.
	    else 
	      b(j)=0.
	    endif
	  enddo	
          call lubksb(ma,n,np,indx,b)
	  do j=1,n
            ma_1(j,i)=b(j)
	  enddo	
	enddo
		
	do i=1,n
	  tmp(i)=vec(i)
	enddo

	do i=1,n
	  vec(i)=0
	  do j=1,n
	    vec(i)=vec(i)+ma_1(i,j)*tmp(j)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine matrix_inverse_doub(ma,n,np,ma_1)
	implicit none

	integer n,np
	double precision ma(np,np),ma_1(np,np),d,b(n),norm
	integer indx(n),i,j

	norm=0d0
	do i=1,n
	  do j=1,n
	    norm=max(norm,abs(ma(i,j)))
	  enddo
	enddo
	    
	do i=1,n
	  do j=1,n
	    ma(i,j)=ma(i,j)/norm
	  enddo
	enddo

        call ludcmp_doub(ma,n,np,indx,d)

	do i=1,n
	  do j=1,n
	    if (j.eq.i) then
	      b(j)=1d0
	    else 
	      b(j)=0d0
	    endif
	  enddo	

          call lubksb_doub(ma,n,np,indx,b)

	  do j=1,n
            ma_1(j,i)=b(j)
	  enddo	
	enddo

	do i=1,n
	  do j=1,n
	    ma_1(i,j)=ma_1(i,j)/norm
	  enddo
	enddo
		
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine matrix_inverse(ma,n,ma_1)
	implicit none

	integer n
	real ma(n,n),ma_1(n,n),d,b(n)
	integer indx(n),i,j

        call ludcmp(ma,n,n,indx,d)

	do i=1,n
	  do j=1,n
	    if (j.eq.i) then
	      b(j)=1.
	    else 
	      b(j)=0.
	    endif
	  enddo	

          call lubksb(ma,n,n,indx,b)

	  do j=1,n
            ma_1(j,i)=b(j)
	  enddo	
	enddo
		
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine find_slope_2D(np,n,arr,aa,bb,cc)
	implicit none

	integer np,n
	real arr(np,3),aa,bb,cc
	integer i,j
	real c(3,3),vec(3),c_1(3,3)

	do i=1,3
	  do j=1,3
	    c(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  c(1,1)=c(1,1)+1.
	  c(1,2)=c(1,2)+arr(i,1)
	  c(1,3)=c(1,3)+arr(i,2)
	  vec(1)=vec(1)+arr(i,3)		  

	  c(2,1)=c(2,1)+arr(i,1)
	  c(2,2)=c(2,2)+arr(i,1)*arr(i,1)
	  c(2,3)=c(2,3)+arr(i,2)*arr(i,1)
	  vec(2)=vec(2)+arr(i,3)*arr(i,1)		  

	  c(3,1)=c(3,1)+arr(i,2)
	  c(3,2)=c(3,2)+arr(i,1)*arr(i,2)
	  c(3,3)=c(3,3)+arr(i,2)*arr(i,2)
	  vec(3)=vec(3)+arr(i,3)*arr(i,2)		  
	enddo

        call matrix_inverse(c,3,c_1)

	aa=0.
	bb=0.
	cc=0.
	do i=1,3
	  aa=aa+c_1(1,i)*vec(i)
	  bb=bb+c_1(2,i)*vec(i)
	  cc=cc+c_1(3,i)*vec(i)
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function fit_func(x,y,nc,nx,n)
	implicit none

	real x,y,fit_func
	integer nc,nx,n
	integer px,py

	px=mod(n-1,nx)
	py=(n-1)/nx

	fit_func=x**px*y**py
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function func_val(x,y,nc,nx,c)
	implicit none

	real x,y,func_val
	integer nc,nx
	real c(nc),fit_func
	integer i

	func_val=0.
	do i=1,nc
	   func_val=func_val+fit_func(x,y,nc,nx,i)*c(i)
	enddo
	   
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_2D(np,n,arr,nc,nx,c)
	implicit none

c  fit a 2D function with f(x,y)=\Sigma_{i=1}^{nc}func_i(x,y)c_i
	
	integer np,n,nc,nx
	real arr(np,3),c(nc)
	real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
	real x,y,f,smax
	real fit_func
	external fit_func
	integer i,j,u,v

	do i=1,nc
	  do j=1,nc
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	do i=1,n
	  x=arr(i,1)
	  y=arr(i,2)
	  f=arr(i,3)

	  do j=1,nc
	    temp(j)=fit_func(x,y,nc,nx,j)
	  enddo
	
	  do u=1,nc
	    do v=1,nc
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)
	    enddo
	    vec(u)=vec(u)+f*temp(u)	  
	  enddo
    
	enddo

        call matrix_inverse(coe,nc,coe_1)

	do j=1,nc
	  c(j)=0.
	  do i=1,nc
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_poly_2D_2(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y+c4x^2+c5xy+c6y^2
	integer np,n
	real arr(np,3),c(6)
	real coe(6,6),coe_1(6,6),vec(6),temp(6)
	real x,y,f,smax
	integer i,j,u,v

	do i=1,6
	  do j=1,6
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	smax=0.
	do i=1,n
	  smax=max(smax,abs(arr(i,1)))
	  smax=max(smax,abs(arr(i,2)))
	enddo
	smax=1./smax

	do i=1,n
	  x=arr(i,1)*smax
	  y=arr(i,2)*smax
	  f=arr(i,3)

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y
	  temp(4)=x*x
	  temp(5)=x*y
	  temp(6)=y*y

	  do u=1,6
	    do v=1,6
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)
	    enddo
	    vec(u)=vec(u)+f*temp(u)	  
	  enddo
    
	enddo

        call matrix_inverse(coe,6,coe_1)

	do j=1,6
	  c(j)=0.
	  do i=1,6
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	c(2)=c(2)*smax
	c(3)=c(3)*smax

	c(4)=c(4)*smax*smax
	c(5)=c(5)*smax*smax
	c(6)=c(6)*smax*smax


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_poly_2D_3(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y+c4x^2+c5xy+c6y^2+c7x^3+c8x^2y+c9xy^2+c10y^3
	integer np,n,nc
	parameter (nc=10)
	
	real arr(np,3),c(nc)
	real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
	real x,y,f,smax,smax2,smax3
	integer i,j,u,v

	do i=1,nc
	  do j=1,nc
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	smax=0.
	do i=1,n
	  smax=max(smax,abs(arr(i,1)))
	  smax=max(smax,abs(arr(i,2)))
	enddo
	smax=1./smax

	do i=1,n
	  x=arr(i,1)*smax
	  y=arr(i,2)*smax
	  f=arr(i,3)

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y

	  temp(4)=x*x
	  temp(5)=x*y
	  temp(6)=y*y

	  temp(7)=x**3
	  temp(8)=x*x*y
	  temp(9)=x*y*y
	  temp(10)=y**3

	  do u=1,nc
	    do v=1,nc
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)
	    enddo
	    vec(u)=vec(u)+f*temp(u)	  
	  enddo
    
	enddo

        call matrix_inverse(coe,nc,coe_1)

	do j=1,nc
	  c(j)=0.
	  do i=1,nc
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	c(2)=c(2)*smax
	c(3)=c(3)*smax

	smax2=smax*smax

	c(4)=c(4)*smax2
	c(5)=c(5)*smax2
	c(6)=c(6)*smax2

	smax3=smax2*smax

	c(7)=c(7)*smax3
	c(8)=c(8)*smax3
	c(9)=c(9)*smax3
	c(10)=c(10)*smax3


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_poly_2D_4(np,n,arr,c)
	implicit none

c fit a 2D function with f=c1+c2x+c3y+c4x^2+c5xy+c6y^2+c7x^3+c8x^2y+c9xy^2+c10y^3
c     +c11x^4+c12x^3y+c13x^2y^2+c14xy^3+c15y^4
	integer np,n,nc
	parameter (nc=15)
	
	real arr(np,3),c(nc)
	real coe(nc,nc),coe_1(nc,nc),vec(nc),temp(nc)
	real x,y,f,smax,smax2,smax3,smax4,x2,y2,xy
	integer i,j,u,v

	do i=1,nc
	  do j=1,nc
	    coe(i,j)=0.
	  enddo
	  vec(i)=0.
	enddo

	smax=0.
	do i=1,n
	  smax=max(smax,abs(arr(i,1)))
	  smax=max(smax,abs(arr(i,2)))
	enddo
	smax=1./smax

	do i=1,n
	  x=arr(i,1)*smax
	  y=arr(i,2)*smax
	  f=arr(i,3)

	  x2=x*x
	  y2=y*y
	  xy=x*y

	  temp(1)=1.
	  temp(2)=x
	  temp(3)=y
	  
	  temp(4)=x2
	  temp(5)=xy
	  temp(6)=y2

	  temp(7)=x2*x
	  temp(8)=x*xy
	  temp(9)=xy*y
	  temp(10)=y2*y

	  temp(11)=x2*x2
	  temp(12)=x2*xy
	  temp(13)=xy*xy
	  temp(14)=xy*y2
	  temp(15)=y2*y2

	  do u=1,nc
	    do v=1,nc
	      coe(u,v)=coe(u,v)+temp(u)*temp(v)
	    enddo
	    vec(u)=vec(u)+f*temp(u)	  
	  enddo
    
	enddo

        call matrix_inverse(coe,nc,coe_1)

	do j=1,nc
	  c(j)=0.
	  do i=1,nc
	    c(j)=c(j)+coe_1(j,i)*vec(i)
	  enddo
	enddo

	c(2)=c(2)*smax
	c(3)=c(3)*smax

	smax2=smax*smax

	c(4)=c(4)*smax2
	c(5)=c(5)*smax2
	c(6)=c(6)*smax2

	smax3=smax2*smax

	c(7)=c(7)*smax3
	c(8)=c(8)*smax3
	c(9)=c(9)*smax3
	c(10)=c(10)*smax3

	smax4=smax2*smax2

	c(11)=c(11)*smax4
	c(12)=c(12)*smax4
	c(13)=c(13)*smax4
	c(14)=c(14)*smax4
	c(15)=c(15)*smax4


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_linear_2D(n,np,x,xt,coe)
	implicit none

	integer n,np
	double precision x(np,2),xt(np,2),coe(2,3)
	double precision matx(3,3),bvec1(3),bvec2(3),vec(3)
	double precision matx_1(3,3)
	integer i,j,u,v
	

        do u=1,3
          do v=1,3
	    matx(u,v)=0d0
	  enddo
	  bvec1(u)=0d0
	  bvec2(u)=0d0
        enddo

	do i=1,n
c    xt(i,1)=coe(1,1)*x(i,1)+coe(1,2)*x(i,2)+coe(1,3)
c    xt(i,2)=coe(2,1)*x(i,1)+coe(2,2)*x(i,2)+coe(2,3)

	  vec(1)=x(i,1)
	  vec(2)=x(i,2)
	  vec(3)=1d0	  

	  do u=1,3
            do v=1,3
	      matx(u,v)=matx(u,v)+vec(u)*vec(v)
	    enddo
	    bvec1(u)=bvec1(u)+xt(i,1)*vec(u)
	    bvec2(u)=bvec2(u)+xt(i,2)*vec(u)
	  enddo
	enddo	    	

	call matrix_inverse_doub(matx,3,3,matx_1)

	do i=1,3
	  coe(1,i)=0d0
	  coe(2,i)=0d0
	  do j=1,3
	    coe(1,i)=coe(1,i)+matx_1(i,j)*bvec1(j)
	    coe(2,i)=coe(2,i)+matx_1(i,j)*bvec2(j)
	  enddo
	enddo

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine plot_comparison(np,n,x,y,de,nbin,imethod
     .,func,filename)
	implicit none
	
	integer n,np,nbin,i,j,is,imethod,nb
	real func
	external func
	real x(np),y(np),de(np),yy(np),dd(np)
	character filename*(*),fname*200
	real xmin,xmax,dx,xbin(nbin)
	common /bound_pass/ xmin,xmax

	real arr(nbin,4),c,dc

        fname=trim(filename)//'.dat'
        open(unit=50,file=fname)
        rewind 50

	dx=(xmax-xmin)/nbin

	nb=5

	do i=1,nbin
	  arr(i,1)=xmin+dx*(i-0.5)
	  is=0	
	  do j=1,n
	    if (x(j).ge.xmin+dx*(i-1).and.x(j).le.xmin+dx*i) then
	      is=is+1
	      yy(is)=y(j)
	      dd(is)=de(j)
	    endif
	  enddo
	  call statis(is,np,yy,dd,nb,c,dc,imethod)
          arr(i,2)=c
          arr(i,3)=dc
	  arr(i,4)=is

c          imagename='chi2.fits'
c         call make_plots(is,np,yy,dd,nb,c,dc*4.,imagename)
	  write(*,*) i,(arr(i,j),j=1,4)
	  write(50,*) i,(arr(i,j),j=1,4)
	enddo

	close(50)

        fname=trim(filename)//'.fits'
	call map_function(func,nbin,nbin,arr,500,500,fname)

        return
        END 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine pick_anormaly(n,np,arr,mark)
	implicit none
	
	integer n,np
	real arr(np),tmp(np),med,sig
	integer mark(np),i,num

	num=0
	do i=1,n
	  if (mark(i).eq.1) then
	    num=num+1
	    tmp(num)=arr(i)
	  endif
	enddo

	call sort(num,np,tmp)

	med=tmp(num/2)
	sig=0.5*(tmp(num*5/6)-tmp(num/6))
	
	do i=1,n
	  if (mark(i).eq.1) then
	    if (abs(arr(i)-med).gt.3.*sig) mark(i)=0
	  endif
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine fit_map(n,np,ns,n1,n2,image,posi,npl,coe)
	implicit none

	integer n,np,ns,npl,n1,n2
	real image(np,ns,ns),posi(np,2),coe(npl,ns,ns)
	real mt(npl,npl),bv(npl),x,y,vec(npl)
	integer i,j,k,px,py,order,nn,u,v

	
	do i=n1,n2
	  do j=n1,n2

            do u=1,npl
	      do v=1,npl
	        mt(u,v)=0.
	      enddo
	      bv(u)=0.
	    enddo
	    
	    do k=1,n
	      x=posi(k,1)
	      y=posi(k,2)

	      px=0
	      py=-1
	      order=-1
	      do nn=1,npl
	        if (py.eq.order) then
	          order=order+1
	          px=order
	          py=0
	        else
	          px=px-1
	          py=py+1
	        endif
	        vec(nn)=x**px*y**py
	      enddo

	      do u=1,npl
	        do v=1,npl
	          mt(u,v)=mt(u,v)+vec(u)*vec(v)
	        enddo
	        bv(u)=bv(u)+image(k,i,j)*vec(u)
	      enddo
	    enddo

	    call matrix_inverse2(mt,npl,npl,bv)
	    do u=1,npl
	      coe(u,i,j)=bv(u)
	    enddo

	  enddo
	enddo	    	
	
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_model(ns,n1,n2,posi,npl,coe,model)
	implicit none

	integer ns,npl,n1,n2
	real posi(2),coe(npl,ns,ns)
	integer i,j,px,py,order,nn,u,v
	real model(ns,ns)

	do i=n1,n2
	  do j=n1,n2
	    model(i,j)=0.	    
	    px=0
	    py=-1
	    order=-1
	    do nn=1,npl
	      if (py.eq.order) then
	        order=order+1
	        px=order
	        py=0
	      else
	        px=px-1
	        py=py+1
	      endif
	      model(i,j)=model(i,j)+posi(1)**px*posi(2)**py*coe(nn,i,j)
	    enddo
	  enddo
	enddo	    		
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_chi2(ns,n1,n2,image,model,chi2)
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
        subroutine get_star_property(n,star,p)
        implicit none

        integer n
        real star(n,n),p,peak,thresh,area
        integer i,j
	real pi
	parameter (pi=3.1415926) 
        
	peak=star(1,1)
	do i=1,n
	  do j=1,n
	    if (star(i,j).gt.peak) peak=star(i,j)
	  enddo
	enddo

        thresh=exp(-1.)*peak

        area=0.
        do i=1,n
          do j=1,n
	    if (star(i,j).ge.thresh) area=area+1.
	  enddo
	enddo
	
	p=sqrt(area/pi)
        
        return
        end      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ana_chi2(n,map1,map2,p)
        implicit none

        integer n
        real map1(n,n),map2(n,n),p,flux
        integer i,j,n1,n2

        n1=n/4
        n2=(n/4)*3
        p=0.
        flux=0.
        do i=n1,n2
          do j=n1,n2
            flux=flux+(map1(i,j)+map2(i,j))*0.5
            p=p+(map1(i,j)-map2(i,j))**2
          enddo
        enddo
        p=p/flux

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
