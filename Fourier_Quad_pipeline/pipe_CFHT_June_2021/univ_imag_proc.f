	subroutine mark_source(n,stamp,weight,sig,source_thresh
     .,core_thresh,boundx,boundy,total_flux,total_area,peak
     .,half_light_flux,half_light_area,flag,radius,xp,yp)
	implicit none

	integer n,flag
	real stamp(n,n),sig,source_thresh,core_thresh,radius,r2max,r_2
	integer weight(n,n),mark(n,n),half_light_area,xp,yp
	real thresh1,thresh2,total_flux,peak,half_light_flux,thresh
	integer cc,r2min,r2,i,j,u,v,changed,ix,iy,total_area
	integer xmin,xmax,ymin,ymax,r2_thresh,s,boundx(2),boundy(2)

	thresh1=core_thresh*sig
	thresh2=source_thresh*sig

	cc=n/2+1
	r2min=n*n*2
	ix=0
	iy=0

	r2_thresh=(n/8)**2

	do i=1,n
	  do j=1,n
	    if (stamp(i,j).ge.thresh1.and.weight(i,j).gt.0) then
	      mark(i,j)=1
	      r2=(i-cc)**2+(j-cc)**2
	      if (r2.lt.r2min) then
	       	r2min=r2
	        ix=i
	        iy=j
	      endif
	    elseif (weight(i,j).eq.0) then
	      mark(i,j)=1
	    else
	      mark(i,j)=0	
	    endif
	  enddo
	enddo

	if (r2min.gt.r2_thresh) then
	  flag=-2
	  return
	endif

	mark(ix,iy)=2
        changed=1
        xmin=ix
        xmax=xmin
	ymin=iy
        ymax=ymin

	total_flux=stamp(ix,iy)
	total_area=1
	peak=stamp(ix,iy)
	xp=ix
	yp=iy

	do while (changed.eq.1)	 
          changed=0
          do i=xmin,xmax
            do j=ymin,ymax
	      if (mark(i,j).eq.2) then
	        do u=max(i-1,1),min(i+1,n)
	          do v=max(j-1,1),min(j+1,n)
	            if (mark(u,v).lt.2.and.stamp(u,v).ge.thresh2
     ..and.weight(u,v).gt.0) then
	              mark(u,v)=2
	              changed=1
	              xmin=min(xmin,u)
	              ymin=min(ymin,v)	
	              xmax=max(xmax,u)
	              ymax=max(ymax,v)	
		      if (stamp(u,v).gt.peak) then
	                peak=stamp(u,v)
	                xp=u
	                yp=v
	              endif
		      total_flux=total_flux+stamp(u,v)
		      total_area=total_area+1		              
	            elseif (weight(u,v).eq.0) then
		      flag=-1
	              return
	 	    endif	 
	          enddo
	        enddo
	      endif
	    enddo
	  enddo
	enddo

	half_light_flux=0.
	half_light_area=0
	thresh=peak*0.5
	
	r2max=0.
	do i=1,n
	  do j=1,n
	    if (mark(i,j).eq.2) then
	      r_2=(i-xp)**2+(j-yp)**2
	      r2max=max(r2max,r_2)
	      if (stamp(i,j).ge.thresh) then
	        half_light_flux=half_light_flux+stamp(i,j)
	        half_light_area=half_light_area+1
	      endif
	    endif
	  enddo
	enddo	     

	radius=sqrt(r2max)
        total_flux=total_flux/sig
        half_light_flux=half_light_flux/sig
	peak=peak/sig

	boundx(1)=xmin
	boundx(2)=xmax
	boundy(1)=ymin
	boundy(2)=ymax

        changed=1
        do while (changed.eq.1)	 
          changed=0	
	  do i=1,n
	    do j=1,n
	      if (mark(i,j).ne.1) cycle
	      mark(i,j)=-1
	      do u=max(i-1,1),min(i+1,n)
	        do v=max(j-1,1),min(j+1,n)
	          if (mark(u,v).eq.0.and.stamp(u,v).ge.thresh2) then
	            mark(u,v)=1
	          endif	 
	        enddo
	      enddo
	      changed=1
	    enddo
	  enddo
	enddo

	do i=1,n
	  do j=1,n
	    if (mark(i,j).ne.-1) cycle
	    mark(i,j)=1
	    do u=max(i-2,1),min(i+2,n)
	      do v=max(j-2,1),min(j+2,n)
	        if (mark(u,v).eq.0) then
	          mark(u,v)=1
	        elseif (mark(u,v).eq.2) then
	          flag=-1
	          return
	        endif	 
	      enddo
	    enddo
          enddo
        enddo

	do i=1,n
	  do j=1,n
	    if (mark(i,j).eq.1) then
	      mark(i,j)=0
	      weight(i,j)=0
	    endif
          enddo
        enddo

	flag=10

	s=2
	do while (s.le.10)
	  do i=1,n
	    do j=1,n
	      if (mark(i,j).eq.s) then
	        do u=max(i-1,1),min(i+1,n)
	          do v=max(j-1,1),min(j+1,n)
	            if (mark(u,v).eq.0) then
	              mark(u,v)=s+1
	              if (weight(u,v).eq.0) then
	                flag=mark(u,v)-2
		        return
	              endif
	            endif
	          enddo
	        enddo
	      endif
            enddo	
          enddo
          s=s+1
	enddo

		

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mark_noise(n,stamp,weight,sig,source_thresh
     .,core_thresh)
	implicit none

	integer n
	real stamp(n,n),sig,source_thresh,core_thresh
	integer mark(n,n),weight(n,n)
	real thresh1,thresh2
	integer i,j,u,v,changed

	thresh1=core_thresh*sig
	thresh2=source_thresh*sig

	do i=1,n
	  do j=1,n
	    if (stamp(i,j).ge.thresh1.and.weight(i,j).gt.0) then
	      mark(i,j)=1
	    elseif (weight(i,j).eq.0) then
	      mark(i,j)=1
	    else
	      mark(i,j)=0	
	    endif
	  enddo
	enddo

        changed=1
        do while (changed.eq.1)	 
          changed=0	
	  do i=1,n
	    do j=1,n
	      if (mark(i,j).ne.1) cycle
	      mark(i,j)=-1
	      do u=max(i-1,1),min(i+1,n)
	        do v=max(j-1,1),min(j+1,n)
	          if (mark(u,v).eq.0.and.stamp(u,v).ge.thresh2) then
	            mark(u,v)=1
	          endif	 
	        enddo
	      enddo
	      changed=1
	    enddo
	  enddo
	enddo

	do i=1,n
	  do j=1,n
	    if (mark(i,j).ne.-1) cycle
	    mark(i,j)=1
	    do u=max(i-2,1),min(i+2,n)
	      do v=max(j-2,1),min(j+2,n)
	        if (mark(u,v).eq.0) mark(u,v)=1
	      enddo
	    enddo
          enddo
        enddo

	do i=1,n
	  do j=1,n
	    if (mark(i,j).eq.1) weight(i,j)=0
          enddo
        enddo


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine flatten_stamp(ns,nl,stamp,weight,ierror)
	implicit none

	integer ns,nl,ierror
	real stamp(nl,nl)
	integer weight(nl,nl)

	real arr(nl*nl,3),aa,bb,cc
	integer i,j,np,d1,d2

	ierror=0
	d1=(nl-ns)/2
	d2=nl-d1+1
	np=0
	do i=1,nl
	  do j=1,nl
	    if ((i.le.d1.or.i.ge.d2.or.j.le.d1.or.j.ge.d2)
     ..and.weight(i,j).eq.1) then
	      np=np+1
	      arr(np,1)=i
	      arr(np,2)=j
	      arr(np,3)=stamp(i,j)	
	    endif
	  enddo
	enddo
	
	if (np.le.(nl*nl-ns*ns)*0.3) then
	  ierror=-1
	  return 
	else  
	  call find_slope_2D(nl*nl,np,arr,aa,bb,cc)
	  do i=1,nl
	    do j=1,nl
              stamp(i,j)=stamp(i,j)-aa-bb*i-cc*j
	    enddo
 	  enddo
	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine decorate_stamp(ns,sig,weights,stamp)
	implicit none

	integer ns
	real stamp(ns,ns),sig,gasdev
	integer weights(ns,ns) 
	integer i,j

	do i=1,ns
	  do j=1,ns
	    if (weights(i,j).eq.0) then
              stamp(i,j)=gasdev()*sig
	    endif
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image33(nx,ny,map)
	implicit none

	integer nx,ny
	real map(nx,ny),temp(nx,ny),f(3,3)
	integer i,j,u,v

	do i=2,nx-1
	  do j=2,ny-1
	    do u=1,3
	      do v=1,3
	        f(u,v)=map(i+u-2,j+v-2)
	      enddo
	    enddo
	    call smooth_grid33(f)
	    temp(i,j)=f(2,2)
	    if (i.eq.2) then
	      temp(1,j)=f(1,2)
	    endif
	    if (i.eq.nx-1) then
	      temp(nx,j)=f(3,2)
	    endif
	    if (j.eq.2) then
	      temp(i,1)=f(2,1)
	    endif
	    if (j.eq.ny-1) then
	      temp(i,ny)=f(2,3)
	    endif
	    if (i.eq.2.and.j.eq.2) then
	      temp(1,1)=f(1,1)
	    endif
	    if (i.eq.2.and.j.eq.ny-1) then
	      temp(1,ny)=f(1,3)
	    endif
	    if (i.eq.nx-1.and.j.eq.2) then
	      temp(nx,1)=f(3,1)
	    endif
	    if (i.eq.nx-1.and.j.eq.ny-1) then
	      temp(nx,ny)=f(3,3)
	    endif
	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=temp(i,j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid33(f)
	implicit none

c  f(x)=a1+a2*x+a3*y+a4*x*2+a5*x*y+a6*y*y


	real f(3,3)
	real vec(6),matx(6,6),x(3),y(3)
	real matx_1(6,6),a(6)

	integer i,j,u,v,k

	integer first
        SAVE first,matx_1
        DATA first /0/
	
        do i=1,3	
	  x(i)=-2+i	
	  y(i)=-2+i	
	enddo

	if (first.eq.0) then
 
	  first=1

	  do u=1,6
	    do v=1,6	
	      matx(u,v)=0.
	    enddo
	  enddo

 	  do i=1,3
	    do j=1,3
	      vec(1)=1.
	      vec(2)=x(i)
	      vec(3)=y(j)
	      vec(4)=x(i)*x(i)
	      vec(5)=x(i)*y(j)
	      vec(6)=y(j)*y(j)
	      do u=1,6
	        do v=1,6	
	          matx(u,v)=matx(u,v)+vec(u)*vec(v)
		enddo
	      enddo
	    enddo
	  enddo
	
   	  call matrix_inverse(matx,6,matx_1)
	endif

	do i=1,6
	  vec(i)=0.
	enddo

	do i=1,3
	  do j=1,3
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,3
	  do j=1,3
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image55(nx,ny,map,ord)
	implicit none

	integer nx,ny,ord
	real map(nx,ny),temp(nx,ny),f(5,5)
	integer i,j,u,v

	do i=3,nx-2
	  do j=3,ny-2
	    do u=1,5
	      do v=1,5
	        f(u,v)=map(i+u-3,j+v-3)
	      enddo
	    enddo

	    if (ord.eq.1) then
	      call smooth_grid55(f)
	    else
	      call smooth_grid55_3rd_order(f)
	    endif

	    temp(i,j)=f(3,3)
	    if (i.eq.3) then
	      temp(1,j)=f(1,3)
	      temp(2,j)=f(2,3)
	    endif
	    if (i.eq.nx-2) then
	      temp(nx,j)=f(5,3)
	      temp(nx-1,j)=f(4,3)
	    endif
	    if (j.eq.3) then
	      temp(i,1)=f(3,1)
	      temp(i,2)=f(3,2)
	    endif
	    if (j.eq.ny-2) then
	      temp(i,ny)=f(3,5)
	      temp(i,ny-1)=f(3,4)
	    endif
	    if (i.eq.3.and.j.eq.3) then
	      temp(1,1)=f(1,1)
	      temp(1,2)=f(1,2)
	      temp(2,1)=f(2,1)
	      temp(2,2)=f(2,2)
	    endif
	    if (i.eq.3.and.j.eq.ny-2) then
	      temp(1,ny)=f(1,5)
	      temp(1,ny-1)=f(1,4)
	      temp(2,ny)=f(2,5)
	      temp(2,ny-1)=f(2,4)
	    endif
	    if (i.eq.nx-2.and.j.eq.3) then
	      temp(nx,1)=f(5,1)
	      temp(nx-1,1)=f(4,1)
	      temp(nx,2)=f(5,2)
	      temp(nx-1,2)=f(4,2)
	    endif
	    if (i.eq.nx-2.and.j.eq.ny-2) then
	      temp(nx,ny)=f(5,5)
	      temp(nx-1,ny)=f(4,5)
	      temp(nx,ny-1)=f(5,4)
	      temp(nx-1,ny-1)=f(4,4)
	    endif
	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=temp(i,j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid55_3rd_order(f)
	implicit none

c  f(x)=a1+a2*x+a3*y+a4*x*2+a5*x*y+a6*y*y+a7*x^3+a8*x^2*y+a9*x*y^2+a10*y^3


	real f(5,5)
	real vec(10),matx(10,10),x(5),y(5)
	real matx_1(10,10),a(10)

	integer i,j,u,v,k

	integer first
        SAVE first,matx_1
        DATA first /0/
	
        do i=1,5	
	  x(i)=-3+i	
	  y(i)=-3+i	
	enddo

	if (first.eq.0) then
 
	  first=1

	  do u=1,10
	    do v=1,10	
	      matx(u,v)=0.
	    enddo
	  enddo

 	  do i=1,5
	    do j=1,5
	      vec(1)=1.
	      vec(2)=x(i)
	      vec(3)=y(j)
	      vec(4)=x(i)*x(i)
	      vec(5)=x(i)*y(j)
	      vec(6)=y(j)*y(j)
	      vec(7)=x(i)*x(i)*x(i)
	      vec(8)=x(i)*x(i)*y(j)
	      vec(9)=x(i)*y(j)*y(j)
	      vec(10)=y(j)*y(j)*y(j)
	      do u=1,10
	        do v=1,10	
	          matx(u,v)=matx(u,v)+vec(u)*vec(v)
		enddo
	      enddo
	    enddo
	  enddo
	
   	  call matrix_inverse(matx,10,matx_1)
	endif

	do i=1,10
	  vec(i)=0.
	enddo

	do i=1,5
	  do j=1,5
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	    vec(7)=vec(7)+f(i,j)*x(i)*x(i)*x(i)	    
	    vec(8)=vec(8)+f(i,j)*x(i)*x(i)*y(j)	    
	    vec(9)=vec(9)+f(i,j)*x(i)*y(j)*y(j)	    
	    vec(10)=vec(10)+f(i,j)*y(j)*y(j)*y(j)	    
	  enddo
	enddo

	do i=1,10
	  a(i)=0.
	  do j=1,10
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,5
	  do j=1,5
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)
     .+a(7)*x(i)*x(i)*x(i)+a(8)*x(i)*x(i)*y(j)	    
     .+a(9)*x(i)*y(j)*y(j)+a(10)*y(j)*y(j)*y(j)	    
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid55(f)
	implicit none

c  f(x)=a1+a2*x+a3*y+a4*x*2+a5*x*y+a6*y*y

	real f(5,5)
	real vec(6),matx(6,6),x(5),y(5)
	real matx_1(6,6),a(6)

	integer i,j

	integer first
        SAVE first,matx_1
        DATA first /0/
	

	if (first.eq.0) then
 
	  first=1
          matx(1,1)=25. 
          matx(1,2)=0.    
          matx(1,3)=0.    
          matx(1,4)=50.    
          matx(1,5)=0.    
          matx(1,6)=50.    

          matx(2,1)=0. 
          matx(2,2)=50.    
          matx(2,3)=0.    
          matx(2,4)=0.    
          matx(2,5)=0.    
          matx(2,6)=0.    

          matx(3,1)=0. 
          matx(3,2)=0.    
          matx(3,3)=50.    
          matx(3,4)=0.    
          matx(3,5)=0.    
          matx(3,6)=0.    

          matx(4,1)=50. 
          matx(4,2)=0.    
          matx(4,3)=0.    
          matx(4,4)=170.    
          matx(4,5)=0.    
          matx(4,6)=100.    

          matx(5,1)=0. 
          matx(5,2)=0.    
          matx(5,3)=0.    
          matx(5,4)=0.    
          matx(5,5)=100.    
          matx(5,6)=0.    

          matx(6,1)=50. 
          matx(6,2)=0.    
          matx(6,3)=0.    
          matx(6,4)=100.    
          matx(6,5)=0.    
          matx(6,6)=170.    
	
   	  call matrix_inverse(matx,6,matx_1)
	endif

	do i=1,5	
	  x(i)=-3+i	
	  y(i)=-3+i	
	enddo

	do i=1,6
	  vec(i)=0.
	enddo

	do i=1,5
	  do j=1,5
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*x(i)	    
	    vec(3)=vec(3)+f(i,j)*y(j)	    
	    vec(4)=vec(4)+f(i,j)*x(i)*x(i)	    
	    vec(5)=vec(5)+f(i,j)*x(i)*y(j)	    
	    vec(6)=vec(6)+f(i,j)*y(j)*y(j)	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1(i,j)*vec(j)
	  enddo
	enddo

	do i=1,5
	  do j=1,5
	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
     .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)	    
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image55_hole(nx,ny,map)
	implicit none

	integer nx,ny
	real map(nx,ny),temp(nx,ny),f(-2:2,-2:2),fc
	integer i,j,u,v,cx,cy,xh,yh
	
	cx=nx/2+1
	cy=ny/2+1

	do i=1,nx
	  do j=1,ny


	    do u=-2,2
	      do v=-2,2
	        f(u,v)
     .=map(mod(i+u-1+nx,nx)+1,mod(j+v-1+ny,ny)+1)
	      enddo
	    enddo

	    xh=cx-i
	    yh=cy-j
	    if (abs(xh).gt.2.or.abs(yh).gt.2) then
	      xh=2
	      yh=2
	    endif
	    call smooth_grid55_with_hole(f,xh,yh,fc)

	    temp(i,j)=fc

	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=temp(i,j)
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_image55_hole_ln(nx,ny,map)
	implicit none

	integer nx,ny
	real map(nx,ny),temp(nx,ny),f(-2:2,-2:2),fc
	integer i,j,u,v,cx,cy,xh,yh
	real delta,mapmax,mapmin
	
	mapmin=map(1,1)
	mapmax=mapmin

	do i=1,nx
	  do j=1,ny
	    mapmin=min(mapmin,map(i,j))
	    mapmax=max(mapmax,map(i,j))
	  enddo
	enddo	
	
	delta=0.0001*(mapmax-mapmin)

	do i=1,nx
	  do j=1,ny
	    map(i,j)=log(map(i,j)+delta)
	  enddo
	enddo	

	cx=nx/2+1
	cy=ny/2+1

	do i=1,nx
	  do j=1,ny


	    do u=-2,2
	      do v=-2,2
	        f(u,v)
     .=map(mod(i+u-1+nx,nx)+1,mod(j+v-1+ny,ny)+1)
	      enddo
	    enddo

	    xh=cx-i
	    yh=cy-j
	    if (abs(xh).gt.2.or.abs(yh).gt.2) then
	      xh=2
	      yh=2
	    endif
	    call smooth_grid55_with_hole(f,xh,yh,fc)

	    temp(i,j)=fc

	  enddo
	enddo
	
	do i=1,nx
	  do j=1,ny
	    map(i,j)=exp(temp(i,j))-delta
	  enddo
	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine smooth_grid55_with_hole(f,xh,yh,fc)
	implicit none

	real fc,f(-2:2,-2:2)
	integer xh,yh

	real vec(6),matx(6,6),matx_1(6,6),a(6)
	integer i,j,u,v,hi,hj

	real matx_1_basis(-2:2,-2:2,6,6)
	integer first
        SAVE first,matx_1_basis
        DATA first /0/
	

	if (first.eq.0) then
 
	  first=1
	

	  do hi=-2,2	
	    do hj=-2,2

              do u=1,6
	        do v=1,6 
                  matx(u,v)=0.
	        enddo
	      enddo

	      do i=-2,2
	        do j=-2,2
	          if (abs(i).eq.2.and.abs(j).eq.2) cycle
	          if (i.eq.hi.and.j.eq.hj) cycle
	          vec(1)=1.
	          vec(2)=i
	          vec(3)=j
	          vec(4)=i*i
	          vec(5)=i*j
	          vec(6)=j*j
	          do u=1,6
	            do v=1,6
	              matx(u,v)=matx(u,v)+vec(u)*vec(v)
	            enddo
	          enddo
	        enddo
	      enddo	
       	      call matrix_inverse(matx,6,matx_1)

              do u=1,6
	        do v=1,6 
                  matx_1_basis(hi,hj,u,v)=matx_1(u,v)
	        enddo
	      enddo

	    enddo
	  enddo

	endif

	do i=1,6
	  vec(i)=0.
	enddo

	do i=-2,2
	  do j=-2,2
            if (abs(i).eq.2.and.abs(j).eq.2) cycle
	    if (i.eq.xh.and.j.eq.yh) cycle
	    vec(1)=vec(1)+f(i,j)	    
	    vec(2)=vec(2)+f(i,j)*i	    
	    vec(3)=vec(3)+f(i,j)*j	    
	    vec(4)=vec(4)+f(i,j)*i*i	    
	    vec(5)=vec(5)+f(i,j)*i*j	    
	    vec(6)=vec(6)+f(i,j)*j*j	
	  enddo
	enddo

	do i=1,6
	  a(i)=0.
	  do j=1,6
	    a(i)=a(i)+matx_1_basis(xh,yh,i,j)*vec(j)
	  enddo
	enddo

	fc=a(1)

c	do i=1,5
c	  do j=1,5
c	    f(i,j)=a(1)+a(2)*x(i)+a(3)*y(j)+a(4)*x(i)*x(i)
c    .+a(5)*x(i)*y(j)+a(6)*y(j)*y(j)	    
c	  enddo
c	enddo
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine process_powers(n,sourcep,noisep)
	implicit none

	integer n,offc,cc
	real sourcep(n,n),noisep(n,n)
	integer i,j
	real temp

	do i=1,n
	  do j=1,n
	    sourcep(i,j)=sourcep(i,j)-noisep(i,j)
	  enddo
	enddo

	temp=0.
	do i=2,n-1
	  temp=temp+sourcep(i,1)+sourcep(i,n)+sourcep(1,i)+sourcep(n,i)
	enddo	
			
	temp=temp/(4.*(n-2.))

	do i=1,n
	  do j=1,n
	    sourcep(i,j)=sourcep(i,j)-temp
	  enddo
	enddo

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine regularize_power(nx,ny,power,star_smooth)
	implicit none

	integer nx,ny,star_smooth
	real power(nx,ny),temp
	integer i,j,cx,cy

	cx=nx/2+1
	cy=ny/2+1

	if (star_smooth.ge.1) then

	  temp=1./power(cx,cy)

	  do i=1,nx
	    do j=1,ny
	      power(i,j)=power(i,j)*temp
	    enddo
	  enddo

	else

	  temp=4./(power(cx+1,cy)+power(cx-1,cy)
     .+power(cx,cy+1)+power(cx,cy-1))

	  do i=1,nx
	    do j=1,ny
	      power(i,j)=power(i,j)*temp
	    enddo
	  enddo
 	  power(cx,cy)=1.

	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_power_area(nx,ny,power,area,thresh_ratio)
	implicit none
	
	integer nx,ny,area
	real thresh_ratio
	real power(nx,ny),thresh,de,rs_2,flx,temp,dx,dy,dx2,dy2,r
	integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
	integer u,v,x,y,tempi

	cx=nx/2+1
	cy=ny/2+1

	thresh=power(cx,cy)*thresh_ratio

	do i=1,nx
	  do j=1,ny
	    mark(i,j)=0
	  enddo
	enddo

	area0=0
	mark(cx,cy)=1
	area=1
	stack(area,1)=cx
	stack(area,2)=cy
	
	do while (area.gt.area0)
	  tempi=area
	  do i=area0+1,tempi
	    x=stack(i,1)
	    y=stack(i,2)
	    do u=max(x-1,1),min(x+1,nx)
	      do v=max(y-1,1),min(y+1,ny)
	        if (mark(u,v).eq.0.and.power(u,v).ge.thresh) then
	          mark(u,v)=1
	          area=area+1
	          stack(area,1)=u
	          stack(area,2)=v
	        endif
	      enddo
	    enddo
	  enddo
	  area0=tempi
	enddo

	area=(area-1)/2

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_power_e(nx,ny,power,e,thresh_ratio)
	implicit none
	
	integer nx,ny,area
	real e(2),norm,thresh_ratio,power(nx,ny),thresh
	integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
	integer u,v,x,y,tempi
	real kx,ky

	cx=nx/2+1
	cy=ny/2+1

	thresh=power(cx,cy)*thresh_ratio

	e(1)=0.
	e(2)=0.
	norm=0.

	do i=1,nx
	  do j=1,ny
	    mark(i,j)=0
	  enddo
	enddo

	area0=0
	mark(cx,cy)=1
	area=1
	stack(area,1)=cx
	stack(area,2)=cy
	
	do while (area.gt.area0)
	  tempi=area
	  do i=area0+1,tempi
	    x=stack(i,1)
	    y=stack(i,2)
	    do u=max(x-1,1),min(x+1,nx)
	      do v=max(y-1,1),min(y+1,ny)
	        if (mark(u,v).eq.0.and.power(u,v).ge.thresh) then
	          mark(u,v)=1
	          area=area+1
	          stack(area,1)=u
	          stack(area,2)=v
	
	          kx=u-cx
	          ky=v-cy
	          e(1)=e(1)+power(u,v)*(kx*kx-ky*ky)
	          e(2)=e(2)+power(u,v)*2.*kx*ky
	          norm=norm+power(u,v)*(kx*kx+ky*ky)

	        endif
	      enddo
	    enddo
	  enddo
	  area0=tempi
	enddo

	e(1)=e(1)/norm
	e(2)=e(2)/norm

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_power_all(nx,ny,power,e,size,thresh_ratio)
	implicit none
	
	integer nx,ny,area
	real e(2),norm,thresh_ratio,power(nx,ny),thresh
	integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
	integer u,v,x,y,tempi
	real kx,ky,size

	cx=nx/2+1
	cy=ny/2+1

	thresh=power(cx,cy)*thresh_ratio

	e(1)=0.
	e(2)=0.
	norm=0.
	size=0.
	
	do i=1,nx
	  do j=1,ny
	    mark(i,j)=0
	  enddo
	enddo

	area0=0
	mark(cx,cy)=1
	area=1
	stack(area,1)=cx
	stack(area,2)=cy
	
	do while (area.gt.area0)
	  tempi=area
	  do i=area0+1,tempi
	    x=stack(i,1)
	    y=stack(i,2)
	    do u=max(x-1,1),min(x+1,nx)
	      do v=max(y-1,1),min(y+1,ny)
	        if (mark(u,v).eq.0.and.power(u,v).ge.thresh) then
	          mark(u,v)=1
	          area=area+1
	          stack(area,1)=u
	          stack(area,2)=v
	
	          kx=u-cx
	          ky=v-cy
	          e(1)=e(1)+power(u,v)*(kx*kx-ky*ky)
	          e(2)=e(2)+power(u,v)*2.*kx*ky
	          norm=norm+power(u,v)*(kx*kx+ky*ky)
		  size=size+power(u,v)
		  
	        endif
	      enddo
	    enddo
	  enddo
	  area0=tempi
	enddo

	e(1)=e(1)/norm
	e(2)=e(2)/norm
	size=norm/size
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_power_shape(nx,ny,power,e,thresh_ratio)
	implicit none
	
	integer nx,ny,area
	real e,e1,e2,norm,thresh_ratio
	real power(nx,ny),thresh,de,rs_2,flx,temp,dx,dy,dx2,dy2,r
	integer i,j,cx,cy,mark(nx,ny),stack(nx*ny,2),area0
	integer u,v,x,y,tempi

	real kx,ky

	cx=nx/2+1
	cy=ny/2+1

	thresh=power(cx,cy)*thresh_ratio

	e1=0.
	e2=0.
	norm=0.

	do i=1,nx
	  do j=1,ny
	    mark(i,j)=0
	  enddo
	enddo

	area0=0
	mark(cx,cy)=1
	area=1
	stack(area,1)=cx
	stack(area,2)=cy
	
	do while (area.gt.area0)
	  tempi=area
	  do i=area0+1,tempi
	    x=stack(i,1)
	    y=stack(i,2)
	    do u=max(x-1,1),min(x+1,nx)
	      do v=max(y-1,1),min(y+1,ny)
	        if (mark(u,v).eq.0.and.power(u,v).ge.thresh) then
	          mark(u,v)=1
	          area=area+1
	          stack(area,1)=u
	          stack(area,2)=v
	
	          kx=u-cx
	          ky=v-cy
c	          e1=e1+power(u,v)*(kx*kx-ky*ky)
c	          e2=e2+power(u,v)*2.*kx*ky
c	          norm=norm+power(u,v)*(kx*kx+ky*ky)
	          e1=e1+(kx*kx-ky*ky)
	          e2=e2+2.*kx*ky
	          norm=norm+(kx*kx+ky*ky)

	        endif
	      enddo
	    enddo
	  enddo
	  area0=tempi
	enddo

	e1=e1/norm
	e2=e2/norm
	e=e1*e1+e2*e2	

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine get_sig_med(nx,ny,image,sig,med)
	implicit none

	integer nx,ny
	real image(nx,ny),sig,med

	integer npp,margin
	parameter (npp=1000)
	parameter (margin=2)

	integer i,ix,iy
	real ran1,pix(npp)

	do i=1,npp
	  ix=int(ran1()*(nx-2.*margin)+margin)
	  iy=int(ran1()*(ny-2.*margin)+margin)
	  pix(i)=image(ix,iy)
	enddo

	call sort(npp,npp,pix)
	
	sig=0.5*(pix(5*npp/6)-pix(npp/6))
	med=pix(npp/2)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine get_entropy(nx,ny,image,sig,med,r,entropy)
	implicit none

	integer nx,ny,r
	real image(nx,ny),sig,med,entropy(nx,ny)
	real a20,a40,a60,a80,p(5),a(5),val
	parameter (a20=0.25)
	parameter (a40=0.52)
	parameter (a60=0.85)
	parameter (a80=1.28)
	
	integer i,j,ii,jj,ix,iy,k

	a(1)=a20*sig
	a(2)=a40*sig
	a(3)=a60*sig
	a(4)=a80*sig
	a(5)=100.*sig


	do i=1,nx
	  do j=1,ny

	    do k=1,5
	      p(k)=0.
	    enddo

	    do ii=i-r,i+r
              ix=mod(ii-1,nx)+1
	      do jj=j-r,j+r
                iy=mod(jj-1,ny)+1
	        val=abs(image(ix,iy)-med)		
   	        k=1
	        do while (val.gt.a(k).and.k.le.4)
	          k=k+1
	        enddo
	        p(k)=p(k)+1.	
	      enddo
	    enddo

	    entropy(i,j)=0.
	    do k=1,5
	      if (p(k).gt.0) entropy(i,j)=entropy(i,j)+p(k)*log(p(k))
	    enddo	  

	  enddo
	enddo
	
	
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	subroutine remove_continuous(nx,ny,npx,npy,map,func,ord)
	implicit none

	integer nx,ny,npx,npy,ord
	real map(npx,npy),maps(nx,ny)
	real func
	external func

	integer i,j

	do i=1,nx
	  do j=1,ny
	    maps(i,j)=func(map(i,j),1)
	  enddo
	enddo

	if (ord.eq.5) then
 	  call smooth_image55(nx,ny,maps,2)
	elseif (ord.eq.4) then
 	  call smooth_image55(nx,ny,maps,1)
	else
	  call smooth_image33(nx,ny,maps)
	endif

	do i=1,nx
	  do j=1,ny
	    maps(i,j)=func(maps(i,j),-1)
	  enddo
	enddo

	do i=1,nx
	  do j=1,ny
	    map(i,j)=map(i,j)-maps(i,j)
	  enddo
	enddo


	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function iden(x,n)
	implicit none

	integer n
	real x,iden

	iden=x

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function loga(x,n)
	implicit none

	integer n
	real x,loga
	
	if (n.eq.1) then
	  if (x.lt.-30.) then
	    loga=-log(-2.*x)
	  else
	    loga=log(x+sqrt(x**2+1.))
	  endif
	else
	  loga=0.5*(exp(x)-exp(-x))
	endif

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine draw_shear_expo(n,map,ncp,nc,pc
     .,nsp,ns,sk,intensity,thickness)
	 implicit none
         
	 integer n,ncp,nc,ns,nsp
         real map(n,n),pc(ncp,4),sk(nsp,5),e(ns),tt(ns)
         real intensity,thickness,xmin,xmax,ymin,ymax
         real ratiox,ratioy,dx,dy,margin,dd,tmp
         real kmin,kmax,emax,ratio,cos2t,cost,sint
         real x1,x2,y1,y2
         parameter (margin=0.05)
         integer i,j

         do i=1,n
           do j=1,n
             map(i,j)=0.
           enddo
         enddo
        
         xmin=1e10
         xmax=-xmin
         ymin=1e10
         ymax=-ymin
         
         do i=1,nc
           xmin=min(pc(i,1),xmin)
           xmin=min(pc(i,3),xmin)
           xmax=max(pc(i,1),xmax)
           xmax=max(pc(i,3),xmax)
           ymin=min(pc(i,2),ymin)
           ymin=min(pc(i,4),ymin)
           ymax=max(pc(i,2),ymax)
           ymax=max(pc(i,4),ymax)
         enddo   

         dx=(xmax-xmin)*margin
         dy=(ymax-ymin)*margin
         xmin=xmin-dx
         xmax=xmax+dx
         ymin=ymin-dy
         ymax=ymax+dy         
         ratiox=(n-1.)/(xmax-xmin)
         ratioy=(n-1.)/(ymax-ymin)

         do i=1,nc
           x1=(pc(i,1)-xmin)*ratiox+1.
           x2=(pc(i,3)-xmin)*ratiox+1.
           y1=(pc(i,2)-ymin)*ratioy+1.
           y2=(pc(i,4)-ymin)*ratioy+1.
           call draw_rectangle(n,n,map,x1,y1
     .,x2,y2,intensity,thickness)
         enddo   


         do i=1,ns
           sk(i,1)=(sk(i,1)-xmin)*ratiox+1.
           sk(i,2)=(sk(i,2)-ymin)*ratioy+1.
         enddo

         tmp=ns
         dd=n/(sqrt(tmp))*0.25

         kmin=10000.
         kmax=-kmin
         emax=0.
         do i=1,ns
           kmin=min(kmin,sk(i,3))
           kmax=max(kmax,sk(i,3))
           e(i)=sqrt(sk(i,4)**2+sk(i,5)**2)
           tt(i)=e(i)
         enddo

         call sort(ns,ns,tt)

         tmp=max(abs(kmax),abs(kmin))*0.01
         ratio=intensity*0.7/(kmax-kmin+tmp)
         
         do i=1,ns
           tmp=(sk(i,3)-kmin)*ratio+1.
           call draw_box_fill(n,n,map,sk(i,1)-dd,sk(i,2)-dd
     .,sk(i,1)+dd,sk(i,2)+dd,tmp)
         enddo
         
         ratio=dd/tt(ns*3/4)
         do i=1,ns
           if (e(i).le.0.) cycle 
           cos2t=sk(i,4)/e(i)
           cost=sqrt((1.+cos2t)*0.5)
           sint=sqrt((1.-cos2t)*0.5)
           if (sk(i,5).gt.0) sint=-sint
           x1=sk(i,1)-e(i)*ratio*cost
           x2=sk(i,1)+e(i)*ratio*cost
           y1=sk(i,2)-e(i)*ratio*sint
           y2=sk(i,2)+e(i)*ratio*sint
           call draw_line(n,n,map,x1,y1,x2,y2,intensity,0.)
         enddo

         call reverse_color(n,n,map)
         
	 return
	 end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine map_function(func,np,n,x,nx,ny,imagename)
        implicit none

c Plot monotonic function & its data

        integer np,n,nx,ny
        real x(np,4),func
	character imagename*(*)
	external func

	real ratio,ymin,ymax,xmin,xmax,dx,dy,x1,x2,y1,y2,xc,yc,sig
	integer nm,i
	parameter (nm=100)
        real map(nx,ny),model(nm,2)
        real yratio,yoffset,xratio,xoffset

        xmax=x(1,1)
	xmin=x(1,1)
        do i=1,n
	  if (x(i,4).gt.0) then
            if (x(i,1).gt.xmax) xmax=x(i,1)
            if (x(i,1).lt.xmin) xmin=x(i,1)
	  endif
        enddo
	dx=xmax-xmin
	  
        ymax=x(1,2)
	ymin=x(1,2)
        do i=1,n
	  if (x(i,4).gt.0) then
c            if (x(i,2)+x(i,3).gt.ymax) ymax=x(i,2)+x(i,3)
c            if (x(i,2)-x(i,3).lt.ymin) ymin=x(i,2)-x(i,3)
            if (x(i,2).gt.ymax) ymax=x(i,2)
            if (x(i,2).lt.ymin) ymin=x(i,2)
	  endif
        enddo
	ymax=max(ymax,func(xmax))
	ymin=min(ymin,func(xmin))

	dy=ymax-ymin

	xmax=xmax+dx*0.1
	xmin=xmin-dx*0.1
	
	ymax=ymax+dy*0.1
	ymin=ymin-dy*0.1

	dx=(xmax-xmin)/(nm-1.)	
	do i=1,nm
	  model(i,1)=xmin+dx*(i-1.)
	  model(i,2)=func(model(i,1))
	enddo

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do i=1,n	
	  if (x(i,4).gt.0) then
            xc=(x(i,1)-xoffset)*xratio+nx/2
            yc=(x(i,2)-yoffset)*yratio+ny/2
            sig=x(i,3)*yratio
            call draw_error(nx,ny,map,xc,yc,sig,180.,1.)
	  endif
        enddo	

        do i=1,nm-1
          x1=(model(i,1)-xoffset)*xratio+nx/2
          x2=(model(i+1,1)-xoffset)*xratio+nx/2
          y1=(model(i,2)-yoffset)*yratio+ny/2
          y2=(model(i+1,2)-yoffset)*yratio+ny/2
          call draw_dashed_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
        enddo	

        x1=(0.-xoffset)*xratio+nx/2
        x2=x1
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(0.-yoffset)*yratio+ny/2
        y2=y1
        call draw_line(nx,ny,map,x1,y1,x2,y2,255.,1.)      
	 
        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  
	   

        return
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine data_binning(np,n,x,a,b,dx,imagename)
	implicit none

	integer np,n,nbin
	real x(np),a,b,dx

	character imagename*(*)
	integer NMAX
	parameter (NMAX=200)
	real num(NMAX),xbin(NMAX)

        real ymin,ymax,xmin,xmax,ratio
	
	integer nx,ny
	parameter (nx=500)
	parameter (ny=500)
        real map(nx,ny)

        real yratio,yoffset,xratio,xoffset
        real intensity,thickness

        integer i,ibin
	real x1,x2,y1,y2

	xmin=a
	xmax=b

	nbin=int((xmax-xmin)/dx+1)
	do ibin=1,nbin+1
	  xbin(ibin)=xmin+dx*(ibin-1.)
	enddo

	ymin=0
	ymax=0
	do ibin=1,nbin
	  num(ibin)=0.
          do i=1,n
	    if (x(i).ge.xbin(ibin).and.x(i).lt.xbin(ibin+1)) then
	      num(ibin)=num(ibin)+1
	    endif 
	    if (x(i).gt.xmax) then
	      pause 'need to increase box size to:'
	      write(*,*) x(i)
	    endif
          enddo
	  ymax=max(ymax,num(ibin))
	enddo

	ymax=ymax*1.2

	xmin=a
	xmax=b

	ratio=0.95	   

        xratio=nx/(xmax-xmin)*ratio
        yratio=ny/(ymax-ymin)*ratio
  	   
        xoffset=(xmax+xmin)*0.5
        yoffset=(ymax+ymin)*0.5
     
        call clean_map(nx,ny,map)
        intensity=100.
        thickness=0.1
       
        x1=(xmin-xoffset)*xratio+nx/2
        x2=(xmax-xoffset)*xratio+nx/2
        y1=(ymin-yoffset)*yratio+ny/2
        y2=(ymax-yoffset)*yratio+ny/2
        call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      

        do ibin=1,nbin
          x1=(xbin(ibin)-xoffset)*xratio+nx/2
          x2=(xbin(ibin+1)-xoffset)*xratio+nx/2
          y1=(ymin-yoffset)*yratio+ny/2
          y2=(num(ibin)-yoffset)*yratio+ny/2
          call draw_rectangle(nx,ny,map,x1,y1,x2,y2,255.,1.)      
        enddo	

        call reverse_color(nx,ny,map)
        call writeimage(imagename,nx,ny,nx,ny,map)	  

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
        subroutine fit_line(np,n,x,a,b)
        implicit none

        integer np,n
        real x(np,2),a,b
	integer i

	real sumy,sumx,sumxx,sumxy

	sumy=0
	sumx=0
	sumxx=0
	sumxy=0

	do i=1,n
	  sumy=sumy+x(i,2)
	  sumx=sumx+x(i,1)
   	  sumxx=sumxx+x(i,1)*x(i,1)
 	  sumxy=sumxy+x(i,2)*x(i,2)
	enddo  	  
	
	b=(n*sumxy-sumx*sumy)/(n*sumxx-sumx*sumx)
	a=(sumy-b*sumx)/n

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fit_line2(nbin,x,y,sig,a,b)
        implicit none

        integer nbin
        real x(nbin),y(nbin),sig(nbin),a,b
	integer i

	real sumy,sumx,sumxx,sumxy,temp,sumn

	sumy=0
	sumx=0
	sumxx=0
	sumxy=0
	sumn=0.

	do i=1,nbin
	  temp=sig(i)**(-2)
	  sumn=sumn+temp
	  sumy=sumy+y(i)*temp
	  sumx=sumx+x(i)*temp
   	  sumxx=sumxx+x(i)*x(i)*temp
 	  sumxy=sumxy+x(i)*y(i)*temp
	enddo  	  
	
	b=(sumn*sumxy-sumx*sumy)/(sumn*sumxx-sumx*sumx)
	a=(sumy-b*sumx)/sumn

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc			
       subroutine clean_map(nx,ny,map)
       implicit none
	   
       integer nx,ny,i,j
       real map(nx,ny)

       do i=1,nx
         do j=1,ny
           map(i,j)=0.
         enddo
       enddo
	   
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 subroutine draw_rectangle(nx,ny,map,x1,y1,x2,y2,intensity
     .,thickness)
	 implicit none
	
	 integer nx,ny
        real map(nx,ny)
        real x1,y1,x2,y2,intensity,thickness

	 call draw_line(nx,ny,map,x1,y1,x1,y2,intensity,thickness)
	 call draw_line(nx,ny,map,x2,y1,x2,y2,intensity,thickness)
	 call draw_line(nx,ny,map,x1,y1,x2,y1,intensity,thickness)
	 call draw_line(nx,ny,map,x1,y2,x2,y2,intensity,thickness)

	 return
	 end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
	 subroutine reverse_color(nx,ny,map)
	 implicit none

	 integer nx,ny,i,j
	 real map(nx,ny)
	 
	 do i=1,nx
	   do j=1,ny
	     map(i,j)=255.-map(i,j)
	   enddo
        enddo
	 
	 return
	 end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
         subroutine draw_error_ln(nx,ny,map,x,y,sig1,sig2
     .,intensity,thickness)
         implicit none
	   
  	 integer nx,ny
         real map(nx,ny)
         real x,y,sig1,sig2,intensity,thickness
	 real x1,x2,y1,y2

	 x1=x
	 x2=x
	 y1=y-sig2
	 y2=y+sig1
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 x1=x-3.*thickness
	 x2=x+3.*thickness
	 y1=y
	 y2=y
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y+sig1
	 y2=y+sig1
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y-sig2
	 y2=y-sig2
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine draw_error(nx,ny,map,x,y,sig,intensity,thickness)
         implicit none
	   
  	 integer nx,ny
         real map(nx,ny)
         real x,y,sig,intensity,thickness
	 real x1,x2,y1,y2

	 x1=x
	 x2=x
	 y1=y-sig
	 y2=y+sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 x1=x-3.*thickness
	 x2=x+3.*thickness
	 y1=y
	 y2=y
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y+sig
	 y2=y+sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)

	 y1=y-sig
	 y2=y-sig
	   
	 call draw_line(nx,ny,map,x1,y1,x2,y2,intensity,thickness)
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 subroutine draw_line(nx,ny,map,x1,y1,x2,y2,intensity
     .,thickness)
         implicit none
	   
  	  integer nx,ny
         real map(nx,ny)
         real x1,y1,x2,y2,intensity,thickness
	   
	 real r,dr,cosx,cosy,x,y	   
         integer ndots,i
	   
         r=sqrt((x2-x1)**2+(y2-y1)**2)	   	  

	   
         if (r.gt.0.) then
	   
           cosx=(x2-x1)/r
           cosy=(y2-y1)/r	
           ndots=int(r)	   
           dr=r/ndots	    
	      
           x=x1
           y=y1
           call draw_dot(nx,ny,map,x,y,intensity,thickness)
           do i=1,ndots
   	     x=x+dr*cosx
   	     y=y+dr*cosy
   	     call draw_dot(nx,ny,map,x,y,intensity,thickness)
 	   enddo	   	  	   
         endif	   
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine draw_dashed_line(nx,ny,map,x1,y1,x2,y2
     .,intensity,thickness)
         implicit none
	   
  	  integer nx,ny
         real map(nx,ny)
         real x1,y1,x2,y2,intensity,thickness
	   
	  real r,dr,cosx,cosy,x,y	   
         integer ndots,i
	   
         r=sqrt((x2-x1)**2+(y2-y1)**2)	   	  

	   
         if (r.gt.0.) then
	   
           cosx=(x2-x1)/r
           cosy=(y2-y1)/r	
           ndots=int(r)	   
           dr=r/ndots	    
	      
           x=x1
           y=y1
           call draw_dot(nx,ny,map,x,y,intensity,thickness)
           do i=1,ndots
   	      x=x+dr*cosx
   	      y=y+dr*cosy
	      if (mod(i,20).le.10) then
     	        call draw_dot(nx,ny,map,x,y,intensity,thickness)
	      endif
 	    enddo	   	  	   
         endif	   
	   
         return
         end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   subroutine draw_dot(nx,ny,map,x,y,intensity,thickness)
	   implicit none
	   
	   integer nx,ny
	   real map(nx,ny)
	   real x,y,intensity,thickness
	   
	   integer ix,iy,i,j
	   
	   ix=int(x+0.5)
	   iy=int(y+0.5)
	   
           do i=int(ix-thickness+0.5),int(ix+thickness+0.5)
             do j=int(iy-thickness+0.5),int(iy+thickness+0.5)
		if (i.ge.1.and.i.le.nx.and.j.ge.1.and.j.le.ny) 
     .map(i,j)=intensity
	     enddo
	   enddo
	   	  
	   
	   return
	   end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         subroutine draw_box_fill(nx,ny,map,x1,y1,x2,y2,intensity)
         implicit none

         integer nx,ny
         real map(nx,ny),x1,x2,y1,y2,intensity

	 integer ix1,ix2,iy1,iy2,i,j
	   
	 ix1=max(1,int(x1+0.5))
         ix2=min(nx,int(x2+0.5))
	 iy1=max(1,int(y1+0.5))
         iy2=min(ny,int(y2+0.5))
         
         do i=ix1,ix2
           do j=iy1,iy2
             map(i,j)=intensity
	   enddo
	 enddo
     
         

         return
         end

