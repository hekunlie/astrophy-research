	function chi2_MPI(n,np,x,de,nbin,c)
	use mpi
	implicit none

        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs
	integer ierr

	integer np,n,nbin
	real x(np),de(np),c,temp,chi2_MPI
	integer NMAX
	parameter (NMAX=200)
	real xbin(NMAX)
	common /chi2_mpi_pass/ xbin
	real diff(NMAX),summ(NMAX),diff_t(NMAX),summ_t(NMAX)
	real y,abs_y
	integer ibin,i

	do ibin=1,NMAX
	  diff(ibin)=0.
	  summ(ibin)=0.
	enddo
	
	if (my_id.ne.0) then
	  do i=1,n
	    y=x(i)-c*de(i)
	    abs_y=abs(y)
	    if (abs_y.gt.xbin(nbin)) then
	      ibin=nbin+1
	    else	    
	      ibin=1
	      do while (abs_y.gt.xbin(ibin).and.ibin.lt.nbin)
	        ibin=ibin+1
	      enddo
	    endif
            summ(ibin)=summ(ibin)+1.
            if (y.gt.0) diff(ibin)=diff(ibin)+1.
            if (y.lt.0) diff(ibin)=diff(ibin)-1.
	  enddo	  
	endif

	call MPI_Reduce(diff,diff_t,NMAX,mpi_float,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)
	call MPI_Reduce(summ,summ_t,NMAX,mpi_float,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)

	if (my_id.eq.0) then
	  temp=0.
	  do ibin=1,nbin+1
	    temp=temp+diff_t(ibin)**2/(2.*summ_t(ibin))	
	  enddo	  	  
        endif
  	call MPI_Bcast(temp,1,mpi_float,0,MPI_COMM_WORLD,ierr)
	
	chi2_MPI=temp
	
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine statis_MPI(n,np,x,de,nbin,c,dc)
	use mpi
	implicit none
	
        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs

	integer ierr
	
	integer np,n,nbin
	real x(np),de(np),xx(n)
	real c,dc
	integer NMAX
	parameter (NMAX=200)
	real cbin(2,NMAX)

	integer dn
	real xbin(NMAX),xbin_tmp(NMAX)
	common /chi2_mpi_pass/ xbin

	real c1,c2,cc,temp,chi2_MPI
	real v1,v2,vc,thresh,a(3)
	integer change
	integer i,j,ibin

	real xb,db,x2b

	do i=1,NMAX
	   xbin(i)=0.
	   xbin_tmp(i)=0.
	enddo
	
	if (my_id.ne.0) then
	  dn=n/(nbin+1)
	  do i=1,n
	    xx(i)=abs(x(i))
	  enddo
	  call sort(n,n,xx)		
  	  do i=1,nbin
  	    xbin(i)=xx(dn*i)
	  enddo		
	endif	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_Reduce(xbin,xbin_tmp,NMAX,mpi_float,MPI_SUM,0
     .,MPI_COMM_WORLD,ierr)
	if (my_id.eq.0) then
	  do i=1,NMAX
	    xbin(i)=xbin_tmp(i)/(num_procs-1.)
c	    write(*,*) 'bin:',i,xbin(i)
	  enddo
c	  xbin(nbin)=xbin_tmp(nbin)/(num_procs-1.)
c	  do i=1,nbin-1
c	     xbin(i)=(xbin(nbin)*i)/nbin
c	     write(*,*) 'bin:',i,xbin(i)
c	  enddo   
        endif
  	call MPI_Bcast(xbin,NMAX,mpi_float,0,MPI_COMM_WORLD,ierr)
	
	c=0.
	dc=0.2	
	c1=c-dc
	c2=c+dc	
	change=1	
	thresh=20.

	do while (change.eq.1)
	  change=0
	  cc=(c1+c2)*0.5
	  
	  v1=chi2_MPI(n,np,x,de,nbin,(c1+cc)*0.5)
	  v2=chi2_MPI(n,np,x,de,nbin,(c2+cc)*0.5)
	  vc=chi2_MPI(n,np,x,de,nbin,cc)	
	  if (v1.gt.vc+thresh) then
	    c1=(c1+cc)*0.5
	    change=1
	  endif
	  if (v2.gt.vc+thresh) then
	    c2=(c2+cc)*0.5
	    change=1
	  endif
	enddo

	
	dc=(c2-c1)/(NMAX-1.)

	do i=1,NMAX
	  cbin(1,i)=c1+(i-1.)*dc
	  cbin(2,i)=chi2_MPI(n,np,x,de,nbin,cbin(1,i))
c	  if (my_id.eq.0) write(*,*) i,cbin(1,i),cbin(2,i)
	enddo

	call simple_quadratic_fitting(NMAX,NMAX,cbin,a)

	c=-a(2)/2./a(1)
	dc=1./sqrt(2.*a(1))

c	if (my_id.eq.0) write(*,*) c,dc
	
c	write(*,*) 'answer:',c,dc
c	pause


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function chi2(n,np,x,de,nbin,c)
	implicit none

	integer np,n,nbin
	real x(np),de(np),c,chi2

	integer NMAX
	parameter (NMAX=200)
	real xbin(NMAX)
	common /chi2_pass/ xbin

	real diff(NMAX),summ(NMAX),summ_ave
	common /test_pass/ diff,summ

	real y,abs_y
	integer ibin,i

	do ibin=1,nbin+1
	  diff(ibin)=0.
	  summ(ibin)=0.
	enddo

	do i=1,n
	  y=x(i)-c*de(i)
	  abs_y=abs(y)
	  if (abs_y.gt.xbin(nbin)) then
	    ibin=nbin+1
	  else	    
	    ibin=1
	    do while (abs_y.gt.xbin(ibin).and.ibin.lt.nbin)
	      ibin=ibin+1
	    enddo
	  endif
          summ(ibin)=summ(ibin)+1.
          if (y.gt.0) diff(ibin)=diff(ibin)+1.
          if (y.lt.0) diff(ibin)=diff(ibin)-1.
	enddo	  
	
	chi2=0.

	summ_ave=n/(nbin+1.)

	do ibin=1,nbin+1
	  chi2=chi2+diff(ibin)**2/(2.*summ(ibin))	
	enddo	  	  

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine statis(n,np,x,de,nbin,c,dc,imethod)
	implicit none

	integer np,n,nbin,imethod
	real x(np),de(np),xx(n)
	real c,dc
	integer NMAX
	parameter (NMAX=200)
	real cbin(2,NMAX)

	integer dn
	real xbin(NMAX)
	common /chi2_pass/ xbin

	real c1,c2,cc,temp,chi2
	real v1,v2,vc,thresh,a(3)
	integer change
	integer i,j,ibin

	real diff(NMAX),summ(NMAX)
	common /test_pass/ diff,summ

	real xb,db,x2b

	if (imethod.eq.0) then
	  xb=0.
	  db=0.
	  x2b=0.
	  do i=1,n
	    xb=xb+x(i)
	    db=db+de(i)
	    x2b=x2b+x(i)*x(i)  	
c	    write(*,*) i,x(i),de(i)
	  enddo
c	  write(*,*) n,xb,db,x2b
          c=xb/db
          dc=sqrt(x2b/db/db)
	  return
	endif

	c=0.
	dc=0.2

	dn=n/(nbin+1)
	do i=1,n
	  xx(i)=abs(x(i))
	enddo
	call sort(n,n,xx)		
	do i=1,nbin
	  xbin(i)=xx(dn*i)
	enddo		

	c1=c-dc
	c2=c+dc	
	change=1	
	thresh=20.

	do while (change.eq.1)
	  change=0
	  cc=(c1+c2)*0.5	
	  v1=chi2(n,np,x,de,nbin,(c1+cc)*0.5)
	  v2=chi2(n,np,x,de,nbin,(c2+cc)*0.5)
	  vc=chi2(n,np,x,de,nbin,cc)	
	  if (v1.gt.vc+thresh) then
	    c1=(c1+cc)*0.5
	    change=1
	  endif
	  if (v2.gt.vc+thresh) then
	    c2=(c2+cc)*0.5
	    change=1
	  endif
	enddo

	
	dc=(c2-c1)/(NMAX-1.)

	do i=1,NMAX
	  cbin(1,i)=c1+(i-1.)*dc
	  cbin(2,i)=chi2(n,np,x,de,nbin,cbin(1,i))
c	  write(*,*) i,cbin(1,i),cbin(2,i)
	enddo

	call simple_quadratic_fitting(NMAX,NMAX,cbin,a)

	c=-a(2)/2./a(1)
	dc=1./sqrt(2.*a(1))

c	write(*,*) 'answer:',c,dc
c	pause


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


