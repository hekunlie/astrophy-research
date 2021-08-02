        subroutine mpi_distribute(num_job,job,message)
        use mpi
        implicit none

        integer my_id,num_procs
        common /MPIpar/ my_id,num_procs

        external job        
        character*(*) message

	integer num_job,i,j,k,complete
	integer source,tag,ierr,status(mpi_status_size)

        
	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

	i=0
	if (my_id.eq.0) then
	   i=1
	   j=num_job
	endif
	complete=0
	do while (complete.eq.0)
	  if (my_id.ne.0) then
	    call MPI_SEND(i,1,mpi_int,0,0,MPI_COMM_WORLD,ierr)
	    call MPI_RECV(i,1,mpi_int,0,0,MPI_COMM_WORLD,status,ierr)
	    if (i.eq.0) then
	      complete=1
	    else
	      call job(i)
	    endif
	  else
	    call MPI_RECV(k,1,mpi_int,MPI_ANY_SOURCE
     .,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
	    tag=status(MPI_TAG)
	    source=status(MPI_SOURCE)
	    call MPI_SEND(i,1,mpi_int,source,tag,MPI_COMM_WORLD,ierr)
	    if (k.gt.0) then
	      j=j-1
	    endif
	    write(*,*) source,i,j,trim(message)
	    if (i.ne.0) i=i+1
	    if (i.gt.num_job) i=0
	    if (j.eq.0) complete=1
	  endif	  
	enddo            

	call MPI_BARRIER(MPI_COMM_WORLD, ierr ) ! synchronize all nodes 

        
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


