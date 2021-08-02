        subroutine map_shear_catalog()
        implicit none
        include 'para.inc'

	character*(strl) IMAGE_FILE(NMAX_EXPO,NMAX_CHIP)
	character*(strl) DIR_OUTPUT(NMAX_EXPO),ASTROMETRY_CAT(NMAX_EXPO)
	character*(strl) SOURCE_CAT(NMAX_EXPO),DIR_PSF(NMAX_EXPO)	
	common /filename_pass/ IMAGE_FILE,DIR_OUTPUT,DIR_PSF
     .,ASTROMETRY_CAT,SOURCE_CAT
	integer N_EXPO,N_CHIP(NMAX_EXPO)
	common /file_stat_pass/ N_CHIP,N_EXPO
        
        character*(strl) expo,PREFIX,filename,catfile,catlist(NMAX_EXPO)
        character*(strl) dir_out
        
	double precision CRVAL(2),ra,dec
        integer num_gal(0:359,-90:90),num_expo(0:359,-90:90)
        integer num_cat(0:359,-90:90)
        integer ir,id,i,j,iexpo,ierror,ncat,icat,found

        do i=0,359
          do j=-90,90
            num_gal(i,j)=0
            num_expo(i,j)=0
            num_cat(i,j)=0
          enddo
        enddo

        ncat=0
        
        do iexpo=1,N_EXPO           

          call get_expo_name(IMAGE_FILE(iexpo,1),expo)        
  	  filename=trim(DIR_OUTPUT(iexpo))//'/result/'
     .//trim(expo)//'_all.cat'           
          open(unit=10,file=filename,status='old',iostat=ierror)
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            pause
          endif
          rewind 10
          read(10,*)
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) CRVAL(1),CRVAL(2)
            if (ierror.lt.0) cycle
            ir=int(CRVAL(1)+0.5)
            if (ir.eq.360) ir=0  
            id=int(CRVAL(2)+90.5)-90
            num_gal(ir,id)=num_gal(ir,id)+1
          enddo  
          close(10)

           
          call get_PREFIX(IMAGE_FILE(iexpo,1),PREFIX)
          filename=trim(DIR_OUTPUT(iexpo))//'/astrometry/'
     .//trim(PREFIX)//'.head'
          open(unit=10,file=filename,status='old',iostat=ierror)
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            pause
          endif
          rewind 10
          read(10,*) CRVAL(1),CRVAL(2)
          close(10)
        
          ir=int(CRVAL(1)+0.5)
          if (ir.eq.360) ir=0  
          id=int(CRVAL(2)+90.5)-90
          num_expo(ir,id)=num_expo(ir,id)+1

          
	  catfile=SOURCE_CAT(iexpo)
	  if (GALAXY_CAT_folder_provided.eq.1) then
	    call generate_gal_cat_file_name(CRVAL,catfile)
	  endif

          found=0
          do icat=1,ncat
            if (trim(catlist(icat)).eq.trim(catfile)) then
              found=1
              exit
            endif
          enddo

          if (found.eq.1) then
            cycle
          else
            ncat=ncat+1
            catlist(ncat)=catfile
          endif
             
          open(unit=10,file=catfile,status='old',iostat=ierror)
          if (ierror.ne.0) then
            write(*,*) trim(catfile),' is missing!'
	    pause 
          endif   
          rewind 10
          read(10,*)         
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) ra,dec           
 	    if (ierror.lt.0) cycle
            ir=int(ra+0.5)
            if (ir.eq.360) ir=0  
            id=int(dec+90.5)-90
            num_cat(ir,id)=num_cat(ir,id)+1
          enddo
          close(10)
          
          write(*,*) iexpo,ncat,trim(DIR_OUTPUT(iexpo))

        enddo

        call get_parent_path(DIR_OUTPUT(1),dir_out)

        filename=trim(dir_out)//'/num_cat.dat'
        open(unit=10,file=filename,status='replace')
        filename=trim(dir_out)//'/num_gal.dat'
        open(unit=20,file=filename,status='replace')
        filename=trim(dir_out)//'/num_expo.dat'
        open(unit=30,file=filename,status='replace')
        rewind 10
        rewind 20
        rewind 30
        write(10,*) '# ra   dec   num_cat  '
        write(20,*) '# ra   dec   num_gal  '
        write(30,*) '# ra   dec   num_expo  '
        do i=0,359
          do j=-90,90              
            if (num_cat(i,j).gt.0)
     .write(10,*) i,j,num_cat(i,j)
            if (num_gal(i,j).gt.0)
     .write(20,*) i,j,num_gal(i,j)
            if (num_expo(i,j).gt.0)
     .write(30,*) i,j,num_expo(i,j)
          enddo
        enddo
        close(10)
        close(20)
        close(30)

        
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine get_parent_path(dir,p_dir)
        implicit none
        
	character*(*) dir,p_dir
	integer i,p_slash,n

	p_slash=0

        n=len(trim(dir))
	do i=n,1,-1
	  if (dir(i:i).eq.'/'.and.p_slash.eq.0) then
            p_slash=i
            exit
          endif
	enddo
	
	p_dir=dir(1:p_slash-1)

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
