        subroutine combine_expo_catalog(nchip,IMAGE_FILE
     .,DIR_OUTPUT,iexpo)
        implicit none
        include 'para.inc'

c     This subroutine is for combining the shear catalog with the external source catalog.
        
        integer nchip,iexpo
	character*(strl) IMAGE_FILE(NMAX_CHIP),DIR_OUTPUT        
        character*(strl) expo,PREFIX,filename
        character*500 cat_content,cat_list1,cat_list2
        
	integer ichip,ierror,u,i,m,n
	real cat(npara),g1c,g2c

c	integer eye_check(NMAX_EXPO)
c	common /eye_pass/ eye_check

        real anomaly(NMAX_EXPO)
        common /anomaly_pass/ anomaly

        call get_expo_name(IMAGE_FILE(1),expo)        
	filename=trim(DIR_OUTPUT)//'/result/'//trim(expo)//'_all.cat'   
        
        open(unit=20,file=filename,status='replace',iostat=ierror)
        rewind 20

        do ichip=1,nchip
           
          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
	  filename
     .=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//'_orig.cat'        
          open(unit=15,file=filename,status='old',iostat=ierror)
          rewind 15
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            pause
          endif
          read(15,'(A)') cat_list2
          read(15,'(A)',iostat=ierror)  cat_content     
          if (ierror.lt.0) then
            close(15) 
            cycle
          else
            close(15)
            exit
          endif
        enddo

        
        n=0
        m=0
	do ichip=1,nchip

          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)

	  filename
     .=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//'_shear.dat'        
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            pause
          endif
          read(10,'(A)') cat_list1 

	  filename
     .=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//'_orig.cat'        
          open(unit=15,file=filename,status='old',iostat=ierror)
          rewind 15
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            pause
          endif
          read(15,*)

          if (ichip.eq.1) then
            write(20,*) trim(cat_list2),' ichip ',trim(cat_list1)
c            if (eye_check(iexpo).ne.1) then
            if (anomaly(iexpo).gt.chi2_thresh) then
              close(10)
              close(15)
              close(20)
              write(*,*) trim(DIR_OUTPUT)//'/'//trim(expo),n,m
              return
            endif
          endif
          
          do while (ierror.ge.0)

            read(10,*,iostat=ierror) (cat(u),u=1,ih2)
            if (ierror.lt.0) cycle
            read(15,'(A)') cat_content 
            
            if (cat(i_imax).ge.ns.or.cat(i_jmax).ge.ns) then
              m=m+1
              cycle
            endif                               
            n=n+1
              
            g1c=cat(igf1)+g1_c
            g2c=cat(igf2)+g2_c

            cat(ig1)=cat(ig1)-g1c*cat(ide)+g1c*cat(ih1)+g2c*cat(ih2)
            cat(ig2)=cat(ig2)-g2c*cat(ide)+g1c*cat(ih2)-g2c*cat(ih1)

            write(20,*) trim(cat_content),ichip,(cat(u),u=1,ih2)
              
          enddo
          close(10)
          close(15)
          
        enddo
        write(*,*) trim(DIR_OUTPUT),n,m


	close(20)


	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
