	real pi,arc_convert
	parameter (pi=3.1415926)
	parameter (arc_convert=pi/180.)

	integer npx,npy,strl
	parameter (npx=3000)
	parameter (npy=5000)
	parameter (strl=150)
	
	integer ASTROMETRY_trivial,PROCESS_stage
	parameter (ASTROMETRY_trivial=0)
	parameter (PROCESS_stage=11*13)

	integer GAIA_folder_provided
        parameter (GAIA_folder_provided=0)

        integer GALAXY_CAT_folder_provided
        parameter (GALAXY_CAT_folder_provided=0)

	integer ext_cat,ext_PSF,CCD_split,nct,ncx
	parameter (ext_cat=1)
	parameter (ext_PSF=0)
	parameter (CCD_split=2)
	parameter (nct=12)
	parameter (ncx=3)

	integer expo_global
	parameter (expo_global=0)

	integer nstar_min,psf_order,npo,npl,nstar_min_local
	parameter (psf_order=9)
	parameter (npo=(psf_order+2)*(psf_order+1)/2)
	parameter (nstar_min=npo*3/2)
	parameter (npl=6)
	parameter (nstar_min_local=npl*2)

	integer ns,chip_margin,ns_2,nl,nl_2,flag_thresh
	parameter (ns=48)
	parameter (chip_margin=8)
	parameter (ns_2=ns/2)
	parameter (nl_2=ns_2+chip_margin)
	parameter (nl=nl_2*2)
	parameter (flag_thresh=3)

	integer len_g,len_s,nstar_max,ngal_max,npara
	parameter (len_g=40)
	parameter (len_s=15)
	parameter (ngal_max=10000)
	parameter (nstar_max=1000)
	parameter (npara=25)

	integer npd
	parameter (npd=18)	

c 1+3+4+1+5+6+1+7+8+1+9+10+1

	integer blocksize
	parameter (blocksize=200)

	real dz_thresh
	parameter (dz_thresh=0.1)
	
	real source_thresh,core_thresh
	parameter (source_thresh=2.)	
	parameter (core_thresh=4.)


	integer area_max,area_thresh
	parameter (area_max=ns*ns)
	parameter (area_thresh=6)

	integer gal_smooth,star_smooth
	parameter (gal_smooth=2)
	parameter (star_smooth=2)
	
	real SNR_thresh_star,saturation_thresh,SNR_PSF
	parameter (SNR_thresh_star=80)
	parameter (SNR_PSF=100)
	parameter (saturation_thresh=25000)

	

	integer iclass,igf1,igf2,isnr,ira,idec,ih_flux
	integer i_imax,iexpo_num,ig1,ig2,ide,ih1,ih2,ipeak,iflag
	integer i_jmax,ih_area,iflux2,iflux_alt,istar,isig
	
	parameter (isig=4)
	parameter (istar=5)
	parameter (ipeak=5)
	parameter (i_imax=6)
	parameter (i_jmax=7)
	parameter (ih_flux=8)
	parameter (ih_area=9)
	parameter (iflag=10)
	parameter (iclass=12)
	parameter (ira=iclass+1)
	parameter (idec=iclass+2)
	parameter (igf1=iclass+3)
	parameter (igf2=iclass+4)
	parameter (ig1=iclass+5)
	parameter (ig2=iclass+6)
	parameter (ide=iclass+7)
	parameter (ih1=iclass+8)
	parameter (ih2=iclass+9)
	parameter (iexpo_num=ih2+1)
	parameter (isnr=iclass-1)
	parameter (iflux_alt=iclass)
	parameter (iflux2=isnr)

	integer NMAX_EXPO,NMAX_CHIP,num_procs_min
	parameter (NMAX_EXPO=2000)
	parameter (NMAX_CHIP=40)
	parameter (num_procs_min=40)

	integer nchip_max,ndim
	parameter (nchip_max=36)
	parameter (ndim=npd*2+nchip_max*2+2)
	
	integer ngfieldmax
	parameter (ngfieldmax=100000)

	real g1_c,g2_c
	parameter (g1_c=0.)
	parameter (g2_c=0.)

	integer np_anmly
	parameter (np_anmly=40)

	real chi2_thresh
	parameter (chi2_thresh=0.00024)