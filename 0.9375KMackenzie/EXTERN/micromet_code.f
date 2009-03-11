c micromet_code.f


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  subroutine MICROMETCODE(nx,ny,xmn,ymn,deltax,deltay,
     &  iyear_init,imonth_init,iday_init,xhour_init,dt,undef,
     &  ifill,iobsint,dn,iter,curve_len_scale,slopewt,curvewt,
     &  topo_land,curvature,terrain_slope,slope_az,topoflag,snow_d,
     &  Tair_grid,rh_grid,uwind_grid,vwind_grid,Qsi_grid,prec_grid,
     &  i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,i_prec_flag,
     &  windspd_grid,winddir_grid,windspd_flag,winddir_flag,sprec,
     &  windspd_min,Qli_grid,i_longwave_flag,vegtype,forest_LAI,
     &  iyear,imonth,iday,xhour,lapse_rate_user_flag,
     &  iprecip_lapse_rate_user_flag,use_shortwave_obs,
     &  use_longwave_obs,use_sfc_pressure_obs,sfc_pressure,
     &  calc_subcanopy_met,vegsnowdepth,gap_frac,cloud_frac_factor,
     &  barnes_lg_domain,k_stn,xlat_grid,nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,Tair_orig,rh_orig,windspd_orig,winddir_orig,prec_orig,
     &  nveg_geotop,LAIw,LAIs)

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell
      integer nstns_orig ! number of input values

      double precision xstn_orig(nstns_max)     ! input stn x coords
      double precision ystn_orig(nstns_max)     ! input stn y coords
      real Tair_orig(nstns_max)     ! input values
      real rh_orig(nstns_max)       ! input values
      real winddir_orig(nstns_max)  ! input values
      real windspd_orig(nstns_max)  ! input values
      real prec_orig(nstns_max)     ! input values
      real elev_orig(nstns_max)     ! station elevation
	  
      real dn                  ! average observation spacing
      real topo_land(nx_max,ny_max) ! grid topography
      real topo(nx_max,ny_max) ! grid topography
      real xlat_grid(nx_max,ny_max) ! lat (dec deg) of cell centers

      real Tair_grid(nx_max,ny_max)   ! output values
      real rh_grid(nx_max,ny_max)     ! output values
      real uwind_grid(nx_max,ny_max)  ! output, E-W wind component
      real vwind_grid(nx_max,ny_max)  ! output, N-S wind component
      real windspd_grid(nx_max,ny_max)
      real winddir_grid(nx_max,ny_max)
      real Qsi_grid(nx_max,ny_max)    ! output
      real Qli_grid(nx_max,ny_max)    ! output
      real prec_grid(nx_max,ny_max)   ! output
      real sprec(nx_max,ny_max)   ! output
      real sfc_pressure(nx_max,ny_max)

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour
      real dt                    ! model time step, in seconds
      integer iter               ! model iteration
      integer iyear_init     ! model start year
      integer imonth_init    ! model start month
      integer iday_init      ! model start day
      real xhour_init        ! model start hour
      integer J_day          ! model Julian day, actually day-of-year

      real undef       ! undefined value
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      real curvature(nx_max,ny_max)     ! topographic curvature
      real slope_az(nx_max,ny_max)      ! azimuth of topographic slope
      real terrain_slope(nx_max,ny_max) ! terrain slope
      real vegtype(nx_max,ny_max)
	  real snow_d(nx_max,ny_max)
      real vegsnowdepth(nvegtypes)
	  real gap_frac(nvegtypes)

      real curve_len_scale   ! length scale for curvature calculation
      real slopewt           ! wind model slope weight
      real curvewt           ! wind model curvature weight
	  real topoflag

      integer i_tair_flag,i_rh_flag,i_wind_flag,i_solar_flag,
     &  i_prec_flag,i_longwave_flag,
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag,n_stns_used

      real windspd_flag,winddir_flag,windspd_min,calc_subcanopy_met,
     &  T_lapse_rate,Td_lapse_rate,precip_lapse_rate,
     &  use_shortwave_obs,use_longwave_obs,use_sfc_pressure_obs,
     &  run_enbal,run_snowpack,cloud_frac_factor,
     &  barnes_lg_domain

      real forest_LAI(nvegtypes)
	  real LAIw(nvegtypes)
	  real LAIs(nvegtypes)

      integer k_stn(nx_max,ny_max,nstns_max)
      
      integer i,j,nveg_geotop
	  
	  n_stns_used=nstns_orig
	  	    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Check max_values
	  if(nveg_geotop.ne.nvegtypes) then
	    write(6,*) 'nveg_geotop:',nveg_geotop,' nveg_liston:',
     &   nvegtypes
        write(6,*) 'INCONSISTENCY ERROR, set them equal'
		read(*,*)
      endif
	  
	  if(nx.ne.nx_max) then
	    write(6,*) 'nx:',nx,' nx_max:',nx_max
        write(6,*) 'INCONSISTENCY ERROR, set them equal'
		read(*,*)
      endif	  

	  if(ny.ne.ny_max) then
	    write(6,*) 'ny:',ny,' ny_max:',ny_max
        write(6,*) 'INCONSISTENCY ERROR, set them equal'
		read(*,*)
      endif	  

	  if(nstns_orig.ne.nstns_max) then
	    write(6,*) 'nstns_orig:',nstns_orig,' nstns_max:',nstns_max
        write(6,*) 'INCONSISTENCY ERROR, set them equal'
		read(*,*)
      endif	  	 

c Topography including or not snow depth
      if (topoflag.eq.1.0) then
        do i=1,nx
          do j=1,ny
			topo(i,j) = topo_land(i,j) + snow_d(i,j)
          enddo
        enddo
      elseif (topoflag.eq.0.0) then
        do i=1,nx
          do j=1,ny
            topo(i,j) = topo_land(i,j)
          enddo
        enddo
      endif

	   	  
c Calculate what the current simulation date should be.
      call get_model_time(iyear_init,imonth_init,iday_init,
     &  xhour_init,iter,dt,iyear,imonth,iday,xhour,J_day)


c Read in the observations for this time step, and build an array of
c   valid observations to be interpolated.
      call get_obs_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,
     &  windspd_orig,winddir_orig,prec_orig,iter)

c      write(6,*) nx,' ',ny,' ',nx_max,' ',ny_max
c	  read(*,*)
c	  do i=1,nx
c	    do j=1,ny
c		  write(6,*) 'i:',i,' j:',j,' ',vegtype(i,j)
c		enddo
c	  enddo
c	  read(*,*)
	  
	 	 
c Make the topographic calculations required by the wind and solar
c   radiation models.  These calculations are not fixed in time
c   because as the snow depth evolves it modifies the "topography".
      call topo_data(nx,ny,deltax,deltay,topo,
     &  curvature,terrain_slope,slope_az,curve_len_scale,undef)

c Calculate the temperature and dew-point lapse rates to be used in
c   the interpoations.
      call get_lapse_rates(imonth,iday,T_lapse_rate,
     &  Td_lapse_rate,xlat_grid(1,1),lapse_rate_user_flag,
     &  precip_lapse_rate,iprecip_lapse_rate_user_flag)

c Calculate the forest lai for each of the five forest types, and
c   for this day of the simulation (the lai varies seasonally for
c   the case of deciduous trees).
      call get_lai(J_day,forest_LAI,LAIw,LAIs)

c TEMPERATURE.
      if (i_tair_flag.eq.1) then
c       print *,'   solving for temperature'
       call temperature(nx,ny,deltax,deltay,xmn,ymn,
     &    nstns_orig,xstn_orig,ystn_orig,Tair_orig,dn,Tair_grid,
     &    undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,
     &    topo,T_lapse_rate,barnes_lg_domain,n_stns_used,k_stn)
c       do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,Tair_grid(i,j),topo(i,j)
c           enddo
c        enddo
c		write(6,*) 'TEMPERATURE'
c        read(*,*)   
      endif

c RELATIVE HUMIDITY.
      if (i_rh_flag.eq.1) then
c       print *,'   solving for relative humidity'
        call relative_humidity(nx,ny,deltax,deltay,xmn,ymn,
     &    nstns_orig,xstn_orig,ystn_orig,rh_orig,dn,rh_grid,undef,
     &    ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,topo,
     &    Tair_orig,Tair_grid,Td_lapse_rate,barnes_lg_domain,
     &    n_stns_used,k_stn)
c        do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,rh_grid(i,j)
c           enddo
c        enddo
c        write(6,*) 'RH'
c        read(*,*) 
      endif

c WIND SPEED AND DIRECTION.
      if (i_wind_flag.eq.1) then
c       print *,'   solving for wind speed and direction'
        call wind(nx,ny,deltax,deltay,xmn,ymn,windspd_orig,
     &    nstns_orig,xstn_orig,ystn_orig,dn,undef,ifill,
     &    iobsint,iyear,imonth,iday,xhour,elev_orig,
     &    winddir_orig,uwind_grid,vwind_grid,slopewt,curvewt,
     &    curvature,slope_az,terrain_slope,windspd_grid,
     &    winddir_grid,windspd_flag,winddir_flag,windspd_min,
     &    vegtype,forest_LAI,calc_subcanopy_met,vegsnowdepth,
     &    barnes_lg_domain,n_stns_used,k_stn,topo)
c        do i=1,nx
c          do j=1,ny
c                write(6,*) i,j,winddir_grid(i,j),uwind_grid(i,j),
c     &                     vwind_grid(i,j),topo(i,j)
c           enddo
c        enddo
c        write(6,*) 'WIND'   
c		read(*,*) 
      endif

c SOLAR RADIATION.
      
      
      if (i_solar_flag.eq.1) then
c       print *,'   solving for solar radiation'
        call solar(nx,ny,xhour,J_day,topo,rh_grid,Tair_grid,
     &    xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,
     &    forest_LAI,T_lapse_rate,Td_lapse_rate,
     &    calc_subcanopy_met,gap_frac,cloud_frac_factor,undef)
c        do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,Qsi_grid(i,j),topo(i,j)
c           enddo
c        enddo
c        write(6,*) 'SHORTWAVE'   
c        read(*,*)     

c If requested, modify the model output to account for shortwave
c   radiation observations.
        if (use_shortwave_obs.eq.1.0) then
          if (barnes_lg_domain.eq.1.0) then
            print *,'The model is not configured to assimilate'
            print *,'  solar data with barnes_lg_domain = 1.0.'
            stop
          endif
          call shortwave_data(nx,ny,deltax,deltay,xmn,ymn,
     &      iyear,imonth,iday,xhour,undef,Qsi_grid,iter)
        endif
      endif
 
c INCOMING LONGWAVE RADIATION.
      if (i_longwave_flag.eq.1) then
c        print *,'   solving for incoming longwave radiation'
        call longwave(nx,ny,rh_grid,Tair_grid,Qli_grid,topo,
     &    vegtype,forest_LAI,T_lapse_rate,Td_lapse_rate,
     &    calc_subcanopy_met,cloud_frac_factor,undef)
c        do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,Qli_grid(i,j)
c           enddo
c        enddo
c        write(6,*) 'LONGWAVE'   
c        read(*,*)       

c If requested, modify the model output to account for longwave
c   radiation observations.
        if (use_longwave_obs.eq.1.0) then
          if (barnes_lg_domain.eq.1.0) then
            print *,'The model is not configured to assimilate'
            print *,'  longwave data with barnes_lg_domain = 1.0.'
            stop
          endif
          call longwave_data(nx,ny,deltax,deltay,xmn,ymn,
     &      iyear,imonth,iday,xhour,undef,Qli_grid,iter)
        endif
      endif

c PRECIPITATION.
      if (i_prec_flag.eq.1) then
c        print *,'   solving for precipitation'
        call precipitation(nx,ny,deltax,deltay,xmn,ymn,
     &    nstns_orig,xstn_orig,ystn_orig,prec_orig,dn,prec_grid,
     &    undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,
     &    topo,Tair_grid,sprec,iter,
     &    precip_lapse_rate,barnes_lg_domain,n_stns_used,k_stn)
c        do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,prec_grid(i,j)
c           enddo
c        enddo
c        write(6,*) 'PRECIPITATION'   
c        read(*,*)  
      endif

c SURFACE PRESSURE.
c Surface pressure is used in EnBal and SnowMass.  If needed for
c   this SnowModel simulation, calculate the distribution here.
	  call pressure(nx,ny,topo,sfc_pressure,undef)

c If requested, modify the model output to account for surface
c   pressure observations.
	  if (use_sfc_pressure_obs.eq.1.0) then
		if (barnes_lg_domain.eq.1.0) then
		  print *,'The model is not configured to assimilate'
		  print *,'  pressure data with barnes_lg_domain = 1.0.'
		  stop
		endif
		call sfc_pressure_data(nx,ny,deltax,deltay,xmn,ymn,
     &      iyear,imonth,iday,xhour,undef,sfc_pressure,iter)
	  endif
c        do i=1,nx
c           do j=1,ny
c                write(6,*) i,j,sfc_pressure(i,j)
c           enddo
c        enddo
c        write(6,*) 'PRESSURE'   
c        read(*,*)            
	  
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine precipitation(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns_orig,xstn_orig,ystn_orig,prec_orig,dn,prec_grid,
     &  undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,
     &  topo,Tair_grid,sprec,iter,
     &  precip_lapse_rate,barnes_lg_domain,n_stns_used,k_stn)

c Interpolate the observed precipitation values to the grid.  Also
c   interpolate the station elevations to a reference surface.  Use
c   a precipitation "lapse rate", or adjustment factor to define
c   the precipitation on the actual elevation grid.  The reason the
c   interpolated station elevations are used as the topographic
c   reference surface (instead of something like sea level), is
c   because the precipitation adjustment factor is a non-linear 
c   function of elevation difference.
 
c The adjustment factor that is used comes from: Thornton, P. E.,
c   S. W. Running, and M. A. White, 1997: Generating surfaces of
c   daily meteorological variables over large regions of complex
c   terrain.  J. Hydrology, 190, 214-251.

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real prec(nstns_max) ! input values
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation
      real prec_orig(nstns_max) ! input values

      real dn                           ! average observation spacing
      real topo(nx_max,ny_max)          ! grid topography
      real prec_grid(nx_max,ny_max)     ! output values
      real Tair_grid(nx_max,ny_max)     ! input values
      real sprec(nx_max,ny_max)         ! output values
      real topo_ref_grid(nx_max,ny_max) ! reference surface

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real delta_topo,alfa,Tf,precip_lapse_rate_m,precip_lapse_rate,
     &  barnes_lg_domain
      integer i,j,iter,n_stns_used
      integer k_stn(nx_max,ny_max,nstns_max)

c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,prec_orig,prec)

c Use the barnes oi scheme to interpolate the station elevation data
c   to the grid, so that it can be used as a topographic reference
c   surface.
	  
      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,elev,dn,topo_ref_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c Use the barnes oi scheme to interpolate the station data to
c   the grid.

      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,prec,dn,prec_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c Convert the precipitation "lapse rate" (km-1) to m-1.
      precip_lapse_rate_m = precip_lapse_rate / 1000.0

c Convert the gridded station data to the actual gridded elevations.
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref_grid(i,j)
c Don't let the elevation difference be greater than some number
c  (like 1800 meters gives a factor of 4.4).  If it is too large
c   you get huge precipitation adjustment, a divide by zero, or
c   even negative adjustments for high elevations).
		  delta_topo = min(delta_topo,1800.0)
		  alfa = precip_lapse_rate_m * delta_topo
		  prec_grid(i,j) = prec_grid(i,j) * (1.0 + alfa)/(1.0 - alfa)
        enddo
      enddo

c Convert the precipitation values from mm to m swe.  Also, make
c   sure the interpolation has not created any negetive
c   precipitation values.
      do j=1,ny
        do i=1,nx
		  prec_grid(i,j) = prec_grid(i,j) / 1000.0
		  prec_grid(i,j) = max(0.0,prec_grid(i,j))
        enddo
      enddo

c Use the temperature distribution to define whether this
c   precipitation is falling as rain or snow (generate a
c   snow-precipitation array following the air temperature
c   threshold defined by Auer (1974) = 2.0 C).  Note here that,
c   if you ever want it, rain_prec = prec_grid - sprec.  This
c   snow precipitation (sprec) is assumed to be in meters
c   snow-water-equivalent per time step.
      Tf = 273.16
      do j=1,ny
        do i=1,nx
		  if (Tair_grid(i,j).lt.2.0+Tf) then
			sprec(i,j) = prec_grid(i,j)
		  else
		    sprec(i,j) = 0.0
		  endif
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine longwave(nx,ny,rh_grid,Tair_grid,Qli_grid,topo,
     &  vegtype,forest_LAI,T_lapse_rate,Td_lapse_rate,
     &  calc_subcanopy_met,cloud_frac_factor,undef)

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real Tair_grid(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)
      real Qli_grid(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real vegtype(nx_max,ny_max)

      real T_lapse_rate,Td_lapse_rate,es,e,emiss_cloud,
     &  A,B,C,Tf,Stef_Boltz,cloud_frac,E1,X1,Y1,Z1,E2,X2,Y2,Z2,
     &  Xs,Ys,Zs,forest_frac,E3,X3,Y3,Z3,alfa,calc_subcanopy_met,
     &  cloud_frac_factor,undef

      real forest_LAI(nvegtypes)

      integer i,j,nveg

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

c Over ice.
c     A = 6.1115 * 100.0
c     B = 22.452
c     C = 272.55

c Define the freezing temperature to be used to convert from C to K.
      Tf = 273.16

c Define the Stefan Boltzmann constant.
      Stef_Boltz = 5.6696e-8

c Constants required for Iziomon et al. (2003).
      E1 = 200.0
      X1 = 0.35
      Y1 = 0.100
      Z1 = 0.224

      E2 = 1500.0
      X2 = 0.43
      Y2 = 0.115
      Z2 = 0.320

c Assume the X and Y coefficients increase linearly to 3000 m,
c   and adjust Z to create a best fit to the CLPX data.
      E3 = 3000.0
      X3 = 0.51
      Y3 = 0.130
      Z3 = 1.100

      do j=1,ny
        do i=1,nx
          
c Compute the cloud fraction.
           call get_cloudfrac(Tair_grid(i,j),rh_grid(i,j),topo(i,j),
     &      cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

c Calculate the vapor pressure.
           es = A * exp((B * (Tair_grid(i,j) - Tf))/
     &      (C + (Tair_grid(i,j) - Tf)))
           e = es * rh_grid(i,j) / 100.0

c Compute Qli following Iziomon et al. (2003).
           if (topo(i,j).lt.E1) then
            Xs = X1
            Ys = Y1
            Zs = Z1
           elseif (topo(i,j).gt.E2) then
            Xs = X3
            Ys = Y3
            Zs = Z3
           else
            Xs = X1 + (topo(i,j) - E1) * (X3 - X1)/(E3 - E1)
            Ys = Y1 + (topo(i,j) - E1) * (Y3 - Y1)/(E3 - E1)
            Zs = Z1 + (topo(i,j) - E1) * (Z3 - Z1)/(E3 - E1)
           endif

           alfa = 1.083
           emiss_cloud = alfa *
     &      (1.0 - Xs * exp((- Ys) * e/Tair_grid(i,j))) *
     &      (1.0 + Zs * cloud_frac**2)

           Qli_grid(i,j) = emiss_cloud * Stef_Boltz * Tair_grid(i,j)**4

c Modify the incoming longwave radiation for the forest canopy.
           if (vegtype(i,j).le.(real(nvegtypes))) then
            if (calc_subcanopy_met.eq.1.0) then

c Define the forest-canopy par.
              nveg = nint(vegtype(i,j))
              if (forest_LAI(nveg).lt.0.2) then
                forest_frac = 0.5 * forest_LAI(nveg)
              else
                forest_frac =
     &            min(1.0,0.55 + 0.29 * log(forest_LAI(nveg)))
              endif
              Qli_grid(i,j) = Qli_grid(i,j) * (1.0 - forest_frac) +
     &          (Stef_Boltz * Tair_grid(i,j)**4) * forest_frac
            endif
		   endif
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solar(nx,ny,xhour,J_day,topo,rh_grid,Tair_grid,
     &  xlat_grid,Qsi_grid,slope_az,terrain_slope,dt,vegtype,
     &  forest_LAI,T_lapse_rate,Td_lapse_rate,
     &  calc_subcanopy_met,gap_frac,cloud_frac_factor,undef)

c First take the surface gridded fields of Tair and RH, and
c   calculate Td for the topographic surface.  Then use Tair and
c   Td, and the associated lapse rates, and calculate Tair and Td
c   for the 700 mb level.  Use these surfaces to calculate RH at
c   700 mb, and convert these values to cloud fraction following
c   Walcek, C. J., 1994: Cloud cover and its relationship to
c   relative humidity during a spring midlatitude cyclone.  Mon.
c   Wea. Rev., 122, 1021-1035.

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real topo(nx_max,ny_max)
      real Tair_grid(nx_max,ny_max)
      real rh_grid(nx_max,ny_max)
      real Qsi_grid(nx_max,ny_max)
      real slope_az(nx_max,ny_max)
      real terrain_slope(nx_max,ny_max)
      real vegtype(nx_max,ny_max)
      real xlat_grid(nx_max,ny_max)

      real xhour                 ! model decimal hour
      real xxhour
      integer J_day              ! model day of year
      integer i,j,ihrs_day,ihour,nveg
      real dt,cloud_frac,Qsi_tmp,Qsi_sum,trans_veg,
     &  T_lapse_rate,Td_lapse_rate,calc_subcanopy_met,
     &  cloud_frac_factor,undef

      real forest_LAI(nvegtypes)
	  real gap_frac(nvegtypes)

      ihrs_day = 24

      do j=1,ny
        do i=1,nx
            
c Compute the cloud fraction.
          call get_cloudfrac(Tair_grid(i,j),rh_grid(i,j),topo(i,j),
     &      cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

c Compute the incoming solar radiation.  The solar_rad subroutine
c   calculates the instantaneous incoming solar radiation, so
c   if the time step is very long, account for this by calculating
c   the incoming solar radiation every 3 hours and then taking the
c   average.
          if (dt.le.10800.0) then
            call solar_rad(Qsi_grid(i,j),J_day,xlat_grid(i,j),
     &        cloud_frac,xhour,slope_az(i,j),terrain_slope(i,j))
          
c            if(i.le.2) then
c                if(j.le.2) then
c                    write(*,*) i,j,Qsi_grid(i,j)
c                    read(*,*)
c                endif
c            endif
                
          elseif (dt.eq.86400.0) then
            Qsi_sum = 0.0
            do ihour=3,ihrs_day,3
              xxhour = real(ihour)
              call solar_rad(Qsi_tmp,J_day,xlat_grid(i,j),
     &          cloud_frac,xxhour,slope_az(i,j),terrain_slope(i,j))
                Qsi_sum = Qsi_sum + Qsi_tmp
            enddo
            Qsi_grid(i,j) = Qsi_sum / (real(ihrs_day)/3.0)
c            write(6,*) i,j,Qsi_grid(i,j)
          else
            print *,'The model may not do what you want with this dt'
            stop
          endif

c Modify the incoming solar radiation for the forest canopy.
          if (vegtype(i,j).le.(real(nvegtypes))) then
            if (calc_subcanopy_met.eq.1.0) then

c Define the forest-canopy transmissivity.  0.71 provided a
c   best-fit to the observations, when averaged over the two years
c   of hourly data.
              nveg = nint(vegtype(i,j))
              trans_veg = exp((- 0.71) * forest_LAI(nveg))

c Account for any gaps in the forest canopy that will allow
c   direct incoming solar radiation to reach the snow surface.
              trans_veg = gap_frac(nveg) * (1.0 - trans_veg) + trans_veg

              Qsi_grid(i,j) = trans_veg * Qsi_grid(i,j)
            endif
          endif

        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_cloudfrac(Tair_grid,rh_grid,topo,
     &  cloud_frac,T_lapse_rate,Td_lapse_rate,cloud_frac_factor)

      implicit none

      real Td_lapse_rate,topo_ref,delta_topo,A,B,C,e,es,dx,
     &  T_lapse_rate,press_ratio,f_max,one_minus_RHe,f_1,
     &  Td_700,Tair_700,rh_700,cloud_frac,Tair_grid,rh_grid,
     &  topo,Td_grid,Tf,cloud_frac_factor

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

c Over ice.
c     A = 6.1115 * 100.0
c     B = 22.452
c     C = 272.55

c Define the freezing temperature to be used to convert from C to K.
      Tf = 273.16

c Assume that 700 mb is equivalent to 3000 m in a standard
c   atmosphere.
      topo_ref = 3000.0

c Define the ratio of 700 mb level pressure to the surface pressure
c   (~1000 mb).
      press_ratio = 0.7

c Assume dx = 80.0 km, for Walcek (1994).
      dx = 80.0

c Walcek coefficients.
      f_max = 78.0 + 80.0/15.5
      one_minus_RHe = 0.196 + (0.76-80.0/2834.0) * (1.0 - press_ratio)
      f_1 = f_max * (press_ratio - 0.1) / 0.6 / 100.0

c Convert the gridded topo-surface RH to Td.
      es = A * exp((B * (Tair_grid - Tf))/(C + (Tair_grid - Tf)))
      e = es * max(10.0,rh_grid) / 100.0
      Td_grid = C * log(e/A) / (B - log(e/A)) + Tf

c Convert the topo-surface temperature values to 700 mb values.
      delta_topo = topo - topo_ref
      Td_700 = Td_grid + Td_lapse_rate * delta_topo
      Tair_700 = Tair_grid + T_lapse_rate * delta_topo

c Convert each Td to a gridded relative humidity (0-1).
      e = A * exp((B * (Td_700 - Tf))/(C + (Td_700 - Tf)))
      es = A * exp((B * (Tair_700 - Tf))/(C + (Tair_700 - Tf)))
      rh_700 = e/es
      rh_700 = min(1.0,rh_700)
      rh_700 = max(0.0,rh_700)

c Use this RH at 700 mb to define the cloud fraction (0-1).
      cloud_frac = f_1 * exp((rh_700 - 1.0)/one_minus_RHe)

c If the user wants to, reduce the calculate cloud fraction by the
c   cloud_frac_factor.
      cloud_frac = cloud_frac_factor * cloud_frac

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solar_rad(Qsi,J_day,xlat,
     &  cloud_frac,xxhour,slope_az,terrain_slope)

      implicit none

      integer J_day

      real solar_const,days_yr,Trop_Can,solstice,pi,deg2rad,
     &  Qsi_direct,Qsi_diffuse,cos_i,cos_Z,Qsi,xlat,sin_z,xxhour,
     &  cloud_frac,slope_az,terrain_slope,sol_dec,hr_angl,
     &  trans_direct,trans_diffuse,Qsi_trans_dir,Qsi_trans_dif,
     &  sun_azimuth,slope_az_S0

c Required constants.
      solar_const = 1370.
      days_yr = 365.25
      Trop_Can = 0.41
      solstice = 173.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

c COMPUTE THE BASIC SOLAR RADIATION PAR.

c Compute the solar declination angle (radians).
      sol_dec = Trop_Can *
     &  cos(2.*pi * (real(J_day) - solstice)/days_yr)
      
c Compute the sun's hour angle (radians).
      hr_angl = (xxhour * 15.0 - 180.0) * deg2rad

c Compute cos_Z.  Note that the sin of the solar elevation angle,
c   sin_alfa, is equal to the cosine of the solar zenith angle,
c   cos_Z.
      cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + 
     &  cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl)
      cos_Z = max(0.0,cos_Z)

c Account for clouds, water vapor, pollution, etc.
      trans_direct = (0.6 + 0.2 * cos_Z) * (1.0 - cloud_frac)
      trans_diffuse = (0.3 + 0.1 * cos_Z) * cloud_frac

c Compute the solar radiation transmitted through the atmosphere.
      Qsi_trans_dir = solar_const * trans_direct
      Qsi_trans_dif = solar_const * trans_diffuse

c COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.

c The sine of the solar zenith angle.
      sin_Z = sqrt(1.0 - cos_Z*cos_Z)

c Azimuth of the sun, with south having zero azimuth.
      sun_azimuth = 
     &  asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))

c Make the corrections so that the angles below the local horizon
c   are still measured from the normal to the slope.
      if (hr_angl.lt.0.0) then
        if (hr_angl.lt.sun_azimuth) sun_azimuth = - pi - sun_azimuth
      elseif (hr_angl.gt.0.0) then
        if (hr_angl.gt.sun_azimuth) sun_azimuth = pi - sun_azimuth
      endif

c Build, from the variable with north having zero azimuth, a 
c   slope_azimuth value with south having zero azimuth.
      if (slope_az.ge.180.0) then
        slope_az_S0 = slope_az - 180.0
      else
        slope_az_S0 = slope_az + 180.0
      endif

c Compute the angle between the normal to the slope and the angle
c   at which the direct solar radiation impinges on the sloping
c   terrain (radians).
      cos_i = cos(terrain_slope * deg2rad) * cos_Z + 
     &  sin(terrain_slope * deg2rad) * sin_Z * 
     &  cos(sun_azimuth - slope_az_S0 * deg2rad)

c Adjust the topographic correction due to local slope so that
c   the correction is zero if the sun is below the local horizon 
c   (i.e., the slope is in the shade) or if the sun is below the
c   global horizon.
      if (cos_i.lt.0.0) cos_i = 0.0
      if (cos_Z.le.0.0) cos_i = 0.0

c Adjust the solar radiation for slope, etc.
      Qsi_direct = cos_i * Qsi_trans_dir
      Qsi_diffuse = cos_Z * Qsi_trans_dif

c Combine the direct and diffuse solar components.
      Qsi = Qsi_direct + Qsi_diffuse

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wind(nx,ny,deltax,deltay,xmn,ymn,windspd_orig,
     &  nstns_orig,xstn_orig,ystn_orig,dn,undef,ifill,
     &  iobsint,iyear,imonth,iday,xhour,elev_orig,
     &  winddir_orig,uwind_grid,vwind_grid,slopewt,curvewt,
     &  curvature,slope_az,terrain_slope,windspd_grid,
     &  winddir_grid,windspd_flag,winddir_flag,windspd_min,
     &  vegtype,forest_LAI,calc_subcanopy_met,vegsnowdepth,
     &  barnes_lg_domain,n_stns_used,k_stn,topo)

c This program takes the station wind speed and direction, converts
c   them to u and v components, interpolates u and v to a grid,
c   converts the gridded values to speed and direction, and then
c   runs a simple wind model that adjusts those speeds and
c   directions according to topographic slope and curvature
c   relationships.  The resulting speeds and directions are
c   converted to u and v components and passed back to the main
c   program to be written to a file.  (All of the conversion
c   between u-v and speed-dir is done because of the problems
c   with interpolating over the 360/0 direction line.)

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation

      real windspd_orig(nstns_max) ! input values
      real winddir_orig(nstns_max) ! input values

      real speed(nstns_max)  ! input values
      real dir(nstns_max)    ! input values
      real u(nstns_max)      ! u component of wind
      real v(nstns_max)      ! v component of wind

      real dn                  ! average observation spacing
      real uwind_grid(nx_max,ny_max)  ! output values
      real vwind_grid(nx_max,ny_max)  ! output values
      real u_grid(nx_max,ny_max)  ! temporary u wind component
      real v_grid(nx_max,ny_max)  ! temporary v wind component
      real winddir_grid(nx_max,ny_max)  ! temporary wind direction
      real windspd_grid(nx_max,ny_max)  ! temporary wind speed

      real curvature(nx_max,ny_max)     ! topographic curvature
      real slope_az(nx_max,ny_max)      ! azimuth of topographic slope
      real terrain_slope(nx_max,ny_max) ! terrain slope
      real vegtype(nx_max,ny_max)
      real vegsnowdepth(nvegtypes)
      real topo(nx_max,ny_max)

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real pi,deg2rad,rad2deg,slopewt,curvewt
      integer i,j,k,n_stns_used
      integer k_stn(nx_max,ny_max,nstns_max)      
      real windspd_flag,winddir_flag,u_sum,v_sum,windspd_min,
     &  calc_subcanopy_met,barnes_lg_domain

      real forest_LAI(nvegtypes)

	 
c Define the required constants.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

c Filter through the original input data, and eliminate any
c   missing values (i.e., make sure each wind direction is paired
c   up with a wind speed.
      call get_good_values2(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,windspd_orig,winddir_orig,
     &  speed,dir)

c Convert these station data to u and v wind components.
      do k=1,nstns
        speed(k) = max(windspd_min,speed(k))
        u(k) = (- speed(k)) * sin(deg2rad * dir(k))
        v(k) = (- speed(k)) * cos(deg2rad * dir(k))
      enddo

c Use the barnes oi scheme to interpolate the station data to
c   the grid.
c U component.

      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,u,dn,u_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c V component.

      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,v,dn,v_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c Convert these u and v components to speed and directions.
      do j=1,ny
        do i=1,nx
          winddir_grid(i,j) = 270.0 -
     &      rad2deg*atan2(v_grid(i,j),u_grid(i,j))
          if (winddir_grid(i,j).ge.360.0)
     &      winddir_grid(i,j) = winddir_grid(i,j)-360.0
          windspd_grid(i,j) = sqrt(u_grid(i,j)**2 + v_grid(i,j)**2)
        enddo
      enddo

c Modify the wind speed and direction according to simple
c   wind-topography relationships.
      call topo_mod_winds(nx,ny,winddir_grid,slopewt,curvewt,
     &  windspd_grid,uwind_grid,vwind_grid,curvature,slope_az,
     &  terrain_slope,vegtype,forest_LAI,calc_subcanopy_met,
     &  vegsnowdepth,topo,undef)

c Avoid problems of zero (low) winds (for example, turbulence
c   theory, log wind profile, etc., says that we must have some
c   wind.  Thus, some equations blow up when the wind speed gets
c   very small).
      do j=1,ny
        do i=1,nx
          if (windspd_grid(i,j).lt.windspd_min) then
            windspd_grid(i,j) = windspd_min
            uwind_grid(i,j) = (- windspd_grid(i,j)) *
     &        sin(deg2rad*winddir_grid(i,j))
            vwind_grid(i,j) = (- windspd_grid(i,j)) *
     &        cos(deg2rad*winddir_grid(i,j))
          endif
        enddo
      enddo

c Find the maximum wind speed in the domain, and the
c   domain-averaged wind direction.
      windspd_flag = 0.0
      u_sum = 0.0
      v_sum = 0.0
      do j=1,ny
        do i=1,nx
          windspd_flag = max(windspd_flag,windspd_grid(i,j))
          u_sum = u_sum + uwind_grid(i,j)
          v_sum = v_sum + vwind_grid(i,j)
        enddo
      enddo
      u_sum = u_sum / real(nx*ny)
      v_sum = v_sum / real(nx*ny)
      winddir_flag = 270.0 - rad2deg*atan2(v_sum,u_sum)
      if (winddir_flag.ge.360.0) winddir_flag = winddir_flag-360.0

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine relative_humidity(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns_orig,xstn_orig,ystn_orig,rh_orig,dn,rh_grid,undef,
     &  ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,topo,
     &  Tair_orig,Tair_grid,Td_lapse_rate,barnes_lg_domain,
     &  n_stns_used,k_stn)

c This procedure follows: Kunkel, K. E., 1989: Simple procedures for
c   extrapolation of humidity variables in the mountainous Western
c   United States. J. Climate, 2, 656-669.

c First convert stn relative humidity to dew-point temperature.  Use
c   the Td lapse rate to take the stn Td to sea level.  Interpolate
c   the stn Td to the grid.  Use the Td lapse rate to take the sea
c   level grid to the actual elevations.  Convert each Td to
c   relative humidity.

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation

      real Tair_orig(nstns_max)  ! input values
      real rh_orig(nstns_max)    ! input values

      real Tair(nstns_max)  ! input values
      real Td(nstns_max)    ! input values
      real rh(nstns_max)    ! input values

      real dn                  ! average observation spacing
      real topo(nx_max,ny_max) ! grid topography
      real Tair_grid(nx_max,ny_max)  ! output values
      real Td_grid(nx_max,ny_max)    ! output values
      real rh_grid(nx_max,ny_max)    ! output values

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real Td_lapse_rate,topo_ref,delta_topo,A,B,C,e,es,Tf,
     &  barnes_lg_domain
      integer i,j,k,n_stns_used
      integer k_stn(nx_max,ny_max,nstns_max)
	  
c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

c Over ice.
c     A = 6.1115 * 100.0
c     B = 22.452
c     C = 272.55

c Define the freezing temperature to be used to convert from C to K.
      Tf = 273.16

c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values2(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,Tair_orig,rh_orig,
     &  Tair,rh)

c Convert the stn relative humidity to Td.
      do k=1,nstns

c Saturation vapor pressure at temperature, T.
        es = A * exp((B * (Tair(k) - Tf))/(C + (Tair(k) - Tf)))

c Dew point temperature for a given temperature and relative humidity.
        e = es * max(10.0,rh(k)) / 100.0
        Td(k) = C * log(e/A) / (B - log(e/A)) + Tf

      enddo

c Define the topographic reference surface.
      topo_ref = 0.0

c Convert the station data to sea level values.
      do k=1,nstns
        delta_topo = topo_ref - elev(k)
        Td(k) = Td(k) + Td_lapse_rate * delta_topo
      enddo

c Use the barnes oi scheme to interpolate the station data to
c   the grid.

      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,Td,dn,Td_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c Convert these grid values back to the actual gridded elevations.
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref
          Td_grid(i,j) = Td_grid(i,j) + Td_lapse_rate * delta_topo
        enddo
      enddo

c Convert each Td to a gridded relative humidity.
      do j=1,ny
        do i=1,nx
          e = A * exp((B * (Td_grid(i,j) - Tf)) /
     &      (C + (Td_grid(i,j) - Tf)))
          es = A * exp((B * (Tair_grid(i,j) - Tf)) /
     &     (C + (Tair_grid(i,j) - Tf)))
          rh_grid(i,j) = 100.0 * e/es

c Make sure the interpolation processes has not created any values
c   above 100 and below 0.
          rh_grid(i,j) = min(100.0,rh_grid(i,j))
          rh_grid(i,j) = max(0.0,rh_grid(i,j))
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine temperature(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns_orig,xstn_orig,ystn_orig,Tair_orig,dn,Tair_grid,
     &  undef,ifill,iobsint,iyear,imonth,iday,xhour,elev_orig,
     &  topo,T_lapse_rate,barnes_lg_domain,n_stns_used,k_stn)

c The lapse rate used depends on the month of the year, and is
c   defined by: Kunkel, K. E., 1989: Simple procedures for
c   extrapolation of humidity variables in the mountainous Western
c   United States. J. Climate, 2, 656-669.

c First adjust the stn temperatures to a common level (sea level),
c   assuming this lapse rate.  Then interpolate the temperatures
c   to the model grid.  Then use the topography data and lapse
c   rate to adjust the gridded temperatures values back to the
c   actual elevation.

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real Tair(nstns_max) ! input values
      real elev(nstns_max) ! station elevation
      real undef           ! undefined value

      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! station elevation
      real Tair_orig(nstns_max) ! input values

      real dn                  ! average observation spacing
      real topo(nx_max,ny_max) ! grid topography
      real Tair_grid(nx_max,ny_max) ! output values

      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      real T_lapse_rate,topo_ref,delta_topo,barnes_lg_domain
      integer i,j,k,n_stns_used
      integer k_stn(nx_max,ny_max,nstns_max)
	  
c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,Tair_orig,Tair)

c Define the topographic reference surface.
      topo_ref = 0.0

c Convert the station data to sea level values.
      do k=1,nstns
        delta_topo = topo_ref - elev(k)
        Tair(k) = Tair(k) + T_lapse_rate * delta_topo
      enddo

c Use the barnes oi scheme to interpolate the station data to
c   the grid.

      call interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,Tair,dn,Tair_grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

c Convert these grid values back to the actual gridded elevations.
      do j=1,ny
        do i=1,nx
          delta_topo = topo(i,j) - topo_ref
          Tair_grid(i,j) = Tair_grid(i,j) + T_lapse_rate * delta_topo
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
      subroutine topo_mod_winds(nx,ny,winddir_grid,slopewt,curvewt,
     &  windspd_grid,uwind_grid,vwind_grid,curvature,slope_az,
     &  terrain_slope,vegtype,forest_LAI,calc_subcanopy_met,
     &  vegsnowdepth,topo,undef)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,nveg

      real pi,deg2rad,rad2deg,slopewt,curvewt,dirdiff,windwt,
     &  wslope_max,beta,veg_ht,a,canopy_windwt,calc_subcanopy_met,undef
      real forest_LAI(nvegtypes)

      real curvature(nx_max,ny_max)
      real windspd_grid(nx_max,ny_max)
      real winddir_grid(nx_max,ny_max)
      real uwind_grid(nx_max,ny_max)
      real vwind_grid(nx_max,ny_max)
      real wind_slope(nx_max,ny_max)
      real slope_az(nx_max,ny_max)
      real terrain_slope(nx_max,ny_max)
      real vegtype(nx_max,ny_max)
      real topo(nx_max,ny_max)
      real vegsnowdepth(nvegtypes)

c Compute the wind modification factor which is a function of
c   topography and wind direction following Liston and Sturm (1998).

c Define the required constants.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0
      rad2deg = 180.0 / pi

c Compute the slope in the direction of the wind.
      do i=1,nx
        do j=1,ny
          wind_slope(i,j) = deg2rad * terrain_slope(i,j) *
     &      cos(deg2rad * (winddir_grid(i,j) - slope_az(i,j)))
        enddo
      enddo

c Scale the wind slope such that the max abs(wind slope) has a value
c   of abs(0.5).  Include a 1 mm slope in slope_max to prevent
c   divisions by zero in flat terrain where the slope is zero.
      wslope_max = 0.0 + 0.001
      do j=1,ny
        do i=1,nx
          wslope_max = max(wslope_max,abs(wind_slope(i,j)))
        enddo
      enddo
      do j=1,ny
        do i=1,nx
          wind_slope(i,j) = wind_slope(i,j) / (2.0 * wslope_max)
        enddo
      enddo

c Calculate the wind speed and direction adjustments.  The
c   curvature and wind_slope values range between -0.5 and +0.5.
c   Valid slopewt and curvewt values are between 0 and 1, with
c   values of 0.5 giving approximately equal weight to slope and
c   curvature.  I suggest that slopewt and curvewt be set such
c   that slopewt + curvewt = 1.0.  This will limit the total
c   wind weight to between 0.5 and 1.5 (but this is not required).
      do i=1,nx
        do j=1,ny
		
c Compute the wind weighting factor.
          windwt = 1.0 + slopewt * wind_slope(i,j) +
     &      curvewt * curvature(i,j)

c Generate the terrain-modified wind speed.
          windspd_grid(i,j) = windwt * windspd_grid(i,j)

c Further modify the wind speed to account for forest canopies.
          if (vegtype(i,j).le.(real(nvegtypes))) then
            if (calc_subcanopy_met.eq.1.0) then
              nveg = nint(vegtype(i,j))

c Define the canopy wind-weighting factor.  Assume z=0.6*canopy_ht,
c   and the canopy_ht equals the vegetation snow-holding depth.
              beta = 0.9
              veg_ht = vegsnowdepth(nveg)
              a = beta * forest_LAI(nveg)
              canopy_windwt = exp((- a)*(1.0 - (0.6*veg_ht)/veg_ht))
              windspd_grid(i,j) = canopy_windwt * windspd_grid(i,j)

            endif
          endif

c Modify the wind direction according to Ryan (1977).  Note that it
c   is critical that "dirdiff" handles the cases where the slope
c   azimuth and the wind direction are on different sides of the
c   360-0 line.
          if (slope_az(i,j).gt.270.0.and.
     &      winddir_grid(i,j).lt.90.0) then
            dirdiff = slope_az(i,j) - winddir_grid(i,j) - 360.0
          elseif (slope_az(i,j).lt.90.0.and.
     &      winddir_grid(i,j).gt.270.0) then
            dirdiff = slope_az(i,j) - winddir_grid(i,j) + 360.0
          else
            dirdiff = slope_az(i,j) - winddir_grid(i,j)
          endif
          if (abs(dirdiff).le.90.0) then
            winddir_grid(i,j) = winddir_grid(i,j) - 0.5 *
     &        min(wind_slope(i,j)*rad2deg,45.0) *
     &        sin(deg2rad * (2.0 * dirdiff))
            if (winddir_grid(i,j).gt.360.0) then
              winddir_grid(i,j) = winddir_grid(i,j) - 360.0
            elseif (winddir_grid(i,j).lt.0.0) then
              winddir_grid(i,j) = winddir_grid(i,j) + 360.0
            endif
          endif

c Extract the u and v wind components.
          uwind_grid(i,j) = (- windspd_grid(i,j)) *
     &      sin(deg2rad*winddir_grid(i,j))
          vwind_grid(i,j) = (- windspd_grid(i,j)) *
     &      cos(deg2rad*winddir_grid(i,j))
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine topo_data(nx,ny,deltax,deltay,topo,
     &  curvature,terrain_slope,slope_az,curve_len_scale,undef)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny,inc,h,k

      real pi,rad2deg,deltax,deltay,deltaxy,curve_len_scale,curve_max,
     &  undef,topo1,topo2,topo3,topo4,topo5,topo6,topo7,topo8

      real topo(nx_max,ny_max)
      real dzdx(nx_max,ny_max)
      real dzdy(nx_max,ny_max)
      real curvature(nx_max,ny_max)
      real slope_az(nx_max,ny_max)
      real terrain_slope(nx_max,ny_max)

	  
c Compute the topographic information required to run the wind
c   model.

c Deal with the model running at a point, or along single or double
c   lines.
      if (nx.le.2  .or.  ny.le.2) then
        do i=1,nx
          do j=1,ny
            curvature(i,j) = 0.0
            terrain_slope(i,j) = 0.0
            slope_az(i,j) = 0.0
          enddo
        enddo
      else

c Define the required constants.
        pi = 2.0 * acos(0.0)
        rad2deg = 180.0 / pi

c CURVATURE CALCULATIONS.

c Compute the average grid increment.
        deltaxy = 0.5 * (deltax + deltay)

c Convert the length scale to an appropriate grid increment.
        inc = max(1,nint(curve_len_scale/deltaxy))

c Compute the curvature.
        do i=1,nx
          do j=1,ny
                          
            k=inc
            do while(i-k.lt.1 .or. j-k.lt.1)
                k=k-1
            enddo
            do while(topo(i-k,j-k).eq.undef)
                k=k-1
            enddo
            topo1=topo(i-k,j-k)
            
            k=inc
            do while(i+k.gt.nx .or. j+k.gt.ny)
                k=k-1
            enddo
            do while(topo(i+k,j+k).eq.undef)
                k=k-1
            enddo
            topo2=topo(i+k,j+k)
            
            k=inc
            do while(i+k.gt.nx .or. j-k.lt.1)
                k=k-1
            enddo
            do while(topo(i+k,j-k).eq.undef)
                k=k-1
            enddo
            topo3=topo(i+k,j-k)
            
            k=inc
            do while(i-k.lt.1 .or. j+k.gt.ny)
                k=k-1
            enddo
            do while(topo(i-k,j+k).eq.undef)
                k=k-1
            enddo
            topo4=topo(i-k,j+k)
            
            k=inc
            do while(i+k.gt.nx)
                k=k-1
            enddo
            do while(topo(i+k,j).eq.undef)
                k=k-1
            enddo
            topo5=topo(i+k,j)
            
            k=inc
            do while(i-k.lt.1)
                k=k-1
            enddo
            do while(topo(i-k,j).eq.undef)
                k=k-1
            enddo
            topo6=topo(i-k,j)
            
            k=inc
            do while(j+k.gt.ny)
                k=k-1
            enddo
            do while(topo(i,j+k).eq.undef)
                k=k-1
            enddo
            topo7=topo(i,j+k)
            
            k=inc
            do while(j-k.lt.1)
                k=k-1
            enddo
            do while(topo(i,j-k).eq.undef)
                k=k-1
            enddo
            topo8=topo(i,j-k)

            curvature(i,j) = (4.0 * topo(i,j) -
     &        topo1 - topo2 - topo3 - topo4) /
     &        (sqrt(2.0) * 16.0 * real(inc) * deltaxy) +
     &        (4.0 * topo(i,j) -
     &        topo5 - topo6 - topo7 - topo8) /
     &        (16.0 * real(inc) * deltaxy)

          enddo
        enddo

c Scale the curvature such that the max abs(curvature) has a value
c   of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
c   divisions by zero in flat terrain where the curvature is zero.
        curve_max = 0.0 + 0.001
        do j=1,ny
          do i=1,nx
            curve_max = max(curve_max,abs(curvature(i,j)))
          enddo
        enddo
        do j=1,ny
          do i=1,nx
            curvature(i,j) = curvature(i,j) / (2.0 * curve_max)
          enddo
        enddo

c SLOPE CALCULATIONS.

c Find dzdx.
        do j=1,ny
          k=1
          do while(topo(k,j).eq.undef .and. k.lt.nx)
              k=k+1
          enddo
          if(k.lt.nx) then
            h=nx
            do while(topo(h,j).eq.undef)
                h=h-1
            enddo
            if(h.eq.k) then
                dzdx(h,j)=0
            else
                dzdx(k,j) = (topo(k+1,j) - topo(k,j)) / deltax
                if((k+1).le.(h-1)) then
                    do i=k+1,h-1
                        dzdx(i,j) = (topo(i+1,j) - topo(i-1,j)) / 
     &                              (2.0 * deltax)
                    enddo
                endif
                dzdx(h,j) = (topo(h,j) - topo(h-1,j)) / deltax
            endif
          endif
        enddo

c Find dzdy.
        do i=1,nx
          k=1
          do while(topo(i,k).eq.undef .and. k.lt.ny)
              k=k+1
          enddo
          if(k.lt.ny) then
            h=ny
            do while(topo(i,h).eq.undef)
                h=h-1
            enddo
            if(h.eq.k) then
                dzdy(i,h)=0
            else
                dzdy(i,k) = (topo(i,k+1) - topo(i,k)) / deltay
                if((k+1).le.(h-1)) then
                    do j=k+1,h-1
                        dzdy(i,j) = (topo(i,j+1) - topo(i,j-1)) / 
     &                              (2.0 * deltay)
                    enddo
                endif
                dzdy(i,h) = (topo(i,h) - topo(i,h-1)) / deltay
            endif
          endif
        enddo

c Calculate the terrain slope and azimuth.
        do i=1,nx
          do j=1,ny

c Some compilers will not allow dzdx and dzdy to both be 0.0 in
c   the atan2 computation.
c           if (abs(dzdx(i,j)).lt.1e-10) dzdx(i,j) = 1e-10
            if (abs(dzdy(i,j)).lt.1e-10) dzdy(i,j) = 1e-10

c Compute the slope azimuth, making sure that north has zero
c   azimuth.  Also note that for the Ryan wind rotation, the
c   azimuth values must range from 0 to 360.
            slope_az(i,j) = rad2deg *
     &        (3.0 / 2.0 * pi - atan2(dzdy(i,j),dzdx(i,j)))
            if (slope_az(i,j).ge.360.0) slope_az(i,j) =
     &        slope_az(i,j) - 360.0

c Compute the slope of the terrain.
            terrain_slope(i,j) = rad2deg *
     &        atan(sqrt(dzdx(i,j)*dzdx(i,j) + dzdy(i,j)*dzdy(i,j)))

          enddo
        enddo

      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine interpolate(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,var,dn,grid,undef,ifill,iobsint,
     &  iyear,imonth,iday,xhour,barnes_lg_domain,n_stns_used,
     &  k_stn)

      implicit none

      include 'snowmodel.inc'

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn  ! center x coords of lower left grid cell
      double precision ymn  ! center y coords of lower left grid cell

      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values

      double precision xstn_tmp(nstns_max) ! input stn x coords
      double precision ystn_tmp(nstns_max) ! input stn y coords
      real var_tmp(nstns_max)  ! input values

      integer nstns        ! number of input values, all good
      real undef           ! undefined value
      real dn                  ! average observation spacing
      real grid(nx_max,ny_max) ! output values

      integer i,j      ! col, row counters
      integer ifill    ! flag (=1) forces a value in every cell
      integer iobsint  ! flag (=1) use dn value from .par file

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour

      integer k_stn(nx_max,ny_max,nstns_max)      
	  integer k,n_stns_used
      real barnes_lg_domain
	  

c Use the barnes oi scheme to grid the station data. If there is
c   only a single station, distribute those data uniformly over
c   the grid.  In the event that there are no valid observations,
c   send an error message and stop (although this should have been
c   caught as part of a preprocessor step).

c The interpolation can be done two different ways:
c   First, barnes_oi does the interpolation by processing all of
c     the available station data for each model grid cell.
c   Second, barnes_oi_ij does the interpolation by processing only
c     the "n_stns_used" number of stations for each model grid cell.
c   For small domains, with relatively few met stations (100's),
c   the first way is best.  For large domains (like the
c   United States, Globe, Pan-Arctic, North America, Greenland)
c   and many met stations (like 1000's), the second approach is the
c   most efficient.  But, the second approach carries the following
c   restrictions: 1) there can be no missing data for the fields of
c   interest; 2) there can be no missing stations (all stations
c   must exist throughout the simulation period); and 3) the
c   station met file must list the stations in the same order for
c   all time steps.  In addition, the code limits the number of
c   nearest stations used to be 5 or less.
	  	  
      if (nstns.ge.2) then
        call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)

        if (barnes_lg_domain.eq.1.0) then

          do j=1,ny
            do i=1,nx

c Use that nearest station list to extract the station information
c   to be used in the interpolation.
              do k=1,n_stns_used
                xstn_tmp(k) = xstn(k_stn(i,j,k))
                ystn_tmp(k) = ystn(k_stn(i,j,k))
                var_tmp(k) = var(k_stn(i,j,k))
              enddo

c Do the interpolation for this model grid cell.
              call barnes_oi_ij(deltax,deltay,xmn,ymn,
     &          n_stns_used,xstn_tmp,ystn_tmp,var_tmp,dn,grid,
     &          undef,ifill,i,j)

            enddo
          enddo

        else

          call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,
     &      nstns,xstn,ystn,var,dn,grid,undef,ifill)

        endif

      elseif (nstns.eq.1) then

        call single_stn(nx,ny,nstns,var,grid)

      else
	  
        print *,'found no valid obs data at this time step'
        print *,'  model time =', iyear,imonth,iday,xhour
        stop

      endif
	  
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var)

      implicit none

      include 'snowmodel.inc'

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! input stn elevation
      real var(nstns_max)  ! input values
      real undef           ! undefined value
      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! input stn elevation
      real var_orig(nstns_max)  ! input values

      integer k

c Before performing the interpolation, sort through the data and
c   toss out any undefined values.
      nstns = 0
      do k=1,nstns_orig
        if (var_orig(k).ne.undef) then
          nstns = nstns + 1
          xstn(nstns) = xstn_orig(k)
          ystn(nstns) = ystn_orig(k)
          var(nstns) = var_orig(k)
          elev(nstns) = elev_orig(k)
        endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_good_values2(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,var1_orig,var2_orig,
     &  var1,var2)

c Account for the special case where you must have two coincident
c   values to do the interpolation, like Tair and rh to interpolate
c   rh, and wind speed and dir to interpolate the winds.

      implicit none

      include 'snowmodel.inc'

      integer nstns        ! number of input values, all good
      integer nstns_orig   ! number of input values
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real elev(nstns_max) ! input stn elevation
      real undef           ! undefined value
      double precision xstn_orig(nstns_max) ! input stn x coords
      double precision ystn_orig(nstns_max) ! input stn y coords
      real elev_orig(nstns_max) ! input stn elevation

      real var1_orig(nstns_max)  ! input values
      real var2_orig(nstns_max)  ! input values
      real var1(nstns_max)  ! input values
      real var2(nstns_max)  ! input values

      integer k

c Before performing the interpolation, sort through the data and
c   toss out any undefined values.
      nstns = 0
      do k=1,nstns_orig
        if (var1_orig(k).ne.undef .and. var2_orig(k).ne.undef) then
          nstns = nstns + 1
          xstn(nstns) = xstn_orig(k)
          ystn(nstns) = ystn_orig(k)
          var1(nstns) = var1_orig(k)
          var2(nstns) = var2_orig(k)
          elev(nstns) = elev_orig(k)
        endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_obs_data(nstns_orig,Tair_orig,rh_orig,xstn_orig,
     &  ystn_orig,elev_orig,iyear,imonth,iday,xhour,undef,
     &  windspd_orig,winddir_orig,prec_orig,iter)

      implicit none
	  
      include 'snowmodel.inc'	  
 
      integer iyr,imo,idy      ! year, month, and day of data
      real xhr                 ! decimal hour
      integer idstn            ! station id number

      integer k,nstns_orig,iter
      integer iyear,imonth,iday

      real Tair_orig(nstns_max),rh_orig(nstns_max)
      real winddir_orig(nstns_max),windspd_orig(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)
      real elev_orig(nstns_max),xhour,prec_orig(nstns_max)
      real undef               ! undefined value

c MicroMet assumes that the air temperature comes in as deg C, but
c   all computations must be done in K.  Check for this and make
c   an appropriate adjustment.

	  do k=1,nstns_orig
        
		if (Tair_orig(k).lt.150.0 .and. Tair_orig(k).ne.undef)
     &    Tair_orig(k) = Tair_orig(k) + 273.16
	 
	  enddo

	  return
	  end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_model_time(iyear_init,imonth_init,iday_init,
     &  xhour_init,iter,dt,iyear,imonth,iday,xhour,J_day)

      implicit none

      integer iyear,imonth,iday  ! model year, month, and day
      real xhour                 ! model decimal hour
      real dt                    ! model time step, in seconds
      integer iter               ! model iteration
      integer iyear_init         ! model start year
      integer imonth_init        ! model start month
      integer iday_init          ! model start day
      real xhour_init            ! model decimal start hour
      integer J_day              ! model day of year

c Misc. variables.
      real xmin,xmin_frac,xhour_tmp,xday
      integer imin,ihour,last

      integer lastday(12)
      data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

c Convert the simulation time to the exact year, month, day, and
c   decimal hour.

c Number of minutes.
      xmin = ((real(iter) - 1.0) * dt) / 60.0 + xhour_init * 60.0
      imin = mod(int(xmin),60)
      xmin_frac = real(imin) / 60.0

c Model integration time in decimal hours.
      xhour_tmp = xmin / 60.0
      ihour = mod(int(xhour_tmp),24)
      xhour = real(ihour) + xmin_frac

c Number of days.
      xday = xhour_tmp / 24.0
      iday = iday_init + int(xday)

c Month and year, while accounting for leap-years.
      imonth = imonth_init
      iyear = iyear_init
      
 20   continue
      last = lastday(imonth)
      if (imonth.eq.2 .and. mod(iyear,4).eq.0
     &  .and. (mod(iyear,100).ne.0 .or. mod(iyear,1000).eq.0)) then
        last = last + 1
      endif
      if (iday.gt.last) then
        iday = iday - last
        imonth = imonth + 1
        if (imonth.gt.12) then
          imonth = 1
          iyear = iyear + 1
        endif
        go to 20
      endif

c Calculate the day of year (1...365,366) corresponding to the date
c   iyear-imonth-iday.
      J_day = iday
     &  + min(1,max(0,imonth-1))*31
     &  + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4))))
     &  + min(1,max(0,imonth-3))*31
     &  + min(1,max(0,imonth-4))*30
     &  + min(1,max(0,imonth-5))*31
     &  + min(1,max(0,imonth-6))*30
     &  + min(1,max(0,imonth-7))*31
     &  + min(1,max(0,imonth-8))*31
     &  + min(1,max(0,imonth-9))*30
     &  + min(1,max(0,imonth-10))*31
     &  + min(1,max(0,imonth-11))*30
     &  + min(1,max(0,imonth-12))*31

c     print *, iyear,imonth,iday,xhour,J_day

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_lapse_rates(imonth,iday,T_lapse_rate,
     &  Td_lapse_rate,xlat,lapse_rate_user_flag,
     &  precip_lapse_rate,iprecip_lapse_rate_user_flag)

      implicit none

      integer imonth,iday,mbefore,mafter,k,months,
     &  lapse_rate_user_flag,iprecip_lapse_rate_user_flag
      real weight,T_lapse_rate,Td_lapse_rate,A,B,C,xlat,
     &  precip_lapse_rate
      parameter (months=12)
      integer lastday(months)
      data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

c The lapse rate units are in deg_C km-1.  They are converted to
c   negative_deg_C m-1 below.
      real lapse_rate(months)
      real lapse_rate_nohem(months)
      real lapse_rate_sohem(months)
      real lapse_rate_user(months)
      data lapse_rate_nohem /4.4,5.9,7.1,7.8,8.1,8.2,
     &                       8.1,8.1,7.7,6.8,5.5,4.7/
      data lapse_rate_sohem /8.1,8.1,7.7,6.8,5.5,4.7,
     &                       4.4,5.9,7.1,7.8,8.1,8.2/

c If you want to use the 'user' array, put your monthly values in
c   here and set lapse_rate_user_flag = 1 in the .par file.
      data lapse_rate_user /4.4,5.9,7.1,7.8,8.1,8.2,
     &                      8.1,8.1,7.7,6.8,5.5,4.7/

c The initial vapor pressure coeffs are in units of km-1.
      real am(months)
      real am_nohem(months)
      real am_sohem(months)
      real am_user(months)
      data am_nohem /0.41,0.42,0.40,0.39,0.38,0.36,
     &               0.33,0.33,0.36,0.37,0.40,0.40/
      data am_sohem /0.33,0.33,0.36,0.37,0.40,0.40,
     &               0.41,0.42,0.40,0.39,0.38,0.36/

c If you want to use the 'user' array, put your monthly values in
c   here and set lapse_rate_user_flag = 1 in the .par file.
      data am_user /0.41,0.42,0.40,0.39,0.38,0.36,
     &              0.33,0.33,0.36,0.37,0.40,0.40/

c The precipitation lapse rate units are in km-1.
      real prec_lapse_rate(months)
      real precip_lapse_rate_nohem(months)
      real precip_lapse_rate_sohem(months)
      real precip_lapse_rate_user(months)
      data precip_lapse_rate_nohem /0.35,0.35,0.35,0.30,0.25,0.20,
     &                              0.20,0.20,0.20,0.25,0.30,0.35/
      data precip_lapse_rate_sohem /0.20,0.20,0.20,0.25,0.30,0.35,
     &                              0.35,0.35,0.35,0.30,0.25,0.20/

c If you want to use the 'user' array, put your monthly values in
c   here and set iprecip_lapse_rate_user_flag = 1 in the .par file.
      data precip_lapse_rate_user /0.35,0.35,0.35,0.30,0.25,0.20,
     &                             0.20,0.20,0.20,0.25,0.30,0.35/

c Air and dewpoint temperature.
      do k=1,months
        if (lapse_rate_user_flag.eq.0) then
          if (xlat.lt.0.0) then
            lapse_rate(k) = lapse_rate_sohem(k)
            am(k) = am_sohem(k)
          else
            lapse_rate(k) = lapse_rate_nohem(k)
            am(k) = am_nohem(k)
          endif
        elseif (lapse_rate_user_flag.eq.1) then
          lapse_rate(k) = lapse_rate_user(k)
          am(k) = am_user(k)
        endif
      enddo

c Precipitation.
      do k=1,months
        if (iprecip_lapse_rate_user_flag.eq.0) then
          if (xlat.lt.0.0) then
            prec_lapse_rate(k) = precip_lapse_rate_sohem(k)
          else
            prec_lapse_rate(k) = precip_lapse_rate_nohem(k)
          endif
        elseif (iprecip_lapse_rate_user_flag.eq.1) then
          prec_lapse_rate(k) = precip_lapse_rate_user(k)
        endif
      enddo

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.
      A = 6.1121 * 100.0
      B = 17.502
      C = 240.97

c Over ice.
c     A = 6.1115 * 100.0
c     B = 22.452
c     C = 272.55

c Find the month before and after the day in question.
      if (iday.le.15) then
        mbefore = imonth - 1
        if (mbefore.eq.0) mbefore = 12
        mafter = imonth
        weight = (real(lastday(mbefore)) - 15. + real(iday)) /
     &    real(lastday(mbefore))
      else
        mbefore = imonth
        mafter = imonth + 1
        if (mafter.eq.13) mafter = 1
        weight = (real(iday) - 15.) / real(lastday(mbefore))
      endif

c Define the temperature lapse rate (deg C/m).
      T_lapse_rate = (- (weight * lapse_rate(mafter) +
     &  (1. - weight) * lapse_rate(mbefore))) / 1000.0

c Define the dew-point temperature lapse rate (deg C/m).
      Td_lapse_rate = (- ((weight * am(mafter) +
     &  (1. - weight) * am(mbefore)) * C)) / (B * 1000.0)

c Define the precipitation lapse rate (km-1).
      precip_lapse_rate = weight * prec_lapse_rate(mafter) +
     &  (1. - weight) * prec_lapse_rate(mbefore)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)

      implicit none

      integer nx,ny,nstns
      real deltax,deltay,dn
      real dn_max           ! the max obs spacing, dn_r
      real dn_min           ! dn_r, for large n
      integer iobsint       ! flag (=1) use dn value from .par file

c Calculate an appropriate filtered wavelength value.  First
c   calculate dn for the case of severely nonuniform data, and
c   then for the case where there is just about a station for
c   every grid cell.  Then assume that the average of these two
c   is a reasonable value to use in the interpolation.
        dn_max = sqrt(deltax*real(nx) * deltay*real(ny)) *
     &    ((1.0 + sqrt(real(nstns))) / (real(nstns) - 1.0))
        dn_min = sqrt((deltax*real(nx) * deltay*real(ny)) /
     &    real(nstns))

        if (iobsint.eq.1) then
c         dn = dn
        else
          dn = 0.5 * (dn_min + dn_max)
        endif

c       print *,'You are using an average obs spacing of',dn
c       print *,'  the program indicates a min, max range of',
c    &    dn_min,dn_max

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine barnes_oi(nx,ny,deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,var,dn,grid,undef,ifill)

c This is an implementation of the Barnes objective analysis scheme
c   as described in:
c
c   Koch, S. E., M. DesJardins, and P. J. Kocin, 1983: An
c   interactive Barnes objective map analysis scheme for use with
c   satellite and conventional data. J. Climate and Applied
c   Meteorology, 22(9), 1487-1503.

      implicit none

      include 'snowmodel.inc'

      real gamma
      parameter (gamma=0.2)
      real pi

      integer nx       ! number of x output values
      integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn !center x coords of lower left grid cell
      double precision ymn !center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values
      integer nflag        ! determines if output will be undef value
      real undef           ! undefined value

      real dn                  ! average observation spacing
      real grid(nx_max,ny_max) ! output values

      integer i,j      ! col, row counters
      integer mm,nn    ! station counters
      integer ifill    ! flag (=1) forces a value in every cell

      double precision xg,yg !temporary x and y coords of current cell
      real w1,w2       ! weights for Gauss-weighted average
      real wtot1,wtot2 ! sum of weights
      real ftot1,ftot2 ! accumulators for values, corrections
      real dsq         ! delx**2 + dely**2
      double precision xa,ya       ! x, y coords of current station
      double precision xb,yb       ! x, y coords of current station
      real dvar(nstns_max)   ! estimated error

      real xkappa_1    ! Gauss constant for first pass
      real xkappa_2    ! Gauss constant for second pass
      real rmax_1      ! maximum scanning radii, for first
      real rmax_2      ! and second passes
      real anum_1      ! numerator, beyond scanning radius,
      real anum_2      ! for first and second passes

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Compute the first and second pass values of the scaling parameter
c   and the maximum scanning radius used by the Barnes scheme.
c   Values above this maximum will use a 1/r**2 weighting.  Here I
c   have assumed a gamma value of 0.2.

c First-round values, Eqn (13).
      pi = 2.0 * acos(0.0)
      xkappa_1 = 5.052 * (2.0*dn/pi)**2

c Define the maximum scanning radius to have weight defined by
c   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
c Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
      rmax_1 = xkappa_1 * 30.0 * log(10.0)
      anum_1 = 1.0e-30 * rmax_1

c Second-round values, Eqn (4).
      xkappa_2 = gamma * xkappa_1
      rmax_2 = rmax_1 * gamma
      anum_2 = 1.0e-30 * rmax_2

c Scan each input data point and construct estimated error, dvar, at
c   that point.
      do 222 nn=1,nstns

        xa = xstn(nn)
        ya = ystn(nn)
        wtot1 = 0.0
        ftot1 = 0.0

        do 111 mm=1,nstns

          xb = xstn(mm)
          yb = ystn(mm)
          dsq = (xb - xa)**2 + (yb - ya)**2

          if (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)

          else

c Assume a 1/r**2 weight.
            w1 = anum_1/dsq

          endif

          wtot1 = wtot1 + w1
          ftot1 = ftot1 + w1 * var(mm)

  111   continue    ! end loop on sites m

        if (wtot1.eq.0.0) print *,'stn wt totals zero'

        dvar(nn) = var(nn) - ftot1/wtot1

  222 continue        ! end prediction loop on sites nn

c Grid-prediction loop.  Generate the estimate using first set of
c   weights, and correct using error estimates, dvar, and second
c   set of weights.

      do 666 j=1,ny
      do 555 i=1,nx

c xcoords of grid nodes at index i,j
c ycoords of grid nodes at index i,j
        xg = xmn + deltax * (real(i) - 1.0)
        yg = ymn + deltay * (real(j) - 1.0)

c Scan each input data point.
        ftot1 = 0.0
        wtot1 = 0.0
        ftot2 = 0.0
        wtot2 = 0.0
        nflag = 0

        do 333 nn=1,nstns
           
          xa = xstn(nn)
          ya = ystn(nn)
          dsq = (xg - xa)**2 + (yg - ya)**2

          if (dsq.le.rmax_2) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = exp((- dsq)/xkappa_2)

          elseif (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = anum_2/dsq

          else

c Assume a 1/r**2 weight.
            w1 = anum_1/dsq
            nflag = nflag + 1
c With anum_2/dsq.
            w2 = gamma * w1

          endif

          wtot1 = wtot1 + w1
          wtot2 = wtot2 + w2
          ftot1 = ftot1 + w1 * var(nn)
          ftot2 = ftot2 + w2 * dvar(nn)
           
  333   continue    ! end loop on data sites nn

        if (wtot1.eq.0.0 .or. wtot2.eq.0.0) print *,'wts total zero'

        if (ifill.eq.1) then
          grid(i,j) = ftot1/wtot1 + ftot2/wtot2
        else
          if (nflag.lt.nstns) then
            grid(i,j) = ftot1/wtot1 + ftot2/wtot2
          else
            grid(i,j) = undef
          endif
        endif

  555 continue         ! end loop on cols i
  666 continue         ! end loop on rows j

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine barnes_oi_ij(deltax,deltay,xmn,ymn,
     &  nstns,xstn,ystn,var,dn,grid,
     &  undef,ifill,i,j)

c This is an implementation of the Barnes objective analysis scheme
c   as described in:
c
c   Koch, S. E., M. DesJardins, and P. J. Kocin, 1983: An
c   interactive Barnes objective map analysis scheme for use with
c   satellite and conventional data. J. Climate and Applied
c   Meteorology, 22(9), 1487-1503.

      implicit none

      include 'snowmodel.inc'

      real gamma
      parameter (gamma=0.2)
      real pi

c     integer nx       ! number of x output values
c     integer ny       ! number of y output values
      real deltax      ! grid increment in x
      real deltay      ! grid increment in y
      double precision xmn !center x coords of lower left grid cell
      double precision ymn !center y coords of lower left grid cell

      integer nstns        ! number of input values, all good
      double precision xstn(nstns_max) ! input stn x coords
      double precision ystn(nstns_max) ! input stn y coords
      real var(nstns_max)  ! input values
      integer nflag        ! determines if output will be undef value
      real undef           ! undefined value

      real dn                  ! average observation spacing
      real grid(nx_max,ny_max) ! output values

      integer i,j      ! col, row counters
      integer mm,nn    ! station counters
      integer ifill    ! flag (=1) forces a value in every cell

      double precision xg,yg !temporary x and y coords of current cell
      real w1,w2       ! weights for Gauss-weighted average
      real wtot1,wtot2 ! sum of weights
      real ftot1,ftot2 ! accumulators for values, corrections
      real dsq         ! delx**2 + dely**2
      double precision xa,ya       ! x, y coords of current station
      double precision xb,yb       ! x, y coords of current station
      real dvar(nstns_max)   ! estimated error

      real xkappa_1    ! Gauss constant for first pass
      real xkappa_2    ! Gauss constant for second pass
      real rmax_1      ! maximum scanning radii, for first
      real rmax_2      ! and second passes
      real anum_1      ! numerator, beyond scanning radius,
      real anum_2      ! for first and second passes

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Compute the first and second pass values of the scaling parameter
c   and the maximum scanning radius used by the Barnes scheme.
c   Values above this maximum will use a 1/r**2 weighting.  Here I
c   have assumed a gamma value of 0.2.

c First-round values, Eqn (13).
      pi = 2.0 * acos(0.0)
      xkappa_1 = 5.052 * (2.0*dn/pi)**2

c Define the maximum scanning radius to have weight defined by
c   wt = 1.0 x 10**(-30) = exp(-rmax_1/xkappa_1)
c Also scale the 1/r**2 wts so that when dsq = rmax, the wts match.
      rmax_1 = xkappa_1 * 30.0 * log(10.0)
      anum_1 = 1.0e-30 * rmax_1

c Second-round values, Eqn (4).
      xkappa_2 = gamma * xkappa_1
      rmax_2 = rmax_1 * gamma
      anum_2 = 1.0e-30 * rmax_2

c Scan each input data point and construct estimated error, dvar, at
c   that point.
      do 222 nn=1,nstns

        xa = xstn(nn)
        ya = ystn(nn)
        wtot1 = 0.0
        ftot1 = 0.0

        do 111 mm=1,nstns

          xb = xstn(mm)
          yb = ystn(mm)
          dsq = (xb - xa)**2 + (yb - ya)**2

          if (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)

          else

c Assume a 1/r**2 weight.
            w1 = anum_1/dsq

          endif

          wtot1 = wtot1 + w1
          ftot1 = ftot1 + w1 * var(mm)

  111   continue    ! end loop on sites m

        if (wtot1.eq.0.0) print *,'stn wt totals zero'

        dvar(nn) = var(nn) - ftot1/wtot1

  222 continue        ! end prediction loop on sites nn

c Grid-prediction loop.  Generate the estimate using first set of
c   weights, and correct using error estimates, dvar, and second
c   set of weights.

c     do 666 j=1,ny
c     do 555 i=1,nx

c xcoords of grid nodes at index i,j
c ycoords of grid nodes at index i,j
        xg = xmn + deltax * (real(i) - 1.0)
        yg = ymn + deltay * (real(j) - 1.0)

c Scan each input data point.
        ftot1 = 0.0
        wtot1 = 0.0
        ftot2 = 0.0
        wtot2 = 0.0
        nflag = 0

        do 333 nn=1,nstns
           
          xa = xstn(nn)
          ya = ystn(nn)
          dsq = (xg - xa)**2 + (yg - ya)**2

          if (dsq.le.rmax_2) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = exp((- dsq)/xkappa_2)

          elseif (dsq.le.rmax_1) then

            w1 = exp((- dsq)/xkappa_1)
            w2 = anum_2/dsq

          else

c Assume a 1/r**2 weight.
            w1 = anum_1/dsq
            nflag = nflag + 1
c With anum_2/dsq.
            w2 = gamma * w1

          endif

          wtot1 = wtot1 + w1
          wtot2 = wtot2 + w2
          ftot1 = ftot1 + w1 * var(nn)
          ftot2 = ftot2 + w2 * dvar(nn)
           
  333   continue    ! end loop on data sites nn

        if (wtot1.eq.0.0 .or. wtot2.eq.0.0) print *,'wts total zero'

        if (ifill.eq.1) then
          grid(i,j) = ftot1/wtot1 + ftot2/wtot2
        else
          if (nflag.lt.nstns) then
            grid(i,j) = ftot1/wtot1 + ftot2/wtot2
          else
            grid(i,j) = undef
          endif
        endif

c 555 continue         ! end loop on cols i
c 666 continue         ! end loop on rows j

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine single_stn(nx,ny,nstns,var,grid)

      implicit none

      include 'snowmodel.inc'

      integer nstns    ! number of input values, all good
      integer nx       ! number of x output values
      integer ny       ! number of y output values

      real var(nstns_max)      ! input values
      real grid(nx_max,ny_max) ! output values
      integer i,j              ! col, row counters

c Assign the station value to every grid cell.
      do j=1,ny
        do i=1,nx
          grid(i,j) = var(nstns)
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pressure(nx,ny,topo,sfc_pressure,undef)

      implicit none

      include 'snowmodel.inc'

      integer i,j,nx,ny

      real topo(nx_max,ny_max),sfc_pressure(nx_max,ny_max)
      real one_atmos,scale_ht,undef

      one_atmos = 101300.0
      scale_ht = 8500.0

c Compute the average station pressure (in Pa).
      do j=1,ny
        do i=1,nx
          sfc_pressure(i,j) = one_atmos * exp((- topo(i,j))/scale_ht)
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine shortwave_data(nx,ny,deltax,deltay,xmn,ymn,
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

c This program takes observations as discrete points, compares
c   those observations with a gridded model representation of those
c   observations, at the corresponding grid cells, computes a
c   difference between the observations and modeled values, fits
c   a gridded surface through those differences, and adds the
c   difference grid to the model grid.  Thus, correcting the model
c   outputs with the observations.

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx_max,ny_max),delta_var_grid(nx_max,ny_max)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,
     &  idstn_orig,nstns_orig,i,j

c Open the observation data file.
      if (iter.eq.1) open (unit=71,file='extra_met/shortwave.dat')

c Read the data describing the time, location, and variable values
c   for each station, at this time step.  Here I have assumed that
c   the data file is in the 'non-single-station' format (with a
c   station count listed at the begining at each new time step).
      read(71,*) nstns_orig
      do k=1,nstns_orig
        read(71,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

c Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

c If there are no observational data at this time step, use the
c   modeled values without any modification.  If there are some
c   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

c For incoming shortwave, incoming longwave, and surface pressure,
c   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

      open (74,file='extra_met/shortwave_grid.gdat',
     &  form='unformatted',access='direct',recl=4*nx*ny)
      write (74,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
      close(74)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine longwave_data(nx,ny,deltax,deltay,xmn,ymn,
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

c This program takes observations as discrete points, compares
c   those observations with a gridded model representation of those
c   observations, at the corresponding grid cells, computes a
c   difference between the observations and modeled values, fits
c   a gridded surface through those differences, and adds the
c   difference grid to the model grid.  Thus, correcting the model
c   outputs with the observations.

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx_max,ny_max),delta_var_grid(nx_max,ny_max)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,
     &  idstn_orig,nstns_orig,i,j

c Open the observation data file.
      if (iter.eq.1) open (unit=72,file='extra_met/longwave.dat')

c Read the data describing the time, location, and variable values
c   for each station, at this time step.  Here I have assumed that
c   the data file is in the 'non-single-station' format (with a
c   station count listed at the begining at each new time step).
      read(72,*) nstns_orig
      do k=1,nstns_orig
        read(72,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

c Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

c If there are no observational data at this time step, use the
c   modeled values without any modification.  If there are some
c   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

c For incoming shortwave, incoming longwave, and surface pressure,
c   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

c     open (75,file='extra_met/longwave_grid.gdat',
c    &  form='unformatted',access='direct',recl=4*nx*ny)
c     write (75,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
c     close(75)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sfc_pressure_data(nx,ny,deltax,deltay,xmn,ymn,
     &  iyear,imonth,iday,xhour,undef,var_grid,iter)

c This program takes observations as discrete points, compares
c   those observations with a gridded model representation of those
c   observations, at the corresponding grid cells, computes a
c   difference between the observations and modeled values, fits
c   a gridded surface through those differences, and adds the
c   difference grid to the model grid.  Thus, correcting the model
c   outputs with the observations.

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,xhour,xhr,elev(nstns_max),undef
      real var_grid(nx_max,ny_max),delta_var_grid(nx_max,ny_max)
      real var_obs(nstns_max),elev_orig(nstns_max),var_orig(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)
      double precision xstn_orig(nstns_max),ystn_orig(nstns_max)

      integer nx,ny,k,nstns,iter,iyr,imo,idy,iyear,imonth,iday,
     &  idstn_orig,nstns_orig,i,j

c Open the observation data file.
      if (iter.eq.1) open (unit=73,file='extra_met/sfc_pressure.dat')

c Read the data describing the time, location, and variable values
c   for each station, at this time step.  Here I have assumed that
c   the data file is in the 'non-single-station' format (with a
c   station count listed at the begining at each new time step).
      read(73,*) nstns_orig
      do k=1,nstns_orig
        read(73,*) iyr,imo,idy,xhr,idstn_orig,xstn_orig(k),
     &    ystn_orig(k),elev_orig(k),var_orig(k)
      enddo

c Compare the observation time with the model time.
      if (iyr.ne.iyear .or. imo.ne.imonth .or. idy.ne.iday
     &  .or. xhr.ne.xhour) then
        print *,'model time does not match obs data input time'
        print *,'  model =', iyear,imonth,iday,xhour
        print *,'  obs   =', iyr,imo,idy,xhr
        stop
      endif

c Filter through the original input data, and eliminate any
c   missing values.
      call get_good_values1(nstns_orig,xstn_orig,ystn_orig,
     &  elev_orig,undef,nstns,xstn,ystn,elev,var_orig,var_obs)

c If there are no observational data at this time step, use the
c   modeled values without any modification.  If there are some
c   good data, do the correction/data assimilation.
      if (nstns.gt.0) then
        call DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,
     &    nstns,var_obs,delta_var_grid,var_grid)
      endif

c For incoming shortwave, incoming longwave, and surface pressure,
c   make sure no negetive numbers have been produced.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = max(0.0,var_grid(i,j))
        enddo
      enddo

c     open (76,file='extra_met/sfc_pressure_grid.gdat',
c    &  form='unformatted',access='direct',recl=4*nx*ny)
c     write (76,rec=iter) ((delta_var_grid(i,j),i=1,nx),j=1,ny)
c     close(76)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine DATA_ASSIM(nx,ny,deltax,deltay,xmn,ymn,xstn,ystn,
     &  nstns,var_obs,delta_var_grid,var_grid)

      implicit none

      include 'snowmodel.inc'

      real deltax,deltay,undef,dn
      real var_grid(nx_max,ny_max),delta_var_grid(nx_max,ny_max)
      real var_model(nstns_max),var_obs(nstns_max),
     &  delta_var(nstns_max)

      double precision xmn,ymn
      double precision xstn(nstns_max),ystn(nstns_max)

      integer ii(nstns_max),jj(nstns_max)
      integer nx,ny,i,j,ifill,iobsint,k,nstns

c Convert the x and y locations to (ii,jj) locations.
      do k=1,nstns
        ii(k) = 1 + nint((xstn(k) - xmn) / deltax)
        jj(k) = 1 + nint((ystn(k) - ymn) / deltay)
      enddo

c Extract the modeled data at the appropriate grid cells.
      do k=1,nstns
        var_model(k) = var_grid(ii(k),jj(k))
      enddo

c Calculate the difference between the modeled variable and the
c   observation at each point/grid cell.
      do k=1,nstns
        delta_var(k) = var_obs(k) - var_model(k)
      enddo

c Now that I have the differences calculated at each observation
c   point, interpolate them over the simulation domain.  Use the
c   barnes oi scheme to create the distribution. If there is
c   only a single station, distribute those data uniformly over
c   the domain.  Make sure that ifill=1, and then undef is not
c   really used (so it does not have to be the same as defined in
c   the .par file).
      undef = -9999.0
      ifill = 1
      iobsint = 0

      if (nstns.ge.2) then
        call get_dn(nx,ny,deltax,deltay,nstns,dn,iobsint)
        call barnes_oi(nx,ny,deltax,deltay,xmn,ymn,
     &    nstns,xstn,ystn,delta_var,dn,delta_var_grid,undef,ifill)
      elseif (nstns.eq.1) then
        call single_stn(nx,ny,nstns,delta_var,delta_var_grid)
      endif

c Use the gridded delta surface to correct the modeled variable.
      do j=1,ny
        do i=1,nx
          var_grid(i,j) = var_grid(i,j) + delta_var_grid(i,j)
        enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_lai(J_day,forest_LAI,vlai_winter,vlai_summer)

      implicit none

      include 'snowmodel.inc'

      integer J_day,n

      real vlai_summer(nvegtypes),vlai_winter(nvegtypes),
     &  forest_LAI(nvegtypes)

      real pi,daysinyr,tmax,tmin,peak_jday,dtseason,vtrans,tseason,
     &  fseason

c Note: A maximum forest LAI of 5.0 will give almost zero (like
c   10 W m^2) incoming solar under the canopy.  Values for Fraser
c   Experimental Forest in Colorado are 2-3 (LSOS site = 1.8,
c   Kelly's/Gus' site = 2.3).

c Calculate a seasonally varying temperature, assuming a max and
c   min temperature and a cos distribution peaking in mid July
c   (J_day = 200).  Then use this to define the seasonal lai
c   variation.
      pi = 2.0 * acos(0.0)
      daysinyr = 366.0
      tmax = 298.0
      tmin = 273.0
      peak_jday = 200.0

      dtseason = tmax - tmin
      vtrans = tmin + dtseason / 2.0

      tseason = vtrans + dtseason / 2.0 *
     &  cos(2.0 * pi / daysinyr * (real(J_day) - peak_jday))

      fseason = 0.0016 * (tmax - tseason)**2

      do n=1,nvegtypes
        forest_LAI(n) = (1.0 - fseason) * vlai_summer(n) +
     &    fseason * vlai_winter(n)
      enddo

c      print *,J_day,(forest_LAI(n),n=1,nvegtypes)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

