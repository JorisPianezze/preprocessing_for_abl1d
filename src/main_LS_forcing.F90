PROGRAM main
   !!======================================================================
   !!                     ***  PROGRAM main  ***
   !!
   !! ** Purpose : compute geostrophic wind or horizontal pressure gradient
   !!              from ECMWF atmospheric model vertical levels altitude, 
   !!              temperature and humidity 3D fields
   !! 
   !!======================================================================
   !! History : 2016-10  (F. Lemarié)  Original code
   !!   
   !!----------------------------------------------------------------------
   USE module_io       ! I/O routines
   USE module_grid     ! compute input and output grids 
   !!
   IMPLICIT NONE
   !!----------------------------------------------------------------------
   !! 
   !!  
   !! 
   !!----------------------------------------------------------------------
   !
   INTEGER                                  :: ji,jj,jk,kt, nhym, nhyi
   INTEGER                                  :: jpka                  ! number of vertical levels for input and target grids 
   INTEGER                                  :: jpi , jpj             ! number of grid points in x and y directions     
   INTEGER                                  :: status             
   INTEGER                                  :: jptime,ctrl
   INTEGER                                  :: ioerr,ncid
   INTEGER, ALLOCATABLE, DIMENSION(:,:  )   :: ind  
   INTEGER, PARAMETER                       :: stdout  = 6
   !!
   REAL(8)                                  :: cff,tv
   !!
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: lev                ! A coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: A_w                ! A coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: B_w                ! B coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: B_r                ! B coefficients to reconstruct ECMWF grid   
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: A_r                ! A coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: ph,lon,lat   
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:  ) :: e3t,ghw                ! thickness of vertical layers in target grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: tmp1d, tmp_fullw, tmp_fullm ! temporary/working 1D arrays
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: humi
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: temp
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: uhpg
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: vhpg
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: slp, zsurf
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: dx,dy,ff_t,tmask,tmask2
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: FX,FE,wrkx,wrke
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: dZx,dZe
   !!
   CHARACTER(len=500)                       :: file_u,file_v,file_hpg, file_geos ! ECMWF files containing wind components 
   CHARACTER(len=500)                       :: file_t,file_q, file_m             ! ECMWF files containing tracers and mask 
   CHARACTER(len=500)                       :: file_z,file_p,cn_dir              ! ECMWF files containing surface geopot and pressure
   CHARACTER(len=500)                       :: out_file
   CHARACTER(len=500)                       :: slp_smth_method ! Method used to smooth the slp : choice between 'hanning' or 'laplacian'
   CHARACTER(len=500)                       :: namelistf
   CHARACTER(len=500)                       :: argument   
   CHARACTER(len= 20),DIMENSION(4)          :: dimnames
   CHARACTER(len= 20),DIMENSION(4)          :: varnames
   CHARACTER(6)                             :: mask_var           ! name of mask variable in file_m file  
   !!
   LOGICAL                                  :: ln_read_zsurf      ! read surface geopotential or not
   LOGICAL                                  :: ln_read_mask       ! read land-sea mask or not
   LOGICAL                                  :: ln_perio_latbc     ! use periodic BC along the domain latitudinal edges (for global data) or use zero-gradient BC (for regional data)
   LOGICAL                                  :: ln_c1d             ! output only a single column in output file
   LOGICAL                                  :: ln_hpg_frc         ! compute horizontal pressure gradient
   LOGICAL                                  :: ln_geo_wnd         ! compute goestrophic wind components 
   LOGICAL                                  :: ln_slp_smth        ! apply gibbs oscillation filetring on mean sea level pressure
   LOGICAL                                  :: ln_slp_log         ! log(sea-level pressure) or sea-level pressure
   LOGICAL                                  :: ln_lsm_land        ! if T mask is 1 over land and 0 over ocean if F it is the other way around
   LOGICAL                                  :: ln_impose_z1       ! impose the altitude of the first level in target grid
   INTEGER                                  :: ptemp_method       ! way to compute potential temperature 
                                                                  ! = 0  (absolute temperature)
                                                                  ! = 1  (potential temperature with local ref pressure)
                                                                  ! = 2  (potential temperature with global ref pressure on temperature perturbation)
                                                                  ! = 3  (potential temperature with global ref pressure)
                                                                  ! = 4  (potential with a constant global ref pressure)
   INTEGER                                  :: slp_smth_iter      ! Number of iterations of the smoothing filter  
   !!
   REAL(8), PARAMETER     :: grav  =   9.80665
   REAL(8), PARAMETER     :: Rd    = 287.058   
   REAL(8), PARAMETER     :: zvir  =   0.609133 
   REAL(8), PARAMETER     :: omega =   7.292116e-05
   REAL(8), PARAMETER     :: rad   =   3.141592653589793 / 180.
   REAL(8), PARAMETER     :: rt    = 6371229.
   
   !!---------------------------------------------------------------------
   !! List of variables read in the namelist file 
   NAMELIST/nml_opt/    ptemp_method, ln_slp_log, ln_slp_smth, slp_smth_method, slp_smth_iter, ln_read_mask, ln_perio_latbc, & 
                        ln_hpg_frc, ln_geo_wnd, ln_c1d, ln_read_zsurf, ln_lsm_land
   NAMELIST/nml_fld/    cn_dir, file_u, file_v, file_t,               &
      &                 file_q, file_z, file_p, file_hpg, file_geos,  &
      &                 file_m, mask_var
   !!

   WRITE(stdout,*) '                                                           '
   WRITE(stdout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(stdout,*) '  ~~~               get_atm_LSfrc.exe                 ~~~  '
   WRITE(stdout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

   !! get the namelist file name
   CALL get_command_argument(1, argument, ctrl, status)
   !
   SELECT CASE(status)
   CASE(0)
      namelistf = trim(argument)
   CASE(-1)
      WRITE(stdout,*) "### Error: file name too long"
      STOP
   CASE DEFAULT
      namelistf = 'namelist_abl_tools'
   END SELECT
   !!---------------------------------------------------------------------

   !!---------------------------------------------------------------------
   !! read namelist variables
   ctrl = 0
   OPEN(50, file=namelistf, status='old', form='formatted', access='sequential', iostat=ioerr)
   IF (ioerr /= 0) ctrl = ctrl + 1   
   READ(50,nml_opt, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1   
   READ(50,nml_fld, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   
   IF (ctrl > 0) then
      WRITE(stdout,*) "### E R R O R while reading namelist file '",trim(namelistf),"'"
      WRITE(stdout,*) " ctrl = ",ctrl
      STOP
   ELSE 
      WRITE(stdout,*) " Namelist file ",trim(namelistf), " OK "    
   END IF

   IF(ln_hpg_frc) THEN
      WRITE(stdout,*) "### E R R O R in namelist variable "
      WRITE(stdout,*) "ln_hpg_frc not yet implemented in CROCO, please use ln_geo_wnd=.true. and ln_hpg_frc=.false."
      STOP
      out_file = trim(cn_dir)//'/'//trim(file_hpg)
   ELSE IF(ln_geo_wnd) THEN
      out_file = trim(cn_dir)//'/'//trim(file_geos)   
   ELSE
      WRITE(stdout,*) "### E R R O R in namelist variable "
      WRITE(stdout,*) "either ln_hpg_frc or ln_geo_wnd should be set to True"
      STOP       
   END IF

   IF ( ln_slp_smth .and. ( ( trim(slp_smth_method) .ne. 'laplacian' ) .and. ( trim(slp_smth_method) .ne. 'hanning' )  ) ) THEN
      WRITE(stdout,*) " ### E R R O R please define slp_smth_method = 'laplacian' or 'hanning' if you want to smooth slp field"
      STOP
   ELSE
      WRITE(stdout,*) " Smoothing of SLP activated with ", slp_smth_iter, "iteration(s) of ", trim(slp_smth_method), " filter"
   END IF
   
   !!---------------------------------------------------------------------   


   !!---------------------------------------------------------------------
   ! check files content
   ctrl = 0
   CALL Read_Ncdf_dim('lev'  ,trim(cn_dir)//'/'//trim(file_t),jpka  )
   !                            geop_surf  temp  humi  mslp
   varnames = [character(len=6) :: 'z',    't',  'q', 'msl']
   !
   IF ( jpka .le. 1 ) ctrl = ctrl + 1 
   IF (ln_read_zsurf) THEN
     IF( .not. VAR_EXISTENCE( trim(varnames(1)) , trim(cn_dir)//'/'//trim(file_z) ) ) ctrl = ctrl + 1
   END IF
   IF (  .not. VAR_EXISTENCE( trim(varnames(2)) , trim(cn_dir)//'/'//trim(file_t) ) ) ctrl = ctrl + 1 
   IF (  .not. VAR_EXISTENCE( trim(varnames(3)) , trim(cn_dir)//'/'//trim(file_q) ) ) ctrl = ctrl + 1   
   IF (  .not. VAR_EXISTENCE( trim(varnames(4)) , trim(cn_dir)//'/'//trim(file_p) ) ) ctrl = ctrl + 1
       
   IF ( ctrl > 0 ) THEN  
      WRITE(stdout,*) "### E R R O R while reading ECMWF atmospheric files "
      STOP
   ELSE
      WRITE(stdout,*) " ECMWF atmospheric files OK "        
   END IF
   !!---------------------------------------------------------------------


   !!---------------------------------------------------------------------
   !! read the dimensions for the input files
   CALL Read_Ncdf_dim ( 'time', trim(cn_dir)//'/'//trim(file_t), jptime  ) 
   CALL Read_Ncdf_dim ( 'lon' , trim(cn_dir)//'/'//trim(file_t), jpi     )     
   CALL Read_Ncdf_dim ( 'lat' , trim(cn_dir)//'/'//trim(file_t), jpj     )     
   CALL Read_Ncdf_dim ( 'nhym', trim(cn_dir)//'/'//trim(file_t), nhym    )
   CALL Read_Ncdf_dim ( 'nhyi', trim(cn_dir)//'/'//trim(file_t), nhyi    )
   !WRITE(stdout,*) "jpka, jptime, jpi, jpj, nhym, nhyi: ", jpka, jptime, jpi, jpj, nhym, nhyi
   !
   !!---------------------------------------------------------------------


   !!---------------------------------------------------------------------
   !! allocate arrays  
   ALLOCATE( A_w    (               0:jpka   ) )         
   ALLOCATE( B_w    (               0:jpka   ) )     
   ALLOCATE( B_r    (               1:jpka   ) )
   ALLOCATE( A_r    (               1:jpka   ) )   

   ALLOCATE( e3t    ( 1:jpi, 1:jpj, 1:jpka   ) )
   ALLOCATE( ghw    ( 1:jpi, 1:jpj, 0:jpka   ) ) 
   ALLOCATE( slp    ( 1:jpi, 1:jpj           ) )
   ALLOCATE( zsurf  ( 1:jpi, 1:jpj           ) )
   ALLOCATE( temp   ( 1:jpi, 1:jpj, 1:jpka, 1) )
   ALLOCATE( humi   ( 1:jpi, 1:jpj, 1:jpka, 1) )   
   ALLOCATE( uhpg   ( 1:jpi, 1:jpj, 1:jpka, 1) )
   ALLOCATE( vhpg   ( 1:jpi, 1:jpj, 1:jpka, 1) )   
   ALLOCATE( ph     (               0:jpka   ) ) 
   ALLOCATE( dx     ( 1:jpi, 1:jpj           ) )
   ALLOCATE( dy     ( 1:jpi, 1:jpj           ) )
   ALLOCATE( FX     ( 1:jpi, 1:jpj           ) ) 
   ALLOCATE( FE     ( 1:jpi, 1:jpj           ) ) 
   ALLOCATE( dzx    ( 0:jpi, 1:jpj           ) ) 
   ALLOCATE( dze    ( 1:jpi, 0:jpj           ) ) 
   ALLOCATE( wrkx   ( 1:jpi, 1:jpj           ) )      
   ALLOCATE( wrke   ( 1:jpi, 1:jpj           ) )  
   ALLOCATE( tmask  ( 1:jpi, 1:jpj           ) ) 
   ALLOCATE( tmask2 ( 1:jpi, 1:jpj           ) ) 
   IF (jpka.NE.nhym) THEN
     ALLOCATE( tmp_fullw(1:nhyi) )
     ALLOCATE( tmp_fullm(1:nhym) )
   END IF

   !!---------------------------------------------------------------------
   !! Read the mask and remove some closed seas 
   IF (ln_read_mask) THEN
     CALL init_atm_mask( jpi, jpj, trim(cn_dir)//'/'//trim(file_m), trim(mask_var), ln_lsm_land, tmask)
   ELSE
     tmask(:,:) = 1.
   END IF
   tmask2 (:,:) = 1.
   !!


   !! Read the static A and B coefficients for the ECMWF vertical grid
   IF (jpka.EQ.nhym) THEN
     CALL Read_Ncdf_var ( 'hyai', trim(cn_dir)//'/'//trim(file_t), A_w )
     CALL Read_Ncdf_var ( 'hybi', trim(cn_dir)//'/'//trim(file_t), B_w )   
     CALL Read_Ncdf_var ( 'hyam', trim(cn_dir)//'/'//trim(file_t), A_r ) 
     CALL Read_Ncdf_var ( 'hybm', trim(cn_dir)//'/'//trim(file_t), B_r ) 
   ELSE
     CALL Read_Ncdf_var ( 'hyai', trim(cn_dir)//'/'//trim(file_t), tmp_fullw )
     A_w(0:jpka) = tmp_fullw(nhyi-(jpka+1)+1:nhyi)
     CALL Read_Ncdf_var ( 'hybi', trim(cn_dir)//'/'//trim(file_t), tmp_fullw )
     B_w(0:jpka) = tmp_fullw(nhyi-(jpka+1)+1:nhyi)
     CALL Read_Ncdf_var ( 'hyam', trim(cn_dir)//'/'//trim(file_t), tmp_fullm )
     A_r(1:jpka) = tmp_fullm(nhym-jpka+1:nhym)
     CALL Read_Ncdf_var ( 'hybm', trim(cn_dir)//'/'//trim(file_t), tmp_fullm )
     B_r(1:jpka) = tmp_fullm(nhym-jpka+1:nhym)
   END IF

   ALLOCATE( lat(1:jpj), lon(1:jpi), ff_t(1:jpi,1:jpj), lev(1:jpka) )
   CALL Read_Ncdf_var ( 'lon' , trim(cn_dir)//'/'//trim(file_t), lon )      !<-- longitude
   CALL Read_Ncdf_var ( 'lat' , trim(cn_dir)//'/'//trim(file_t), lat )      !<-- latitude  
   CALL Read_Ncdf_var ( 'lev' , trim(cn_dir)//'/'//trim(file_t), lev )

   !!---------------------------------------------------------------------
   !++ Compute Coriolis frequency at cell centers
   !
   DO jj = 1, jpj
      DO ji = 1, jpi
         ff_t(ji,jj) = 2. * omega * SIN( rad * lat(jj) )
      END DO
   END DO

   !!---------------------------------------------------------------------
   !++ Compute dx and dy at cell centers
   !    
   DO jj = 2, jpj-1
      DO ji = 1, jpi-1
         dx(ji,jj) = rt * rad * ( lon(ji+1) - lon(ji) ) * COS( rad * lat(jj) )
      END DO
   END DO 
   dx(jpi,1:jpj) = dx(jpi-1,1:jpj)
   dx(1:jpi,jpj) = dx(1:jpi,jpj-1)
   dx(1:jpi,  1) = dx(1:jpi,  2)
   !++ 
   DO jj = 1, jpj-1
      DO ji = 1, jpi   
         dy(ji,jj) = rt * rad * abs( lat(jj+1) - lat(jj) )
      END DO
   END DO   
   dy(1:jpi,jpj) = dy(1:jpi,jpj-1) 


   !!---------------------------------------------------------------------   
   !! create output file
   !!
   status = nf90_create ( trim(out_file), NF90_WRITE, ncid )
   status = nf90_close  ( ncid )
   !
   CALL Write_Ncdf_dim ( 'lon' , trim(out_file), jpi    )
   CALL Write_Ncdf_dim ( 'lat' , trim(out_file), jpj    )
   CALL Write_Ncdf_dim ( 'lev' , trim(out_file), jpka   )
   CALL Write_Ncdf_dim ( 'nhym', trim(out_file), jpka   )
   CALL Write_Ncdf_dim ( 'nhyi', trim(out_file), jpka+1 )
   CALL Write_Ncdf_dim ( 'time', trim(out_file), 0      )
   !
   CALL Write_Ncdf_var ( 'lon' , 'lon' , trim(out_file), lon, 'double' )
   CALL Write_Ncdf_var ( 'lat' , 'lat' , trim(out_file), lat, 'double' )
   CALL Write_Ncdf_var ( 'lev' , 'lev' , trim(out_file), lev, 'double' )
   CALL Write_Ncdf_var ( 'hyai', 'nhyi', trim(out_file), A_w, 'double' )
   CALL Write_Ncdf_var ( 'hybi', 'nhyi', trim(out_file), B_w, 'double' )
   CALL Write_Ncdf_var ( 'hyam', 'nhym', trim(out_file), A_r, 'double' )
   CALL Write_Ncdf_var ( 'hybm', 'nhym', trim(out_file), B_r, 'double' )

   !!---------------------------------------------------------------------   
   !! Initialize the name of the dimensions for geostrophic winds in the output file
   !!
   dimnames(1) = 'lon'
   dimnames(2) = 'lat'       
   dimnames(3) = 'lev'
   dimnames(4) = 'time'
   CALL Write_Ncdf_var ( 'tmask', dimnames(1:2), trim(out_file), tmask, 'float' )

   !!---------------------------------------------------------------------
   ! Read time variable
   ALLOCATE(tmp1d (1:jptime))    
   CALL Read_Ncdf_var ( 'time', trim(cn_dir)//'/'//trim(file_t), tmp1d )    
   !!---------------------------------------------------------------------     

   DO kt=1,jptime  
   !
      WRITE(stdout,*) ' Treat time = ',kt,'/',jptime
      !
      CALL Write_Ncdf_var( 'time', dimnames(4:4), trim(out_file), tmp1d(kt:kt), kt, 'double' )
      !
      IF( kt == 1 ) THEN
         CALL Duplicate_lon_lat_time( trim(cn_dir)//'/'//trim(file_t), out_file )
         CALL Duplicate_lev_hyb     ( trim(cn_dir)//'/'//trim(file_t), out_file ) 
      ENDIF      
      !
      CALL Read_Ncdf_var ( varnames(4) , trim(cn_dir)//'/'//trim(file_p), slp,  kt )      !<-- (log of) surface pressure / msl
      !
      IF (ln_slp_log) THEN
        DO jj = 1, jpj
           DO ji = 1, jpi      
             slp(ji,jj) = exp( slp(ji,jj) )
           END DO
        END DO
      ENDIF
      !
      IF (ln_read_zsurf) THEN
        CALL Read_Ncdf_var ( varnames(1) , trim(cn_dir)//'/'//trim(file_z), zsurf,  kt )      !<-- surface geopotential 
      ELSE
        zsurf(:,:) = 0.
      END IF
      !    
      CALL Read_Ncdf_var ( varnames(2) , trim(cn_dir)//'/'//trim(file_t), temp,  kt    )      !<-- temperature   
      CALL Read_Ncdf_var ( varnames(3) , trim(cn_dir)//'/'//trim(file_q), humi,  kt    )      !<-- humidity        
      !
      ! Smoothing of surface fields to remove gibbs oscillations (must be done on both fields or none of them)
      IF( ln_slp_smth ) CALL smooth_field ( jpi, jpj, slp (:,:), tmask (:,:), slp_smth_method, slp_smth_iter )
      !
      ! Compute horizontal gradient of slp in x-direction (FX = dslp / dx)
      ! ++ Add initialization
      FX(:,:) = 0.0
      !
      DO jj = 1, jpj
         DO ji = 2, jpi-1     
            ! ++ JORIS ++
            ! Ajout dun test sur le mask
            IF ( tmask(ji,jj) == 1.0 .AND. tmask(ji+1,jj) == 1.0 ) THEN
              cff            = 2. / ( dx(ji+1,jj) + dx(ji,jj) )
              FX( ji ,jj )   = cff * ( slp( ji+1,jj ) - slp( ji,jj ) )
            ELSE IF ( tmask(ji,jj) == 0.0 .OR. (tmask(ji+1,jj) == 0.0 .AND. tmask(ji-1,jj) == 0.0) ) THEN
              FX( ji ,jj )   = 0.0
            ELSE IF ( tmask(ji,jj) == 1.0 .AND. tmask(ji+1,jj) == 0.0 ) THEN
              FX( ji ,jj )   = FX( ji-1 ,jj )
            END IF
            ! -- JORIS --
         END DO
      END DO 

      IF (ln_perio_latbc) THEN
        ! apply periodicity       
        DO jj = 1, jpj
           cff             = 2./( dx(1,jj)+dx(jpi,jj) )
           FX( jpi ,jj )   = cff*( slp( 1,jj ) - slp( jpi,jj ) )            
        END DO    
      ELSE
        ! apply no-gradient
        DO jj = 1, jpj
           FX( jpi ,jj )   = FX( jpi-1 ,jj )
           FX(   1 ,jj )   = FX(     2 ,jj )
        END DO
      ENDIF    

      CALL Write_Ncdf_var( 'dslpdx', dimnames(1:2), trim(out_file), FX*tmask, 'float' )

      ! Compute horizontal gradient of slp in y-direction (FE = dslp / dy)
      ! ++ Add initialization
      FE(:,:) = 0.0
      !
      DO jj = 2, jpj-1
         DO ji = 1, jpi     
            ! ++ JORIS ++ 
            ! Ajout dun test sur le mask
            IF ( tmask(ji,jj) == 1.0 .AND. tmask(ji,jj+1) == 1.0 ) THEN
              cff            = 2./( dy(ji,jj+1)+dy(ji,jj) )
              FE( ji ,jj )   = cff*( slp( ji,jj+1 ) - slp( ji,jj ) )
            ELSE IF ( tmask(ji,jj) == 0.0 .OR. (tmask(ji,jj+1) == 0.0 .AND. tmask(ji,jj-1) == 0.0) ) THEN
              FE( ji ,jj )   = 0.0
           ELSE IF ( tmask(ji,jj) == 1.0 .AND. tmask(ji,jj+1) == 0.0 ) THEN
              FE( ji ,jj )   = FE( ji ,jj-1 )
            END IF
            ! -- JORIS --
         END DO
      END DO   
      ! apply no-gradient     
      DO ji = 1, jpi
         FE( ji ,jpj )   = FE( ji ,jpj-1 )          
         FE( ji ,  1 )   = FE( ji ,    2 )          
      END DO              

      CALL Write_Ncdf_var( 'dslpdy', dimnames(1:2), trim(out_file), FE*tmask, 'float' )

      ! Compute the altitude at layer interfaces
      ghw(:,:,1:jpka)  = 0.
      ghw(:,:,  jpka)  = zsurf(:,:) * (1. / grav)
      DO jj = 1, jpj
         DO ji = 1, jpi

            DO jk=0,jpka
               ph(jk) = A_w( jk ) + B_w( jk ) * slp( ji, jj )    !<-- Pa
            END DO

            !ph(0) = 0.1
            IF ( nhym .EQ. jpka) ph(0) = 1.

            DO jk=jpka,1,-1
               tv      = temp( ji, jj, jk, 1 ) * ( 1. + zvir*humi( ji, jj, jk, 1 ) )  !<-- Virtual temperature
               e3t ( ji, jj, jk   ) =  (1./grav)*( Rd * tv * log( ph( jk ) / ph( jk-1 ) ) )
               ghw ( ji, jj, jk-1 ) =  e3t( ji, jj, jk ) + ghw( ji, jj, jk )
            END DO  

         END DO
      END DO          
      !
      !++ Compute the geostrophic winds
      !
      !////////////   
      DO jk=1,jpka
      !////////////           

         ! Compute horizontal gradient of altitude in x-direction along the coordinate dZx = (dz / dx)s
         DO jj = 1, jpj
            DO ji = 1, jpi-1  
               cff          = 1./(dx(ji+1,jj)+dx(ji,jj))
               dZx(ji,jj)   = cff*(  (ghw( ji+1, jj, jk-1) - ghw( ji, jj, jk-1))   &
                  &                 +(ghw( ji+1, jj, jk  ) - ghw( ji, jj, jk  )) )
            END DO
         END DO   

         IF (ln_perio_latbc) THEN
           ! apply periodicity
           ji=jpi
           DO jj = 1, jpj
              cff          = 1./(dx(1,jj)+dx(jpi,jj)) 
              dZx(ji,jj)   = cff*( (ghw( 1, jj, jk-1) - ghw( jpi, jj, jk-1))   &
                 &                +(ghw( 1, jj, jk  ) - ghw( jpi, jj, jk  )) )
              dZx( 0,jj)   = dZx(jpi,jj)                
           END DO 
         ELSE
           ! apply no-gradient
           DO jj = 1, jpj
              dZx(jpi,jj)     = dZx(jpi-1,jj)
              dZx(  0,jj)     = dZx(    1,jj)            
           END DO         
         END IF
         
         ! Compute horizontal gradient of altitude in y-direction along the coordinate dZy = (dz / dy)s
         DO jj = 1, jpj-1
            DO ji = 1, jpi     
               cff          = 1./(dy(ji,jj)+dy(ji,jj+1))
               dZe(ji,jj)   = cff*(  (ghw( ji, jj+1, jk-1) - ghw( ji, jj, jk-1))  &
                  &                 +(ghw( ji, jj+1, jk  ) - ghw( ji, jj, jk  )) )   
            END DO
         END DO 
         ! apply no-gradient 
         DO ji = 1, jpi
            dZe(ji,jpj)     = dZe(ji,jpj-1)
            dZe(ji,  0)     = dZe(ji,    1)            
         END DO              
         
         ! Compute horizontal pressure gradient in x-direction along the coordinate wrkX = (dp/dx)s
         DO jj = 1, jpj
            DO ji = 2, jpi       
               cff = slp(ji,jj)*(B_w(jk)-B_w(jk-1))+(A_w(jk)-A_w(jk-1))
               wrkX(ji,jj) = B_r(jk)*0.5*(FX(ji,jj)+FX(ji-1,jj))  &
                  &               *(ghw(ji,jj,jk)-ghw(ji,jj,jk-1))/cff
            END DO
         END DO

         IF (ln_perio_latbc) THEN
           ! apply periodicity
           ji = 1
           DO jj = 1, jpj
                 cff = slp(ji,jj)*(B_w(jk)-B_w(jk-1))+(A_w(jk)-A_w(jk-1))
                 wrkX(ji,jj) = B_r(jk)*0.5*(FX(ji,jj)+FX(jpi,jj))  &
                    &               *(ghw(ji,jj,jk)-ghw(ji,jj,jk-1))/cff            
           END DO         
         ELSE
           ! apply no gradient  
           DO jj = 1, jpj
              wrkX(1,jj) = wrkX(2,jj)
           END DO     
         END IF 
      
         ! Compute horizontal pressure gradient in y-direction along the coordinate wrkE = (dp/dy)s
         DO jj = 2, jpj
            DO ji = 1, jpi       
               cff = slp(ji,jj)*(B_w(jk)-B_w(jk-1))+(A_w(jk)-A_w(jk-1))
               wrkE(ji,jj) = B_r(jk)*0.5*(FE(ji,jj)+FE(ji,jj-1))  &
                  &               *(ghw(ji,jj,jk)-ghw(ji,jj,jk-1))/cff
            END DO
         END DO    
         ! apply no gradient  
         jj = 1
         DO ji = 1, jpi       
            wrkE(ji,1) = wrkE(ji,2) 
         END DO   

         !+++ Finalize pressure gradient computation
         IF(ln_hpg_frc) THEN
            DO jj=1,jpj
               DO ji=1,jpi                     
               uhpg(ji,jj,jk,1) = - grav*( wrkX(ji,jj) - 0.5*( dZx(ji,jj)+dZx(ji-1,jj) ) ) 
               vhpg(ji,jj,jk,1) =   grav*( wrkE(ji,jj) - 0.5*( dZe(ji,jj)+dZe(ji,jj-1) ) )
               END DO   
            END DO   
         ELSE
            DO jj=1,jpj               
               DO ji=1,jpi    
               cff = grav/ff_t(ji,jj)
               if(abs(ff_t(ji,jj)) < 3.e-5) cff = 0.    
               ! minus sign for uhpg because y-derivatives are inverted
               vhpg(ji,jj,jk,1) = - cff*( wrkX(ji,jj) - 0.5*( dZx(ji,jj)+dZx(ji-1,jj) ) ) 
               uhpg(ji,jj,jk,1) = - cff*( wrkE(ji,jj) - 0.5*( dZe(ji,jj)+dZe(ji,jj-1) ) )
               END DO   
            END DO          
         ENDIF
      !////////////        
      END DO
      !////////////        

      CALL Write_Ncdf_var ( 'uhpg', dimnames(1:4), trim(out_file), uhpg, kt, 'float' )
      CALL Write_Ncdf_var ( 'vhpg', dimnames(1:4), trim(out_file), vhpg, kt, 'float' )
      ! DEBUG:
      !uhpg(:,:,:,1) = ghw(:,:,:)
      !vhpg(:,:,:,1) = e3t(:,:,:)
      !CALL Write_Ncdf_var ( 'ghw', dimnames(1:4), trim(out_file), uhpg, kt, 'float' )
      !CALL Write_Ncdf_var ( 'e3t', dimnames(1:4), trim(out_file), vhpg, kt, 'float' )
   !
   END DO ! kt
   !
   !
   !DEALLOCATE(zsurf,slp,temp,humi,uhpg,vhpg)
   DEALLOCATE(zsurf,slp,temp,uhpg,vhpg)
   IF (jpka.NE.nhym) DEALLOCATE(tmp_fullw,tmp_fullm)
   !
   STOP
   !
END PROGRAM main
