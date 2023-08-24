PROGRAM main
   !!======================================================================
   !!                     ***  PROGRAM main  ***
   !!
   !! ** Purpose : Vertical interpolation of ECMWF dataset on a given fixed 
   !!              vertical grid
   !! 
   !!======================================================================
   !! History : 2016-10  (F. Lemarié)  Original code
   !!   
   !!----------------------------------------------------------------------
   USE module_io       ! I/O routines
   USE module_interp   ! vertical interpolation routines
   USE module_grid     ! compute input and output grids 
   !!
   IMPLICIT NONE
   !!----------------------------------------------------------------------
   !! 
   !!  
   !! 
   !!----------------------------------------------------------------------
   !
   INTEGER                                  :: ji,jj,jk,kt, jk_in, nhym, nhyi
   INTEGER                                  :: jpka_in, jpka      ! number of vertical levels for input and target grids 
   INTEGER                                  :: jpi , jpj          ! number of grid points in x and y directions     
   INTEGER                                  :: iloc, jloc         ! grid indexes for c1d case
   INTEGER                                  :: status             
   INTEGER                                  :: jptime,ctrl
   INTEGER                                  :: ioerr
   INTEGER, ALLOCATABLE, DIMENSION(:,:  )   :: ind  
   INTEGER, PARAMETER                       :: stdout  = 6
   INTEGER, PARAMETER                       :: jp_weno = 1
   INTEGER, PARAMETER                       :: jp_spln = 2   
   !!
   REAL(8)                                  :: hc,hmax,theta_s,z1 ! parameters related to the target vertical grid
   REAL(8)                                  :: cff
   !!  
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: A_w                ! A coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: A_wa               ! A coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: B_w                ! B coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: B_wa               ! B coefficients to reconstruct ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: tmp1d, tmp_fullw, tmp_fullm ! temporary/working 1D arrays
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: e3t,e3w            ! thickness of vertical layers in target grid
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: ght,ghw            ! altitude of vertical grid points 
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: e3_bak
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:  ) :: ghw_in             ! altitude of cell interfaces  of ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:  ) :: e3t_in             ! thickness of vertical layers in ECMWF grid
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: humi
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: temp
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: varout, varc1d 
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: slp, zsurf
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: tmask              ! land-sea mask  
   !!
   CHARACTER(len=500)                       :: file_u,file_v,file_hpg,file_geos  ! ECMWF files containing wind components 
   CHARACTER(len=500)                       :: file_t,file_q, file_m             ! ECMWF files containing tracers and mask 
   CHARACTER(len=500)                       :: file_z,file_p,cn_dir,file_in      ! ECMWF files containing surface geopot and pressure
   CHARACTER(len=500)                       :: grd_file, abl_file, drwn_file, out_file
   CHARACTER(len=500)                       :: slp_smth_method    ! Method used to smooth the slp : choice between 'hanning' or 'laplacian'
   CHARACTER(len=500)                       :: namelistf,stmp
   CHARACTER(len=500)                       :: argument, var_file
   CHARACTER(len= 20),DIMENSION(4)          :: dimnames
   CHARACTER(len= 20),DIMENSION(9)          :: varnames, outnames
   CHARACTER(len=500),DIMENSION(9)          :: filnames
   CHARACTER(6)                             :: mask_var           ! name of mask variable in file_m file            
   CHARACTER(6)                             :: var_name
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

   !!---------------------------------------------------------------------
   !! List of variables read in the namelist file 
   NAMELIST/nml_dom/    jpka, hmax, theta_s, hc, ln_impose_z1, z1   
   NAMELIST/nml_out/    grd_file, abl_file, drwn_file, var_name
   NAMELIST/nml_c1d/    iloc, jloc
   NAMELIST/nml_opt/    ptemp_method, ln_slp_log, ln_slp_smth, slp_smth_method, slp_smth_iter, ln_read_mask, ln_perio_latbc, &
                        ln_hpg_frc, ln_geo_wnd, ln_c1d, ln_read_zsurf, ln_lsm_land
   NAMELIST/nml_fld/    cn_dir, file_u, file_v, file_t,               &
      &                 file_q, file_z, file_p, file_hpg, file_geos,  &
      &                 file_m, mask_var

   WRITE(stdout,*) '                                                           '
   WRITE(stdout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(stdout,*) '  ~~~              vinterp_atm_frc.exe                ~~~  '
   WRITE(stdout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      
   !!
   !! get the namelist file name
   CALL get_command_argument(1, argument, ctrl, status)
   !
   SELECT CASE(status)
   CASE(0)
      namelistf = trim(argument)
   CASE(-1)
      WRITE(stdout,*) " ### Error: file name too long"
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
   READ(50,nml_dom, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1   
   READ(50,nml_opt, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   READ(50,nml_fld, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   READ(50,nml_out, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1  
   IF( ln_c1d ) then
      print*,'c1d is activated'
      READ(50,nml_c1d, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   END IF

   IF (ctrl > 0) then
      WRITE(stdout,*) " ### E R R O R while reading namelist file '",trim(namelistf),"'"
      WRITE(stdout,*) " ctrl = ",ctrl
      STOP
   ELSE 
      WRITE(stdout,*) " Namelist file ",trim(namelistf)," OK "    
   END IF
   IF( ln_hpg_frc .AND. ln_geo_wnd ) THEN
      WRITE(stdout,*) " ### E R R O R conflicting options "
      WRITE(stdout,*) " ln_hpg_frc and ln_geo_wnd can not both be set to True"
      STOP   
   END IF
   
   SELECT CASE (ptemp_method) 
   CASE(0)
      WRITE(stdout,*) " Absolute temperature option is activated"
   CASE(1)   
      WRITE(stdout,*) " Potential temperature option with  local reference pressure is activated"   
   CASE(2)    
      WRITE(stdout,*) " Potential temperature option with global reference pressure on temperature perturbation is activated"   
   CASE(3)    
      WRITE(stdout,*) " Potential temperature option with global reference pressure is activated"
   CASE(4)    
      WRITE(stdout,*) " Potential temperature option with global reference pressure is activated"  
   END SELECT
   
   IF ( ln_slp_smth .and. ( ( trim(slp_smth_method) .ne. 'laplacian' ) .and. ( trim(slp_smth_method) .ne. 'hanning' )  ) ) THEN
      WRITE(stdout,*) " ### E R R O R please define slp_smth_method = 'laplacian' or 'hanning' if you want to smooth slp field"
      STOP
   ELSE
      WRITE(stdout,*) " Smoothing of SLP activated with ", slp_smth_iter, "iteration(s) of ", trim(slp_smth_method), " filter"
   END IF
   
   IF(ln_hpg_frc) THEN
      WRITE(stdout,*) "### E R R O R in namelist variable "
      WRITE(stdout,*) "ln_hpg_frc not yet implemented in CROCO, please use ln_geo_wnd=.true. and ln_hpg_frc=.false."
      STOP
      WRITE(stdout,*) " Large-scale pressure gradient will be interpolated"
   ENDIF
   IF(ln_geo_wnd) THEN
      WRITE(stdout,*) " Geostrophic winds will be interpolated"
   ENDIF
   !!---------------------------------------------------------------------   
       
   !!-------------------------------------------------------------------------------------   
   !! list of variables to treat
   !!
   varnames = [character(len=6)   ::    'z',    't',    'q',    'u',    'v',  'msl',  'lsm',    'uhpg',   'vhpg' ]
   outnames = [character(len=6)   ::     '', 'tair', 'humi', 'uwnd', 'vwnd',     '',     '',    'uhpg',   'vhpg' ]
   filnames = [character(len=500) :: file_z, file_t, file_q, file_u, file_v, file_p, file_m, file_geos, file_hpg ]
 
   !!---------------------------------------------------------------------
   ! check files content
   ctrl = 0
   CALL Read_Ncdf_dim('lev'  ,trim(cn_dir)//'/'//trim(file_t),jpka_in  )
   !
   IF (ln_read_zsurf) THEN
     !IF( .not. VAR_EXISTENCE( trim(varnames(1)) , trim(cn_dir)//'/'//trim(file_z) )      &
     !   & .or. jpka_in /= 1 ) ctrl = ctrl + 1
     IF( .not. VAR_EXISTENCE( trim(varnames(1)) , trim(cn_dir)//'/'//trim(file_z) ) ) ctrl = ctrl + 1 
   END IF
   IF ( .not. VAR_EXISTENCE( trim(varnames(2)) , trim(cn_dir)//'/'//trim(file_t) ) ) ctrl = ctrl + 1 
   IF ( .not. VAR_EXISTENCE( trim(varnames(3)) , trim(cn_dir)//'/'//trim(file_q) ) ) ctrl = ctrl + 1   
   IF ( .not. VAR_EXISTENCE( trim(varnames(4)) , trim(cn_dir)//'/'//trim(file_u) ) ) ctrl = ctrl + 1
   IF ( .not. VAR_EXISTENCE( trim(varnames(5)) , trim(cn_dir)//'/'//trim(file_v) ) ) ctrl = ctrl + 1 
   IF ( .not. VAR_EXISTENCE( trim(varnames(6)) , trim(cn_dir)//'/'//trim(file_p) ) ) ctrl = ctrl + 1 
   IF (ln_read_mask) THEN
     IF ( .not. VAR_EXISTENCE( trim(varnames(7)) , trim(cn_dir)//'/'//trim(file_m) ) ) ctrl = ctrl + 1 
     varnames(7) = TRIM(mask_var)
   END IF
   IF(ln_hpg_frc) THEN
      IF ( .not. VAR_EXISTENCE( trim(varnames(8)) , trim(cn_dir)//'/'//trim(file_hpg) ) ) ctrl = ctrl + 1 
      IF ( .not. VAR_EXISTENCE( trim(varnames(9)) , trim(cn_dir)//'/'//trim(file_hpg) ) ) ctrl = ctrl + 1    
   END IF
   IF(ln_geo_wnd) THEN
      IF ( .not. VAR_EXISTENCE( trim(varnames(8)) , trim(cn_dir)//'/'//trim(file_geos) ) ) ctrl = ctrl + 1 
      IF ( .not. VAR_EXISTENCE( trim(varnames(9)) , trim(cn_dir)//'/'//trim(file_geos) ) ) ctrl = ctrl + 1    
   END IF
   !       
   IF ( ctrl > 0 ) THEN  
      WRITE(stdout,*) " ### E R R O R while reading ECMWF atmospheric files "
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
   WRITE(stdout,*) " jpka_in, jptime, jpi, jpj, nhym, nhyi: ", jpka_in, jptime, jpi, jpj, nhym, nhyi
   !
   !!---------------------------------------------------------------------


   !!---------------------------------------------------------------------
   !! allocate arrays  
   ALLOCATE( A_w    (               0:jpka_in) )         
   ALLOCATE( B_w    (               0:jpka_in) )
   ALLOCATE( e3t_in ( 1:jpi, 1:jpj, 1:jpka_in) )
   ALLOCATE( ghw_in ( 1:jpi, 1:jpj, 0:jpka_in) ) 
   ALLOCATE( slp    ( 1:jpi, 1:jpj              ) )
   ALLOCATE( zsurf  ( 1:jpi, 1:jpj              ) )
   ALLOCATE( temp   ( 1:jpi, 1:jpj, 1:jpka_in, 1) )
   ALLOCATE( humi   ( 1:jpi, 1:jpj, 1:jpka_in, 1) )       
   ALLOCATE( varout ( 1:jpi, 1:jpj, 1:jpka+1 , 1) )
   ALLOCATE( ind    ( 1:jpi, 1:jpj              ) )
   ALLOCATE( e3_bak ( 1:jpi, 1:jpj              ) )  
   ALLOCATE( ght    (               1:jpka+1    ) )
   ALLOCATE( ghw    (               1:jpka+1    ) )   
   ALLOCATE( e3t    (               1:jpka+1    ) )
   ALLOCATE( e3w    (               1:jpka+1    ) )      
   ALLOCATE( tmask  ( 1:jpi, 1:jpj              ) )
   IF( ln_c1d ) ALLOCATE( varc1d( 1:3, 1:3, 1:jpka+1 , 1 ) )
   IF( ln_c1d ) varc1d( 1:3, 1:3, 1:jpka+1 , 1 ) = 0.
   IF (jpka_in.NE.nhym) THEN
     ALLOCATE( tmp_fullw(1:nhyi) )
     ALLOCATE( tmp_fullm(1:nhym) )
   END IF
   !
   varout(:,:,:,1)    = 0.
      
   !!---------------------------------------------------------------------
   !! Read the mask and remove some closed seas 
   IF (ln_read_mask) THEN
     CALL init_atm_mask(jpi,jpj,trim(cn_dir)//'/'//trim(file_m),trim(mask_var),ln_lsm_land,tmask)
   ELSE
     tmask(:,:) = 1.
   END IF
   !!   
   
   !!---------------------------------------------------------------------
   !! Compute the altitude and layer thickness of the target grid
   CALL init_target_grid ( jpka, ght, ghw, e3t, e3w, hmax, hc, theta_s,   &
      &                                              ln_impose_z1, z1  )

   !! Write the grid file for the target grid
   CALL Write_Grid_File  ( jpka, ght, ghw, e3t, e3w, trim(cn_dir)//'/'//trim(grd_file) )

   !! Read the static A and B coefficients for the ECMWF vertical grid
   IF (jpka_in.EQ.nhym) THEN
     CALL Read_Ncdf_var    ( 'hyai', trim(cn_dir)//'/'//trim(file_t), A_w )
     CALL Read_Ncdf_var    ( 'hybi', trim(cn_dir)//'/'//trim(file_t), B_w )   
   ELSE
     CALL Read_Ncdf_var    ( 'hyai', trim(cn_dir)//'/'//trim(file_t), tmp_fullw )
     A_w(0:jpka_in) = tmp_fullw(nhyi-(jpka_in+1)+1:nhyi)
     CALL Read_Ncdf_var    ( 'hybi', trim(cn_dir)//'/'//trim(file_t), tmp_fullw )
     B_w(0:jpka_in) = tmp_fullw(nhyi-(jpka_in+1)+1:nhyi)
   END IF


   !!---------------------------------------------------------------------   
   !! create output file
   !!
   IF (Len_Trim(var_name) == 0) THEN
     out_file = trim(cn_dir)//'/'//trim(abl_file)
   ELSE
     out_file = trim(cn_dir)//'/'//trim(var_name)//'_'//trim(abl_file)
   END IF
   IF(ln_c1d) THEN
      CALL Init_output_File_c1d ( jpi, jpj, jpka, trim(cn_dir)//'/'//trim(file_t), out_file, tmask(:,:), iloc, jloc )
   ELSE 
      CALL Init_output_File ( jpi, jpj, jpka, trim(cn_dir)//'/'//trim(file_t), out_file, tmask(:,:) )
   END IF

   !!---------------------------------------------------------------------   
   !! Initialize the name of the dimensions for the result of the interpolation
   !!
   dimnames(1) = 'lon'
   dimnames(2) = 'lat'       
   dimnames(3) = 'jpka'
   dimnames(4) = 'time'
   !CALL Write_Ncdf_var( 'tmask', dimnames(1:2), trim(out_file), tmask, 'float' )

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
         CALL add_globatt_real( out_file, "jpka"   , REAL(jpka) )
         CALL add_globatt_real( out_file, "hmax"   , hmax       )
         CALL add_globatt_real( out_file, "theta_s", theta_s    )
         CALL add_globatt_real( out_file, "hc"     , hc         )
         IF (ln_impose_z1) CALL add_globatt_real( out_file, "z1", z1 )
      ENDIF
      !
      ! read SLP
      CALL Read_Ncdf_var ( varnames(6) , trim(cn_dir)//'/'//trim(file_p), slp,  kt )      !<-- (log of) surface pressure
      IF (ln_slp_log) THEN
        DO jj = 1, jpj
           DO ji = 1, jpi      
             slp( ji ,jj )   = exp( slp(ji ,jj ) )
           END DO
        END DO
      ENDIF
      !
      ! read ZSURF
      IF (ln_read_zsurf) THEN
        CALL Read_Ncdf_var ( varnames(1) , trim(cn_dir)//'/'//trim(file_z), zsurf,  kt )      !<-- surface geopotential 
      ELSE
        zsurf(:,:) = 0.
      END IF
      !    
      ! Smoothing of SLP and ZSURF to remove gibbs oscillations (must be done on both fields or none of them)
      IF( ln_slp_smth ) CALL smooth_field( jpi, jpj, slp(:,:), tmask(:,:), slp_smth_method, slp_smth_iter )
      IF (ln_read_zsurf.AND.ln_slp_smth) CALL smooth_field( jpi, jpj, zsurf(:,:), tmask(:,:), slp_smth_method, slp_smth_iter )
      !IF( ln_slp_smth ) CALL DTV_Filter( jpi, jpj, slp(:,:), tmask(:,:), 25, kt )  !<-- not yet robust enough
      !
      ! read TEMP and HUMI
      CALL Read_Ncdf_var ( varnames(2) , trim(cn_dir)//'/'//trim(file_t), temp, kt    )      !<-- temperature   
      CALL Read_Ncdf_var ( varnames(3) , trim(cn_dir)//'/'//trim(file_q), humi, kt    )      !<-- humidity        
      !
      ! Reconstruct the ERA-Interim vertical grid in terms of altitude 
      ghw_in(:,:,1:jpka_in)  = 0.
      ghw_in(:,:,  jpka_in)  = zsurf(:,:) * (1. / grav)
      CALL get_atm_grid( jpi, jpj, jpka_in, slp, temp,  &
         &                          humi, A_w, B_w, e3t_in, ghw_in )
      !
      ! Compute potential temperature
      CALL get_pot_temp    ( jpi, jpj, jpka_in, slp, temp, A_w, B_w,   &
         &                   tmask(:,:), ptemp_method, humi(:,:,jpka_in,1), 0.5*ghw_in(:,:,jpka_in-1) )

      ! Flip the vertical axis to go from k=0 at the bottom to k=N_in at the top of the atmosphere
      CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, e3t_in          )       
      CALL flip_vert_dim ( 0, jpka_in, jpi, jpj, ghw_in          ) 
      CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, temp(:,:,:,1) )       
      CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, humi(:,:,:,1) ) 
      
      ! Correct the layer thickness to match hmax  
      DO jj = 1, jpj
         DO ji = 1, jpi          
            cff = 0.
            DO jk=1,jpka_in
               cff = cff + e3t_in( ji, jj, jk )
               IF ( cff > hmax ) THEN
                  jk_in = jk
                  EXIT
               END IF 
            END DO  
            ind    ( ji, jj       ) = jk_in    
            e3_bak ( ji, jj       ) = e3t_in ( ji, jj, jk_in )  ! store the value of the original layer thickness 
            e3t_in ( ji, jj, jk_in) = e3t_in ( ji, jj, jk_in ) - ( cff - hmax )
         END DO  
      END DO 
      !
      IF (Len_Trim(var_name) == 0) THEN
        !
        ! Interpolation of potential TEMP
        CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
            &                         temp, e3t_in, e3_bak, e3t, varout, jp_weno )  
        varout(:,:,1,1) = varout(:,:,2,1)
        IF(ln_c1d) THEN
           DO jj=1,3
              DO ji=1,3
                 varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
              END DO
           END DO
           CALL Write_Ncdf_var ('tair', dimnames(1:4), out_file, varc1d, kt, 'float' )
        ELSE
           CALL Write_Ncdf_var ('tair', dimnames(1:4), out_file, varout, kt, 'float' ) 
        ENDIF
        !
        ! Interpolation of Q      
        CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
            &                        humi, e3t_in, e3_bak, e3t, varout, jp_weno  )
        ! FL: dirty loop (possible issue in boundary conditions for WENO scheme)
        DO jj = 1, jpj
           DO ji = 1, jpi
              DO jk = 2, jpka+1
                 varout(ji,jj,jk,1) = MAX(varout(ji,jj,jk,1),1.E-08)                       !<-- negative values in ECMWF
              END DO
           END DO
        END DO

        varout(:,:,1,1) = varout(:,:,2,1) 

        IF(ln_c1d) THEN
           DO jj=1,3
              DO ji=1,3
                 varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
              END DO
           END DO
           CALL Write_Ncdf_var  ('humi', dimnames(1:4), out_file, varc1d, kt, 'float' )
        ELSE
           CALL Write_Ncdf_var ('humi', dimnames(1:4), out_file, varout, kt, 'float' )  
        ENDIF

        !
        ! Interpolate large-scale HPG or geostrophic wind
        !
        IF ( ln_hpg_frc .OR. ln_geo_wnd ) THEN   
           !
           IF( ln_hpg_frc ) THEN 
              file_in = trim(file_hpg )
           ELSE
              file_in = trim(file_geos)
           ENDIF

           ! Read geostrophic wind
           CALL Read_Ncdf_var ( varnames(8) , trim(cn_dir)//'/'//trim(file_in), temp, kt )        
           CALL Read_Ncdf_var ( varnames(9) , trim(cn_dir)//'/'//trim(file_in), humi, kt ) 
           CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, temp( :,:,:,1 ) )       
           CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, humi( :,:,:,1 ) ) 

           ! Interpolation of geostrophic U
           CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
               &                         temp, e3t_in, e3_bak, e3t, varout, jp_spln  )      
           varout(:,:,1,1) = varout(:,:,2,1)
           
           IF(ln_c1d) THEN
              DO jj=1,3
                 DO ji=1,3
                    varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
                 END DO
              END DO
              CALL Write_Ncdf_var ( 'uhpg', dimnames(1:4), out_file, varc1d, kt, 'float' )
           ELSE
              CALL Write_Ncdf_var ( 'uhpg', dimnames(1:4), out_file, varout, kt, 'float' )
           ENDIF
           !
           ! Interpolation of geostrophic V
           CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
               &                         humi, e3t_in, e3_bak, e3t, varout, jp_spln  )       
           varout(:,:,1,1) = varout(:,:,2,1)

           IF(ln_c1d) THEN
              DO jj=1,3
                 DO ji=1,3
                    varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
                 END DO
              END DO
              CALL Write_Ncdf_var ( 'vhpg', dimnames(1:4), out_file, varc1d, kt, 'float' )
           ELSE
              CALL Write_Ncdf_var ( 'vhpg', dimnames(1:4), out_file, varout, kt, 'float' )
           ENDIF
        
        END IF       
        !
        ! Interpolation of total winds
        !      
        ! Read wind
        CALL Read_Ncdf_var ( varnames(4) , trim(cn_dir)//'/'//trim(file_u), temp, kt )        
        CALL Read_Ncdf_var ( varnames(5) , trim(cn_dir)//'/'//trim(file_v), humi, kt ) 
        CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, temp( :,:,:,1 ) )       
        CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, humi( :,:,:,1 ) ) 
        !      
        ! Interpolation of total U
        CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
            &                         temp, e3t_in, e3_bak, e3t, varout, jp_spln  )      
        varout(:,:,1,1) = varout(:,:,2,1)

        IF(ln_c1d) THEN
           DO jj=1,3
              DO ji=1,3
                 varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
              END DO
           END DO
           CALL Write_Ncdf_var ( 'uwnd', dimnames(1:4), out_file, varc1d, kt, 'float' )
        ELSE
           CALL Write_Ncdf_var ( 'uwnd', dimnames(1:4), out_file, varout, kt, 'float' )
        ENDIF
        !
        ! Interpolation of total V
        CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
            &                         humi, e3t_in, e3_bak, e3t, varout, jp_spln  )       
        varout(:,:,1,1) = varout(:,:,2,1)

        IF(ln_c1d) THEN
           DO jj=1,3
              DO ji=1,3
                 varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
              END DO
           END DO
           CALL Write_Ncdf_var ( 'vwnd', dimnames(1:4), out_file, varc1d, kt, 'float' )
        ELSE
           CALL Write_Ncdf_var ( 'vwnd', dimnames(1:4), out_file, varout, kt, 'float' )
        ENDIF

      ELSE
        !
        ctrl     =  minval(pack([(ji,ji=1,size(outnames))],outnames==var_name))
        var_file = filnames(ctrl)
        !
        !
        ! Interpolation of var_name
        !      
        ! Read var_name
        CALL Read_Ncdf_var ( varnames(ctrl) , trim(cn_dir)//'/'//trim(var_file), temp, kt )
        CALL flip_vert_dim ( 1, jpka_in, jpi, jpj, temp( :,:,:,1 ) )       
        !      
        ! Interpolation of var_name
        CALL zinterp ( jpi, jpj, jpka, jpka_in, ind,  &
            &                         temp, e3t_in, e3_bak, e3t, varout, jp_spln  )      
        varout(:,:,1,1) = varout(:,:,2,1)

        IF(ln_c1d) THEN
           DO jj=1,3
              DO ji=1,3
                 varc1d( ji, jj , 1:jpka+1 , 1 ) = varout(iloc,jloc, 1:jpka+1 , 1 )
              END DO
           END DO
           CALL Write_Ncdf_var ( outnames(ctrl), dimnames(1:4), out_file, varc1d, kt, 'float' )
        ELSE
           CALL Write_Ncdf_var ( outnames(ctrl), dimnames(1:4), out_file, varout, kt, 'float' )
        ENDIF
        !
      END IF
   !
   END DO ! kt
   !
   WRITE(*,*) 'Altiude of IFS grid : ghw_in(10,10,:) = ', ghw_in(10,10,:)
   WRITE(*,*) 'Altiude of IFS grid : e3t_in(10,10,:) = ', e3t_in(10,10,:)
   !
   DEALLOCATE(zsurf,slp,temp,humi,varout)     
   IF (jpka_in.NE.nhym) DEALLOCATE(tmp_fullw,tmp_fullm)
   !
   STOP
   !
END PROGRAM main







Function to_upper (str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(3), Intent(In) :: str
    Character(3)             :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

End Function to_upper
