PROGRAM main
   !!======================================================================
   !!                     ***  PROGRAM main  ***
   !!
   !! ** Purpose : horizontal extrapolation ("drowning") of ECMWF dataset
   !!
   !! 
   !!======================================================================
   !! History : 2016-10  (F. Lemarié)  Original code largely inspired by SOSIE (L. Brodeau & J.M. Molines)
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
   INTEGER                                  :: ji,jj,jk,kt,kv,jjp1,jjm1,jip1,jim1
   INTEGER                                  :: jpka,jpvar            ! number of vertical levels for input and target grids 
   INTEGER                                  :: jpi , jpj             ! number of grid points in x and y directions                 
   INTEGER                                  :: jptime,ctrl,niter,status,ncid
   INTEGER                                  :: ioerr
   INTEGER, PARAMETER                       :: stdout   = 6
   INTEGER, PARAMETER                       :: max_iter = 50
   !!
   REAL(8)                                  :: cff,cff1,cff2
   !!  
   REAL(8), ALLOCATABLE, DIMENSION(:,:,:,:) :: varin
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: tmask
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: mask_coast,maskv              ! land-sea mask   
   REAL(8), ALLOCATABLE, DIMENSION(:,:    ) :: varold,varnew   
   REAL(8), ALLOCATABLE, DIMENSION(:      ) :: tmp1d
   !!
   CHARACTER(len=500)                       :: file_u,file_v,file_hpg,file_geos ! ERAi files containing wind components 
   CHARACTER(len=500)                       :: file_t,file_q, file_m   ! ERAi files containing tracers and mask 
   CHARACTER(len=500)                       :: file_z,file_p,cn_dir    ! ERAi files containing surface geopot and pressure
   CHARACTER(len=500)                       :: abl_file, grd_file, drwn_file, out_file
   CHARACTER(len=500)                       :: slp_smth_method ! Method used to smooth the slp : choice between 'hanning' or 'laplacian'
   CHARACTER(len=500)                       :: namelistf              
   CHARACTER(len=500)                       :: argument   
   CHARACTER(len= 20), DIMENSION(4)         :: dimnames
   CHARACTER(:), ALLOCATABLE, DIMENSION(:)  :: varnames
   CHARACTER(6)                             :: mask_var            ! name of mask variable in file_m file  
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
   INTEGER                                  :: ptemp_method       ! way to compute potential temperature 
   INTEGER                                  :: slp_smth_iter      ! Number of iterations of the smoothing filter       


   !!---------------------------------------------------------------------
   !! List of variables read in the namelist file 
   NAMELIST/nml_out/    grd_file, abl_file, drwn_file, var_name
   NAMELIST/nml_opt/    ptemp_method, ln_slp_log, ln_slp_smth, slp_smth_method, slp_smth_iter, ln_read_mask, ln_perio_latbc, &
      &                 ln_hpg_frc, ln_geo_wnd, ln_c1d, ln_read_zsurf, ln_lsm_land
   NAMELIST/nml_fld/    cn_dir, file_u, file_v, file_t,              &
      &                 file_q, file_z, file_p, file_hpg, file_geos, &
      &                 file_m, mask_var

   WRITE(stdout,*) '                                                           '
   WRITE(stdout,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(stdout,*) '  ~~~               drown_atm_frc.exe                 ~~~  '
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
   READ(50,nml_opt, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   READ(50,nml_fld, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   READ(50,nml_out, iostat=ioerr); IF (ioerr /= 0) ctrl = ctrl + 1
   
   IF (ctrl > 0) then
      WRITE(stdout,*) "### E R R O R while reading namelist file '",trim(namelistf),"'"
      WRITE(stdout,*) " ctrl = ",ctrl
      STOP
   ELSE 
      WRITE(stdout,*) " Namelist file ",trim(namelistf)," OK "    
   END IF
   
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
   ENDIF

   !!-------------------------------------------------------------------------------------   
   !! list of variables to treat
   !!
   IF (Len_Trim(var_name) == 0) THEN
     jpvar    = 6
     allocate(character(len=20) :: varnames(jpvar) )
     varnames = [character(len=8) :: 'uwnd', 'vwnd', 'tair', 'humi', 'uhpg', 'vhpg' ]
     IF( .NOT. Var_Existence( 'uhpg' , trim(cn_dir)//'/'//trim(abl_file) ) ) jpvar = jpvar - 1
     IF( .NOT. Var_Existence( 'vhpg' , trim(cn_dir)//'/'//trim(abl_file) ) ) jpvar = jpvar - 1 
   ELSE
     jpvar     = 1 
     allocate(character(len=20) :: varnames(jpvar) )
     varnames  = var_name
     drwn_file = trim(var_name)//'_'//trim(drwn_file)
   END IF 
   WRITE(stdout,*) " Number of variables to treat : ",jpvar    
   !!-------------------------------------------------------------------------------------    
   !! read the dimensions for the input files
   CALL Read_Ncdf_dim ( 'jpka' , trim(cn_dir)//'/'//trim(abl_file), jpka    ) 
   CALL Read_Ncdf_dim ( 'time' , trim(cn_dir)//'/'//trim(abl_file), jptime  ) 
   CALL Read_Ncdf_dim ( 'lon'  , trim(cn_dir)//'/'//trim(abl_file), jpi     )     
   CALL Read_Ncdf_dim ( 'lat'  , trim(cn_dir)//'/'//trim(abl_file), jpj     )        
   WRITE(stdout,*) " jpka, jptime, jpi, jpj: ", jpka, jptime, jpi, jpj
   !
   !!---------------------------------------------------------------------


   !!---------------------------------------------------------------------
   !! allocate arrays  
   ALLOCATE(varin (1:jpi,1:jpj,1:jpka,1))
   ALLOCATE(tmask (1:jpi,1:jpj    ))
   ALLOCATE(maskv      (1:jpi,1:jpj )) 
   ALLOCATE(mask_coast (1:jpi,1:jpj ))      
   ALLOCATE(varold     (1:jpi,1:jpj ))  
   ALLOCATE(varnew     (1:jpi,1:jpj ))     

   !!---------------------------------------------------------------------
   !! Read the mask and remove some closed seas 
   IF (ln_read_mask) THEN
     CALL init_atm_mask(jpi,jpj,trim(cn_dir)//'/'//trim(file_m),trim(mask_var),ln_lsm_land,tmask)
   ELSE
     tmask(:,:) = 1.
   END IF
   !!   



   !!---------------------------------------------------------------------   
   !! create output file
   !!
   out_file = trim(cn_dir)//'/'//trim(drwn_file)
   CALL Init_output_File ( jpi, jpj, jpka-1, trim(cn_dir)//'/'//trim(abl_file), out_file, tmask(:,:) )

   !!---------------------------------------------------------------------   
   !! Initialize the name of the dimensions for the result of the drowning
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

   !===========
   DO kv = 1,jpvar
   !===========
      DO kt = 1,jptime
      !===========

         IF( kv == 1 ) THEN
           CALL Write_Ncdf_var( 'time', dimnames(4:4), trim(out_file), tmp1d(kt:kt), kt, 'double' )
         ENDIF
         WRITE(stdout,*) ' Treat variable ',trim(varnames(kv)),' at time =', kt, '/', jptime
         CALL Read_Ncdf_var ( trim(varnames(kv)) , trim(cn_dir)//'/'//trim(abl_file), varin,  kt )         

         !===========
         DO jk = 1,jpka
         !===========

            ! Read variable to treat
            maskv  (1:jpi,1:jpj) = tmask(1:jpi,1:jpj )
            varold (1:jpi,1:jpj) = varin(1:jpi,1:jpj,jk,1)
            varnew (1:jpi,1:jpj) = varin(1:jpi,1:jpj,jk,1)            
            mask_coast(1:jpi,1:jpj) = 0.

            !$$$$$$$$$$$$$$$$
            DO niter = 1, max_iter
            !$$$$$$$$$$$$$$$$                 

               varold(1:jpi,1:jpj) = varnew (1:jpi,1:jpj)

               IF ( .NOT. (ANY(maskv == 0))  ) THEN
                  EXIT
               END IF

               !+++++++++++++++++++++++++++++++++++
               ! Build mask_coast
               DO jj = 1,jpj
                  DO ji = 1,jpi

                     IF (ln_perio_latbc) THEN
                       jip1 = ji + 1; if(ji==jpi) jip1 = 1
                       jim1 = ji - 1; if(ji==  1) jim1 = jpi
                     ELSE
                       jip1 = ji + 1; if(ji==jpi) jip1 = jpi-1
                       jim1 = ji - 1; if(ji==  1) jim1 = 2
                     END IF
                     jjp1 = jj + 1; if(jj==jpj) jjp1 = jpj-1
                     jjm1 = jj - 1; if(jj==  1) jjm1 = 2

                     cff  = (1.-maskv(ji,jj))
                     cff1 = maskv(ji,jjp1)+maskv(jip1,jj)+maskv(jim1,jj)+maskv(ji,jjm1) 
                     mask_coast(ji,jj) = cff * cff1 
                     IF( mask_coast(ji,jj) > 0. ) mask_coast(ji,jj) = 1.

                  END DO
               END DO   
               !+++++++++++++++++++++++++++++++++++          

               DO jj=1,jpj
                  DO ji=1,jpi

                     IF ( mask_coast(ji,jj) == 1. ) THEN

                        IF (ln_perio_latbc) THEN
                          jip1 = ji + 1; if(ji==jpi) jip1 = 1
                          jim1 = ji - 1; if(ji==  1) jim1 = jpi
                        ELSE
                          jip1 = ji + 1; if(ji==jpi) jip1 = jpi-1
                          jim1 = ji - 1; if(ji==  1) jim1 = 2
                        END IF
                        jjp1 = jj + 1; if(jj==jpj) jjp1 = jpj-1
                        jjm1 = jj - 1; if(jj==  1) jjm1 = 2

                        cff  = maskv(jim1,jjm1)+maskv(jim1,jj)+maskv(jim1,jjp1) &
                           & + maskv(ji  ,jjm1)+               maskv(ji  ,jjp1) &
                           & + maskv(jip1,jjm1)+maskv(jip1,jj)+maskv(jip1,jjp1)

                        varnew(ji,jj) = (1./cff)*(                                   &
                           &                    varold(jim1,jjm1)*maskv(jim1,jjm1)   &
                           &              +     varold(jim1,jj  )*maskv(jim1,jj  )   &
                           &              +     varold(jim1,jjp1)*maskv(jim1,jjp1)   &
                           &              +     varold(ji  ,jjm1)*maskv(ji  ,jjm1)   &
                           &              +     varold(ji  ,jjp1)*maskv(ji  ,jjp1)   &                
                           &              +     varold(jip1,jjm1)*maskv(jip1,jjm1)   &
                           &              +     varold(jip1,jj  )*maskv(jip1,jj  )   &
                           &              +     varold(jip1,jjp1)*maskv(jip1,jjp1)   )

                     END IF

                  END DO
               END DO   
               !+++++++++++++++++++++++++++++++++++
               maskv(1:jpi,1:jpj) = maskv(1:jpi,1:jpj) + mask_coast(1:jpi,1:jpj)
               !+++++++++++++++++++++++++++++++++++                   
           !$$$$$$$$$$$$$$$$
           END DO
           !$$$$$$$$$$$$$$$$

           ! DEBUG
           !IF (( kv == 1 ).AND.( kt == 1 ).AND.( jk == 1 )) THEN
           !   CALL Write_Ncdf_var( 'maskv', dimnames(1:2), trim(out_file), maskv, 'float' )
           !END IF

           ! SMOOTHING AFTER DROWNING -> only over continent: bug ?
           !maskv (1:jpi,1:jpj     ) = 1. - tmask(1:jpi,1:jpj  )
           !CALL smooth_field( jpi, jpj, varnew, maskv, slp_smth_method, slp_smth_iter )    
           !CALL smooth_field( jpi, jpj, varnew, maskv, slp_smth_method, slp_smth_iter )

           varin(1:jpi,1:jpj,jk,1) = varnew(1:jpi,1:jpj)

        !===========
        END DO ! jk
        !===========

        ! MASK REMAINING BAD (MISSING/UNITIALIZED) VALUES (IN ANTARCTIC ONLY)
        !WHERE (varin < -1.e+10) varin = -9.e+33

        CALL Write_Ncdf_var ( trim(varnames(kv)), dimnames(1:4), trim(cn_dir)//'/'//drwn_file, varin, kt, 'float' )

      !===========
      END DO ! time
      !===========

   !===========
   END DO ! variable
   !===========
   
   STOP
   !
END PROGRAM main
