!---------------------------------------------------------------------------------
! This code is NetCDF IO routine
!        
! Author: Jung H. C. (dhtyphoon@naver.com)
! History
!       - 2018-03-08: initial version
!    
!---------------------------------------------------------------------------------



      MODULE NETCDF_IO

      USE netcdf

      IMPLICIT NONE

      CONTAINS

!!==============================================================================

      SUBROUTINE NF_CREATE

      USE netcdf_var

      IMPLICIT NONE

      !!NF90_CREATE
      ncstat=NF90_CREATE(TRIM(nc_fname), NF90_64BIT_OFFSET, nc_file_id)
      IF ( ncstat /= NF90_NOERR ) THEN
         print*,'NF90_CREATE Error'
      ENDIF

      END SUBROUTINE NF_CREATE

!!==============================================================================
      SUBROUTINE NF_OPEN

      USE netcdf_var

      IMPLICIT NONE

      !!NF90_CREATE
      ncstat=NF90_OPEN(TRIM(nc_fname), NF90_64BIT_OFFSET, nc_file_id)
      IF ( ncstat /= NF90_NOERR ) THEN
          print*,'NF90_OPEN Error'
      ENDIF

      END SUBROUTINE NF_OPEN

!!==============================================================================

      SUBROUTINE DEF_DIM

      USE netcdf_var

      IMPLICIT NONE

      !!NF90_DEF_DIM
      ncstat=NF90_DEF_DIM(nc_file_id, 'xaxis_1', nx, nx_id)
      IF ( ncstat /= NF90_NOERR ) THEN
         print*,'NF90_DEF_DIM xaxis_1 Error'
      ENDIF
      ncstat=NF90_DEF_DIM(nc_file_id, 'yaxis_1', ny, ny_id)
      IF ( ncstat /= NF90_NOERR ) THEN
         print*,'NF90_DEF_DIM yaxis_1 Error'
      ENDIF
      ncstat=NF90_DEF_DIM(nc_file_id, 'depth', nz, nz_id)
      IF ( ncstat /= NF90_NOERR ) THEN
         print*,'NF90_DEF_DIM zaxis_1 Error'
      ENDIF
      ncstat=NF90_DEF_DIM(nc_file_id, 'time', NF90_UNLIMITED, nt_id)
      IF ( ncstat /= NF90_NOERR ) THEN
         print*,'NF90_DEF_DIM time Error'
      ENDIF

      nc_dims_id(1)=nx_id
      nc_dims_id(2)=ny_id
      nc_dims_id(3)=nz_id
      nc_dims_id(4)=nt_id

      END SUBROUTINE DEF_DIM

!!=============================================================================


      SUBROUTINE DEF_VAR(nc_var_id)

      USE netcdf_var

      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nc_var_id

      !!NF90_DEF_VAR
      ncstat=NF90_DEF_VAR(nc_file_id, TRIM(varname_in), NF90_double,&
                          nc_dims_id, nc_var_id)
      IF (ncstat/=NF90_NOERR) THEN
          print*,'NF90_DEF_VAR ', TRIM(varname_in), ' Error'
      ENDIF


      END SUBROUTINE DEF_VAR


!!=============================================================================

      SUBROUTINE DEF_TIME

      USE netcdf_var

      IMPLICIT NONE

       !!NF90_DEF_VAR
       ncstat=NF90_DEF_VAR(nc_file_id, 'time', NF90_double,&
                           nt_id, time_id)
       IF (ncstat/=NF90_NOERR) THEN
          print*,'DEF_TIME 1 Error'
       ENDIF

       ncstat=NF90_PUT_ATT(nc_file_id, time_id, 'units', TRIM(times))
       IF (ncstat/=NF90_NOERR) THEN
           print*,'DEF_TIME 2 ERROR'
       ENDIF


      END SUBROUTINE DEF_TIME

!!===============================================================================
     
      SUBROUTINE DEF_DEPTH

      USE netcdf_var

      IMPLICIT NONE

      !!NF90_DEF_VAR
      ncstat=NF90_DEF_VAR(nc_file_id, 'depth', NF90_double,&
                          nz_id, depth_id)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'DEF_Z1 Error'
      ENDIF

      ncstat=NF90_PUT_ATT(nc_file_id, depth_id, 'units', 'm')
      IF (ncstat/=NF90_NOERR) THEN
          print*,'DEF_Z2 ERROR'
      ENDIF

      END SUBROUTINE DEF_DEPTH

!!===============================================================================

      SUBROUTINE PUT_ATT(nc_var_id)

      USE netcdf_var

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nc_var_id

      ncstat=NF90_PUT_ATT(nc_file_id, nc_var_id, 'long_name', TRIM(long_name))
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_ATT long_name ', long_name,'  ERROR'
      ENDIF

      ncstat=NF90_PUT_ATT(nc_file_id, nc_var_id, 'title', TRIM(short_name))
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_ATT title ', short_name,'  ERROR'
      ENDIF

      ncstat=NF90_PUT_ATT(nc_file_id, nc_var_id, 'short_name', TRIM(short_name))
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_ATT short_name ', short_name,'  ERROR'
      ENDIF

      ncstat=NF90_PUT_ATT(nc_file_id, nc_var_id, 'units', TRIM(units))
      IF (ncstat/=NF90_NOERR) THEN
        print*,'NF90_PUT_ATT units ', units,'  ERROR'
      ENDIF

      END SUBROUTINE PUT_ATT

!!================================================================================

      SUBROUTINE END_DEF

      USE netcdf_var

      IMPLICIT NONE

      ncstat=NF90_ENDDEF(nc_file_id)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_ENDDEF Error'
      ENDIF

      END SUBROUTINE END_DEF

!!================================================================================

      SUBROUTINE PUT_VAR(nc_var_id, var_data, nlev, itt)

      USE netcdf_var

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nc_var_id
      INTEGER, INTENT(IN) :: nlev
      INTEGER, INTENT(IN) :: itt
      REAL, DIMENSION(:,:,:), INTENT(IN) :: var_data
      INTEGER, DIMENSION(4) :: start, edges

      REAL, ALLOCATABLE, DIMENSION(:) :: dum

      if (allocated(dum) ) deallocate(dum);allocate(dum(nlev));dum=0.0

      dum(1:nlev)=var_data(1,1,1:nlev)
      start(1) = 1; edges(1) = 1
      start(2) = 1; edges(2) = 1
      start(3) = 1; edges(3) = nlev
      start(4) = itt; edges(4) = 1

      ncstat=NF90_PUT_VAR(nc_file_id, nc_var_id, dum, start, edges)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_VAR Error'
      ENDIF

      END SUBROUTINE PUT_VAR

!!================================================================================

      SUBROUTINE PUT_TIME(itt, var_data)

      USE netcdf_var

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: itt
      REAL, INTENT(IN) :: var_data

      INTEGER, DIMENSION(4) :: start, edge
      REAL, DIMENSION(4) :: dump

      start(1)=itt
      edge(1)=1
      dump(1)=var_data

      ncstat=NF90_PUT_VAR(nc_file_id, time_id, dump, start, edge)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_TIME Error'
      ENDIF

      END SUBROUTINE PUT_TIME

!!================================================================================

      SUBROUTINE PUT_DEPTH(var_data)

      USE netcdf_var

      IMPLICIT NONE

      REAL, DIMENSION(:), INTENT(IN) :: var_data

!      INTEGER, DIMENSION(4) :: start, edge
!      REAL, DIMENSION(4) :: dump

!      start(1)=1
!      edge(1)=
!      dump(1)=var_data

      ncstat=NF90_PUT_VAR(nc_file_id, depth_id, var_data)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_PUT_DEPTH Error'
      ENDIF

      END SUBROUTINE PUT_DEPTH
!!================================================================================

      SUBROUTINE GET_VAR4D(varname, var_data)
     
      USE netcdf_var

      IMPLICIT NONE

      INTEGER             :: priv_var_id
      CHARACTER(LEN=32), INTENT(IN)         :: varname
      REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: var_data
   
      ncstat=NF90_INQ_VARID(nc_file_id, TRIM(varname), priv_var_id)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_INQ_VARID Error',varname
      ENDIF

      varname_in=varname

      ncstat=NF90_GET_VAR(nc_file_id, priv_var_id, var_data)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_GET_VAR Error'
      ENDIF

      END SUBROUTINE GET_VAR4D


!!================================================================================

      SUBROUTINE GET_VAR3D(varname, var_data)

      USE netcdf_var

      IMPLICIT NONE

      INTEGER             :: priv_var_id
      CHARACTER(LEN=32), INTENT(IN)         :: varname
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: var_data

      ncstat=NF90_INQ_VARID(nc_file_id, TRIM(varname), priv_var_id)
      IF (ncstat/=NF90_NOERR) THEN
        print*,'NF90_INQ_VARID Error',varname
      ENDIF

      varname_in=varname

      ncstat=NF90_GET_VAR(nc_file_id, priv_var_id, var_data)
      IF (ncstat/=NF90_NOERR) THEN
        print*,'NF90_GET_VAR Error'
      ENDIF

      END SUBROUTINE GET_VAR3D

!!================================================================================      
      SUBROUTINE NF_CLOSE

      USE netcdf_var

      IMPLICIT NONE

      ncstat=NF90_CLOSE(nc_file_id)
      IF (ncstat/=NF90_NOERR) THEN
         print*,'NF90_CLOSE Error'
      ENDIF
      
      END SUBROUTINE NF_CLOSE


      END MODULE NETCDF_IO

