!---------------------------------------------------------------------------------
! This code is NetCDF IO routine
!        
! Author: Jung H. C. (dhtyphoon@naver.com)
! History
!       - 2018-03-08: initial version
!    
!---------------------------------------------------------------------------------

      MODULE netcdf_var

      IMPLICIT NONE

      CHARACTER(LEN=100) :: nc_fname
      INTEGER :: nc_file_id

      INTEGER :: nx, ny, nz
      INTEGER :: nx_id, ny_id, nz_id, nt_id
      INTEGER :: time_id, depth_id
      INTEGER, DIMENSION(4) :: nc_dims_id

      CHARACTER(LEN=32) :: varname_in
      CHARACTER(LEN=100) :: long_name, short_name, units, times

      INTEGER :: ncstat

      END MODULE netcdf_var
