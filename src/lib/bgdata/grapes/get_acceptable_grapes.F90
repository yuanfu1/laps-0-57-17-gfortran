SUBROUTINE get_acceptable_grapes(bgpath,i4time_needed,nfiles, &
           i4times_found,filenames)

  USE time_utils

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: bgpath
  INTEGER,INTENT(IN)            :: i4time_needed
  INTEGER,INTENT(OUT)           :: nfiles
  INTEGER,INTENT(OUT)           :: i4times_found(2)
  CHARACTER(LEN=256)            :: filenames(4) ! plus 1 dim. by Wei-Ting (130312) to put previous time
  CHARACTER(LEN=256)            :: thisfile

  ! Local vars
  INTEGER  :: delta_t, delta_t_1, delta_t_2, istatus
  INTEGER, PARAMETER  :: max_files = 500
  INTEGER, PARAMETER  :: dtmax = 86400
  INTEGER, PARAMETER  :: dtmin = -86400
  CHARACTER(LEN=11),PARAMETER   :: filter = "grapesinput"
  CHARACTER(LEN=24)             :: grptime
  INTEGER                       :: i4time
  INTEGER                       :: fname_len

  ! Local variables
  INTEGER                       :: total_files, grapes_files, i
  CHARACTER(LEN=256)            :: all_files(max_files)

  delta_t_1 = dtmax
  delta_t_2 = dtmin  
  CALL get_file_names(bgpath,total_files,all_files,max_files,istatus)
  
  grapes_files = 0
  DO i=1,total_files
    IF (all_files(i)(LEN_TRIM(bgpath)+2:(LEN_TRIM(bgpath)+12)) .EQ. filter) THEN
      grapes_files = grapes_files + 1
      WRITE(*,*) '  Found GRAPES file: ', TRIM(all_files(i))
      all_files(grapes_files) = all_files(i)
    ELSE IF (all_files(i)(LEN_TRIM(bgpath)+2:(LEN_TRIM(bgpath)+9)) .EQ. "namelist") THEN
      WRITE(*,*) '  Found Namelist file: ', TRIM(all_files(i))
      all_files(max_files) = all_files(i) ! store namelist file in the last slot of all_files for later use
    END IF
    WRITE(*,*) '  Found file: ', TRIM(all_files(i)(LEN_TRIM(bgpath)+2:LEN_TRIM(all_files(i))))
    WRITE(*,*) '  Path length: ', LEN_TRIM(bgpath), LEN_TRIM(all_files(i))
  END DO

  DO i=1,grapes_files
    WRITE(*,*) 'GRAPES file: ', TRIM(all_files(i))

    thisfile = all_files(i)
    fname_len = LEN_TRIM(all_files(i))

    ! WRF file format is assumed to be:    wrfout_d01_YYYY-MM-DD_HH:MM:SS
    ! GRAPES file format is assumed to be: grapesinputYYYYMMDDHH
    PRINT*, '  Checking file: ', thisfile(fname_len-20:fname_len-10), ' ', &
      thisfile(fname_len-9:fname_len-6), ' ', thisfile(fname_len-5:fname_len)

    IF (thisfile(fname_len-20:fname_len-10) .EQ. filter) THEN
      ! wrftime = thisfile(fname_len-18:fname_len) // ".0000"
      ! Convert to GRAPES time format: YYYYMMDDHH -> YYYY-MM-DD_HH:MM:SS.0000
      grptime = thisfile(fname_len-9:fname_len-6)//'-'// &
        thisfile(fname_len-5:fname_len-4)//'-'// &
        thisfile(fname_len-3:fname_len-2)//'_'// &
        thisfile(fname_len-1:fname_len) // ":00:00.0000"
      CALL mm5_to_lapstime(grptime,i4time)
      
      PRINT*,'GRAPES file time string: ', grptime, ' -> i4time: ', i4time, i4time_needed

      IF (i4time .EQ. i4time_needed) THEN
        print *, "File matches time needed: ",i, TRIM(thisfile)
        nfiles = 1
        i4times_found(1) = i4time
        filenames(1)     = thisfile
        IF (i .GT. 1) THEN
          filenames(2)  = all_files(i-1) ! Modified by Wei-Ting (130326) to get previous time
        ENDIF
        filenames(4) = all_files(max_files) ! the namelist file is stored in the last slot of all_files
        RETURN
      ELSE
       ! delta_t = i4time_needed - i4time
        delta_t = i4time - i4time_needed
        IF (delta_t .GT. 0) THEN
          IF (delta_t .LT. delta_t_1) THEN
            delta_t_1 = delta_t
            i4times_found(1) = i4time
            filenames(1) = thisfile
          ENDIF
        ELSE
          IF (delta_t .GT. delta_t_2) THEN
            delta_t_2 = delta_t
            i4times_found(2) = i4time
            filenames(2) = thisfile
            IF (i .GT. 1) THEN
               filenames(3) = all_files(i-1) ! added by Wei-Ting (130326) to get previous time
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  filenames(4) = all_files(max_files)
  
  ! If we are here, then we did not find an exact
  ! match.  We need to check if we got appropriate
  ! bounding files
 
  IF ( (i4times_found(1)-i4time_needed .LT. dtmax).AND. & ! add -i4time_needed by Wei-Ting (130326)
       (i4times_found(2)-i4time_needed .GT. dtmin) ) THEN ! add -i4time_needed by Wei-Ting (130326)
     nfiles = 2
  ELSE
     nfiles = 0
  ENDIF
  RETURN
END SUBROUTINE get_acceptable_grapes