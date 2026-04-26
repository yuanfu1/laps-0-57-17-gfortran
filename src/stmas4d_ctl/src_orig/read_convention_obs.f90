!doc====================================================================
!doc This routine is designed to read all conventional observation data
!doc and is a private routine for observations module.
!
! History:
!  Creation: Dec. 2009 by Yuanfu Xie
!doc====================================================================

SUBROUTINE read_convetion_obs(source,obs,num_obs)

  CHARACTER, INTENT(IN) :: source      ! data source to be used
  INTEGER,  INTENT(OUT) :: num_obs     ! Total conventional obs read
  REAL,     INTENT(OUT) :: obs(:,:)    ! Conventional bservations read
                                       ! First dimension: observation,
                                       ! location/time;
                                       ! Second dimension: num obs.

  SELECT CASE(source)

  CASE ('LAPS', 'laps')

    CALL get_convention_laps_obs(obs,num_obs)

  CASE DEFAULTS
    WRITE(6,*) 'read_convention_obs: undefined data source: ',source
    STOP
  END SELECT

  RETURN

END SUBROUTINE read_convention_obs
