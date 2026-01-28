c Subroutine HMS converts TIME in seconds into days, hours, mins,secs. 
      SUBROUTINE DHMS(TIME, NDAYS, NHOURS, MINS, SECS) 
      REAL*8 TIME, TIME2, SECS 
      INTEGER NDAYS, NHOURS, MINS 
      TIME2  = DABS(TIME)
      NDAYS  = INT(TIME2 / 86400.0)
      TIME2  = TIME2 - 86400.0 * NDAYS 
      NHOURS = INT(TIME2 / 3600.0) 
      TIME2  = TIME2 - 3600.0 * NHOURS 
      MINS   = INT(TIME2 / 60.0) 
      SECS   = TIME2 - 60.0 * MINS 
      END
