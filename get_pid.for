c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      GET_PID.FOR    (FEG  29 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Get PID process of a program.
c
c------------------------------------------------------------------------------
c
      subroutine get_pid (pid)
c
      implicit none
c
      INTEGER*4 getpid, getuid, getgid
      INTEGER pid, uid, gid
c
      pid = getpid()
c      uid = getuid()
c      gid = getgid()
c
c      print *, "The current process ID is ", pid
c      print *, "Your numerical user ID is ", uid
c      print *, "Your numerical group ID is ", gid
c
      end
c
c------------------------------------------------------------------------------
c
