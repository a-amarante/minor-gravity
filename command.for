c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMMAND.FOR    (FEG  28 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Get overall CPU % usage.
c
c------------------------------------------------------------------------------
c
      subroutine command ( string, ex )
c
      implicit none
c
c      character*8 string / 'ls s*' /
      character*1000 string
      integer ex
      INTEGER*4 status, system
c
      status = system( string(1:ex) )
      if ( status .ne. 0 ) stop 'system: error'
c
      return
      end
c
c------------------------------------------------------------------------------
c
