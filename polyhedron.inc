c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      POLYHEDRON.INC    (12 July 2013)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author:
c
c Andre Amarante Luiz (A. Amarante) - a.amarante.br@gmail.com
c
c Parameters that you may want to alter at some point:
c
c novmax  = maximum number of vertices
c nopmax  = maximum number of faces
c noed = maximum number of vertices per face
c CMAX  = maximum number of dump interval
c NMESS = maximum number of messages in message.in
c NFILES = maximum number of files that can be open at the same time
c nocen = maximum number of centroids
c TINY  = a small number
c
      integer novmax, nopmax, noed, CMAX, NMESS, NFILES, nocen
      real*8 TINY
c
      parameter (novmax = 200000)
      parameter (nopmax = 200000)
      parameter (noed = 10)
      parameter (CMAX = 1000000)
      parameter (NMESS = 130)
      parameter (NFILES = 1)
      parameter (nocen = 100000)
      parameter (TINY = 4.D-15)
c
c------------------------------------------------------------------------------
c
c Constants:
c
c DR = conversion factor from degrees to radians
c gc = gravitational constant in km^3/s^2/kg
c
      real*8 pi,DR
c
      parameter (pi = 3.141592653589793d0)
      parameter (DR = pi / 180.d0)
c      parameter (gc = 6.6719842229637d-20)
