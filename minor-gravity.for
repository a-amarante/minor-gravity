c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MINOR-GRAVITY.FOR    (UEMS   20 Ago 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Minor-Gravity package computes numerically the irregular gravity field around
c minor bodies (such as Asteroids and Comets), through polyhedron or mascons
c techniques.
c
c******************************************************************************
c
c Copyright notice for SEG distribution:
c
c Copyright (c) 2012 by the Society of Exploration Geophysicists.
c For more information, go to http://software.seg.org/2012/0001 .
c You must read and accept usage terms at:
c http://software.seg.org/disclaimer.txt before use.
c
c-------------------------------------------------------------------------------
c
c	PROGRAMM "polyhedron.f"
c
c	D.Tsoulis                        Thessaloniki, June 2010
c	
c
c	This programm computes the potential, its first
c	and second derivatives of a homogenous polyhedron
c	according to Petrovic (J of G, 1996). The triple
c	integrals of V, Vi and Vij (i,j=1..3) are transformed
c	twice by means of the divergence theorem of Gauss. The
c	transission from volume integrals into line integrals
c	is accomplished in two steps as follows:
c
c	                      GAUSS
c	1. Volume Integral  ----------> Surface Integral
c	2. Surface Integral ----------> Line Integral
c
c
c	Literature:
c	-----------
c	1. Petrovic S. (1996): Determination of the potential of
c	     homogeneous polyhedral bodies using line integrals, 
c	     Journal of Geodesy 71, 44 - 52.
c	2. Werner R.A. and D.J. Scheeres (1997): Exterior gravitation
c	     of a polyhedron derived and compared with harmonic and
c	     mascon gravitation representations of asteroid 4769
c	     Castalia, Celestial Mechanics and Dynamical Astronomy 65,
c	     313 - 344.
c
c-------------------------------------------------------------------------------
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      character*80 infile(6),outfile(1),dumpfile(2),mem(NMESS)
      integer*4 nograd,nflush,opflag,l0,nstored,l1,l2
c
      integer nov,nop,noc,noe,k,lmem(NMESS),opt(7)
      real*8 xp,yp,zp,omg,factor
c
      integer i, j, l
      integer novl,nopl,noel,kl
      integer i0,ilay
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 xc(nocen),yc(nocen),zc(nocen),mc(nocen),dens(nocen)
      real*8 tmp1, tmp2, s_1, s2, s_3, s_5
      real*8 tmpxx, tmpyy, tmpzz, tmpxy, tmpyz, tmpxz
      real*8 xl(novmax),yl(novmax),zl(novmax)
      real*8 vx,vy,vz,del,inc
      real*8 a(3),v,vij,C,grav,phi,theta,omega,sumlapl
      real*8 dx,dy,dz
      real*8 eixa,eixb,eixc,test,a2,b2,c2
      dimension vij(3,3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
      real*8 xi,xf,yi,yf,zi,zf,incx,incy,incz
      real*8 treal,tuser,tsys
      real*8 gc
      real*8 cubx,cuby,cubz
      real*8 xo,yo,zo
c
c-------------------------------------------------------------------------------
c
c chama a subrotina mio_in para verificar os nomes e as existências dos arquivos de entrada (.in), saída (.out) e despejo (.dmp), criar os arquivos de saída e despejo e prosseguir com a integração no ponto onde ela foi interrompida
      call mio_in (infile,outfile,dumpfile,mem,lmem,nov,nop,noc,
     %  noe,k,dens,omg,xv,yv,zv,xc,yc,zc,mc,opt,opflag,nograd,
     %  nflush,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,incx,incy,
     %  incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
c
      opflag = opflag + 1
c
      call timestamp2( )
      if (opflag.gt.0) then
        call mio_elapse (treal,tuser,tsys,82,mem,lmem)
        call mio_remaining (treal,opflag,nograd,mem,lmem)
      end if
      call mio_sizeout (opt(1),opflag,nograd,nflush,outfile,
     %  mem,lmem)
      call mio_host ( )
c
      if (opt(4).eq.1) then
c
c abre o arquivo dos pontos em que se pretende calcular o potencial e suas derivadas (13)
 450    open  (10+3, file=infile(3), status='old', err=450)
c
c posiciona a função read para ler o ponto inicial do próximo intervalo a partir do qual se pretendia calcular o potencial e suas derivadas quando a integração foi interrompida
c        opflag = opflag + 1
        do 20 l0=0,opflag-1
          read(10+3,*) xp,yp,zp
  20    continue
      else
        if (incx.eq.0.d0.and.xi.ne.xf) call mio_err (6,mem(1),lmem(1),
     %    mem(97),lmem(97),' ',1,infile(4),80)
        if (incy.eq.0.d0.and.yi.ne.yf) call mio_err (6,mem(1),lmem(1),
     %    mem(98),lmem(98),' ',1,infile(4),80)
        if (incz.eq.0.d0.and.zi.ne.zf) call mio_err (6,mem(1),lmem(1),
     %    mem(99),lmem(99),' ',1,infile(4),80)
      end if
c
c lê os pontos em que se pretende calcular o potencial e suas derivadas iniciando sempre do sucessor do último ponto do último intervalo em que a integração parou
      nstored = 0
c
      if (eixa.eq.0.d0) then
        a2 = 0.d0
      else
        a2 = 1.d0/eixa
      end if
      if (eixb.eq.0.d0) then
        b2 = 0.d0
      else
        b2 = 1.d0/eixb
      end if
      if (eixc.eq.0.d0) then
        c2 = 0.d0
      else
        c2 = 1.d0/eixc
      end if
c
      do 2013 l0=opflag,nograd-1
	if (opt(4).eq.1) read(10+3,*,end=2020) xp,yp,zp
c
        l1 = -1
        if (opt(4).eq.0) xp = xi
        if (opt(4).eq.1) xf = xp
        if (opt(4).eq.1) xi = xp
        do while (xp.le.xf+TINY)
          if (opt(4).eq.0) yp = yi
          if (opt(4).eq.1) yf = yp
          if (opt(4).eq.1) yi = yp
          do while (yp.le.yf+TINY)
            if (opt(4).eq.0) zp = zi
            if (opt(4).eq.1) zf = zp
            if (opt(4).eq.1) zi = zp
            do while (zp.le.zf+TINY)
c
            l1 = l1 + 1
            if (opt(4).eq.0) l2 = l1
            if (opt(4).eq.1) l2 = l0
c
        test = xp*xp*a2 + yp*yp*b2 + zp*zp*c2 -1.d0
c        if (a2.eq.0.d0.and.b2.eq.0.d0.and.c2.eq.0.d0) test = 1.1d0
        if ((test.gt.0.d0.or.test.eq.-1.d0).and.l2.ge.opflag) then
c no programa original o potencial e suas derivadas são sempre calculados na origem, isto é, no ponto de coordenadas (0,0,0). Logo, para calcularmos o potencial e suas derivadas em qualquer outro ponto P devemos ter nossa origem nesse novo ponto P e para isso basta transladarmos todas as coordenadas dos pontos do poliedro, subtraindo-as das coordenadas do ponto P
c
        a(1) = 0.d0
        a(2) = 0.d0
        a(3) = 0.d0
        v    = 0.d0
        vij(1,1)= 0.d0
        vij(2,2)= 0.d0
        vij(3,3)= 0.d0
        vij(1,2)= 0.d0
        vij(1,3)= 0.d0
        vij(2,3)= 0.d0
c
c polyhedron method
        if (opt(3).eq.0) then
          if (opt(2).eq.1) then
            call polyhedron (nov,nop,xv,yv,zv,noe,k,
     %        dens(opt(2)),xp,yp,zp,vx,vy,vz,v,
     %        vij,gc)
            a(1) = a(1)  -  vx
            a(2) = a(2)  -  vy
            a(3) = a(3)  -  vz
          else
            do i0 = 1, nop
              l = 1
c compute new vertices
              do ilay = 1, opt(2)
                inc = del * ilay
                do j = 1, noe(i0)
                  l = l + 1
                  xl(l) = xv(k(i0,j)) * inc
                  yl(l) = yv(k(i0,j)) * inc
                  zl(l) = zv(k(i0,j)) * inc
                enddo
              enddo
              novl = l
c compute new faces
              l = 1
c compute top face of pyramid
              do j = 1, noe(i0)
                kl(l,j) = j + 1
              enddo
              noel(l) = noe(i0)
c compute lateral faces of pyramid
              do j = 1, noe(i0)-1
                l = l + 1
                kl(l,1) = 1
                kl(l,2) = j + 2
                kl(l,3) = j + 1
                noel(l) = 3
              enddo
c compute last lateral face of pyramid
              l = l + 1
              kl(l,1) = 1
              kl(l,2) = 2
              kl(l,3) = noe(i0) + 1
              noel(l) = 3
c compute pyramid's potential
              nopl = l
              call polyhedron (novl,nopl,xl,yl,zl,noel,kl,
     %          dens(1),xp,yp,zp,vx,vy,vz,v,vij,gc)
              a(1) = a(1)  -  vx
              a(2) = a(2)  -  vy
              a(3) = a(3)  -  vz
c compute faces of pyramidal frustums
              do ilay = 2, opt(2)
                l = 1
c compute bottom face (reverse)
                do j = 1, noe(i0)
                  kl(l,j) = (ilay-2)*noe(i0) + (noe(i0)-j) + 2
                enddo
                noel(l) = noe(i0)
c compute lateral faces
                do j = 1, noe(i0)-1
                  l = l + 1
                  kl(l,1) = (ilay-2)*noe(i0) + j + 1
                  kl(l,2) = kl(l,1) + 1
                  kl(l,3) = ((ilay-1)*noe(i0) + j) + 2
                  kl(l,4) = kl(l,3) - 1
                  noel(l) = 4
                enddo
c compute last lateral face
                l = l + 1
                kl(l,1) = (ilay-1)*noe(i0) + 1
                kl(l,2) = kl(l,1) - (noe(i0)-1)
                kl(l,3) = ilay*noe(i0) - (noe(i0)-1) + 1
                kl(l,4) = ilay*noe(i0) + 1
                noel(l) = 4
c compute top face (no reverse)
                l = l + 1
                do j = 1, noe(i0)
                  kl(l,j) = (ilay-1)*noe(i0) + j + 1
                enddo
                noel(l) = noe(i0)
c compute pyramidal frustum's potential
                nopl = l
                call polyhedron (novl,nopl,xl,yl,zl,noel,kl,
     %            dens(ilay),xp,yp,zp,vx,vy,vz,v,vij,gc)
                a(1) = a(1)  -  vx
                a(2) = a(2)  -  vy
                a(3) = a(3)  -  vz
              enddo
            enddo
          end if
c          v = -v
c          vij(1,1) = -vij(1,1)
c          vij(1,2) = -vij(1,2)
c          vij(1,3) = -vij(1,3)
c          vij(2,1) = -vij(2,1)
c          vij(2,2) = -vij(2,2)
c          vij(2,3) = -vij(2,3)
c          vij(3,1) = -vij(3,1)
c          vij(3,2) = -vij(3,2)
c          vij(3,3) = -vij(3,3)
c
c mascon method
        else
          do j = 1, noc
            dx = xp - xc(j)
            dy = yp - yc(j)
            dz = zp - zc(j)
            s2 = dx*dx + dy*dy + dz*dz
            s_1 = 1.d0 / sqrt(s2)
            s_3 = s_1 * s_1 * s_1
            s_5 = s_3 * s_1 * s_1
            s_5 = 3.d0 * s_5
            tmp1 = s_3 * mc(j)
            tmp2 = s_1 * mc(j)
            tmpxx = s_5 * dx * dx
            tmpyy = s_5 * dy * dy
            tmpzz = s_5 * dz * dz
            tmpxy = s_5 * dx * dy
            tmpyz = s_5 * dy * dz
            tmpxz = s_5 * dx * dz
            a(1) = a(1)  +  tmp1 * dx
            a(2) = a(2)  +  tmp1 * dy
            a(3) = a(3)  +  tmp1 * dz
c            v = v - tmp2
c            vij(1,1) = vij(1,1) + mc(j) * (s_3 - tmpxx)
c            vij(2,2) = vij(2,2) + mc(j) * (s_3 - tmpyy)
c            vij(3,3) = vij(3,3) + mc(j) * (s_3 - tmpzz)
c            vij(1,2) = vij(1,2) - mc(j) * tmpxy
c            vij(1,3) = vij(1,3) - mc(j) * tmpxz
c            vij(2,3) = vij(2,3) - mc(j) * tmpxy
            v = v + tmp2
            vij(1,1) = vij(1,1) + mc(j) * (tmpxx - s_3)
            vij(2,2) = vij(2,2) + mc(j) * (tmpyy - s_3)
            vij(3,3) = vij(3,3) + mc(j) * (tmpzz - s_3)
            vij(1,2) = vij(1,2) + mc(j) * tmpxy
            vij(1,3) = vij(1,3) + mc(j) * tmpxz
            vij(2,3) = vij(2,3) + mc(j) * tmpxy
          end do
          a(1) = -a(1)
          a(2) = -a(2)
          a(3) = -a(3)
          vij(2,1) = vij(1,2)
          vij(3,1) = vij(1,3)
          vij(3,2) = vij(2,3)
        end if
c
        vx = a(1)
        vy = a(2)
        vz = a(3)
c
        xo=xp
        yo=yp
        zo=zp
c        if (opt(4).eq.0) then
c          xo=zp
c          yo=yp
c          zo=xp
c        end if
c
c Output data for all bodies
        call mio_out (opt,nstored,l2,xo,yo,zo,v,vx,vy,vz,vij,
     %    omg,nflush,nograd,outfile,dumpfile,mem,lmem,nov,
     %    nop,noc,dens,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,
     %    zf,incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
c
        end if
c
              if (l1.ge.nograd-1) goto 2222
              if (zi.eq.zf) exit
              zp = zp + dabs(incz)
            end do
            if (yi.eq.yf) exit
            yp = yp + dabs(incy)
          end do
          if (xi.eq.xf) exit
          xp = xp + dabs(incx)
        end do
c
        if (opt(4).eq.0) goto 2222
c
c fim da leitura dos pontos em que se pretende calcular o potencial e suas derivadas
2013  continue
c
c fecha o arquivo de entrada dos pontos em que se pretende calcular o potencial e suas derivadas (13)
2020  if (opt(4).eq.1) close (10+3)
c
      nograd = l0
      l1     = l0
2222  nograd = l1
c
      if (nstored.ne.0) then
        xo=xp
        yo=yp
        zo=zp
c        if (opt(4).eq.0) then
c          xo=zp
c          yo=yp
c          zo=xp
c        end if
        call mio_out (opt,nstored-1,l1-1,xo,yo,zo,v,vx,vy,vz,vij,
     %    omg,nflush,nograd,outfile,dumpfile,mem,lmem,nov,nop,noc,
     %    dens,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,incx,incy,
     %    incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
      end if
c
      call timestamp2( )
c
c termina a execução do programa
      write (*,'(a)') mem(10)(1:lmem(10))
      stop
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_IN.FOR    (ErikSoft   4 May 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 12 July 2013)
c
c Reads names, masses, coordinates and velocities of all the bodies,
c and integration parameters for the MERCURY integrator package. 
c If DUMPFILE(4) exists, the routine assumes this is a continuation of
c an old integration, and reads all the data from the dump files instead
c of the input files.
c
c N.B. All coordinates are with respect to the central body!!
c ===
c
c------------------------------------------------------------------------------
c
      subroutine mio_in (infile,outfile,dumpfile,mem,lmem,nov,
     %  nop,noc,noe,k,dens,omg,xv,yv,zv,xc,yc,zc,mc,opt,opflag,
     %  nograd,nflush,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,
     %  incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      character*80 infile(6),outfile(1),dumpfile(2),mem(NMESS)
      integer*4 nograd,nflush,opflag
      integer nov,nop,noc,noe,k,lmem(NMESS),opt(7)
      real*8 xv(novmax),yv(novmax),zv(novmax),omg,factor
      real*8 xc(nocen),yc(nocen),zc(nocen),mc(nocen),dens(nocen)
      real*8 eixa,eixb,eixc
      dimension noe(nopmax)
      dimension k(nopmax,noed)
      real*8 xi,xf,yi,yf,zi,zf,incx,incy,incz
      real*8 treal,tuser,tsys
      real*8 gc
c
c Local
      character*3 c3
      character*1 c1
      character*80 filename,c80
      logical test,oldflag
      integer j,i,l,i0,lim(2,1000),nsub,lineno
      character*15000 string
      real*8 T
      real*8 mcen,rcen,vc,xct,yct,zct
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1b(3),T2(3),TP(3)
      real*8 J0(3,3),JC(3,3),eigV(3,3),eigVa(3)
      real*8 massd,C20,C22,eA,eB,eC
      integer kkk
      real*8 dot,nx,ny,nz
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,len
      integer signal1(nopmax),signal2(nopmax)
      real*8 face1(nopmax),face2(nopmax)
      real*8 cubx,cuby,cubz
      dimension kkk(nopmax,noed)
      integer giulia
      real*8 mcenb,xctb,yctb,zctb
c
c------------------------------------------------------------------------------
c
      treal = 0.d0
      tuser = 0.d0
      tsys  = 0.d0
c
      do j = 1, 80
        filename(j:j) = ' '
      end do
      do j = 1, 6
        infile(j)   = filename
      end do
      do j = 1, 1
        outfile(j)  = filename
      end do
      do j = 1, 2
        dumpfile(j) = filename
      end do
c
      opt(3) = 0
      opt(4) = 0
      opt(5) = 0
      opt(6) = 0
      opt(7) = 0
c
      call system('mkdir -p out')
c      call system('mkdir -p dmp')
c
c Read in output messages
c      inquire (file='message.in', exist=test)
c      if (.not.test) then
c        write (*,'(/,2a)') ' ERROR: This file is needed to start',
c     %    ' the integration:  message.in'
c        stop
c      end if
c      open (16, file='message.in', status='old')
c  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
c      goto 10
c  20  close (16)
c
c Read in filenames and check for duplicate filenames
      inquire (file='files.in', exist=test)
c      if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88),
c     %  ' ',1,'files.in',8)
      if (.not.test) then
        write (*,'(/,2a)') ' ERROR: This file is needed to start',
     %    ' the integration:  files.in'
        stop
      end if
c
      open (15, file='files.in', status='old')
      do j = 1, 5
 400    read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 400
        if (j.eq.5) then
          infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        end if
      end do
      close(15)
c
c Read in output messages
      inquire (file=infile(5), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(5)
        stop
      end if
      open (16, file=infile(5), status='old')
  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
      goto 10
  20  close (16)
c
      open (15, file='files.in', status='old')
c
c Input files
      do j = 1, 6
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do i = 1, j - 1
          if (infile(j).eq.infile(i)) call mio_err (6,mem(1),lmem(1),
     %      mem(3),lmem(3),infile(j),80,mem(4),lmem(4))
        end do
      end do
c
c Output files
      do j = 1, 1
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        outfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do i = 1, j - 1
          if (outfile(j).eq.outfile(i)) call mio_err (6,mem(1),
     %      lmem(1),mem(3),lmem(3),outfile(j),80,mem(4),lmem(4))
        end do
        do i = 1, 6
          if (outfile(j).eq.infile(i)) call mio_err (6,mem(1),lmem(1),
     %      mem(3),lmem(3),outfile(j),80,mem(4),lmem(4))
        end do
      end do
c
c Dump files
      do j = 1, 2
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        dumpfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do i = 1, j - 1
          if (dumpfile(j).eq.dumpfile(i)) call mio_err (6,mem(1),
     %      lmem(1),mem(3),lmem(3),dumpfile(j),80,mem(4),lmem(4))
        end do
        do i = 1, 4
          if (dumpfile(j).eq.infile(i)) call mio_err (6,mem(1),
     %      lmem(1),mem(3),lmem(3),dumpfile(j),80,mem(4),lmem(4))
        end do
        do i = 1, 1
          if (dumpfile(j).eq.outfile(i)) call mio_err (6,mem(1),
     %      lmem(1),mem(3),lmem(3),dumpfile(j),80,mem(4),lmem(4))
        end do
      end do
      close (15)
c
c Find out if this is an old integration (i.e. does the restart file exist)
      inquire (file=dumpfile(2), exist=oldflag)
c
c------------------------------------------------------------------------------
c
c  READ  IN  INTEGRATION  PARAMETERS
c
c Check if the file containing integration parameters exists, and open it
      filename = infile(4)
      if (oldflag) filename = dumpfile(1)
      inquire (file=filename, exist=test)
      if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),lmem(2),
     %  ' ',1,filename,80)
  30  open  (14, file=filename, status='old', err=30)
c
c Read integration parameters
      lineno = 0
      do j = 1, 26
  40    lineno = lineno + 1
        read (14,'(a15000)') string
        if (string(1:1).eq.')') goto 40
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 40
        c80 = string(lim(1,nsub):lim(2,nsub))
c        if (j.eq.1) read (c80,*,err=661) nov
c        if (j.eq.2) read (c80,*,err=661) nop
        if (j.eq.1) read (c80,*,err=661) nograd
        if (j.eq.2) read (c80,*,err=661) nflush
        if (j.eq.3) then
          c1 = c80(1:1)
          if(c1.eq.'l'.or.c1.eq.'L') then
            opt(1) = 1
          else if (j.eq.3.and.(c1.eq.'m'.or.c1.eq.'M')) then
            opt(1) = 2
          else if (j.eq.3.and.(c1.eq.'h'.or.c1.eq.'H')) then
            opt(1) = 3
          else
            goto 661
          end if
        end if
        if (j.eq.4.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(3) = 1
        if (j.eq.5.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(6) = 1
        if (j.eq.6.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(7) = 1
        if (j.eq.7) read (c80,*,err=661) cubx
        if (j.eq.8) read (c80,*,err=661) cuby
        if (j.eq.9) read (c80,*,err=661) cubz
        if (j.eq.10.and.(c80(1:1).eq.'y'.or.c80(1:1).eq.'Y')) opt(4) = 1
        if (j.eq.11) read (c80,*,err=661)  xi
        if (j.eq.12) read (c80,*,err=661)  xf
        if (j.eq.13) read (c80,*,err=661) yi
        if (j.eq.14) read (c80,*,err=661) yf
        if (j.eq.15) read (c80,*,err=661) zi
        if (j.eq.16) read (c80,*,err=661) zf
        if (j.eq.17) read (c80,*,err=661) incx
        if (j.eq.18) read (c80,*,err=661) incy
        if (j.eq.19) read (c80,*,err=661) incz
        if (j.eq.20) then
          l = nsub
          i0 = 1
          c80 = string(lim(1,l):lim(2,l))
          do while (c80(1:7).ne.'deepest')
            read (c80,*,err=661) dens(i0)
            l  = l  - 1
            i0 = i0 + 1
            c80 = string(lim(1,l):lim(2,l))
          end do
          opt(2) = i0 - 1
        end if
        if (j.eq.21) read (c80,*,err=661) T
        if (j.eq.22) read (c80,*,err=661) factor
        if (j.eq.23) read (c80,*,err=661) eixa
        if (j.eq.24) read (c80,*,err=661) eixb
        if (j.eq.25) read (c80,*,err=661) eixc
        if (j.eq.26) read (c80,*,err=661) gc
      end do
c
c tira o módulo dos parâmetros para o caso de ser digitado algum valor negativo
c      nov = abs(nov)
c      nop = abs(nop)
      nograd = abs(nograd)
      nflush = abs(nflush)
      do l = 1, opt(2)
        dens(l) = dabs(dens(l))
      end do
      omg = 2.0d0*pi/(T*3600.0d0)
      factor = dabs(factor)
      eixa = eixa * eixa
      eixb = eixb * eixb
      eixc = eixc * eixc
      cubx = dabs(cubx)
      cuby = dabs(cuby)
      cubz = dabs(cubz)
c      if (opt(6).eq.1.or.opt(7).eq.1) opt(3) = 1
      if (opt(6).eq.1) opt(3) = 1
      close (14)
c
c Check if the file containing integration parameters exists, and open it
      if (opt(4).eq.1) then
        filename = infile(3)
        if (oldflag) filename = infile(3)
        inquire (file=filename, exist=test)
        if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),lmem(2),
     %    ' ',1,filename,80)
      end if
c
c------------------------------------------------------------------------------
c
      write(6,'(3x,a34,a23)') 'Reading Polyhedron information and',
     %  ' computing centroids...'
c
c      if (nov.gt.novmax) call mio_err (6,mem(1),lmem(1),mem(16),
c     %  lmem(16),' ',1,mem(14),lmem(14))
c      if (nop.gt.nopmax) call mio_err (6,mem(1),lmem(1),mem(17),
c     %  lmem(17),' ',1,mem(14),lmem(14))
c
      call readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,
     %  noe,k,norm,wfac,factor)
c
c      call cross_ast (nov,nop,noe,k,xv,yv,zv,face1)
c
c      call masc_mass (nov,nop,xv,yv,zv,noe,k,dens(opt(2)),
c     %   mc,mcen,vc)
c      write(*,'(2(1x,1p,e35.25))') vc,mcen
c
c      call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %  mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c
c      giulia = 0
c      do j = 1, opt(2) - 1
c        do i = j + 1, opt(2)
c          if (dens(j).ne.dens(i)) giulia = 1
c        end do
c      end do
c
      if (opt(7).eq.0) then
c        if (giulia.eq.0) then
          call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,k,
     %      norm,wfac,T0,T1b,T2,TP)
          vc = T0
          mcen = dens(1) * vc
          call compcenpolyhedron (dens(1),mcen,T0,T1b,T2,TP,
     %      xct,yct,zct,J0)
          noc = 1
c        else if (giulia.eq.1) then
c          if (opt(6).eq.0) then
c            call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %        mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c          else if (opt(6).eq.1) then
c            call masc_layer2 (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %        mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,
c     %        cuby,cubz)
c          end if
c        end if
      else if (opt(7).eq.1) then
        write(6,'(3x,a35,a15)') 'Reading Polyhedron information from',
     %    ' cube file...'
c        call readbinpol2b (infile,mem,lmem,noc,mcenb,vc,J0,
c     %    xc,yc,zc,mc,xct,yct,zct,mcen)
        call readbinpol2c (infile,mem,lmem,noc,mcen,vc,J0,
     %    xc,yc,zc,mc,xct,yct,zct,T1b,T2,TP)
c        xctb = 0.d0
c        yctb = 0.d0
c        zctb = 0.d0
        do j = 1, noc
c          xctb = xctb + mc(j) * xc(j)
c          yctb = yctb + mc(j) * yc(j)
c          zctb = zctb + mc(j) * zc(j)
          mc(j) = mc(j) * gc
        end do
c        xctb = xctb / mcenb
c        yctb = yctb / mcenb
c        zctb = zctb / mcenb
      end if
c
c      write(*,'(2(1x,1p,e35.25),i)') vc,mcen,noc
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
c computes total center of mass of the polyhedron
c      if (opt(7).eq.0.and.giulia.eq.1) then
c      xct = 0.d0
c      yct = 0.d0
c      zct = 0.d0
c      do j = 1, noc
c        xct = xct + mc(j) * xc(j)
c        yct = yct + mc(j) * yc(j)
c        zct = zct + mc(j) * zc(j)
c        mc(j) = mc(j) * gc
c      end do
c      xct = xct / mcen
c      yct = yct / mcen
c      zct = zct / mcen
c      end if
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c
c      call compVolumeIntegrals (nov,nop,nop,xv,yv,zv,noe,
c     %  k,norm,wfac,T0,T1b,T2,TP)
c      write(*,'(2(1x,1p,e35.25))') vc,T0
c      call compcenpolyhedron (dens(opt(2)),mcen,T0,T1b,
c     %  T2,TP,xct,yct,zct,J0)
c      write(*,'(3(1x,1p,e35.25))') xct,yct,zct
c      write(*,'(1x,1p,e35.25)') T1b(1)
c      write(*,'(1x,1p,e35.25)') T1b(2)
c      write(*,'(1x,1p,e35.25)') T1b(3)
c      write(*,'(1x,1p,e35.25)') T2(1)
c      write(*,'(1x,1p,e35.25)') T2(2)
c      write(*,'(1x,1p,e35.25)') T2(3)
c      write(*,'(1x,1p,e35.25)') TP(1)
c      write(*,'(1x,1p,e35.25)') TP(2)
c      write(*,'(1x,1p,e35.25)') TP(3)
c      write(*,'(3(1x,1p,e35.25))') J0(1,1),J0(1,2),J0(1,3)
c      write(*,'(3(1x,1p,e35.25))') J0(2,1),J0(2,2),J0(2,3)
c      write(*,'(3(1x,1p,e35.25))') J0(3,1),J0(3,2),J0(3,3)
c
      call masc_trans (J0,xct,yct,zct,mcen,JC)
c
      rcen = ((3.0/(4.0*pi))*vc)**(1.0/3.0)
c
      C20 = -1.d0/(2.d0*mcen)*(2.d0*JC(3,3)-JC(1,1)-JC(2,2))
      C22 =  1.d0/(4.d0*mcen)*(JC(2,2)-JC(1,1))
      massd = (JC(2,2)-JC(1,1))/(JC(3,3)-JC(1,1))
c
      if (.not.oldflag) then
        write(*,'(/,a,i7)') mem(111)(1:lmem(111)),nov
        write(*,'(a,i7)') mem(112)(1:lmem(112)),nop
c        if ((opt(7).eq.0.and.giulia.eq.1).or.opt(7).eq.1
c     %    .or.opt(3).eq.0) then
        if (opt(3).eq.0.or.opt(7).eq.1) then
          write(*,'(a,i7)') mem(113)(1:lmem(113)),noc
          write(*,'(a,1p,e22.15)') mem(58)(1:lmem(58)),vc
          write(*,'(a,1p,e22.15)') mem(57)(1:lmem(57)),rcen
        end if
c        if (opt(7).eq.1) then
c        write(*,'(a,1p,e22.15)') mem(56)(1:lmem(56)),mcenb
c        write(*,'(2a,3(1p,e22.15,a))') mem(54)(1:lmem(54)),
c     %    mem(53)(1:lmem(53)),xctb,mem(51)(1:lmem(51)),
c     %    yctb,mem(51)(1:lmem(51)),zctb,mem(52)(1:lmem(52))
c        else
        write(*,'(a,1p,e22.15)') mem(56)(1:lmem(56)),mcen
        write(*,'(2a,3(1p,e22.15,a))') mem(54)(1:lmem(54)),
     %    mem(53)(1:lmem(53)),xct,mem(51)(1:lmem(51)),
     %    yct,mem(51)(1:lmem(51)),zct,mem(52)(1:lmem(52))
c        end if
c        end if
        write(*,'(a)') mem(55)(1:lmem(55))
        write(*,'(3(1x,1p,e22.15))') JC(1,1),JC(1,2),JC(1,3)
        write(*,'(3(1x,1p,e22.15))') JC(2,1),JC(2,2),JC(2,3)
        write(*,'(3(1x,1p,e22.15))') JC(3,1),JC(3,2),JC(3,3)
        write(*,'(a)') mem(101)(1:lmem(101))
        write(*,'(3(1x,1p,e22.15))') JC(1,1)/mcen,JC(2,2)/mcen,
     %    JC(3,3)/mcen
        write(*,'(a)') mem(104)(1:lmem(104))
        write(*,'(2(a,1p,e22.15))') mem(106)(1:lmem(106)),
     %    C20,mem(107)(1:lmem(107)),C22
        write(*,'(a,1p,e22.15)') mem(105)(1:lmem(105)),
     %    massd
      endif
c
      call EigenVectors (6,mem,lmem,JC,eigV,eigVa,oldflag)
c
      eA = sqrt(5.d0*(eigVa(2)+eigVa(3)-eigVa(1))/(2.d0*mcen))
      eB = sqrt(5.d0*(eigVa(1)+eigVa(3)-eigVa(2))/(2.d0*mcen))
      eC = sqrt(5.d0*(eigVa(1)+eigVa(2)-eigVa(3))/(2.d0*mcen))
c
      if (.not.oldflag) then
        write(*,'(3(a,1p,e22.15))') mem(102)(1:lmem(102)),
     %    eA,mem(103)(1:lmem(103)),eB,mem(103)(1:lmem(103)),eC
c
        write(*,'(a)') mem(108)(1:lmem(108))
        write(*,'(3(1x,1p,e22.15))') T1b(1),T1b(2),T1b(3)
        write(*,'(3(1x,1p,e22.15))') T2(1),T2(2),T2(3)
        write(*,'(3(1x,1p,e22.15))') TP(1),TP(2),TP(3)
      endif
c
c signal1 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal1(j) = -1
c        if (dot.gt.(0.d0)) signal1(j) =  1
c      end do
c
c      call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
c      call masc_rot (nov,xv,yv,zv,xct,yct,zct,eigV)
c
      call masc_rot1 (nov,xv,yv,zv,xct,yct,zct)
      call cross_ast (nov,nop,noe,k,xv,yv,zv,face1)
      call masc_rot2 (nov,xv,yv,zv,eigV)
      call cross_ast (nov,nop,noe,k,xv,yv,zv,face2)
c
c      call reorientation (face1,face2,nov,nop,noe,kkk)
      call reorientation2 (face1,face2,nov,nop,noe,k)
c
      if (opt(3).eq.1.and.opt(7).eq.0) then
c        call masc_layer3 (nov,nop,xv,yv,zv,noe,k,dens(1),
c     %  mcenb,vc,xc,yc,zc,mem,lmem,noc,cubx,cuby,cubz,mc)
        call masc_layer4b (nov,nop,xv,yv,zv,noe,k,dens(1),
     %  mcenb,vc,xc,yc,zc,mem,lmem,noc,cubx,cuby,cubz,mc)
        do j = 1, noc
          mc(j) = mc(j) * gc
        end do
        if (.not.oldflag) then
          write(*,'(a,i7)') mem(113)(1:lmem(113)),noc
          write(*,'(a,1p,e22.15)') mem(58)(1:lmem(58)),vc
          rcen = ((3.0/(4.0*pi))*vc)**(1.0/3.0)
          write(*,'(a,1p,e22.15)') mem(57)(1:lmem(57)),rcen
        end if
      end if
c
c      call masc_cen (nop,xv,yv,zv,noe,k,xc,yc,zc)
c
c signal2 comparation
c compute normal faces
c      do j = 1, nop
c        dx1 = xv(k(i,2)) - xv(k(i,1))
c        dy1 = yv(k(i,2)) - yv(k(i,1))
c        dz1 = zv(k(i,2)) - zv(k(i,1))
c        dx2 = xv(k(i,3)) - xv(k(i,1))
c        dy2 = yv(k(i,3)) - yv(k(i,1))
c        dz2 = zv(k(i,3)) - zv(k(i,1))
c        nx = dy1 * dz2 - dy2 * dz1
c        ny = dz1 * dx2 - dz2 * dx1
c        nz = dx1 * dy2 - dx2 * dy1
c        len = dsqrt(nx * nx + ny * ny + nz * nz)
c        nx = nx / len
c        ny = ny / len
c        nz = nz / len
c dot product
c        dot = xv(k(j,1))*nx+yv(k(j,1))*ny+zv(k(j,1))*nz
c        if (dot.lt.(0.d0)) signal2(j) = -1
c        if (dot.gt.(0.d0)) signal2(j) =  1
c      end do
c      do j = 1, nop
c        do l = 1, noe(j)
c          kkk(j,l) = k(j,l)
c        end do
c      end do
c      do j = 1, nop
c        if (signal1(j).ne.signal2(j)) then
c          write(*,*) j,signal1(j),signal2(j)
c          do l = 1, noe(j)
c            k(j,l) = kkk(j,noe(j)-l+1)
c          end do
c        end if
c      end do
c
c      if (opt(3).eq.1) then
c      if (opt(7).eq.0.and.giulia.eq.0) then
c        if (opt(6).eq.0) then
c          call masc_layer (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0)
c        else if (opt(6).eq.1) then
c          call masc_layer2 (nov,nop,xv,yv,zv,noe,k,dens,opt(2),
c     %      mc,mcen,vc,xc,yc,zc,mem,lmem,noc,T1b,T2,TP,J0,cubx,
c     %      cuby,cubz)
c        end if
c        xct = 0.d0
c        yct = 0.d0
c        zct = 0.d0
c        do j = 1, noc
c          xct = xct + mc(j) * xc(j)
c          yct = yct + mc(j) * yc(j)
c          zct = zct + mc(j) * zc(j)
c          mc(j) = mc(j) * gc
c        end do
c        xct = xct / mcen
c        yct = yct / mcen
c        zct = zct / mcen
c        rcen = ((3.0/(4.0*PI))*vc)**(1.0/3.0)
c        if (.not.oldflag) then
c        write(*,'(a,i7)') mem(113)(1:lmem(113)),noc
c        write(*,'(a,1p,e22.15)') mem(58)(1:lmem(58)),vc
c        write(*,'(a,1p,e22.15)') mem(57)(1:lmem(57)),rcen
c        write(*,'(a,1p,e22.15)') mem(56)(1:lmem(56)),mcen
c        write(*,'(2a,3(1p,e22.15,a))') mem(54)(1:lmem(54)),
c     %    mem(53)(1:lmem(53)),xct,mem(51)(1:lmem(51)),
c     %    yct,mem(51)(1:lmem(51)),zct,mem(52)(1:lmem(52))
c        end if
c      else if (opt(7).eq.0.and.giulia.eq.1) then
c        call masc_rot (noc,xc,yc,zc,xct,yct,zct,eigV)
c      end if
c      end if
c
c      if (.not.oldflag) then
c        write(*,'(a)') mem(108)(1:lmem(108))
c        write(*,'(3(1x,1p,e22.15))') T1b(1),T1b(2),T1b(3)
c        write(*,'(3(1x,1p,e22.15))') T2(1),T2(2),T2(3)
c        write(*,'(3(1x,1p,e22.15))') TP(1),TP(2),TP(3)
c      endif
c
      write(*,'(/)')
c
c------------------------------------------------------------------------------
c
c  IF  CONTINUING  AN  OLD  INTEGRATION
c
      if (oldflag) then
c Read in energy and angular momentum variables, and convert to internal units
 330    open (32, file=dumpfile(2), status='old', err=330)
        read (32,'(a150)') string
        call mio_spl (150,string,nsub,lim)
c se a subrotina mio_spl não encontrar duas substrings no vetor de caracteres string, então é por que o programa terminou a escrita antes de chegar no final da primeira linha (NÃO escreveu "EOL" no final da primeira linha) do arquivo restart.dmp; daí utiliza o arquivo restart.tmp para ler o índice do último ponto de despejo que se pretendia calcular o potencial antes da integração ser interrompida
        if(nsub.lt.2) then
 331      open (42, file='restart.tmp', status='old', err=331)
          read (42,'(a150)') string
          call mio_spl (150,string,nsub,lim)
          c80 = string(lim(1,1):lim(2,1))
          read (c80,*,err=667) opflag
          read (42,*) treal,tuser,tsys
          close (42)
        else
          c80 = string(lim(1,1):lim(2,1))
          read (c80,*,err=666) opflag
          read (32,*) treal,tuser,tsys
        end if
        close (32)
c abre o arquivo de saída de dados (21) para escrever a mensagem 9 do arquivo message.in nele e na tela
c 433    open  (20+1, file=outfile(1), status='old', access='append',
c     %    err=433)
c        write (20+1,'(/,a,i10,/)') mem(7)(1:lmem(7)),opflag
c        close (20+1)
c termina a execução do programa caso opflag for igual a nograd-1
        if(opflag.eq.nograd-1) then
          write (*,'(a)') mem(10)(1:lmem(10))
          stop
        end if
        write (6,'(/,a,i10,/)') mem(9)(1:lmem(9)),opflag+1
c caso os nomes dos arquivos de saída forem diferentes no arquivo files.in a cada reinício da integração, então cria um novo arquivo de saída diferente a cada reinício da integração
        do j = 1, 1
          inquire (file=outfile(j), exist=test)
          if (j.eq.1.and.(.not.test)) then
            write(6,'(3(1x),a,1x,a6,1x,a,/)') 'WARNING: the file',
     %        outfile(j),'does not exist! Creating it...'
 431        open  (20+j, file=outfile(j), status='new', err=431)
            close (20+j)
          end if
        end do
c
c------------------------------------------------------------------------------
c
c  IF  STARTING  A  NEW  INTEGRATION
c
      else
        opflag = -1
c
c Check that element and close-encounter files don't exist, and create them
        do j = 1, 1
          inquire (file=outfile(j), exist=test)
          if (test) call mio_err (6,mem(1),lmem(1),mem(5),lmem(5),
     %      ' ',1,outfile(j),80)
 430      open  (20+j, file=outfile(j), status='new', err=430)
          close (20+j)
        end do
c
c Check that dump files don't exist, and then create them
        do j = 1, 2
          inquire (file=dumpfile(j), exist=test)
          if (test) call mio_err (6,mem(1),lmem(1),mem(5),lmem(5),
     %      ' ',1,dumpfile(j),80)
 450      open  (30+j, file=dumpfile(j), status='new', err=450)
          close (30+j)
        end do
c chama a subrotina mio_dump para escrever pela primeira vez nos arquivos de despejo (.dmp e .tmp)
        call mio_dump (mem,lmem,nov,nop,nograd,nflush,opt,dens,omg,
     %    opflag,dumpfile,outfile,noc,factor,eixa,eixb,eixc,xi,xf,
     %    yi,yf,zi,zf,incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
        write (*,'(a)') mem(11)(1:lmem(11))
      end if
c
      return
c
c Error reading from the input file containing integration parameters
 661  write (c3,'(i3)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c3,3,
     %  mem(7),lmem(7))
c
c Error reading from the input file for Big or Small bodies
 666  call mio_err (6,mem(1),lmem(1),mem(8),lmem(8),'restart.dmp',11,
     %  ' ',1)
c
c Error reading epoch of Big bodies
 667  call mio_err (6,mem(1),lmem(1),mem(8),lmem(8),'restart.tmp',11,
     %  ' ',1)
c
c------------------------------------------------------------------------------
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_ERR.FOR    (ErikSoft  6 December 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Writes out an error message and terminates Mercury.
c
c------------------------------------------------------------------------------
c
      subroutine mio_err (unit,s1,ls1,s2,ls2,s3,ls3,s4,ls4)
c
      implicit none
c
c Input/Output
      integer unit,ls1,ls2,ls3,ls4
      character*80 s1,s2,s3,s4
c
c------------------------------------------------------------------------------
c
      write (*,'(/,1a)') ' ERROR: Programme terminated.'
c
      write (unit,'(/,3a,/,2a)') s1(1:ls1),s2(1:ls2),s3(1:ls3),
     %  ' ',s4(1:ls4)
      stop
c
c------------------------------------------------------------------------------
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_SPL.FOR    (ErikSoft  14 November 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Given a character string STRING, of length LEN bytes, the routine finds 
c the beginnings and ends of NSUB substrings present in the original, and 
c delimited by spaces. The positions of the extremes of each substring are 
c returned in the array DELIMIT.
c Substrings are those which are separated by spaces or the = symbol.
c
c------------------------------------------------------------------------------
c
      subroutine mio_spl (len,string,nsub,delimit)
c
      implicit none
c
c Input/Output
      integer len,nsub,delimit(2,1000)
      character*1 string(len)
c
c Local
      integer j,k
      character*1 c
c
c------------------------------------------------------------------------------
c
      nsub = 0
      j = 0
      c = ' '
      delimit(1,1) = -1
c
c Find the start of string
  10  j = j + 1
      if (j.gt.len) goto 99
      c = string(j)
      if (c.eq.' '.or.c.eq.'=') goto 10
c
c Find the end of string
      k = j
  20  k = k + 1
      if (k.gt.len) goto 30
      c = string(k)
      if (c.ne.' '.and.c.ne.'=') goto 20
c
c Store details for this string
  30  nsub = nsub + 1
      delimit(1,nsub) = j
      delimit(2,nsub) = k - 1
c
      if (k.lt.len) then
        j = k
        goto 10
      end if
c
  99  continue
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_OUT.FOR    (ErikSoft   13 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 12 July 2013)
c
c Writes output variables for each object to an output file. Each variable
c is scaled between the minimum and maximum possible values and then
c written in a compressed format using ASCII characters.
c The output variables are:
c  r = the radial distance
c  theta = polar angle
c  phi = azimuthal angle
c  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
c                             kinetic energies. (Note that 0 < fv < 1).
c  vtheta = polar angle of velocity vector
c  vphi = azimuthal angle of the velocity vector
c
c If this is the first output (OPFLAG = -1), or the first output since the 
c number of the objects or their masses have changed (OPFLAG = 1), then 
c the names, masses and spin components of all the objects are also output.
c
c N.B. Each object's distance must lie between RCEN < R < RMAX
c ===  
c
c------------------------------------------------------------------------------
c
      subroutine mio_out (opt,nstored,l,xp,yp,zp,v,vx,vy,vz,
     %  vij,omg,nflush,nograd,outfile,dumpfile,mem,lmem,nov,
     %  nop,noc,dens,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,
     %  incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer*4 nstored,l,nflush,nograd,noc
      real*8 xp,yp,zp,v,vx,vy,vz,vij(3,3),omg
      character*80 outfile(1),dumpfile(2),mem(NMESS)
      integer nov,nop,lmem(NMESS),opt(7)
      real*8 dens(nocen),factor
      real*8 eixa,eixb,eixc
      real*8 xi,xf,yi,yf,zi,zf,incx,incy,incz
      real*8 treal,tuser,tsys
      real*8 gc
      real*8 cubx,cuby,cubz
c
c Local
      integer*4 k
      integer len, nchar
      character*80 header
      character*104 c(CMAX)
      character*8 mio_fl2c,mio_re2c
      character*8 fout
c
      character*13 fout2
c
c------------------------------------------------------------------------------
c
      save c
c
c Create the format list, FOUT, used when outputting the orbital elements
      if (opt(1).eq.1) nchar = 3
      if (opt(1).eq.2) nchar = 5
      if (opt(1).eq.3) nchar = 8
c
c Store details of each new close-encounter minimum
      nstored = nstored + 1
      c(nstored)(1:8)   = mio_fl2c (xp)
c
      c(nstored)(9:16)  = mio_fl2c (yp)
c
      c(nstored)(17:24) = mio_fl2c (zp)
c
      c(nstored)(25:32) = mio_fl2c (v)
      c(nstored)(25+nchar-1:25+nchar-1) = c(nstored)(32:32)
c
      c(nstored)(25+nchar:32+nchar)     = mio_fl2c (vx)
      c(nstored)(25+2*nchar-1:25+2*nchar-1) =
     %  c(nstored)(32+nchar:32+nchar)
c
      c(nstored)(25+2*nchar:32+2*nchar) = mio_fl2c (vy)
      c(nstored)(25+3*nchar-1:25+3*nchar-1) =
     %  c(nstored)(32+2*nchar:32+2*nchar)
c
      c(nstored)(25+3*nchar:32+3*nchar) = mio_fl2c (vz)
      c(nstored)(25+4*nchar-1:25+4*nchar-1) =
     %  c(nstored)(32+3*nchar:32+3*nchar)
c
      c(nstored)(25+4*nchar:32+4*nchar) = mio_fl2c (vij(1,1))
      c(nstored)(25+5*nchar-1:25+5*nchar-1) =
     %  c(nstored)(32+4*nchar:32+4*nchar)
c
      c(nstored)(25+5*nchar:32+5*nchar) = mio_fl2c (vij(2,2))
      c(nstored)(25+6*nchar-1:25+6*nchar-1) =
     %  c(nstored)(32+5*nchar:32+5*nchar)
c
      c(nstored)(25+6*nchar:32+6*nchar) = mio_fl2c (vij(3,3))
      c(nstored)(25+7*nchar-1:25+7*nchar-1) =
     %  c(nstored)(32+6*nchar:32+6*nchar)
c
      c(nstored)(25+7*nchar:32+7*nchar) = mio_fl2c (vij(1,2))
      c(nstored)(25+8*nchar-1:25+8*nchar-1) =
     %  c(nstored)(32+7*nchar:32+7*nchar)
c
      c(nstored)(25+8*nchar:32+8*nchar) = mio_fl2c (vij(1,3))
      c(nstored)(25+9*nchar-1:25+9*nchar-1) =
     %  c(nstored)(32+8*nchar:32+8*nchar)
c
      c(nstored)(25+9*nchar:32+9*nchar) = mio_fl2c (vij(2,3))
      c(nstored)(25+10*nchar-1:25+10*nchar-1) =
     %  c(nstored)(32+9*nchar:32+9*nchar)
c
c If required, output the stored close encounter details
      if (nstored.ge.nflush.or.l.ge.nograd-1) then
        len = 25+10*nchar-1
        fout(1:8) = '(a   ,$)'
        if (len.lt.100) write (fout(3:4),'(i2)') len
        if (len.ge.100) write (fout(3:5),'(i3)') len
c
  10    open (21, file=outfile(1), status='old', access='append',err=10)
c
c Compose a header line with time, number of objects and relevant parameters
        header(1:4) = mio_re2c (dble(nstored), 0.d0, 2517630975.99d0)
        header(5:12) = mio_fl2c (omg)
        fout2(1:13) = '(a1,i1,a12,$)'
        write (21,fout2) char(12),opt(1),header(1:12)
        do k = 1, nstored
          write (21,fout) c(k)(1:len)
        end do
        close (21)
        nstored = 0
c chama a subrotina mio_dump para escrever nos arquivos de despejo (.dmp e .tmp)
        call mio_dump (mem,lmem,nov,nop,nograd,nflush,opt,dens,omg,l,
     %    dumpfile,outfile,noc,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,
     %    incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_RE2C.FOR    (ErikSoft  27 June 1999)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts a REAL*8 variable X, where XMIN <= X < XMAX, into an ASCII string
c of 8 characters, using the new format compression: 
c
c X is first converted to base 224, and then each base 224 digit is converted 
c to an ASCII character, such that 0 -> character 32, 1 -> character 33...
c and 223 -> character 255.
c
c ASCII characters 0 - 31 (CTRL characters) are not used, because they
c cause problems when using some operating systems.
c
c------------------------------------------------------------------------------
c
      function mio_re2c (x,xmin,xmax)
c
      implicit none
c
c Input/output
      real*8 x,xmin,xmax
      character*8 mio_re2c
c
c Local
      integer j
      real*8 y,z
c
c------------------------------------------------------------------------------
c
      mio_re2c(1:8) = '        '
      y = (x - xmin) / (xmax - xmin)
c
      if (y.ge.1) then
        do j = 1, 8
          mio_re2c(j:j) = char(255)
        end do
      else if (y.gt.0) then
        z = y
        do j = 1, 8
          z = mod(z, 1.d0) * 224.d0
          mio_re2c(j:j) = char(int(z) + 32)
        end do
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_FL2C.FOR    (ErikSoft  1 July 1998)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts a (floating point) REAL*8 variable X, into a CHARACTER*8 ASCII 
c string, using the new format compression:
c
c X is first converted to base 224, and then each base 224 digit is converted 
c to an ASCII character, such that 0 -> character 32, 1 -> character 33...
c and 223 -> character 255.
c The first 7 characters in the string are used to store the mantissa, and the
c eighth character is used for the exponent.
c
c ASCII characters 0 - 31 (CTRL characters) are not used, because they
c cause problems when using some operating systems.
c
c N.B. X must lie in the range -1.e112 < X < 1.e112
c ===
c
c------------------------------------------------------------------------------
c
      function mio_fl2c (x)
c
      implicit none
c
c Input/Output
      real*8 x
      character*8 mio_fl2c
c
c Local
      integer ex
      real*8 ax,y
      character*8 mio_re2c
c
c------------------------------------------------------------------------------
c
      if (x.eq.0) then
        y = .5d0
      else
        ax = abs(x)
        ex = int(log10(ax))
        if (ax.ge.1) ex = ex + 1
        y = ax*(10.d0**(-ex))
        if (y.eq.1) then
          y = y * .1d0
          ex = ex + 1
        end if
        y = sign(y,x) *.5d0 + .5d0
      end if
c
      mio_fl2c(1:8) = mio_re2c (y, 0.d0, 1.d0)
      ex = ex + 112
      if (ex.gt.223) ex = 223
      if (ex.lt.0) ex = 0
      mio_fl2c(8:8) = char(ex+32)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_DUMP.FOR    (ErikSoft   21 February 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 13 July 2013)
c
c Writes masses, coordinates, velocities etc. of all objects, and integration
c parameters, to dump files. Also updates a restart file containing other
c variables used internally by MERCURY.
c
c------------------------------------------------------------------------------
c
      subroutine mio_dump (mem,lmem,nov,nop,nograd,nflush,opt,dens,omg,
     %  l,dumpfile,outfile,noc,factor,eixa,eixb,eixc,xi,xf,yi,yf,zi,zf,
     %  incx,incy,incz,treal,tuser,tsys,gc,cubx,cuby,cubz)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer*4 l,nflush,nograd,noc
      character*80 dumpfile(2),mem(NMESS),outfile(1)
      integer nov,nop,lmem(NMESS),opt(7)
      real*8 omg,dens(noc),factor
      real*8 eixa,eixb,eixc
      real*8 xi,xf,yi,yf,zi,zf,incx,incy,incz
      real*8 treal,tuser,tsys
      real*8 gc
      real*8 cubx,cuby,cubz
c
c Local
      integer idp,k
      real*4 dtime,etime,t(2)
      real*8 treal2,tuser2,tsys2,treal0,tuser0,tsys0,treal3,tuser3,tsys3
c
c------------------------------------------------------------------------------
c
      save treal0,tuser0,tsys0

      treal2 = dtime ( t )
      tuser2 = t(1)
      tsys2  = t(2)

      if (opt(5).eq.0) then
        treal0 = 0.d0
        tuser0 = 0.d0
        tsys0  = 0.d0
        opt(5) = 1
      end if
      treal3 = etime ( t )
      tuser3 = t(1)
      tsys3  = t(2)
      treal = treal + treal3 - treal0
      tuser = tuser + tuser3 - tuser0
      tsys  = tsys  + tsys3  - tsys0
      treal0= treal3
      tuser0= tuser3
      tsys0 = tsys3
c
c Dump to temporary files (idp=1) and real dump files (idp=2)
      do idp = 1, 2
c
c Dump the integration parameters
  40    if (idp.eq.1) open (31,file='polyhedron.tmp',status='unknown',
     %    err=40)
  45    if (idp.eq.2) open (31, file=dumpfile(1), status='old', err=45)
c
c Important parameters
        write (31,'(a)') mem(20)(1:lmem(20))
        write (31,'(a)') mem(21)(1:lmem(21))
        write (31,'(a)') mem(36)(1:lmem(36))
        write (31,'(a)') mem(22)(1:lmem(22))
        write (31,'(a)') mem(23)(1:lmem(23))
        write (31,'(a)') mem(22)(1:lmem(22))
c        write (31,*) mem(24)(1:lmem(24)),nov
c        write (31,*) mem(25)(1:lmem(25)),nop
        write (31,*) mem(26)(1:lmem(26)),nograd
        write (31,*) mem(27)(1:lmem(27)),nflush
c
c Integration options
        write (31,'(a)') mem(22)(1:lmem(22))
        write (31,'(a)') mem(28)(1:lmem(28))
        write (31,'(a)') mem(22)(1:lmem(22))
        if (opt(1).eq.1) then
          write (31,'(2a)') mem(29)(1:lmem(29)),mem(33)(1:lmem(33))
        else if (opt(1).eq.3) then
          write (31,'(2a)') mem(29)(1:lmem(29)),mem(35)(1:lmem(35))
        else
          write (31,'(2a)') mem(29)(1:lmem(29)),mem(34)(1:lmem(34))
        end if
        if (opt(3).eq.0) then
          write (31,'(2a)') mem(37)(1:lmem(37)),mem(12)(1:lmem(12))
        else
          write (31,'(2a)') mem(37)(1:lmem(37)),mem(13)(1:lmem(13))
        end if
        if (opt(6).eq.0) then
          write (31,'(2a)') mem(116)(1:lmem(116)),mem(12)(1:lmem(12))
        else
          write (31,'(2a)') mem(116)(1:lmem(116)),mem(13)(1:lmem(13))
        end if
        if (opt(7).eq.0) then
          write (31,'(2a)') mem(117)(1:lmem(117)),mem(12)(1:lmem(12))
        else
          write (31,'(2a)') mem(117)(1:lmem(117)),mem(13)(1:lmem(13))
        end if
        write (31,'(a,1p1e22.15)') mem(118)(1:lmem(118)),cubx
        write (31,'(a,1p1e22.15)') mem(119)(1:lmem(119)),cuby
        write (31,'(a,1p1e22.15)') mem(120)(1:lmem(120)),cubz
        if (opt(4).eq.0) then
          write (31,'(2a)') mem(71)(1:lmem(71)),mem(12)(1:lmem(12))
        else
          write (31,'(2a)') mem(71)(1:lmem(71)),mem(13)(1:lmem(13))
        end if
        write (31,'(a,1p1e22.15)') mem(72)(1:lmem(72)),xi
        write (31,'(a,1p1e22.15)') mem(73)(1:lmem(73)),xf
        write (31,'(a,1p1e22.15)') mem(74)(1:lmem(74)),yi
        write (31,'(a,1p1e22.15)') mem(75)(1:lmem(75)),yf
        write (31,'(a,1p1e22.15)') mem(76)(1:lmem(76)),zi
        write (31,'(a,1p1e22.15)') mem(77)(1:lmem(77)),zf
        write (31,'(a,1p1e22.15)') mem(78)(1:lmem(78)),incx
        write (31,'(a,1p1e22.15)') mem(79)(1:lmem(79)),incy
        write (31,'(a,1p1e22.15)') mem(80)(1:lmem(80)),incz
c
c Infrequently-changed parameters
        write (31,'(a)') mem(22)(1:lmem(22))
        write (31,'(a)') mem(30)(1:lmem(30))
        write (31,'(a)') mem(22)(1:lmem(22))
        write (31,'(a)',advance='no') mem(38)(1:lmem(38))
        do k = 1, opt(2)-1
          write (31,'(1x,1p1e22.15)',advance='no')
     %      dens(opt(2)-k+1)
        enddo
        write (31,'(1x,1p1e22.15)',advance='yes') dens(1)
        write (31,*) mem(32)(1:lmem(32)),2.0d0*pi/(omg*3600.0d0)
        write (31,'(a,1p1e22.15)') mem(39)(1:lmem(39)),factor
        write (31,'(a,1p1e22.15)') mem(40)(1:lmem(40)),dsqrt(eixa)
        write (31,'(a,1p1e22.15)') mem(41)(1:lmem(41)),dsqrt(eixb)
        write (31,'(a,1p1e22.15)') mem(42)(1:lmem(42)),dsqrt(eixc)
        write (31,'(a)') mem(22)(1:lmem(22))
        write (31,'(a,1p1e22.15)') mem(43)(1:lmem(43)),gc
        close (31)
c
c Create new version of the restart file
  60    if (idp.eq.1) open (32, file='restart.tmp', status='unknown',
     %    err=60)
  65    if (idp.eq.2) open (32, file=dumpfile(2), status='old', err=65)
        write (32,'(1x,i10,10x,a3)') l,'EOL'
        write (32,*) treal
        write (32,*) tuser
        write (32,*) tsys
        close (32)
      end do
c
      if (l.ne.-1) then
        call mio_elapse (treal2,tuser2,tsys2,83,mem,lmem)
        call mio_elapse (treal,tuser,tsys,82,mem,lmem)
        call mio_remaining (treal,l+1,nograd,mem,lmem)
        call mio_sizeout (opt(1),l,nograd,nflush,outfile,mem,lmem)
        call mio_pid ( )
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READPOLYHEDRON.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read vertices, faces and compute normal faces of a polyhedron.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine readPolyhedron (infile,mem,lmem,nov,nop,xv,yv,zv,noe,k,
     %  norm,w,factor)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(6)
      integer nov,nop,noe,k
      real*8 xv(novmax),yv(novmax),zv(novmax)
      real*8 norm(nopmax,3),w(nopmax),factor
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      logical test
      integer i,j,l
      character*80 filename
      real*8 dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len
      integer lim(2,1000),nsub
      character*15000 string
      character*80 c80
      character*5 c5
      integer lineno
      integer flag
c
c------------------------------------------------------------------------------
c
      do j = 1, 80
        filename(j:j) = ' '
      end do
c
      do j = 1, 2
        filename = infile(j)
        inquire (file=filename, exist=test)
        if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),
     %    lmem(2),' ',1,filename,80)

  33    open  (10+j, file=filename, status='old', err=33)
c
        lineno = 0
        if (j.eq.1) then
c lê os vértices do poliedro
            i = 0
c          do 21 i = 1, nov
  40        lineno = lineno + 1
c            read (10+j,'(a15000)') string
            read (10+j,'(a15000)',end=2107) string
            if (string(1:1).eq.')') goto 40
            call mio_spl (15000,string,nsub,lim)
            if (lim(1,1).eq.-1) goto 40
            i = i + 1
            if (i.gt.novmax) call mio_err (6,mem(1),lmem(1),mem(16),
     %        lmem(16),' ',1,mem(14),lmem(14))
            c80 = string(lim(1,1):lim(2,1))
            read (c80,*,err=661) xv(i)
            c80 = string(lim(1,2):lim(2,2))
            read (c80,*,err=661) yv(i)
            c80 = string(lim(1,3):lim(2,3))
            read (c80,*,err=661) zv(i)
c            read(10+j,*) xv(i),yv(i),zv(i)
            xv(i) = xv(i) * factor
            yv(i) = yv(i) * factor
            zv(i) = zv(i) * factor
c  21      continue
            goto 40
 2107       nov = i
        end if
        if (j.eq.2) then
c lê as faces e o número de vértices por face do poliedro
            flag = 0
            i = 0
c          do 32 i=1,nop
  41        lineno = lineno + 1
c            read (10+j,'(a15000)') string
            read (10+j,'(a15000)',end=2108) string
            if (string(1:1).eq.')') goto 41
            call mio_spl (15000,string,nsub,lim)
            if (lim(1,1).eq.-1) goto 41
c            c80 = string(lim(1,1):lim(2,1))
c            read (c80,*,err=662) noe(i)
            i = i + 1
            if (i.gt.nopmax) call mio_err (6,mem(1),lmem(1),mem(17),
     %        lmem(17),' ',1,mem(14),lmem(14))
            noe(i) = nsub
            do l = 1, noe(i)
c              c80 = string(lim(1,l+1):lim(2,l+1))
              c80 = string(lim(1,l):lim(2,l))
              read (c80,*,err=662) k(i,l)
c              if (k(i,l).eq.0) goto 663
              if (k(i,l).eq.0) flag = 1
            end do
c lê a primeira coluna do arquivo 15 onde é guardado o número inteiro (de até 3 algarismos (i3), por exemplo, 120) de vértices por face no vetor noe(i). OBS: advance='no' faz com que read NÃO avance para a leitura da próxima coluna, isto é, faz com que o cursor fique posicionado após a leitura da primeira coluna
c            read(10+j,'(i3)',advance='no') noe(i)
c lê as colunas restantes das linhas onde são guardados as posições dos vértices das faces na matriz k(i,j)
c            read(10+j,*) (k(i,l),l=1,noe(i))
c  32      continue
            goto 41
 2108       nop = i
c
            if (flag.eq.1) then
              do i = 1, nop
                do l = 1, noe(i)
                  k(i,l) = k(i,l) + 1
                end do
              end do
            end if
        end if
        close (10+j)
      end do
c
c compute face normal and offset w from first 3 vertices
      do i = 1, nop
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        norm(i,1) = nx / len
        norm(i,2) = ny / len
        norm(i,3) = nz / len
        w(i) = - norm(i,1) * xv(k(i,1))
     %         - norm(i,2) * yv(k(i,1))
     %         - norm(i,3) * zv(k(i,1))
      end do
c
c------------------------------------------------------------------------------
c
      return
c
c Error reading from the input file containing integration parameters
 661  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c5,5,
     %  mem(18),lmem(18))
c
 662  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c5,5,
     %  mem(19),lmem(19))
c
 663  write (c5,'(i5)') lineno
      call mio_err (6,mem(1),lmem(1),mem(110),lmem(110),c5,5,
     %  mem(19),lmem(19))
c
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_MASS.FOR    (8 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal face of a polyhedron.
c Also computes the total mass and the total volumen of a polyhedron.
c
c------------------------------------------------------------------------------
c
      subroutine masc_mass (nov,nop,xv,yv,zv,noe,k,d,m,mt,svt)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 m(nocen),d,mt,svt
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
      integer i,j
      real*8 v,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,sv
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,nx,ny,nz,len,dot
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      svt= 0.0
      mt = 0.0
c
      do i = 1, nop
c
c compute normal face
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        nx = nx / len
        ny = ny / len
        nz = nz / len
c dot product
        dot = xv(k(i,1))*nx+yv(k(i,1))*ny+zv(k(i,1))*nz
c
c compute volume of tetrahedron
        sv = 0.0
c
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          v  = (x4-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))+
     %         ((y4-y1)*((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)))+
     %         (z4-z1)*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
c
          v  = sign(v,dot)
          sv = sv + v
        end do
        v    = sv / 6.d0
        m(i) = d * v
        mt   = mt + m(i)
        svt = svt + v
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_CEN.FOR    (8 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the center of mass of each tetrahedron face of a polyhedron.
c
c------------------------------------------------------------------------------
c
      subroutine masc_cen (nov,nop,xv,yv,zv,noe,k,xc,yc,zc)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen)
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
      integer i,j,l
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
      l  = 0
c
      do i = 1, nop
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          l = l + 1
c
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          xc(l) = (x1 + x2 + x3 + x4) * 0.25d0
          yc(l) = (y1 + y2 + y3 + y4) * 0.25d0
          zc(l) = (z1 + z2 + z3 + z4) * 0.25d0
        end do
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c
c------------------------------------------------------------------------------
c
      subroutine masc_layer (nov,nop,xv,yv,zv,noe,k,d,lay,m,mt,vt,
     %  xc,yc,zc,mem,lmem,l0,T1t,T2t,TPt,Jt)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,lay,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),m(nocen),d(nocen),mt,vt
      real*8 T1t(3),T2t(3),TPt(3),Jt(3,3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,l,ilay,novl,nopl
      integer noel,kl
      real*8 del,inc,nx,ny,nz
      real*8 xl(nov),yl(nov),zl(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 fvt,fmt,vol,vol1,vol2
      real*8 J0(3,3)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      del = 1.d0 / lay
      vt = 0.d0
      mt = 0.d0
      l0 = 0
c
      T1t(1) = 0.0
      T1t(2) = 0.0
      T1t(3) = 0.0
      T2t(1) = 0.0
      T2t(2) = 0.0
      T2t(3) = 0.0
      TPt(1) = 0.0
      TPt(2) = 0.0
      TPt(3) = 0.0
c central point (origin)
      xl(1) = 0.0
      yl(1) = 0.0
      zl(1) = 0.0
c
      do i = 1, nop
        fvt = 0.d0
        fmt = 0.d0
c
c compute new vertices
        l = 1
        do ilay = 1, lay
          inc = del * ilay
          do j = 1, noe(i)
            l = l + 1
            xl(l) = xv(k(i,j)) * inc
            yl(l) = yv(k(i,j)) * inc
            zl(l) = zv(k(i,j)) * inc
          enddo
        enddo
        novl = l
c
c compute new faces
        l = 1
c compute top face of pyramid
        do j = 1, noe(i)
          kl(l,j) = j + 1
        enddo
        noel(l) = noe(i)
c compute lateral faces of pyramid
        do j = 1, noe(i)-1
          l = l + 1
          kl(l,1) = 1
          kl(l,2) = j + 2
          kl(l,3) = j + 1
          noel(l) = 3
        enddo
c compute last lateral face of pyramid
        l = l + 1
        kl(l,1) = 1
        kl(l,2) = 2
        kl(l,3) = noe(i) + 1
        noel(l) = 3
c
c compute volumens
        nopl = l
c compute normal faces and w vector
        do j = 1, nopl
          call compnormalface (nov,nop,j,xl,yl,zl,kl,nx,ny,nz)
          norm(j,1) = nx
          norm(j,2) = ny
          norm(j,3) = nz
          wfac(j) = - norm(j,1) * xl(kl(j,1))
     %              - norm(j,2) * yl(kl(j,1))
     %              - norm(j,3) * zl(kl(j,1))
        enddo
c compute volumen of a irregular pyramid (first tetrahedron face)
        call compVolumeIntegrals (nov,nop,nopl,xl,yl,zl,noel,kl,
     %    norm,wfac,T0,T1,T2,TP)
        fvt = fvt + T0
        T1t(1) = T1t(1) + T1(1)
        T1t(2) = T1t(2) + T1(2)
        T1t(3) = T1t(3) + T1(3)
        T2t(1) = T2t(1) + T2(1)
        T2t(2) = T2t(2) + T2(2)
        T2t(3) = T2t(3) + T2(3)
        TPt(1) = TPt(1) + TP(1)
        TPt(2) = TPt(2) + TP(2)
        TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
        l0 = l0 + 1
        m(l0) = d(1) * T0
        fmt = fmt + m(l0)
        call compcenpolyhedron (d(1),m(l0),T0,T1,T2,TP,
     %    xc(l0),yc(l0),zc(l0),J0)
        Jt(1,1) = Jt(1,1) + J0(1,1)
        Jt(1,2) = Jt(1,2) + J0(1,2)
        Jt(1,3) = Jt(1,3) + J0(1,3)
        Jt(2,1) = Jt(2,1) + J0(2,1)
        Jt(2,2) = Jt(2,2) + J0(2,2)
        Jt(2,3) = Jt(2,3) + J0(2,3)
        Jt(3,1) = Jt(3,1) + J0(3,1)
        Jt(3,2) = Jt(3,2) + J0(3,2)
        Jt(3,3) = Jt(3,3) + J0(3,3)
c
c compute faces of pyramidal frustums
        do ilay = 2, lay
          l0 = l0 + 1
          l = 1
c compute bottom face (reverse)
          do j = 1, noe(i)
            kl(l,j) = (ilay-2)*noe(i) + (noe(i)-j) + 2
          enddo
          noel(l) = noe(i)
c compute lateral faces
          do j = 1, noe(i)-1
            l = l + 1
            kl(l,1) = (ilay-2)*noe(i) + j + 1
            kl(l,2) = kl(l,1) + 1
            kl(l,3) = ((ilay-1)*noe(i) + j) + 2
            kl(l,4) = kl(l,3) - 1
            noel(l) = 4
          enddo
c compute last lateral face
          l = l + 1
          kl(l,1) = (ilay-1)*noe(i) + 1
          kl(l,2) = kl(l,1) - (noe(i)-1)
          kl(l,3) = ilay*noe(i) - (noe(i)-1) + 1
          kl(l,4) = ilay*noe(i) + 1
          noel(l) = 4
c compute top face (no reverse)
          l = l + 1
          do j = 1, noe(i)
            kl(l,j) = (ilay-1)*noe(i) + j + 1
          enddo
          noel(l) = noe(i)
c
c compute volumens
          nopl = l
c compute normal faces and w vector
          do j = 1, nopl
            call compnormalface (nov,nop,j,xl,yl,zl,kl,nx,ny,nz)
            norm(j,1) = nx
            norm(j,2) = ny
            norm(j,3) = nz
            wfac(j) = - norm(j,1) * xl(kl(j,1))
     %                - norm(j,2) * yl(kl(j,1))
     %                - norm(j,3) * zl(kl(j,1))
          enddo
c compute volumens of a irregular pyramidal frustums
          call compVolumeIntegrals (nov,nop,nopl,xl,yl,zl,noel,kl,
     %      norm,wfac,T0,T1,T2,TP)
          fvt = fvt + T0
          T1t(1) = T1t(1) + T1(1)
          T1t(2) = T1t(2) + T1(2)
          T1t(3) = T1t(3) + T1(3)
          T2t(1) = T2t(1) + T2(1)
          T2t(2) = T2t(2) + T2(2)
          T2t(3) = T2t(3) + T2(3)
          TPt(1) = TPt(1) + TP(1)
          TPt(2) = TPt(2) + TP(2)
          TPt(3) = TPt(3) + TP(3)
c          write(*,*) T0,d(lay)*T0
c          call compvoltetrahedron (nov,nop,5,xl,yl,zl,noel,kl,vol)
c          vol2 = vol
c          call compvoltetrahedron (nov,nop,1,xl,yl,zl,noel,kl,vol)
c          vol1 = -vol
c          write(*,*) vol1,vol2,vol2-vol1
c          stop
c compute pieces of masses
          m(l0) = d(ilay) * T0
          fmt = fmt + m(l0)
c compute center of masses of irregular pyramidal frustums
          call compcenpolyhedron (d(ilay),m(l0),T0,T1,
     %      T2,TP,xc(l0),yc(l0),zc(l0),J0)
          Jt(1,1) = Jt(1,1) + J0(1,1)
          Jt(1,2) = Jt(1,2) + J0(1,2)
          Jt(1,3) = Jt(1,3) + J0(1,3)
          Jt(2,1) = Jt(2,1) + J0(2,1)
          Jt(2,2) = Jt(2,2) + J0(2,2)
          Jt(2,3) = Jt(2,3) + J0(2,3)
          Jt(3,1) = Jt(3,1) + J0(3,1)
          Jt(3,2) = Jt(3,2) + J0(3,2)
          Jt(3,3) = Jt(3,3) + J0(3,3)
        enddo
        vt = vt + fvt
        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(15),
     %    lmem(15),' ',1,mem(14),lmem(14))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPNORMALFACE.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Compute normal face.
c
c------------------------------------------------------------------------------
c
      subroutine compnormalface (nov,nop,i,xv,yv,zv,k,nx,ny,nz)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,i,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 nx,ny,nz
      dimension k(nopmax,noed)
c
c Local
c      integer i
      real*8 dx1,dy1,dz1,dx2,dy2,dz2,len
c
c------------------------------------------------------------------------------
c
      nx = 0.d0
      ny = 0.d0
      nz = 0.d0
c
c      do i = 1, nop
        dx1 = xv(k(i,2)) - xv(k(i,1))
        dy1 = yv(k(i,2)) - yv(k(i,1))
        dz1 = zv(k(i,2)) - zv(k(i,1))
        dx2 = xv(k(i,3)) - xv(k(i,2))
        dy2 = yv(k(i,3)) - yv(k(i,2))
        dz2 = zv(k(i,3)) - zv(k(i,2))
        nx = dy1 * dz2 - dy2 * dz1
        ny = dz1 * dx2 - dz2 * dx1
        nz = dx1 * dy2 - dx2 * dy1
        len = dsqrt(nx * nx + ny * ny + nz * nz)
        nx = nx / len
        ny = ny / len
        nz = nz / len
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPVOLTETRAHEDRON.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Compute volumen of a pyramid.
c
c------------------------------------------------------------------------------
c
      subroutine compvoltetrahedron (nov,nop,i,xv,yv,zv,noe,k,vol)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,i,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 vol
      dimension noe(nop)
      dimension k(nop,noed)
c
c Local
c      integer i,j
      integer j
      real*8 x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
      real*8 v,sv,nx,ny,nz,dot
c
c------------------------------------------------------------------------------
c
      x1 = 0.0
      y1 = 0.0
      z1 = 0.0
c
c      do i = 1, nop
c
c compute normal face
        call compnormalface (nov,nop,i,xv,yv,zv,k,nx,ny,nz)
        dot = xv(k(i,1))*nx+yv(k(i,1))*ny+zv(k(i,1))*nz
c
c compute volume of tetrahedron
        sv = 0.0
c
        x2 = xv(k(i,1))
        y2 = yv(k(i,1))
        z2 = zv(k(i,1))
c
        j = 2
        do while (j.le.noe(i))
          x3 = xv(k(i,j))
          y3 = yv(k(i,j))
          z3 = zv(k(i,j))
c
          j  =  j  +  1
          x4 = xv(k(i,j))
          y4 = yv(k(i,j))
          z4 = zv(k(i,j))
c
          v  = (x4-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))+
     %         ((y4-y1)*((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)))+
     %         (z4-z1)*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
c
          v  = sign(v,dot)
          sv = sv + v
        end do
        vol = sv / 6.d0
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPCENPOLYHEDRON.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes centroid of a polyhedron. Also computes intertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine compcenpolyhedron (d,m,T0,T1,T2,TP,xc,yc,zc,J)
c
      implicit none
c
c Input/Output
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 xc,yc,zc,d,m
      real*8 J(3,3)
c Local
      real*8 r(3)
c
c------------------------------------------------------------------------------
c
c as coordenadas do centro de massa R_cm de um corpo são dadas pela integral volumétrica R_cm=1/M*int(r)dm, onde M é a massa total do corpo, r é o vetor posição e dm=ro*dv é o elemento de massa. Se a massa está distribuída de forma homogênea, a densidade será constante, assim, fazendo uso das relações dm=ro*dv e int(dm)=ro*int(dv) => M=ro*V as coordenadas do centro de massa são dadas por R_cm=ro*int(r)dv/(ro*V) => R_cm=int(r)dv/V. OBS: quando utilizamos o volume o centro de massa independe da densidade. Explicitando as componentes do vetor r ficamos com X_cm=int(x)dxdydz/V, Y_cm=int(y)dxdydz/V e Z_cm=int(z)dxdydz/V
c compute center of mass
c coordenada X_cm=int(x)dxdydz/V do centro de massa
      r(1) = T1(1) / T0
      xc   = r(1)
c coordenada Y_cm=int(y)dxdydz/V do centro de massa
      r(2) = T1(2) / T0
      yc   = r(2)
c coordenada Z_cm=int(z)dxdydz/V do centro de massa
      r(3) = T1(3) / T0
      zc   = r(3)
c
c compute inertia tensor
c de modo geral o elemento I_ij do tensor de inércia é dado pela integral volumétrica I_ij=int(r^2*dk_ij-x_i*x_j)ro*dv; onde r=x^2+y^2+z^2, dk_ij é o delta de Kronecker (1 se i=j, 0 se i!=j), ro é densidade volumétrica, dv=dxdydz é o elemento de volume e x_i com i=1,2,3 são as coordenadas dos vértices, ou seja, x_1=x, x_2=y e x_3=z, o mesmo vale para x_j
c momento de inércia I_xx=ro*int(y^2+z^2)dxdydz
      J(1,1) = d * (T2(2) + T2(3))
c momento de inércia I_yy=ro*int(z^2+x^2)dxdydz
      J(2,2) = d * (T2(3) + T2(1))
c momento de inércia I_zz=ro*int(x^2+y^2)dxdydz
      J(3,3) = d * (T2(1) + T2(2))
c momento de inércia I_xy=-ro*int(xy)dxdydz=I_yx
      J(1,2) =-d * TP(1)
      J(2,1) = J(1,2)
c momento de inércia I_yz=-ro*int(yz)dxdydz=I_zy
      J(2,3) =-d * TP(2)
      J(3,2) = J(2,3)
c momento de inércia I_zx=-ro*int(zx)dxdydz=I_xz
      J(3,1) =-d * TP(3)
      J(1,3) = J(3,1)
c teorema do eixo paralelo (Huygens-Steiner) I_ij=I_ij_cm+M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm) => I_ij_cm=I_ij-M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm)
c translate inertia tensor to center of mass
c      J(1,1) = J(1,1) - m * (r(2)*r(2) + r(3)*r(3))
c      J(2,2) = J(2,2) - m * (r(3)*r(3) + r(1)*r(1))
c      J(3,3) = J(3,3) - m * (r(1)*r(1) + r(2)*r(2))
c
c      J(1,2) = J(1,2) + m * r(1) * r(2)
c      J(2,1) = J(1,2)
c      J(2,3) = J(2,3) + m * r(2) * r(3)
c      J(3,2) = J(2,3)
c      J(3,1) = J(3,1) + m * r(3) * r(1)
c      J(1,3) = J(3,1)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPVOLUMEINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes volume integrals of each face. In general computes volumen of any
c polyhedron. Also computes intertia's axes.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compVolumeIntegrals (nov,nop,nopl,xv,yv,zv,noe,k,norm,
     %  wfac,T0,T1,T2,TP)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,nopl,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
      integer A,B,C
      real*8 nx, ny, nz, TEST, w
      real*8 Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca
c
c------------------------------------------------------------------------------
c
      T0    = 0.d0
      T1(1) = 0.d0
      T1(2) = 0.d0
      T1(3) = 0.d0
      T2(1) = 0.d0
      T2(2) = 0.d0
      T2(3) = 0.d0
      TP(1) = 0.d0
      TP(2) = 0.d0
      TP(3) = 0.d0
c
      do i = 1, nopl
        nx = dabs(norm(i,1))
        ny = dabs(norm(i,2))
        nz = dabs(norm(i,3))
        if (nx.gt.ny.and.nx.gt.nz) then
          C = 0
        else if (ny.gt.nz) then
          C = 1
        else
          C = 2
        end if
c
        A = mod(C + 1, 3)
        B = mod(A + 1, 3)
        A = A + 1
        B = B + 1
        C = C + 1
c
        Fa   = 0.d0
        Fb   = 0.d0
        Fc   = 0.d0
        Faa  = 0.d0
        Fbb  = 0.d0
        Fcc  = 0.d0
        Faaa = 0.d0
        Fbbb = 0.d0
        Fccc = 0.d0
        Faab = 0.d0
        Fbbc = 0.d0
        Fcca = 0.d0
c
        w = wfac(i)
c
        call compFaceIntegrals(nov,nop,i,xv,yv,zv,noe,k,norm,w,A,B,C,
     %    Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca)
c
        if (A.eq.1) then
          TEST = Fa
        else if (B.eq.1) then
          TEST = Fb
        else
          TEST = Fc
        end if
c
        T0 = T0 + norm(i,1) * TEST
c
        T1(A) = T1(A) + norm(i,A) * Faa
        T1(B) = T1(B) + norm(i,B) * Fbb
        T1(C) = T1(C) + norm(i,C) * Fcc
        T2(A) = T2(A) + norm(i,A) * Faaa
        T2(B) = T2(B) + norm(i,B) * Fbbb
        T2(C) = T2(C) + norm(i,C) * Fccc
        TP(A) = TP(A) + norm(i,A) * Faab
        TP(B) = TP(B) + norm(i,B) * Fbbc
        TP(C) = TP(C) + norm(i,C) * Fcca
c
      end do
c
      T1(1) = T1(1) / 2.d0
      T1(2) = T1(2) / 2.d0
      T1(3) = T1(3) / 2.d0
      T2(1) = T2(1) / 3.d0
      T2(2) = T2(2) / 3.d0
      T2(3) = T2(3) / 3.d0
      TP(1) = TP(1) / 2.d0
      TP(2) = TP(2) / 2.d0
      TP(3) = TP(3) / 2.d0
c
      do i = 1, 3
        if (dabs(T1(i)).le.TINY) T1(i) = 0.d0 
      end do
      do i = 1, 3
        if (dabs(T2(i)).le.TINY) T2(i) = 0.d0 
      end do
      do i = 1, 3
        if (dabs(TP(i)).le.TINY) TP(i) = 0.d0 
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPFACEINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes face integrals of each face.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compFaceIntegrals (nov,nop,nf,xv,yv,zv,noe,k,n,w,A,B,C,
     %  Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,nf,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 n(nopmax,3),w
      integer A,B,C
      real*8 Fa,Fb,Fc,Faa,Fbb,Fcc,Faaa,Fbbb,Fccc,Faab,Fbbc,Fcca
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      real*8 k1, k2b, k3, k4
      real*8 P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb
      real*8 sqr,cube
c
c------------------------------------------------------------------------------
c
      call compProjectionIntegrals(nov,nop,nf,xv,yv,zv,noe,k,A,B,
     %  P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)
c
      k1 = 1.d0 / n(nf,C)
      k2b = k1 * k1
      k3 = k2b * k1
      k4 = k3 * k1
c
      Fa = k1 * Pa
      Fb = k1 * Pb
      Fc = -k2b * (n(nf,A)*Pa + n(nf,B)*Pb + w*P1)
c
      Faa = k1 * Paa
      Fbb = k1 * Pbb
      Fcc = k3 * (sqr(n(nf,A))*Paa + 2.d0*n(nf,A)*n(nf,B)*Pab +
     %  sqr(n(nf,B))*Pbb + w*(2.d0*(n(nf,A)*Pa + n(nf,B)*Pb) + w*P1))
c
      Faaa = k1 * Paaa
      Fbbb = k1 * Pbbb
      Fccc = -k4 * (cube(n(nf,A))*Paaa +
     %  3.d0*sqr(n(nf,A))*n(nf,B)*Paab +
     %  3.d0*n(nf,A)*sqr(n(nf,B))*Pabb + cube(n(nf,B))*Pbbb +
     %  3.d0*w*(sqr(n(nf,A))*Paa + 2.d0*n(nf,A)*n(nf,B)*Pab +
     %  sqr(n(nf,B))*Pbb) + w*w*(3.d0*(n(nf,A)*Pa + n(nf,B)*Pb) + w*P1))
c
      Faab = k1 * Paab
      Fbbc = -k2b * (n(nf,A)*Pabb + n(nf,B)*Pbbb + w*Pbb)
      Fcca = k3 * (sqr(n(nf,A))*Paaa + 2.d0*n(nf,A)*n(nf,B)*Paab +
     %  sqr(n(nf,B))*Pabb + w*(2.d0*(n(nf,A)*Paa + n(nf,B)*Pab) + w*Pa))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      COMPPROJECTIONINTEGRALS.FOR    (15 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes projection integrals of each face.
c
c Adapted by A. Amarante (Fortran 77)
c Brian Mirtich, "Fast and Accurate Computation of Polyhedral Mass Properties,
c " journal of graphics tools, volume 1, number 1, 1996.
c
c------------------------------------------------------------------------------
c
      subroutine compProjectionIntegrals(nov,nop,nf,xv,yv,zv,noe,k,A,B,
     %  P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,nf,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      integer A,B
      real*8 P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      real*8 a0, a1, da
      real*8 b0, b1, db
      real*8 a0_2, a0_3, a0_4, b0_2, b0_3, b0_4
      real*8 a1_2, a1_3, b1_2, b1_3
      real*8 C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb
      real*8 Cab, Kab, Caab, Kaab, Cabb, Kabb
      integer i
c
c------------------------------------------------------------------------------
c
      P1   = 0.d0
      Pa   = 0.d0
      Pb   = 0.d0
      Paa  = 0.d0
      Pab  = 0.d0
      Pbb  = 0.d0
      Paaa = 0.d0
      Paab = 0.d0
      Pabb = 0.d0
      Pbbb = 0.d0
c
      do i = 0, noe(nf)-1
c
        if (A.eq.1) then
          a0 = xv(k(nf,i+1))
          a1 = xv(k(nf,mod((i+1),noe(nf))+1))
        else if (A.eq.2) then
          a0 = yv(k(nf,i+1))
          a1 = yv(k(nf,mod((i+1),noe(nf))+1))
        else if (A.eq.3) then
          a0 = zv(k(nf,i+1))
          a1 = zv(k(nf,mod((i+1),noe(nf))+1))
        end if
c
        if (B.eq.1) then
          b0 = xv(k(nf,i+1))
          b1 = xv(k(nf,mod((i+1),noe(nf))+1))
        else if (B.eq.2) then
          b0 = yv(k(nf,i+1))
          b1 = yv(k(nf,mod((i+1),noe(nf))+1))
        else if (B.eq.3) then
          b0 = zv(k(nf,i+1))
          b1 = zv(k(nf,mod((i+1),noe(nf))+1))
        end if
c
        da = a1 - a0
        db = b1 - b0
c
        a0_2 = a0 * a0
        a0_3 = a0_2 * a0
        a0_4 = a0_3 * a0
c
        b0_2 = b0 * b0
        b0_3 = b0_2 * b0
        b0_4 = b0_3 * b0
c
        a1_2 = a1 * a1
        a1_3 = a1_2 * a1
c
        b1_2 = b1 * b1
        b1_3 = b1_2 * b1
c
        C1 = a1 + a0
c
        Ca = a1*C1 + a0_2
        Caa = a1*Ca + a0_3
        Caaa = a1*Caa + a0_4
c
        Cb = b1*(b1 + b0) + b0_2
        Cbb = b1*Cb + b0_3
        Cbbb = b1*Cbb + b0_4
c
        Cab = 3.d0*a1_2 + 2.d0*a1*a0 + a0_2
        Kab = a1_2 + 2.d0*a1*a0 + 3.d0*a0_2
c
        Caab = a0*Cab + 4.d0*a1_3
        Kaab = a1*Kab + 4.d0*a0_3
c
        Cabb = 4.d0*b1_3 + 3.d0*b1_2*b0 + 2.d0*b1*b0_2 + b0_3
c
        Kabb = b1_3 + 2.d0*b1_2*b0 + 3.d0*b1*b0_2 + 4.d0*b0_3
c
        P1   = P1   + db*C1
        Pa   = Pa   + db*Ca
        Paa  = Paa  + db*Caa
        Paaa = Paaa + db*Caaa
        Pb   = Pb   + da*Cb
        Pbb  = Pbb  + da*Cbb
        Pbbb = Pbbb + da*Cbbb
        Pab  = Pab  + db*(b1*Cab  + b0*Kab)
        Paab = Paab + db*(b1*Caab + b0*Kaab)
        Pabb = Pabb + da*(a1*Cabb + a0*Kabb)
      end do
c
      P1   = P1 / 2.0d0
      Pa   = Pa / 6.0d0
      Paa  = Paa / 12.0d0
      Paaa = Paaa / 20.0d0
      Pb   = Pb / (-6.0d0)
      Pbb  = Pbb / (-12.0d0)
      Pbbb = Pbbb / (-20.0d0)
      Pab  = Pab / 24.0d0
      Paab = Paab / 60.0d0
      Pabb = Pabb / (-60.0d0)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      SQR.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c A squared function.
c
c------------------------------------------------------------------------------
c
      function sqr (x)
c
      implicit none
c
c Input/Output
      real*8 x,sqr
c
c------------------------------------------------------------------------------
c
      sqr = x * x
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CUBE.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c A cubic function.
c
c------------------------------------------------------------------------------
c
      function cube (x)
c
      implicit none
c
c Input/Output
      real*8 x,cube
c
c------------------------------------------------------------------------------
c
      cube = x * x * x
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_TRANS.FOR    (18 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translate inertia tensor to the center of mass.
c
c------------------------------------------------------------------------------
c
      subroutine masc_trans (J0,xc,yc,zc,m,J)
c
      implicit none
c
c Input/Output
      real*8 J0(3,3),J(3,3)
      real*8 xc,yc,zc,m
c
c------------------------------------------------------------------------------
c
c teorema do eixo paralelo (Huygens-Steiner) I_ij=I_ij_cm+M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm) => I_ij_cm=I_ij-M*((R_cm)^2*dk_ij-X_i_cm*X_j_cm)
c translate inertia tensor to center of mass
      J(1,1) = J0(1,1) - m * (yc*yc + zc*zc)
      J(2,2) = J0(2,2) - m * (zc*zc + xc*xc)
      J(3,3) = J0(3,3) - m * (xc*xc + yc*yc)
c
      J(1,2) = J0(1,2) + m * xc * yc
      J(2,1) = J(1,2)
      J(2,3) = J0(2,3) + m * yc * zc
      J(3,2) = J(2,3)
      J(3,1) = J0(3,1) + m * zc * xc
      J(1,3) = J(3,1)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      EIGENVECTORS.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes normalized eigenvectors for inertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine EigenVectors (unit,mem,lmem,J,eigVecs,eigVals,oldflag)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 J(3,3),eigVecs(3,3),eigVals(3)
      logical oldflag
c
c Local
      integer iEig
      real*8 Ixx,Iyy,Izz,Ixy,Iyz,Izx,Iyx,Izy,Ixz
      real*8 CC3,CC2,CC1,CC0
      real*8 CubicCoeffs(4)
      real*8 sqr
      real*8 eigeval(3)
      real*8 JI(3,3)
      integer i0, j0, i1, j1, iDrop
      integer iMaxDet, i1MaxDet, jMaxDet, j1MaxDet
      real*8 det, maxDet
      real*8 result
      real*8 x, y, z, lat, wlon
      real*8 raio
c
c------------------------------------------------------------------------------
c
      Ixx = J(1,1)
      Iyy = J(2,2)
      Izz = J(3,3)
      Ixy = J(1,2)
      Iyz = J(2,3)
      Izx = J(3,1)
      Iyx = Ixy
      Izy = Iyz
      Ixz = Izx
c
      CC3 = -1.0
      CC2 = Ixx + Iyy + Izz
      CC1 = sqr(Ixy) + sqr(Iyz) + sqr(Izx)
     %  - (Ixx*Iyy + Iyy*Izz + Izz*Ixx)
      CC0 = (2.0*Ixy*Iyz*Izx) + (Ixx*Iyy*Izz)
     %  - (Izz*sqr(Ixy) + Ixx*sqr(Iyz) + Iyy*sqr(Izx))
c      write(*,'(4(1x,1p,e35.25))') CC3,CC2,CC1,CC0
c      write(*,'(1x,1p,e35.25)') J(1,1)
c      write(*,'(1x,1p,e35.25)') J(1,2)
c      write(*,'(1x,1p,e35.25)') J(1,3)
c      write(*,'(1x,1p,e35.25)') J(2,1)
c      write(*,'(1x,1p,e35.25)') J(2,2)
c      write(*,'(1x,1p,e35.25)') J(2,3)
c      write(*,'(1x,1p,e35.25)') J(3,1)
c      write(*,'(1x,1p,e35.25)') J(3,2)
c      write(*,'(1x,1p,e35.25)') J(3,3)
c
      CubicCoeffs(4) = CC3
      CubicCoeffs(3) = CC2
      CubicCoeffs(2) = CC1
      CubicCoeffs(1) = CC0
c
      do iEig = 1, 3
        eigVals(iEig) = J(iEig,iEig)
      end do
c
      call CubicRoots(unit,mem,lmem,CubicCoeffs,eigVals)
c
      call evalcubic (eigVals(1),CubicCoeffs,eigeval(1))
      call evalcubic (eigVals(2),CubicCoeffs,eigeval(2))
      call evalcubic (eigVals(3),CubicCoeffs,eigeval(3))
c
      if (.not.oldflag) then
        write(unit,'(1x,a)') mem(59)(1:lmem(59))
        do j0 = 1, 3
          write(unit,'(4x,1p,e35.25,1x,a,e35.25,0p,a)')
     %      eigVals(j0),mem(53)(1:lmem(53)),eigeval(j0),
     %      mem(52)(1:lmem(52))
        enddo
      endif
c
c------------------------------------------------------------------------------
c
      do iEig = 1, 3
c
c        write(*,*) 'Matrix for eigenvalue ',iEig,' = ',eigVals(iEig)
c
        do j0 = 1, 3
          do i0 = 1, 3
            JI(j0,i0) = J(j0,i0)
          enddo
          JI(j0,j0) = JI(j0,j0) - eigVals(iEig)
c          do i0 = 1, 3
c            write(*,'(1x,1p,e19.12,0p)') JI(j0,i0)
c          enddo
c          write(*,'(/)')
        enddo
c
        maxDet = 0.0
        do i0 = 1, 2
          do j0 = 1, 2
            do i1 = i0+1, 3
              do j1 = j0+1, 3
                det = dabs(JI(j0,i1)*JI(j1,i0)-JI(j0,i0)*JI(j1,i1))
                if (det.gt.maxDet) then
                  maxDet = det
                  iMaxDet = i0
                  i1MaxDet = i1
                  jMaxDet = j0
                  j1MaxDet = j1
                end if
              enddo
            enddo
          enddo
        enddo
c
        if (maxDet.eq.0.d0) then
          call mio_err (6,mem(1),lmem(1),mem(67),lmem(67),
     %    ' ',1,' ',1)
        end if
c
        i0 = iMaxDet
        i1 = i1MaxDet
        j0 = jMaxDet
        j1 = j1MaxDet
c
        iDrop = 6 - (i0 + i1)
        det = JI(j0,i0)*JI(j1,i1) - JI(j0,i1)*JI(j1,i0)
        eigVecs(iEig,i1) = -(JI(j0,i0)*JI(j1,iDrop) -
     %    JI(j0,iDrop)*JI(j1,i0))/det
        eigVecs(iEig,i0) =  (JI(j0,i1)*JI(j1,iDrop) -
     %    JI(j0,iDrop)*JI(j1,i1))/det
c
        eigVecs(iEig,iDrop) = 1.0 / dsqrt( 1.0 +
     %    sqr(eigVecs(iEig,i0)) + sqr(eigVecs(iEig,i1)) )
        eigVecs(iEig,i1) = eigVecs(iEig,i1) *
     %    eigVecs(iEig,iDrop)
        eigVecs(iEig,i0) = eigVecs(iEig,i0) *
     %    eigVecs(iEig,iDrop)
c
        if (.not.oldflag) then
          write(unit,'(1x,a9,i1,a8)',advance='no') 'Eigenvec ',iEig,
     %      ' sol''n: '
        endif
c
        do j0 = 1, 3
          result = 0.0
          do i0 = 1, 3
            result = result + JI(j0,i0) * eigVecs(iEig,i0)
          enddo
          if (.not.oldflag) then
            write(unit,'(1x,1p,e19.12)',advance='no') result
          endif
        enddo
        if (.not.oldflag) then
          write(unit,'(1x)',advance='yes')
        endif
      enddo
c
c------------------------------------------------------------------------------
c
c os autovetores já estão normalizados |1|
      if (.not.oldflag) then
        write(unit,'(a)') mem(60)(1:lmem(60))
      endif
      do iEig = 1, 3
        if (.not.oldflag) then
          write(unit,'(3(1x,1p,e22.15))') eigVecs(1,iEig),
     %      eigVecs(2,iEig), eigVecs(3,iEig)
        endif
      enddo
c
      if (.not.oldflag) then
        write(unit,'(a)') mem(68)(1:lmem(68))
      endif
      do iEig = 1, 3
        x = eigVecs(iEig,1)
        y = eigVecs(iEig,2)
        z = eigVecs(iEig,3)
c
c        if (z.lt.-1.d0) then
c          lat = -999.0
c        else if (z.gt.1.d0) then
c          lat = 999.0
c        else
c          lat = asin(z)
c        end if
c        lat = lat * 1.d0 / DR
c
c        if (x.eq.0.d0.and.y.eq.0.d0) then
c          wlon = 0.0
c        else if (x.eq.0.d0) then
c          if (y.lt.0.d0) then
c            wlon =  90.0
c          else
c            wlon = -90.0
c          end if
c        else
c          wlon = atan2( -y, x)
c        end if
c        wlon = wlon * 1.d0 / DR
c        if (wlon.lt.0.d0) wlon = wlon + 360.0
c
        raio= dsqrt(x*x+y*y+z*z)
        lat = asin(z/raio)
        lat = lat * 1.d0 / DR
        wlon= atan2( y, x)
        wlon = wlon * 1.d0 / DR
        do while (wlon.lt.0.d0)
          wlon = wlon + 360.0
        enddo

c
        if (.not.oldflag) then
          write(unit,'(2(2x,f15.6),a)') lat,wlon,mem(69)(1:lmem(69))
        endif
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CUBICROOTS.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes normalized eigenvectors for inertia tensor.
c
c------------------------------------------------------------------------------
c
      subroutine cubicroots (unit,mem,lmem,CubicCoeffs,Roots)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 Roots(3),CubicCoeffs(4)
c
c Local
      real*8 inflectx,determx,loxmincub,hixmaxcub
      real*8 sqr
      real*8 arg1,arg2
c
c------------------------------------------------------------------------------
c
      inflectx = CubicCoeffs(3) / 3.0
      determx = sqr(inflectx) + (CubicCoeffs(2) / 3.0)
c
      if (determx < 0.0) then
        call printfcc(unit,mem,lmem,CubicCoeffs)
        call mio_err (6,mem(1),lmem(1),mem(62),lmem(62),
     %  ' ',1,' ',1)
      end if
c
      determx = dsqrt(determx)
c
      if (determx.eq.0.0) then
        Roots(1) = inflectx
        Roots(2) = inflectx
        Roots(3) = inflectx
        return
      end if
c
      loxmincub = inflectx - determx
      hixmaxcub = inflectx + determx
c
c      write(*,'(1x,a,e35.25)') 'inflectx  = ',inflectx
c      write(*,'(1x,a,e35.25)') 'determx   = ',determx
c      write(*,'(1x,a,e35.25)') 'loxmincub = ',loxmincub
c      write(*,'(1x,a,e35.25)') 'hixmaxcub = ',hixmaxcub
c      write(*,'(1x,a,4(e35.25),/)') 'cubcoeffs = ',
c     %  CubicCoeffs(4),CubicCoeffs(3),CubicCoeffs(2),CubicCoeffs(1)
c
      arg1 = loxmincub
      arg2 = loxmincub-2.0*determx
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(1))
      arg1 = loxmincub
      arg2 = hixmaxcub
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(2))
      arg1 = hixmaxcub+2.0*determx
      arg2 = hixmaxcub
      call SearchCubicRoot(unit,mem,lmem,CubicCoeffs,arg1,arg2,Roots(3))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      PRINTFCC.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print cubic coefficients.
c
c------------------------------------------------------------------------------
c
      subroutine printfcc (unit,mem,lmem,CubicCoeffs)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 CubicCoeffs(4)
c
c------------------------------------------------------------------------------
c
c      write(6,'(1x,a,4(/,2x,e35.25))') mem(61)(1:lmem(61)),
c     %  CubicCoeffs(1),CubicCoeffs(2),CubicCoeffs(3),CubicCoeffs(4)
      write(unit,'(1x,a,4(/,2x,e35.25))') mem(61)(1:lmem(61)),
     %  CubicCoeffs(1),CubicCoeffs(2),CubicCoeffs(3),CubicCoeffs(4)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      SEARCHCUBICROOT.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Search for cubic roots.
c
c------------------------------------------------------------------------------
c
      subroutine SearchCubicRoot(unit,mem,lmem,CubicCoeffs,lox,hix,root)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer unit,lmem(NMESS)
      character*80 mem(NMESS)
      real*8 CubicCoeffs(4),lox,hix,root
c
c Local
      integer iter
      real*8 lo,hi,mid,midx
c
c------------------------------------------------------------------------------
c
      midx = (lox + hix) / 2.0
      iter = 0
c
      call evalcubic (lox,CubicCoeffs,lo)
      if (lo.eq.0.d0) then
        root = lox
        return
      end if
      call evalcubic (hix,CubicCoeffs,hi)
      if (hi.eq.0.d0) then
        root = hix
        return
      end if
c
      if (lo.gt.0.d0.or.hi.lt.0.d0) then
        write(6,'(1x,2(a,e35.25),/,1x,2(a,e35.25),/)')
     %    mem(63)(1:lmem(63)),lox,',',lo,
     %    mem(64)(1:lmem(64)),hix,',',hi
        call printfcc(unit,mem,lmem,CubicCoeffs)
        call mio_err (6,mem(1),lmem(1),mem(65),lmem(65),
     %  ' ',1,' ',1)
        root = midx
        return
      end if
c
      do while (lox.ne.midx.and.midx.ne.hix)
        iter = iter + 1
        if (iter.gt.1000) then
          write(6,'(1x,2(a,e35.25),/,1x,2(a,e35.25),/)')
     %      mem(63)(1:lmem(63)),lox,',',lo,
     %      mem(64)(1:lmem(64)),hix,',',hi
          call printfcc(unit,mem,lmem,CubicCoeffs)
          call mio_err (6,mem(1),lmem(1),mem(66),lmem(66),
     %    ' ',1,' ',1)
          root = midx
          return
        end if
c
        call evalcubic (midx,CubicCoeffs,mid)
        if (mid.eq.0.d0) exit
        if (mid.gt.0.d0) then
          hix = midx
        else
          lox = midx
        end if
        midx = (lox + hix) / 2.0
      end do
c
      root = midx
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      EVALCUBIC.FOR    (16 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes a cubic fuction.
c
c------------------------------------------------------------------------------
c
      subroutine evalcubic (X,CubicCoeffs,r)
c
      implicit none
c
c Input/Output
      real*8 X,CubicCoeffs(4),r
c
c------------------------------------------------------------------------------
c
      r = CubicCoeffs(1)+(X)*(CubicCoeffs(2)+(X)*
     %  (CubicCoeffs(3)+(X)*CubicCoeffs(4)))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT.FOR    (19 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translates vertices to the center of mass and rotates to the main inertia axes.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot (nov,xv,yv,zv,xc,yc,zc,T)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc,yc,zc,T(3,3)
c
c Local
      integer i,j,k
      real*8 x(2,3),xcm(3),TT(3,3)
c
c------------------------------------------------------------------------------
c
      xcm(1) = xc
      xcm(2) = yc
      xcm(3) = zc
c
c computes transposed matrix
c      do i = 1, 3
c        do j = 1, 3
c          TT(i,j) = T(j,i)
c        enddo
c      enddo

c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
        do j = 1, 3
          x(1,j) = x(1,j) - xcm(j)
        enddo
c
c rotate
        do j = 1, 3
          x(2,j) = 0.d0
          do k = 1, 3
c            x(2,j) = x(2,j) + TT(j,k) * x(1,k)
            x(2,j) = x(2,j) + T(j,k) * x(1,k)
          enddo
        enddo
        xv(i) = x(2,1)
        yv(i) = x(2,2)
        zv(i) = x(2,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      POLYHEDRON.FOR    (FEG   7 March 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c------------------------------------------------------------------------------
c
      subroutine polyhedron (nov,nop,x2,y2,z2,noe,k,d,xp,yp,zp,
     %  vx,vy,vz,v,vij,gc)
c
c------------------------------------------------------------------------------
c
c******************************************************************************
c CHANGES IN VERSION 3 (23 July 2013)
c******************************************************************************
c
c Author:
c
c Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c-------------------------------------------------------------------------------
c
c Copyright notice for SEG distribution:
c
c Copyright (c) 2012 by the Society of Exploration Geophysicists.
c For more information, go to http://software.seg.org/2012/0001 .
c You must read and accept usage terms at:
c http://software.seg.org/disclaimer.txt before use.
c
c-------------------------------------------------------------------------------


c	PROGRAMM "polyhedron.f"
c
c	D.Tsoulis                        Thessaloniki, June 2010
c	
c
c	This programm computes the potential, its first
c	and second derivatives of a homogenous polyhedron
c	according to Petrovic (J of G, 1996). The triple
c	integrals of V, Vi and Vij (i,j=1..3) are transformed
c	twice by means of the divergence theorem of Gauss. The
c	transission from volume integrals into line integrals
c	is accomplished in two steps as follows:
c
c	                      GAUSS
c	1. Volume Integral  ----------> Surface Integral
c	2. Surface Integral ----------> Line Integral
c
c
c	Literature:
c	-----------
c	1. Petrovic S. (1996): Determination of the potential of
c	     homogeneous polyhedral bodies using line integrals, 
c	     Journal of Geodesy 71, 44 - 52.
c	2. Werner R.A. and D.J. Scheeres (1997): Exterior gravitation
c	     of a polyhedron derived and compared with harmonic and
c	     mascon gravitation representations of asteroid 4769
c	     Castalia, Celestial Mechanics and Dynamical Astronomy 65,
c	     313 - 344.
c
c	One has to calculate:
c	--------------------
c
c	1. h(P)       : The distance of point P from the plane p (PP')
c	2. s(P)       : -1, 1
c	3. h(PQ)      : The distance of P' from line segment pq (P'P'')
c 	4. s(PQ)      : -1, 1
c	5. cos(Np,ei) : Np is the plane normal to plane p
c	6. cos(npq,ei): npq is the line normal to line segment pq
c	7. LNpq, ANpq : (Transcendental Functions of the
c	                coordinates of the body's vertices)
c
c
c	ALL of the above quantities are computed automatically
c	using only the following three files:
c	1. file 'xyz' which contains the coordinates x(i), y(i) and
c	   z(i) of the i vertices of the polyhedron
c	2. file 'dat' is a line whith (nop) elements, which declare
c	   the number of edges building each of the (nop) planes,
c	   in the same order that they appear in file 'topology'
c	3. file 'topology' where the number of vertices belonging
c 	   at each facet is given, in such an order that the plane
c	   normal is always pointing outside the body
c
c	EXAMPLES FOR THESE FILES
c	------------------------
c
c	file 'xyz':
c
c	  10  10  20
c       30  10  20
c       30  20  20
c       10  20  20
c       10  10  10
c	  30  10  10
c	  30  20  10
c	  10  20  10
c
c	file 'topology':
c
c	  1 2 4 
c	  2 3 4
c	  6 5 8 7
c	  8 5 1 4
c	  6 7 3 2
c	  5 6 2 1
c	  7 8 4 3
c
c	file 'dat':
c
c	  3 3 4 4 4 4 4
c
c
c

c
c	IMPORTANT NOTICE
c	----------------
c
c	The origin of the reference coordinate system (ei) is 
c	situated at P, i.e.:  
c
c	                      Xp = Yp = Zp = 0
c
c	Thus, the coordinates of all vertices of the polyhedron 
c	in file 'xyz' refer to P.
c
c


c
c
c	parameter       stands for
c       --------------------------
c	nop         :   Number of planes
c	nov         :   Number of vertices
c	noe(i)      :   Number of edges (varying from facet to
c			facet according to the information given
c			in file 'dat'. As starting value for 
c			parameter (noe) we choose a conveniently
c			big value (noed), sufficient for our matrices 
c			definitions
c
c
c
c     Parameters 'nop' and 'nov' have to be edited by the user
c
c
c


	implicit real*8 (a-h,o-z)
c
        include 'polyhedron.inc'
c      include 'polyhedron.inc'
c
c Input/Output
c      character*80 infile(4),outfile(1),dumpfile(2),mem(NMESS)
c      integer nov,nop,opt(1),lmem(NMESS)
c      integer*4 nograd,nflush,opflag,l
c      real*8 d0,d,T,omg,x2(novmax),y2(novmax),z2(novmax),xp,yp,zp,J
        integer nov,nop
c
        real*8 d,x2(nov),y2(nov),z2(nov),xp,yp,zp
c
        real*8 gc
c
c	Edit number of planes (nop) and number of vertices (nov)
c
c
c	parameter(nop=7,nov=8,noed=100)


	dimension e(3,3)
c
	dimension noe(nop)
c
	dimension x(novmax),y(novmax),z(novmax),h(nopmax),s(nopmax),
     %  cs(nopmax,3)
	dimension xo(nopmax),yo(nopmax),zo(nopmax)
	dimension xx(nopmax),yy(nopmax),zz(nopmax)
	dimension xxx(nopmax,noed),yyy(nopmax,noed),zzz(nopmax,noed)
c
	dimension k(nopmax,noed)
c
	dimension g(nopmax,noed,3)
	dimension plnorm(nopmax,3),segnorm(nopmax,noed,3)
	dimension dcnp(nopmax,3),dcnpq(nopmax,noed,3)
	dimension hh(nopmax,noed),ss(nopmax,noed)
	dimension cspq(nopmax,noed,3)
	dimension vij(3,3)

	snorm(xipq,yipq,zipq)=dsqrt(xipq*xipq+yipq*yipq+zipq*zipq)
	fln(aa,bb,cc,dd)=dlog((aa+bb)/(cc+dd))
	fan(ee,f,gg,t,yf,w)=datan((ee*t)/(f*w))-datan((ee*gg)/(f*yf))

c
c	open(10,file='xyzposnew',form='formatted',status='unknown')
c	open(11,file='topoaut',form='formatted',status='unknown')
c	open(12,file='dataut',form='formatted',status='unknown')
c      open(13,file='polaut2_out',form='formatted',access='append')

c chama a subrotina mio_in para verificar os nomes e as existências dos arquivos de entrada (.in), saída (.out) e despejo (.dmp), criar os arquivos de saída e despejo e prosseguir com a integração no ponto onde ela foi interrompida
c      call mio_in (infile,outfile,dumpfile,mem,lmem,nov,nop,
c     %  nograd,nflush,opt,d0,T,x2,y2,z2,noe,k,opflag)

c
c	read(12,*) (noe(j),j=1,nop)



c converte a densidade d de g/cm^3 para kg/km^3. OBS: os vértices devem estar em quilômetros
c        d=d0*(1000.0d0*100.0d0)**3/1000.0d0
c encontra a velocidade de rotação omg do asteróide tendo o seu período de rotação T
c        omg=2.0d0*pi/(T*3600.0d0)

c abre o arquivo dos pontos em que se pretende calcular o potencial e suas derivadas (13)
c 450  open  (10+3, file=infile(3), status='old', err=450)
c posiciona a função read para ler o ponto inicial do próximo intervalo a partir do qual se pretendia calcular o potencial e suas derivadas quando a integração foi interrompida
c      opflag = opflag + 1
c      do 20 l=0,opflag-1
c        read(10+3,*) xp,yp,zp
c  20  continue

c
c	do 20 i=1,nov
c	read(10,*) x(i),y(i),z(i)
c20	continue


c	do 30 i=1,nop
c	read(11,*) (k(i,j),j=1,noe(i))
c30	continue


c
c	Vector basis (ei) situated at P
c

	e(1,1)=1
	e(1,2)=0
	e(1,3)=0
	e(2,1)=0
	e(2,2)=1
	e(2,3)=0
	e(3,1)=0
	e(3,2)=0
	e(3,3)=1


c
c	Constants
c
c
c	gc=6.67259e-08
c	d=2.67
c	pi=4.*atan(1.)

c lê os pontos em que se pretende calcular o potencial e suas derivadas iniciando sempre do sucessor do último ponto do último intervalo em que a integração parou
c      nstored = 0
c      do 2013 l=opflag,nograd-1
c	read(10+3,*,end=2020) xp,yp,zp

c no programa original o potencial e suas derivadas são sempre calculados na origem, isto é, no ponto de coordenadas (0,0,0). Logo, para calcularmos o potencial e suas derivadas em qualquer outro ponto P devemos ter nossa origem nesse novo ponto P e para isso basta transladarmos todas as coordenadas dos pontos do poliedro, subtraindo-as das coordenadas do ponto P
	do 2014 i=1,nov
	x(i)=x2(i)-xp
        y(i)=y2(i)-yp
        z(i)=z2(i)-zp
2014	continue
c
c
c	Building vectors Gij // Equation (17)
c
c

	v=0
	vx=0
	vy=0
	vz=0
	do 100 i=1,nop
c	do 80 i=1,nop
	do 70 j=1,noe(i)

	if(j.eq.noe(i)) then 

	g(i,j,1)=x(k(i,1))-x(k(i,j))
	g(i,j,2)=y(k(i,1))-y(k(i,j))
	g(i,j,3)=z(k(i,1))-z(k(i,j))


	else


	g(i,j,1)=x(k(i,j+1))-x(k(i,j))
	g(i,j,2)=y(k(i,j+1))-y(k(i,j))
	g(i,j,3)=z(k(i,j+1))-z(k(i,j))

	end if


70	continue

c
c
c	Computation of plane Normals (Np) // Equation (18)
c
c

c	do 80 i=1,nop
	call cross (g(i,1,1),g(i,1,2),g(i,1,3),g(i,2,1),g(i,2,2),
     $g(i,2,3),plnorm(i,1),plnorm(i,2),plnorm(i,3))
	qnorm=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	plnorm(i,1)=plnorm(i,1)/qnorm
	plnorm(i,2)=plnorm(i,2)/qnorm
	plnorm(i,3)=plnorm(i,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We computed a normed vector, i.e. Np's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c

	dcnp(i,1)=plnorm(i,1)
	dcnp(i,2)=plnorm(i,2)
	dcnp(i,3)=plnorm(i,3)


c80	continue



c
c
c	Computation of line segment normals (npq) // Equation (19)
c
c

c	do 90 i=1,nop
	do 90 j=1,noe(i)
	call cross (g(i,j,1),g(i,j,2),g(i,j,3),plnorm(i,1),
     $plnorm(i,2),plnorm(i,3),segnorm(i,j,1),
     $segnorm(i,j,2),segnorm(i,j,3))
	qnorm=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	segnorm(i,j,1)=segnorm(i,j,1)/qnorm
	segnorm(i,j,2)=segnorm(i,j,2)/qnorm
	segnorm(i,j,3)=segnorm(i,j,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We are computing a normed vector, i.e. npq's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c
	dcnpq(i,j,1)=segnorm(i,j,1)
	dcnpq(i,j,2)=segnorm(i,j,2)
	dcnpq(i,j,3)=segnorm(i,j,3)


90	continue



c
c
c	cos(Np,ei), cos(npq,ei)
c
c
c	cs(i,j) = dcnp(i,j)
c
c	and
c
c	cspq(i,j,m) = dcnpq(i,j,m), m = 1,2,3
c
c

c	do 23 i=1,nop
	qnorm1=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	do 23 j=1,3
	qnorm2=snorm(e(j,1),e(j,2),e(j,3))
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),e(j,1),e(j,2),
     $e(j,3),dotnpi)
	cs(i,j)=dotnpi/(qnorm1*qnorm2)


23	continue
	
c	do 24 i=1,nop
	do 24 j=1,noe(i)
 	do 24 m=1,3
	qnorm1=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	qnorm2=snorm(e(m,1),e(m,2),e(m,3))
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),e(m,1),
     $e(m,2),e(m,3),dotnpq)
	cspq(i,j,m)=dotnpq/(qnorm1*qnorm2)


24	continue

c
c
c	Distances of P from planes Sp ( h(i) ) 
c	by subroutine 'planedist'
c
c	xo(i), yo(i) and zo(i) are the intersections of Sp(i)
c	with axes x, y and z respectively
c
c	and 
c
c	plane normal orientation ( s(i) ) // Equation (20)
c
c

c	do 60 i=1,nop
	call planedist (x,y,z,k(i,1),k(i,2),k(i,3),h(i),xo(i),yo(i),
     $zo(i))


c60	continue



c	do 41 i=1,nop
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),-x(k(i,1)),
     $-y(k(i,1)),-z(k(i,1)),plnrmor)
	if (plnrmor.gt.0) then
	s(i)=-1
	else if (plnrmor.lt.0) then
	s(i)=1
	else if (plnrmor.eq.0) then
	s(i)=0
	end if


c41	continue





c
c
c	Computation of the projections of P on the 'p' planes
c	of the polyhedron. [P'(x',y',z')]
c
c	xx(i), yy(i), zz(i) // Equation (21)
c
c

c	do 31 i=1,nop

	if (xo(i).ge.0.and.yo(i).ge.0.and.zo(i).ge.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (xo(i).lt.0.and.yo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.yo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0.and.zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	end if


c31	continue




c
c	segment normal orientation	
c
c	senrmor = npq (dot) Gpq(ij-->P') // Equation (22)
c

c	do 42 i=1,nop
	do 42 j=1,noe(i)
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),
     $xx(i)-x(k(i,j)),yy(i)-y(k(i,j)),zz(i)-z(k(i,j)),senrmor)
	if (abs(senrmor).lt.1e-10) then
	ss(i,j)=0
	go to 42
	else if (senrmor.gt.0) then
	ss(i,j)=-1
	else if (senrmor.lt.0) then
	ss(i,j)=1
	else if (senrmor.eq.0) then
	ss(i,j)=0
	end if


42	continue



c
c
c	Distances of P' from line segments Gpq (hpq) 
c
c	hh(i,j)
c
c	Computation of the projections of P' on each
c	line segment Gpq, i.e. P''(x'',y''',z''')
c	by calling subroutine 'segmproj'
c
c	xxx(i,j), yyy(i,j), zzz(i,j) // Equations (23) - (25)
c
c

c	do 32 i=1,nop
	do 32 j=1,noe(i)

	if (j.eq.noe(i)) then

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,1)),y(k(i,1)),z(k(i,1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg
	go to 131

	end if

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg

131	continue




32	continue



c
c
c	Computation of V and Vi
c
c

c	v=0
c	vx=0
c	vy=0
c	vz=0

c	do 100 i=1,nop

	sum1=0
	sum2=0
	singarea=0
	singsegm=0
	singvert=0
	
	do 200 j=1,noe(i)




	thisvert=0

	gpqnorm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

c	
c	check if P' lays inside the polygon Np(i)
c

	if (ss(i,j).eq.1) then
	singarea=singarea+1
	end if

c
c	check if P' lays inside the segment Gpq 
c
c	or
c
c	on the vertex (i,j)
c
c	Condition controlled is  
c
c	ss(i,j) = 0
c

	if (ss(i,j).eq.0) then



	 if(j.eq.noe(i)) then 

	 gx1=xx(i)-x(k(i,1))
	 gx2=xx(i)-x(k(i,j))
	 gy1=yy(i)-y(k(i,1))
	 gy2=yy(i)-y(k(i,j))
	 gz1=zz(i)-z(k(i,1))
	 gz2=zz(i)-z(k(i,j))

	 else  

	 gx2=xx(i)-x(k(i,j))
	 gx1=xx(i)-x(k(i,j+1))
	 gy2=yy(i)-y(k(i,j))
	 gy1=yy(i)-y(k(i,j+1))
	 gz2=zz(i)-z(k(i,j))
	 gz1=zz(i)-z(k(i,j+1))

	 end if

	 enorm1=snorm(gx1,gy1,gz1)
	 enorm2=snorm(gx2,gy2,gz2)




c
c	check for the Vertex (i,j) 
c

	  if (enorm1.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.noe(i)) then

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,1,1),g(i,1,2),
     $g(i,1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,j+1,1),g(i,j+1,2),g(i,j+1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,j+1,1),g(i,j+1,2),
     $g(i,j+1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if


	  if (enorm2.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.1) then

	   g1norm=snorm(g(i,noe(i),1),g(i,noe(i),2),g(i,noe(i),3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,noe(i),1),-g(i,noe(i),2),-g(i,noe(i),3),
     $g(i,1,1),g(i,1,2),g(i,1,3),gdot)
     
	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j-1,1),g(i,j-1,2),g(i,j-1,3))
	   g2norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

	call dot(-g(i,j-1,1),-g(i,j-1,2),-g(i,j-1,3),g(i,j,1),g(i,j,2),
     $g(i,j,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if




	      if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

	      singsegm=singsegm+1

	      end if




	end if









c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from field point P.
c
c						



	if (j.eq.noe(i)) then

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,1)),yyy(i,j)-y(k(i,1)),
     $zzz(i,j)-z(k(i,1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,1)),y(k(i,1)),z(k(i,1)))

	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1968
	end if


	else 

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,j+1)),yyy(i,j)-y(k(i,j+1)),
     $zzz(i,j)-z(k(i,j+1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1968
	end if



	end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1972
	end if

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if


	go to 1972






1968	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	else if (d1.eq.d2) then

	s1=-d1
	s2=d2
	r1=-f1
	r2=f2

	end if


1972	continue


	if (h(i).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1973
	   end if


	if ((s1+s2).lt.1e-10.and.(r1+r2).lt.1e-10) then
	c1=0
	c2=0
	go to 1973
	end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	
	go to 1973

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(i,j).eq.0) then
	c2=0
	else
	c2=fan(h(i),hh(i,j),s1,s2,r1,r2)
	end if

1973	continue

	sum1=sum1+ss(i,j)*hh(i,j)*c1
	sum2=sum2+ss(i,j)*c2



200	continue


c
c	If P' lays inside the polygon Ap then
c	ALL of ss(i,j) = 1. In that case sing = noe(i)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(i)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(i)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(i) 
c	where
c
c	theta = arccos{ (g1 (dot) g2) / (|g1||g2|) }, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c
c	theta has already been computed (see above)
c


	if (singarea.eq.noe(i)) then
	area=2*pi*h(i)
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi*h(i)
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta*h(i)
	else
	vertex=0
	end if


c     
c     Final computation of Equations (10) and (11)
c


	v=v+s(i)*h(i)*(sum1+h(i)*sum2-area-edge-vertex)
	vx=vx+cs(i,1)*(sum1+h(i)*sum2-area-edge-vertex)
	vy=vy+cs(i,2)*(sum1+h(i)*sum2-area-edge-vertex)
	vz=vz+cs(i,3)*(sum1+h(i)*sum2-area-edge-vertex)




100	continue

c
c	v=(v*gc*d)/2
c	vx=abs(vx*gc*d)
c	vy=abs(vy*gc*d)
c	vz=abs(vz*gc*d)
	v=(v*abs(gc*d))/2
	vx=vx*abs(gc*d)
	vy=vy*abs(gc*d)
	vz=vz*abs(gc*d)




c
c
c	Computation of Vij
c
c
	do 901 i=1,3
	do 901 j=1,3
	vij(i,j)=0
901	continue

	do 400 i=1,3
	do 500 j=1,3


		do 600 m=1,nop

		sum1=0
		sum2=0

		singarea=0
		singsegm=0
		singvert=0



		do 700 n=1,noe(m)

		thisvert=0

		gpqnorm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))

c	
c		check if P' lays inside the polygon Np(m)
c

		if (ss(m,n).eq.1) then
		singarea=singarea+1
		end if



c
c		check if P' lays inside the segment Gpq
c
c		or
c
c		on the vertex (m,n)
c
c		Condition controlled is  
c
c		ss(m,n) = 0
c

		if (ss(m,n).eq.0) then

		if (n.eq.noe(m)) then 

		gx1=xx(m)-x(k(m,1))
		gx2=xx(m)-x(k(m,n))
		gy1=yy(m)-y(k(m,1))
		gy2=yy(m)-y(k(m,n))
		gz1=zz(m)-z(k(m,1))
		gz2=zz(m)-z(k(m,n))

		else  

		gx2=xx(m)-x(k(m,n))
		gx1=xx(m)-x(k(m,n+1))
		gy2=yy(m)-y(k(m,n))
		gy1=yy(m)-y(k(m,n+1))
		gz2=zz(m)-z(k(m,n))
		gz1=zz(m)-z(k(m,n+1))

		end if

		enorm1=snorm(gx1,gy1,gz1)
		enorm2=snorm(gx2,gy2,gz2)



c
c		check if P' lays on the Vertex (m,n) 
c

		if (enorm1.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.noe(m)) then

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))

	


	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,1,1),g(m,1,2),
     $g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,n+1,1),g(m,n+1,2),g(m,n+1,3))




	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,n+1,1),g(m,n+1,2),
     $g(m,n+1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if



		if (enorm2.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.1) then

	g1norm=snorm(g(m,noe(m),1),g(m,noe(m),2),g(m,noe(m),3))
	g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))


	call dot(-g(m,noe(m),1),-g(m,noe(m),2),-g(m,noe(m),3),
     $g(m,1,1),g(m,1,2),g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

	g1norm=snorm(g(m,n-1,1),g(m,n-1,2),g(m,n-1,3))
	g2norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))


	call dot(-g(m,n-1,1),-g(m,n-1,2),-g(m,n-1,3),g(m,n,1),g(m,n,2),
     $g(m,n,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if

		if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

		singsegm=singsegm+1

		end if

		end if





c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from filed point P.
c
c



		if (n.eq.noe(m)) then

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,1)),yyy(m,n)-
     $y(k(m,1)),
     $zzz(m,n)-z(k(m,1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,1)),y(k(m,1)),z(k(m,1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1958
	end if

		else

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,n+1)),
     $yyy(m,n)-y(k(m,n+1)),zzz(m,n)-z(k(m,n+1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,n+1)),y(k(m,n+1)),z(k(m,n+1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1958
	end if

		end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1982
	end if

c
c	1-D Coordinate system orientation
c

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if

	go to 1982

1958	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	end if


1982	continue

	if (h(m).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1983
	   end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	go to 1983

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(m,n).eq.0) then
	c2=0
	else
	c2=fan(h(m),hh(m,n),s1,s2,r1,r2)
	end if

1983	continue




		sum1=sum1+cspq(m,n,j)*c1
		sum2=sum2+ss(m,n)*c2



700	continue



c
c	If P' lays inside the polygon Ap then
c	ALL of ss(m,n) = 1. In that case sing = noe(m)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(m)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(m)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(m) 
c	where
c
c	theta = arccos{(g1 (dot) g2) / (|g1||g2|)}, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c


	if (singarea.eq.noe(m)) then
	area=2*pi
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta
	else
	vertex=0
	end if


c     
c     Final computation of Equations (12)
c

	
	vij(i,j)=vij(i,j)+cs(m,i)*(sum1+s(m)*cs(m,j)*(sum2-area
     $-edge-vertex))



600	continue

500	continue
400	continue


	do 900 i=1,3
	do 900 j=1,3
	vij(i,j)=vij(i,j)*gc*d
900	continue


c
c	write(13,*) 'V = ',v
c	write(13,*) 'Vx = ',vx
c	write(13,*) 'Vy = ',vy
c	write(13,*) 'Vz = ',vz
c	write(13,*) 'Vxx = ',vij(1,1)
c	write(13,*) 'Vyy = ',vij(2,2)
c	write(13,*) 'Vzz = ',vij(3,3)
c	write(13,*) 'Vxy = ',vij(1,2)
c	write(13,*) 'Vxz = ',vij(1,3)
c	write(13,*) 'Vyz = ',vij(2,3)

c cálculo do pseudo potencial (D. J. SCHEERES AND S. J. OSTRO - Orbits Close to Asteroid 4769 Castalia (equações 18 e 19))
c        J = -omg*omg*(xp*xp + yp*yp)/2.0d0 -v

c
c Output data for all bodies
c      call mio_out (opt,nstored,l,xp,yp,zp,v,vx,vy,vz,
c     %  vij,omg,nflush,nograd,outfile,dumpfile,mem,lmem,nov,
c     %  nop,d0,T)

c fim da leitura dos pontos em que se pretende calcular o potencial e suas derivadas
c2013  continue
c fecha o arquivo de entrada dos pontos em que se pretende calcular o potencial e suas derivadas (13)
c2020    close (10+3)
c termina a execução do programa
c      write (*,'(a)') mem(10)(1:lmem(10))
c	stop
c1001	format(3i4)
c1002	format(7(4i2)/)

        return
	end


c
c
c	subroutines
c
c
	subroutine planedist (x,y,z,kp1,kp2,kp3,hp,xax,yax,zax)
c
c	computes the distance of P(xp,yp,zp) from 
c	the plane a*x+b*y+c*z+d=0
c	for our case (xp=yp=zp=0) it is
c	h = d / sqrt(a^2+b^2+c^2)
c
	implicit real*8 (a-h,o-z)
	dimension x(8),y(8),z(8)
	a=(y(kp2)-y(kp1))*(z(kp3)-z(kp1))-
     $(y(kp3)-y(kp1))*(z(kp2)-z(kp1))
	b=(z(kp2)-z(kp1))*(x(kp3)-x(kp1))-
     $(z(kp3)-z(kp1))*(x(kp2)-x(kp1))
	c=(x(kp2)-x(kp1))*(y(kp3)-y(kp1))-
     $(x(kp3)-x(kp1))*(y(kp2)-y(kp1))
	d=-a*x(kp1)-b*y(kp1)-c*z(kp1)
	hp=abs(d/dsqrt(a*a+b*b+c*c))

	if (z(kp1).eq.z(kp2).and.z(kp2).eq.z(kp3)) then
	xax=0
	yax=0
	zax=-d/c
	else if (y(kp1).eq.y(kp2).and.y(kp2).eq.y(kp3)) then
	xax=0
	yax=-d/b
	zax=0
	else if (x(kp1).eq.x(kp2).and.x(kp2).eq.x(kp3)) then
	xax=-d/a
	yax=0
	zax=0
	end if

	if (a.eq.0) then
	xax=0
	else
	xax=-d/a
	end if

	if (b.eq.0) then
	yax=0
	else
	yax=-d/b
	end if

	if (c.eq.0) then
	zax=0
	else
	zax=-d/c
	end if



	return
	end


	subroutine segmproj(x1,y1,z1,x2,y2,z2,xo,yo,zo,xs,ys,zs)
c
c	Computes the projections of P'(xx(i),yy(i),zz(i)) from
c	each line segment lpq, i.e. the coordinate values
c	P''(xxx(i,j),yyy(i,j),zzz(i,j)). P'' is considered to lie
c	on the intersection of 3 planes (see notes), namely
c	plane XoX1X2, plane XX1X2 and the plane perpendicular to
c	vector X1X2. Building the respective plane equations, we get
c	a 3x3 system of equations whose solution gives the 
c	vector xxx. The equations (23) - (25) build this system.
c
	implicit real*8 (a-h,o-z)
	r1=x2-x1
	r2=y2-y1
	r3=z2-z1
	e1=x1-xo
	e2=y1-yo
	e3=z1-zo
	a1=x2-x1
	b1=y2-y1
	c1=z2-z1
	d1=a1*xo+b1*yo+c1*zo
 	call cross (e1,e2,e3,r1,r2,r3,a2,b2,c2)	
	call cross (a2,b2,c2,r1,r2,r3,a3,b3,c3)
	call dot (a2,b2,c2,xo,yo,zo,d2)
	call dot (a3,b3,c3,x1,y1,z1,d3)
	det=a1*b2*c3+b1*c2*a3+c1*a2*b3-a3*b2*c1-b3*c2*a1-c3*a2*b1
	dx=d1*b2*c3+b1*c2*d3+c1*d2*b3-d3*b2*c1-b3*c2*d1-c3*d2*b1
	dy=a1*d2*c3+d1*c2*a3+c1*a2*d3-a3*d2*c1-d3*c2*a1-c3*a2*d1
	dz=a1*b2*d3+b1*d2*a3+d1*a2*b3-a3*b2*d1-b3*d2*a1-d3*a2*b1
	xs=dx/det
	ys=dy/det
	zs=dz/det
	return
	end




	subroutine cross (a1,a2,a3,b1,b2,b3,cr1,cr2,cr3)
c
c	Cross product A x B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cr1=a2*b3-a3*b2
	cr2=a3*b1-a1*b3
	cr3=a1*b2-a2*b1
	return
	end


	subroutine dot (a1,a2,a3,b1,b2,b3,cdot)
c
c	Dot product A . B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cdot=a1*b1+a2*b2+a3*b3
	return
	end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_ELAPSE.FOR    (FEG   30 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print out elapsed time of Mercury integration.
c
c------------------------------------------------------------------------------
c
      subroutine mio_elapse (treal,tuser,tsys,nmem,mem,lmem)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/output
      integer lmem(NMESS),nmem
      real*8 treal,tuser,tsys
      character*80 mem(NMESS)
c
c Local
      integer tday,thour,tmin
      real*8 tsec
      character*1000 string
c
c------------------------------------------------------------------------------
c
      write(6,'(a)') mem(nmem)(1:lmem(nmem))
c
      call dhms(treal,tday,thour,tmin,tsec)
      write(6,'(a,i5,a,2(i2,a),f6.3,a)') mem(84)(1:lmem(84)),tday,
     %  mem(87)(1:lmem(87)),thour,mem(88)(1:lmem(88)),tmin,
     %  mem(89)(1:lmem(89)),tsec,mem(90)(1:lmem(90))
c
      call dhms(tuser,tday,thour,tmin,tsec)
      write(6,'(a,i5,a,2(i2,a),f6.3,a)') mem(85)(1:lmem(85)),tday,
     %  mem(87)(1:lmem(87)),thour,mem(88)(1:lmem(88)),tmin,
     %  mem(89)(1:lmem(89)),tsec,mem(90)(1:lmem(90))
c
      call dhms(tsys,tday,thour,tmin,tsec)
      write(6,'(a,i5,a,2(i2,a),f6.3,a)') mem(86)(1:lmem(86)),tday,
     %  mem(87)(1:lmem(87)),thour,mem(88)(1:lmem(88)),tmin,
     %  mem(89)(1:lmem(89)),tsec,mem(90)(1:lmem(90))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_REMAINING.FOR    (FEG   30 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print out integration time remaining.
c
c------------------------------------------------------------------------------
c
      subroutine mio_remaining (treal,opflag,nograd,mem,lmem)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/output
      integer lmem(NMESS)
      integer*4 opflag,nograd
      real*8 treal
      character*80 mem(NMESS)
c
c Local
      real*8 trem
      integer tday,thour,tmin
      real*8 tsec
      character*1000 string
      integer size
c
c------------------------------------------------------------------------------
c
      if (opflag.ne.0) then
        trem = treal * (abs(nograd-opflag) / (1.d0 * opflag))
        call dhms(trem,tday,thour,tmin,tsec)
        write(6,'(a,i5,a,2(i2,a),f6.3,a)') mem(91)(1:lmem(91)),tday,
     %    mem(87)(1:lmem(87)),thour,mem(88)(1:lmem(88)),tmin,
     %    mem(89)(1:lmem(89)),tsec,mem(90)(1:lmem(90))
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_SIZEOUT.FOR    (FEG   30 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print out size of Mercury outfiles.
c
c------------------------------------------------------------------------------
c
      subroutine mio_sizeout (opt,opflag,nograd,nflush,outfile,mem,lmem)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/output
      integer opt,lmem(NMESS)
      integer*4 opflag,nograd,nflush
      character*80 mem(NMESS),outfile(1)
c
c Local
      character*1000 string
      integer j,nchar
      integer*8 size
c
c------------------------------------------------------------------------------
c
      if (opflag.eq.0) then
        if (opt.eq.1) nchar = 3
        if (opt.eq.2) nchar = 5
        if (opt.eq.3) nchar = 8
c
        size = 14 * (int ((1.d0 * nograd) / nflush))
        size = size + (25+10*nchar-1) * nograd
c
        write(6,'(a,i15,1x,2(a),a80)') mem(92)(1:lmem(92)),size,
     %    mem(94)(1:lmem(94)),mem(95)(1:lmem(95)),outfile(1)
      else
        write(6,'(a)') mem(93)(1:lmem(93))
        do j = 1, 1
          string(1:6) = 'ls -l '
          string(7:87) = outfile(j)
          string(88:128) = '| awk ''{printf $5 " bytes -> " $9 "\n"}'''
          call command ( string, 128 )
        end do
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_PID.FOR    (FEG   30 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print out PID process of Mercury.
c
c------------------------------------------------------------------------------
c
      subroutine mio_pid ( )
c
      implicit none
c
c Local
      integer pid
      real*8 nbr
      integer ex
      character*1000 string
c
c------------------------------------------------------------------------------
c
      call get_pid (pid)
      string(1:6) = 'ps -p '
      nbr = pid
      if (nbr.ne.0) then
        ex = int(log10(nbr))
      else
        ex = 0
      end if
      ex = ex + 1
      write (string(7:7+ex-1),'(i0)') pid
      string(7+ex:43+ex) = ' -o %cpu,%mem,cputime,etime,bsdstart,'
      string(44+ex:66+ex) = 'pid,comm,state,ni,euser'
c      string(1:43) = 'ps -C mercury_guara4 -o %cpu,%mem,bsdstart,'
c      string(44:84) = 'bsdtime,c,comm,command,cp,cputime,egroup,'
c      string(85:125) = 'etime,euid,euser,gid,lstart,ni,pcpu,pgid,'
c      string(126:167) = 'pid,pmem,ppid,rss,ruid,size,start,sz,time,'
c      string(168:209) = 'tname,vsz,args,nice,rgroup,ruser,tty,user,'
c      string(210:251) = 'addr,atime,flags,pri,rgid,sid,state,stime,'
c      string(252:265) = 'uid,wchan,xpid'
      call command ( string, 66+ex )
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_HOST.FOR    (FEG   22 Ago 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Print out hostname information.
c
c------------------------------------------------------------------------------
c
      subroutine mio_host ( )
c
      implicit none
c
c Local
      character*1000 string
c
c------------------------------------------------------------------------------
c
      string(1:20)    = 'echo "$(hostname -I '
      string(21:54)   = '| awk ''{printf $1 " " $2 " "}''; '
      string(55:66)   = 'hostname -i '
      string(67:93)   = '| awk ''{printf $1 " "}''; '
      string(94:105)   = 'hostname -s '
      string(106:132)  = '| awk ''{printf $1 " "}''; '
      string(133:144) = 'hostname -d '
      string(145:172) = '| awk ''{printf $1 "\n"}'')"'

      call command ( string, 172)
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      CROSS_AST.FOR    (4 October 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Computes cross product of faces.
c
c------------------------------------------------------------------------------
c
      subroutine cross_ast (nov,nop,noe,k,xv,yv,zv,face)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 face(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j
      real*8 AB(3,noed),norm(3)
      real*8 prod,temp,temp1,temp2,xc(3)
c
c------------------------------------------------------------------------------
c
      do i = 1, nop
        do j = 1, noe(i)-1
          AB(1,j) = xv(k(i,j+1))-xv(k(i,j))
          AB(2,j) = yv(k(i,j+1))-yv(k(i,j))
          AB(3,j) = zv(k(i,j+1))-zv(k(i,j))
        end do
        norm(1) = AB(2,1)*AB(3,2)-AB(3,1)*AB(2,2)
        norm(2) = AB(3,1)*AB(1,2)-AB(1,1)*AB(3,2)
        norm(3) = AB(1,1)*AB(2,2)-AB(2,1)*AB(1,2)
c        prod = xv(k(i,1))*norm(1)+yv(k(i,1))*norm(2)+
c     %    zv(k(i,1))*norm(3)
        xc(1)   = (xv(k(i,1))+xv(k(i,2))+xv(k(i,3)))/3.d0
        xc(2)   = (yv(k(i,1))+yv(k(i,2))+yv(k(i,3)))/3.d0
        xc(3)   = (zv(k(i,1))+zv(k(i,2))+zv(k(i,3)))/3.d0
        prod = xc(1)*norm(1)+xc(2)*norm(2)+xc(3)*norm(3)
c        temp1 = dsqrt(xv(k(i,1))*xv(k(i,1))+yv(k(i,1))*yv(k(i,1))+
c     %    zv(k(i,1))*zv(k(i,1)))
c        temp1 = dsqrt(xc(1)*xc(1)+xc(2)*xc(2)+xc(3)*xc(3))
c        temp2 = dsqrt(norm(1)*norm(1)+norm(2)*norm(2)+norm(3)*norm(3))
c        temp = prod/(temp1*temp2)
c        temp = acos(temp)
c        temp = temp / DR
c        write(*,*) i,temp
        face(i) = sign(1.d0,prod)
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      REORIENTATION.FOR    (4 October 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Reorientate faces.
c
c------------------------------------------------------------------------------
c
      subroutine reorientation (face1,face2,nov,nop,noe,k)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 face1(nop),face2(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j,l,n
      integer k0
      dimension k0(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      n = 0
      do i = 1, nop
        if (face1(i).ne.face2(i)) n = n + 1
      end do
c
      if (n.eq.nop) then
        do i = 1, nop
          l = 0
          do j = noe(i), 1, -1
            l = l + 1
            k0(i,l) = k(i,j)
          end do
          do j = 1, noe(i)
            k(i,j) = k0(i,j)
          end do
        end do
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT1.FOR    (7 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Translates vertices to the center of mass.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot1 (nov,xv,yv,zv,xc,yc,zc)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc,yc,zc
c
c Local
      integer i,j,k
      real*8 x(2,3),xcm(3)
c
c------------------------------------------------------------------------------
c
      xcm(1) = xc
      xcm(2) = yc
      xcm(3) = zc
c
c computes transposed matrix
c      do i = 1, 3
c        do j = 1, 3
c          TT(i,j) = T(j,i)
c        enddo
c      enddo

c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
        do j = 1, 3
          x(1,j) = x(1,j) - xcm(j)
        enddo
c
        xv(i) = x(1,1)
        yv(i) = x(1,2)
        zv(i) = x(1,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_ROT2.FOR    (7 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Rotates to the main inertia axes.
c
c------------------------------------------------------------------------------
c
      subroutine masc_rot2 (nov,xv,yv,zv,T)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 T(3,3)
c
c Local
      integer i,j,k
      real*8 x(2,3),TT(3,3)
c
c------------------------------------------------------------------------------
c
      do i = 1, nov
c
c translate
        x(1,1) = xv(i)
        x(1,2) = yv(i)
        x(1,3) = zv(i)
c
c rotate
        do j = 1, 3
          x(2,j) = 0.d0
          do k = 1, 3
c            x(2,j) = x(2,j) + TT(j,k) * x(1,k)
            x(2,j) = x(2,j) + T(j,k) * x(1,k)
          enddo
        enddo
        xv(i) = x(2,1)
        yv(i) = x(2,2)
        zv(i) = x(2,3)
      enddo
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      REORIENTATION2.FOR    (9 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Reorientate faces.
c
c------------------------------------------------------------------------------
c
      subroutine reorientation2 (face1,face2,nov,nop,noe,k)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer nov,nop,noe,k
      real*8 face1(nop),face2(nop)
      dimension noe(nop)
      dimension k(nopmax,noed)
c Local
      integer i,j,l
      integer k0
      dimension k0(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      do i = 1, nop
        if (face1(i)*face2(i).eq.-1) then
          l = 0
          do j = noe(i), 1, -1
            l = l + 1
            k0(i,l) = k(i,j)
          end do
          do j = 1, noe(i)
            k(i,j) = k0(i,j)
          end do
        end if
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER2.FOR    (9 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer2 (nov,nop,xv,yv,zv,noe,k,d,lay,m,mt,vt,
     %  xc,yc,zc,mem,lmem,l0,T1t,T2t,TPt,Jt,cubx0,cuby0,cubz0)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,lay,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),m(nocen),d(nocen),mt,vt
      real*8 T1t(3),T2t(3),TPt(3),Jt(3,3)
      real*8 cubx0,cuby0,cubz0
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,l,ilay,novl,nopl
      integer noel,kl
      real*8 del,inc,nx,ny,nz
      real*8 xl(nov),yl(nov),zl(nov)
      real*8 norm(nopmax,3),wfac(nopmax)
      real*8 T0,T1(3),T2(3),TP(3)
      real*8 fvt,fmt,vol,vol1,vol2
      real*8 J0(3,3)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor
      integer astflag1,astflag2
      real*8 xvc(8),yvc(8),zvc(8)
      real*8 cubx,cuby,cubz
      real*8 xl2(nov),yl2(nov),zl2(nov)
      dimension noel(nopmax)
      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
      del = 1.d0 / lay
      vt = 0.d0
      mt = 0.d0
      l0 = 0
c
      T1t(1) = 0.0
      T1t(2) = 0.0
      T1t(3) = 0.0
      T2t(1) = 0.0
      T2t(2) = 0.0
      T2t(3) = 0.0
      TPt(1) = 0.0
      TPt(2) = 0.0
      TPt(3) = 0.0
c central point (origin)
      xl(1) = 0.0
      yl(1) = 0.0
      zl(1) = 0.0
c
c      cubx = 1.d0 / cubx0
c      cuby = 1.d0 / cuby0
c      cubz = 1.d0 / cubz0
      cubx = cubx0
      cuby = cuby0
      cubz = cubz0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
      do ilay = 1, lay
        inc = del * ilay
c
        do i = 1, nov
c
c compute new vertices
          xl(i) = xv(i) * inc
          yl(i) = yv(i) * inc
          zl(i) = zv(i) * inc
          if (ilay.gt.1) then
            xl2(i) = xv(i) * del * (ilay-1)
            yl2(i) = yv(i) * del * (ilay-1)
            zl2(i) = zv(i) * del * (ilay-1)
          end if
c
c          if (i.eq.1) then
c            xmin=xl(i)
c            xmax=xl(i)
c            ymin=yl(i)
c            ymax=yl(i)
c            zmin=zl(i)
c            zmax=zl(i)
c          end if
c
c          if(xl(i).lt.xmin) xmin = xl(i)
c          if(xl(i).gt.xmax) xmax = xl(i)
c          if(yl(i).lt.ymin) ymin = yl(i)
c          if(yl(i).gt.ymax) ymax = yl(i)
c          if(zl(i).lt.zmin) zmin = zl(i)
c          if(zl(i).gt.zmax) zmax = zl(i)
        end do
c
        fvt = 0.d0
        fmt = 0.d0
c
        zstep = zmin
        do while (zstep.lt.zmax)
          ystep = ymin
          do while (ystep.lt.ymax)
            xstep = xmin
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
              call tetrah (nov,nop,xl,yl,zl,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c
              if (astflag1.eq.1) then
                if (ilay.le.1) then
                  astflag2 = 0
                else
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
                  call tetrah (nov,nop,xl2,yl2,zl2,noe,k,
     %              xcoor,ycoor,zcoor,astflag2)
                end if
                if (astflag2.eq.0) then
c compute cubic vertices
                xvc(1) = xstep
                yvc(1) = ystep
                zvc(1) = zstep
                xvc(2) = xstep + cubx
                yvc(2) = ystep
                zvc(2) = zstep
                xvc(3) = xstep + cubx
                yvc(3) = ystep + cuby
                zvc(3) = zstep
                xvc(4) = xstep
                yvc(4) = ystep + cuby
                zvc(4) = zstep
                xvc(5) = xstep
                yvc(5) = ystep
                zvc(5) = zstep + cubz
                xvc(6) = xstep + cubx
                yvc(6) = ystep
                zvc(6) = zstep + cubz
                xvc(7) = xstep + cubx
                yvc(7) = ystep + cuby
                zvc(7) = zstep + cubz
                xvc(8) = xstep
                yvc(8) = ystep + cuby
                zvc(8) = zstep + cubz
c
c compute cubic's faces
                noel(1) = 4
                kl(1,1) = 4
                kl(1,2) = 3
                kl(1,3) = 2
                kl(1,4) = 1
c
                noel(2) = 4
                kl(2,1) = 1
                kl(2,2) = 2
                kl(2,3) = 6
                kl(2,4) = 5
c
                noel(3) = 4
                kl(3,1) = 2
                kl(3,2) = 3
                kl(3,3) = 7
                kl(3,4) = 6
c
                noel(4) = 4
                kl(4,1) = 3
                kl(4,2) = 4
                kl(4,3) = 8
                kl(4,4) = 7
c
                noel(5) = 4
                kl(5,1) = 4
                kl(5,2) = 1
                kl(5,3) = 5
                kl(5,4) = 8
c
                noel(6) = 4
                kl(6,1) = 5
                kl(6,2) = 6
                kl(6,3) = 7
                kl(6,4) = 8
c
c compute volumens
                nopl = 6
c compute normal faces and w vector
                do j = 1, nopl
                  call compnormalface (nov,nop,j,xvc,yvc,zvc,kl,
     %              nx,ny,nz)
                  norm(j,1) = nx
                  norm(j,2) = ny
                  norm(j,3) = nz
                  wfac(j) = - norm(j,1) * xvc(kl(j,1))
     %                      - norm(j,2) * yvc(kl(j,1))
     %                      - norm(j,3) * zvc(kl(j,1))
                enddo
c compute volumen of a cube
                call compVolumeIntegrals (nov,nop,nopl,xvc,yvc,zvc,
     %            noel,kl,norm,wfac,T0,T1,T2,TP)
                fvt = fvt + T0
                T1t(1) = T1t(1) + T1(1)
                T1t(2) = T1t(2) + T1(2)
                T1t(3) = T1t(3) + T1(3)
                T2t(1) = T2t(1) + T2(1)
                T2t(2) = T2t(2) + T2(2)
                T2t(3) = T2t(3) + T2(3)
                TPt(1) = TPt(1) + TP(1)
                TPt(2) = TPt(2) + TP(2)
                TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
                l0 = l0 + 1
                m(l0) = d(ilay) * T0
                fmt = fmt + m(l0)
c compute center of masses of a cube
                call compcenpolyhedron (d(ilay),m(l0),T0,T1,
     %            T2,TP,xc(l0),yc(l0),zc(l0),J0)
                Jt(1,1) = Jt(1,1) + J0(1,1)
                Jt(1,2) = Jt(1,2) + J0(1,2)
                Jt(1,3) = Jt(1,3) + J0(1,3)
                Jt(2,1) = Jt(2,1) + J0(2,1)
                Jt(2,2) = Jt(2,2) + J0(2,2)
                Jt(2,3) = Jt(2,3) + J0(2,3)
                Jt(3,1) = Jt(3,1) + J0(3,1)
                Jt(3,2) = Jt(3,2) + J0(3,2)
                Jt(3,3) = Jt(3,3) + J0(3,3)
              end if
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
        vt = vt + fvt
        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(15),
     %    lmem(15),' ',1,mem(14),lmem(14))
      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READBINPOL2B.FOR    (16 October 2017)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read computed center of mass from cube file.
c
c------------------------------------------------------------------------------
c
      subroutine readbinpol2b (infile,mem,lmem,noc,mcen,vc,J0,
     %  xc,yc,zc,mc)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(6)
      real*8 xc(noc),yc(noc),zc(noc),mc(noc)
      real*8 mcen,vc,J0(3,3)
      integer noc
c
c Local
      integer lim(2,1000),nsub
      integer j,itmp
      character*15000 string
      character*80 c80
      integer lineno
      character*7 c7
      logical test
c
c------------------------------------------------------------------------------
c
      inquire (file=infile(6), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(6)
        stop
      end if
      open (10, file=infile(6), status='old')
      lineno = 0
      write (*,'(a13,1x,a)') 'Reading file:',infile(6)
c
      open (10, file=infile(6), status='old', access='sequential')
c
  40  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 40
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 40
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) noc
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) mcen
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) vc
c
  50  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 50
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 50
c
      do j = 1, 3
  60    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 60
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 60
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) J0(j,1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) J0(j,2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) J0(j,3)
      end do
c
      do j = 1, noc
  70    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 70
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 70
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xc(j)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yc(j)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zc(j)
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mc(j)
      end do
      close (10)
c
 667  continue
      if(lineno.eq.0) then
        write(6,'(/,3(1x),2(a,1x),/)') 'Error:',
     %    'No center mass point was found.'
        stop
      end if      
c
c termina a execução do programa
      write (6,'(a)') 'Process completed successfully!'
c
c------------------------------------------------------------------------------
c
      return
c
c------------------------------------------------------------------------------
c
c Error reading from the input file containing integration parameters
 661  write (c7,'(i7)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c7,7,
     %  mem(7),lmem(7))
c
c------------------------------------------------------------------------------
c
      end
c
c------------------------------------------------------------------------------
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      TETRAH.FOR    (6 September 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Given a point and a polyhedron, check if the point is inside (iflag=1) or
c outside (iflag=0) the polyhedron using the tetrahedron method.
c
c------------------------------------------------------------------------------
c
      subroutine tetrah (nov,nop,xv,yv,zv,noe,k,x,y,z,iflag)
c
      implicit none
      include 'polyhedron.inc'
      real*8 EPS
      parameter (EPS = 0.00001)
c
c Input/Output
      integer nov,nop,noe,k,iflag
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 x,y,z
      dimension noe(nop)
      dimension k(nopmax,noed)
c
c Local
      integer i,j,flag,cross
      integer inout(5)
      real*8 A(4,4),D(5),Dt,b(4),Px,Py,Pz
      real*8 GETDET
c
c------------------------------------------------------------------------------
c
      Px = x
      Py = y
      Pz = z
c
1234  continue
      cross = 0
      do i = 1, nop
        A(1,1) = 0.d0
        A(1,2) = 0.d0
        A(1,3) = 0.d0
        A(1,4) = 1.d0
        do j = 1, noe(i)
          A(j+1,1) = xv(k(i,j))
          A(j+1,2) = yv(k(i,j))
          A(j+1,3) = zv(k(i,j))
          A(j+1,4) = 1.d0
        end do
        D(1) = GETDET(A,4)
        if (abs(D(1)).le.TINY) cycle
        inout(1) = int(sign(1.d0,D(1)))
c
        do j = 1, noe(i)
          if (j.eq.noe(i)) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(5) = GETDET(A,4)
        if (abs(D(5)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(5) = int(sign(1.d0,D(5)))
c        b(4) = D(5)/D(1)
c
        do j = 1, noe(i)
          if (j.eq.noe(i)-1) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(4) = GETDET(A,4)
        if (abs(D(4)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(4) = int(sign(1.d0,D(4)))
c        b(3) = D(4)/D(1)
c
        do j = 1, noe(i)
          if (j.eq.noe(i)-2) then
            A(j+1,1) = Px
            A(j+1,2) = Py
            A(j+1,3) = Pz
            A(j+1,4) = 1.d0
          else
            A(j+1,1) = xv(k(i,j))
            A(j+1,2) = yv(k(i,j))
            A(j+1,3) = zv(k(i,j))
            A(j+1,4) = 1.d0
          end if
        end do
        D(3) = GETDET(A,4)
        if (abs(D(3)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(3) = int(sign(1.d0,D(3)))
c        b(2) = D(3)/D(1)
c
        A(1,1) = Px
        A(1,2) = Py
        A(1,3) = Pz
        A(1,4) = 1.d0
        do j = 1, noe(i)
          A(j+1,1) = xv(k(i,j))
          A(j+1,2) = yv(k(i,j))
          A(j+1,3) = zv(k(i,j))
          A(j+1,4) = 1.d0
        end do
        D(2) = GETDET(A,4)
        if (abs(D(2)).le.TINY) then
          Px = Px + EPS
          Py = Py + EPS
          Pz = Pz + EPS
          goto 1234
        end if
        inout(2) = int(sign(1.d0,D(2)))
c        b(1) = D(2)/D(1)
c
c        Dt = D(2)+D(3)+D(4)+D(5)
c
c        write(*,*) i,D(1),inout(1)
c        write(*,*) D(2),D(3),D(4),D(5)
c        write(*,*) inout(2),inout(3),inout(4),inout(5)
c        write(*,*) D(1),Dt
c        write(*,*) b(1),b(2),b(3),b(4)
c
        flag = 1
        do j = 2, 5
          if (inout(1)*inout(j).lt.0) flag = 0
        end do
c        if (flag.eq.1) exit
        if (flag.eq.1) cross = cross + 1
c
        if (cross.eq.1) exit
      end do
c
      iflag = 1
      if (mod(cross,2).eq.0) iflag = 0
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      GETDET.FOR    (6 September 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
C A general purpose function written in FORTRAN77 to calculate determinant of a
C square matrix
C Passed parameters:
C A = the matrix
C N = dimension of the square matrix
C A modification of a code originally written by Ashwith J. Rego, available from
C http://www.dreamincode.net/code/snippet1273.htm
C Modified by Syeilendra Pramuditya, available from http://wp.me/p61TQ-zb
C Last modified on January 13, 2011
c
c------------------------------------------------------------------------------
c
      function GETDET (A,N)
c
      IMPLICIT REAL*8 (A-H,O-Z)
c
c Input/Output
      REAL*8 A(N,N)
      real*8 GETDET
c
c Local
      REAL*8 ELEM(N,N)
      REAL*8 M, TEMP
      INTEGER I, J, K, L
      LOGICAL DETEXISTS
c
c------------------------------------------------------------------------------
c
      DO I=1,N
        DO J=1,N
          ELEM(I,J)=A(I,J)
        END DO
      END DO
      DETEXISTS = .TRUE.
      L = 1
!CONVERT TO UPPER TRIANGULAR FORM
      DO K = 1, N-1
        IF (DABS(ELEM(K,K)).LE.1.0D-20) THEN
          DETEXISTS = .FALSE.
          DO I = K+1, N
            IF (ELEM(I,K).NE.0.0) THEN
              DO J = 1, N
                TEMP = ELEM(I,J)
                ELEM(I,J)= ELEM(K,J)
                ELEM(K,J) = TEMP
              END DO
              DETEXISTS = .TRUE.
              L=-L
              EXIT
            END IF
          END DO
          IF (DETEXISTS .EQV. .FALSE.) THEN
            GETDET = 0
            RETURN
          END IF
        END IF
        DO J = K+1, N
          M = ELEM(J,K)/ELEM(K,K)
          DO I = K+1, N
            ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
          END DO
        END DO
      END DO
!CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
      GETDET = L
      DO I = 1, N
        GETDET = GETDET * ELEM(I,I)
      END DO
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      READBINPOL2C.FOR    (28 August 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Read computed center of mass from cube file.
c
c------------------------------------------------------------------------------
c
      subroutine readbinpol2c (infile,mem,lmem,noc,mcen,vcb,J0,
     %  xc,yc,zc,mc,xct,yct,zct,T1b,T2,TP)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      character*80 infile(6)
      real*8 xc(noc),yc(noc),zc(noc),mc(noc)
      real*8 mcen,vc,J0(3,3),mcenb
      integer noc
      real*8 xct,yct,zct
      real*8 T1b(3),T2(3),TP(3)
c
c Local
      integer lim(2,1000),nsub
      integer j,itmp
      character*15000 string
      character*80 c80
      integer lineno
      character*7 c7
      logical test
      real*8 xctb,yctb,zctb,vcb,cubx,cuby,cubz
c
c------------------------------------------------------------------------------
c
      inquire (file=infile(6), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(6)
        stop
      end if
      open (10, file=infile(6), status='old')
      lineno = 0
      write (*,'(a13,1x,a)') 'Reading file:',infile(6)
c
      open (10, file=infile(6), status='old', access='sequential')
c
  40  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 40
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 40
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) noc
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) mcenb
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) vcb
      c80 = string(lim(1,4):lim(2,4))
      read (c80,*,err=661,end=667) xctb
      c80 = string(lim(1,5):lim(2,5))
      read (c80,*,err=661,end=667) yctb
      c80 = string(lim(1,6):lim(2,6))
      read (c80,*,err=661,end=667) zctb
c
  50  read (10,'(a15000)',end=667) string
      lineno = lineno + 1
      if (string(1:1).eq.')') goto 50
      call mio_spl (15000,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 50
      c80 = string(lim(1,1):lim(2,1))
      read (c80,*,err=661,end=667) cubx
      c80 = string(lim(1,2):lim(2,2))
      read (c80,*,err=661,end=667) cuby
      c80 = string(lim(1,3):lim(2,3))
      read (c80,*,err=661,end=667) cubz
c
  80    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 80
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 80
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xct
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yct
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zct
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mcen
        c80 = string(lim(1,5):lim(2,5))
        read (c80,*,err=661,end=667) vc
c
      do j = 1, 3
  60    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 60
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 60
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) J0(j,1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) J0(j,2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) J0(j,3)
      end do
c
  61    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 61
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 61
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T1b(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T1b(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T1b(3)
  62    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 62
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 62
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) T2(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) T2(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) T2(3)
  63    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 63
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 63
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) TP(1)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) TP(2)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) TP(3)
c
      do j = 1, noc
  70    read (10,'(a15000)',end=667) string
        lineno = lineno + 1
        if (string(1:1).eq.')') goto 70
        call mio_spl (15000,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 70
        c80 = string(lim(1,1):lim(2,1))
        read (c80,*,err=661,end=667) xc(j)
        c80 = string(lim(1,2):lim(2,2))
        read (c80,*,err=661,end=667) yc(j)
        c80 = string(lim(1,3):lim(2,3))
        read (c80,*,err=661,end=667) zc(j)
        c80 = string(lim(1,4):lim(2,4))
        read (c80,*,err=661,end=667) mc(j)
      end do
      close (10)
c
 667  continue
      if(lineno.eq.0) then
        write(6,'(/,3(1x),2(a,1x),/)') 'Error:',
     %    'No center mass point was found.'
        stop
      end if      
c
c termina a execução do programa
      write (6,'(a)') 'Process completed successfully!'
c
c------------------------------------------------------------------------------
c
      return
c
c------------------------------------------------------------------------------
c
c Error reading from the input file containing integration parameters
 661  write (c7,'(i7)') lineno
      call mio_err (6,mem(1),lmem(1),mem(6),lmem(6),c7,7,
     %  mem(7),lmem(7))
c
c------------------------------------------------------------------------------
c
      end
c
c------------------------------------------------------------------------------
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER3.FOR    (27 Ago 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer3 (nov,nop,xv,yv,zv,noe,k,d,
     %  fmt,fvt,xc,yc,zc,mem,lmem,l0,cubx,cuby,cubz,mc)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),d,mc(nocen)
      real*8 cubx,cuby,cubz
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
c      integer i,j,l,ilay,novl,nopl
c      integer noel,kl
c      real*8 del,inc,nx,ny,nz
c      real*8 xl(nov),yl(nov),zl(nov)
c      real*8 norm(nopmax,3),wfac(nopmax)
c      real*8 T0,T1(3),T2(3),TP(3)
c      real*8 fvt,fmt,vol,vol1,vol2
c      real*8 J0(3,3)
      real*8 fvt,fmt
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor,mass,vol
      integer astflag1
c      integer astflag1,astflag2
c      real*8 xvc(8),yvc(8),zvc(8)
c      real*8 cubx,cuby,cubz
c      real*8 xl2(nov),yl2(nov),zl2(nov)
c      dimension noel(nopmax)
c      dimension kl(nopmax,noed)
c
c------------------------------------------------------------------------------
c
c      del = 1.d0 / lay
c      vt = 0.d0
c      mt = 0.d0
      l0 = 0
c
c      T1t(1) = 0.0
c      T1t(2) = 0.0
c      T1t(3) = 0.0
c      T2t(1) = 0.0
c      T2t(2) = 0.0
c      T2t(3) = 0.0
c      TPt(1) = 0.0
c      TPt(2) = 0.0
c      TPt(3) = 0.0
c central point (origin)
c      xl(1) = 0.0
c      yl(1) = 0.0
c      zl(1) = 0.0
c
c      cubx = 1.d0 / cubx0
c      cuby = 1.d0 / cuby0
c      cubz = 1.d0 / cubz0
c      cubx = cubx0
c      cuby = cuby0
c      cubz = cubz0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
c      do ilay = 1, lay
c        inc = del * ilay
c
c        do i = 1, nov
c
c compute new vertices
c          xl(i) = xv(i) * inc
c          yl(i) = yv(i) * inc
c          zl(i) = zv(i) * inc
c          if (ilay.gt.1) then
c            xl2(i) = xv(i) * del * (ilay-1)
c            yl2(i) = yv(i) * del * (ilay-1)
c            zl2(i) = zv(i) * del * (ilay-1)
c          end if
c
c          if (i.eq.1) then
c            xmin=xl(i)
c            xmax=xl(i)
c            ymin=yl(i)
c            ymax=yl(i)
c            zmin=zl(i)
c            zmax=zl(i)
c          end if
c
c          if(xl(i).lt.xmin) xmin = xl(i)
c          if(xl(i).gt.xmax) xmax = xl(i)
c          if(yl(i).lt.ymin) ymin = yl(i)
c          if(yl(i).gt.ymax) ymax = yl(i)
c          if(zl(i).lt.zmin) zmin = zl(i)
c          if(zl(i).gt.zmax) zmax = zl(i)
c        end do
c
        fvt = 0.d0
        fmt = 0.d0
c
        zstep = zmin
        do while (zstep.lt.zmax)
          ystep = ymin
          do while (ystep.lt.ymax)
            xstep = xmin
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c
              if (astflag1.eq.1) then
c                if (ilay.le.1) then
c                  astflag2 = 0
c                else
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c                  call tetrah (nov,nop,xl2,yl2,zl2,noe,k,
c     %              xcoor,ycoor,zcoor,astflag2)
c                end if
c                if (astflag2.eq.0) then
c compute cubic vertices
c                xvc(1) = xstep
c                yvc(1) = ystep
c                zvc(1) = zstep
c                xvc(2) = xstep + cubx
c                yvc(2) = ystep
c                zvc(2) = zstep
c                xvc(3) = xstep + cubx
c                yvc(3) = ystep + cuby
c                zvc(3) = zstep
c                xvc(4) = xstep
c                yvc(4) = ystep + cuby
c                zvc(4) = zstep
c                xvc(5) = xstep
c                yvc(5) = ystep
c                zvc(5) = zstep + cubz
c                xvc(6) = xstep + cubx
c                yvc(6) = ystep
c                zvc(6) = zstep + cubz
c                xvc(7) = xstep + cubx
c                yvc(7) = ystep + cuby
c                zvc(7) = zstep + cubz
c                xvc(8) = xstep
c                yvc(8) = ystep + cuby
c                zvc(8) = zstep + cubz
c
c compute cubic's faces
c                noel(1) = 4
c                kl(1,1) = 4
c                kl(1,2) = 3
c                kl(1,3) = 2
c                kl(1,4) = 1
c
c                noel(2) = 4
c                kl(2,1) = 1
c                kl(2,2) = 2
c                kl(2,3) = 6
c                kl(2,4) = 5
c
c                noel(3) = 4
c                kl(3,1) = 2
c                kl(3,2) = 3
c                kl(3,3) = 7
c                kl(3,4) = 6
c
c                noel(4) = 4
c                kl(4,1) = 3
c                kl(4,2) = 4
c                kl(4,3) = 8
c                kl(4,4) = 7
c
c                noel(5) = 4
c                kl(5,1) = 4
c                kl(5,2) = 1
c                kl(5,3) = 5
c                kl(5,4) = 8
c
c                noel(6) = 4
c                kl(6,1) = 5
c                kl(6,2) = 6
c                kl(6,3) = 7
c                kl(6,4) = 8
c
c compute volumens
c                nopl = 6
c compute normal faces and w vector
c                do j = 1, nopl
c                  call compnormalface (nov,nop,j,xvc,yvc,zvc,kl,
c     %              nx,ny,nz)
c                  norm(j,1) = nx
c                  norm(j,2) = ny
c                  norm(j,3) = nz
c                  wfac(j) = - norm(j,1) * xvc(kl(j,1))
c     %                      - norm(j,2) * yvc(kl(j,1))
c     %                      - norm(j,3) * zvc(kl(j,1))
c                enddo
c compute volumen of a cube
c                call compVolumeIntegrals (nov,nop,nopl,xvc,yvc,zvc,
c     %            noel,kl,norm,wfac,T0,T1,T2,TP)
c                fvt = fvt + T0
c                T1t(1) = T1t(1) + T1(1)
c                T1t(2) = T1t(2) + T1(2)
c                T1t(3) = T1t(3) + T1(3)
c                T2t(1) = T2t(1) + T2(1)
c                T2t(2) = T2t(2) + T2(2)
c                T2t(3) = T2t(3) + T2(3)
c                TPt(1) = TPt(1) + TP(1)
c                TPt(2) = TPt(2) + TP(2)
c                TPt(3) = TPt(3) + TP(3)
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
c                m(l0) = d(ilay) * T0
c                fmt = fmt + m(l0)
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
c compute center of masses of a cube
c                call compcenpolyhedron (d(ilay),m(l0),T0,T1,
c     %            T2,TP,xc(l0),yc(l0),zc(l0),J0)
c                Jt(1,1) = Jt(1,1) + J0(1,1)
c                Jt(1,2) = Jt(1,2) + J0(1,2)
c                Jt(1,3) = Jt(1,3) + J0(1,3)
c                Jt(2,1) = Jt(2,1) + J0(2,1)
c                Jt(2,2) = Jt(2,2) + J0(2,2)
c                Jt(2,3) = Jt(2,3) + J0(2,3)
c                Jt(3,1) = Jt(3,1) + J0(3,1)
c                Jt(3,2) = Jt(3,2) + J0(3,2)
c                Jt(3,3) = Jt(3,3) + J0(3,3)
c              end if
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c        vt = vt + fvt
c        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(24),
     %    lmem(24),' ',1,mem(2),lmem(2))
c
c compute mascons' masses
c        do i = 1, l0
c          m(i) = mcen / l0
c        enddo
c
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MASC_LAYER4B.FOR    (28 Ago 2020)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Calculates the pieces of masses of each pyramidal frustum layer of a polyhedron.
c Also calculates the volumens, centroids and inertia tensor of each pyramidal
c frustum layer.
c

c------------------------------------------------------------------------------
c
      subroutine masc_layer4b (nov,nop,xv,yv,zv,noe,k,d,
     %  fmt,fvt,xc,yc,zc,mem,lmem,l0,cubx,cuby,cubz,mc)
c
      implicit none
      include 'polyhedron.inc'
c
c Input/Output
      integer lmem(NMESS)
      character*80 mem(NMESS)
      integer nov,nop,noe,k,l0
      real*8 xv(nov),yv(nov),zv(nov)
      real*8 xc(nocen),yc(nocen),zc(nocen),d,mc(nocen)
      real*8 cubx,cuby,cubz,eixa,eixb,eixc
      dimension noe(nopmax)
      dimension k(nopmax,noed)
c
c Local
      integer i
      real*8 fvt,fmt
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
      real*8 xstep,ystep,zstep,xcoor,ycoor,zcoor,mass,vol
      integer astflag1
c      real*8 a2,b2,c2,test
c
c------------------------------------------------------------------------------
c
      l0 = 0
c
c compute xmin, xmax, ymin, ymax, zmin e zmax
      do i = 1, nov
        if (i.eq.1) then
          xmin=xv(i)
          xmax=xv(i)
          ymin=yv(i)
          ymax=yv(i)
          zmin=zv(i)
          zmax=zv(i)
        end if
c
        if(xv(i).lt.xmin) xmin = xv(i)
        if(xv(i).gt.xmax) xmax = xv(i)
        if(yv(i).lt.ymin) ymin = yv(i)
        if(yv(i).gt.ymax) ymax = yv(i)
        if(zv(i).lt.zmin) zmin = zv(i)
        if(zv(i).gt.zmax) zmax = zv(i)
      end do
c
        fvt = 0.d0
        fmt = 0.d0
c
c      if (eixa.eq.0.d0) then
c        a2 = 0.d0
c      else
c        a2 = 1.d0/eixa
c      end if
c      if (eixb.eq.0.d0) then
c        b2 = 0.d0
c      else
c        b2 = 1.d0/eixb
c      end if
c      if (eixc.eq.0.d0) then
c        c2 = 0.d0
c      else
c        c2 = 1.d0/eixc
c      end if
c
c 1º Octante
c        write(6,'(a21)') 'Starting 1 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c 2º Octante
c        write(6,'(a21)') 'Starting 2 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep + cubz
        end do
c
c 3º Octante
c        write(6,'(a21)') 'Starting 3 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep + cubz
        end do
c
c 4º Octante
c        write(6,'(a21)') 'Starting 4 Octant...'
        zstep = 0.d0
        do while (zstep.lt.zmax)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep + 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep + cubz
        end do
c
c------------------------------------------------------------------------
c
c 5º Octante
c        write(6,'(a21)') 'Starting 5 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep - cubz
        end do
c
c 6º Octante
c        write(6,'(a21)') 'Starting 6 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.lt.ymax)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep + 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep + cuby
          end do
          zstep = zstep - cubz
        end do
c
c 7º Octante
c        write(6,'(a21)') 'Starting 7 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.gt.xmin)
              xcoor = xstep - 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep - cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep - cubz
        end do
c
c 8º Octante
c        write(6,'(a21)') 'Starting 8 Octant...'
        zstep = 0.d0
        do while (zstep.gt.zmin)
          ystep = 0.d0
          do while (ystep.gt.ymin)
            xstep = 0.d0
            do while (xstep.lt.xmax)
              xcoor = xstep + 0.5d0 * cubx
              ycoor = ystep - 0.5d0 * cuby
              zcoor = zstep - 0.5d0 * cubz
c check if the point is inside (iflag=1) or outside (iflag=0) the polyhedron using the tetrahedron method
c              astflag1 = 1
c              test = xcoor*xcoor*a2+ycoor*ycoor*b2+zcoor*zcoor*c2 -1.d0
c              if (test.gt.0.d0.or.test.eq.-1.d0) then
              call tetrah (nov,nop,xv,yv,zv,noe,k,xcoor,ycoor,zcoor,
     %          astflag1)
c              endif
c
              if (astflag1.eq.1) then
c compute pieces of masses
                l0 = l0 + 1
                vol = cubx * cuby * cubz
                mass = d * vol
                mc(l0) = mass
                fmt = fmt + mass
                fvt = fvt + vol
                xc(l0) = xcoor
                yc(l0) = ycoor
                zc(l0) = zcoor
              end if
              xstep = xstep + cubx
            end do
            ystep = ystep - cuby
          end do
          zstep = zstep - cubz
        end do
c
c        vt = vt + fvt
c        mt = mt + fmt
        if (l0.gt.nocen) call mio_err (6,mem(1),lmem(1),mem(24),
     %    lmem(24),' ',1,mem(2),lmem(2))
c
c compute mascons' masses
c        do i = 1, l0
c          m(i) = mcen / l0
c        enddo
c
c      end do
c
c------------------------------------------------------------------------------
c
      return
      end
c
