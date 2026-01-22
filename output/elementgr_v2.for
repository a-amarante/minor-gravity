c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      ELEMENTGR.FOR    (FEG   26 April 2016)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c******************************************************************************
c CHANGES IN VERSION 2 (17 Dec 2025)
c******************************************************************************
c
c Author: Andre Amarante (A. Amarante) - andre.amarante@unesp.br
c
c Makes output files containing potentials from data created by Polyhedron.
c
c As alterações foram feitas através dos marcadores ##Am##.
c
c------------------------------------------------------------------------------
c
      implicit none
c ##A1##
      include 'polyhedron.inc'
c
      integer nchar,allflag,firstflag,ninfile,nel,iel(30)
      integer*4 lprevious,nline,j,m,nout,nlout,line_grid
      integer*8 line_num,nc,ne
      integer*4 nstored,ltmp,nstoredmax
      logical test
      character*80 mem(NMESS),cc
      character*104 c(CMAX)
      integer distyle,lenhead,lmem(NMESS)
      character*250 string,fout,header,infile(50)
      integer nmaster,nopen,nwait,nsub,lim(2,100)
      real*8 xi,xf,yi,yf,zi,zf,omg,xp,yp,zp
      character*1 check,c1
      character*250 mheader
      integer lenmhead
      integer itmp,i,k,l,n,precision,lenin
      character*8 master_id(CMAX)
      integer unit(CMAX),master_unit(CMAX)
      character*125 str1,str2
      character*6 fin
      real*8 el(30,CMAX)
      real*8 mio_c2re, mio_c2fl
      character*5 fomg
c ##A1##
c
c      integer itmp,i,j,k,l,iback(CMAX),precision,lenin
c      integer nmaster,nopen,nwait,nbig,nsml,nbod,nsub,lim(2,100)
c      integer year,month,timestyle,line_num,lenhead,lmem(NMESS)
c      integer nchar,algor,centre,allflag,firstflag,ninfile,nel,iel(22)
c      integer nbod1,nbig1,unit(CMAX),code(CMAX),master_unit(CMAX)
c      real*8 time,teval,t0,t1,tprevious,rmax,rcen,rfac,rhocgs,temp
c      real*8 mcen,jcen(3),el(22,CMAX),s(3),is(CMAX),ns(CMAX),a(CMAX)
c      real*8 mio_c2re, mio_c2fl,fr,theta,phi,fv,vtheta,vphi,gm
c      real*8 x(3,CMAX),v(3,CMAX),xh(3,CMAX),vh(3,CMAX),m(CMAX)
c      logical test
c      character*250 string,fout,header,infile(50)
c      character*80 mem(NMESS),cc,c(CMAX)
c      character*8 master_id(CMAX),id(CMAX)
c      character*5 fin
c      character*1 check,style,type,c1
c      character*2 c2
c
c------------------------------------------------------------------------------
c
      allflag = 0
c ##A2##
      lprevious = 0
c ##A2##
c      rhocgs = AU * AU * AU * K2 / MSUN
c
c Read in output messages
c      inquire (file='message.in', exist=test)
c      if (.not.test) then
c        write (*,'(/,2a)') ' ERROR: This file is needed to continue: ',
c     %    ' message.in'
c        stop
c      end if
c      open (14, file='message.in', status='old')
c  10  continue
c        read (14,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
c      goto 10
c  20  close (14)
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
          infile(1)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        end if
      end do
      close(15)
c
c Read in output messages
      inquire (file=infile(1), exist=test)
      if (.not.test) then
        write (*,'(/,3a)') ' ERROR: This file is needed to start',
     %    ' the integration:  ',infile(1)
        stop
      end if
      open (16, file=infile(1), status='old')
  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
      goto 10
  20  close (16)
c
      do k = 1, 50
        do j = 1, 250
          infile(k)(j:j) = ' '
        end do
      end do
c
c Open file containing parameters for this programme
      inquire (file='element.in', exist=test)
      if (test) then
        open (10, file='element.in', status='old')
      else
c ##A3##
        call mio_err (6,mem(1),lmem(1),mem(2),lmem(2),' ',1,
     %    'element.in',10)
c ##A3##
      end if
c
c Read number of input files
  30  read (10,'(a250)') string
      if (string(1:1).eq.')') goto 30
      call mio_spl (250,string,nsub,lim)
      if (lim(1,1).eq.-1) goto 30
      read (string(lim(1,nsub):lim(2,nsub)),*) ninfile
c
c Make sure all the input files exist
      do j = 1, ninfile
  40    read (10,'(a250)') string
        if (string(1:1).eq.')') goto 40
        call mio_spl (250,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 40
        infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        inquire (file=infile(j), exist=test)
c ##A6##
        if (.not.test) call mio_err (6,mem(1),lmem(1),mem(2),
     %    lmem(2),' ',1,infile(j),250)
c ##A6##
      end do
c
c What type elements does the user want?
c      centre = 0
c  45  read (10,'(a250)') string
c      if (string(1:1).eq.')') goto 45
c      call mio_spl (250,string,nsub,lim)
c      c2 = string(lim(1,nsub):(lim(1,nsub)+1))
c      if (c2.eq.'ce'.or.c2.eq.'CE'.or.c2.eq.'Ce') then
c        centre = 0
c      else if (c2.eq.'ba'.or.c2.eq.'BA'.or.c2.eq.'Ba') then
c        centre = 1
c      else if (c2.eq.'ja'.or.c2.eq.'JA'.or.c2.eq.'Ja') then
c        centre = 2
c      else
c        call mio_err (6,mem(81),lmem(81),mem(107),lmem(107),' ',1,
c     %    '       Check element.in',23)
c      end if
c
c Read parameters used by this programme
c ##A7##
      do j = 1, 6
  52    read (10,'(a250)') string
        if (string(1:1).eq.')') goto 52
        call mio_spl (250,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 52
        if (j.eq.1) read (string(lim(1,nsub):lim(2,nsub)),*) xi
        if (j.eq.2) read (string(lim(1,nsub):lim(2,nsub)),*) xf
        if (j.eq.3) read (string(lim(1,nsub):lim(2,nsub)),*) yi
        if (j.eq.4) read (string(lim(1,nsub):lim(2,nsub)),*) yf
        if (j.eq.5) read (string(lim(1,nsub):lim(2,nsub)),*) zi
        if (j.eq.6) read (string(lim(1,nsub):lim(2,nsub)),*) zf
      end do
      distyle = 1
      do j = 1, 3
  50    read (10,'(a250)') string
        if (string(1:1).eq.')') goto 50
        call mio_spl (250,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 50
        c1 = string(lim(1,nsub):lim(2,nsub))
        if (j.eq.1) read (string(lim(1,nsub):lim(2,nsub)),*) nout
        if (j.eq.2) call m_format (string,distyle,nel,iel,fout,header,
     %    lenhead,mheader,lenmhead)
        if (j.eq.3) read (string(lim(1,nsub):lim(2,nsub)),*) nlout
      end do
      nout = abs(nout)
      nlout = abs(nlout)
c ##A7##
c
c Read in the names of the objects for which orbital elements are required
      nopen = 0
      nwait = 0
      nmaster = 0
  60  continue
        read (10,'(a250)',end=70) string
        call mio_spl (250,string,nsub,lim)
        if (string(1:1).eq.')'.or.lim(1,1).eq.-1) goto 60
c
c Either open an aei file for this object or put it on the waiting list
        nmaster = nmaster + 1
        itmp = min(7,lim(2,1)-lim(1,1))
        master_id(nmaster)='        '
        master_id(nmaster)(1:itmp+1) = string(lim(1,1):lim(1,1)+itmp)
        if (nopen.lt.NFILES) then
          nopen = nopen + 1
          master_unit(nmaster) = 10 + nopen
c ##A15##
          call mio_aei (master_id(nmaster),'.dat',master_unit(nmaster),
     %      header,lenhead,mem,lmem,mheader,lenmhead)
c ##A15##
        else
          nwait = nwait + 1
          master_unit(nmaster) = -2
        end if
      goto 60
c
  70  continue
c If no objects are listed in ELEMENT.IN assume that all objects are required
      if (nopen.eq.0) allflag = 1
      close (10)
c
c ##A20##
      if (allflag.eq.1) then
        nmaster = nmaster + 1
        write(str1,*) nmaster
        str2 = 'v'
        string = str2//str1
        call mio_spl (250,string,nsub,lim)
        string(1:(lim(2,1)-lim(1,1)+1)+(lim(2,2)-lim(1,2)+1))=
     %  string(lim(1,1):lim(2,1))//string(lim(1,2):lim(2,2))
        call mio_spl (250,string,nsub,lim)
        itmp = min(7,lim(2,1)-lim(1,1))
        master_id(nmaster)='        '
        master_id(nmaster)(1:itmp+1) = string(lim(1,1):lim(1,1)+itmp)
        if (nopen.lt.NFILES) then
          nopen = nopen + 1
          master_unit(nmaster) = 10 + nopen
          call mio_aei (master_id(nmaster),'.dat',master_unit(nmaster),
     %      header,lenhead,mem,lmem,mheader,lenmhead)
        else
          nwait = nwait + 1
          master_unit(nmaster) = -2
        end if
      end if
      do j = 1, nmaster
        unit(j) = master_unit(j)
      end do
      master_unit(nmaster+1) = -3
      unit(nmaster+1) = master_unit(nmaster+1)
      line_grid = 0
      m = 0
      n = 1
c ##A20##
c
c------------------------------------------------------------------------------
c
c  LOOP  OVER  EACH  INPUT  FILE  CONTAINING  INTEGRATION  DATA
c
c  90  continue
      firstflag = 0
c ##A21##
      do i = 1, ninfile
        line_num = 0
        write (*,'(a13,1x,a)') 'Reading file:',infile(i)
        open (10, file=infile(i), status='old', access='sequential')
c
c Loop over each time slice
 100    continue
        nc = 0
        line_num = line_num + 1
        read (10,'(a1)',end=900,eor=900,err=666,advance='no') check
c ##A21##
c        line_num = line_num - 1
c        backspace 10
c
c Check if this is an old style input file
c        if (ichar(check).eq.12.and.(style.eq.'0'.or.style.eq.'1'.or.
c     %    style.eq.'2'.or.style.eq.'3'.or.style.eq.'4')) then
c          write (*,'(/,2a)') ' ERROR: This is an old style data file',
c     %      '        Try running m_elem5.for instead.'
c          stop
c        end if
        if (ichar(check).ne.12) goto 666
c
c------------------------------------------------------------------------------
c
c  IF  SPECIAL  INPUT,  READ  TIME,  PARAMETERS,  NAMES,  MASSES  ETC.
c
c        if (type.eq.'a') then
          line_num = line_num + 1
c ##A22##
          read (10,'(i1)',end=666,eor=666,err=666,advance='no')
     %    precision
c          backspace 10
          if (precision.eq.1) then
            nchar = 3
          else if (precision.eq.2) then
                 nchar = 5
               else if (precision.eq.3) then
                      nchar = 8
                    else
                      goto 666
                    end if
          fomg(1:5) = '(12a)'
          lenin = 12
          read (10,fomg,end=666,eor=666,err=666,advance='no')
     %    (cc(nc:nc),nc=1,lenin)
          line_num = line_num + nc - 1
c ##A22##
c
c Decompress the time, number of objects, central mass and J components etc.
c          time = mio_c2fl (cc(1:8))
c          nbig = int(.5d0 + mio_c2re(cc(9:16), 0.d0, 11239424.d0, 3))
c          nsml = int(.5d0 + mio_c2re(cc(12:19),0.d0, 11239424.d0, 3))
c          mcen = mio_c2fl (cc(15:22))
c          jcen(1) = mio_c2fl (cc(23:30))
c          jcen(2) = mio_c2fl (cc(31:38))
c          jcen(3) = mio_c2fl (cc(39:46))
c          rcen = mio_c2fl (cc(47:54))
c          rmax = mio_c2fl (cc(55:62))
c          rfac = log10 (rmax / rcen)
c ##A23##
          nstored = int(.5d0 + mio_c2re(cc(1:4), 0.d0,
     %    2517630976.d0, 4))
          omg = mio_c2fl (cc(5:12),8)
c ##A23##
c Read in strings containing compressed data for each object
c          do j = 1, nbig + nsml
c            line_num = line_num + 1
c            read (10,'(a)',err=666) c(j)(1:51)
c          end do
c
c Create input format list
c ##A28##
c          if (precision.eq.1) nchar = 3
c          if (precision.eq.2) nchar = 5
c          if (precision.eq.3) nchar = 8
          lenin = 25+10*nchar-1
c          fin(1:6) = '(a000)'
c          if (lenin.lt.100) write (fin(4:5),'(i2)') lenin
c          if (lenin.ge.100) write (fin(3:5),'(i3)') lenin
c          lenin = 4  +  13 * nchar
          nstoredmax = 0
          do j = 1, nstored
c            line_num = line_num + 1
c            read (10,fin,end=911,err=666) c(j)(5:lenin)
            do nc = 1, lenin
              read (10,'(a1)',end=667,eor=667,err=666,advance='no')
     %        c(j)(nc:nc)
              if (ichar(c(j)(nc:nc)).eq.12) then
                line_num = line_num - 1
                backspace 10
                do ne = 1, line_num + nc
                  read (10,'(a1)',advance='no') c1
                end do
                line_num = line_num + 1
                goto 666
              end if
            end do
            nstoredmax = nstoredmax + 1
            line_num = line_num + nc - 1
          end do
c ##A28##
c
c For each object decompress its name, code number, mass, spin and density
c          do j = 1, nbig + nsml
c            k = int(.5d0 + mio_c2re(c(j)(1:8),0.d0,11239424.d0,3))
c            id(k) = c(j)(4:11)
c            el(18,k) = mio_c2fl (c(j)(12:19))
c            s(1) = mio_c2fl (c(j)(20:27))
c            s(2) = mio_c2fl (c(j)(28:35))
c            s(3) = mio_c2fl (c(j)(36:43))
c            el(21,k) = mio_c2fl (c(j)(44:51))
c ##A29##
          nline = 0
          do j = 1, nstoredmax
            xp = mio_c2fl (c(j)(1:8),8)
            yp = mio_c2fl (c(j)(9:16),8)
            zp = mio_c2fl (c(j)(17:24),8)
            if(.not.((xp.ge.xi.and.xp.le.xf).and.
     %      (yp.ge.yi.and.yp.le.yf).and.
     %      (zp.ge.zi.and.zp.le.zf))) then
              goto 1001
            end if
            nline = nline + 1
            el(1,nline) = xp
            el(2,nline) = yp
            el(3,nline) = zp
            el(4,nline) = mio_c2fl (c(j)(25:25+nchar-1),nchar)
            el(4,nline) = -el(4,nline)
            el(5,nline) = mio_c2fl (c(j)(25+nchar:25+2*nchar-1),nchar)
            el(6,nline) = mio_c2fl (c(j)(25+2*nchar:25+3*nchar-1),nchar)
            el(7,nline) = mio_c2fl (c(j)(25+3*nchar:25+4*nchar-1),nchar)
            el(8,nline) = mio_c2fl (c(j)(25+4*nchar:25+5*nchar-1),nchar)
            el(9,nline) = mio_c2fl (c(j)(25+5*nchar:25+6*nchar-1),nchar)
            el(10,nline) = mio_c2fl (c(j)(25+6*nchar:25+7*nchar-1),
     %      nchar)
            el(11,nline) = mio_c2fl (c(j)(25+7*nchar:25+8*nchar-1),
     %      nchar)
            el(12,nline) = mio_c2fl (c(j)(25+8*nchar:25+9*nchar-1),
     %      nchar)
            el(13,nline) = mio_c2fl (c(j)(25+9*nchar:25+10*nchar-1),
     %      nchar)
            el(14,nline) = -omg*omg*(el(1,nline)*el(1,nline)+
     %      el(2,nline)*el(2,nline))/2.0d0 + el(4,nline)
            el(15,nline) = dsqrt(el(5,nline)*el(5,nline)+
     %      el(6,nline)*el(6,nline)+el(7,nline)*el(7,nline))
            el(16,nline)=acos(el(5,nline)/el(15,nline))/DR
            el(17,nline)=acos(el(6,nline)/el(15,nline))/DR
            el(18,nline)=acos(el(7,nline)/el(15,nline))/DR
            el(19,nline)=el(8,nline)+el(9,nline)+el(10,nline)
            el(20,nline) = omg*omg*(el(1,nline)*el(1,nline)+
     %      el(2,nline)*el(2,nline))/2.0d0 + el(4,nline)
            el(21,nline)=-omg*(-el(2,nline)*el(5,nline)+
     %      el(1,nline)*el(6,nline))
            el(22,nline) = omg*omg*el(1,nline)+el(5,nline)
            el(23,nline) = omg*omg*el(2,nline)+el(6,nline)
            el(24,nline) = el(7,nline)
            el(25,nline) = dsqrt(el(22,nline)*el(22,nline)+
     %      el(23,nline)*el(23,nline)+el(24,nline)*el(24,nline))
            el(26,nline) = dsqrt(xp*xp+yp*yp+zp*zp)
            el(27,nline) = -omg*omg*(el(1,nline)*el(1,nline)+
     %      el(2,nline)*el(2,nline))/2.0d0
            el(28,nline) = omg*omg*el(1,nline)
            el(29,nline) = omg*omg*el(2,nline)
            el(30,nline) = dsqrt(el(28,nline)*el(28,nline)+
     %      el(29,nline)*el(29,nline))
 1001       continue
          end do
          if (nline.eq.0) goto 100
          line_grid = line_grid + nline
c ##A29##
c
c Calculate spin rate and longitude & inclination of spin vector
c            temp = sqrt(s(1)*s(1) + s(2)*s(2) + s(3)*s(3))
c            if (temp.gt.0) then
c              call mce_spin (1.d0,el(18,k)*K2,temp*K2,el(21,k)*
c     %            rhocgs,el(20,k))
c              temp = s(3) / temp
c              if (abs(temp).lt.1) then
c                is(k) = acos (temp)
c                ns(k) = atan2 (s(1), -s(2))
c              else
c                if (temp.gt.0) is(k) = 0.d0
c                if (temp.lt.0) is(k) = PI
c                ns(k) = 0.d0
c              end if
c            else
c              el(20,k) = 0.d0
c              is(k) = 0.d0
c              ns(k) = 0.d0
c            end if
c
c Find the object on the master list
c            unit(k) = 0
c            do l = 1, nmaster
c              if (id(k).eq.master_id(l)) unit(k) = master_unit(l)
c              unit(l) = master_unit(l)
c            end do
c
c If object is not on the master list, add it to the list now
c            if (unit(k).eq.0) then
c              nmaster = nmaster + 1
c              master_id(nmaster) = id(k)
c
c Either open an aei file for this object or put it on the waiting list
c              if (allflag.eq.1) then
c                if (nopen.lt.NFILES) then
c                  nopen = nopen + 1
c                  master_unit(nmaster) = 10 + nopen
c                  call mio_aei (master_id(nmaster),'.aei',
c     %              master_unit(nmaster),header,lenhead,mem,lmem)
c                else
c                  nwait = nwait + 1
c                  master_unit(nmaster) = -2
c                end if
c              else
c                master_unit(nmaster) = -1
c              end if
c              unit(k) = master_unit(nmaster)
c            end if
c          end do
c
c------------------------------------------------------------------------------
c
c  IF  NORMAL  INPUT,  READ  COMPRESSED  ORBITAL  VARIABLES  FOR  ALL  OBJECTS
c
c        else if (type.eq.'b') then
c          line_num = line_num + 1
c          read (10,'(3x,a14)',err=666) cc(1:14)
c
c Decompress the time and the number of objects
c          time = mio_c2fl (cc(1:8))
c          nbig = int(.5d0 + mio_c2re(cc(9:16),  0.d0, 11239424.d0, 3))
c          nsml = int(.5d0 + mio_c2re(cc(12:19), 0.d0, 11239424.d0, 3))
c          nbod = nbig + nsml
c          if (firstflag.eq.0) t0 = time
c
c Read in strings containing compressed data for each object
c          do j = 1, nbod
c            line_num = line_num + 1
c            read (10,fin,err=666) c(j)(1:lenin)
c          end do
c
c Look for objects for which orbital elements are required
c          m(1) = mcen * K2
c          do j = 1, nbod
c            code(j) = int(.5d0 + mio_c2re(c(j)(1:8), 0.d0,
c     %        11239424.d0, 3))
c            if (code(j).gt.CMAX) then
c              write (*,'(/,2a)') mem(81)(1:lmem(81)),
c     %          mem(90)(1:lmem(90))
c              stop
c            end if
c
c Decompress orbital variables for each object
c            l = j + 1
c            m(l) = el(18,code(j)) * K2
c            fr     = mio_c2re (c(j)(4:11), 0.d0, rfac,  nchar)
c            theta  = mio_c2re (c(j)(4+  nchar:11+  nchar), 0.d0, PI,
c     %               nchar)
c            phi    = mio_c2re (c(j)(4+2*nchar:11+2*nchar), 0.d0, TWOPI,
c     %               nchar)
c            fv     = mio_c2re (c(j)(4+3*nchar:11+3*nchar), 0.d0, 1.d0,
c     %               nchar)
c            vtheta = mio_c2re (c(j)(4+4*nchar:11+4*nchar), 0.d0, PI,
c     %               nchar)
c            vphi   = mio_c2re (c(j)(4+5*nchar:11+5*nchar), 0.d0, TWOPI,
c     %               nchar)
c            call mco_ov2x (rcen,rmax,m(1),m(l),fr,theta,phi,fv,
c     %        vtheta,vphi,x(1,l),x(2,l),x(3,l),v(1,l),v(2,l),v(3,l))
c            el(16,code(j)) = sqrt(x(1,l)*x(1,l) + x(2,l)*x(2,l)
c     %                     + x(3,l)*x(3,l))
c          end do
c
c Convert to barycentric, Jacobi or close-binary coordinates if desired
c          nbod1 = nbod + 1
c          nbig1 = nbig + 1
c          call mco_iden (jcen,nbod1,nbig1,temp,m,x,v,xh,vh)
c          if (centre.eq.1) call mco_h2b (jcen,nbod1,nbig1,temp,m,xh,vh,
c     %      x,v)
c          if (centre.eq.2) call mco_h2j (jcen,nbod1,nbig1,temp,m,xh,vh,
c     %      x,v)
c          if (centre.eq.0.and.algor.eq.11) call mco_h2cb (jcen,nbod1,
c     %      nbig1,temp,m,xh,vh,x,v)
c
c Put Cartesian coordinates into element arrays
c          do j = 1, nbod
c            k = code(j)
c            l = j + 1
c            el(10,k) = x(1,l)
c            el(11,k) = x(2,l)
c            el(12,k) = x(3,l)
c            el(13,k) = v(1,l)
c            el(14,k) = v(2,l)
c            el(15,k) = v(3,l)
c
c Convert to Keplerian orbital elements
c            gm = (mcen + el(18,k)) * K2
c            call mco_x2el (gm,el(10,k),el(11,k),el(12,k),el(13,k),
c     %        el(14,k),el(15,k),el(8,k),el(2,k),el(3,k),el(7,k),
c     %        el(5,k),el(6,k))
c            el(1,k) = el(8,k) / (1.d0 - el(2,k))
c            el(9,k) = el(1,k) * (1.d0 + el(2,k))
c            el(4,k) = mod(el(7,k) - el(5,k) + TWOPI, TWOPI)
c Calculate true anomaly
c            if (el(2,k).eq.0) then
c              el(17,k) = el(6,k)
c            else
c              temp = (el(8,k)*(1.d0 + el(2,k))/el(16,k) - 1.d0) /el(2,k)
c              temp = sign (min(abs(temp), 1.d0), temp)
c              el(17,k) = acos(temp)
c              if (sin(el(6,k)).lt.0) el(17,k) = TWOPI - el(17,k)
c            end if
c Calculate obliquity
c            el(19,k) = acos (cos(el(3,k))*cos(is(k))
c     %        + sin(el(3,k))*sin(is(k))*cos(ns(k) - el(5,k)))
c
c Convert angular elements from radians to degrees
c            do l = 3, 7
c              el(l,k) = mod(el(l,k) / DR, 360.d0)
c            end do
c            el(17,k) = el(17,k) / DR
c            el(19,k) = el(19,k) / DR
c          end do
c
c Convert time to desired format
c          if (timestyle.eq.0) t1 = time
c          if (timestyle.eq.1) call mio_jd_y (time,year,month,t1)
c          if (timestyle.eq.2) t1 = time - t0
c          if (timestyle.eq.3) t1 = (time - t0) / 365.25d0
c
c If output is required at this epoch, write elements to appropriate files
c          if (firstflag.eq.0.or.abs(time-tprevious).ge.teval) then
c            firstflag = 1
c            tprevious = time
c
c Write required elements to the appropriate aei file
c ##A30##
            ltmp = line_grid - nline
            do k = 1, nline
              ltmp = ltmp + 1
              if (firstflag.eq.0.or.(ltmp-lprevious).ge.nout) then
                if (firstflag.eq.0) lprevious = ltmp-1
                if (firstflag.eq.1) lprevious = ltmp
                firstflag = 1
                do while (unit(n).eq.-1)
                  n = n + 1
                end do
                if (unit(n).eq.-2) then
c Close aei files
                  do j = 1, nopen
                    close (10+j)
                  end do
                  nopen = 0
c
c If some objects remain on waiting list, read through input files again
                  if (nwait.gt.0) then
                    do j = 1, nmaster
                      if (master_unit(j).ge.10) master_unit(j) = -1
                      if (master_unit(j).eq.-2.and.nopen.lt.NFILES) then
                        nopen = nopen + 1
                        nwait = nwait - 1
                        master_unit(j) = 10 + nopen
                        call mio_aei (master_id(j),'.dat',
     %  master_unit(j),header,lenhead,mem,lmem,mheader,lenmhead)
                        if (master_unit(j).eq.-1) then
                          nopen = nopen - 1
                          n = n + 1
                        end if
                      end if
                    end do
                  end if
                  do j = 1, nmaster
                    unit(j) = master_unit(j)
                  end do
                end if
                if (unit(n).eq.-3) then
c Close aei files
                  do j = 1, nopen
                    close (10+j)
                  end do
                  nopen = 0
 1234             nmaster = nmaster + 1
                  write(str1,*) nmaster
                  str2 = 'v'
                  string = str2//str1
                  call mio_spl (250,string,nsub,lim)
                  string(1:(lim(2,1)-lim(1,1)+1)+(lim(2,2)-lim(1,2)+1))=
     %            string(lim(1,1):lim(2,1))//string(lim(1,2):lim(2,2))
                  call mio_spl (250,string,nsub,lim)
                  itmp = min(7,lim(2,1)-lim(1,1))
                  master_id(nmaster)='        '
                  master_id(nmaster)(1:itmp+1) =
     %            string(lim(1,1):lim(1,1)+itmp)
                  if (nopen.lt.NFILES) then
                    nopen = nopen + 1
                    master_unit(nmaster) = 10 + nopen
                    call mio_aei (master_id(nmaster),'.dat',
     %  master_unit(nmaster),header,lenhead,mem,lmem,mheader,lenmhead)
                    if (master_unit(nmaster).eq.-1) then
                     nopen = nopen - 1
                     n = n + 1
                     goto 1234
                    end if
                  else
                    nwait = nwait + 1
                    master_unit(nmaster) = -2
                  end if
                  do j = 1, nmaster
                    unit(j) = master_unit(j)
                  end do
                  master_unit(nmaster+1) = -3
                  unit(nmaster+1) = master_unit(nmaster+1)
                end if
                m = m + 1
                if (unit(n).ge.10) then
                  write (unit(n),fout) (el(iel(l),k),l=1,nel)
                end if
                if (nlout.ne.0) then
                  if (mod(m,nlout).eq.0) then
                    n = n + 1
                    m = 0
                  end if
                end if
              end if
            end do
c ##A30##
c          end if
c
c------------------------------------------------------------------------------
c
c  IF  TYPE  IS  NOT  'a'  OR  'b',  THE  INPUT  FILE  IS  CORRUPTED
c
c        else
c          goto 666
c        end if
c
c Move on to the next time slice
        goto 100
c
c If input file is corrupted, try to continue from next uncorrupted time slice
c ##A18##
 667    continue
        nc = nc - 1
 666    continue
        write (*,'(2a,/,a,i11)') mem(46)(1:lmem(46)),
     %    infile(i)(1:60),mem(49)(1:lmem(49)),line_num + nc + 1
        if (nc.ne.0) line_num = line_num - 1
        line_num = line_num + nc
c ##A18##
        c1 = ' '
        do while (ichar(c1).ne.12)
          line_num = line_num + 1
          read (10,'(a1)',end=900,eor=900,advance='no') c1
        end do
        line_num = line_num - 1
        backspace 10
        do ne = 1, line_num
          read (10,'(a1)',advance='no') c1
        end do
c
c ##ERROR Chambers##
c para procurar o próximo pedaço de tempo não corrompido no mesmo arquivo de saída e NÃO no próximo arquivo de saída
        goto 100
c ##ERROR Chambers##
c
c Move on to the next file containing integration data
 900    continue
        close (10)
      end do
c
c Close aei files
c      do j = 1, nopen
c        close (10+j)
c      end do
c      nopen = 0
c
c If some objects remain on waiting list, read through input files again
c      if (nwait.gt.0) then
c        do j = 1, nmaster
c          if (master_unit(j).ge.10) master_unit(j) = -1
c          if (master_unit(j).eq.-2.and.nopen.lt.NFILES) then
c            nopen = nopen + 1
c            nwait = nwait - 1
c            master_unit(j) = 10 + nopen
c            call mio_aei (master_id(j),'.aei',master_unit(j),header,
c     %        lenhead,mem,lmem)
c          end if
c        end do
c        goto 90
c      end if
c
c------------------------------------------------------------------------------
c
c  CREATE  A  SUMMARY  OF  FINAL  MASSES  AND  ELEMENTS
c
c      open (10, file='element.out', status='unknown')
c      rewind 10
c
c      if (timestyle.eq.0.or.timestyle.eq.2) then
c        write (10,'(/,a,f18.5,/)') ' Time (days): ',t1
c      else if (timestyle.eq.1) then
c        write (10,'(/,a,i10,1x,i2,1x,f8.5,/)') ' Date: ',year,month,t1
c      else if (timestyle.eq.3) then
c        write (10,'(/,a,f18.7,/)') ' Time (years): ',t1
c      end if
c      write (10,'(2a,/)') '              a        e       i      mass',
c     %  '    Rot/day  Obl'
c
c Sort surviving objects in order of increasing semi-major axis
c      do j = 1, nbod
c        k = code(j)
c        a(j) = el(1,k)
c      end do
c      call mxx_sort (nbod,a,iback)
c
c Write values of a, e, i and m for surviving objects in an output file
c      do j = 1, nbod
c        k = code(iback(j))
c        write (10,213) id(k),el(1,k),el(2,k),el(3,k),el(18,k),el(20,k),
c     %      el(19,k)
c      end do
c
c------------------------------------------------------------------------------
c
c Format statements
c 213  format (1x,a8,1x,f8.4,1x,f7.5,1x,f7.3,1p,e11.4,0p,1x,f6.3,1x,f6.2)
c
c ##A31##
      if(line_grid.eq.0) then
        write(6,'(/,3(1x),2(a,1x),/)') 'No point potential was found.',
     %  'See the file element.in'
      end if      
c termina a execução do programa
      write (*,'(a)') mem(48)(1:lmem(48))
      stop
c ##A31##
c
      end
c
c ##A4##
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
c ##A5##
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
      integer len,nsub,delimit(2,100)
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
c ##A8##
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      M_FORMAT.FOR    (ErikSoft   31 January 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 14 July 2013)
c
c Makes an output format list and file header for the orbital-element files
c created by M_ELEM3.FOR
c Also identifies which orbital elements will be output for each object.
c
c------------------------------------------------------------------------------
c
      subroutine m_format (string,distyle,nel,iel,fout,header,lenhead,
     %  mheader,lenmhead)
c
c ##A9##
      implicit none
c
c Input/Output
      integer distyle,nel,iel(30),lenhead
      character*250 string,header,fout
      character*250 mheader
      integer lenmhead
c
c Local
c ##Ubuntu-22-LTS##
c      integer i,j,pos,nsub,lim(2,20),formflag,lenfout,f1,f2,itmp
      integer i,j,pos,nsub,lim(2,100),formflag,lenfout,f1,f2,itmp
c ##Ubuntu-22-LTS##
      character*3 elcode(30)
      character*4 elhead(30)
      integer k,posc
      character*80 c80
c ##A9##
c
c------------------------------------------------------------------------------
c
c ##A10##
      data elcode/ 'x','y','z','v','vx','vy','vz','vxx','vyy','vzz',
     %  'vxy','vxz','vyz','j','a','phi','the','ome','lap','ene','pow',
     %  'vtx','vty','vtz','at','r','c','cax','cay','ca'/
      data elhead/ '  X ','  Y ','  Z ','  V ',' Vx ',' Vy ',' Vz ',
     %  ' Vxx',' Vyy',' Vzz',' Vxy',' Vxz',' Vyz','  J ','Acel',
     %  ' Phi','Thet','Omeg','Lapl','Enrg','PowG',' Vtx',' Vty',
     %  ' Vtz','AccT','Radi','VC','VCx','VCy','CA'/
c ##A10##
c
c Initialize header to a blank string
      do i = 1, 250
        header(i:i) = ' '
      end do
c
c Create part of the format list and header for the required time style
c      if (timestyle.eq.0.or.timestyle.eq.2) then
c ##A11##
        fout(1:3) = '(1x'
        lenfout = 3
        header(1:1) = ' '
        lenhead = 1
c ##A11##
c      else if (timestyle.eq.1) then
c        fout(1:21) = '(1x,i10,1x,i2,1x,f8.5'
c        lenfout = 21
c        header(1:23) = '    Year/Month/Day     '
c        lenhead = 23
c      else if (timestyle.eq.3) then
c        fout(1:9) = '(1x,f18.7'
c        lenfout = 9
c        header(1:19) = '    Time (years)   '
c        lenhead = 19
c      end if
c
c ##A12##
      mheader(1:42) = 'distance (vertice unit) / time (seconds) /'
      mheader(43:58) = ' angle (degrees)'
      lenmhead = 58
c ##A12##
c
c Identify the required elements
      call mio_spl (250,string,nsub,lim)
c ##A13##
      do i = 1, nsub
        iel(i) = 0
      end do
      do i = 1, nsub
        do j = 1, 30
          do k = lim(1,i), lim(2,i)
            if (.not.(LGT(string(k:k),'9'))) then
              posc = k
              goto 10
            end if
          end do
  10      do k = 1, 80
            c80(k:k) = ' '
          end do
          c80(1:(posc-lim(1,i))) = string(lim(1,i):(posc-1))
          if (LLE(c80(1:3),elcode(j)(1:3)).and.
     %      LGE(c80(1:3),elcode(j)(1:3))) iel(i) = j
        end do
        if (iel(i).eq.0) then
          write(6,'(/,3(1x),a,1x,i2,1x,2(a,1x),/)') 'Error code',i,
     %  'selected is invalid.','See the file element.in'
          stop
        end if
c ##A13##
      end do
      nel = nsub
c
c For each element, see whether normal or exponential notation is required
      do i = 1, nsub
        formflag = 0
c ##A14##
        do k = lim(1,i), lim(2,i)
          if (.not.(LGT(string(k:k),'9'))) then
            posc = k
            goto 20
          end if
        end do
  20    do j = posc, lim(2,i)
          if (formflag.eq.0) pos = j
          if (string(j:j).eq.'.') formflag = 1
          if (string(j:j).eq.'e') formflag = 2
        end do
c
c Create the rest of the format list and header
        if (formflag.eq.1) then
          read (string(posc:pos-1),*) f1
          read (string(pos+1:lim(2,i)),*) f2
          write (fout(lenfout+1:lenfout+10),'(a10)') ',1x,f  .  '
          write (fout(lenfout+6:lenfout+7),'(i2)') f1
          write (fout(lenfout+9:lenfout+10),'(i2)') f2
          lenfout = lenfout + 10
        else if (formflag.eq.2) then
          read (string(posc:pos-1),*) f1
          read (string(pos+1:lim(2,i)-1),*) f2
          write (fout(lenfout+1:lenfout+16),'(a16)') ',1x,1p,e  .  ,0p'
          write (fout(lenfout+9:lenfout+10),'(i2)') f1
          write (fout(lenfout+12:lenfout+13),'(i2)') f2
c ##A14##
          lenfout = lenfout + 16
        end if
        itmp = (f1 - 4) / 2
        header(lenhead+itmp+2:lenhead+itmp+5) = elhead(iel(i))
        lenhead = lenhead + f1 + 1
      end do
c
      lenfout = lenfout + 1
      fout(lenfout:lenfout) = ')'
c
c------------------------------------------------------------------------------
c
      return
      end
c
c ##A16##
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_AEI.FOR    (ErikSoft   31 January 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 19 July 2013)
c
c Creates a filename and opens a file to store aei information for an object.
c The filename is based on the name of the object.
c
c------------------------------------------------------------------------------
c
      subroutine mio_aei (id,extn,unitnum,header,lenhead,mem,lmem,
     %  mheader,lenmhead)
c
      implicit none
c ##A17##
      include 'polyhedron.inc'
c
c Input/Output
      integer unitnum,lenhead,lmem(NMESS)
      character*4 extn
      character*8 id
      character*250 header
      character*80 mem(NMESS)
      character*250 mheader
      integer lenmhead
c ##A17##
c
c Local
c ##Ubuntu-22-LTS##
c      integer j,k,itmp,nsub,lim(2,4)
      integer j,k,itmp,nsub,lim(2,100)
c ##Ubuntu-22-LTS##
      logical test
      character*1 bad(5)
      character*250 filename
c
c------------------------------------------------------------------------------
c
      data bad/ '*', '/', '.', ':', '&'/
c
c Create a filename based on the object's name
      call mio_spl (8,id,nsub,lim)
      itmp = min(7,lim(2,1)-lim(1,1))
      filename(1:itmp+1) = id(1:itmp+1)
      filename(itmp+2:itmp+5) = extn
      do j = itmp + 6, 250
        filename(j:j) = ' '
      end do
c
c Check for inappropriate characters in the filename
      do j = 1, itmp + 1
        do k = 1, 5
          if (filename(j:j).eq.bad(k)) filename(j:j) = '_'
        end do
      end do
c
c If the file exists already, give a warning and don't overwrite it
      inquire (file=filename, exist=test)
      if (test) then
c ##A18##
        write (*,'(/,3a)') mem(46)(1:lmem(46)),mem(5)(1:lmem(5)),
     %    filename(1:80)
c ##A18##
        unitnum = -1
      else
        open (unitnum, file=filename, status='new')
c ##A19##
c        write (unitnum, '(/,30x,a,//,a)') mheader(1:lenmhead),
c     %  header(1:lenhead)
        write (unitnum, '(a1,/,a1,30x,a,/,a1,/,a1,a)') '#','#',
     %    mheader(1:lenmhead),'#','#',header(1:lenhead)
c ##A19##
      end if
c
c------------------------------------------------------------------------------
c
      return
      end
c
c ##A24##
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_C2FL.FOR    (ErikSoft   5 June 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers (adapted by Andre Amarante - 16 July 2013)
c
c CHARACTER*8 ASCII string into a REAL*8 variable.
c
c N.B. X will lie in the range -1.e112 < X < 1.e112
c ===
c
c------------------------------------------------------------------------------
c
      function mio_c2fl (c,nchar)
c
      implicit none
c
c Input/Output
      real*8 mio_c2fl
      character*8 c
c ##A25##
      integer nchar
c ##A25##
c
c Local
      real*8 x,mio_c2re
      integer ex
c
c------------------------------------------------------------------------------
c
c ##A26##
      x = mio_c2re (c(1:nchar), 0.d0, 1.d0, nchar-1)
      x = x * 2.d0 - 1.d0
      ex = mod(ichar(c(nchar:nchar)) + 256, 256) - 32 - 112
c ##A26##
      mio_c2fl = x * (10.d0**dble(ex))
c
c------------------------------------------------------------------------------
c
      return
      end
c
c ##A27##
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      MIO_C2RE.FOR    (ErikSoft   5 June 2001)
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c Author: John E. Chambers
c
c Converts an ASCII string into a REAL*8 variable X, where XMIN <= X < XMAX,
c using the new format compression:
c
c X is assumed to be made up of NCHAR base-224 digits, each one represented
c by a character in the ASCII string. Each digit is given by the ASCII
c number of the character minus 32.
c The first 32 ASCII characters (CTRL characters) are avoided, because they
c cause problems when using some operating systems.
c
c------------------------------------------------------------------------------
c
      function mio_c2re (c,xmin,xmax,nchar)
c
      implicit none
c
c Input/output
      integer nchar
      real*8 xmin,xmax,mio_c2re
      character*8 c
c
c Local
      integer j
      real*8 y
c
c------------------------------------------------------------------------------
c
      y = 0
      do j = nchar, 1, -1
        y = (y + dble(mod(ichar(c(j:j)) + 256, 256) - 32)) / 224.d0
      end do
c
      mio_c2re = xmin  +  y * (xmax - xmin)
c
c------------------------------------------------------------------------------
c
      return
      end
c
