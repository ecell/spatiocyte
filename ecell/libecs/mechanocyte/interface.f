      MODULE TEMP
      USE iolibsw

      real(c_double),bind(C)::tnext      
      real(8) PIP2m,PIP3m,PIP3a,PI3Km,PTENm
      real(8) tdump(100),idump
      real(c_double),bind(C)::tstop
      integer,dimension(1)::isuff
      character(len=20)::fnum

      CONTAINS

      subroutine extractnumeric(ipf,idot) bind(C)
      use iso_c_binding
      implicit none
      integer(C_INT)::idot
      character(kind=c_char, len=1), dimension(20), intent(in)::ipf
      character (len=20) :: newipf
      newipf = " "
      loop_string: do i=1, 20
      if ( ipf(i) == c_null_char ) then
         exit loop_string
      else
         newipf (i:i) = ipf(i)
      end if
      end do loop_string
      fnum=newipf(idot+1:idot+3)
      isuff=(iachar(fnum(1:1))-48)*100+
     1      (iachar(fnum(2:2))-48)*10+
     2       iachar(fnum(3:3))-48
      print *,'isuff=',isuff
      end subroutine extractnumeric

      subroutine openfile(ipf,idiv,idot) bind(C)
      use iso_c_binding
      implicit none
      integer(C_INT)::idiv,idot
      character(kind=c_char, len=1), dimension(20), intent(in)::ipf
      character (len=20) :: newipf
      character (len=20)::datafile
      integer :: i,iframe
      newipf = " "
      loop_string: do i=1, 20
      if ( ipf(i) == c_null_char ) then
         exit loop_string
      else
         newipf (i:i) = ipf(i)
      end if
      end do loop_string
      open(unit=idiv,file=newipf,status='old')
      call iordfile(iframe,idiv)
      close(idiv)
      datafile=newipf
      datafile(idot+1:idot+3)='dat'
      open(31,file=datafile,status='new')
      end subroutine openfile
  
      subroutine wrfile(idiv,isve,ipf,idot,delt,logInt) bind(C)
      use iso_c_binding
      implicit none
      integer(C_INT),bind(C)::idiv,isve,idot
      character(kind=c_char, len=1), dimension(20), intent(in)::ipf
      character(len=20) :: opf
      real(c_double)::delt,logInt
      open(unit=idiv,file='test_out')
      call iowrfile(0,idiv)
      close(idiv)
      print *,fnum
      print *,isuff
      if (isve.eq.1) then
            isuff=isuff+1
            print *,isuff
            opf = " "
            loop_string: do i=1, 20
            if ( ipf(i) == c_null_char ) then
               exit loop_string
            else
               opf(i:i) = ipf(i)
            end if
            end do loop_string
            write(fnum,90)isuff
            opf(idot+1:idot+3)=fnum
            open(unit=21,file=opf,status='unknown')
            call iowrfile(0,21)
            close(21)      
!
           if (time_.ge.tstop*0.999) then
               print *,'time,tstop',time_,tstop
               close(31)
               stop
           else
               idump=idump+1
               if (logInt.eq.0) then
                  tnext=tdump(idump)
               else
                  tnext=tnext+logInt
               endif
           endif
        endif
   90 format(i3.3)
      end subroutine wrfile
  

      subroutine timedump(delt,logInt) bind(C)
      use iso_c_binding
      implicit none
      real(c_double)::delt,logInt
      real(8) ndump
!--if delt is zero look for time of dumps in file tdump
      if (delt.eq.zero) then
         open(44,file='tdump',status='old')
         read(44,*) ndump
         do i=1,ndump
            read(44,*) tdump(i)
         enddo
         tstop=tdump(ndump)
c--look for current dump position
         if (time_.gt.tstop) then
            print *,'time>tstop',time_,tstop
            stop
         else
            do i=1,ndump
               if (tdump(i).gt.time_) then
                  idump=i
                  goto 20
               endif
            enddo
         endif
   20    print *,'ndump,idump,tstop',ndump,idump,tstop
         tnext=tdump(idump)
      else
         tstop=time_+delt
         tnext=time_+logInt
      endif
      print *,'tnext,tstop',tnext,tstop
      end subroutine timedump

      subroutine initsvec() bind(C)
      USE iolibsw
      use iso_c_binding
      integer in
      real(8) x,y
      do isn=1,ns
         do lvn=1,3
            x = (hvec(1,lvn,isn))
            y = (hvec(2,lvn,isn))
            !if(x.lt.5d-6) svec(12,lvn,isn) = 1.444605d12 !PIP2m
            if(y.lt.2.5d-6) then
               svec(12,lvn,isn) = 1.277d13 !PIP2m
               svec(9,lvn,isn) = 0.627d12 !PTENm
            endif
            if(x.lt.-6d-6) then
               svec(8,lvn,isn) = 3.8546d13 !PI3Km
            endif
         enddo
      enddo
      end subroutine initsvec

      subroutine initarea(area) bind(C)
      USE iolibsw
      use iso_c_binding
      implicit none
      real(c_double)::area
      real(8) surfintv,surfintd,surfinte
      real(8),dimension(3,NSM)::cnode=1d0
      call gosurfintn(cnode,surfintv,surfintd,surfinte)
      area = surfintv+surfintd+surfinte
      end subroutine initarea

      subroutine getSurfaceMolecules(kom, cnt)
      USE iolibsw
      real(8) cnt,surfintv,surfintd,surfinte
      real(8) cnode(3,NSM)

      do i=1,ns
         cnode(1,i)=svec(kom,1,i)
         cnode(2,i)=svec(kom,2,i)
         cnode(3,i)=svec(kom,3,i)
      enddo
      call gosurfintn(cnode,surfintv,surfintd,surfinte)
      cnt = surfintv+surfintd+surfinte
      end subroutine getSurfaceMolecules

      subroutine clchm(area) bind(C)
      use iolibsw
      use iso_c_binding
      real(c_double)::area
      integer in,iPIP2,iPIP3,iPIP3a,iPI3K,iPTEN
      real(8) PIPtc,PTENtc,PI3Ktc
      real(8) PIP2,PI3K,PTEN,maxConc,new,val
      real(8) k1,k2,k3,k4,k5,k6,k7,k8,k9,k10
      iPIP2 = 12
      iPIP3 = 11
      iPIP3a = 10
      iPTEN = 9
      iPI3K = 8

      k1 = 4d-2
      k2 = 2d-14
      k3 = 2d-13
      k4 = 4d-14
      k5 = 2d-14
      k6 = 5d-14
      k7 = 5d-14
      k8 = 0.09
      k9 = 0.02
      k10 = 0.0001
      maxConc = 8.7816d13
      PIPtc = 1.089171d13
      PTENtc = 7.62666d12
      PI3Ktc = 1.44957d13

      call getSurfaceMolecules(iPIP2,val)
      PIP2 = val
      call getSurfaceMolecules(iPIP3,val)
      PIP2 = PIP2 + val
      call getSurfaceMolecules(iPIP3a,val)
      PIP2 = PIP2 + val
      PIP2 = PIPtc-PIP2/area 

      call getSurfaceMolecules(iPI3K,PI3K)
      PI3K = PI3Ktc-PI3K/area

      call getSurfaceMolecules(iPTEN,PTEN)
      PTEN = PTENtc-PTEN/area
c      print *,"PTEN:",int(PTEN*area),"PI3K:",int(PI3K*area),
c     1        "PIP2:",int(PIP2*area)

      do isn=1,ns
         do lvn=1,3
            PIP2m = svec(iPIP2,lvn,isn)
            PIP3m = svec(iPIP3,lvn,isn)
            PIP3a = svec(iPIP3a,lvn,isn)
            PTENm = svec(iPTEN,lvn,isn)
            PI3Km = svec(iPI3K,lvn,isn)

            sdot(iPIP2,lvn,isn) = k1*PIP2-k5*PIP2m*PI3Km+k6*PIP3m*PTENm+
     1                            k7*PIP3a*PTENm-k10*PIP2m
            sdkr(iPIP2,lvn,isn) = -k5*PI3Km-k10
            new = PIP2m + sdot(iPIP2,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP2,lvn,isn) = 0
              sdkr(iPIP2,lvn,isn) = 0
            endif

            sdot(iPTEN,lvn,isn) = k2*PTEN*PIP2m-k8*PTENm
            sdkr(iPTEN,lvn,isn) = -k8
            new = PTENm + sdot(iPTEN,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPTEN,lvn,isn) = 0
              sdkr(iPTEN,lvn,isn) = 0
            endif

            sdot(iPIP3a,lvn,isn) = -k3*PIP3a*PI3K+k4*PIP3m*PIP3m-
     1                             k7*PIP3a*PTENm-k9*PIP3a
            sdkr(iPIP3a,lvn,isn) = -k3*PI3k-k7*PTENm-k9
            new = PIP3a + sdot(iPIP3a,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP3a,lvn,isn) = 0
              sdkr(iPIP3a,lvn,isn) = 0
            endif

            sdot(iPIP3,lvn,isn) = k3*PIP3a*PI3K-k4*PIP3m*PIP3m+
     1                             k5*PIP2m*PI3Km-k6*PIP3m*PTENm-
     2                             k9*PIP3m
            sdkr(iPIP3,lvn,isn) = -2*k4*PIP3m-k6*PTENm-k9
            new = PIP3m + sdot(iPIP3,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPIP3,lvn,isn) = 0
              sdkr(iPIP3,lvn,isn) = 0
            endif

            sdot(iPI3K,lvn,isn) = k3*PIP3a*PI3K-k9*PI3Km
            sdkr(iPI3K,lvn,isn) = -k9
            new = PI3Km + sdot(iPI3K,lvn,isn)*tstp
            if (new.gt.maxConc) then
              sdot(iPI3K,lvn,isn) = 0
              sdkr(iPI3K,lvn,isn) = 0
            endif
         enddo
      enddo
      end subroutine clchm

      subroutine chksurmol() bind(C)
      use iso_c_binding
      iPIP2 = 12
      iPIP3 = 11
      iPIP3a = 10
      iPTEN = 9
      iPI3K = 8
      call getSurfaceMolecules(12,PIP2m)
      call getSurfaceMolecules(11,PIP3m)
      call getSurfaceMolecules(10,PIP3a)
      call getSurfaceMolecules(9,PTENm)
      call getSurfaceMolecules(8,PI3Km)
      print *,"time:",time_
      print *,"PIP2m:",int(PIP2m),"PI3Km:",int(PI3Km),"PIP3m:",
     1           int(PIP3m),"PIP3a:",int(PIP3a),"PTENm:",int(PTENm)
      end subroutine chksurmol


      END MODULE TEMP
