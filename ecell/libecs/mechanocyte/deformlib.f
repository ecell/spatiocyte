      MODULE deformlib
      USE iolibsw
      USE MOLIBSW
      USE AVLIBSW
      USE CMLIBSW
      
      real(c_double),bind(C)::tnext 
      real(8) tdump(400)
      integer idump
      real(c_double),bind(C)::tstop
      integer,dimension(1)::isuff
      character(len=20)::fnum

      CONTAINS
      
      subroutine testcut() bind(C)
      use iso_c_binding
      character ipf*20, opf*20, datafile*20
      character fnum*3, suff*3, pref*3
      real(8) epsl, avdt, cmdt(12),volint, tdump(400),hmax,rmean
      real(8) surfintv,surfintd,surfinte   
      real(8) delt,tsurf,r0,a0,aratio,yq,xq,al,vmax,cnwmin
      real(8) alphamin,ratiomn,rmax,rmin,dxy,dz,rxy,alpha
      real(8) const1,const2,const3,const4,const5,const6
      real(8),dimension(3,NSM)::cnode=1d0

!      goto 88
      call getarg(1,ipf)
      idot=index(ipf,'.')
      if (idot.eq.0) then
         write(*,*)'invalid inputfile name'
         stop
      endif
!--extract numeric field
      fnum=ipf(idot+1:idot+3)
      isuff=(iachar(fnum(1:1))-48)*100+
     1      (iachar(fnum(2:2))-48)*10+
     2       iachar(fnum(3:3))-48
      print *,'isuff=',isuff
!
!--open,read and close the input file
!
      
      open(unit=11,file=ipf,status='old')
      call iordfile(iframe,11)
      close(11)
!--setup data log file
      datafile=ipf
      datafile(idot+1:idot+3)='dat'
      open(31,file=datafile,status='new')
!
!--figure out time dumps
!
      delt=dscp(3,idfrm)
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
         tnext=time_+delt
      endif
      print *,'tnext,tstop',tnext,tstop
!  
!   88 continue
      idphi=1 !phase rub
      idvis=1 !viscosity
      idvfx=1 !fixed velocity at boundaries
      idgam=1 !surface tension
      idpsi=1 !network contractility
      idsfr=1 !surface forces
      idbfr=0 !body forces
      iddrg=0 !drag at boundaries
      idtrc=0 !boundary tractions
      idhyc=0 !boundary solvent permeability
      call clphi
      call clvis
      call clpsi
      call govolint(cnode,volint)
      call gosurfintn(cnode,surfintv,surfintd,surfinte)
      tsurf=surfintv+surfintd+surfinte
      !print *,volint
      const1=volint*3d0
      const2=4d0*3.141592d0
      const3=1d0/3d0
      r0=(const1/const2)**const3
      !r0=(volint*3.0/(4.0*3.141592))**(1./3.)
      !print *,r0
      a0=4.0*3.141592d0*r0**2
      aratio=tsurf/a0
      !print *,tsurf,a0,aratio
      call clgam(aratio)
!
      do iq=1,nq
         vfixv(1,iq)=1d0
         vfixv(2,iq)=1d0
         vfixv(3,iq)=1d0
      enddo
      if (time_.le.1d-5) then
         do is=1,ns
            svec(4,1,is)=max(svec(4,1,is),1d-4)
         enddo
         do il=1,nl
            cxprm(4,il)=1d0
            cxval(4,il)=2d-1
!           if (time.lt.5d1) then
!              cxval(4,il)=2d-1*(1d0+time/1d2)
!           else
!              cxval(4,il)=3d-1
!           endif
            iq=iqol(il)
            yq=0d0
            xq=0d0
            do isn=1,4
               is=isoq(isn,iq)
               xq=xq+hvec(1,1,is)
               yq=yq+hvec(2,1,is)
               !print *,'xq ',il,xq
            enddo
            al=atan(abs(yq)/(abs(xq)+1d-20))*180d0/3.141592d0
            !print *,'al ',al
            if (xq.lt.0d0.or.al.gt.45d0) then
!           if (xq.lt.0d0) then
               cxprm(4,il)=0d0
!              if (al.lt.75d0) then
!                 vfixv(1,iq)=0d0
!                 vfixv(2,iq)=0d0
!              endif
            endif
            evec(4,il)=cxprm(4,il)
         enddo
      else
         do il=1,nl
            cxprm(4,il)=evec(4,il)
            cxval(4,il)=2d-1
!           if (time.lt.5d1) then
!              cxval(4,il)=2d-1*(1d0+time/1d2)
!           else
!              cxval(4,il)=3d-1
!           endif
         enddo
      endif
      !do i=1,nl
          !print *,'cxprm ',cxprm(4,i)
      !enddo
!
      call vvec1
      call clsfr
!     call clbfr
!
      epsl=1d-4
      idebug=1
      kount=0
  10  call modriver(icyc,epsl,idebug)
      kount=kount+1
      open(unit=12,file='test_out')
      call iowrfile(0,12)
      close(12)
      if (icyc.gt.10) goto 10
!     stop
!
      idebug=1
  100 continue
         vmax=maxval(sqrt(hvec(7,1:3,1:ns)**2+hvec(8,1:3,1:ns)**2+
     1                                        hvec(9,1:3,1:ns)**2))
         cnwmin=minval(svec(1,1:3,1:ns))
         print *,'vmax=',vmax,'  cnwmin=',cnwmin
         call avtstep(avdt)
         if (avdt.lt.1d-9) then
            print *,'small avdt!',avdt
            stop
         endif
         call clchm
         call vvec1
         call BC
         call cmstpsiz(cmdt,1d-2)
!        tstp=cmdt(2)
!        tstp=min(cmdt(2),cmdt(4))
!        tstp=min(cmdt(2),2d-2)
         print *,'chemistry step,advection step',tstp,avdt
!        print *,'cmdt(1:4)',cmdt(1:4)
!        tstp=min(avdt,tstp)
         tstp=avdt
         if (tstp.ge.(tnext-time_)) then
            tstp=tnext-time_
            isve=1
         else
            isve=0
         endif
         print *,'advt,tstp',avdt,tstp
         !print *,'aaaaaaaaaaaaaaaaaaaaaaaaaaaaa',svec(4,1,1:ns)
         !call dfdriver('lagrangian')
         !do i =1,ns
         !  print *,i,'  ',svec(4,1,i)
         !enddo
         print *,'min 4',minval(svec(4,1,1:ns))
         !call avgridmo('lagrangian')
         !call avfield(1,'nw')
!        call avfield(2,'ca')
!        cnode=svec(1,:,:)
!        call govolint(cnode,volint)
!        print *,'total network',volint
!
         call clphi
         call clvis
         call clpsi
         call vvec1
         call clsfr
!        call clbfr
         call clgam(aratio)
         do 
           call modriver(icyc,epsl,idebug)
           if (icyc.lt.10) exit
           open(unit=12,file='test_out')
           call iowrfile(0,12)
           close(12)
         enddo
!        open(unit=12,file='test_out')
!        call iowrfile(0,12)
!        close(12)
!        stop
!        call govolint(cnode,volint)
!        print *,'min zmid/ztop',minval(hvec(3,2,1:ns)/hvec(3,3,1:ns))
!        print *,'zmin',minval(hvec(3,3,1:ns))
         alphamin=vbig
         ratiomn=vbig
         rmean=0d0
         rmax=0d0
         rmin=vbig
         do il=1,nl
           is=isol(1,il)
           dxy=sqrt((hvec(1,2,is)-hvec(1,1,is))**2+
     1              (hvec(2,2,is)-hvec(2,1,is))**2)
           dz=hvec(3,2,is)-hvec(3,1,is)
           rxy=sqrt(hvec(1,1,is)**2+hvec(2,1,is)**2)
           rmean=rmean+rxy
           rmax=max(rmax,rxy)
           rmin=min(rmin,rxy)
           alpha=atan(dz/dxy)
           alphamin=min(alpha,alphamin)
           ratiomn=min(hvec(3,2,is)/hvec(3,3,is),ratiomn)
         enddo
         rmean=rmean/nl
         !print *,'real alphamin: ',alphamin
         print *,'alphamin',alphamin*180d0/3.141592d0,'ratiomn=',ratiomn
         time_=time_+tstp
         open(unit=12,file='test_out')
         call iowrfile(0,12)
         close(12)
         hmax=maxval(hvec(3,3,1:ns))
         call govolint(cnode,volint)
         call gosurfintn(cnode,surfintv,surfintd,surfinte)
         tsurf=surfintv+surfintd+surfinte
         const4 = volint*3d0
         const5 = 4d0*3.141592d0
         const6 = 1d0/3d0
         !print *,tsurf,const4,const5,const6
         r0=(const4/const5)**(const6)
         !r0=(volint*3.0/(4.0*3.141592))**(1./3.)
         a0=4d0*3.141592d0*r0**2
         !print *,volint
         !print *,r0
         !print *,a0
         aratio=tsurf/a0
         print *,'i,volume, aratio, time_',i,volint,aratio, time_
         write(31,*) time_, rmin,rmax,hmax,volint,aratio
         if (isve.eq.1) then
            isuff=isuff+1
            opf=ipf
            write(fnum,90)isuff
            print *,'isuff=',isuff
            opf(idot+1:idot+3)=fnum
            open(unit=21,file=opf,status='unknown')
            call iowrfile(0,21)
            close(21)
!
!--check if we have reached the end
!
            if (time_.ge.tstop*0.999) then
               print *,'time_,tstop',time_,tstop
               close(31)
               stop
            else
               idump=idump+1
               tnext=tdump(idump)
               print *,'tnext',tnext
            endif
         endif
         goto 100
!
!
   90 format(i3.3)
  131 format(4(a2,1pe9.2))
      stop

      end subroutine testcut

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
      
      subroutine timedump(delt) bind(C)
      use iso_c_binding
      implicit none
      real(c_double)::delt
      integer ndump
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
         tnext=time_+delt
      endif
      print *,'tnext,tstop',tnext,tstop
      end subroutine timedump

      subroutine clphi() bind(C)
      implicit none
      integer isn
      integer lvn
      do isn=1,ns
         do lvn=1,3
            phi(lvn,isn)=1d11*svec(1,lvn,isn)
         enddo
      enddo
      return
      end subroutine clphi

      subroutine clvis() bind(C)    
      do isn=1,ns
         do lvn=1,3
            vis(lvn,isn)=1d6*svec(1,lvn,isn)
         enddo
      enddo
      return
      end subroutine clvis

      subroutine clpsi() bind(C)
      do isn=1,ns
         do lvn=1,3
            psi(lvn,isn)=-1d5*svec(1,lvn,isn)
         enddo
      enddo
      return
      end subroutine clpsi

      subroutine clgam(afac) bind(C)
      implicit none
      real(c_double)::afac
      real(8) afac3
      integer iq,il
      !print *,afac
      afac3=afac**3
      !print *,afac3
      do iq=1,nq
         gamv(iq)=1d-2*afac3
         gamd(iq)=1d-2*afac3
      enddo
      do il=1,nl
         game(il)=1d-2*afac3
      enddo
      return
      end subroutine clgam
      
      subroutine vvec1() bind(C)
      do iq=1,nq
         v4=0d0
         d4=0d0
         do isn=1,4
            is=isoq(isn,iq)
            v4=v4+svec(4,1,is)
            d4=d4+svec(1,3,is)*svec(2,3,is)
         enddo
         vvec(1,iq)=0.25*v4*1d-2
         dvec(1,iq)=0.25*d4
      enddo
      do il=1,nl
         e6=0d0
         do isn=1,2
            is=isol(isn,il)
            do lv=1,3
               e6=e6+svec(1,lv,is)*svec(2,lv,is)
            enddo
         enddo
         evec(1,il)=e6/6d0
      enddo
      return
      end subroutine vvec1

      subroutine clsfr() bind(C)
!idsfr=+1
      do iq=1,nq
         sfrv(iq)=-4d2*vvec(1,iq)
         sfrd(iq)=1d1*dvec(1,iq)
      enddo
      do il=1,nl
         sfre(il)=1d1*evec(1,il)
      enddo
!     do il=1,nl
!        iq=iqol(il)
!        sfrv(iq)=0d0
!     enddo
!idsfr=-1
!     do is=1,ns
!        sfrn(1,is)=-5d0*svec(4,1,is)
!        sfrn(3,is)=1d1*svec(2,3,is)*svec(1,3,is)
!     enddo
!     do il=1,nl
!        is=isol(1,il)
!        sfrn(2,is)=1d1*svec(2,2,is)*svec(1,2,is)
!        sfrn(1,is)=1d1*svec(2,1,is)*svec(1,1,is)
!     enddo
      return
      end subroutine clsfr

      subroutine openunit12() bind(C)
      use iso_c_binding
      open(unit=12,file='test_out')
      call iowrfile(0,12)
      close(12)
      end subroutine openunit12
      
      subroutine clchm() bind(C)
      tautheta=2d0
      theta0=1d-3
      taum=1d0
      do isn=1,ns
         do lvn=1,3
            thetaeq=theta0*(svec(2,lvn,isn)+1d0)
            svec(3,lvn,isn)=thetaeq
            sdot(1,lvn,isn)=(thetaeq-svec(1,lvn,isn))/
     1                      (tautheta/svec(2,lvn,isn))
            sdkr(1,lvn,isn)=-1d0/(tautheta/svec(2,lvn,isn))
            sdot(2,lvn,isn)=(5d-2-svec(2,lvn,isn))/taum
            sdkr(2,lvn,isn)=-1d0/taum
            sdot(4,lvn,isn)=(1d-4-svec(4,lvn,isn))/(2d0*taum)
            sdkr(4,lvn,isn)=-1d0/(2d0*taum)
         enddo
      enddo
      do il=1,nl
!        if (time.lt.5d1) then
!           cxval(4,il)=2d-1*(1d0+time/1d2)
!        else
!           cxval(4,il)=3d-1
!        endif
         if (cxprm(4,il).eq.0d0.and.cxprm(4,ilol(1,il)).eq.0d0) then
            isn=isol(1,il)
            do lvn=1,3
               sdot(1,lvn,isn)=1.0d0*sdot(1,lvn,isn)
               sdkr(1,lvn,isn)=1.0d0*sdkr(1,lvn,isn)
            enddo
!           print *,'decay time',1/sdkr(1,1,isn)
         endif
      enddo
!         
      end subroutine clchm

      subroutine bc() bind(C)
      do iq=1,nq
         vxsrc(2,iq)=vvec(1,iq)
      enddo
      return
      end subroutine bc

      subroutine wrfile(isve,ipf,idot,rmin,rmax,hmax,volint, 
     1                  aratio) bind(C)
      implicit none
      integer(C_INT)::isve,idot,rmin,rmax,hmax,volint,aratio
      character(kind=c_char, len=1), dimension(20), intent(in)::ipf
      character(len=20) :: opf
      write(31,*) time_, rmin,rmax,hmax,volint,aratio
      if (isve.eq.1) then
            isuff=isuff+1
            opf = " "
            loop_string: do i=1, 20
            if ( ipf(i) == c_null_char ) then
               exit loop_string
            else
               opf(i:i) = ipf(i)
            end if
            end do loop_string
            write(fnum,90)isuff
            print *,'isuff=',isuff
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
               tnext=tdump(idump)
               print *,'tnext',tnext
           endif
        endif
   90 format(i3.3)
  131 format(4(a2,1pe9.2))
      end subroutine wrfile
      
      subroutine checkcxprm() bind(C)
      use iolibsw 
      do i =1,nl
         print *,"cxprm ",i,cxprm(4,i)
      enddo
      return
      end subroutine checkcxprm

      subroutine precisearatio(volint,tsurf,aratio) bind(C)
      real(c_double)::volint
      real(c_double)::tsurf
      real(c_double)::aratio
      print *,tsurf
      print *,volint
      r0=(volint*3.0/(4.0*3.141592))**(1./3.)
      a0=4.0*3.141592*r0**2
      print *,r0
      print *,a0
      aratio=tsurf/a0
      print *,aratio
      return
      end subroutine precisearatio




      END MODULE DEFORMLIB

