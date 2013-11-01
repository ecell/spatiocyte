      MODULE CMLIBSW
!
! library of diffusion-reaction routines.
!
! General idea (VOLUME diffusion of concentration c)
!     Qd=<dHjx*dHkx>|iq
!     Qm=<Hj*Hk>|iq
!
!     Qs=stiffness matrix
!       =(1-dt*<d(cdot)/dc>|iq)*Qm+dt*Diff*Qd
!
!     vldq=load vector
!         =Qm*[(1-<d(cdot)/dc>|iq)*c(:)+dt*c(:)]
!
!  boundary contributions:
!     L=<SjSk>*permeability
!     vld_boundary=dt*(source+c_ext*perm)*<S>
!
!  And we solve:
!      (Qs+L's)*c_new(:)= vld's
!
      USE IOLIBSW
      USE MXLIBSW
      USE ISO_C_BINDING
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  these are component of the stiffness and mass matrices that can be 
!  reused as long as the mesh does not move. they are therefore saved.
!--volume diffusion matrices
      real(8),save::Qd(3,4,3,4,NSM),Qm(3,4,3,4,NSM) 
!--surface diffusion matrices
      real(8),save::Qdv(4,4,NSM),Qmv(4,4,NSM)!ventral
      real(8),save::Qdd(4,4,NSM),Qmd(4,4,NSM)!dorsal
      real(8),save::Qde(3,2,3,2,NSM),Qme(3,2,3,2,NSM)!edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      CONTAINS!module subprograms follow
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfdriver(oldopt) bind(C)
!
! This is the main driver for diffusion/reaction. 
! One call evolves all the components of svec 
!
      implicit none
!
      character(kind=c_char,len=1),dimension(18),intent(in)::oldopt !choices are "eul" or "lag"
      character(len=10)::opt
      integer,save::ichk=-1
!
      real(8) con(3*NSM),difcon ! node concentrations, diff. coeff.
!
!--volume diffusion arrays:
      real(8) Qs(3,4,3,4,NSM) !interior stiffness matrix
!               (level,stack,level,stack,element)
!   ventral,dorsal,edge contributions to stiffness matrix:    
      real(8) Lv(4,4,NSM),Ld(4,4,NSM),Le(3,2,3,2,NLM)
      real(8) vldq(3,4,NSM) !interior load vector
!                 (level,stack,element)
!   ventral,dorsal,edge contributions to load vector:
      real(8) vldv(4,NSM),vldd(4,NSM),vlde(3,2,NLM) 
!
!--surface diffusion arrays
!   ventral, dorsal, edge stiffness matrices
      real(8) Qsv(4,4,NSM),Qsd(4,4,NSM),Qse(3,2,3,2,NLM) 
!   contact line contribution to stiffness matrix
      real(8) Lc(2,2,NLM)
!   ventral, dorsal, edge load vectors
      real(8) vldqv(4,NSM),vldqd(4,NSM),vldqe(3,2,NLM)
!   contact line contribution to load vector:
      real(8) vldc(2,NLM)
!
      real(8) C(NSM*3)! stiffness matrix conditioner
!
      integer kom, nskind, js,jl,inode,il,k
      opt=" "
      loop_string: do i=1, 10
      if ( oldopt(i) == c_null_char ) then
         exit loop_string
      else
         opt (i:i) = oldopt(i)
      end if
      end do loop_string

      !print *,opt
!
!--do two checks to see if further calculation is necessary
      if(idsdc.le.0)return!no diffusion descriptors
      if(tstp.le.vtiny)return!no time interval  
!
!--we may have to compute grid-determined mass and stiffness 
!   matrix components because this is eithe the first call to 
!   dfdriver, or the grid is moving (non-eulerian case)
      if (ichk.lt.0) then
        call dfqdfmxs!init Qd and Qm
        call dfqdfmxs2!init Qde,Qdv,Qdd, and Qme,Qmv,Qmd
      endif
      ichk=+1
      if (index(opt,'eul').le.0) then
         call dfqdfmxs!init Qd and Qm
         call dfqdfmxs2!init Qde,Qdv,Qdd, and Qme,Qmv,Qmd
      endif
!
!--start loop to diffuse the components of SVEC
      do kom=1,NKS
         difcon=dscp(kom,idsdc)!diffusion coef. of svec(kom,*,*)
         nskind=nint(dscp(kom,idspk))
!--check diffusion coeff
         if (difcon.lt.vtiny) then 
            call cmdriver(kom)!no diffusion case--pure chemistry
         else
            con=0d0! initialize concentration vector
!--nskind=1, volume variable, do volume diffusion
            !print *,'nskind',nskind
            if (nskind.eq.1) then
               call dfbdymat(Lv,Ld,Le,kom)!assemble bdy permeability mx
               call dfqesmat(Qs,kom)!assemble Qs=stiffness mx
               call cgscndmx(C,Qs,Lv,Ld,Le)!assemble preconditioning mx
               call dfbsrlod(vldv,vldd,vlde,kom)!assemble boundary loads
               call dfmaslod(vldq,kom)!assemble element loads
!--load the concentration vector 
               do js=1,ns
                  do jl=1,3
                     inode=(js-1)*3+jl
                     con(inode)=svec(kom,jl,js)
                     !print *,js,'     ',con(inode)
                  enddo
               enddo
!
!--Solve for con(in) at t=t+dt [existing field will be overwritten!] 
!
!--solve by conjugate gradient method
               call cgssolve(con,Qs,Lv,Ld,Le,C,vldq,vldv,vldd,vlde,20,k)
!--solve by LU decomposition method
!              call lussolve(con,Qs,Lv,Ld,Le,vldq,vldv,vldd,vlde)
!--replace svec(kom,:,:) with new concentrations
               do js=1,ns
                  do jl=1,3
                     inode=(js-1)*3+jl
                     if (con(inode).lt.0.0) then
                       print *,'negative concentration,level,stack,kom',
     1                          con(inode),jl,js,kom
                     endif
                     svec(kom,jl,js)=max(con(inode),vtiny)
                     !print *,js,'   ', svec(kom,jl,js)
                  enddo
               enddo
!
!--nskind=2, surface variable, do surface diffusion
            elseif (nskind.eq.2) then
               call dfqesmat2(Qsv,Qsd,Qse,kom)!assemble stiffness mx
               call cgscndmx2(C,Qsv,Qsd,Qse)!preconditioning matrix
               call dfmaslod2(vldqv,vldqd,vldqe,kom)!surface elmt loads
!--load ventral, dorsal, and edge nodes into concentration vector
               do js=1,ns
                  inode=js
                  con(inode)=svec(kom,1,js)
               enddo
               do js=1,ns
                  inode=js+ns
                  con(inode)=svec(kom,3,js)
               enddo
               do il=1,nl
                  js=isol(1,il)
                  inode=il+2*ns
                  con(inode)=svec(kom,2,js)
               enddo
!--Solve for con(in) at t=t+dt [existing field will be overwritten!] 
!
!--solve by conjugate gradient
               call cgssolve2(con,Qsv,Qsd,Qse,C,vldqv,vldqd,vldqe,20)
!--solve by LU decomposition:
!              call lussolve2(con,Qsv,Qsd,Qse,vldqv,vldqd,vldqe)
!
               do inode=1,2*ns+nl
                  if (con(inode).lt.0.0) then
!                     print *,'negative surface concentration,inode,kom',
!     1                        con(inode),inode,kom
                  endif
               enddo 
!--replace svec(kom,:,:) with new concentrations
               do js=1,ns
                  inode=js
                  svec(kom,1,js)=max(con(inode),vtiny)
               enddo
               do js=1,ns
                  inode=js+ns
                  svec(kom,3,js)=max(con(inode),vtiny)
               enddo
               do il=1,nl
                  js=isol(1,il)
                  inode=il+2*ns
                  svec(kom,2,js)=max(con(inode),vtiny)
               enddo
!--nskind=3, bottom surface variable, do surface diffusion
            elseif (nskind.eq.3) then
               call dfbdymat3(Lc,kom)!assemble bdy permeability mx
               call dfqesmat3(Qsv,kom)!assemble stiffness mx
               call cgscndmx3(C,Qsv,Lc)!preconditioning matrix
               call dfbsrlod3(vldc,kom)!assemble boundary loads
               call dfmaslod3(vldqv,kom)!surface elmt loads
!--load ventral into concentration vector
               do js=1,ns
                  inode=js
                  con(inode)=svec(kom,1,js)
                  !print *,inode,con(inode)
               enddo
!--Solve for con(in) at t=t+dt [existing field will be overwritten!] 
!
!--solve by conjugate gradient
               call cgssolve3(con,Qsv,Lc,C,vldqv,vldc,20)
!--solve by LU decomposition:
!              call lussolve3(con,Qsv,Lc,vldqv,vldc)
!
               do inode=1,ns
                  if (con(inode).lt.0.0) then
                     print *,'negative surface concentration,inode,kom',
     1                        con(inode),inode,kom
                  endif
               enddo 
!--replace svec(kom,:,:) with new concentrations
               do js=1,ns
                  inode=js
                  svec(kom,1,js)=max(con(inode),vtiny)
                  !print *,'nskind3','  ',js,'  ',svec(kom,1,js)
               enddo
            endif
         endif
      enddo
      return
      end subroutine dfdriver!(opt)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cmdriver(kom)
!
! this subroutine will do the chemistry on svec node field kom
! DO NOT CALL at the same time as dfdriver as this would do 
! the chemistry twice, unless there are no diffusion descriptors
!
      implicit none
!
      integer kom!the component of svec to do chemistry on
      real(8) sfac
      integer nskind,is,ilv,istack,iq,il
!
      nskind=nint(dscp(kom,idspk))
!--nskind=1: volume variable, do volume chemistry
      if (nskind.eq.1) then
         do is=1,ns
            do ilv=1,3
               sfac=1d0-sdkr(kom,ilv,is)*tstp
               if (sfac.gt.vtiny) then
                  svec(kom,ilv,is)=(svec(kom,ilv,is)*sfac+
     1                              sdot(kom,ilv,is)*tstp)/sfac
               else
                 print *,'time step too big!,stack,level,kom',is,ilv,kom
                 print *,'1/sdkr,tstp',1d0/sdkr(kom,ilv,is),tstp
               endif
            enddo
         enddo
!--nskind=2: surface variable, do surface chemistry
      elseif(nskind.eq.2) then
         do iq=1,nq
            do is=1,4
               istack=isoq(is,iq)
!--do ventral and dorsal nodes
               do ilv=1,3,2
                  sfac=1d0-sdkr(kom,ilv,istack)*tstp
                  if (sfac.gt.vtiny) then
                     svec(kom,ilv,istack)=(svec(kom,ilv,istack)*sfac+
     1                                   sdot(kom,ilv,istack)*tstp)/sfac
                  else
                     print *,'time step too big!,stack,level,kom',
     1                        is,ilv,kom
                     print *,'1/sdkr,tstp',1d0/sdkr(kom,ilv,istack),tstp
                  endif
               enddo
            enddo
         enddo
!--do middle edge nodes
         do il=1,nl
            istack=isol(1,il)
            sfac=1d0-sdkr(kom,2,istack)*tstp
            if (sfac.gt.vtiny) then
                svec(kom,2,istack)=(svec(kom,2,istack)*sfac+
     1                              sdot(kom,2,istack)*tstp)/sfac
            else
               print *,'time step too big!edge,kom',il,kom
               print *,'1/sdkr,tstp',1d0/sdkr(kom,2,istack),tstp
            endif
         enddo
      endif
      return
      end subroutine cmdriver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! begin stiffness matrix assembly routines for diffusion operators
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfqdfmxs
!
!  this subroutine computes grid and only grid-dependent VOLUME
!  diffusion matrices.
!   Qd(jl,js,kl,ks,iq)=dVol*<(dHjx))*(dHkx)+(dHjy)*(dHky)+
!                                          +(dHjz)*(dHkz)>|iq
!   Qm(jl,js,ks,kl,iq)=dVol*<H(j)*H(k)>|iq 
!
      implicit none
!
      integer iq,isg,lvg,js,jl,ks,kl
      real(8) dvol,Fqd,Fqm
!
      Qd=0d0
      Qm=0d0
!--loop over elements
      do iq=1,nq
!--loop over Gauss points stacks and levels
         do isg=1,4
            do lvg=1,3
               dvol=det(lvg,isg,iq)*wgp(lvg)! weighted vol. of GP
c--loop over j and k nodes by stack and level
               do ks=1,4
                  do kl=1,3
                     do js=1,4
                        do jl=1,3
             Fqd=dHg(1,jl,js,lvg,isg,iq)*dHg(1,kl,ks,lvg,isg,iq)+
     1           dHg(2,jl,js,lvg,isg,iq)*dHg(2,kl,ks,lvg,isg,iq)+
     2           dHg(3,jl,js,lvg,isg,iq)*dHg(3,kl,ks,lvg,isg,iq)
             Qd(jl,js,kl,ks,iq)=Qd(jl,js,kl,ks,iq)+Fqd*dvol
             Fqm=H(jl,js,lvg,isg)*H(kl,ks,lvg,isg)
             Qm(jl,js,kl,ks,iq)=Qm(jl,js,kl,ks,iq)+Fqm*dvol
                        enddo
                     enddo
                  enddo
               enddo
            enddo 
         enddo
      enddo
      return
      end subroutine dfqdfmxs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfqdfmxs2
!
!  this subroutine computes grid and only grid-dependent SURFACE
!  diffusion matrices.
!  Qdv(js,ks,iq)=dArea*<(dS4vjx))*(dS4vkx)+(dS4vjy)*(dS4vky)
!                                         +(dS4vjz)*(dS4vkz)>|iq
!  Qdd(js,ks,iq)=dArea*<(dS4djx))*(dS4dkx)+(dS4djy)*(dS4dky)
!                                         +(dS4djz)*(dS4dkz)>|iq
!  Qde(jl,js,kl,ks,iq)=dArea*<(dS6djx))*(dS6dkx)+(dS6djy)*(dS6dky)
!                                              +(dS6djz)*(dS6dkz)>|il
!  Qmv(js,ks,iq)=dArea*<(S4j)*(S4k)>|iq
!  Qmd(js,ks,iq)=dArea*<(S4j)*(S4k)>|iq
!  Qme(jl,js,kl,ks,iq)=dArea*<(S6j)*(S6k)>|il
!
      implicit none
!
      integer iq,il,isg,lvg,js,jl,ks,kl
      real(8) dAg,dAvg,dAdg,Fqd,Fqm,Fqdd,Fqdv,Fqme,Fqmv
!
!--ventral/dorsal surfaces 
      Qdv=0d0
      Qmv=0d0
      Qdd=0d0
      Qmd=0d0
!--loop over elements
      do iq=1,nq
!--loop over 4 GP 
         do isg=1,4
           dAvg=dAv(0,isg,iq)
           dAdg=dAd(0,isg,iq)
            do ks=1,4
               do js=1,4
                  Fqdv=dS4vg(1,js,isg,iq)*dS4vg(1,ks,isg,iq)+
     1                 dS4vg(2,js,isg,iq)*dS4vg(2,ks,isg,iq)+
     2                 dS4vg(3,js,isg,iq)*dS4vg(3,ks,isg,iq)
                  Qdv(js,ks,iq)=Qdv(js,ks,iq)+Fqdv*dAvg
                  Fqdd=dS4dg(1,js,isg,iq)*dS4dg(1,ks,isg,iq)+
     1                 dS4dg(2,js,isg,iq)*dS4dg(2,ks,isg,iq)+
     2                 dS4dg(3,js,isg,iq)*dS4dg(3,ks,isg,iq)
                  Qdd(js,ks,iq)=Qdd(js,ks,iq)+Fqdd*dAdg
                  Fqm=S4(js,isg)*S4(ks,isg)
                  Qmv(js,ks,iq)=Qmv(js,ks,iq)+Fqm*dAvg
                  Qmd(js,ks,iq)=Qmd(js,ks,iq)+Fqm*dAdg
               enddo
            enddo 
         enddo
      enddo
!--edge surfaces 
      Qde=0d0
      Qme=0d0
!--loop over edges
      do il=1,nl
!--six GPs in two stacks, 3 levels
         do isg=1,2
            do lvg=1,3
               dAg=wgp(lvg)*dAe(0,lvg,isg,il)
               do ks=1,2
                  do kl=1,3
                     do js=1,2
                        do jl=1,3
             Fqd=dS6g(1,jl,js,lvg,isg,il)*dS6g(1,kl,ks,lvg,isg,il)+
     1           dS6g(2,jl,js,lvg,isg,il)*dS6g(2,kl,ks,lvg,isg,il)+
     2           dS6g(3,jl,js,lvg,isg,il)*dS6g(3,kl,ks,lvg,isg,il)
             Qde(jl,js,kl,ks,il)=Qde(jl,js,kl,ks,il)+Fqd*dAg
             Fqm=S6(jl,js,lvg,isg)*S6(kl,ks,lvg,isg)
             Qme(jl,js,kl,ks,il)=Qme(jl,js,kl,ks,il)+Fqm*dAg
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine dfqdfmxs2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfqesmat(Qs,kom)
!
! This subroutine assembles the interior contribution to the stiffness 
! matrix for VOLUME diffusion. 
!
      implicit none
!
! stiffness matrix
      real(8) Qs(3,4,3,4,NSM)!stiffness matrix
!               (jlevel,jstack,klevel,kstack,element)
      integer kom !index of svec field
!
      real(8) DC,DK,QDCOF,QMCOF
      integer iq,isn,istack,lvn,ks,kl,js,jl
!--node weight parameters 1/12*(1/6,2/3,1/6) for averaging 
!   a field over the element:
      real(8), dimension(3), parameter :: wn=(/0.013888888889,
     1                                         0.055555555556,
     2                                         0.013888888889/)
!
      DC=dscp(kom,idsdc)!load diffusion const of svec(kom,*)
      QDCOF=DC*tstp!diffusion const x time step.
!--loop over elements
      do iq=1,nq
         DK=0d0 !decay rate = d(sdot)/d(svec)
!--average decay rate over element
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               DK=DK+wn(lvn)*sdkr(kom,lvn,istack)
            enddo
         enddo
         QMCOF=1.0-DK*tstp!1.0-decay rate x time step.
!--double loop over all element nodes
         do ks=1,4
            do kl=1,3
               do js=1,4
                  do jl=1,3
                     Qs(jl,js,kl,ks,iq)=QMCOF*Qm(jl,js,kl,ks,iq)+
     1                                  QDCOF*Qd(jl,js,kl,ks,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine dfqesmat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfqesmat2(Qsv,Qsd,Qse,kom)
!
! This subroutine assembles the stiffness matrices for SURFACE 
! diffusion. 
!
      implicit none
!
! ventral, dorsal, and edge stiffness matrices:
      real(8) Qsv(4,4,NSM),Qsd(4,4,NSM) !(jstack,kstack,element)
      real(8) Qse(3,2,3,2,NLM)!(jlevel,jstack,klevel,kstack,edge)
      integer kom !svec field index
!
      real(8) DC, DK,DKd,DKv,QDCOF,QMCOF,QMCOFv,QMCOFd
      integer il,iq,isn,istack,lvn,ks,kl,js,jl
!
!--node weight parameters 1/6*(1/6,2/3,1/6) for averaging a field
!   over an edge surface.
      real(8), dimension(3), parameter :: wn=(/0.027777777778,
     1                                         0.111111111111,
     2                                         0.027777777778/)
!
      DC=dscp(kom,idsdc)!load diffusion const of svec(kom,*)
      QDCOF=DC*tstp!diffusion const x time step.
!--ventral/dorsal
      do iq=1,nq
         DKv=0d0
         DKd=0d0
!--average decay rates (d(sdot)/d(s)) over quad. surface
         do isn=1,4
            istack=isoq(isn,iq)
            DKv=DKv+0.25*sdkr(kom,1,istack)
            DKd=DKd+0.25*sdkr(kom,3,istack)
         enddo
         QMCOFv=1d0-DKv*tstp!1.0-decay rate x time step.
         QMCOFd=1d0-DKd*tstp!1.0-decay rate x time step.
         do ks=1,4
            do js=1,4
               Qsv(js,ks,iq)=QMCOFv*Qmv(js,ks,iq)+QDCOF*Qdv(js,ks,iq)
               Qsd(js,ks,iq)=QMCOFd*Qmd(js,ks,iq)+QDCOF*Qdd(js,ks,iq)
            enddo
         enddo
      enddo
!--edges
      do il=1,nl
         DK=0d0
!--average decay rates (d(sdot)/d(s)) over edge surface
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               DK=DK+wn(lvn)*sdkr(kom,lvn,istack)
            enddo
         enddo
         QMCOF=1.0-DK*tstp!1.0+decay rate x time step.
         do ks=1,2
            do kl=1,3
               do js=1,2
                  do jl=1,3
                     Qse(jl,js,kl,ks,il)=QMCOF*Qme(jl,js,kl,ks,il)+
     1                                   QDCOF*Qde(jl,js,kl,ks,il)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine dfqesmat2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfqesmat3(Qsv,kom)
!
! This subroutine assembles the stiffness matrices for BOTTOM SURFACE 
! diffusion. 
!
      implicit none
!
! ventral stiffness matrices:
      real(8) Qsv(4,4,NSM) !(jstack,kstack,element)
      integer kom !svec field index
!
      real(8) DC, DK,DKv,QDCOF,QMCOF,QMCOFv
      integer il,iq,isn,istack,lvn,ks,kl,js,jl
!
      DC=dscp(kom,idsdc)!load diffusion const of svec(kom,*)
      !print *,'DC',DC
      QDCOF=DC*tstp!diffusion const x time step.
      !print *,'QDCOF',QDCOF
!--ventral/dorsal
      do iq=1,nq
         DKv=0d0
!--average decay rates (d(sdot)/d(s)) over quad. surface
         do isn=1,4
            istack=isoq(isn,iq)
            DKv=DKv+0.25*sdkr(kom,1,istack)
            !print *,'DKv',DKv
         enddo
         QMCOFv=1d0-DKv*tstp!1.0-decay rate x time step.
         !print *,'QMCOFv',QMCOFv
         do ks=1,4
            do js=1,4
               Qsv(js,ks,iq)=QMCOFv*Qmv(js,ks,iq)+QDCOF*Qdv(js,ks,iq)
               !print *,'Qsv',Qsv(js,ks,iq)
            enddo
         enddo
      enddo
      return
      end subroutine dfqesmat3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfbdymat(Lv,Ld,Le,kom)
!
! this subroutine computes the ventral, dorsal, and edge surface 
! contributions Lv, Ld, and Le to the stiffness matrix for VOLUME 
! diffusion.  These are only non-zero in the setting of Dirichlet 
! boundary condition for which an external value is set along with
! a permeability.
!    xprm(kom,*)=bdy permeability  (volume/area/time) 
!
      implicit none
!
      real(8) Lv(4,4,NSM), Ld(4,4,NSM) !ventral and dorsal mx
!                (jstack,kstack,element)
      real(8) Le(3,2,3,2,NLM) !edge mx
!                (jlevel,jstack,klevel,kstack,edge)
      integer kom
!
      real(8) Fee, Fed, Fev, facg, facgv, facgd, fack,fackd,fackv
      integer il,isg,lvg,ks,kl,js,jl,iq
!
!--dorsal and ventral surfaces
      Ld=0d0
      Lv=0d0
!--loop over elements top and bottom surfaces
      do iq=1,nq
         Fev=vxprm(kom,iq)*tstp!ventral perm. at iq for svec(kom,:,:)
         Fed=dxprm(kom,iq)*tstp!dorsal permeability at iq for svec(kom
!--loop over Gauss points, dA elementary area of GP
         do isg=1,4
            facgd=Fed*dAd(0,isg,iq)
            facgv=Fev*dAv(0,isg,iq)
!--now loop 4x4 pairs of nodes from ventral and dorsal surfaces
            do ks=1,4
               fackd=facgd*S4(ks,isg)
               fackv=facgv*S4(ks,isg)
               do js=1,4
                  Ld(js,ks,iq)=Ld(js,ks,iq)+fackd*S4(js,isg)
                  Lv(js,ks,iq)=Lv(js,ks,iq)+fackv*S4(js,isg)
               enddo
            enddo
         enddo
      enddo
!--edge surfaces
      Le=0d0
      do il=1,nl
         Fee=exprm(kom,il)*tstp
!--loop over 2 stacks x 3 levels of GPs, dA elementary area of GP
         do isg=1,2
            do lvg=1,3
               facg=wgp(lvg)*dAe(0,lvg,isg,il)*Fee
!--loop over 6x6 nodes of the edge surface
               do ks=1,2
                  do kl=1,3
                     fack=facg*S6(kl,ks,lvg,isg)
                     do js=1,2
                        do jl=1,3
                           Le(jl,js,kl,ks,il)=Le(jl,js,kl,ks,il)+
     1                                        fack*S6(jl,js,lvg,isg)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end subroutine dfbdymat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfbdymat3(Lc,kom)
!
! this subroutine computes the contact line
! contribution Lc, to the stiffness matrix for surface
! diffusion.  These are only non-zero in the setting of Dirichlet 
! boundary condition for which an external value is set along with
! a permeability.
!    xprm(kom,*)=bdy permeability  (volume/area/time) 
!
      implicit none
!
      real(8) Lc(2,2,NLM)
      integer kom
!
      real(8) Fc, Ladd
      integer il,isg,lvg,ks,kl,js,jl,iq
!
!--run around the contact line with one GP integration (hence 0.25=0.5*0.5)
      Lc=0d0
      do il=1,nl
         Fc=cxprm(kom,il)*tstp
         !print *,'Fc',il,kom,Fc,cxprm(kom,il),tstp
         Ladd=0.25d0*cnn(0,il)*Fc
         Lc(1,1,il)=Ladd
         Lc(2,1,il)=Ladd
         Lc(1,2,il)=Ladd
         Lc(2,2,il)=Ladd
      enddo
!
      return
      end subroutine dfbdymat3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   end stiffness matrix assembly routines for diffusion operators
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   begin load vector assembly routines for diffusion operators
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine dfmaslod(vldq,kom)
!
! this subroutine computes the interior contribution to the load 
! vector for VOLUME diffusion (backward-Euler diffusion step)
!       c_rhs=c_old+dt*c_dot-dt*(dc_dot/dc)*c_old
!
! It makes use of the mass matrix: Qm(jl,js,kl,ks,iq)=<Hj*Hk>
!
      implicit none
      real(8) vldq(3,4,NSM)! load vector (node lev,node stack,element)
      integer kom! component of svec to be evolved
!
      real(8) c_rhs,DK
      integer iq,istack,lvn,isn,js,jl
!
!--node weight parameters 1/12*(1/6,2/3,1/6) for element averaging
      real(8), dimension(3), parameter :: wn=(/0.013888888889,
     1                                         0.055555555556,
     2                                         0.013888888889/)
!
      vldq=0d0
!--loop over elements
      do iq=1,nq
!--average to a single decay rate for each element
         DK=0d0
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               DK=DK+wn(lvn)*sdkr(kom,lvn,istack)
            enddo
         enddo
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               c_rhs=svec(kom,lvn,istack)*(1d0-DK*tstp)+
     1               sdot(kom,lvn,istack)*tstp
               do js=1,4
                  do jl=1,3
                     vldq(jl,js,iq)=vldq(jl,js,iq)+
     1                              c_rhs*Qm(jl,js,lvn,isn,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo   
      return
      end subroutine dfmaslod!(vldq,kom)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfmaslod2(vldqv,vldqd,vldqe,kom)
!
! this subroutine computes the ventral, dorsal and edge contributions 
! vldqv, vldqd, vldqe, to the load vector for SURFACE diffusion 
!                                  (backward-Euler diffusion step)
!       c_rhs=c_old+dt*c_dot-dt*(dc_dot/dc)*c_old
!
! It makes use of the mass matrices: Qmv/d(js,ks,iq)=<S4j*S4k> and
!                                    Qe(jl,js,kl,ks,il)=<S6j*S6k>
!
      implicit none
!
      real(8) vldqv(4,NSM),vldqd(4,NSM)!ventral and dorsal loads
      real(8) vldqe(3,2,NLM)!edge loads
      integer kom!component of svec
!
      real(8) c_rhs,c_rhsv,c_rhsd,DK,DKv,DKd
      integer il,iq,istack,lvn,isn,js,jl
!
!--node weight parameters 1/6*(1/6,2/3,1/6)
      real(8), dimension(3), parameter :: wn=(/0.027777777778,
     1                                         0.111111111111,
     2                                         0.027777777778/)
!
!--ventral and dorsal loads 
      vldqv=0d0
      vldqd=0d0
!--loop over all element ventral and dorsal surfaces
      do iq=1,nq
!--compute a single decay for each surface 
         DKv=0d0
         DKd=0d0
         do isn=1,4
            istack=isoq(isn,iq)
            DKv=DKv+0.25*sdkr(kom,1,istack)
            DKd=DKd+0.25*sdkr(kom,3,istack)
         enddo
         do isn=1,4
            istack=isoq(isn,iq)
            c_rhsv=svec(kom,1,istack)*(1d0-DKv*tstp)+
     1             sdot(kom,1,istack)*tstp
            c_rhsd=svec(kom,3,istack)*(1d0-DKd*tstp)+
     1             sdot(kom,3,istack)*tstp
! each node contributes to 4 nodes
            do js=1,4
               vldqv(js,iq)=vldqv(js,iq)+c_rhsv*Qmv(js,isn,iq)
               vldqd(js,iq)=vldqd(js,iq)+c_rhsd*Qmd(js,isn,iq)
            enddo
         enddo
      enddo   
!--edge loads
      vldqe=0d0
      do il=1,nl
!--compute a single decay for each edge surfaces 
         DK=0d0
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               DK=DK+wn(lvn)*sdkr(kom,lvn,istack)
            enddo
         enddo
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               c_rhs=svec(kom,lvn,istack)*(1d0-DK*tstp)+
     1               sdot(kom,lvn,istack)*tstp
               do js=1,2
                  do jl=1,3
                     vldqe(jl,js,il)=vldqe(jl,js,il)+
     1                               c_rhs*Qme(jl,js,lvn,isn,il)
                  enddo
               enddo
            enddo
         enddo
      enddo  
      return
      end subroutine dfmaslod2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfmaslod3(vldqv,kom)
!
! this subroutine computes the ventral contributions 
! vldqv to the load vector for BOTTOM SURFACE diffusion 
!                                  (backward-Euler diffusion step)
!       c_rhs=c_old+dt*c_dot-dt*(dc_dot/dc)*c_old
!
! It makes use of the mass matrix: Qmv/d(js,ks,iq)=<S4j*S4k> 
!
      implicit none
!
      real(8) vldqv(4,NSM)!ventral loads
      integer kom!component of svec
!
      real(8) c_rhs,c_rhsv,DK,DKv
      integer il,iq,istack,lvn,isn,js,jl
!
!--ventral loads 
      vldqv=0d0
!--loop over all element ventral surfaces
      do iq=1,nq
!--compute a single decay for each surface 
         DKv=0d0
         do isn=1,4
            istack=isoq(isn,iq)
            DKv=DKv+0.25*sdkr(kom,1,istack)
         enddo
         do isn=1,4
            istack=isoq(isn,iq)
            c_rhsv=svec(kom,1,istack)*(1d0-DKv*tstp)+
     1             sdot(kom,1,istack)*tstp
! each node contributes to 4 nodes
            do js=1,4
               vldqv(js,iq)=vldqv(js,iq)+c_rhsv*Qmv(js,isn,iq)
            enddo
         enddo
      enddo   
      return
      end subroutine dfmaslod3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfbsrlod(vldv,vldd,vlde,kom)
!
! boundary contributions to node loads for VOLUME diffusion.
! 1) Neumann boundary condition -> source term 
! 2) Dirichlet boundary condition -> permeability/external value terms
!    xsrc(kom,*)=bdy source strength in (mass/area/time) 
!    xprm(kom,*)=bdy permeability  (volume/area/time) 
!    xval(kom,*)=external value of svec 
!
      implicit none
!
      real(8) vldd(4,NSM),vldv(4,NSM)!dorsal and ventral loads
      real(8) vlde(3,2,NLM)!edge load
      integer kom
!
      real(8) Fld, Fldv, Fldd, facg, facgd, facgv
      integer il, isg, lvg, isn, lvn, iq
!
!--ventral and dorsal surfaces
      vldd=0d0
      vldv=0d0
!--loop over elements ventral and dorsal sides
      do iq=1,nq
         Fldv=(vxsrc(kom,iq)+vxval(kom,iq)*vxprm(kom,iq))*tstp
         Fldd=(dxsrc(kom,iq)+dxval(kom,iq)*dxprm(kom,iq))*tstp
!--loop over GPs
         do isg=1,4
            facgv=dAv(0,isg,iq)*Fldv
            facgd=dAd(0,isg,iq)*Fldd
!--add contribution to four nodes
            do isn=1,4
               vldv(isn,iq)=vldv(isn,iq)+S4(isn,isg)*facgv
               vldd(isn,iq)=vldd(isn,iq)+S4(isn,isg)*facgd
            enddo
         enddo
      enddo
!--edge surfaces
      vlde=0d0
      do il=1,nl
         Fld=(exsrc(kom,il)+exval(kom,il)*exprm(kom,il))*tstp
!--loop over GPs
         do isg=1,2
            do lvg=1,3
               facg=dAe(0,lvg,isg,il)*Fld*wgp(lvg)
!--add conntribution to 6 nodes
               do isn=1,2
                  do lvn=1,3
                     vlde(lvn,isn,il)=vlde(lvn,isn,il)+
     1                                facg*S6(lvn,isn,lvg,isg)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine dfbsrlod
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dfbsrlod3(vldc,kom)
!
! boundary (contact line) contributions to node loads for BOTTOM SURFACE 
! diffusion.
! 1) Neumann boundary condition -> source term 
! 2) Dirichlet boundary condition -> permeability/external value terms
!    cxsrc(kom,*)=bdy source strength in (mass/area/time) 
!    cxprm(kom,*)=bdy permeability  (volume/area/time) 
!    cxval(kom,*)=external value of svec 
!
      implicit none
!
      real(8) vldc(2,NLM)!contact line loads
      integer kom
!
      real(8) Flc, vadd
      integer il 
!
!--contact line boundary 
!--one GP integration
      vldc=0d0
      do il=1,nl
         Flc=(cxsrc(kom,il)+cxval(kom,il)*cxprm(kom,il))*tstp
         vadd=Flc*0.5d0*cnn(0,il)
         vldc(1,il)=vldc(1,il)+vadd
         vldc(2,il)=vldc(2,il)+vadd
      enddo
!
      return
      end subroutine dfbsrlod3
! 
! end assembly routines for diffusion operator
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cmstpsiz(tstep,accu) bind(C)
!
! This subroutine generates a timestep recommendation for
! svec evolution based on second derivative estimates
! i.e. sdot/sddot
!
      implicit none
! 
      real(c_double)::tstep(NKN)
      real(c_double)::accu
!
      real(8),save::svec0(NKN,3,NSM)=0d0
      real(8),save::svec1(NKN,3,NSM)=0d0
      real(8),save::svec2(NKN,3,NSM)=0d0
      real(8),save::t2=0d0,t1=0d0,t0=0d0
      integer,save::kount=0
!
      integer is,ilv,k,il,nskind,ismin3
      real(8) d2sdt2, dsdt, s0, step
      !do i=1,12
      !print *,tstep(i)
      !enddo
!
! shift 1 -> 2; 0 -> 1; current -> 0
      t2=t1
      t1=t0
      t0=time_
      svec2=svec1
      svec1=svec0
      svec0=svec
      kount=kount+1
      if (kount.lt.3) then 
         tstep=max(tstp,tiny)
         return
      endif
      tstep=big
! loop over svec fields
      do k=1,NKN
         nskind=nint(dscp(k,idspk))
         if (nskind.eq.1) then !volume field
            do is=1,ns
               do ilv=1,3
!--second derivative estimated at t0
                  d2sdt2=abs(2*svec2(k,ilv,is)/((t2-t1)*(t2-t0))+
     1                       2*svec1(k,ilv,is)/((t1-t0)*(t1-t2))+
     2                       2*svec0(k,ilv,is)/((t0-t1)*(t0-t2)))
!                 dsdt=abs(svec2(k,ilv,is)*(t0-t1)/((t2-t1)*(t2-t0))+
!    1                     svec1(k,ilv,is)*(t0-t2)/((t1-t0)*(t1-t2))+
!    2                   svec0(k,ilv,is)*(2*t0-t1-t2)/((t0-t1)*(t0-t2)))
                  s0=svec0(k,ilv,is)
                  step=sqrt((2d0*(s0+tiny)*accu)/max(d2sdt2,vtiny))
                  tstep(k)=min(tstep(k),step)
               enddo
            enddo
         elseif(nskind.eq.2) then !surface field
            do is=1,ns
               do ilv=1,3,2! ventral and dorsal nodes
                  d2sdt2=abs(2*svec2(k,ilv,is)/((t2-t1)*(t2-t0))+
     1                       2*svec1(k,ilv,is)/((t1-t0)*(t1-t2))+
     2                       2*svec(k,ilv,is)/((t0-t1)*(t0-t2)))
!                 dsdt=abs(svec2(k,ilv,is)*(t0-t1)/((t2-t1)*(t2-t0))+
!    1                     svec1(k,ilv,is)*(t0-t2)/((t1-t0)*(t1-t2))+
!    2                    svec(k,ilv,is)*(2*t0-t1-t2)/((t0-t1)*(t0-t2)))
                  s0=svec(k,ilv,is)
                  step=sqrt((2d0*(s0+tiny)*accu+tiny)/max(d2sdt2,vtiny))
                  tstep(k)=min(tstep(k),step)
               enddo
            enddo
            do il=1,nl ! pure (middle) edge nodes
               is=isol(1,il)
               d2sdt2=abs(2*svec2(k,2,is)/((t2-t1)*(t2-t0))+
     1                    2*svec1(k,2,is)/((t1-t0)*(t1-t2))+
     2                    2*svec(k,2,is)/((t0-t1)*(t0-t2)))
!              dsdt=abs(svec2(k,2,is)*(t0-t1)/((t2-t1)*(t2-t0))+
!    1                  svec1(k,2,is)*(t0-t2)/((t1-t0)*(t1-t2))+
!    2                  svec(k,2,is)*(2*t0-t1-t2)/((t0-t1)*(t0-t2)))
               s0=svec(k,2,is)
               step=sqrt((2d0*(s0+tiny)*accu+tiny)/max(d2sdt2,vtiny))
               tstep(k)=min(tstep(k),step)
            enddo
         elseif(nskind.eq.3) then !bottom surface field
            do is=1,ns
               ilv=1 !ventral nodes only
               d2sdt2=abs(2*svec2(k,ilv,is)/((t2-t1)*(t2-t0))+
     1                    2*svec1(k,ilv,is)/((t1-t0)*(t1-t2))+
     2                    2*svec(k,ilv,is)/((t0-t1)*(t0-t2)))
!                 dsdt=abs(svec2(k,ilv,is)*(t0-t1)/((t2-t1)*(t2-t0))+
!    1                     svec1(k,ilv,is)*(t0-t2)/((t1-t0)*(t1-t2))+
!    2                    svec(k,ilv,is)*(2*t0-t1-t2)/((t0-t1)*(t0-t2)))
               s0=svec(k,ilv,is)
               step=sqrt((2d0*(s0+tiny)*accu+tiny)/max(d2sdt2,vtiny))
!              svec(5,ilv,is)=step
!              print *,'is,step',is,step
               if (step.lt.tstep(k)) then
!                 ismin3=is
                  tstep(k)=min(tstep(k),step)
               endif
            enddo
         endif
      enddo
      !print *,'ismin3',ismin3
      !print *,'cmstpz: tstep(1:4)',tstep(1:4)
      tstep=min(tstep,tstp*1.1)
      tstep=max(tstep,tstp*0.9)
      return
      end subroutine cmstpsiz  
! 
      end MODULE CMLIBSW! diffusion-reaction library.
