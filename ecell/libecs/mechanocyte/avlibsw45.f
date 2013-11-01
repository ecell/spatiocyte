      MODULE AVLIBSW
      USE ISO_C_BINDING
!
! library of advection/mesh motion routines.
!
      USE IOLIBSW
      USE lupack
!
      real(8),save::xzr(3,3,NSM)!initial node positions
      real(8),save::xnw(3,3,NSM)!node positions after
                                !lagrangian motion by vnw
      real(8),save::xca(3,3,NSM)!node positions after
                                !lagrangian motion by vvol
      real(8),save::xfn(3,3,NSM)!final node positions
      real(8),save::utg1(3,3,NSM)!useful unit tangent vectors
      real(8),save::utg2(3,3,NSM)!useful unit tangent vectors
      real(8),save::dilfacv(3,NSM)!volume dilation coefficient for node
      real(8),save::dilfaca(3,NSM)!area dilation coefficient for node
      logical,save::avnw
      logical,save::avca
      integer,save::clfix(NSM)
!
      contains!the following subroutines  used in advection operations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
      subroutine avgridmo(oldopt) bind(C)
!
! This subroutine displaces the grid
!
      real(8) volzr(3,NSM), areazr(3,NSM)
      real(8) dt
      integer nadv,iflag,imx(2)
      character(kind=c_char,len=1),dimension(18),intent(in)::oldopt !choices are "eul" or "lag"
      character(len=10)::opt
!
      opt=" "
      loop_string: do i=1, 10
      if ( oldopt(i) == c_null_char ) then
         exit loop_string
      else
         opt (i:i) = oldopt(i)
      end if
      end do loop_string
      dt=tstp
!
! load initial node positions and update according to velocities
! recall hvec(1:3,:,:) = x-y-z positions
!        hvec(7:9,:,:) = vnx-vny-vnz network velocities
!        hvec(4:6,:,:) = vvrlx-vvrly-vvrlz volume vel. wrt to network
!
      xzr=hvec(1:3,:,:)
      xnw=xzr+dt*hvec(7:9,:,:)
      xca=xzr+dt*(hvec(7:9,:,:)+hvec(4:6,:,:))
      avnw=.false.
      avca=.false.
!     print *,'xzr(3,3,547),hvec(9,3,547)',xzr(3,3,547),hvec(9,3,547)
!     print *,'xnw(3,3,547)',xnw(3,3,547)
!
!--compute volume and area dilation for network motion
!--(volume motion dilation is zero)
!
!--calculate and save initial volumes and areas
      call godriver(xzr)
      volzr=voln
      areazr=arean
      call godriver(xnw)
      do is=1,ns
         dilfacv(1,is)=volzr(1,is)/voln(1,is)
         dilfacv(2,is)=volzr(2,is)/voln(2,is)
         dilfacv(3,is)=volzr(3,is)/voln(3,is)
         dilfaca(1,is)=areazr(1,is)/arean(1,is)
         dilfaca(3,is)=areazr(3,is)/arean(3,is)
         if (ilos(1,is).ne.0) then !an edge stack
            dilfaca(2,is)=areazr(2,is)/arean(2,is)
         else
            dilfaca(2,is)=0d0
         endif
      enddo
!     print *,'max dilfacv',maxval(dilfacv(:,1:ns))
      imx=maxloc(dilfacv(:,1:ns))
!     print *,'level,stack',imx
!     print *,'volzr,voln',volzr(imx(1),imx(2)),voln(imx(1),imx(2))
!     print *,'hvec(7:9,2,imx(2))',hvec(7:9,2,imx(2))
!     print *,'hvec(7:9,3,imx(2))',hvec(7:9,3,imx(2))
!
      if (index(opt,'eul').gt.0) then !mesh does not move.
         xfn=xzr(:,:,:) !return mesh to original position
      else                            !grid motion is allowed
         xfn=xnw(:,:,:) !mesh follows the network
!
! enforce substratum node positions while trying to preserve volume
! this is done because the penalty method for BC leads to very small
! but non-zero velocities.
         do is=1,ns 
            xfn(3,3,is)=xfn(3,3,is)-xfn(3,1,is)
            xfn(3,2,is)=xfn(3,2,is)-xfn(3,1,is)
            xfn(3,1,is)=0d0 
         enddo
!        print *,'xfn(3,3,547)',xfn(3,3,547)
!
         call edgadv(xfn,nadv)
!        print *,'post-edgadv xfn(3,3,547)',xfn(3,3,547)
!     deal with very small contact angles
!
         call edgslide(iflagv,iflagd,xfn)
!        print *,'post-edgslide xfn(3,3,547)',xfn(3,3,547)
         if (iflagd+iflagv.ne.0) call godriver(xfn)
         if (iflagd.ne.0) call edgadv(xfn,nadv)
!        print *,'post-edgadv2 xfn(3,3,547)',xfn(3,3,547)
!
         call rezdriver(xfn)!rezone the grid
         call godriver(xfn)
         irezfix=0
!
         do is=1,ns
            if (abs(xfn(3,1,is)).ge.1d-9) then
               print *,'xfn(3,1,is),is',xfn(3,1,is),is
               print *,'hvec(3,1,is),dt',hvec(3,1,is),dt
               stop
            endif
            if (xfn(3,3,is).lt.1d-9) then
              print *,'xfn(3,3,is),is',xfn(3,3,is),is
              print *,'xfn(3,3,is),is',xfn(3,3,is),is
              print *,'hvec(3,3,is),hvec(9,3,is),dt',
     1                 hvec(3,3,is),hvec(9,3,is),dt
              stop
            endif
            do il=1,3
               hvec(1,il,is)=xfn(1,il,is)
               hvec(2,il,is)=xfn(2,il,is)
               hvec(3,il,is)=xfn(3,il,is)
            enddo
         enddo
!
      endif
!
      return
      end subroutine avgridmo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avrz(xav,xrz,kqrz,xirz,etarz,zetarz)
!
! This subroutine finds the element and intrinsic coordinates
! with which to interpolate the node values of the final
! mesh with respect to the advected mesh 
!
      implicit none
!
      real(8) xav(3,3,NSM) ! coordinates of advected mesh
      real(8) xrz(3,3,NSM) ! coordinates of rezoned mesh
      integer kqrz(3,NSM)  ! interpolating element for rezoned nodes
      real(8) xirz(3,NSM)  ! xi value for rezoned nodes
      real(8) etarz(3,NSM) ! eta value for rezoned nodes
      real(8) zetarz(3,NSM)! zeta value for rezoned nodes
!
      real(8) xyz(3)
      real(8) xi,eta,zeta
      integer list(NSM)
      integer is,lv,kq,icorn
!
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
!--middle node (level=2)
!--prepare the list of neighboring elements
            list(1:kqos(is))=lqos(iqos(is):iqos(is)+kqos(is)-1)
            xyz=xrz(:,2,is)
            call find_element(xyz,xav,kqos(is),list,kq,xi,eta,zeta)
            if (kq.eq.0) then ! we look over the entire mesh
               call find_element(xyz,xav,nq,idxarr,kq,xi,eta,zeta)
               if (kq.eq.0) then
                  print *,'no interpolating element found for stack',
     1                    is,' level 2'
                  stop
               endif
            endif
            kqrz(2,is)=kq
            xirz(2,is)=xi
            etarz(2,is)=eta
            zetarz(2,is)=zeta
!
!--dorsal and ventral nodes (level=1 and 3)
            do lv=1,3,2
               xyz=xrz(:,lv,is)
               call find_surface(xyz,xav,kq,xi,eta,zeta)
               if (kq.eq.0) then
                  print *,'no interpolating element found for stack',
     1                    is,' level 1'
                  stop
               endif
               kqrz(lv,is)=kq
               xirz(lv,is)=xi
               etarz(lv,is)=eta
               zetarz(lv,is)=zeta
            enddo
         else !--edge stack case
!--contact line first (ventral node)
            xyz=xrz(:,1,is)
            call find_cl(xyz,xav,kq,xi,eta,zeta)
            kqrz(1,is)=kq
            xirz(1,is)=xi
            etarz(1,is)=eta
            zetarz(1,is)=zeta
!--middle and dorsal node
            do lv=2,3
               xyz=xrz(:,lv,is)
               call find_surface(xyz,xav,kq,xi,eta,zeta)
               kqrz(lv,is)=kq
               xirz(lv,is)=xi
               etarz(lv,is)=eta
               zetarz(lv,is)=zeta
            enddo
         endif
      enddo
!
      return
      end subroutine avrz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avfield(kom,oldopt) bind(C)
!
      implicit none
!
      integer(C_INT)::kom !svec field to advect
      character(kind=c_char,len=1),dimension(18),intent(in)::oldopt !choices are "eul" or "lag"
      character(len=2)::opt
!
      integer kqnw(3,NSM)  
      real(8) xinw(3,NSM), etanw(3,NSM), zetanw(3,NSM)
      integer kqca(3,NSM)  
      real(8) xica(3,NSM), etaca(3,NSM), zetaca(3,NSM)
!
      real(8) svtmp(3,NSM), svnew(3,NSM)
      integer nskind, is
      opt=" "
      loop_string: do i=1, 10
      if ( oldopt(i) == c_null_char ) then
         exit loop_string
      else
         opt (i:i) = oldopt(i)
      end if
      end do loop_string
!
      nskind=nint(dscp(kom,idspk))
      if (index(opt,'nw').gt.0) then !network advection
!--check if already computed the interpolating coordinates, if not, do it
         if (.not.avnw) then
            call avrz(xnw,xfn,kqnw,xinw,etanw,zetanw)
            avnw=.true.
         endif
!--apply dilation coefficient
         if (nskind.eq.1) then !this is a volume variable
            svtmp=svec(kom,:,:)*dilfacv
         elseif (nskind.eq.2) then
            svtmp=svec(kom,:,:)*dilfaca
         else 
            print *,'unrecognized variable type'
            print *,'kom=',kom,'   nskind=',nskind
            stop
         endif
!--interpolate
         call avinter(nskind,svtmp,kqnw,xinw,etanw,zetanw,svnew)
      elseif (index(opt,'ca').gt.0) then !bulk advection
         if (.not.avca) then
            call avrz(xca,xfn,kqca,xica,etaca,zetaca)
            avca=.true.
         endif
         if (nskind.eq.1.or.nskind.eq.2) then
            svtmp=svec(kom,:,:)
         else 
            print *,'unrecognized variable type'
            print *,'kom=',kom,'   nskind=',nskind
            stop
         endif
!--interpolate
         call avinter(nskind,svtmp,kqca,xica,etaca,zetaca,svnew)
      endif
      svec(kom,:,:)=svnew
!
      return
      end subroutine avfield
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avinter(nskind,svtmp,kqs,xis,etas,zetas,svnew)
!
      implicit none
!
      integer nskind
      real(8) svtmp(3,NSM)
      integer kqs(3,NSM)
      real(8) xis(3,NSM), etas(3,NSM), zetas(3,NSM)
      real(8) svnew(3,NSM)
!
      real(8) xi,eta,zeta
      real(8) dsdxi,dsdeta,dsdzeta
      real(8) snodes(3,4), s
      integer is,lv,i,l,kq,isn,il
!
      if (nskind.eq.1) then !volume variable
         do is=1,ns
            do lv=1,3
               kq=kqs(lv,is)
               xi=xis(lv,is)
               eta=etas(lv,is)
               zeta=zetas(lv,is)
               do i=1,4
                  isn=isoq(i,kq)
                  do l=1,3
                     snodes(l,i)=svtmp(l,isn)
                  enddo
               enddo
               call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
               svnew(lv,is)=s
            enddo
         enddo
      elseif (nskind.eq.2) then !surface variable
         do is=1,ns  !loop over dorsal and ventral nodes
            do lv=1,3,2
               kq=kqs(lv,is)
               xi=xis(lv,is)
               eta=etas(lv,is)
               zeta=zetas(lv,is)
               do i=1,4
                  isn=isoq(i,kq)
                  do l=1,3
                     snodes(l,i)=svtmp(l,isn)
                  enddo
               enddo
               call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
               svnew(lv,is)=s
            enddo
         enddo
         do il=1,nl  !loop over mid-point edge nodes
            is=isol(1,il)
            kq=kqs(2,is)
            xi=xis(2,is)
            eta=etas(2,is)
            zeta=zetas(2,is)
            do i=1,4
               isn=isoq(i,kq)
               do l=1,3
                  snodes(l,i)=svtmp(l,isn)
               enddo
            enddo
            call gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
            svnew(2,is)=s
         enddo
      else
         print *,'unable to advect nskind=',nskind
         stop
      endif
!
      return
      end subroutine avinter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine edgadv(xnew,nadv)
!
! this subroutine advances the edge nodes along the substrate
!
      implicit none
!
      real(8) xnew(3,3,NSM)
      integer nadv
! 
      real(8) xold(3,3,NSM)
      real(8) x1,y1,x2,y2,x3,y3,z2,z3,z2new,dx1,dy1
      real(8) dxpm,dypm,xn,yn
      real(8) zeta1, zeta2, H1, H2, H3
      integer il, is, isp, ism
!
      xold=xnew !backup the node positions
      nadv=0
      clfix=0
      do il=1,nl
         is=isol(1,il)
         z2=xold(3,2,is) 
         z3=xold(3,3,is) 
         if (z3.lt.0d0) then
! we are screwed because the dorsal point is below the substrate
            print *,'edgadv: dorsal point below substrate!',is,z3
            stop
         elseif (z2.lt.0d0) then !we are also screwed
            print *,'edgadv: middle point below substrate!',is,z2
            stop
         elseif (z3.gt.4d0*z2) then  !this should be the usual case
            x1=xold(1,1,is)         
            y1=xold(2,1,is) 
            x2=xold(1,2,is)         
            y2=xold(2,2,is) 
            x3=xold(1,3,is)         
            y3=xold(2,3,is) 
            zeta1=1d0/(1d0-0.5d0*z3/z2)
            H1=0.5d0*zeta1*(zeta1-1d0)
            H2=1d0-zeta1**2
            H3=0.5d0*zeta1*(zeta1+1d0)
!--displacement
            dx1=(H1*x1+H2*x2+H3*x3)-x1
            dy1=(H1*y1+H2*y2+H3*y3)-y1
!--check if the motion is outward or inward
            if ((dx1*cnn(1,il)+dy1*cnn(2,il)).ge.0d0) then
               xnew(1,1,is)=dx1+x1
               xnew(2,1,is)=dy1+y1
               print *,'edgadv: moved outward!',is,zeta1
               clfix(is)=1
            else
               xnew(1,1,is)=dx1+x1
               xnew(2,1,is)=dy1+y1
               print *,'edgadv: moved inward!',is,zeta1
               clfix(is)=1
            endif
!--in all cases move level 2 point such that z2 back to z3/4
            z2new=0.25d0*z3
            call findzeta(0d0,z2,z3,z2new,zeta2)
            H1=0.5d0*zeta2*(zeta2-1d0)
            H2=1d0-zeta2**2
            H3=0.5d0*zeta2*(zeta2+1d0)
            xnew(1,2,is)=H1*x1+H2*x2+H3*x3
            xnew(2,2,is)=H1*y1+H2*y2+H3*y3
            xnew(3,2,is)=H2*z2+H3*z3
            irezfix(is)=1
            nadv=nadv+1
         endif
      enddo
!         
      return 
      end subroutine edgadv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine findzeta(z1,z2,z3,z,zeta)
!
! Given z1 < z2 < z3 node positions (zeta=-1,0,1), find the intrinsic 
! coord zeta that gives position z
!
      implicit none
!
      real(8) z1, z2, z3, z, zeta
!
      real(8) tol
      parameter(tol=1d-6)
      real(8) zetad, zetau, zi
      integer i
!
      if (z.lt.z1.or.z.gt.z3.or.z2.lt.z1.or.z2.gt.z3) then 
         print *,'ERROR in findzeta',z1,z2,z3,z
         stop
      endif
!
      if (z.lt.z2) then
         zetad=-1d0
         zetau=0d0
      elseif (z.gt.z2) then
         zetad=0d0
         zetau=1d0
      else !z=z2 case
         zeta=0d0
         return
      endif
!
      do i=1,25
         zeta=0.5d0*(zetad+zetau)
         zi=0.5d0*zeta*(zeta-1d0)*z1+(1d0-zeta**2)*z2+
     1      0.5d0*zeta*(zeta+1d0)*z3
         if (abs(z-zi).lt.tol) then !close enough
            exit 
         elseif (zi.lt.z) then !zi too small, increase zetad
            zetad=zeta
         elseif (zi.gt.z) then !zi too big, decrease zetau
            zetau=zeta
         endif
      enddo
      if (i.ge.25) then
         print *,'failure to converge in findzeta' 
         print *,'z1,z2,z3,z,zi,zeta',z1,z2,z3,z,zi
         print *,'zeta,abs(z-zi)',zeta,abs(z-zi)
      endif
!
      return
      end subroutine findzeta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezdriver(xn)
!
! this subroutine is the driver for the rezoning of the mesh
!
      implicit none
!
      real(8) xn(3,3,NSM),cnode(3,NSM)
!
      real(8) reztol
      parameter(reztol=1d-3)
      integer icmx
      parameter(icmx=500)
      integer icyc,is,lv,iflag
      real(8) c_newton,dratio,volint,wtot1,wtot2,dwtot,dd
!
      cnode=1d0
!
!--rezone contact line (ventral edge nodes)
!
      do icyc=1,icmx
         call rezcl(xn,dd)
!        print *,'rezcl: dd, icyc',dd,icyc
         if (dd.lt.0.1d0*reztol.and.icyc.gt.3) exit
         call goquick(xn)
      enddo
      print *,'rezcl: dd, icyc',dd,icyc
      if (icyc.ge.icmx) print *,'rezcl warning',dd
      call godriver(xn)
!     call govolint(cnode,volint)
!     print *,'done rezcl, volint=',volint
!
!--rezone the middle and dorsal edge nodes
!
      call rezedg(xn)
!     print *,'post-rezedg xn(3,3,547)',xn(3,3,547)
      call godriver(xn)
!     call govolint(cnode,volint)
!     print *,'done rezedg, volint=',volint
!
!     deal with very small contact angles
!
!     call edgslide(iflag,xn)
!     if (iflag.ne.0) call godriver(xn)
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call godriver(xn)
        call goquick(xn)
        call rezbot(xn,c_newton,dratio,wtot2)
!       call godriver(xn)
!       call govolint(cnode,volint)
!       print *,'volint,wtot2,dratio,icyc',volint,wtot2,dratio,icyc
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
      print *,'rezbot: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) print *,'rezdriver warning',dratio
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call godriver(xn)
        call goquick(xn)
        call rezback(xn,c_newton,dratio,wtot2)
!       call godriver(xn)
!       call govolint(cnode,volint)
!       print *,'volint,wtot2,dratio,icyc',volint,wtot2,dratio,icyc
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
!     print *,'post-rezback xn(3,3,547)',xn(3,3,547)
      print *,'rezback: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) print *,'rezdriver warning',dratio
!
      c_newton=0.5d0 !under-relaxation 
      wtot1=0d0
      do icyc=1,icmx
        call goquick(xn)
        call rezone(xn,c_newton,dratio,wtot2)
!       call godriver(xn)
!       call govolint(cnode,volint)
!       print *,'volint,wtot2,dratio,icyc',volint,wtot2,dratio,icyc
        dwtot=abs(wtot2-wtot1)/wtot2
        if (dwtot.lt.reztol/nq.and.icyc.gt.3) exit
        wtot1=wtot2
      enddo
      print *,'rezone: wtot2, icyc',wtot2,icyc
      if (icyc.ge.icmx) print *,'rezdriver warning',dratio
!
      return
      end subroutine rezdriver 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezcl(xn,dd)
!
! This subroutine moves each ventral edge nodes toward the "midpoint"
! between its two neighbors:
! Let the original edge point be 0, the original neighbors A and B
! The target edge point O' is given by the intersection of
! 1) the substratum plane z=0
! 2) the volume conservation plane
!    (x-x0)*vnn_x+(y-y0)*vnn_y+(z-z0)*vnn_z=0
! 1) and 2) (z0=0) give the line (x-x0)*vnn_x+(y-y0)*vnn_y=0; z=0
! 3) the median between A and B
!       [x-(xA+xB)/2](xB-xA)+[y-(yA+yB)/2](yB-yA)=0
!
!  we move half-way to the target so xnew=(x+x0)/2 etc.
!
! User specified:
!   clfix(NLM), when non zero point is not moved
!   IMPORTANT: at the end of the routine clfix is zeroed out.
!
      implicit none
!
      real(8) xn(3,3,NSM) !in old, out new node x-y coordinates
      real(8) dd
!
      real(8) dxtg(2,NLM) !displ. of 2nd node of il in tg volume plane
      real(8) dxcl(2,NLM) !displ. of 2nd node of il in contact line
      real(8) x,y,xA,yA,xB,yB,dx,dy,det2,xM,yM,rhs1,rhs2,xt,yt,sint
      real(8) dAB,dxA,dyA,dxB,dyB,dA2,dB2,disp
      integer il,isA,is,isB,lv
!
      do il=1,nl !loop over edges
         isA=isol(1,il)
         is=isol(2,il)
         isB=isol(2,ilol(2,il))
!        if (is.ne.578) then
         if (is.ne.383) then
!        if (is.ne.965) then
!        if (clfix(is).eq.0) then
            x=xn(1,1,is)
            y=xn(2,1,is)
            xA=xn(1,1,isA)
            yA=xn(2,1,isA)
            xB=xn(1,1,isB)
            yB=xn(2,1,isB)
            dx=xB-xA
            dy=yB-yA
            dAB=sqrt(dx**2+dy**2)
            sint=sqrt(vnn(1,1,is)**2+vnn(2,1,is)**2)
            xM=0.5*(xA+xB)
            yM=0.5*(yA+yB)
!
!--work on the tangent plane displacement
            if (sint.lt.5d-2) then
! the tangent plane and the substratum essentialy coincide,
! (angle <= 3deg) we solve"
! (xnew-x)(yA-yB)+(ynew-y)(xB-xA)=0
! (xnew-xM)(xB-xA)+(ynew-yM)(yB-yA)=0
               det2=-dy**2-dx**2
               rhs1=-x*dy+y*dx
               rhs2=xM*dx+yM*dy
               xt=(rhs1*dy-rhs2*dx)/det2
               yt=(-rhs2*dy-rhs1*dx)/det2
!              print *,'xA,xB,x,xM,xnew',xA,xB,x,xM,0.5d0*(x+xt)
!              print *,'yA,yB,y,yM,ynew',yA,yB,y,yM,0.5d0*(y+yt)
               dxtg(1,il)=xt-x
               dxtg(2,il)=yt-y
!              print *,'il,n,sint',il,n,sint
            else
! solve:
!       xnew*(xB-xA)+ynew*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]
! xnew*vnn(1,1,is)+ynew*vnn(2,1,is)=x*vnn(1,1,is)+y*vnn(2,1,is)
               det2=dx*vnn(2,1,is)-dy*vnn(1,1,is)
!check that the plane/substrate intersection and the median are not //
               if (abs(det2).gt.2d-2*dAB) then
                  rhs1=xM*dx+yM*dy
                  rhs2=x*vnn(1,1,is)+y*vnn(2,1,is)
                  xt=(rhs1*vnn(2,1,is)-rhs2*dy)/det2
                  yt=(rhs2*dx-rhs1*vnn(1,1,is))/det2
                  dxtg(1,il)=xt-x
                  dxtg(2,il)=yt-y
               else
                  print *,'warning in rezcl: '
                  print *,'il,is,dAB,det2',il,is,dAB,det2
                  print *,'node number=',(is-1)*3+1
                  print *,'vn,vny,vnz',vnn(1:3,1,is)
                  print *,'dx,dy',dx,dy
                  print *,'sint',sint
                  print *,'return to continue'
                  read *
               endif
            endif
!
!--work on the contact line displacement
            dxA=x-xA
            dyA=y-yA
            dxB=x-xB
            dyB=y-yB
            dA2=dxA**2+dyA**2
            dB2=dxB**2+dyB**2
            if (dA2.lt.dB2) then !we move along the B line
!--compute intersection of B line with median
!   xt*(xB-xA)+yt*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]       
!   xt*(yB-y)+yt*(x-xB)=x*(yB-y)+y*(x-xB)
               det2=dx*dxB+dy*dyB
               if (abs(det2).lt.1d-16*dB2) then
                  print *,'error in rezcl'
                  print *,'det2,dB2',det2,dB2
                  print *,'is,il',is,il
                  stop
               endif
               rhs1=xM*dx+yM*dy
               rhs2=x*(yB-y)+y*(x-xB)
               xt=(rhs1*dxB-rhs2*dy)/det2
               yt=(rhs2*dx+rhs1*dyB)/det2
               dxcl(1,il)=xt-xn(1,1,is)
               dxcl(2,il)=yt-xn(2,1,is)
            elseif(dB2.lt.dA2) then!we move along the A line
!--compute intersection of B line with median
!   xt*(xB-xA)+yt*(yB-yA)=0.5*[(xA+xB)*(xB-xA)+(yA+yB)*(yB-yA)]       
!   xt*(yA-y)+yt*(x-xA)=x*(yA-y)+y*(x-xA)
               det2=dx*dxA+dy*dyA
               if (abs(det2).lt.1d-16*dA2) then
                  print *,'error in rezcl'
                  print *,'det2,dA2',det2,dA2
                  print *,'is,il',is,il
                  stop
               endif
               rhs1=xM*dx+yM*dy
               rhs2=x*(yA-y)+y*(x-xA)
               xt=(rhs1*dxA-rhs2*dy)/det2
               yt=(rhs2*dx+rhs1*dyA)/det2
               dxcl(1,il)=xt-x
               dxcl(2,il)=yt-y
            else
               dxcl(1,il)=0d0
               dxcl(2,il)=0d0
            endif
         else
            dxtg(1,il)=0d0
            dxtg(2,il)=0d0
            dxcl(1,il)=0d0
            dxcl(2,il)=0d0
         endif
      enddo
!
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=0d0
!           hvec(8,lv,is)=0d0
!           hvec(9,lv,is)=0d0
!        enddo
!     enddo
!     do il=1,nl
!        is=isol(2,il)
!        hvec(7,1,is)=0.2d0*(0.9d0*dxtg(1,il)+0.1d0*dxcl(1,il))
!        hvec(8,1,is)=0.2d0*(0.9d0*dxtg(2,il)+0.1d0*dxcl(2,il))
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!
! implement displacements
      dd=0d0
      do il=1,nl
         isA=isol(1,il)
         is=isol(2,il)
         isB=isol(2,ilol(2,il))
         dAB=sqrt((xn(1,1,isA)-xn(1,1,isB))**2+
     1            (xn(2,1,isA)-xn(2,1,isB))**2)
         dx=0.2d0*(0.9d0*dxtg(1,il)+0.1d0*dxcl(1,il))
         dy=0.2d0*(0.9d0*dxtg(2,il)+0.1d0*dxcl(2,il))
         dd=max(dd,sqrt(dx**2+dy**2)/dAB)
         xn(1,1,is)=xn(1,1,is)+dx
         xn(2,1,is)=xn(2,1,is)+dy
      enddo
!
      return
      end subroutine rezcl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezedg(xnew)
!
! this subroutine rezones the middle and dorsal edge nodes
!
!     the dorsal and middle nodes are shifted along their
!     tangent planes so as to lie in the plane perpendicular
!     to the substrate z=0 and containing contact lione normal.
!     The displacement keeps z constant
!
      implicit none
!
      real(8) xnew(3,3,NSM)
!
      real(8) xold(3,3,NSM)
      real(8) det2,rhs1,rhs2
      real(8) xv,yv,xo,yo,zo,zn,tnx,tny,tnz
      real(8) eiznx,eizny,eiznz,teizx,teizy,teizz,eizx,eizy,eizz
      real(8) xA,yA,zA,xB,yB,zB
      real(8) xe,ye,clnx,clny
      integer il,is

!--back up copy of xnew
      xold=xnew
!--loop over edge stacks
      do il=1,nl
         is=isol(1,il)
!--edge stack ventral coord (z=0)
         xe=xold(1,1,is)
         ye=xold(2,1,is)
!--normal to the contact line at isol(1,il)
         clnx=cnn(1,il)
         clny=cnn(2,il)
!--we slide the dorsal edge node along its tangent plane to the
!  plane perpendicular that contains the contact line normal.
!  We require that z remain constant
         xA=xold(1,3,is)
         yA=xold(2,3,is)
         zA=xold(3,3,is)
!--normal to the tangent plane
         tnx=vnn(1,3,is)
         tny=vnn(2,3,is)
         tnz=vnn(3,3,is)
! new position xB, yB, zB given by:
! clny*(xB-xe)-clnx*(yB-ye)=0 (in cl normal perpendicular plane) 
! (xB-xA)*tnx+(yB-yA)*tny+(zB-zA)*tnz=0 (in tangent plane)
! zB=zA
         det2=clny*tny+clnx*tnx
         if (abs(det2).gt.vtiny) then
            rhs1=clny*xe-clnx*ye
            rhs2=xA*tnx+yA*tny
            xnew(1,3,is)=(rhs1*tny+rhs2*clnx)/det2
            xnew(2,3,is)=(rhs2*clny-rhs1*tnx)/det2
            xnew(3,3,is)=zA
         endif
!
!--do the same for the middle edge node
         xA=xold(1,2,is)
         yA=xold(2,2,is)
         zA=xold(3,2,is)
!--normal to the tangent plane
         tnx=vnn(1,2,is)
         tny=vnn(2,2,is)
         tnz=vnn(3,2,is)
! new position xB, yB, zB given by:
! clny*(xB-xe)-clnx*(yB-ye)=0 (in cl normal perpendicular plane) 
! (xB-xA)*tnx+(yB-yA)*tny+(zB-zA)*tnz=0 (in tangent plane)
! zB=zA
         det2=clny*tny+clnx*tnx
         if (abs(det2).gt.vtiny) then
            rhs1=clny*xe-clnx*ye
            rhs2=xA*tnx+yA*tny
            xnew(1,2,is)=(rhs1*tny+rhs2*clnx)/det2
            xnew(2,2,is)=(rhs2*clny-rhs1*tnx)/det2
            xnew(3,2,is)=zA
         endif
      enddo
!
      return
      end subroutine rezedg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine edgslide(iflagv,iflagd,xnew)
!
      implicit none
!
      integer iflagv,iflagd
      real(8) xnew(3,3,NSM)
!
      real(8) a_out, a_in
      parameter(a_out=5d0*3.141592d0/180d0) !5 degrees
      parameter(a_in=45d0*3.141592d0/180d0) !45 degrees
      real(8) dxy12,dz12,dx12,dy12,dxy13,dz13,dx1,dy1,dx2,dy2,dx3,dy3
      real(8) dx13,dy13,ddxy13
      real(8) dV1,dV2,dV3,dV21,alpha,alpha3,tanalph,a3min,ddotn
      real(8) Bmat(3,3),rhs(3),Bvec(3)
      real(8) d2,d3,clnx,clny,c2,c3,z2,z3,det
      integer il,is
!
      iflagv=0
      do il=1,nl
         is=isol(1,il)
         dx12=xnew(1,2,is)-xnew(1,1,is)
         dy12=xnew(2,2,is)-xnew(2,1,is)
         ddotn=dx12*cnn(1,il)+dy12*cnn(2,il)
         dxy12=sqrt(dx12**2+dy12**2)
         dz12=xnew(3,2,is)-xnew(3,1,is)
         alpha=atan(dz12/max(dxy12,vtiny))
         if (alpha.lt.a_out.and.ddotn.gt.0d0) then
!           print *,'old alpha',alpha*180d0/3.141592
            if (abs(dx12).gt.1d-6*dxy12) then
               dx1=0.5d0*(dxy12**2-(dz12/tan(a_out))**2)/
     1                     (dx12+dy12**2/dx12)
               dy1=dx1*dy12/dx12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=-(dV1/(vnn(0,2,is)))/
     1              (vnn(1,2,is)+vnn(2,2,is)*dy12/dx12)
               dy2=dx2*dy12/dx12
            else
               dx1=0d0
               dy1=0.5d0*(dxy12**2-(dz12/tan(a_out))**2)/dy12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=0d0
               dy2=-(dV1/(vnn(0,2,is)))/vnn(2,2,is)
            endif
            dV2=vnn(0,2,is)*(vnn(1,2,is)*dx2+vnn(2,2,is)*dy2)
            xnew(1,1,is)=xnew(1,1,is)+dx1
            xnew(2,1,is)=xnew(2,1,is)+dy1
            xnew(1,2,is)=xnew(1,2,is)+dx2
            xnew(2,2,is)=xnew(2,2,is)+dy2
            iflagv=iflagv+1
!           print *,'slid bottom edge node is',is,alpha
!           print *,'vnn(0:3,1,is)',vnn(0:3,1,is)
!           print *,'dV1,dV2',dV1,dV2
!           print *,'dx1,dx2,dy1,dy2',dx1,dx2,dy1,dy2
!           dx12=xnew(1,2,is)-xnew(1,1,is)
!           dy12=xnew(2,2,is)-xnew(2,1,is)
!           dxy12=sqrt(dx12**2+dy12**2)
!           dz12=xnew(3,2,is)-xnew(3,1,is)
!           alpha=atan(dz12/max(dxy12,vtiny))
!           print *,'new alpha',alpha*180d0/3.141592
         elseif (alpha.lt.a_in.and.ddotn.lt.0d0) then
!           print *,'old alpha',alpha*180d0/3.141592
            if (abs(dx12).gt.1d-6*dxy12) then
               dx1=0.5d0*(dxy12**2-(dz12/tan(a_in))**2)/
     1                     (dx12+dy12**2/dx12)
               dy1=dx1*dy12/dx12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=-(dV1/(vnn(0,2,is)))/
     1              (vnn(1,2,is)+vnn(2,2,is)*dy12/dx12)
               dy2=dx2*dy12/dx12
            else
               dx1=0d0
               dy1=0.5d0*(dxy12**2-(dz12/tan(a_in))**2)/dy12
               dV1=vnn(0,1,is)*(vnn(1,1,is)*dx1+vnn(2,1,is)*dy1)
               dx2=0d0
               dy2=-(dV1/(vnn(0,2,is)))/vnn(2,2,is)
            endif
            dV2=vnn(0,2,is)*(vnn(1,2,is)*dx2+vnn(2,2,is)*dy2)
            xnew(1,1,is)=xnew(1,1,is)+dx1
            xnew(2,1,is)=xnew(2,1,is)+dy1
            xnew(1,2,is)=xnew(1,2,is)+dx2
            xnew(2,2,is)=xnew(2,2,is)+dy2
            iflagv=iflagv+1
!           print *,'slid in bottom edge node is',is,alpha
!           print *,'vnn(0:3,1,is)',vnn(0:3,1,is)
!           print *,'dV1,dV2',dV1,dV2
!           print *,'dx1,dx2,dy1,dy2',dx1,dx2,dy1,dy2
!           dx12=xnew(1,2,is)-xnew(1,1,is)
!           dy12=xnew(2,2,is)-xnew(2,1,is)
!           dxy12=sqrt(dx12**2+dy12**2)
!           dz12=xnew(3,2,is)-xnew(3,1,is)
!           alpha=atan(dz12/max(dxy12,vtiny))
!           print *,'new alpha',alpha*180d0/3.141592
!           print *,'return to continue'
!           read  *
         endif   
      enddo
!
      if (iflagv.ne.0) then
         print *,iflagv,' ventral edge nodes were slid'
         call goquick(xnew)
      endif
!
       iflagd=0
!      a3min=vbig
!      do il=1,nl
!         is=isol(1,il)
!         dx13=xnew(1,3,is)-xnew(1,1,is)
!         dy13=xnew(2,3,is)-xnew(2,1,is)
!         ddxy13=dx13**2+dy13**2
!         dxy13=sqrt(ddxy13)
!         z3=xnew(3,3,is)
!         tanalph=z3/max(dxy12,vtiny)
!         alpha3=atan(tanalph)
!         a3min=min(alpha3,a3min)
!         if (alpha3.lt.acrit) then
!! solve:
!! nx*dx+ny*dy+nz*dz=0 (in volume tangent plane)
!! dy13*dx-dx13*dy=0 (in 13 direction on the xy plane)
!! -dx13*z3*dx-y13*z3*dy+ddxy13*dz=(tancrit-tan)*dxy13**3/2
!            Bmat(1,1)=vnn(1,3,is)
!            Bmat(1,2)=vnn(2,3,is)
!            Bmat(1,3)=vnn(3,3,is)
!            rhs(1)=0d0
!            Bmat(2,1)=dy13
!            Bmat(2,2)=-dx13
!            Bmat(2,3)=0d0
!            rhs(2)=0d0
!            Bmat(3,1)=-dx13*z3
!            Bmat(3,2)=-dy13*z3
!            Bmat(3,3)=ddxy13
!            rhs(3)=(tan(2d0*acrit)-tanalph)*ddxy13*dxy13
!            call mat3solv(Bmat,rhs,Bvec)
!            print *,'is',is
!            print *,'xold',xnew(1:3,3,is)
!            xnew(1,3,is)=xnew(1,3,is)+Bvec(1)
!            xnew(2,3,is)=xnew(2,3,is)+Bvec(2)
!            xnew(3,3,is)=xnew(3,3,is)+Bvec(3)
!            print *,'xnew',xnew(1:3,3,is)
!            iflagd=iflagd+1
!         endif   
!      enddo
!      if (iflagd.ne.0) then
!         print *,iflagd,' dorsal edge nodes were slid'
!         call goquick(xnew)
!      endif
!      print *,'edgslide: a3min',a3min*180d0/3.141592
      
!
      return
      end subroutine edgslide
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezbot(xn,c_newton,dratio,wtot)
!
! This subroutine rezones the nodes
! so as to minimize the Winslow functional
!
      implicit none
!
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
      real(8) dratio
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dS1,dS2
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) B(11,11),BLU(11,11),dis(11),vrhs(11)
      integer ipvt(11), isng
      real(8) z2min,z3min
      integer is2min,is3min
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Win=0d0
      Wx=0d0
      Wy=0d0
      Wz=0d0
      Wxx=0d0
      Wyy=0d0
      Wzz=0d0
      Wxy=0d0
      Wxz=0d0
      Wyz=0d0
      dxn=0d0
!
      dz1=0d0
      dz2=0d0
      do iq=1,nq !outer loop over dorsal quads
         do isg=1,4 !loop over the gauss points of iq
            dx1=dxSv(1,isg,iq)
            dx2=dxSv(2,isg,iq)
            dy1=dySv(1,isg,iq)
            dy2=dySv(2,isg,iq)
            call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,
     1                   g11,g22,F)
            Jcb=dAv(0,isg,iq)
            dJ=1d0/Jcb
            Wg=F*dJ
            Wtot=Wtot+Wg
            do isn=1,4
               istack=isoq(isn,iq)
               dS1=dS4(1,isn,isg)
               dS2=dS4(2,isn,isg)
             call dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                        Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
               Wx(1,istack)=Wx(1,istack)+Wxi
               Wy(1,istack)=Wy(1,istack)+Wyi
               Wxx(1,istack)=Wxx(1,istack)+Wxxi
               Wyy(1,istack)=Wyy(1,istack)+Wyyi
               Wxy(1,istack)=Wxy(1,istack)+Wxyi
            enddo
         enddo
      enddo
!
!     print *,'rezbot: Wtot',Wtot
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation)
!
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--ventral nodes
            gradWtg1=-Wx(1,is)
            gradWtg2=-Wy(1,is)
            HessWtg11=Wxx(1,is)
            HessWtg22=Wyy(1,is)
            HessWtg12=Wxy(1,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,1,is)=c_newton*dtg1
            dxn(2,1,is)=c_newton*dtg2
         else!an edge stack
            dxn(1,1,is)=0d0
            dxn(2,1,is)=0d0
         endif 
      enddo
!
! for debugging might want a dump
!     if (mod(kount,50).eq.0) then
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=dxn(1,lv,is)
!           hvec(8,lv,is)=dxn(2,lv,is)
!           hvec(9,lv,is)=dxn(3,lv,is)
!        enddo
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!     read *
!     hvec(7:9,:,1:ns)=0d0
!     endif
!     kount=kount+1
!
!--update the node positions while limiting the displacements
!
      dratio=0d0
      dmean=0d0
      dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxn(1,1,is)**2+dxn(2,1,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dscale=abs(voln(1,is))**(1./3.)
            dratio=max(dratio,dist/dscale)
            if (dist.gt.deps*dscale) then
               fac=deps*dist/dscale
               dxn(1,1,is)=dxn(1,1,is)*fac
               dxn(2,1,is)=dxn(2,1,is)*fac
            endif
         endif
      enddo
      dmean=dmean/ns
!     print *,'dmax,dmean',dmax,dmean
!
      do is=1,ns
         do lv=1,3
            xn(1,lv,is)=xn(1,lv,is)+dxn(1,lv,is)
            xn(2,lv,is)=xn(2,lv,is)+dxn(2,lv,is)
         enddo
      enddo
!
      return
      end subroutine rezbot
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezback(xn,c_newton,dratio,wtot)
!
! This subroutine rezones the nodes
! so as to minimize the Winslow functional
!
      implicit none
!
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
      real(8) dratio
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dS1,dS2
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) B(11,11),BLU(11,11),dis(11),vrhs(11)
      integer ipvt(11), isng
      real(8) z2min,z3min
      integer is2min,is3min
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Win=0d0
      Wx=0d0
      Wy=0d0
      Wz=0d0
      Wxx=0d0
      Wyy=0d0
      Wzz=0d0
      Wxy=0d0
      Wxz=0d0
      Wyz=0d0
      dxn=0d0
!
      do iq=1,nq !outer loop over dorsal quads
         do isg=1,4 !loop over the gauss points of iq
            dx1=dxSd(1,isg,iq)
            dx2=dxSd(2,isg,iq)
            dy1=dySd(1,isg,iq)
            dy2=dySd(2,isg,iq)
            dz1=dzSd(1,isg,iq)
            dz2=dzSd(2,isg,iq)
            call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,
     1                   g11,g22,F)
            Jcb=dAd(0,isg,iq)
            dJ=1d0/Jcb
            Wg=F*dJ
            Wtot=Wtot+Wg
            do isn=1,4
               istack=isoq(isn,iq)
               dS1=dS4(1,isn,isg)
               dS2=dS4(2,isn,isg)
             call dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                        Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
               Wx(3,istack)=Wx(3,istack)+Wxi
               Wy(3,istack)=Wy(3,istack)+Wyi
               Wz(3,istack)=Wz(3,istack)+Wzi
               Wxx(3,istack)=Wxx(3,istack)+Wxxi
               Wyy(3,istack)=Wyy(3,istack)+Wyyi
               Wzz(3,istack)=Wzz(3,istack)+Wzzi
               Wxy(3,istack)=Wxy(3,istack)+Wxyi
               Wxz(3,istack)=Wxz(3,istack)+Wxzi
               Wyz(3,istack)=Wyz(3,istack)+Wyzi
            enddo
         enddo
      enddo
!
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               dx1=dxSe(1,lvg,isg,il)
               dx2=dxSe(2,lvg,isg,il)
               dy1=dySe(1,lvg,isg,il)
               dy2=dySe(2,lvg,isg,il)
               dz1=dzSe(1,lvg,isg,il)
               dz2=dzSe(2,lvg,isg,il)
               call Winslow(dx1,dx2,dy1,dy2,dz1,dz2,
     1                      g11,g22,F)
               Jcb=dAe(0,lvg,isg,il)
               dJ=1d0/Jcb
               Wg=F*dJ*wgp(lvg)
               Wtot=Wtot+Wg
               do isn=1,2
                  istack=isol(isn,il)
                  do lvn=1,3
                     dS1=dS6(1,lvn,isn,lvg,isg)
                     dS2=dS6(2,lvn,isn,lvg,isg)
                     call dWin(dx1,dx2,dy1,dy2,dz1,dz2,
     1                         g11,g22,F,Jcb,dJ,dS1,dS2,
     3                         Wxi,Wyi,Wzi,
     3                         Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
                     Wx(lvn,istack)=Wx(lvn,istack)+Wxi*wgp(lvg)
                     Wy(lvn,istack)=Wy(lvn,istack)+Wyi*wgp(lvg)
                     Wz(lvn,istack)=Wz(lvn,istack)+Wzi*wgp(lvg)
                     Wxx(lvn,istack)=Wxx(lvn,istack)+Wxxi*wgp(lvg)
                     Wyy(lvn,istack)=Wyy(lvn,istack)+Wyyi*wgp(lvg)
                     Wzz(lvn,istack)=Wzz(lvn,istack)+Wzzi*wgp(lvg)
                     Wxy(lvn,istack)=Wxy(lvn,istack)+Wxyi*wgp(lvg)
                     Wxz(lvn,istack)=Wxz(lvn,istack)+Wxzi*wgp(lvg)
                     Wyz(lvn,istack)=Wyz(lvn,istack)+Wyzi*wgp(lvg)
                  enddo
               enddo
            enddo
         enddo
      enddo

!     print *,'rezback: Wtot',Wtot
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation)
!
      z2min=vbig
      z3min=vbig
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--ventral and middle nodes do not move here
            dxn(1,1,is)=0d0
            dxn(2,1,is)=0d0
            dxn(3,1,is)=0d0
            dxn(1,2,is)=0d0
            dxn(2,2,is)=0d0
            dxn(3,2,is)=0d0
!--dorsal nodes
            u1x=utg1(1,3,is)
            u1y=utg1(2,3,is)
            u1z=utg1(3,3,is)
            u2x=utg2(1,3,is)
            u2y=utg2(2,3,is)
            u2z=utg2(3,3,is)
            gradWtg1=-Wx(3,is)*u1x-Wy(3,is)*u1y-Wz(3,is)*u1z
            gradWtg2=-Wx(3,is)*u2x-Wy(3,is)*u2y-Wz(3,is)*u2z
            HessWtg11=+u1x*u1x*Wxx(3,is)+u1y*u1y*Wyy(3,is)
     1                +u1z*u1z*Wzz(3,is)
     2           +2d0*(u1x*u1y*Wxy(3,is)+u1x*u1z*Wxz(3,is)+
     3                 u1y*u1z*Wyz(3,is))
            HessWtg22=+u2x*u2x*Wxx(3,is)+u2y*u2y*Wyy(3,is)
     1                +u2z*u2z*Wzz(3,is)
     2           +2d0*(u2x*u2y*Wxy(3,is)+u2x*u2z*Wxz(3,is)+
     3                 u2y*u2z*Wyz(3,is))
            HessWtg12=+u1x*u2x*Wxx(3,is)+u1y*u2y*Wyy(3,is)
     1                +u1z*u2z*Wzz(3,is)+
     2                (u1x*u2y+u1y*u2x)*Wxy(3,is)+
     3                (u1x*u2z+u1z*u2x)*Wxz(3,is)+
     4                (u1y*u2z+u1z*u2y)*Wyz(3,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,3,is)=c_newton*(dtg1*u1x+dtg2*u2x)
            dxn(2,3,is)=c_newton*(dtg1*u1y+dtg2*u2y)
            dxn(3,3,is)=c_newton*(dtg1*u1z+dtg2*u2z)
         else !an edge stack
!-ventral node pinned
            dxn(1,1,is)=0d0
            dxn(2,1,is)=0d0
            dxn(3,1,is)=0d0
!-middle and dorsal node move in the intersection
! of the tangent plane and the CL normal plane
! with the constraint that dz3/z3=dz2/z2
!
!--Lagrange multiplier method
!
            il=ilos(1,is)
            B=0d0
            B(1,1)=Wxx(2,is); B(1,2)=Wxy(2,is); B(1,3)=Wxz(2,is)
            B(2,1)=Wxy(2,is); B(2,2)=Wyy(2,is); B(2,3)=Wyz(2,is)
            B(3,1)=Wxz(2,is); B(3,2)=Wyz(2,is); B(3,3)=Wzz(2,is)
            B(1,4)=vnn(1,2,is);B(2,4)=vnn(2,2,is);B(3,4)=vnn(3,2,is)
            B(4,1)=vnn(1,2,is);B(4,2)=vnn(2,2,is);B(4,3)=vnn(3,2,is)
            B(1,5)=cnn(2,il);B(2,5)=-cnn(1,il)
            B(5,1)=cnn(2,il);B(5,2)=-cnn(1,il)
            B(6,6)=Wxx(3,is); B(6,7)=Wxy(3,is); B(6,8)=Wxz(3,is)
            B(7,6)=Wxy(3,is); B(7,7)=Wyy(3,is); B(7,8)=Wyz(3,is)
            B(8,6)=Wxz(3,is); B(8,7)=Wyz(3,is); B(8,8)=Wzz(3,is)
            B(6,9)=vnn(1,3,is);B(7,9)=vnn(2,3,is);B(8,9)=vnn(3,3,is)
            B(9,6)=vnn(1,3,is);B(9,7)=vnn(2,3,is);B(9,8)=vnn(3,3,is)
            B(6,10)=cnn(2,il);B(7,10)=-cnn(1,il)
            B(10,6)=cnn(2,il);B(10,7)=-cnn(1,il)
            B(11,3)=xn(3,3,is);B(11,8)=-xn(3,2,is)
            B(3,11)=xn(3,3,is);B(8,11)=-xn(3,2,is)
            call lufa(B,11,BLU,isng,ipvt)
            if (isng.ne.0) print *,
     1         'ill conditioned matrix in rezone, isng=',isng
            vrhs=0d0
            vrhs(1)=-Wx(2,is); vrhs(2)=-Wy(2,is); vrhs(3)=-Wz(2,is)
            vrhs(6)=-Wx(3,is); vrhs(7)=-Wy(3,is); vrhs(8)=-Wz(3,is)
            call lusl(BLU,vrhs,isng,ipvt,11,dis)
            dxn(1,2,is)=c_newton*dis(1)
            dxn(2,2,is)=c_newton*dis(2)
            dxn(3,2,is)=c_newton*dis(3)
            dxn(1,3,is)=c_newton*dis(6)
            dxn(2,3,is)=c_newton*dis(7)
            dxn(3,3,is)=c_newton*dis(8)
!           print *,'is,il',is,il
!           print *,'dxn(1:3,3,is)',dxn(1:3,3,is)
!           print *,'dxn(1:3,2,is)',dxn(1:3,2,is)
!           print *,'dxn3*vnn3',dxn(1,3,is)*vnn(1,3,is)+
!    1                          dxn(2,3,is)*vnn(2,3,is)+
!    2                          dxn(3,3,is)*vnn(3,3,is)
!           print *,'dxn2*vnn2',dxn(1,2,is)*vnn(1,2,is)+
!    1                          dxn(2,2,is)*vnn(2,2,is)+
!    2                          dxn(3,2,is)*vnn(3,2,is)
!           print *,'dxn3*cl',dxn(1,3,is)*cnn(2,il)-
!    1                        dxn(2,3,is)*cnn(1,il)
!           print *,'dxn2*cl',dxn(1,2,is)*cnn(2,il)-
!    1                        dxn(2,2,is)*cnn(1,il)
            if (xn(3,3,is).lt.z3min) then
               is3min=is
               z3min=xn(3,3,is)
            endif
            if (xn(3,2,is).lt.z2min) then
               is2min=is
               z2min=xn(3,2,is)
            endif
         endif 
      enddo
!     print *,'is3min,z,dz',is3min,xn(3,3,is3min),dxn(3,3,is3min)
!     print *,'is3min,z0,dzv',is3min,hvec(3,3,is3min),
!    1                        hvec(9,3,is3min)*tstp
!     print *,'is2min,z,dz',is2min,xn(3,2,is2min),dxn(3,2,is2min)
!
! for debugging might want a dump
!     if (mod(kount,50).eq.0) then
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=dxn(1,lv,is)
!           hvec(8,lv,is)=dxn(2,lv,is)
!           hvec(9,lv,is)=dxn(3,lv,is)
!        enddo
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!     stop
!     hvec(7:9,:,1:ns)=0d0
!     endif
!     kount=kount+1
!
!--update the node positions while limiting the displacements
!
!     print *,'dxn(3,2,547)',xn(3,2,547),dxn(3,2,547)
!     print *,'dxn(3,3,547)',xn(3,3,547),dxn(3,3,547)
      dratio=0d0
      dmean=0d0
      dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dscale=abs(voln(3,is))**(1./3.)
            dratio=max(dratio,dist/dscale)
            if (dist.gt.deps*dscale) then
               fac=deps*dscale/dist
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
            endif
         else !an edge stack
            dscale=min(abs(voln(3,is))**(1./3.),xn(3,2,is),xn(3,3,is))
            dist=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dratio=max(dratio,dist/dscale)
            if (dist.gt.(deps*dscale)) then
               fac=deps*dscale/dist
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
               dxn(1,2,is)=dxn(1,2,is)*fac
               dxn(2,2,is)=dxn(2,2,is)*fac
               dxn(3,2,is)=dxn(3,2,is)*fac
            endif
         endif
      enddo
      dmean=dmean/(3*ns)
!     print *,'dmax,dmean',dmax,dmean
!
      do is=1,ns
         do lv=1,3
            xn(1,lv,is)=xn(1,lv,is)+dxn(1,lv,is)
            xn(2,lv,is)=xn(2,lv,is)+dxn(2,lv,is)
            xn(3,lv,is)=xn(3,lv,is)+dxn(3,lv,is)
         enddo
      enddo
!
      return
      end subroutine rezback
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezone(xn,c_newton,dratio,wtot)
!
! This subroutine rezones the nodes
! so as to minimize the Winslow functional
!
      implicit none
!
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
      real(8) dratio
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dH1,dH2,dH3
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ2, dJ3
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) B(11,11),BLU(11,11),dis(11),vrhs(11)
      integer ipvt(11), isng
      real(8) z2min,z3min
      integer is2min,is3min
      real(8) dmax,dmean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Win=0d0
      Wx=0d0
      Wy=0d0
      Wz=0d0
      Wxx=0d0
      Wyy=0d0
      Wzz=0d0
      Wxy=0d0
      Wxz=0d0
      Wyz=0d0
      dxn=0d0
!
      do iq=1,nq !outer loop over iq quadrilateral
         do isg=1,4 !loop over the gauss points of iq
            do lvg=1,3
               dx1=dxV(1,lvg,isg,iq)
               dx2=dxV(2,lvg,isg,iq)
               dx3=dxV(3,lvg,isg,iq)
               dy1=dyV(1,lvg,isg,iq)
               dy2=dyV(2,lvg,isg,iq)
               dy3=dyV(3,lvg,isg,iq)
               dz1=dzV(1,lvg,isg,iq)
               dz2=dzV(2,lvg,isg,iq)
               dz3=dzV(3,lvg,isg,iq)
               call TTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                  g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,Wg)
               Win(iq)=Win(iq)+wgp(lvg)*Wg
               Wtot=Wtot+wgp(lvg)*Wg
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     dH1=dH(1,lvn,isn,lvg,isg)
                     dH2=dH(2,lvn,isn,lvg,isg)
                     dH3=dH(3,lvn,isn,lvg,isg)
                     call dTTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                         g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,
     2                         dH1,dH2,dH3,
     3                         Wxi,Wyi,Wzi,
     4                         Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
                     Wx(lvn,istack)=Wx(lvn,istack)+Wxi*wgp(lvg)
                     Wy(lvn,istack)=Wy(lvn,istack)+Wyi*wgp(lvg)
                     Wz(lvn,istack)=Wz(lvn,istack)+Wzi*wgp(lvg)
                     Wxx(lvn,istack)=Wxx(lvn,istack)+Wxxi*wgp(lvg)
                     Wyy(lvn,istack)=Wyy(lvn,istack)+Wyyi*wgp(lvg)
                     Wzz(lvn,istack)=Wzz(lvn,istack)+Wzzi*wgp(lvg)
                     Wxy(lvn,istack)=Wxy(lvn,istack)+Wxyi*wgp(lvg)
                     Wxz(lvn,istack)=Wxz(lvn,istack)+Wxzi*wgp(lvg)
                     Wyz(lvn,istack)=Wyz(lvn,istack)+Wyzi*wgp(lvg)
                  enddo
               enddo
            enddo
         enddo
      enddo
!     print *,'Wtot',Wtot
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation)
!
      z2min=vbig
      z3min=vbig
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--ventral nodes are restricted to move in the z=0 plane
            gradWtg1=-Wx(1,is)
            gradWtg2=-Wy(1,is)
            HessWtg11=Wxx(1,is)
            HessWtg22=Wyy(1,is)
            HessWtg12=Wxy(1,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,1,is)=c_newton*dtg1
            dxn(2,1,is)=c_newton*dtg2
            dxn(1,1,is)=0d0
            dxn(2,1,is)=0d0
            dxn(3,1,is)=0d0
!--middle nodes
            gradW(1)=-Wx(2,is)
            gradW(2)=-Wy(2,is)
            gradW(3)=-Wz(2,is)
            HessW(1,1)=Wxx(2,is)
            HessW(2,1)=Wxy(2,is)
            HessW(3,1)=Wxz(2,is)
            HessW(1,2)=Wxy(2,is)
            HessW(2,2)=Wyy(2,is)
            HessW(3,2)=Wyz(2,is)
            HessW(1,3)=Wxz(2,is)
            HessW(2,3)=Wyz(2,is)
            HessW(3,3)=Wzz(2,is)
            call mat3solv(HessW,gradW,dvec)
            dxn(1,2,is)=c_newton*dvec(1)
            dxn(2,2,is)=c_newton*dvec(2)
            dxn(3,2,is)=c_newton*dvec(3)
!--dorsal nodes
            u1x=utg1(1,3,is)
            u1y=utg1(2,3,is)
            u1z=utg1(3,3,is)
            u2x=utg2(1,3,is)
            u2y=utg2(2,3,is)
            u2z=utg2(3,3,is)
            gradWtg1=-Wx(3,is)*u1x-Wy(3,is)*u1y-Wz(3,is)*u1z
            gradWtg2=-Wx(3,is)*u2x-Wy(3,is)*u2y-Wz(3,is)*u2z
            HessWtg11=+u1x*u1x*Wxx(3,is)+u1y*u1y*Wyy(3,is)
     1                +u1z*u1z*Wzz(3,is)
     2           +2d0*(u1x*u1y*Wxy(3,is)+u1x*u1z*Wxz(3,is)+
     3                 u1y*u1z*Wyz(3,is))
            HessWtg22=+u2x*u2x*Wxx(3,is)+u2y*u2y*Wyy(3,is)
     1                +u2z*u2z*Wzz(3,is)
     2           +2d0*(u2x*u2y*Wxy(3,is)+u2x*u2z*Wxz(3,is)+
     3                 u2y*u2z*Wyz(3,is))
            HessWtg12=+u1x*u2x*Wxx(3,is)+u1y*u2y*Wyy(3,is)
     1                +u1z*u2z*Wzz(3,is)+
     2                (u1x*u2y+u1y*u2x)*Wxy(3,is)+
     3                (u1x*u2z+u1z*u2x)*Wxz(3,is)+
     4                (u1y*u2z+u1z*u2y)*Wyz(3,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,3,is)=c_newton*(dtg1*u1x+dtg2*u2x)
            dxn(2,3,is)=c_newton*(dtg1*u1y+dtg2*u2y)
            dxn(3,3,is)=c_newton*(dtg1*u1z+dtg2*u2z)
            dxn(1,3,is)=0d0
            dxn(2,3,is)=0d0
            dxn(3,3,is)=0d0
         else !an edge stack
!-ventral node pinned
            dxn(1,1,is)=0d0
            dxn(2,1,is)=0d0
            dxn(3,1,is)=0d0
!-middle and dorsal node move in the intersection
! of the tangent plane and the CL normal plane
! with the constraint that dz3/z3=dz2/z2
!
!--Lagrange multiplier method
!
            il=ilos(1,is)
            B=0d0
            B(1,1)=Wxx(2,is); B(1,2)=Wxy(2,is); B(1,3)=Wxz(2,is)
            B(2,1)=Wxy(2,is); B(2,2)=Wyy(2,is); B(2,3)=Wyz(2,is)
            B(3,1)=Wxz(2,is); B(3,2)=Wyz(2,is); B(3,3)=Wzz(2,is)
            B(1,4)=vnn(1,2,is);B(2,4)=vnn(2,2,is);B(3,4)=vnn(3,2,is)
            B(4,1)=vnn(1,2,is);B(4,2)=vnn(2,2,is);B(4,3)=vnn(3,2,is)
            B(1,5)=cnn(2,il);B(2,5)=-cnn(1,il)
            B(5,1)=cnn(2,il);B(5,2)=-cnn(1,il)
            B(6,6)=Wxx(3,is); B(6,7)=Wxy(3,is); B(6,8)=Wxz(3,is)
            B(7,6)=Wxy(3,is); B(7,7)=Wyy(3,is); B(7,8)=Wyz(3,is)
            B(8,6)=Wxz(3,is); B(8,7)=Wyz(3,is); B(8,8)=Wzz(3,is)
            B(6,9)=vnn(1,3,is);B(7,9)=vnn(2,3,is);B(8,9)=vnn(3,3,is)
            B(9,6)=vnn(1,3,is);B(9,7)=vnn(2,3,is);B(9,8)=vnn(3,3,is)
            B(6,10)=cnn(2,il);B(7,10)=-cnn(1,il)
            B(10,6)=cnn(2,il);B(10,7)=-cnn(1,il)
            B(11,3)=xn(3,3,is);B(11,8)=-xn(3,2,is)
            B(3,11)=xn(3,3,is);B(8,11)=-xn(3,2,is)
            call lufa(B,11,BLU,isng,ipvt)
            if (isng.ne.0) print *,
     1         'ill conditioned matrix in rezone, isng=',isng
            vrhs=0d0
            vrhs(1)=-Wx(2,is); vrhs(2)=-Wy(2,is); vrhs(3)=-Wz(2,is)
            vrhs(6)=-Wx(3,is); vrhs(7)=-Wy(3,is); vrhs(8)=-Wz(3,is)
            call lusl(BLU,vrhs,isng,ipvt,11,dis)
            dxn(1,2,is)=c_newton*dis(1)
            dxn(2,2,is)=c_newton*dis(2)
            dxn(3,2,is)=c_newton*dis(3)
            dxn(1,3,is)=c_newton*dis(6)
            dxn(2,3,is)=c_newton*dis(7)
            dxn(3,3,is)=c_newton*dis(8)
            dxn(1,2,is)=0d0
            dxn(2,2,is)=0d0
            dxn(3,2,is)=0d0
            dxn(1,3,is)=0d0
            dxn(2,3,is)=0d0
            dxn(3,3,is)=0d0
!           print *,'is,il',is,il
!           print *,'dxn(1:3,3,is)',dxn(1:3,3,is)
!           print *,'dxn(1:3,2,is)',dxn(1:3,2,is)
!           print *,'dxn3*vnn3',dxn(1,3,is)*vnn(1,3,is)+
!    1                          dxn(2,3,is)*vnn(2,3,is)+
!    2                          dxn(3,3,is)*vnn(3,3,is)
!           print *,'dxn2*vnn2',dxn(1,2,is)*vnn(1,2,is)+
!    1                          dxn(2,2,is)*vnn(2,2,is)+
!    2                          dxn(3,2,is)*vnn(3,2,is)
!           print *,'dxn3*cl',dxn(1,3,is)*cnn(2,il)-
!    1                        dxn(2,3,is)*cnn(1,il)
!           print *,'dxn2*cl',dxn(1,2,is)*cnn(2,il)-
!    1                        dxn(2,2,is)*cnn(1,il)
            if (xn(3,3,is).lt.z3min) then
               is3min=is
               z3min=xn(3,3,is)
            endif
            if (xn(3,2,is).lt.z2min) then
               is2min=is
               z2min=xn(3,2,is)
            endif
         endif 
      enddo
!     print *,'is3min,z,dz',is3min,xn(3,3,is3min),dxn(3,3,is3min)
!     print *,'is3min,z0,dzv',is3min,hvec(3,3,is3min),
!    1                        hvec(9,3,is3min)*tstp
!     print *,'is2min,z,dz',is2min,xn(3,2,is2min),dxn(3,2,is2min)
!
! for debugging might want a dump
!     if (mod(kount,50).eq.0) then
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=dxn(1,lv,is)
!           hvec(8,lv,is)=dxn(2,lv,is)
!           hvec(9,lv,is)=dxn(3,lv,is)
!        enddo
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!     stop
!     hvec(7:9,:,1:ns)=0d0
!     endif
!     kount=kount+1
!
!--update the node positions while limiting the displacements
!
      dratio=0d0
      dmean=0d0
      dmax=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            do lv=1,3
              dist=sqrt(dxn(1,lv,is)**2+dxn(2,lv,is)**2+dxn(3,lv,is)**2)
              dmax=max(dmax,dist)
              dmean=dist+dmean
               dscale=abs(voln(lv,is))**(1./3.)
               dratio=max(dratio,dist/dscale)
               if (dist.gt.deps*dscale) then
                  fac=deps*dscale/dist
                  dxn(1,lv,is)=dxn(1,lv,is)*fac
                  dxn(2,lv,is)=dxn(2,lv,is)*fac
                  dxn(3,lv,is)=dxn(3,lv,is)*fac
               endif
            enddo
         else !an edge stack
            dscale=min(abs(voln(3,is))**(1./3.),xn(3,2,is),xn(3,3,is))
            dist=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dmax=max(dmax,dist)
            dmean=dist+dmean
            dratio=max(dratio,dist/dscale)
            if (dist.gt.(deps*dscale)) then
               fac=deps*dscale/dist
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
               dxn(1,2,is)=dxn(1,2,is)*fac
               dxn(2,2,is)=dxn(2,2,is)*fac
               dxn(3,2,is)=dxn(3,2,is)*fac
            endif
         endif
      enddo
      dmean=dmean/(3*ns)
!     print *,'dmax,dmean',dmax,dmean
!
      do is=1,ns
         do lv=1,3
            xn(1,lv,is)=xn(1,lv,is)+dxn(1,lv,is)
            xn(2,lv,is)=xn(2,lv,is)+dxn(2,lv,is)
            xn(3,lv,is)=xn(3,lv,is)+dxn(3,lv,is)
         enddo
      enddo
!
      return
      end subroutine rezone
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Winslow(dx1,dx2,dy1,dy2,dz1,dz2,
     1                   g11,g22,F)
!
      implicit none
      real(8) dx1,dx2,dy1,dy2,dz1,dz2
      real(8) g11,g22,F
!
      g11=dx1*dx1+dy1*dy1+dz1*dz1
      g22=dx2*dx2+dy2*dy2+dz2*dz2
      F=g11+g22
      return
      end subroutine Winslow
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dWin(dx1,dx2,dy1,dy2,dz1,dz2,g11,g22,F,Jcb,dJ,dS1,dS2,
     3                Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
!
      implicit none
      real(8) dx1,dx2,dy1,dy2,dz1,dz2
      real(8) g11,g22,F,Jcb,dJ,dS1,dS2
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
!
      real(8) g11x,g11y,g11z,g22x,g22y,g22z
      real(8) g11xx,g11yy,g11zz,g22xx,g22yy,g22zz
!
      real(8) Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Nx,Ny,Nz,Nx_y,Nx_z,Ny_x,Ny_z,Nz_x,Nz_y
      real(8) dJ2,dJ3,Jcbx,Jcby,Jcbz,Jcbxx,Jcbyy,Jcbzz,Jcbxy,Jcbxz,Jcbyz
!
!
!               1st derivatives of metric
!--diagonal terms
      g11x=2d0*dS1*dx1; g11y=2d0*dS1*dy1; g11z=2d0*dS1*dz1
      g22x=2d0*dS2*dx2; g22y=2d0*dS2*dy2; g22z=2d0*dS2*dz2
!
!              2nd derivatives of metric
!--diagnonal terms
      g11xx=2d0*dS1*dS1; g11yy=g11xx; g11zz=g11xx
      g22xx=2d0*dS2*dS2; g22yy=g22xx; g22zz=g22xx
!--note: xy cross terms are zero
!
!             gradient of F
!
      Fx=g11x+g22x
      Fy=g11y+g22y
      Fz=g11z+g22z
!
!             Hessian of F
!
      Fxx=g11xx+g22xx
      Fyy=g11yy+g22yy
      Fzz=g11zz+g22zz
!--note: xy cross terms are zero
!
!           gradient of Jacobian
!
      Nx=dy1*dz2-dy2*dz1
      Ny=dz1*dx2-dz2*dx1
      Nz=dx1*dy2-dx2*dy1
!
      Nx_y=dS1*dz2-dS2*dz1
      Nx_z=dy1*dS2-dy2*dS1
      Ny_x=dz1*dS2-dz2*dS1
      Ny_z=dS1*dx2-dS2*dx1
      Nz_x=dS1*dy2-dS2*dy1
      Nz_y=dx1*dS2-dx2*dS1
!
! Nx_yy=Nx_yz=0 etc.
!
      Jcbx=dJ*(Ny*Ny_x+Nz*Nz_x)
      Jcby=dJ*(Nx*Nx_y+Nz*Nz_y)
      Jcbz=dJ*(Nx*Nx_z+Ny*Ny_z)
!
!     Hessian of Jacobian
!
      Jcbxx=dJ*(-Jcbx*Jcbx+Ny_x*Ny_x+Nz_x*Nz_x)
      Jcbyy=dJ*(-Jcby*Jcby+Nx_y*Nx_y+Nz_y*Nz_y)
      Jcbzz=dJ*(-Jcbz*Jcbz+Nx_z*Nx_z+Ny_z*Ny_z)
      Jcbxy=dJ*(-Jcbx*Jcby+Nz_x*Nz_y)
      Jcbxz=dJ*(-Jcbx*Jcbz+Ny_x*Ny_z)
      Jcbyz=dJ*(-Jcby*Jcbz+Nx_y*Nx_z)
!
      dJ2=dJ*dJ
      Wxi=(Fx*Jcb-F*Jcbx)*dJ2
      Wyi=(Fy*Jcb-F*Jcby)*dJ2
      Wzi=(Fz*Jcb-F*Jcbz)*dJ2
      dJ3=dJ*dJ2
      Wxxi=(Fxx*Jcb*Jcb-2d0*Jcbx*(Fx*Jcb-F*Jcbx)-F*Jcb*Jcbxx)*dJ3
      Wyyi=(Fyy*Jcb*Jcb-2d0*Jcby*(Fy*Jcb-F*Jcby)-F*Jcb*Jcbyy)*dJ3
      Wzzi=(Fzz*Jcb*Jcb-2d0*Jcbz*(Fz*Jcb-F*Jcbz)-F*Jcb*Jcbzz)*dJ3
      Wxyi=(-(Fx*Jcby+Fy*Jcbx)*Jcb+2d0*F*Jcbx*Jcby-F*Jcbxy*Jcb)*dJ3
      Wxzi=(-(Fx*Jcbz+Fz*Jcbx)*Jcb+2d0*F*Jcbx*Jcbz-F*Jcbxz*Jcb)*dJ3
      Wyzi=(-(Fy*Jcbz+Fz*Jcby)*Jcb+2d0*F*Jcby*Jcbz-F*Jcbyz*Jcb)*dJ3
!
      return
      end subroutine dWin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine TTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1               g11,g22,g33,g12,g13,g23,
     2               F,Jcb,dJ2,dJ3,Wg)
!
      implicit none
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23
      real(8) F,Jcb,dJ2,dJ3,Wg
!
      real(8) a(3,3)
!
      g11=dx1*dx1+dy1*dy1+dz1*dz1
      g22=dx2*dx2+dy2*dy2+dz2*dz2
      g33=dx3*dx3+dy3*dy3+dz3*dz3
      g12=dx1*dx2+dy1*dy2+dz1*dz2
      g13=dx1*dx3+dy1*dy3+dz1*dz3
      g23=dx2*dx3+dy2*dy3+dz2*dz3
      F=g11*g22+g11*g33+g22*g33-g12*g12-g13*g13-g23*g23
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcb=det3x3(a)
      dJ2=1d0/(Jcb*Jcb)
      dJ3=1d0/(Jcb*Jcb*Jcb)
      Wg=F/Jcb
      return
      end subroutine TTM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dTTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,
     2                dH1,dH2,dH3,
     3                Wxi,Wyi,Wzi,
     4                Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
!
      implicit none
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3
      real(8) dH1,dH2,dH3
      real(8) Wxi,Wyi,Wzi
      real(8) Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
!
      real(8) Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcbx,Jcby,Jcbz
!
      real(8) g11x,g11y,g11z,g22x,g22y,g22z,g33x,g33y,g33z
      real(8) g12x,g12y,g12z,g13x,g13y,g13z,g23x,g23y,g23z
      real(8) g11xx,g11yy,g11zz,g22xx,g22yy,g22zz,g33xx,g33yy,g33zz
      real(8) g12xx,g12yy,g12zz,g13xx,g13yy,g13zz,g23xx,g23yy,g23zz
      real(8) a(3,3)
!
!
!               1st derivatives of metric
!--diagonal terms
      g11x=2d0*dH1*dx1; g11y=2d0*dH1*dy1; g11z=2d0*dH1*dz1
      g22x=2d0*dH2*dx2; g22y=2d0*dH2*dy2; g22z=2d0*dH2*dz2
      g33x=2d0*dH3*dx3; g33y=2d0*dH3*dy3; g33z=2d0*dH3*dz3
!
!--off diagonal terms
      g12x=dH1*dx2+dH2*dx1; g12y=dH1*dy2+dH2*dy1; g12z=dH1*dz2+dH2*dz1
      g13x=dH1*dx3+dH3*dx1; g13y=dH1*dy3+dH3*dy1; g13z=dH1*dz3+dH3*dz1
      g23x=dH2*dx3+dH3*dx2; g23y=dH2*dy3+dH3*dy2; g23z=dH2*dz3+dH3*dz2
!
!              2nd derivatives of metric
!--diagnonal terms
      g11xx=2d0*dH1*dH1; g11yy=g11xx; g11zz=g11xx
      g22xx=2d0*dH2*dH2; g22yy=g22xx; g22zz=g22xx
      g33xx=2d0*dH3*dH3; g33yy=g33xx; g33zz=g33xx
!
!--off diagonal terms
      g12xx=2d0*dH1*dH2; g12yy=g12xx; g12zz=g12xx
      g13xx=2d0*dH1*dH3; g13yy=g13xx; g13zz=g13xx
      g23xx=2d0*dH2*dH3; g23yy=g23xx; g23zz=g23xx
!--note: cross terms are zero
!
!             gradient of F
!
      Fx=g11x*g22+g11*g22x+g11x*g33+g11*g33x+g22x*g33+g22*g33x
     1                 -2d0*g12x*g12-2d0*g13x*g13-2d0*g23x*g23
      Fy=g11y*g22+g11*g22y+g11y*g33+g11*g33y+g22y*g33+g22*g33y
     1                 -2d0*g12y*g12-2d0*g13y*g13-2d0*g23y*g23
      Fz=g11z*g22+g11*g22z+g11z*g33+g11*g33z+g22z*g33+g22*g33z
     1                 -2d0*g12z*g12-2d0*g13z*g13-2d0*g23z*g23
!
!             Hessian of F
!--diagonal terms
      Fxx=g11xx*g22+2d0*g11x*g22x+g11*g22xx+
     1    g11xx*g33+2d0*g11x*g33x+g11*g33xx+   
     2    g22xx*g33+2d0*g22x*g33x+g22*g33xx
     3    -2d0*g12x*g12x-2d0*g12xx*g12
     4    -2d0*g13x*g13x-2d0*g13xx*g13
     5    -2d0*g23x*g23x-2d0*g23xx*g23
      Fyy=g11yy*g22+2d0*g11y*g22y+g11*g22yy+
     1    g11yy*g33+2d0*g11y*g33y+g11*g33yy+   
     2    g22yy*g33+2d0*g22y*g33y+g22*g33yy
     3    -2d0*g12y*g12y-2d0*g12yy*g12
     4    -2d0*g13y*g13y-2d0*g13yy*g13
     5    -2d0*g23y*g23y-2d0*g23yy*g23
      Fzz=g11zz*g22+2d0*g11z*g22z+g11*g22zz+
     1    g11zz*g33+2d0*g11z*g33z+g11*g33zz+   
     2    g22zz*g33+2d0*g22z*g33z+g22*g33zz
     3    -2d0*g12z*g12z-2d0*g12zz*g12
     4    -2d0*g13z*g13z-2d0*g13zz*g13
     5    -2d0*g23z*g23z-2d0*g23zz*g23
!
!--off diagonal terms
      Fxy=g11x*g22y+g11y*g22x+g11x*g33y+g11y*g33x+g22x*g33y+g22y*g33x
     1    -2d0*g12x*g12y-2d0*g13x*g13y-2d0*g23x*g23y 
      Fxz=g11x*g22z+g11z*g22x+g11x*g33z+g11z*g33x+g22x*g33z+g22z*g33x
     1    -2d0*g12x*g12z-2d0*g13x*g13z-2d0*g23x*g23z 
      Fyz=g11y*g22z+g11z*g22y+g11y*g33z+g11z*g33y+g22y*g33z+g22z*g33y
     1    -2d0*g12y*g12z-2d0*g13y*g13z-2d0*g23y*g23z 
!
!           gradient of Jacobian
!
!--Jcbx
      a(1,1)=dH1; a(1,2)=dH2; a(1,3)=dH3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcbx=det3x3(a)
!--Jcby
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dH1; a(2,2)=dH2; a(2,3)=dH3
      a(3,1)=dz1; a(3,2)=dz2; a(3,3)=dz3
      Jcby=det3x3(a)
!--Jcbz
      a(1,1)=dx1; a(1,2)=dx2; a(1,3)=dx3
      a(2,1)=dy1; a(2,2)=dy2; a(2,3)=dy3
      a(3,1)=dH1; a(3,2)=dH2; a(3,3)=dH3
      Jcbz=det3x3(a)
!
!     Hessian of Jacobian=0
!
      Wxi=(Fx*Jcb-F*Jcbx)*dJ2
      Wyi=(Fy*Jcb-F*Jcby)*dJ2
      Wzi=(Fz*Jcb-F*Jcbz)*dJ2
      Wxxi=(Fxx*Jcb*Jcb-2d0*Jcbx*(Fx*Jcb-F*Jcbx))*dJ3
      Wyyi=(Fyy*Jcb*Jcb-2d0*Jcby*(Fy*Jcb-F*Jcby))*dJ3
      Wzzi=(Fzz*Jcb*Jcb-2d0*Jcbz*(Fz*Jcb-F*Jcbz))*dJ3
      Wxyi=(Fxy*Jcb*Jcb-(Fx*Jcby+Fy*Jcbx)*Jcb+F*2d0*Jcbx*Jcby)*dJ3
      Wxzi=(Fxz*Jcb*Jcb-(Fx*Jcbz+Fz*Jcbx)*Jcb+F*2d0*Jcbx*Jcbz)*dJ3
      Wyzi=(Fyz*Jcb*Jcb-(Fy*Jcbz+Fz*Jcby)*Jcb+F*2d0*Jcby*Jcbz)*dJ3
!
      return
      end subroutine dTTM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      function det3x3(a)
      implicit none
      real(8) det3x3
      real(8) a(3,3)
      det3x3=+a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     1       -a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))
     2       +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
      end function det3x3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mat3solv(Bmat,rhs,Bvec)
!
! this subroutine solves a 3x3 linear system with
! coefficient Bmat and right hand side rhs and places
! the result in Bvec
!
      implicit none
!
      real(8) Bmat(3,3),rhs(3),Bvec(3)
!
      real(8) cofac(3,3), Inv(3,3)
      real(8) ddet
      integer i, j
!
            
!--compute cofactor matrix
      cofac(1,1)=+Bmat(2,2)*Bmat(3,3)-Bmat(3,2)*Bmat(2,3)
      cofac(2,1)=-Bmat(1,2)*Bmat(3,3)+Bmat(3,2)*Bmat(1,3)
      cofac(3,1)=+Bmat(1,2)*Bmat(2,3)-Bmat(2,2)*Bmat(1,3)
      cofac(1,2)=-Bmat(2,1)*Bmat(3,3)+Bmat(3,1)*Bmat(2,3)
      cofac(2,2)=+Bmat(1,1)*Bmat(3,3)-Bmat(3,1)*Bmat(1,3)
      cofac(3,2)=-Bmat(1,1)*Bmat(2,3)+Bmat(2,1)*Bmat(1,3)
      cofac(1,3)=+Bmat(2,1)*Bmat(3,2)-Bmat(3,1)*Bmat(2,2)
      cofac(2,3)=-Bmat(1,1)*Bmat(3,2)+Bmat(3,1)*Bmat(1,2)
      cofac(3,3)=+Bmat(1,1)*Bmat(2,2)-Bmat(2,1)*Bmat(1,2)
!--determinant
      ddet=1d0/(Bmat(1,1)*cofac(1,1)+Bmat(2,1)*cofac(2,1)+
     1                               Bmat(3,1)*cofac(3,1))
      if (abs(ddet).gt.vbig) then
         print *,'mat3solv: determinant nearly 0!',ddet
      endif
!
      do j=1,3
         do i=1,3
            Inv(i,j)=ddet*cofac(j,i)
         enddo
      enddo
!
      Bvec=0d0
      do j=1,3
         do i=1,3
            Bvec(j)=Bvec(j)+Inv(i,j)*rhs(i)
         enddo
      enddo
! 
      return
      end subroutine mat3solv         
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine avtstep(avtstp) bind(C)
!
! this subroutine computes a recommended advection time step:
! for node x,y,z and neighboring node xn,yn,xn at distance dn
! and joined by unit vector un
!  dt=0.05*dn/abs(v.un)
!
      implicit none
!
      real(c_double)::avtstp
!
      integer isp(4), ism(4)
      data isp/2,3,4,1/, ism/4,1,2,3/
      integer iq,is,istack,istackn,istackp
      real(8) x,y,z,vx,vy,vz,xn,yn,zn,dn,unx,uny,unz
      real(8) z2,vz2,z3,vz3
!
      avtstp=vbig
      !print *,'1',avtstp
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            istackn=isoq(ism(is),iq)
            istackp=isoq(isp(is),iq)
!
!--ventral node
            x=hvec(1,1,istack)   
            y=hvec(2,1,istack)   
            z=hvec(3,1,istack)   
            vx=hvec(7,1,istack)   
            vy=hvec(8,1,istack)   
            vz=hvec(9,1,istack)   
!--the previous ventral node in the element
            xn=hvec(1,1,istackn)   
            yn=hvec(2,1,istackn)   
            zn=hvec(3,1,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            !print *,'dn',dn
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            !print *,'check xyz',x,y,z
            !print *,'check v',vx,vy,vz
            !print *,'check un',unx,uny,unz
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'check abs',vtiny+abs(vx*unx+vy*uny+vz*unz)
            !print *,'check this',dn/(vtiny+abs(vx*unx+vy*uny+vz*unz))
            !print *,'2',avtstp
!--the next ventral node in the element
            xn=hvec(1,1,istackp)   
            yn=hvec(2,1,istackp)   
            zn=hvec(3,1,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'3',avtstp
!--the node above (level 2) in the element
            xn=hvec(1,2,istack)   
            yn=hvec(2,2,istack)   
            zn=hvec(3,2,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'4',avtstp
!
!--middle node
            x=hvec(1,2,istack)   
            y=hvec(2,2,istack)   
            z=hvec(3,2,istack)   
            vx=hvec(7,2,istack)   
            vy=hvec(8,2,istack)   
            vz=hvec(9,2,istack)   
!--the previous middle node in the element
            xn=hvec(1,2,istackn)   
            yn=hvec(2,2,istackn)   
            zn=hvec(3,2,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'5',avtstp
!--the next middle node in the element
            xn=hvec(1,2,istackp)   
            yn=hvec(2,2,istackp)   
            zn=hvec(3,2,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'6',avtstp
!--the node below (level 1) in the element
            xn=hvec(1,1,istack)   
            yn=hvec(2,1,istack)   
            zn=hvec(3,1,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'7',avtstp
!--the node above (level 3) in the element
            xn=hvec(1,3,istack)   
            yn=hvec(2,3,istack)   
            zn=hvec(3,3,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'8',avtstp
!
!--dorsal node
            x=hvec(1,3,istack)   
            y=hvec(2,3,istack)   
            z=hvec(3,3,istack)   
            vx=hvec(7,3,istack)   
            vy=hvec(8,3,istack)   
            vz=hvec(9,3,istack)   
!--the previous dorsal node in the element
            xn=hvec(1,3,istackn)   
            yn=hvec(2,3,istackn)   
            zn=hvec(3,3,istackn)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'9',avtstp
!--the next dorsal node in the element
            xn=hvec(1,3,istackp)   
            yn=hvec(2,3,istackp)   
            zn=hvec(3,3,istackp)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'10',avtstp
!--the node below (level 2) in the element
            xn=hvec(1,2,istack)   
            yn=hvec(2,2,istack)   
            zn=hvec(3,2,istack)   
            dn=sqrt((xn-x)**2+(yn-y)**2+(zn-z)**2)
            unx=(xn-x)/dn
            uny=(yn-y)/dn
            unx=(yn-y)/dn
            avtstp=min(avtstp,dn/(vtiny+abs(vx*unx+vy*uny+vz*unz)))
            !print *,'11',avtstp
         enddo
      enddo
      avtstp=0.05d0*avtstp
      !print *,'12',avtstp
!
!--make sure that the middle and top nodes will not descend toward
!  the substratum too fast.
      do istack=1,ns
         z2=hvec(3,2,istack)
         vz2=hvec(9,2,istack)
         if (z2.gt.0d0.and.vz2.lt.0d0) then
            avtstp=min(avtstp,-0.05d0*z2/vz2)
            !print *,'13',avtstp
         endif
         z3=hvec(3,3,istack)
         vz3=hvec(9,3,istack)
         if (z3.gt.0d0.and.vz3.lt.0d0) then
            avtstp=min(avtstp,-0.05d0*z3/vz3)
            !print *,'14',avtstp
         endif
      enddo
!
      return
      end subroutine avtstep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine goquick(nxyz)
!
! This subroutine is an amalgam of goisomap and govnn
! that computes new derivatives of extrinsic wrt intrinsic coord.
! and a new volume normal at the boundaries.
! It does NOT compute new Jacobians, volumes, etc.
!
      implicit none
!
      real(8) nxyz(3,3,NSM)
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8) fac1x,fac2x,fac3x,fac1y,fac2y,fac3y,fac1z,fac2z,fac3z
      real(8) dVx,dVy,dvZ,vnorm,vepsilon
      real(8) tgx,tgy,tgz,dtgnorm
      integer iq,isn,lvn,istack,isg,lvg,j,il
!
      dxV=0d0
      dyV=0d0
      dzV=0d0
      vnn=0d0
!--loop through the elements
      do iq=1,nq
!--loop over Gauss points in cell
         do isg=1,4
            do lvg=1,3
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     do j=1,3
                        dxV(j,lvg,isg,iq)=dxV(j,lvg,isg,iq)+
     1             nxyz(1,lvn,istack)*dH(j,lvn,isn,lvg,isg)
                        dyV(j,lvg,isg,iq)=dyV(j,lvg,isg,iq)+
     1             nxyz(2,lvn,istack)*dH(j,lvn,isn,lvg,isg)
                        dzV(j,lvg,isg,iq)=dzV(j,lvg,isg,iq)+
     1             nxyz(3,lvn,istack)*dH(j,lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
               fac1x=(+dyV(2,lvg,isg,iq)*dzV(3,lvg,isg,iq)
     1                -dyV(3,lvg,isg,iq)*dzV(2,lvg,isg,iq))*wgp(lvg)
               fac2x=(-dyV(1,lvg,isg,iq)*dzV(3,lvg,isg,iq)
     1                +dyV(3,lvg,isg,iq)*dzV(1,lvg,isg,iq))*wgp(lvg)
               fac3x=(+dyV(1,lvg,isg,iq)*dzV(2,lvg,isg,iq)
     1                -dyV(2,lvg,isg,iq)*dzV(1,lvg,isg,iq))*wgp(lvg)
               fac1y=(-dxV(2,lvg,isg,iq)*dzV(3,lvg,isg,iq)
     1                +dxV(3,lvg,isg,iq)*dzV(2,lvg,isg,iq))*wgp(lvg)
               fac2y=(+dxV(1,lvg,isg,iq)*dzV(3,lvg,isg,iq)
     1                -dxV(3,lvg,isg,iq)*dzV(1,lvg,isg,iq))*wgp(lvg)
               fac3y=(-dxV(1,lvg,isg,iq)*dzV(2,lvg,isg,iq)
     1                +dxV(2,lvg,isg,iq)*dzV(1,lvg,isg,iq))*wgp(lvg)
               fac1z=(+dxV(2,lvg,isg,iq)*dyV(3,lvg,isg,iq)
     1                -dxV(3,lvg,isg,iq)*dyV(2,lvg,isg,iq))*wgp(lvg)
               fac2z=(-dxV(1,lvg,isg,iq)*dyV(3,lvg,isg,iq)
     1                +dxV(3,lvg,isg,iq)*dyV(1,lvg,isg,iq))*wgp(lvg)
               fac3z=(+dxV(1,lvg,isg,iq)*dyV(2,lvg,isg,iq)
     1                -dxV(2,lvg,isg,iq)*dyV(1,lvg,isg,iq))*wgp(lvg)
c--contribution of this GP to 12 nodes
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3 
                     dVx=dH(1,lvn,isn,lvg,isg)*fac1x+
     1                   dH(2,lvn,isn,lvg,isg)*fac2x+
     2                   dH(3,lvn,isn,lvg,isg)*fac3x
                     dVy=dH(1,lvn,isn,lvg,isg)*fac1y+
     1                   dH(2,lvn,isn,lvg,isg)*fac2y+
     2                   dH(3,lvn,isn,lvg,isg)*fac3y
                     dVz=dH(1,lvn,isn,lvg,isg)*fac1z+
     1                   dH(2,lvn,isn,lvg,isg)*fac2z+
     2                   dH(3,lvn,isn,lvg,isg)*fac3z
                     vnn(1,lvn,istack)=vnn(1,lvn,istack)+dVx
                     vnn(2,lvn,istack)=vnn(2,lvn,istack)+dVy
                     vnn(3,lvn,istack)=vnn(3,lvn,istack)+dVz
                  enddo
               enddo
            enddo
         enddo
      enddo
      
!
!--normalize the vnns
!
      do isn=1,ns
         do lvn=1,3
            vnorm=sqrt(vnn(1,lvn,isn)**2+vnn(2,lvn,isn)**2
     1                                  +vnn(3,lvn,isn)**2)
            vnn(0,lvn,isn)=vnorm
            vepsilon=1d-9*abs(voln(lvn,isn))**(2./3.)
            vnn(1,lvn,isn)=vnn(1,lvn,isn)/(vnorm+vepsilon)
            vnn(2,lvn,isn)=vnn(2,lvn,isn)/(vnorm+vepsilon)
            vnn(3,lvn,isn)=vnn(3,lvn,isn)/(vnorm+vepsilon)
         enddo
      enddo
!
!--compute some useful unit vectors at boundary nodes
!
      do isn=1,ns
         if (ilos(1,isn).eq.0) then!not an edge stack
            do lvn=1,3,2 !dorsal and surface nodes
! pick which one of the x and y axis most perpendicular to
               if (abs(vnn(1,lvn,isn)).lt.abs(vnn(2,lvn,isn))) then
! utg1 = vnn x e_x
                  tgy=vnn(3,lvn,isn)
                  tgz=-vnn(2,lvn,isn)
                  dtgnorm=1d0/(sqrt(tgy*tgy+tgz*tgz))
                  utg1(1,lvn,isn)=0d0
                  utg1(2,lvn,isn)=tgy*dtgnorm
                  utg1(3,lvn,isn)=tgz*dtgnorm
               else
! utg1 = vnn x e_y
                  tgx=-vnn(3,lvn,isn)
                  tgz=vnn(1,lvn,isn)
                  dtgnorm=1d0/(sqrt(tgx*tgx+tgz*tgz))
                  utg1(1,lvn,isn)=tgx*dtgnorm
                  utg1(2,lvn,isn)=0d0
                  utg1(3,lvn,isn)=tgz*dtgnorm
               endif
! utg2 = vnn x utg1
               utg2(1,lvn,isn)=+vnn(2,lvn,isn)*utg1(3,lvn,isn)
     1                         -vnn(3,lvn,isn)*utg1(2,lvn,isn)
               utg2(2,lvn,isn)=-vnn(1,lvn,isn)*utg1(3,lvn,isn)
     1                         +vnn(3,lvn,isn)*utg1(1,lvn,isn)
               utg2(3,lvn,isn)=+vnn(1,lvn,isn)*utg1(2,lvn,isn)
     1                         -vnn(2,lvn,isn)*utg1(1,lvn,isn)
            enddo
         else !it is an edge stack
            do lvn=1,3 !all three nodes on boundary
               tgx=vnn(2,lvn,isn)
               tgy=-vnn(1,lvn,isn)
               dtgnorm=1d0/(sqrt(tgx*tgx+tgy*tgy))
! utg1 = vnn x e_z
               utg1(1,lvn,isn)=tgx*dtgnorm
               utg1(2,lvn,isn)=tgy*dtgnorm
               utg1(3,lvn,isn)=0d0
            enddo
         endif
      enddo
!
      return
      end subroutine goquick
!       
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            OLD ROUTINES
!
!
      subroutine rezbot2(xn,c_newton,dratio,wtot)
!
! This subroutine rezones the interior ventral nodes
! so as to minimize the 2D Winslow functional
!
      implicit none
!
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
      real(8) dratio
!
      real(8) deps !maximum node motion parameter
      parameter(deps=5d-1)
      real(8) dxb(2,NSM)
!
      real(8) dx1,dx2,dy1,dy2
      real(8) g11,g22 !metric coeff.
      real(8) dS_1,dS_2
      real(8) F,Fx,Fy,Fxx,Fyy,Fxy
      real(8) Jcb, Jcbx,Jcby,dJ2, dJ3
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(NSM),Wy(NSM) !grad of Winslow functional
      real(8) Wxx(NSM),Wyy(NSM),Wxy(NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wxxi,Wyyi,Wxyi
      real(8) HessWtg11,HessWtg12,HessWtg22 
      real(8) gradWtg1, gradWtg2, ddet, dx, dy
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
      dxb=0d0
!
!--Initialize the W's  
      Wtot=0d0
      Win=0d0
      Wx=0d0
      Wy=0d0
      Wxx=0d0
      Wyy=0d0
      Wxy=0d0
!
      do iq=1,nq !outer loop over iq quadrilateral
         do isg=1,4 !loop over the gauss points of iq
            dx1=dxSv(1,isg,iq)
            dx2=dxSv(2,isg,iq)
            dy1=dySv(1,isg,iq)
            dy2=dySv(2,isg,iq)
            call TTMbot(dx1,dx2,dy1,dy2,
     1                  g11,g22,F,Jcb,dJ2,dJ3,Wg)
            Win(iq)=Win(iq)+Wg
            Wtot=Wtot+Wg
            do isn=1,4
               istack=isoq(isn,iq)
               dS_1=dS4(1,isn,isg)
               dS_2=dS4(2,isn,isg)
               call dTTMbot(dx1,dx2,dy1,dy2,
     1                      g11,g22,F,Jcb,dJ2,dJ3,
     2                      dS_1,dS_2,
     3                      Wxi,Wyi,Wxxi,Wyyi,Wxyi)
               Wx(istack)=Wx(istack)+Wxi
               Wy(istack)=Wy(istack)+Wyi
               Wxx(istack)=Wxx(istack)+Wxxi
               Wyy(istack)=Wyy(istack)+Wyyi
               Wxy(istack)=Wxy(istack)+Wxyi
            enddo
         enddo
      enddo
      print *,'Wtot',Wtot
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation) 
!
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
            gradWtg1=-Wx(is)
            gradWtg2=-Wy(is)
            HessWtg11=Wxx(is)
            HessWtg22=Wyy(is)
            HessWtg12=Wxy(is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dx=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dy=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxb(1,is)=c_newton*dx
            dxb(2,is)=c_newton*dy
         endif
      enddo
! for debugging might want a dump
!     if (mod(kount,50).eq.0) then
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=0d0
!           hvec(8,lv,is)=0d0
!           hvec(9,lv,is)=0d0
!        enddo
!        hvec(7,1,is)=dxb(1,is)
!        hvec(8,1,is)=dxb(2,is)
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!     hvec(7:9,:,1:ns)=0d0
!     endif
!     kount=kount+1
!
!--update the node positions while limiting the displacements
!
      dratio=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            dist=sqrt(dxb(1,is)**2+dxb(2,is)**2)
            dscale=abs(voln(lv,is))**(1./3.)
            dratio=max(dratio,dist/dscale)
            if (dist.gt.deps*dscale) then
               fac=deps*dscale/dist
               dxb(1,is)=dxb(1,is)*fac
               dxb(2,is)=dxb(2,is)*fac
            endif
         endif
      enddo
!
      do is=1,ns
         xn(1,1,is)=xn(1,1,is)+dxb(1,is)
         xn(2,1,is)=xn(2,1,is)+dxb(2,is)
      enddo
!
      return
      end subroutine rezbot2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine TTMbot(dx1,dx2,dy1,dy2,
     1                  g11,g22,F,Jcb,dJ2,dJ3,Wg)
!
      implicit none
      real(8) dx1,dx2,dy1,dy2
      real(8) g11,g22
      real(8) F,Jcb,dJ2,dJ3,Wg
!
      g11=dx1*dx1+dy1*dy1
      g22=dx2*dx2+dy2*dy2
      F=g11+g22
      Jcb=dx1*dy2-dx2*dy1
      dJ2=1d0/(Jcb*Jcb)
      dJ3=1d0/(Jcb*Jcb*Jcb)
      Wg=F/Jcb
!
      return
      end subroutine TTMbot
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine dTTMbot(dx1,dx2,dy1,dy2,
     1                   g11,g22,F,Jcb,dJ2,dJ3,
     2                   dS_1,dS_2,
     3                   Wxi,Wyi,Wxxi,Wyyi,Wxyi)
!
      implicit none
      real(8) dx1,dx2,dy1,dy2
      real(8) g11,g22,F,Jcb,dJ2,dJ3
      real(8) dS_1,dS_2
      real(8) Wxi,Wyi,Wxxi,Wyyi,Wxyi
!
      real(8) Fx,Fy,Fz,Fxx,Fyy,Fxy
      real(8) Jcbx,Jcby
!
      real(8) g11x,g11y,g22x,g22y
      real(8) g12x,g12y,g12z,g13x,g13y,g13z,g23x,g23y,g23z
      real(8) g11xx,g11yy,g22xx,g22yy
      real(8) g12xx,g12yy
!
!
! 1st derivatives of metric
      g11x=2d0*dS_1*dx1; g11y=2d0*dS_1*dy1
      g22x=2d0*dS_2*dx2; g22y=2d0*dS_2*dy2
!
! 2nd derivatives of metric
      g11xx=2d0*dS_1*dS_1; g11yy=g11xx
      g22xx=2d0*dS_2*dS_2; g22yy=g22xx
!--note: cross terms are zero
!
! gradient of F
      Fx=g11x+g22x
      Fy=g11y+g22y
!
! Hessian of F
!--diagonal terms
      Fxx=g11xx+g22xx
      Fyy=g11yy+g22yy
      Fxy=0d0
!
!           gradient of Jacobian
!
      Jcbx=dS_1*dy2-dy1*dS_2
      Jcby=dx1*dS_2-dx2*dS_1
!     Hessian of Jacobian=0
!
      Wxi=(Fx*Jcb-F*Jcbx)*dJ2
      Wyi=(Fy*Jcb-F*Jcby)*dJ2
      Wxxi=(Fxx*Jcb*Jcb-2d0*Jcbx*(Fx*Jcb-F*Jcbx))*dJ3
      Wyyi=(Fyy*Jcb*Jcb-2d0*Jcby*(Fy*Jcb-F*Jcby))*dJ3
      Wxyi=(-(Fx*Jcby+Fy*Jcbx)*Jcb+F*2d0*Jcbx*Jcby)*dJ3
!
      return
      end subroutine dTTMbot

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rezone_old(xn,c_newton,dratio,wtot)
!
! This subroutine rezones the nodes
! so as to minimize the Winslow functional
!
      implicit none
!
      real(8) xn(3,3,NSM) !node coordinates (xyz,level,stack)
      real(8) c_newton
      real(8) dratio
!
      real(8) deps !maximum node motion parameter
      parameter(deps=1d-1)
      real(8) dxn(3,3,NSM)
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) g11,g22,g33,g12,g13,g23 !metric coeff.
      real(8) dH1,dH2,dH3
      real(8) F,Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxy,Fxz,Fyz
      real(8) Jcb, Jcbx,Jcby,Jcbz, dJ2, dJ3
      real(8) Win(NSM),Wtot,Wg
      real(8) Wx(3,NSM),Wy(3,NSM),Wz(3,NSM) !grad of Winslow functional
      real(8) Wxx(3,NSM),Wyy(3,NSM),Wzz(3,NSM) !Hessian of Winslow
      real(8) Wxy(3,NSM),Wxz(3,NSM),Wyz(3,NSM) !Hessian of Winslow
      real(8) Wxi,Wyi,Wzi,Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi
      real(8) HessW(3,3), gradW(3), dvec(3),gradWn
      real(8) HessWtg11,HessWtg12,HessWtg22
      real(8) gradWtg1, gradWtg2, dtg1,dtg2,ddet
      real(8) u1x,u1y,u1z,u2x,u2y,u2z
      real(8) b(4,4),blu(4,4),dis(4),vrhs(4)
      integer ipvt(4), isng
      real(8) b9(9,9),blu9(9,9),dis9(9),vrhs9(9)
      integer ipvt9(9), isng9
      real(8) z2mean,z3mean
      integer,save::kount=0
!
      real(8) tnx,tny,tnz,dn,fac,dist,dscale
      integer iq,isg,lvg,isn,lvn,is,lv,istack,il,lve
!
!--Initialize the W's
      Wtot=0d0
      Win=0d0
      Wx=0d0
      Wy=0d0
      Wz=0d0
      Wxx=0d0
      Wyy=0d0
      Wzz=0d0
      Wxy=0d0
      Wxz=0d0
      Wyz=0d0
!
      do iq=1,nq !outer loop over iq quadrilateral
         do isg=1,4 !loop over the gauss points of iq
            do lvg=1,3
               dx1=dxV(1,lvg,isg,iq)
               dx2=dxV(2,lvg,isg,iq)
               dx3=dxV(3,lvg,isg,iq)
               dy1=dyV(1,lvg,isg,iq)
               dy2=dyV(2,lvg,isg,iq)
               dy3=dyV(3,lvg,isg,iq)
               dz1=dzV(1,lvg,isg,iq)
               dz2=dzV(2,lvg,isg,iq)
               dz3=dzV(3,lvg,isg,iq)
               call TTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                  g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,Wg)
               Win(iq)=Win(iq)+wgp(lvg)*Wg
               Wtot=Wtot+wgp(lvg)*Wg
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     dH1=dH(1,lvn,isn,lvg,isg)
                     dH2=dH(2,lvn,isn,lvg,isg)
                     dH3=dH(3,lvn,isn,lvg,isg)
                     call dTTM(dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3,
     1                         g11,g22,g33,g12,g13,g23,F,Jcb,dJ2,dJ3,
     2                         dH1,dH2,dH3,
     3                         Wxi,Wyi,Wzi,
     4                         Wxxi,Wyyi,Wzzi,Wxyi,Wxzi,Wyzi)
                     Wx(lvn,istack)=Wx(lvn,istack)+Wxi*wgp(lvg)
                     Wy(lvn,istack)=Wy(lvn,istack)+Wyi*wgp(lvg)
                     Wz(lvn,istack)=Wz(lvn,istack)+Wzi*wgp(lvg)
                     Wxx(lvn,istack)=Wxx(lvn,istack)+Wxxi*wgp(lvg)
                     Wyy(lvn,istack)=Wyy(lvn,istack)+Wyyi*wgp(lvg)
                     Wzz(lvn,istack)=Wzz(lvn,istack)+Wzzi*wgp(lvg)
                     Wxy(lvn,istack)=Wxy(lvn,istack)+Wxyi*wgp(lvg)
                     Wxz(lvn,istack)=Wxz(lvn,istack)+Wxzi*wgp(lvg)
                     Wyz(lvn,istack)=Wyz(lvn,istack)+Wyzi*wgp(lvg)
                  enddo
               enddo
            enddo
         enddo
      enddo
      print *,'Wtot',Wtot
!     stop
!
!-- We now compute the node motions needed to minimize the winslow
!   functional (but with 50% under-relaxation)
!
      z2mean=0d0
      z3mean=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then!not an edge stack
!--ventral nodes are restricted to move in the z=0 plane
            gradWtg1=-Wx(1,is)
            gradWtg2=-Wy(1,is)
            HessWtg11=Wxx(1,is)
            HessWtg22=Wyy(1,is)
            HessWtg12=Wxy(1,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,1,is)=c_newton*dtg1
            dxn(2,1,is)=c_newton*dtg2
            dxn(3,1,is)=0d0
!--middle nodes
            gradW(1)=-Wx(2,is)
            gradW(2)=-Wy(2,is)
            gradW(3)=-Wz(2,is)
            HessW(1,1)=Wxx(2,is)
            HessW(2,1)=Wxy(2,is)
            HessW(3,1)=Wxz(2,is)
            HessW(1,2)=Wxy(2,is)
            HessW(2,2)=Wyy(2,is)
            HessW(3,2)=Wyz(2,is)
            HessW(1,3)=Wxz(2,is)
            HessW(2,3)=Wyz(2,is)
            HessW(3,3)=Wzz(2,is)
            call mat3solv(HessW,gradW,dvec)
            dxn(1,2,is)=c_newton*dvec(1)
            dxn(2,2,is)=c_newton*dvec(2)
            dxn(3,2,is)=c_newton*dvec(3)
!--dorsal nodes
            u1x=utg1(1,3,is)
            u1y=utg1(2,3,is)
            u1z=utg1(3,3,is)
            u2x=utg2(1,3,is)
            u2y=utg2(2,3,is)
            u2z=utg2(3,3,is)
!           gradWn=sqrt(Wx(3,is)**2+Wy(3,is)**2+Wz(3,is)**2)
            gradWtg1=-Wx(3,is)*u1x-Wy(3,is)*u1y-Wz(3,is)*u1z
            gradWtg2=-Wx(3,is)*u2x-Wy(3,is)*u2y-Wz(3,is)*u2z
            HessWtg11=+u1x*u1x*Wxx(3,is)+u1y*u1y*Wyy(3,is)
     1                +u1z*u1z*Wzz(3,is)
     2                +2d0*(u1x*u1y*Wxy(3,is)+u1x*u1z*Wxz(3,is)+
     3                    u1y*u1z*Wyz(3,is))
            HessWtg22=+u2x*u2x*Wxx(3,is)+u2y*u2y*Wyy(3,is)
     1                +u2z*u2z*Wzz(3,is)
     2                +2d0*(u2x*u2y*Wxy(3,is)+u2x*u2z*Wxz(3,is)+
     3                    u2y*u2z*Wyz(3,is))
            HessWtg12=+u1x*u2x*Wxx(3,is)+u1y*u2y*Wyy(3,is)
     1                +u1z*u2z*Wzz(3,is)+
     2                (u1x*u2y+u1y*u2x)*Wxy(3,is)+
     3                (u1x*u2z+u1z*u2x)*Wxz(3,is)+
     4                (u1y*u2z+u1z*u2y)*Wyz(3,is)
            ddet=1d0/(HessWtg11*HessWtg22-HessWtg12*HessWtg12)
            dtg1=(gradWtg1*HessWtg22-gradWtg2*HessWtg12)*ddet
            dtg2=(gradWtg2*HessWtg11-gradWtg1*HessWtg12)*ddet
            dxn(1,3,is)=c_newton*(dtg1*u1x+dtg2*u2x)
            dxn(2,3,is)=c_newton*(dtg1*u1y+dtg2*u2y)
            dxn(3,3,is)=c_newton*(dtg1*u1z+dtg2*u2z)
!
!--Lagrange multiplier method
!
!              b=0d0
!              b(1,1)=Wxx(3,is)
!              b(2,1)=Wxy(3,is)
!              b(3,1)=Wxz(3,is)
!              b(1,2)=Wxy(3,is)
!              b(2,2)=Wyy(3,is)
!              b(3,2)=Wyz(3,is)
!              b(1,3)=Wxz(3,is)
!              b(2,3)=Wyz(3,is)
!              b(3,3)=Wzz(3,is)
!              b(4,1)=vnn(1,3,is)
!              b(4,2)=vnn(2,3,is)
!              b(4,3)=vnn(3,3,is)
!              b(1,4)=vnn(1,3,is)
!              b(2,4)=vnn(2,3,is)
!              b(3,4)=vnn(3,3,is)
!              call lufa(b,4,blu,isng,ipvt)
!              if (isng.ne.0) print *,
!    1            'ill conditioned matrix in rezone, isng=',isng
!              vrhs=0d0
!              vrhs(1)=-Wx(3,is)
!              vrhs(2)=-Wy(3,is)
!              vrhs(3)=-Wz(3,is)
!              call lusl(blu,vrhs,isng,ipvt,4,dis)
!              print *,'3,is',3,is
!              print *,'1: dxn,dis',dxn(1,3,is),0.5*dis(1)
!              print *,'2: dxn,dis',dxn(2,3,is),0.5*dis(2)
!              print *,'3: dxn,dis',dxn(3,3,is),0.5*dis(3)
!              print *,'dot products',dxn(1,3,is)*vnn(1,3,is)+
!    1                                dxn(2,3,is)*vnn(2,3,is)+
!    2                                dxn(3,3,is)*vnn(3,3,is),
!    3 0.5*(dis(1)*vnn(1,3,is)+dis(2)*vnn(2,3,is)+dis(3)*vnn(3,3,is))
!              print *,'gradWn-gradHdotn',gradWn-
!    1     abs(Wx(3,is)*vnn(1,3,is)+Wy(3,is)*vnn(2,3,is)+
!    1                                Wz(3,is)*vnn(3,3,is))
         else !it is an edge stack
            lve=3
            do lv=1,lve
               u1x=utg1(1,lv,is)
               u1y=utg1(2,lv,is)
               u1z=utg1(3,lv,is)
               gradWtg1=-Wx(lv,is)*u1x-Wy(lv,is)*u1y-Wz(lv,is)*u1z
               HessWtg11=+u1x*u1x*Wxx(lv,is)+u1y*u1y*Wyy(lv,is)
     1                   +u1z*u1z*Wzz(lv,is)
     2              +2d0*(u1x*u1y*Wxy(lv,is)+u1x*u1z*Wxz(lv,is)+
     3                    u1y*u1z*Wyz(lv,is))
               dtg1=gradWtg1/HessWtg11
               dxn(1,lv,is)=c_newton*dtg1*u1x
               dxn(2,lv,is)=c_newton*dtg1*u1y
               dxn(3,lv,is)=c_newton*dtg1*u1z !this should be zero
!              if (lv.eq.1) print *,'e: dxn(3,1,is),is',dxn(3,1,is),is
            enddo
            if (lve.eq.1) then
!
!--Lagrange multiplier method
!
            b9=0d0
            b9(1,1)=Wxx(2,is); b9(1,2)=Wxy(2,is); b9(1,3)=Wxz(2,is)
            b9(2,1)=Wxy(2,is); b9(2,2)=Wyy(2,is); b9(2,3)=Wyz(2,is)
            b9(3,1)=Wxz(2,is); b9(3,2)=Wyz(2,is); b9(3,3)=Wzz(2,is)
            b9(1,4)=vnn(1,2,is);b9(2,4)=vnn(2,2,is);b9(3,4)=vnn(3,2,is)
            b9(4,1)=vnn(1,2,is);b9(4,2)=vnn(2,2,is);b9(4,3)=vnn(3,2,is)
            b9(5,5)=Wxx(3,is); b9(5,6)=Wxy(3,is); b9(5,7)=Wxz(3,is)
            b9(6,5)=Wxy(3,is); b9(6,6)=Wyy(3,is); b9(6,7)=Wyz(3,is)
            b9(7,5)=Wxz(3,is); b9(7,6)=Wyz(3,is); b9(7,7)=Wzz(3,is)
            b9(5,8)=vnn(1,3,is);b9(6,8)=vnn(2,3,is);b9(7,8)=vnn(3,3,is)
            b9(8,5)=vnn(1,3,is);b9(8,6)=vnn(2,3,is);b9(8,7)=vnn(3,3,is)
            b9(9,3)=xn(3,3,is);b9(9,7)=-xn(3,2,is)
            b9(3,9)=xn(3,3,is);b9(7,9)=-xn(3,2,is)
            call lufa(b9,9,blu9,isng9,ipvt9)
            if (isng9.ne.0) print *,
     1         'ill conditioned matrix in rezone, isng=',isng
            vrhs9=0d0
            vrhs9(1)=-Wx(2,is); vrhs9(2)=-Wy(2,is); vrhs9(3)=-Wz(2,is)
            vrhs9(5)=-Wx(3,is); vrhs9(6)=-Wy(3,is); vrhs9(7)=-Wz(3,is)
            call lusl(blu9,vrhs9,isng9,ipvt9,9,dis9)
            dxn(1,2,is)=c_newton*dis9(1)
            dxn(2,2,is)=c_newton*dis9(2)
            dxn(3,2,is)=c_newton*dis9(3)
            dxn(1,3,is)=c_newton*dis9(5)
            dxn(2,3,is)=c_newton*dis9(6)
            dxn(3,3,is)=c_newton*dis9(7)
            dxn(3,3,is)=c_newton*2d0*dis9(3)
!           print *,'dz3/z3,dz2/d2',
!    1               dxn(3,3,is)/xn(3,3,is),dxn(3,2,is)/xn(3,2,is)
!           print *,'dot products 3 and 2',
!    1      dis9(5)*vnn(1,3,is)+dis9(6)*vnn(2,3,is)+dis9(7)*vnn(3,3,is),
!    2      dis9(1)*vnn(1,2,is)+dis9(2)*vnn(2,2,is)+dis9(3)*vnn(3,2,is)
!           print *,'z2,z3',xn(3,2,is),xn(3,3,is)
!           z2mean=z2mean+xn(3,2,is)
!           z3mean=z3mean+xn(3,3,is)
            endif
         endif
      enddo
!     print *,'z2mean,z3mean',z2mean/nl,z3mean/nl
!
! for debugging might want a dump
!     if (mod(kount,50).eq.0) then
!     do is=1,ns
!        do lv=1,3
!           hvec(1,lv,is)=xn(1,lv,is)
!           hvec(2,lv,is)=xn(2,lv,is)
!           hvec(3,lv,is)=xn(3,lv,is)
!           hvec(4,lv,is)=vnn(1,lv,is)
!           hvec(5,lv,is)=vnn(2,lv,is)
!           hvec(6,lv,is)=vnn(3,lv,is)
!           hvec(7,lv,is)=dxn(1,lv,is)
!           hvec(8,lv,is)=dxn(2,lv,is)
!           hvec(9,lv,is)=dxn(3,lv,is)
!        enddo
!     enddo
!     open(66,file='av.dump')
!     call iowrfile(99,66)
!     close(66)
!     hvec(7:9,:,1:ns)=0d0
!     endif
!     kount=kount+1
!
!--update the node positions while limiting the displacements
!
      dratio=0d0
      do is=1,ns
         if (ilos(1,is).eq.0) then !not an edge stack
            do lv=1,3
               dist=sqrt(dxn(1,lv,is)**2+
     1                   dxn(2,lv,is)**2+dxn(3,lv,is)**2)
               dscale=abs(voln(lv,is))**(1./3.)
               dratio=max(dratio,dist/dscale)
               if (dist.gt.deps*dscale) then
                  fac=deps*dist/dscale
                  dxn(1,lv,is)=dxn(1,lv,is)*fac
                  dxn(2,lv,is)=dxn(2,lv,is)*fac
                  dxn(3,lv,is)=dxn(3,lv,is)*fac
               endif
            enddo
         else !an edge stack, restrict motion more
            dscale=abs(voln(3,is))**(1./3.)
            dist=sqrt(dxn(1,3,is)**2+dxn(2,3,is)**2+dxn(3,3,is)**2)
            dratio=max(dratio,dist/dscale)
            if (dist.gt.(0.2d0*deps*dscale)) then
               fac=0.2d0*deps*dist/dscale
               dxn(1,3,is)=dxn(1,3,is)*fac
               dxn(2,3,is)=dxn(2,3,is)*fac
               dxn(3,3,is)=dxn(3,3,is)*fac
               dxn(1,2,is)=dxn(1,2,is)*fac
               dxn(2,2,is)=dxn(2,2,is)*fac
               dxn(3,2,is)=dxn(3,2,is)*fac
            endif
            dscale=abs(voln(1,is))**(1./3.)
            dist=sqrt(dxn(1,1,is)**2+dxn(2,1,is)**2+dxn(3,1,is)**2)
            dratio=max(dratio,dist/dscale)
            if (dist.gt.(0.2d0*deps*dscale)) then
               fac=0.2d0*deps*dist/dscale
               dxn(1,1,is)=dxn(1,1,is)*fac
               dxn(2,1,is)=dxn(2,1,is)*fac
               dxn(3,1,is)=dxn(3,1,is)*fac
            endif
         endif
      enddo
!
      do is=1,ns
         do lv=1,3
            xn(1,lv,is)=xn(1,lv,is)+dxn(1,lv,is)
            xn(2,lv,is)=xn(2,lv,is)+dxn(2,lv,is)
            xn(3,lv,is)=xn(3,lv,is)+dxn(3,lv,is)
         enddo
      enddo
!
      return
      end subroutine rezone_old
!
      END MODULE AVLIBSW!end library of advection routines.

