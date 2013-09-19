      MODULE mxlibsw
!
      use iolibsw
!
      CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BEGIN DIRECT MATRIX INVERSION (LU DECOMPOSITION) ROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine lussolve(wv,Q,Lv,Ld,Le,vldq,vldv,vldd,vlde)
!
! This will solve: (Q+Lv+Ld+Le)*wv=vldq+vldv+vldd_vlde for wv
!  by LU decomposition. This is VOLUME SCALAR diffusion.
! This routine is the homolog to cgssolve
!
      use lupack
!
      implicit none
!
      real(8) wv(3*NSM)!concentration
      real(8) Q(3,4,3,4,NSM)!element stiffness matrix
      real(8) Lv(4,4,NSM),Ld(4,4,NSM)!ventral/dorsal boundary stiffness
      real(8) Le(3,2,3,2,NLM)!edge boundary stiffness matrix
      real(8) vldq(3,4,NSM)!element load vector
      real(8) vldv(4,NSM),vldd(4,NSM)!ventral/dorsal boundary loads
      real(8) vlde(3,2,NLM)!edge boundary load vector
!
!  Kstif=(3ns x 3ns) stiffness matrix, Kluf in LU decompositon
!  con=(3ns) concentration vector, vload=(3ns) load vector
!  Kinv, Kcheck for debugging only 
!  ipvt(3ns) pivot vector
!
      real(8), allocatable, dimension(:,:) :: Kstif, Kluf
      real(8), allocatable, dimension(:,:) :: Kinv, Kcheck
      real(8), allocatable, dimension(:) :: con, vload, vload2
      integer, allocatable, dimension(:) :: ipvt
      integer mrank, iq, il, isng, idebug, itrouble
      integer js, jstack, jlv, jnode, j
      integer is, istack, ilv, inode, i
      real(8) pKstif, pKinv, vtest
!
      idebug=0
      mrank=3*ns
      allocate (Kstif(mrank,mrank),Kluf(mrank,mrank))
      allocate (con(mrank),vload(mrank))
      allocate (ipvt(mrank))
      if (idebug.eq.1) then
         allocate (Kinv(mrank,mrank),Kcheck(mrank,mrank),vload2(mrank)) 
      endif
!
!--assemble stiffness matrix
!
      Kstif=0d0
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            do jlv=1,3
               jnode=(jstack-1)*3+jlv
               do is=1,4
                  istack=isoq(is,iq)
                  do ilv=1,3
                     inode=(istack-1)*3+ilv
                     Kstif(inode,jnode)=Kstif(inode,jnode)
     1                                  +Q(ilv,is,jlv,js,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
!--ventral contributions (level=1)
      jlv=1
      ilv=1
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            jnode=(jstack-1)*3+jlv
            do is=1,4
               istack=isoq(is,iq)
               inode=(istack-1)*3+ilv
               Kstif(inode,jnode)=Kstif(inode,jnode)+Lv(is,js,iq)
            enddo
         enddo
      enddo
!--dorsal contributions (level=3)
      jlv=3
      ilv=3
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            jnode=(jstack-1)*3+jlv
            do is=1,4
               istack=isoq(is,iq)
               inode=(istack-1)*3+ilv
               Kstif(inode,jnode)=Kstif(inode,jnode)+Ld(is,js,iq)
            enddo
         enddo
      enddo
!--edge contributions
      do il=1,nl
         do js=1,2
            jstack=isol(js,il)
            do jlv=1,3
               jnode=(jstack-1)*3+jlv
               do is=1,2
                  istack=isol(is,il)
                  do ilv=1,3
                     inode=(istack-1)*3+ilv
                     Kstif(inode,jnode)=Kstif(inode,jnode)
     1                                  +Le(ilv,is,jlv,js,il)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--assemble load vector
!
      vload=0d0
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            do ilv=1,3
               inode=(istack-1)*3+ilv
               vload(inode)=vload(inode)+vldq(ilv,is,iq)
            enddo
         enddo
      enddo
!--ventral contributions (level=1)
      ilv=1
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            inode=(istack-1)*3+ilv
            vload(inode)=vload(inode)+vldv(is,iq)
         enddo
      enddo
!--dorsal contribution (level=3)
      ilv=3
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            inode=(istack-1)*3+ilv
            vload(inode)=vload(inode)+vldd(is,iq)
         enddo
      enddo
!--edge contributions
      do il=1,nl
         do is=1,2
            istack=isol(is,il)
            do ilv=1,3
               inode=(istack-1)*3+ilv
               vload(inode)=vload(inode)+vlde(ilv,is,il)
            enddo
         enddo
      enddo
!
!--do LU decomposition
      call lufa(Kstif,mrank,Kluf,isng,ipvt)
      if (isng.ne.0) then
         print *,'ill conditioned stiffness matrix in lussolve, isng=',
     1            isng
         return
      endif
!--solve for concentrations
      if (idebug.eq.1) print *,'calling lusl to solve for scalar'
      call lusl(Kluf,vload,isng,ipvt,mrank,con)
!
!--concentration vector to return
      do i=1,3*ns
         wv(i)=con(i)
      enddo
!
!--we are done, except if debugging flag on
!
      itrouble=0
      if (idebug.eq.1) then
!-- get the inverse matrix
         print *,'calling luiv for inverse matrix'
         call luiv(Kluf,isng,ipvt,mrank,Kinv)
!--get the P-1 norms for the matrices and compute the condition number
         call lup1(Kstif,mrank,pKstif)
         call lup1(Kinv,mrank,pKinv)
         print *,'condition number',pKstif*pKinv
         vload2=matmul(Kstif,con)
         do i=1,mrank
            vtest=abs(vload(i)-vload2(i))/abs(vload(i))
            if (vtest.ge.1d-6) then
               print *,'vload,vload2',vload(i),vload2(i)
               itrouble=1
            endif
         enddo
         Kcheck=matmul(Kstif,Kinv)
         do j=1,mrank
            do i=1,mrank
               if (i.eq.j) then
                  if (abs(1d0-Kcheck(i,j)).ge.1d-6) then
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               else
                  if (abs(Kcheck(i,j)).ge.1d-9) then 
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               endif
            enddo
         enddo
         deallocate(Kinv,Kcheck,vload2)
         print *,'lussolve: itrouble=',itrouble
      endif
!
      deallocate(Kstif,Kluf,vload,con)
!
      return
      end subroutine lussolve   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine lussolve2(wv,Qv,Qd,Qe,vldqv,vldqd,vldqe)
!
! This will solve: (Qv+Qd+Qe)*wv=vldqv+vldqd+vldqe for wv
!  by LU decomposition. This is SURFACE SCALAR diffusion.
! This routine is the homolog to cgssolve2
!
      use lupack
!
      implicit none
!
      real(8) wv(3*NSM)!concentration
      real(8) Qv(4,4,NSM),Qd(4,4,NSM)!ventral/dorsal stiffness matrices
      real(8) Qe(3,2,3,2,NLM)!edge surface stiffness matrices
      real(8) vldqd(4,NSM),vldqv(4,NSM)!ventral/dorsal load vector
      real(8) vldqe(3,2,NLM)!edge surface load vector
!
!  Kstif=(2ns+nl x 2ns+nl) stiffness matrix, Kluf in LU decompositon
!  con=(2ns+nl) concentration vectors, vload=(2ns+nl) load vector
!  Kinv, Kcheck for debugging only 
!  ipvt(2ns+nl) pivot vector
!
      real(8), allocatable, dimension(:,:) :: Kstif, Kluf
      real(8), allocatable, dimension(:,:) :: Kinv, Kcheck
      real(8), allocatable, dimension(:) :: con, vload, vload2
      integer, allocatable, dimension(:) :: ipvt
      integer mrank, iq, il, isng, idebug, itrouble
      integer js, jstack, jlv, jnode, j
      integer is, istack, ilv, inode, i
      integer ne(3,2), js1,js2,il2
      real(8) pKstif, pKinv, vtest
!
      idebug=0
      mrank=2*ns+nl
      allocate (Kstif(mrank,mrank),Kluf(mrank,mrank))
      allocate (con(mrank),vload(mrank))
      allocate (ipvt(mrank))
      if (idebug.eq.1) then
         allocate (Kinv(mrank,mrank),Kcheck(mrank,mrank),vload2(mrank)) 
      endif
!
!--assemble stiffness matrix
!
      Kstif=0d0
!--ventral contributions 
      do iq=1,nq
         do js=1,4
            jnode=isoq(js,iq)
            do is=1,4
               inode=isoq(is,iq)
               Kstif(inode,jnode)=Kstif(inode,jnode)+Qv(is,js,iq)
            enddo
         enddo
      enddo
!--dorsal contributions 
      do iq=1,nq
         do js=1,4
            jnode=isoq(js,iq)+ns
            do is=1,4
               inode=isoq(is,iq)+ns
               Kstif(inode,jnode)=Kstif(inode,jnode)+Qd(is,js,iq)
            enddo
         enddo
      enddo
!--edge contributions
      do il=1,nl
         js1=isol(1,il)    
         js2=isol(2,il)
         il2=ilol(2,il)
!--six nodes on the edge surface
         ne(1,1)=js1 !first stack, ventral
         ne(3,1)=js1+ns !first stack, dorsal
         ne(1,2)=js2 !second stack, ventral
         ne(3,2)=js2+ns !second stack, dorsal
         ne(2,1)=il+2*ns !first stack, middle
         ne(2,2)=il2+2*ns !second stack, middle
         do js=1,2
            do jlv=1,3
               jnode=ne(jlv,js)
               do is=1,2
                  do ilv=1,3
                     inode=ne(ilv,is)
                     Kstif(inode,jnode)=Kstif(inode,jnode)
     1                                  +Qe(ilv,is,jlv,js,il)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--assemble load vector
!
      vload=0d0
!--ventral contributions 
      do iq=1,nq
         do is=1,4
            inode=isoq(is,iq)
            vload(inode)=vload(inode)+vldqv(is,iq)
         enddo
      enddo
!--dorsal contribution 
      do iq=1,nq
         do is=1,4
            inode=isoq(is,iq)+ns
            vload(inode)=vload(inode)+vldqd(is,iq)
         enddo
      enddo
!--edge contributions
      do il=1,nl
         js1=isol(1,il)    
         js2=isol(2,il)
         il2=ilol(2,il)
!--six nodes on the edge surface
         ne(1,1)=js1 !first stack, ventral
         ne(3,1)=js1+ns !first stack, dorsal
         ne(1,2)=js2 !second stack, ventral
         ne(3,2)=js2+ns !second stack, dorsal
         ne(2,1)=il+2*ns !first stack, middle
         ne(2,2)=il2+2*ns !second stack, middle
         do is=1,2
            do ilv=1,3
               inode=ne(ilv,is)
               vload(inode)=vload(inode)+vldqe(ilv,is,il)
            enddo
         enddo
      enddo
!
!--do LU decomposition
      call lufa(Kstif,mrank,Kluf,isng,ipvt)
      if (isng.gt.0) then
         print *,'ill conditioned stiffness matrix in lussolve2, isng=',
     1            isng
         return
      endif
!--solve for concentrations
      call lusl(Kluf,vload,isng,ipvt,mrank,con)
!
!--assemble concentration vector
      do inode=1,2*ns+nl
         wv(inode)=con(inode)
      enddo
!
!--we are done, except if debugging flag on
!
      itrouble=0
      if (idebug.eq.1) then
!-- get the inverse matrix
         call luiv(Kluf,isng,ipvt,mrank,Kinv)
!--get the P-1 norms for the matrices and compute the condition number
         call lup1(Kstif,mrank,pKstif)
         call lup1(Kinv,mrank,pKinv)
         print *,'condition number',pKstif*pKinv
         vload2=matmul(Kstif,con)
         do i=1,mrank
            vtest=abs(vload(i)-vload2(i))/abs(vload(i))
            if (vtest.ge.1d-6) then
               print *,'vload,vload2',vload(i),vload2(i)
               itrouble=1
            endif
         enddo
         Kcheck=matmul(Kstif,Kinv)
         do j=1,mrank
            do i=1,mrank
               if (i.eq.j) then
                  if (abs(1d0-Kcheck(i,j)).ge.1d-6) then
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               else
                  if (abs(Kcheck(i,j)).ge.1d-9) then 
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               endif
            enddo
         enddo
         deallocate(Kinv,Kcheck,vload2)
         print *,'lussolve2: itrouble=',itrouble
      endif
!
      deallocate(Kstif,Kluf,vload,con)
!
      return
      end subroutine lussolve2   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine lussolve3(wv,Qv,Lc,vldqv,vldc)
!
! This will solve: (Qv+Lc)*wv=vldqv+vldc for wv
! by LU decomposition. This is BOTTOM SURFACE SCALAR diffusion.
! This routine is the homolog to cgssolve3
!
      use lupack
!
      implicit none
!
      real(8) wv(3*NSM)!concentration
      real(8) Qv(4,4,NSM)!ventral stiffness matrices
      real(8) Lc(2,2,NLM)!contact line stiffness matrix
      real(8) vldqv(4,NSM)!ventral load vector
      real(8) vldc(2,NLM)!contact line boundary load vector
!
!  Kstif=(ns x ns) stiffness matrix, Kluf in LU decompositon
!  con=(ns) concentration vectors, vload=(ns) load vector
!  Kinv, Kcheck for debugging only 
!  ipvt(ns) pivot vector
!
      real(8), allocatable, dimension(:,:) :: Kstif, Kluf
      real(8), allocatable, dimension(:,:) :: Kinv, Kcheck
      real(8), allocatable, dimension(:) :: con, vload, vload2
      integer, allocatable, dimension(:) :: ipvt
      integer mrank, iq, il, isng, idebug, itrouble
      integer js, jstack, jlv, jnode, j
      integer is, istack, ilv, inode, i
      integer ne(3,2), js1,js2,il2
      real(8) pKstif, pKinv, vtest
!
      idebug=0
      mrank=ns
      allocate (Kstif(mrank,mrank),Kluf(mrank,mrank))
      allocate (con(mrank),vload(mrank))
      allocate (ipvt(mrank))
      if (idebug.eq.1) then
         allocate (Kinv(mrank,mrank),Kcheck(mrank,mrank),vload2(mrank)) 
      endif
!
!--assemble stiffness matrix
!
      Kstif=0d0
!--ventral contributions 
      do iq=1,nq
         do js=1,4
            jnode=isoq(js,iq)
            do is=1,4
               inode=isoq(is,iq)
               Kstif(inode,jnode)=Kstif(inode,jnode)+Qv(is,js,iq)
            enddo
         enddo
      enddo
!--contact line contributions
      do il=1,nl
         do js=1,2
            jnode=isol(js,il) !contact line node is ventral
            do is=1,2
               inode=isol(is,il) !contact line node is ventral
               Kstif(inode,jnode)=Kstif(inode,jnode)+Lc(is,js,il)
            enddo
         enddo
      enddo
!
!--assemble load vector
!
      vload=0d0
!--ventral contributions 
      do iq=1,nq
         do is=1,4
            inode=isoq(is,iq)
            vload(inode)=vload(inode)+vldqv(is,iq)
         enddo
      enddo
!--contact line contributions
      do il=1,nl
         do is=1,2
            inode=isol(is,il)
            vload(inode)=vload(inode)+vldc(is,il)
         enddo
      enddo
!
!--do LU decomposition
      call lufa(Kstif,mrank,Kluf,isng,ipvt)
      if (isng.gt.0) then
         print *,'ill conditioned stiffness matrix in lussolve2, isng=',
     1            isng
         return
      endif
!--solve for concentrations
      call lusl(Kluf,vload,isng,ipvt,mrank,con)
!
!--assemble concentration vector
      do inode=1,ns
         wv(inode)=con(inode)
      enddo
!
!--we are done, except if debugging flag on
!
      itrouble=0
      if (idebug.eq.1) then
!-- get the inverse matrix
         call luiv(Kluf,isng,ipvt,mrank,Kinv)
!--get the P-1 norms for the matrices and compute the condition number
         call lup1(Kstif,mrank,pKstif)
         call lup1(Kinv,mrank,pKinv)
         print *,'condition number',pKstif*pKinv
         vload2=matmul(Kstif,con)
         do i=1,mrank
            vtest=abs(vload(i)-vload2(i))/abs(vload(i))
            if (vtest.ge.1d-6) then
               print *,'vload,vload2',vload(i),vload2(i)
               itrouble=1
            endif
         enddo
         Kcheck=matmul(Kstif,Kinv)
         do j=1,mrank
            do i=1,mrank
               if (i.eq.j) then
                  if (abs(1d0-Kcheck(i,j)).ge.1d-6) then
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               else
                  if (abs(Kcheck(i,j)).ge.1d-9) then 
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               endif
            enddo
         enddo
         deallocate(Kinv,Kcheck,vload2)
         print *,'lussolve2: itrouble=',itrouble
      endif
!
      deallocate(Kstif,Kluf,vload,con)
!
      return
      end subroutine lussolve3   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine luvsolve(wv,Q,Lv,Ld,Le,vldq,vldv,vldd,vlde)
!
! This will solve: (Q+Lv+Ld+Le)*wv=vldq+vldv+vldd_vlde for wv
!  by LU decomposition. wv is a VECTOR field.
! This routine is the homolog to cgvsolve
!
      use lupack
!
      implicit none
 
!
      real(8) wv(3,3*NSM)!vector field
      real(8) Q(3,3,4,3,3,4,NSM)!element stiffness matrix
      real(8) Lv(3,4,3,4,NSM),Ld(3,4,3,4,NSM)!ventral/dorsal stiffness
      real(8) Le(3,3,2,3,3,2,NLM)!edge boundary stiffness
      real(8) vldq(3,3,4,NSM)!element load vector
      real(8) vldv(3,4,NSM), vldd(3,4,NSM)!ventral/dorsal boundary load
      real(8) vlde(3,3,2,NLM)!edge boundary load vector
!
      integer mrank, iq, il, isng, idebug, itrouble
      integer js, jstack, jlv, jnode, jx, j
      integer is, istack, ilv, inode, ix, i
      real(8) pKstif, pKinv, vtest
!
!  Kstif=(9ns x 9ns) stiffness matrix, Kluf in LU decompositon
!  vel=(9ns) velocity vector, vload=(9ns) load vector
!  Kinv, Kcheck, vload2 for debugging only 
!  ipvt(9ns) pivot vector
!
      real(8), allocatable, dimension(:,:) :: Kstif, Kluf
      real(8), allocatable, dimension(:,:) :: Kinv, Kcheck
      real(8), allocatable, dimension(:) :: vel, vload, vload2
      integer, allocatable, dimension(:) :: ipvt
!
      idebug=1
      mrank=9*ns
      allocate (Kstif(mrank,mrank),Kluf(mrank,mrank))
      allocate (vel(mrank),vload(mrank)) 
      allocate (ipvt(mrank))
      if (idebug.eq.1) then
         allocate (Kinv(mrank,mrank),Kcheck(mrank,mrank),vload2(mrank))
      endif
!
!--assemble stiffness matrix
      if (idebug.eq.1) print *,'assembling stiffness matrix'
!
!--interior contributions
      Kstif=0d0
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            do jlv=1,3
               jnode=(jstack-1)*3+jlv
               do jx=1,3
                  j=(jnode-1)*3+jx
                  do is=1,4
                     istack=isoq(is,iq)
                     do ilv=1,3
                        inode=(istack-1)*3+ilv
                        do ix=1,3
                           i=(inode-1)*3+ix
                           Kstif(i,j)=Kstif(i,j)+
     1                                Q(ix,ilv,is,jx,jlv,js,iq)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!--ventral contributions (level=1)
      jlv=1
      ilv=1
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            jnode=(jstack-1)*3+jlv
            do jx=1,3
               j=(jnode-1)*3+jx
               do is=1,4
                  istack=isoq(is,iq)
                  inode=(istack-1)*3+ilv
                  do ix=1,3
                     i=(inode-1)*3+ix
                     Kstif(i,j)=Kstif(i,j)+Lv(ix,is,jx,js,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
!--dorsal contributions (level=3)
      jlv=3
      ilv=3
      do iq=1,nq
         do js=1,4
            jstack=isoq(js,iq)
            jnode=(jstack-1)*3+jlv
            do jx=1,3
               j=(jnode-1)*3+jx
               do is=1,4
                  istack=isoq(is,iq)
                  inode=(istack-1)*3+ilv
                  do ix=1,3
                     i=(inode-1)*3+ix
                     Kstif(i,j)=Kstif(i,j)+Ld(ix,is,jx,js,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
!--edge contributions
      do il=1,nl
         do js=1,2
            jstack=isol(js,il)
            do jlv=1,3
               jnode=(jstack-1)*3+jlv
               do jx=1,3
                  j=(jnode-1)*3+jx
                  do is=1,2
                     istack=isol(is,il)
                     do ilv=1,3
                        inode=(istack-1)*3+ilv
                        do ix=1,3
                           i=(inode-1)*3+ix
                           Kstif(i,j)=Kstif(i,j)+
     1                                Le(ix,ilv,is,jx,jlv,js,il)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--assemble load vector
!
      if (idebug.eq.1) print *,'assembling load vector'
!--interior contributions
      vload=0d0
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            do ilv=1,3
               inode=(istack-1)*3+ilv
               do ix=1,3
                  i=(inode-1)*3+ix
                  vload(i)=vload(i)+vldq(ix,ilv,is,iq)
               enddo
            enddo
         enddo
      enddo
!--ventral contributions (level=1)
      ilv=1
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            inode=(istack-1)*3+ilv
            do ix=1,3
               i=(inode-1)*3+ix
               vload(i)=vload(i)+vldv(ix,is,iq)
            enddo
         enddo
      enddo
!--dorsal contribution (level=3)
      ilv=3
      do iq=1,nq
         do is=1,4
            istack=isoq(is,iq)
            inode=(istack-1)*3+ilv
            do ix=1,3
               i=(inode-1)*3+ix
               vload(i)=vload(i)+vldd(ix,is,iq)
            enddo
         enddo
      enddo
c--edge contributions
      do il=1,nl
         do is=1,2
            istack=isol(is,il)
            do ilv=1,3
               inode=(istack-1)*3+ilv
               do ix=1,3
                  i=(inode-1)*3+ix
                  vload(i)=vload(i)+vlde(ix,ilv,is,il)
               enddo
            enddo
         enddo
      enddo
!
!--do LU decomposition
!
      if (idebug.eq.1) print *,'begin LU decomp. of stiffness matrix'
      call lufa(Kstif,mrank,Kluf,isng,ipvt)
      if (isng.gt.0) then
         print *,'ill conditioned stiffness matrix in luvsolve, isng=',
     1            isng
         return
      endif
!--solve for velocities
      if (idebug.eq.1) print *,'solving for velocities'
      call lusl(Kluf,vload,isng,ipvt,mrank,vel)
!
!--assemble velocity vector
      do is=1,ns
         do ilv=1,3
            inode=(is-1)*3+ilv
            do ix=1,3
               i=(inode-1)*3+ix
               wv(ix,inode)=vel(i)
            enddo
         enddo
      enddo
!
!--we are done, except if debugging flag on
!
      if (idebug.eq.1) then
         itrouble=0
!-- get the inverse matrix
         call luiv(Kluf,isng,ipvt,mrank,Kinv)
!--get the P-1 norms for the matrices and compute the condition number
         call lup1(Kstif,mrank,pKstif)
         call lup1(Kinv,mrank,pKinv)
         print *,'condition number',pKstif*pKinv
         vload2=matmul(Kstif,vel)
         do i=1,mrank
            vtest=abs(vload(i)-vload2(i))/abs(vload(i))
            if (vtest.ge.1d-6) then
               print *,'vload,vload2',vload(i),vload2(i)
               itrouble=1
            endif
         enddo
         Kcheck=matmul(Kstif,Kinv)
         do j=1,mrank
            do i=1,mrank
               if (i.eq.j) then
                  if (abs(1d0-Kcheck(i,j)).ge.1d-6) then
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               else
                  if (abs(Kcheck(i,j)).ge.1d-9) then 
                     print *,'i,j,Kcheck',i,j,Kcheck(i,j)
                     itrouble=1
                  endif
               endif
            enddo
         enddo
         deallocate(Kinv,Kcheck,vload2)
         print *,'luvsolve: itrouble=',itrouble
      endif
!
      deallocate(Kstif,Kluf,vload,vel)
!
      return
      end subroutine luvsolve   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   END OF LU DECOMPOSITION ROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  BEGIN CONJUGATE GRADIENT ROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  begin CG routines for VOLUME SCALAR diffusion
!
!!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! 
!
      subroutine cgssolve(wv,Q,Lv,Ld,Le,C,vldq,vldv,vldd,vlde,mcyc,icyc)
! 
! This will solve: (Q+Lv+Ld+Le)*wv=vldq+vldv+vldd_vlde for wv
!  by the conjugate gradient method. This is VOLUME SCALAR diffusion.
! This routine is the homolog to lussolve.
!
      real(8) wv(3*NSM)!concentration
      real(8) Q(3,4,3,4,NSM)!element stiffness matrix
      real(8) Lv(4,4,NSM),Ld(4,4,NSM)!ventral/dorsal boundary stiffness
      real(8) Le(3,2,3,2,NLM)!edge boundary stiffness matrix
      real(8) vldq(3,4,NSM)!element load vector
      real(8) vldv(4,NSM),vldd(4,NSM)!ventral/dorsal boundary loads
      real(8) vlde(3,2,NLM)!edge boundary load vector
      real(8) C(3*NSM)! conditioning matrix
      integer mcyc!maximum number of CG iterations
      integer icyc !actual number of iterations
!
      real(8) gv(3*NSM),dv(3*NSM),fv(3*NSM)
      real(8) del0,del1,beta,ddotf,tau
      integer inode,i
!
!--compute the initial gradient field.
      call cgsmxmul(gv,Q,Lv,Ld,Le,wv)
      call cgsasmbl(gv,vldq,vldv,vldd,vlde)!gv(in)=gv(in)-vld(in)
!--gv is now the residual e.g. gv=Q*wv-vld 
!
      do inode=1,3*ns
         dv(inode)=-gv(inode)*C(inode)!init the search direction field.
      enddo
!--dv is now the pre-conditioned residual
!
      del0=0.0
      do inode=1,3*ns
         del0=del0-dv(inode)*gv(inode)!negative inner prod of dv*gv
      enddo
!
      icyc=0!initialize the cycle count
   50 icyc=icyc+1!return here to start a new cycle.
!     print *,'del0,icyc,mcyc',del0,icyc,mcyc
      if(del0.le.0.0.or.icyc.gt.mcyc)return!solution reached
!
!--compute the congugate of the search direction field.
      call cgsmxmul(fv,Q,Lv,Ld,Le,dv)
!
      ddotf=0.0
      do inode=1,3*ns
         ddotf=ddotf+dv(inode)*fv(inode)!compute the energy norm; ddotf
      enddo
!     print *,'ddotf=',ddotf
!
      if(ddotf.le.vtiny)return!exact solution,(normal but rare).
!
      tau=del0/ddotf!compute tau
!     print *,'del0,ddotf,tau,icyc=',del0,ddotf,tau,icyc
      do inode=1,3*ns
         wv(inode)=wv(inode)+tau*dv(inode)!new solution vector
         gv(inode)=gv(inode)+tau*fv(inode)!new gradient field
      enddo
!
      do inode=1,3*ns
         fv(inode)=gv(inode)*C(inode)!precondition the new gradient field
      enddo
!
      del1=0.0
      do inode=1,3*ns
         del1=del1+gv(inode)*fv(inode)!compute inner product of gv&fv
      enddo
!
      beta=del1/del0
      do inode=1,3*ns
         dv(inode)=-fv(inode)+beta*dv(inode)!new search direction field.
      enddo
!
      del0=del1!reinitialize del0
      go to 50! start a new cycle.
      end subroutine cgssolve
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgsmxmul(f,Q,Lv,Ld,Le,w)
!
! This subroutine does computes the matrix multiplication:
!    f = (Q+Lv+Ld+Le)*w
!  this is for SCALAR VOLUME diffusion
!
      implicit none
!
      real(8) f(3*NSM),w(3*NSM)
      real(8) Q(3,4,3,4,NSM)
      real(8) Lv(4,4,NSM), Ld(4,4,NSM), Le(3,2,3,2,NLM)
!
! Temporary sums: vq = Q*wv; vv = Lv*wv, etc.
      real(8) vq(3,4,NSM), vv(4,NSM), vd(4,NSM), ve(3,2,NLM) 
      real(8) wn,wi 
      integer iq, isn,lvn,is,ilv,il,nstack, inode,istack
!
      vq=0d0
!--loop over elements
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               wn=w(inode)
               do is=1,4
                  do ilv=1,3
                     vq(ilv,is,iq)=vq(ilv,is,iq)-
     1                             Q(ilv,is,lvn,isn,iq)*wn
!--this is the same since Q is symmetric
!                    istack=isoq(is,iq)
!                    inode=(istack-1)*3+ilv
!                    wi=w(inode)
!                    vq(lvn,isn,iq)=vq(lvn,isn,iq)-
!    1                              Q(ilv,is,lvn,isn,iq)*wi
                  enddo
               enddo
            enddo
         enddo
      enddo
!--loop over ventral/dorsal boundaries
      vv=0d0
      vd=0d0
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=(nstack-1)*3+1
            wn=w(inode)
            do is=1,4
               vv(is,iq)=vv(is,iq)-Lv(is,isn,iq)*wn
            enddo
            inode=(nstack-1)*3+3
            wn=w(inode)
            do is=1,4
               vd(is,iq)=vd(is,iq)-Ld(is,isn,iq)*wn
            enddo
         enddo
      enddo 
!--loop over edge boundaries
      ve=0d0
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               wn=w(inode)
               do is=1,2
                  do ilv=1,3
                     ve(ilv,is,il)=ve(ilv,is,il)-
     1                             Le(ilv,is,lvn,isn,il)*wn
                  enddo
               enddo
            enddo
         enddo
      enddo
      f=0d0
      call cgsasmbl(f,vq,vv,vd,ve)!assemble f by adding the temp sums
      return
      end subroutine cgsmxmul
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgsasmbl(vn,vq,vv,vd,ve)
!
!  this computes vn=-(vq+vv+vd+ve)
!  this is for SCALAR VOLUME diffusion
!
      implicit none
!
      real(8) vn(3*NSM)
      real(8) vq(3,4,NSM), vv(4,NSM), vd(4,NSM), ve(3,2,NLM)
!
      integer iq,il,isn,nstack,inode,lvn
!
!--element contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               vn(inode)=vn(inode)-vq(lvn,isn,iq)
            enddo
         enddo
      enddo
!--ventral/dorsal contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=(nstack-1)*3+1
            vn(inode)=vn(inode)-vv(isn,iq)
            inode=(nstack-1)*3+3
            vn(inode)=vn(inode)-vd(isn,iq)
         enddo
      enddo
!--edge contributions
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               vn(inode)=vn(inode)-ve(lvn,isn,il)
            enddo
         enddo
      enddo
      return
      end subroutine cgsasmbl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgscndmx(C,Q,Lv,Ld,Le)
!
! compute the conditoning matrix C used by the conjugate
! gradient solver for SCALAR VOLUME diffusion
!
!   C=diag(Q+Lv+Ld+Le)^-1
!
      implicit none
!
      real(8) C(3*NSM),Q(3,4,3,4,NSM)
      real(8) Lv(4,4,NSM),Ld(4,4,NSM),Le(3,2,3,2,NLM)
!
      integer iq,isn,istack,lvn,inode,il
!
      C=0d0
!
!--interior contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               inode=(istack-1)*3+lvn
               C(inode)=C(inode)+Q(lvn,isn,lvn,isn,iq)
            enddo
         enddo
      enddo
!
!--boundary contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
!--ventral nodes
            inode=(istack-1)*3+1
            C(inode)=C(inode)+Lv(isn,isn,iq)
!--dorsal nodes
            inode=(istack-1)*3+3
            C(inode)=C(inode)+Ld(isn,isn,iq)
         enddo
      enddo
!--edge nodes
      do il=1,nl
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               inode=(istack-1)*3+lvn
               C(inode)=C(inode)+Le(lvn,isn,lvn,isn,il)
            enddo
         enddo
      enddo
!--invert (diagonal) matrix
      do inode=1,3*ns
         C(inode)=1d0/C(inode)
      enddo
      end subroutine cgscndmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  end CG routines for VOLUME SCALAR diffusion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  begin CG routines for SURFACE SCALAR diffusion
!
!!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! 
!
      subroutine cgssolve2(wv,Qv,Qd,Qe,C,vldqv,vldqd,vldqe,mcyc)
!
! This will solve: (Qv+Qd+Qe)*wv=vldqv+vldqd+vldqe for wv by the 
! conjugate gradient method. This is for SURFACE SCALAR diffusion.
! This routine is the homolog to lussolve2
!
      implicit none
!
      real(8) wv(3*NSM)!concentration
      real(8) Qv(4,4,NSM),Qd(4,4,NSM)!ventral/dorsal stiffness matrices
      real(8) Qe(3,2,3,2,NLM)!edge surface stiffness matrices
      real(8) vldqd(4,NSM),vldqv(4,NSM)!ventral/dorsal load vector
      real(8) vldqe(3,2,NLM)!edge surface load vector
      real(8) C(3*NSM)! conditioning matrix
      integer mcyc!maximum number of CG iterations
!
      real(8) gv(3*NSM),dv(3*NSM),fv(3*NSM)
      real(8) del0,del1,beta,ddotf,tau
      integer icyc,inode
!
!--compute the initial gradient field.
      call cgsmxmul2(gv,Qv,Qd,Qe,wv)
      call cgsasmbl2(gv,vldqv,vldqd,vldqe)!gv(in)=gv(in)-vld(in)
!--gv is now the residual e.g. gv=Q*wv-vld 
!
      do inode=1,2*ns+nl
         dv(inode)=-gv(inode)*C(inode)!init the search direction field.
      enddo
!--dv is now the pre-conditioned residual
!
      del0=0.0
      do inode=1,2*ns+nl
         del0=del0-dv(inode)*gv(inode)!negative inner prod of dv*gv
      enddo
!     print *,'del0=',del0
!
      icyc=0!initialize the cycle count
   50 icyc=icyc+1!return here to start a new cycle.
!     print *,'del0,icyc,mcyc',del0,icyc,mcyc
      if(del0.le.0.0.or.icyc.gt.mcyc)return!solution reached
!
!--compute the congugate of the search direction field.
      call cgsmxmul2(fv,Qv,Qd,Qe,dv)
      ddotf=0.0
      do inode=1,2*ns+nl
         ddotf=ddotf+dv(inode)*fv(inode)!compute the energy norm; ddotf
      enddo
!     print *,'ddotf=',ddotf
!
      if(ddotf.le.vtiny*vtiny)return!exact solution,(normal but rare).
!
      tau=del0/ddotf!compute tau
!     print *,'del0,ddotf,tau=',del0,ddotf,tau
      do inode=1,2*ns+nl
         wv(inode)=wv(inode)+tau*dv(inode)!new solution vector
         gv(inode)=gv(inode)+tau*fv(inode)!new gradient field
      enddo
!
      do inode=1,2*ns+nl
         fv(inode)=gv(inode)*C(inode)!precondition the new gradient field
      enddo
!
      del1=0.0
      do inode=1,2*ns+nl
         del1=del1+gv(inode)*fv(inode)!compute inner product of gv&fv
      enddo
!
      beta=del1/del0
      do inode=1,2*ns+nl
         dv(inode)=-fv(inode)+beta*dv(inode)!new search direction field.
      enddo
!
      del0=del1!reinitialize del0
      go to 50! start a new cycle.
      end subroutine cgssolve2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgsmxmul2(f,Qv,Qd,Qe,w)
!
! This subroutine does computes the matrix multiplication:
!    f = (Qv+Qd+Qe)*w
!  this is for SCALAR SURFACE diffusion
!
      implicit none
!
      real(8) f(3*NSM),w(3*NSM)
      real(8) Qe(3,2,3,2,NLM),Qv(4,4,NSM),Qd(4,4,NSM)
!
! Temporary sums: vv = Qv*wv, etc.
      real(8) vv(4,NSM), vd(4,NSM), ve(3,2,NLM)
!
      real(8) wn 
      integer iq, isn,lvn,is,ilv,il,il2,nstack, inode
!
      vv=0d0
!--loop over ventral surfaces
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack
            wn=w(inode)
            do is=1,4
               vv(is,iq)=vv(is,iq)-Qv(is,isn,iq)*wn
            enddo
         enddo
      enddo
      vd=0d0
!--loop over dorsal surfaces
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack+ns
            wn=w(inode)
            do is=1,4
               vd(is,iq)=vd(is,iq)-Qd(is,isn,iq)*wn
            enddo
         enddo
      enddo
!--loop over edges surfaces
      ve=0d0
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
!--ventral node
            inode=nstack
            wn=w(inode)
            do is=1,2
               do ilv=1,3
                  ve(ilv,is,il)=ve(ilv,is,il)-Qe(ilv,is,1,isn,il)*wn
               enddo
            enddo
!--dorsal node
            inode=nstack+ns
            wn=w(inode)
            do is=1,2
               do ilv=1,3
                  ve(ilv,is,il)=ve(ilv,is,il)-Qe(ilv,is,3,isn,il)*wn
               enddo
            enddo
         enddo
!--middle nodes 1st, isol(1,il)
         inode=il+2*ns
         wn=w(inode)
         do is=1,2
            do ilv=1,3
               ve(ilv,is,il)=ve(ilv,is,il)-Qe(ilv,is,2,1,il)*wn
            enddo
         enddo
         il2=ilol(2,il)
         inode=il2+2*ns
         wn=w(inode)
         do is=1,2
            do ilv=1,3
               ve(ilv,is,il)=ve(ilv,is,il)-Qe(ilv,is,2,2,il)*wn
            enddo
         enddo
      enddo
      f=0d0
!
      call cgsasmbl2(f,vv,vd,ve)!assemble f by adding the temp sums
!
      return
      end subroutine cgsmxmul2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgsasmbl2(vn,vv,vd,ve)
!
!  this computes vn=-(vv+vd+ve)
!  this is for SCALAR SURFACE diffusion
!
      implicit none
!
      real(8) vn(3*NSM)
      real(8) vv(4,NSM), vd(4,NSM), ve(3,2,NLM)
!
      integer iq,il,il2,isn,nstack,inode,lvn
!
!--ventral element contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack
            vn(inode)=vn(inode)-vv(isn,iq)
         enddo
      enddo
!--dorsal element contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack+ns
            vn(inode)=vn(inode)-vd(isn,iq)
         enddo
      enddo
!--edge contributions
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
!--bottom node
            inode=nstack
            vn(inode)=vn(inode)-ve(1,isn,il)
!--top node
            inode=nstack+ns
            vn(inode)=vn(inode)-ve(3,isn,il)
         enddo
!--middle nodes
         inode=il+2*ns
         vn(inode)=vn(inode)-ve(2,1,il)
         il2=ilol(2,il)
         inode=il2+2*ns
         vn(inode)=vn(inode)-ve(2,2,il)
      enddo
!
      return
      end subroutine cgsasmbl2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgscndmx2(C,Qsv,Qsd,Qse)
!
      implicit none
!
! compute the conditoning matrix C used by the conjugate
! gradient solver for SCALAR SURFACE diffusion
!
!   C=diag(Qsv+Qsd+Qse)^-1
!
      real(8) C(3*NSM)
      real(8) Qsv(4,4,NSM),Qsd(4,4,NSM), Qse(3,2,3,2,NLM)
!
      integer iq,isn,istack,lvn,inode,il,il2
!
      C=0d0
!
!--ventral contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
            inode=istack
            C(inode)=C(inode)+Qsv(isn,isn,iq)
         enddo
      enddo
!--dorsal contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
            inode=istack+ns
            C(inode)=C(inode)+Qsd(isn,isn,iq)
         enddo
      enddo
!--edge contributions 
      do il=1,nl
         do isn=1,2
            istack=isol(isn,il)
!--ventral node
            inode=istack
            C(inode)=C(inode)+Qse(1,isn,1,isn,il)
!--dorsal node
            inode=istack+ns
            C(inode)=C(inode)+Qse(3,isn,3,isn,il)
         enddo
!--middle node - 1st stack
         inode=il+2*ns
         C(inode)=C(inode)+Qse(2,1,2,1,il)
!--middle node - second stack
         il2=ilol(2,il)
         inode=il2+2*ns
         C(inode)=C(inode)+Qse(2,2,2,2,il)
      enddo
!--invert (diagonal) matrix
      do inode=1,2*ns+nl
         C(inode)=1d0/C(inode)
      enddo
!
      return
      end subroutine cgscndmx2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  begin CG routines for BOTTOM SURFACE SCALAR diffusion
!
!!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! 
!
      subroutine cgssolve3(wv,Qv,Lc,C,vldqv,vldc,mcyc)
!
! This will solve: (Qv+Lc)*wv=vldqv+vldc for wv by the 
! conjugate gradient method. This is for BOTTOM SURFACE SCALAR diffusion.
! This routine is the homolog to lussolve3
!
      implicit none
!
      real(8) wv(3*NSM)!concentration
      real(8) Qv(4,4,NSM)!ventral stiffness matrices
      real(8) Lc(2,2,NLM)! contact line stiffness matrices
      real(8) vldqv(4,NSM)!ventral load vector
      real(8) vldc(2,NLM)!contact line boundary load vector
      real(8) C(3*NSM)! conditioning matrix
      integer mcyc!maximum number of CG iterations
!
      real(8) gv(3*NSM),dv(3*NSM),fv(3*NSM)
      real(8) del0,del1,beta,ddotf,tau
      integer icyc,inode
!
!--compute the initial gradient field.
      call cgsmxmul3(gv,Qv,Lc,wv)
      call cgsasmbl3(gv,vldqv,vldc)!gv(in)=gv(in)-vld(in)
!--gv is now the residual e.g. gv=Q*wv-vld 
!
      do inode=1,ns
         dv(inode)=-gv(inode)*C(inode)!init the search direction field.
      enddo
!--dv is now the pre-conditioned residual
!
      del0=0.0
      do inode=1,ns
         del0=del0-dv(inode)*gv(inode)!negative inner prod of dv*gv
      enddo
!     print *,'del0=',del0
!
      icyc=0!initialize the cycle count
   50 icyc=icyc+1!return here to start a new cycle.
!     print *,'del0,icyc,mcyc',del0,icyc,mcyc
      if(del0.le.0.0.or.icyc.gt.mcyc)return!solution reached
!
!--compute the congugate of the search direction field.
      call cgsmxmul3(fv,Qv,Lc,dv)
      ddotf=0.0
      do inode=1,ns
         ddotf=ddotf+dv(inode)*fv(inode)!compute the energy norm; ddotf
      enddo
!     print *,'ddotf=',ddotf
!
      if(ddotf.le.vtiny*vtiny)return!exact solution,(normal but rare).
!
      tau=del0/ddotf!compute tau
!     print *,'del0,ddotf,tau=',del0,ddotf,tau
      do inode=1,ns
         wv(inode)=wv(inode)+tau*dv(inode)!new solution vector
         gv(inode)=gv(inode)+tau*fv(inode)!new gradient field
      enddo
!
      do inode=1,ns
         fv(inode)=gv(inode)*C(inode)!precondition the new gradient field
      enddo
!
      del1=0.0
      do inode=1,ns
         del1=del1+gv(inode)*fv(inode)!compute inner product of gv&fv
      enddo
!
      beta=del1/del0
      do inode=1,ns
         dv(inode)=-fv(inode)+beta*dv(inode)!new search direction field.
      enddo
!
      del0=del1!reinitialize del0
      go to 50! start a new cycle.
!
      end subroutine cgssolve3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgsmxmul3(f,Qv,Lc,w)
!
! This subroutine does computes the matrix multiplication:
!    f = (Qv+Lc)*w
!  this is for BOTTOM SCALAR SURFACE diffusion
!
      implicit none
!
      real(8) f(3*NSM),w(3*NSM)
      real(8) Qv(4,4,NSM)
      real(8) Lc(2,2,NLM)
!
! Temporary sums: vv = Qv*wv, etc.
      real(8) vv(4,NSM), vc(2,NLM)
!
      real(8) wn 
      integer iq, isn,lvn,is,ilv,il,il2,nstack, inode
!
      vv=0d0
!--loop over ventral surfaces
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack
            wn=w(inode)
            do is=1,4
               vv(is,iq)=vv(is,iq)-Qv(is,isn,iq)*wn
            enddo
         enddo
      enddo
!--loop over contact line
      vc=0d0
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
            inode=nstack
            wn=w(inode)
            do is=1,2
               vc(is,il)=vc(is,il)-Lc(is,isn,il)*wn
            enddo
         enddo
      enddo
      f=0d0
!
      call cgsasmbl3(f,vv,vc)!assemble f by adding the temp sums
!
      return
      end subroutine cgsmxmul3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgsasmbl3(vn,vv,vc)
!
!  this computes vn=-(vv+vc)
!  this is for SCALAR SURFACE diffusion
!
      implicit none
!
      real(8) vn(3*NSM)
      real(8) vv(4,NSM)
      real(8) vc(2,NLM)
!
      integer iq,il,il2,isn,nstack,inode,lvn
!
!--ventral element contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=nstack
            vn(inode)=vn(inode)-vv(isn,iq)
         enddo
      enddo
!--contact line contributions
      do il=1,nl
         do isn=1,2
            inode=isol(isn,il)
            vn(inode)=vn(inode)-vc(isn,il)
         enddo
      enddo
!
      return
      end subroutine cgsasmbl3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgscndmx3(C,Qsv,Lc)
!
      implicit none
!
! compute the conditoning matrix C used by the conjugate
! gradient solver for SCALAR SURFACE diffusion
!
!   C=diag(Qsv+Lc)^-1
!
      real(8) C(3*NSM)
      real(8) Qsv(4,4,NSM)
      real(8) Lc(2,2,NLM)
!
      integer iq,isn,istack,lvn,inode,il,il2
!
      C=0d0
!
!--ventral contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
            inode=istack
            C(inode)=C(inode)+Qsv(isn,isn,iq)
         enddo
      enddo
!--contact line contributions
      do il=1,nl
         do isn=1,2
            inode=isol(isn,il) !contact line nodes are ventral
            C(inode)=C(inode)+Lc(isn,isn,il)
         enddo
      enddo
!--invert (diagonal) matrix
      do inode=1,ns
         C(inode)=1d0/C(inode)
      enddo
!
      return
      end subroutine cgscndmx3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  end CG routines for BOTTOM SURFACE SCALAR diffusion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  begin CG routines for VECTORS 
!
!!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! !!!!!! 
!
      subroutine cgvsolve(wv,Q,Lv,Ld,Le,C,vldq,vldv,vldd,vlde,mcyc,icyc)
! 
! This will solve: (Q+Lv+Ld+Le)*wv=vldq+vldv+vldd+vlde for wv
!  (where wv is a VECTOR field) by the conjugate gradient method. 
! This routine is the homolog to luvsolve.
!
      implicit none
!
      real(8) wv(3,3*NSM)!velocity field
      real(8) Q(3,3,4,3,3,4,NSM)!element stiffness matrix
      real(8) Lv(3,4,3,4,NSM),Ld(3,4,3,4,NSM)!ventral/dorsal stiffness
      real(8) Le(3,3,2,3,3,2,NLM)!edge boundary stiffness matrix
      real(8) vldq(3,3,4,NSM)!element load vector
      real(8) vldv(3,4,NSM), vldd(3,4,NSM)!ventral/dorsal boundary loads
      real(8) vlde(3,3,2,NLM)!edge boundary load vector
      real(8) C(3,3,3*NSM)!conditioning matrix
      integer mcyc !maximum number of CG iterations
      integer icyc !actual number of CG iterations
!
      real(8) del0,del1,beta,ddotf,tau,gvji
      integer inode,jx,ix
      integer isn,lvn
      real(8) gv(3,3*NSM),dv(3,3*NSM),fv(3,3*NSM)
!
      icyc=0!init the cycle count.
!
!--compute the residual gv = (Q+L) wv - vld
      call cgvmxmul(gv,Q,Lv,Ld,Le,wv)
      call cgvasmbl(gv,vldq,vldv,vldd,vlde)
!
!--use the residual gv to compute the pre-conditioned residual dv
      dv=0d0
      do inode=1,3*ns
         do jx=1,3
            gvji=gv(jx,inode)
            do ix=1,3
               dv(ix,inode)=dv(ix,inode)-C(ix,jx,inode)*gvji
            enddo
         enddo
      enddo
!
      del0=0.0
      do inode=1,3*ns
         do ix=1,3
            del0=del0-dv(ix,inode)*gv(ix,inode)!negative inner prod of dv*gv
         enddo
      enddo
!
!--cycle starts here
   50 continue
      icyc=icyc+1!increase the cycle count.
!     print *,'icyc,del0',icyc,del0
      if(del0.le.0.0.or.icyc.gt.mcyc) then
!        print *,'cgvsolve: del0/icyc return del0,ddotf',del0,ddotf
         return!primary exit mode
      endif
!
!--compute the conjugate of the search direction field.
      call cgvmxmul(fv,Q,Lv,Ld,Le,dv)
!
      ddotf=0.0!
      do inode=1,3*ns
         do ix=1,3
            ddotf=ddotf+dv(ix,inode)*fv(ix,inode)
         enddo
      enddo
!     print *,'ddotf',ddotf
!
!--check if exact solution reached.
      if(ddotf.le.vtiny) then
!       print *,'cgvsolve: ddotf return del0,ddotf',del0,ddotf
        return!(a rare but normal exit)
      endif
!
      tau=del0/ddotf!compute tau
      do inode=1,3*ns
         do ix=1,3
            wv(ix,inode)=wv(ix,inode)+tau*dv(ix,inode)!new solution
            gv(ix,inode)=gv(ix,inode)+tau*fv(ix,inode)!new gradient
         enddo
      enddo
!     print *,'del0,icyc,mcyc',del0,icyc,mcyc
!     print *,'ddotf,tau',ddotf,tau
!
      fv=0d0
      do inode=1,3*ns
         do jx=1,3
            gvji=gv(jx,inode)
            do ix=1,3
               fv(ix,inode)=fv(ix,inode)+C(ix,jx,inode)*gvji
            enddo
         enddo
      enddo
!
      del1=0.0
      do inode=1,3*ns
         do ix=1,3
            del1=del1+gv(ix,inode)*fv(ix,inode)!inner product gv*fv
         enddo
      enddo
!
!--get new search direction
      beta=del1/del0
      do inode=1,3*ns
         do ix=1,3
            dv(ix,inode)=-fv(ix,inode)+beta*dv(ix,inode)
         enddo
      enddo
!
      del0=del1!reinitialize del0
      go to 50!start a new cycle.
!
      end subroutine cgvsolve
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgvmxmul(f,Q,Lv,Ld,Le,w)
!
! This subroutine does computes the matrix multiplication:
!    f = (Q+Lv+Ld+Le)*w
!  this is for SCALAR VOLUME diffusion
!
      implicit none
!
      real(8) f(3,3*NSM)
      real(8) Q(3,3,4,3,3,4,NSM)
      real(8) Lv(3,4,3,4,NSM), Ld(3,4,3,4,NSM), Le(3,3,2,3,3,2,NLM)
      real(8) w(3,3*NSM)
!
! Temporary sums: vq= Q*wv, vv=Lv*wv, etc.
      real(8) vq(3,3,4,NSM),vv(3,4,NSM),vd(3,4,NSM),ve(3,3,2,NLM)
      real(8) wni
      integer iq,isn,nstack,lvn,inode,inx,js,jlv,jx,il
!
!--element contributions
      vq=0d0
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               do inx=1,3
                  wni=w(inx,inode)
                  do js=1,4
                     do jlv=1,3
                        do jx=1,3
                           vq(jx,jlv,js,iq)=vq(jx,jlv,js,iq)
     1                          -Q(jx,jlv,js,inx,lvn,isn,iq)*wni
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!--ventral/dorsal boundary contributions
      vv=0d0
      vd=0d0
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=(nstack-1)*3+1
            do inx=1,3
               wni=w(inx,inode)
               do js=1,4
                  do jx=1,3
                     vv(jx,js,iq)=vv(jx,js,iq)-
     1                            Lv(jx,js,inx,isn,iq)*wni
                  enddo
               enddo
            enddo
            inode=(nstack-1)*3+3
            do inx=1,3
               wni=w(inx,inode)
               do js=1,4
                  do jx=1,3
                     vd(jx,js,iq)=vd(jx,js,iq)-
     1                            Ld(jx,js,inx,isn,iq)*wni
                  enddo
               enddo
            enddo
         enddo
      enddo
!--edge boundary contributions
      ve=0d0
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               do inx=1,3
                  wni=w(inx,inode)
                  do js=1,2
                     do jlv=1,3
                        do jx=1,3
                           ve(jx,jlv,js,il)=ve(jx,jlv,js,il)-
     1                               Le(inx,lvn,isn,jx,jlv,js,il)*wni
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      f=0d0
      call cgvasmbl(f,vq,vv,vd,ve)! f=-(vq+vv+vd+ve)
!
      return
      end subroutine cgvmxmul
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine cgvasmbl(vn,vq,vv,vd,ve)
!
!  this computes vn=-(vq+vv+vd+ve)
!  this is for VECTOR FIELDS
!
      implicit none
      real(8) vn(3,3*NSM)
      real(8) vq(3,3,4,NSM),vv(3,4,NSM),vd(3,4,NSM),ve(3,3,2,NLM)
      integer iq,isn,nstack,lvn,inode,inx,il
!
!--element contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               do inx=1,3
                  vn(inx,inode)=vn(inx,inode)-vq(inx,lvn,isn,iq)
               enddo
            enddo
         enddo
      enddo
!--ventral/dorsal contributions
      do iq=1,nq
         do isn=1,4
            nstack=isoq(isn,iq)
            inode=(nstack-1)*3+1
            do inx=1,3
               vn(inx,inode)=vn(inx,inode)-vv(inx,isn,iq)
            enddo
            inode=(nstack-1)*3+3
            do inx=1,3
               vn(inx,inode)=vn(inx,inode)-vd(inx,isn,iq)
            enddo
         enddo
      enddo
c--edge contributions
      do il=1,nl
         do isn=1,2
            nstack=isol(isn,il)
            do lvn=1,3
               inode=(nstack-1)*3+lvn
               do inx=1,3
                  vn(inx,inode)=vn(inx,inode)-ve(inx,lvn,isn,il)
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine cgvasmbl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine cgvcndmx(C,Q,Lv,Ld,Le)
!
! compute the conditoning matrix C used by the conjugate
! gradient solver for VECTOR FIELDS diffusion
!
!   C=diag(Q+Lv+Ld+Le)^-1 
!     where diag really means (3x3) diagonal matrices
!
      implicit none
!
      real(8) C(3,3,3*NSM)
      real(8) Q(3,3,4,3,3,4,NSM) 
      real(8) Lv(3,4,3,4,NSM), Ld(3,4,3,4,NSM), Le(3,3,2,3,3,2,NLM)
!
      real(8) cofac(3,3),ddet
      integer iq,isn,lvn,istack,inode,ix,jx,il
!
! initialize conditioning matrix
      C=0d0
!
!--interior contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               inode=(istack-1)*3+lvn
               do ix=1,3
                  do jx=1,3
                     C(jx,ix,inode)=C(jx,ix,inode)+
     1                              Q(jx,lvn,isn,ix,lvn,isn,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--boundary contributions
      do iq=1,nq
         do isn=1,4
            istack=isoq(isn,iq)
!--ventral nodes
            inode=(istack-1)*3+1
            do ix=1,3
               do jx=1,3
                  C(jx,ix,inode)=C(jx,ix,inode)+Lv(jx,isn,ix,isn,iq)
               enddo
            enddo
!--dorsal nodes
            inode=(istack-1)*3+3
            do ix=1,3
               do jx=1,3
                  C(jx,ix,inode)=C(jx,ix,inode)+Ld(jx,isn,ix,isn,iq)
               enddo
            enddo
         enddo
      enddo
!--edge nodes
      do il=1,nl
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               inode=(istack-1)*3+lvn
               do ix=1,3
                  do jx=1,3
                     C(jx,ix,inode)=C(jx,ix,inode)+
     1                              Le(jx,lvn,isn,ix,lvn,isn,il)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--now invert
!
      do inode=1,3*ns
!--cofactor matrix
         cofac(1,1)=+C(2,2,inode)*C(3,3,inode)-C(3,2,inode)*C(2,3,inode)
         cofac(2,1)=-C(1,2,inode)*C(3,3,inode)+C(3,2,inode)*C(1,3,inode)
         cofac(3,1)=+C(1,2,inode)*C(2,3,inode)-C(2,2,inode)*C(1,3,inode)
         cofac(1,2)=-C(2,1,inode)*C(3,3,inode)+C(3,1,inode)*C(2,3,inode)
         cofac(2,2)=+C(1,1,inode)*C(3,3,inode)-C(3,1,inode)*C(1,3,inode)
         cofac(3,2)=-C(1,1,inode)*C(2,3,inode)+C(2,1,inode)*C(1,3,inode)
         cofac(1,3)=+C(2,1,inode)*C(3,2,inode)-C(3,1,inode)*C(2,2,inode)
         cofac(2,3)=-C(1,1,inode)*C(3,2,inode)+C(3,1,inode)*C(1,2,inode)
         cofac(3,3)=+C(1,1,inode)*C(2,2,inode)-C(2,1,inode)*C(1,2,inode)
!--determinant 
         ddet=1d0/(C(1,1,inode)*cofac(1,1)+C(2,1,inode)*cofac(2,1)+
     1                                     C(3,1,inode)*cofac(3,1))
!--invert conditioning matrix
         do jx=1,3
            do ix=1,3
               C(ix,jx,inode)=ddet*cofac(jx,ix)
            enddo
         enddo
      enddo
!
      return
      end subroutine cgvcndmx
!
      END MODULE MXLIBSW!this ends the list of contained subprograms
