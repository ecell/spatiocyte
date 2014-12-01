      program rifgmv
!
! This program reads a sw dump and generates gmv source file
!
      use iolibsw
!
      integer nmat(3,NSM)!node materials
      real(8) sc
      parameter(sc=1d4)!scaling
      character ipfile*66
!
!--get name of source file from the command line and read it
!
      call getarg(1,ipfile)
      open(unit=11,file=ipfile,status='old')
      call iordfile(num,11)
      close(11)
      print *,'number of stacks ns=',ns
!
!--total number of nodes 3xstacks
      nn=ns*3
!
      open (unit=12,file='rif.gmv')
      write (12,'(A14)')'gmvinput ascii'
!
      write (12,*) ' '
      write (12,*) 'nodev ',nn
!
!--write node positions
      do is=1,ns
         do j=1,3
            write(12,*) hvec(1,j,is)*sc,hvec(2,j,is)*sc,hvec(3,j,is)*sc
         enddo
      enddo 
      write (12,*) ' '
!
! cells are defined by their face topologies
!
      write (12,*) 'cells ',nq
      do iq=1,nq
         write (12,*) ' '
         write(12,*) 'general 6'
! cells have two 4-node faces and four 6-node faces
         write(12,*) '4 4 6 6 6 6'
         is1=isoq(1,iq)
         is2=isoq(2,iq)
         is3=isoq(3,iq)
         is4=isoq(4,iq)
         n11=(is1-1)*3+1;n21=(is1-1)*3+2;n31=(is1-1)*3+3
         n12=(is2-1)*3+1;n22=(is2-1)*3+2;n32=(is2-1)*3+3
         n13=(is3-1)*3+1;n23=(is3-1)*3+2;n33=(is3-1)*3+3
         n14=(is4-1)*3+1;n24=(is4-1)*3+2;n34=(is4-1)*3+3
         write(12,*) n11,n14,n13,n12, !ventral face
     1               n31,n32,n33,n34, !dorsal face
     2               n12,n22,n32,n31,n21,n11,
     3               n13,n23,n33,n32,n22,n12,
     2               n14,n24,n34,n33,n23,n13,
     2               n11,n21,n31,n34,n24,n14
      enddo
!
! node materials
! interior nodes have material=1, surface nodes have material=2
! by default all level=1 and 3 are surface
!
      do is=1,ns
         nmat(1,is)=2
         nmat(2,is)=1
         nmat(3,is)=2
      enddo
! now scan for those middle nodes that are on the edge
      do il=1,nl
         is=isol(1,il)
         nmat(2,is)=2
      enddo
      write (12,*) ' '
      write (12,*) 'material 11 1' ! materials(11), node centered(1)
      write (12,*) 'interior_node' !1
      write (12,*) 'surface_node'  !2
      write (12,*) 'ventral'       !3
      write (12,*) 'dorsal'        !4
      write (12,*) 'edge_1'        !5
      write (12,*) 'edge_2'        !6
      write (12,*) 'edge_3'        !7
      write (12,*) 'edge_4'        !8
      write (12,*) 'edge_5'        !9
      write (12,*) 'edge_6'       !10
      write (12,*) 'substrate'    !11
      write(12,*) ((nmat(jl,js),jl=1,3),js=1,ns)
      write (12,*) ' '
!
!--the default node velocity is the network velocity
!
      write (12,*) 'velocity 1' !node centered(1)
      write(12,*) ((hvec(7,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ((hvec(8,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ((hvec(9,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ' '
!
!--now do variables
!
      write(12,*)'variable'
      write(12,*) ' '
!--pressure is always there
      write(12,*)'pressure 1'
      write(12,*) ((hvec(10,jl,js),jl=1,3),js=1,ns)
      write(12,*) ' '
!
!--look at svec for variables in use
      do k=1,NKS
         nskind=nint(dscp(k,idspk))
         if (nskind.eq.1.or.nskind.eq.2.or.nskind.eq.3) then
            write(12,*) names(k,idspk)//'1'
            write(12,*) ((svec(k,jl,js),jl=1,3),js=1,ns)
            write(12,*) ' '
         endif
      enddo
      write(12,*) 'endvars '
      write(12,*) ' '
!
!--other velocities written as vectors:
!
      write(12,*) 'vectors'
      write(12,*) ' '
      write(12,*) 'vv_relative 1 3 1' !node centered(1),xyz(3),named(1)
      write(12,*) 'vvrx'
      write(12,*) 'vvry'
      write(12,*) 'vvrz'
      write(12,*) ((hvec(4,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ((hvec(5,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ((hvec(6,jl,js)*sc,jl=1,3),js=1,ns)
      write(12,*) ' '
      write(12,*) 'vol_vel 1 3 1' !node centered(1),xyz(3),named(1)
      write(12,*) 'vvx'
      write(12,*) 'vvy'
      write(12,*) 'vvz'
      write(12,*) (((hvec(4,jl,js)+
     1               hvec(7,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*) (((hvec(5,jl,js)+
     1               hvec(8,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*) (((hvec(6,jl,js)+
     1               hvec(9,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*) ' '
      write(12,*) 'endvect'
      write(12,*) ' '
!
!--write tracers (these are useful to directly 
!                 examine variables at a node)
      write(12,*) 'tracers ',nn
      write(12,*) ((hvec(1,j,is)*sc,j=1,3),is=1,ns)
      write(12,*) ((hvec(2,j,is)*sc,j=1,3),is=1,ns)
      write(12,*) ((hvec(3,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vnx'
      write(12,*) ((hvec(7,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vny'
      write(12,*) ((hvec(8,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vnz'
      write(12,*) ((hvec(9,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vvrelx'
      write(12,*) ((hvec(4,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vvrely'
      write(12,*) ((hvec(5,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vvrelz'
      write(12,*) ((hvec(6,j,is)*sc,j=1,3),is=1,ns)
      write(12,*)'vvx'
      write(12,*) (((hvec(4,jl,js)+
     1               hvec(7,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*)'vvy'
      write(12,*) (((hvec(5,jl,js)+
     1               hvec(8,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*)'vvz'
      write(12,*) (((hvec(6,jl,js)+
     1               hvec(9,jl,js))*sc,jl=1,3),js=1,ns)
      write(12,*)'pressure'
      write(12,*) ((hvec(10,j,is),j=1,3),is=1,ns)
      do k=1,NKS
         nskind=nint(dscp(k,idspk))
         if (nskind.eq.1.or.nskind.eq.2.or.nskind.eq.3) then
            write(12,*) names(k,idspk)
            write(12,*) ((svec(k,jl,js),jl=1,3),js=1,ns)
         endif
      enddo
      write(12,*) ' '
      write(12,*) 'endtrace'
      write(12,*) ' '
      write(12,*) 'surface ',nl+2*nq
!--ventral surface
      do iq=1,nq
         is1=isoq(1,iq)
         is2=isoq(2,iq)
         is3=isoq(3,iq)
         is4=isoq(4,iq)
         write(12,*)4,(is1-1)*3+1,(is2-1)*3+1,(is3-1)*3+1,(is4-1)*3+1
      enddo
!--dorsal surface
      do iq=1,nq
         is1=isoq(1,iq)
         is2=isoq(2,iq)
         is3=isoq(3,iq)
         is4=isoq(4,iq)
         write(12,*)4,(is1-1)*3+3,(is2-1)*3+3,(is3-1)*3+3,(is4-1)*3+3
      enddo
!edge surfaces
      do il=1,nl
         is1=isol(1,il)
         is2=isol(2,il)
         n11=(is1-1)*3+1;n21=(is1-1)*3+2;n31=(is1-1)*3+3
         n12=(is2-1)*3+1;n22=(is2-1)*3+2;n32=(is2-1)*3+3
         write(12,*) 6,n12,n22,n32,n31,n21,n11
      enddo
      write(12,*) ' '
      write(12,*) 'surfmats '
! the ventral surface has material 3, the dorsal surface has material 4
! the edge surfaces have their material defined by ibol+4
      write(12,*) (3,iq=1,nq),(4,iq=1,nq),(ibol(il)+4,il=1,nl)
      write(12,*) ' '
      write(12,*) 'surfvars'
      write(12,*) ' '
      write(12,*) 'surf1'
      write(12,*) (vvec(1,iq),iq=1,nq),(dvec(1,iq),iq=1,nq),
     1                                 (evec(1,il),il=1,nl)
      write(12,*) 'surf2'
      write(12,*) (vvec(2,iq),iq=1,nq),(dvec(2,iq),iq=1,nq),
     1                                 (evec(2,il),il=1,nl)
      write(12,*) ' '
      write(12,*) 'endsvars'
      write(12,*) ' '
      write(12,*) 'polygons'
! all the polygons are are triangles with
! one of the vertices the center of gravity of 4 nodes
! and the other two vertices two of those nodes.
! ventral
      do iq=1,nq
         is1=isoq(1,iq)
         is2=isoq(2,iq)
         is3=isoq(3,iq)
         is4=isoq(4,iq)
         x1=hvec(1,1,is1)*sc;y1=hvec(2,1,is1)*sc;z1=hvec(3,1,is1)*sc
         x2=hvec(1,1,is2)*sc;y2=hvec(2,1,is2)*sc;z2=hvec(3,1,is2)*sc
         x3=hvec(1,1,is3)*sc;y3=hvec(2,1,is3)*sc;z3=hvec(3,1,is3)*sc
         x4=hvec(1,1,is4)*sc;y4=hvec(2,1,is4)*sc;z4=hvec(3,1,is4)*sc
         xcg=0.25*(x1+x2+x3+x4)
         ycg=0.25*(y1+y2+y3+y4)
         zcg=0.25*(z1+z2+z3+z4)
         write(12,*) 3,3,x1,x2,xcg,y1,y2,ycg,z1,z2,zcg
         write(12,*) 3,3,x2,x3,xcg,y2,y3,ycg,z2,z3,zcg
         write(12,*) 3,3,x3,x4,xcg,y3,y4,ycg,z3,z4,zcg
         write(12,*) 3,3,x4,x1,xcg,y4,y1,ycg,z4,z1,zcg
      enddo
! dorsal
      do iq=1,nq
         is1=isoq(1,iq)
         is2=isoq(2,iq)
         is3=isoq(3,iq)
         is4=isoq(4,iq)
         x1=hvec(1,3,is1)*sc;y1=hvec(2,3,is1)*sc;z1=hvec(3,3,is1)*sc
         x2=hvec(1,3,is2)*sc;y2=hvec(2,3,is2)*sc;z2=hvec(3,3,is2)*sc
         x3=hvec(1,3,is3)*sc;y3=hvec(2,3,is3)*sc;z3=hvec(3,3,is3)*sc
         x4=hvec(1,3,is4)*sc;y4=hvec(2,3,is4)*sc;z4=hvec(3,3,is4)*sc
         xcg=0.25*(x1+x2+x3+x4)
         ycg=0.25*(y1+y2+y3+y4)
         zcg=0.25*(z1+z2+z3+z4)
         write(12,*) 4,3,x1,x2,xcg,y1,y2,ycg,z1,z2,zcg
         write(12,*) 4,3,x2,x3,xcg,y2,y3,ycg,z2,z3,zcg
         write(12,*) 4,3,x3,x4,xcg,y3,y4,ycg,z3,z4,zcg
         write(12,*) 4,3,x4,x1,xcg,y4,y1,ycg,z4,z1,zcg
      enddo
! edges face has six nodes that are subdivided in lower 4 and upper 4
      do il=1,nl
         is1=isol(1,il)
         is2=isol(2,il)
         x11=hvec(1,1,is1)*sc;y11=hvec(2,1,is1)*sc;z11=hvec(3,1,is1)*sc
         x21=hvec(1,2,is1)*sc;y21=hvec(2,2,is1)*sc;z21=hvec(3,2,is1)*sc
         x31=hvec(1,3,is1)*sc;y31=hvec(2,3,is1)*sc;z31=hvec(3,3,is1)*sc
         x12=hvec(1,1,is2)*sc;y12=hvec(2,1,is2)*sc;z12=hvec(3,1,is2)*sc
         x22=hvec(1,2,is2)*sc;y22=hvec(2,2,is2)*sc;z22=hvec(3,2,is2)*sc
         x32=hvec(1,3,is2)*sc;y32=hvec(2,3,is2)*sc;z32=hvec(3,3,is2)*sc
         x1u=(3./8.)*x31+(3./4.)*x21+(-1./8.)*x11
         y1u=(3./8.)*y31+(3./4.)*y21+(-1./8.)*y11
         z1u=(3./8.)*z31+(3./4.)*z21+(-1./8.)*z11
         x2u=(3./8.)*x32+(3./4.)*x22+(-1./8.)*x12
         y2u=(3./8.)*y32+(3./4.)*y22+(-1./8.)*y12
         z2u=(3./8.)*z32+(3./4.)*z22+(-1./8.)*z12
         x1l=(3./8.)*x11+(3./4.)*x21+(-1./8.)*x31
         y1l=(3./8.)*y11+(3./4.)*y21+(-1./8.)*y31
         z1l=(3./8.)*z11+(3./4.)*z21+(-1./8.)*z31
         x2l=(3./8.)*x12+(3./4.)*x22+(-1./8.)*x32
         y2l=(3./8.)*y12+(3./4.)*y22+(-1./8.)*y32
         z2l=(3./8.)*z12+(3./4.)*z22+(-1./8.)*z32
         xcgl=(3./16.)*(x11+x12)+(3./8.)*(x21+x22)+(-1./16.)*(x31+x32)
         ycgl=(3./16.)*(y11+y12)+(3./8.)*(y21+y22)+(-1./16.)*(y31+y32)
         zcgl=(3./16.)*(z11+z12)+(3./8.)*(z21+z22)+(-1./16.)*(z31+z32)
         xcgu=(-1./16.)*(x11+x12)+(3./8.)*(x21+x22)+(3./16.)*(x31+x32)
         ycgu=(-1./16.)*(y11+y12)+(3./8.)*(y21+y22)+(3./16.)*(y31+y32)
         zcgu=(-1./16.)*(z11+z12)+(3./8.)*(z21+z22)+(3./16.)*(z31+z32)
         write(12,*) ibol(il)+4,3,x11,x12,xcgl,y11,y12,ycgl,z11,z12,zcgl
         write(12,*) ibol(il)+4,3,x11,x1l,xcgl,y11,y1l,ycgl,z11,z1l,zcgl
         write(12,*) ibol(il)+4,3,x12,x2l,xcgl,y12,y2l,ycgl,z12,z2l,zcgl
         write(12,*) ibol(il)+4,3,x21,x1l,xcgl,y21,y1l,ycgl,z21,z1l,zcgl
         write(12,*) ibol(il)+4,3,x22,x2l,xcgl,y22,y2l,ycgl,z22,z2l,zcgl
         write(12,*) ibol(il)+4,3,x21,x22,xcgl,y21,y22,ycgl,z21,z22,zcgl
         write(12,*) ibol(il)+4,3,x21,x22,xcgu,y21,y22,ycgu,z21,z22,zcgu
         write(12,*) ibol(il)+4,3,x21,x1u,xcgu,y21,y1u,ycgu,z21,z1u,zcgu
         write(12,*) ibol(il)+4,3,x22,x2u,xcgu,y22,y2u,ycgu,z22,z2u,zcgu
         write(12,*) ibol(il)+4,3,x31,x1u,xcgu,y31,y1u,ycgu,z31,z1u,zcgu
         write(12,*) ibol(il)+4,3,x32,x2u,xcgu,y32,y2u,ycgu,z32,z2u,zcgu
         write(12,*) ibol(il)+4,3,x31,x32,xcgu,y31,y32,ycgu,z31,z32,zcgu
      enddo
!--place big  square
      xmx=maxval(hvec(1,1:3,1:ns))*sc
      xmn=minval(hvec(1,1:3,1:ns))*sc
      ymx=maxval(hvec(2,1:3,1:ns))*sc
      ymn=minval(hvec(2,1:3,1:ns))*sc
!     xmx=2d2
!     xmn=-2d2
!     ymx=2d2
!     ymn=-2d2
      dy=ymx-ymn
      dx=xmx-xmn
      if (dx.gt.dy) then
         x1=xmn-0.5*dx
         x2=xmx+0.5*dx
         ymean=0.5*(ymx+ymn)
         y1=ymean-dx
         y2=ymean+dx
         write(12,*) 11,4,x1,x2,x2,x1,y1,y1,y2,y2,0d0,0d0,0d0,0d0
      else
         y1=ymn-0.5*dy
         y2=ymx+0.5*dy
         xmean=0.5*(xmx+xmn)
         x1=xmean-dy
         x2=xmean+dy
         write(12,*) 11,4,x1,x2,x2,x1,y1,y1,y2,y2,0d0,0d0,0d0,0d0
      endif
      write(12,*) 'endpoly'
      write(12,*) ' '
      write(12,*) 'probtime ',dscp(01,idfrm)
      write(12,*) ' '
      write (12,*) 'endgmv'
      close(12)
      stop
      end

