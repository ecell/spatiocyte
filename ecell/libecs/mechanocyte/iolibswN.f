      MODULE iolibsw
      USE iso_c_binding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module takes care of the I/O and of all the preliminary 
! mesh calculations that must be accomplished before using the FEM.
! It also contains a number of general purpose subroutines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!  Basic mesh structure is laminar with one element layer
!  Elements have 12 nodes, with 4 nodes placed between the
!  top and bottom quadrilaterals of the classic 8-node brick.
!  Nodes are referenced by their global stack number and 
!  and by their local level number (1,2, or 3).  Therefore
!  each stack encompasses three nodes and each element is
!  made up of four stacks.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer,parameter::NKN=12! #components in hvec
      integer,parameter::NKS=12! #components in svec
      integer,parameter::NKL=4!  #components surface vectors.
      integer,parameter::NKD=12! #components in descriptor-vectors.
      integer,parameter::NSM=2048!max storage for stacks/hexahedrons
      integer,parameter::NLM=512!max storage for edges 
      integer,parameter::NDM=64!max storage for descriptor vectors
!
! fixed parameters used for setting arithmetic roundoff control
      real(8),parameter::zero=0d0,tiny=1.0d-16,vtiny=1.0d-32
      real(8),parameter::one=1d0,big=1.0d+16,vbig=1.0d+32
      real(8),parameter::half=0.5d0,quarter=0.25d0,two=2d0 
!
!-- Three Gauss point quadrature weight parameters 5/9,8/9,5/9
      real(8), dimension(3), parameter :: wgp=(/0.55555555555556,
     1                                          0.88888888888889,
     2                                          0.55555555555556/)
!-- kronecker delta
      integer, dimension(0:3,0:3), 
     1    parameter::kd=reshape((/1,0,0,0,0,1,0,0,
     2                            0,0,1,0,0,0,0,1/),(/4,4/))
!-- index array
      integer, dimension(NSM), parameter::idxarr=(/(i,i=1,NSM)/)
!
!--read from the input file
      character(len=66)::mshtitle!(1st line of input file).
      character(len=66),dimension(NDM)::heads!=descriptor-headings.
      character(len=11),dimension(NKD,NDM)::names!descriptor-names.
      real(c_double),dimension(NKN,3,NSM),bind(C)::hvec=0!node hard vectors
      real(8),dimension(NKS,3,NSM)::svec=0!node soft vectors
      real(8),dimension(NKL,NSM)::dvec=0!dorsal surface vectors.
      real(8),dimension(NKL,NSM)::vvec=0!ventral surface vectors.
      real(8),dimension(NKL,NLM)::evec=0!edge vectors.
      real(c_double),dimension(NKD,NDM),bind(C)::dscp=0!descriptive vectors 
      integer(C_INT),dimension(4,NSM),bind(C)::isoq=0!index of ith stack of Qj
      integer,dimension(2,NLM)::isol=0!index of ith stack of Lj
      integer,dimension(NLM)::ibol=0!index of edge bdy of Lj
      integer(C_INT),bind(C)::ns=0!number of stacks (must be less than NSM)
      integer(C_INT),bind(C)::nq=0!number of hexahedrons (must be less than NSM)
      integer::nl=0!number of edges (must be less than NLM)
      integer::nd=0!number of descriptor vectors (must be less than NDM)
      real(c_double),bind(C)::time_=0!the current value of the time_
      real(c_double),bind(C)::tstp=0!the current value of the time_ step
!
! COMPUTED QUANTITIES:
!
!  Geometry:
      integer,dimension(NSM)::iloq=0!index of ith-L of Q_j
      integer,dimension(2,NLM)::ilol=0!index of neighbors of L_j
      integer,dimension(2,NSM)::ilos=0!index of edge neighbors stack j
      integer,dimension(NLM)::iqol=0!index of Q interior to L_j
      integer,dimension(NSM)::isos=0!index of stack interior to stack
      integer(C_INT),dimension(NSM),bind(C)::kqos=0!number of elements containing stack 
      integer(C_INT),dimension(NSM),bind(C)::iqos=0!starting point in lqos list for stack
      integer(C_INT),dimension(4*NSM),bind(C)::lqos=0!list of elements in which stacks 
                                      !belong:
                                      !lqos(iqos(is):iqos(is)+kqos(is)-1)
!
!  SHAPE FUNCTIONS
      real(8),dimension(3,4,3,4)::H=0! 12-node brick shape function
!--node level; node stack; gauss point level; gauss point stack 
!
      real(8),dimension(3,3,4,3,4)::dH=0! nablaH wrt intrisic space coords
!--coord; node level; node stack; gauss point level; gauss point stack 
!
      real(8),dimension(4,4)::S4=0! 4-node surface shape function
!--node stack, gauss point stack
!
      real(8),dimension(2,4,4)::dS4=0! nablaS4 wrt intrisic surface coords
!--xi,eta coord; node stack; gauss point stack
!
      real(8),dimension(3,2,3,2)::S6=0! 6-node surface shape function
!--node level; node stack; gauss point level; gauss point stack
!
      real(8),dimension(2,3,2,3,2)::dS6=0! nablaS6 wrt intrisic coords
!--xi, zeta coord; node level; node stack; gp level; gp stack
!
      real(8),dimension(3,3,4,3,4,NSM)::dHg=0!dH_i wrt to x_j @ gp_k 
!--x-y-z; node level; node stack; gauss point level; gauss point stack
!
      real(8),dimension(3,4,NSM)::det=0!Jacobian det at gp_i of Q_j
!--Gauss point level, Gauss point stack, quadrilateral
!
      real(8),dimension(0:3,4,NSM)::dAd=0!out dorsal surface area vector
      real(8),dimension(0:3,4,NSM)::dAv=0!out ventral surface area vector
!--norm,x-y-z components.; gauss point stack; quadrilateral
      real(8),dimension(0:3,3,2,NLM)::dAe=0!out edge surface area vector
!--norm,x-y-z components.; gp level; gauss point stack; edge element
!
      real(8),dimension(3,4,4,NSM)::dS4dg=0!, dS4 dorsal wrt to x_j
      real(8),dimension(3,4,4,NSM)::dS4vg=0!, dS4 ventral wrt to x_j
!--x-y-z component; node stack, Gauss point stack, quadrilateral
      real(8),dimension(3,3,2,3,2,NLM)::dS6g=0!, dS6 edge wrt to x_j
!--x-y-z; node level, stack; Gauss point level, stack, edge element
!
! der. of x,y,z wrt to xi and eta (or zeta) at GP on surfaces
      real(8),dimension(2,4,NSM)::dxSv=0,dySv=0,dzSv=0
      real(8),dimension(2,4,NSM)::dxSd=0,dySd=0,dzSd=0
! xi-eta, GP stack, quadrilateral
      real(8),dimension(2,3,2,NLM)::dxSe=0,dySe=0,dzSe=0
! xi-zeta, GP level, GP stack, edge element
!
      real(8),dimension(0:3,3,NSM)::snn=0 !unit outward normal at node
! norm,x,y,z,level,stack, if norm=0 interior node, or no curvature.
!
! der. of x,y,z wrt to xi eta and zeta at GP for each element
      real(8),dimension(3,3,4,NSM)::dxV=0,dyV=0,dzV=0
! xi-eta-zeta, GP level, GP stack, volume element
!
      real(8),dimension(0:3,3,NSM)::vnn=0 !normalized volume gradient
! norm,x,y,z,level,stack, 
!
      real(8),dimension(3,NSM)::voln=0 !volume associated with each node
!
      real(8),dimension(3,NSM)::arean=0 !area associated with each node
!
      real(8),dimension(0:2,NLM)::cnn=0 !contact line normal at isol(1,il)
! line length, nx, ny (nz=0 because in plane).
!
      real(8),dimension(3,3,NSM)::vtg1=0!1st unit volume tangent vector
      real(8),dimension(3,3,NSM)::vtg2=0!2nd unit volume tangent vector
!
      real(8),dimension(NLM)::ca=0!contact angle at node isol(1,il)
!
! idxxx = position of relevant parameters in descriptors dscp(:,idxxx)
      integer ::idmsh=0!index of 'MESH' descriptors
      integer(C_INT),bind(C)::idfrm=0!index of 'FRAME' descriptors
      integer ::idhpk=0!index of 'NODE HARD VECTOR PACK' descriptors
      integer ::idspk=0!index of 'NODE SOFT VECTOR PACK' descriptors
      integer ::idbpk=0!index of 'BOUNDARY VECTOR PACK' descriptors
      integer ::idnfx=0!index of 'VELOCITY CONSTRAINT'descriptors
      integer ::idgmo=0!index of 'GRID MOTION' descriptors
      integer(C_INT),bind(C)::idvfx=0!index of 'VELOCITY CONSTRAINTS'
      integer(C_INT),bind(C)::idvis=0!index of 'NETWORK VISCOSITY' descriptors
      integer(C_INT),bind(C)::idphi=0!index of 'PHASE RUB' descriptors
      integer(C_INT),bind(C)::iddrg=0!index of 'DRAG' descriptors
      integer(C_INT),bind(C)::idpsi=0!index of 'CONTRACT' descriptors
      integer(C_INT),bind(C)::idbfr=0!index of 'BODY FOR' descriptors
      integer(C_INT),bind(C)::idhyc=0!index of 'HYDRAULIC COND' descriptors
      integer ::idosm=0!index of 'OSMOTIC' descriptors
      integer(C_INT),bind(C)::idgam=0!index of 'SURFACE TEN' descriptors
      integer(C_INT),bind(C)::idtrc=0!index of 'RY TRAC' descriptors
      integer(C_INT),bind(C)::idsfr=0!index of 'SURFACE FOR' descriptors
      integer ::idexv=0!index of 'EDGE EXTERNAL VALUE' descriptors
      integer ::idexp=0!index of 'EDGE EXTERNAL PERMEABILITY' descript
      integer ::idexs=0!index of 'EDGE EXTERNAL SOURCE' descriptors
      integer ::idprm=0!index of 'EDGE EXTERNAL PERMEABILITY' descript
      integer ::idscm=0!index of 'SVEC CHEMISTRY' descriptors
      integer ::idsdc=0!index of 'SVEC DIFFUSION' descriptors
      integer ::iderg=0!index of 'ENERGY DISSIPATION' descriptors
!
! arrays for constitutive laws
!
!
!  chemistry: svec component, node level, node stack
      real(8),dimension(NKS,3,NSM)::sdot=0!d/dt(svec(i,*))from chemistry.
      real(8),dimension(NKS,3,NSM)::sdkr=0!d(sdot(i,*))/d(svec(i,*)).
!
! diffusion boundary conditions: 
!   first index -> corresponding svec field
!   second index -> element number (for ventral and dorsal) or
!                   edge number (for edge surface)
!--Neumann:
      real(8),dimension(NKS,NSM)::vxsrc=0!source at ventral surface
      real(8),dimension(NKS,NSM)::dxsrc=0!source at dorsal surface
      real(8),dimension(NKS,NLM)::exsrc=0!source at edge surface
      real(8),dimension(NKS,NLM)::cxsrc=0!source at contact line (surf dif)
!--Dirichlet:
      real(8),dimension(NKS,NSM)::vxval=0!ext value at ventral surface
      real(8),dimension(NKS,NSM)::dxval=0!ext value at dorsal surface
      real(8),dimension(NKS,NLM)::exval=0!ext value at edge surface
      real(8),dimension(NKS,NLM)::cxval=0!ext value at contact line (surf dif)
      real(8),dimension(NKS,NSM)::vxprm=0!permeability at ventral surf
      real(8),dimension(NKS,NSM)::dxprm=0!permeability at dorsal surf
      real(8),dimension(NKS,NLM)::exprm=0!permeability at edge surf
      real(8),dimension(NKS,NLM)::cxprm=0!permeability at contact line (surf dif)
!
!  Momentum
      real(8),dimension(3,NSM)::vis=0!network shear viscosity
      real(8),dimension(3,NSM)::lam=0!network pure dilation viscosity
      real(8),dimension(3,NSM)::phi=0!network-solvent friction coefficent
      real(8),dimension(3,NSM)::psi=0!network contractility
      real(8),dimension(3,3,NSM)::bfr=0!network bodyforce vector
!--surface tension
      real(8),dimension(NSM)::gamv=0!ventral
      real(8),dimension(NSM)::gamd=0!dorsal
      real(8),dimension(NLM)::game=0!edge
!--fixed boundary velocities (xyz component, element or edge index)
      real(8),dimension(3,NSM)::vfixv=0!ventral
      real(8),dimension(3,NSM)::vfixd=0!dorsal
      real(8),dimension(3,NLM)::vfixe=0!edge
! external velocity components x,y,z at boundary (used if vfix ne 0)
      real(8),dimension(3,NSM)::vbv=0!ventral
      real(8),dimension(3,NSM)::vbd=0!dorsal
      real(8),dimension(3,NLM)::vbe=0!edge
!--surface friction coefficients (x/y/z/tangential/normal, element or edge)
      real(8),dimension(5,NSM)::drgv=0!ventral 
      real(8),dimension(5,NSM)::drgd=0!dorsal 
      real(8),dimension(5,NLM)::drge=0!edge 
! external boundary tractions
      real(8),dimension(3,NSM)::trcv=0!ventral
      real(8),dimension(3,NSM)::trcd=0!dorsal
      real(8),dimension(3,NLM)::trce=0!edge
!--surface forces
      real(8),dimension(NSM)::sfrv=0!ventral 
      real(8),dimension(NSM)::sfrd=0!dorsal 
      real(8),dimension(NLM)::sfre=0!edge 
      real(8),dimension(3,NSM)::sfrn=0!node-defined
!--solvent boundary permeability 
      real(8),dimension(NSM)::hycv=0!ventral
      real(8),dimension(NSM)::hycd=0!dorsal
      real(8),dimension(NLM)::hyce=0!edge
!--relative solvent osmotic activity
      real(8),dimension(NSM)::osmv=0!ventral
      real(8),dimension(NSM)::osmd=0!dorsal
      real(8),dimension(NLM)::osme=0!edge
!
! Arrays for grid motion
!
      integer,dimension(NSM)::irezfix=0! stack motion =0 allowed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      CONTAINS!the following are the subroutines of the io module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine iowrfile(ifrm,idiv)
!
! Archives to disk a simulation file
!
      character (len=11):: score= '***** *****'
      write(idiv,200)ifrm,mshtitle!write the frame# and title
      dscp(1,idfrm)=time_
      dscp(2,idfrm)=tstp
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nd,'CONTROL PARAMETERS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nd
         write(idiv,200)i,heads(i)!heading for the i-th discriptor
         write(idiv,210)(names(k,i),k=01,06)
         write(idiv,220)( dscp(k,i),k=01,06)
         write(idiv,210)(names(k,i),k=07,12)
         write(idiv,220)( dscp(k,i),k=07,12)
      enddo
!
!--write POINTERS of EDGE TOPOLOGY
!   each edge has two stacks of 3 nodes associated with it
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nl,'POINTERS OF EDGE TOPOLOGY '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nl,3
         write(idiv,230)(j,'=',isol(1,j),isol(2,j),ibol(j),j=i,i+2)
      enddo
!
!--write POINTERS of ELEMENT TOPOLOGY
!   each element has 4 stacks of 3 nodes associated with it
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nq,'POINTERS OF ELEMENT TOPOLOGY '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nq,2
         write(idiv,235)('*****',j,'=',(isoq(k,j),k=1,4),j=i,i+1)
      enddo
!
!--write DATA DEFINED at SURFACES 
!   each surface (edge, ventral,or dorsal) has 4 fields available 
!   for storage of data.
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nl,'EDGE SURFACE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nl
         write(idiv,240)i,'=',(evec(k,i),k=01,04)
      enddo
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nq,'VENTRAL SURFACE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nq
         write(idiv,240)i,'=',(vvec(k,i),k=01,04)
      enddo
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nq,'DORSAL SURFACE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nq
         write(idiv,240)i,'=',(dvec(k,i),k=01,04)
      enddo
!
c--write DATA DEFINED ON VERTEX NODES
!   each node is identified by its global stack number and
!   its level (1,2, or 3). At each node 12 fields are
!   "hard data" (e.g. xyz position) and 12 fields are
!   "soft data" (e.g. chemical concentration)
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)ns,'NODE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,ns
         do l=1,3
            write(idiv,240)i,'=',(hvec(k,l,i),k=01,04)
            write(idiv,250)      (hvec(k,l,i),k=05,08)
            write(idiv,250)      (hvec(k,l,i),k=09,12)
            write(idiv,240)i,'=',(svec(k,l,i),k=01,04)
            write(idiv,250)      (svec(k,l,i),k=05,08)
            write(idiv,250)      (svec(k,l,i),k=09,12)
         enddo
      enddo
      write(idiv,210)score,score,score,score,score,score
      call iobuffer(10,idiv)!flush the io buffer
!
!--FORMAT STATEMENTS
  200 format(1x,i4,1x,a66)!formate for io of title+headings
  210 format(6(1x,a11))!format for names
  220 format(6(1x,1pe11.4))!format for dscp values 
  230 format(3(1x,i4,a1,3(1x,i5)))!format;for L pointers
  235 format(2(1x,a5,i5,a1,4(1x,i5)))!format;for Q pointers
  240 format(i5,a1,2x,4(1x,1pe15.7))
  250 format(5x,1x,2x,4(1x,1pe15.7))
!
      return
      end subroutine iowrfile
!
      subroutine iobuffer(nflush,idiv)
      integer nflush,idiv,i
      real(8) rn
      do i=1,nflush
         rn=6*(i-1)
         write(idiv,30)rn+1,rn+2,rn+3,rn+4,rn+5,rn+6
      enddo
   30 format(6(1x,1pe11.4))
      return
      end subroutine iobuffer!(nflush,idiv)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine iordfile(ifrm,idiv)
!
!--This surbroutines reads a FEM shallow water simulation file
!
      character (len=66) label
      character (len=11):: score = '***** *****'
!
!     print *,'in iordfile'
      read(idiv,200)ifrm,mshtitle!read the frame# and title
      print *,'read title'
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nd,label! read ND 
      print *,'read label'
      read(idiv,210)score,score,score,score,score,score
      do i=1,nd!
         read(idiv,200)id,heads(i)!heading for the i-th descriptor
         read(idiv,210)(names(k,i),k=01,06)
         read(idiv,220)( dscp(k,i),k=01,06)
         read(idiv,210)(names(k,i),k=07,12)
         read(idiv,220)( dscp(k,i),k=07,12)
      enddo
      print *,'read DESCRIPTORS'
!--read POINTERS of EDGE TOPOLOGY
!   each edge has two stacks of 3 nodes associated with it
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nl,label!READ ne = number of EDGEs
      read(idiv,210)score,score,score,score,score,score
      do i=1,nl,3
         read(idiv,230)(isol(1,j),isol(2,j),ibol(j),j=i,i+2)
      enddo
      print *,'read POINTERS of EDGE TOPOLOGY'
!
!--read POINTERS of ELEMENT TOPOLOGY
!   each element has 4 stacks of 3 nodes associated with it
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nq,label!read nq=number of quadrilaterials
      read(idiv,210)score,score,score,score,score,score
      do i=1,nq,2
         read(idiv,235)((isoq(k,j),k=1,4),j=i,i+1)
      enddo
      print *,'read POINTERS of ELEMENT TOPOLOGY'
!
!--read DATA DEFINED at surfaces
!   each surface (edge, ventral,or dorsal) has 4 fields available 
!   for storage of data.
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nl,label!'EDGE SURFACE DATA VECTORS '
      read(idiv,210)score,score,score,score,score,score
      do i=1,nl
         read(idiv,250)(evec(k,i),k=01,04)
      enddo
      print *,'read EDGE SURFACE DATA VECTORS'
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nq,label!'VENTRAL SURFACE DATA VECTORS '
      read(idiv,210)score,score,score,score,score,score
      do i=1,nq
         read(idiv,250)(vvec(k,i),k=01,04)
      enddo
      print *,'read VENTRAL SURFACE DATA VECTORS'
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nq,label!'DORSAL SURFACE DATA VECTORS '
      read(idiv,210)score,score,score,score,score,score
      do i=1,nq
         read(idiv,250)(dvec(k,i),k=01,04)
      enddo
      print *,'read DORSAL SURFACE DATA VECTORS'
!
!--read DATA DEFINED ON NODES
!   each node is identified by its global stack number and
!   its level (1,2, or 3). At each node 12 fields are
!   "hard data" (e.g. xyz position) and 12 fields are
!   "soft data" (e.g. chemical concentration)
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)ns,label!number of stacks
      read(idiv,210)score,score,score,score,score,score
      do i=1,ns
!-- loop over the 3 levels of the stack
         do l=1,3
            read(idiv,250)(hvec(k,l,i),k=01,04)
            read(idiv,250)(hvec(k,l,i),k=05,08)
            read(idiv,250)(hvec(k,l,i),k=09,12)
            read(idiv,250)(svec(k,l,i),k=01,04)
            read(idiv,250)(svec(k,l,i),k=05,08)
            read(idiv,250)(svec(k,l,i),k=09,12)
         enddo
      enddo
!
      call iddscp !id the descriptors
      call topomsh !compute mesh topology dependent information
      call ioshapef!initialize shape functions
      call godriver(hvec(1:3,:,:))!computes geometry dependent stuff
!
      ns=nint(dscp(01,idmsh)) !number of stacks
      if (ns.gt.NSM) then
         print *,'too many stacks: ns, NSM',ns,NSM
         stop
      endif
      nq=nint(dscp(02,idmsh)) !number of elements
      if (nq.gt.NSM) then
         print *,'too many elements: nq, NSM',nq,NSM
         stop
      endif
      nl=nint(dscp(03,idmsh)) !number of edge surfaces
      if (nl.gt.NLM) then
         print *,'too many edges: nl, NLM',nl,NLM
         stop
      endif
      nd=nint(dscp(04,idmsh)) !number of descriptors
      if (nd.gt.NDM) then
         print *,'too many descriptors: nd, NDM',nd,NDM
         stop
      endif
      time_=dscp(01,idfrm)!current value of time_ 
      tstp=max(dscp(02,idfrm),tiny)!starting value of tstp ne.0
!
!--FORMAT STATEMENTS
  200 format(1x,i4,1x,a66)!formate for io of title+headings
  210 format(6(1x,a11))!format for names
  220 format(6(1x,1pe11.4))!format for dscp values 
  230 format(3(1x,4x,1x,3(1x,i5)))!format;for L pointers
  235 format(2(1x,5x,5x,1x,4(1x,i5)))!format;for Q pointers
  250 format(5x,1x,2x,4(1x,1pe15.7))!evec, hvec, and svec data
!
      return
      end subroutine iordfile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine iddscp
!
! this subroutine  figures out which descriptors are present and where
      call iogetidx(idmsh,'MESH')
      call iogetidx(idfrm,'FRAME')
      call iogetidx(idhpk,'NODE HARD VECTOR PACK')
      call iogetidx(idspk,'NODE SOFT VECTOR PACK')
      call iogetidx(idbpk,'BOUNDARY VECTOR PACK')
      call iogetidx(idnfx,'VELOCITY CONSTRAINT')
      call iogetidx(idgmo,'GRID MOTION')
      call iogetidx(idvfx,'VELOCITY CONSTRAINTS')
      call iogetidx(idvis,'NETWORK VISCOSITY')
      call iogetidx(idpsi,'CONTRACT')
      call iogetidx(idphi,'PHASE RUB')
      call iogetidx(idtrc,'BOUNDARY TRACTION')
      call iogetidx(idgam,'SURFACE TENSION')
      call iogetidx(idsvs,'SURFACE VISCOSITY')
      call iogetidx(iddrg,'SURFACE DRAG')
      call iogetidx(idsfr,'SURFACE FOR')
      call iogetidx(idhyc,'HYDRAULIC COND')
      call iogetidx(idosm,'OSMOTIC')
      call iogetidx(iderg,'ENERGY DISSIPATION')
      call iogetidx(idsdc,'SVEC DIFFUSION')
      call iogetidx(idscm,'SVEC CHEMISTRY')
      call iogetidx(idexv,'EDGE EXTERNAL VALUE')
      call iogetidx(idexs,'EDGE EXTERNAL SOURCE')
      call iogetidx(idexp,'EDGE EXTERNAL PERMEABILITY')
      call iogetidx(idldc,'LVEC DIFFUSION')
      call iogetidx(idxlv,'LVEC EXTERNAL')
      call iogetidx(idlvs,'LVEC EDGE SOURCE')
      call iogetidx(idlvp,'LVEC EDGE PERM')
      call iogetidx(idlcm,'LVEC CHEMISTRY')
      call iogetidx(idgfr,'MEMBRANE NETWORK INTERACTION')
!
      if(idmsh.le.0)then;write(*,*)'IDMSH.LE.ZERO!';stop;endif!FATAL ER
      if(idfrm.le.0)then;write(*,*)'IDFRM.LE.ZERO!';stop;endif!FATAL ER
      end subroutine iddscp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine iogetidx(idx,string)
      character(len=*) string
! Search a list of character strings; heads(j) for j=1,2,.. nh,
! for the first occurance of a pattern. If heads(k) is the result
! set idx=k. If no match is found set idx=0.
      idx=0
      do i=1,NDM
        if (INDEX(heads(i),string).gt.0) idx=i
      enddo
      return
      end subroutine iogetidx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine topomsh
! this subroutine populates arrays determined by the topology of the mesh
!
      call golostop!init ilos(*,*)=edges to which stack belongs
      call gololtop!init ilol(*,*)=neighboring edges of an edge
      call goloqoltop!init iloq(*) and iqol(*) (element-edge arrays)
      call gosostop!init isos(*)
      call goqostop!init kqos(*),iqos(*),lqos(*)
!
      end subroutine topomsh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine golostop
!
! this subroutine compute the edges boundaries (if any) that
! a stack belongs to. (otherwise, it is 0)
!
      implicit none
!
      integer il,is1,is2
!
      do il=1,nl
         is1=isol(1,il)
         is2=isol(2,il)
         ilos(2,is1)=il
         ilos(1,is2)=il
      enddo
      return
!
      end subroutine golostop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gololtop
!
! this subroutine populates the array ilol
!  ilol(1,il)=index of edge-elemt clockwise of il-th edge elemt
!  ilol(2,il)=index of edge-elemt anti-clockwise of il-th edge elemt
!
      implicit none
!
      integer il,jl,is1,is2
!
      ilol=0
      do il=1,nl
         is1=isol(1,il)
         is2=isol(2,il)
         ilol(1,il)=0
         ilol(2,il)=0
         do jl=nl,1,-1
            if(is1.eq.isol(2,jl))ilol(1,il)=jl
            if(is2.eq.isol(1,jl))ilol(2,il)=jl
         enddo
      enddo
      return
      end subroutine gololtop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine goloqoltop
!
! This subroutine populates the arrays iqol and iloq: for each edge,
! it determines the index of the interior element, for each element
! it determines the edge (if any)
! It also reorders isoq for that element so that 
! isoq(1,iqol(il)) =  isol(1,il) 
! isoq(2,iqol(il)) =  isol(2,il)
!
      implicit none
!
      integer il,iq,ism,isp,is1,is2,is3,is4,ispp
      iqol=0
      do il=1,nl
         ism=isol(1,il)
         isp=isol(2,il)
         ispp=isol(2,ilol(2,il))
! for each edge, determines the two stacks, and check if
! the pair matches with any element
         do iq=1,nq
            is1=isoq(1,iq)
            is2=isoq(2,iq)
            is3=isoq(3,iq)
            is4=isoq(4,iq)
            if (is1.eq.ism.and.is2.eq.isp) then
               iqol(il)=iq
               iloq(iq)=il
               if (is3.eq.ispp) then
                  print *,'WARNING: free stack at edge'
                  print *,'Please run dcorn'
                  stop
               endif
               exit
            elseif(is2.eq.ism.and.is3.eq.isp) then
               iqol(il)=iq
               iloq(iq)=il
               if (is4.eq.ispp) then
                  print *,'WARNING: free stack at edge'
                  print *,'Please run dcorn'
                  stop
               endif
               isoq(1,iq)=is2
               isoq(2,iq)=is3
               isoq(3,iq)=is4
               isoq(4,iq)=is1
               exit
            elseif(is3.eq.ism.and.is4.eq.isp) then
               iqol(il)=iq
               iloq(iq)=il
               if (is1.eq.ispp) then
                  print *,'WARNING: free stack at edge'
                  print *,'Please run dcorn'
                  stop
               endif
               isoq(1,iq)=is3
               isoq(2,iq)=is4
               isoq(3,iq)=is1
               isoq(4,iq)=is2
               exit
            elseif(is4.eq.ism.and.is1.eq.isp) then
               iqol(il)=iq
               iloq(iq)=il
               if (is2.eq.ispp) then
                  print *,'WARNING: free stack at edge'
                  print *,'Please run dcorn'
                  stop
               endif
               isoq(1,iq)=is4
               isoq(2,iq)=is1
               isoq(3,iq)=is2
               isoq(4,iq)=is3
               exit
            endif
         enddo
      enddo
!
      return
      end subroutine goloqoltop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosostop
!
! this subroutine populates the array isos
! isos(istack)=interior stack to a boundary stack 
!
      implicit none
!
      integer il,is,is1,is2,iq,i,istack,isp,ism
!
      isos=0
      do il=1,nl
         iq=iqol(il)
         istack=isoq(1,iq) !should be the same as isol(1,il)
         isos(istack)=isoq(4,iq)
      enddo
!
      return
      end subroutine gosostop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine goqostop
!
! this subroutine computes the "qos" arrays
!  kqos= # of elements to which stack s belongs
!  iqos= starting index in lqos
!  lqos= list of elements
!  i.e. the elements to which stack is belongs are given by:
!                         lqos(iqos(is):iqos(is)+kqos(is)-1)
!
      implicit none
!
      integer istack,iq,isq,k,i!,z,zz
!
      k=1
      do istack=1,ns !loop over stacks
         iqos(istack)=k!starting point in the list for istack
         do iq=1,nq  !loop over elements
            do isq=1,4
               if (isoq(isq,iq).eq.istack) then
                  lqos(k)=iq
                  kqos(istack)=kqos(istack)+1 !advance counter
                  k=k+1 !advance list counter
                  exit !move to the next quadrilateral
               endif
            enddo
         enddo
      enddo
      !do i=1,4
      !z=iqos(isoq(i,0))
      !zz=kqos(z)
     
      !print*,'sngsklgnskfbsjkdfsdklfnsldkf',zz
      !enddo
!     
      return
      end subroutine goqostop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine godriver(nxyz)
!
! This subroutine computes mesh geometry-dependent quantities
! assuming node positions nxyz and assuming that the topology
! is as defined by isoq, isol, etc.
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
! compute element-specific gradient of shape functions and associated stuff
      call goisomap(nxyz)
! compute surface-specific shape functions and ancillary data
      call gobdisom(nxyz)
! compute outward normal at each surface node
      call gosnn(nxyz)
! compute volume gradient at each node
      call govnn(nxyz)
! compute contact line normals
      call gocnn(nxyz)
! compute orthogonal pairs of (volume) tangent vectors for each surface node
      call govtg(nxyz)
! compute contact angles
      call goca(nxyz)
!
      return
      end subroutine godriver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Begin shape function routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine ioshapef
!
! This subroutine computes the values of shape functions at Gauss points
!
!  Element shape function at Gauss points
!     H(3,4,3,4): node (level, stack); gp (level, stack)
!  Element shape function gradient at GP (wrt to intrisic coord xi,eta,zeta)    
!     dH(3,3,4,3,4): coord.xi-eta-zeta; node (level,stack); gp (level,stack)
!  Ventral/Dorsal quadrilateral surface shape function at Gauss points
!     S4(4,4): node stack, surface Gauss point
!  Gradient of ventral/dorsal surface shape function at Gauss points
!     dS4(2,4,4) !node stack; comp.; surface Gauss point
!  Edge surface shape function at Gauss points
!     S6(3,2,3,2) !node level, stack, surface GP level, stack
!  Gradient of edge surface shape function at Gauss points
!     dS6(2,3,2,3,2)!node level;stack;comp.;GP level;stack
!
      implicit none
! Element nodes intrisic coordinates
      integer N(3,3,4) !xi,eta,zeta coord; level in stack; stack#
! Element Gauss points intrisic coordinates
      real(8) GP(3,3,4)!xi,eta,zeta coord; level in stack; stack#
! Ventral/Dorsal surface nodes intrisic coordinates
      integer N4(2,4) !xi,eta coord; stack#
! Ventral/Dorsal surface Gauss points intrisic coordinates
      real(8) GP4(2,4)!xi,eta coord; stack#
! Edge surface nodes intrisic coordinates
      integer N6(2,3,2) !xi,zeta; level in stack, stack #
! Edge surface Gauss points intrisic coordinates
      real(8) GP6(2,3,2)!xi,zeta; level in stack, stack #
!
      integer isg,lvg,isn,lvn
!     real(8), parameter::zgp=1.0 
!     real(8), parameter::xgp=1.0
!     real(8), parameter::zgm=-1.0
!     real(8), parameter::xgm=-1.0
      real(8), parameter::zgp=0.77459666924 !sqrt(3/5)
      real(8), parameter::xgp=0.57735026919 !1/sqrt(3)
      real(8), parameter::zgm=-0.77459666924 
      real(8), parameter::xgm=-0.57735026919
      data N/-1,-1,-1,-1,-1,0,-1,-1,+1,
     1       +1,-1,-1,+1,-1,0,+1,-1,+1,
     2       +1,+1,-1,+1,+1,0,+1,+1,+1,
     3       -1,+1,-1,-1,+1,0,-1,+1,+1/
      data N4/-1,-1,+1,-1,+1,+1,-1,+1/
      data N6/-1,-1,-1,0,-1,+1,
     1        +1,-1,+1,0,+1,+1/
      data GP/xgm,xgm,zgm,xgm,xgm,0d0,xgm,xgm,zgp,
     1        xgp,xgm,zgm,xgp,xgm,0d0,xgp,xgm,zgp,
     2        xgp,xgp,zgm,xgp,xgp,0d0,xgp,xgp,zgp,
     3        xgm,xgp,zgm,xgm,xgp,0d0,xgm,xgp,zgp/
      data GP4/xgm,xgm,xgp,xgm,xgp,xgp,xgm,xgp/
      data GP6/xgm,zgm,xgm,0,xgm,zgp,
     1         xgp,zgm,xgp,0,xgp,zgp/
      real(8) Hgp,dH1,dH2,dH3
!
!--shape functions evaluated at GPs of quadrilateral surface elements
!   loop by Gauss point stack, contributing node stack
      do isg=1,4
         do isn=1,4
!--working in the plane zeta=1
            call getH(N4(1,isn),N4(2,isn),1,GP4(1,isg),GP4(2,isg),1d0,
     2                Hgp,dH1,dH2,dH3)
            S4(isn,isg)=Hgp
            dS4(1,isn,isg)=dH1
            dS4(2,isn,isg)=dH2
         enddo
      enddo
!--shape functions evaluated at GPs of 6-node polygon surface elements
!   loop by Gauss point stack, GP level, contributing node stack
      do isg=1,2
         do lvg=1,3
            do isn=1,2
               do lvn=1,3
!--working in the plane eta=1
                  call getH(N6(1,lvn,isn),1,N6(2,lvn,isn),
     1                      GP6(1,lvg,isg),1d0,GP6(2,lvg,isg),
     2                      Hgp,dH1,dH2,dH3)
                  S6(lvn,isn,lvg,isg)=Hgp
                  dS6(1,lvn,isn,lvg,isg)=dH1
                  dS6(2,lvn,isn,lvg,isg)=dH3 !3 because 2 is eta deriv.
               enddo
            enddo
         enddo
      enddo
!--shape functions evaluated at GPs of finite element volumes
!   loop by Gauss point stack, GP level, contributing node stack
      do isg=1,4
         do lvg=1,3
            do isn=1,4
               do lvn=1,3
                  call getH(N(1,lvn,isn),N(2,lvn,isn),N(3,lvn,isn),
     1                      GP(1,lvg,isg),GP(2,lvg,isg),GP(3,lvg,isg),
     2                      Hgp,dH1,dH2,dH3)
                  H(lvn,isn,lvg,isg)=Hgp
                  dH(1,lvn,isn,lvg,isg)=dH1
                  dH(2,lvn,isn,lvg,isg)=dH2
                  dH(3,lvn,isn,lvg,isg)=dH3
!                 print *,'lvn,isn,lvg,isg',lvn,isn,lvg,isg
!                 print *,'Hgp,dH1,dH2,dH3',Hgp,dH1,dH2,dH3
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine ioshapef
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getH(nxi,neta,nzeta,xi,eta,zeta,
     1                Hshape,dHdxi,dHdeta,dHdzeta)
!
! This subroutine computes the value of the shape function generated by
! a node at arbitrary intrinsic coordinates
!
      implicit none
!
      integer nxi,neta,nzeta !node intrinsic coordinates
      real(8) xi,eta,zeta    !intrinsic coordinates where H needed
      real(8) Hshape,dHdxi,dHdeta,dHdzeta !shape function 
!
      real(8) xin,etan,zetan
!
      xin=dble(nxi)
      etan=dble(neta)
      zetan=dble(nzeta)
      if (nzeta.ne.0) then !top or bottom nodes
         Hshape=0.125d0*(1d0+xin*xi)*(1d0+etan*eta)*zeta*(zeta+zetan)
         dHdxi=0.125d0*xin*(1d0+etan*eta)*zeta*(zeta+zetan)
         dHdeta=0.125d0*(1d0+xin*xi)*etan*zeta*(zeta+zetan)
         dHdzeta=0.125d0*(1d0+xin*xi)*(1d0+etan*eta)*(2d0*zeta+zetan)
      else !middle nodes
         Hshape=0.25d0*(1d0+xin*xi)*(1d0+etan*eta)*(1d0-zeta*zeta)
         dHdxi=0.25d0*xin*(1d0+etan*eta)*(1d0-zeta*zeta)
         dHdeta=0.25d0*(1d0+xin*xi)*etan*(1d0-zeta*zeta)
         dHdzeta=-0.5d0*(1d0+xin*xi)*(1d0+etan*eta)*zeta
      endif
!
      return
      end subroutine getH
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getS4(nxi,neta,xi,eta,
     1                 Sshape,dSdxi,dSdeta,d2Sdxi2,d2Sdeta2,d2Sdxideta)
!
! This subroutine computes the value of the 4-node surface shape 
! function generated by a node at arbitrary intrinsic coordinates
!
      implicit none
!
      integer nxi,neta !node intrinsic coordinates
      real(8) xi,eta   !intrinsic coordinates where H needed
      real(8) Sshape,dSdxi,dSdeta !shape function 
      real(8) d2Sdxi2,d2Sdeta2,d2Sdxideta
!
      real(8) xin,etan
!
      xin=dble(nxi)
      etan=dble(neta)
      Sshape=0.25d0*(1d0+xin*xi)*(1d0+etan*eta)
      dSdxi=0.25d0*xin*(1d0+etan*eta)
      dSdeta=0.25d0*(1d0+xin*xi)*etan
      d2Sdxi2=0d0
      d2Sdeta2=0d0
      d2Sdxideta=0.25*xin*etan
!
      return
      end subroutine getS4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getS6(nxi,nzeta,xi,zeta,
     1              Sshape,dSdxi,dSdzeta,d2Sdxi2,d2Sdzeta2,d2Sdxidzeta)
!
! This subroutine computes the value of the 6-node surface shape 
! function generated by a node at arbitrary intrinsic coordinates
!
      implicit none
!
      integer nxi,nzeta !node intrinsic coordinates
      real(8) xi,zeta    !intrinsic coordinates where H needed
      real(8) Sshape,dSdxi,dSdzeta !shape function 
      real(8) d2Sdxi2,d2Sdzeta2,d2Sdxidzeta
!
      real(8) xin,zetan
!
      xin=dble(nxi)
      zetan=dble(nzeta)
      if (nzeta.ne.0) then !top or bottom nodes
         Sshape=0.25d0*(1d0+xin*xi)*zeta*(zeta+zetan)
         dSdxi=0.25d0*xin*zeta*(zeta+zetan)
         dSdzeta=0.25d0*(1d0+xin*xi)*(2d0*zeta+zetan)
         d2Sdxi2=0d0
         d2Sdzeta2=0.5d0*(1d0+xin*xi)
         d2Sdxidzeta=0.25d0*xin*(2d0*zeta+zetan)
      else !middle nodes
         Sshape=0.5d0*(1d0+xin*xi)*(1d0-zeta*zeta)
         dSdxi=0.5d0*xin*(1d0-zeta*zeta)
         dSdzeta=-(1d0+xin*xi)*zeta
         d2Sdxi2=0d0
         d2Sdzeta2=-(1d0+xin*xi)
         d2Sdxidzeta=-xin*zeta
      endif
!
      return
      end subroutine getS6
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getL2(nxi,xi,Lshape,dLdxi,d2Ldxi2)
!
! This subroutine computes the value of the 2-node line shape 
! function generated by a node at arbitrary intrinsic coordinates
!
      implicit none
!
      integer nxi !node intrinsic coordinates
      real(8) xi  !intrinsic coordinates where H needed
      real(8) Lshape,dLdxi,d2Ldxi2 !shape function 
!
      real(8) xin
!
      xin=dble(nxi)
      Lshape=0.5d0*(1d0+xin*xi)
      dLdxi=0.5d0*xin
      d2Ldxi2=0d0
!
      return
      end subroutine getL2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getL3(nzeta,zeta,Lshape,dLdzeta,d2Ldzeta2)
!
! This subroutine computes the value of the 3-node line shape 
! function generated by a node at arbitrary intrinsic coordinates
!
      implicit none
!
      integer nzeta !node intrinsic coordinates
      real(8) zeta    !intrinsic coordinates where H needed
      real(8) Lshape,dLdzeta,d2Ldzeta2 !shape function 
!
      real(8) zetan
!
      zetan=dble(nzeta)
      if (nzeta.ne.0) then !top or bottom nodes
         Lshape=0.5d0*zeta*(zeta+zetan)
         dLdzeta=0.5d0*(2d0*zeta+zetan)
         d2Ldzeta2=1d0
      else !middle nodes
         Lshape=1d0-zeta*zeta
         dLdzeta=-2d0*zeta
         d2Ldzeta2=-2d0
      endif
!
      return
      end subroutine getL3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gets(snodes,xi,eta,zeta,s,dsdxi,dsdeta,dsdzeta)
!
! Given values of scalar s defined at element nodes this subroutine 
! computes s and its intrinsic derivatives at arbitrary intrinsic 
! coordinates xi,eta,zeta
!
      implicit none
!
      real(8) snodes(3,4)  !s at node (level,stack)
      real(8) xi,eta,zeta  !input intrinsic coordinates
      real(8) s,dsdxi,dsdeta,dsdzeta  !s and deriv. at coord.
!
      integer, dimension(4), parameter :: nxi=(/-1,+1,+1,-1/)
      integer, dimension(4), parameter :: neta=(/-1,-1,+1,+1/)
      integer, dimension(3), parameter :: nzeta=(/-1,0,+1/)
!
      real(8) sn
      real(8) Hshape,dHdxi,dHdeta,dHdzeta
      integer isn,lvn
!
      s=0d0
      dsdxi=0d0
      dsdeta=0d0
      dsdzeta=0d0
      do isn=1,4
         do lvn=1,3
            call getH(nxi(isn),neta(isn),nzeta(lvn),xi,eta,zeta,
     1                Hshape,dHdxi,dHdeta,dHdzeta)
            sn=snodes(lvn,isn)
            s=s+Hshape*sn
            dsdxi=dsdxi+dHdxi*sn
            dsdeta=dsdeta+dHdeta*sn
            dsdzeta=dsdzeta+dHdzeta*sn
         enddo
      enddo
!
      return
      end subroutine gets
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getv12(vnodes,c,v,dvdc)
!
! Given values of vector v defined at element nodes this subroutine 
! computes v and its intrinsic derivatives at arbitrary intrinsic 
! coordinates xi,eta,zeta
!
      implicit none
!
      real(8) vnodes(3,3,4)  !v at node (component,level,stack)
      real(8) c(3)  !input intrinsic coordinates
      real(8) v(3)  !interpolated vector
      real(8) dvdc(3,3)  ! vector derivatives wrt to intrinsic coord
                         ! (component, intrinsic coord)
!
      integer, dimension(4), parameter :: nxi=(/-1,+1,+1,-1/)
      integer, dimension(4), parameter :: neta=(/-1,-1,+1,+1/)
      integer, dimension(3), parameter :: nzeta=(/-1,0,+1/)
!
      real(8) Hshape,dHdxi,dHdeta,dHdzeta
      real(8) vn(3)
      integer isn,lvn
!
      v=0d0
      dvdc=0d0
      do isn=1,4
         do lvn=1,3
            call getH(nxi(isn),neta(isn),nzeta(lvn),c(1),c(2),c(3),
     1                Hshape,dHdxi,dHdeta,dHdzeta)
            vn=vnodes(:,lvn,isn)
            v=v+Hshape*vn
            dvdc(:,1)=dvdc(:,1)+dHdxi*vn
            dvdc(:,2)=dvdc(:,2)+dHdeta*vn
            dvdc(:,3)=dvdc(:,3)+dHdzeta*vn
         enddo
      enddo
!
      return
      end subroutine getv12
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getv4(vnodes,xi,eta,
     1                 v,dvdxi,dvdeta,d2vdxi2,d2vdeta2,d2vdxideta)
!
! Given values of vector v defined at element nodes this subroutine 
! computes v and its intrinsic derivatives at arbitrary intrinsic 
! coordinates xi,eta
!
      implicit none
!
      real(8) vnodes(3,4)  !v at nodes (component,stack)
      real(8) xi,eta  !input intrinsic coordinates
      real(8) v(3)  !interpolated vector
      real(8) dvdxi(3)  ! vector derivative wrt to xi
      real(8) dvdeta(3) ! vector derivative wrt to eta
      real(8) d2vdxi2(3)    ! 2nd derivative wrt to xi (=0)
      real(8) d2vdeta2(3)    ! 2nd derivative wrt to eta (=0)
      real(8) d2vdxideta(3)  ! 2nd cross-derivative (non zero)
!
      integer, dimension(4), parameter :: nxi=(/-1,+1,+1,-1/)
      integer, dimension(4), parameter :: neta=(/-1,-1,+1,+1/)
!
      real(8) Sshape,dSdxi,dSdeta,dSdzeta
      real(8) d2Sdxi2,d2Sdeta2,d2Sdxideta
      real(8) vn(3)
      integer isn
!
      v=0d0
      dvdxi=0d0
      dvdeta=0d0
      d2vdxi2=0d0
      d2vdeta2=0d0
      d2vdxideta=0d0
      do isn=1,4
         call getS4(nxi(isn),neta(isn),xi,eta,
     1              Sshape,dSdxi,dSdeta,d2Sdxi2,d2Sdeta2,d2Sdxideta)
         vn=vnodes(:,isn)
         v=v+Sshape*vn
         dvdxi=dvdxi+dSdxi*vn
         dvdeta=dvdeta+dSdeta*vn
         d2vdxideta=d2vdxideta+d2Sdxideta*vn
      enddo
!
      return
      end subroutine getv4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getv6(vnodes,xi,zeta,
     1                 v,dvdxi,dvdzeta,d2vdxi2,d2vdzeta2,d2vdxidzeta)
!
! Given values of vector v defined at element nodes this subroutine 
! computes v and its intrinsic derivatives at arbitrary intrinsic 
! coordinates xi,zeta
!
      implicit none
!
      real(8) vnodes(3,3,2)  !v at nodes (component,level,stack)
      real(8) xi,zeta  !input intrinsic coordinates
      real(8) v(3)  !interpolated vector
      real(8) dvdxi(3) ! vector derivative wrt to xi
      real(8) dvdzeta(3) ! vector derivative wrt to zeta
      real(8) d2vdxi2(3) ! 2nd derivative wrt to xi (=0)
      real(8) d2vdzeta2(3) ! 2nd derivative wrt to zeta
      real(8) d2vdxidzeta(3) ! 2nd cross-derivative 
!
      integer, dimension(4), parameter :: nxi=(/-1,+1,+1,-1/)
      integer, dimension(3), parameter :: nzeta=(/-1,0,+1/)
!
      real(8) Sshape,dSdxi,dSdzeta
      real(8) d2Sdxi2,d2Sdzeta2,d2Sdxidzeta
      real(8) vn(3)
      integer isn,lvn
!
      v=0d0
      dvdxi=0d0
      dvdzeta=0d0
      d2vdxi2=0d0
      d2vdzeta2=0d0
      d2vdxidzeta=0d0
      do isn=1,2
         do lvn=1,3
            call getS6(nxi(isn),nzeta(lvn),xi,zeta,
     1             Sshape,dSdxi,dSdzeta,d2Sdxi2,d2Sdzeta2,d2Sdxidzeta)
            vn=vnodes(:,lvn,isn)
            v=v+Sshape*vn
            dvdxi=dvdxi+dSdxi*vn
            dvdzeta=dvdzeta+dSdzeta*vn
            d2vdzeta2=d2vdzeta2+d2Sdzeta2*vn
            d2vdxidzeta=d2vdxidzeta+d2Sdxidzeta*vn
         enddo
      enddo
!
      return
      end subroutine getv6
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getv2(vnodes,xi,v,dvdxi,d2vdxi2)
!
! Given values of vector v defined at line element nodes this 
! subroutine computes v and its intrinsic derivatives at arbitrary 
! intrinsic coordinates xi
!
      implicit none
!
      real(8) vnodes(3,2)  !v at nodes (component,stack)
      real(8) xi !input intrinsic coordinates
      real(8) v(3)  !interpolated vector
      real(8) dvdxi(3)  ! vector derivative wrt to xi
      real(8) d2vdxi2(3)    ! 2nd derivative wrt to xi (=0)
!
      integer, dimension(2), parameter :: nxi=(/-1,+1/)
!
      real(8) Lshape,dLdxi,d2Ldxi2
      real(8) vn(3)
      integer isn
!
      v=0d0
      dvdxi=0d0
      d2vdxi2=0d0
      do isn=1,2
         call getL2(nxi(isn),xi,Lshape,dLdxi,d2Ldxi2)
         vn=vnodes(:,isn)
         v=v+Lshape*vn
         dvdxi=dvdxi+dLdxi*vn
      enddo
!
      return
      end subroutine getv2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getv3(vnodes,zeta,v,dvdzeta,d2vdzeta2)
!
! Given values of vector v defined at line element nodes this 
! subroutine computes v and its intrinsic derivatives at arbitrary 
! intrinsic coordinates zeta
!
      implicit none
!
      real(8) vnodes(3,3)  !v at nodes (component,level)
      real(8) zeta  !input intrinsic coordinate
      real(8) v(3)  !interpolated vector
      real(8) dvdzeta(3) ! vector derivative wrt to zeta
      real(8) d2vdzeta2(3) ! 2nd derivative wrt to zeta
!
      integer, dimension(3), parameter :: nzeta=(/-1,0,+1/)
!
      real(8) Lshape,dLdzeta,d2Ldzeta2
      real(8) vn(3)
      integer lvn
!
      v=0d0
      dvdzeta=0d0
      d2vdzeta2=0d0
      do lvn=1,3
         call getL3(nzeta(lvn),zeta,Lshape,dLdzeta,d2Ldzeta2)
         vn=vnodes(:,lvn)
         v=v+Lshape*vn
         dvdzeta=dvdzeta+dLdzeta*vn
         d2vdzeta2=d2vdzeta2+d2Ldzeta2*vn
      enddo
!
      return
      end subroutine getv3
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getcoord(xyz,xyzn,xi,eta,zeta,iflag)
!
! Given the node positions xyzn of an element, this routine
! compute the xi,eta,zeta corresponding to a supplied xyz positions.
! this uses Newton's method.
! 
!     implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,3,4)   ! xyz coord. of 12 element nodes
      real(8) xi,eta,zeta ! input=initial xi,eta,zeta guess
                          ! output=actual xi,eta,zeta
      integer iflag !=-1 far outside the box, 
                    !=0  no convergence in 20 steps
                    !=+1 close to box, good convergence
!
      real(8) c(3)
      real(8) dc(3)
      real(8) xyzi(3)
      real(8) dxyz(3)
      real(8) dxyzdc(3,3)
      real(8) dcdxyz(3,3)
      integer i
!      
      c(1)=xi
      c(2)=eta
      c(3)=zeta
      iflag=0
      do i=1,20
         call getv12(xyzn,c,xyzi,dxyzdc)
         dxyz=xyz-xyzi
         call Mat3Inv(dxyzdc,dcdxyz)
         dc=matmul(dcdxyz,dxyz)
!        print *,'c',c
         c=c+dc
!        print *,'dc',dc
         if (maxval(abs(c)).gt.2d0) then
!--too far out
            iflag=-1
            exit
         endif
         if (maxval(abs(dc)).lt.1d-9) then
!--we have converged
            iflag=+1
            exit
         endif
      enddo
      xi=c(1)
      eta=c(2)
      zeta=c(3)
!
!--if there was no convergence
      if (iflag.eq.0) then
          call getcoord2(xyz,xyzn,xi,eta,zeta,iflag)
          print *,'getcoord: warning getcoord2 was called!'
      endif
!
      return 
      end subroutine getcoord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getcoord2(xyz,xyzn,xi,eta,zeta,iflag)
!
! Given the node positions xyzn of an element, this routine
! compute the xi,eta,zeta corresponding to a supplied xyz positions.
! The method used is a variant of bisection in intrinsic space.
! 
!     implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,3,4)   ! xyz coord. of 12 element nodes
      real(8) xi,eta,zeta ! output=actual xi,eta,zeta
      integer iflag !=-1 far outside the box, 
                    !=+1 close to box
!
      real(8) DD, DDmin, DDmax, cfac
      real(8) xin(4), etan(4)
      data xin/-1d0,+1d0,+1d0,-1d0/
      data etan/-1d0,-1d0,+1d0,+1d0/
      real(8) c(3), cc(3), dc(3,8)
      data dc/-1d0,-1d0,-1d0,+1d0,-1d0,-1d0,
     1        +1d0,+1d0,-1d0,-1d0,+1d0,-1d0,
     2        -1d0,-1d0,+1d0,+1d0,-1d0,+1d0,
     1        +1d0,+1d0,+1d0,-1d0,+1d0,+1d0/
      real(8) xyzi(3), dxyzdc(3,3)
!
!--determine which of the corner nodes is the closest
!
      DDmin=vbig
      DDmax=0d0
      do is=1,4
         do lv=1,3,2
            DD=(xyz(1)-xyzn(1,lv,is))**2+(xyz(2)-xyzn(2,lv,is))**2
     1                                  +(xyz(3)-xyzn(3,lv,is))**2
            if (DD.lt.DDmin) then
               DDmin=DD
               cc(1)=xin(is)
               cc(2)=etan(is)
               cc(3)=lv-2d0
            endif
            DDmax=max(DD,DDmax)
         enddo
      enddo
!
!--check if one of the corner nodes is it.
      if (DDmin.lt.DDmax*1d-4) then
         c=cc
         return
      endif
!
      cfac=1d0
      do i=1,10
         DDmin=vbig
         do j=1,8 !loop around the 8 corners
            c=cc+cfac*dc(:,j)
            call getv12(xyzn,c,xyzi,dxyzdc)
            DD=(xyz(1)-xyzi(1))**2+(xyz(2)-xyzi(2))**2
     1                            +(xyz(3)-xyzi(3))**2
            if (DD.lt.DDmin) then
               DDmin=DD
               cc=c !update center
            endif
         enddo
         cfac=0.5d0*cfac !decrease size of cube by 2
      enddo
!
      return
      end subroutine getcoord2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getS4coord(xyz,xyzn,xi,eta,DD,iflag)
!
! Given the node positions xyzn of an S4 surface element, this routine
! compute the xi,eta corresponding to the minimum distance to a supplied 
! xyz positions.
! 
      implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,4)     ! xyz coord. of 4 nodes of element
      real(8) xi,eta        ! input=initial xi,eta guess
                            ! output=actual xi,eta
      real(8) DD            ! distance to surface
      integer iflag !=-1 far outside the box, 
                    !=0  no convergence in 20 steps
                    !=+1 close to box, good convergence
!
      real(8) v(3), dvdxi(3), dvdeta(3)
      real(8) d2vdxi2(3), d2vdeta2(3), d2vdxideta(3)
      real(8) dx, dy, dz
      real(8) dDDdxi, dDDdeta, d2DDdxi2, d2DDdeta2, d2DDdxideta
      real(8) det, dxi, deta
      real(8) xyz2(3,2), xi2, DD2
      integer i
!      
      iflag=0
      do i=1,20
         call getv4(xyzn,xi,eta,
     1              v,dvdxi,dvdeta,d2vdxi2,d2vdeta2,d2vdxideta)
         dx=v(1)-xyz(1)
         dy=v(2)-xyz(2)
         dz=v(3)-xyz(3)
         DD=dx**2+dy**2+dz**2
         dDDdxi=2d0*(dx*dvdxi(1)+dy*dvdxi(2)+dz*dvdxi(3))
         dDDdeta=2d0*(dx*dvdeta(1)+dy*dvdeta(2)+dz*dvdeta(3))
! note that we use here d2vdxi2=d2vdeta2=0
         d2DDdxi2=2d0*(dvdxi(1)**2+dvdxi(2)**2+dvdxi(3)**2)
         d2DDdeta2=2d0*(dvdeta(1)**2+dvdeta(2)**2+dvdeta(3)**2)
         d2DDdxideta=2d0*(dx*d2vdxideta(1)+dvdxi(1)*dvdeta(1)+
     2                    dy*d2vdxideta(2)+dvdxi(2)*dvdeta(2)+
     4                    dz*d2vdxideta(3)+dvdxi(3)*dvdeta(3))
         det=d2DDdxi2*d2DDdeta2-d2DDdxideta**2
         if (abs(det).gt.vtiny) then
            dxi=(dDDdeta*d2DDdxideta-dDDdxi*d2DDdeta2)/det
            deta=(dDDdxi*d2DDdxideta-dDDdeta*d2DDdxi2)/det
         else
            print *,'vanishing determinant in getS4coord',det
            stop
         endif
         xi=xi+dxi
         eta=eta+deta
         if (max(abs(xi),abs(eta)).gt.2d0) then
!--too far out 
            iflag=-1
            exit
         endif
         if (max(abs(dxi),abs(deta)).lt.1d-4) then
!--we have converged
            iflag=+1
            exit
         endif
      enddo
!
      if (max(abs(xi),abs(eta)).gt.1d0) then
         DD=vbig
!--we have to check the 4 borders of the quadrilateral 
         xyz2(:,1)=xyzn(:,1)
         xyz2(:,2)=xyzn(:,2)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=xi2; eta=-1d0
         endif
!
         xyz2(:,1)=xyzn(:,2)
         xyz2(:,2)=xyzn(:,3)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=+1d0; eta=xi2
         endif
!
         xyz2(:,1)=xyzn(:,4)
         xyz2(:,2)=xyzn(:,3)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=xi2; eta=+1d0
         endif
!
         xyz2(:,1)=xyzn(:,1)
         xyz2(:,2)=xyzn(:,4)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=-1d0; eta=xi2
         endif
      endif
!
      return 
      end subroutine getS4coord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getS6coord(xyz,xyzn,xi,zeta,DD,iflag)
!
! Given the node positions xyzn of an S6 surface element, this routine
! compute the xi,zeta corresponding to the minimum distance to a supplied 
! xyz positions.
! 
      implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,3,2)   ! xyz coord. of 6 nodes of element
      real(8) xi,zeta       ! input=initial xi,zeta guess
                            ! output=actual xi,zeta
      real(8) DD            ! distance to surface
      real(8) xyzS(3)       ! nearest location on surface
      integer iflag !=-1 far outside the box, 
                    !=0  no convergence in 20 steps
                    !=+1 close to box, good convergence
!
      real(8) v(3), dvdxi(3), dvdzeta(3)
      real(8) d2vdxi2(3), d2vdzeta2(3), d2vdxidzeta(3)
      real(8) dx, dy, dz
      real(8) dDDdxi, dDDdzeta, d2DDdxi2, d2DDdzeta2, d2DDdxidzeta
      real(8) det, dxi, dzeta, fnewton
      real(8) xyz2(3,2),xi2,DD2,xyz3(3,3),zeta3,DD3
      integer i
!      
      iflag=0
      do i=1,100
         call getv6(xyzn,xi,zeta,
     1              v,dvdxi,dvdzeta,d2vdxi2,d2vdzeta2,d2vdxidzeta)
         dx=v(1)-xyz(1)
         dy=v(2)-xyz(2)
         dz=v(3)-xyz(3)
         DD=dx**2+dy**2+dz**2
         dDDdxi=2d0*(dx*dvdxi(1)+dy*dvdxi(2)+dz*dvdxi(3))
         dDDdzeta=2d0*(dx*dvdzeta(1)+dy*dvdzeta(2)+dz*dvdzeta(3))
! note that we use here d2vdxi2=0
         d2DDdxi2=2d0*(dvdxi(1)**2+dvdxi(2)**2+dvdxi(3)**2)
         d2DDdzeta2=2d0*(dvdzeta(1)**2+dx*d2vdzeta2(1)+
     1                   dvdzeta(2)**2+dy*d2vdzeta2(2)+
     2                   dvdzeta(3)**2+dz*d2vdzeta2(3))
         d2DDdxidzeta=2d0*(dx*d2vdxidzeta(1)+dvdxi(1)*dvdzeta(1)+
     1                     dy*d2vdxidzeta(2)+dvdxi(2)*dvdzeta(2)+
     2                     dz*d2vdxidzeta(3)+dvdxi(3)*dvdzeta(3))
         det=d2DDdxi2*d2DDdzeta2-d2DDdxidzeta**2
         if (abs(det).gt.vtiny) then
            dxi=0.5d0*(dDDdzeta*d2DDdxidzeta-dDDdxi*d2DDdzeta2)/det
            dzeta=0.5d0*(dDDdxi*d2DDdxidzeta-dDDdzeta*d2DDdxi2)/det
         else
            print *,'vanishing determinant in getS6coord',det
            stop
         endif
!        print *,'i,DD',i,DD
!        print *,'xi,dxi,zeta,dzeta',xi,dxi,zeta,dzeta
!        print *,'xyz',xyz(:)
!        print *,'v',v(:)
         if (abs(dxi).gt.0.5) then
            fnewton=0.5/abs(dxi)
            dxi=fnewton*dxi
            dzeta=fnewton*dzeta
         endif
         if (abs(dzeta).gt.0.5) then
            fnewton=0.5/abs(dzeta)
            dxi=fnewton*dxi
            dzeta=fnewton*dzeta
         endif
         xi=xi+dxi
         zeta=zeta+dzeta
         if (max(abs(xi),abs(zeta)).gt.3.0) then
!--too far out 
            iflag=-1
            exit
         endif
         if (max(abs(dxi),abs(dzeta)).lt.1d-4) then
!--we have converged
            iflag=+1
            exit
         endif
      enddo
!
      if (max(abs(xi),abs(zeta)).gt.1d0) then
!--we have to check the borders of the surface
         DD=vbig
!--bottom
         xyz2(:,1)=xyzn(:,1,1)
         xyz2(:,2)=xyzn(:,1,2)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=xi2; zeta=-1d0
         endif
!--top
         xyz2(:,1)=xyzn(:,3,1)
         xyz2(:,2)=xyzn(:,3,2)
         xi2=0d0
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            xi=xi2; zeta=+1d0
         endif
!--side 1
         xyz3(:,:)=xyzn(:,:,1)
         zeta3=0d0
         call getL3coord(xyz,xyz3,zeta3,DD3,iflag)
         if (DD3.lt.DD) then
            DD=DD3
            xi=-1d0; zeta=zeta3
         endif
!--side 2
         xyz3(:,:)=xyzn(:,:,2)
         zeta3=0d0
         call getL3coord(xyz,xyz3,zeta3,DD3,iflag)
         if (DD3.lt.DD) then
            DD=DD3
            xi=+1d0; zeta=zeta3
         endif
      endif
         
      return 
      end subroutine getS6coord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getL2coord(xyz,xyzn,xi,DD,iflag)
!
! Given the node positions xyzn of an L2 line element, this routine
! compute the xi corresponding to the minimum distance to a supplied 
! xyz positions.
! 
      implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,2)     ! xyz coord. of 2 nodes of element
      real(8) xi            ! input=initial xi guess
                            ! output=actual xi
      real(8) DD            ! distance to line
      integer iflag !=0  no convergence in 20 steps
                    !=+1 close to box, good convergence
!
      real(8) v(3), dvdxi(3), d2vdxi2(3)
      real(8) dx, dy, dz
      real(8) dDDdxi, d2DDdxi2
      real(8) dxi
      integer i
!      
      iflag=0
      do i=1,20
         call getv2(xyzn,xi,v,dvdxi,d2vdxi2)
         dx=v(1)-xyz(1)
         dy=v(2)-xyz(2)
         dz=v(3)-xyz(3)
         DD=dx**2+dy**2+dz**2
         dDDdxi=2d0*(dx*dvdxi(1)+dy*dvdxi(2)+dz*dvdxi(3))
! note that we use here d2vdxi2=0
         d2DDdxi2=2d0*(dvdxi(1)**2+dvdxi(2)**2+dvdxi(3)**2)
         if (abs(d2DDdxi2).gt.vtiny) then
            dxi=-dDDdxi/d2DDdxi2
         else
            print *,'vanishing d2DDdxi2 in getL2coord',d2DDdxi2
            stop
         endif
         xi=xi+dxi
         if (abs(xi).gt.3.0) exit !--too far out 
         if (abs(dxi).lt.1d-4) then
!--we have converged
            iflag=+1
            exit
         endif
      enddo
!
      if (xi.lt.-1d0) then
         xi=-1d0
         dx=xyzn(1,1)-xyz(1)
         dy=xyzn(2,1)-xyz(2)
         dz=xyzn(3,1)-xyz(3)
         DD=dx**2+dy**2+dz**2
      elseif (xi.gt.+1d0) then
         xi=+1d0
         dx=xyzn(1,2)-xyz(1)
         dy=xyzn(2,2)-xyz(2)
         dz=xyzn(3,2)-xyz(3)
         DD=dx**2+dy**2+dz**2
      endif
!
      return 
      end subroutine getL2coord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine getL3coord(xyz,xyzn,zeta,DD,iflag)
!
! Given the node positions xyzn of an L3 line element, this routine
! compute the zeta corresponding to the minimum distance to a supplied 
! xyz positions.
! 
      implicit none
!
      real(8) xyz(3)        ! xyz coordinate of desired point
      real(8) xyzn(3,3)     ! xyz coord. of 3 nodes of element
      real(8) zeta          ! input=initial zeta guess
                            ! output=actual zeta
      real(8) DD            ! distance to surface
      integer iflag !=0  no convergence in 20 steps
                    !=+1 close to box, good convergence
!
      real(8) v(3), dvdzeta(3), d2vdzeta2(3)
      real(8) dx, dy, dz
      real(8) dDDdzeta, d2DDdzeta2
      real(8) dzeta
      integer i
!      
      iflag=0
      do i=1,20
         call getv3(xyzn,zeta,v,dvdzeta,d2vdzeta2)
         dx=v(1)-xyz(1)
         dy=v(2)-xyz(2)
         dz=v(3)-xyz(3)
         DD=dx**2+dy**2+dz**2
         dDDdzeta=2d0*(dx*dvdzeta(1)+
     1                 dy*dvdzeta(2)+v(3)*dz*dvdzeta(3))
         d2DDdzeta2=2d0*(dvdzeta(1)**2+dx*d2vdzeta2(1)+
     2                   dvdzeta(2)**2+dy*d2vdzeta2(2)+
     4                   dvdzeta(3)**2+dz*d2vdzeta2(3))
         if (abs(d2DDdzeta2).gt.vtiny) then
            dzeta=-dDDdzeta/d2DDdzeta2
         else
            print *,'vanishing d2DDdzeta2 in getS3coord',d2DDdzeta2
            stop
         endif
         zeta=zeta+dzeta
         if (abs(zeta).gt.3.0) exit !--too far out 
         if (abs(dzeta).lt.1d-4) then
!--we have converged
            iflag=+1
            exit
         endif
      enddo
!
      if (zeta.lt.-1d0) then
         zeta=-1d0
         dx=xyzn(1,1)-xyz(1)
         dy=xyzn(2,1)-xyz(2)
         dz=xyzn(3,1)-xyz(3)
         DD=dx**2+dy**2+dz**2
      elseif (zeta.gt.+1d0) then
         zeta=+1d0
         dx=xyzn(1,3)-xyz(1)
         dy=xyzn(2,3)-xyz(2)
         dz=xyzn(3,3)-xyz(3)
         DD=dx**2+dy**2+dz**2
      endif
!
      return 
      end subroutine getL3coord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine find_element(xyz,nxyz,qkount,qlist,kq,xi,eta,zeta)
!
!     given a point with coordinates x,y,z, the subroutine searches the
!     mesh for the element kq that contains it. if kq=0, no element
!     was found
!
      implicit none
!
      real(8) xyz(3) !point position
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
      integer qkount ! number of elements to search in list
      integer qlist(NSM) !list of elements to search
      integer kq !the element that contains the point
      real(8) xi,eta,zeta !the intrinsic coordinates of the point in kq 
!
      real(8) xyzel(3,3,4)
      integer k, iq, isn, istack, lvn, jx, iflag
!
      kq=0
      do k=1,qkount !loop over items in list
         iq=qlist(k)
         do isn=1,4 !load node positions of element iq
            istack=isoq(isn,iq)
            do lvn=1,3
               do jx=1,3
                  xyzel(jx,lvn,isn)=nxyz(jx,lvn,istack)
               enddo
            enddo
         enddo
!--initial guesses are for the center of the element
         xi=0d0
         eta=0d0
         zeta=0d0
         call getcoord(xyz,xyzel,xi,eta,zeta,iflag)
         if (iflag.eq.1) then
            if (abs(xi).le.1d0.and.abs(eta).le.1d0
     1                        .and.abs(zeta).le.1d0) then
               kq=iq !xyz inside iq
               exit !we are done
            endif
         endif
      enddo 
!
      return
      end subroutine find_element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine find_surface(xyz,nxyz,kq,xiq,etaq,zetaq)
!
!  given a point with coordinates x,y,z, the subroutine searches the
!  mesh for the surface nearest to it.
!
      implicit none
!
      real(8) xyz(3) !point position
      real(8) nxyz(3,3,NSM)
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
      integer kq !the element that contains the point
      real(8) xiq,etaq,zetaq !the intrinsic coordinates of the point
!
      real(8) areamin, DDmin, DDmin2, DD
      integer is, lvn, ismin, lvmin, il 
!
!-- find the closest surface node
!
      DDmin=vbig
!--ventral and dorsal nodes
      do is=1,ns
         do lvn=1,3,2
            DD=(xyz(1)-nxyz(1,lvn,is))**2+
     1         (xyz(2)-nxyz(2,lvn,is))**2+(xyz(3)-nxyz(3,lvn,is))**2
            if (DD.lt.DDmin) then
               DDmin2=DDmin
               DDmin=DD
               ismin=is
               lvmin=lvn
            endif
         enddo
      enddo
!-edge nodes
      do il=1,nl
         is=isol(1,il)
         DD=(xyz(1)-nxyz(1,2,is))**2+
     1      (xyz(2)-nxyz(2,2,is))**2+(xyz(3)-nxyz(3,2,is))**2
         if (DD.lt.DDmin) then
            DDmin2=DDmin
            DDmin=DD
            ismin=is
            lvmin=2
         endif
      enddo
!
      if (ilos(1,ismin).eq.0) then !interior stack, level 1 or 3
         call find_quad(xyz,nxyz(:,lvmin,:),ismin,kq,xiq,etaq,DD)
         zetaq=lvmin-2d0
      else !edge stack
         call find_edge(xyz,nxyz,ismin,lvmin,kq,xiq,etaq,zetaq,DD)
      endif
!         
      if (DD.gt.DDmin.and.DD.gt.DDmin2*1d-3) then
         print *,'error in find_surface'
         print *,'DD,DDmin',DD,DDmin
         print *,'ismin,lvmin',ismin,lvmin
         print *,'kq=',kq
         print *,'xiq,etaq,zetaq',xiq,etaq,zetaq

!        stop
      endif
!         
      areamin=arean(lvmin,ismin)
      if (DDmin.gt.2d0*areamin) then
         print *,'DDmin too large!',DDmin
         print *,'ismin,lvmin',ismin,lvmin
         print *,'areamin',areamin
         stop
      endif
!
      return
      end subroutine find_surface
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine find_quad(xyz,lxyz,ismin,kq,xiq,etaq,DDmin)
!
!  given a point with coordinates x,y,z, and the nearest edge node ismin,
!  lvmin=3, the subroutine searches the mesh for the quad surface element 
!  kq nearest to it, returns the intrinsic coordinates xiq,etaq, of the 
!  position on the surface closest to the point, and the distance DD. 
!
      implicit none
!
      real(8) xyz(3) !point position
      real(8) lxyz(3,NSM) ! xyz node positions at a given level;  stacks 1-NSM
      integer ismin !the stack nearest to the point (input)
      integer kq !the element which is nearest to the point (output)
      real(8) xiq,etaq !the intrinsic coordinates of the point in kq 
      real(8) DDmin !minimum distance square to point of interest
!
      real(8) xyz4(3,4), xi,eta,DD
      real(8) xi4(4),eta4(4)
      data xi4/-1d0,1d0,1d0,-1d0/
      data eta4/-1d0,-1d0,1d0,1d0/
      integer i, iq, isn, istack, isnmin, iflag
!
      DDmin=vbig
      do i=0,kqos(ismin)-1 
         iq=lqos(iqos(ismin)+i) !loop over elements containing stack ismin
!--load 4 nodes of surface
         do isn=1,4 
            istack=isoq(isn,iq)
            xyz4(:,isn)=lxyz(:,istack)
            if (istack.eq.ismin) isnmin=isn
         enddo
         xi=xi4(isnmin); eta=eta4(isnmin) !initial guesses
         call getS4coord(xyz,xyz4,xi,eta,DD,iflag)
         if (DD.lt.DDmin) then
            DDmin=DD
            kq=iq
            xiq=xi; etaq=eta
         endif
      enddo 
!
      if (abs(xiq).gt.1d0.or.abs(etaq).gt.1d0) then
         print *,'error in find_quad: ismin,lvmin',ismin
         print *,'DDmin,xiq,etaq,zetaq',DDmin,xiq,etaq
         stop
      endif
!
      return
      end subroutine find_quad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine find_edge(xyz,nxyz,ismin,lvmin,kq,xiq,etaq,zetaq,DDmin)
!
!  given a point with coordinates x,y,z, and the nearest edge node ismin,
!  lvmin, the subroutine searches the mesh for the edge element kq nearest 
!  to it, returns the intrinsic coordinates xiq,etaq,zetaq of the position
!  on the surface closest to the point, and the distance DDmin. 
!
      implicit none
!
      real(8) xyz(3) !point position
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
      integer ismin,lvmin !stack and level of closest node
      integer kq !the element that contains the point
      real(8) xiq,etaq,zetaq !the intrinsic coordinates of the point 
      real(8) DDmin
!
      real(8) xi,eta,zeta
      real(8) xyz6(3,3,2), xyz4(3,4), DD
      integer il, iq, lvn, istack, isn, i
      integer iflag
!
!--the closest edge stack belongs to 2 elements, find which
!   one contributes the closest surface to the point
!
      DDmin=vbig
!
      do i=1,2 ! note for i=1, xi of stack is +1, for i=2, xi=-1
         il=ilos(i,ismin)
         iq=iqol(il)
!--edge (frontal surface)
         do isn=1,2
            istack=isol(isn,il)
            xyz6(:,:,isn)=nxyz(:,:,istack)
         enddo
         xi=3d0-2d0*i; zeta=lvmin-2d0
         call getS6coord(xyz,xyz6,xi,zeta,DD,iflag)
         if (DD.lt.DDmin) then
            DDmin=DD
            kq=iq
            xiq=xi; etaq=-1d0; zetaq=zeta
         endif
!
!--ventral surface
         do isn=1,4 
            istack=isoq(isn,iq)
            xyz4(:,isn)=nxyz(:,1,istack)
         enddo
         xi=3d0-2d0*i; eta=-1d0
         call getS4coord(xyz,xyz4,xi,eta,DD,iflag)
         if (DD.lt.DDmin) then
            DDmin=DD
            kq=iq
            xiq=xi; etaq=eta; zetaq=-1d0
         endif
!
!--dorsal surface
         do isn=1,4 
            istack=isoq(isn,iq)
            xyz4(:,isn)=nxyz(:,3,istack)
         enddo
         xi=3d0-2d0*i; eta=-1d0
         call getS4coord(xyz,xyz4,xi,eta,DD,iflag)
         if (DD.lt.DDmin) then
            DDmin=DD
            kq=iq
            xiq=xi; etaq=eta; zetaq=+1d0
         endif
      enddo
!
      if (abs(xiq).gt.1d0.or.abs(etaq).gt.1d0
     1                   .or.abs(zetaq).gt.1d0) then
         print *,'error in find_edge: ismin,lvmin',ismin,lvmin
         print *,'DDmin,xiq,etaq,zetaq',DDmin,xiq,etaq,zetaq
         stop
      endif
!
      return
      end subroutine find_edge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine find_cl(xyz,nxyz,kq,xiq,etaq,zetaq)
!
!  given a point with coordinates x,y,z, the subroutine searches the
!  mesh for the contact line nearest to it.
      implicit none
!
      real(8) xyz(3) !point position
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
      integer kq !the element that contains the point
      real(8) xiq,etaq,zetaq !the intrinsic coordinates of the point in kq 
!
      real(8) xyz2(3,2), xi2, DD2, DDmin, DD
      integer il, is, ismin, i, iflag
!--find closest node on contact line
      DDmin=vbig
      do il=1,nl
         is=isol(1,il)
         DD=(xyz(1)-nxyz(1,1,is))**2+(xyz(2)-nxyz(2,1,is))**2+
     1                               (xyz(3)-nxyz(3,1,is))**2
         if (DD.lt.DDmin) then
            DDmin=DD
            ismin=is
         endif
      enddo
!
!--node belongs to two line elements
      DD=vbig
      do i=1,2
         il=ilos(i,ismin)
         xyz2(:,1)=nxyz(:,1,isol(1,il))
         xyz2(:,2)=nxyz(:,1,isol(2,il))
         xi2=3d0-2d0*i
         call getL2coord(xyz,xyz2,xi2,DD2,iflag)
         if (DD2.lt.DD) then
            DD=DD2
            kq=iqol(il)
            xiq=xi2; etaq=-1d0; zetaq=-1d0
         endif
      enddo
!
      if (DD.gt.DDmin.or.abs(xiq).gt.1d0) then 
         print *,'error in find_cl'
         print *,'DD,DDmin',DD,DDmin
         print *,'ismin=',ismin
         print *,'kq,xiq',kq,xiq
         stop
      endif
!
      return
      end subroutine find_cl
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine Mat3Inv(Mat,Inv)
!
! Given a 3x3 matrix Mat, the subroutine returns
! the inverse matrix Inv
!
      implicit none
!
      real(8) Mat(3,3)
      real(8) Inv(3,3)
!
      real(8) cofac(3,3)
      real(8) ddet
      integer i, j
!
!--compute cofactor matrix
      cofac(1,1)=+Mat(2,2)*Mat(3,3)-Mat(3,2)*Mat(2,3)
      cofac(2,1)=-Mat(1,2)*Mat(3,3)+Mat(3,2)*Mat(1,3)
      cofac(3,1)=+Mat(1,2)*Mat(2,3)-Mat(2,2)*Mat(1,3)
      cofac(1,2)=-Mat(2,1)*Mat(3,3)+Mat(3,1)*Mat(2,3)
      cofac(2,2)=+Mat(1,1)*Mat(3,3)-Mat(3,1)*Mat(1,3)
      cofac(3,2)=-Mat(1,1)*Mat(2,3)+Mat(2,1)*Mat(1,3)
      cofac(1,3)=+Mat(2,1)*Mat(3,2)-Mat(3,1)*Mat(2,2)
      cofac(2,3)=-Mat(1,1)*Mat(3,2)+Mat(3,1)*Mat(1,2)
      cofac(3,3)=+Mat(1,1)*Mat(2,2)-Mat(2,1)*Mat(1,2)
!--determinant
      ddet=1d0/(Mat(1,1)*cofac(1,1)+Mat(2,1)*cofac(2,1)+
     1                              Mat(3,1)*cofac(3,1))
!     print *,'ddet',ddet
      if (abs(ddet).gt.vbig) then
         print *,'Mat3Inv: determinant nearly 0!',ddet
      endif
!
      do j=1,3
         do i=1,3
            Inv(i,j)=ddet*cofac(j,i)
         enddo
      enddo
!
      return
      end subroutine Mat3Inv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine goisomap(nxyz)
!
! This subroutine computes shape function quantities for each element
! at each GP and differentiated wrt to extrinsic xyz coordinates      
!
!  det(3,4,NSM) Jacobian: GP level, GP stack, element
!  dHg(3,3,4,3,4,NSM) Shape function gradient:
!       x-y-z; node level; node stack; GP level; GP stack; element
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      logical,save::firstwarn=.true.
      real(8) x(3,4), y(3,4), z(3,4), dx(3), dy(3), dz(3)
      real(8) cofacx(3), cofacy(3), cofacz(3), ddet, chkdet
      real(8) xg,yg,zg,dh1,dh2,dh3
      integer iq,isn,lvn,istack,isg,lvg,j,iw,lw
!
      chkdet=vbig
!--loop through the elements
      do iq=1,nq
!--store node coordinates: 4 stacks, 3 levels = 12 nodes
         do isn=1,4
            istack=isoq(isn,iq)
            do lvn=1,3
               x(lvn,isn)=nxyz(1,lvn,istack)
               y(lvn,isn)=nxyz(2,lvn,istack)
               z(lvn,isn)=nxyz(3,lvn,istack)
            enddo
         enddo
!--loop over Gauss points in cell
         do isg=1,4
            do lvg=1,3
!--zero out Jacobian entries (e.g. dx(2)=x*(d_2 H))
               dx=0.0d0
               dy=0.0d0
               dz=0.0d0
               do isn=1,4
                  do lvn=1,3
                     do j=1,3
                        dx(j)=dx(j)+x(lvn,isn)*dH(j,lvn,isn,lvg,isg)
                        dy(j)=dy(j)+y(lvn,isn)*dH(j,lvn,isn,lvg,isg)
                        dz(j)=dz(j)+z(lvn,isn)*dH(j,lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
               do j=1,3
                  dxV(j,lvg,isg,iq)=dx(j)
                  dyV(j,lvg,isg,iq)=dy(j)
                  dzV(j,lvg,isg,iq)=dz(j)
               enddo
!--cofactor matrix
               cofacx(1)=+dy(2)*dz(3)-dy(3)*dz(2)
               cofacx(2)=-dy(1)*dz(3)+dy(3)*dz(1)
               cofacx(3)=+dy(1)*dz(2)-dy(2)*dz(1)
               cofacy(1)=-dx(2)*dz(3)+dx(3)*dz(2)
               cofacy(2)=+dx(1)*dz(3)-dx(3)*dz(1)
               cofacy(3)=-dx(1)*dz(2)+dx(2)*dz(1)
               cofacz(1)=+dx(2)*dy(3)-dx(3)*dy(2)
               cofacz(2)=-dx(1)*dy(3)+dx(3)*dy(1)
               cofacz(3)=+dx(1)*dy(2)-dx(2)*dy(1)
!--determinant of Jacobian
               det(lvg,isg,iq)=+dx(1)*cofacx(1)
     1                         +dx(2)*cofacx(2)
     2                         +dx(3)*cofacx(3)
               if (det(lvg,isg,iq).le.vtiny) then
                  print*,'NEGATIVE Q-ELEMENT IN GOISOMAP',iq
                  print *,'lvg,isg,det=',lvg,isg,det(lvg,isg,iq)
                  if (firstwarn) then
                     do iw=1,ns
                        do lw=1,3
                           hvec(1,lw,iw)=nxyz(1,lw,iw)
                           hvec(2,lw,iw)=nxyz(2,lw,iw)
                           hvec(3,lw,iw)=nxyz(3,lw,iw)
                        enddo
                     enddo
                     open(66,file='isomap.dump')
                     call iowrfile(99,66)
                     close(66)
                     print *,'enter to continue'
                     read *
                     firstwarn=.false.
                  endif
               endif 
               ddet=1.0/max(det(lvg,isg,iq),vtiny)
!--for each node function, derivative wrt to x y z coord at GP
               do isn=1,4
                  do lvn=1,3
                     dHg(1,lvn,isn,lvg,isg,iq)=ddet*(
     1                      cofacx(1)*dH(1,lvn,isn,lvg,isg)+
     2                      cofacx(2)*dH(2,lvn,isn,lvg,isg)+
     3                      cofacx(3)*dH(3,lvn,isn,lvg,isg))
                     dHg(2,lvn,isn,lvg,isg,iq)=ddet*(
     1                      cofacy(1)*dH(1,lvn,isn,lvg,isg)+
     2                      cofacy(2)*dH(2,lvn,isn,lvg,isg)+
     3                      cofacy(3)*dH(3,lvn,isn,lvg,isg))
                     dHg(3,lvn,isn,lvg,isg,iq)=ddet*(
     1                      cofacz(1)*dH(1,lvn,isn,lvg,isg)+
     2                      cofacz(2)*dH(2,lvn,isn,lvg,isg)+
     3                      cofacz(3)*dH(3,lvn,isn,lvg,isg))
                  enddo
               enddo
            enddo
         enddo
      enddo 
!
      return
      end subroutine goisomap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gobdisom(nxyz)
!
! This subroutine computes boundary surface shape function-related 
! quantities at Gauss points
!
! derivatives of 3-space extrinsic x,y,z wrt to intrinsic 
! 2-space xi and eta (or zeta) at GP on ventral, dorsal, edge surfaces
!  dxSv(2,4,NSM),dySv(2,4,NSM),dzSv(2,4,NSM) xi-eta; GP stack; element
!  dxSd(2,4,NSM),dySd(2,4,NSM),dzSd(2,4,NSM) xi-eta; GP stack; element
!  dxSe(2,3,2,NLM),dySe(2,3,2,NLM),dzSe(2,3,2,NLM)
!                            xi-zeta; GP level; GP stack; edge 
!
! (dxS(1),dyS(1),dzS(1))^(dxS(2),dyS(2),dzS(2))
!                     = elementary surface vector at Gauss point
!  dAv(0:3,4,NSM) ventral: norm-x-y-z ; GP stack; element
!  dAd(0:3,4,NSM) dorsal:  norm-x-y-z ; GP stack; element
!  dAe(0:3,3,2,NLM) edge:
!             norm,x-y-z components; GP level; GP stack; edge 
!
! 2-dimensional surface gradients wrt to extrinsic xyz coordinates
!  dS4vg(3,4,4,NSM) ventral: xyz; node stack; GP stack; element
!  dS4dg(3,4,4,NSM) dorsal:  xyz; node stack; GP stack; element
!  dS6g(3,3,2,3,2,NLM) edge: 
!       xyz; node level, node stack; GP level, GP stack, edge 
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8) x4(4), y4(4), z4(4)
      real(8) x6(3,2), y6(3,2), z6(3,2)
      real(8) dx(2), dy(2), dz(2)
      real(8) e1(3), e2(3), e3(3), J2D(2,2), Jinv(2,2)
      real(8) fac, dAnorm, dSe1, dSe2, e1norm, e3norm, ddetJ,e3dS
      integer iq,il,isn,istack,isg,j,lvn,lvg
!
!--quadrilateral surfaces 
      do iq=1,nq
!--load nodes for ventral surface (level=1)
         do isn=1,4
            istack=isoq(isn,iq)
            x4(isn)=nxyz(1,1,istack)
            y4(isn)=nxyz(2,1,istack)
            z4(isn)=nxyz(3,1,istack)
         enddo
!--loop over Gauss points on surface 
         do isg=1,4
!--zero out tangents entries 
            dx=0.0d0
            dy=0.0d0
            dz=0.0d0
            do isn=1,4
               do j=1,2
                  dx(j)=dx(j)+x4(isn)*dS4(j,isn,isg)
                  dy(j)=dy(j)+y4(isn)*dS4(j,isn,isg)
                  dz(j)=dz(j)+z4(isn)*dS4(j,isn,isg)
               enddo
            enddo
            do j=1,2
               dxSv(j,isg,iq)=dx(j)
               dySv(j,isg,iq)=dy(j)
               dzSv(j,isg,iq)=dz(j)
            enddo
!--outward elementary ventral surface area vector
            dAv(1,isg,iq)=-dy(1)*dz(2)+dy(2)*dz(1)
            dAv(2,isg,iq)=+dx(1)*dz(2)-dx(2)*dz(1)
            dAv(3,isg,iq)=-dx(1)*dy(2)+dx(2)*dy(1)
            dAv(0,isg,iq)=dsqrt(dAv(1,isg,iq)**2+dAv(2,isg,iq)**2+
     1                                           dAv(3,isg,iq)**2)
!--compute units vectors:
!--first is d/dxi
            e1norm=1d0/dsqrt(dx(1)**2+dy(1)**2+dz(1)**2)
            e1(1)=dx(1)*e1norm
            e1(2)=dy(1)*e1norm
            e1(3)=dz(1)*e1norm
!--third is d/deta x d/dxi
            e3norm=1d0/dAv(0,isg,iq)
            e3(1)=dAv(1,isg,iq)*e3norm
            e3(2)=dAv(2,isg,iq)*e3norm
            e3(3)=dAv(3,isg,iq)*e3norm
!--second unit vector is -e3xe1
            e2(1)=-e3(2)*e1(3)+e3(3)*e1(2)
            e2(2)=+e3(1)*e1(3)-e3(3)*e1(1)
            e2(3)=-e3(1)*e1(2)+e3(2)*e1(1)
!--compute the 2D Jacobian matrix=ei_j dj_k (in the e1 e2 coords), invert
            J2D(1,1)=e1(1)*dx(1)+e1(2)*dy(1)+e1(3)*dz(1)
            J2D(2,1)=e2(1)*dx(1)+e2(2)*dy(1)+e2(3)*dz(1)
            J2D(1,2)=e1(1)*dx(2)+e1(2)*dy(2)+e1(3)*dz(2)
            J2D(2,2)=e2(1)*dx(2)+e2(2)*dy(2)+e2(3)*dz(2)
            ddetJ=1d0/(J2D(1,1)*J2D(2,2)-J2D(1,2)*J2D(2,1))
            Jinv(1,1)=+J2D(2,2)*ddetj
            Jinv(2,1)=-J2D(1,2)*ddetj
            Jinv(1,2)=-J2D(1,2)*ddetj
            Jinv(2,2)=+J2D(1,1)*ddetj
            do isn=1,4
!--compute dS4 in e1 e2 coordinates
               dSe1=dS4(1,isn,isg)*Jinv(1,1)+dS4(2,isn,isg)*Jinv(2,1)
               dSe2=dS4(1,isn,isg)*Jinv(1,2)+dS4(2,isn,isg)*Jinv(2,2)
!--return to original coordinates
               dS4vg(1,isn,isg,iq)=dSe1*e1(1)+dSe2*e2(1)
               dS4vg(2,isn,isg,iq)=dSe1*e1(2)+dSe2*e2(2)
               dS4vg(3,isn,isg,iq)=dSe1*e1(3)+dSe2*e2(3)
            enddo
         enddo
      enddo
!
!--dorsal surfaces:
      do iq=1,nq
!--load nodes for dorsal surface (level=3)
         do isn=1,4
            istack=isoq(isn,iq)
            x4(isn)=nxyz(1,3,istack)
            y4(isn)=nxyz(2,3,istack)
            z4(isn)=nxyz(3,3,istack)
         enddo
!--loop over Gauss points in cell
         do isg=1,4
!--zero out tangents entries 
            dx=0.0d0
            dy=0.0d0
            dz=0.0d0
            do isn=1,4
               do j=1,2
                  dx(j)=dx(j)+x4(isn)*dS4(j,isn,isg)
                  dy(j)=dy(j)+y4(isn)*dS4(j,isn,isg)
                  dz(j)=dz(j)+z4(isn)*dS4(j,isn,isg)
               enddo
            enddo
            do j=1,2
               dxSd(j,isg,iq)=dx(j)
               dySd(j,isg,iq)=dy(j)
               dzSd(j,isg,iq)=dz(j)
            enddo
!--outward elementary dorsal surface area vector
            dAd(1,isg,iq)=+dy(1)*dz(2)-dy(2)*dz(1)
            dAd(2,isg,iq)=-dx(1)*dz(2)+dx(2)*dz(1)
            dAd(3,isg,iq)=+dx(1)*dy(2)-dx(2)*dy(1)
            dAd(0,isg,iq)=dsqrt(dAd(1,isg,iq)**2+dAd(2,isg,iq)**2+
     1                                           dAd(3,isg,iq)**2)
!--compute units vectors:
!--first is d/dxi
            e1norm=1d0/dsqrt(dx(1)**2+dy(1)**2+dz(1)**2)
            e1(1)=dx(1)*e1norm
            e1(2)=dy(1)*e1norm
            e1(3)=dz(1)*e1norm
!--third is d/dxi x d/deta
            e3norm=1d0/dAd(0,isg,iq)
            e3(1)=dAd(1,isg,iq)*e3norm
            e3(2)=dAd(2,isg,iq)*e3norm
            e3(3)=dAd(3,isg,iq)*e3norm
!--second unit vector is e3xe1
            e2(1)=+e3(2)*e1(3)-e3(3)*e1(2)
            e2(2)=-e3(1)*e1(3)+e3(3)*e1(1)
            e2(3)=+e3(1)*e1(2)-e3(2)*e1(1)
!--compute the 2D Jacobian matrix=ei_j dj_k (in the e1 e2 coords), invert
            J2D(1,1)=e1(1)*dx(1)+e1(2)*dy(1)+e1(3)*dz(1)
            J2D(2,1)=e2(1)*dx(1)+e2(2)*dy(1)+e2(3)*dz(1)
            J2D(1,2)=e1(1)*dx(2)+e1(2)*dy(2)+e1(3)*dz(2)
            J2D(2,2)=e2(1)*dx(2)+e2(2)*dy(2)+e2(3)*dz(2)
            ddetJ=1d0/(J2D(1,1)*J2D(2,2)-J2D(1,2)*J2D(2,1))
            Jinv(1,1)=+J2D(2,2)*ddetj
            Jinv(2,1)=-J2D(1,2)*ddetj
            Jinv(1,2)=-J2D(1,2)*ddetj
            Jinv(2,2)=+J2D(1,1)*ddetj
            do isn=1,4
!--compute dS4 in e1 e2 coordinates
               dSe1=dS4(1,isn,isg)*Jinv(1,1)+dS4(2,isn,isg)*Jinv(2,1)
               dSe2=dS4(1,isn,isg)*Jinv(1,2)+dS4(2,isn,isg)*Jinv(2,2)
!--return to original coordinates
               dS4dg(1,isn,isg,iq)=dSe1*e1(1)+dSe2*e2(1)
               dS4dg(2,isn,isg,iq)=dSe1*e1(2)+dSe2*e2(2)
               dS4dg(3,isn,isg,iq)=dSe1*e1(3)+dSe2*e2(3)
            enddo
         enddo
      enddo
!
!--do edge surface elements
      do il=1,nl
!--load edge surface nodes
         do isn=1,2
            istack=isol(isn,il)
            do lvn=1,3
               x6(lvn,isn)=nxyz(1,lvn,istack)
               y6(lvn,isn)=nxyz(2,lvn,istack)
               z6(lvn,isn)=nxyz(3,lvn,istack)
            enddo
         enddo
         do isg=1,2
            do lvg=1,3
!--zero out tangents entries 
               dx=0.0d0
               dy=0.0d0
               dz=0.0d0
               do isn=1,2
                  do lvn=1,3
                     do j=1,2
                        dx(j)=dx(j)+x6(lvn,isn)*dS6(j,lvn,isn,lvg,isg)
                        dy(j)=dy(j)+y6(lvn,isn)*dS6(j,lvn,isn,lvg,isg)
                        dz(j)=dz(j)+z6(lvn,isn)*dS6(j,lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
               do j=1,2
                  dxSe(j,lvg,isg,il)=dx(j)
                  dySe(j,lvg,isg,il)=dy(j)
                  dzSe(j,lvg,isg,il)=dz(j)
               enddo
!--outward elementary edge surface area vector
               dAe(1,lvg,isg,il)=+dy(1)*dz(2)-dy(2)*dz(1)
               dAe(2,lvg,isg,il)=-dx(1)*dz(2)+dx(2)*dz(1)
               dAe(3,lvg,isg,il)=+dx(1)*dy(2)-dx(2)*dy(1)
               dAe(0,lvg,isg,il)=dsqrt(dAe(1,lvg,isg,il)**2+
     1                                 dAe(2,lvg,isg,il)**2+
     1                                 dAe(3,lvg,isg,il)**2)
!--compute units vectors:
!--first is d/dxi
               e1norm=1d0/dsqrt(dx(1)**2+dy(1)**2+dz(1)**2)
               e1(1)=dx(1)*e1norm
               e1(2)=dy(1)*e1norm
               e1(3)=dz(1)*e1norm
!--third is d/dxi x d/dzeta
               e3norm=1d0/dAe(0,lvg,isg,il)
               e3(1)=dAe(1,lvg,isg,il)*e3norm
               e3(2)=dAe(2,lvg,isg,il)*e3norm
               e3(3)=dAe(3,lvg,isg,il)*e3norm
!--second unit vector is e3xe1
               e2(1)=+e3(2)*e1(3)-e3(3)*e1(2)
               e2(2)=-e3(1)*e1(3)+e3(3)*e1(1)
               e2(3)=+e3(1)*e1(2)-e3(2)*e1(1)
!--compute the 2D Jacobian matrix=ei_j dj_k (in the e1 e2 coords), invert
               J2D(1,1)=e1(1)*dx(1)+e1(2)*dy(1)+e1(3)*dz(1)
               J2D(2,1)=e2(1)*dx(1)+e2(2)*dy(1)+e2(3)*dz(1)
               J2D(1,2)=e1(1)*dx(2)+e1(2)*dy(2)+e1(3)*dz(2)
               J2D(2,2)=e2(1)*dx(2)+e2(2)*dy(2)+e2(3)*dz(2)
               ddetJ=1d0/(J2D(1,1)*J2D(2,2)-J2D(1,2)*J2D(2,1))
               Jinv(1,1)=+J2D(2,2)*ddetj
               Jinv(2,1)=-J2D(1,2)*ddetj
               Jinv(1,2)=-J2D(1,2)*ddetj
               Jinv(2,2)=+J2D(1,1)*ddetj
               do isn=1,2
                  do lvn=1,3
!--compute dS6 in e1 e2 coordinates
                     dSe1=dS6(1,lvn,isn,lvg,isg)*Jinv(1,1)+
     1                    dS6(2,lvn,isn,lvg,isg)*Jinv(2,1)
                     dSe2=dS6(1,lvn,isn,lvg,isg)*Jinv(1,2)+
     1                    dS6(2,lvn,isn,lvg,isg)*Jinv(2,2)
!--return to original coordinates
                     dS6g(1,lvn,isn,lvg,isg,il)=dSe1*e1(1)+dSe2*e2(1)
                     dS6g(2,lvn,isn,lvg,isg,il)=dSe1*e1(2)+dSe2*e2(2)
                     dS6g(3,lvn,isn,lvg,isg,il)=dSe1*e1(3)+dSe2*e2(3)
                  enddo
               enddo
            enddo
         enddo
      enddo   
      end subroutine gobdisom
!
!  end of shape function routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosnn(nxyz)
! 
! this subroutine computes the outward normal at each surface nodes
! the method used is to find the gradient of the area with respect
! to motion of a node (see surface tension routine in molib for more
! details) and to define that as the normal direction if the curvature 
! is zero, we use the normals computed from the nearest GP instead
!
!   snn(0:3,3,NSM) norm of gradA-x-y-z; node level; node stack
!      if norm=0 interior node, or no curvature.
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8) snngp(3,3,NSM),wnn(3,NSM)
! snngp is weighted average GP normals, and wnn is the sum of weights
      real(8) dAvx,dAvy,dAvz,ddAvn
      real(8) dAdx,dAdy,dAdz,ddAdn
      real(8) dAex,dAey,dAez,ddAen
      real(8) dAx,dAy,dAz
      real(8) fac1x,fac1y,fac1z,fac2x,fac2y,fac2z
      real(8) xg,yg,zg, dng
      real(8) snnn, snngpn, sense, dwnn
      integer iq, isg, isn, il, lvg, lvn, istack, is, ix
!
      arean=0d0
      snn=0d0
      snngp=0d0 
      wnn=0d0
!--loop over ventral surfaces 
      do iq=1,nq
!--loop over 4 Gauss points 
         do isg=1,4
            dAvx=dAv(1,isg,iq)
            dAvy=dAv(2,isg,iq)
            dAvz=dAv(3,isg,iq)
            ddAvn=1d0/dAv(0,isg,iq)
            fac1x=ddAvn*(dAvz*dySv(2,isg,iq)-dAvy*dzSv(2,isg,iq))
            fac1y=ddAvn*(dAvx*dzSv(2,isg,iq)-dAvz*dxSv(2,isg,iq))
            fac1z=ddAvn*(dAvy*dxSv(2,isg,iq)-dAvx*dySv(2,isg,iq))
            fac2x=ddAvn*(dAvz*dySv(1,isg,iq)-dAvy*dzSv(1,isg,iq))
            fac2y=ddAvn*(dAvx*dzSv(1,isg,iq)-dAvz*dxSv(1,isg,iq))
            fac2z=ddAvn*(dAvy*dxSv(1,isg,iq)-dAvx*dySv(1,isg,iq))
c--contribution of this GP to 4 nodes
            xg=0d0
            yg=0d0
            zg=0d0
            do isn=1,4
               dAx=dS4(1,isn,isg)*fac1x-dS4(2,isn,isg)*fac2x
               dAy=dS4(1,isn,isg)*fac1y-dS4(2,isn,isg)*fac2y
               dAz=dS4(1,isn,isg)*fac1z-dS4(2,isn,isg)*fac2z
               istack=isoq(isn,iq)
               snn(1,1,istack)=snn(1,1,istack)-dAx
               snn(2,1,istack)=snn(2,1,istack)-dAy
               snn(3,1,istack)=snn(3,1,istack)-dAz
               xg=xg+nxyz(1,1,istack)*S4(isn,isg)
               yg=yg+nxyz(2,1,istack)*S4(isn,isg)
               zg=zg+nxyz(3,1,istack)*S4(isn,isg)
               arean(1,istack)=arean(1,istack)+dAv(0,isg,iq)*S4(isn,isg)
            enddo
            istack=isoq(isg,iq)
! weight is inversely prop. to distance
            dng=1d0/sqrt((nxyz(1,1,istack)-xg)**2+
     1           (nxyz(2,1,istack)-yg)**2+(nxyz(3,1,istack)-zg)**2)
            wnn(1,istack)=wnn(1,istack)+dng
            snngp(1,1,istack)=snngp(1,1,istack)+dAvx*dng
            snngp(2,1,istack)=snngp(2,1,istack)+dAvy*dng
            snngp(3,1,istack)=snngp(3,1,istack)+dAvz*dng
         enddo
      enddo
!--loop over dorsal surfaces 
      do iq=1,nq
         do isg=1,4
            dAdx=dAd(1,isg,iq)
            dAdy=dAd(2,isg,iq)
            dAdz=dAd(3,isg,iq)
            ddAdn=1d0/dAd(0,isg,iq)
            fac1x=ddAdn*(dAdz*dySd(2,isg,iq)-dAdy*dzSd(2,isg,iq))
            fac1y=ddAdn*(dAdx*dzSd(2,isg,iq)-dAdz*dxSd(2,isg,iq))
            fac1z=ddAdn*(dAdy*dxSd(2,isg,iq)-dAdx*dySd(2,isg,iq))
            fac2x=ddAdn*(dAdz*dySd(1,isg,iq)-dAdy*dzSd(1,isg,iq))
            fac2y=ddAdn*(dAdx*dzSd(1,isg,iq)-dAdz*dxSd(1,isg,iq))
            fac2z=ddAdn*(dAdy*dxSd(1,isg,iq)-dAdx*dySd(1,isg,iq))
!--contribution of this GP to 4 nodes
            xg=0d0
            yg=0d0
            zg=0d0
            do isn=1,4
               dAx=dS4(1,isn,isg)*fac1x-dS4(2,isn,isg)*fac2x
               dAy=dS4(1,isn,isg)*fac1y-dS4(2,isn,isg)*fac2y
               dAz=dS4(1,isn,isg)*fac1z-dS4(2,isn,isg)*fac2z
!--dA sign different from above because dorsal node
!--ordering different from ventral
               istack=isoq(isn,iq)
               snn(1,3,istack)=snn(1,3,istack)+dAx
               snn(2,3,istack)=snn(2,3,istack)+dAy
               snn(3,3,istack)=snn(3,3,istack)+dAz
               xg=xg+nxyz(1,3,istack)*S4(isn,isg)
               yg=yg+nxyz(2,3,istack)*S4(isn,isg)
               zg=zg+nxyz(3,3,istack)*S4(isn,isg)
               arean(3,istack)=arean(3,istack)+dAd(0,isg,iq)*S4(isn,isg)
            enddo
            istack=isoq(isg,iq)
            dng=1d0/sqrt((nxyz(1,3,istack)-xg)**2+
     1           (nxyz(2,3,istack)-yg)**2+(nxyz(3,3,istack)-zg)**2)
            wnn(3,istack)=wnn(3,istack)+dng
            snngp(1,3,istack)=snngp(1,3,istack)+dAdx*dng
            snngp(2,3,istack)=snngp(2,3,istack)+dAdy*dng
            snngp(3,3,istack)=snngp(3,3,istack)+dAdz*dng
         enddo
      enddo
!--loop over edge surfaces 
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               dAex=dAe(1,lvg,isg,il)
               dAey=dAe(2,lvg,isg,il)
               dAez=dAe(3,lvg,isg,il)
               ddAen=wgp(lvg)/dAe(0,lvg,isg,il)
               fac1x=ddAen*(dAez*dySe(2,lvg,isg,il)-
     1                      dAey*dzSe(2,lvg,isg,il))
               fac1y=ddAen*(dAex*dzSe(2,lvg,isg,il)-
     1                      dAez*dxSe(2,lvg,isg,il))
               fac1z=ddAen*(dAey*dxSe(2,lvg,isg,il)-
     1                      dAex*dySe(2,lvg,isg,il))
               fac2x=ddAen*(dAez*dySe(1,lvg,isg,il)-
     1                      dAey*dzSe(1,lvg,isg,il))
               fac2y=ddAen*(dAex*dzSe(1,lvg,isg,il)-
     1                      dAez*dxSe(1,lvg,isg,il))
               fac2z=ddAen*(dAey*dxSe(1,lvg,isg,il)-
     1                      dAex*dySe(1,lvg,isg,il))
!--contribution of this GP to 6 nodes
               xg=0d0
               yg=0d0
               zg=0d0
               do isn=1,2
                  istack=isol(isn,il)
                  do lvn=1,3
                     dAx=dS6(1,lvn,isn,lvg,isg)*fac1x-
     1                   dS6(2,lvn,isn,lvg,isg)*fac2x
                     dAy=dS6(1,lvn,isn,lvg,isg)*fac1y-
     1                   dS6(2,lvn,isn,lvg,isg)*fac2y
                     dAz=dS6(1,lvn,isn,lvg,isg)*fac1z-
     1                   dS6(2,lvn,isn,lvg,isg)*fac2z
                     snn(1,lvn,istack)=snn(1,lvn,istack)+dAx
                     snn(2,lvn,istack)=snn(2,lvn,istack)+dAy
                     snn(3,lvn,istack)=snn(3,lvn,istack)+dAz
                     xg=xg+nxyz(1,lvn,istack)*S6(lvn,isn,lvg,isg)
                     yg=yg+nxyz(2,lvn,istack)*S6(lvn,isn,lvg,isg)
                     zg=zg+nxyz(3,lvn,istack)*S6(lvn,isn,lvg,isg)
                     arean(lvn,istack)=arean(lvn,istack)+
     1                     dAe(0,lvg,isg,il)*S6(lvn,isn,lvg,isg)
                  enddo
               enddo
               istack=isol(isg,il)
               dng=1d0/sqrt((nxyz(1,2,istack)-xg)**2+
     1              (nxyz(2,2,istack)-yg)**2+(nxyz(3,2,istack)-zg)**2)
               wnn(2,istack)=wnn(2,istack)+dng
               snngp(1,lvg,istack)=snngp(1,lvg,istack)+dAex*dng
               snngp(2,lvg,istack)=snngp(2,lvg,istack)+dAey*dng
               snngp(3,lvg,istack)=snngp(3,lvg,istack)+dAez*dng
            enddo
         enddo
      enddo
!--normalize the weighted average per GP
      do is=1,ns
         do lvn=1,3
            dwnn=1d0/(vtiny+wnn(lvn,istack))
            do ix=1,3
               snngp(ix,lvn,is)=snngp(ix,lvn,is)*dwnn
            enddo
         enddo
      enddo
!--loop over ventral nodes
      do is=1,ns
         snnn=sqrt(snn(1,1,is)**2+snn(2,1,is)**2+snn(3,1,is)**2)
         snngpn=sqrt(snngp(1,1,is)**2+snngp(2,1,is)**2+snngp(3,1,is)**2)
! snn is a surface gradient with units [L]
! snngp is a surface with units [L]**2
         if (snnn.lt.1d-9*sqrt(snngpn)) then
!--0 mean curvature, compute normal per snngp
            snn(0,1,is)=0d0
            snn(1,1,is)=snngp(1,1,is)/snngpn
            snn(2,1,is)=snngp(2,1,is)/snngpn
            snn(3,1,is)=snngp(3,1,is)/snngpn
         else
! snn (grad A) is not necessarily outward normal (inward for curv. <0)
            sense=snn(1,1,is)*snngp(1,1,is)+snn(2,1,is)*snngp(2,1,is)
     1                                     +snn(3,1,is)*snngp(3,1,is)
            if (sense.ge.0d0) then
               snn(0,1,is)=snnn
               snn(1,1,is)=snn(1,1,is)/snnn
               snn(2,1,is)=snn(2,1,is)/snnn
               snn(3,1,is)=snn(3,1,is)/snnn
            else
               snn(0,1,is)=-snnn
               snn(1,1,is)=-snn(1,1,is)/snnn
               snn(2,1,is)=-snn(2,1,is)/snnn
               snn(3,1,is)=-snn(3,1,is)/snnn
            endif
         endif
      enddo
!--loop over dorsal nodes
      do is=1,ns
         snnn=sqrt(snn(1,3,is)**2+snn(2,3,is)**2+snn(3,3,is)**2)
         snngpn=sqrt(snngp(1,3,is)**2+snngp(2,3,is)**2+snngp(3,3,is)**2)
         if (snnn.lt.1d-9*sqrt(snngpn)) then
            snn(0,3,is)=0d0
            snn(1,3,is)=snngp(1,3,is)/snngpn
            snn(2,3,is)=snngp(2,3,is)/snngpn
            snn(3,3,is)=snngp(3,3,is)/snngpn
         else
            sense=snn(1,3,is)*snngp(1,3,is)+snn(2,3,is)*snngp(2,3,is)
     1                                     +snn(3,3,is)*snngp(3,3,is)
            if (sense.ge.0d0) then
               snn(0,3,is)=snnn
               snn(1,3,is)=snn(1,3,is)/snnn
               snn(2,3,is)=snn(2,3,is)/snnn
               snn(3,3,is)=snn(3,3,is)/snnn
            else
               snn(0,3,is)=-snnn
               snn(1,3,is)=-snn(1,3,is)/snnn
               snn(2,3,is)=-snn(2,3,is)/snnn
               snn(3,3,is)=-snn(3,3,is)/snnn
            endif
         endif
      enddo
!--loop over edge nodes (lv=2 only)
      do il=1,nl
         is=isol(1,il)
         snnn=sqrt(snn(1,2,is)**2+snn(2,2,is)**2+snn(3,2,is)**2)
         snngpn=sqrt(snngp(1,2,is)**2+snngp(2,2,is)**2+snngp(3,2,is)**2)
         if (snnn.lt.1d-9*sqrt(snngpn)) then
            snn(0,2,is)=0d0
            snn(1,2,is)=snngp(1,2,is)/snngpn
            snn(2,2,is)=snngp(2,2,is)/snngpn
            snn(3,2,is)=snngp(3,2,is)/snngpn
         else
            sense=snn(1,2,is)*snngp(1,2,is)+snn(2,2,is)*snngp(2,2,is)
     1                                     +snn(3,2,is)*snngp(3,2,is)
            if (sense.ge.0d0) then
               snn(0,2,is)=snnn
               snn(1,2,is)=snn(1,2,is)/snnn
               snn(2,2,is)=snn(2,2,is)/snnn
               snn(3,2,is)=snn(3,2,is)/snnn
            else
               snn(0,2,is)=-snnn
               snn(1,2,is)=-snn(1,2,is)/snnn
               snn(2,2,is)=-snn(2,2,is)/snnn
               snn(3,2,is)=-snn(3,2,is)/snnn
            endif
         endif
      enddo
      return
      end subroutine gosnn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine govnn(nxyz)
! 
! this subroutine computes the volume gradient at each node
!
!   vnn(0:3,3,NSM) norm of gradV-x-y-z; node level; node stack
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8) dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3
      real(8) fac1x,fac2x,fac3x,fac1y,fac2y,fac3y,fac1z,fac2z,fac3z
      real(8) dVx,dVy,dVz,volgp,vnorm,vepsilon
      integer iq,isg,lvg,isn,lvn,istack 
!
      vnn=0d0
      voln=0d0
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over Gauss points
            do lvg=1,3
               dx1=dxV(1,lvg,isg,iq) !derivative of x wrt to xi @GP
               dx2=dxV(2,lvg,isg,iq)
               dx3=dxV(3,lvg,isg,iq)
               dy1=dyV(1,lvg,isg,iq)
               dy2=dyV(2,lvg,isg,iq)
               dy3=dyV(3,lvg,isg,iq)
               dz1=dzV(1,lvg,isg,iq)
               dz2=dzV(2,lvg,isg,iq)
               dz3=dzV(3,lvg,isg,iq)
               fac1x=(+dy2*dz3-dy3*dz2)*wgp(lvg)
               fac2x=(-dy1*dz3+dy3*dz1)*wgp(lvg)
               fac3x=(+dy1*dz2-dy2*dz1)*wgp(lvg)
               fac1y=(-dx2*dz3+dx3*dz2)*wgp(lvg)
               fac2y=(+dx1*dz3-dx3*dz1)*wgp(lvg)
               fac3y=(-dx1*dz2+dx2*dz1)*wgp(lvg)
               fac1z=(+dx2*dy3-dx3*dy2)*wgp(lvg)
               fac2z=(-dx1*dy3+dx3*dy1)*wgp(lvg)
               fac3z=(+dx1*dy2-dx2*dy1)*wgp(lvg)
               volgp=det(lvg,isg,iq)*wgp(lvg)
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
                     voln(lvn,istack)=voln(lvn,istack)+
     1                                H(lvn,isn,lvg,isg)*volgp
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--normalize
!
      do isn=1,ns
         do lvn=1,3
            vnorm=sqrt(vnn(1,lvn,isn)**2+vnn(2,lvn,isn)**2
     1                                  +vnn(3,lvn,isn)**2)
            vnn(0,lvn,isn)=vnorm
            vepsilon=1d-9*abs(voln(lvn,isn))**(2./3.)
            if (voln(lvn,isn).le.0d0) then
               print *,'govnn warning!: voln < 0'
               print *,'lvn,isn,voln',lvn,isn,voln(lvn,isn)
            endif
            vnn(1,lvn,isn)=vnn(1,lvn,isn)/(vnorm+vepsilon)
            vnn(2,lvn,isn)=vnn(2,lvn,isn)/(vnorm+vepsilon)
            vnn(3,lvn,isn)=vnn(3,lvn,isn)/(vnorm+vepsilon)
         enddo
      enddo 
!
      return
      end subroutine govnn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gocnn(nxyz)
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8) xm,ym,xp,yp,nx,ny,xn,yn
      integer il,is,ism,isp
! 
! this subroutine computes the contact line normal at each node
!
      do il=1,nl
         is=isol(1,il)
         ism=isol(1,ilol(1,il))
         isp=isol(2,il)
         xm=hvec(1,1,ism)
         ym=hvec(2,1,ism)
         xp=hvec(1,1,isp)
         yp=hvec(2,1,isp)
         nx=yp-ym
         ny=xm-xp
         cnn(1,il)=nx/sqrt(nx**2+ny**2) !xnormal at isol(1,il)
         cnn(2,il)=ny/sqrt(nx**2+ny**2) !ynormal at isol(1,il)
         xn=hvec(1,1,is)
         yn=hvec(2,1,is)
         cnn(0,il)=sqrt((xn-xp)**2+(yn-yp)**2) !length of contact line il
      enddo
!
      return
      end subroutine gocnn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine govtg(nxyz)
!
!  for every surface nodes, this subroutine computes a pair
!  of orthogonal volume tangent vectors 
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      integer isn,lvn
      real(8) tgx, tgy, tgz,dtgnorm
!
      do isn=1,ns
         if (ilos(1,isn).eq.0) then!not an edge stack
            do lvn=1,3,2 !dorsal and surface nodes
! pick which one of the x and y axis most perpendicular to
               if (abs(vnn(1,lvn,isn)).lt.abs(vnn(2,lvn,isn))) then
! vtg1 = vnn x e_x
                  tgy=vnn(3,lvn,isn)
                  tgz=-vnn(2,lvn,isn)
                  dtgnorm=1d0/(sqrt(tgy*tgy+tgz*tgz))
                  vtg1(1,lvn,isn)=0d0
                  vtg1(2,lvn,isn)=tgy*dtgnorm
                  vtg1(3,lvn,isn)=tgz*dtgnorm
               else
! vtg1 = vnn x e_y
                  tgx=-vnn(3,lvn,isn)
                  tgz=vnn(1,lvn,isn)
                  dtgnorm=1d0/(sqrt(tgx*tgx+tgz*tgz))
                  vtg1(1,lvn,isn)=tgx*dtgnorm
                  vtg1(2,lvn,isn)=0d0
                  vtg1(3,lvn,isn)=tgz*dtgnorm
               endif
! vtg2 = vnn x utg1
               vtg2(1,lvn,isn)=+vnn(2,lvn,isn)*vtg1(3,lvn,isn)
     1                         -vnn(3,lvn,isn)*vtg1(2,lvn,isn)
               vtg2(2,lvn,isn)=-vnn(1,lvn,isn)*vtg1(3,lvn,isn)
     1                         +vnn(3,lvn,isn)*vtg1(1,lvn,isn)
               vtg2(3,lvn,isn)=+vnn(1,lvn,isn)*vtg1(2,lvn,isn)
     1                         -vnn(2,lvn,isn)*vtg1(1,lvn,isn)
            enddo
         else !it is an edge stack
            do lvn=1,3 !all three nodes on boundary
               tgx=vnn(2,lvn,isn)
               tgy=-vnn(1,lvn,isn)
               dtgnorm=1d0/(sqrt(tgx*tgx+tgy*tgy))
! vtg1 = vnn x e_z
               vtg1(1,lvn,isn)=tgx*dtgnorm
               vtg1(2,lvn,isn)=tgy*dtgnorm
               vtg1(3,lvn,isn)=0d0
            enddo
         endif
      enddo
!
      return
      end subroutine govtg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine goca(nxyz)
!
! This subroutine computes the contact angle for each contact line node
!
      implicit none
!
      real(8) nxyz(3,3,NSM) 
!             xyz positions of nodes at levels 1,2,3; stacks 1-NSM
!
      real(8),parameter::pi=3.141592d0
      real(8) vnodes(3,3)  !v at nodes (component,level)
      real(8) tg(3) !tangent vector
      real(8) v(3),  d2vdzeta2(3) 
      real(8) tgdotcnn,tgxy,tgz
      integer il,is 
!
      do il=1,nl
         is=isol(1,il)
         vnodes(:,:)=nxyz(:,:,is)
         call getv3(vnodes,-1d0,v,tg,d2vdzeta2)
!--dot product with outward normal
         tgdotcnn=tg(1)*cnn(1,il)+tg(2)*cnn(2,il)
         tgxy=sqrt(tg(1)**2+tg(2)**2)
         tgz=tg(3)
         if (tgdotcnn.ge.0d0) then
            if (tgz.ge.0d0) then
               ca(il)=0.5d0*pi+atan(tgxy/(tgz+vtiny))
            else
               ca(il)=pi+atan(abs(tgz)/(tgxy+vtiny))
            endif
         else
            if (tgz.ge.0d0) then
               ca(il)=atan(tgz/(tgxy+vtiny))
            else
               ca(il)=-atan(abs(tgz)/(tgxy+vtiny))
            endif
         endif
      enddo
!
      return
      end subroutine goca               
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine govolint(cnode,volint)
!
!--the routine computes the volume integral volint of scalar field cnode
!
      implicit none
      real(8) cnode(3,NSM),volint
      real(8) fac
      integer iq,isg,lvg,isn,ksn,lvn
!
      volint=0.0d0
!--loop of over elements
      do iq=1,nq
!--loop over 12 Gauss points
         do isg=1,4
            do lvg=1,3
!--fac=Jacobian * weight of GP
               fac=det(lvg,isg,iq)*wgp(lvg)
               do isn=1,4
                  ksn=isoq(isn,iq)
                  do lvn=1,3
                     volint=volint+fac*cnode(lvn,ksn)*H(lvn,isn,lvg,isg)
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine govolint
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosurfintn(cnode,surfintv,surfintd,surfinte)
!
!--this subroutine computes the dorsal, ventral, and edge surface 
!   integral of a node field cnode.
!   
      real(8) cnode(3,NSM)
      real(8) surfintd,surfintv,surfinte
      real(8) dAdn, dAvn, dAen
!
      surfintd=0.0d0
      surfintv=0.0d0
!--loop over elements for dorsal and ventral surfaces
      do iq=1,nq
!--loop over (surface) Gauss points
         do isg=1,4
!--compute elementary surface area for each GP
            dAvn=dAv(0,isg,iq)
            dAdn=dAd(0,isg,iq)
!--loop over 4 stacks of element; 1 is ventral, 3 is dorsal
            do isn=1,4
               ksn=isoq(isn,iq)
               surfintv=surfintv+cnode(1,ksn)*S4(isn,isg)*dAvn
               surfintd=surfintd+cnode(3,ksn)*S4(isn,isg)*dAdn
            enddo
         enddo
      enddo
!--loop over edges
      surfinte=0.0d0
      do il=1,nl
!--loop over edge GP
         do isg=1,2
            do lvg=1,3
!--compute weighted elementary area for each GP
               dAen=wgp(lvg)*dAe(0,lvg,isg,il)
!--loop over edge nodes for each contribution
               do isn=1,2
                  ksn=isol(isn,il)
                  do lvn=1,3
                     surfinte=surfinte+
     1                        cnode(lvn,ksn)*S6(lvn,isn,lvg,isg)*dAen
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine gosurfintn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosurfintd(csurf,surfint)
!  this subroutine computes the surface integral of a dorsal surface field
      real(8) csurf(*),surfint
!
      surfint=0.0d0
      do iq=1,nq
         do isg=1,4
            surfint=surfint+csurf(iq)*dAd(0,isg,iq)
         enddo
      enddo
      return
      end subroutine gosurfintd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosurfintv(csurf,surfint)
!  this subroutine computes the surface integral of a ventral surface field
      real(8) csurf(*),surfint
!
      surfint=0.0d0
      do iq=1,nq
         do isg=1,4
            surfint=surfint+csurf(iq)*dAv(0,isg,iq)
         enddo
      enddo
      return
      end subroutine gosurfintv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine gosurfinte(csurf,surfint)
!  this subroutine computes the surface integral of an edge surface field
      real(8) csurf(*),surfint
!
      surfint=0d0
!--loop over edges
      do il=1,nl
!--loop over two GP stacks
         do isg=1,2
!--loop over three GP levels
            do lvg=1,3
               areae=areae+wgp(lvg)*csurf(il)*dAe(0,lvg,isg,il)
            enddo
         enddo
      enddo
      return
      end subroutine gosurfinte
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine oldread(idiv,lvec,nvec)
!
!  this subroutine reads the old standard 2D RIF format
!  and returns nd,nq,ns,nl, names and dscp as in shallow water
!  but also nvec (->hvec and ->svec) and lvec (->evec). 
!
      real(8) nvec(12,NSM)! vertex-node vectors (N_j)
      real(8) lvec(12,NLM)! line-node  vectors.
      character (len=66) label
      character (len=11):: score = '***** *****'
!
      nvec=0d0
      lvec=0d0
!
      read(idiv,200)ifrm,mshtitle!read the frame# and title
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nd,label! read ND
      read(idiv,210)score,score,score,score,score,score
      do i=1,nd!
         read(idiv,200)id,heads(i)!heading for the i-th descriptor
         read(idiv,210)(names(k,i),k=01,06)
         read(idiv,220)( dscp(k,i),k=01,06)
         read(idiv,210)(names(k,i),k=07,12)
         read(idiv,220)( dscp(k,i),k=07,12)
      enddo
!
! read POINTERS of EDGE TOPOLOGY
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nl,label!READ NL = number of EDGEs
      read(idiv,210)score,score,score,score,score,score
      do i=1,nl,3
         read(idiv,230)(isol(1,j),isol(2,j),ibol(j),j=i,i+2)
      enddo
!
!**** read POINTERS of QUADRILATERIAL TOPOLOGY
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nq,label!read nq=number of quadrilaterials
      read(idiv,210)score,score,score,score,score,score
      do i=1,nq,2
         read(idiv,235)((isoq(k,j),k=1,4),j=i,i+1)
      enddo 
!
!**** read DATA DEFINED at EDGE NODES
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)nl,label!'EDGE NODE DATA VECTORS '
      read(idiv,210)score,score,score,score,score,score
      do i=1,nl
         read(idiv,240)(lvec(k,i),k=01,04)
         read(idiv,250)(lvec(k,i),k=05,08)
         read(idiv,250)(lvec(k,i),k=09,12)
      enddo
!
!**** read DATA DEFINED ON VERTEX NODES
      read(idiv,210)score,score,score,score,score,score
      read(idiv,200)ns,label!number of vertex nodes
      read(idiv,210)score,score,score,score,score,score
      do i=1,ns
         read(idiv,240)  (nvec(k,i),k=01,04)
         read(idiv,250)  (nvec(k,i),k=05,08)
         read(idiv,250)  (nvec(k,i),k=09,12)
      enddo
      call iddscp !id the descriptors
!
!**** FORMAT STATEMENTS
  200 format(1x,i4,1x,a66)!formate for io of title+headings
  210 format(6(1x,a11))!format for names
  220 format(6(1x,1pe11.4))!format for dscp values
  230 format(3(1x,4x,1x,3(1x,i5)))!format;for L pointers
  235 format(2(1x,5x,5x,1x,4(1x,i5)))!format;for Q pointers
  240 format(5x,1x,2x,4(1x,1pe15.7))!1st line of lvec and nvec
  250 format(5x,1x,2x,4(1x,1pe15.7))!2nd line of lvec and nvec
      return
      end subroutine oldread
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine oldwrite(ifrm,idiv)
!
!  this subroutine writes the old standard 2D RIF format
!  using the bottom of the mesh
!  but also nvec (->hvec and ->svec) and lvec (->evec). 
!
      real(8) nvec(12,NSM)! vertex-node vectors (N_j)
      real(8) lvec(12,NLM)! line-node  vectors.
      character (len=11):: score= '***** *****'
!
      lvec=0d0
      nvec=0d0
      nvec(1:2,1:ns)=hvec(1:2,1,1:ns)
!
      write(idiv,200)ifrm,mshtitle!write the frame# and title
      dscp(1,idfrm)=time_
      dscp(2,idfrm)=tstp
!
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nd,'CONTROL PARAMETERS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nd
         write(idiv,200)i,heads(i)!heading for the i-th discriptor
         write(idiv,210)(names(k,i),k=01,06)
         write(idiv,220)( dscp(k,i),k=01,06)
         write(idiv,210)(names(k,i),k=07,12)
         write(idiv,220)( dscp(k,i),k=07,12)
      enddo
!
!**** write POINTERS of EDGE TOPOLOGY
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nl,'POINTERS OF EDGE TOPOLOGY '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nl,3
         write(idiv,230)(j,'=',isol(1,j),isol(2,j),ibol(j),j=i,i+2)
      enddo
!
!**** write POINTERS of QUADRILATERIAL TOPOLOGY
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nq,'POINTERS OF QUADRILATERIAL TOPOLOGY '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nq,2
         write(idiv,235)('*****',j,'=',(isoq(k,j),k=1,4),j=i,i+1)
      enddo 
!
!**** write DATA DEFINED at EDGE NODES
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)nl,'EDGE NODE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,nl
         write(idiv,240)i,'=',(lvec(k,i),k=01,04)
         write(idiv,250)      (lvec(k,i),k=05,08)
         write(idiv,250)      (lvec(k,i),k=09,12)
      enddo
!
!**** write DATA DEFINED ON VERTEX NODES
      write(idiv,210)score,score,score,score,score,score
      write(idiv,200)ns,'VERTEX NODE DATA VECTORS '
      write(idiv,210)score,score,score,score,score,score
      do i=1,ns
         write(idiv,240)i,'=',(nvec(k,i),k=01,04)
         write(idiv,250)      (nvec(k,i),k=05,08)
         write(idiv,250)      (nvec(k,i),k=09,12)
      enddo 
      write(idiv,210)score,score,score,score,score,score
      call iobuffer(10,idiv)!flush the io buffer
!
!**** FORMAT STATEMENTS
  200 format(1x,i4,1x,a66)!formate for io of title+headings
  210 format(6(1x,a11))!format for names
  220 format(6(1x,1pe11.4))!format for dscp values
  230 format(3(1x,i4,a1,3(1x,i5)))!format;for L pointers
  235 format(2(1x,a5,i5,a1,4(1x,i5)))!format;for Q pointers
  240 format(i5,a1,2x,4(1x,1pe15.7))!1st line of lvec and nvec
  250 format(5x,1x,2x,4(1x,1pe15.7))!1st line of lvec and nvec
!
      return
      end subroutine oldwrite!(ifrm,idiv)

      END MODULE IOLIBSW!this ends the list of contained subprograms
