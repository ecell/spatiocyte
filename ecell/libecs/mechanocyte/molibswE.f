      MODULE MOLIBSW
      USE ISO_C_BINDING
!
! library of momentum-solver routines
!
! General idea:
!
! 1) Solve bulk cytoplasmic momentum equation for v_n: "velocity eq."
!
!      div [nu(grad v_n + T(grad v_n))] = grad P + div Psi
!
!    or in FEM form:
!                 Qv = sum nu<dHjx*dHky>|iq
!                 uldq = Psi loads
!                 vldq = uldq + Pressure loads
!    so that
!                 Qv * vnw = vldq  gives vnw
!
!
! 2) Use computed v_n to solve the solvent 
!    momentum eq. + incompressibility for pressure P: "pressure eq."
!
!        div [(grad P) / phi] = div v_n
!
! for numerical stability, we use:
!
!       div [v_n - (grad P_new)/phi] = rlax*(P_new-P_old)
!
!        where rlax ~ (div v_n)/P
!
!      in FEM form:
!        Qp = (1/phi)*<dHjx*dHkx>|iq + rlax*Hj*Hk
!        pldv = Pold * Hj +dHj . v_n
!
!  We iterate between 1) and 2) until we have converged.
!
!
      USE IOLIBSW
      USE MXLIBSW
!
!-- velocity equation arrays
      real(8) Qv(3,3,4,3,3,4,NSM) !element velocity stiffness matrix
!         vj_xyz; j level; j stack; vk_xyz; k level; k stack; iq element
      real(8) Lvv(3,4,3,4,NSM), Ldv(3,4,3,4,NSM)!boundary stiffness
      real(8) Lev(3,3,2,3,3,2,NLM)!edge boundary velocity stiffness
      real(8) Cv(3,3,3*NSM) !preconditioner for velocity stiffness matrix
      real(8) uldq(3,3,4,NSM)!non pressure element velocity loads
      real(8) vldq(3,3,4,NSM)!element velocity loads (include pressure)
      real(8) vldv(3,4,NSM), vldd(3,4,NSM)!boundary velocity loads
      real(8) vlde(3,3,2,NLM)!edge boundary velocity loads
      real(8) vnw(3,3*NSM)!network velocity
!
!--pressure equation
      real(8) Qp(3,4,3,4,NSM) !pressure stiffness matrix
!             j level; j stack; k level; k stack; iq element
      real(8) Lvp(4,4,NSM),Ldp(4,4,NSM)! boundary pressure stiffness
      real(8) Lep(3,2,3,2,NLM)! edge boundary pressure stiffness
      real(8) Cp(NSM*3) !preconditioner for pressure stiffness matrix
      real(8) pldq(3,4,NSM) !element pressure loads
      real(8) pldv(4,NSM), pldd(4,NSM) !boundary pressure loads
      real(8) plde(3,2,NLM)!edge boundary pressure loads
      real(8) pon(3*NSM)!pressure
!
      real(8) rlv(3,3,NSM)!relative velocity
      real(8) ftn(3,3,NSM)!loads due to surface tension
!
      contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine modriver(icyc,epsl,idebug) bind(C)
!
!  when called, this routine iterates between the velocity
!  equation and the pressure equation until changes in pressure
!  are small, and we have converged to a solution for the network
!  velocity.
!
      implicit none
!
      integer(C_INT)::icyc!number of iterations
      real(c_double)::epsl!tolerance for pdel (p change during iteration)
      integer(C_INT)::idebug!debugging flag
!
      integer icycmx!max iterations
      parameter(icycmx=200)
      integer ipanic,kpanic! panic flag, counter
      real(8) cpen!penalty coefficient used for boundary vel. cond.
      real(8) rlax!relaxation coefficient for pressure diffusion
      real(8),save::rlax_old=0d0
      real(8) divvsum,divvav1,divvav0!average div v estimates
      real(8) vdel, pdel !relative velocity and pressure change
      real(8) vdel_old, pdel_old !relative velocity and pressure change
      real rnum(NSM*3) !random number array for pressure perturbation
      integer nccg!max conjugate gradient iterations
      parameter(nccg=50)
      integer inode,isn,ilv,k
      real(8) ppert,pmax,pmin,vmax,vmin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!             PRELIMINARY VELOCITY EQUATION COMPUTATIONS
!
      if (idvis.le.0) then
         write(*,*)'MODRIVER ERROR; idvis=0'
         stop
      endif
!--set penalty coeff. used for vel. boundary conditions.
      call mocstpen(cpen)
!
!--initialize stiffness matrices and load vectors
      Qv=0d0; Lvv=0d0; Ldv=0d0; Lev=0d0 !velocity stiffness
      uldq=0d0; vldv=0d0; vldd=0d0; vlde=0d0 !velocity loads
!
!--assemble element velocity stiffness matrix (viscous term)
      call movisvmx(Qv,cpen) !shear contribution
      call modilvmx(Qv,cpen) !"pure dilation" contribution
!
!--assemble boundary velocity stiffness matrices
      if(idvfx.gt.0) call mofixvmx(Lvv,Ldv,Lev,cpen)!vel. constraints
      if(iddrg.gt.0) call modrgvmx(Lvv,Ldv,Lev,cpen)!bdy. vel. drags 
!
!--compute velocity conditioning matrix, Cv 
      call cgvcndmx(Cv,Qv,Lvv,Ldv,Lev)
!
!--assemble element velocity load vector
      if(idpsi.gt.0)call moctrvld(uldq)!apply contractility loads
      if(idbfr.gt.0)call mobodvld(uldq)!body force loads
      if(idsfr.gt.0)call mosfrvld(uldq)!apply surface force loads
      if(idsfr.lt.0)call mosfrvld2(uldq)!apply surface force loads
!
!--assemble boundary velocity load vectors
      if(idvfx.gt.0)call moxvnvld(vldv,vldd,vlde,Lvv,Ldv,Lev)!ext. vel. 
      if(idtrc.gt.0)call motrcvld(vldv,vldd,vlde)!bdy traction loads
!-surface tension velocity vector loads:
      if(idgam.gt.0)call mostnvld_mrc(vldv,vldd,vlde)!Marc Herant method
!     if(idgam.gt.0)call mostnvld_mxd(vldv,vldd,vlde)!Micah Dembo method
!
      if (idebug.ge.2) then 
         print *,'modriver:cpen',cpen
         print *,'modriver: idvfx',idvfx
         print *,'modriver: iddrg',iddrg
         print *,'modriver: idpsi',idpsi
         print *,'modriver: idbfr',idbfr
         print *,'modriver: idsfr',idsfr
         print *,'modriver: idtrc',idtrc
         print *,'modriver: idgam',idgam
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         begin calculation of rlax 
!
! we perturb the pressure and examine the resulting change in div v_n
!           rlax=delta(div v_n)/delta(P)
!
      ipanic=0  
   10 continue !this is where we restart in case of a PANIC
      vdel_old=0d0
      pdel_old=0d0
      kpanic=0
!
      pmax=maxval(hvec(10,1:3,1:ns))
      pmin=minval(hvec(10,1:3,1:ns))
      if (idebug.ge.2) print *,'prelim. pmax,pmin',pmax,pmin
!--1% perturbation of the pressure or 1d-5 whichever is more
      ppert=1d-5+(abs(pmax)+abs(pmin))*1d-2
      call random_seed()
      call random_number(rnum)
!--load pressures and network velocities, and pertub pressure
      do isn=1,ns
         do ilv=1,3
            inode=(isn-1)*3+ilv
            pon(inode)=hvec(10,ilv,isn)+ppert*rnum(inode)
            vnw(1,inode)=hvec(7,ilv,isn)
            vnw(2,inode)=hvec(8,ilv,isn)
            vnw(3,inode)=hvec(9,ilv,isn)
         enddo
      enddo
!--put non-pressure velocity loads uldq in vldq 
!                    and add the pressure loads
      vldq=uldq
      call moprsvld(vldq,pon)
!
!--solve for network velocities
!   conjugate gradient solver:
      call cgvsolve(vnw,Qv,Lvv,Ldv,Lev,Cv,vldq,vldv,vldd,vlde,2*nccg,k)
      if (idebug.ge.1) print *,'back from cgvsolve, k=',k
!   LU decomposition solver:
!     call luvsolve(vnw,Qv,Lvv,Ldv,Lev,vldq,vldv,vldd,vlde)
!
!--get velocity divergence
      call modivv(divvsum,divvav1,vnw)
!     print *,'divvsum,divvav1',divvsum,divvav1
!
!--reload unperturbed pressures and original velocities
      do isn=1,ns
         do ilv=1,3
            inode=(isn-1)*3+ilv
            pon(inode)=hvec(10,ilv,isn)
!           print *,'vnwxyz',vnw(1,inode),vnw(2,inode),vnw(3,inode)
            vnw(1,inode)=hvec(7,ilv,isn)
            vnw(2,inode)=hvec(8,ilv,isn)
            vnw(3,inode)=hvec(9,ilv,isn)
         enddo
      enddo
!--put non-pressure velocity loads uldq in vldq 
!                    and add the pressure loads
      vldq=uldq
      call moprsvld(vldq,pon)
!
!--solve for network velocities: 
!          (Qv+Lvv+Ldv+Lev)*vnw = vldq+vldv+vldd+vlde
!   conjugate gradient solver:
      call cgvsolve(vnw,Qv,Lvv,Ldv,Lev,Cv,vldq,vldv,vldd,vlde,2*nccg,k)
      if (idebug.ge.2) print *,'back from cgvsolve, k=',k
!   LU decomposition solver:
!     call luvsolve(vnw,Qv,Lvv,Ldv,Lev,vldq,vldv,vldd,vlde)
!
!--get velocity divergence
      call modivv(divvsum,divvav0,vnw)
!     print *,'divvsum,divvav0',divvsum,divvav0
!
!--get rlax from the change in velocity divergences 
!                    from the pressure pertubations 
      if (idphi.eq.0) then !Stokes flow, conservative rlax
         rlax=1d3*(1d0+ipanic)*abs(divvav1-divvav0)/ppert
      else
         rlax=1d2*(1d0+ipanic)*abs(divvav1-divvav0)/ppert
      endif
!--rlax cannot decrease to less than 70% of previous rlax
      rlax=max(rlax,rlax_old*0.7d0)
      rlax_old=rlax
      if (idebug.ge.1) then
         print *,'divvav1,divvav0,ppert',divvav1,divvav0,ppert
         print *,'rlax=',rlax
      endif
!
!               end of rlax computations
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!              PRELIMINARY PRESSURE EQUATION COMPUTATIONS
!
!--do some checking on the inverse permeability (resistance) phi
!  we set the max for Stokes limit.
      call mophichk
!
!--initialize pressure stiffness matrices 
      Qp=0d0; Lvp=0d0; Ldp=0d0; Lep=0d0 
!
!--assemble stiffness matrices for pressure equation
      call mophipmx(Qp) !pressure diffusion matrix
      call momaspmx(Qp,rlax)!pressure mass matrix contribution
!--boundary hydraulic conductivity
      if(idhyc.gt.0)call mohycpmx(Lvp,Ldp,Lep)
      if (idebug.ge.2) then
         print *,'modriver: idphi=',idphi
         print *,'modriver: idhyc=',idhyc
      endif
!
!--compute Pressure preconditioning matrix,
      call cgscndmx(Cp,Qp,Lvp,Ldp,Lep)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ITERATIVE LOOP BETWEEN PRESSURE AND VELOCITY EQUATIONS BEGINS HERE
!
!  the only FEM arrays that need to be  recomputated during 
!  the iterations are:
!       the pressure loads
!       the pressure contribution to the velocity load
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      icyc=0!initialize pressure diffusion cycle counter
  200 icyc=icyc+1!advance the PD-cycle counter
!

!--initialize pressure load vectors
      pldq=0d0; pldv=0d0; pldd=0d0; plde=0d0 !pressure loads
!
!--assemble pressure loads
      call momaspld(pldq,pon,rlax)!add mass loads 
      call modivpld(pldq,vnw)!add velocity divergence load
!--osmotic pressure load
      if(idhyc.gt.0)call mohycpld(pldv,pldd,plde)
!
!--compute new pressure field:
!            (Qp+Lvp+Ldp+Lep)*pon = pldq+pldv+pldd+plde
!--conjugate gradient method
      call cgssolve(pon,Qp,Lvp,Ldp,Lep,Cp,pldq,pldv,pldd,plde,nccg,k)
      if (idebug.ge.2) print *,'back from cgssolve, k=',k
!--LU decomposition method
!     call lussolve(pon,Qp,Lvp,Ldp,Lep,pldq,pldv,pldd,plde)
!
!--check for trouble
      pmin=minval(pon(1:3*ns))
      pmax=maxval(pon(1:3*ns))
      if (pmin.lt.-big.or.pmax.gt.big) then
         print *,'modriver PANIC in pressure computation: pmin,pmax',
     1            pmin,pmax
         call mopanic(ipanic)
         goto 10
      endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  UPDATE VELOCITY EQUATION:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--put non-pressure velocity loads uldq in vldq 
!                    and add the pressure loads
      vldq=uldq 
      call moprsvld(vldq,pon)

!--solve for network velocities: 
!--conjugate gradient method
      call cgvsolve(vnw,Qv,Lvv,Ldv,Lev,Cv,vldq,vldv,vldd,vlde,2*nccg,k)
      if (idebug.ge.2) print *,'back from cgvsolve, k=',k
!--LU decomposition method
!     call luvsolve(vnw,Qv,Lvv,Ldv,Lev,vldq,vldv,vldd,vlde)
!
!--check for trouble
      vmin=minval(vnw(:,1:3*ns))
      vmax=maxval(vnw(:,1:3*ns))
      if (pmin.lt.-big.or.pmax.gt.big) then
         print *,'modriver PANIC in velocity computation: pmin,pmax',
     1            pmin,pmax
         call mopanic(ipanic)
         goto 10
      endif
!
      if (idebug.ge.2) then
         call modivv(divvsum,divvav1,vnw)
         pmin=minval(pon(1:3*ns))
         pmax=maxval(pon(1:3*ns))
         print *,'modriver: divvsum,pmin,pmax',divvsum,pmin,pmax
      endif
!
!--check convergence 
      call mochk(pdel,pon,vdel,vnw)
      if (idebug.ge.1) print *,'modriver: pdel,vdel',pdel,vdel,icyc
!
! we require better than 1% progress in velocity or pressure or 
! we increment the panic counter
      if ((((pdel_old-pdel)/pdel).lt.0.001).and.
     1    (((vdel_old-vdel)/vdel).lt.0.001)) then
         kpanic=kpanic+1
         if (idebug.ge.1.and.kpanic.ne.1) print*,'kpanic=',kpanic
      else
         kpanic=0
      endif
      if (kpanic.gt.10) then
         print *,'modriver PANIC due to bad convergence'
         print *,'pdel,pdel_old',pdel,pdel_old
         print *,'vdel,vdel_old',vdel,vdel_old
         call mopanic(ipanic)
         goto 10
      endif
      pdel_old=pdel
      vdel_old=vdel
         
!
      if(icyc.lt.5)go to 200!icyc too small 
!--if we have nor converged, start new iteration except if icyc>icycmx
      if((pdel.gt.epsl.or.vdel.gt.epsl).and.icyc.le.icycmx) go to 200
!
      if(icyc.gt.icycmx)write(*,300)icyc,nccg,pdel,rlax!convergence failure
!
!   END OF ITERATIONS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--begin cleanup phase by solving for rlv=-GRAD(pon)/phi.
      call morlvpol(rlv,pon)
!
! update hvec's with newly computed fields.
      do isn=1,ns
         do ilv=1,3
            inode=(isn-1)*3+ilv
            hvec(4,ilv,isn)=rlv(1,ilv,isn)
            hvec(5,ilv,isn)=rlv(2,ilv,isn)
            hvec(6,ilv,isn)=rlv(3,ilv,isn)
            hvec(7,ilv,isn)=vnw(1,inode)
            hvec(8,ilv,isn)=vnw(2,inode)
            hvec(9,ilv,isn)=vnw(3,inode)
            hvec(10,ilv,isn)=pon(inode)
         enddo
      enddo
!
      return
!
!--format statements
  300 format('icyc=',i4,' nccg=',i4,' |pdel|=',1pe10.3,' rlax=',1pe10.3)
!
      end subroutine modriver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   BEGIN VELOCITY EQUATION SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mocstpen(cpen)
!
! this subroutine computes the penalty parameter cpen (dimension 
! viscosity)that will be used to enforce dirichlet boundary 
! conditions for the velocity. 
!
      implicit none
!     
      real(8) cpen
!
      real(8) visav! average viscosity through the mesh
      integer istack, ilv
!
      visav=0.0
      do istack=1,ns
         do ilv=1,3
            visav=visav+vis(ilv,istack)
         enddo
      enddo
      visav=visav/(3*ns)
!
      cpen=(1.0e+07)*visav
!
      return
      end subroutine mocstpen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Begin subroutines for velocity stiffness matrices 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine movisvmx(Qv,cpen)
!
!  Assemble Qv, the FEM stiffness matrix operator for
!    div [vis (grad(*)+T(grad(*)))]
!
!  USER-SPECIFIED (constitutive law):
!  vis(level,stack)=coefficent of "shear" viscosity for each node
!
! Note that the viscosity must not be so large as to make the 
! stiffness matrix ill-conditioned. For this reason, vis=min(vis,cpen).
!
      implicit none
!
      real(8) Qv(3,3,4,3,3,4,NSM)
!         vj_xyz; j level; j stack; vk_xyz; k level; k stack; iq element
      real(8) cpen! the penalty coefficient
!
      real(8) visgp, vfac, vfacx,vfacy,vfacz 
      real(8) delxx,delxy,delxz,delyx,delyy,delyz,delzx,delzy,delzz
      integer iq,isg,lvg,isn,istack,lvn,ks,kl,kx,js,jl,jx
!
!--loop over elements
      do iq=1,nq
!--loop over GPs stacks and levels
         do isg=1,4
            do lvg=1,3
!--compute viscosity at GP
               visgp=0d0
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     visgp=visgp+H(lvn,isn,lvg,isg)*vis(lvn,istack)
                  enddo
               enddo
               vfac=min(visgp,cpen)*det(lvg,isg,iq)*wgp(lvg)
!--loop over k and j vertices by stack, level, and component
               do ks=1,4
                  do kl=1,3
                     vfacx=vfac*dHg(1,kl,ks,lvg,isg,iq)
                     vfacy=vfac*dHg(2,kl,ks,lvg,isg,iq)
                     vfacz=vfac*dHg(3,kl,ks,lvg,isg,iq)
                     do js=1,4
                        do jl=1,3
!--compute all the possible dH_j*dH_k combinations
                           delxx=dHg(1,jl,js,lvg,isg,iq)*vfacx
                           delxy=dHg(1,jl,js,lvg,isg,iq)*vfacy
                           delxz=dHg(1,jl,js,lvg,isg,iq)*vfacz
                           delyx=dHg(2,jl,js,lvg,isg,iq)*vfacx
                           delyy=dHg(2,jl,js,lvg,isg,iq)*vfacy
                           delyz=dHg(2,jl,js,lvg,isg,iq)*vfacz
                           delzx=dHg(3,jl,js,lvg,isg,iq)*vfacx
                           delzy=dHg(3,jl,js,lvg,isg,iq)*vfacy
                           delzz=dHg(3,jl,js,lvg,isg,iq)*vfacz
                 Qv(1,jl,js,1,kl,ks,iq)=Qv(1,jl,js,1,kl,ks,iq)
     1                                        +2d0*delxx+delyy+delzz
                 Qv(2,jl,js,1,kl,ks,iq)=Qv(2,jl,js,1,kl,ks,iq)+delyx
                 Qv(3,jl,js,1,kl,ks,iq)=Qv(3,jl,js,1,kl,ks,iq)+delzx
                 Qv(1,jl,js,2,kl,ks,iq)=Qv(1,jl,js,2,kl,ks,iq)+delxy
                 Qv(2,jl,js,2,kl,ks,iq)=Qv(2,jl,js,2,kl,ks,iq)
     1                                        +delxx+2d0*delyy+delzz
                 Qv(3,jl,js,2,kl,ks,iq)=Qv(3,jl,js,2,kl,ks,iq)+delzy
                 Qv(1,jl,js,3,kl,ks,iq)=Qv(1,jl,js,3,kl,ks,iq)+delxz
                 Qv(2,jl,js,3,kl,ks,iq)=Qv(2,jl,js,3,kl,ks,iq)+delyz
                 Qv(3,jl,js,3,kl,ks,iq)=Qv(2,jl,js,2,kl,ks,iq)
     1                                        +delxx+delyy+2d0*delzz
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine movisvmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine modilvmx(Qv,cpen)
!
!  Assemble Qv, the FEM stiffness matrix operator for
!    div [lam (grad(*)]
!
!  USER-SPECIFIED (constitutive law):
!  lam(level,stack)="pure dilation" Lame viscosity coefficent at node
!
! Note that the viscosity must not be so large as to make the 
! stiffness matrix ill-conditioned. For this reason, lam=min(lam,cpen).
!
      implicit none
!
      real(8) Qv(3,3,4,3,3,4,NSM)
!         vj_xyz; j level; j stack; vk_xyz; k level; k stack; iq element
      real(8) cpen! the penalty coefficient
!
      real(8) lamgp, vfac, vfacx,vfacy,vfacz 
      real(8) delxx,delxy,delxz,delyx,delyy,delyz,delzx,delzy,delzz
      integer iq,isg,lvg,isn,istack,lvn,ks,kl,kx,js,jl,jx
!
      if (maxval(lam(1:3,1:ns)).le.vtiny) return
!
!--loop over elements
      do iq=1,nq
!--loop over GPs stacks and levels
         do isg=1,4
            do lvg=1,3
!--compute viscosity at GP
               lamgp=0d0
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     lamgp=lamgp+H(lvn,isn,lvg,isg)*lam(lvn,istack)
                  enddo
               enddo
               vfac=min(lamgp,cpen)*det(lvg,isg,iq)*wgp(lvg)
!--loop over k and j vertices by stack, level, and component
               do ks=1,4
                  do kl=1,3
                     vfacx=vfac*dHg(1,kl,ks,lvg,isg,iq)
                     vfacy=vfac*dHg(2,kl,ks,lvg,isg,iq)
                     vfacz=vfac*dHg(3,kl,ks,lvg,isg,iq)
                     do js=1,4
                        do jl=1,3
!--compute all the possible dH_j*dH_k combinations
                           delxx=dHg(1,jl,js,lvg,isg,iq)*vfacx
                           delxy=dHg(1,jl,js,lvg,isg,iq)*vfacy
                           delxz=dHg(1,jl,js,lvg,isg,iq)*vfacz
                           delyx=dHg(2,jl,js,lvg,isg,iq)*vfacx
                           delyy=dHg(2,jl,js,lvg,isg,iq)*vfacy
                           delyz=dHg(2,jl,js,lvg,isg,iq)*vfacz
                           delzx=dHg(3,jl,js,lvg,isg,iq)*vfacx
                           delzy=dHg(3,jl,js,lvg,isg,iq)*vfacy
                           delzz=dHg(3,jl,js,lvg,isg,iq)*vfacz
                 Qv(1,jl,js,1,kl,ks,iq)=Qv(1,jl,js,1,kl,ks,iq)+delxx
                 Qv(2,jl,js,1,kl,ks,iq)=Qv(2,jl,js,1,kl,ks,iq)+delyx
                 Qv(3,jl,js,1,kl,ks,iq)=Qv(3,jl,js,1,kl,ks,iq)+delzx
                 Qv(1,jl,js,2,kl,ks,iq)=Qv(1,jl,js,2,kl,ks,iq)+delxy
                 Qv(2,jl,js,2,kl,ks,iq)=Qv(2,jl,js,2,kl,ks,iq)+delyy
                 Qv(3,jl,js,2,kl,ks,iq)=Qv(3,jl,js,2,kl,ks,iq)+delzy
                 Qv(1,jl,js,3,kl,ks,iq)=Qv(1,jl,js,3,kl,ks,iq)+delxz
                 Qv(2,jl,js,3,kl,ks,iq)=Qv(2,jl,js,3,kl,ks,iq)+delyz
                 Qv(3,jl,js,3,kl,ks,iq)=Qv(2,jl,js,2,kl,ks,iq)+delzz
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine modilvmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine mofixvmx(Lvv,Ldv,Lev,cpen)
!
! This routine applies boundary velocity constraints to the 
! boundary velocity stiffness matrices
! The idea is that in the  Q * vnw = vload equation
! we introduce a diagonal Q term that's large enough to force
! the corresponding vnw to near zero. If a non zero vnw is
! desired, a matching term must be added to vload (this
! is done in subroutine moxvnvld).
!
! Note that cpen must be large enough so that it is much 
! bigger than all the other terms in the stiffness matrix
! but not so big that it makes the matrix ill-conditioned...
!
!--fixed boundary velocities (xyz component, element or edge index)
!    set to +1 for constraint, 0 for no constraint
!     real(8),dimension(3,NSM)::vfixv=0!ventral
!     real(8),dimension(3,NSM)::vfixd=0!dorsal
!     real(8),dimension(3,NLM)::vfixe=0!edge
!
! At this time, the option to fix normal and/or tangential velocities
! does not exist, and one should use the drag subroutine modrgvmx to
! achieve that.
!
      implicit none
!
      real(8) Lvv(3,4,3,4,NSM), Ldv(3,4,3,4,NSM)!ventral & dorsal
      real(8) Lev(3,3,2,3,3,2,NLM)!edge boundary stiffness matrix
      real(8) cpen! the penalty coefficient
!
      real(8) facg, facv, facd, facvg, facvv, facvd
      integer il, isg, lvg, ks, kl, kx, iq
!
!--ventral and dorsal surfaces
      !print *,'vfixv ',vfixv(1,1)
      !print *,'vfixv ',vfixv(2,2)
      do iq=1,nq 
         do isg=1,4 !loop over GPs
            facv=dAv(0,isg,iq)*cpen
            facd=dAd(0,isg,iq)*cpen
            do ks=1,4 !loop over 4 nodes
               facvv=facv*S4(ks,isg)
               facvd=facd*S4(ks,isg)
               do kx=1,3
!--because we only fix vx, vy, or vz, L is diagonal
                  Lvv(kx,ks,kx,ks,iq)=Lvv(kx,ks,kx,ks,iq)
     1                                +vfixv(kx,iq)*facvv
                  Ldv(kx,ks,kx,ks,iq)=Ldv(kx,ks,kx,ks,iq)
     1                                +vfixd(kx,iq)*facvd
               enddo
            enddo
         enddo
      enddo
!
!--edge surfaces
      do il=1,nl 
         do isg=1,2
            do lvg=1,3
               facg=dAe(0,lvg,isg,il)*wgp(lvg)*cpen
               do ks=1,2
                  do kl=1,3
                     facvg=facg*S6(kl,ks,lvg,isg)
                     do kx=1,3
!--because we only fix vx, vy, or vz, L is diagonal
             Lev(kx,kl,ks,kx,kl,ks,il)=Lev(kx,kl,ks,kx,kl,ks,il)
     1                                 +vfixe(kx,il)*facvg
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine mofixvmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine modrgvmx(Lvv,Ldv,Lev,cpen)
!
! This routine applies boundary velocity drags to the  boundary 
! velocity stiffness matrices
!
! Note that the drag must not be so large as to make the matrix
! ill conditioned. For this reason, drag=min(drag,cpen).
!
! USER-SPECIFIED (constitutive law)
! surface friction coefficients (normal/tangential,ventral dorsal edge)
!     real(8),dimension(5,NSM)::drgv=0!ventral
!     real(8),dimension(5,NSM)::drgd=0!dorsal
!     real(8),dimension(5,NLM)::drge=0!edge
!
      implicit none
!
      real(8) Lvv(3,4,3,4,NSM), Ldv(3,4,3,4,NSM)!ventral & dorsal
      real(8) Lev(3,3,2,3,3,2,NLM)!edge boundary stiffness matrix
      real(8) cpen! the penalty coefficient
!
      real(8) drgx(3),drgn,drgt,dAHk,dAHjk,snk,snj,snjk,dAen
      integer iq,il,isg,lvg,ks,kl,kx,kstack,js,jl,jx,jstack
!
!--ventral boundary
      do iq=1,nq
         drgx(1:3)=min(drgv(1:3,iq),cpen) !xyz drag
         drgn=min(drgv(4,iq),cpen) 
         drgt=min(drgv(5,iq),cpen)
         do isg=1,4 !loop over 4 GPs
            do ks=1,4
               kstack=isoq(ks,iq)
               dAHk=dAv(0,isg,iq)*S4(ks,isg)
               do kx=1,3
!--kx component of normal at ventral (level=1) node ks
                  snk=snn(kx,1,kstack) 
                  do js=1,4
                     jstack=isoq(js,iq)
                     dAHjk=dAHk*S4(js,isg)
                     Lvv(kx,js,kx,ks,iq)=Lvv(kx,js,kx,ks,iq)+
     1                                   dAHjk*drgx(kx) !xyz drag
                     do jx=1,3
!--jx component of normal at ventral (level=1) node js
                        snj=snn(jx,1,jstack)
                        snjk=snj*snk !element of the two-vector nn
!--vnormal= v.nn, vtangent= v - v.nn
                        Lvv(jx,js,kx,ks,iq)=Lvv(jx,js,kx,ks,iq)+
     1                         dAHjk*(drgt*(kd(jx,kx)-snjk)+drgn*snjk)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--dorsal boundary
      do iq=1,nq
         drgx(1:3)=min(drgd(1:3,iq),cpen) !xyz drag
         drgn=min(drgd(4,iq),cpen)
         drgt=min(drgd(5,iq),cpen)
         do isg=1,4
            do ks=1,4
               kstack=isoq(ks,iq)
               dAHk=dAd(0,isg,iq)*S4(ks,isg)
               do kx=1,3
!--kx component of normal at dorsal (level=3) node ks
                  snk=snn(kx,3,kstack)
                  do js=1,4
                     jstack=isoq(js,iq)
                     dAHjk=dAHk*S4(js,isg)
                     Ldv(kx,js,kx,ks,iq)=Ldv(kx,js,kx,ks,iq)+
     1                                   dAHjk*drgx(kx) !xyz drag
                     do jx=1,3
!--jx component of normal at dorsal (level=3) node js
                        snj=snn(jx,3,jstack)
                        snjk=snj*snk
!--vnormal= v.nn, vtangent= v - v.nn
                        Ldv(jx,js,kx,ks,iq)=Ldv(jx,js,kx,ks,iq)+
     1                         dAHjk*(drgt*(kd(jx,kx)-snjk)+drgn*snjk)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--edge boundary
      do il=1,nl
         drgx(1:3)=min(drge(1:3,iq),cpen) !xyz drag
         drgn=min(drge(4,il),cpen)
         drgt=min(drge(5,il),cpen)
         do isg=1,2 !loop over 2 stacks x 3 levels = 6 GPs
            do lvg=1,3
               dAen=wgp(lvg)*dAe(0,lvg,isg,il)
               do ks=1,2
                  kstack=isol(ks,il)
                  do kl=1,3
                     dAHk=dAen*S6(kl,ks,lvg,isg)
                     do kx=1,3
!--kx component of normal at edge node kl, ks
                        snk=snn(kx,kl,kstack)
                        do js=1,2
                           jstack=isol(js,il)
                           do jl=1,3
                              dAHjk=dAHk*S6(jl,js,lvg,isg)
            Lev(kx,jl,js,kx,kl,ks,il)=Lev(kx,jl,js,kx,kl,ks,il)+
     1                                            dAHjk*drgx(kx) !xyz drag
                              do jx=1,3
!--jx component of normal at edge node jl, js
                                 snj=snn(jx,jl,jstack)
                                 snjk=snj*snk !j,k component of nn
!--vnormal= v.nn, vtangent= v - v.nn
            Lev(jx,jl,js,kx,kl,ks,il)=Lev(jx,jl,js,kx,kl,ks,il)+
     1                           dAHjk*(drgt*(kd(jx,kx)-snjk)+drgn*snjk)
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine modrgvmx               
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  End subroutines for velocity stiffness matrices 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Begin subroutine for velocity load vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine moctrvld(vldq)
!
! This subroutine computes the contribution to the velocity
! load vector vldq from an isotropic network contractility.
!
!  USER SPECIFIED (constitutive law):
!     psi(ilv,istack)=contractile stress at each node (<0 = swelling).
!   
      implicit none
!
      real(8) vldq(3,3,4,NSM)!xyz; level; stack; element
!
      real(8) psig, facg
      integer iq, isg, lvg, is, istack, ilv, js, jl, jx
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 4 Gauss points stacks
            do lvg=1,3 !loop over 3 GPs levels
!--compute contractility at Gauss point
               psig=0d0
               do is=1,4
                  istack=isoq(is,iq)
                  do ilv=1,3
                     psig=psig+psi(ilv,istack)*H(ilv,is,lvg,isg)
                  enddo
               enddo
               facg=wgp(lvg)*det(lvg,isg,iq)*psig
!--compute contribution of node jl,js of iq to jx component of load
               do js=1,4
                  do jl=1,3
                     do jx=1,3
                        vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)-
     1                                    facg*dHg(jx,jl,js,lvg,isg,iq)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return 
!
      end subroutine moctrvld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      subroutine mobodvld(vldq)
!
! This subroutine computes the contribution to the velocity
! load vector vldq from an external body force.
!
!  USER SPECIFIED (constitutive law):
!   bfr(3,ilv,istack)=force (density) vector at each node
!
      implicit none
!
      real(8) vldq(3,3,4,NSM)!velocity load vector
!
      real(8) fg(3)
      integer iq, isg, lvg, isn, istack, lvn, ix
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over Gauss points
            do lvg=1,3
! compute the force at the GP
               fg=0d0
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     do ix=1,3
                        fg(ix)=fg(ix)+
     1                         bfr(ix,lvn,istack)*H(lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
! weigh by volume associate with GP
               do ix=1,3
                  fg(ix)=fg(ix)*det(lvg,isg,iq)*wgp(lvg)
               enddo
! apply load to the nodes
               do isn=1,4
                  do lvn=1,3
                     do ix=1,3
                        vldq(ix,lvn,isn,iq)=vldq(ix,lvn,isn,iq)+
     1                                      fg(ix)*H(lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine mobodvld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mosfrvld(vldq)
!
! This subroutine computes the contribution to the velocity
! load vector vldq from a surface disjoining stress.
!
!  USER SPECIFIED (constitutive law):
!     disjoining stresses (if <0 it is a contractile stress).
!     real(8) sfrv(NSM) ventral
!     real(8) sfrd(NSM) dorsal
!     real(8) sfre(NLM) edge
!
      implicit none
!
      real(8) vldq(3,3,4,NSM)!velocity load vector
!
      real(8) dsfr(3,3,4)
      real(8) dvol, volg, ddAvn, ddAdn, ddAen, facg 
      real(8) unx, uny, unz, sfrxx, sfrxy, sfrxz, sfryy, sfryz, sfrzz
      integer iq,il,lvg,isg,js,jl,jx,isp,ism,is,istackg
!
! NOTE: the code could be made more compact, but would then be less
!       transparent and much more difficult to understand.
!
!--loop over ventral surfaces
      do iq=1,nq
         volg=0d0
         dsfr=0d0
         do isg=1,4
            ddAvn=1d0/dAv(0,isg,iq)
!--outward unit normal surface GP
            unx=dAv(1,isg,iq)*ddAvn
            uny=dAv(2,isg,iq)*ddAvn
            unz=dAv(3,isg,iq)*ddAvn
            dvol=wgp(1)*det(1,isg,iq)!volume of associated interior GP
            volg=volg+dvol!keep a running sum for normalization later
            facg=dvol*dAv(0,isg,iq)*sfrv(iq)
!--construct the stress tensor
            sfrxx=facg*unx*unx
            sfrxy=facg*unx*uny
            sfrxz=facg*unx*unz
            sfryy=facg*uny*uny
            sfryz=facg*uny*unz
            sfrzz=facg*unz*unz
!--compute (nabla H).(sfr)
            do js=1,4
               do jl=1,3
                  dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxx+
     2                          dHg(2,jl,js,1,isg,iq)*sfrxy+
     3                          dHg(3,jl,js,1,isg,iq)*sfrxz
                  dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxy+
     2                          dHg(2,jl,js,1,isg,iq)*sfryy+
     3                          dHg(3,jl,js,1,isg,iq)*sfryz
                  dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxz+
     2                          dHg(2,jl,js,1,isg,iq)*sfryz+
     3                          dHg(3,jl,js,1,isg,iq)*sfrzz
               enddo
            enddo
         enddo
!--the terms are weighted by volume, normalize
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                                 dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
!--loop over dorsal surfaces (same as ventral)
      do iq=1,nq
         volg=0d0
         dsfr=0d0
         do isg=1,4
            ddAdn=1d0/dAd(0,isg,iq)
            unx=dAd(1,isg,iq)*ddAdn
            uny=dAd(2,isg,iq)*ddAdn
            unz=dAd(3,isg,iq)*ddAdn
            dvol=wgp(3)*det(3,isg,iq)
            volg=volg+dvol
            facg=dvol*dAd(0,isg,iq)*sfrd(iq)
            sfrxx=facg*unx*unx
            sfrxy=facg*unx*uny
            sfrxz=facg*unx*unz
            sfryy=facg*uny*uny
            sfryz=facg*uny*unz
            sfrzz=facg*unz*unz
            do js=1,4
               do jl=1,3
                  dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxx+
     2                          dHg(2,jl,js,3,isg,iq)*sfrxy+
     3                          dHg(3,jl,js,3,isg,iq)*sfrxz
                  dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxy+
     2                          dHg(2,jl,js,3,isg,iq)*sfryy+
     3                          dHg(3,jl,js,3,isg,iq)*sfryz
                  dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxz+
     2                          dHg(2,jl,js,3,isg,iq)*sfryz+
     3                          dHg(3,jl,js,3,isg,iq)*sfrzz
               enddo
            enddo
         enddo
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                              dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
!--loop over edges
      do il=1,nl
         ism=isol(1,il) !determine the two stacks
         isp=isol(2,il)
         iq=iqol(il)
         volg=0d0
         dsfr=0d0
         do isg=1,4 !loop over interior GPs
            istackg=isoq(isg,iq)
!--check whether the GP belongs to a stack close to the edge
            if (istackg.eq.ism.or.istackg.eq.isp) then
               if(istackg.eq.ism) is=1
               if(istackg.eq.isp) is=2
               do lvg=1,3 !loop over levels in the stack
!--unit normal at corresponding surface GP
                  ddAen=1d0/dAe(0,lvg,is,il)
                  unx=dAe(1,lvg,is,il)*ddAen
                  uny=dAe(2,lvg,is,il)*ddAen
                  unz=dAe(3,lvg,is,il)*ddAen
!--the volume of the associated interior GP
                  dvol=wgp(lvg)*det(lvg,isg,iq)
                  volg=volg+dvol
                  facg=dvol*wgp(lvg)*dAe(0,lvg,is,il)*sfre(il)
!--compute the stress tensor
                  sfrxx=facg*unx*unx
                  sfrxy=facg*unx*uny
                  sfrxz=facg*unx*unz
                  sfryy=facg*uny*uny
                  sfryz=facg*uny*unz
                  sfrzz=facg*unz*unz
!                 print *,'sfrxx',sfrxx
!--compute (nabla H).(sfr)
                  do js=1,4
                     do jl=1,3
                        dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxx+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfrxy+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfrxz
                        dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxy+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfryy+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfryz
                        dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxz+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfryz+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfrzz
                     enddo
                  enddo
               enddo
            endif
         enddo
!--the terms are weighted by volume, normalize
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                              dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
      return 
      end subroutine mosfrvld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mosfrvld2(vldq)
!
! This subroutine computes the contribution to the velocity
! load vector vldq from a surface node-defined disjoining stress.
!
!  USER SPECIFIED (constitutive law):
!     disjoining stresses (if <0 it is a contractile stress).
!     real(8) sfrn(3,NSM) 
!
      implicit none
!
      real(8) vldq(3,3,4,NSM)!velocity load vector
!
      real(8) dsfr(3,3,4)
      real(8) dvol, volg, ddAvn, ddAdn, ddAen, facg, sfrg
      real(8) unx, uny, unz, sfrxx, sfrxy, sfrxz, sfryy, sfryz, sfrzz
      integer iq,il,lvg,isg,js,jl,jx,isp,ism,is,istackg,jstack
!
! NOTE: the code could be made more compact, but would then be less
!       transparent and much more difficult to understand.
!
!--loop over ventral surfaces
      do iq=1,nq
         volg=0d0
         dsfr=0d0
         do isg=1,4
            ddAvn=1d0/dAv(0,isg,iq)
!--outward unit normal surface GP
            unx=dAv(1,isg,iq)*ddAvn
            uny=dAv(2,isg,iq)*ddAvn
            unz=dAv(3,isg,iq)*ddAvn
            dvol=wgp(1)*det(1,isg,iq)!volume of associated interior GP
            volg=volg+dvol!keep a running sum for normalization later
            sfrg=0d0
            do js=1,4
               jstack=isoq(js,iq)
               sfrg=sfrg+S4(js,isg)*sfrn(1,jstack)
            enddo
            facg=dvol*dAv(0,isg,iq)*sfrg
!--construct the stress tensor
            sfrxx=facg*unx*unx
            sfrxy=facg*unx*uny
            sfrxz=facg*unx*unz
            sfryy=facg*uny*uny
            sfryz=facg*uny*unz
            sfrzz=facg*unz*unz
!--compute (nabla H).(sfr)
            do js=1,4
               do jl=1,3
                  dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxx+
     2                          dHg(2,jl,js,1,isg,iq)*sfrxy+
     3                          dHg(3,jl,js,1,isg,iq)*sfrxz
                  dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxy+
     2                          dHg(2,jl,js,1,isg,iq)*sfryy+
     3                          dHg(3,jl,js,1,isg,iq)*sfryz
                  dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                          dHg(1,jl,js,1,isg,iq)*sfrxz+
     2                          dHg(2,jl,js,1,isg,iq)*sfryz+
     3                          dHg(3,jl,js,1,isg,iq)*sfrzz
               enddo
            enddo
         enddo
!--the terms are weighted by volume, normalize
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                                 dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
!--loop over dorsal surfaces (same as ventral)
      do iq=1,nq
         volg=0d0
         dsfr=0d0
         do isg=1,4
            ddAdn=1d0/dAd(0,isg,iq)
            unx=dAd(1,isg,iq)*ddAdn
            uny=dAd(2,isg,iq)*ddAdn
            unz=dAd(3,isg,iq)*ddAdn
            dvol=wgp(3)*det(3,isg,iq)
            volg=volg+dvol
            sfrg=0d0
            do js=1,4
               jstack=isoq(js,iq)
               sfrg=sfrg+S4(js,isg)*sfrn(3,jstack)
            enddo
            facg=dvol*dAd(0,isg,iq)*sfrg
            sfrxx=facg*unx*unx
            sfrxy=facg*unx*uny
            sfrxz=facg*unx*unz
            sfryy=facg*uny*uny
            sfryz=facg*uny*unz
            sfrzz=facg*unz*unz
            do js=1,4
               do jl=1,3
                  dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxx+
     2                          dHg(2,jl,js,3,isg,iq)*sfrxy+
     3                          dHg(3,jl,js,3,isg,iq)*sfrxz
                  dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxy+
     2                          dHg(2,jl,js,3,isg,iq)*sfryy+
     3                          dHg(3,jl,js,3,isg,iq)*sfryz
                  dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                          dHg(1,jl,js,3,isg,iq)*sfrxz+
     2                          dHg(2,jl,js,3,isg,iq)*sfryz+
     3                          dHg(3,jl,js,3,isg,iq)*sfrzz
               enddo
            enddo
         enddo
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                              dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
!--loop over edges
      do il=1,nl
         ism=isol(1,il) !determine the two stacks
         isp=isol(2,il)
         iq=iqol(il)
         volg=0d0
         dsfr=0d0
         do isg=1,4 !loop over interior GPs
            istackg=isoq(isg,iq)
!--check whether the GP belongs to a stack close to the edge
            if (istackg.eq.ism.or.istackg.eq.isp) then
               if(istackg.eq.ism) is=1
               if(istackg.eq.isp) is=2
               do lvg=1,3 !loop over levels in the stack
!--unit normal at corresponding surface GP
                  ddAen=1d0/dAe(0,lvg,is,il)
                  unx=dAe(1,lvg,is,il)*ddAen
                  uny=dAe(2,lvg,is,il)*ddAen
                  unz=dAe(3,lvg,is,il)*ddAen
!--the volume of the associated interior GP
                  dvol=wgp(lvg)*det(lvg,isg,iq)
                  volg=volg+dvol
                  sfrg=0d0
                  do js=1,2
                     jstack=isol(js,il)
                     do jl=1,3
                        sfrg=sfrg+S6(jl,js,lvg,isg)*sfrn(jl,jstack)
                     enddo
                  enddo
                  facg=dvol*wgp(lvg)*dAe(0,lvg,is,il)*sfrg
!--compute the stress tensor
                  sfrxx=facg*unx*unx
                  sfrxy=facg*unx*uny
                  sfrxz=facg*unx*unz
                  sfryy=facg*uny*uny
                  sfryz=facg*uny*unz
                  sfrzz=facg*unz*unz
!                 print *,'sfrxx',sfrxx
!--compute (nabla H).(sfr)
                  do js=1,4
                     do jl=1,3
                        dsfr(1,jl,js)=dsfr(1,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxx+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfrxy+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfrxz
                        dsfr(2,jl,js)=dsfr(2,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxy+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfryy+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfryz
                        dsfr(3,jl,js)=dsfr(3,jl,js)+
     1                                dHg(1,jl,js,lvg,isg,iq)*sfrxz+
     2                                dHg(2,jl,js,lvg,isg,iq)*sfryz+
     3                                dHg(3,jl,js,lvg,isg,iq)*sfrzz
                     enddo
                  enddo
               enddo
            endif
         enddo
!--the terms are weighted by volume, normalize
         do js=1,4
            do jl=1,3
               do jx=1,3
                  vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                              dsfr(jx,jl,js)/volg
               enddo
            enddo
         enddo
      enddo
!
      return 
      end subroutine mosfrvld2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine moxvnvld(vldv,vldd,vlde,Lvv,Ldv,Lev)
! 
! This routine works in conjuction with subroutine applies mofixvmx
! to fix velocities at the boundaries.
! The idea is that in the  Q * vnw = vload equation
! a matching term is introduced into vload to enforce a certain vnw.
!
! external velocity components x,y,z at boundary (set by user)
!   real(8) vbv(3,NSM)! ventral
!   real(8) vbd(3,NSM)! dorsal
!   real(8) vbe(3,NLM)! edge
!
      implicit none
!
!--boundary velocity load vectors
      real(8) vldv(3,4,NSM),vldd(3,4,NSM),vlde(3,3,2,NLM)
!--boundary velocity stiffness matrices
      real(8) Lvv(3,4,3,4,NSM), Ldv(3,4,3,4,NSM), Lev(3,3,2,3,3,2,NLM)
!
      real(8) vbei,vbvi,vbdi
      integer il, iq, isn, ilv, ix, jsn, jlv, jx
!
!--loop over ventral and dorsal surfaces
      do iq=1,nq
         do isn=1,4
            do ix=1,3
               do jsn=1,4
                  do jx=1,3
                     vldv(ix,isn,iq)=vldv(ix,isn,iq)+
     1                               vbv(jx,iq)*Lvv(jx,jsn,ix,isn,iq)
                     vldd(ix,isn,iq)=vldd(ix,isn,iq)+
     1                               vbd(jx,iq)*Ldv(jx,jsn,ix,isn,iq)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!--loop over edge surfaces
      do il=1,nl
         do isn=1,2
            do ilv=1,3
               do ix=1,3
                  do jsn=1,2
                     do jlv=1,3
                        do jx=1,3
                  vlde(ix,ilv,isn,il)=vlde(ix,ilv,isn,il)+
     1               vbe(jx,il)*Lev(jx,jlv,jsn,ix,ilv,isn,il)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine moxvnvld!(vldl,Lv)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine motrcvld(vldv,vldd,vlde)
!
! this subroutine computes the contribution to the boundary velocity
! load vectors from boundary tractions.
!
! external boundary tractions (specified by user)
!     real(8)trcv(3,NSM)!ventral
!     real(8)trcd(3,NSM)!dorsal
!     real(8)trce(3,NLM)!edge
!
      implicit none
!
      real(8) vldv(3,4,NSM),vldd(3,4,NSM),vlde(3,3,2,NLM)
!
      real(8) dAen
      integer iq,isg,il,lvg,js,jl,jx
!
!--ventral and dorsal tractions
      do iq=1,nq!loop overs surfaces
         do isg=1,4!loop over Gauss points
            do js=1,4!loop over nodes
               do jx=1,3
                  vldv(jx,js,iq)=vldv(jx,js,iq)+
     1                           trcv(jx,iq)*dAv(0,isg,iq)*S4(js,isg)
                  vldd(jx,js,iq)=vldd(jx,js,iq)+
     1                           trcd(jx,iq)*dAv(0,isg,iq)*S4(js,isg)
               enddo
            enddo
         enddo
      enddo
!
!--edge tractions
      do il=1,nl!loop over edge surfaces
         do isg=1,2
            do lvg=1,3
               dAen=dAe(0,lvg,isg,il)*wgp(lvg)
               do js=1,2
                  do jl=1,3
                     do jx=1,3
                        vlde(jx,lvg,isg,il)=vlde(jx,lvg,isg,il)
     1                      +trce(jx,il)*dAen*S6(jl,js,lvg,isg)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end subroutine motrcvld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mostnvld_mrc(vldv,vldd,vlde)
!
! This subroutine computes the surface tension contribution to
! the velocity loads per Marc Herant.
!
! The method uses the principle of virtual work:
!    given that the motion of a node by dx changes the area of a
!    a surface boundary element by dA, the work associated with
!    this motion is gamma*dA
! The force load associated with this work is therefore -gamma grad A
!    where grad A is the gradient of the area of the surface element
!    wrt to motion of the node.
!
! Let x,y,z be the standard extrinsic coordinates, and xi and eta
! the FEM surface intrinsic coordinates. One can show somewhat
! laboriously that:
!                      (dx/deta)         (dx/dxi)    (Nx)
!   grad A_GP = [dS(1) (dy/deta) - dS(2) (dy/dxi)] ^ (Ny)
!                      (dz/deta)         (dz/dxi)    (Nz)
!
!   where A_GP is the area associated with a Gauss point,
!   dS is the gradient of the shape function from the node wrt 
!   to xi(1) and eta(2) at that Gauss point, dx/deta... and dx/dxi
!   are the tangent vectors at this Gauss point, and
!   Nx... the normal vector at this Gauss point.
!
!   Then grad A = sum_GP w_GP grap A_GP  (sum over surface GPs)
!    where w_GP are the Gauss point weights
!
!  USER-SPECIFIED (constitutive law)
!     real(8) gamv(NSM)  ventral tension
!     real(8) gamd(NSM)  dorsal tension
!     real(8) game(NLM)  edge tension
!
      implicit none
! 
!--boundary (ventral, dorsal, edge) velocity loads
      real(8) vldv(3,4,NSM),vldd(3,4,NSM),vlde(3,3,2,NLM)
!
      real(8) dAvx,dAvy,dAvz,ddAvn
      real(8) dAdx,dAdy,dAdz,ddAdn
      real(8) dAex,dAey,dAez,ddAen
      real(8) felem(3,NSM)
      real(8) fac1x,fac1y,fac1z,fac2x,fac2y,fac2z,fx,fy,fz
      integer iq,isg,isn,ix,il,lvg,lvn,istack,is,ilv
!
!--ventral surface tension
      ftn=0d0 !the loads for surface tension are kept for debugging
      do iq=1,nq
         do isg=1,4
            dAvx=dAv(1,isg,iq)
            dAvy=dAv(2,isg,iq)
            dAvz=dAv(3,isg,iq)
            ddAvn=gamv(iq)/dAv(0,isg,iq)
            !print *,gamv(iq),dAv(0,isg,iq)
            fac1x=ddAvn*(dAvz*dySv(2,isg,iq)-dAvy*dzSv(2,isg,iq))
            fac1y=ddAvn*(dAvx*dzSv(2,isg,iq)-dAvz*dxSv(2,isg,iq))
            fac1z=ddAvn*(dAvy*dxSv(2,isg,iq)-dAvx*dySv(2,isg,iq))
            fac2x=ddAvn*(dAvz*dySv(1,isg,iq)-dAvy*dzSv(1,isg,iq))
            fac2y=ddAvn*(dAvx*dzSv(1,isg,iq)-dAvz*dxSv(1,isg,iq))
            fac2z=ddAvn*(dAvy*dxSv(1,isg,iq)-dAvx*dySv(1,isg,iq))
!--contribution of this GP to 4 nodes
            do isn=1,4
               fx=dS4(1,isn,isg)*fac1x-dS4(2,isn,isg)*fac2x
               fy=dS4(1,isn,isg)*fac1y-dS4(2,isn,isg)*fac2y
               fz=dS4(1,isn,isg)*fac1z-dS4(2,isn,isg)*fac2z
               vldv(1,isn,iq)=vldv(1,isn,iq)+fx
               vldv(2,isn,iq)=vldv(2,isn,iq)+fy
               vldv(3,isn,iq)=vldv(3,isn,iq)+fz
               !print *, vldv(2,isn,iq)
               istack=isoq(isn,iq)
               ftn(1,1,istack)=ftn(1,1,istack)+fx
               ftn(2,1,istack)=ftn(2,1,istack)+fy
               ftn(3,1,istack)=ftn(3,1,istack)+fz
            enddo
         enddo
      enddo
!
!--dorsal surface tension
      do iq=1,nq
         do isg=1,4
            dAdx=dAd(1,isg,iq)
            dAdy=dAd(2,isg,iq)
            dAdz=dAd(3,isg,iq)
            ddAdn=gamd(iq)/dAd(0,isg,iq)
            fac1x=ddAdn*(dAdz*dySd(2,isg,iq)-dAdy*dzSd(2,isg,iq))
            fac1y=ddAdn*(dAdx*dzSd(2,isg,iq)-dAdz*dxSd(2,isg,iq))
            fac1z=ddAdn*(dAdy*dxSd(2,isg,iq)-dAdx*dySd(2,isg,iq))
            fac2x=ddAdn*(dAdz*dySd(1,isg,iq)-dAdy*dzSd(1,isg,iq))
            fac2y=ddAdn*(dAdx*dzSd(1,isg,iq)-dAdz*dxSd(1,isg,iq))
            fac2z=ddAdn*(dAdy*dxSd(1,isg,iq)-dAdx*dySd(1,isg,iq))
!--contribution of this GP to 4 nodes
            do isn=1,4
               fx=dS4(1,isn,isg)*fac1x-dS4(2,isn,isg)*fac2x
               fy=dS4(1,isn,isg)*fac1y-dS4(2,isn,isg)*fac2y
               fz=dS4(1,isn,isg)*fac1z-dS4(2,isn,isg)*fac2z
!--fx sign different from above because dorsal node
!--ordering different from ventral
               vldd(1,isn,iq)=vldd(1,isn,iq)-fx
               vldd(2,isn,iq)=vldd(2,isn,iq)-fy
               vldd(3,isn,iq)=vldd(3,isn,iq)-fz
               istack=isoq(isn,iq)
               ftn(1,3,istack)=ftn(1,3,istack)-fx
               ftn(2,3,istack)=ftn(2,3,istack)-fy
               ftn(3,3,istack)=ftn(3,3,istack)-fz
            enddo
         enddo
      enddo
!
!--edge surface tension
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               dAex=dAe(1,lvg,isg,il)
               dAey=dAe(2,lvg,isg,il)
               dAez=dAe(3,lvg,isg,il)
               ddAen=wgp(lvg)*game(il)/dAe(0,lvg,isg,il)
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
!--contribution of this GP to 4 nodes
               do isn=1,2
                  istack=isol(isn,il)
                  do lvn=1,3
                     fx=dS6(1,lvn,isn,lvg,isg)*fac1x-
     1                  dS6(2,lvn,isn,lvg,isg)*fac2x
                     fy=dS6(1,lvn,isn,lvg,isg)*fac1y-
     1                  dS6(2,lvn,isn,lvg,isg)*fac2y
                     fz=dS6(1,lvn,isn,lvg,isg)*fac1z-
     1                  dS6(2,lvn,isn,lvg,isg)*fac2z
                     vlde(1,lvn,isn,il)=vlde(1,lvn,isn,il)-fx
                     vlde(2,lvn,isn,il)=vlde(2,lvn,isn,il)-fy
                     vlde(3,lvn,isn,il)=vlde(3,lvn,isn,il)-fz
                     ftn(1,lvn,istack)=ftn(1,lvn,istack)-fx
                     ftn(2,lvn,istack)=ftn(2,lvn,istack)-fy
                     ftn(3,lvn,istack)=ftn(3,lvn,istack)-fz
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      end subroutine mostnvld_mrc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mostnvld_mxd(vldv,vldd,vlde)
!
! This subroutine computes the surface tension contribution to
! the velocity loads per Micah Dembo
!
!  USER-SPECIFIED (constitutive law)
!     real(8) gamv(NSM)  ventral tension
!     real(8) gamd(NSM)  dorsal tension
!     real(8) game(NLM)  edge tension
!
      implicit none
!
      real(8) vldv(3,4,NSM),vldd(3,4,NSM),vlde(3,3,2,NLM)
!
      real(8) facv,facd,face,dAvn,dAdn,dAen
      integer iq,isg,isn,ix,il,lvg,lvn,istack,is
!
      ftn=0d0
!--ventral and dorsal
      do iq=1,nq
         do isg=1,4
            facv=dAv(0,isg,iq)*gamv(iq)
            facd=dAd(0,isg,iq)*gamd(iq)
            do isn=1,4
               istack=isoq(isn,iq)
               do ix=1,3
                  vldv(ix,isn,iq)=vldv(ix,isn,iq)-
     1                            dS4vg(ix,isn,isg,iq)*facv
        ftn(ix,1,istack)=ftn(ix,1,istack)-dS4vg(ix,isn,isg,iq)*facv
                  vldd(ix,isn,iq)=vldd(ix,isn,iq)-
     1                            dS4dg(ix,isn,isg,iq)*facd
        ftn(ix,3,istack)=ftn(ix,3,istack)-dS4dg(ix,isn,isg,iq)*facd
               enddo
            enddo
         enddo
      enddo
!
!--edges
      do il=1,nl
         iq=iqol(il)
         do isg=1,2
            do lvg=1,3
               face=wgp(lvg)*dAe(0,lvg,isg,il)*game(il)
               do isn=1,2
                  istack=isol(isn,il)
                  do lvn=1,3
                     do ix=1,3
                        vlde(ix,lvn,isn,il)=vlde(ix,lvn,isn,il)-
     1                                face*dS6g(ix,lvn,isn,lvg,isg,iq)
        ftn(ix,lvn,istack)=ftn(ix,lvn,istack)-
     1                     face*dS6g(ix,lvn,isn,lvg,isg,iq)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine mostnvld_mxd
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine moprsvld(vldq,pon)
!
!  This subroutine conputes the pressure field contribution
!  to the velocity load which is grad P
!
      implicit none
!
      real(8) vldq(3,3,4,NSM)!velocity load vector
      real(8) pon(NSM*3)!pressure field
!
      real(8) pig,wig
      integer iq,isg,lvg,is,istack,ilv,inode,js,jl,jx
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 12 Gauss points
            do lvg=1,3
!--compute pressure at Gauss point
               pig=0d0
               do is=1,4
                  istack=isoq(is,iq)
                  do ilv=1,3
                     inode=3*(istack-1)+ilv
                     pig=pig+pon(inode)*H(ilv,is,lvg,isg)
                  enddo
               enddo
               wig=pig*det(lvg,isg,iq)*wgp(lvg)
!
!--loop p contribution over nodes
               do js=1,4 
                  do jl=1,3
                     do jx=1,3 !gradient component in jx direction
                        vldq(jx,jl,js,iq)=vldq(jx,jl,js,iq)+
     1                                    wig*dHg(jx,jl,js,lvg,isg,iq)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return 
      end subroutine moprsvld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  End subroutines for velocity load vectors
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!   END VELOCITY EQUATION SUBROUTINES
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
      subroutine modivv(divvsum,divvav,vnw)
!
! this subroutine computes the divergence of the network
! velocity field vnw and reports its volume weighted average divave
!
      implicit none
!
      real(8) divvsum !sum divergence
      real(8) divvav !average absolute divergence
      real(8) vnw(3,3*NSM) !network velocity
!
      real(8) divg,volsum,volg
      integer iq,isg,lvg,isn,nstack,lvn,inode,ix
!
      divvav=0d0
      divvsum=0d0
      volsum=0d0
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over Gauss points
            do lvg=1,3
               volg=det(lvg,isg,iq)*wgp(lvg)
!--for each node, add up contribution of the GP to divergence
               divg=0d0
               do isn=1,4 
                  nstack=isoq(isn,iq)
                  do lvn=1,3
                     inode=(nstack-1)*3+lvn
                     do ix=1,3
                        divg=divg+vnw(ix,inode)*
     1                            dHg(ix,lvn,isn,lvg,isg,iq)
                     enddo
                  enddo
               enddo
               divvav=divvav+abs(divg)*volg
               divvsum=divvsum+divg*volg
               volsum=volsum+volg
            enddo
         enddo
      enddo
      divvav=divvav/volsum
      divvsum=divvsum/volsum
!
      return
      end subroutine modivv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   BEGIN PRESSURE EQUATION SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mophichk
!
!   this routine computes the maximum inverse permeability (resistance) 
!   phimx(level,stack) of the  solvent wrt to the network that is
!   allowable for network stability.
!   phimx=viscosity/L^2  where L is the length scale of an element
!
!   if idphi=0 we set phi=phimx for Stokes flow
!   if idphi=1 we check that phi<phimx and issue a warning if warranted
!
!  USER-SPECIFIED (constitutive law unless Stokes limit):
!  phi(level,stack)=network inverse permeability
!
      implicit none
!
      logical,save::firstwarn=.true.
      real(8) phimx(3,NSM)
      real(8) dvol, chk
      integer iq,is,istack,ilv,nodemx(2)
!
      phimx=0d0
! for each node, find the corresponding Gauss point associated with
! the biggest volume 
      do iq=1,nq!loop over elements
         do is=1,4
            istack=isoq(is,iq)
            do ilv=1,3
               dvol=wgp(ilv)*det(ilv,is,iq)
               phimx(ilv,istack)=max(phimx(ilv,istack),dvol)
            enddo
         enddo
      enddo
! now loop again to get it right
      do istack=1,ns
       do ilv=1,3
          phimx(ilv,istack)=vis(ilv,istack)/phimx(ilv,istack)**(2./3.)
       enddo
      enddo
!
      if(idphi.eq.1) then
         chk=maxval(phi(1:3,1:ns)/phimx(1:3,1:ns))
         if (chk.gt.1d0.and.firstwarn) then
            print *,
     1     'molib: mophichk: WARNING! excessive resistance phi: chk=',
     2                                                          chk
            nodemx=maxloc(phi(1:3,1:ns)/phimx(1:3,1:ns))
            print *,'stack,level',nodemx
            print *,'viscosity',vis(nodemx(1),nodemx(2))
            print *,'phi',phi(nodemx(1),nodemx(2))
            print *,'phimx',phimx(nodemx(1),nodemx(2))
            print *,'press any key to continue'
            read *
            firstwarn=.false.
         endif
      else
         phi=phimx
      endif
!
      return
      end subroutine mophichk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Begin subroutines for pressure stiffness matrices 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mophipmx(Qp)
!
!  Assemble Qp, the FEM stiffness matrix operator for
!    div [grad(P)/phi] used in the pressure equation
! Qp(jl,js,kl,ks,iq)=<(dHjx))*(dHkx)+(dHjy)*(dHky)+(dHjz)*(dHkz)>|iq
!
!  USER-SPECIFIED (constitutive law unless Stokes limit):
!  phi(level,stack)=network inverse permeability
!
      implicit none
!
      real(8) Qp(3,4,3,4,NSM)!j level, stack; k level, stack, element
!
      real(8) phig,wig,dQ
      integer iq,isg,lvg,isn,nstack,lvn,ks,kl,js,jl
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 12 Gauss points
            do lvg=1,3
!--compute phi at the Gauss point
               phig=0d0
               do isn=1,4
                  nstack=isoq(isn,iq)
                  do lvn=1,3
                     phig=phig+phi(lvn,nstack)*H(lvn,isn,lvg,isg)
                  enddo
               enddo
               if (phig.le.0d0) then
                  print *,'mophipmx warning! phig<0'
                  print *,'iq,isg,lvg,phig',iq,isg,lvg,phig
!--reset to minimum phig in element
                  phig=vbig
                  do isn=1,4
                     nstack=isoq(isn,iq)
                     do lvn=1,3
                        phig=min(phig,phi(lvn,nstack))
                     enddo
                  enddo
                  if (phig.le.0d0) then
                     print *,'phig<0,iq',phig,iq
                     stop
                  endif
               endif
               wig=det(lvg,isg,iq)*wgp(lvg)/phig
               if (wig.lt.0d0) then
                  print *,'wig,det(lvg,isg,iq),phig',
     1                     wig,det(lvg,isg,iq),phig
                  print *,'iq,isg,lvg',iq,isg,lvg
                  stop
               endif
!
!--compute the stiffness matrix by 
!--looping over all pairs of nodes in the element
               do ks=1,4
                  do kl=1,3
                     do js=1,4
                        do jl=1,3
             dQ=dHg(1,jl,js,lvg,isg,iq)*dHg(1,kl,ks,lvg,isg,iq)+
     1          dHg(2,jl,js,lvg,isg,iq)*dHg(2,kl,ks,lvg,isg,iq)+
     2          dHg(3,jl,js,lvg,isg,iq)*dHg(3,kl,ks,lvg,isg,iq)
             Qp(jl,js,kl,ks,iq)=Qp(jl,js,kl,ks,iq)+dQ*wig
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine mophipmx!(Qp)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine momaspmx(Qp,rlax)
!
! Add the rlax contribution from the mass matrix to the pressure 
! stiffness matrix : rlax*<Hj*Hk>|iq
!
      implicit none
!
      real(8) rlax !pressure relaxation const
      real(8) Qp(3,4,3,4,NSM) !pressure stiffness matrix
!
      real(8) dvol,Hk
      integer iq,isg,lvg,ks,kl,js,jl
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 12 Gauss points
            do lvg=1,3
               dvol=rlax*det(lvg,isg,iq)*wgp(lvg)
               do ks=1,4
                  do kl=1,3
                     Hk=H(kl,ks,lvg,isg)*dvol
                     do js=1,4
                        do jl=1,3
           Qp(jl,js,kl,ks,iq)=Qp(jl,js,kl,ks,iq)+H(jl,js,lvg,isg)*Hk
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine momaspmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mohycpmx(Lvp,Ldp,Lep)
!
! This subroutine adds a hydraulic conductivity term
! to the boundary pressure stiffness matrices
! it works in conjunction with mohycpld which takes of the loads.
!
! USER SPECIFIED
!--solvent boundary permeability
!     real(8) hycv(NSM)!ventral
!     real(8) hycd(NSM)!dorsal
!     real(8) hyce(NLM)!edge
!
      implicit none
!
      real(8) Lvp(4,4,NSM), Ldp(4,4,NSM)!ventral/dorsal stiffness
      real(8) Lep(3,2,3,2,NLM)!edge boundary pressure stiffness
!
      real(8) facg,facgv,facgd,fack,fackv,fackd
      integer iq,il,isg,lvg,ks,kl,js,jl
!
      do iq=1,nq
         do isg=1,4
            facgv=dAv(0,isg,iq)*hycv(iq)
            facgd=dAd(0,isg,iq)*hycd(iq)
            do ks=1,4
               fackv=facgv*S4(ks,isg)
               fackd=facgd*S4(ks,isg)
               do js=1,4
                  Lvp(js,ks,iq)=Lvp(js,ks,iq)+fackv*S4(js,isg)
                  Ldp(js,ks,iq)=Ldp(js,ks,iq)+fackd*S4(js,isg)
               enddo
            enddo
         enddo
      enddo
!
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               facg=wgp(lvg)*dAe(0,lvg,isg,il)*hyce(il)
               do ks=1,2
                  do kl=1,3
                     fack=facg*S6(kl,ks,lvg,isg)
                     do js=1,2
                        do jl=1,3
                           Lep(jl,js,kl,ks,il)=Lep(jl,js,kl,ks,il)+
     1                                        fack*S6(jl,js,lvg,isg)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
!
      enddo
!
      return
      end subroutine mohycpmx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  End subroutines for pressure stiffness matrices 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Begin subroutine for pressure load vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine momaspld(pldq,pon,rlax)
!
! add the pressure mass matrix contribution to the pressure load
!
      real(8) pldq(3,4,NSM) !element pressure load
      real(8) pon(3*NSM) !pressure field
      real(8) rlax !pressure diffusion relaxation coefficient
!
      real(8) pig, wig
      integer iq,isg,lvg,is,istack,ilv,inode,js,jl
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 12 Gauss points
            do lvg=1,3
!
!--compute pressure at Gauss point
               pig=0d0
               do is=1,4
                  istack=isoq(is,iq)
                  do ilv=1,3
                     inode=3*(istack-1)+ilv
                     pig=pig+pon(inode)*H(ilv,is,lvg,isg)
                  enddo
               enddo
               wig=rlax*pig*det(lvg,isg,iq)*wgp(lvg)
!--add the pressure term to the pressure load
               do js=1,4
                  do jl=1,3
                     pldq(jl,js,iq)=pldq(jl,js,iq)+wig*H(jl,js,lvg,isg)
                  enddo
               enddo
            enddo 
         enddo
      enddo
!
      return
      end subroutine momaspld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine modivpld(pldq,vnw)
!
! this subroutine computes the load due to the network velocity
! divergence appearing in the pressure equation
!
      implicit none
!
      real(8) pldq(3,4,NSM) !pressure load
      real(8) vnw(3,3*NSM) !network velocity
!
      real(8) divg
      integer iq,isg,lvg,isn,istack,lvn,inode
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over 12 Gauss points
            do lvg=1,3
!--compute velocity divergence at GP
               divg=0d0
               do isn=1,4
                  istack=isoq(isn,iq)
                  do lvn=1,3
                     inode=(istack-1)*3+lvn
                     divg=divg+dHg(1,lvn,isn,lvg,isg,iq)*vnw(1,inode)
     1                        +dHg(2,lvn,isn,lvg,isg,iq)*vnw(2,inode)
     2                        +dHg(3,lvn,isn,lvg,isg,iq)*vnw(3,inode)
                  enddo
               enddo
!--weigh by GP volume
               divg=divg*det(lvg,isg,iq)*wgp(lvg)
               do isn=1,4
                  do lvn=1,3
                     pldq(lvn,isn,iq)=pldq(lvn,isn,iq)-
     1                                H(lvn,isn,lvg,isg)*divg
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine modivpld
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mohycpld(pldv,pldd,plde)
!
! this subroutine adds the pressure loads due to 
! membrane hydraulic conductivity couple with differential
! osmotic activity. It works in conjunction with mohycpmx
!
! USER SPECIFIED
!--solvent boundary permeability
!     real(8) hycv(NSM)!ventral
!     real(8) hycd(NSM)!dorsal
!     real(8) hyce(NLM)!edge
!--solvent boundary osmotic activity
!     real(8)osmv(NSM)!ventral
!     real(8)osmd(NSM)!dorsal
!     real(8)osme(NLM)!edge
!
      implicit none
!
      real(8) pldv(4,NSM), pldd(4,NSM)!ventral/dorsal loads
      real(8) plde(3,2,NLM)! edge pressure load
!
      real(8) facg,facgv,facgd
      integer iq,il,isg,lvg,isn,lvn
!
!--ventral/dorsal boundary
      do iq=1,nq
         do isg=1,4
            facgv=dAv(0,isg,iq)*hycv(iq)*osmv(iq)
            facgd=dAd(0,isg,iq)*hycd(iq)*osmd(iq)
            do isn=1,4
               pldv(isn,iq)=pldv(isn,iq)+facgv*S4(isn,isg)
               pldd(isn,iq)=pldd(isn,iq)+facgd*S4(isn,isg)
            enddo
         enddo
      enddo
!
!--edge boundary
      do il=1,nl
         do isg=1,2
            do lvg=1,3
               facg=wgp(lvg)*dAe(0,lvg,isg,il)*hyce(il)*osme(il)
               do isn=1,2
                  do lvn=1,3
                     plde(lvn,isn,il)=plde(lvn,isn,il)+
     1                                facg*S6(lvn,isn,lvg,isg)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      return
      end subroutine mohycpld!(pldl)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  End subroutine for pressure load vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   END PRESSURE EQUATION SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mochk(pdel,pnew,vdel,vnw)
!
!  check convergence
!
      implicit none
!
      real(8) pdel !relative pressure change
      real(8) pnew(3*NSM) !current pressure
      real(8) vdel !relative velocity change
      real(8) vnw(3,3*NSM) !current velocity
!
      real(8) pmax,pmin !max and min pressure
      real(8),save::pold(3*NSM)=0.0 !prior pressure
      real(8) vmax !max velocity norm
      real(8) vnew(3*NSM) !current velocity norm
      real(8),save::vold(3*NSM)=0.0 !prior velocity norm
!       
      pmax=maxval(pnew(1:3*ns))
      pmin=minval(pnew(1:3*ns))
      pdel=maxval(abs(pnew(1:3*ns)-pold(1:3*ns)))
      pold=pnew
      pdel=pdel/(1d-6+abs(pmax)+abs(pmin))
!
      vnew=sqrt(vnw(1,:)**2+vnw(2,:)**2+vnw(3,:)**2)
      vmax=maxval(vnew)
      vdel=maxval(abs(vnew(1:3*ns)-vold(1:3*ns)))
      vold=vnew
      !print *,'oldvdel',vdel
      vdel=vdel/(1d-16+vmax)
      !print *,'pmax,pmin,plde',pmax,pmin,pdel
      !print *,'vnew,vmax,vdel',vmax,vdel
!
      return
      end subroutine mochk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine mopanic(ipanic)
!
      integer ipanic
!
      integer isn,ilv
!
      if (ipanic.ge.100) then
         print *,'dumping in mopanic.dump'
         open(66,file='mopanic.dump')
         call iowrfile(99,66)
         stop
      else
         do isn=1,ns
            do ilv=1,3
               hvec(7,ilv,isn)=0d0
               hvec(8,ilv,isn)=0d0
               hvec(9,ilv,isn)=0d0
               hvec(10,ilv,isn)=0d0
            enddo
         enddo
         ipanic=2*ipanic+1
      endif
!
      return
      end subroutine mopanic
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine morlvpol(rlv,pon)
!
! this subroutine computes the relative volume velocity:
!     v_v - v_n = - grad P/phi
!
!   where v_v = theta_s * v_s + theta_n * v_n
!
!  USER-SPECIFIED (constitutive law unless Stokes limit):
!  phi(level,stack) = network inverse permeability (resistance)
!
      real(8) rlv(3,3,NSM) !relative velocity
      real(8) pon(3*NSM) !current pressure
!
      real(8) gradp(3,3,4) !pressure gradient at GP
      real(8) xg(3) !coordinates of GP
      real(8) won(3,NSM)!a weighting array
      integer nhyc(3,NSM)! hydraulic conductivity node flag
      real(8) dg,ddg
      integer iq,isg,igstack,lvg,isn,lvn,inode,ix,istack,is,ilv
!
      won=0d0
      rlv=0d0
!
      do iq=1,nq !loop over elements
         do isg=1,4 !loop over Gauss points
            igstack=isoq(isg,iq)
            do lvg=1,3
!--compute the gradient of the pressure and the position at the GP
               gradp=0d0
               xg=0d0
               do isn=1,4
                  nstack=isoq(isn,iq)
                  do lvn=1,3
                     inode=(nstack-1)*3+lvn
                     do ix=1,3
                        gradp(ix,lvg,isg)=gradp(ix,lvg,isg)+
     1                            pon(inode)*dHg(ix,lvn,isn,lvg,isg,iq)
                        xg(ix)=xg(ix)+
     1                         hvec(ix,lvn,nstack)*H(lvn,isn,lvg,isg)
                     enddo
                  enddo
               enddo
!--dg=distance of GP point to its "associated" node
               dg=sqrt((xg(1)-hvec(1,lvg,igstack))**2+
     1                 (xg(2)-hvec(2,lvg,igstack))**2+
     1                 (xg(3)-hvec(3,lvg,igstack))**2)
               ddg=1d0/(max(vtiny,dg))
!--add the gradient to the GP associated node with weigh 1/distance
               do ix=1,3
                  rlv(ix,lvg,igstack)=rlv(ix,lvg,igstack)+
     1                                ddg*gradp(ix,lvg,isg)
               enddo 
!--keep track of the sum of weights for the node
               won(lvg,igstack)=won(lvg,igstack)+ddg
            enddo
         enddo
      enddo
!
!--set base approximation to the relative volume velocity
      do is=1,ns !loop over stacks
         do ilv=1,3 !loop over levels
            fac=-1d0/(phi(ilv,is)*won(ilv,is))
            do ix=1,3
               rlv(ix,ilv,is)=fac*rlv(ix,ilv,is)
            enddo
         enddo
      enddo
!
!--loop over all surface elements to determine those
!  that have permeability to the solvent
      nhyc=0
      if (idhyc.ne.0) then
!--ventral
         do iq=1,nq
            if (hycv(iq).gt.vtiny) then
               do is=1,4
                  nhyc(1,isoq(is,iq))=1
               enddo
            endif
         enddo
!--dorsal
         do iq=1,nq
            if (hycd(iq).gt.vtiny) then
               do is=1,4
                  nhyc(3,isoq(is,iq))=1
               enddo
            endif
         enddo
!--edge
         do il=1,nl
            if (hyce(il).gt.vtiny) then
               do is=1,2
                  istack=isol(is,il)
                  do ilv=1,3
                     nhyc(ilv,istack)=1
                  enddo
               enddo
            endif
         enddo
      endif
!
!--loop over all surface nodes; if nhyc = 0 then
!  we have zero permeability and the normal component
!  of the relative volume velocity needs to be zeroed out.
      do istack=1,ns !loop over stacks
!--ventral node
         if (nhyc(1,istack).eq.0) then
            rlvdotn=rlv(1,1,istack)*snn(1,1,istack)+
     1              rlv(2,1,istack)*snn(2,1,istack)+
     2              rlv(3,1,istack)*snn(3,1,istack)
            rlv(1,1,istack)=rlv(1,1,istack)-rlvdotn*snn(1,1,istack)
            rlv(2,1,istack)=rlv(2,1,istack)-rlvdotn*snn(2,1,istack)
            rlv(3,1,istack)=rlv(3,1,istack)-rlvdotn*snn(3,1,istack)
         endif
!--dorsal node
         if (nhyc(3,istack).eq.0) then
            rlvdotn=rlv(1,3,istack)*snn(1,3,istack)+
     1              rlv(2,3,istack)*snn(2,3,istack)+
     2              rlv(3,3,istack)*snn(3,3,istack)
            rlv(1,3,istack)=rlv(1,3,istack)-rlvdotn*snn(1,3,istack)
            rlv(2,3,istack)=rlv(2,3,istack)-rlvdotn*snn(2,3,istack)
            rlv(3,3,istack)=rlv(3,3,istack)-rlvdotn*snn(3,3,istack)
         endif
      enddo
!
!--edge middle nodes
      do il=1,nl
         istack=isol(1,il)
         if (nhyc(2,istack).eq.0) then
            rlvdotn=rlv(1,2,istack)*snn(1,2,istack)+
     1              rlv(2,2,istack)*snn(2,2,istack)+
     2              rlv(3,2,istack)*snn(3,2,istack)
            rlv(1,2,istack)=rlv(1,2,istack)-rlvdotn*snn(1,2,istack)
            rlv(2,2,istack)=rlv(2,2,istack)-rlvdotn*snn(2,2,istack)
            rlv(3,2,istack)=rlv(3,2,istack)-rlvdotn*snn(3,2,istack)
         endif
      enddo   
!          
      return
      end subroutine morlvpol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
      END MODULE MOLIBSW
