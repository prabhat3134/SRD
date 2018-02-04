      program sediment2d

C Colloids in solvent, with possibility to add polymers. The solvent
C is represented by 'particles' streaming continuously and coarse
C grained onto the cells of a regular lattice. All solvent particles
C interact with the colloid via a steep potential (vide infra).
C In each cell, after a  number of streaming (MD) steps, a collision
C is taking place, reorienting the velocities of the 'particles',
C but conserving energy, momentum, and mass (.i.e. hydrodynamics).
C
C 2D VERSION: flows / external forces are in y-direction.
C
C WALL VERSION: walls may be introduced in x- and/or y-direction.
C               (stick boundary walls in this version)
C 
C Interactions:
C
C Vcc = 4*V0*((sigcc/r)**48 - (sigcc/r)**24 + 1/4 )
C       + V0y*sigcc/r*exp(-kappa*(r-sigcc))
C Vpp = V0p*exp(-3r**2/4Rg**2)                      [may be changed later]
C Vss = 0 (except for collision steps)
C Vcp = 4*V0*((sigcp/r)**12 - (sigcp/r)**6 + 1/4 )  [may be changed later]
C Vcs = 4*V0*((sigcs/r)**12 - (sigcs/r)**6 + 1/4 )
C
C where V0 is about 2.5kT, V0y depends on the Yukawa screened coulomb
C interactions, V0p is of order kT, sigcc is the hard-core
C colloid diameter, kappa the inverse screening length,
C Rg is the polymer gyration radius, 
C and sigcp is the colloid-polymer interaction range (equal
C to sigcc/2+Rg). All hard interactions are cut off after 2**(1/n)
C sigma. The solvent may be gravity driven in the z-direction (effect
C of a pressure gradient), or the colloids may experience an external
C field, or both (to constrain the centre-of-mass).
C
C Author: Johan T. Padding
C Date:   24 November 2004

      integer N,Ndepl,Ncol,Ncolup
      integer Lx,Ly,Lxdepl,Lydepl,Lxcol,Lycol
      integer itherm,idens,iflow,ihydro,iorder
      integer icol,isort,ncell,ncelldepl,ncellcol,i,j
      integer totsteps,isave,icont,step,ssample,sframe,steller
      integer vsample,vsampletrue,vframe,vframetrue,vteller,vtellertrue
      integer fsample,fframe,fteller
      integer sqtsample,sqtframe,sqtteller,grsample,rbin,ybin,ipdb
      logical yxdens,fixcom,fixsolv,fixdepl,fixcol
      logical movie,wallx,wally
      double precision dt,Mdepl,Mcol,rlist,rcutcc,rcutccy,rcutcp,
     :                 rcutcs,rcutpp
      double precision sigcc,sigcp,sigcs,Rg,V0,V0y,V0p,rminsq,friccol,
     :                 fricdepl
      double precision kappa,omega,time,direct,pi,avdens,Foverlap
      double precision vflow,Etot,K,temp,Vpot,nu,ptot(2)
      double precision grav,vavdepl
      double precision sumvavdepl, dummy
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue
      common / block1b / rstore
      double precision,pointer,dimension(:,:) :: rstore
      common / block1c / fstore
      double precision,pointer,dimension(:,:) :: fstore
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      common / block2a / nearwall
      logical,pointer,dimension(:) :: nearwall
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) ::
     :                  Fdrx,Fdry,Fdrxsq,Fdrysq
      common / block5 / Rtemp,Vtemp,listtemp
      double precision,pointer,dimension(:) :: Rtemp,Vtemp,listtemp
#ifdef USEOPENMP
      integer omp_get_max_threads, maxthreads
      common / block5a / head_priv
      double precision,pointer,dimension(:,:) :: head_priv
#endif
      double precision store(1000)
      double precision sac(1000),scount(1000)
      double precision vact(2,1000),vacttrue(2,1000),vactcount(1000)
      double precision vactcounttrue(1000)
      double precision vacr(4,100),vacrcount(2,100),vsed(2),vsedtrue(2)
      double precision fact(2,10000),factcount(10000)
      double precision sqtx(10,100),sqty(10,100),
     :                 sqtcount(100)
      double precision dcount,dbin
      double precision grbin(100),delgr,totgr
      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision drsamp,rx,ry
      double precision maxcputime,time1,time2

#ifdef USEOPENMP
      maxthreads = omp_get_max_threads()
      write(*,*) "OpenMP max threads:", maxthreads
#endif
      call CPU_TIME (time1)
      call readcontrol(Ncol,Ncolup,Mcol,sigcc,kappa,V0,V0y,rcutccy,
Csed
     : sigcs,friccol,fricdepl,
Csed
     : Ndepl,Mdepl,Rg,V0p,N,icol,dt,totsteps,maxcputime,isave,icont,
     : rlist,Lx,Ly,wallx,wally,iflow,vflow,itherm,
     : ihydro,grav,omega,fixcom,fixsolv,fixdepl,fixcol,
     : fsample,vsample,vsampletrue,
     : ssample,sqtsample,idens,dbin,grsample,delgr,iorder,ipdb,movie)
      
      sigcp = 0.5d0*sigcc + Rg
      rcutcp = 2.d0**(1.d0/6.d0)*sigcp
      rcutcs = 2.d0**(1.d0/6.d0)*sigcs
      rcutcc = 2.d0**(1.d0/24.d0)*sigcc
      rcutpp = 3.5d0*Rg
      pi = 4.d0*datan(1.d0)
      avdens = dble(N)/
     :  (dble(Lx*Ly)-pi*(sigcs**2)*dble(Ncol))
      Foverlap = sigcs*avdens
      isort = INT(1.d0/dt) !sort when particles have moved one cell
      isort = isort - mod(isort,icol) !sort on integer number of icol
      if (isort.eq.0) isort = 1 !in case of a large dt

      if (V0p.ne.0.d0 .and. rlist.lt.rcutpp) stop 'rlist too small'

      Lxdepl = int(dble(Lx)/rlist)
      if (Lxdepl.lt.3) Lxdepl = 3
      Lydepl = int(dble(Ly)/rlist)
      if (Lydepl.lt.3) Lydepl = 3
      ncelldepl = Lxdepl*Lydepl

      if (V0y.eq.0.d0 .or. (V0y.ne.0.d0 .and. rcutcc.gt.rcutccy)) then
         Lxcol = int(dble(Lx)/rcutcc)
         if (Lxcol.lt.3) Lxcol = 3
         Lycol = int(dble(Ly)/rcutcc)
         if (Lycol.lt.3) Lycol = 3
      else
         Lxcol = int(dble(Lx)/rcutccy)
         if (Lxcol.lt.3) Lxcol = 3
         Lycol = int(dble(Ly)/rcutccy)
         if (Lycol.lt.3) Lycol = 3
      endif
      ncellcol = Lxcol*Lycol

      ncell = (Lx+1)*(Ly+1)
      yxdens = .false.
      if (idens.lt.0) then
         idens = -idens
         yxdens = .true.
      endif

#ifdef USEOPENMP
      allocate(head_priv(ncell,maxthreads))
#endif
      allocate(listtemp(N),Rtemp(3*N),Vtemp(3*N))

      if (Ncol.gt.0) then
         rminsq = 0.25d0*(dble(ncell)/dble(Ncol))
      else
         rminsq = 0.d0
      endif

      if (N.ne.0) then
         allocate(R(2*N),V(2*N),F(2*N))
         allocate(list(N))
      else
         allocate(R(2),V(2),F(2))
         allocate(list(1))
         do j=1,2
            F(j) = 0.d0
         enddo
         list(1) = 0
      endif
      if (Ndepl.ne.0) then
         allocate(Rdepl(2*Ndepl),Vdepl(2*Ndepl),Fdepl(2*Ndepl))
         allocate(listdepl(Ndepl),biglistdepl(Ndepl))
      else
         allocate(Rdepl(2),Vdepl(2),Fdepl(2))
         allocate(listdepl(1),biglistdepl(1))
         do j=1,2
            Fdepl(j) = 0.d0
         enddo
         listdepl(1) = 0
      endif
      if (Ncol.ne.0) then
         allocate(Rcol(2*Ncol),Vcol(2*Ncol),Fcol(2*Ncol))
         allocate(Fdrx(2*Ncol),Fdry(2*Ncol))
         allocate(Fdrxsq(2*Ncol),Fdrysq(2*Ncol))
         allocate(Rcolprev(2*Ncol))
         allocate(vstore(1000,2*Ncol),vstoretrue(1000,2*Ncol))
         allocate(rstore(2*Ncol,1000),fstore(10000,2*Ncol))
         allocate(listcol(Ncol))
      else
         allocate(Rcol(2),Vcol(2),Fcol(2))
         allocate(Rcolprev(2))
         allocate(Fdrx(2),Fdry(2))
         allocate(Fdrxsq(2),Fdrysq(2))
         allocate(vstore(1000,2),vstoretrue(1000,2))
         allocate(rstore(2,1000),fstore(10000,2))
         allocate(listcol(1))
         do j = 1,2
            Fcol(j) = 0.d0
         enddo
         do j = 1,2
            Fdrx(j) = 0.d0
            Fdry(j) = 0.d0
            Fdrxsq(j) = 0.d0
            Fdrysq(j) = 0.d0
         enddo
         do j=1,2
            do i = 1,1000
               vstore(i,j) = 0.d0
               vstoretrue(i,j) = 0.d0
            enddo
         enddo
         listcol(1) = 0
      endif
      allocate(head(ncell),nearwall(ncell))
      allocate(headdepl(ncell))
      allocate(bigheaddepl(ncelldepl))
      allocate(headcol(ncellcol))
      allocate(map(4*ncelldepl))
      allocate(mapcol(4*ncellcol))
      if (idens.ne.0) then
         rbin = INT(0.5d0*dble(Lx)/dbin-0.5d0)
         ybin = INT(dble(Ly)/dbin-0.5d0)
         allocate(dens(ybin+1,rbin+1),denssolv(ybin+1,rbin+1))
         allocate(flowfieldr(ybin+1,rbin+1),flowfieldy(ybin+1,rbin+1))
         allocate(flowhit(ybin+1,rbin+1))
      endif     
      !write(*,*) 'Current step: ',step,randy()
      call loadranstate
     
      call vel_pos_init(N,Lx,Ly)
      if (icont.eq.1) then
         call readconf(step,N,Ndepl,Ncol)
         call loadcorr(drsamp,dcount,sac,scount,sframe,steller,
     :        fframe,fteller,fact,factcount,
     :        vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,vactcount,
     :        vactcounttrue,vframe,vframetrue,vteller,vtellertrue,
     :        grbin,totgr,sumvavdepl,sqtx,sqty,sqtcount,
     :        sqtframe,sqtteller)
      else
         sumvavdepl = 0.d0
         if (ssample.ne.0) then
            call setupsac(sac,scount,sframe,steller)
         endif
         if (sqtsample.ne.0) then
            call setupsqt(sqtx,sqty,sqtcount,sqtframe,sqtteller)
         endif
         if (fsample.ne.0) then
            call setupfac(fact,factcount,fframe,fteller)
         endif
         if (vsample.ne.0) then
            call setupvac(vsed,vacr,vacrcount,vact,vactcount,
     :                    vframe,vteller,Ncol)
         endif
         if (vsampletrue.ne.0) then
            call setupvactrue(vsedtrue,vacttrue,vactcounttrue,
     :                    vframetrue,vtellertrue,Ncol)
         endif
         if (grsample.ne.0) then
            call setupgr(grbin,totgr)
         endif
         call setupdrift(Ncol,drsamp)
         if (idens.ne.0) then
            call setupdens(dcount,rbin,ybin)
         endif
      endif


      call maps(Lxdepl,Lydepl,Lxcol,Lycol)

      call determinenearwall(Lx,Ly,wallx,wally)

      call cells(N,Ndepl,Ncol,Lx,Ly,Lxdepl,Lydepl,
     :   Lxcol,Lycol,
     :   icol,step,rx,ry,V0p,wallx,wally,isort,1)

      time = dble(step)*dt
      if (dsin(omega*time).ge.0.d0) then
         direct = 1.d0
      else
         direct = -1.d0
      endif

      call force(N,Ndepl,Ncol,Ncolup,Lx,Ly,wallx,wally,
     :  Lxdepl,Lydepl,Lxcol,Lycol,
     :  sigcc,kappa,sigcp,sigcs,Rg,rlist,rcutcc,rcutccy,rcutcp,rcutcs,
     :  rcutpp,Mdepl,V0,V0y,V0p,Vpot,grav,direct,rx,ry,
     :  fixcom,fixsolv,fixdepl,fixcol,Foverlap)

      drsamp = drsamp+1.d0
      write(*,*) 'Simulation of solvent flow around colloids'
      write(*,*)
     :   '   Step      Temp      Etot      Ekin      Epot',
     :   '       P-x       P-y'

10    if (step.lt.totsteps) then
         step = step+1
         call movea(N,Ndepl,Ncol,dt,Mdepl,Mcol,step,ssample,
     :              ihydro,fricdepl,friccol)
         call periodic(N,Ndepl,Ncol,Lx,Ly,wallx,wally,Rg,sigcc,
     :                 iflow,vflow,dt,ihydro)
         call cells(N,Ndepl,Ncol,Lx,Ly,Lxdepl,Lydepl,
     :      Lxcol,Lycol,
     :      icol,step,rx,ry,V0p,wallx,wally,isort,1)
         time = dble(step)*dt
         if (dsin(omega*time).ge.0.d0) then
            direct = 1.d0
         else
            direct = -1.d0
         endif
         call force(N,Ndepl,Ncol,Ncolup,Lx,Ly,wallx,wally,
     :    Lxdepl,Lydepl,Lxcol,Lycol,
     :    sigcc,kappa,sigcp,sigcs,Rg,rlist,rcutcc,rcutccy,rcutcp,rcutcs,
     :    rcutpp,Mdepl,V0,V0y,V0p,Vpot,grav,direct,rx,ry,
     :    fixcom,fixsolv,fixdepl,fixcol,Foverlap)
         drsamp = drsamp+1.d0
         if (ihydro.ne.-1) then
            call moveb(N,Ndepl,Ncol,dt,Mdepl,Mcol,vavdepl)
         else
            vavdepl = 0.d0
         endif
         sumvavdepl = sumvavdepl+vavdepl
         if (mod(step,icol).eq.0 .and. ihydro.ne.-1) then
            call collide(N,Ndepl,Ncol,Mdepl,Mcol,Lx,Ly,
     :                   wallx,wally,avdens,
     :                   ihydro,iflow,vflow,rx,ry,dt,step,ssample,
     :                   fixcom,fixsolv,fixdepl,itherm,K,temp,icol)
         endif
         if (ihydro.lt.-1) call anderson(N)

         if (idens.ne.0 .and. mod(step,idens).eq.0) then
            call densitymap(dcount,rbin,ybin,dbin,
     :                      N,Ndepl,Lx,Ly,yxdens)
         endif
         if (ssample.ne.0 .and. mod(step,ssample).eq.0) then
            call stressac(store,sac,scount,sframe,steller,
     :                    N,Ndepl,Mdepl)
         endif
         if (sqtsample.ne.0 .and. mod(step,sqtsample).eq.0) then
            call calcsqt(sqtx,sqty,sqtcount,sqtframe,sqtteller,
     :                   Ncol,Lx,Ly)
         endif
         if (fsample.ne.0 .and. mod(step,fsample).eq.0) then
            call calcfact(fact,factcount,fframe,fteller,Ncol)
         endif
         if (vsample.ne.0 .and. mod(step,vsample).eq.0) then
            call velac(vsed,vacr,vacrcount,sigcs,vact,vactcount,
     :                 vframe,vteller,Ncol,Lx,Ly,wallx,wally,
     :                 ihydro,dt,vsample)
         endif
         if (vsampletrue.ne.0 .and. mod(step,vsampletrue).eq.0) then
            call velactrue(vsedtrue,vacttrue,vactcounttrue,
     :                 vframetrue,vtellertrue,Ncol,ihydro)
         endif
         if (grsample.ne.0 .and. mod(step,grsample).eq.0) then
            call radial(grbin,totgr,delgr,Ncol,Lx,Ly,wallx,wally)
         endif
         if (iorder.ne.0 .and. mod(step,iorder).eq.0) then
            call order(rminsq,Ncol,Ncolup,Lx,step,dt,direct,grav)
         endif
         if (idens.ne.0 .and. mod(step,10*idens).eq.0) then
            call updatedens(dcount,rbin,ybin,dbin,yxdens)
            call updatedrift(Ncol,drsamp)
         endif
         if (ssample.ne.0 .and. mod(step,10*ssample).eq.0) then
            if (ihydro.ne.-1) then
               call momentum(N,Ndepl,Ncol,Mdepl,Mcol,ptot)
               if (mod(step,icol).ne.0) then !separate kin. calculation
                  call kinetic(N,Lx,Ly,ihydro,vflow,K,temp)
               endif
            else
               K = 0.d0
               temp = 1.d0
               ptot(1) = 0.d0
               ptot(2) = 0.d0
            endif
            call updatesac(sac,scount,ssample,Lx,Ly,
     :                     N,Ndepl,dt,Mdepl,nu)
            Etot = K+Vpot
            write(*,'(I8,7(1X,E9.3))') 
     :       step,temp,Etot,K,Vpot,(ptot(j),j=1,2)
         else if (ssample.eq.0 .and. mod(step,1000).eq.0) then
            if (ihydro.ne.-1) then
               call momentum(N,Ndepl,Ncol,Mdepl,Mcol,ptot)
               if (mod(step,icol).ne.0) then !separate kin. calculation
                  call kinetic(N,Lx,Ly,ihydro,vflow,K,temp)
               endif
            else
               K = 0.d0
               temp = 1.d0
               ptot(1) = 0.d0
               ptot(2) = 0.d0
            endif
            Etot = K+Vpot
            write(*,'(I8,7(1X,E9.3))') 
     :       step,temp,Etot,K,Vpot,(ptot(j),j=1,2)
         endif
         if (sqtsample.ne.0 .and. mod(step,10*sqtsample).eq.0) then
            call updatesqt(sqtx,sqty,sqtcount,
     :                     sqtsample,dt,Lx,Ly)
         endif
         if (fsample.ne.0 .and. mod(step,50*fsample).eq.0) then
            call updatefact(fact,factcount,fsample,dt)
         endif
         if (vsample.ne.0 .and. mod(step,10*vsample).eq.0) then
            call updatevac(vsed,vacr,vacrcount,sigcs,vact,vactcount,
     :                     vteller,vsample,dt,Ncol)
         endif
         if (vsampletrue.ne.0 .and. mod(step,10*vsampletrue).eq.0) then
            call updatevactrue(vsedtrue,vacttrue,vactcounttrue,
     :                     vtellertrue,vsampletrue,dt,Ncol)
         endif
         if (grsample.ne.0 .and. mod(step,10*grsample).eq.0) then
            call updategr(grbin,totgr,delgr,Ncol,Lx,Ly)
         endif
         if (ipdb.ne.0 .and. mod(step,ipdb).eq.0) then
            call pdb(movie,Lx,Ly,Ncol,Ncolup,step)
         endif
         if (mod(step,isave).eq.0) then
            call saveconf(step,N,Ndepl,Ncol)
            call savecorr(drsamp,dcount,sac,scount,sframe,steller,
     :             fframe,fteller,fact,factcount,
     :             vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,
     :             vactcount,vactcounttrue,vframe,vframetrue,
     :             vteller,vtellertrue,grbin,totgr,sumvavdepl,
     :             sqtx,sqty,sqtcount,sqtframe,sqtteller)
            call saveranstate
         endif
	 ! for velocity averaging-our algorithm
         if (step > 500000) then
         call cells(N,Ndepl,Ncol,Lx,Ly,Lxdepl,Lydepl,
     :      Lxcol,Lycol,
     :      icol,step,rx,ry,V0p,wallx,wally,isort,0)
              call vel_avg(Lx,Ly)
         end if

         call CPU_TIME (time2)
         if (time2-time1.gt.maxCPUtime) then
            write(*,*) 'Maximum cput time consumed'
            goto 20
         endif
         goto 10
      endif

 20   if (idens.ne.0) then
         call updatedens(dcount,rbin,ybin,dbin,yxdens)
         deallocate(dens,denssolv,flowfieldr,flowfieldy,flowhit)
      endif
      call updatedrift(Ncol,drsamp)
      if (ssample.ne.0) then
         call updatesac(sac,scount,ssample,Lx,Ly,
     :                  N,Ndepl,dt,Mdepl,nu)
      endif
      if (sqtsample.ne.0) then
         call updatesqt(sqtx,sqty,sqtcount,
     :                  sqtsample,dt,Lx,Ly)
      endif
      if (fsample.ne.0) then
         call updatefact(fact,factcount,fsample,dt)
      endif
      if (vsample.ne.0) then
         call updatevac(vsed,vacr,vacrcount,sigcs,vact,vactcount,
     :                  vteller,vsample,dt,Ncol)
      endif
      if (vsampletrue.ne.0) then
         call updatevactrue(vsedtrue,vacttrue,vactcounttrue,
     :                  vtellertrue,vsampletrue,dt,Ncol)
      endif
      if (grsample.ne.0) then
         call updategr(grbin,totgr,delgr,Ncol,Lx,Ly)
      endif
      call saveranstate

C      deallocate(R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol)
C      deallocate(Rcolprev)
C      deallocate(list,head,nearwall,vstore,vstoretrue,rstore,fstore)
C      deallocate(listdepl,biglistdepl,listcol)
C      deallocate(headdepl,bigheaddepl,headcol)
C      deallocate(Fdrx,Fdry,Fdrxsq,Fdrysq)

      end


      subroutine readcontrol(Ncol,Ncolup,Mcol,sigcc,kappa,V0,V0y,
     : rcutccy,
Csed
     : sigcs,friccol,fricdepl,
Csed
     : Ndepl,Mdepl,Rg,V0p,N,icol,dt,totsteps,maxcputime,isave,icont,
     : rlist,Lx,Ly,wallx,wally,iflow,vflow,
     : itherm,ihydro,grav,omega,fixcom,fixsolv,fixdepl,fixcol,fsample,
     : vsample,vsampletrue,ssample,sqtsample,idens,dbin,grsample,delgr,
     : iorder,ipdb,movie)

      integer N,Ndepl,Ncol,Ncolup
      integer Lx,Ly
      integer icol,totsteps,isave,icont,iwallx,iwally
      integer vsample,vsampletrue,ssample,idens,itherm,iflow,ihydro
      integer fsample,sqtsample,iorder,ipdb
      double precision  dt,Mdepl,Mcol,vflow,dbin,grav,omega,frac
      double precision  sigcc,kappa,V0,V0y,rcutccy,Rg,V0p,rlist,delgr
Csed
      double precision  sigcs,friccol,fricdepl,maxcputime
Csed
      logical fixcom,fixsolv,fixdepl,fixcol,movie,wallx,wally
      integer ifixsolv,ifixdepl,ifixcol,grsample

      open(unit=10,file='control.dat',form='formatted')
      read(10,*)
      read(10,*) Ncol
      read(10,*) frac
      read(10,*) Mcol
      read(10,*) sigcc
      read(10,*) V0
      read(10,*) V0y
      read(10,*) kappa
      read(10,*) rcutccy
      read(10,*) friccol
      read(10,*)
      read(10,*) Ndepl
      read(10,*) Mdepl
      read(10,*) Rg
      read(10,*) V0p
      read(10,*) fricdepl
      read(10,*)
      read(10,*) N
Csed
      read(10,*) sigcs
Csed
      read(10,*) icol
      read(10,*)
      read(10,*) dt
      read(10,*) totsteps
      read(10,*) maxcputime
      read(10,*) isave
      read(10,*) icont
      read(10,*) rlist
      read(10,*) Lx
      read(10,*) Ly
      read(10,*) iwallx
      read(10,*) iwally
      read(10,*) iflow
      read(10,*) vflow
      read(10,*) itherm
      read(10,*) ihydro
      read(10,*) grav
      read(10,*) omega
      read(10,*) ifixsolv
      read(10,*) ifixdepl
      read(10,*) ifixcol
      read(10,*)
      read(10,*) fsample
      read(10,*) vsampletrue
      read(10,*) vsample
      read(10,*) ssample
      read(10,*) sqtsample
      read(10,*) idens
      read(10,*) dbin
      read(10,*) grsample
      read(10,*) delgr
      read(10,*) iorder
      read(10,*) ipdb
      close(unit=10)

      Ncolup = INT(frac*dble(Ncol))
Csed
      if (grav.lt.0.d0) then
         grav = Mcol*grav
      endif
Csed
      wallx   = .false.
      wally   = .false.
      fixcom  = .false.
      fixsolv = .false.
      fixdepl = .false.
      fixcol  = .false.
      if (iwallx  .eq.1) wallx   = .true.
      if (iwally  .eq.1) wally   = .true.
      if (ifixsolv.eq.1) fixsolv = .true.
      if (ifixsolv.eq.-1) fixcom = .true.
      if (ifixdepl.eq.1) fixdepl = .true.
      if (ifixcol .eq.1) fixcol  = .true.

      if (fixsolv .and. grav.gt.0.d0)
     :    stop 'Cannot fix solvent AND have field on solvent'
      if (dble(Lx).lt.2.d0*rlist) stop 'Lx < 2 rlist'
      if (dble(Ly).lt.2.d0*rlist) stop 'Ly < 2 rlist'

      if (ssample.ne.0 .and. mod(ssample,icol).ne.0)
     :   stop 'ssample must be integer multiple of icol'

      movie = .false.
      if (ipdb.lt.0) then
         ipdb = -ipdb
         movie = .true.
      endif

      return
      end


      subroutine readconf(step,N,Ndepl,Ncol)

      integer step
      integer N,NN,Ncol,NNcol,Ndepl,NNdepl
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      open(unit=11,file='pos.dat',action ='read')
      read(11) step
      read(11) NN,NNdepl,NNcol
      if (NN.ne.N) stop 'Error in # solvent particles in config.dat'
      if (NNdepl.ne.Ndepl) stop 'Error in # depletants in config.dat'
      if (NNcol.ne.Ncol) stop 'Error in # colloids in config.dat'
      read(11) R,V,Rdepl,Vdepl,Rcol,Vcol
      close(unit=11)

      return
      end


      subroutine saveconf(step,N,Ndepl,Ncol)

      integer step
      integer N,Ndepl,Ncol
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      open(unit=11,file='config.dat',form='unformatted')
      write(11) step
      write(11) N,Ndepl,Ncol
      write(11) R,V,Rdepl,Vdepl,Rcol,Vcol
      close(unit=11)

      return
      end


      subroutine movea(N,Ndepl,Ncol,dt,Mdepl,Mcol,step,ssample,
     :                 ihydro,fricdepl,friccol)

      integer N,Ndepl,Ncol,i,i2,step,ssample,ihydro
      double precision  dt,dt2,dtsq2,Mdepl,Mdeplinv,Mcol,Mcolinv
      double precision  xprev,yprev,dx,dy
      double precision  vx,vy,fx,fy
      double precision  gauss,drdepl,drcol,fricdepl,friccol,dummy
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block3 / sxyfr,sxyrot,sxykin
      double precision sxyfr,sxyrot,sxykin

      if (ihydro.ne.-1) then

      dt2 = dt / 2.d0
      dtsq2 = dt * dt2
      Mdeplinv = 1.d0/Mdepl
      Mcolinv  = 1.d0/Mcol

      if (ssample.ne.0 .and. mod(step,ssample).eq.0) then
         sxykin = 0.d0
C$OMP parallel
C$OMP& default(shared)
C$OMP& private(i,i3,xprev,yprev,vx,vy,fx,fy,dx,dy)
C$OMP& reduction(+:sxykin)
C$OMP do
         do i = 1,N
            i2 = 2*(i-1)
            xprev = R(i2+1)
            yprev = R(i2+2)
            vx    = V(i2+1)
            vy    = V(i2+2)
            fx    = F(i2+1)
            fy    = F(i2+2)
            dx = dt*vx + dtsq2*fx
            dy = dt*vy + dtsq2*fy
            R(i2+1) = xprev+dx
            R(i2+2) = yprev+dy
            V(i2+1) = vx + dt2*fx
            V(i2+2) = vy + dt2*fy
            sxykin = sxykin + vy/dt*
     :               (dble(int(xprev+dx))-dble(int(xprev)))
         enddo
C$OMP end do
C$OMP end parallel
         do i = 1,Ndepl
            i2 = 2*(i-1)
            xprev = Rdepl(i2+1)
            yprev = Rdepl(i2+2)
            vx    = Vdepl(i2+1)
            vy    = Vdepl(i2+2)
            fx    = Fdepl(i2+1)
            fy    = Fdepl(i2+2)
            dx = dt*vx + dtsq2*Mdeplinv*fx
            dy = dt*vy + dtsq2*Mdeplinv*fy
            Rdepl(i2+1) = xprev+dx
            Rdepl(i2+2) = yprev+dy
            sxykin = sxykin + vy/dt*
     :               (dble(int(xprev+dx))-dble(int(xprev)))
            Vdepl(i2+1) = vx + dt2*Mdeplinv*fx
            Vdepl(i2+2) = vy + dt2*Mdeplinv*fy
         enddo
      else
C$OMP parallel
C$OMP& default(shared)
C$OMP& private(i,vx,fx)
C$OMP do
         do i = 1,2*N
            vx = V(i)
            fx = F(i)
            V(i) = vx + dt2*fx
            R(i) = R(i) + dt*vx + dtsq2*fx
         enddo
C$OMP end do
C$OMP end parallel
         do i = 1,2*Ndepl
            vx = Vdepl(i)
            fx = Fdepl(i)
            Vdepl(i) = vx + dt2*Mdeplinv*fx
            Rdepl(i) = Rdepl(i) + dt*vx + dtsq2*Mdeplinv*fx
         enddo
      endif

      do i = 1,2*Ncol
         vx = Vcol(i)
         fx = Fcol(i)
         Vcol(i) = vx + dt2*Mcolinv*fx
         Rcol(i) = Rcol(i) + dt*vx + dtsq2*Mcolinv*fx
      enddo

      else !pure Brownian dynamics

      Mdeplinv = dt/fricdepl
      Mcolinv  = dt/friccol
      drdepl = dsqrt(2.d0*Mdeplinv)
      drcol  = dsqrt(2.d0*Mcolinv)
      dx = 0.d0
      dy = 0.d0
      do i = 1,Ndepl
         fx = drdepl*gauss(dummy)
         fy = drdepl*gauss(dummy)
         dx = dx+Mdepl*fx
         dy = dy+Mdepl*fy
         Rdepl(2*i-1) = Rdepl(2*i-1) + Mdeplinv*Fdepl(2*i-1) + fx
         Rdepl(2*i)   = Rdepl(2*i)   + Mdeplinv*Fdepl(2*i)   + fy
      enddo
      do i = 1,Ncol
         fx = drcol*gauss(dummy)
         fy = drcol*gauss(dummy)
         dx = dx+Mcol*fx
         dy = dy+Mcol*fy
         Rcol(2*i-1) = Rcol(2*i-1) + Mcolinv*Fcol(2*i-1) + fx
         Rcol(2*i)   = Rcol(2*i)   + Mcolinv*Fcol(2*i)   + fy
      enddo
      dx = dx/(Mdepl*dble(Ndepl)+Mcol*dble(Ncol))
      dy = dy/(Mdepl*dble(Ndepl)+Mcol*dble(Ncol))
      do i = 1,Ndepl
         Rdepl(2*i-1) = Rdepl(2*i-1) - dx
         Rdepl(2*i)   = Rdepl(2*i)   - dy
      enddo
      do i = 1,Ncol
         Rcol(2*i-1) = Rcol(2*i-1) - dx
         Rcol(2*i)   = Rcol(2*i)   - dy
      enddo

      endif

      return
      end


      subroutine moveb(N,Ndepl,Ncol,dt,Mdepl,Mcol,vavdepl)

      integer N,Ndepl,Ncol,i
      double precision  dt,dt2,Mdepl,Mdeplinv,Mcol,Mcolinv,Mtot
      double precision  vavdepl
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      dt2 = dt / 2.d0
      Mdeplinv = 1.d0/Mdepl
      Mcolinv = 1.d0/Mcol

C$OMP parallel
C$OMP& default(shared)
C$OMP& private(i)
C$OMP do
      do i = 1,2*N
         V(i) = V(i) + dt2*F(i)
      enddo
C$OMP end do
C$OMP end parallel
      do i = 1,2*Ndepl
         Vdepl(i) = Vdepl(i) + dt2*Mdeplinv*Fdepl(i)
      enddo
      vavdepl = 0.d0
      if (Ndepl.gt.0) then
         do i = 1,Ndepl
            vavdepl = vavdepl+Vdepl(2*i)
         enddo
         vavdepl = vavdepl/dble(Ndepl)
      endif
      do i = 1,2*Ncol
         Vcol(i) = Vcol(i) + dt2*Mcolinv*Fcol(i)
      enddo

      return
      end


      subroutine periodic(N,Ndepl,Ncol,Lx,Ly,wallx,wally,Rg,sigcc,
     :                    iflow,vflow,dt,ihydro)

      integer N,Ndepl,Ncol,ihydro,i,i2
      integer Lx,Ly,iflow,Lybound
      logical wallx,wally
      double precision  vflow,dt,Rg,sigcc,hsigcc,Rg2,gauss,dummy,frac
      double precision  rLx,rLy,rLx2,rLy2,rx,ry,vx,vy
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      rLx = dble(Lx)
      rLy = dble(Ly)
      hsigcc = 0.5d0*sigcc
      Rg2 = 2.d0*Rg
      rLx2 = 2.d0*rLx
      rLy2 = 2.d0*rLy
      
      if (ihydro.ne.-1) then

C$OMP parallel
C$OMP& default(shared)
C$OMP& private(i,i2,rx,ry,vx,vy,frac)
C$OMP do
      do i = 1,N
         i2 = 2*i-1
         rx = R(i2)
         ry = R(i2+1)
         vx = V(i2)
         vy = V(i2+1)
         if (wallx) then
            if (rx.lt.0.d0) then
               frac = rx/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               R(i2)   = -rx
               R(i2+1) = ry - frac*vy
               V(i2)   = -vx
               V(i2+1) = -vy
            else if (rx.ge.rLx) then
               frac = (rx-rLx)/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               R(i2)   = rLx2-rx
               R(i2+1) = ry - frac*vy
               V(i2)   = -vx
               V(i2+1) = -vy
            endif
         else !.not.wallx
            if (rx.lt.0.d0) then
               R(i2) = rx+rLx
            else if (rx.ge.rLx) then
               R(i2) = rx-rLx
            endif
         endif !wallx
         if (wally) then
            if (ry.lt.0.d0) then
               frac = ry/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               R(i2)   = rx - frac*vx
               R(i2+1) = -ry
               V(i2)   = -vx
               V(i2+1) = -vy
            else if (ry.ge.rLy) then
               frac = (ry-rLy)/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               R(i2)   = rx - frac*vx
               R(i2+1) = rLy2-ry
               V(i2)   = -vx
               V(i2+1) = -vy
            endif
         else !.not.wally
            if (ry.lt.0.d0) then
               R(i2+1) = ry+rLy
            else if (ry.ge.rLy) then
               R(i2+1) = ry-rLy
            endif
         endif !wally
C Correct out of bound positions due to multiple bounce-back
         if (wallx) then
            rx = R(i2)
            if (rx.lt.0.d0) then
               R(i2) = 0.d0
            else if (rx.ge.rLx) then
               R(i2) = rLx-1.d-6
            endif
         endif
      enddo !i=1,N
C$OMP end do
C$OMP end parallel

      do i = 1,Ndepl
         i2 = 2*i-1 
         rx = Rdepl(i2) 
         ry = Rdepl(i2+1) 
         vx = Vdepl(i2)
         vy = Vdepl(i2+1)
         if (wallx) then
            if (rx.lt.Rg) then
               frac = (rx-Rg)/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rdepl(i2)   = Rg2-rx
               Rdepl(i2+1) = ry - frac*vy
               Vdepl(i2)   = -vx
               Vdepl(i2+1) = -vy
            else if (rx.ge.rLx-Rg) then
               frac = (rx+Rg-rLx)/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rdepl(i2)   = rLx2-Rg2-rx
               Rdepl(i2+1) = ry - frac*vy
               Vdepl(i2)   = -vx
               Vdepl(i2+1) = -vy
            endif
         else !.not.wallx
            if (rx.lt.0.d0) then
               Rdepl(i2) = rx+rLx
            else if (rx.ge.rLx) then
               Rdepl(i2) = rx-rLx
            endif
         endif !wallx
         if (wally) then
            if (ry.lt.Rg) then
               frac = (ry-Rg)/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rdepl(i2)   = rx - frac*vx
               Rdepl(i2+1) = Rg2-ry
               Vdepl(i2)   = -vx
               Vdepl(i2+1) = -vy
            else if (ry.ge.rLy-Rg) then
               frac = (ry+Rg-rLy)/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rdepl(i2)   = rx - frac*vx
               Rdepl(i2+1) = rLy2-Rg2-ry
               Vdepl(i2)   = -vx
               Vdepl(i2+1) = -vy
            endif
         else !.not.wally
            if (ry.lt.0.d0) then
               Rdepl(i2+1) = ry+rLy
            else if (ry.ge.rLy) then
               Rdepl(i2+1) = ry-rLy
            endif
         endif !wally
C Correct out of bound positions due to multiple bounce-back
         if (wallx) then
            rx = Rdepl(i2)
            if (rx.lt.Rg) then
               Rdepl(i2) = Rg
            else if (rx.ge.rLx-Rg) then
               Rdepl(i2) = rLx-Rg-1.d-6
            endif
         endif
      enddo !i=1,Ndepl

      do i = 1,Ncol
         i2 = 2*i-1 
         rx = Rcol(i2) 
         ry = Rcol(i2+1)  
         vx = Vcol(i2)
         vy = Vcol(i2+1)
         if (wallx) then
            if (rx.lt.hsigcc) then
               frac = (rx-hsigcc)/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rcol(i2)   = sigcc-rx
               Rcol(i2+1) = ry - frac*vy
               Vcol(i2)   = -vx
               Vcol(i2+1) = -vy
            else if (rx.ge.rLx-hsigcc) then
               frac = (rx+hsigcc-rLx)/(vx*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rcol(i2)   = rLx2-sigcc-rx
               Rcol(i2+1) = ry - frac*vy
               Vcol(i2)   = -vx
               Vcol(i2+1) = -vy
            endif
         else !.not.wallx
            if (rx.lt.0.d0) then
               Rcol(i2) = rx+rLx
            else if (rx.ge.rLx) then
               Rcol(i2) = rx-rLx
            endif
         endif !wallx
         if (wally) then
            if (ry.lt.hsigcc) then
               frac = (ry-hsigcc)/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rcol(i2)   = rx - frac*vx
               Rcol(i2+1) = sigcc-ry
               Vcol(i2)   = -vx
               Vcol(i2+1) = -vy
            else if (ry.ge.rLy-hsigcc) then
               frac = (ry+hsigcc-rLy)/(vy*dt)
               if (frac.lt.0.d0) frac = 0.d0
               if (frac.gt.1.d0) frac = 1.d0
               frac = frac*2.d0*dt
               Rcol(i2)   = rx - frac*vx
               Rcol(i2+1) = rLy2-sigcc-ry
               Vcol(i2)   = -vx
               Vcol(i2+1) = -vy
            endif
         else !.not.wally
            if (ry.lt.0.d0) then
               Rcol(i2+1) = ry+rLy
            else if (ry.ge.rLy) then
               Rcol(i2+1) = ry-rLy
            endif
         endif !wally
C Correct out of bound positions due to multiple bounce-back
         if (wallx) then
            rx = Rcol(i2)
            if (rx.lt.hsigcc) then
               Rcol(i2) = hsigcc
            else if (rx.ge.rLx-hsigcc) then
               Rcol(i2) = rLx-hsigcc-1.d-6
            endif
         endif
      enddo !i=1,Ncol

C Inflow boundary conditions for the solvent (and solvent only!)
      if (iflow.eq.1) then
         Lybound = INT(Ly/8) - 1
         if (Lybound.lt.0) Lybound = 0
         do i = 1, N
            i2 = 2*(i-1)
            if (INT(R(i2+2)).le.Lybound) then
               V(i2+1) = gauss(dummy)
               V(i2+2) = vflow+gauss(dummy)
            endif
         enddo
      endif

      else !ihydro.eq.-1, slip boundaries

      do i = 1,Ncol
         i2 = 2*i-1 
         rx = Rcol(i2) 
         ry = Rcol(i2+1)  
         if (wallx) then
            if (rx.lt.hsigcc) then
               Rcol(i2)   = sigcc-rx
            else if (rx.ge.rLx-hsigcc) then
               Rcol(i2)   = rLx2-sigcc-rx
            endif
         else !.not.wallx
            if (rx.lt.0.d0) then
               Rcol(i2) = rx+rLx
            else if (rx.ge.rLx) then
               Rcol(i2) = rx-rLx
            endif
         endif !wallx
         if (wally) then
            if (ry.lt.hsigcc) then
               Rcol(i2+1) = sigcc-ry
            else if (ry.ge.rLy-hsigcc) then
               Rcol(i2+1) = rLy2-sigcc-ry
            endif
         else !.not.wally
            if (ry.lt.0.d0) then
               Rcol(i2+1) = ry+rLy
            else if (ry.ge.rLy) then
               Rcol(i2+1) = ry-rLy
            endif
         endif !wally
      enddo

      endif

      return
      end


      subroutine maps(Lxdepl,Lydepl,Lxcol,Lycol)

      integer Lxdepl, Lydepl, Lxcol, Lycol
      integer ix, iy, imap, icellf
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol

         do 40 iy = 1, Lydepl

            do 30 ix = 1, Lxdepl

               imap = (icellf(ix,iy,Lxdepl,Lydepl)-1)*4

               MAP( IMAP + 1  ) = ICELLf(IX + 1,IY    ,Lxdepl,Lydepl)
               MAP( IMAP + 2  ) = ICELLf(IX + 1,IY + 1,Lxdepl,Lydepl)
               MAP( IMAP + 3  ) = ICELLf(IX    ,IY + 1,Lxdepl,Lydepl)
               MAP( IMAP + 4  ) = ICELLf(IX - 1,IY + 1,Lxdepl,Lydepl)

30          continue

40       continue

         do 140 iy = 1, Lycol

            do 130 ix = 1, Lxcol

               imap = (icellf(ix,iy,Lxcol,Lycol)-1)*4

               MAPcol( IMAP + 1  ) = ICELLf(IX + 1,IY    ,Lxcol,Lycol)
               MAPcol( IMAP + 2  ) = ICELLf(IX + 1,IY + 1,Lxcol,Lycol)
               MAPcol( IMAP + 3  ) = ICELLf(IX    ,IY + 1,Lxcol,Lycol)
               MAPcol( IMAP + 4  ) = ICELLf(IX - 1,IY + 1,Lxcol,Lycol)

130         continue

140      continue

      return
      end


      integer function icellf(ix,iy,Lx,Ly)

      integer ix,iy,Lx,Ly

      icellf = 1 + mod(ix-1+Lx,Lx)
     :           + mod(iy-1+Ly,Ly) *Lx

      return
      end


      subroutine determinenearwall(Lx,Ly,wallx,wally)

      integer Lx,Ly,ncell,icell,ix,iy
      logical wallx,wally
      common / block2a / nearwall
      logical,pointer,dimension(:) :: nearwall

      ncell = (Lx+1)*(Ly+1)
      do icell = 1,ncell
         nearwall(icell) = .false.
      enddo
      if (wallx) then
            do iy = 0,Ly
               do ix = 0,Lx,Lx
                  icell = 1 + ix + iy*(Lx+1)
                  nearwall(icell) = .true.
               enddo
            enddo
      endif
      if (wally) then
            do iy = 0,Ly,Ly
               do ix = 0,Lx
                  icell = 1 + ix + iy*(Lx+1)
                  nearwall(icell) = .true.
               enddo
            enddo
      endif

      return
      end


      subroutine cells(N,Ndepl,Ncol,Lx,Ly,Lxdepl,Lydepl,Lxcol,Lycol,
     :                 icol,step,rx,ry,V0p,wallx,wally,isort,rgs)

      integer N,Ndepl,Ncol,i,i2,icell,ncell,ncelldepl,ncellcol,rgs
      integer icol,step,isort,j,j2,m
      integer Lx,Ly,Lxdepl,Lydepl,Lxcol,Lycol,ix,iy
      double precision  randy,rx,ry,V0p,rLxinv,rLyinv
      double precision  x,y
      logical wallx, wally, wallz, empty
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      common / block5 / Rtemp,Vtemp,listtemp
      double precision,pointer,dimension(:) :: Rtemp,Vtemp,listtemp
#ifdef USEOPENMP
      common / block5a / head_priv
      double precision,pointer,dimension(:,:) :: head_priv
      integer omp_get_max_threads, maxthreads
      integer omp_get_num_threads, omp_get_thread_num
      integer ithread, this_thread, numthreads
#endif

      ncell = (Lx+1)*(Ly+1)
      do i = 1, ncell
         head(i)  = 0
         headdepl(i) = 0
      enddo
      rx = 0.d0
      ry = 0.d0

      if (mod(step,icol).eq.0) then
C        random displacements [0,1] to retain Gaussian invariance
         rx = randy()
         ry = randy()
      endif
      if (rgs == 0) then
	      rx = 0.d0
	      ry = 0.d0
      end if

#ifdef USEOPENMP
      if (mod(step,isort).eq.0) then !sort SRD coordinates every isort
C$OMP parallel 
C$OMP& default(shared)
C$OMP& private(this_thread,icell,i,i2,j,x,y,ix,iy)
         numthreads = omp_get_num_threads()
         this_thread = omp_get_thread_num() + 1
         do icell = 1, ncell
            head_priv(icell,this_thread) = 0
         enddo
C$OMP do
         do i = 1,N
            i2 = 2*(i-1)
            x = R(i2+1)
            y = R(i2+2)
            Rtemp(i2+1) = x
            Rtemp(i2+2) = y
            do m = 1,2
               Vtemp(i2+m) = V(i2+m)
            enddo
            ix = INT( x + rx )
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = INT( y + ry )
            if (.not.wally .and. iy.eq.Ly) iy = 0
            icell = 1 + ix + iy*(Lx+1)
            listtemp(i) = head_priv(icell,this_thread)
            head_priv(icell,this_thread) = i
         enddo
C$OMP end do
C$OMP end parallel
         j = 0
         j2 = -2
         do icell = 1, ncell
            empty = .true.
            ithread = 0
 10         ithread = ithread+1
            if (ithread.gt.numthreads) goto 30
            i = head_priv(icell,ithread)
 20         if (i.eq.0) goto 10
            j = j+1
            j2 = j2+2
            i2 = 2*(i-1)
            do m = 1,2
               R(j2+m) = Rtemp(i2+m)
               V(j2+m) = Vtemp(i2+m)
            enddo
            if (empty) then
               head(icell) = j
               empty = .false.
            endif
            list(j) = j+1
            i = listtemp(i)
            goto 20
 30         if (.not.empty) list(j) = 0
         enddo
      else !do not sort SRD coordinates
C$OMP parallel 
C$OMP& default(shared)
C$OMP& private(this_thread,icell,i,i2,ix,iy,empty,ithread,j)
         numthreads = omp_get_num_threads()
         this_thread = omp_get_thread_num() + 1
         do icell = 1, ncell
            head_priv(icell,this_thread) = 0
         enddo
C$OMP do
         do i = 1, N
            i2 = 2*(i-1)
            ix = INT( R(i2+1) + rx )
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = INT( R(i2+2) + ry )
            if (.not.wally .and. iy.eq.Ly) iy = 0
            icell = 1 + ix + iy*(Lx+1)
            list(i) = head_priv(icell,this_thread)
            head_priv(icell,this_thread) = i
         enddo
C$OMP end do
C$OMP barrier
C$OMP do
         do icell = 1, ncell
            empty = .true.
            do ithread = 1,numthreads
               i = head_priv(icell,ithread)
               if (i.ne.0) then
                  if (empty) then
                     head(icell) = i
                     empty = .false.
                  else
                     list(j) = i
                  endif
 70               j = i
                  i = list(i)
                  if (i.ne.0) goto 70
               endif
            enddo !ithread
         enddo
C$OMP end do
C$OMP end parallel
      endif
#else
      if (mod(step,isort).eq.0) then !sort SRD coordinates every isort
         do i = 1,N
            i2 = 2*(i-1)
            x = R(i2+1)
            y = R(i2+2)
            Rtemp(i2+1) = x
            Rtemp(i2+2) = y
            do m = 1,2
               Vtemp(i2+m) = V(i2+m)
            enddo
            ix = INT( x + rx )
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = INT( y + ry )
            if (.not.wally .and. iy.eq.Ly) iy = 0
            icell = 1 + ix + iy*(Lx+1)
            listtemp(i) = head(icell)
            head(icell) = i
         enddo
         j = 0
         j2 = -2
         do icell = 1, ncell
            empty = .true.
            i = head(icell)
 120        if (i.ne.0) then
               j = j+1
               j2 = j2+2
               i2 = 2*(i-1)
               do m = 1,2
                  R(j2+m) = Rtemp(i2+m)
                  V(j2+m) = Vtemp(i2+m)
               enddo
               if (empty) then
                  head(icell) = j
                  empty = .false.
               endif
               list(j) = j+1
               i = listtemp(i)
               goto 120
            endif
            if (.not.empty) list(j) = 0
         enddo
      else
         do i = 1, N
            i2 = 2*i
            ix = INT( R(i2-1) + rx )
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = INT( R(i2)   + ry )
            if (.not.wally .and. iy.eq.Ly) iy = 0
            icell = 1 + ix + iy*(Lx+1)
            list(i) = head(icell)
            head(icell) = i
         enddo
      endif
#endif
      do i = 1, Ndepl
         i2 = 2*i
         ix = INT( Rdepl(i2-1) + rx )
         if (.not.wallx .and. ix.eq.Lx) ix = 0
         iy = INT( Rdepl(i2) + ry )
         if (.not.wally .and. iy.eq.Ly) iy = 0
         icell = 1 + ix + iy*(Lx+1)
         listdepl(i) = headdepl(icell)
         headdepl(icell) = i
      enddo

      ncelldepl = Lxdepl*Lydepl
      do i = 1, ncelldepl
         bigheaddepl(i) = 0
      enddo
      if (V0p.eq.0.d0) then
         do i = 1,Ndepl
            biglistdepl(i) = 0
         enddo
      else
         rLxinv = dble(Lxdepl)/dble(Lx)
         rLyinv = dble(Lydepl)/dble(Ly)
         do i = 1,Ndepl
            i2 = 2*i
            ix = INT( rLxinv*Rdepl(i2-1) )
            if (ix.eq.Lxdepl) ix = 0
            iy = INT( rLyinv*Rdepl(i2) )
            if (iy.eq.Lydepl) iy = 0
            icell = 1 + ix + iy*Lxdepl
            biglistdepl(i) = bigheaddepl(icell)
            bigheaddepl(icell) = i
         enddo
      endif

      ncellcol = Lxcol*Lycol
      do i = 1, ncellcol
         headcol(i) = 0
      enddo
      rLxinv = dble(Lxcol)/dble(Lx)
      rLyinv = dble(Lycol)/dble(Ly)
      do i = 1, Ncol
         i2 = 2*i
         ix = INT( rLxinv*Rcol(i2-1) )
         if (ix.eq.Lxcol) ix = 0
         iy = INT( rLyinv*Rcol(i2) )
         if (iy.eq.Lycol) iy = 0
         icell = 1 + ix + iy*Lxcol
         listcol(i) = headcol(icell)
         headcol(icell) = i
      enddo

      return
      end


      subroutine force(N,Ndepl,Ncol,Ncolup,Lx,Ly,wallx,wally,
     :  Lxdepl,Lydepl,Lxcol,Lycol,
     :  sigcc,kappa,sigcp,sigcs,Rg,rlist,rcutcc,rcutccy,rcutcp,rcutcs,
     :  rcutpp,Mdepl,V0,V0y,V0p,Vpot,grav,direct,rx,ry,
     :  fixcom,fixsolv,fixdepl,fixcol,Foverlap)

      integer N,Ndepl,Ncol,Ncolup
      integer Lx,Ly,Lxdepl,Lydepl,Lxcol,Lycol
      integer ixmin,ixmax,iymin,iymax,ix,iy
      integer icell,icellcol,jcell,jcell0,nabor,ncelldepl,ncellcol
      integer i,i2,j,j2
      double precision  sigcc,kappa,sigcp,sigcs,Rg,rcutcc,rcutccy,
     :                  rcutcp,rcutcs
      double precision  rcutpp,rlist,rlistsq,grav,direct,Foverlap
      double precision  sigccsq,sigcssq,sigcpsq,Rgsq,Rgsqinv
      double precision  rcutccsq,rcutccysq,rcutcpsq,rcutcssq,rcutppsq
      double precision  Mdepl,Mtotinv,V0,V0y,V0p,Vpot
      double precision  xi,yi,xij,yij,rij,rijsq
      double precision  invrij,rx,ry
      double precision  rLx,rLy,rLx2,rLy2,dV,r6,r12,r24,r48,f0
      double precision  fxtot,fytot,fxtotdepl,fytotdepl
      double precision  fxi,fyi
      logical fixcom,fixsolv,fixdepl,fixcol,wallx,wally

      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      common / block3 / sxyfr,sxyrot,sxykin
      double precision sxyfr,sxyrot,sxykin
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) :: Fdrx,Fdry,Fdrxsq,Fdrysq

C ### Initialisation

      sigccsq = sigcc**2
      sigcpsq = sigcp**2
      sigcssq = sigcs**2
      Rgsq    = Rg**2
      Rgsqinv = 1.5d0/Rgsq
      rlistsq = rlist**2
      rcutccsq = rcutcc**2
      rcutccysq = rcutccy**2
      rcutcpsq = rcutcp**2
      rcutcssq = rcutcs**2
      rcutppsq = rcutpp**2
      rLx = dble(Lx)
      rLy = dble(Ly)
      rLx2 = rLx/2.d0
      rLy2 = rLy/2.d0
      Vpot = 0.d0
      sxyfr = 0.d0
      sxyrot = 0.d0
      do i = 1,2*N
         F(i) = 0.d0
      enddo
      do i = 1,2*Ndepl
         Fdepl(i) = 0.d0
      enddo
      do i = 1,2*Ncol
         Fcol(i) = 0.d0
      enddo
      if (grav.gt.0.d0) then !external field on solvent
         do i = 1,N
            F(2*i) = direct*grav
         enddo
         if (fixcom) then !opposing field on colloids
            do i = 1,Ncol
               Fcol(2*i) = -direct*grav*dble(N)/dble(Ncol)
            enddo
         endif
      endif
      if (grav.lt.0.d0) then !external field on colloids
Csed
C         do i = 1,Ncol/2
C            Fcol(2*i) = -direct*grav
C         enddo
C         do i = Ncol/2+1,Ncol
C            Fcol(2*i) = direct*grav
C         enddo
         do i = 1,Ncolup
            Fcol(2*i) = -direct*grav
         enddo
         do i = Ncolup+1,Ncol
            Fcol(2*i) = direct*grav
         enddo
         if (fixcom) then !opposing field on solvent
            do i = 1,N
               F(2*i) = dble(2*Ncolup-Ncol)*direct*grav/dble(N)
            enddo
         endif
Csed
      endif
      fxtot = 0.d0 !total force on solvent
      fytot = 0.d0
      fxtotdepl = 0.d0 !total force on depletant
      fytotdepl = 0.d0

      ncellcol = Lxcol*Lycol
      do 9000 icellcol = 1, ncellcol

         i = headcol(icellcol)
 6000    if (i.gt.0) then
            i2 = 2*(i-1)
            xi = Rcol(i2+1)
            yi = Rcol(i2+2)
            fxi = Fcol(i2+1) 
            fyi = Fcol(i2+2) 

C ### Colloid - colloid interactions

            if (fixcol) goto 1000

            j = listcol(i)

 6100       if (j.gt.0) then
               j2 = 2*(j-1)
               xij = xi - Rcol(j2+1)
               yij = yi - Rcol(j2+2)
               if (.not.wallx) then
                  if (xij.gt.rLx2) then
                     xij = xij-rLx
                  elseif (xij.le.-rLx2) then
                     xij = xij+rLx
                  endif
               endif
               if (.not.wally) then
                  if (yij.gt.rLy2) then
                     yij = yij-rLy
                  elseif (yij.le.-rLy2) then
                     yij = yij+rLy
                  endif
               endif
               rijsq = xij**2 + yij**2
               rij = dsqrt(rijsq)
               if (rijsq.lt.rcutccsq) then
                  r24 = (sigccsq/rijsq)**12
                  r48 = r24**2
                  dV = 4.d0*V0*(r48-r24+0.25d0)
                  Vpot = Vpot+dV
                  f0 = 96.d0*V0*(2.d0*r48-r24)/rijsq
                  fxi = fxi+f0*xij
                  fyi = fyi+f0*yij
                  Fcol(j2+1) = Fcol(j2+1)-f0*xij
                  Fcol(j2+2) = Fcol(j2+2)-f0*yij
                  sxyfr = sxyfr + f0*xij*yij
                  if (rij.lt.2.10d0*sigcs) then !compensate overlap f.
                     f0 = Foverlap*dsqrt(4.d0-(rij/(1.05d0*sigcs))**2)
                     fxi = fxi+f0*xij
                     fyi = fyi+f0*yij
                     Fcol(j2+1) = Fcol(j2+1)-f0*xij
                     Fcol(j2+2) = Fcol(j2+2)-f0*yij
                     sxyfr = sxyfr + f0*xij*yij
                  endif
               endif
               if (V0y.ne.0.d0 .and. rijsq.lt.rcutccysq) then
                  invrij = 1.d0/rij
                  dV = V0y*sigcc*invrij*dexp(-kappa*(rij-sigcc))
                  Vpot = Vpot+dV
                  f0 = dV*(kappa+invrij)*invrij
                  fxi = fxi+f0*xij
                  fyi = fyi+f0*yij
                  Fcol(j2+1) = Fcol(j2+1)-f0*xij
                  Fcol(j2+2) = Fcol(j2+2)-f0*yij
                  sxyfr = sxyfr + f0*xij*yij
               endif
               j = listcol(j)
               goto 6100
            endif

            jcell0 = 4*(icellcol-1)

            do 6300 nabor = 1,4

               jcell = mapcol(jcell0+nabor)
               j = headcol(jcell)
 
 6200          if (j.gt.0) then
                  j2 = 2*(j-1)
                  xij = xi - Rcol(j2+1)
                  yij = yi - Rcol(j2+2)
                  if (.not.wallx) then
                     if (xij.gt.rLx2) then
                        xij = xij-rLx
                     elseif (xij.le.-rLx2) then
                        xij = xij+rLx
                     endif
                  endif
                  if (.not.wally) then
                     if (yij.gt.rLy2) then
                        yij = yij-rLy
                     elseif (yij.le.-rLy2) then
                        yij = yij+rLy
                     endif
                  endif
                  rijsq = xij**2 + yij**2
                  rij = dsqrt(rijsq)
                  if (rijsq.lt.rcutccsq) then
                     r24 = (sigccsq/rijsq)**12
                     r48 = r24**2
                     dV = 4.d0*V0*(r48-r24+0.25d0)
                     Vpot = Vpot+dV
                     f0 = 96.d0*V0*(2.d0*r48-r24)/rijsq
                     fxi = fxi+f0*xij
                     fyi = fyi+f0*yij
                     Fcol(j2+1) = Fcol(j2+1)-f0*xij
                     Fcol(j2+2) = Fcol(j2+2)-f0*yij
                     sxyfr = sxyfr + f0*xij*yij
                     if (rij.lt.2.10d0*sigcs) then !compensate overlap f
                        f0=Foverlap*dsqrt(4.d0-(rij/(1.05d0*sigcs))**2)
                        fxi = fxi+f0*xij
                        fyi = fyi+f0*yij
                        Fcol(j2+1) = Fcol(j2+1)-f0*xij
                        Fcol(j2+2) = Fcol(j2+2)-f0*yij
                        sxyfr = sxyfr + f0*xij*yij
                     endif
                  endif
                  if (V0y.ne.0.d0 .and. rijsq.lt.rcutccysq) then
                     invrij = 1.d0/rij
                     dV = V0y*sigcc*invrij*dexp(-kappa*(rij-sigcc))
                     Vpot = Vpot+dV
                     f0 = dV*(kappa+invrij)*invrij
                     fxi = fxi+f0*xij
                     fyi = fyi+f0*yij
                     Fcol(j2+1) = Fcol(j2+1)-f0*xij
                     Fcol(j2+2) = Fcol(j2+2)-f0*yij
                     sxyfr = sxyfr + f0*xij*yij
                  endif
                  j = listcol(j)
                  goto 6200
               endif

 6300       continue !next nabor

C ### Colloid - depletant interactions

 1000       if (Ndepl.eq.0) goto 1100
            ixmin = NINT(xi-rcutcp+rx-0.5d0)
            ixmax = NINT(xi+rcutcp+rx-0.5d0)
            if (wallx) then
               if (ixmin.lt.0) ixmin = 0
               if (ixmax.ge.Lx) ixmax = Lx
            else
               if (ixmin.lt.0) ixmin = ixmin+Lx
               if (ixmax.ge.Lx) ixmax = ixmax-Lx
            endif
            iymin = NINT(yi-rcutcp+ry-0.5d0)
            iymax = NINT(yi+rcutcp+ry-0.5d0)
            if (wally) then
               if (iymin.lt.0) iymin = 0
               if (iymax.ge.Ly) iymax = Ly
            else
               if (iymin.lt.0) iymin = iymin+Ly
               if (iymax.ge.Ly) iymax = iymax-Ly
            endif

            ix = ixmin-1
 1010       ix = ix+1
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = iymin-1
 1020       iy = iy+1
            if (.not.wally .and. iy.eq.Ly) iy = 0

            icell = 1 + ix + iy*(Lx+1)
            j = headdepl(icell)
 1040       if (j.gt.0) then
               j2 = 2*(j-1)
               xij = xi - Rdepl(j2+1)
               yij = yi - Rdepl(j2+2)
               if (.not.wallx) then
                  if (xij.gt.rLx2) then
                     xij = xij-rLx
                  elseif (xij.le.-rLx2) then
                     xij = xij+rLx
                  endif
               endif
               if (.not.wally) then
                  if (yij.gt.rLy2) then
                     yij = yij-rLy
                  elseif (yij.le.-rLy2) then
                     yij = yij+rLy
                  endif
               endif
               rijsq = xij**2 + yij**2
               if (rijsq.lt.rcutcpsq) then
                  r6 = (sigcpsq/rijsq)**3
                  r12 = r6**2
                  dV = 4.d0*V0*(r12-r6+0.25d0)
                  Vpot = Vpot+dV
                  f0 = 24.d0*V0*(2.d0*r12-r6)/rijsq
                  fxi = fxi+f0*xij
                  fyi = fyi+f0*yij
                  Fdrx(i) = Fdrx(i)+f0*xij
                  Fdry(i) = Fdry(i)+f0*yij
                  Fdrxsq(i) = Fdrxsq(i)+f0*f0*xij*xij
                  Fdrysq(i) = Fdrysq(i)+f0*f0*yij*yij
                  Fdepl(j2+1) = Fdepl(j2+1)-f0*xij
                  Fdepl(j2+2) = Fdepl(j2+2)-f0*yij
                  fxtotdepl = fxtotdepl - f0*xij
                  fytotdepl = fytotdepl - f0*yij
                  sxyfr = sxyfr + f0*xij*yij
               endif
               j = listdepl(j)
               goto 1040
            endif !j>0

            if (iy.ne.iymax) goto 1020
            if (ix.ne.ixmax) goto 1010

C ### Colloid - solvent interactions

1100        if (N.eq.0) goto 1190
            ixmin = NINT(xi-rcutcs+rx-0.5d0)
            ixmax = NINT(xi+rcutcs+rx-0.5d0)
            if (wallx) then
               if (ixmin.lt.0) ixmin = 0
               if (ixmax.ge.Lx) ixmax = Lx
            else
               if (ixmin.lt.0) ixmin = ixmin+Lx
               if (ixmax.ge.Lx) ixmax = ixmax-Lx
            endif
            iymin = NINT(yi-rcutcs+ry-0.5d0)
            iymax = NINT(yi+rcutcs+ry-0.5d0)
            if (wally) then
               if (iymin.lt.0) iymin = 0
               if (iymax.ge.Ly) iymax = Ly
            else
               if (iymin.lt.0) iymin = iymin+Ly
               if (iymax.ge.Ly) iymax = iymax-Ly
            endif

            ix = ixmin-1
 1110       ix = ix+1
            if (.not.wallx .and. ix.eq.Lx) ix = 0
            iy = iymin-1
 1120       iy = iy+1
            if (.not.wally .and. iy.eq.Ly) iy = 0

            icell = 1 + ix + iy*(Lx+1)

            j = head(icell)
 1140       if (j.gt.0) then
               j2 = 2*(j-1)
               xij = xi - R(j2+1)
               yij = yi - R(j2+2)
               if (.not.wallx) then
                  if (xij.gt.rLx2) then
                     xij = xij-rLx
                  elseif (xij.le.-rLx2) then
                     xij = xij+rLx
                  endif
               endif
               if (.not.wally) then
                  if (yij.gt.rLy2) then
                     yij = yij-rLy
                  elseif (yij.le.-rLy2) then
                     yij = yij+rLy
                  endif
               endif
               rijsq = xij**2 + yij**2
               if (rijsq.lt.rcutcssq) then
                  r6 = (sigcssq/rijsq)**3
                  r12 = r6**2
                  dV = 4.d0*V0*(r12-r6+0.25d0)
                  Vpot = Vpot+dV
                  f0 = 24.d0*V0*(2.d0*r12-r6)/rijsq
                  fxi = fxi+f0*xij
                  fyi = fyi+f0*yij
                  Fdrx(Ncol+i) = Fdrx(Ncol+i)+f0*xij
                  Fdry(Ncol+i) = Fdry(Ncol+i)+f0*yij
                  Fdrxsq(Ncol+i) = Fdrxsq(Ncol+i)+f0*f0*xij*xij
                  Fdrysq(Ncol+i) = Fdrysq(Ncol+i)+f0*f0*yij*yij
                  F(j2+1) = F(j2+1)-f0*xij
                  F(j2+2) = F(j2+2)-f0*yij
                  fxtot = fxtot - f0*xij
                  fytot = fytot - f0*yij
                  sxyfr = sxyfr + f0*xij*yij
               endif
               j = list(j)
               goto 1140
            endif !j>0

            if (iy.ne.iymax) goto 1120
            if (ix.ne.ixmax) goto 1110

 1190       Fcol(i2+1) = fxi
            Fcol(i2+2) = fyi

            i = listcol(i)
            goto 6000 !next colloid in current icellcol
         endif

 9000 continue !next icellcol

C ### Depletant - depletant interactions

      if (Ndepl.eq.0 .or. V0p.eq.0) goto 1400
      ncelldepl = Lxdepl*Lydepl
      do 5000 icell = 1,ncelldepl
         i = bigheaddepl(icell)
 1900    if (i.gt.0) then
            i2 = 2*(i-1)
            xi = Rdepl(i2+1)
            yi = Rdepl(i2+2)
            fxi = Fdepl(i2+1)
            fyi = Fdepl(i2+2)
            j = biglistdepl(i)
 2000       if (j.gt.0) then !consider this pair
               j2 = 2*(j-1)
               xij = xi - Rdepl(j2+1)
               yij = yi - Rdepl(j2+2)
               if (.not.wallx) then
                  if (xij.gt.rLx2) then
                     xij = xij-rLx
                  elseif (xij.le.-rLx2) then
                     xij = xij+rLx
                  endif
               endif
               if (.not.wally) then
                  if (yij.gt.rLy2) then
                     yij = yij-rLy
                  elseif (yij.le.-rLy2) then
                     yij = yij+rLy
                  endif
               endif
               rijsq = xij**2 + yij**2
               if (rijsq.lt.rcutppsq) then
                  dV = V0p*dexp(-0.75d0*rijsq/Rgsq)
                  Vpot = Vpot+dV
                  f0 = dV*Rgsqinv
                  fxi = fxi+f0*xij
                  fyi = fyi+f0*yij
                  Fdepl(j2+1) = Fdepl(j2+1)-f0*xij
                  Fdepl(j2+2) = Fdepl(j2+2)-f0*yij
                  sxyfr = sxyfr + f0*xij*yij
               endif
               j = biglistdepl(j)
               goto 2000
            endif !last j in current cell reached

            jcell0 = 4*(icell-1)
            do 4000 nabor = 1,4
               jcell = map(jcell0+nabor) 
               j = bigheaddepl(jcell)
 3000          if (j.ne.0) then
                  j2 = 2*(j-1)
                  xij = xi - Rdepl(j2+1)
                  yij = yi - Rdepl(j2+2)
                  if (.not.wallx) then
                     if (xij.gt.rLx2) then
                        xij = xij-rLx 
                     elseif (xij.le.-rLx2) then
                        xij = xij+rLx
                     endif
                  endif
                  if (.not.wally) then
                     if (yij.gt.rLy2) then
                        yij = yij-rLy 
                     elseif (yij.le.-rLy2) then
                        yij = yij+rLy
                     endif
                  endif
                  rijsq = xij**2 + yij**2
                  if (rijsq.lt.rcutppsq) then
                     dV = V0p*dexp(-0.75*rijsq/Rgsq)
                     Vpot = Vpot+dV
                     f0 = dV*Rgsqinv
                     fxi = fxi+f0*xij
                     fyi = fyi+f0*yij
                     Fdepl(j2+1) = Fdepl(j2+1)-f0*xij
                     Fdepl(j2+2) = Fdepl(j2+2)-f0*yij
                     sxyfr = sxyfr + f0*xij*yij
                  endif
                  j = biglistdepl(j)
                  goto 3000
               endif !last j in current cell reached
 4000       continue !next neighbouring cell
            Fdepl(i2+1) = fxi
            Fdepl(i2+2) = fyi
            i = biglistdepl(i)
            goto 1900
         endif !last i in current cell reached
 5000 continue !next cell

C ### Depletant - solvent interactions by SRD interactions

C ### Solvent - solvent interaction by SRD interactions

C ### Subtract nett forces to keep various centres-of-mass fixed

 1400 if (fixdepl) then !fix depletant c.o.m.
         if (.not.wallx) then
            do i = 1,Ndepl
               Fdepl(2*i-1) = Fdepl(2*i-1) - fxtotdepl/dble(Ndepl)
            enddo
         endif
         if (.not.wally) then
            do i = 1,Ndepl
               Fdepl(2*i)   = Fdepl(2*i)   - fytotdepl/dble(Ndepl)
            enddo
         endif
      endif
      if (fixsolv) then !fix solvent c.o.m.
         if (.not.wallx) then
            do i = 1,N
               F(2*i-1) = F(2*i-1) - fxtot/dble(N)
            enddo
         endif
         if (.not.wally) then
            do i = 1,N
               F(2*i)   = F(2*i)   - fytot/dble(N)
            enddo
         endif
      endif
      if (fixcol) then !fix EACH colloid individually
C        Calculate total force on all colloids and
C        remove individual forces on the colloids  
         fxi = 0.d0
         fyi = 0.d0
         do i = 1, Ncol
            fxi = fxi + Fcol(2*i-1) 
            Fcol(2*i-1) = 0.d0
            fyi = fyi + Fcol(2*i) 
            Fcol(2*i)   = 0.d0
         enddo
C        Redistribute total force on colloids over all other particles
C        to keep TOTAL (solvent plus polymers) c.o.m. fixed
C         Mtotinv = 1.d0/(dble(N)+Mdepl*dble(Ndepl))
C         do i = 1,N
C            F(2*i-1) = F(2*i-1)+Mtotinv*fxi
C            F(2*i)   = F(2*i)  +Mtotinv*fyi
C         enddo
C         do i = 1,Ndepl
C            Fdepl(2*i-1) = Fdepl(2*i-1)+Mdepl*Mtotinv*fxi
C            Fdepl(2*i)   = Fdepl(2*i)  +Mdepl*Mtotinv*fyi
C         enddo
      endif

      return
      end


      subroutine momentum(N,Ndepl,Ncol,Mdepl,Mcol,ptot)

      integer N,Ndepl,Ncol,i,i2,j
      double precision Mdepl,Mcol,ptot(2)
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      do j = 1,2
         ptot(j) = 0.d0
      enddo
      do i = 1,N
         i2 = 2*(i-1)
         do j = 1,2
            ptot(j) = ptot(j)+V(i2+j)
         enddo
      enddo
      do i = 1,Ndepl
         i2 = 2*(i-1)
         do j = 1,2
            ptot(j) = ptot(j)+Mdepl*Vdepl(i2+j)
         enddo
      enddo
      do i = 1,Ncol
         i2 = 2*(i-1)
         do j = 1,2
            ptot(j) = ptot(j)+Mcol*Vcol(i2+j)
         enddo
      enddo

      return
      end


      subroutine kinetic(N,Lx,Ly,ihydro,vflow,kin,temp)

C calculate kinetic energy of solvent

      integer N,Lx,Ly,ncell,ihydro,i,i2,j,icell
      double precision  kin,temp,free,vflow,Mtot,Vav(2)
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol

      kin = 0.d0
      ncell = (Lx+1)*(Ly+1)

      if (ihydro.eq.0) then

         free = 2.d0*dble(N-1)
         do i = 1,N
            kin = kin + V(2*i-1)**2 + (V(2*i)-vflow)**2
         enddo

      else

      free = 2.d0*dble(N-ncell) !each cell represents one constraint
      do 2000 icell = 1,ncell
C calculate average velocity of this cell
         do j = 1,2
            Vav(j) = 0.d0
         enddo
         Mtot = 0.d0
         i = head(icell)
 1000    if (i.gt.0) then
            i2 = 2*(i-1)
            do j = 1,2
               kin = kin + (V(i2+j))**2
               Vav(j) = Vav(j)+V(i2+j)
            enddo
            Mtot = Mtot + 1.d0
            i = list(i)
            goto 1000
         endif
C calculate kinetic energy relative to average velocity
         if (Mtot.ne.0.d0) then
            do j = 1,2
               Vav(j) = Vav(j)/Mtot
               kin = kin - Mtot*(Vav(j)**2)
            enddo
         else !no solvent in this cell, increase 'free' by 2
            free = free+2.d0
         endif
 2000 continue

      endif

      temp = kin/free
      kin = 0.5d0*kin

      return
      end


      subroutine collide(N,Ndepl,Ncol,Mdepl,Mcol,Lx,Ly,
     :                   wallx,wally,avdens,
     :                   ihydro,iflow,vflow,rx,ry,dt,step,ssample,
     :                   fixcom,fixsolv,fixdepl,itherm,kin,temp,icol)

      integer N,Ndepl,Ncol,i,i2,j,k,icell,itherm,icell2,istart
      integer Lx,Ly,ihydro,iflow,step,ssample,icol,ncell
      double precision  Mdepl,Mcol,det,gauss,dummy,vflow,rx,ry,dt,
     :                  kin,temp,free,avdens
      double precision  pi,az,randy,vv,vzav,fact
      double precision  vx
      double precision  vtemp(2),vadd(2)
      double precision,pointer,dimension(:) :: Vprev,Vprevdepl
      logical,pointer,dimension(:) :: done
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      common / block2a / nearwall
      logical,pointer,dimension(:) :: nearwall
      common / block3 / sxyfr,sxyrot,sxykin
      double precision sxyfr,sxyrot,sxykin
      double precision Vav(2),Mtot,Vtotsolv(2),Vtotdepl(2)
      double precision,pointer,dimension(:) :: Vavsolv
      double precision,pointer,dimension(:) :: randy1
      double precision,pointer,dimension(:,:) :: gauss1
      logical fixcom,fixsolv,fixdepl,wallx,wally

      if (ssample.ne.0 .and. mod(step,ssample).eq.0) then
         allocate(Vprev(2*N),Vprevdepl(2*Ndepl))
         do i = 1,2*N
            Vprev(i) = V(i)
         enddo
         do i = 1,2*Ndepl
            Vprevdepl(i) = Vdepl(i)
         enddo
      endif

      kin = 0.d0
      do i = 1,2
         Vtotsolv(i) = 0.d0
         Vtotdepl(i) = 0.d0
      enddo
      pi = 4.d0*datan(1.d0)
      ncell = (Lx+1)*(Ly+1)
      if (ihydro.eq.0) then
         free = 2.d0*dble(N-1)
      else
         free = 2.d0*dble(N-ncell)
      endif
      allocate(Vavsolv(2*ncell))
      allocate(randy1(ncell))
      allocate(gauss1(2,ncell))

C Predetermine random numbers for use in parallel loop
      do icell = 1,ncell
         randy1(icell) = randy()
         if (nearwall(icell)) then
            do j = 1,2
               gauss1(j,icell) = gauss(dummy)
            enddo
         endif
      enddo 

C$OMP parallel
C$OMP& default(shared)
C$OMP& private(icell,icell2,j,Vav,Mtot,i,i2,vx)
C$OMP& private(az,vadd,vtemp)
C$OMP& reduction(+:Vtotsolv,Vtotdepl,kin,free)
C$OMP do
      do 2000 icell = 1,ncell
         icell2 = 2*(icell-1)
         do j = 1,2
            Vav(j) = 0.d0
            Vavsolv(icell2+j) = 0.d0
         enddo
C calculate average velocity
         Mtot = 0.d0
         i = head(icell)
 1000    if (i.gt.0) then
            i2 = 2*(i-1)
            do j = 1,2
               vx = V(i2+j)
               Vavsolv(icell2+j) = Vavsolv(icell2+j) + vx
               Vav(j) = Vav(j) + vx
            enddo
            Mtot = Mtot + 1.d0
            i = list(i)
            goto 1000
         endif
         if (Mtot.eq.0.d0) then
            free = free+2.d0
            goto 999
         endif
         do j = 1,2
            Vavsolv(icell2+j) = Vavsolv(icell2+j)/Mtot
         enddo
C Add gaussian noise due to finite temperature of walls
         if (nearwall(icell) .and. Mtot.lt.avdens) then
            do j = 1,2
               Vav(j) = Vav(j) +
     :            dsqrt(avdens-Mtot)*gauss1(j,icell)
            enddo
            Mtot = avdens
         endif
 999     i = headdepl(icell)
 1001    if (i.gt.0) then
            i2 = 2*(i-1)
            do j = 1,2
               Vav(j) = Vav(j) + Mdepl*Vdepl(i2+j)
            enddo
            Mtot = Mtot + Mdepl
            i = listdepl(i)
            goto 1001
         endif
         if (Mtot.eq.0.d0) goto 2000
         do j = 1,2
            Vav(j) = Vav(j)/Mtot
         enddo
C choose random axis
         az = randy1(icell)
C rotate (90 degrees for now)
C Rotate each velocity relative to average (c.o.m.) velocity in cell
         i = head(icell)
10       if (i.gt.0) then
            i2 = 2*(i-1)
            do j = 1,2
               vtemp(j) = V(i2+j)-Vav(j)
            enddo
            if (az.le.0.5d0) then !-90 degrees
               vadd(1) = -vtemp(2)
               vadd(2) = vtemp(1)
            else !+90 degrees
               vadd(1) = vtemp(2)
               vadd(2) = -vtemp(1)
            endif
            do j = 1,2
               V(i2+j) = vadd(j) + Vav(j)
               Vtotsolv(j) = Vtotsolv(j)+V(i2+j)
            enddo
            if (ihydro.eq.0) then !kinetic energy from abs. vel.
               kin = kin+V(i2+1)**2 + (V(i2+2)-vflow)**2
            else
               do j = 1,2
                  kin = kin + vadd(j)**2
               enddo
            endif
            i = list(i)
            goto 10
         endif !i>0
         i = headdepl(icell)
11       if (i.gt.0) then
            i2 = 2*(i-1)
            do j = 1,2
               vtemp(j) = Vdepl(i2+j)-Vav(j)
            enddo
            if (az.le.0.5d0) then !-90 degrees
               vadd(1) = -vtemp(2)
               vadd(2) = vtemp(1)
            else !+90 degrees
               vadd(1) = vtemp(2)
               vadd(2) = -vtemp(1)
            endif
            do j = 1,2
               Vdepl(i2+j) = vadd(j)+Vav(j)
               Vtotdepl(j) = Vtotdepl(j)+Vdepl(i2+j)
            enddo
            i = listdepl(i)
            goto 11
         endif !i>0
 2000 continue !icell
C$OMP end do
C$OMP end parallel
      kin = 0.5d0*kin

      deallocate(randy1)
      deallocate(gauss1)

C Remove net momentum of solvent if c.o.m. needs to be fixed.
      if (fixsolv) then
         do j = 1,2
            Vtotsolv(j) = Vtotsolv(j)/dble(N)
         enddo
         do i = 1,N
            i2 = 2*(i-1)
            if (.not.wallx) V(i2+1) = V(i2+1)-Vtotsolv(1)
            if (.not.wally) V(i2+2) = V(i2+2)-Vtotsolv(2)
         enddo
         if (ihydro.eq.0) then
            if (.not.wallx) kin = kin-0.5d0*dble(N)*(Vtotsolv(1)**2)
            if (.not.wally) kin = kin-0.5d0*dble(N)*(Vtotsolv(2)**2)
         endif
      else if (fixcom) then !fix total c.o.m.
         do i = 1,Ncol
            i2 = 2*(i-1)
            do j = 1,2
               Vtotsolv(j) = Vtotsolv(j) + Mcol*Vcol(i2+j)
            enddo
         enddo
         do j = 1,2
            Vtotsolv(j) = Vtotsolv(j)/(dble(N)+Mcol*dble(Ncol))
         enddo
         do i = 1,N
            i2 = 2*(i-1)
            if (.not.wallx) V(i2+1) = V(i2+1)-Vtotsolv(1)
            if (.not.wally) V(i2+2) = V(i2+2)-Vtotsolv(2)
         enddo
         do i = 1,Ncol
            i2 = 2*(i-1)
            if (.not.wallx) Vcol(i2+1) = Vcol(i2+1)-Vtotsolv(1)
            if (.not.wally) Vcol(i2+2) = Vcol(i2+2)-Vtotsolv(2)
         enddo
         if (ihydro.eq.0) then
            if (.not.wallx) kin = kin-0.5d0*dble(N)*(Vtotsolv(1)**2)
            if (.not.wally) kin = kin-0.5d0*dble(N)*(Vtotsolv(2)**2)
         endif
      endif !fixsolv/fixcom
C Remove net momentum of depletant if c.o.m. needs to be fixed
      if (fixdepl) then
         do j = 1,2
            Vtotdepl(j) = Vtotdepl(j)/dble(Ndepl)
         enddo
         do i = 1,Ndepl
            i2 = 2*(i-1)
            if (.not.wallx) Vdepl(i2+1) = Vdepl(i2+1)-Vtotdepl(1)
            if (.not.wally) Vdepl(i2+2) = Vdepl(i2+2)-Vtotdepl(2)
         enddo
      endif

C Hydrodynamics in the solvent is turned off if velocities of
C the solvent particles are interchanged randomly.
      if (mod(N,2).ne.0) stop 'N must be even'
      if (ihydro.eq.0) then
         allocate(done(N))
         do i = 1, N
            done(i) = .false.
         enddo
         do 100 istart = 1, N, 20
            do 200 i = istart, istart+19
               if (istart.gt.N) goto 300
               if (done(i)) goto 200
30             j = istart+INT(20*randy())
               if (j.gt.N) goto 30
               if (done(j)) goto 30
               vtemp(1) = V(2*i-1)
               vtemp(2) = V(2*i)
               V(2*i-1) = V(2*j-1)
               V(2*i)   = V(2*j)
               V(2*j-1) = vtemp(1)
               V(2*j)   = vtemp(2)
               done(i)  = .true.
               done(j)  = .true.
200         continue
100      continue
300      continue
         deallocate(done)
      endif

C Alternatively, hydrodynamics in the solvent is turned off if
C velocities are redrawn from a Maxwellian distribution (bath)
C every |ihydro| collision intervals. Note: this also acts as
C a thermostat, so separate thermostating will not be needed.
 
      if (ihydro.le.-2 .and. mod(step,abs(icol*ihydro)).eq.0) then
         call anderson(N)
      endif

C Thermostat the solvent if desired
      temp = 2.d0*kin/free
      if (itherm.eq.1) then
         fact = dsqrt(1.d0/temp)
         if (ihydro.eq.0) then !scale abs. vel.
            do i = 1, N
               V(2*i-1) = fact*V(2*i-1)
               V(2*i) = vflow+fact*(V(2*i)-vflow)
            enddo
         else !scale rel. vel.
C$OMP parallel
C$OMP& default(shared)
C$OMP& private(icell,icell2,vtemp,i,i2,j)
C$OMP do
            do 4000 icell = 1,ncell
               icell2 = 2*(icell-1)
               do j = 1,2
                  vtemp(j) = Vavsolv(icell2+j)
               enddo
               i = head(icell)
 3000          if (i.gt.0) then
                  i2 = 2*(i-1)
                  do j = 1,2
                     V(i2+j) = vtemp(j) + fact*(V(i2+j)-vtemp(j))
                  enddo
                  i = list(i)
                  goto 3000
               endif
 4000       continue
C$OMP end do
C$OMP end parallel
         endif
      endif !itherm=1
      deallocate(Vavsolv)

C Calculate rotational part of stress tensor due to shifting of cells
      if (ssample.ne.0 .and. mod(step,ssample).eq.0) then
         sxyrot = 0.d0
         do i = 1,N
            i2 = 2*(i-1)
            sxyrot = sxyrot + (V(i2+2)-Vprev(i2+2))/dt*
     :                     (dble(int(R(i2+1)))-dble(int(R(i2+1)+rx)))
         enddo
         do i = 1,Ndepl
            i2 = 2*(i-1)
            sxyrot = sxyrot + Mdepl*(Vdepl(i2+2)-Vprevdepl(i2+2))/dt*
     :              (dble(int(Rdepl(i2+1)))-dble(int(Rdepl(i2+1)+rx)))
         enddo
         deallocate(Vprev,Vprevdepl)
      endif

      return
      end
      
      subroutine vel_pos_init(N,Lx,Ly)
	implicit none 
      !Writing a function to initialize the positions and velocities of all the particles 
 
      integer i,j,N,Lx,Ly
      double precision randy,gauss,rx,ry
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      
      do i = 1,N
	rx = randy()
	R(2*i-1) = rx*Lx
	ry = randy()
	R(2*i) = ry*Ly
	V(2*i-1) = gauss(2)
	V(2*i) = gauss(2)
      end do  
      return
      end


      subroutine anderson(N)

C Hydrodynamics in the solvent is turned off if
C velocities are redrawn from a Maxwellian distribution (bath)
C every |ihydro| collision intervals. Note: this also acts as
C a thermostat, so separate thermostating will not be needed.

      integer N,i,j
      double precision vtemp(2),gauss,dummy
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      do j = 1,2
         vtemp(j) = 0.d0
      enddo
      do i = 1,N
         do j = 1,2
            V(2*(i-1)+j) = gauss(dummy)
            vtemp(j) = vtemp(j) + V(2*(i-1)+j)
         enddo
      enddo
      do j = 1,2
         vtemp(j) = vtemp(j)/dble(N)
      enddo
      do i = 1,N
         do j = 1,3
            V(2*(i-1)+j) = V(2*(i-1)+j) - vtemp(j)
         enddo
      enddo

      return
      end
