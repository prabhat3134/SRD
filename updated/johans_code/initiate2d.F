      program initiate

C Generates initial file for solvent flow around colloid
C simulation program laning.exe.
C
C Author:
C Johan T. Padding
C Date:
C 27 May 2004

      integer N,Ndepl,Ncol,i,icol,icol2
      integer Lx,Ly,step,ihydro
      double precision  Mdepl,Mcol,sigcc,sigcs,kappa,Rg,lambda
      double precision  V0,V0y,rcutccy,V0p,vflow
      logical fixcol,wallx,wally

      common / block1 / R,V,Rdepl,Vdepl,Rcol,Vcol
      double precision,pointer,dimension(:) :: R,V,Rdepl,Vdepl,Rcol,Vcol

      call readcontrol(Ncol,Mcol,sigcc,sigcs,kappa,V0,V0y,rcutccy,
     :   Ndepl,Mdepl,Rg,V0p,N,Lx,Ly,wallx,wally,vflow,fixcol,ihydro)

      write(*,*) 'N=',N
      write(*,*) 'Ndepl=',Ndepl
      write(*,*) 'Ncol=',Ncol

      if (N.ne.0) then
         allocate(R(2*N),V(2*N))
      else
         allocate(R(2),V(2))
         do i = 1,2
            R(i) = 0.d0
            V(i) = 0.d0
         enddo
      endif
      if (Ndepl.ne.0) then
         allocate(Rdepl(2*Ndepl),Vdepl(2*Ndepl))
      else
         allocate(Rdepl(2),Vdepl(2))
         do i = 1,2
            Rdepl(i) = 0.d0
            Vdepl(i) = 0.d0
         enddo
      endif
      if (Ncol.ne.0) then
         allocate(Rcol(2*Ncol),Vcol(2*Ncol))
      else
         allocate(Rcol(2),Vcol(2))
         do i = 1,2
            Rcol(i) = 0.d0
            Vcol(i) = 0.d0
         enddo
      endif

      if (Ncol.le.3) then
       do icol = 1,Ncol
         icol2 = 2*(icol-1)
         write(*,*)
         write(*,*) 'Colloid number ',icol
         write(*,*) 'Enter x-coordinate'
         read(*,*) Rcol(icol2+1)
         write(*,*) 'Enter y-coordinate'
         read(*,*) Rcol(icol2+2)
       enddo
      endif

      call loadranstate
      write(*,*) 'Generating positions'
      call generatepos(Ncol,sigcc,sigcs,kappa,V0,V0y,rcutccy,
     :                 Ndepl,Rg,V0p,N,Lx,Ly,wallx,wally)
      if (ihydro.ne.-1) then
         write(*,*) 'Generating velocities'
         call generatevel(N,Ncol,Mcol,Ndepl,Mdepl,vflow,fixcol)
      else
         do i = 1,2*N
            V(i) = 0.d0
         enddo
         do i = 1,2*Ndepl
            Vdepl(i) = 0.d0
         enddo
         do i = 1,2*Ncol
            Vcol(i) = 0.d0
         enddo
      endif
      step = 0
      write(*,*) 'Saving'
      call saveconf(step,N,Ndepl,Ncol)
      call saveranstate

      end


      subroutine readcontrol(Ncol,Mcol,sigcc,sigcs,kappa,V0,V0y,rcutccy,
     :   Ndepl,Mdepl,Rg,V0p,N,Lx,Ly,wallx,wally,vflow,fixcol,ihydro)

      integer N,Ndepl,Ncol
      integer Lx,Ly,iwallx,iwally,ihydro,ifixcol
      double precision Mdepl,Mcol,sigcc,sigcs,Rg
      double precision V0,V0y,V0p,vflow,kappa,rcutccy
      logical fixcol,wallx,wally

      open(unit=10,file='control.dat',form='formatted')
      read(10,*) 
      read(10,*) Ncol
      read(10,*) !frac
      read(10,*) Mcol
      read(10,*) sigcc
      read(10,*) V0
      read(10,*) V0y
      read(10,*) kappa
      read(10,*) rcutccy
      read(10,*) !friccol
      read(10,*)
      read(10,*) Ndepl
      read(10,*) Mdepl
      read(10,*) Rg
      read(10,*) V0p
      read(10,*) !fricdepl
      read(10,*)
      read(10,*) N
      read(10,*) sigcs
      read(10,*) !icol
      read(10,*)
      read(10,*) !dt
      read(10,*) !totsteps
      read(10,*) !maxcputime
      read(10,*) !isave
      read(10,*) !icont
      read(10,*) !rlist
      read(10,*) Lx
      read(10,*) Ly
      read(10,*) iwallx
      read(10,*) iwally
      read(10,*) !iflow
      read(10,*) vflow
      read(10,*) !itherm
      read(10,*) ihydro
      read(10,*) !grav
      read(10,*) !omega
      read(10,*) !ifixsolv
      read(10,*) !ifixdepl
      read(10,*) ifixcol
      read(10,*)
      read(10,*) !fsample
      read(10,*) !vsampletrue
      read(10,*) !vsample
      read(10,*) !ssample
      read(10,*) !sqtsample
      read(10,*) !idens
      read(10,*) !dbin
      read(10,*) !grsample
      read(10,*) !delgr
      read(10,*) !iorder
      read(10,*) !ipdb
      close(unit=10)

      fixcol  = .false.
      if (ifixcol .eq.1) fixcol  = .true.
      wallx = .false.
      if (iwallx.eq.1) wallx = .true.
      wally = .false.
      if (iwally.eq.1) wally = .true.

      return
      end


      subroutine generatepos(Ncol,sigcc,sigcs,kappa,V0,V0y,rcutccy,
     :                       Ndepl,Rg,V0p,N,Lx,Ly,wallx,wally)

      integer N,Ndepl,Ncol
      integer Lx,Ly,i,icol,icol2,i2,idepl
      double precision sigcc,V0,V0y,Rg,V0p,Rgsq,kappa,rcutccy
      double precision xi,yi,xij,yij,rij,rijsq,rcut,rr,randy
      double precision rLx,rLy,rLx2,rLy2,r6,r12,r24,r48
      double precision Vpot,prob,sigcs,sigcp,rcutcc,rcutcp,rcutcs,rcutpp
      logical wallx,wally
      common / block1 / R,V,Rdepl,Vdepl,Rcol,Vcol
      double precision,pointer,dimension(:) :: R,V,Rdepl,Vdepl,Rcol,Vcol
      
      Rgsq = Rg**2
      sigcp = 0.5d0*sigcc + Rg
      rcutcc = 2.d0**(1.d0/24.d0)*sigcc
      rcutcp = 2.d0**(1.d0/6.d0)*sigcp
      rcutcs = 2.d0**(1.d0/6.d0)*sigcs
      rcutpp = 3.5d0*Rg

      rLx = dble(Lx)
      rLy = dble(Ly)
      rLx2 = rLx/2.d0
      rLy2 = rLy/2.d0

C Do a Boltzmann fill of the box: first put in all polymers, given
C the positions of the colloids, then put in the solvent, given the
C positions of the colloids AND polymers (not entirely correct, but
C good enough).
C ### Colloids:
      if (Ncol.le.3) goto 8 !explicitly given
      i = 1
5     xi = randy()*dble(Lx)
      yi = randy()*dble(Ly)
      if (wallx) then
         if (xi.lt.0.5d0*sigcc) goto 5
         if (xi.gt.rLx-0.5d0*sigcc) goto 5
      endif
      if (wally) then
         if (yi.lt.0.5d0*sigcc) goto 5
         if (yi.gt.rLy-0.5d0*sigcc) goto 5
      endif
      Vpot = 0.d0
      do icol = 1,i-1
         icol2 = 2*(icol-1)
         xij = xi - Rcol(icol2+1)
         yij = yi - Rcol(icol2+2)
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
         rijsq = (xij**2 + yij**2)
         if (rijsq.lt.rcutcc**2) then
            r24 = sigcc**24/rijsq**12
            r48 = r24**2
            Vpot = Vpot+4.d0*V0*(r48-r24+0.25d0)
         endif
         if (V0y.ne.0.d0 .and. rijsq.lt.rcutccy**2) then
            rij = dsqrt(rijsq)
            Vpot = Vpot+V0y*dexp(-kappa*(rij-sigcc))*sigcc/rij
         endif
      enddo !icol
      prob = dexp(-Vpot)
      if (prob.lt.0.0001d0) goto 5
      rr = randy()
      if (rr.lt.prob) then
         i2 = 2*(i-1)
         Rcol(i2+1) = xi 
         Rcol(i2+2) = yi 
         i = i+1
      endif
      write(*,*) 'colloid ',i
      if (i.le.Ncol) goto 5
C ### Depletant:
8     i = 1
      if (Ndepl.eq.0) goto 15
10    xi = randy()*dble(Lx)
      yi = randy()*dble(Ly)
      if (wallx) then
         if (xi.lt.Rg) goto 10
         if (xi.gt.rLx-Rg) goto 10
      endif
      if (wally) then
         if (yi.lt.Rg) goto 10
         if (yi.gt.rLy-Rg) goto 10
      endif
      if (Ndepl.eq.1) then
         xi = 0.d0
         yi = 0.d0
      endif
      Vpot = 0.d0
      do icol = 1,Ncol
         icol2 = 2*(icol-1)
         xij = xi - Rcol(icol2+1)
         yij = yi - Rcol(icol2+2)
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
         rijsq = (xij**2 + yij**2)
         if (rijsq.lt.rcutcp**2) then
            r6 = sigcp**6/rijsq**3
            r12 = r6**2
            Vpot = Vpot+(r12-r6+0.25d0)
         endif
      enddo !icol
      Vpot = 4.d0*V0*Vpot
      if (V0p.eq.0.d0) goto 11
      do idepl = 1, i-1
         i2 = 2*(idepl-1)
         xij = xi - Rdepl(i2+1)
         yij = yi - Rdepl(i2+2)
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
         rijsq = (xij**2 + yij**2)
         if (rijsq.lt.rcutpp**2) then
            Vpot = Vpot + V0p*dexp(-0.75d0*rijsq/Rgsq)
         endif
      enddo
11    prob = dexp(-Vpot) !probability of accepting this config
      if (prob.lt.1.d0) goto 10
      rr = randy()
      if (rr.lt.prob) then
         i2 = 2*(i-1)
         Rdepl(i2+1) = xi 
         Rdepl(i2+2) = yi 
         i = i+1
      endif
      if (mod(i,10).eq.0) write(*,*) 'depletant ',i
      if (i.le.Ndepl) goto 10
C ### Solvent:
      i = 1
15    if (N.eq.0) goto 30
20    xi = randy()*dble(Lx)
      yi = randy()*dble(Ly)
      Vpot = 0.d0
      do icol = 1,Ncol
         icol2 = 2*(icol-1)
         xij = xi - Rcol(icol2+1)
         yij = yi - Rcol(icol2+2)
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
         rijsq = (xij**2 + yij**2)
         if (rijsq.lt.rcutcs**2) then
            r6 = sigcs**6/rijsq**3
            r12 = r6**2
            Vpot = Vpot+(r12-r6+0.25d0)
         endif
      enddo !icol
      Vpot = 4.d0*V0*Vpot
      prob = dexp(-Vpot) !probability of accepting this config
      if (prob.lt.1.d0) goto 20
      rr = randy()
      if (rr.lt.prob) then
         i2 = 2*(i-1)
         R(i2+1) = xi 
         R(i2+2) = yi 
         i = i+1
      endif
      if (mod(i,1000).eq.0) write(*,*) 'solvent ',i
      if (i.le.N) goto 20
30    continue

      return
      end


      subroutine generatevel(N,Ncol,Mcol,Ndepl,Mdepl,vflow,
     :                       fixcol)

      integer N,Ndepl,Ncol,i,i2,j
      double precision  Mdepl,Mcol,Minv,vflow,gauss,dummy
      double precision  totv(2),totvdepl(2),totvcol(2),vv
      logical fixcol
      common / block1 / R,V,Rdepl,Vdepl,Rcol,Vcol
      double precision,pointer,dimension(:) :: R,V,Rdepl,Vdepl,Rcol,Vcol

C Assign Maxwellian distributed velocities and calculate total momentum
C for solvent, depletant, and colloids separately

      do j = 1,2
         totv(j) = 0.d0
         totvdepl(j) = 0.d0
         totvcol(j) = 0.d0
      enddo
      do i = 1,N
         i2 = 2*(i-1)
         do j = 1,2
            vv = gauss(dummy)
            V(i2+j) = vv
            totv(j) = totv(j)+vv
            if (j.eq.2) V(i2+j) = V(i2+j) + vflow !add average flowfield
         enddo
      enddo
      Minv = 1.d0/dsqrt(Mdepl)
      do i = 1, Ndepl
         i2 = 2*(i-1)
         do j = 1,2
            vv = Minv*gauss(dummy)
            Vdepl(i2+j) = vv
            totvdepl(j) = totvdepl(j)+vv
            if (j.eq.2) Vdepl(i2+j) = Vdepl(i2+j) + vflow !add av. flowf
         enddo
      enddo
      if (fixcol) then
         do i = 1, 2*Ncol
            Vcol(i) = 0.d0
         enddo
      else
         Minv = 1.d0/dsqrt(Mcol)
         do i = 1, Ncol
            i2 = 2*(i-1)
            do j = 1,2
               vv = Minv*gauss(dummy)
               Vcol(i2+j) = vv
               totvcol(j) = totvcol(j)+vv
            enddo
         enddo
      endif
C Subtract total momentum (relative to flowfield)
      do j = 1,2
         totv(j) = totv(j)/dble(N)
         totvdepl(j) = totvdepl(j)/dble(Ndepl)
         totvcol(j) = totvcol(j)/dble(Ncol)
      enddo
      do i = 1,N
         i2 = 2*(i-1)
         do j = 1,2
            V(i2+j) = V(i2+j) - totv(j)
         enddo
      enddo
      do i = 1,Ndepl
         i2 = 2*(i-1)
         do j = 1,2
            Vdepl(i2+j) = Vdepl(i2+j) - totvdepl(j)
         enddo
      enddo
      do i = 1,Ncol
         i2 = 2*(i-1)
         do j = 1,2
            Vcol(i2+j) = Vcol(i2+j) - totvcol(j)
         enddo
      enddo

      return
      end


      subroutine saveconf(step,N,Ndepl,Ncol)

      integer step
      integer N,Ndepl,Ncol
      common / block1 / R,V,Rdepl,Vdepl,Rcol,Vcol
      double precision,pointer,dimension(:) :: R,V,Rdepl,Vdepl,Rcol,Vcol

      open(unit=11,file='config.dat',form='unformatted')
      write(11) step
      write(11) N,Ndepl,Ncol
      write(11) R,V,Rdepl,Vdepl,Rcol,Vcol
      close(unit=11)

      return
      end

