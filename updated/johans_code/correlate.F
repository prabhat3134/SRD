      subroutine loadcorr(drsamp,dcount,sac,scount,sframe,steller,
     :            fframe,fteller,fact,factcount,
     :            vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,
     :            vactcount,vactcounttrue,vframe,vframetrue,vteller,
     :            vtellertrue,grbin,totgr,sumvavdepl,
     :            sqtx,sqty,sqtcount,sqtframe,sqtteller)

      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue
      common / block1b / rstore
      double precision,pointer,dimension(:,:) :: rstore
      common / block1c / fstore
      double precision,pointer,dimension(:,:) :: fstore
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) ::
     :                  Fdrx,Fdry,Fdrxsq,Fdrysq
      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision drsamp,dcount,sac(1000),scount(1000)
      double precision sqtx(10,100),sqty(10,100),sqtcount(100)
      double precision fact(2,10000),factcount(10000)
      double precision vact(2,1000),vacttrue(2,1000),vactcount(1000)
      double precision vactcounttrue(1000)
      double precision vsed(2),vsedtrue(2),vacr(4,100),vacrcount(2,100)
      double precision grbin(100),totgr,sumvavdepl
      integer sframe,steller,vframe,vframetrue,vteller,vtellertrue
      integer fframe,fteller,sqtframe,sqtteller
      logical ex
         
      inquire(file='corr.dat',exist=ex)
      if (ex) then
         open(unit=14,file='corr.dat',form='unformatted')
         read(14) Fdrx,Fdry,Fdrxsq,Fdrysq,drsamp,dcount,
     :            dens,denssolv,flowfieldr,flowfieldy,flowhit,
     :            sac,scount,sframe,steller,Rcolprev,
     :            fframe,fteller,fact,factcount,fstore,
     :            vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,
     :            vstore,vstoretrue,vactcount,vactcounttrue,
     :            vframe,vframetrue,vteller,vtellertrue,
     :            grbin,totgr,sumvavdepl,rstore,sqtx,sqty,
     :            sqtcount,sqtframe,sqtteller
         close(unit=14)
      else
         stop 'continuation of run not possible: no corr.dat'
      endif

      return
      end
      
      
      subroutine savecorr(drsamp,dcount,sac,scount,sframe,steller,
     :            fframe,fteller,fact,factcount,
     :            vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,
     :            vactcount,vactcounttrue,vframe,vframetrue,
     :            vteller,vtellertrue,grbin,totgr,sumvavdepl,
     :            sqtx,sqty,sqtcount,sqtframe,sqtteller)

 
      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue
      common / block1b / rstore
      double precision,pointer,dimension(:,:) :: rstore
      common / block1c / fstore
      double precision,pointer,dimension(:,:) :: fstore
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) ::
     :                  Fdrx,Fdry,Fdrxsq,Fdrysq
      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision drsamp,dcount,sac(1000),scount(1000)
      double precision sqtx(10,100),sqty(10,100),sqtcount(100)
      double precision fact(2,10000),factcount(10000)
      double precision vact(2,1000),vacttrue(2,1000),vactcount(1000)
      double precision vactcounttrue(1000)
      double precision vsed(2),vsedtrue(2),vacr(4,100),vacrcount(2,100)
      double precision grbin(100),totgr,sumvavdepl
      integer sframe,steller,vframe,vframetrue,vteller,vtellertrue
      integer fframe,fteller,sqtframe,sqtteller

      open(unit=14,file='corr.dat',form='unformatted',status='unknown')
      write(14) Fdrx,Fdry,Fdrxsq,Fdrysq,drsamp,dcount,
     :            dens,denssolv,flowfieldr,flowfieldy,flowhit,
     :            sac,scount,sframe,steller,Rcolprev,
     :            fframe,fteller,fact,factcount,fstore,
     :            vsed,vsedtrue,vacr,vacrcount,vact,vacttrue,
     :            vstore,vstoretrue,vactcount,vactcounttrue,
     :            vframe,vframetrue,vteller,vtellertrue,
     :            grbin,totgr,sumvavdepl,rstore,sqtx,sqty,
     :            sqtcount,sqtframe,sqtteller
      close(unit=14)

      return
      end


      subroutine setupsac(sac,scount,sframe,steller)
 
      double precision sac(1000),scount(1000)
      integer sframe,steller,i

      do i = 1,1000
         sac(i) = 0.d0
         scount(i) = 0.d0
      enddo
      sframe = 0
      steller = 0

      return
      end


      subroutine updatesac
     :   (sac,scount,ssample,Lx,Ly,N,Ndepl,dt,Mdepl,nu)
 
      double precision sac(1000),scount(1000)
      integer ssample,Lx,Ly,N,Ndepl,icorr
      double precision nu,dum,dt,Mdepl

      open (unit=12,file='modulus.dat',status='unknown')
      do icorr = 1,1000
         if (scount(icorr).ne.0.d0) then
            dum = sac(icorr)/(dble(Lx*Ly)*scount(icorr))
            write(12,*) dble((icorr-1)*ssample)*dt,dum
            if (icorr.eq.1) then
               nu = 0.5d0*dum*dble(ssample)*dt
            else
               nu = nu+dum*dble(ssample)*dt
            endif
         endif
      enddo
      close(unit=12)
C      nu = nu*dble(Lx*Ly*Lz)/(dble(N)+Mdepl*dble(Ndepl))
      nu = nu*dble(Lx*Ly)/dble(N)
 
      return
      end
 
 
      subroutine stressac(store,sac,scount,sframe,steller,
     :                    N,Ndepl,Mdepl)

      double precision store(1000)
      double precision sac(1000),scount(1000)
      integer sframe,steller,N,Ndepl,lastfra,i,i2,icorr,iframe
      double precision Mdepl,sxy,sxz,syz,si
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block3 / sxyfr,sxyrot,sxykin
      double precision sxyfr,sxyrot,sxykin
 
      steller = steller + 1
      lastfra = sframe
      sframe  = sframe + 1
      if (sframe.gt.1000) sframe = 1
      sxy = sxykin + sxyrot
      store(sframe) = sxy+sxyfr
         si = store(sframe)
         iframe = sframe+1
         do icorr = 1,1000
            iframe = iframe-1
            if (iframe.lt.1) then
               if (steller.le.1000) goto 25
               iframe = 1000
            endif
            sac(icorr) = sac(icorr) + si*store(iframe)
            scount(icorr) = scount(icorr) + 1.d0
         enddo
25       continue
 
      return
      end


      subroutine setupsqt(sqtx,sqty,sqtcount,sqtframe,sqtteller)
 
      double precision sqtx(10,100),sqty(10,100),sqtcount(100)
      integer sqtframe,sqtteller,i,j

      do i = 1,100
         do j = 1,10
           sqtx(j,i) = 0.d0
           sqty(j,i) = 0.d0
         enddo
         sqtcount(i) = 0.d0
      enddo
      sqtframe = 0
      sqtteller = 0

      return
      end


      subroutine updatesqt(sqtx,sqty,sqtcount,
     :                     sqtsample,dt,Lx,Ly)
 
      double precision sqtx(10,100),sqty(10,100),sqtcount(100)
      integer sqtsample,icorr
      integer Lx,Ly,iq
      double precision dt,pi,riq,wavex,wavey,qx,qy

      pi = 4.d0*datan(1.d0)
      wavex = 2.d0*pi/dble(Lx)
      wavey = 2.d0*pi/dble(Ly)
      open (unit=11,file='wavevectors.dat',status='unknown')
      do iq = 1,10
         riq = dble(iq)
         qx = riq*wavex
         qy = riq*wavey
         write(11,*) qx,qy
      enddo
      close (unit=11)

      open (unit=12,file='sqtx.dat',status='unknown')
      open (unit=13,file='sqty.dat',status='unknown')
      do icorr = 1,100
         if (sqtcount(icorr).ne.0.d0) then
            write(12,'(11(E10.4,1X))') dble((icorr-1)*sqtsample)*dt,
     :         (sqtx(iq,icorr)/sqtcount(icorr),iq=1,10)
            write(13,'(11(E10.4,1X))') dble((icorr-1)*sqtsample)*dt,
     :         (sqty(iq,icorr)/sqtcount(icorr),iq=1,10)
         endif
      enddo
      close(unit=12)
      close(unit=13)
 
      return
      end
 

      subroutine calcsqt(sqtx,sqty,sqtcount,sqtframe,sqtteller,
     :                   Ncol,Lx,Ly)

      double precision sqtx(10,100),sqty(10,100),sqtcount(100)
      integer sqtframe,sqtteller,Ncol,lastfra,i,i2,j,j2,icorr,iframe
      integer Lx,Ly,iq
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1b / rstore
      double precision,pointer,dimension(:,:) :: rstore
      double precision pi,riq,wavex,wavey,xi,yi,dx,dy
      double precision qx(10),qy(10)

      pi = 4.d0*datan(1.d0)
      wavex = 2.d0*pi/dble(Lx)
      wavey = 2.d0*pi/dble(Ly)
      do iq = 1,10
         riq = dble(iq)
         qx(iq) = riq*wavex
         qy(iq) = riq*wavey
      enddo

      sqtteller = sqtteller + 1
      lastfra = sqtframe
      sqtframe  = sqtframe + 1
      if (sqtframe.gt.100) sqtframe = 1
      do i = 1,2*Ncol
         rstore(i,sqtframe) = Rcol(i)
      enddo
      iframe = sqtframe+1
      do icorr = 1,100
         iframe = iframe-1
         if (iframe.lt.1) then
            if (sqtteller.le.100) goto 25
            iframe = 100
         endif
         do i = 1,Ncol
            i2 = 2*(i-1)
            xi = rstore(i2+1,sqtframe)
            yi = rstore(i2+2,sqtframe)
            do j = 1,Ncol
               j2 = 2*(j-1)
               dx = xi - rstore(j2+1,iframe)
               dy = yi - rstore(j2+2,iframe)
               do iq = 1,10
                sqtx(iq,icorr) = sqtx(iq,icorr) + dcos(qx(iq)*dx)
                sqty(iq,icorr) = sqty(iq,icorr) + dcos(qy(iq)*dy)
               enddo
            enddo 
         enddo
         sqtcount(icorr) = sqtcount(icorr) + dble(Ncol)
      enddo
25    continue
 
      return
      end


      subroutine setupfac(fact,factcount,fframe,fteller)

      double precision fact(2,10000),factcount(10000)
      integer fframe,fteller,Ncol,i,j

      do i = 1,10000
         factcount(i) = 0.d0
         do j = 1,2
            fact(j,i) = 0.d0
         enddo
      enddo
      fframe = 0
      fteller = 0
          
      return
      end


      subroutine updatefact(fact,factcount,fsample,dt)

      double precision fact(2,10000),factcount(10000),dt,dum1,dum2
      integer fsample,icorr
     
      open (unit=12,file='fact.dat',status='unknown')
      do icorr = 1,10000
         if (factcount(icorr).ne.0.d0) then
            dum1 = fact(1,icorr)/factcount(icorr)
            dum2 = fact(2,icorr)/factcount(icorr)
            write(12,*) dble((icorr-1)*fsample)*dt,dum1,dum2
         endif
      enddo
      close(unit=12)

      return
      end


      subroutine calcfact(fact,factcount,fframe,fteller,Ncol)

      double precision fact(2,10000),factcount(10000)
      integer fframe,fteller,Ncol,lastfra,i2,iframe,icorr,i
       common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1c / fstore
      double precision,pointer,dimension(:,:) :: fstore
      double precision fx,fy
 
      fteller = fteller+1
      lastfra = fframe
      fframe = fframe + 1
      if (fframe.gt.10000) fframe = 1

      do i = 1,Ncol
         i2 = 2*(i-1)
         fx = Fcol(i2+1)
         fy = Fcol(i2+2)
         fstore(fframe,i2+1) = fx
         fstore(fframe,i2+2) = fy
         iframe = fframe+1
         do icorr = 1,10000
            iframe = iframe-1
            if (iframe.lt.1) then
               if (fteller.le.10000) goto 25
               iframe = 10000
            endif
            fact(1,icorr) = fact(1,icorr) + fx*fstore(iframe,i2+1)  
            fact(2,icorr) = fact(2,icorr) + fy*fstore(iframe,i2+2)
            factcount(icorr) = factcount(icorr) + 1.d0
         enddo
 25      continue
      enddo

      return
      end


      subroutine setupvac(vsed,vacr,vacrcount,vact,vactcount,
     :                    vframe,vteller,Ncol)
 
      double precision vact(2,1000),vactcount(1000)
      double precision vsed(2),vacr(4,100),vacrcount(2,100)
      integer vframe,vteller,Ncol,i,j
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue

      do i = 1,1000
         vactcount(i) = 0.d0
         do j = 1,2
            vact(j,i) = 0.d0
         enddo
      enddo
      do i = 1,100
         do j = 1,2
            vacrcount(j,i) = 0.d0
         enddo
         do j = 1,4
            vacr(j,i) = 0.d0
         enddo
      enddo
      vsed(1) = 0.d0
      vsed(2) = 0.d0
      vframe = 0
      vteller = 0
      do i = 1,2*Ncol
         Rcolprev(i) = Rcol(i)
      enddo

      return
      end


      subroutine setupvactrue(vsedtrue,vacttrue,vactcounttrue,
     :                    vframetrue,vtellertrue,Ncol)
 
      double precision vacttrue(2,1000),vactcounttrue(1000),vsedtrue(2)
      integer vframetrue,vtellertrue,Ncol,i,j

      do i = 1,1000
         vactcounttrue(i) = 0.d0
         do j = 1,2
            vacttrue(j,i) = 0.d0
         enddo
      enddo
      vsedtrue(1) = 0.d0
      vsedtrue(2) = 0.d0
      vframetrue = 0
      vtellertrue = 0

      return
      end


      subroutine updatevac(vsed,vacr,vacrcount,sigcs,vact,vactcount,
     :                     vteller,vsample,dt,Ncol)
 
      double precision vact(2,1000),vactcount(1000)
      double precision vsed(2),vacr(4,100),vacrcount(2,100),sigcs
      integer vteller,vsample,Ncol,icorr
      double precision dum1,dum2,dum3,dum4,dt,diam,avx,avy

      diam = 2.d0*sigcs

      open (unit=12,file='vsed.dat',status='unknown')
      if (vteller.ne.0) then
         avx = vsed(1)/dble(Ncol*vteller)
         avy = vsed(2)/dble(Ncol*vteller)
         write(12,*) 'Average x-velocity = ',avx
         write(12,*) 'Average y-velocity = ',avy
      else
         avx = 0.d0
         avy = 0.d0
         write(12,*) 'Warning: vteller is still zero'
      endif
      close(unit=12)

      open (unit=12,file='vact.dat',status='unknown')
      do icorr = 1,1000
         if (vactcount(icorr).ne.0.d0) then
            dum1 = vact(1,icorr)/vactcount(icorr) - avx**2
            dum2 = vact(2,icorr)/vactcount(icorr) - avy**2
            write(12,*) dble((icorr-1)*vsample)*dt,dum1,dum2
         endif
      enddo
      close(unit=12)

      open (unit=12,file='vacr.dat',status='unknown')
      do icorr = 1,100
         if (vacrcount(1,icorr).ne.0.d0) then
            dum1 = vacr(1,icorr)/vacrcount(1,icorr) - avx**2
            dum3 = vacr(3,icorr)/vacrcount(1,icorr) - avy**2
         else
            dum1 = 0.d0
            dum3 = 0.d0
         endif
         if (vacrcount(2,icorr).ne.0.d0) then
            dum2 = vacr(2,icorr)/vacrcount(2,icorr) - avx**2
            dum4 = vacr(4,icorr)/vacrcount(2,icorr) - avy**2
         else
            dum2 = 0.d0
            dum4 = 0.d0
         endif
         write(12,'(5D12.4,1X)') dble(icorr-1)*diam,dum1,dum2,dum3,dum4
      enddo
      close(unit=12)
 
      return
      end
 
 
      subroutine updatevactrue(vsedtrue,vacttrue,vactcounttrue,
     :                     vtellertrue,vsampletrue,dt,Ncol)
 
      double precision vacttrue(2,1000),vactcounttrue(1000)
      double precision vsedtrue(2)
      integer vtellertrue,vsampletrue,Ncol,icorr
      double precision dum1,dum2,dt,avxtrue,avytrue

      open (unit=12,file='vsedtrue.dat',status='unknown')
      if (vtellertrue.ne.0) then
         avxtrue = vsedtrue(1)/dble(Ncol*vtellertrue)
         avytrue = vsedtrue(2)/dble(Ncol*vtellertrue)
         write(12,*) 'Average x-velocity = ',avxtrue
         write(12,*) 'Average y-velocity = ',avytrue
      else
         avxtrue = 0.d0
         avytrue = 0.d0
         write(12,*) 'Warning: vteller is still zero'
      endif
      close(unit=12)

      open (unit=12,file='vacttrue.dat',status='unknown')
      do icorr = 1,1000
         if (vactcounttrue(icorr).ne.0.d0) then
            dum1 = vacttrue(1,icorr)/vactcounttrue(icorr) - avxtrue**2
            dum2 = vacttrue(2,icorr)/vactcounttrue(icorr) - avytrue**2
            write(12,*) dble((icorr-1)*vsampletrue)*dt,dum1,dum2
         endif
      enddo
      close(unit=12)

      return
      end


      subroutine velac(vsed,vacr,vacrcount,sigcs,vact,vactcount,
     :                 vframe,vteller,Ncol,Lx,Ly,wallx,wally,
     :                 ihydro,dt,vsample)

      double precision vact(2,1000),vactcount(1000)
      double precision vsed(2),vacr(4,100),vacrcount(2,100),sigcs,invdt
      integer vframe,vteller,Ncol,lastfra,i,i2,j,j2,icorr,iframe
      integer Lx,Ly,xbin,ybin,ihydro,vsample
      logical wallx,wally
      double precision dt,vx,vy,vxj,vyj,xi,yi,xij,yij,rLx,rLy,rLx2,rLy2,
     :                 invdiam
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue
 
      rLx = dble(Lx)
      rLy = dble(Ly)
      rLx2 = 0.5d0*rLx
      rLy2 = 0.5d0*rLy
      invdiam = 0.5d0/sigcs
      invdt = 1.d0/(dt*dble(vsample))

      vteller = vteller + 1
      lastfra = vframe
      vframe  = vframe + 1
      if (vframe.gt.1000) vframe = 1

      do i = 1,Ncol
         xij = Rcol(2*i-1) - Rcolprev(2*i-1)
         if (.not.wallx) then
            if (xij.gt.rLx2) xij = xij-rLx
            if (xij.lt.-rLx2) xij = xij+rLx
         endif
         vstore(vframe,2*i-1) = xij*invdt
         yij = Rcol(2*i) - Rcolprev(2*i)
         if (.not.wally) then
            if (yij.gt.rLy2) yij = yij-rLy
            if (yij.lt.-rLy2) yij = yij+rLy
         endif
         vstore(vframe,2*i) = yij*invdt
      enddo
      do i = 1,2*Ncol
         Rcolprev(i) = Rcol(i)
      enddo

      do i = 1,Ncol
         i2 = 2*(i-1)
         vx = vstore(vframe,i2+1)
         vy = vstore(vframe,i2+2)
         vsed(1) = vsed(1) + vx
         vsed(2) = vsed(2) + vy
         iframe = vframe+1
         do icorr = 1,1000
            iframe = iframe-1
            if (iframe.lt.1) then
               if (vteller.le.1000) goto 25
               iframe = 1000
            endif
            vact(1,icorr) = vact(1,icorr) + vx*vstore(iframe,i2+1)
            vact(2,icorr) = vact(2,icorr) + vy*vstore(iframe,i2+2)
            vactcount(icorr) = vactcount(icorr) + 1.d0
         enddo
25       continue
      enddo

      do i = 1,Ncol
         i2 = 2*(i-1)
         vx = vstore(vframe,i2+1)
         vy = vstore(vframe,i2+2)
         xi = Rcol(i2+1)
         yi = Rcol(i2+2)
         do j = i,Ncol
            j2 = 2*(j-1)
            xij = xi-Rcol(j2+1)
            yij = yi-Rcol(j2+2)
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
            xij = dabs(xij)
            yij = dabs(yij)
            xbin = NINT(xij*invdiam)+1
            ybin = NINT(yij*invdiam)+1
            if (xbin.eq.1 .or. ybin.eq.1) then
               vxj = vstore(vframe,j2+1)
               vyj = vstore(vframe,j2+2)
               if (ybin.eq.1 .and. xbin.le.100) then
                  vacr(1,xbin) = vacr(1,xbin) + vx*vxj
                  vacr(3,xbin) = vacr(3,xbin) + vy*vyj
                  vacrcount(1,xbin) = vacrcount(1,xbin) + 1.d0
               endif
               if (xbin.eq.1 .and. ybin.le.100) then
                  vacr(2,ybin) = vacr(2,ybin) + vx*vxj
                  vacr(4,ybin) = vacr(4,ybin) + vy*vyj
                  vacrcount(2,ybin) = vacrcount(2,ybin) + 1.d0
               endif
            endif
         enddo 
      enddo
 
      return
      end


      subroutine velactrue(vsedtrue,vacttrue,vactcounttrue,
     :                 vframetrue,vtellertrue,Ncol,ihydro)

      double precision vacttrue(2,1000),vactcounttrue(1000),vsedtrue(2)
      integer vframetrue,vtellertrue,Ncol,lastfra,i,i2,j2,icorr,iframe
      integer ihydro
      double precision vxtrue,vytrue
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block1a / Rcolprev,vstore,vstoretrue
      double precision,pointer,dimension(:) :: Rcolprev
      double precision,pointer,dimension(:,:) :: vstore,vstoretrue
 
      vtellertrue = vtellertrue + 1
      lastfra = vframetrue
      vframetrue  = vframetrue + 1
      if (vframetrue.gt.1000) vframetrue = 1

      if (ihydro.ne.-1) then
         do i = 1,2*Ncol
            vstoretrue(vframetrue,i) = Vcol(i)
         enddo
      else !no true velocity in BD
         do i = 1,2*Ncol
            vstoretrue(vframetrue,i) = 0.d0
         enddo
      endif

      do i = 1,Ncol
         i2 = 2*(i-1)
         vxtrue = vstoretrue(vframetrue,i2+1)
         vytrue = vstoretrue(vframetrue,i2+2)
         vsedtrue(1) = vsedtrue(1) + vxtrue
         vsedtrue(2) = vsedtrue(2) + vytrue
         iframe = vframetrue+1
         do icorr = 1,1000
            iframe = iframe-1
            if (iframe.lt.1) then
               if (vtellertrue.le.1000) goto 25
               iframe = 1000
            endif
            vacttrue(1,icorr) = vacttrue(1,icorr) 
     :                        + vxtrue*vstoretrue(iframe,i2+1)
            vacttrue(2,icorr) = vacttrue(2,icorr) 
     :                        + vytrue*vstoretrue(iframe,i2+2)
            vactcounttrue(icorr) = vactcounttrue(icorr) + 1.d0
         enddo
25       continue
      enddo

      return
      end


      subroutine setupdens(dcount,rbin,ybin)

      integer rbin,ybin,ir,iy
      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision dcount,dbin,avr,avy,norm,pi
      logical yxdens

      dcount = 0.d0
      do ir = 1, rbin+1
         do iy = 1, ybin+1
            dens(iy,ir) = 0.d0
            denssolv(iy,ir) = 0.d0
         enddo
      enddo

      return
      end


      subroutine updatedens(dcount,rbin,ybin,dbin,yxdens)
 
      integer rbin,ybin,irbin,iybin
      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision dcount,dbin,avr,avy,norm,pi
      logical yxdens

      pi = 4.d0*datan(1.d0)

      open(unit=10,file='dens.dat')
      open(unit=11,file='denssolv.dat')
      do iybin = 1,ybin+1
         avy = dble(iybin-1)*dbin
         do irbin = 1,rbin+1
           avr = dble(irbin-1)*dbin
           norm = 2.d0*dbin**2
           if (iybin.eq.1) norm = norm*0.5d0
           norm = norm*dcount
           write(10,'(3(E12.5,1X))') avy,avr,dens(iybin,irbin)/norm
           write(11,'(3(E12.5,1X))') avy,avr,denssolv(iybin,irbin)/norm
         enddo
         write(10,*)
         write(11,*)
      enddo
      close(10)
      close(11)

      open(unit=12,file='flowfield.dat')
      do iybin = 1,ybin+1
         avy  = dble(iybin-1)*dbin
         do irbin = 1,rbin+1
            avr = dble(irbin-1)*dbin
            norm = flowhit(iybin,irbin)
            if (norm.ne.0.d0) then
               write(12,'(4(E12.5,1X))') avy,avr,
     :         flowfieldy(iybin,irbin)/norm,
     :         flowfieldr(iybin,irbin)/norm
            endif
         enddo
      enddo
      close(12)

      return
      end


      subroutine densitymap(dcount,rbin,ybin,dbin,
     :                      N,Ndepl,Lx,Ly,yxdens)

      common / blockdens / dens,denssolv,flowfieldr,flowfieldy,flowhit
      double precision,pointer,dimension(:,:) :: dens,denssolv
      double precision,pointer,dimension(:,:) :: flowfieldr,flowfieldy,
     :                                           flowhit
      double precision dcount,dbin,rLx,rLy,rr,x,y,vx,vy
      integer Ndepl,N,rbin,ybin,Lx,Ly,i,i2,irbin,iybin
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      logical yxdens

      dcount = dcount+1.d0
      rLx = dble(Lx)
      rLy = dble(Ly)
      do i = 1, Ndepl
         i2 = 2*(i-1)
         x = Rdepl(i2+1)
         y = Rdepl(i2+2)
         rr = x
         irbin = nint(rr/dbin)+1
         iybin = nint(y/dbin)+1
         if (irbin.le.(rbin+1) .and. iybin.le.(ybin+1)) then
            dens(iybin,irbin) = dens(iybin,irbin)+1.d0
         endif
      enddo
      do i = 1, N
         i2 = 2*(i-1)
         x = R(i2+1)
         y = R(i2+2)
         vx = V(i2+1)
         vy = V(i2+2)
         rr = x
         irbin = nint(rr/dbin)+1
         iybin = nint(y/dbin)+1
         if (irbin.le.(rbin+1) .and. iybin.le.(ybin+1)) then
            denssolv(iybin,irbin) = denssolv(iybin,irbin)+1.d0
            flowfieldr(iybin,irbin) = flowfieldr(iybin,irbin)+vx
            flowfieldy(iybin,irbin) = flowfieldy(iybin,irbin)+vy
            flowhit(iybin,irbin) = flowhit(iybin,irbin)+1.d0
         endif
      enddo
         
      return
      end


      subroutine setupdrift(Ncol,drsamp)

      integer Ncol,i
      double precision drsamp
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) ::
     :                  Fdrx,Fdry,Fdrxsq,Fdrysq

      drsamp = 0.d0
      do i = 1,2*Ncol
         Fdrx(i) = 0.d0
         Fdry(i) = 0.d0
         Fdrxsq(i) = 0.d0
         Fdrysq(i) = 0.d0
      enddo

      return
      end


      subroutine updatedrift(Ncol,drsamp)

      integer Ncol
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block4 / Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision,pointer,dimension(:) ::
     :                  Fdrx,Fdry,Fdrxsq,Fdrysq
      double precision drsamp, del(2)
      double precision fx,fy
      integer i,ic,m

      open(unit=11,file='drift.dat')
      write(11,*) 'Depletant drift force'
      do ic = 1, Ncol
         del(1) = Fdrxsq(ic)/drsamp-(Fdrx(ic)/drsamp)**2
         del(2) = Fdrysq(ic)/drsamp-(Fdry(ic)/drsamp)**2
         do i = 1,2
            if (del(i).gt.0.d0) then
               del(i) = dsqrt(del(i))
            else
               del(i) = 0.d0
            endif
         enddo
         fx = Fdrx(ic)
         fy = Fdry(ic)
         write(11,'(6(E10.4,1X))') (Rcol(2*(ic-1)+m),m=1,2),
     :      fx/drsamp,fy/drsamp,
     :      (del(m),m=1,2)
      enddo
      write(11,*) 'Solvent drift force'
      do ic = Ncol+1, 2*Ncol
         del(1) = Fdrxsq(ic)/drsamp-(Fdrx(ic)/drsamp)**2
         del(2) = Fdrysq(ic)/drsamp-(Fdry(ic)/drsamp)**2
         do i = 1,2
            if (del(i).gt.0.d0) then
               del(i) = dsqrt(del(i))
            else
               del(i) = 0.d0
            endif
         enddo
         fx = Fdrx(ic)
         fy = Fdry(ic)
         write(11,'(6(E10.4,1X))') (Rcol(2*(ic-Ncol-1)+m),m=1,2),
     :      fx/drsamp,fy/drsamp,
     :      (del(m),m=1,2)
      enddo
      close(11)

      return
      end


      subroutine setupgr(grbin,totgr)

      double precision grbin(100),totgr
      integer i

      do i = 1,100
         grbin(i) = 0.d0
      enddo
      totgr = 0.d0

      return
      end


      subroutine radial(grbin,totgr,delgr,Ncol,Lx,Ly,wallx,wally)

      double precision grbin(100),totgr,delgr
      integer Ncol,Lx,Ly,ibin,i,i2,j,j2
      double precision xi,yi,xij,yij,rijsq,rij,rmaxsq
      double precision rLx,rLy,rLx2,rLy2
      logical wallx,wally
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      rmaxsq = (100.d0*delgr)**2
      rLx = dble(Lx)
      rLy = dble(Ly)
      rLx2 = 0.5d0*rLx
      rLy2 = 0.5d0*rLy
      totgr = totgr+1.d0

      do i = 2, Ncol
         i2 = 2*(i-1)
         xi = Rcol(i2+1)
         yi = Rcol(i2+2)
         do j = 1,i-1
            j2 = 2*(j-1)
            xij = xi-Rcol(j2+1)
            yij = yi-Rcol(j2+2)
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
            if (rijsq.lt.rmaxsq) then
               rij = dsqrt(rijsq)
               ibin = INT(rij/delgr)+1
               grbin(ibin) = grbin(ibin)+1.d0
            endif
         enddo
      enddo
      if (Ncol.eq.2) then
         open(unit=100,file='displacement.dat',position='append')
         write(100,'(2(D20.10,1X))') xij,yij 
         close(unit=100)
      endif

      return
      end


      subroutine updategr(grbin,totgr,delgr,Ncol,Lx,Ly)

      double precision grbin(100),totgr,delgr
      integer Ncol,Lx,Ly,i
      double precision vol,norm,pi,r

      pi = 4.d0*datan(1.d0)
      norm = dble(Lx*Ly) /
     :    (0.5d0*dble(Ncol*(Ncol-1))*pi*totgr*delgr**2)

      open(unit=13,file='gr.dat')
      do i = 1,100
         r = (dble(i)-0.5d0)*delgr
         vol = dble(i**2) - dble((i-1)**2)
         write(13,'(2(D10.4,1X))') r,grbin(i)*norm/vol
      enddo
      close(13)

      return
      end


      subroutine order(rminsq,Ncol,Ncolup,Lx,step,dt,direct,grav)

      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      integer step,Lx,Ncol,Ncolup,i,i2,j,j2
      double precision dt,time,direct,grav,phitot,rminsq
      double precision rijsq,xij,xi,rLx,rLx2

      time = dt*dble(step)
      rLx = dble(Lx)
      rLx2 = rLx/2.d0
      phitot = 0.d0
      do 10 i = 1,Ncolup
         i2 = 2*(i-1)
         xi = Rcol(i2+1)
         do j = Ncolup+1,Ncol
            j2 = 2*(j-1)
            xij = xi - Rcol(j2+1)
            if (xij.gt.rLx2) then
               xij = xij-rLx
            elseif (xij.le.-rLx2) then
               xij = xij+rLx
            endif
            rijsq = xij**2
            if (rijsq.le.rminsq) goto 10 !at least one is close
         enddo
         phitot = phitot+1.d0
10    continue
      do 20 i = Ncolup+1,Ncol
         i2 = 2*(i-1)
         xi = Rcol(i2+1)
         do j = 1,Ncolup
            j2 = 2*(j-1)
            xij = xi - Rcol(j2+1)
            if (xij.gt.rLx2) then
               xij = xij-rLx
            elseif (xij.le.-rLx2) then
               xij = xij+rLx
            endif
            rijsq = xij**2
            if (rijsq.le.rminsq) goto 20 !at least one is close
         enddo
         phitot = phitot+1.d0
20    continue
      open(unit=15,file='order.dat',position='append')
      write(15,'(3D20.10)') time,phitot/dble(Ncol),-direct*grav
      close(15)

      return
      end
      
      subroutine vel_avg(Lx,Ly)

      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      common / block2 / list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer,pointer,dimension(:) :: list,listdepl,biglistdepl,listcol,
     :                  head,headdepl,bigheaddepl,headcol,map,mapcol
      integer Lx,Ly,i,j,j1,k,ipar,plane_pos
      integer, save :: t_count=0,file_count =0
      double precision,save :: vx_avg(102) =0.0d0, vy_avg(102)=0.0d0
      double precision, save :: tot_part(102) = 0.0d0 
      double precision weight
      character(len=100) :: fname
      
      t_count = t_count + 1
      DO i=0,Lx
      DO j1=0,Ly
	ipar = head(1+i+j1*(Lx+1))					
	do while (ipar/=0)
                j = 2*ipar-1 
                plane_pos = i
                weight = R(j)-plane_pos
                k = i+1
		vx_avg(k)   = vx_avg(k)   + V(j)*(1-weight)
		vx_avg(k+1) = vx_avg(k+1) + V(j)*weight
		vy_avg(k)   = vy_avg(k)   + V(j+1)*(1-weight)
		vy_avg(k+1) = vy_avg(k+1) + V(j+1)*weight
		tot_part(k) = tot_part(k) + (1-weight)
		tot_part(k+1) = tot_part(k+1) + weight
		ipar = list(ipar)
	end do
      END DO
      END DO
      if (modulo(t_count,100000)==0) then
	file_count = file_count+1
	! Averaging all the velocities over the respective particles and time
	vx_avg = vx_avg/tot_part
	vy_avg = vy_avg/tot_part
	
	write (fname, "(A12,I0,A4)") "Poise_vx_vy_",file_count,".dat"         
	open (unit=20,file=trim(fname),action="write",status="replace")
	
	DO i=1,Lx+1
		write(20,*) (i*1.0-1),vx_avg(i),vy_avg(i)
	END DO
        close(20)	
	
	vx_avg = 0.0d0
	vy_avg = 0.0d0
	tot_part=0.0d0
	t_count = 0	
        end if    		
      return
      end


