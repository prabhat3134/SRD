      subroutine pdb(movie,Lx,Ly,Ncol,Ncolup,step)

      logical movie
      integer Lx,Ly,Ncol,Ncolup,step,i,i2,m
      common / block1 / R,V,F,Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol
      double precision,pointer,dimension(:) :: R,V,F,
     :                        Rdepl,Vdepl,Fdepl,Rcol,Vcol,Fcol

      if (movie) then
         open(unit=99,file='config.pdb',position='append')
      else
         open(unit=99,file='config.pdb')
      endif

      write(99,900) 'REMARK TIMESTEP',step
      write(99,1000) 'CRYST1',dble(Lx),dble(Ly),1.d0,
     :   90.0,90.0,90.0,'P 1',1
 900  format (A15,1X,I10)
1000  format (A6,3F9.3,3F7.2,1X,A3,8X,I3)

      do i = 1,Ncolup
         i2 = 2*(i-1)
         write(99,1020) 'HETATM',i,'C',1,'UNK',1,
     :      (Rcol(i2+m),m=1,2),0.d0,1.0,1.0
      enddo
      do i = Ncolup+1,Ncol
         i2 = 2*(i-1)
         write(99,1020) 'HETATM',i,'S',1,'UNK',1,
     :      (Rcol(i2+m),m=1,2),0.d0,1.0,1.0
      enddo
1020  format (A6,I5,2X,A1,I2,1X,A3,2X,I4,4X,3F8.2,2F6.2)

      write(99,1070) 'END'
1070  format(A3)

      close(unit=99)

      return
      end
