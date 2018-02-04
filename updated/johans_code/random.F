      double precision function randy()

c     start_doc
c     Name
c       randy
c       Adapted by Johan Padding for use with double precision.
c
c     Usage
c       integer i, j, k, l
c       double precision random, randy
c       call rstart (i, j, k, l)     ! seed the generator once and for all
c       random = randy()
c
c     Description
c       Universal portable random number generator.  Uniformly
c       distributed random numbers in [0, 1).
c       Ref. G Marsaglia, A Zaman and W W Tsang, Statistics and
c       Probability Letters 8, 35 (1990)
c
c       Correct performance of this implementation was tested using the
c       test driver of table 3 of the reference.
c
c     Parameters
c       i, j, k, l   integer seeds: i, j, k: 1-178 and not all 1
c                                   l:       0-168
c
c     Memo
c       random number generator
c
c     Imported
c
c     Author
c       fmp nov 90, hacked from above reference
c
c     end_doc

      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      double precision u(97)
      common / crandy / u, c, cd, cm, ip, jp

      randy = u(ip) - u(jp)
      if (randy .lt. 0.d0) randy = randy + 1.d0
      u(ip) = randy

      ip = ip - 1
      if (ip .eq. 0) ip = 97
      jp = jp - 1
      if (jp .eq. 0) jp = 97

      c = c - cd
      if (c .lt. 0.0) c = c + cm

      randy = randy - c
      if (randy .lt. 0.d0) randy = randy + 1.d0
      return
      end


      subroutine rstart (iii, jjj, kkk, lll)

      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      double precision u(97)
      common / crandy / u, c, cd, cm, ip, jp

      i = iii
      j = jjj
      k = kkk
      l = lll

      do 2 ii = 1, 97
         s = 0.0
         t = 0.5
         do 3 jj = 1, 24
            m = mod (mod (i * j, 179) * k, 179)
            i = j
            j = k
            k = m
            l = mod (53 * l + 1, 169)
            if (mod (l * m, 64) .ge. 32) s = s + t
            t = 0.5 * t
 3       continue
         u(ii) = s
 2    continue

      c  =   362436.0 / 16777216.0
      cd =  7654321.0 / 16777216.0
      cm = 16777213.0 / 16777216.0

      ip = 97
      jp = 33

      return
      end


      subroutine loadranstate

      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      double precision u(97)
      common / crandy / u, c, cd, cm, ip, jp

      logical ex
   
      inquire(file='ranstate',exist=ex)
      if (ex) then
         open(unit=40,file='ranstate',status='old',form='unformatted')
         read(40) u,c,cd,cm,ip,jp
         close(40)
      else
         iii = 13
         jjj = 121
         kkk = 45
         lll = 43
         call rstart(iii,jjj,kkk,lll)
      endif

      return
      end


      subroutine saveranstate

      implicit double precision (a-h, o-z)
      implicit integer (i-n)
      double precision u(97)
      common / crandy / u, c, cd, cm, ip, jp

      open(unit=40,file='ranstate',status='unknown',form='unformatted')
      write(40) u,c,cd,cm,ip,jp
      close(40)

      return
      end


        DOUBLE PRECISION FUNCTION GAUSS ( DUMMY )

C    *******************************************************************
C    ** RANDOM VARIATE FROM THE STANDARD NORMAL DISTRIBUTION.         **
C    **                                                               **
C    ** THE DISTRIBUTION IS GAUSSIAN WITH ZERO MEAN AND UNIT VARIANCE.**
C    **                                                               **
C    ** REFERENCE:                                                    **
C    **                                                               **
C    ** KNUTH D, THE ART OF COMPUTER PROGRAMMING, (2ND EDITION        **
C    **    ADDISON-WESLEY), 1978                                      **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** REAL FUNCTION randy()                                         **
C    **    RETURNS A UNIFORM RANDOM VARIATE ON THE RANGE ZERO TO ONE  **
C    *******************************************************************

        DOUBLE PRECISION A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138d0, A3 = 0.252408784d0 )
        PARAMETER ( A5 = 0.076542912d0, A7 = 0.008355968d0 )
        PARAMETER ( A9 = 0.029899776d0                   )

        DOUBLE PRECISION SUM, R, R2
        DOUBLE PRECISION randy, DUMMY
        INTEGER     I

C    *******************************************************************

        SUM = 0.d0

        DO 10 I = 1, 12

           SUM = SUM + randy()

10      CONTINUE

        R  = ( SUM - 6.d0 ) / 4.d0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 )
     :          * R

        RETURN
        END


        SUBROUTINE BIVAR ( VAR1, VAR2, COVAR, G1, G2 )

C    ** TWO NUMBERS G1 AND G2 ARE SAMPLED FROM A BIVARIATE GAUSSIAN
C    ** DISTRIBUTION WITH SD's VAR1 AND VAR2 AND CORRELATION COEFFICIENT
C    ** (COVARIANCE) COVAR

        DOUBLE PRECISION VAR1, VAR2, COVAR, G1, G2
        DOUBLE PRECISION GAUSS1, GAUSS2, DUMMY, GAUSS

        GAUSS1 = GAUSS ( DUMMY )
        GAUSS2 = GAUSS ( DUMMY )

        G1 = VAR1 * GAUSS1
        G2 = VAR2 * ( COVAR*GAUSS1 + DSQRT(1.D0-COVAR*COVAR) * GAUSS2 )

        RETURN
        END

       

