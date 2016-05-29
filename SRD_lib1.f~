module SRD_lib1

implicit none

INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828

 contains
!*******************************************
!	To generate random seed
subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end subroutine init_random_seed

!*****************************************************************************
! random generator
FUNCTION ran()
		IMPLICIT NONE
			INTEGER, PARAMETER :: K4B=selected_int_kind(9)			
			REAL :: ran
		!“Minimal” random number generator of Park and Miller combined with a Marsaglia shift
		!sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
		!values). This fully portable, scalar generator has the “traditional” (not Fortran 90) calling
		!sequence with a random deviate as the returned function value: call with idum a negative
		!integer to initialize; thereafter, do not alter idum except to reinitialize. The period of this
		!generator is about 3.1 × 10^18.
		INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
		REAL, SAVE :: am
		INTEGER(K4B), SAVE :: ix=-1,iy=-1,k, idum =-1

		if (idum <= 0 .or. iy < 0) then 			!Initialize.
			am=nearest(1.0,-1.0)/IM
			iy=ior(ieor(888889999,abs(idum)),1)
			ix=ieor(777755555,abs(idum))
			idum=abs(idum)+1 				!Set idum positive.		
		end if

		ix=ieor(ix,ishft(ix,13)) 			!Marsaglia shift sequence with period 232 − 1.
		ix=ieor(ix,ishft(ix,-17))	
		ix=ieor(ix,ishft(ix,5))
		k=iy/IQ 					!Park-Miller sequence by Schrage’s method,
		iy=IA*(iy-k*IQ)-IR*k 

		if (iy < 0) iy=iy+IM
		ran=am*ior(iand(IM,ieor(ix,iy)),1) 	!Combine the two generators with masking to
END FUNCTION ran


!***********************************************************************************
! For eliminating particle within a box
subroutine box_eliminate(rx,ry,rx1,ry1,l0,b0,width)
implicit none

real(kind=dp), dimension(:) :: rx,ry
integer ::  l0, b0, width
logical :: l1(size(rx)), l2(size(rx))
real(kind=dp), dimension(:) , allocatable::rx1,ry1

l1 = ry<b0+width .and. ry>b0 .and. rx<l0+width .and. rx>l0
l2 = .not. l1
rx1 = pack(rx,l2)
ry1 = pack(ry,l2)

end subroutine box_eliminate

!******************************************
! Getting normal distribution from uniform distribution

subroutine random_normal(vel,np,av,std)
implicit none
integer  :: np, i
real(kind=dp) :: vel(np), v1(np), v2(np)
real :: av,std

DO i=1,np
v1(i) = ran()
v2(i) = ran()
END DO

! Box Muller Transformation
vel = std*sqrt(-2*log(v1))*cos(2*pi*v2)+av

end subroutine random_normal

!******************************************
! getting weibull distribution for velocity
subroutine random_weibull(vel,np,kbt)
implicit none
real(kind=dp) :: vel(np)
integer :: np,i
real :: a,b,v,Pv, vmax, vprob, Pmax, kbt, P
i=1
a = sqrt(2*kbt)
b = 2.0
vprob = sqrt(kbt)		! most probable speed
vmax = 5*vprob			! maximum speed possible
Pmax = 1/sqrt(kbt*e)		!maximum probability

do while (i<=np)	
	v = ran()*vmax
	Pv = ( b/(a**b) ) * v**(b-1) * EXP( -1*(v/a)**b )	
	P = ran()*Pmax	
	if (Pv>P) then
		vel(i)=v
		i=i+1
	end if	
end do
end subroutine random_weibull

!******************************************
! Intialize the domain
subroutine initialize(rx,ry,vx,vy, Ly, Lx, np, av, std)
implicit none
INTEGER :: np, Ly, Lx, i
REAL :: av, std
REAL(kind=dp) :: rx(np), ry(np), vx(np), vy(np)

do i=1, np
rx(i) = ran()*Lx
ry(i) = ran()*Ly
end do

call random_normal(vx,np,av,std)
call random_normal(vy,np,av,std)

end subroutine initialize

!*******************************************
! Cell partition for dividing particles into grid cells
subroutine partition(rx,ry,head,list,Lx,Ly,np)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx, ry
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: np, Lx, Ly, xindex, yindex,i
head = 0
list = 0

do i=1,np
	yindex = ceiling(ry(i))	
	xindex = ceiling(rx(i))		
	if (yindex==0) yindex=1 
	if (xindex==0) xindex=1
	! linked list algorithm
	if (head(yindex,xindex)==0) then
		head(yindex,xindex)=i
	else 
		list(i) = head(yindex,xindex)
		head(yindex,xindex)=i
	end if	
end do
end subroutine partition

!*******************************************
! Collision Scheme
subroutine collision(vx, vy, temp, tempx, tempy, head, list, Lx, Ly, alpha)
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: i,j,count1, ipar, jpar, Lx, Ly, out_unit, k
real	:: r, alpha, vx_temp, vy_temp, alpha1
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), temp(Ly,Lx), tempx(Ly,Lx), tempy(Ly,Lx)


vxcom=0.0
vycom=0.0
tempx=0.0
tempy=0.0
temp =0.0

do i=1,Ly	
	do j=1,Lx
		count1 = 0
				
		r = ran()
		
		if (r >=0.5) then
			alpha1 = -alpha
		else 
			alpha1 = alpha
		end if			
		
		ipar = head(i,j)					
		do while (ipar/=0)			
			vxcom(i,j) = vxcom(i,j) + vx(ipar)
			vycom(i,j) = vycom(i,j) + vy(ipar)
			count1 = count1 + 1
			ipar = list(ipar)
		end do
		if (count1/=0) then
			vxcom(i,j) = vxcom(i,j)/count1
			vycom(i,j) = vycom(i,j)/count1
		end if		
		
		! Rotation of velocity
		ipar = head(i,j)
		do while (ipar/=0)
			!tempx(i,j) = tempx(i,j) + ((vx(ipar)-vxcom(i,j))**2)/2 
			!tempy(i,j) = tempy(i,j) + ((vy(ipar)-vycom(i,j))**2)/2
			temp(i,j) = temp(i,j) + (vx(ipar)**2 + vy(ipar)**2) 	
			vx_temp  = vx(ipar)-vxcom(i,j) 
			vy_temp  = vy(ipar)-vycom(i,j) 			
			vx(ipar) = vxcom(i,j) + cos(alpha1)*vx_temp + sin(alpha1)*vy_temp
			vy(ipar) = vycom(i,j) - sin(alpha1)*vx_temp + cos(alpha1)*vy_temp			
			ipar = list(ipar)			
		end do	
		temp(i,j) = temp(i,j)/count1
		!tempx(i,j) = tempx(i,j)/(count1-1)
		!tempy(i,j) = tempy(i,j)/(count1-1) 				
	end do
	
end do		
	
end subroutine collision

!****************************************************************
! Another subroutine for thermal boundary conditions
subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt, np)
implicit none
integer, parameter :: np1=1000		! Expected no. of particles crossing boundary
real(kind=dp), dimension(:) :: rx, ry, rx1, ry1, vx, vy
real :: g, dt_c, kbt, av, std, t_app, t_dec
integer :: Ly, i, np, j
logical :: Wall(np)
real(kind=dp) :: vxwall(np1), vywall(np1)

std   = sqrt(kbt)
j = 0
av =0.0

Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)
call random_normal(vxwall,np1,av,std)
call random_weibull(vywall,np1,kbt)

do i = 1,np
	if (Wall(i)) then	
		if ((ry1(i)-ry(i))<0.0) then
			j = j+1
			t_app = ry(i)/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = vywall(j)
	      	ry1(i) = 0.0 + vy(i)*t_dec	
		else 
			j = j+1
			t_app = (Ly-ry(i))/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vywall(j)
	      	ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i)  = vx(i)-g*dt_c
		rx1(i) = rx(i) + vx(i)*t_app + (g*t_app**2)/2
		!rx1(i)  = rx(i) + vx(i)*t_app     
		rx1(i) = rx1(i) + vxwall(j)*t_dec + (g*t_dec**2)/2
      	vx(i)  = vxwall(j) + g*t_dec
      	!rx1(i)  = rx1(i) + vx(i)*t_dec       	       				
	end if
	if (j==np1) then		! if more particles than expected are crossing, generate new numbers
		call random_normal(vxwall,np1,av,std)
		call random_weibull(vywall,np1,kbt)
		j=0	
	end if
end do

end subroutine thermal_wall

!***********************************************************************
! Writing data to files for analysis
subroutine average(rx,ry,vx,vy,head,list, density, Lx ,Ly)
implicit none
real(kind=dp), dimension(:) :: rx, ry, vx, vy
integer , dimension(:,:)    :: head
integer, dimension (:)      :: list
integer :: Lx, Ly, i, j, den=0, ipar
real(kind=dp) :: density(Ly,Lx), vx_com(Ly,Lx), vy_com(Ly,Lx)
density(Ly,Lx)=0.0
vx_com(Ly,Lx)=0.0
vy_com(Ly,Lx)=0.0

do j=1, Lx
	do i=1, Ly
		ipar = head(i,j)		
		do while (ipar/=0)
			vx_com(i,j) = vx_com(i,j) + vx(ipar)
			vy_com(i,j) = vy_com(i,j) + vy(ipar)
			den = den +1
			ipar = list(ipar)
		end do
		if (den/=0) then		
		density(i,j) = den
		vx_com(i,j)  = vx_com(i,j)/den	
		vy_com(i,j)  = vy_com(i,j)/den
		end if	
	end do
end do

end subroutine average
!***********************************************************************
! Subroutine for bounce back boundary condition at wall
subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt,np)
implicit none
real(kind=dp), dimension(np)  :: rx, ry
real(kind=dp), dimension(np)  :: rx1, ry1, vx, vy
real :: g, dt_c, kbt, t_app, t_dec
integer :: Ly, i, np
logical :: Wall(size(ry))

Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)  			  		

do i = 1,np
	if (Wall(i)) then			
		if ((ry1(i)-ry(i))<0.0) then			
			t_app = ry(i)/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vy(i)
	      	ry1(i) = 0.0 + vy(i)*t_dec	
		else 			
			t_app = (Ly-ry(i))/abs(vy(i))
	      	t_dec = dt_c - t_app
	      	vy(i) = -vy(i)
	      	ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i)  = vx(i)-g*(dt_c-t_app)
		rx1(i)  = rx(i) + vx(i)*t_app      
      		vx(i)  = -vx(i) + g*t_dec
      		rx1(i)  = rx1(i) + vx(i)*t_dec  		
	end if	
end do

end subroutine bounce_wall



end module SRD_lib1
