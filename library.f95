module SRD_library
use omp_lib
implicit none
INTEGER, PARAMETER :: dp = selected_real_kind(15, 307), Gama = 10, m=1, a0=1
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828d0
! Grid Parameters
INTEGER, PARAMETER :: Ly = 160, Lx = 1600, np = Ly*Lx*Gama, half_plane = 1
REAL(kind=dp), PARAMETER ::  alpha = pi/2.0d0, kbT = 1.0d0, dt_c =0.1d0
! Forcing 
REAL(kind=dp) :: avg=0.0d0, std=sqrt(kbT/(m*1.0d0)), f_b = 5.0d-4
! time values
INTEGER :: tmax=200000, t_avg = 500000, avg_interval=1, ensemble_num = 100000
! RGS, streaming
INTEGER :: random_grid_shift = 1, verlet = 1, grid_up_down
! Thermostat
INTEGER :: mb_scaling = 0, MC_scaling = 1, mbs_freq=50
REAL(kind=dp)    :: force(2), mu_tot, MC_strength = 0.25d0
LOGICAL :: xy(2)=[.TRUE., .FALSE.], temperature = .TRUE., wall_thermal = .FALSE.
! File naming 
CHARACTER(len=100) :: file_name='Poiseuille_flow', data_path='./'     !file_name of size 20
! cylinder parameters
INTEGER :: obst = 1, grid_check(0:Ly+1,Lx)=0 
REAL(kind=dp) :: rad = 10d0, xp = Lx/4.0d0, yp = Ly/2.0d0
REAL(kind=dp),ALLOCATABLE :: theta_intersect(:)   
LOGICAL, ALLOCATABLE ::  obst_par(:)
!! We should not define  very large static arrays typically of np. 
!! Use dynamic allocation instead for such arrays and keep stack size unlimited.
!! This precaution is important when using openmp and can be ignored without openmp.
contains
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
		INTEGER(K4B), SAVE :: ix=-1,iy=-1,k, idum =-2, call_count=0
		IF (call_count == 0) then 
			call clockx(idum)
			idum = -modulo(idum,10000)
		END IF
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
		call_count = call_count + 1
END FUNCTION ran


!***********************************************************************************
! For eliminating particle within a box
subroutine box_eliminate(x_dummy, y_dummy)
implicit none
integer :: i
real(kind=dp), dimension(:) :: x_dummy, y_dummy
logical, allocatable :: l1(:), l2(:)

ALLOCATE(l1(np), l2(np))

i = 1
! l1 = ry<b0+width .and. ry>b0 .and. rx<l0+width .and. rx>l0
l1 = (x_dummy - xp)**2 + (y_dummy - yp)**2 < rad**2
l2 = .not. l1

do while (i <= size(x_dummy))
	if (l1(i)) then
		x_dummy(i) = ran()*Lx
		y_dummy(i) = ran()*Ly

		if ((x_dummy(i) - xp)**2 + (y_dummy(i) - yp)**2 > rad**2) then 
			i = i+1
		end if
	else 
		i = i+1
	end if
end do

DEALLOCATE(l1,l2)
end subroutine box_eliminate

!******************************************
! Getting normal distribution from uniform distribution
subroutine random_normal(vel,p_count,av0,std0)
implicit none
integer  :: p_count, i
real(kind=dp) :: vel(p_count), v1(p_count), v2(p_count)
real(kind=dp) :: av0,std0

DO i=1,p_count
v1(i) = ran()
v2(i) = ran()
END DO

! Box Muller Transformation
vel = std0*sqrt(-2*log(v1))*cos(2*pi*v2)+av0

end subroutine random_normal

!******************************************
! getting weibull distribution for velocity
subroutine random_weibull(vel,p_count)
implicit none
integer :: p_count,i
real(kind=dp) :: vel(p_count)
real :: a,b,v,Pv, vmax, vprob, Pmax, P
i=1
a = sqrt(2*kbT)
b = 2.0
vprob = sqrt(kbT)		! most probable speed
vmax = 5*vprob			! maximum speed possible
Pmax = 1/sqrt(kbT*e)		!maximum probability

do while (i<=p_count)	
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
subroutine initialize(x_dummy, y_dummy, rx,ry,vx,vy, head, list)
implicit none
INTEGER ::  i, j, ipar, p_count, head(:,:), list(:)
REAL(kind=dp) ::  block(3), vp_max 
REAL(kind = dp) :: x_dummy(:), y_dummy(:), mu_kin, mu_col
REAL(kind=dp) :: rx(:), ry(:), vx(:), vy(:), vxcom(Ly,Lx),vycom(Ly,Lx)

mu_kin = (gama*kbT*dt_c)*(gama/((gama - 1.0 + exp(-1.0*gama))*(1.0-cos(2.0*alpha)))- 0.5)
mu_col = ((1.0-cos(alpha))/(12.0*dt_c))*(gama - 1.0 + exp(-1.0*gama))
mu_tot = mu_kin + mu_col
vp_max = (Gama * Ly**2.0 *f_b)/(8.0*mu_tot)
!write(*,*) vp_max, mu_tot
!stop
block = [xp,yp,rad]
do i=1, np
	x_dummy(i) = ran()*Lx
	y_dummy(i) = ran()*Ly
end do
IF (obst == 1) THEN 
	ALLOCATE(obst_par(np))
	ALLOCATE(theta_intersect(np))
	call box_eliminate(x_dummy, y_dummy)
	call gridcheck(grid_check, block)
END IF
rx = x_dummy
ry = y_dummy

call random_normal(vx,np,avg,std)
call random_normal(vy,np,avg,std)

call partition(rx,ry,head,list)
vxcom = 0.0d0
vycom = 0.0d0
DO j=1,Lx
DO i=1,Ly
	p_count = 0
	ipar = head(i,j)					
	do while (ipar/=0)			
		vxcom(i,j) = vxcom(i,j) + vx(ipar)
		vycom(i,j) = vycom(i,j) + vy(ipar)
		p_count = p_count + 1
		ipar = list(ipar)
	end do
	if (p_count/=0) then
		vxcom(i,j) = vxcom(i,j)/p_count
		vycom(i,j) = vycom(i,j)/p_count
	end if	

	ipar = head(i,j)
	do while (ipar/=0)
		vx(ipar) = vx(ipar) - vxcom(i,j) 
		vy(ipar) = vy(ipar) - vycom(i,j)
		ipar = list(ipar)
	end do
END DO
END DO
end subroutine initialize

!*********************************************
! Streaming step in SRD
subroutine streaming(rx, ry, rx1, ry1, vx, vy)
implicit none
integer :: i
real(kind = dp),dimension(:) :: rx, ry, rx1, ry1, vx, vy

IF (verlet == 2) THEN	 	! LEAPFROG algorithm
        !$OMP PARALLEL DO
        DO i=1,np
                vx(i)  = vx(i) + f_b*dt_c
                rx1(i) = rx(i) + vx(i)*dt_c  
                ry1(i) = ry(i) + vy(i)*dt_c
        END DO
        !$OMP END PARALLEL DO
ELSE IF (verlet ==1) THEN	!Verlet Algorithm
        !$OMP PARALLEL DO
        DO i=1,np
                rx1(i) = rx(i) + vx(i)*dt_c + (f_b*dt_c**2)/2.0 
                ry1(i) = ry(i) + vy(i)*dt_c
                vx(i)  = vx(i) + f_b*dt_c
        END DO
        !$OMP END PARALLEL DO
ELSE				!EULER algorithm
        !$OMP PARALLEL DO
        DO i=1,np
                rx1(i) = rx(i) + vx(i)*dt_c 
                ry1(i) = ry(i) + vy(i)*dt_c
                vx(i)  = vx(i) + f_b*dt_c
        END DO
        !$OMP END PARALLEL DO
END IF
IF (obst==1) obst_par  = (rx1 - xp)**2 + (ry1 - yp)**2 < rad**2

end subroutine streaming
!*************************************************
! Find the intersection angles, and, the time required to intersect, for particles which cross the cylinder boundary
! Enters the module every time step
! Using linear interpolation between initial and final positions, which can be subjected to change
subroutine par_in_cyl( rx, ry, rx1, ry1, vx, vy )
implicit none

integer :: i
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1, vx, vy
real(kind =dp) :: dx, dy, dr, de, disc, sol1(2), sol2(2), th, sol(2) 
real(kind = dp) :: exp1(2), exp2(2) 

theta_intersect = 0.0d0
! finding the solution of y - mx -c = 0 and x^2 + y^2 = r^2 using the quadratic formula 
! Solving the intersection angle for every particle inside the cylinder (using the list obst_par(i)) 

!$OMP PARALLEL IF(np>100000)
!$OMP DO PRIVATE(i, dx, dy, dr, de, disc, sol1, sol2, sol, exp1, exp2) SCHEDULE(guided)
do i=1,np
	if(obst_par(i)) then
! Calculating the slope, of the line joining the initial and final particle positions. 
! for reference, check: 	http://mathworld.wolfram.com/Circle-LineIntersection.html	
! The line is shifted to origin. 
		dx = rx1(i) - rx(i)
		dy = ry1(i) - ry(i)
		dr = SQRT(dx**2 + dy**2)
		de = (rx(i) -xp)*(ry1(i)-yp) - (rx1(i)-xp)*(ry(i) -yp)
		disc = (rad**2)*(dr**2)-(de**2) 
		
		if (disc > 0) then 
			sol1(1) = (de*dy + SIGN(1.0d0,dy)*dx*SQRT(disc))/(dr**2)
			sol2(1) = (de*dy - SIGN(1.0d0,dy)*dx*SQRT(disc))/(dr**2)
			sol1(2) = (-de*dx + ABS(dy)*SQRT(disc))/(dr**2)
			sol2(2) = (-de*dx - ABS(dy)*SQRT(disc))/(dr**2)
		else
			write(*,*) "Negative Discriminant: particle doesn't collide",rx(i), ry(i), rx1(i), ry1(i)
		endif

		exp1(1) = (sol1(1) + xp - rx1(i))*(sol1(1) + xp - rx(i)) 
		exp2(1) = (sol2(1) + xp - rx1(i))*(sol2(1) + xp - rx(i))
		exp1(2) = (sol1(2) + yp - ry1(i))*(sol1(2) + yp - ry(i)) 
		exp2(2) = (sol2(2) + yp - ry1(i))*(sol2(2) + yp - ry(i)) 
		
		if ( exp1(1) < 0.0 .and. exp1(2) < 0.0 ) then
			sol = sol1
		elseif ( exp2(1) < 0.0 .and. exp2(2) < 0.0 ) then 
			sol = sol2
		else
			write(*,*) "intersection algorithm fails",rx(i), ry(i), rx1(i), ry1(i)
		endif 

		th = atan2(sol(2), sol(1))
		theta_intersect(i) = merge(th, th + 2*pi, th >= 0)
	
	endif
end do
!$OMP END DO
!$OMP END PARALLEL


call pos_thermal_update(rx,ry,rx1,ry1,vx,vy)

end subroutine par_in_cyl
!*********************************************
! based on the angle theta_intersect of the particles inside the cylinder, place them back to the right position
! Velocity in radial coordinate is [radial, tangential] and the Euler algorithm is implemented by default
subroutine pos_bounce_update(rx, ry, rx1, ry1, vx, vy)

implicit none
integer :: i, j
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1,  vx, vy
real(kind=dp) :: t_app, t_sep, jacobian(2,2), th, vr(2), vxy(2), force_r(2)

force = 0.0d0

!$OMP PARALLEL IF(np>100000)
!$OMP DO PRIVATE(th, jacobian, t_app, t_sep, vxy, vr, force_r ) SCHEDULE(guided) REDUCTION(+:force)
DO i=1,np 
	if(obst_par(i)) then 
		th = theta_intersect(i)
		jacobian = reshape([cos(th),- sin(th), sin(th), cos(th)],shape(jacobian))
		t_app = abs(yp+rad*sin(th) - ry(i))/abs(vy(i))
		t_sep = dt_c - t_app
	
		if (t_sep < 0)  write(*,*) 'Cylinder collision: negative time', rx(i),ry(i),rx1(i),ry1(i),vx(i),vy(i),th 
		! x velocity changes in approach -> reversal -> x velocity changes in separation	
		vx(i) = vx(i) - f_b*dt_c
		
		! velocity of approach at the surface of the cylinder
		vx(i) = vx(i) + f_b*t_app
		
		! coordinate transformation for reversal in dirn of normal and tangential velocity and propogation
		vxy = [vx(i), vy(i)]
		vr  = MATMUL(jacobian,vxy)
		vr  = -vr
		vxy = MATMUL(transpose(jacobian), vr)
		
		! Force calculation by change in momentum of each striking particle
		force_r = -2*vr
		force   = force +   MATMUL(transpose(jacobian), force_r)

		! velocity of separation at the surface of the cylinder (reversal)
		vx(i) = vxy(1)
		vy(i) = vxy(2)
		
		! propogation from the point of intersection
		rx1(i) = xp+ rad*cos(th) + vx(i)*t_sep
		ry1(i) = yp+ rad*sin(th) + vy(i)*t_sep
		
		! updating velocity for the next time step
		vx(i)  = vx(i) + f_b*t_sep

	endif
end do
!$OMP END DO
!$OMP END PARALLEL
end subroutine pos_bounce_update

!*********************************************
! based on the angle theta_intersect of the particles inside the cylinder, place them back to the right position
! Velocity in radial coordinate is [radial, tangential] and the Euler algorithm is implemented by default
subroutine pos_thermal_update(rx, ry, rx1, ry1, vx, vy)

implicit none
integer :: i, j
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1, vx, vy
real(kind=dp) :: t_app, t_sep, jacobian(2,2), th, v_polar_app(2), v_polar_sep(2), vxy(2), force_polar(2), vxwall(1), vywall(1)
force = 0.0d0

!$OMP PARALLEL IF(np>100000)
!$OMP DO PRIVATE(th, jacobian, t_app, t_sep, vxy, v_polar_app, v_polar_sep, force_polar, vxwall, vywall ) SCHEDULE(guided) REDUCTION(+:force)
do i=1, np
	if(obst_par(i)) then 
		th = theta_intersect(i)
		jacobian = reshape([cos(th),- sin(th), sin(th), cos(th)],shape(jacobian))
		t_app = abs(yp+rad*sin(th) - ry(i))/abs(vy(i))
		t_sep = dt_c - t_app
	
		if (t_sep < 0)  write(*,*) 'Cylinder collision: negative time', rx(i),ry(i),rx1(i),ry1(i),vx(i),vy(i),th 
		! x velocity changes in approach -> reversal -> x velocity changes in separation	
		vx(i) = vx(i) - f_b*dt_c
		! velocity of approach at the surface of the cylinder
		 vx(i) = vx(i) + f_b*t_app
		
		! coordinate transformation for reversal in dirn of normal and tangential velocity and propogation
		vxy = [vx(i), vy(i)]
		v_polar_app  = MATMUL(jacobian,vxy)
		call random_weibull(vywall,1)
		call random_normal(vxwall,1,avg,std)
		v_polar_sep  = [vywall(1), vxwall(1)] 
		vxy = MATMUL(transpose(jacobian), v_polar_sep)
		
		! Force calculation by change in momentum of each striking particle
		force_polar = -1.0*( v_polar_sep - v_polar_app)
		force   = force +   MATMUL(transpose(jacobian), force_polar)

		! velocity of separation at the surface of the cylinder (reversal)
		vx(i) = vxy(1)
		vy(i) = vxy(2)
		
		! propogation from the point of intersection
		rx1(i) = xp+ rad*cos(th) + vx(i)*t_sep 
		ry1(i) = yp+ rad*sin(th) + vy(i)*t_sep
		
		! updating velocity for the next time step
		vx(i)  = vx(i) + f_b*t_sep
	endif
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine pos_thermal_update
!*************************************************
! implementing the periodic boundary condition
! Periodic boundary condition in x and y as per the condition given by xy
subroutine periodic_xy( rx1, ry1)
implicit none
real(kind = dp) :: rx1(:), ry1(:) 
INTEGER :: i

if (xy(1)) then
        !$OMP PARALLEL DO
        DO i=1,np
                IF (rx1(i) < 0.0) THEN
                        rx1(i) = rx1(i) + Lx
                ELSE IF (rx1(i) > Lx*1.0) THEN
                        rx1(i) = rx1(i) - Lx
                END IF
        END DO
        !$OMP END PARALLEL DO
end if

if (xy(2)) then
	!$OMP PARALLEL DO
        DO i=1,np
                IF (ry1(i) < 0.0) THEN
                        ry1(i) = ry1(i) + Ly
                ELSE IF (ry1(i) > Ly*1.0) THEN
                        ry1(i) = ry1(i) - Ly
                END IF
        END DO
        !$OMP END PARALLEL DO
end if

end subroutine periodic_xy
!********************************************************
! Grid check, gives which grids are overlapping with obstacles for further use in collision
subroutine gridcheck(grid_check, block)
implicit none
INTEGER :: grid_check(:,:), i, j, x, y
REAL(kind=dp) :: block(:), xc, yc, l, b, R, rtmp

grid_check(0:Ly+1,:) = 0
if (SIZE(block)==3) then	! for cylinder [center x, center y, radius]
	xc = block(1)
	yc = block(2)
	R  = block(3)
	DO i = CEILING(yc-R)+1,CEILING(yc+R)+1	! i is in y direction
	DO j = CEILING(xc-R),CEILING(xc+R)	! j is in x direction
		rtmp = SQRT((j-xc)**2 + (i-1-yc)**2)
		if (rtmp < R) then
			grid_check(i,j) = 1
			cycle
		end if
		rtmp = SQRT((j-xc)**2 + (i-2-yc)**2)
		if (rtmp < R) then
			grid_check(i,j) = 1
			cycle
		end if
		rtmp = SQRT((j-1-xc)**2 + (i-1-yc)**2)
		if (rtmp < R) then
			grid_check(i,j) = 1
			cycle
		end if
		rtmp = SQRT((j-1-xc)**2 + (i-2-yc)**2)
		if (rtmp < R) then
			grid_check(i,j) = 1
			cycle
		end if
	END DO
	END DO
else
	xc = block(1)
	yc = block(2)
	l  = block(3)
	b  = block(4)
	! for including square or rectangular block [base x, base y, length, breadth]
end if
end subroutine gridcheck

!*******************************************
! Cell partition for dividing particles into grid cells
! OPENMP not implemented as it severely reduces the performance
subroutine partition(rx1,ry1,head,list)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx1, ry1
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: xindex, yindex,i
head = 0
list = 0

do i=1,np
	yindex = ceiling(ry1(i))	
	xindex = ceiling(rx1(i))		
	! linked list algorithm
	if (head(yindex,xindex)==0) then
		head(yindex,xindex)=i
	else 
		list(i) = head(yindex,xindex)
		head(yindex,xindex)=i
	end if	
end do
end subroutine partition

!**********************************************************************************
! Cell partition for dividing particles into grid cells with random grid shift
! OPENMP not implemented as it severely reduces the performance
subroutine partition_rgs(rx1,ry1,head1,list1)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx1, ry1
integer, dimension(:) :: list1
integer, dimension(:,:)    ::   head1
integer :: xindex, yindex, i 
real(kind = dp) :: x_rand, y_rand, xindex_temp, yindex_temp

head1(0:Ly+1,:) = 0		!for negative indics, use indexing fully and not the default one
list1(:) = 0
x_rand = ran()-0.5
y_rand = ran()-0.5
grid_up_down=merge(0,1,y_rand<0)

do i=1,np
	if (xy(2)) then
		yindex_temp = mod(ry1(i)+y_rand, Ly*1.0) 
		yindex_temp = merge(yindex_temp,yindex_temp+Ly,yindex_temp>=0)
		yindex = ceiling(yindex_temp)
	else
		yindex = ceiling(ry1(i) + y_rand) 
	end if

	xindex_temp = mod(rx1(i) + x_rand, Lx*1.0)			! generate random number		
	xindex_temp = merge(xindex_temp,xindex_temp+Lx,xindex_temp>=0)
	xindex = ceiling(xindex_temp)

	! linked list algorithm
	if (head1(yindex,xindex)==0) then
		head1(yindex,xindex)=i
	else 
		list1(i) = head1(yindex,xindex)
		head1(yindex,xindex)=i
	end if

end do
end subroutine partition_rgs

!***********************************************************************************************
!********** MBS Velocity Scaling by sterling and logarithmic approach **************************
!***********************************************************************************************
function vel_scale_ln(Ek,p_count) result(vel_scale)
implicit none

integer :: p_count, f, dummy, fs
real(kind=dp) :: vel_scale, Ek, prob, Ek_hat, Emax, Eprob, p, pmax, E1, gf2, log_temp

dummy = 0
f = 2*(p_count-1)					! Total spatial DOF of all particles
fs = (f/2.0)-1						 
! Sterling approximation
gf2 = fs*log(fs*1.0)-fs + 0.5*log(2*pi*fs)+ log(1.0+1.0/(12*fs) + 1.0/(288*(fs**2)))	! logarithmic of gamma(f/2)

Eprob = 0.5*(f-2)*kbT
log_temp = ((f/2.0)*log(Eprob/kbT))-log(Eprob)-(Eprob/kbT)-gf2
!pmax =  (((f-2)/(2*e))**((f/2.0)-1))/(kbT*gf2)
pmax = exp(log_temp)
Emax = 10*Eprob
do while (dummy==0)
	E1 = ran()*Emax
	log_temp = ((f/2.0)*log(E1/kbT))-log(E1)-(E1/kbT)-gf2
	prob = exp(log_temp)
	p = ran()*pmax

	if (prob > p) then
		Ek_hat = E1
		dummy = 1
	end if
end do
vel_scale = sqrt(Ek_hat/Ek)


end function vel_scale_ln

!*****************************************************************************************************************************
!********************************* Monte Carlo velocity scaling **************************************************************
!*****************************************************************************************************************************
function vel_scale_MC(Ek, p_count) result(vel_scale)
implicit none
integer :: p_count
real(kind=dp) :: a, Ek, psi, rand_temp, inv_psi, S, A_cap, d, pmin, vel_scale 

! physical dimension of the system
d = 2.0d0                                                  

! Generate a random between psi belonging to [1,1+c] 
a = 1.0d0
psi = a + MC_strength*ran()                          ! scale the random number to a value between 1 and 1+c

! Generate random number to choose vel_scale based on the range 0.5 to 1
rand_temp = ran()
inv_psi= 1.0d0/psi
S = merge(psi, inv_psi, rand_temp > 0.5d0) 

! Calculate A
A_cap = ((S)**(d*(p_count-1)))*exp(-1.0d0*(Ek/kbT)*(S**2.0d0 - 1.0d0))

! Choosing to perform the scaling
rand_temp = ran()
pmin = min(1.0d0,A_cap)
vel_scale = merge(S, 1.0d0, rand_temp < pmin)

end function vel_scale_MC
!*****************************************************************************************************************************
! Collision Scheme
subroutine collision(vx, vy, head, list )
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer ,save :: mbs_iter = 0
integer :: i,j,count1, ipar, jpar, k, deficit,virtual(Ly,Lx)
real(kind=dp)	:: r, vx_temp, vy_temp, alpha1,  rel_vel_scale
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), Ek(Ly,Lx), a(2)

jpar = 0
vxcom=0.0
vycom=0.0
Ek = 0.0
virtual = 0
mbs_iter = mbs_iter + 1
!$OMP PARALLEL DO PRIVATE(i,count1,r, alpha1, ipar, deficit, vx_temp,vy_temp,rel_vel_scale) REDUCTION(+:jpar)
do j=1,Lx	
	do i=1,Ly
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

		deficit = 0
		! adding ghost particles for cells overlapping the obstacle
		if (grid_check(i,j)==1 .and. (count1 < Gama) .and. (count1 > 0)) then
			deficit = Gama - count1
			call random_normal(a,2,0.0d0,sqrt(deficit*kbT))	
			vxcom(i,j) = vxcom(i,j) + a(1)
			vycom(i,j) = vycom(i,j) + a(2)
			virtual(i,j) = 1
		end if
		if (count1/=0) then
			vxcom(i,j) = vxcom(i,j)/count1
			vycom(i,j) = vycom(i,j)/count1
		end if	
		jpar = jpar + count1
		! Calculation of the energy
		ipar = head(i,j)
		do while (ipar/=0)
			Ek(i,j) = Ek(i,j) + 0.5*((vx(ipar) -vxcom(i,j))**2 +  (vy(ipar) - vycom(i,j))**2)
			ipar = list(ipar)
		end do
		
		rel_vel_scale = 1.0
		!Obtain the scale factor for MB scaling 
		if (mb_scaling == 1 .and. mod(mbs_iter, mbs_freq) == 0 .and. count1 >= 3 .and. virtual(i,j)==0) then
			rel_vel_scale = vel_scale_ln(Ek(i,j),count1) 
		else if (MC_scaling == 1 .and. count1>=2 .and. virtual(i,j)==0) then
			rel_vel_scale = vel_scale_MC(Ek(i,j),count1)
		endif
		! Performing the collision step with the relevant velocity scale
		! Rotation of velocity
		ipar = head(i,j)
		do while (ipar/=0) 				
			vx_temp  = vx(ipar)-vxcom(i,j) 
			vy_temp  = vy(ipar)-vycom(i,j) 			
			vx(ipar) = vxcom(i,j) + rel_vel_scale*(1.0*cos(alpha1)*vx_temp + sin(alpha1)*vy_temp)
			vy(ipar) = vycom(i,j) + rel_vel_scale*(-1.0*sin(alpha1)*vx_temp + cos(alpha1)*vy_temp)
			ipar = list(ipar)			
		end do
	end do
end do	
!$OMP END PARALLEL DO

IF (jpar /= np) write(*,*) "Particles lost during collision"
end subroutine collision

!***********************************************************************************************
!********************************* Collision for Random Grid shift *****************************
!***********************************************************************************************
subroutine collision_rgs(vx, vy, head1, list1)
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list1
integer, dimension(:,:)    ::   head1
integer, save :: mbs_iter=0
integer :: i,j,count1, ipar, jpar, k, deficit, virtual(0:Ly+1,Lx)
real(kind=dp)	:: r, vx_temp, vy_temp, alpha1,  rel_vel_scale
real(kind=dp) :: vxcom(0:Ly+1,Lx), vycom(0:Ly+1,Lx), Ek(0:Ly+1,Lx), a(2)

jpar = 0
vxcom(0:Ly+1,:)=0.0
vycom(0:Ly+1,:)=0.0
Ek(0:Ly+1,:) = 0.0
virtual(0:Ly+1,:) = 0
mbs_iter = mbs_iter + 1

!$OMP PARALLEL DO PRIVATE(i,j,count1,r, alpha1, ipar, deficit, a, vx_temp,vy_temp,rel_vel_scale)REDUCTION(+:jpar)
do j=1,Lx
	do i=0,Ly+1
		count1 = 0
				
		r = ran()

		if (r >=0.5) then
			alpha1 = -alpha
		else 
			alpha1 = alpha
		end if			
		
		ipar = head1(i,j)					
		do while (ipar/=0)			
			vxcom(i,j) = vxcom(i,j) + vx(ipar)
			vycom(i,j) = vycom(i,j) + vy(ipar)
			count1 = count1 + 1
			ipar = list1(ipar)
		end do

		deficit = 0
		! adding ghost particles for wall
		if (.NOT.(xy(2)) .and. (i==0+grid_up_down .or. i== Ly+grid_up_down .or. grid_check(i,j)==1) .and. (count1 < Gama) .and. (count1 > 0)) then
			deficit = Gama - count1
			call random_normal(a,2,0.0d0,sqrt(deficit*kbT))	
			vxcom(i,j) = vxcom(i,j) + a(1)
			vycom(i,j) = vycom(i,j) + a(2)
			virtual(i,j) = 1
		end if
		
		if (count1/=0) then
			vxcom(i,j) = vxcom(i,j)/(count1 + deficit)
			vycom(i,j) = vycom(i,j)/(count1 + deficit)
		end if
		jpar = jpar + count1
		! Calculation of the energy
		ipar = head1(i,j)
		do while (ipar/=0)
			Ek(i,j) = Ek(i,j) + 0.5*((vx(ipar) -vxcom(i,j))**2 +  (vy(ipar) - vycom(i,j))**2)
			ipar = list1(ipar)
		end do
		
		rel_vel_scale = 1.0
		!Obtain the scale factor for MB scaling 
		if (mb_scaling == 1 .and. mod(mbs_iter, mbs_freq) == 0 .and. count1 >= 3 .and. virtual(i,j)==0 ) then
			rel_vel_scale = vel_scale_ln(Ek(i,j),count1) 
		else if (MC_scaling == 1 .and. count1>=2 .and. virtual(i,j)==0) then
			rel_vel_scale = vel_scale_MC(Ek(i,j),count1)
		endif
		
		! Performing the collision step with the relevant velocity scale
		! Rotation of velocity
		ipar = head1(i,j)
		do while (ipar/=0) 				
			vx_temp  = vx(ipar)-vxcom(i,j) 
			vy_temp  = vy(ipar)-vycom(i,j) 			
			vx(ipar) = vxcom(i,j) + rel_vel_scale*(1.0*cos(alpha1)*vx_temp + sin(alpha1)*vy_temp)
			vy(ipar) = vycom(i,j) + rel_vel_scale*(-1.0*sin(alpha1)*vx_temp + cos(alpha1)*vy_temp)
			ipar = list1(ipar)			
		end do
	end do
end do	
!$OMP END PARALLEL DO
IF (jpar /= np) write(*,*) "Particles lost during collision", jpar, np
end subroutine collision_rgs

!****************************************************************
! Subroutine for thermal boundary conditions for a cylinder
subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy)
implicit none
integer, parameter :: np1=1000		! Expected no. of particles crossing boundary
real(kind=dp), dimension(:) :: rx, ry, rx1, ry1, vx, vy
real(kind=dp) ::  t_app, t_dec
integer :: i, j
!logical :: Wall(np)
real(kind=dp) :: vxwall(1), vywall(1)

j = 0
!Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)
!call random_normal(vxwall,np1,avg,std)
!call random_weibull(vywall,np1)

!$OMP PARALLEL IF(np>100000) 
!$OMP DO PRIVATE(t_app,t_dec,vxwall,vywall) SCHEDULE(guided)
do i = 1,np
	if (ry1(i) > (Ly*1.0) .or. ry1(i) < 0.0) then	
		if ((ry1(i)-ry(i))<0.0) then
			t_app = ry(i)/abs(vy(i))
			t_dec = dt_c - t_app
                        call random_weibull(vywall,1)
			vy(i) = vywall(1)
			ry1(i) = 0.0 + vy(i)*t_dec	
		else 
			t_app = (Ly-ry(i))/abs(vy(i))
			t_dec = dt_c - t_app
                        call random_weibull(vywall,1)
			vy(i) = -vywall(1)
			ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i)  = vx(i)-f_b*dt_c
		IF (verlet == 2) THEN
			vx(i)  = vx(i) + f_b*(dt_c+t_app)/2
			rx1(i) = rx(i) + vx(i)*t_app
                        call random_normal(vxwall,1,avg,std)
			vx(i)  = vxwall(1) 
			vx(i)  = vx(i) + f_b*t_dec/2
			rx1(i) = rx1(i) + vx(i)*t_dec
		ELSE IF (verlet == 1) THEN
			rx1(i) = rx(i) + vx(i)*t_app + (f_b*t_app**2)/2 
                        call random_normal(vxwall,1,avg,std)
			vx(i)  = vxwall(1)
			rx1(i) = rx1(i) + vx(i)*t_dec + (f_b*t_dec**2)/2 
			vx(i)  = vx(i) + f_b*t_dec
		ELSE
			rx1(i) = rx(i) + vx(i)*t_app 
                        call random_normal(vxwall,1,avg,std)
			vx(i)  = vxwall(1)
			rx1(i) = rx1(i) + vx(i)*t_dec 
			vx(i)  = vx(i) + f_b*dt_c
		END IF
	end if
!	if (j==np1) then		! if more particles than expected are crossing, generate new numbers
!		call random_normal(vxwall,np1,avg,std)
!		call random_weibull(vywall,np1)
!		j=0	
!	end if
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine thermal_wall

!***********************************************************************
!***********************************************************************
! Writing data to files for analysis
subroutine v_avg(vx,vy, rx, ry, head,list)
implicit none
INTEGER,PARAMETER :: grid_points = half_plane*Ly + 1
real(kind=dp), dimension(:) :: vx, vy, rx, ry
integer   :: head(:,:), list(:)
integer, save ::  t_count = 0, file_count = 0
integer :: i, j, k, ipar, p_count, out_unit, y_low, y_up, density(Ly,Lx)
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), temp_com(Ly,Lx),  weight, plane_pos
real(kind=dp), save :: vx1(Ly,Lx)=0.0d0, vy1(Ly,Lx)=0.0d0, forcing(2) = 0.0d0, tx(Ly,Lx) = 0.0d0, dens(Ly,Lx) = 0.0,&
vx_avg(grid_points)=0.0, vy_avg(grid_points)=0.0,  tot_part_count(grid_points)=0.0	! every grid points considered from 0 to Ly
logical :: Lexist
CHARACTER(LEN=200)  :: fname

vxcom   = 0.0
vycom   = 0.0
out_unit = 20
temp_com = 0.0
density   = 0.0

t_count = t_count+1
!$OMP PARALLEL DO PRIVATE(i,p_count,ipar,plane_pos,weight,k)
do j = 1,Lx    			  			
do i = 1,Ly
	p_count = 0
	ipar = head(i,j)					
	do while (ipar/=0)			
		vxcom(i,j) = vxcom(i,j) + vx(ipar)
		vycom(i,j) = vycom(i,j) + vy(ipar)
		p_count = p_count + 1
		IF (half_plane==2) THEN
			plane_pos = i-0.5
			weight  = (ry(ipar)-plane_pos)/0.5
			If (weight < 0.0) weight = 1 + weight
			k = merge(i*half_plane, i*half_plane-1, ry(ipar) .gt. plane_pos)
		ELSE
			plane_pos = i-1
			weight = ry(ipar)-plane_pos
			k = i
		END IF
		vx_avg(k)   = vx_avg(k)   + vx(ipar)*(1-weight)
		vx_avg(k+1) = vx_avg(k+1) + vx(ipar)*weight
		vy_avg(k)   = vy_avg(k)   + vy(ipar)*(1-weight)
		vy_avg(k+1) = vy_avg(k+1) + vy(ipar)*weight
		tot_part_count(k) = tot_part_count(k) + (1-weight)
		tot_part_count(k+1) = tot_part_count(k+1) + weight
		ipar = list(ipar)
	end do
	density(i,j) = p_count
	if(temperature .and. p_count > 1) then
		ipar = head(i,j)    			  		
		do while (ipar/=0)
			temp_com(i,j) = temp_com(i,j) + 0.5*((vx(ipar) - vxcom(i,j))**2 + (vy(ipar) - vycom(i,j))**2) 					
			ipar = list(ipar)		
		end do
		temp_com(i,j) = temp_com(i,j)/(p_count - 1)
	end if
end do
end do 
!$OMP END PARALLEL DO
	
! write(*,*) temp_com(1,1)
! 	write (fname, "(A<LEN(trim(file_name))>,A10)") trim(file_name),"_energy"                            
! 	inquire(file = trim(data_path)//trim(fname)//'.dat', exist = Lexist)  	
! 	if (Lexist) then	
! 		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',status="old",action="write",position="append")
!     	else
! 		open(unit = out_unit, file=trim(data_path)//trim(fname)//'.dat', status = "new", action = "write")
! 	end if	
!     	write (out_unit,"(<Ly>F10.5)") temp_com(:,1)
!     	close(out_unit)

vx1 = vx1 + vxcom
vy1 = vy1 + vycom
dens = dens + density    			  	
forcing = forcing + force
if (temperature) tx = tx + temp_com

if (modulo(t_count,ensemble_num)==0) then
	file_count = file_count+1
    	
	! Averaging all the velocities over the respective particles and time
	vx_avg = vx_avg/tot_part_count
	vy_avg = vy_avg/tot_part_count
	
	vx1 = vx1/dens
	vy1 = vy1/dens
	forcing = forcing/t_count
	tx = tx/t_count
	dens = dens/t_count

	write (fname, "(A<LEN(trim(file_name))>,A10)") trim(file_name),"_drag_lift"   !this format of writing is valid only in ifort, not gfortran                  
	inquire(file = trim(data_path)//trim(fname)//'.dat', exist = Lexist)  	
	if (Lexist) then	
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',status="old",action="write",position="append")
    	else
		open(unit = out_unit, file=trim(data_path)//trim(fname)//'.dat', status = "new", action = "write")
	end if	
    	write (out_unit,*) forcing		  		
    	close(out_unit)

	write (fname, "(A<LEN(trim(file_name))>,A7,I0)") trim(file_name),"_vx_vy_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
	
	DO i=1,half_plane*Ly+1
		write(out_unit,*) (i*1.0-1)/(half_plane*1.0),vx_avg(i),vy_avg(i)
	END DO
!	write (out_unit,"(<half_plane*Ly+1>F10.5)") vx_avg	!could also directly use the size(vx_avg) in format string
!	write (out_unit,"(<half_plane*Ly+1>F10.5)") vy_avg
!	do i = 1,Ly
!		write (out_unit,"(<Lx>F10.5)") vx1(i,:)
!	end do
!	do i = 1,Ly
!		write (out_unit,"(<Lx>F10.5)") vy1(i,:)
!	end do
	close(out_unit)	
	
	write (fname, "(A<LEN(trim(file_name))>,A5,I0)") trim(file_name),"_rho_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
    	
	do i = 1,Ly
		write (out_unit,"(<Lx>F10.5)") dens(i,:)
	end do
	close(out_unit)	

	if (temperature) then
		write (fname, "(A<LEN(trim(file_name))>,A6,I0)") trim(file_name),"_temp_",file_count                             
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
		do i = 1,Ly
			write (out_unit,"(<Lx>F10.5)") tx(i,:)
		end do
		close(out_unit)
	end if
	vx1 = 0.0
	vy1 = 0.0
	vx_avg = 0.0
	vy_avg = 0.0
	tot_part_count=0.0
	tx  = 0.0
	dens= 0.0
	t_count = 0	
	forcing = 0.0
end if    		

end subroutine v_avg
!***********************************************************************
! Subroutine for bounce back boundary condition at wall
subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy)
implicit none
real(kind=dp), dimension(:)  :: rx, ry, rx1, ry1, vx, vy
real(kind=dp) :: t_app, t_dec
integer :: i
!logical :: Wall(np)

!Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)  			  		

!$OMP PARALLEL IF(np>100000)
!$OMP DO PRIVATE(t_app,t_dec) SCHEDULE(guided)
do i = 1,np
	if (ry1(i) > (Ly*1.0) .or. ry1(i) < 0.0)   then			
		if ((ry1(i)-ry(i))<0.0) then			
			t_app  = ry(i)/abs(vy(i))
			t_dec  = dt_c - t_app
			vy(i)  = -vy(i)
			ry1(i) = 0.0 + vy(i)*t_dec	
		else 			
			t_app  = (Ly-ry(i))/abs(vy(i))
			t_dec  = dt_c - t_app
			vy(i)  = -vy(i)
			ry1(i) = Ly + vy(i)*t_dec	
		end if		
		vx(i) = vx(i)-f_b*dt_c
		IF (verlet == 2) THEN
			vx(i)  = vx(i) + f_b*(dt_c+t_app)/2
			rx1(i) = rx(i) + vx(i)*t_app
			vx(i)  = -vx(i) 
			vx(i)  = vx(i) + f_b*t_dec/2
			rx1(i) = rx1(i) + vx(i)*t_dec
		ELSE IF (verlet == 1) THEN
			rx1(i) = rx(i) + vx(i)*t_app + (f_b*t_app**2)/2
			vx(i)  = vx(i) + f_b*t_app
			vx(i)  = -vx(i)
			rx1(i) = rx1(i) + vx(i)*t_dec + (f_b*t_dec**2)/2 
			vx(i)  = vx(i) + f_b*t_dec
		ELSE
			rx1(i) = rx(i) + vx(i)*t_app 
			vx(i)  = -vx(i)
			rx1(i) = rx1(i) + vx(i)*t_dec 
			vx(i)  = vx(i) + f_b*dt_c
		END IF
	end if	
end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine bounce_wall

!**********************************************************************************
! for writing file which gives the details of all parameters used in the simulation
SUBROUTINE param_file()
implicit none

OPEN(UNIT = 10, FILE="Parameters.txt",ACTION = "write",STATUS="replace")

write(10,*)"Parameters used in the current simulation to which the data belongs"
write(10,*) "Particle density: ",Gama
write(10,*) "alpha rotation: ",alpha
write(10,*) "kbT: ",kbT
write(10,*) "Domain size Ly,Lx:",[Ly,Lx]
write(10,*)
write(10,*) "Time step used: ", dt_c
write(10,*) "Total time: ", tmax
write(10,*) "Total iterations: ",tmax/dt_c
write(10,*) "Force applied:",f_b
write(10,*)
write(10,*) "Streaming Algorithm:  ",merge(merge('Leapfrog','  Verlet',verlet/=1),'   Euler',verlet/=0)
write(10,*) "Periodicity in x: ",merge('Yes',' No',xy(1))
write(10,*) "Periodicity in y: ",merge('Yes',' No',xy(2))
write(10,*) "Random Grid Shift applied: ",merge('Yes',' No',random_grid_shift==1)
write(10,*) "Wall Boundary Condition: ",merge('Thermal wall',' Bounce wall',wall_thermal)
write(10,*)
write(10,*) "MC applied: ",merge('Yes',' No',MC_scaling==1)
IF (MC_scaling ==1) write(10,*) "Strength of MC: ",MC_strength
write(10,*) "MBS applied: ",merge('Yes',' No',mb_scaling==1)
IF (mb_scaling == 1) write(10,*) "MBS applied in number of iterations: ",mbs_freq
write(10,*) "Averaging starts from iteration: ",t_avg
write(10,*) "Interval between samples considered for average: ",avg_interval
write(10,*) "Number of samples considered for average : ",ensemble_num
write(10,*) "Mid plane of every cell used for averaging: ", merge('Yes',' No',half_plane == 2)
write(10,*)
write(10,*) "Analytical Viscosity as per given condition:  ",mu_tot 
write(10,*) "Maximum Poiseuille velocity (analytical): ", (Gama * Ly**2.0 *f_b)/(8.0*mu_tot)
write(10,*)
IF (obst == 1) write(10,*) "Cylinder centre position: ",xp,yp," and its radius: ",rad
close(10)

END SUBROUTINE param_file

end module SRD_library
