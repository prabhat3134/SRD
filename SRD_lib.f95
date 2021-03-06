module SRD_library

implicit none


INTEGER, PARAMETER :: dp = selected_real_kind(15, 307), long = selected_int_kind(range(1)*2)
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828
REAL, PARAMETER ::  kbT = 1.0, dt_c = 1.0, alpha = pi/2
INTEGER, PARAMETER :: Ly = 50, Lx = 50, Gama = 10, m=1, a0=1, ensemble_num = 50000, half_plane =2, freq = 100
REAL, PARAMETER :: rad = 10, xp = Lx/4.0, yp = Ly/2.0				! cylinder parameters
INTEGER, PARAMETER :: random_grid_shift = 1, mb_scaling = 1, obst = 0, verlet = 1, MPC_AT = 1
INTEGER :: grid_check(Ly+2,Lx)=0
LOGICAL, PARAMETER :: xy(2)=[.TRUE., .FALSE.], temperature = .TRUE., wall_thermal = .TRUE.
REAL(kind=dp)  :: force(2), mu_tot
CHARACTER(len=100) :: file_name='Poise_verlet_f50_g1e5', data_path='./'	!file_name of size 20
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
logical :: l1(size(x_dummy)), l2(size(x_dummy))

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
subroutine random_weibull(vel,np)
implicit none
real(kind=dp) :: vel(np)
integer :: np,i
real :: a,b,v,Pv, vmax, vprob, Pmax, P
i=1
a = sqrt(2*kbT)
b = 2.0
vprob = sqrt(kbT)		! most probable speed
vmax = 5*vprob			! maximum speed possible
Pmax = 1/sqrt(kbT*e)		!maximum probability

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
subroutine initialize(x_dummy, y_dummy, rx,ry,vx,vy, np, av, std, fb)
implicit none
INTEGER :: np, i
REAL :: av, std, block(3), vp_max, fb 
REAL(kind = dp) :: x_dummy(:), y_dummy(:), mu_kin, mu_col
REAL(kind=dp) :: rx(np), ry(np), vx(np), vy(np)

mu_kin = (gama*kbT*dt_c)*(gama/((gama - 1.0 + exp(-1.0*gama))*(1.0-cos(2.0*alpha)))- 0.5)
mu_col = ((1.0-cos(alpha))/(12.0*dt_c))*(gama - 1.0 + exp(-1.0*gama))
mu_tot = mu_kin + mu_col
vp_max = (gama * Ly**2.0 *fb)/(8.0*mu_tot)
block = [xp,yp,rad]
do i=1, np
	x_dummy(i) = ran()*Lx
	y_dummy(i) = ran()*Ly
end do

IF (obst == 1) call box_eliminate(x_dummy, y_dummy)
rx = x_dummy
ry = y_dummy

call random_normal(vx,np,av,std)
call random_normal(vy,np,av,std)
IF (obst ==1) call gridcheck(grid_check, block)

! impose a Poiseuille profile in the start
vx = vx + (0.6*4.0*vp_max*(Ly - ry)*ry)/(Ly**2.0)

end subroutine initialize

!*********************************************
! Streaming step in SRD
subroutine streaming(rx, ry, rx1, ry1, vx, vy, np, l1, g)
implicit none
integer :: np,i
real(kind = dp),dimension(:) :: rx, ry, rx1, ry1, vx, vy
real :: g
logical :: l1(:)

IF (verlet == 2) THEN	 	! LEAPFROG algorithm
	vx  = vx + g*dt_c
	rx1 = rx + vx*dt_c
ELSE IF (verlet ==1) THEN	!Verlet Algorithm
	rx1 = rx + vx*dt_c + (g*dt_c**2)/2.0 
	vx  = vx + g*dt_c
ELSE				!EULER algorithm
	rx1 = rx + vx*dt_c 
	vx  = vx + g*dt_c
END IF
ry1 = ry + vy*dt_c

IF (obst==1) l1  = (rx1 - xp)**2 + (ry1 - yp)**2 < rad**2

!l1 = (rx1 >= xp - 2*rad) && (rx1 <= xp + 2*rad) && (ry1 >= yp - 2*rad) && (ry1 <= yp + 2*rad)

end subroutine streaming
!*************************************************
! Dind the intersection angles, and, the time required to intersect, for particles which cross the cylinder boundary
! Enters the module every time step
subroutine par_in_cyl(l1, rx, ry, rx1, ry1,vx,vy, theta_intersect, g)

implicit none
integer :: i
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1, vx, vy, theta_intersect
logical, dimension(:) :: l1
real(kind =dp) :: dx, dy, dr, de, disc, x_sol1, x_sol2, y_sol1, y_sol2, th, x_sol, y_sol 
real(kind = dp) :: exprx1, exprx2, expry1, expry2 
real :: g

theta_intersect = 0.0
y_sol = 0.0
x_sol = 0.0
! finding the solution of y - mx -c = 0 and x^2 + y^2 = r^2 using the quadratic formula 
! Solving the intersection angle for every particle inside the cylinder (using the list l1(i)) 
do i=1,size(rx1)
	if(l1(i)) then
! Calculating the slope, of the line joining the initial and final particle positions. 
! for reference, check: 	http://mathworld.wolfram.com/Circle-LineIntersection.html	
! The line is shifted to origin. 
		dx = rx1(i) - rx(i)
		dy = ry1(i) - ry(i)
		dr = SQRT(dx**2 + dy**2)
		de = (rx(i) -xp)*(ry1(i)-yp) - (rx1(i)-xp)*(ry(i) -yp)
		disc = (rad**2)*(dr**2)-(de**2) 
		
		if (disc > 0) then 
			x_sol1 = (de*dy + SIGN(1.0d0,dy)*dx*SQRT(disc))/(dr**2)
			x_sol2 = (de*dy - SIGN(1.0d0,dy)*dx*SQRT(disc))/(dr**2)
			y_sol1 = (-de*dx + ABS(dy)*SQRT(disc))/(dr**2)
			y_sol2 = (-de*dx - ABS(dy)*SQRT(disc))/(dr**2)
		else
			write(*,*) "Negative Discriminant: particle doesn't collide",rx(i), ry(i), rx1(i), ry1(i)
		endif

		exprx1 = (x_sol1 + xp - rx1(i))*(x_sol1 + xp - rx(i)) 
		exprx2 = (x_sol2 + xp - rx1(i))*(x_sol2 + xp - rx(i))
		expry1 = (y_sol1 + yp - ry1(i))*(y_sol1 + yp - ry(i)) 
		expry2 = (y_sol2 + yp - ry1(i))*(y_sol2 + yp - ry(i)) 
		
		if ( exprx1 < 0.0 .and. expry1 < 0.0 ) then
			x_sol = x_sol1
			y_sol = y_sol1
		elseif ( exprx2 < 0.0 .and. expry2 < 0.0 ) then 
			x_sol = x_sol2
			y_sol = y_sol2
		else
			write(*,*) "intersection algorithm fails",rx(i), ry(i), rx1(i), ry1(i)
		endif 

		th = atan2(y_sol, x_sol)
		theta_intersect(i) = merge(th, th + 2*pi, th >= 0)
	
	endif
end do

call pos_thermal_update(l1,theta_intersect,rx,ry,rx1,ry1,vx,vy, g)

end subroutine par_in_cyl
!*********************************************
! based on the angle theta_intersect of the particles inside the cylinder, place them back to the right position
subroutine pos_bounce_update(l1, theta_intersect,rx, ry, rx1, ry1, vx, vy, g)

implicit none
integer :: i, j
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1, theta_intersect, vx, vy
real(kind=dp) :: t_app, t_sep, jacobian(2,2), th, vr(2), vxy(2), force_r(2)
real :: g
logical, dimension(:) :: l1
force = 0.0
do i=1, size(rx1)
	if(l1(i)) then 
		th = theta_intersect(i)
		jacobian = reshape([cos(th),- sin(th), sin(th), cos(th)],shape(jacobian))
		t_app = abs(yp+rad*sin(th) - ry(i))/abs(vy(i))
		t_sep = dt_c - t_app
	
		if (t_sep < 0)  write(*,*) 'Cylinder collision: negative time'		
		! x velocity changes in approach -> reversal -> x velocity changes in separation	
		vx(i) = vx(i) - g*dt_c
		
		! velocity of approach at the surface of the cylinder
		! vx(i) = vx(i) + g*t_app
		
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
		vx(i)  = vx(i) + g*dt_c

	endif
end do
end subroutine pos_bounce_update

!*********************************************
! based on the angle theta_intersect of the particles inside the cylinder, place them back to the right position
subroutine pos_thermal_update(l1, theta_intersect,rx, ry, rx1, ry1, vx, vy, g)

implicit none
integer, parameter :: np1=500
integer :: i, j
real(kind=dp), dimension(:) :: rx,ry, rx1, ry1, theta_intersect, vx, vy
real(kind=dp) :: t_app, t_sep, jacobian(2,2), th, v_polar_app(2), v_polar_sep(2), vxy(2), force_polar(2), vr_wall(np1),vt_wall(np1)
real :: g, av, std
logical, dimension(:) :: l1
av = 0.0
force = 0.0
std = sqrt(kbT)
j = 0
call random_normal(vt_wall,np1,av,std)
call random_weibull(vr_wall,np1)

do i=1, size(rx1)
	if(l1(i)) then 
		j = j+1
		th = theta_intersect(i)
		jacobian = reshape([cos(th),- sin(th), sin(th), cos(th)],shape(jacobian))
		t_app = abs(yp+rad*sin(th) - ry(i))/abs(vy(i))
		t_sep = dt_c - t_app
	
		if (t_sep < 0)  write(*,*) 'Cylinder collision: negative time'		
		! x velocity changes in approach -> reversal -> x velocity changes in separation	
		vx(i) = vx(i) - g*dt_c
		! velocity of approach at the surface of the cylinder
		! vx(i) = vx(i) + g*t_app
		
		! coordinate transformation for reversal in dirn of normal and tangential velocity and propogation
		vxy = [vx(i), vy(i)]
		v_polar_app  = MATMUL(jacobian,vxy)
		v_polar_sep  = [vr_wall(j), vt_wall(j)] 
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
		vx(i)  = vx(i) + g*dt_c

		if (j==np1) then		! if more particles than expected are crossing, generate new numbers
			call random_normal(vt_wall,np1,av,std)
			call random_weibull(vr_wall,np1)
			j=0	
		end if
	endif
end do

end subroutine pos_thermal_update
!*************************************************
! implementing the periodic boundary condition
! Periodic boundary condition in x and y as per the condition given by xy
subroutine periodic_xy( rx1, ry1, np)
implicit none
integer :: np
logical :: l1(np)
real(kind = dp) :: rx1(np), ry1(np)

if (xy(1)) then
	l1 = (rx1 <= 0.0)
	rx1 = mod(rx1, Lx*1.0)
	where (l1)
		rx1 = rx1 + Lx
	end where
end if

if (xy(2)) then
	l1 = (ry1 <= 0.0)
	ry1 = mod(ry1, Ly*1.0)
	where (l1)
		ry1 = ry1 + Ly
	end where
end if
end subroutine periodic_xy
!********************************************************
! Grid check, gives which grids are overlapping with obstacles for further use in collision
subroutine gridcheck(grid_check, block)
implicit none
INTEGER :: grid_check(Ly+2,Lx), i, j, x, y
REAL :: block(:), xc, yc, l, b, R, rtmp

grid_check = 0
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
subroutine partition(rx1,ry1,head,list,Lx,Ly,np)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx1, ry1
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer :: np, Lx, Ly, xindex, yindex,i
head = 0
list = 0

do i=1,np
	yindex = ceiling(ry1(i))	
	xindex = ceiling(rx1(i))		
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

!**********************************************************************************
! Cell partition for dividing particles into grid cells with random grid shift
subroutine partition_rgs(rx1,ry1,head1,list1,Lx,Ly,np)    	! assumed no grid shifting
implicit none
real(kind=dp), dimension(:) :: rx1, ry1
integer, dimension(:) :: list1
integer, dimension(:,:)    ::   head1
integer :: np, Lx, Ly, xindex, yindex,i
real(kind = dp) :: x_rand, y_rand, xindex_temp, yindex_temp

head1 = 0		! Head starts from [1, Ly+1]
list1 = 0
x_rand = ran()	! random shift from [0,1]
y_rand = ran()
do i=1,np
	if (xy(2)) then
		yindex_temp = mod(ry1(i)+y_rand, Ly*1.0) 
		yindex = ceiling(yindex_temp) 
	else
		yindex = ceiling(ry1(i) + y_rand) 		
	end if

	xindex_temp = mod(rx1(i) + x_rand, Lx*1.0)						
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
function vel_scale_ln(Ek,np) result(vel_scale)
implicit none

integer :: np, f, dummy, fs
real(kind=dp) :: vel_scale, Ek, prob, Ek_hat, Emax, Eprob, p, pmax, E1, gf2, log_temp

dummy = 0
f = 2*(np-1)					! Total spatial DOF of all particles
fs = (f/2.0)-1						 
! Sterling approximation
gf2 = fs*log(fs*1.0)-fs + 0.5*log(2*pi*fs)+ log(1.0+1.0/(12*fs) + 1.0/(288*(fs**2)))	! logarithmic of gamma(f/2)

Eprob = 0.5*(f-2)*kbT
log_temp = ((f/2.0)*log(Eprob/kbT))-log(Eprob)-(Eprob/kbT)-gf2
!pmax =  (((f-2)/(2*e))**((f/2.0)-1))/(kbT*gf2)
pmax = exp(log_temp)
Emax = 10*Eprob
!write(*,*) ek,np,f,gf2,eprob,pmax,emax
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
!************************************************************************************************************
! Collision Scheme
subroutine collision(vx, vy, temp, temp_com, tempy, head, list, Lx, Ly )
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list
integer, dimension(:,:)    ::   head
integer, save :: mbs_iter = 0
integer :: i,j,count1, ipar, jpar, Lx, Ly, out_unit, k, part_iter
real(kind=dp)	:: r, vx_temp, vy_temp, alpha1,  rel_vel_scale, vx_rel, vy_rel, AT_vx(Gama*Ly*Lx), AT_vy(Gama*Ly*Lx), AT_mean(2)
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), temp(Ly,Lx), temp_com(Ly,Lx), tempy(Ly,Lx), Ek(Ly,Lx)

IF (MPC_AT ==1) THEN
	call random_normal(AT_vx,Gama*Ly*Lx,0.0,sqrt(kbT))
	call random_normal(AT_vy,Gama*Ly*Lx,0.0,sqrt(kbT))
END IF

part_iter = 0
out_unit=20
vxcom=0.0
vycom=0.0
temp_com=0.0
tempy=0.0
temp =0.0
Ek = 0.0
mbs_iter = mbs_iter + 1
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
		
		temp_com(i,j) =  (0.5*(vxcom(i,j)**2 +  vycom(i,j)**2))
		! Rotation of velocity
		IF (MPC_AT == 1) THEN
			AT_mean(1) = SUM(AT_vx(part_iter+1:part_iter + count1))/count1
			AT_mean(2) = SUM(AT_vy(part_iter+1:part_iter + count1))/count1
			ipar = head(i,j)
			do while (ipar/=0) 	
				part_iter = part_iter + 1
				vx(ipar) = vxcom(i,j) + AT_vx(part_iter) - AT_mean(1) 
				vy(ipar) = vycom(i,j) + AT_vy(part_iter) - AT_mean(2)		
				ipar = list(ipar)			
			end do	
		ELSE
			ipar = head(i,j)
			do while (ipar/=0) 	
				temp(i,j) = temp(i,j) + (0.5*(vx(ipar)**2 -vxcom(i,j)**2 +  vy(ipar)**2 - vycom(i,j)**2))**2 			
				vx_temp  = vx(ipar)-vxcom(i,j) 
				vy_temp  = vy(ipar)-vycom(i,j) 			
				vx(ipar) = vxcom(i,j) + cos(alpha1)*vx_temp + sin(alpha1)*vy_temp
				vy(ipar) = vycom(i,j) - sin(alpha1)*vx_temp + cos(alpha1)*vy_temp		
				Ek(i,j) = Ek(i,j) + 0.5*((vx(ipar) -vxcom(i,j))**2 +  (vy(ipar) - vycom(i,j))**2)  			
				ipar = list(ipar)			
				part_iter = part_iter + 1
			end do	
			temp(i,j) = sqrt(temp(i,j)/count1)
			
			! obtaining the velocity scale factor
			if (mb_scaling == 1 .and. mod(mbs_iter,freq) == 0 .and. count1 >= 3 ) then 
				! count1 >=3: This condition satisfies for the minimum number of particles required for MB scaling
				rel_vel_scale = vel_scale_ln(Ek(i,j),count1) 					! Implementing MB scaling 
					
				! scaling the relative velocities
				ipar = head(i,j)			
				
				do while (ipar/=0)  
					vx_rel = vx(ipar)-vxcom(i,j)
					vy_rel = vy(ipar)-vycom(i,j)						! Relative velocities after collision
					vx(ipar)  = vxcom(i,j) + rel_vel_scale*(vx_rel) 			
					vy(ipar)  = vycom(i,j) + rel_vel_scale*(vy_rel) 			! Add scaled relative velocities 		
					ipar = list(ipar)							! after collision using MBS scheme
				end do	
			end if
		END IF
	end do
	!write(out_unit,*), vxcom(i,:)
end do	

IF (part_iter .lt. Gama*Ly*Lx) write(*,*) "Particles missing from domain"
end subroutine collision

!***********************************************************************************************
!********************************* Collision for Random Grid shift *****************************
!***********************************************************************************************
subroutine collision_rgs(vx, vy, head1, list1, Lx, Ly, np)
implicit none
real(kind=dp), dimension(:) :: vx, vy
integer, dimension(:) :: list1
integer, dimension(:,:)    ::   head1
integer, save :: mbs_iter=0
integer :: i,j,count1, ipar, jpar, Lx, Ly, out_unit, k, deficit, virtual(Ly+1,Lx), part_iter
real(kind=dp)	:: r, vx_temp, vy_temp, alpha1,  rel_vel_scale, vx_rel, vy_rel, AT_vx(Gama*Ly*Lx), AT_vy(Gama*Ly*Lx), AT_mean(2)
real(kind=dp) :: vxcom(Ly+1,Lx), vycom(Ly+1,Lx), Ek(Ly+1,Lx), a(2)
integer ::  np

IF (MPC_AT ==1) THEN
	call random_normal(AT_vx,Gama*Ly*Lx,0.0,sqrt(kbT))
	call random_normal(AT_vy,Gama*Ly*Lx,0.0,sqrt(kbT))
END IF
part_iter = 0
out_unit=20
vxcom=0.0
vycom=0.0
Ek = 0.0
virtual = 0
mbs_iter = mbs_iter + 1
do i=1,Ly+1	
	do j=1,Lx
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
		if ((i==1 .or. i== Ly+1 .or. grid_check(i,j)==1) .and. (count1 < Gama) .and. (count1 > 0)) then
			deficit = Gama - count1
			call random_normal(a,2,0.0,sqrt(deficit*kbT))	
			vxcom(i,j) = vxcom(i,j) + a(1)
			vycom(i,j) = vycom(i,j) + a(2)
			virtual(i,j) = 1
		end if
		
		if (count1/=0) then
			vxcom(i,j) = vxcom(i,j)/(count1 + deficit)
			vycom(i,j) = vycom(i,j)/(count1 + deficit)
		end if	
		
		! Rotation of velocity
		IF (MPC_AT == 1) THEN
			AT_mean(1) = SUM(AT_vx(part_iter+1:part_iter + count1))/count1
			AT_mean(2) = SUM(AT_vy(part_iter+1:part_iter + count1))/count1
			ipar = head1(i,j)
			do while (ipar/=0) 	
				part_iter = part_iter + 1
				vx(ipar) = vxcom(i,j) + AT_vx(part_iter) - AT_mean(1) 
				vy(ipar) = vycom(i,j) + AT_vy(part_iter) - AT_mean(2)		
				ipar = list1(ipar)			
			end do	
		ELSE
			ipar = head1(i,j)
			do while (ipar/=0) 				
				vx_temp  = vx(ipar)-vxcom(i,j) 
				vy_temp  = vy(ipar)-vycom(i,j) 			
				vx(ipar) = vxcom(i,j) + cos(alpha1)*vx_temp + sin(alpha1)*vy_temp
				vy(ipar) = vycom(i,j) - sin(alpha1)*vx_temp + cos(alpha1)*vy_temp
				Ek(i,j) = Ek(i,j) + 0.5*((vx(ipar) -vxcom(i,j))**2 +  (vy(ipar) - vycom(i,j))**2) 					
				ipar = list1(ipar)			
				part_iter = part_iter + 1
			end do	
			!temp(i,j) = sqrt(temp(i,j)/count1)
			
			! obtaining the velocity scale factor
			if (mb_scaling == 1 .and. mod(mbs_iter,freq) == 0) then 
				! count1 >=3: This condition satisfies for the minimum number of particles required for MB scaling
				if (count1 >= 3 .and. virtual(i,j)==0 ) then
					
					rel_vel_scale = vel_scale_ln(Ek(i,j),count1) 					! Implementing MB scaling 

					! scaling the relative velocities
					ipar = head1(i,j)			
					do while (ipar/=0)  
						vx_rel = vx(ipar)-vxcom(i,j)
						vy_rel = vy(ipar)-vycom(i,j)						! Relative velocities after collision
						vx(ipar)  = vxcom(i,j) + rel_vel_scale*(vx_rel) 			
						vy(ipar)  = vycom(i,j) + rel_vel_scale*(vy_rel) 			! Add scaled relative velocities
						ipar = list1(ipar)							! after collision using MBS scheme
					end do	
				end if	
			end if	
		END IF
	end do
end do	
IF (part_iter .lt. Gama*Ly*Lx) write(*,*) "Particles missing from domain"
end subroutine collision_rgs

!****************************************************************
! Subroutine for thermal boundary conditions for a cylinder
subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, g, np)
implicit none
integer, parameter :: np1=1000		! Expected no. of particles crossing boundary
real(kind=dp), dimension(:) :: rx, ry, rx1, ry1, vx, vy
real :: g, av, std, t_app, t_dec
integer :: Ly, i, np, j
logical :: Wall(np)
real(kind=dp) :: vxwall(np1), vywall(np1)

std   = sqrt(kbT)
j = 0
av =0.0

Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)
call random_normal(vxwall,np1,av,std)
call random_weibull(vywall,np1)

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
		IF (verlet == 2) THEN
			vx(i)  = vx(i) + g*(dt_c+t_app)/2
			rx1(i) = rx(i) + vx(i)*t_app
			vx(i)  = vxwall(j) 
			vx(i)  = vx(i) + g*t_dec/2
			rx1(i) = rx1(i) + vx(i)*t_dec
		ELSE IF (verlet == 1) THEN
			rx1(i) = rx(i) + vx(i)*t_app + (g*t_app**2)/2 
			vx(i)  = vxwall(j)
			rx1(i) = rx1(i) + vx(i)*t_dec + (g*t_dec**2)/2 
			vx(i)  = vx(i) + g*t_dec
		ELSE
			rx1(i) = rx(i) + vx(i)*t_app 
			vx(i)  = vxwall(j)
			rx1(i) = rx1(i) + vx(i)*t_dec 
			vx(i)  = vx(i) + g*dt_c
		END IF
	end if
	if (j==np1) then		! if more particles than expected are crossing, generate new numbers
		call random_normal(vxwall,np1,av,std)
		call random_weibull(vywall,np1)
		j=0	
	end if
end do

end subroutine thermal_wall

!***********************************************************************
!***********************************************************************
! Writing data to files for analysis
subroutine v_avg(vx,vy, rx, ry, head,list)
implicit none
INTEGER,PARAMETER :: grid_points = half_plane*Ly + 1
real(kind=dp), dimension(:) :: vx, vy, rx, ry
integer , dimension(:,:)    :: head
integer, dimension (:)      :: list
integer, save ::  t_count = 0, file_count = 0
integer :: i, j, k, ipar, p_count, out_unit, y_low, y_up
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), temp_com(Ly,Lx), density(Ly,Lx), weight, plane_pos
real(kind=dp), save :: vx1(Ly,Lx)=0.0, vy1(Ly,Lx)=0.0, forcing(2) = 0.0, tx(Ly,Lx) = 0.0, dens(Ly,Lx) = 0.0,&
vx_avg(grid_points)=0.0, vy_avg(grid_points)=0.0,  tot_part_count(grid_points)=0.0	! every grid points considered from 0 to Ly
logical :: Lexist
CHARACTER(LEN=200)  :: fname

vxcom   = 0.0
vycom   = 0.0
out_unit = 20
temp_com = 0.0
density   = 0.0

t_count = t_count+1
do i = 1,Ly
do j = 1,Lx    			  			
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
	if (p_count/=0) then
		vxcom(i,j) = vxcom(i,j)/p_count
		vycom(i,j) = vycom(i,j)/p_count
	end if	
	
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
dens = dens + density    			  	
!forcing = forcing + force
if (temperature) tx = tx + temp_com

if (modulo(t_count,ensemble_num)==0) then
	file_count = file_count+1
    	
	! Averaging all the velocities over the respective particles and time
	vx_avg = vx_avg/tot_part_count
	vy_avg = vy_avg/tot_part_count

	forcing = forcing/t_count
	tx = tx/t_count
	dens = dens/t_count
	write (fname, "(A<LEN(trim(file_name))>,A10)") trim(file_name),"_drag_lift"                            
	inquire(file = trim(data_path)//trim(fname)//'.dat', exist = Lexist)  	
	if (Lexist) then	
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',status="old",action="write",position="append")
    	else
		open(unit = out_unit, file=trim(data_path)//trim(fname)//'.dat', status = "new", action = "write")
	end if	
    	write (out_unit,*) forcing		  		
    	close(out_unit)

	write (fname, "(A<LEN(trim(file_name))>,A7,I0.2)") trim(file_name),"_vx_vy_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")

	DO i=1,half_plane*Ly+1
	write(out_unit,*) vx_avg(i),vy_avg(i)
	END DO
!	write (out_unit,"(<half_plane*Ly+1>F10.5)") vx_avg	!could also directly use the size(vx_avg) in format string
!	write (out_unit,"(<half_plane*Ly+1>F10.5)") vy_avg
	close(out_unit)	
	
	write (fname, "(A<LEN(trim(file_name))>,A5,I0.2)") trim(file_name),"_rho_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
    	
	do i = 1,Ly
		write (out_unit,"(<Lx>F10.5)") dens(i,:)
	end do
	close(out_unit)	

	if (temperature) then
		write (fname, "(A<LEN(trim(file_name))>,A6,I0.2)") trim(file_name),"_temp_",file_count                             
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
		do i = 1,Ly
			write (out_unit,"(<Lx>F10.5)") tx(i,:)
		end do
		close(out_unit)
	end if
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
subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, g, np)
implicit none
real(kind=dp), dimension(np)  :: rx, ry
real(kind=dp), dimension(np)  :: rx1, ry1, vx, vy
real :: g, t_app, t_dec
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
		vx(i) = vx(i)-g*dt_c
		IF (verlet == 2) THEN
			vx(i)  = vx(i) + g*(dt_c+t_app)/2
			rx1(i) = rx(i) + vx(i)*t_app
			vx(i)  = -vx(i) 
			vx(i)  = vx(i) + g*t_dec/2
			rx1(i) = rx1(i) + vx(i)*t_dec
		ELSE IF (verlet == 1) THEN
			rx1(i) = rx(i) + vx(i)*t_app + (g*t_app**2)/2
			vx(i)  = vx(i) + g*t_app
			vx(i)  = -vx(i)
			rx1(i) = rx1(i) + vx(i)*t_dec + (g*t_dec**2)/2 
			vx(i)  = vx(i) + g*t_dec
		ELSE
			rx1(i) = rx(i) + vx(i)*t_app 
			vx(i)  = -vx(i)
			rx1(i) = rx1(i) + vx(i)*t_dec 
			vx(i)  = vx(i) + g*dt_c
		END IF
	end if	
end do

end subroutine bounce_wall

!**********************************************************************************
! for writing file which gives the details of all parameters used in the simulation
SUBROUTINE param_file(tmax, t_avg, g, avg_interval)
implicit none
REAL  :: g
INTEGER :: tmax, t_avg, avg_interval

OPEN(UNIT = 10, FILE="Parameters.txt",ACTION = "write",STATUS="replace")
Write(10,*) "Simulation used in the algorithm: ", merge('MPC-AT','   SRD', MPC_AT == 1)
write(10,*)"Parameters used in the current simulation to which the data belongs"
write(10,*) "Particle density: ",Gama
write(10,*) "alpha rotation: ",alpha
write(10,*) "kbT: ",kbT
write(10,*) "Domain size Ly,Lx:",[Ly,Lx]
write(10,*)
write(10,*) "Time step used: ", dt_c
write(10,*) "Total time: ", tmax
write(10,*) "Total iterations: ",tmax/dt_c
write(10,*) "Force applied:",g
write(10,*)
write(10,*) "Streaming Algorithm:  ",merge(merge('Leapfrog','  Verlet',verlet/=1),'   Euler',verlet/=0)
write(10,*) "Periodicity in x: ",merge('Yes',' No',xy(1))
write(10,*) "Periodicity in y: ",merge('Yes',' No',xy(2))
write(10,*) "Random Grid Shift applied: ",merge('Yes',' No',random_grid_shift==1)
write(10,*) "Wall Boundary Condition: ",merge('Thermal wall',' Bounce wall',wall_thermal)
write(10,*)
write(10,*) "MBS applied: ",merge('Yes',' No',mb_scaling==1)
write(10,*) "MBS applied in number of iterations: ",freq
write(10,*) "Averaging starts from iteration: ",t_avg
write(10,*) "Interval between samples considered for average: ",avg_interval
write(10,*) "Number of samples considered for average : ",ensemble_num
write(10,*) "Mid plane of every cell used for averaging: ", merge('Yes',' No',half_plane == 2)
write(10,*)
write(10,*) "Analytical Viscosity as per given condition:  ",mu_tot 
close(10)

END SUBROUTINE param_file

end module SRD_library
