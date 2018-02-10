module SRD_library
use omp_lib
implicit none
INTEGER, PARAMETER :: dp = selected_real_kind(15, 307), Gama = 10, m=1, a0=1
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828d0, aspect_ratio = 0.10d0
! Grid Parameters
INTEGER, PARAMETER :: Ly = 32, Lx = 32, np = Ly*Lx*Gama, half_plane = 1
REAL(kind=dp), PARAMETER :: alpha = pi/2.0d0, kbT = 1.0d0, dt_c = 1.0d0
! Forcing 
REAL(kind=dp) :: avg=0.0d0, std=sqrt(kbT/(m*1.0d0)), f_b = 0.0d-4
! time values
INTEGER :: tmax = 5e4, t_avg = 3e4, avg_interval=1, ensemble_num = 5e3
! RGS, streaming
INTEGER :: random_grid_shift = 0, verlet = 1, grid_up_down, grid_check(0:Ly+1,Lx)=0 
! Thermostat
INTEGER :: mb_scaling = 0, MC_scaling = 0, mbs_freq = 20
REAL(kind=dp) :: force(2), mu_tot, MC_strength = 0.25d0
LOGICAL :: xy(2)=[.TRUE., .TRUE.], temperature = .TRUE., wall_thermal = .FALSE.
LOGICAL :: R_P = .FALSE., slip = .TRUE.
REAL(kind=dp) :: RP_ratio = 3.0d0 
! Scrambling Boundary Condition for periodic X-Y condition
LOGICAL :: scrambling = .TRUE.
INTEGER :: scramble_exponent = 3
REAL(kind=dp) :: scramble_B = 5.0, scramble_U0 = 0.2d0
! File naming 
INTEGER :: wall_oscillatory = 0
LOGICAL :: image = .FALSE., dynamic_density = .FALSE. ,Init_unity_temp = .FALSE., write_poise_vel = .FALSE.
CHARACTER(len=100) :: file_name='equi', data_path='./' 
! cylinder parameters
! obst_shape = 1 (for cylinder), 2 (for square)
INTEGER :: obst = 0, obst_shape = 1
REAL(kind=dp) :: rad = Ly*(aspect_ratio/2.0d0), xp = Lx/4.0d0, yp = Ly/2.0d0
REAL(kind=dp) :: obst_x = Lx/4.0d0, obst_y = Ly/2.0d0, obst_breadth = Ly*aspect_ratio, obst_length = 400 
REAL(kind=dp),ALLOCATABLE :: theta_intersect(:)   
LOGICAL, ALLOCATABLE ::  obst_par(:)
! Analysis of particle dynamics in equilibrium system, assuming periodic in both direction
LOGICAL :: MSD = .FALSE., STRUCTURE_FUNC = .FALSE., trans_vel_SF = .FALSE.
INTEGER,PARAMETER :: tau = 1000, total_t = 50000, t_start = 50000
double complex,allocatable :: rho_k(:,:,:), S_k(:,:,:)
real,allocatable :: S_factor(:,:,:), MSD_xy(:,:,:), par_temp(:), MSD_value(:), par_periodic(:,:)
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
			idum = -mclock()
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
! For eliminating particles within a box
subroutine box_eliminate(x_dummy, y_dummy)
implicit none
integer :: i
real(kind=dp), dimension(:) :: x_dummy, y_dummy
logical, allocatable :: l1(:), l2(:)

ALLOCATE(l1(np), l2(np))

i = 1
if (obst_shape == 1) then
	l1 = (x_dummy - xp)**2 + (y_dummy - yp)**2 < rad**2
else if (obst_shape == 2) then 
	l1 = abs(x_dummy - obst_x) < obst_length/2.0d0 .and. abs(y_dummy - obst_y) < obst_breadth/2.0d0 
endif
l2 = .not. l1

do while (i <= size(x_dummy))
	if (l1(i)) then
		x_dummy(i) = ran()*Lx
		y_dummy(i) = ran()*Ly
		
		if ( (obst_shape == 1) .and. ((x_dummy(i) - xp)**2 + (y_dummy(i) - yp)**2 > rad**2)) then 
			i = i+1
		else if ( (obst_shape == 2) .and. (abs(x_dummy(i) - obst_x) > obst_length/2.0d0 .or. abs(y_dummy(i) - obst_y) > obst_breadth/2.0d0)) then
			i= i+1
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
!******************************************************
! MSD and Structure factor calculation
subroutine equilibrium(rx, ry, vx, vy)
implicit none
integer :: i,j,k, kx, ky
real(kind=dp) :: rx(:), ry(:), vx(:), vy(:), k1, k2
integer, save :: equi_index=0

equi_index = equi_index + 1
IF (equi_index .gt. total_t) return
IF (STRUCTURE_FUNC) THEN
	do kx = 0,(Lx/2)-1
	do ky = 0,(Ly/2)-1
		k1 = kx*2*pi/(Lx*1.0d0)
		k2 = ky*2*pi/(Ly*1.0d0)
		par_temp = k1*rx + k2*ry 
                IF(trans_vel_SF) THEN
                        rho_k(kx, ky, equi_index) = cmplx(sum( vy*cos(par_temp)),-sum( vy*sin(par_temp)))
                ELSE
                        rho_k(kx, ky, equi_index) = cmplx(sum(cos(par_temp)),-sum(sin(par_temp)))
                END IF
	end do
	end do
END IF
IF (MSD) THEN
	MSD_xy(:,1,equi_index) = rx + par_periodic(:,1)*Lx
	MSD_xy(:,2,equi_index) = ry + par_periodic(:,2)*Ly
END IF
end subroutine equilibrium
!*******************************************************
! Actual calculation of MSD and Structure factor
subroutine equi_cal()
implicit none
integer :: i,j,k
IF (STRUCTURE_FUNC) THEN
	!$OMP PARALLEL DO PRIVATE(i,j)
	DO i = 0,tau
		S_k(0:(Lx/2)-1, 0:(Ly/2)-1, i) = cmplx(0.0d0, 0.0d0)
		DO j = 1, total_t-i
			S_k(:,:,i) = S_k(:,:,i) + CONJG(rho_k(:,:,j+i))*rho_k(:,:,j)
		END DO
		S_k(:,:,i) = sum(S_k,3)/((total_t-i)*1.0d0)
		S_factor(:,:,i) = ABS( S_k(:,:,i) )/(np*1.0d0)
	END DO
	!$OMP END PARALLEL DO
	open(unit=20,file='DSF.dat',action='write',status = 'replace')
	DO j=0,tau
		write(20,"(3F15.5)") S_factor(1,0,j),  S_factor(0,1,j), S_factor(1,1,j)
	END DO
	close(20)
END IF
IF (MSD) THEN
	!$OMP PARALLEL DO PRIVATE(i,j,par_temp)
	DO i = 1,tau
		par_temp = 0.0d0
		DO j = 1,total_t-i
			par_temp = par_temp + (MSD_xy(:,1,j+i)-MSD_xy(:,1,j))**2 + (MSD_xy(:,2,j+i)-MSD_xy(:,2,j))**2
		END DO
		MSD_value(i) = SUM(par_temp)/((total_t-i)*np*1.0d0)
	END DO
	!$OMP END PARALLEL DO
	open(unit=20,file='MSD.dat',action='write',status = 'replace')
	DO j=1,tau
		write(20,"(I0,F10.5)") j, MSD_value(j) 
	END DO
	close(20)
END IF

end subroutine equi_cal

!******************************************
! Initialize the domain for a Riemann problem

subroutine Riemann_initial ( rx, ry ,vx, vy )
implicit none
REAL(kind=dp),dimension(:) :: rx, ry, vx, vy
REAL(kind=dp) :: xc, xr, temp
INTEGER :: i, k1, k2

xc = Lx/2.0d0
xr = xc/2.0d0
temp = sqrt( RP_ratio )
! for rho_ratio
k1 = CEILING( 0.5*np/(RP_ratio + 1.0d0))
k2 = k1 + CEILING( RP_ratio*np/(RP_ratio + 1))

! for uniform rho
!k1 = np/4 
!k2 = k1 + np/2 

DO i = 1,np
        rx(i) = ran()
END DO
rx(1:k1) = rx(1:k1) * xr
rx(k1+1: k2) = rx( k1+1: k2)* 2*xr + xr 
rx( k2+1:np) = rx(k2+1:np)*xr + xc + xr

!call random_normal(vx(1:k1), k1, avg, std )
!call random_normal(vx(k1+1: k2),(k2-k1) ,avg, temp )
!call random_normal(vx(k2+1: np),(np-k2) ,avg, std )

!call random_normal(vy(1:k1), k1, avg, std )
!call random_normal(vy(k1+1: k2),(k2-k1) ,avg, temp )
!call random_normal(vy(k2+1: np),(np-k2) ,avg, std )

end subroutine Riemann_initial

!******************************************
! write velocity data for the shock front for a Riemann problem
subroutine RP_shock_front( rx, ry, vx, vy, head, list, iter )
implicit none
INTEGER :: head(:,:), list(:), iter
INTEGER :: i,j,ipar, t1, t2
REAL(kind=dp), DIMENSION(:) :: rx, ry, vx, vy
REAL(kind=dp) :: U_shock = 1.70d0
CHARACTER(len=100) :: fname

t1 = 1500+ int(U_shock*iter) - 25
t2 = 1500+ int(U_shock*iter) + 15
write (fname, "(A12,I0)") "RP_vel_stat_",iter                             
open (unit=20, file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
DO j=t1,t2
DO i=1,Ly
        ipar = head(i,j)
        DO WHILE( ipar/=0)
                write(20, '(4F10.5)') rx(ipar), ry(ipar), vx(ipar), vy(ipar)
                ipar = list(ipar)
        END DO
END DO
END DO
close(20)
end subroutine RP_shock_front


!******************************************
! Intialize the domain
subroutine initialize(x_dummy, y_dummy, rx,ry,vx,vy, head, list)
implicit none
INTEGER ::  i, j, ipar, p_count, head(:,:), list(:)
REAL(kind=dp) ::  vp_max, vel_scale, temp
REAL(kind = dp) :: x_dummy(:), y_dummy(:), mu_kin, mu_col, RE
REAL(kind=dp) :: rx(:), ry(:), vx(:), vy(:), vxcom(Ly,Lx),vycom(Ly,Lx), Ek(Ly,Lx)

mu_kin = (gama*kbT*dt_c)*(gama/((gama - 1.0 + exp(-1.0*gama))*(1.0-cos(2.0*alpha)))- 0.5)
mu_col = ((1.0-cos(alpha))/(12.0*dt_c))*(gama - 1.0 + exp(-1.0*gama))
mu_tot = mu_kin + mu_col
vp_max = (Gama * Ly**2.0 *f_b)/(8.0*mu_tot)
RE = Gama*vp_max*merge(2*rad, Ly*1.0d0, obst ==1 ) /mu_tot
!write(*,*) vp_max, mu_tot, RE
!stop
do i=1, np
	x_dummy(i) = ran()*Lx
	y_dummy(i) = ran()*Ly
end do
IF (obst == 1) THEN 
	ALLOCATE(obst_par(np))
	ALLOCATE(theta_intersect(np))
	call box_eliminate(x_dummy, y_dummy)
	call gridcheck(grid_check)
END IF

rx = x_dummy
ry = y_dummy

call random_normal(vx,np,avg,std)
call random_normal(vy,np,avg,std)

IF (R_P) call Riemann_initial( rx, ry, vx, vy )

call partition(rx,ry,head,list)
vxcom = 0.0d0
vycom = 0.0d0
Ek = 0.0d0
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
		Ek(i,j)  = Ek(i,j) + 0.5*( vx(ipar)**2 + vy(ipar)**2 )
		ipar = list(ipar)
	end do
        IF (Init_unity_temp) THEN
                if (p_count .gt. 1) then
                        temp = Ek(i,j)/( (p_count - 1)*1.0d0)
                        vel_scale = sqrt(kbT/temp)
                        ipar = head(i,j)
                        do while (ipar/=0)
                                vx(ipar) = vx(ipar) * vel_scale
                                vy(ipar) = vy(ipar) * vel_scale
                                ipar = list(ipar)
                        end do
                end if
        END IF
END DO
END DO

IF (MSD) allocate( MSD_xy(np,2,0:total_t), par_periodic(np,2), par_temp(np), MSD_value(tau))
IF (STRUCTURE_FUNC) THEN
	allocate(rho_k(0:(Lx/2)-1,0:(Ly/2)-1, total_t), S_k(0:(Lx/2)-1,0:(Ly/2)-1, 0:tau))
	allocate(S_factor(0:(Lx/2)-1,0:(Ly/2)-1, 0:tau), par_temp(np))
END IF
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
IF (obst == 1) THEN
	IF (obst_shape==1) THEN
		obst_par  = (rx1 - xp)**2 + (ry1 - yp)**2 < rad**2
	ELSE
		obst_par = (abs(rx1 - obst_x) < obst_length/2.0d0 .and. abs(ry1 - obst_y) < obst_breadth/2.0d0)
	END IF
END IF

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
subroutine periodic_xy( rx1, ry1, vx, vy)
implicit none
REAL(kind = dp) :: rx1(:), ry1(:), vx(:), vy(:), dout
INTEGER :: i
!$OMP PARALLEL
if (xy(1)) then
    !$OMP DO
    DO i=1,np
        IF (rx1(i) < 0.0) THEN
            rx1(i) = rx1(i) + Lx
            IF (MSD) par_periodic(i,1) = par_periodic(i,1) - 1.0d0
        ELSE IF (rx1(i) > Lx*1.0) THEN
            rx1(i) = rx1(i) - Lx
            IF (MSD) par_periodic(i,1) = par_periodic(i,1) + 1.0d0
        END IF
    END DO
    !$OMP END DO
end if

if (xy(2)) then
    !$OMP DO
    DO i=1,np
        IF (ry1(i) < 0.0) THEN
            ry1(i) = ry1(i) + Ly
            IF (MSD) par_periodic(i,2) = par_periodic(i,2) - 1.0d0
        ELSE IF (ry1(i) > Ly*1.0) THEN
            ry1(i) = ry1(i) - Ly
            IF (MSD) par_periodic(i,2) = par_periodic(i,2) + 1.0d0
        END IF
    END DO
    !$OMP END DO
end if

if ( scrambling ) then
    !$OMP DO PRIVATE(dout)
	DO i=1,np
		IF ( rx1(i) <= scramble_B ) THEN
			IF ( ry1(i) >= scramble_B .AND. ry1(i) <= Ly-scramble_B) THEN
				!region 1 
                dout = rx1(i)
                call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 2)
			ELSE IF ( ry1(i) < scramble_B) THEN
				!region 8 
                dout = rx1(i)
				IF ( ry1(i) >= dout ) THEN
                    dout = rx1(i)
                    call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 2)
				ELSE 
					dout = ry1(i)
                    call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 1)
				END IF
			ELSE
			!region 2 

			ENDIF
		ELSE IF ( rx1(i) >= Lx-scramble_B ) THEN
			IF ( ry1(i) >= scramble_B .AND. ry1(i) <= Ly-scramble_B) THEN
			!region 5
				dout = Lx - rx1(i)
                call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 2)
			ELSE IF ( ry1(i) < scramble_B) THEN
			!region 6

			ELSE
			!region 4

			ENDIF

		ELSE
			IF ( ry1(i) >= Ly-scramble_B) THEN
			! region 3
				dout = Ly - ry1(i)
                call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 1)
			ELSE IF ( ry1(i) <= scramble_B) THEN
			! region 7
				dout = ry1(i)
                call scramble_pos_vel( rx1(i), ry1(i), vx(i), vy(i), dout, 1)
			ENDIF
		END IF
	END DO
    !$OMP END DO
end if
!$OMP END PARALLEL 
end subroutine periodic_xy
!********************************************************
subroutine scramble_pos_vel(rx, ry, vx, vy, dout, flag)
implicit none
INTEGER :: flag
REAL(kind=dp) :: rx, ry, vx, vy, dout
REAL(kind=dp) :: prob_disp, normal_vel(1)

prob_disp = (1.0 - dout/scramble_B)**scramble_exponent
IF ( prob_disp >= ran() ) THEN
    IF ( flag == 1 ) THEN   ! horizontal displacement
        rx = (Lx - 2.0*dout)*ran() + dout 
    ELSE
        ry = (Ly - 2.0*dout)*ran() + dout 
    ENDIF
    call random_normal( normal_vel, 1, scramble_U0, std) 
    vx = normal_vel(1)
    call random_normal( normal_vel, 1, avg, std) 
    vy = normal_vel(1)
END IF
end subroutine scramble_pos_vel

!********************************************************
! Grid check, gives which grids are overlapping with obstacles for further use in collision
subroutine gridcheck(grid_check)
implicit none
INTEGER :: grid_check(:,:), i, j, x, y
REAL(kind=dp) ::  xc, yc, l, b, R, rtmp

!grid_check(0:Ly+1,:) = 0	! This works for ifort compiler
grid_check(:,:) = 0		! This works for gfortran compiler
if (obst_shape==1) then	! for cylinder [center x, center y, radius]
	xc = xp
	yc = yp
	R  = rad
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
	xc = obst_x 
	yc = obst_y
	l  = obst_length
	b  = obst_breadth
	DO i = CEILING(yc- b/2.0d0)+1,CEILING(yc+b/2.0d0)	! i is in y direction
	DO j = CEILING(xc- l/2.0d0)+1,CEILING(xc+ l/2.0d0)	! j is in x direction
		grid_check(i,j) = 1
	END DO
	END DO
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

!head1(0:Ly+1,:) = 0		! This works for ifort
head1(:,:) = 0			! This kind of initialization works for gfortran (for negative indices, use indexing fully and not the default one)
list1(:) = 0
x_rand = ran()-0.5
y_rand = ran()-0.5
grid_up_down = merge(0,1,y_rand<0)

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
			Vx_temp  = vx(ipar)-vxcom(i,j) 
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
integer, save ::  t_count = 0, file_count = 0, steady_count=0
integer :: i, j, k, ipar, p_count, out_unit, density(Ly,Lx),den_t(4)
real(kind=dp) :: vxcom(Ly,Lx), vycom(Ly,Lx), temp_com(Ly,Lx),  weight, plane_pos, vxc, vyc
real(kind=dp), save :: vx1(Ly,Lx)=0.0d0, vy1(Ly,Lx)=0.0d0, forcing(2) = 0.0d0, tx(Ly,Lx) = 0.0d0, dens(Ly,Lx) = 0.0d0,&
vx_avg(grid_points)=0.0d0, vy_avg(grid_points)=0.0d0,  tot_part_count(grid_points)=0.0d0, converge(Ly,Lx,2)=0.0d0 ! every grid points considered from 0 to Ly
logical :: Lexist
CHARACTER(LEN=200)  :: fname, fmt_spec
integer,save :: interval = 0, den_count(4) = 0, png_write_count=0
real(kind=dp),save :: den_interval(Ly,Lx,4) = 0.0d0

vxcom   = 0.0
vycom   = 0.0
out_unit = 20
temp_com = 0.0
density   = 0
den_t = (/1,10,25,50/)
t_count = t_count+1
!$OMP PARALLEL DO PRIVATE(i,j,p_count,vxc,vyc,ipar,plane_pos,weight,k)
do j = 1,Lx    			  			
do i = 1,Ly
	p_count = 0
	ipar = head(i,j)					
	do while (ipar/=0)			
		vxcom(i,j) = vxcom(i,j) + vx(ipar)
		vycom(i,j) = vycom(i,j) + vy(ipar)
		p_count = p_count + 1
		IF ( write_poise_vel ) THEN
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
		END IF
		ipar = list(ipar)
	end do
	density(i,j) = p_count
	if(temperature .and. p_count > 1) then
		vxc = vxcom(i,j)/p_count
		vyc = vycom(i,j)/p_count
		ipar = head(i,j)    			  		
		do while (ipar/=0)
			temp_com(i,j) = temp_com(i,j) + 0.5*((vx(ipar) - vxc)**2 + (vy(ipar) - vyc)**2) 					
			ipar = list(ipar)		
		end do
		temp_com(i,j) = temp_com(i,j)/(p_count - 1)
	end if
end do
end do 
!$OMP END PARALLEL DO


IF (dynamic_density) THEN
	den_count(1) = den_count(1) + 1
	write (fname, "(A14,I0)") "Dyn_density_1_",den_count(1)                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
	write (fmt_spec,'(a,I0,a)') '( ', Lx ,'F10.5)'
	DO j=1,Ly
		write (out_unit, fmt = fmt_spec) density(j,:)*1.0d0 
	END DO
	close(out_unit)
        IF (temperature) THEN
                write (fname, "(A11,I0)") "Dyn_temp_1_",den_count(1)                             
                open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
                write (fmt_spec,'(a,I0,a)') '( ', Lx ,'F10.5)'
                DO j=1,Ly
                        write (out_unit, fmt = fmt_spec) temp_com(j,:)*1.0d0 
                END DO
                close(out_unit)
        END IF
END IF
IF (image)THEN
	DO i=1,4
		den_interval(:,:,i) = den_interval(:,:,i) + density
	END DO
	interval = interval + 1
	DO i=1,4
		IF (modulo( interval , den_t(i) )==0 ) THEN
			den_count(i) = den_count(i) + 1
			write (fname, "(A10,I0,A1,I0)") "shock_rho_",den_t(i),"_",den_count(i)                             
			open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
			 write (fmt_spec,'(a,I0,a)') '( ', Lx ,'F10.5)'
			DO j=1,Ly
				write (out_unit, fmt = fmt_spec ) den_interval(j,1:Lx,i)/(den_t(i)*1.0d0)
			END DO
			close(out_unit)
			write (fname, "(A10,I0,A1,I0)") "shock_col_",den_t(i),"_",den_count(i)                             
			open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
			DO j=1,Lx
				write (out_unit,"(4F10.5)") (den_interval(k,j,i)/(den_t(i)*1.0d0),k=20,80,20)
			END DO
			close(out_unit)
			den_interval(:,:,i) = 0.0d0
		END IF
	END DO
	IF (interval == den_t(4)) then
		interval = 0
		! call the script to plot the .png files & delete the data files
		write (fname, "(I0,A1,I0,A1,I0,A1,I0, A1,I0)")den_count(1)," ",den_count(2)," ",den_count(3)," ",den_count(4)," ",png_write_count                              
		call system('./png_script.sh '//trim(fname))
		png_write_count = png_write_count + 1 
	ENDIF
END IF

vx1 = vx1 + vxcom
vy1 = vy1 + vycom
dens = dens + density    			  	
forcing = forcing + force
if (temperature) tx = tx + temp_com

!IF (t_count == 1000) THEN
!	vxcom = vx1/dens
!	vycom = vy1/dens
!	IF (obst == 1) THEN
!		DO j=1,Lx
!		DO i=1,Ly
!			IF (isnan(vxcom(i,j))) vxcom(i,j) = 0.0d0
!			IF (isnan(vycom(i,j))) vycom(i,j) = 0.0d0
!		END DO
!		END DO
!	END IF
!	write(*,*) sum( abs( vxcom - converge(:,:,1) ) + abs( vycom - converge(:,:,2) ) )
!	converge(:,:,1) = vxcom
!	converge(:,:,2) = vycom
!END IF

if (modulo(t_count,ensemble_num)==0) then
	file_count = file_count+1
    	
	! Averaging all the velocities over the respective particles and time
	
	vx1 = vx1/dens
	vy1 = vy1/dens
    DO j=1,Lx
    DO i=1,Ly
            IF (isnan(vx1(i,j))) vx1(i,j) = 0.0d0
            IF (isnan(vy1(i,j))) vy1(i,j) = 0.0d0
    END DO
    END DO
	forcing = forcing/t_count
	tx = tx/t_count
	dens = dens/t_count
	write (fmt_spec,'(a,I0,a)') '(A', len_trim(file_name), ', A10)'
	write (fname, fmt = fmt_spec) trim(file_name),"_drag_lift"  
	inquire(file = trim(data_path)//trim(fname)//'.dat', exist = Lexist)  	
	if (Lexist) then	
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',status="old",action="write",position="append")
    	else
		open(unit = out_unit, file=trim(data_path)//trim(fname)//'.dat', status = "new", action = "write")
	end if	
    	write (out_unit,*) forcing		  		
    	close(out_unit)
	
	write (fmt_spec,'(a,I0,a)') '(A', len_trim(file_name), ', A7, I0)'
	write (fname, fmt = fmt_spec ) trim(file_name),"_vx_vy_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
	IF ( write_poise_vel ) THEN
        vx_avg = vx_avg/tot_part_count
        vy_avg = vy_avg/tot_part_count
		DO i=1,half_plane*Ly+1
			write(out_unit,*) (i*1.0-1)/(half_plane*1.0),vx_avg(i),vy_avg(i)
		END DO
	ELSE
		write (fmt_spec,'(a,I0,a)') '( ',Lx,'F10.5)'
		do i = 1,Ly
			write (out_unit, fmt = fmt_spec ) vx1(i,:)
		end do
		do i = 1,Ly
			write (out_unit,fmt = fmt_spec ) vy1(i,:)
		end do
	END IF
	close(out_unit)	
	
	write (fmt_spec,'(a,I0,a)') '(A', len_trim(file_name), ', A5, I0)'
	write (fname, fmt = fmt_spec ) trim(file_name),"_rho_",file_count                             
	open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
	write (fmt_spec,'(a,I0,a)') '( ',Lx,'F10.5)'
	do i = 1,Ly
		write (out_unit, fmt = fmt_spec ) dens(i,:)
	end do
	close(out_unit)	

	if (temperature) then
		write (fmt_spec,'(a,I0,a)') '(A', len_trim(file_name), ', A6, I0)'
		write (fname, fmt = fmt_spec ) trim(file_name),"_temp_",file_count                             
		open (unit=out_unit,file=trim(data_path)//trim(fname)//'.dat',action="write",status="replace")
		write (fmt_spec,'(a,I0,a)') '( ',Lx,'F10.5)'
		do i = 1,Ly
			write (out_unit,fmt = fmt_spec ) tx(i,:)
		end do
		close(out_unit)
	end if
	vx1 = 0.0d0
	vy1 = 0.0d0
	vx_avg = 0.0d0
	vy_avg = 0.0d0
	tot_part_count=0.0d0
	tx  = 0.0d0
	dens= 0.0d0
	t_count = 0	
	forcing = 0.0d0
end if    		

end subroutine v_avg
!***********************************************************************
! Subroutine for bounce back boundary condition at wall
subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy)
implicit none
real(kind=dp), dimension(:)  :: rx, ry, rx1, ry1, vx, vy
real(kind=dp) :: t_app, t_dec, vy_rel
real(kind=dp),save :: lower_wall_i=0.0d0, lower_wall_f=0.0d0, v_low_wall=0.0d0, f_wall=5.0d0, v_wall=0.05d0
integer :: i, wall_freq, iter=0

iter = iter + 1
wall_freq = f_wall/( v_wall * dt_c )	! half wall frequency 
if ( wall_oscillatory == 1 .and. iter .gt. 10000 ) then
        lower_wall_i = lower_wall_f
        lower_wall_f = lower_wall_i + v_wall*dt_c
        v_low_wall  = v_wall
        if (mod(iter, wall_freq) == 0 ) v_wall = -1.0d0*v_wall
        grid_check( 0:int(f_wall), : ) = 0
        grid_check( 0:ceiling(lower_wall_f), : ) = 1	! for collision
else
        lower_wall_f = 0.0d0
        lower_wall_i = 0.0d0
        v_low_wall = 0.0d0
        grid_check( 0,: ) = 1
        grid_check( Ly+1 ,: ) = 1
end if

!$OMP PARALLEL IF(np>100000)
!$OMP DO PRIVATE(t_app,t_dec,vy_rel) SCHEDULE(guided)
do i = 1,np
        if (ry1(i) > (Ly*1.0) .or. ry1(i) < lower_wall_f )   then
                if ( ry1(i) < lower_wall_f ) then
                        vy_rel = vy(i) - v_low_wall
                        if( vy_rel < 0 ) then
                                t_app  = (ry(i) - lower_wall_i) /abs(vy_rel)
                                t_dec  = dt_c - t_app
                                vy(i)  = v_low_wall - vy_rel
if (t_app < 0.0d0) write(*,*) iter, 'negative time-bounce', v_low_wall, vy(i),t_app, vy_rel, ry(i), lower_wall_i
                                ry1(i) = lower_wall_i + t_app*v_low_wall +  vy(i)*t_dec	
                        end if
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
                        vx(i)  = merge( vx(i), -vx(i), slip) 
                        vx(i)  = vx(i) + f_b*t_dec/2
                        rx1(i) = rx1(i) + vx(i)*t_dec
                ELSE IF (verlet == 1) THEN
                        rx1(i) = rx(i) + vx(i)*t_app + (f_b*t_app**2)/2
                        vx(i)  = vx(i) + f_b*t_app
                        vx(i)  = merge( vx(i), -vx(i), slip) 
                        rx1(i) = rx1(i) + vx(i)*t_dec + (f_b*t_dec**2)/2 
                        vx(i)  = vx(i) + f_b*t_dec
                ELSE
                        rx1(i) = rx(i) + vx(i)*t_app 
                        vx(i)  = merge( vx(i), -vx(i), slip) 
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
real(kind=dp) :: Vp,Re

Vp = (Gama * Ly**2.0 *f_b)/(8.0*mu_tot)
Re = (Gama * Vp * Ly )/mu_tot

OPEN(UNIT = 10, FILE="Parameters.txt",ACTION = "write",STATUS="replace")

write(10,*)"Parameters used in the current simulation to which the data belongs"
write(10,*) "Particle density: ",Gama, ", alpha rotation: ",alpha, ", kbT: ",kbT
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
write(10,*) "Maximum Poiseuille velocity (analytical): ", Vp
write(10,*) "Maximum Reynolds number: ", Re
write(10,*)
IF (obst == 1) write(10,*) "Cylinder centre position: ",xp,yp," and its radius: ",rad
IF (wall_oscillatory == 1) write(10,*) 'Oscillating lower wall'
IF (R_P) write(10,*) "Reimann Problem with density ratio: ",RP_ratio
close(10)

END SUBROUTINE param_file

end module SRD_library
