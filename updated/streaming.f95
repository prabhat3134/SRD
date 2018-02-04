module SRD_streaming
use SRD_var
use omp_lib
contains
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
 			IF (MSD) par_periodic(i,1) = par_periodic(i,1) - 1.0d0
                ELSE IF (rx1(i) > Lx*1.0) THEN
                        rx1(i) = rx1(i) - Lx
 			IF (MSD) par_periodic(i,1) = par_periodic(i,1) + 1.0d0
                END IF
        END DO
        !$OMP END PARALLEL DO
end if

if (xy(2)) then
	!$OMP PARALLEL DO
        DO i=1,np
                IF (ry1(i) < 0.0) THEN
                        ry1(i) = ry1(i) + Ly
 			IF (MSD) par_periodic(i,2) = par_periodic(i,2) - 1.0d0
                ELSE IF (ry1(i) > Ly*1.0) THEN
                        ry1(i) = ry1(i) - Ly
 			IF (MSD) par_periodic(i,2) = par_periodic(i,2) + 1.0d0
                END IF
        END DO
        !$OMP END PARALLEL DO
end if

end subroutine periodic_xy

  
!********************************************************************************************************
!******** Find all particles crossing the square and place it back at their respective position *********
!***************** Assumes Euler method for streaming: straight line trajectory *************************
subroutine par_in_box( rx, ry, rx1, ry1, vx, vy, dt_fine , p_force )
implicit none
INTEGER :: wall, i
REAL(kind=dp) :: rx, ry, rx1, ry1, vx, vy
REAL(kind=dp) :: dx, dy, w1, w2, w3, w4, slope, slope_c, p_force(:), dt_fine

w1 = obst_y - obst_breadth/2.0d0
w2 = obst_x + obst_length/2.0d0
w3 = obst_y + obst_breadth/2.0d0
w4 = obst_x - obst_length/2.0d0
dx = rx1  - rx 
dy = ry1  - ry 
if ( rx .gt. w4 .and. rx .lt. w2 ) then
        if ( ry .lt. w1 ) call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 1, p_force , dt_fine )
        if ( ry .gt. w3 ) call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 3, p_force , dt_fine )
        return
elseif ( ry  .gt. w1 .and. ry  .lt. w3 ) then
        if ( rx  .lt. w4 ) call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 4, p_force , dt_fine )
        if ( rx  .gt. w2 ) call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 2, p_force , dt_fine )
        return
else
        slope = atan2( dy, dx)
        if ( rx  .lt. w4 .and. ry  .lt. w1 ) then
                slope_c = atan2( w1-ry , w4- rx  )
                if (slope .gt. slope_c ) then
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 4, p_force , dt_fine )
                        return
                else
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 1, p_force , dt_fine )
                        return
                end if
        end if
        if ( rx  .gt. w2 .and. ry  .lt. w1 ) then
                slope_c = atan2( w1-ry , w2- rx  )
                if (slope .gt. slope_c ) then
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 1, p_force , dt_fine )
                        return
                else
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 2, p_force , dt_fine )
                        return
                end if
        end if
        if ( rx  .gt. w2 .and. ry  .gt. w3 ) then
                slope_c = atan2( w3-ry , w2- rx  )
                if (slope .gt. slope_c ) then
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 2, p_force , dt_fine )
                        return
                else
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 3, p_force , dt_fine )
                        return
                end if
        end if
        if ( rx  .lt. w4 .and. ry  .gt. w3 ) then
                slope_c = atan2( w3-ry , w4- rx  )
                if (slope .gt. slope_c ) then
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 3, p_force , dt_fine )
                        return
                else
                        call box_reflect( rx , ry , rx1 , ry1 , vx , vy , 4, p_force , dt_fine )
                        return
                end if
        end if
end if	
end subroutine par_in_box

!********************************************************************************************************
!******** Reflect all particles from the square and place it back at their respective position *********
!***************** Assumes the bounce wall condition and return forces on the box  **********************

subroutine box_reflect(rx, ry, rx1, ry1, vx, vy, wall, p_force, dt_fine ) 
implicit none
real(kind=dp) :: rx, rx1, ry, ry1, vx, vy, p_force(:)
real(kind=dp) :: t_app, t_dec, ywall, xwall, dt_fine
integer :: wall, vert

vert = mod(wall, 2)
vx = vx - f_b*dt_fine
   SELECT CASE (vert)
      CASE (1)          ! for horizontal wall, y = const
	IF (wall == 1) THEN     !for lower wall
		ywall= (obst_y - obst_breadth/2.0d0)
	ELSE
		ywall= (obst_y + obst_breadth/2.0d0)
	END IF
        t_app = abs( (ywall-ry)/vy ) 
	t_dec = dt_fine - t_app	
        if (t_app - dt_fine < 1e-4 ) then
                rx1 = rx
                ry1 = ry
                vy = -vy
                return
        end if
	if (t_dec .lt. 0) write(*,*) rx,rx1,ry,ry1,vx,vy,t_app,'- obstacle wall -ve time'	
	vy    = -vy
	ry1   = ywall + vy*t_dec
	rx1 = rx + vx*t_app 
	vx = vx + f_b*t_app
	vx = -vx
        p_force = [ -2*vx, -2*vy ]
	rx1 = rx1 + vx*t_dec
	vx  = vx + f_b*t_dec

      CASE (0 ) ! for vertical walls, x= const
	IF (wall == 2) THEN     ! for right wall
		xwall= (obst_x + obst_length/2.0d0)
	ELSE
		xwall= (obst_x - obst_length/2.0d0)
	END IF
        t_app = abs( (xwall-rx)/vx ) 
	t_dec = dt_fine - t_app	
        if (t_app - dt_fine < 1e-4 ) then
                rx1 = rx
                ry1 = ry
                vx = -vx
                return
        end if
	if (t_dec .lt. 0) write(*,*) rx,rx1,ry,ry1,vx,vy,t_app,'- obstacle wall -ve time'	
	ry1   = ry + vy*t_app
	vy    = -vy
	ry1   = ry1 + vy*t_dec
	vx = vx + f_b*t_app
	vx = -vx
        p_force = [ -2*vx, -2*vy ]
	rx1 = xwall + vx * t_dec 
	vx  = vx + f_b*t_dec

   END SELECT

end subroutine box_reflect

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
!
END MODULE SRD_streaming
