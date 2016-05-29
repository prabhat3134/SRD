program test_matlab
		
    		IMPLICIT NONE
    		REAL, PARAMETER ::  e = 2.71828, kbt=1.0
		INTEGER, PARAMETER :: gamma=5,m=1,a0=1  ! Particle parameters
		INTEGER, PARAMETER :: Lx = 10,Ly = 10, l0=50, b0=50, width=10 
		INTEGER, PARAMETER :: np=gamma*Ly*Lx, out_unit=20, dp1 = selected_real_kind(15, 307)
		real(kind=dp1) :: pi=4.D0*DATAN(1.D0)
		REAL :: av=0.0, std=sqrt(kbt/m), alpha, g=0.01, dt_c=1
		INTEGER :: N, tmax=1000, head(Ly,Lx), list(np),t, t_tot, xindex,i,j=0
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), vxwall(5000), vywall(5000)		
		logical :: l1(np)
		t_tot = tmax/dt_c
		alpha = pi/2.0
		
		interface
			subroutine reading (rx,ry,vx,vy, vxwall, vywall, np)	
			real, intent (inout) :: rx,ry,vx,vy, vxwall, vywall, np
			end subroutine swap
			subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy, vxwall, vywall, Ly, dt_c, g, np,j)
			real, intent (inout) :: (rx, ry, rx1, ry1, vx, vy, vxwall, vywall, Ly, dt_c, g, np,j)
			end subroutine swap
			subroutine partition(rx,ry,head,list,Lx,Ly,np)
			real, intent (inout) :: (rx,ry,head,list,Lx,Ly,np)
			end subroutine swap
			subroutine collision(vx, vy, head, list, Lx, Ly, alpha,np)
			real, intent (inout) :: (vx, vy, head, list, Lx, Ly, alpha,np)
			end subroutine swap			
		end interface	
		
		
		call reading(rx,ry,vx,vy, vxwall, vywall, np)		 		
		
		do t=1,t_tot
			vx  = vx + g*dt_c			
			rx1 = rx + vx*dt_c		! Streaming
    			ry1 = ry + vy*dt_c   			
    			call thermal_wall(rx, ry, rx1, ry1, vx, vy, vxwall, vywall, Ly, dt_c, g, np,j)
    			! For periodic boundary condition in x 			
    			
    			l1 = (rx1<=0.0)      			    			  			
    			rx1 = mod(rx1,Lx*1.0)
    			where (l1) 
    				rx1 = rx1+Lx
    			end where     			  			
    			rx = rx1					 			
    			ry = ry1 
    			call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells
    			  	
    			call collision(vx, vy, head, list, Lx, Ly, alpha,np)	! Collision  	    		
    			  			
		end do
		
		open (unit=out_unit,file="../Codes/test_MATLAB.dat",action="write",status="replace")
			do i=1,np
 				write(out_unit,'(F10.3, F10.3, F10.3, F10.3)') rx(i), ry(i), vx(i), vy(i)
 			end do
 		close(out_unit)	
		
		
end program test_matlab































subroutine thermal_wall(rx, ry, rx1, ry1, vx, vy, vxwall, vywall, Ly, dt_c, g, np,j)
implicit none
integer, parameter :: dp1 = selected_real_kind(15, 307)	! Expected no. of particles crossing boundary
real(kind=dp1), dimension(np) :: rx, ry, rx1, ry1, vx, vy
real :: g, dt_c, t_app, t_dec
integer :: Ly, i, np, j
logical :: Wall(size(ry))
real(kind=dp1) :: vxwall(5000), vywall(5000)


Wall = (ry1 > (Ly*1.0) .or. ry1 < 0.0)

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
		vx(i)  = vx(i)-g*(dt_c-t_app)
		rx1(i)  = rx(i) + vx(i)*t_app      
      	vx(i)  = vxwall(j) + g*t_dec
      	rx1(i)  = rx1(i) + vx(i)*t_dec       	      	       				
	end if	
end do

end subroutine thermal_wall

subroutine bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g,np)
implicit none
integer, parameter :: dp1 = selected_real_kind(15, 307)	
real(kind=dp1), dimension(np) :: rx, ry, rx1, ry1, vx, vy
real :: g, dt_c, t_app, t_dec
integer :: Ly, i, np
logical :: Wall(size(ry))
np = size(rx)

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

subroutine reading(rx,ry,vx,vy, vxwall, vywall, np)
implicit none
INTEGER, PARAMETER :: dp1 = selected_real_kind(15, 307)
integer :: np , i, out_unit=2
real(kind=dp1), dimension(np) :: rx,ry,vx,vy
real(kind=dp1), dimension(5000) :: vxwall, vywall


		open (unit=out_unit,file="rx_init.dat",action="read",status="old")
			do i=1,np
				read(out_unit,'(f15.11)'), rx(i)
			end do
			rx = rx*10
 		close(out_unit)
 		open (unit=out_unit,file="ry_init.dat",action="read",status="old")
			do i=1,np
				read(out_unit,'(f15.11)'), ry(i)
			end do
			ry = ry*10
 		close(out_unit)
 		open (unit=out_unit,file="vx_init.dat",action="read",status="old")
			do i=1,np
				read(out_unit,'(f15.12)'), vx(i)
			end do
 		close(out_unit)
 		open (unit=out_unit,file="vy_init.dat",action="read",status="old")
			do i=1,np
				read(out_unit,'(f15.12)'), vy(i)
			end do
 		close(out_unit)
 		open (unit=out_unit,file="vx_wall.dat",action="read",status="old")
			do i=1,5000
				read(out_unit,'(f15.12)'), vxwall(i)
			end do
 		close(out_unit)
 		open (unit=out_unit,file="vy_wall.dat",action="read",status="old")
			do i=1,5000
				read(out_unit,'(f15.12)'), vywall(i)
			end do
 		close(out_unit)
 		
end subroutine reading

subroutine partition(rx,ry,head,list,Lx,Ly,np)    	! assumed no grid shifting
implicit none
INTEGER, PARAMETER :: dp1 = selected_real_kind(15, 307)
real(kind=dp1), dimension(np) :: rx, ry
integer :: np, Lx, Ly, xindex, yindex,i, head(Ly,Lx), list(np)
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


subroutine collision(vx, vy, head, list, Lx, Ly, alpha,np)
implicit none
INTEGER, PARAMETER :: dp1 = selected_real_kind(15, 307)
real(kind=dp1), dimension(np) :: vx, vy
integer :: i,j,k=1,count1=0, ipar, jpar, Lx, Ly, head(Ly,Lx), list(np), np
real	:: r, alpha, vx_temp, vy_temp, num_rand(100)
real(kind=dp1) :: vxcom(Ly,Lx), vycom(Ly,Lx)

vxcom=0.0
vycom=0.0
 open (unit=20,file="coll_rand.dat",action="read",status="old")
		do i=1,100
			read(20,'(f15.9)'), num_rand(i)
		end do			
 close(20)
do j=1,Lx
	do i=1,Ly	
		count1 = 0	
		if (num_rand(k)>=0.5) then
			alpha=-alpha
		end if
		k = k +1		
		
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
			vx_temp  = vx(ipar)-vxcom(i,j) 
			vy_temp  = vy(ipar)-vycom(i,j) 
			if (i==3 .and. j==9) write(*,*) ipar, vx_temp, vy_temp, alpha, num_rand(k-1), k
			!vx(ipar) = vxcom(i,j) + cos(alpha)*vx_temp + sin(alpha)*vy_temp
			!vy(ipar) = vycom(i,j) - sin(alpha)*vx_temp + cos(alpha)*vy_temp
			vx(ipar) = vxcom(i,j) + sin(alpha)*vy_temp
			vy(ipar) = vycom(i,j) - sin(alpha)*vx_temp 
			ipar = list(ipar)		 
		end do					
	end do
end do

end subroutine collision


