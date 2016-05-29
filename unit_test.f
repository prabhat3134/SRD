program unit_test
use SRD_lib
    		IMPLICIT NONE
		INTEGER, PARAMETER :: gamma=5,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 5,Ly = 5, l0=50, b0=50, width=10 
		INTEGER, PARAMETER :: np=gamma*Ly*Lx, out_unit=20, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), alpha=pi/2, g=0.1, dt_c=1
		INTEGER :: N, tmax=2, head(Ly,Lx), list(np),t, t_tot, xindex,i
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), tempx(Ly,Lx), tempy(Ly,Lx), vel
		REAL :: start, finish		
		logical :: l1(np)
		!real, dimension(:) , allocatable::rx1,ry1		! needed when square block is present
		
		t_tot = tmax/dt_c
		
		call cpu_time(start)
		call init_random_seed()
		call initialize(rx,ry,vx,vy, Ly, Lx, np, av, std)
		
		open (unit=out_unit,file="unit_pos.dat",action="write",status="replace")
			do i=1,np
 				write(out_unit,'F10.3, F10.3, F10.3, F10.3') rx(i), ry(i), vx(i), vy(i)
 			end do
 		close(out_unit)	 		
		
		do t=1,t_tot
			vx  = vx + g*dt_c			
			rx1 = rx + vx*dt_c		! Streaming
    			ry1 = ry + vy*dt_c   			
    			call thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt)
    			! For periodic boundary condition in x 			
    			
    			l1 = (rx1<=0.0)      			    			  			
    			rx1 = mod(rx1,Lx*1.0)
    			where (l1) 
    				rx1 = rx1+Lx
    			end where     			  			
    			rx = rx1					 			
    			ry = ry1 		
    			
    			
    			!call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells
    			
    			!write(*,*) 'List for particle at cell(3,2) starts'
    			!xindex = head(3,2)
    			!do while (xindex/=0)
    				!write(*,*) xindex
    				!xindex = list(xindex)
    			!end do    			
 			
    			!call collision(vx, vy, tempx, tempy, head, list, Lx, Ly, alpha)	! Collision  
    			  			
		end do
		
		
end program unit_test

