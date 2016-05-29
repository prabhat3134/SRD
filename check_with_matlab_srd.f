program test_matlab
		
    		use SRD_lib
    		IMPLICIT NONE
		INTEGER, PARAMETER :: gamma=5,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 10,Ly = 10, l0=50, b0=50, width=10 
		INTEGER, PARAMETER :: np=gamma*Ly*Lx, out_unit=20, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), alpha=3*pi/4, g=1e-4, dt_c=1.0, start, finish
		INTEGER :: tmax=209, head(Ly,Lx), list(np), t, t_tot, t_count=0, p_count, i,j, ipar
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np)	
		REAL(kind=dp1) ::	tempx(Ly,Lx), tempy(Ly,Lx), vxcom(Ly)=0.0, vycom(Ly)=0.0, vx1(Ly)=0.0, vy1(Ly)=0.0
		LOGICAL  :: l1(np)
		CHARACTER(LEN=1024)  :: fname
		!real, dimension(:) , allocatable::rx1,ry1		! needed when square block is present
		
		t_tot = tmax/dt_c
		
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
 		
 				 		
		
		do t=1,t_tot
			vx  = vx + g*dt_c			
			rx1 = rx + vx*dt_c		! Streaming
    			ry1 = ry + vy*dt_c  
    			
    				
    			call bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt,np)
    			! For periodic boundary condition in x 
    						
    			
    			l1 = (rx1<=0.0)      			    			  			
    			rx1 = mod(rx1,Lx*1.0)
    			where (l1) 
    				rx1 = rx1+Lx
    			end where     			  			
    			rx = rx1					 			
    			ry = ry1    			
    			
    			call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells
    			
    			if (t==260) then
    				do i=1,Ly
    					do j=1,Lx
    						write(*,*),'Particles in cell ',i,j
    						ipar = head(i,j)
    						do while (ipar/=0)
    							write(*,*), ipar
    							ipar = list(ipar)
    						end do
    					end do
    				end do
 			end if
    			
    			
    			!
 			!open (unit=out_unit,file="../Codes/check_list.txt",action="write",status="replace")
			!do i=1,np
 			!	write(out_unit,*) list(i)
 			!end do
    			  	if (t==tmax) write(*,*), rx(197),ry(197)
    			call collision(vx, vy, tempx, tempy, head, list, Lx, Ly, alpha, t)	! Collision  	    		
    			  if (t==tmax) write(*,*), rx(197),ry(197), vx(197)
    			 
		end do
		
			
		
		
end program test_matlab
