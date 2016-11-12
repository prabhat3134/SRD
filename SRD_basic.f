      	PROGRAM SRD
		use SRD_lib1
    		IMPLICIT NONE
		
		INTEGER, PARAMETER :: Gama=10,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 50,Ly = 50, l0=50, b0=50, width=10, msize = Lx*Ly
		INTEGER, PARAMETER :: np=Gama*Ly*Lx, out_unit=20,  out_unit1=21, out_unit2=22, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), g=1e-3, dt_c=1.0, start, finish, seed
		REAL (kind=dp1) :: alpha,fs = 7.0, twer
		INTEGER :: tmax=200000, head(Ly,Lx), list(np), t, t_tot, t_count=0, p_count, i,j, ipar, rand_for_cell
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), randposx(np), randposy(np), v1(np), v2(np), v3(np), v4(np)
		REAL(kind=dp1) :: temp(Ly,Lx), temp_com(Ly,Lx), tempy(Ly,Lx), vxcom(Ly), vycom(Ly), vx1(Ly), vy1(Ly), vx1_temp(Ly), vy1_temp(Ly)
		!REAL(kind=dp1) :: Tx(Ly,500000), 
		LOGICAL  :: l1(np)
		CHARACTER(LEN=1024)  :: fname
		!real, dimension(:) , allocatable::rx1,ry1		! needed when square block is present
		vxcom(Ly)=0.0
		vycom(Ly)=0.0
		vx1(Ly)=0.0
		vy1(Ly)=0.0
		t_tot = tmax/dt_c
		alpha = pi/2
		call cpu_time(start)
		
		call initialize(rx,ry,vx,vy, Ly, Lx, np, av, std)			                           
    		
		do t=1,t_tot
			
			rx1 = rx + vx*dt_c + (g*dt_c**2)/2	! Streaming
    			ry1 = ry + vy*dt_c  
    			vx  = vx + g*dt_c			  
    			
    			call thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt,np)
    			!call bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, dt_c, g, kbt,np)

    			! For periodic boundary condition in x 	
    			l1 = (rx1<=0.0)      			    			  			
    			rx1 = mod(rx1,Lx*1.0)
    			where (l1) 
    				rx1 = rx1+Lx
    			end where
    			     			
    			rx = rx1					 			
    			ry = ry1   					    			
    			
    			call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells    			
			
			call collision(vx, vy, temp, temp_com, tempy, head, list, Lx, Ly, alpha)	! Collision'			
				
		end do 		

		call cpu_time(finish)
    		print '("CPU Time of code = ",f10.3," seconds.")',finish-start    		

      END PROGRAM SRD	
