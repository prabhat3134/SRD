      	PROGRAM SRD
		use SRD_lib1
    		IMPLICIT NONE
		
		INTEGER, PARAMETER :: Gama=10,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 64,Ly = 64, l0=50, b0=50, width=10, msize = Lx*Ly
		INTEGER, PARAMETER :: np=Gama*Ly*Lx, out_unit=20, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), alpha=pi/2, g=1e-4, dt_c=1.0, start, finish, seed
		INTEGER :: tmax=200000, head(Ly,Lx), list(np), t, t_tot, t_count=0, p_count, i,j, ipar, rand_for_cell
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), randposx(np), randposy(np), v1(np), v2(np), v3(np), v4(np)	
		REAL(kind=dp1) :: temp(Ly,Lx), tempx(Ly,Lx), tempy(Ly,Lx), vxcom(Ly), vycom(Ly), vx1(Ly), vy1(Ly), randcol(msize), vx1_temp(Ly), vy1_temp(Ly)
		!REAL(kind=dp1) :: Tx(Ly,500000), 
		LOGICAL  :: l1(np)
		CHARACTER(LEN=1024)  :: fname
		!real, dimension(:) , allocatable::rx1,ry1		! needed when square block is present
		vxcom(Ly)=0.0
		vycom(Ly)=0.0
		vx1(Ly)=0.0
		vy1(Ly)=0.0
		t_tot = tmax/dt_c
		
		call cpu_time(start)
		
		call initialize(rx,ry,vx,vy, Ly, Lx, np, av, std)	

		!write (fname, "(A14,I0.6)") "../Codes/initial_stats"                             
    			!open (unit=out_unit,file='../Codes/random_stats_testing.dat',action="write",status="replace")
    			!do i = 1,np
    				!write (out_unit,'(F10.5, F10.5, F10.6, F10.6)') rx(i), ry(i), vx(i), vy(i)	
    			  			!write (out_unit,'(F9.6)') vy(i)	  		
    			  		!end do
    			  		!		


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
    			!l1 = (ry1<=0.0)      			    			  			
    			!ry1 = mod(ry1,Ly*1.0)
    			!where (l1) 
    			!	ry1 = ry1+Lx
    			!end where      			
    			rx = rx1					 			
    			ry = ry1   					
    			
    			
    			call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells
    			
			
			call collision(vx, vy, temp, tempx, tempy, head, list, Lx, Ly, alpha)	! Collision

			! To get temperature
			if (mod(t,1000)==0) then
				!write(*,*) "about to get screwed"
				write (fname, "(A14,I0.6)") "../Codes/Tcell",t                             
    			  		open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    			  		do i = 1,Ly
    			  			write (out_unit,*) (temp(i,j), j=1,Lx)	
    			  			!write (out_unit,'(F9.6)') vx(i)	  		
    			  		end do
    			  		close(out_unit)	
			end if
    								
    			!if (mod(t,tmax)==0) then
    						
			!		write (fname, "(A21,I0.6)") "../Codes/thermal_wall",t                             
    			 ! 		open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    			 ! 		do i = 1,np
    			 ! 			write (out_unit,'(F10.5, F10.5, F10.6, F10.6)') rx(i), ry(i), vx(i), vy(i)	
    			 ! 			!write (out_unit,'(F9.6)') vx(i)	  		
    			 ! 		end do
    			 ! 		close(out_unit)	
    			!end if
    			! for averaging the values 
    			!if (t>200000) then 
    			!	t_count = t_count+1
    			!	vxcom   = 0.0
    			!	vycom   = 0.0
    			! 	do  i=1,Ly
    			!  		p_count = 0
    			!  		do j = 1,Lx    			  			
    			!  			ipar = head(i,j)    			  		
    			!  			do while (ipar/=0)
    			!  				vxcom(i) = vxcom(i) + vx(ipar)
    			!  				vycom(i) = vycom(i) + vy(ipar)
    			!  				p_count  = p_count +1
    			!  				ipar     = list(ipar)    			  		
    			!  			end do    			  			
    			!  		end do
    			!  		vxcom(i) = vxcom(i)/p_count
    			!  		vycom(i) = vycom(i)/p_count
    			!  	end do 
    			!  	vx1 = vx1 + vxcom
    			!  	vy1 = vy1 + vycom
    			  	!if (mod(t_count,2500)==0) then
    			  	!	vx1_temp = vx1/t_count
    			  	!	vy1_temp = vy1/t_count    			  		
    			  	!	write (fname, "(A16,I0.6)") "../Codes/vel_poise",t                             
    			  	!	open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    			  	!	do i = 1,Ly
    			  	!		write (out_unit,'(I2.2,F10.5, F10.5)') i, vx1_temp(i), vy1_temp(i)			  		
    			  	!	end do
    			  	!	close(out_unit)					
    			  	!end if    			  	  		  
    			!end if			
		end do
		!close(out_unit)
	! For averaged velocity profile
  		!vx1_temp = vx1/t_count
    	        !vy1_temp = vy1/t_count    			  		
 		!write (fname, "(A26,I0.6)") "../Codes/thermal_vel_poise",t                             
    		!open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    		!do i = 1,Ly
    		!	write (out_unit,'(I2.2,F10.5, F10.5)') i, vx1_temp(i), vy1_temp(i)			  		
    		!end do
    		!close(out_unit)	
 		

		call cpu_time(finish)
    		print '("Time = ",f10.3," seconds.")',finish-start
    		
		!call box_eliminate(rx,ry,rx1,ry1,l0,b0,width)		! for eliminating particle inside the box
		!N = size(rx1)			! Number of particle
		!write(*,*) np,N
		!open (unit=out_unit,file="rx.dat",action="write",status="replace")
 		!write (out_unit,*) rx
 		!open (unit=out_unit,file="ry.dat",action="write",status="replace")
 		!write (out_unit,*) ry
 		!open (unit=out_unit,file="vx.dat",action="write",status="replace")
 		!write (out_unit,*) vx
 		!open (unit=out_unit,file="vy.dat",action="write",status="replace")
 		!write (out_unit,*) vy
 		
 		!call random_weibull(vel,np,kbt)
		!open (unit=out_unit,file="vel.dat",action="write",status="replace")
 		!write (out_unit,*) vel

      END PROGRAM SRD	
