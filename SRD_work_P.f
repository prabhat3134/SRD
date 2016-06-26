      	PROGRAM SRD
		use SRD_lib_P
    		IMPLICIT NONE
		
		INTEGER, PARAMETER :: Gama=10,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 50,Ly = 50, l0=50, b0=50, width=10, msize = Lx*Ly
		INTEGER, PARAMETER :: np=Gama*Ly*Lx, out_unit=20,  out_unit1=21, out_unit2=22, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), g=1e-3, dt_c=1.0, start, finish, seed
		REAL (kind=dp1) :: alpha,fs = 7.0, twer
		INTEGER :: tmax=100000, head(Ly,Lx), list(np), t, t_tot, t_count=0, p_count, i,j, ipar, rand_for_cell
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), randposx(np), randposy(np), v1(np), v2(np), v3(np), v4(np)
		REAL(kind=dp1) :: temp(Ly,Lx), temp_com(Ly,Lx), tempy(Ly,Lx), vxcom(Ly), vycom(Ly), vx1(Ly), vy1(Ly), vx1_temp(Ly), vy1_temp(Ly)		
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
                           
    		open (unit=out_unit,file='./data/mbs_E_fluc_t_1e3.dat',action="write",status="replace")
    		open (unit=out_unit1,file='./data/mbs_E_av_t_1e3.dat',action="write",status="replace")
    			
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

			
			
				write(out_unit,*) sum(temp,2)/Lx
				write(out_unit1,*) sum(temp_com,2)/Lx
				!write(*,*) t
    				!
    			  			!write (out_unit,'(F9.6)') vy(i)	  		
			

			! To get temperature
			!if (mod(t,250)==0 .or. t==2099) then
			!	!write(*,*) "about to get screwed"
			!	write (fname, "(A14,I0.6)") "./data/position",t                             
    			 ! 		open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
			!	do i = 1,np
			!		write (out_unit2,'(F10.5, F10.5, F10.6, F10.6)') rx(i), ry(i), vx(i), vy(i)	
			!	end do
    			  !		do i = 1,Ly
    			  !			write (out_unit,*) (temp(i,j), j=1,Lx)	
    			  !			!write (out_unit,'(F9.6)') vx(i)	  		
    			  !		end do
    			 ! 		close(out_unit2)	
			! end if
    								
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
    			!if (t>100000) then 
    			!	t_count = t_count+1
    			!	vxcom   = 0.0
    			!	vycom   = 0.0
    			 !	do  i=1,Ly
    			  !		p_count = 0
    			  !		do j = 1,Lx    			  			
    			  !			ipar = head(i,j)    			  		
    			  !			do while (ipar/=0)
    			  !				vxcom(i) = vxcom(i) + vx(ipar)
    			  !				vycom(i) = vycom(i) + vy(ipar)
    			  !				p_count  = p_count +1
    			  !				ipar     = list(ipar)    			  		
    			  !			end do    			  			
    			  !		end do
    			  !		vxcom(i) = vxcom(i)/p_count
    			  !		vycom(i) = vycom(i)/p_count
    			  !	end do 
    			  !	vx1 = vx1 + vxcom
    			  !	vy1 = vy1 + vycom
    			 !	if (mod(t_count,2500)==0) then
    			  !		vx1_temp = vx1/t_count
    			  !		vy1_temp = vy1/t_count    			  		
    			  !		write (fname, "(A16,I0.6)") "./data/MBS_thermal_g1e4_vel_p",t                             
    			  !		open (unit=out_unit,file=trim(fname)//'.dat',action="write",status="replace")
    			  !		do i = 1,Ly
    			  !			write (out_unit,'(I2.2,F10.5, F10.5)') i, vx1_temp(i), vy1_temp(i)			  		
    			  !		end do
    			  !		close(out_unit)					
    			  !	end if    			  	  		  
    			!end if			
		end do
		close(out_unit)
		close(out_unit1)
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
    		

      END PROGRAM SRD	
