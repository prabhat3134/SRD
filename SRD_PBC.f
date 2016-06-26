	PROGRAM SRD_PBC      
		use SRD_lib_T
    		IMPLICIT NONE
		
		INTEGER, PARAMETER :: Gama=10,m=1,a0=1  ! Particle parameters
		REAL, PARAMETER  :: kbt=1.0
		INTEGER, PARAMETER :: Lx = 50,Ly = 50, l0=50, b0=50, width=10, msize = Lx*Ly
		INTEGER, PARAMETER :: np=Gama*Ly*Lx, out_unit=20,  out_unit1=21, out_unit2=22, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbt/m), g=1e-4, dt_c=1.0, start, finish, seed
		REAL (kind=dp1) :: alpha,fs = 7.0, twer
		INTEGER :: tmax=90000, head(Ly,Lx), list(np), t, t_tot, t_count=0, p_count, i,j, ipar, rand_for_cell, counter = 0, fc_i, fc_j
		REAL(kind=dp1) :: rx(np), ry(np), vx(np), vy(np), rx1(np), ry1(np), randposx(np), randposy(np), v1(np), v2(np), v3(np), v4(np)
		REAL(kind=dp1) :: temp(Ly,Lx), temp_com(Ly,Lx), tempy(Ly,Lx), vxcom(Ly), vycom(Ly), vx1(Ly), vy1(Ly), vx1_temp(Ly), vy1_temp(Ly)
		REAL(kind=dp1) :: momx_before_col(Ly,Lx), momx_after_col(Ly,Lx), momy_before_col(Ly,Lx), momy_after_col(Ly,Lx)  			
		REAL(kind=dp1) :: ke_before_col(Ly,Lx), ke_after_col(Ly,Lx)
		REAL(kind=dp1) :: momx_diff(Ly,Lx), momy_diff(Ly,Lx), ke_diff(Ly,Lx)
		!REAL(kind=dp1) :: Tx(Ly,500000), 
		LOGICAL  :: l1(np)
		CHARACTER(LEN=1024)  :: fname
		!real, dimension(:) , allocatable::rx1,ry1		! needed when square block is present
		vxcom(Ly)=0.0
		vycom(Ly)=0.0
		vx1(Ly)=0.0
		vy1(Ly)=0.0
		ke_before_col = 0.0
		ke_after_col = 0.0
		momx_before_col = 0.0
		momy_after_col = 0.0
		momx_before_col = 0.0
		momy_after_col = 0.0
		
		t_tot = tmax/dt_c
		alpha = pi/2
		call cpu_time(start)
		
		call initialize(rx,ry,vx,vy, Ly, Lx, np, av, std)			                           
    		
		do t=1,t_tot
			
			rx1 = rx + vx*dt_c + (g*dt_c**2)/2	! Streaming
    			ry1 = ry + vy*dt_c  
    			vx  = vx + g*dt_c			  

    			! Periodic Boundary condition in x
    			l1 = (rx1<=0.0)      
    			rx1 = mod(rx1,Lx*1.0)
    			where (l1) 
    				rx1 = rx1+Lx
    			end where
    			    
			! Periodic Boundary condition in y		
    			l1 = (ry1<=0.0)      
    			ry1 = mod(ry1,Ly*1.0)
    			where (l1) 
    				ry1 = ry1+Ly
    			end where
    			     	
			! Make final changes in the position 		
    			rx = rx1					 			
    			ry = ry1   			

			! Partition the particles into cells
    			call partition(rx,ry,head,list,Lx,Ly,np)	! partition into cells    			
			
			momx_before_col = 0.0
			momy_before_col = 0.0
			momx_after_col = 0.0
			momy_after_col = 0.0
			ke_before_col = 0.0
			ke_after_col = 0.0

			
			! Check for momentum and KE before collision for every cell
			do i=1,Ly
				do j=1,Lx
					ipar = head(i,j)
					do while (ipar/=0)
						momx_before_col(i,j) = momx_before_col(i,j) + vx(ipar) 
						momy_before_col(i,j) = momy_before_col(i,j) + vy(ipar)
						ke_before_col(i,j) = ke_before_col(i,j) + vx(ipar)**2 + vy(ipar)**2
						ipar = list(ipar)
					end do
				end do
			end do			


			! Collision of the particles in the cell
			call collision(vx, vy, temp, temp_com, tempy, head, list, Lx, Ly, alpha)	! Collision'	


			! Check for momentum and KE after collision for every cell
			do i=1,Ly
				do j=1,Lx
					ipar = head(i,j)
					do while (ipar/=0)
						momx_after_col(i,j) = momx_after_col(i,j) + vx(ipar) 
						momy_after_col(i,j) = momy_after_col(i,j) + vy(ipar)
						ke_after_col(i,j) = ke_after_col(i,j) + vx(ipar)**2 + vy(ipar)**2
						ipar = list(ipar)
					end do
				end do
			end do		
				
		end do

		! Difference in KE and x and y momentum before and after collision
			if (mod(t,1000)==0) then
				
				momx_diff = momx_after_col - momx_before_col
				momy_diff = momy_after_col - momy_before_col
				ke_diff = ke_after_col - ke_before_col
																
				counter = counter + 1;
			
				write (fname, "(A14,I0.6)") "./data/momx_diff_t",t                             
    				open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (momx_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    			  	close(out_unit2)	
					
				write (fname, "(A14,I0.6)") "./data/momy_diff_t",t                             
    			  	open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (momy_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    				close(out_unit2)	
					
				write (fname, "(A14,I0.6)") "./data/ke_diff_t",t                             
    				open (unit=out_unit2,file=trim(fname)//'.dat',action="write",status="replace")
				do fc_i = 1,Ly
					write (out_unit2,*) (ke_diff(fc_i,fc_j), fc_j=1,Lx)	
				end do
    				close(out_unit2)	
					
			end if   		

		call cpu_time(finish)
    		print '("CPU Time of code = ",f10.3," seconds.")',finish-start    		

      END PROGRAM SRD_PBC	


		

			

			!write(*,*) t

			
