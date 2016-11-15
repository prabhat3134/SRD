      	PROGRAM SRD
		!use SRD_lib_T
		use SRD_library
		IMPLICIT NONE

		INTEGER, PARAMETER :: l0=50, b0=50, width=10, msize = Lx*Ly
		INTEGER, PARAMETER :: n_total=Gama*Ly*Lx, np = n_total,out_unit=20,  out_unit1=21, out_unit2=22, dp1 = selected_real_kind(15, 307)
		
		REAL :: av=0.0, std=sqrt(kbT/m), g=1e-5, start, finish, seed
		REAL (kind=dp1) :: fs = 7.0, twer, x_dummy(n_total), y_dummy(n_total)
		INTEGER :: tmax=1500000, head(Ly,Lx), list(n_total), head1(Ly+2,Lx)
		INTEGER :: list1(n_total),iter, t_tot, t_count=0, p_count, i,j, ipar, rand_for_cell, counter = 0 		
		INTEGER :: t_avg = 500000, avg_interval=1
		REAL(kind=dp1) :: temp(Ly,Lx), temp_com(Ly,Lx), tempy(Ly,Lx), vxcom(Ly), vycom(Ly), vx1(Ly), vy1(Ly), vx_temp(Ly), vy_temp(Ly)
		LOGICAL :: l1(np), l1_temp(np)
		CHARACTER(LEN=1024)  :: fname
		REAL(kind=dp1)::rx(np),ry(np),rx1(np),ry1(np),vx(np),vy(np), theta_intersect(np)		! needed when square block is present
		vxcom(Ly)=0.0
		vycom(Ly)=0.0
		vx1(Ly)=0.0
		vy1(Ly)=0.0
		
		t_tot = tmax/dt_c
		call cpu_time(start)
		
		call initialize(x_dummy, y_dummy, rx,ry,vx,vy, np, av, std, g)	 
		call param_file(tmax,t_avg,g, avg_interval)
		do iter=1,t_tot			
			call streaming(rx, ry, rx1, ry1, vx, vy, np, l1, g)

			IF (obst==1) call par_in_cyl(l1, rx, ry, rx1, ry1,vx, vy, theta_intersect, g)
			
			if (.NOT. xy(2)) then 
				IF (wall_thermal) THEN
					call thermal_wall(rx, ry, rx1, ry1, vx, vy, Ly, g, np)
				ELSE
					call bounce_wall(rx, ry, rx1, ry1, vx, vy, Ly, g, np)
				END IF
			end if
			call periodic_xy( rx1, ry1,np)
		
			rx = rx1
			ry = ry1

			! Partitioning into cells
			if (random_grid_shift == 0) then  			
				call partition(rx1,ry1,head,list,Lx,Ly,np)
				call collision(vx, vy, temp, temp_com, tempy, head, list, Lx, Ly )	! Collision without random grid shift      
			else 
				call partition_rgs(rx1,ry1,head1,list1,Lx,Ly,np)	
				call collision_rgs(vx, vy, head1, list1, Lx, Ly, np)
			end if    					

			if (mod(iter,10000) == 0) then
				call cpu_time(finish)
				write(*,*) "Time= ",iter,", run time= ",finish-start
			endif

    			! for averaging the values 
			if (MODULO(iter - t_avg, avg_interval) == 0) then
				if (random_grid_shift == 1) then
					call partition(rx,ry,head,list,Lx,Ly,np)
				end if
				call v_avg(vx,vy, rx, ry, head,list) 
 			end if
		end do


		call cpu_time(finish)
    		print '("Time = ",f10.3," seconds.")',finish-start
    		!DEALLOCATE(rx, ry, rx1, ry1, vx, vy,l1,theta_intersect)	
      END PROGRAM SRD	
