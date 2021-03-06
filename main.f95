PROGRAM SRD
	use SRD_library
	IMPLICIT NONE

	INTEGER :: iter, t_tot
	REAL(kind=dp):: start, finish
 	INTEGER, ALLOCATABLE :: head(:,:), list(:), head1(:,:), list1(:)
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:) ::rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy
	ALLOCATE(head(Ly,Lx), list(np), head1(0:Ly+1,Lx), list1(np))
	ALLOCATE(rx(np), ry(np), rx1(np), ry1(np), vx(np), vy(np), x_dummy(np), y_dummy(np))
	t_tot = tmax/dt_c
        call cpu_time(start)
	!$ start = omp_get_wtime()

	call initialize(x_dummy, y_dummy, rx,ry,vx,vy, head, list) 
	call param_file()
	do iter=1,t_tot
		call streaming(rx, ry, rx1, ry1, vx, vy)

		IF (obst==1) call par_in_cyl(rx, ry, rx1, ry1, vx, vy)
		
		if (.NOT. xy(2)) then 
			IF (wall_thermal) THEN
				call thermal_wall(rx, ry, rx1, ry1, vx, vy)
			ELSE
				call bounce_wall(rx, ry, rx1, ry1, vx, vy)
			END IF
		end if
		call periodic_xy( rx1, ry1)
		rx = rx1
		ry = ry1
		! Partitioning into cells
		if (random_grid_shift == 0) then  			
			call partition(rx1, ry1, head, list)
			call collision(vx, vy, head, list )! Collision without random grid shift      
		else 
			call partition_rgs(rx1, ry1, head1, list1)
			call collision_rgs(vx, vy, head1, list1)
		end if    					

		if (mod(iter,10000) == 0) then
                        call cpu_time(finish)
			!$ finish = omp_get_wtime()
			PRINT '("Iteration = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start
		endif

		! for averaging the values 
		if (iter .GT. t_avg .AND. MODULO(iter - t_avg, avg_interval) == 0) then
			if (random_grid_shift == 1) call partition(rx,ry,head,list)
			call v_avg(vx,vy, rx, ry, head,list) 
		end if
	end do
	call cpu_time(finish)
	!$ finish = omp_get_wtime()
	PRINT '("Total Elapsed CPU run time =",F10.2," seconds")', finish-start

	DEALLOCATE(head, list, head1, list1)
	DEALLOCATE(rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy)
	IF(obst == 1 ) THEN
		DEALLOCATE(obst_par)
		DEALLOCATE(theta_intersect)
	END IF
END PROGRAM SRD	
