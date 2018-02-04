PROGRAM SRD
	use SRD_library
	IMPLICIT NONE

	INTEGER :: iter, t_tot, kx, ky, par_iter, tau, i, j
	REAL(kind=dp):: start, finish, k1, k2
 	INTEGER, ALLOCATABLE :: head(:,:), list(:), head1(:,:), list1(:)
	REAL(kind=dp), ALLOCATABLE, DIMENSION(:) ::rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy, par_temp
double complex,allocatable :: rho_k(:,:,:), S_k(:,:,:)
real,allocatable :: S_factor(:,:,:)
	ALLOCATE(head(Ly,Lx), list(np), head1(0:Ly+1,Lx), list1(np))
	ALLOCATE(rx(np), ry(np), rx1(np), ry1(np), vx(np), vy(np), x_dummy(np), y_dummy(np), par_temp(np))
allocate(rho_k(0:24,0:24,30000), S_k(0:24,0:24,5000))
allocate(S_factor(0:24,0:24,5000))
par_iter = 0
tau = 5000
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
		
		if (iter .gt. 200000) then
			par_iter = par_iter + 1
			do kx = 0,24
			do ky = 0,24
				k1 = kx*2*pi/Lx
				k2 = ky*2*pi/Ly
				par_temp = kx*rx + ky*ry 
				rho_k(kx,ky,par_iter) = cmplx(sum(cos(par_temp)),-sum(sin(par_temp)))
			end do
			end do
		end if
	end do

DO i = 1,tau
S_k(0:24,0:24,i) = cmplx(0.0d0, 0.0d0)
DO j = 1,100000-i
	S_k(:,:,i) = S_k(:,:,i) + CONJG(rho_k(:,:,j+i))*rho_k(:,:,j)
END DO
	S_k(:,:,i) = sum(S_k,3)/((20000-i)*1.0d0)
	S_factor(:,:,i) = abs(S_k(:,:,i))/(np*1.0d0)
END DO
open(unit=20,file='DSF.dat',action='write',status = 'replace')
DO j=1,tau
write(20,"(I0,25F10.5)") j,(S_factor(i,i,j), i=0,24) 
END DO
close(20)

	call cpu_time(finish)
	!$ finish = omp_get_wtime()
	PRINT '("Total Elapsed CPU run time =",F10.2," seconds")', finish-start

	DEALLOCATE(head, list, head1, list1)
	DEALLOCATE(rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy)
	IF(obst == 1 ) THEN
		PRINT *, "Deallocation of arrays in obst case"
		DEALLOCATE(obst_par)
		DEALLOCATE(theta_intersect)
	END IF
END PROGRAM SRD	
