PROGRAM SRD
use SRD_library
IMPLICIT NONE

INTEGER :: iter, t_tot, i, j
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
!$ start = omp_get_wtime()
    call streaming(rx, ry, rx1, ry1, vx, vy)
!$ finish = omp_get_wtime()
PRINT '("Streaming = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start

!$ start = omp_get_wtime()
    IF (obst==1) call par_in_cyl(rx, ry, rx1, ry1, vx, vy)
!$ finish = omp_get_wtime()
PRINT '("Par in Cyl = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start
    
    if (.NOT. xy(2)) then 
        IF (wall_thermal) THEN
            call thermal_wall(rx, ry, rx1, ry1, vx, vy)
        ELSE
            call bounce_wall(rx, ry, rx1, ry1, vx, vy)
        END IF
    end if
!$ start = omp_get_wtime()
    call periodic_xy( rx1, ry1, vx, vy )
    rx = rx1
    ry = ry1
!$ finish = omp_get_wtime()
PRINT '("Periodic xy = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start

    ! Partitioning into cells
    if (random_grid_shift == 0) then  			
        call partition(rx1, ry1, head, list)
        call collision(vx, vy, head, list )! Collision without random grid shift      
    else 
!$ start = omp_get_wtime()
        call partition_rgs(rx1, ry1, head1, list1)
!$ finish = omp_get_wtime()
PRINT '("Partition rgs = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start
!$ start = omp_get_wtime()
        call collision_rgs(vx, vy, head1, list1)
!$ finish = omp_get_wtime()
PRINT '("collision rgs = ",I7," CPU run time =",F10.2," seconds")', iter, finish-start
    end if    					

!                if (iter .ge. 50 .and. R_P ) call RP_shock_front( rx, ry, vx, vy, head, list, iter )

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
    
    if (( MSD .or. STRUCTURE_FUNC ) .and. iter .gt. t_start) call equilibrium(rx,ry,vx,vy)
    IF (MSD .and. iter .eq. t_start) THEN
        par_periodic = 0.0d0
        MSD_xy(:,1,0) = rx
        MSD_xy(:,2,0) = ry
    END IF
end do
! To do equilibrium calculation
IF (MSD .or. STRUCTURE_FUNC) call equi_cal()

call cpu_time(finish)
!$ finish = omp_get_wtime()
PRINT '("Total Elapsed CPU run time =",F10.2," seconds")', finish-start

DEALLOCATE(head, list, head1, list1)
DEALLOCATE(rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy)
IF(obst == 1 ) THEN
    DEALLOCATE(obst_par)
    DEALLOCATE(theta_intersect)
END IF
IF (MSD) DEALLOCATE(MSD_xy, par_periodic, par_temp, MSD_value)
IF (STRUCTURE_FUNC) THEN
    DEALLOCATE(rho_k, S_k)
    DEALLOCATE(S_factor, par_temp)
END IF
END PROGRAM SRD	
