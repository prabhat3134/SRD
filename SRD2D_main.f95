PROGRAM SRD
use SRD_library
#define progress(s) OPEN(UNIT = 10, FILE="progress.txt",ACTION = "write",STATUS="old",position="append"); write(10,*) s; CLOSE(10)
IMPLICIT NONE

INTEGER :: iter, iter_start=1, iter_tot, i, j
REAL(kind=dp):: start, finish, remaining
INTEGER, ALLOCATABLE :: head(:,:), list(:), head1(:,:), list1(:)
REAL(kind=dp), ALLOCATABLE, DIMENSION(:) ::rx, ry, rx1, ry1, vx, vy, x_dummy, y_dummy
ALLOCATE(head(Ly,Lx), list(np), head1(0:Ly+1,Lx), list1(np))
ALLOCATE(rx(np), ry(np), rx1(np), ry1(np), vx(np), vy(np), x_dummy(np), y_dummy(np))

iter_tot = tmax/dt_c
    call cpu_time(start)
!$ start = omp_get_wtime()

call initialize(x_dummy, y_dummy, rx,ry,vx,vy, head, list) 
call param_file()
IF ( restart == 1 ) call read_restartFile( rx, ry, vx, vy, iter_start )

do iter = iter_start, iter_tot
    call streaming(rx, ry, rx1, ry1, vx, vy)

    IF (obst==1) call par_in_cyl(rx, ry, rx1, ry1, vx, vy)
    
    if (.NOT. xy(2)) then 
        IF (wall_thermal) THEN
            call thermal_wall(rx, ry, rx1, ry1, vx, vy)
        ELSE
            call bounce_wall(rx, ry, rx1, ry1, vx, vy)
        END IF
    end if
    call periodic_xy( rx1, ry1, vx, vy )
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

!                if (iter .ge. 50 .and. R_P ) call RP_shock_front( rx, ry, vx, vy, head, list, iter )

    if (mod(iter,int(iter_tot/50)) == 0) then
        call cpu_time(finish)
        !$ finish = omp_get_wtime()
        remaining = ( finish - start)*( iter_tot - iter)/( ( iter - iter_start )*1.0 )
        write(progress_txt,'(a,I0,a,F12.5,a,F12.5,a)') "Iter=",iter,", CPU runtime=",finish-start," secs, Est. remain time:", remaining," secs."
        progress(trim( progress_txt ))
    endif

    ! for averaging the values 
    if (iter .GT. t_avg .AND. MODULO(iter - t_avg, avg_interval) == 0) then
        if (random_grid_shift == 1) call partition(rx,ry,head,list)
        call v_avg(vx,vy, rx, ry, head,list) 
    end if

    IF ( restart == 1 .AND. MODULO( iter, restart_iter ) == 0 ) call write_restartFile( rx, ry, vx, vy, iter )
    IF ( ( MSD .or. STRUCTURE_FUNC ) .and. iter .gt. t_start ) call equilibrium(rx,ry,vx,vy)
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
write(progress_txt,*) "Total Elapsed CPU run time =", finish-start," seconds"
progress(trim( progress_txt ))

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
