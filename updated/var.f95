module SRD_var
implicit none

INTEGER, PARAMETER :: dp = selected_real_kind(15, 307), Gama = 10, m=1, a0=1
REAL(kind=dp), PARAMETER :: pi=4.D0*DATAN(1.D0), e = 2.71828d0, aspect_ratio = 0.25
! Grid Parameters
INTEGER, PARAMETER :: Ly = 200, Lx = 2000, np = Ly*Lx*Gama, half_plane = 1
REAL(kind=dp), PARAMETER ::  alpha = pi/2.0d0, kbT = 1.0d0, dt_c = 0.1d0
! Forcing 
REAL(kind=dp) :: avg=0.0d0, std=sqrt(kbT/(m*1.0d0)), f_b = 4.0d-4
! time values
INTEGER :: tmax = 30000, t_avg = 100000, avg_interval=1, ensemble_num = 50000
! RGS, streaming
INTEGER :: random_grid_shift = 1, verlet = 1, grid_up_down, grid_check(0:Ly+1,Lx)=0 
! Thermostat
INTEGER :: mb_scaling = 0, MC_scaling = 1, mbs_freq=1
REAL(kind=dp) :: force(2), mu_tot, MC_strength = 0.25d0
LOGICAL :: xy(2)=[.TRUE., .FALSE.], temperature = .TRUE., wall_thermal = .FALSE., R_P = .FALSE.
! File naming 
LOGICAL :: image = .FALSE., dynamic_density = .FALSE.
 CHARACTER(len=100) :: file_name='cylinder', data_path='./' 
! For oscillating wall
INTEGER :: wall_oscillatory = 0, start_oscillation = 10000
REAL(kind=dp) :: f_wall=5.0d0, v_wall=0.2d0
! cylinder parameters:   obst_shape = 1 (for cylinder), 2 (for square)
INTEGER :: obst = 1, obst_shape = 2
REAL(kind=dp) :: rad = Ly*(aspect_ratio/2.0d0), xp = Lx/4.0d0, yp = Ly/2.0d0
REAL(kind=dp) :: obst_x = Lx/4.0d0, obst_y = Ly/2.0d0, obst_breadth = Ly*aspect_ratio, obst_length = 400 , obst_tfac = 10.0d0
REAL(kind=dp),ALLOCATABLE :: theta_intersect(:)   
LOGICAL, ALLOCATABLE ::  obst_par(:), obst_domain(:)
! Analysis of particle dynamics in equilibrium system, assuming periodic in both direction
LOGICAL :: MSD = .FALSE., STRUCTURE_FUNC = .FALSE.
INTEGER,PARAMETER :: tau = 10000, total_t = 100000, t_start = 200000
double complex,allocatable :: rho_k(:,:,:), S_k(:,:,:)
real,allocatable :: S_factor(:,:,:), MSD_xy(:,:,:), par_temp(:), MSD_value(:), par_periodic(:,:)
!! We should not define  very large static arrays typically of np. 
!! Use dynamic allocation instead for such arrays and keep stack size unlimited.
!! This precaution is important when using openmp and can be ignored without openmp.
end module SRD_var
